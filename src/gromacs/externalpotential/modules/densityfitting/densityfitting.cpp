/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#include "gmxpre.h"

#include "densityfitting.h"

#include <cstdio>


#include <algorithm>
#include <numeric>
#include <ios>

#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/externalpotential/forceplotter.h"
#include "gromacs/externalpotential/externalpotentialIO.h"
#include "gromacs/externalpotential/forceplotter.h"
#include "gromacs/math/quaternion.h"

#include "gromacs/externalpotential/atomgroups/wholemoleculegroup.h"
#include "gromacs/externalpotential/mpi-helper.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/fileio/volumedataio.h"
#include "gromacs/fileio/json.h"
#include "gromacs/math/volumedata/gausstransform.h"
#include "gromacs/math/volumedata/gridreal.h"
#include "gromacs/math/volumedata/gridmeasures.h"
#include "gromacs/fileio/xtcio.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/atoms.h"
#include "emscatteringfactors.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/utility/gmxomp.h"



namespace gmx
{

class volumedata::FastGaussianGriddingForce : public FastGaussianGridding
{
    public:
        using FastGaussianGridding::FastGaussianGridding;
        RVec force(const rvec x, const real weight,  const volumedata::GridReal &div);

};

RVec
volumedata::FastGaussianGriddingForce::force(const rvec x, const real weight, const volumedata::GridReal &densityDerivative)
{
    prepare_2d_grid(x, weight);
    ivec l_min;
    ivec l_max;
    RVec force = {0, 0, 0};

    for (size_t i = XX; i <= ZZ; ++i)
    {
        l_min[i] = grid_index_of_spread_atom_[i] - m_spread_ < 0 ? 0 : grid_index_of_spread_atom_[i]-m_spread_;
        l_max[i] = grid_index_of_spread_atom_[i] + m_spread_ >= grid_->extend()[i] ? grid_->extend()[i]-1 : grid_index_of_spread_atom_[i]+m_spread_;
    }

    std::vector < std::vector < int>> ceilSqrtLUT(m_spread_+1);
    for (int d_z = 0; d_z <= m_spread_; d_z++)
    {
        ceilSqrtLUT[d_z].resize( (int)std::ceil(sqrt(m_spread_*m_spread_ - d_z*d_z))+1);
        for (int d_y = 0; d_y*d_y <= m_spread_*m_spread_ - d_z*d_z; d_y++)
        {
            ceilSqrtLUT[d_z][d_y] = (int)std::ceil(sqrt(m_spread_*m_spread_ - d_z*d_z- d_y*d_y));
        }
    }
    real   spread_zy;
    std::vector<real>::iterator densityDerivativeIterator;
    real * spread_1d_XX;


    RVec     d_x       = densityDerivative.unit_cell_XX();
    RVec     d_y       = densityDerivative.unit_cell_YY();
    RVec     d_z       = densityDerivative.unit_cell_ZZ();
    RVec     shift_z;                                          //< (x-v_0) + dz
    RVec     shift_yz;                                         //< ((x-v_0) + dz) + dy
    RVec     shift_xyz;                                        //< (((x-v_0) + dz) + dy) + dx



    rvec voxel_force; //< The force each voxel contributes

    /*
     * To avoid repeated calculation of the distance from atom to voxel, (x-v)
     * calculate instead for voxel v at i_x,i_y,i_z: (((x-v_0) + dz+...(i_z times)...+dz) + dy+...+dy) + dx...+dx
     */

    /*
     * Use iterators to step though the grid.
     * This relies on the assumption that the grid is stored with z the slowest and x the fasted changing dimension with no padding
     * (x,y,z not being linked to any coordiante system, but short-hand for first, second, third dimension)
     * Loosing generality through this approach, we save substantial time when we don't have to calculate the grid index.
     */

    RVec d_z_offset;
    svmul(l_min[ZZ], d_z, d_z_offset);
    RVec d_y_offset;
    svmul(l_min[YY], d_y, d_y_offset);
    RVec d_x_offset;
    svmul(l_min[XX], d_x, d_x_offset);

    rvec_sub(densityDerivative.gridpoint_coordinate({0, 0, 0}), x, shift_z); // using the grid-origin as reference (v_0-x)
    rvec_inc(shift_z, d_z_offset);                                           // move in "z-slices" until spreading starts


    for (int l_grid_z = l_min[ZZ]; l_grid_z <= l_max[ZZ]; ++l_grid_z)
    {
        shift_yz = shift_z;                                              // start at the beginning of a "y-row"
        rvec_inc(shift_yz, d_y_offset);                                  // move into the "y-row" to the voxel where spreading starts
        auto densityDerivativeGridData   = densityDerivative.access();
        int  d_l_z                       = l_grid_z - grid_index_of_spread_atom_[ZZ];
        int  l_z                         = d_l_z + m_spread_; // the grid index in the 2d spread grid around the atom
        int  globalGridIndexYYStart      = std::max(l_min[YY], grid_index_of_spread_atom_[YY] - ceilSqrtLUT[std::abs(d_l_z)][0]);
        int  globalGridIndexYYEnd        = std::min(l_max[YY], grid_index_of_spread_atom_[YY] + ceilSqrtLUT[std::abs(d_l_z)][0]);
        for (int l_grid_y = globalGridIndexYYStart; l_grid_y <= globalGridIndexYYEnd; ++l_grid_y)
        {
            shift_xyz = shift_yz;            // start at the beginning of an "x-column"
            rvec_inc(shift_xyz, d_x_offset); // move into the "x-row" to the voxel where spreading starts

            int d_l_y        =  l_grid_y - grid_index_of_spread_atom_[YY];
            int l_y          = d_l_y+m_spread_;
            spread_zy                     = spread_2d_[l_z][l_y];

            int globalGridIndexXXStart = std::max(l_min[XX], grid_index_of_spread_atom_[XX] - ceilSqrtLUT[std::abs(d_l_z)][std::abs(d_l_y)]);
            int globalGridIndexXXEnd   = std::min(l_max[XX], grid_index_of_spread_atom_[XX] + ceilSqrtLUT[std::abs(d_l_z)][std::abs(d_l_y)]);
            int localGridIndexXXStart  = globalGridIndexXXStart - grid_index_of_spread_atom_[XX]+m_spread_;
            int numberSpreadVoxelsXX   = globalGridIndexXXEnd-globalGridIndexXXStart;
            densityDerivativeIterator        = densityDerivativeGridData.zy_column_begin(l_grid_z, l_grid_y)+globalGridIndexXXStart;

            spread_1d_XX = &(spread_1d_[XX][localGridIndexXXStart]);

            for (int l_x = 0; l_x <= numberSpreadVoxelsXX; ++l_x)
            {
                if ((*densityDerivativeIterator > GMX_FLOAT_EPS) || (*densityDerivativeIterator < GMX_FLOAT_EPS))
                {

                    /*
                     * The core routine that calcualtes k*V*(x-v)/sigma exp(-(x-v)^2/2*sigma^2) * (rho_div_v)
                     *
                     * *force_density_iterator = k*V/sigma exp(-(x-v)^2/2*sigma^2)
                     * shift_xyz = (x-v)
                     * *densityDerivativeIterator = rho_div_v
                     *
                     */
                    svmul(spread_zy * (*spread_1d_XX) * (*densityDerivativeIterator), shift_xyz, voxel_force);

                    /*
                     * add the force to the total force from all voxels
                     */
                    rvec_inc(force, voxel_force);
                }

                ++spread_1d_XX;
                ++densityDerivativeIterator;
                rvec_inc(shift_xyz, d_x); // next step in grid x-direction
            }
            rvec_inc(shift_yz, d_y);      // next step in grid y-direction
        }
        rvec_inc(shift_z, d_z);           // next step in grid z-direction
    }
    return force;
};

namespace externalpotential
{

std::unique_ptr<ExternalPotential> DensityFitting::create()
{
    return std::unique_ptr<ExternalPotential> (new DensityFitting());
}

real DensityFitting::single_atom_properties(GroupAtom * atom, t_mdatoms * /*mdatoms*/, gmx_localtop_t * /*topology_loc*/, const gmx_mtop_t * /*topology_global*/, const gmx_mtop_atomlookup * atom_lookup)
{
    t_atom * single_atom;
    gmx_mtop_atomnr_to_atom((const gmx_mtop_atomlookup_t)atom_lookup, *(atom->i_global), &single_atom);
    return atomicNumber2EmScatteringFactor(single_atom->atomnumber);
}

void DensityFitting::inv_mul(std::vector<real> &to_invert, const std::vector<real> &multiplier)
{
#pragma omp parallel for num_threads(number_of_threads_) schedule(static, to_invert.size()/number_of_threads_ + 1)
    for (std::size_t i = 0; i < to_invert.size(); ++i)
    {
        to_invert[i] = multiplier[i]/(to_invert[i]);
    }
}

RVec
DensityFitting::shiftedAndOriented(const RVec x)
{
    RVec x_translated;

    rvec_sub(x, centerOfMass_, x_translated);
    orientation_.rotate(x_translated);
    rvec_inc(x_translated, centerOfMass_);

    rvec_inc(x_translated, translation_);
    return x_translated;
}

void DensityFitting::spreadLocalAtoms_(WholeMoleculeGroup * spreadgroup)
{
    simulated_density_->zero();
    std::vector<volumedata::IVec>            minimumUsedGridIndex(number_of_threads_);
    std::vector<volumedata::IVec>            maximumUsedGridIndex(number_of_threads_);

#pragma omp parallel num_threads(number_of_threads_)
    {
        int           thread     = gmx_omp_get_thread_num();
        simulated_density_buffer_[thread]->zero();
        gauss_transform_[thread]->set_grid(std::move(simulated_density_buffer_[thread]));
        GroupIterator beginThreadAtoms = spreadgroup->begin(thread, number_of_threads_);
        GroupIterator endThreadAtoms   = spreadgroup->end(thread, number_of_threads_);
        for (auto atom = beginThreadAtoms; atom != endThreadAtoms; ++atom)
        {
            gauss_transform_[thread]->transform(shiftedAndOriented(*(*atom).xTransformed), *(*atom).properties);
        }
        simulated_density_buffer_[thread] = std::move(gauss_transform_[thread]->finish_and_return_grid());
        minimumUsedGridIndex[thread]      = gauss_transform_[thread]->getMinimumUsedGridIndex();
        maximumUsedGridIndex[thread]      = gauss_transform_[thread]->getMaximumUsedGridIndex();
        volumedata::MrcFile().write("threadlocal" + std::to_string(thread) + ".mrc", *(simulated_density_buffer_[thread]));

    }

    std::vector<std::vector<real>::iterator> contributingThreadLocalVoxelIterators(number_of_threads_);
    std::vector<real>::iterator              simulatedDensityVoxelIterator;

    volumedata::IVec gridStart;
    volumedata::IVec gridEnd;
    for (int i = XX; i <= ZZ; ++i)
    {
        auto ithElementLarger = [i](volumedata::IVec a, volumedata::IVec b){
                return a[i] < b[i];
            };
        gridStart[i] = (*std::min_element(std::begin(minimumUsedGridIndex), std::end(minimumUsedGridIndex), ithElementLarger))[i];
        gridEnd[i]   = (*std::max_element(std::begin(maximumUsedGridIndex), std::end(maximumUsedGridIndex), ithElementLarger))[i];
    }

    // add together all thread-local grids, using only density values, where there is
    // actual spread density
    int  nGridPointsXX        = gridEnd[XX]-gridStart[XX];
    auto simulatedDensityData = simulated_density_->access();

    // store the data accessors for threadlocal grid data in a vector
    std::vector<volumedata::GridDataAccess<real> > simulatedDensityThreadGridData;
    for (int thread = 0; thread < number_of_threads_; ++thread)
    {
        simulatedDensityThreadGridData.emplace_back(simulated_density_buffer_[thread]->extend(), simulated_density_buffer_[thread]->access().data());
    }

    //
    for (int gridIndexZZ = gridStart[ZZ]; gridIndexZZ < gridEnd[ZZ]; ++gridIndexZZ)
    {
        for (int gridIndexYY = gridStart[YY]; gridIndexYY < gridEnd[YY]; ++gridIndexYY)
        {
            contributingThreadLocalVoxelIterators.resize(0);
            // access rows in thread local grids if they contribute to the density
            for (int thread = 0; thread < number_of_threads_; ++thread)
            {
                if ((minimumUsedGridIndex[thread][ZZ] <= gridIndexZZ) && ( gridIndexZZ <= maximumUsedGridIndex[thread][ZZ] ) && (minimumUsedGridIndex[thread][YY] <= gridIndexYY) && (gridIndexYY <= maximumUsedGridIndex[thread][YY]))
                {
                    contributingThreadLocalVoxelIterators.push_back(simulatedDensityThreadGridData[thread].zy_column_begin(gridIndexZZ, gridIndexYY)+gridStart[XX]);
                }
            }
            simulatedDensityVoxelIterator = simulatedDensityData.zy_column_begin(gridIndexZZ, gridIndexYY) + gridStart[XX];
            // step though grid row by row
            for (int gridIndexXX = 0; gridIndexXX < nGridPointsXX; ++gridIndexXX)
            {
                // loop over the local threads, collect voxel contributions and advance all threadlocal iterators to next voxel in row
                for (auto &threadLocalVoxel : contributingThreadLocalVoxelIterators)
                {
                    *simulatedDensityVoxelIterator += *threadLocalVoxel;
                    ++threadLocalVoxel;
                }
                ++simulatedDensityVoxelIterator;
            }
        }
    }
}

real DensityFitting::getTotalScatteringSum_(WholeMoleculeGroup * atomgroup)
{
    real weightssum = atomgroup->local_weights_sum();
    if (mpi_helper() != nullptr)
    {
        mpi_helper()->to_reals_buffer(&weightssum, 1);
        mpi_helper()->sum_allReduce();
    }
    return weightssum;
}

void DensityFitting::setCenterOfMass(WholeMoleculeGroup * atomgroup)
{
    centerOfMass_ = atomgroup->weighted_local_coordinate_sum();
    real weightssum = atomgroup->local_weights_sum();
    if (mpi_helper() != nullptr)
    {
        mpi_helper()->to_reals_buffer(centerOfMass_, 3);
        mpi_helper()->to_reals_buffer(&weightssum, 1);
        mpi_helper()->sum_reduce();
        if (mpi_helper()->isMaster())
        {
            mpi_helper()->from_reals_buffer(centerOfMass_, 3);
            mpi_helper()->from_reals_buffer(&weightssum, 1);
        }
    }
    if (mpi_helper() == nullptr || mpi_helper()->isMaster())
    {
        svmul(1.0/weightssum, centerOfMass_, centerOfMass_);
    }

}

void DensityFitting::sumReduceNormalize_()
{
    if (mpi_helper() != nullptr)
    {
        sum_reduce_simulated_density_();
    }
    if (mpi_helper() == nullptr || mpi_helper()->isMaster())
    {
        simulated_density_->multiply(simulated_density_->grid_cell_volume() * totalScatteringSum_);
    }
}

real DensityFitting::KLDivergenceFromTargetOnMaster(WholeMoleculeGroup * atomgroup)
{
    spreadLocalAtoms_(atomgroup);
    sumReduceNormalize_();
    if (mpi_helper() == nullptr || mpi_helper()->isMaster())
    {
        return volumedata::GridMeasures(*target_density_).getKLSameGrid(*simulated_density_);
    }
    else
    {
        return 0;
    }
}

void
DensityFitting::alignComDensityAtoms()
{
    if (mpi_helper() == nullptr || mpi_helper()->isMaster())
    {
        rvec_sub(target_density_->center_of_mass(), centerOfMass_, translation_);
    }

    if (mpi_helper() != nullptr)
    {
        mpi_helper()->broadcast(translation_, 3);
    }
}

bool
DensityFitting::optimizeTranslation(WholeMoleculeGroup * translationgroup, real &divergenceToCompareTo)
{
    const real step_size                  = 0.02;
    bool       succededReducingDivergence = false;
    // calculate gradient vector
    // normalize to step size
    // move along gradient vector
    // accept translation if new potential better
    auto translationBeforeTrial = translation_;
    spreadLocalAtoms_(translationgroup);
    volumedata::MrcFile().write("test.mrc", *simulated_density_);
    sumReduceNormalize_();

    invertMultiplySimulatedDensity_();
    KLForceCalculation_(translationgroup);
    auto force = translationgroup->local_force_sum();
    if (mpi_helper() != nullptr)
    {
        mpi_helper()->sum_reduce_rvec(force);
    }

    if (mpi_helper() == nullptr || mpi_helper()->isMaster())
    {
        unitv(force, force);
    }

    fprintf(input_output()->output_file(), "#\tMinimization\n"
            "#\tMinimization - New direction:\t[ %g %g %g ]\n"
            "#\tMinimization - Step size:\t%g nm .\n", force[XX], force[YY], force[ZZ], step_size);
    RVec force_scaled;
    if (mpi_helper() == nullptr || mpi_helper()->isMaster())
    {
        svmul(step_size, force, force_scaled);
        rvec_inc(translation_, force_scaled);
        fprintf(input_output()->output_file(), "#\tMinimization - Move attemted to\t[ %g %g %g ]\n", force_scaled[XX], force_scaled[YY], force_scaled[ZZ]);
    }
    if (mpi_helper() != nullptr)
    {
        mpi_helper()->broadcast(&translation_, 1);
    }
    auto new_divergence =  std::abs(KLDivergenceFromTargetOnMaster(translationgroup));
    if (mpi_helper() == nullptr || mpi_helper()->isMaster())
    {

        fprintf(input_output()->output_file(), "#\tMinimization - Divergence before move: %g , after move: %g\n", divergenceToCompareTo, new_divergence);

        if (new_divergence < divergenceToCompareTo)
        {
            succededReducingDivergence   = true;
            divergenceToCompareTo        = new_divergence;
        }
        else
        {
            translation_ = translationBeforeTrial;
        }
    }
    return succededReducingDivergence;
}

bool
DensityFitting::optimizeOrientation(WholeMoleculeGroup * atomgroup, real &divergenceToCompareTo)
{
    const real rotationAngle               = M_PI/128.;
    bool       succeededReducingDivergence = false;

    // calculate gradient vector
    // normalize to step size
    // move along gradient vector
    // accept if new potential better than old
    // move back along gradient vector

    spreadLocalAtoms_(atomgroup);
    sumReduceNormalize_();
    invertMultiplySimulatedDensity_();
    KLForceCalculation_(atomgroup);

    auto torque_direction = atomgroup->local_torque_sum(centerOfMass_); //TODO: check for sign of rotation

    if (mpi_helper() != nullptr)
    {
        mpi_helper()->sum_reduce_rvec(torque_direction);
    }

    if (mpi_helper() == nullptr || mpi_helper()->isMaster())
    {
        unitv(torque_direction, torque_direction);
    }

    fprintf(input_output()->output_file(), "#\tMinimization in Orientation Space \n"
            "#\tMinimization - Torque Vector:\t[ %g %g %g ]\n",
            torque_direction[XX], torque_direction[YY], torque_direction[ZZ]);

    auto orientationBeforeTrial = orientation_;

    if (mpi_helper() == nullptr || mpi_helper()->isMaster())
    {
        orientation_ *= Quaternion(torque_direction, rotationAngle);
        orientation_.normalize();
    }
    if (mpi_helper() != nullptr)
    {
        mpi_helper()->broadcast(&orientation_[0], 4);
    }

    auto new_divergence =  std::abs(KLDivergenceFromTargetOnMaster(atomgroup));

    if (mpi_helper() == nullptr || mpi_helper()->isMaster())
    {

        fprintf(input_output()->output_file(), "#\tMinimization - Divergence before rotation: %g , after rotation: %g\n", divergenceToCompareTo, new_divergence);

        if (new_divergence < divergenceToCompareTo)
        {
            succeededReducingDivergence = true;
            divergenceToCompareTo       = new_divergence;
        }
        else
        {
            orientation_ = orientationBeforeTrial;
            if (mpi_helper() != nullptr)
            {
                mpi_helper()->broadcast(&orientation_[0], 4);
            }
        }
    }

    return succeededReducingDivergence;
}


void DensityFitting::translate_atoms_into_map_(WholeMoleculeGroup * translationgroup)
{

    std::vector<real> upperHalfGrid {
        1.
    };
    std::vector<real> gridPoints {
        0
    };
    Quaternion bestOrientation                  = orientation_;
    RVec       bestTranslation                  = translation_;
    real       bestDivergence                   = std::abs(KLDivergenceFromTargetOnMaster(translationgroup));
    auto       bestDivergenceCurrentOrientation = bestDivergence;
    auto       totalNumberOrientations          = upperHalfGrid.size()*gridPoints.size()*gridPoints.size()*gridPoints.size();
    decltype(totalNumberOrientations) currentOrientationNumber = 0;

    for (auto q0 : upperHalfGrid)
    {
        for (auto q1 : gridPoints)
        {
            for (auto q2 : gridPoints)
            {
                for (auto q3 : gridPoints)
                {

                    fprintf(stderr, "Optimizing Orientation: %3lu / %3lu = ", ++currentOrientationNumber, totalNumberOrientations);

                    orientation_ = Quaternion(Quaternion::QVec {{q0, q1, q2, q3}});
                    orientation_.normalize();
                    fprintf(stderr, "[ %3g %3g %3g %3g ]", orientation_[0], orientation_[1], orientation_[2], orientation_[3]);

                    bestDivergenceCurrentOrientation =  std::abs(KLDivergenceFromTargetOnMaster(translationgroup));

                    while (optimizeTranslation(translationgroup, bestDivergenceCurrentOrientation) || optimizeOrientation(translationgroup, bestDivergenceCurrentOrientation))
                    {
                        fprintf(stderr, "\r\t\t\t\t\t\t\tcurrent fit = %g ; best fit = %g \t\t", bestDivergenceCurrentOrientation, bestDivergence);
                    }

                    if (bestDivergenceCurrentOrientation < bestDivergence)
                    {
                        bestTranslation = translation_;
                        bestOrientation = orientation_;
                        bestDivergence  = bestDivergenceCurrentOrientation;
                    }
                }
            }
        }
    }
    translation_ = bestTranslation;
    orientation_ = bestOrientation;
    orientation_.normalize();
}

void
DensityFitting::initialize_target_density_()
{
    target_density_->normalize();
}

void
DensityFitting::initialize_buffers_()
{
    force_density_.resize(number_of_threads_);
    force_gauss_transform_.resize(number_of_threads_);
    simulated_density_buffer_.resize(number_of_threads_);
    gauss_transform_.resize(number_of_threads_);
// ensure thread local buffers
#pragma omp parallel num_threads(number_of_threads_)
    {
        int           thread     = gmx_omp_get_thread_num();

        auto          forceDensityTmp = std::unique_ptr<volumedata::GridReal>(new volumedata::GridReal());
        forceDensityTmp->copy_grid(*target_density_);

        force_density_[thread] = std::move(forceDensityTmp);

        /* Create a gauss transform object for each density buffer,
         * because the gauss transforms will be carried out simultaneously
         */
        auto forceGaussTransformTmp = std::unique_ptr<volumedata::FastGaussianGriddingForce>(new volumedata::FastGaussianGriddingForce());
        forceGaussTransformTmp->set_sigma(sigma_);
        forceGaussTransformTmp->set_n_sigma(n_sigma_);

        force_gauss_transform_[thread] = std::move(forceGaussTransformTmp);

        auto simulatedDensityTmp = std::unique_ptr<volumedata::GridReal>(new volumedata::GridReal());
        simulatedDensityTmp->copy_grid(*target_density_);
        simulated_density_buffer_[thread] = std::move(simulatedDensityTmp);
        auto gaussTransformTmp = std::unique_ptr<volumedata::FastGaussianGridding>(new volumedata::FastGaussianGridding());
        gauss_transform_[thread] = std::move(gaussTransformTmp);
    }
}

void
DensityFitting::initialize_spreading_()
{
    simulated_density_->copy_grid(*target_density_);
#pragma omp parallel num_threads(number_of_threads_)
    {
        int           thread     = gmx_omp_get_thread_num();
        gauss_transform_[thread]->set_sigma(sigma_);
        gauss_transform_[thread]->set_n_sigma(n_sigma_);
    }
}

void
DensityFitting::initializeKL_(const matrix box, const rvec x[])
{
    if (bWriteXTC_)
    {
        out_  = open_xtc(trajectory_name_.c_str(), "w");
    }
    auto fitatoms = wholemoleculegroup(x, box, 0);
    setCenterOfMass(fitatoms);
    alignComDensityAtoms();

    initialize_target_density_();
    initialize_buffers_();

    totalScatteringSum_ = getTotalScatteringSum_(fitatoms);
    initialize_spreading_();

    spreadLocalAtoms_(fitatoms);
    sumReduceNormalize_();

    reference_density_ = simulated_density_->access().data();

    if (isCenterOfMassCentered_)
    {
        translate_atoms_into_map_(fitatoms);
    }

    absolute_target_divergence_ = volumedata::GridMeasures(*target_density_).getKLSameGrid(*simulated_density_);
    const real timeStepNs = 0.004; // TODO: pass this from simulation input
    optimalDeltaPotentialEnergy_ = fitatoms->num_atoms_global()*std::sqrt((float)every_nth_step_) * timeStepNs * maximumEnergyFluctuationPerAtom_;
    // fitatoms->medianSort();

    fprintf(input_output()->output_file(), "#KL-divergence(target|simulated): %g.\n", absolute_target_divergence_);
    fprintf(input_output()->output_file(), "#---Step-size estimate in DensityFitting potential space---.\n");
    fprintf(input_output()->output_file(), "#   Maximum energy fluctuation per atom: %g.\n", maximumEnergyFluctuationPerAtom_);
    fprintf(input_output()->output_file(), "#   Number of refined atoms: %d .\n", fitatoms->num_atoms_global());
    fprintf(input_output()->output_file(), "#   Fitting every %dth step.\n", every_nth_step_);
    fprintf(input_output()->output_file(), "#   With time step %g.\n", timeStepNs);
    fprintf(input_output()->output_file(), "#   Thus estimating optimal change in refinement potental energy: %g .\n", optimalDeltaPotentialEnergy_);
}

void
DensityFitting::initialize(const matrix box, const rvec x[])
{
    initialize_(box, x);
}

void DensityFitting::sum_reduce_simulated_density_()
{
    mpi_helper()->to_reals_buffer(simulated_density_->access().data().data(), simulated_density_->access().data().size() );
    mpi_helper()->sum_reduce();
    if (mpi_helper()->isMaster())
    {
        mpi_helper()->from_reals_buffer(simulated_density_->access().data().data(), simulated_density_->access().data().size());
    }
}

void DensityFitting::plot_forces(WholeMoleculeGroup * plotatoms)
{

    externalpotential::ForcePlotter plot;
    static int i = 0;
    plot.start_plot_forces("forces" + std::to_string(++i) + ".bild");
    real       max_f    = -1;
    for (auto &atom : *plotatoms)
    {
        if (norm(*atom.force) > max_f)
        {
            max_f = norm(*atom.force);
        }
    }
    fprintf(stderr, "max force:  %g \n", max_f);
    for (auto &atom : *plotatoms)
    {
        plot.plot_force(shiftedAndOriented(*atom.xTransformed), *atom.force, 0, 1.0/max_f);
    }
    plot.stop_plot_forces();

}

void DensityFitting::KLForceCalculation_(WholeMoleculeGroup * fitatoms)
{
#pragma omp parallel num_threads(number_of_threads_) shared(stderr,fitatoms) default(none)
    {
        int           thread             = gmx_omp_get_thread_num();
        auto         &threadlocaldensity = *simulated_density_buffer_[thread];
        threadlocaldensity.access().data() = (*simulated_density_).access().data();
        real          prefactor  = k_/(norm_simulated_*sigma_*sigma_);

        force_gauss_transform_[thread]->set_grid(std::move(force_density_[thread]));
        auto         &threadlocal_force_gauss_transform = force_gauss_transform_[thread];

        GroupIterator beginThreadAtoms = fitatoms->begin(thread, number_of_threads_);
        GroupIterator endThreadAtoms   = fitatoms->end(thread, number_of_threads_);

        for (auto atom = beginThreadAtoms; atom != endThreadAtoms; ++atom)
        {
            /* calculate for atom position x, voxel v, grid cell Volume V, force constant k:
             *
             * sum_v k*V*(x-v)/sigma exp(-(x-v)^2/2*sigma^2) * (rho_div_v)
             *
             * to keep the potential translation invariant, use shifted atoms positions for force calculations.
             * Atoms are translated, such that the potential is minimal with respect to coordinate shifts
             * To speed up the expensive exp(-(x-v)^2/2*sigma^2) part of the force caluclation use fast gaussian gridding
             * Each thread has its own pre-allocated memory for doing the fast gaussian gridding.
             *
             * There is some double work done here, since all atoms have already been spread on a grid to calculate the total spread density.
             * However, keeping the spread grid for all atoms appears to be very memory intense. (natoms * voxelgridsize)
             * The atoms spread weight is k_/(norm_simulated * sigma_^2), so we don't have to do that multiplication later in the force calculation loop
             * simulated_density_->grid_cell_volume()/
             */
            copy_rvec( threadlocal_force_gauss_transform->force(shiftedAndOriented(*(*atom).xTransformed), *((*atom).properties)*prefactor, threadlocaldensity), *(*atom).force);
            orientation_.rotate_backwards(*(*atom).force);
        }
        force_density_[thread] = std::move(threadlocal_force_gauss_transform->finish_and_return_grid());
    }
}

void DensityFitting::invertMultiplySimulatedDensity_()
{
    if (mpi_helper() == nullptr || mpi_helper()->isMaster())
    {
        inv_mul(simulated_density_->access().data(), target_density_->access().data());
    }

    if (mpi_helper() != nullptr)
    {
        mpi_helper()->broadcast(simulated_density_->access().data().data(), simulated_density_->access().data().size());
        mpi_helper()->broadcast(&norm_simulated_, 1);
    }

}

void
DensityFitting::writeTranslatedCoordinates_(WholeMoleculeGroup * atoms, int step)
{
    if (bWriteXTC_)
    {
        std::vector<RVec> x_out(atoms->num_atoms_global());
        int               i = 0;
        for (auto atom : *atoms)
        {
            x_out[i++] = shiftedAndOriented(*atom.xTransformed);
        }
        write_xtc(out_, atoms->num_atoms_global(), step, 0, *(atoms->box()), as_vec_array(x_out.data()), 1000);
    }
}

void
DensityFitting::doPotentialINV_( const matrix box, const rvec x[], const gmx_int64_t /*step*/)
{
    auto fitatoms = wholemoleculegroup(x, box, 0);
#pragma omp parallel num_threads(number_of_threads_)
    {
        std::array<volumedata::GridDataAccess<real>, 3> invertedDensityForcesAccess = {{invertedDensityForces_[XX]->access(), invertedDensityForces_[YY]->access(), invertedDensityForces_[ZZ]->access()}};
        for (auto atom : *fitatoms)
        {
            IVec gridCoordinate     = target_density_->coordinate_to_gridindex_floor_ivec(shiftedAndOriented(*atom.xTransformed));
            clear_rvec(*atom.force);
            if (target_density_->inGrid(gridCoordinate))
            {
                for (size_t dim = XX; dim <= ZZ; dim++)
                {
                    *atom.force[dim] = invertedDensityForcesAccess[dim].at(gridCoordinate);
                }
                orientation_.rotate_backwards(*atom.force);
            }
        }
    }
}

void
DensityFitting::doPotentialCC_( const matrix /*box*/, const rvec /*x*/[], const gmx_int64_t /*step*/)
{

}

void DensityFitting::doPotentialKL_( const matrix box, const rvec x[], const gmx_int64_t step)
{
    auto fitatoms = wholemoleculegroup(x, box, 0);
    setCenterOfMass(fitatoms);

    if ( (step%(every_nth_step_) == 0) && (bWriteXTC_))
    {
        writeTranslatedCoordinates_(fitatoms, step);
    }

    /*
     * Spread all local atoms on a grid
     */
    spreadLocalAtoms_(fitatoms);
    /*
     * Gather the local contributions to the overall spread density on the master node.
     */
    sumReduceNormalize_();
    if (mpi_helper() == nullptr || mpi_helper()->isMaster())
    {
        auto delta_divergence  = volumedata::GridMeasures(*target_density_).getRelativeKLCrossTermSameGrid(*simulated_density_, reference_density_);
        reference_divergence_ -= delta_divergence;
        set_local_potential(k_*reference_divergence_);
        reference_density_             = simulated_density_->access().data();
        exponentialDeltaEnergyAverage_ = 0.3* (k_ * delta_divergence) + 0.7 * exponentialDeltaEnergyAverage_;
        fprintf(input_output()->output_file(), "%8g\t%8g\t%8g\t", reference_divergence_, k_, exponentialDeltaEnergyAverage_);

        if (exponentialDeltaEnergyAverage_ < optimalDeltaPotentialEnergy_)
        {
            k_ *= k_factor_;
        }
        else
        {
            k_ /= k_factor_;
        }
    }
    else
    {
        set_local_potential(0);
    }

    invertMultiplySimulatedDensity_();

    KLForceCalculation_(fitatoms);

    // volumedata::MrcFile simulated_output;
    if (step % (every_nth_step_) == 0)
    {
        // simulated_output.write("simulated.mrc", *simulated_density_);
        plot_forces(fitatoms);
    }
    // fitatoms->parallel_loop(std::bind( &DensityFitting::ForceKernel_KL, this, std::placeholders::_1, std::placeholders::_2));

}

void
DensityFitting::initializeINV_(const matrix /*box*/, const rvec * /*x*/)
{
    //
    for (size_t dimension = 0; dimension <= ZZ; dimension++)
    {
        invertedDensityForces_.emplace_back(std::unique_ptr<volumedata::GridReal>(new volumedata::GridReal()));
        invertedDensityForces_.back()->copy_grid(*target_density_);
        invertedDensityForces_.back()->zero();
        auto spacing                     = invertedDensityForces_.back()->avg_spacing();
        auto invertedDensityForcesAccess = invertedDensityForces_.back()->access();
        auto targetDensityAccess         = target_density_->access();
        // calculate numerical gradient of the map
        for (int i_z = 1; i_z < invertedDensityForces_.back()->extend()[ZZ]-1; i_z++)
        {
            for (int i_y = 1; i_y < invertedDensityForces_.back()->extend()[ZZ]-1; i_y++)
            {
                for (int i_x = 1; i_x < invertedDensityForces_.back()->extend()[ZZ]-1; i_x++)
                {
                    auto ref = targetDensityAccess.at({i_x-1, i_y, i_z});
                    decltype(ref) cumulant = 0.;
                    ivec ii;
                    for (ii[ZZ] = -1; ii[ZZ] <= 1; ++ii[ZZ])
                    {
                        for (ii[YY] = -1; ii[YY] <= 1; ++ii[YY])
                        {
                            for (ii[XX] = -1; ii[XX] <= 1; ++ii[XX])
                            {
                                cumulant += ii[dimension] * targetDensityAccess.at({i_x+ii[XX], i_y+ii[YY], i_z+ii[ZZ]});
                            }
                        }
                    }
                    invertedDensityForcesAccess.at({i_x, i_y, i_z}) = spacing * cumulant / 18.;
                }
            }
        }
    }
}

void
DensityFitting::initializeCC_(const matrix /*box*/, const rvec * /*x*/)
{

}

void DensityFitting::do_potential( const matrix box, const rvec x[], const gmx_int64_t step)
{
    doPotential_(box, x, step);
};

DensityFitting::DensityFitting() : ExternalPotential(),
                                   target_density_(std::unique_ptr<volumedata::GridReal>(new volumedata::GridReal())),
                                   simulated_density_(std::unique_ptr<volumedata::GridReal>(new volumedata::GridReal())),
                                   translation_({0, 0, 0}
                                                ), number_of_threads_ {
    std::max(1, gmx_omp_nthreads_get(emntDefault))
},
exponentialDeltaEnergyAverage_ {
    0.0
}, orientation_ {
    Quaternion::QVec({{1., 0., 0., 0.}})
}
{
}

bool DensityFitting::do_this_step(gmx_int64_t step)
{
    return (step % every_nth_step_ == 0 );
}


void DensityFitting::read_input()
{
    FILE      * inputfile       = input_output()->input_file();
    if (inputfile == nullptr)
    {
        GMX_THROW(InconsistentInputError("Please provide an external potential input file."));
    }
    fseek(inputfile, 0, SEEK_END);
    long fsize = ftell(inputfile);
    fseek(inputfile, 0, SEEK_SET);  //same as rewind(f);

    auto line = (char*) malloc(fsize + 1);;
    fsize = fread(line, fsize, 1, inputfile);
    // TODO: implement JSON scheme for checking input consistency
    json::Object parsed_json {
        std::string(line)
    };

    if (parsed_json.has("k"))
    {
        k_                  = std::stof(parsed_json["k"]);
    }
    else
    {
        GMX_THROW(InconsistentInputError("Force constant needed for densityfitting. Please provide a real number " "k" " in JSON input file."));
    };

    if (parsed_json.has("sigma"))
    {
        sigma_              = std::stof(parsed_json["sigma"]);
    }
    else
    {
        fprintf(stderr, "\n No mobility estimate given, guessing sigma = 0.2 nm . \n");
        sigma_ = 0.2;
    }

    if (parsed_json.has("n_sigma"))
    {
        n_sigma_            = std::stof(parsed_json["n_sigma"]);
    }
    else
    {
        fprintf(stderr, "\n No density spread range provided, guessing n_sigma = 5 . \n");
        n_sigma_ = 5;
    }
    if (parsed_json.has("method"))
    {
        fitMethod_ = parsed_json["method"];
    }
    else
    {
        fitMethod_ = "KL";
    }

    if (fitMethod_ == "KL")
    {
        initialize_  = std::bind(&DensityFitting::initializeKL_, this, std::placeholders::_1, std::placeholders::_2);
        doPotential_ = std::bind(&DensityFitting::doPotentialKL_, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
    }

    if (fitMethod_ == "INV")
    {
        initialize_  = std::bind(&DensityFitting::initializeINV_, this, std::placeholders::_1, std::placeholders::_2);
        doPotential_ = std::bind(&DensityFitting::doPotentialINV_, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
    }


    if (fitMethod_ == "CC")
    {
        initialize_  = std::bind(&DensityFitting::initializeCC_, this, std::placeholders::_1, std::placeholders::_2);
        doPotential_ = std::bind(&DensityFitting::doPotentialCC_, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
    }

    if (parsed_json.has("energy_step"))
    {
        maximumEnergyFluctuationPerAtom_ = std::stof(parsed_json.at("energy_step"));
    }
    else
    {
        maximumEnergyFluctuationPerAtom_ = 1e-4;
    }

    if (parsed_json.has("k_factor"))
    {
        k_factor_ = std::stof(parsed_json.at("k_factor"));
    }
    else
    {
        k_factor_ = 1;
    }
    if (parsed_json.has("center_to_density"))
    {
        isCenterOfMassCentered_ = parsed_json.at("center_to_density") == "true" || parsed_json.at("center_to_density") == "yes";
    }
    else
    {
        isCenterOfMassCentered_ = false;
    }

    if (parsed_json.has("write"))
    {
        bWriteXTC_ = (parsed_json.at("write") == "true") || (parsed_json.at("write") == "yes");
    }
    else
    {
        bWriteXTC_ = false;
    }

    if (parsed_json.has("every_nth_step"))
    {
        every_nth_step_ = std::stoi(parsed_json.at("every_nth_step"));
        k_             *= every_nth_step_;
    }
    else
    {
        every_nth_step_ = 1;
    }

    if (parsed_json.has("target_density"))
    {
        volumedata::MrcFile  target_input_file;
        target_density_name_ = parsed_json["target_density"];
        target_input_file.read(target_density_name_, *target_density_);
        if (parsed_json.has("target_density_used"))
        {
            volumedata::MrcFile target_output_file;
            target_output_file.write(parsed_json["target_density_used"], *target_density_);
        }
        else
        {
            volumedata::MrcFile target_output_file;
            target_output_file.write("target.ccp4", *target_density_);
        }
    }
    else
    {
        GMX_THROW(InconsistentInputError("No target_density name given in input file. A target em-density map is required for refining atom coordinates against it."));
    }


    if (parsed_json.has("trajectory"))
    {
        trajectory_name_ = parsed_json["trajectory"];
    }
    else
    {
        trajectory_name_ = "out.xtc";
    }

    fprintf(input_output()->output_file(), "%s", dumpParsedInput().c_str());
    fprintf(input_output()->output_file(), "KL-div\tk\t%%match\tF_max\tweight\n");


}

std::string
DensityFitting::dumpParsedInput()
{
    std::string result("#input: \n");
    result += "#input:\tk\t\t = " + std::to_string(k_) + "\n";
    result += "#input:\tsigma\t\t = " + std::to_string(sigma_) + "\n";
    result += "#input:\tn_sigma\t\t = " + std::to_string(n_sigma_) + "\n";
    result += "#input:\tbackground_density\t\t = " + std::to_string(background_density_) + "\n";
    result += "#input:\tenergy_step\t\t = " + std::to_string(maximumEnergyFluctuationPerAtom_) + "\n";
    result += "#input:\tk_factor\t\t = " + std::to_string(k_factor_) + "\n";
    result += "#input:\tcenter_to_density\t\t = " + std::to_string(isCenterOfMassCentered_) + "\n";
    result += "#input:\tevery_nth_step\t\t = " + std::to_string(every_nth_step_) + "\n";
    result += "#input:\tmethod\t\t = " + fitMethod_ + "\n";
    result += "#input:\ttarget_density\t\t = " + target_density_name_ + "\n";

    return result;
}

void DensityFitting::broadcast_internal()
{
    if (mpi_helper() != nullptr)
    {
        mpi_helper()->broadcast(&k_, 1);
        mpi_helper()->broadcast(&k_factor_, 1);
        mpi_helper()->broadcast(&every_nth_step_, 1);
        mpi_helper()->broadcast(&isCenterOfMassCentered_, 1);
        mpi_helper()->broadcast(&sigma_, 1);
        mpi_helper()->broadcast(&n_sigma_, 1);
    }
}

void DensityFitting::finish()
{
    if (bWriteXTC_)
    {
        close_xtc(out_);
    }
}

std::string DensityFittingInfo::name                        = "densityfitting";
std::string DensityFittingInfo::shortDescription            = "calculate forces from difference to target density";
const int   DensityFittingInfo::numberIndexGroups           = 1;
const int   DensityFittingInfo::numberWholeMoleculeGroups   = 1;
externalpotential::ModuleCreator DensityFittingInfo::create = DensityFitting::create;

} // namespace externalpotential
} // namespace gmx
