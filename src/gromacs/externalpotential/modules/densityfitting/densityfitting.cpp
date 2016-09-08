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

#include "gromacs/externalpotential/atomgroups/wholemoleculegroup.h"
#include "gromacs/externalpotential/mpi-helper.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/math/volumedata.h"
#include "gromacs/fileio/volumedataio.h"
#include "gromacs/fileio/json.h"
#include "gromacs/math/gausstransform.h"
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
volumedata::FastGaussianGriddingForce::force(const rvec x, const real weight, const volumedata::GridReal &div)
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

    real   spread_zy;
    std::vector<real>::iterator density_ratio_iterator;
    real * spread_1d_XX;


    RVec     d_x       = div.unit_cell_XX();
    RVec     d_y       = div.unit_cell_YY();
    RVec     d_z       = div.unit_cell_ZZ();
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

    rvec_sub(div.gridpoint_coordinate({0, 0, 0}), x, shift_z); // using the grid-origin as reference (v_0-x)
    rvec_inc(shift_z, d_z_offset);                             // move in "z-slices" until spreading starts

    for (int l_grid_z = l_min[ZZ]; l_grid_z <= l_max[ZZ]; ++l_grid_z)
    {
        shift_yz = shift_z;                                              // start at the beginning of a "y-row"
        rvec_inc(shift_yz, d_y_offset);                                  // move into the "y-row" to the voxel where spreading starts

        int d_l_z                  = l_grid_z - grid_index_of_spread_atom_[ZZ];
        int l_z                    = d_l_z + m_spread_; // the grid index in the 2d spread grid around the atom
        int globalGridIndexYYStart = std::max(l_min[YY], grid_index_of_spread_atom_[YY] - (int)std::ceil(sqrt(m_spread_*m_spread_ - d_l_z*d_l_z)));
        int globalGridIndexYYEnd   = std::min(l_max[YY], grid_index_of_spread_atom_[YY] + (int)std::ceil(sqrt(m_spread_*m_spread_ - d_l_z*d_l_z)));
        for (int l_grid_y = globalGridIndexYYStart; l_grid_y <= globalGridIndexYYEnd; ++l_grid_y)
        {
            shift_xyz = shift_yz;            // start at the beginning of an "x-column"
            rvec_inc(shift_xyz, d_x_offset); // move into the "x-row" to the voxel where spreading starts

            int d_l_y        =  l_grid_y - grid_index_of_spread_atom_[YY];
            int l_y          = d_l_y+m_spread_;
            spread_zy                     = spread_2d_[l_z][l_y];

            int globalGridIndexXXStart = std::max(l_min[XX], grid_index_of_spread_atom_[XX] - (int)std::ceil(sqrt(m_spread_*m_spread_ - d_l_z*d_l_z- d_l_y*d_l_y)));
            int globalGridIndexXXEnd   = std::min(l_max[XX], grid_index_of_spread_atom_[XX] + (int)std::ceil(sqrt(m_spread_*m_spread_ - d_l_z*d_l_z- d_l_y*d_l_y)));
            int localGridIndexXXStart  = globalGridIndexXXStart - grid_index_of_spread_atom_[XX]+m_spread_;
            int numberSpreadVoxelsXX   = globalGridIndexXXEnd-globalGridIndexXXStart;
            density_ratio_iterator        = div.zy_column_begin(l_grid_z, l_grid_y)+globalGridIndexXXStart;

            spread_1d_XX = &(spread_1d_[XX][localGridIndexXXStart]);

            for (int l_x = 0; l_x <= numberSpreadVoxelsXX; ++l_x)
            {
                /*
                 * The core routine that calcualtes k*V*(x-v)/sigma exp(-(x-v)^2/2*sigma^2) * (rho_div_v)
                 *
                 * *force_density_iterator = k*V/sigma exp(-(x-v)^2/2*sigma^2)
                 * shift_xyz = (x-v)
                 * *density_ratio_iterator = rho_div_v
                 *
                 */
                svmul(spread_zy * (*spread_1d_XX) * (*density_ratio_iterator), shift_xyz, voxel_force);

                /*
                 * add the force to the total force from all voxels
                 */
                rvec_inc(force, voxel_force);

                ++spread_1d_XX;
                ++density_ratio_iterator;
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

void DensityFitting::SpreadKernel(GroupAtom &atom, const int &thread)
{
    rvec x_translated;
    rvec_add(*atom.xTransformed, translation_, x_translated);
    gauss_transform_[thread]->transform(x_translated, *(atom.properties));
}

void DensityFitting::spread_density_(WholeMoleculeGroup * spreadgroup)
{
    simulated_density_->zero();
    std::vector<volumedata::IVec>            minimumUsedGridIndex(number_of_threads_);
    std::vector<volumedata::IVec>            maximumUsedGridIndex(number_of_threads_);

#pragma omp parallel num_threads(number_of_threads_)
    {
        int           thread     = gmx_omp_get_thread_num();
        simulated_density_buffer_[thread]->zero();
        gauss_transform_[thread]->set_grid(std::move(simulated_density_buffer_[thread]));
        GroupIterator first_atom_for_thread = spreadgroup->begin(thread, number_of_threads_);
        GroupIterator end_of_thread_atoms   = spreadgroup->end(thread, number_of_threads_);
        for (auto atom = first_atom_for_thread; atom != end_of_thread_atoms; ++atom)
        {
            rvec x_translated;
            rvec_add(*(*atom).xTransformed, translation_, x_translated);
            gauss_transform_[thread]->transform(x_translated, *(*atom).properties);
        }
        simulated_density_buffer_[thread] = std::move(gauss_transform_[thread]->finish_and_return_grid());
        minimumUsedGridIndex[thread]      = gauss_transform_[thread]->getMinimumUsedGridIndex();
        maximumUsedGridIndex[thread]      = gauss_transform_[thread]->getMaximumUsedGridIndex();
    }

    std::vector<std::vector<real>::iterator> it_buffer(number_of_threads_);
    std::vector<real>::iterator              it_simulated_density;

    volumedata::IVec gridStart;
    volumedata::IVec gridEnd;
    for (int i = XX; i <= ZZ; ++i)
    {
        gridStart[i] = (*std::min_element(std::begin(minimumUsedGridIndex), std::end(minimumUsedGridIndex), [i](volumedata::IVec a, volumedata::IVec b){return a[i] < b[i]; }))[i];
        gridEnd[i]   = (*std::max_element(std::begin(maximumUsedGridIndex), std::end(maximumUsedGridIndex), [i](volumedata::IVec a, volumedata::IVec b){return a[i] < b[i]; }))[i];
    }

    int nGridPointsXX = gridEnd[XX]-gridStart[XX];
    for (int gridIndexZZ = gridStart[ZZ]; gridIndexZZ < gridEnd[ZZ]; ++gridIndexZZ)
    {
        for (int gridIndexYY = gridStart[YY]; gridIndexYY < gridEnd[YY]; ++gridIndexYY)
        {
            it_buffer.resize(0);
            for (int thread = 0; thread < number_of_threads_; ++thread)
            {
                if ((minimumUsedGridIndex[thread][ZZ] <= gridIndexZZ) && ( gridIndexZZ <= maximumUsedGridIndex[thread][ZZ] ) && (minimumUsedGridIndex[thread][YY] <= gridIndexYY) && (gridIndexYY <= maximumUsedGridIndex[thread][YY]))
                {
                    it_buffer.push_back(simulated_density_buffer_[thread]->zy_column_begin(gridIndexZZ, gridIndexYY)+gridStart[XX]);
                }
            }
            it_simulated_density = simulated_density_->zy_column_begin(gridIndexZZ, gridIndexYY) + gridStart[XX];

            for (int gridIndexXX = 0; gridIndexXX < nGridPointsXX; ++gridIndexXX)
            {
                for (auto voxel = it_buffer.begin(); voxel != it_buffer.end(); ++voxel)
                {
                    *it_simulated_density += **voxel;
                    ++(*voxel);
                }
                ++it_simulated_density;
            }
        }
    }

}

void DensityFitting::translate_atoms_into_map_(const rvec x[], const matrix box)
{
    auto translationgroup            = wholemoleculegroup(x, box, 0);
    RVec center_of_geometry_fitatoms = translationgroup->local_coordinate_sum();
    if (mpi_helper() != nullptr)
    {
        mpi_helper()->to_reals_buffer(center_of_geometry_fitatoms, 3);
        mpi_helper()->sum_reduce();
        if (mpi_helper()->isMaster())
        {
            mpi_helper()->from_reals_buffer(center_of_geometry_fitatoms, 3);
        }
    }

    if (mpi_helper() == nullptr || mpi_helper()->isMaster())
    {
        svmul(1.0/translationgroup->num_atoms_global(), center_of_geometry_fitatoms, center_of_geometry_fitatoms);
        RVec center_of_density = target_density_->center_of_mass();
        rvec_sub(center_of_density, center_of_geometry_fitatoms, translation_);
    }

    if (mpi_helper() != nullptr)
    {
        mpi_helper()->broadcast(translation_, 3);
    }
    const real minimum_step_size = 0.005;
    real       step_size         = 1;
    // calculate gradient vector
    // normalize to step size
    // move along gradient vector
    // check new potential
    // if better, accept, else flip gradient vector, half step-size
    // move back along gradient vector
    // stop if step size too small
    real reference_divergence = 0;
    real new_divergence       = 0;
    while (step_size > minimum_step_size)
    {
        /* calculate the ridig body fitting force on the current atom configuration:
         *  - spread local atoms on grid,
         *  - communicate grid
         *  - get a reference KL divergence ot target density on Master node
         *  - calculate force on each atom
         *  - the total sum of the forces is the rigid body fitting force
         *  - communicate this force
         *
         */
        spread_density_(translationgroup);
        if (mpi_helper() != nullptr)
        {
            sum_reduce_simulated_density_();
        }
        if (mpi_helper() == nullptr || mpi_helper()->isMaster())
        {
            simulated_density_->add_offset(simulated_density_->grid_cell_volume()*background_density_);
            norm_simulated_        = simulated_density_->normalize();
            reference_divergence   = absolute_kl_divergence(target_density_->data(), simulated_density_->data());
        }

        KLForceCalculation_(translationgroup);
        auto force = translationgroup->local_force_sum();

        if (mpi_helper() != nullptr)
        {
            mpi_helper()->to_reals_buffer(force, 3);
            mpi_helper()->sum_reduce();
            if (mpi_helper()->isMaster())
            {
                mpi_helper()->from_reals_buffer(force, 3);
            }
        }

        if (mpi_helper() == nullptr || mpi_helper()->isMaster())
        {
            unitv(force, force);
        }

        int direction = 1;
        fprintf(input_output()->output_file(), "#\tMinimization\n#\tMinimization - New direction:\t[ %g %g %g ]\n#Minimization - Step size:\t%g\n", force[XX], force[YY], force[ZZ], step_size);
        do
        {
            RVec force_scaled;
            if (mpi_helper() == nullptr || mpi_helper()->isMaster())
            {
                svmul(step_size*direction, force, force_scaled);
                rvec_inc(translation_, force_scaled);
                fprintf(input_output()->output_file(), "#\tMinimization - Move attemted to\t[ %g %g %g ]\n", force_scaled[XX], force_scaled[YY], force_scaled[ZZ]);
            }
            if (mpi_helper() != nullptr)
            {
                mpi_helper()->broadcast(&translation_, 1);
            }
            spread_density_(translationgroup);
            if (mpi_helper() != nullptr)
            {
                sum_reduce_simulated_density_();
            }
            if (mpi_helper() == nullptr || mpi_helper()->isMaster())
            {
                simulated_density_->add_offset(simulated_density_->grid_cell_volume()*background_density_);
                norm_simulated_        = simulated_density_->normalize();
                new_divergence         = absolute_kl_divergence(target_density_->data(), simulated_density_->data());
                fprintf(input_output()->output_file(), "#\tMinimization - Divergence before move: %g , after move: %g\n", reference_divergence, new_divergence);
                if (std::abs(new_divergence) < std::abs(reference_divergence))
                {
                    reference_divergence = new_divergence;
                    step_size           *= 2.;
                    break;
                }
                else
                {
                    fprintf(input_output()->output_file(), "#\tMinimization - Moving back now.\n");
                    direction  = -1;
                    step_size /= 2.;
                }
            }

        }
        while (step_size > minimum_step_size);
    }

}



void
DensityFitting::initialize_target_density_()
{
    target_density_->add_offset(target_density_->grid_cell_volume()*background_density_);
    target_density_->normalize();
};
void
DensityFitting::initialize_buffers_()
{
    force_density_.resize(number_of_threads_);
    force_gauss_transform_.resize(number_of_threads_);
    simulated_density_buffer_.resize(number_of_threads_);
    gauss_transform_.resize(number_of_threads_);
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
    simulated_density_->resize();
    for (int thread = 0; thread < std::max(1, gmx_omp_nthreads_get(emntDefault)); ++thread)
    {
        gauss_transform_[thread]->set_sigma(sigma_);
        gauss_transform_[thread]->set_n_sigma(n_sigma_);
    }
}

void
DensityFitting::initialize_reference_density(const rvec x[], const matrix box)
{
    auto spreadgroup = wholemoleculegroup(x, box, 0);
    spread_density_(spreadgroup);
    simulated_density_->add_offset(simulated_density_->grid_cell_volume()*background_density_);
    simulated_density_->normalize();
    reference_density_ = simulated_density_->data();
}

void
DensityFitting::initialize(const matrix box, const rvec x[])
{
    initialize_target_density_();
    initialize_buffers_();
    initialize_spreading_();
    initialize_reference_density(x, box);
    if (isCenterOfMassCentered_)
    {
        translate_atoms_into_map_(x, box);
        reference_density_ = simulated_density_->data();
    }
    absolute_target_divergence_ = std::abs(absolute_kl_divergence(target_density_->data(), simulated_density_->data()));
    fprintf(input_output()->output_file(), "# KL-divergence(target|simulated): %g.\n", absolute_target_divergence_);
    requested_delta_divergence_ = absolute_target_divergence_*(float)every_nth_step_/(float)nsteps_till_fit_;
    // wholemoleculegroup(x, box, 0)->medianSort();
    out_  = open_xtc(trajectory_name_.c_str(), "w");
}

void DensityFitting::sum_reduce_simulated_density_()
{
    mpi_helper()->to_reals_buffer(simulated_density_->data().data(), simulated_density_->data().size() );
    mpi_helper()->sum_reduce();
    if (mpi_helper()->isMaster())
    {
        mpi_helper()->from_reals_buffer(simulated_density_->data().data(), simulated_density_->data().size());
    }
}

void DensityFitting::ForceKernel_KL(GroupAtom &atom, const int &thread)
{
    /* calculate for atom position x, voxel v, grid cell Volume V, force constant k:
     *
     * sum_v k*V*(x-v)/sigma exp(-(x-v)^2/2*sigma^2) * (rho_div_v)
     *
     */

    clear_rvec(*(atom.force));

    /*
     * to keep the potential translation invariant, use shifted atoms positions for force calculations.
     * Atoms are translated, such that the potential is minimal with respect to coordinate shifts
     */
    RVec x_shifted;
    rvec_add(*atom.xTransformed, translation_, x_shifted);

    /*
     * To speed up the expensive exp(-(x-v)^2/2*sigma^2) part of the force caluclation use fast gaussian gridding
     * Each thread has its own pre-allocated memory for doing the fast gaussian gridding.
     *
     * There is some double work done here, since all atoms have already been spread on a grid to calculate the total spread density.
     * However, keeping the spread grid for all atoms appears to be very memory intense. (natoms * voxelgridsize)
     */
    force_gauss_transform_[thread]->set_grid(std::move(force_density_[thread]));
    /*
     * The atoms spread weight is k_/(norm_simulated * sigma_^2), so we don't have to do that multiplication later in the force calculation loop
     * simulated_density_->grid_cell_volume()/
     */
    copy_rvec(
            force_gauss_transform_[thread]->force(x_shifted, *(atom.properties)*k_/(norm_simulated_*sigma_*sigma_), *simulated_density_),
            *(atom.force));

    force_density_[thread] = std::move(force_gauss_transform_[thread]->finish_and_return_grid());

};

void DensityFitting::plot_forces(const rvec x[], const matrix box)
{

    externalpotential::ForcePlotter plot;
    plot.start_plot_forces("forces.bild");
    real max_f    = -1;
    auto fitatoms = wholemoleculegroup(x, box, 0);
    for (auto &atom : *fitatoms)
    {
        if (norm(*(atom.force)) > max_f)
        {
            max_f = norm(*(atom.force));
        }
    }
    rvec x_translated;
    for (auto &atom : *fitatoms)
    {
        rvec_add(*atom.xTransformed, translation_, x_translated);
        plot.plot_force(x_translated, *(atom.force), 0, 1.0/max_f);
    }
    plot.stop_plot_forces();

}

void DensityFitting::KLForceCalculation_(WholeMoleculeGroup * fitatoms)
{
    #pragma omp parallel num_threads(number_of_threads_) shared(fitatoms) default(none)
    {
        int           thread             = gmx_omp_get_thread_num();
        auto         &threadlocaldensity = *simulated_density_buffer_[thread];
        threadlocaldensity.data() = (*simulated_density_).data();
        real          prefactor  = k_/(norm_simulated_*sigma_*sigma_);

        force_gauss_transform_[thread]->set_grid(std::move(force_density_[thread]));
        auto         &threadlocal_force_gauss_transform = force_gauss_transform_[thread];

        GroupIterator first_atom_for_thread = fitatoms->begin(thread, number_of_threads_);
        GroupIterator end_of_thread_atoms   = fitatoms->end(thread, number_of_threads_);

        for (auto atom = first_atom_for_thread; atom != end_of_thread_atoms; ++atom)
        {
            /* calculate for atom position x, voxel v, grid cell Volume V, force constant k:
             *
             * sum_v k*V*(x-v)/sigma exp(-(x-v)^2/2*sigma^2) * (rho_div_v)
             *
             */

            clear_rvec(*(*atom).force);

            /*
             * to keep the potential translation invariant, use shifted atoms positions for force calculations.
             * Atoms are translated, such that the potential is minimal with respect to coordinate shifts
             */
            RVec x_shifted;
            rvec_add(*(*atom).xTransformed, translation_, x_shifted);

            /*
             * To speed up the expensive exp(-(x-v)^2/2*sigma^2) part of the force caluclation use fast gaussian gridding
             * Each thread has its own pre-allocated memory for doing the fast gaussian gridding.
             *
             * There is some double work done here, since all atoms have already been spread on a grid to calculate the total spread density.
             * However, keeping the spread grid for all atoms appears to be very memory intense. (natoms * voxelgridsize)
             */
            /*
             * The atoms spread weight is k_/(norm_simulated * sigma_^2), so we don't have to do that multiplication later in the force calculation loop
             * simulated_density_->grid_cell_volume()/
             */
            copy_rvec( threadlocal_force_gauss_transform->force(x_shifted, *((*atom).properties)*prefactor, threadlocaldensity), *((*atom).force));
        }
        force_density_[thread] = std::move(threadlocal_force_gauss_transform->finish_and_return_grid());
    }
}

void DensityFitting::do_potential( const matrix box, const rvec x[], const gmx_int64_t step)
{
    auto fitatoms = wholemoleculegroup(x, box, 0);

    if (step%(100*every_nth_step_) == 0)
    {
        matrix      box_write;
        copy_mat(box, box_write);
        rvec        x_out[fitatoms->num_atoms_global()];
        int         i = 0;
        for (auto atom : *fitatoms)
        {
            rvec_add(*atom.xTransformed, translation_, x_out[i++]);
        }
        write_xtc(out_, fitatoms->num_atoms_global(), 0, 0, box_write, x_out, 1000);
    }

    /*
     * Spread all local atoms on a grid
     */
    spread_density_(fitatoms);

    /*
     * Gather the local contributions to the overall spread density on the master node.
     */
    if (mpi_helper() != nullptr)
    {
        sum_reduce_simulated_density_();
    }

    if (mpi_helper() == nullptr || mpi_helper()->isMaster())
    {
        simulated_density_->add_offset(simulated_density_->grid_cell_volume()*background_density_);
        norm_simulated_        = simulated_density_->normalize();
        auto delta_divergence  = relative_kl_divergence(target_density_->data(), simulated_density_->data(), reference_density_);
        reference_divergence_ -= delta_divergence;
        set_local_potential(k_*reference_divergence_);
        reference_density_ = simulated_density_->data();
        fprintf(input_output()->output_file(), "%8g\t%8g\t%8g\t", reference_divergence_, k_, -100*reference_divergence_/absolute_target_divergence_);
        if (delta_divergence < requested_delta_divergence_)
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

    if (mpi_helper() == nullptr || mpi_helper()->isMaster())
    {
        inv_mul(simulated_density_->data(), target_density_->data());
    }

    if (mpi_helper() != nullptr)
    {
        mpi_helper()->broadcast(simulated_density_->data().data(), simulated_density_->data().size());
        mpi_helper()->broadcast(&norm_simulated_, 1);
    }

    KLForceCalculation_(fitatoms);
    // fitatoms->parallel_loop(std::bind( &DensityFitting::ForceKernel_KL, this, std::placeholders::_1, std::placeholders::_2));
};

DensityFitting::DensityFitting() : ExternalPotential(),
                                   target_density_(std::unique_ptr<volumedata::GridReal>(new volumedata::GridReal())),
                                   simulated_density_(std::unique_ptr<volumedata::GridReal>(new volumedata::GridReal())),
                                   translation_({0, 0, 0}
                                                ), number_of_threads_ {
    std::max(1, gmx_omp_nthreads_get(emntDefault))
}
{
}

real
DensityFitting::absolute_kl_divergence(std::vector<real> &P, std::vector<real> &Q)
{
    // for numerical stability use a reference density
    if (P.size() != Q.size())
    {
        GMX_THROW(APIError("KL-Divergence calculation requires euqally sized input vectors."));
    }
    auto p           = P.begin();
    auto q           = Q.begin();
    int  size        = Q.size();
    real sum         = 0;
    #pragma omp parallel for num_threads(std::max(1, gmx_omp_nthreads_get(emntDefault))) reduction(+:sum) schedule(static, size/number_of_threads_ + 1)
    for (int i = 0; i < size; ++i)
    {
        if ((p[i] > 0) && (q[i] > 0 ))
        {
            sum += p[i] * log(q[i]);
        }
    }
    return sum;
}


real
DensityFitting::relative_kl_divergence(std::vector<real> &P, std::vector<real> &Q,
                                       std::vector<real> &Q_reference)
{
    // for numerical stability use a reference density
    if (P.size() != Q.size())
    {
        GMX_THROW(APIError("KL-Divergence calculation requires euqally sized input vectors."));
    }
    auto p           = P.begin();
    auto q           = Q.begin();
    auto q_reference = Q_reference.begin();
    int  size        = Q.size();
    real sum         = 0;
    #pragma omp parallel for num_threads(number_of_threads_) reduction(+:sum) schedule(static, size/number_of_threads_ + 1)
    for (int i = 0; i < size; ++i)
    {
        if ((p[i] > 0) && (q[i] > 0 ) && ( q_reference[i] > 0 ))
        {
            sum += p[i] * log(q[i]/(q_reference[i]));
        }
    }

    return sum;
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
    ;


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
    if (parsed_json.has("background_density"))
    {
        background_density_ = std::stof(parsed_json["background_density"]);
    }
    else
    {
        fprintf(stderr, "\n No background density given, guessing background_density = 0.01 nm^(-3). \n");
        background_density_ = 1e-2;
    }

    if (parsed_json.has("nsteps_till_fit"))
    {
        nsteps_till_fit_ = std::stoi(parsed_json.at("nsteps_till_fit"));
    }
    else
    {
        nsteps_till_fit_ = 1000000;
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
            target_input_file.write(parsed_json["target_density_used"], *target_density_);
        }
        else
        {
            target_input_file.write("target.ccp4", *target_density_);
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
    result += "#input:\tnsteps_till_fit\t\t = " + std::to_string(nsteps_till_fit_) + "\n";
    result += "#input:\tk_factor\t\t = " + std::to_string(k_factor_) + "\n";
    result += "#input:\tcenter_to_density\t\t = " + std::to_string(isCenterOfMassCentered_) + "\n";
    result += "#input:\tevery_nth_step\t\t = " + std::to_string(every_nth_step_) + "\n";
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
    close_xtc(out_);
}

std::string DensityFittingInfo::name                        = "densityfitting";
std::string DensityFittingInfo::shortDescription            = "calculate forces from difference to target density";
const int   DensityFittingInfo::numberIndexGroups           = 1;
const int   DensityFittingInfo::numberWholeMoleculeGroups   = 1;
externalpotential::ModuleCreator DensityFittingInfo::create = DensityFitting::create;

} // namespace externalpotential
} // namespace gmx
