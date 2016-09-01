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
#include "gromacs/externalpotential/atomgroups/group.h"
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
    int    l_x_min = l_min[XX]-grid_index_of_spread_atom_[XX]+m_spread_;
    int    n_l_x   = l_max[XX]-l_min[XX];
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
        rvec_inc(shift_yz, d_y_offset);                                  // move in the "y-row" to the voxel where spreading starts

        int l_z = l_grid_z - grid_index_of_spread_atom_[ZZ] + m_spread_; // the grid index in the 2d spread grid around the atom
        for (int l_grid_y = l_min[YY]; l_grid_y <= l_max[YY]; ++l_grid_y)
        {
            shift_xyz = shift_yz;            // start at the beginning of an "x-column"
            rvec_inc(shift_xyz, d_x_offset); // move in the "x-row" to the voxel where spreading starts

            int l_y          = l_grid_y - grid_index_of_spread_atom_[YY]+m_spread_;
            spread_zy                     = spread_2d_[l_z][l_y];
            density_ratio_iterator        = div.zy_column_begin(l_grid_z, l_grid_y)+l_min[XX];

            spread_1d_XX = &(spread_1d_[XX][l_x_min]);

            for (int l_x = 0; l_x <= n_l_x; ++l_x)
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
    return atomicnumber2emscatteringfactor(single_atom->atomnumber);
}

void DensityFitting::inv_mul(std::vector<real> &to_invert, const std::vector<real> &multiplier)
{
#pragma omp parallel for num_threads(number_of_threads_)
    for (std::size_t i = 0; i < to_invert.size(); ++i)
    {
        to_invert[i] = multiplier[i]/(to_invert[i]);
    }
}

void DensityFitting::SpreadKernel(GroupAtom &atom, const int &thread)
{
    rvec x_translated;
    rvec_add(atom.x, translation_, x_translated);
    gauss_transform_[thread]->transform(x_translated, *(atom.properties));
}

void DensityFitting::spread_density_(const rvec x[])
{
    simulated_density_->zero();
#pragma omp parallel for num_threads(number_of_threads_)
    for (int thread = 0; thread < number_of_threads_; ++thread)
    {
        simulated_density_buffer_[thread]->zero();
        gauss_transform_[thread]->set_grid(std::move(simulated_density_buffer_[thread]));

    }

    group(x, 0)->parallel_loop(std::bind( &DensityFitting::SpreadKernel, this, std::placeholders::_1, std::placeholders::_2));

#pragma omp parallel for num_threads(number_of_threads_)
    for (int thread = 0; thread < number_of_threads_; ++thread)
    {
        simulated_density_buffer_[thread] = std::move(gauss_transform_[thread]->finish_and_return_grid());
    }

    std::vector<std::vector<real>::iterator> it_buffer(number_of_threads_);
    std::vector<real>::iterator              it_simulated_density;

    std::vector<volumedata::IVec>            minimumUsedGridIndex;
    std::vector<volumedata::IVec>            maximumUsedGridIndex;
    for (int thread = 0; thread < number_of_threads_; ++thread)
    {
        minimumUsedGridIndex.push_back(gauss_transform_[thread]->getMinimumUsedGridIndex());
        maximumUsedGridIndex.push_back(gauss_transform_[thread]->getMaximumUsedGridIndex());
    }
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

void DensityFitting::translate_atoms_into_map_(const rvec x[])
{
    RVec center_of_geometry_fitatoms = group(x, 0)->local_coordinate_sum();

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
        svmul(1.0/group(x, 0)->num_atoms_global(), center_of_geometry_fitatoms, center_of_geometry_fitatoms);
        RVec center_of_density = target_density_->center_of_mass();
        rvec_sub(center_of_density, center_of_geometry_fitatoms, translation_);
    }

    if (mpi_helper() != nullptr)
    {
        mpi_helper()->broadcast(translation_, 3);
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
    for (int thread = 0; thread < std::max(1, gmx_omp_nthreads_get(emntDefault)); thread++)
    {
        force_density_.push_back((std::unique_ptr<volumedata::GridReal>(new volumedata::GridReal())));
        force_density_[thread]->copy_grid(*target_density_);

        /* Create a gauss transform object for each density buffer,
         * because the gauss transforms will be carried out simultaneously
         */
        force_gauss_transform_.push_back((std::unique_ptr<volumedata::FastGaussianGriddingForce>(new volumedata::FastGaussianGriddingForce())));
        force_gauss_transform_[thread]->set_sigma(sigma_);
        force_gauss_transform_[thread]->set_n_sigma(n_sigma_);
        simulated_density_buffer_.push_back(std::unique_ptr<volumedata::GridReal>(new volumedata::GridReal()));
        simulated_density_buffer_[thread]->copy_grid(*target_density_);
        gauss_transform_.push_back(std::unique_ptr<volumedata::FastGaussianGridding>(new volumedata::FastGaussianGridding()));
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
DensityFitting::initialize_reference_density(const rvec x[])
{
    spread_density_(x);
    simulated_density_->add_offset(simulated_density_->grid_cell_volume()*background_density_);
    simulated_density_->normalize();
    reference_density_ = simulated_density_->data();
}

void
DensityFitting::initialize(const matrix /*box*/, const rvec x[])
{
    initialize_target_density_();
    initialize_buffers_();
    initialize_spreading_();
    if (isCenterOfMassCentered_)
    {
        translate_atoms_into_map_(x);
    }
    initialize_reference_density(x);

    out_  = open_xtc("out.xtc", "a");
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
    rvec_add(atom.x, translation_, x_shifted);

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

void DensityFitting::plot_forces(const rvec x[])
{

    externalpotential::ForcePlotter plot;
    plot.start_plot_forces("forces.bild");
    real max_f = -1;
    for (auto &atom : *group(x, 0))
    {
        if (norm(*(atom.force)) > max_f)
        {
            max_f = norm(*(atom.force));
        }
    }
    rvec x_translated;
    for (auto &atom : *group(x, 0))
    {
        rvec_add(atom.x, translation_, x_translated);
        plot.plot_force(x_translated, *(atom.force), 0, 1.0/max_f);
    }
    plot.stop_plot_forces();

}

void DensityFitting::do_potential( const matrix box, const rvec x[], const gmx_int64_t step)
{

    if (step %10 == 0)
    {
        matrix      box_write;
        copy_mat(box, box_write);
        rvec        x_out[group(x, 0)->num_atoms_global()];
        int         i = 0;
        for (auto atom : *group(x, 0))
        {
            rvec_add(atom.x, translation_, x_out[i++]);
        }
        write_xtc(out_, group(x, 0)->num_atoms_global(), 0, 0, box_write, x_out, 1000);
    }

    /*
     * Spread all local atoms on a grid
     */
    spread_density_(x);

    /*
     * Gather the local contributions to the overall spread density on the master node.
     */
    if (mpi_helper() != nullptr)
    {
        sum_reduce_simulated_density_();
    }

    if (mpi_helper() == nullptr || mpi_helper()->isMaster())
    {
        real delta_divergence;
        simulated_density_->add_offset(simulated_density_->grid_cell_volume()*background_density_);
        norm_simulated_        = simulated_density_->normalize();
        delta_divergence       = relative_kl_divergence(target_density_->data(), simulated_density_->data(), reference_density_, potential_contribution_);
        reference_divergence_ -= delta_divergence;
        set_local_potential(k_*reference_divergence_);
        reference_density_ = simulated_density_->data();
        fprintf(input_output()->output_file(), "%8g\t%8g\t", k_*reference_divergence_, k_);
        if (k_*delta_divergence < .005*every_nth_step_)
        {
            k_ *= k_factor_;
        }

        if (k_*delta_divergence > .005*every_nth_step_)
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

    group(x, 0)->parallel_loop(std::bind( &DensityFitting::ForceKernel_KL, this, std::placeholders::_1, std::placeholders::_2));
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
DensityFitting::relative_kl_divergence(std::vector<real> &P, std::vector<real> &Q,
                                       std::vector<real> &Q_reference, std::vector<real> &buffer)
{
    // for numerical stability use a reference density
    if (P.size() != Q.size())
    {
        GMX_THROW(APIError("KL-Divergence calculation requires euqally sized input vectors."));
    }
    buffer.resize(P.size());
    auto p           = P.begin();
    auto q           = Q.begin();
    auto q_reference = Q_reference.begin();
    int  size        = Q.size();
    real sum         = 0;
    #pragma omp parallel for num_threads(std::max(1, gmx_omp_nthreads_get(emntDefault))) reduction(+:sum)
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
    char      * line            = nullptr;
    std::string file_as_string;
    std::string target_density_name;
    try
    {
        fseek(inputfile, 0, SEEK_END);
        long fsize = ftell(inputfile);
        fseek(inputfile, 0, SEEK_SET);  //same as rewind(f);

        line           = (char*) malloc(fsize + 1);
        fsize          = fread(line, fsize, 1, inputfile);
        file_as_string = std::string(line);
        // TODO: implemetn json scheme for checking input consistency
        json::Object parsed_json(file_as_string);

        k_                  = std::stof(parsed_json["k"]);
        sigma_              = std::stof(parsed_json["sigma"]);
        n_sigma_            = std::stof(parsed_json["n_sigma"]);
        background_density_ = std::stof(parsed_json["background_density"]);
        target_density_name = parsed_json["target_density"];
        std::string k_factor_string;
        try
        {
            k_factor_string         = parsed_json.at("k_factor");
        }
        catch (const std::out_of_range &e)
        {
            k_factor_string = "1";
        }
        k_factor_ = std::stof(k_factor_string);
        try
        {
            isCenterOfMassCentered_ = parsed_json.at("center_to_density") == "true";
        }
        catch (const std::out_of_range &e)
        {
            isCenterOfMassCentered_ = true;
        }

        try
        {
            every_nth_step_ = std::stoi(parsed_json.at("every_nth_step"));
            k_             *= every_nth_step_;
        }
        catch (const std::out_of_range &e)
        {
            every_nth_step_ = 1;
        }

        volumedata::MrcFile         target_input_file;
        target_input_file.read(target_density_name, *target_density_);
    }
    catch (const std::ios_base::failure &e)
    {
        GMX_THROW(gmx::InvalidInputError("Reading input for external potential has failed."));
    }
    fprintf(input_output()->output_file(), "E\tk\tF_max\tweight\n");
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
