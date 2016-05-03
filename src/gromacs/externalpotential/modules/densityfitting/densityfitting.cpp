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
#include <ios>

#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/externalpotential/forceplotter.h"
#include "gromacs/externalpotential/externalpotentialIO.h"
#include "gromacs/externalpotential/forceplotter.h"
#include "gromacs/externalpotential/group.h"
#include "gromacs/externalpotential/mpi-helper.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/math/volumedata.h"
#include "gromacs/fileio/volumedataio.h"
#include "gromacs/fileio/json.h"
#include "gromacs/math/gausstransform.h"
#include "ifgt/Ifgt.h"

namespace gmx
{

namespace externalpotential
{

std::unique_ptr<ExternalPotential> DensityFitting::create()
{
    return std::unique_ptr<ExternalPotential> (new DensityFitting());
}

gmx::AtomProperties * DensityFitting::single_atom_properties(t_mdatoms * mdatoms, gmx_localtop_t * topology_loc)
{
    (void) mdatoms;
    (void) topology_loc;
    return nullptr;
}

void DensityFitting::inv_mul(std::vector<real> &to_invert, const std::vector<real> &multiplier)
{
    auto q = to_invert.begin();
    for (auto const &p : multiplier)
    {
        *q = p/(*q);
        ++q;
    }
}


void DensityFitting::spread_density_(const rvec x[])
{
    simulated_density_->zero();
    rvec x_translated;
    /* compress local atom contributions to density in expansion centers*/
    gauss_transform_->set_grid(std::move(simulated_density_));

    for (auto const &atom : *group(x, 0))
    {
        rvec_add(atom.x, translation_, x_translated);
        gauss_transform_->transform(x_translated, 1);
    }

    simulated_density_ = std::move(gauss_transform_->finish_and_return_grid());

}

void DensityFitting::translate_atoms_into_map_(const matrix box, const rvec x[])
{
    const real minimum_translation              = 1; //< Assuming that 1nm is the smallest reasonable map-shift for a good first approximation placeing the atoms within the density
    const real maximum_number_box_vector_shifts = 4; //< The search range in multiples of box-vector shifts, assuming it is unlikely that a user placed the atoms further away from the box than this number
    real       ref_mean                         = simulated_density_->mean();
    real       this_mean                        = simulated_density_->mean();
    RVec       best_translation                 = { 0, 0, 0 };
    RVec       translation_XX                   = box[XX];
    RVec       translation_YY                   = box[YY];
    RVec       translation_ZZ                   = box[ZZ];
    RVec       last_best_translation            = {0, 0, 0};

    for (real factor = maximum_number_box_vector_shifts; norm(translation_XX) > minimum_translation && norm(translation_YY) > minimum_translation && norm(translation_ZZ) > minimum_translation; factor /= 2)
    {
        for (int shift_XX = -1; shift_XX <= 1; ++shift_XX)
        {
            svmul(shift_XX*factor, box[XX], translation_XX);

            for (int shift_YY = -1; shift_YY <= 1; ++shift_YY)
            {
                svmul(shift_YY*factor, box[YY], translation_YY);

                for (int shift_ZZ = -1; shift_ZZ <= 1; ++shift_ZZ)
                {
                    svmul(shift_ZZ*factor, box[ZZ], translation_ZZ);
                    rvec_add(last_best_translation, translation_XX, translation_);
                    rvec_inc(translation_, translation_YY);
                    rvec_inc(translation_, translation_ZZ);

                    spread_density_(x);
                    this_mean = simulated_density_->mean();

                    if ( (this_mean > ref_mean))
                    {
                        ref_mean = this_mean;
                        copy_rvec(translation_, best_translation);
                        fprintf(stderr, " \n reference mean : %g ", ref_mean);
                    }
                }
            }
        }
        copy_rvec(best_translation, last_best_translation);
    }
    copy_rvec(best_translation, translation_);

}

void DensityFitting::minimize_map_potential_through_translation_(const matrix box, const rvec x[])
{
    spread_density_(x);
    const real minimum_translation = 0.01; //< Assuming that approx 1/100th of a typical bond length is the smallest reasonable map-shift

    simulated_density_->add_offset(simulated_density_->grid_cell_volume()*background_density_);
    simulated_density_->normalize();
    reference_density_  = simulated_density_->data();
    RVec last_best_translation = {0, 0, 0};
    RVec best_translation      = {0, 0, 0};

    RVec translation_XX = box[XX];
    RVec translation_YY = box[YY];
    RVec translation_ZZ = box[ZZ];

    real current_potential;

    for (real factor = 0.5; norm(translation_XX) > minimum_translation && norm(translation_YY) > minimum_translation && norm(translation_ZZ) > minimum_translation; factor /= 2)
    {
        for (int shift_XX = -1; shift_XX <= 1; ++shift_XX)
        {
            svmul(shift_XX*factor, box[XX], translation_XX);

            for (int shift_YY = -1; shift_YY <= 1; ++shift_YY)
            {
                svmul(shift_YY*factor, box[YY], translation_YY);

                for (int shift_ZZ = -1; shift_ZZ <= 1; ++shift_ZZ)
                {
                    svmul(shift_ZZ*factor, box[ZZ], translation_ZZ);
                    rvec_add(last_best_translation, translation_XX, translation_);
                    rvec_inc(translation_, translation_YY);
                    rvec_inc(translation_, translation_ZZ);

                    spread_density_(x);
                    simulated_density_->add_offset(simulated_density_->grid_cell_volume()*background_density_);
                    simulated_density_->normalize();
                    current_potential = potential(target_density_->data(), simulated_density_->data());
                    if  (current_potential <= 1e-10)
                    {
                        reference_density_  = simulated_density_->data();
                        relative_potential_ = 0;
                        copy_rvec(translation_, best_translation);
                        fprintf(stderr, "\nBest potential : %g ", current_potential);
                    }
                }
            }
        }
        copy_rvec(best_translation, last_best_translation);
    }
    copy_rvec(best_translation, translation_);
}

void DensityFitting::translation_removal(const matrix box, const rvec x[])
{

/* Perform naive grid search for map translation.
 *
 * try coordinate translations in 27 box vector combination directions, scaled by factor
 * if a coordinate translation generates a map that fits better than the previous reference map
 * make this the new reference map, store the translation as best
 *
 * after all 27 translations tried, use the previous best translation as starting point
 * try translations half as large as before
 *
 * finish, if map translations become smaler than MD-relevant length-scales
 *
 * due to normalisation and off-setting, zero density might fit better than a density in the wrong spot
 * ( ---/\ better fits ----- than /\--- )
 * as a consequence, starting from very far, with very bad initial guesses
 * will result in sending the coordinates as far as possible from the target map.
 */
    translate_atoms_into_map_(box, x);
    minimize_map_potential_through_translation_(box, x);

    volumedata::MrcFile simulated_output;
    simulated_output.write("simulated_initial.ccp4", *simulated_density_);

}

void DensityFitting::initialize(const matrix box, const rvec x[])
{

    target_density_->add_offset(target_density_->grid_cell_volume()*background_density_);
    target_density_->normalize();

    simulated_density_->copy_grid(*target_density_);
    simulated_density_->resize();

    for (int i = 0; i < std::max(1, gmx_omp_nthreads_get(emntDefault)); i++)
    {
        force_density_.push_back((std::unique_ptr<volumedata::GridReal>(new volumedata::GridReal())));
        force_density_[i]->copy_grid(*target_density_);
        force_gauss_transform_.push_back((std::unique_ptr<volumedata::FastGaussianGridding>(new volumedata::FastGaussianGridding())));
        force_gauss_transform_[i]->set_sigma(sigma_);
        force_gauss_transform_[i]->set_n_sigma(n_sigma_);
    }

    gauss_transform_->set_sigma(sigma_);
    gauss_transform_->set_n_sigma(n_sigma_);
    translation_ = {0, 0, 0};
    spread_density_(x);
    simulated_density_->add_offset(simulated_density_->grid_cell_volume()*background_density_);
    simulated_density_->normalize();
    reference_density_  = simulated_density_->data();
    relative_potential_ = 0;
    translation_removal(box, x);
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
    clear_rvec(atom.force);
    RVec x_shifted;

    rvec_add(atom.x, translation_, x_shifted);
    force_density_[thread]->zero();
    force_gauss_transform_[thread]->set_grid(std::move(force_density_[thread]));
    force_gauss_transform_[thread]->transform(x_shifted, k_*simulated_density_->grid_cell_volume()/sigma_);
    force_density_[thread] = std::move(force_gauss_transform_[thread]->finish_and_return_grid());

    RVec shift_xyz, shift_yz, shift_z;
    RVec d_x, d_y, d_z;
    rvec_sub(simulated_density_->gridpoint_coordinate(0), x_shifted, shift_xyz);
    shift_yz = shift_xyz;
    shift_z  = shift_xyz;
    d_x      = simulated_density_->unit_cell_XX();
    d_y      = simulated_density_->unit_cell_YY();
    d_z      = simulated_density_->unit_cell_ZZ();
    rvec force;
    std::vector<real>::iterator simulated_density_iterator = simulated_density_->data().begin();
    std::vector<real>::iterator force_density_iterator     = force_density_[thread]->data().begin();

    for (int i_ZZ = 0; i_ZZ < simulated_density_->extend()[ZZ]; i_ZZ++)
    {
        shift_yz = shift_z;
        for (int i_YY = 0; i_YY < simulated_density_->extend()[YY]; i_YY++)
        {
            shift_xyz = shift_yz;
            for (int i_XX = 0; i_XX < simulated_density_->extend()[XX]; i_XX++)
            {

                svmul(*force_density_iterator * *simulated_density_iterator, shift_xyz, force);

                rvec_inc(atom.force, force);
                ++simulated_density_iterator;
                ++force_density_iterator;
                rvec_inc(shift_xyz, d_x);
            }
            rvec_inc(shift_yz, d_y);

        }
        rvec_inc(shift_z, d_z);

    }


/*
    for (int i_ZZ = 0; i_ZZ < simulated_density_->extend()[ZZ]; i_ZZ++)
    {
        shift_yz = shift_z;
        for (int i_YY = 0; i_YY < simulated_density_->extend()[YY]; i_YY++)
        {
            shift_xyz = shift_yz;
            for (int i_XX = 0; i_XX < simulated_density_->extend()[XX]; i_XX++)
            {

                svmul(*simulated_density_iterator * (*force_density_iterator), shift_xyz, force);
   ++simulated_density_iterator;
   ++force_density_iterator;

                rvec_inc(atom.force, force);
                rvec_inc(shift_xyz, d_x);
            }
            rvec_inc(shift_yz, d_y);
        }
        rvec_inc(shift_z, d_z);
    }
 */
    // {
    //     clear_rvec(atom.force);
    //     RVec x_shifted;
    //     rvec_add(atom.x,translation_,x_shifted);
    //     gauss_transform_->set_sigma(sigma_);
    //     gauss_transform_->transform(x_shifted,k_*simulated_density_->grid_cell_volume());
    //
    //     // after doing the spreading shift x relative to grid for easier calculation of distance to voxels
    //     rvec_dec(x_shifted, simulated_density_->translation());
    //
    //     RVec shift_x,shift_y,shift_z;
    //     RVec d_x,d_y,d_z;
    //     d_x=simulated_density_->unit_cell_XX();
    //     d_y=simulated_density_->unit_cell_YY();
    //     d_z=simulated_density_->unit_cell_ZZ();
    //
    //     for (size_t i = 0; i < simulated_density_; i++) {
    //         /* code */
    //     }
    //     for (size_t gridpoint_index = 0; gridpoint_index < simulated_density_->num_gridpoints(); gridpoint_index++)
    //     {
    //         mvmul(impl_->unit_cell_, RVec(i[XX], i[YY], i[ZZ]), result);
    //         rvec_sub(simulated_density_->gridpoint_coordinate(gridpoint_index), x_shifted, d);
    //         svmul(simulated_density_->data()[gridpoint_index]*spread_density/sigma_, d, d);
    //
    //         rvec_add(atom.force, d, atom.force);
    //     }

};

void DensityFitting::plot_forces(const rvec x[])
{

    externalpotential::ForcePlotter plot;
    plot.start_plot_forces("forces.bild");
    real max_f = -1;
    for (auto atom : *group(x, 0))
    {
        if (norm(atom.force) > max_f)
        {
            max_f = norm(atom.force);
        }
    }
    rvec x_translated;
    for (auto &atom : *group(x, 0))
    {
        rvec_add(atom.x, translation_, x_translated);
        plot.plot_force(x_translated, atom.force, 0, 1.0/max_f);
    }
    plot.stop_plot_forces();

}

void DensityFitting::do_potential( const matrix /*box*/, const rvec x[], const gmx_int64_t /*step*/)
{
    spread_density_(x);
    simulated_density_->add_offset(simulated_density_->grid_cell_volume()*background_density_);
    simulated_density_->normalize();

    if (mpi_helper() != nullptr)
    {
        sum_reduce_simulated_density_();
    }

    if (mpi_helper() == nullptr || mpi_helper()->isMaster())
    {
        set_local_potential(potential(target_density_->data(), simulated_density_->data()));
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
    }

    fprintf(stderr, "\nCalculating density fitting forces ... ");
    group(x, 0)->parallel_loop(std::bind( &DensityFitting::ForceKernel_KL, this, std::placeholders::_1, std::placeholders::_2));
    fprintf(stderr, "done.\n");

    plot_forces(x);

};

DensityFitting::DensityFitting() : ExternalPotential(),
                                   target_density_(std::unique_ptr<volumedata::GridReal>(new volumedata::GridReal())),
                                   simulated_density_(std::unique_ptr<volumedata::GridReal>(new volumedata::GridReal())),
                                   gauss_transform_(std::unique_ptr<volumedata::FastGaussianGridding>(new volumedata::FastGaussianGridding()))
{
    ;
}


real DensityFitting::potential(std::vector<real> &P, std::vector<real> &Q)
{
    // for numerical stability use a reference density
    if (P.size() != Q.size())
    {
        GMX_THROW(APIError("KL-Divergence calculation requires euqally sized input vectors."));
    }
    potential_contribution_.resize(P.size());
    auto pot         = potential_contribution_.begin();
    auto q           = Q.begin();
    auto q_reference = reference_density_.begin();
    for (auto p : P)
    {
        if ((*q > 0) && (p > 0))
        {
            *pot = p * log(*q/(*q_reference));
        }
        ++q;
        ++q_reference;
        ++pot;
    }
    std::sort(potential_contribution_.begin(), potential_contribution_.end());

    real result = std::accumulate(potential_contribution_.begin(), potential_contribution_.end(), 0.0);

    relative_potential_ -= k_*simulated_density_->grid_cell_volume()*result;

    // this steps density becomes the reference for the next step
    reference_density_       = simulated_density_->data();

    return relative_potential_;
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

        k_                  = strtof(parsed_json["k"].c_str(), nullptr);
        sigma_              = strtof(parsed_json["sigma"].c_str(), nullptr);
        n_sigma_            = strtof(parsed_json["n_sigma"].c_str(), nullptr);
        background_density_ = strtof(parsed_json["background_density"].c_str(), nullptr);
        target_density_name = parsed_json["target_density"];

        volumedata::MrcFile         target_input_file;
        target_input_file.read(target_density_name, *target_density_);
    }
    catch (const std::ios_base::failure &e)
    {
        GMX_THROW(gmx::InvalidInputError("Reading input for external potential has failed."));
    }
}

void DensityFitting::broadcast_internal()
{
    if (mpi_helper() != nullptr)
    {
        mpi_helper()->broadcast(&k_, 1);
        mpi_helper()->broadcast(&sigma_, 1);
        mpi_helper()->broadcast(&n_sigma_, 1);
    }
}

void DensityFitting::finish()
{

}

std::string DensityFittingInfo::name                        = "densityfitting";
std::string DensityFittingInfo::shortDescription            = "calculate forces from difference to target density";
const int   DensityFittingInfo::numberIndexGroups           = 1;
externalpotential::ModuleCreator DensityFittingInfo::create = DensityFitting::create;

} // namespace externalpotential
} // namespace gmx
