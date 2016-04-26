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
    for (auto p : multiplier)
    {
        *q = p/(*q);
        ++q;
    }
}


void DensityFitting::do_force_plain(const rvec x, rvec force)
{
    clear_rvec(force);
    RVec d;

    // externalpotential::ForcePlotter plot;
    // plot.start_plot_forces("plain_forces.bild");

    for (size_t i = 0; i < simulated_density_->num_gridpoints(); i++)
    {
        rvec_sub(simulated_density_->gridpoint_coordinate(i), x, d);
        svmul(1/(sigma_), d, d);
        svmul(k_*simulated_density_->grid_cell_volume()*simulated_density_->data()[i]*exp(-norm2(d)/2.0), d, d);
        // plot.plot_force(x, d, 0, 1);
        rvec_add(force, d, force);
    }
}

void DensityFitting::do_potential( const matrix /*box*/, const rvec x[], const gmx_int64_t step)
{
    real                pot       = 0;
    volumedata::MrcFile simulated_output;
    simulated_density_->zero();

    gauss_transform_->set_sigma(sigma_); //TODO: real initialisation after broadcast internal
    gauss_transform_->set_n_sigma(n_sigma_);

    /* compress local atom contributions to density in expansion centers*/
    gauss_transform_->set_grid(std::move(simulated_density_));

    for (auto atom : *group(x, 0))
    {
        gauss_transform_->transform(atom.x, 1);
    }

    simulated_density_ = std::move(gauss_transform_->finish_and_return_grid());
    simulated_density_->add_offset(simulated_density_->grid_cell_volume()*background_density_);
    simulated_density_->normalize();

    if (mpi_helper() != nullptr)
    {
        mpi_helper()->to_reals_buffer(simulated_density_->data().data(), simulated_density_->data().size() );
        mpi_helper()->sum_reduce();
        if (mpi_helper()->isMaster())
        {
            mpi_helper()->from_reals_buffer(simulated_density_->data().data(), simulated_density_->data().size());
        }
    }

    if (mpi_helper() == nullptr || mpi_helper()->isMaster())
    {
        if (step == 0)
        {
            reference_density_ = simulated_density_->data();
        }

        pot                      = potential_reference_sum_-k_ * potential(target_density_->data(), simulated_density_->data());

        potential_reference_sum_ = pot;
        reference_density_       = simulated_density_->data();

        inv_mul(simulated_density_->data(), target_density_->data());
        // simulated_output.write("divided.ccp4", *simulated_density_);
    }
    else
    {
        pot = 0;
    }

    externalpotential::ForcePlotter plot;
    plot.start_plot_forces("forces.bild");
    real max_f = -1;
    for (auto atom : *group(x, 0))
    {
        do_force_plain(atom.x, atom.force);
        if (norm(atom.force) > max_f)
        {
            max_f = norm(atom.force);
        }
    }
    for (auto atom : *group(x, 0))
    {
        plot.plot_force(atom.x, atom.force, 0, 1.0/max_f);
    }
    plot.stop_plot_forces();
    fprintf(stderr, "\n\nPOTENTIAL: %g \n\n", pot);
    set_local_potential(pot);

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
    return simulated_density_->grid_cell_volume()*result;
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

        line = (char*) malloc(fsize + 1);
        fread(line, fsize, 1, inputfile);
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
        target_density_->add_offset(target_density_->grid_cell_volume()*background_density_);
        target_density_->normalize();
        simulated_density_->copy_grid(*target_density_);
        simulated_density_->resize();
    }
    catch (const std::ios_base::failure &e)
    {
        GMX_THROW(gmx::InvalidInputError("Cannot read from externalpotential : template inputfile."));
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
