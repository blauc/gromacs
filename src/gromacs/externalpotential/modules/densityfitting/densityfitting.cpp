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

#include "gromacs/externalpotential/externalpotentialIO.h"
#include "gromacs/externalpotential/forceplotter.h"
#include "gromacs/externalpotential/group.h"
#include "gromacs/externalpotential/mpi-helper.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/math/volumedata.h"
#include "gromacs/fileio/volumedataio.h"
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

void DensityFitting::divide(std::vector<real> &result, const std::vector<real> &divisor)
{
    auto q = result.begin();
    for (auto p : divisor)
    {
        *q /= p;
        ++q;
    }
}

void DensityFitting::do_potential( const matrix box, const rvec x[], const gmx_int64_t step)
{
    real pot       = 0;
    real threshold = 0.001;
    /* loop over local atoms, compress map */

    gauss_transform_->set_box(box);
    difference_transform_->set_box(box);

    for (auto atom : *group(x, 0))
    {
        gauss_transform_->compress(atom.x, 1);
    }

    gauss_transform_->sum_reduce();

    if (mpi_helper()->isMaster())
    {
        gauss_transform_->expand_at_ref(*simulated_density_, *target_density_, threshold);
        pot          = -k_ * potential(target_density_->data(), simulated_density_->data());
        divide(simulated_density_->data(), target_density_->data());
        difference_transform_->compress_density(*simulated_density_);
    }
    else
    {
        pot = 0;
    }

    difference_transform_->broadcast_expansion();

    for (auto atom : *group(x, 0))
    {
        difference_transform_->do_force(atom.x, atom.force, gauss_transform_->sigma(), *simulated_density_);
    }

    set_local_potential(pot);

    (void) step;
};

DensityFitting::DensityFitting() : ExternalPotential(),
                                   target_density_(std::unique_ptr<volumedata::GridReal>()),
                                   simulated_density_(std::unique_ptr<volumedata::GridReal>()),
                                   gauss_transform_(std::unique_ptr<Ifgt>()) {};


real DensityFitting::potential(std::vector<real> &P, std::vector<real> &Q)
{
    if (P.size() != Q.size())
    {
        GMX_THROW(APIError("KL-Divergence calculation requires euqally sized input vectors."));
    }
    real result = 0;
    auto q      = Q.begin();
    for (auto p : P)
    {
        result += p * log(*q);
        ++q;
    }
    return result;
}

void DensityFitting::read_input()
{
    FILE      * inputfile       = input_output()->input_file();
    char      * line            = nullptr;
    size_t      len             = 0;
    real        sigma           = 0;
    int         expansion_order = 0;
    std::string target_density_name;
    try
    {
        getline(&line, &len, inputfile);
        k_ = strtof(line, nullptr);

        getline(&line, &len, inputfile);
        sigma = strtof(line, nullptr);

        getline(&line, &len, inputfile);
        expansion_order = strtof(line, nullptr);

        gauss_transform_->init(sigma, expansion_order);
        difference_transform_->init(sigma, expansion_order);

        getline(&line, &len, inputfile);
        target_density_name = std::string(line);
        volumedata::MrcFile         target_input_file;
        gmx::volumedata::GridReal * read_target = new gmx::volumedata::GridReal;
        target_input_file.read(target_density_name, read_target );
        target_density_ = std::unique_ptr<volumedata::GridReal>(read_target);
    }
    catch (const std::ios_base::failure &e)
    {
        GMX_THROW(gmx::InvalidInputError("Cannot read from externalpotential : template inputfile."));
    }

}

void DensityFitting::broadcast_internal()
{
    gauss_transform_->set_mpi_helper(mpi_helper());
    difference_transform_->set_mpi_helper(mpi_helper());
    if (mpi_helper() != nullptr)
    {
        mpi_helper()->broadcast(&k_, 1);
        gauss_transform_->broadcast_internal();
        difference_transform_->broadcast_internal();
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
