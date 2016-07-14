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

#include "template.h"

#include <cstdio>

#include <algorithm>
#include <ios>

#include "gromacs/externalpotential/externalpotentialIO.h"
#include "gromacs/externalpotential/forceplotter.h"
#include "gromacs/externalpotential/atomgroups/group.h"
#include "gromacs/externalpotential/mpi-helper.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/exceptions.h"

namespace gmx
{

namespace externalpotential
{

std::unique_ptr<ExternalPotential> Template::create()
{
    return std::unique_ptr<ExternalPotential> (new Template());
}

gmx::AtomProperties * Template::single_atom_properties(t_mdatoms * mdatoms, gmx_localtop_t * topology_loc)
{
    (void) mdatoms;
    (void) topology_loc;
    return nullptr;
}

void Template::ForceKernel_(GroupAtom &atom, const int &thread)
{
    rvec_sub(atom.x, com2_, *(atom.force));
    svmul(-k_, *(atom.force), *(atom.force));
    potential_[thread]      += k_*distance2(atom.x, com2_)/2.0;
}

void Template::do_potential( const matrix box, const rvec x[], const gmx_int64_t /*step*/)
{
    std::shared_ptr<Group> r1_local = group(x, 0);
    std::shared_ptr<Group> r2_local = group(x, 1);
    RVec                   com1;
    real                   potential = 0;

    com1[XX] = box[XX][XX]*0.25;
    com1[YY] = box[YY][YY]*0.25;
    com1[ZZ] = box[ZZ][ZZ]*0.25;

    for (auto atom : *r1_local)
    {

        rvec_sub(atom.x, com1, *(atom.force));
        svmul(-k_, *(atom.force), *(atom.force));
        potential      += k_*distance2(atom.x, com1)/2.0;
    }

    com2_[XX] = 3*com1[XX];
    com2_[YY] = 3*com1[YY];
    com2_[ZZ] = 3*com1[ZZ];

    potential_.clear();
    potential_.resize(std::max(1, gmx_omp_nthreads_get(emntDefault)));

    r2_local->parallel_loop(std::bind( &Template::ForceKernel_, this, std::placeholders::_1, std::placeholders::_2));

    for (int i = 0; i < std::max(1, gmx_omp_nthreads_get(emntDefault)); i++)
    {
        potential += potential_[i];
    }

    set_local_potential(potential);

};

bool Template::do_this_step(gmx_int64_t /*step*/)
{
    return true;
}

Template::Template() : ExternalPotential() {};

void Template::read_input()
{
    FILE * inputfile = input_output()->input_file();
    char * line      = nullptr;
    size_t len       = 0;
    try
    {
        size_t line_length = getline(&line, &len, inputfile);
        if (line_length < 1)
        {
            throw std::ios_base::failure("Failure");
        }
    }
    catch (const std::ios_base::failure &e)
    {
        GMX_THROW(gmx::InvalidInputError("Cannot read from externalpotential : template inputfile."));
    }

    k_ = strtof(line, nullptr);
}

void Template::broadcast_internal()
{
    if (mpi_helper() != nullptr)
    {
        mpi_helper()->broadcast(&k_, 1);
    }
}

void Template::initialize(const matrix /*box*/, const rvec /*x*/[]){};

void Template::finish()
{

}

std::string TemplateInfo::name                        = "template";
std::string TemplateInfo::shortDescription            = "basic external potential example";
const int   TemplateInfo::numberIndexGroups           = 2;
const int   TemplateInfo::numberWholeMoleculeGroups   = 0;
externalpotential::ModuleCreator TemplateInfo::create = Template::create;

} // namespace externalpotential
} // namespace gmx
