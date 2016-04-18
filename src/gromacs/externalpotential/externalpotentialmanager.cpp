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

#include "externalpotentialmanager.h"

#include <algorithm>

#include "gromacs/externalpotential/externalpotential.h"
#include "gromacs/externalpotential/externalpotentialIO.h"
#include "gromacs/externalpotential/group.h"
#include "gromacs/externalpotential/modules.h"
#include "gromacs/externalpotential/mpi-helper.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/exceptions.h"

namespace gmx
{
namespace externalpotential
{

void Manager::do_potential(const matrix box, const rvec x[], const gmx_int64_t step)
{
    for (auto && it : potentials_)
    {
        it->do_potential( box, x, step);
    }
    return;
};

std::vector<real> Manager::calculate_weights()
{
    std::vector<real> weights;

    // ensure numerical stability by substracting the largest potential from all
    real V_max;  //< the largest of all potentials
    V_max = *std::max_element(V_external_.begin(), V_external_.end());
    real weight_normal;

    for (auto && V : V_external_)
    {
        weights.push_back(exp(-(V-V_max)));
    }
    weight_normal = std::accumulate(weights.begin(), weights.end(), 0.);

    for (auto && w : weights)
    {
        w /= weight_normal;
    }

    return weights;
}

real Manager::add_forces(rvec f[], tensor vir, gmx_int64_t step)
{
    if (potentials_.size() != 0)
    {

        int i = 0;
        for (auto && it : potentials_)
        {
            V_external_[i] = it->sum_reduce_potential_virial();
        }

        std::vector<real> weights = calculate_weights();
        real              V_total = 0;

        for (auto && it : potentials_)
        {
            it->add_forces(f, step, weights[i]);
            it->add_virial(vir, step, weights[i]);
            V_total += weights[i]*V_external_[i];
            ++i;
        }
        return V_total;
    }
    return 0;
};

Manager::Manager(struct ext_pot *external_potential)
{
    registerExternalPotentialModules(&modules_);
    Modules::ModuleProperties curr_module;

    for (int i = 0; i < external_potential->number_external_potentials; ++i)
    {
        try
        {
            curr_module = modules_.module.at(external_potential->inputrec_data[i]->method);
        }
        catch (std::out_of_range)
        {
            GMX_THROW(gmx::InvalidInputError("Unknown Method : "+ std::string(external_potential->inputrec_data[i]->method) + ". Most likely this gromacs version is older than the one to generate the .tpr-file." ));
        }

        if (curr_module.numberIndexGroups != external_potential->inputrec_data[i]->number_index_groups)
        {
            GMX_THROW(gmx::InconsistentInputError("Number of index groups found in tpr file for external potential does not match number of required index groups."));
        }

        potentials_.push_back(curr_module.create());
    }

    V_external_.resize(potentials_.size(), 0);

    return;
};

void Manager::init_groups(struct ext_pot_ir ** ir_data, int n_potentials, bool bParallel)
{
    for (int i_potential = 0; i_potential < n_potentials; ++i_potential)
    {
        for (int i = 0; i < ir_data[i_potential]->number_index_groups; i++)
        {
            std::shared_ptr<Group> group(new Group(ir_data[i_potential]->nat[i], ir_data[i_potential]->ind[i], bParallel));
            potentials_[i_potential]->add_group(std::move(group));
        }
    }
}

void Manager::set_atom_properties( t_mdatoms * mdatom, gmx_localtop_t * topology_loc)
{
    for (auto && potential : potentials_)
    {
        potential->set_atom_properties(mdatom, topology_loc);
    }
}

void Manager::init_mpi(bool bMaster, MPI_Comm mpi_comm_mygroup, int masterrank)
{
    std::shared_ptr<MpiHelper> mpihelper(new MpiHelper(mpi_comm_mygroup, masterrank, bMaster));
    for (auto && potential : potentials_)
    {
        potential->set_mpi_helper(mpihelper);
    }
}

void Manager::init_input_output(struct ext_pot_ir ** ir_data, int n_potentials, const char * basepath)
{
    for (int i = 0; i < n_potentials; ++i)
    {
        std::shared_ptr<ExternalPotentialIO> input_output( new(ExternalPotentialIO));
        input_output->set_input_file(std::string(basepath), std::string(ir_data[i]->inputfilename));
        input_output->set_output_file(std::string(basepath), std::string(ir_data[i]->outputfilename));
        potentials_[i]->set_input_output(std::move(input_output));

        potentials_[i]->read_input();
    }
}

void Manager::broadcast_internal()
{
    for (auto &potential : potentials_)
    {
        potential->broadcast_internal();
    }
}

void Manager::dd_make_local_groups(gmx_ga2la_t  * ga2la)
{
    for (auto && potential : potentials_)
    {
        potential->dd_make_local_groups(ga2la);
    }
};

void Manager::finish()
{
    for (auto &potential : potentials_)
    {
        potential->finish();
    }
};

} //namespace externalpotential
} // namespace gmx
