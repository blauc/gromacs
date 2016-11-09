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
#include <numeric>

#include "gromacs/externalpotential/externalpotential.h"
#include "gromacs/externalpotential/externalpotentialIO.h"
#include "gromacs/externalpotential/atomgroups/group.h"
#include "gromacs/externalpotential/atomgroups/wholemoleculegroup.h"
#include "gromacs/externalpotential/modules.h"
#include "gromacs/externalpotential/mpi-helper.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/math/units.h"
#include "gromacs/topology/mtop_util.h"

namespace gmx
{
namespace externalpotential
{

void Manager::do_potential(const matrix box, const rvec x[], const gmx_int64_t step)
{
    update_whole_molecule_groups(x, box);
    for (auto && it : potentials_)
    {
        if (it->do_this_step(step))
        {
            it->do_potential( box, x, step);
        }
    }
    return;
};

std::vector<real> Manager::calculate_weights(real /*V_forcefield*/, real temperature)
{
    // ensure numerical stability by substracting the largest potential from all
    auto              V_max = *std::max_element(V_external_.begin(), V_external_.end());

    std::vector<real> weights(V_external_.size());
    real              kbT = BOLTZ*temperature;
    auto              boltzmannWeightAfterShift =  [kbT, V_max](real V){
            return exp(-(V-V_max)/kbT);
        };

    std::transform(std::begin(V_external_), std::end(V_external_), std::begin(weights), boltzmannWeightAfterShift);

    auto weight_normal = std::accumulate(weights.begin(), weights.end(), 0.);
    std::for_each(std::begin(weights), std::end(weights), [weight_normal](real &w){w /= weight_normal; });

/* TODO: proper option reading for exp-averaging with system potential energy
    for (int i = 0; i < V_external_.size(); ++i)
    {
        weights[i] *= exp(-(V_external_[i]-V_forcefield)/(BOLTZ*temperature));
    }
 */
    return weights;
}

real
Manager::add_forces(rvec f[], tensor vir, gmx_int64_t step, real forcefield_potential, real temperature)
{
    if (potentials_.size() != 0)
    {
        if (mpi_helper_)
        {
            mpi_helper_->to_reals_buffer(forcefield_potential);
            mpi_helper_->sum_reduce();
            forcefield_potential = mpi_helper_->from_reals_buffer();
        }
        if (step == 0)
        {
            forcefield_potential_reference_ = forcefield_potential;
        }
        auto V = V_external_.begin();
        for (auto && it : potentials_)
        {
            *V = it->sum_reduce_potential_virial();
            ++V;
        }

        std::vector<real> weights = calculate_weights(forcefield_potential-forcefield_potential_reference_, temperature);
        real              V_total = 0;

        auto              weight        = weights.begin();
        auto              V_potential   = V_external_.begin();
        for (auto && potential : potentials_)
        {
            V_total += *weight* *V_potential;
            potential->add_forces(f, step, *weight);
            potential->add_virial(vir, step, *weight);
            ++weight;
            ++V_potential;
        }
        return V_total;
    }
    return 0;
};

Manager::Manager(struct ext_pot *external_potential) : forcefield_potential_reference_ {0}
{
    registerExternalPotentialModules(&modules_);
    Modules::ModuleProperties curr_module_info;

    for (int i = 0; i < external_potential->number_external_potentials; ++i)
    {
        try
        {
            curr_module_info = modules_.module.at(external_potential->inputrec_data[i]->method);
        }
        catch (std::out_of_range)
        {
            GMX_THROW(gmx::InvalidInputError("Unknown Method : "+ std::string(external_potential->inputrec_data[i]->method) + ". Most likely this gromacs version is older than the one to generate the .tpr-file." ));
        }

        if (curr_module_info.numberIndexGroups != external_potential->inputrec_data[i]->number_index_groups)
        {
            GMX_THROW(gmx::InconsistentInputError("Number of index groups found in tpr file for external potential does not match number of required index groups."));
        }

/* TODO: set potentials name and numbers of index groups when creating ; find a better way of defining index groups (refer to index groups by name, instead of number .?)*/
        potentials_.push_back(curr_module_info.create());
    }

    V_external_.resize(potentials_.size(), 0);

    return;
};

void Manager::init_groups(struct ext_pot_ir ** ir_data, bool bParallel)
{
    /* TODO: make only the manager operator on group and call routines affecting groups, externalpotentials only read groups; enable sharing of groups across externalpotentials*/
    for (size_t i_potential = 0; i_potential < potentials_.size(); ++i_potential)
    {
        for (int i_group = 0; i_group < ir_data[i_potential]->number_index_groups; i_group++)
        {
            /* Inititalize a new atom group, then add it to an external potential */
            std::shared_ptr<Group> group(new Group(ir_data[i_potential]->nat[i_group], ir_data[i_potential]->ind[i_group], bParallel));
            potentials_[i_potential]->add_group(group);
        }
    }
}

void Manager::init_whole_molecule_groups(struct ext_pot_ir ** ir_data, const gmx_mtop_t *top_global, const int ePBC, const matrix box, const rvec x[])
{
    /* TODO: make only the manager operator on group and call routines affecting groups, externalpotentials only read groups; enable sharing of groups across externalpotentials
     * This method hsould not need to know ir_data anymore, since all information from there is in the already setup groups
     */

    for (size_t i_potential = 0; i_potential < potentials_.size(); ++i_potential)
    {
        Modules::ModuleProperties curr_module_info = modules_.module.at(ir_data[i_potential]->method);
        for (int i_group = 0; i_group < curr_module_info.numberWholeMoleculeGroups; i_group++)
        {
            if (mpi_helper_)
            {
                whole_groups_.emplace_back(new WholeMoleculeGroup(*(potentials_[i_potential]->group(x, i_group)), mpi_helper_.get(), box, 3, top_global, ePBC) );
            }
            else
            {
                whole_groups_.emplace_back(new WholeMoleculeGroup(*(potentials_[i_potential]->group(x, i_group)), nullptr, box, 3, top_global, ePBC) );
            }
            potentials_[i_potential]->add_wholemoleculegroup(whole_groups_.back().get());
        }
    }
}

void
Manager::update_whole_molecule_groups (const rvec x[], const matrix box)
{
    for (auto &group : whole_groups_)
    {
        group->update_shifts_and_reference(x, box);
    }
}

void Manager::setAtomProperties( t_mdatoms * mdatom, gmx_localtop_t * topology_loc, const gmx_mtop_t * mtop)
{
    gmx_mtop_atomlookup * topology_atomlookup = gmx_mtop_atomlookup_init(mtop);

    for (auto &potential : potentials_)
    {
        potential->set_atom_properties(mdatom, topology_loc, mtop, topology_atomlookup);
    }
}

void Manager::init_mpi(bool bMaster, MPI_Comm mpi_comm_mygroup, int masterrank)
{
    mpi_helper_ = std::unique_ptr<MpiHelper>(new MpiHelper(mpi_comm_mygroup, masterrank, bMaster));
    for (auto &potential : potentials_)
    {
        potential->set_mpi_helper(mpi_helper_.get());
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
    for (auto && group : this->whole_groups_)
    {
        group->set_indices(ga2la);
    }

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

void Manager::initialize(const matrix box, const rvec x[])
{
    for (auto &potential : potentials_)
    {
        potential->initialize(box, x);
    }
}

} //namespace externalpotential
} // namespace gmx
