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

#include "externalpotential.h"
#include "externalpotentialmanager.h"
#include "modules.h"
#include "group.h"
#include "mpi-helper.h"
#include "externalpotentialIO.h"

#include <algorithm>

#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/inputrec.h"
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
    weight_normal = std::accumulate(weights.begin(), weights.end(), 0);

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

Manager::Manager(struct ext_pot *external_potential, bool bMaster, bool bParallel, MPI_Comm mpi_comm_mygroup, int masterrank)
{
    t_ext_pot_ir * ir_data;

    registerExternalPotentialModules(&modules_);

    std::shared_ptr<MpiHelper> mpihelper(new MpiHelper(mpi_comm_mygroup, masterrank, bMaster));

    for (int i = 0; i < external_potential->number_external_potentials; ++i)
    {
        ir_data = external_potential->inputrec_data[i];

        Modules::ModuleProperties curr_module;
        try
        {
            curr_module = modules_.module.at(ir_data->method);
        }
        catch (std::out_of_range)
        {
            GMX_THROW(gmx::InvalidInputError("Unknown Method : "+ std::string(ir_data->method) + ". Most likely this gromacs version is older than the one to generate the .tpr-file.\n" ));
        }

        potentials_.push_back(curr_module.create());

        if (bParallel)
        {
            potentials_.back()->set_mpi_helper(mpihelper);
        }

        if (curr_module.numberIndexGroups != ir_data->number_index_groups)
        {
            GMX_THROW(gmx::InconsistentInputError("Number of index groups found in tpr file for external potential does not match number of required index groups."));
        }

        for (int i = 0; i < curr_module.numberIndexGroups; i++)
        {
            std::unique_ptr<Group> group(new Group(ir_data->nat[i], ir_data->ind[i], bParallel));
            potentials_.back()->add_group(std::move(group));
        }

        if (bMaster)
        {
            std::shared_ptr<ExternalPotentialIO> input_output;
            input_output->set_input_file(std::string(external_potential->basepath), std::string(ir_data->inputfilename));
            input_output->set_output_file(std::string(external_potential->basepath), std::string(ir_data->outputfilename));

            potentials_.back()->set_input_output(std::move(input_output));
            potentials_.back()->read_input();
        }
    }

    V_external_.resize(potentials_.size(), 0);

    return;
};

void Manager::dd_make_local_groups(gmx_ga2la_t  * ga2la)
{
    for (auto && it : potentials_)
    {
        it->dd_make_local_groups(ga2la);
    }
};

} //namespace externalpotential
} // namespace gmx
