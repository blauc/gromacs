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

#include <algorithm>

#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/futil.h"

namespace gmx
{
namespace externalpotential
{

void Manager::do_potential( t_commrec *cr, const matrix box, const rvec x[], const gmx_int64_t step)
{
    for (auto && it : potentials_)
    {
        it->do_potential(cr, box, x, step);
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

real Manager::add_forces(rvec f[], tensor vir, t_commrec *cr, gmx_int64_t step)
{
    if (potentials_.size() != 0)
    {

        int i = 0;
        for (auto && it : potentials_)
        {
            V_external_[i] = it->sum_reduce_potential_virial(cr);
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

Manager::Manager(FILE *fplog, t_inputrec *ir, gmx_mtop_t *mtop, rvec *x, matrix box, t_commrec *cr, const gmx_output_env_t *oenv, unsigned long Flags, bool bVerbose)
{
    t_ext_pot_ir * ir_data;

    FILE          *input_file  = nullptr;
    FILE          *output_file = nullptr;
    std::string    basepath(ir->external_potential->basepath);

    registerExternalPotentialModules(&modules_);

    for (int i = 0; i < ir->external_potential->number_external_potentials; ++i)
    {
        ir_data = ir->external_potential->inputrec_data[i];
        if (MASTER(cr))
        {
            if (!std::string(ir_data->inputfilename).empty())
            {
                input_file  = gmx_ffopen((basepath + "/" + std::string(ir_data->inputfilename)).c_str(), "r");
            }

            if (!std::string(ir_data->outputfilename).empty())
            {
                output_file = gmx_ffopen((basepath + "/" + std::string(ir_data->outputfilename)).c_str(), "w");
            }
        }

        try
        {
            Modules::ModuleProperties curr_module = modules_.module.at(ir_data->method);
            potentials_.push_back(curr_module.create(ir_data, cr, ir,  mtop, x, box, input_file, output_file, fplog, bVerbose, oenv, Flags, curr_module.numberIndexGroups));
        }
        catch (std::out_of_range)
        {
            GMX_THROW(gmx::InvalidInputError("Method " ""+ std::string(ir_data->method) + "" "referenced in the .mdp file was not found registered."
                                             "Most likely the gromacs version used for the mdrun is older than the one to generate the .tpr-file.\n" ));
        }

    }

    V_external_.resize(potentials_.size(), 0);

    return;
};

void Manager::throw_at_input_inconsistency_(t_commrec * cr, t_inputrec * ir, std::string input_file, std::string output_file, int current)
{
    if (current > ir->external_potential->number_external_potentials)
    {
        GMX_THROW(gmx::InconsistentInputError("Number of recognised exernal potentials does not match number of external potentials in mdp file."));
    }
    if (!gmx_fexist( input_file.c_str() ) )
    {
        GMX_THROW(gmx::FileIOError("Cannot open external input file " + input_file + "."));
    }

    FILE * outputfile_p;
    outputfile_p = fopen(output_file.c_str(), "w");
    if (outputfile_p == NULL)
    {
        GMX_THROW(gmx::FileIOError("Cannot open external output file " + output_file + " for writing."));
    }
    else
    {
        fclose(outputfile_p);
    };

    if (PAR(cr) && !DOMAINDECOMP(cr))
    {
        GMX_THROW(gmx::APIError("External potential modules only implemented for domain decomposition ."));
    }
};

void Manager::dd_make_local_groups(gmx_domdec_t *dd)
{
    for (auto && it : potentials_)
    {
        it->dd_make_local_groups(dd);
    }
};

} //namespace externalpotential
} // namespace gmx
