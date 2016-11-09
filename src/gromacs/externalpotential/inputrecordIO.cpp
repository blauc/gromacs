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

#include "inputrecordIO.h"

#include <algorithm>
#include <list>
#include <sstream>
#include <string>
#include <vector>

#include "gromacs/externalpotential/externalpotentialmanager.h"
#include "gromacs/fileio/readinp.h"
#include "gromacs/gmxpreprocess/readir.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/topology/block.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/txtdump.h"

namespace gmx
{
namespace externalpotential
{
std::list<std::string> inputrecordutils::split_string_at_token(const char * names, char token)
{
    std::stringstream        iss(names);
    std::string              single_name {};
    std::list<std::string>   result;
    while (std::getline(iss, single_name, token) )
    {
        result.push_back(single_name);
    }
    return result;
};

char ** inputrecordutils::set_external_potential(int *ninp_p, t_inpfile **inp_p,
                                                 struct ext_pot * ep )
{

    /* variable declarations to use STYPE/CTYPE macro */
    int         ninp;
    const char *tmp;
    t_inpfile  *inp;
    ninp = *ninp_p;
    inp  = *inp_p;
    /* end declarations to use STYPE/CTYPE macro*/

    CTYPE("Where to look for the external potential files");
    snew(ep->basepath, STRLEN);
    STYPE("external-potential-path", ep->basepath, "");
    CTYPE("All registered external potential types. Seperate input for multiple files by " ":" " .");

    std::string              methodname;
    std::list<std::string>   inputfilenames;
    std::list<std::string>   outputfilenames;
    std::list<std::string>   groupnames;
    char                 *** groups_p;
    snew(groups_p, 1);
    char                  ** groups = *groups_p;
    ep->number_external_potentials = 0;
    snew(ep->inputrec_data, 1);
    ep->inputrec_data = nullptr;
    int     i_current_potential = 0;

    Modules modules;
    registerExternalPotentialModules(&modules);

    for (auto && method : modules.module)
    {
        methodname      = method.first;
        inputfilenames  = inputrecordutils::split_string_at_token(get_estr(&ninp, &inp, (methodname +"-input").c_str(), ""), ':');
        outputfilenames = inputrecordutils::split_string_at_token(get_estr(&ninp, &inp, (methodname +"-output").c_str(), ""), ':');
        groupnames      = inputrecordutils::split_string_at_token(get_estr(&ninp, &inp, (methodname +"-groups").c_str(), "system"), ':');

        size_t n_potentials = std::max({inputfilenames.size(), outputfilenames.size(), groupnames.size()});
        ep->number_external_potentials += n_potentials;

        srenew(ep->inputrec_data, ep->number_external_potentials );
        srenew(groups, ep->number_external_potentials );

        if (groupnames.size() != inputfilenames.size() || inputfilenames.size() != outputfilenames.size() )
        {
            fprintf(stderr, "\nNOTE: Number of input files, output files, and groups for external potential " "%s" " does not match.\n", methodname.c_str());
        }

        for (; i_current_potential < ep->number_external_potentials; i_current_potential++)
        {
            snew(ep->inputrec_data[i_current_potential], 1);

            ep->inputrec_data[i_current_potential]->method = strdup(methodname.c_str());

            if (!inputfilenames.empty())
            {
                ep->inputrec_data[i_current_potential]->inputfilename = strdup(inputfilenames.front().c_str());
                if (inputfilenames.size() > 1)
                {
                    inputfilenames.pop_front();
                }
            }
            else
            {
                // leaving outputfilename=nullptr here does not work because do_tpx needs at least one char* for reading/writing
                ep->inputrec_data[i_current_potential]->inputfilename = strdup("\0");
            }

            if (!outputfilenames.empty())
            {
                ep->inputrec_data[i_current_potential]->outputfilename = strdup(outputfilenames.front().c_str());
                if (outputfilenames.size() > 1)
                {
                    outputfilenames.pop_front();
                }
            }
            else
            {
                // leaving outputfilename=nullptr here does not work because do_tpx needs at least one char* for reading/writing
                ep->inputrec_data[i_current_potential]->outputfilename = strdup("\0");
            }

            snew(groups[i_current_potential], STRLEN);
            strcpy(groups[i_current_potential], groupnames.front().c_str());

            if (groupnames.size() > 1)
            {
                groupnames.pop_front();
            }
        }
    }

    /* to use STYPE/CTYPE macro  */
    *ninp_p   = ninp;
    *inp_p    = inp;
    /* end to use STYPE/CTYPE macro  */

    return groups;
};

void inputrecordutils::make_groups(t_ext_pot *ep, char **group_name, struct t_blocka *groups, char **all_group_names)
{
    int                             group_index;
    std::vector<std::string>        curr_group_names;

    t_ext_pot_ir                   *ep_ir;

    Modules modules;
    registerExternalPotentialModules(&modules);

    for (int curr_ep = 0; curr_ep < ep->number_external_potentials; curr_ep++)
    {
        ep_ir                      = ep->inputrec_data[curr_ep];
        curr_group_names           = splitString(std::string(group_name[curr_ep]));
        ep_ir->number_index_groups = curr_group_names.size();
        if (ep_ir->number_index_groups != modules.module.at(ep_ir->method).numberIndexGroups)
        {
            GMX_THROW(gmx::InconsistentInputError("Number of external potential index groups in .mdp file does not match the number of required index groups."));
        }


        snew(ep_ir->nat, ep_ir->number_index_groups);
        snew(ep_ir->ind, ep_ir->number_index_groups);

        for (int curr_group = 0; curr_group < ep_ir->number_index_groups; curr_group++)
        {
            group_index             = search_string( curr_group_names.front().c_str(), groups->nr, all_group_names);
            if (group_index == -1)
            {
                GMX_THROW(gmx::InvalidInputError("Group '"+ curr_group_names.front() + "' referenced in the .mdp file was not found in the index file."
                                                 "Group names must match either [moleculetype] names or custom index group "
                                                 "names, in which case you must supply an index file to the '-n' option "
                                                 "of grompp."));
            }

            curr_group_names.erase(curr_group_names.begin());

            ep_ir->nat[curr_group]  = groups->index[group_index+1] - groups->index[group_index];
            snew(ep_ir->ind[curr_group], ep_ir->nat[curr_group]);
            if (ep_ir->nat[curr_group] > 0)
            {
                snew(ep_ir->ind[curr_group], ep_ir->nat[curr_group]);
                for (int i = 0; i < ep_ir->nat[curr_group]; i++)
                {
                    ep_ir->ind[curr_group][i] = groups->a[groups->index[group_index]+i];
                }
            }

        }
    }
};


void inputrecordutils::pr_externalpotential(FILE *fp, int indent, t_ext_pot * pot)
{

    t_ext_pot_ir * curr_ir;
    pr_str(fp, indent, "external-potential-path", pot->basepath);

    for (int p = 0; p < pot->number_external_potentials; p++)
    {
        curr_ir = pot->inputrec_data[p];
        pr_str(fp, indent, "external-potential-method", curr_ir->method);
        pr_str(fp, indent, "external-potential-inputfile", curr_ir->inputfilename);
        pr_str(fp, indent, "external-potential-outputfile", curr_ir->outputfilename);
    }

}

void inputrecordutils::bc_externalpotential(const t_commrec * cr, struct ext_pot * ep)
{
    if (!MASTER(cr))
    {
        snew(ep->basepath, STRLEN);
    }
    MPI_Bcast(ep->basepath, STRLEN, MPI_BYTE, MASTERRANK(cr), cr->mpi_comm_mygroup);
    MPI_Bcast(&ep->number_external_potentials, sizeof(int), MPI_INT, MASTERRANK(cr), cr->mpi_comm_mygroup);
    if (!MASTER(cr))
    {
        snew(ep->inputrec_data, ep->number_external_potentials);
    }

    struct ext_pot_ir * curr_ir;
    for (int p = 0; p < ep->number_external_potentials; p++)
    {
        if (!MASTER(cr))
        {
            snew(ep->inputrec_data[p], 1);
        }
        curr_ir = ep->inputrec_data[p];

        MPI_Bcast(&(curr_ir->number_index_groups), sizeof(int), MPI_INT, MASTERRANK(cr), cr->mpi_comm_mygroup);

        if (!MASTER(cr))
        {
            snew(curr_ir->inputfilename, STRLEN);
            snew(curr_ir->outputfilename, STRLEN);
            snew(curr_ir->method, STRLEN);
            snew(curr_ir->nat, curr_ir->number_index_groups);
            snew(curr_ir->ind, curr_ir->number_index_groups);
        }

        MPI_Bcast(curr_ir->inputfilename, STRLEN, MPI_BYTE, MASTERRANK(cr), cr->mpi_comm_mygroup);
        MPI_Bcast(curr_ir->outputfilename, STRLEN, MPI_BYTE, MASTERRANK(cr), cr->mpi_comm_mygroup);
        MPI_Bcast(curr_ir->method, STRLEN, MPI_BYTE, MASTERRANK(cr), cr->mpi_comm_mygroup);
        MPI_Bcast(curr_ir->nat, sizeof(int)*curr_ir->number_index_groups, MPI_INT, MASTERRANK(cr), cr->mpi_comm_mygroup);

        for (int i = 0; i < curr_ir->number_index_groups; i++)
        {
            if (!MASTER(cr))
            {
                snew(curr_ir->ind[i], curr_ir->nat[i]);
            }
            MPI_Bcast(curr_ir->ind[i], curr_ir->nat[i], MPI_INT, MASTERRANK(cr), cr->mpi_comm_mygroup);
        }
    }

}

}
} // namespace gmx
