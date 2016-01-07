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

#include <stdio.h>

#include <algorithm>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

#include "gromacs/externalpotential/externalpotential.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/txtdump.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/readinp.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/sim_util.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

#define   block_bc(cr,   d) gmx_bcast(     sizeof(d),     &(d), (cr))
/* Probably the test for (nr) > 0 in the next macro is only needed
 * on BlueGene(/L), where IBM's MPI_Bcast will segfault after
 * dereferencing a null pointer, even when no data is to be transferred. */
#define  nblock_bc(cr, nr, d) { if ((nr) > 0) {gmx_bcast((nr)*sizeof((d)[0]), (d), (cr)); }}
#define    snew_bc(cr, d, nr) { if (!MASTER(cr)) {snew((d), (nr)); }}


namespace externalpotential
{

typedef std::string methodkey_t;

const std::map<methodkey_t, Registry::Info>
Registry::methods = {
    {"template", {&Template::create, "template"}}
};


std::vector<std::string> stringutils::split_input_at_token(int * ninp_p, t_inpfile ** inp_p, std::string mdp_line_lhs, std::string default_rhs, char token)
{
    return split_string_at_token(get_estr(ninp_p, inp_p, mdp_line_lhs.c_str(), default_rhs.c_str()), token);
}

int stringutils::search_string_(const char *s, int ng, char *gn[])
{
    int i;
    for (i = 0; (i < ng); i++)
    {
        if (gmx_strcasecmp(s, gn[i]) == 0)
        {
            return i;
        }
    }
    GMX_THROW(gmx::InvalidInputError("Group"+ std::string(s) + "referenced in the .mdp file was not found in the index file.\n"
                                     "Group names must match either [moleculetype] names or custom index group\n"
                                     "names, in which case you must supply an index file to the '-n' option\n"
                                     "of grompp."));
    return -1;
};

std::vector<std::string> stringutils::split_string_at_token(const char * names, char token)
{
    std::stringstream        iss(names);
    std::string              single_name {};
    std::vector<std::string> result;
    while (std::getline(iss, single_name, token) )
    {
        result.push_back(single_name);
    }
    return result;
};


std::vector<std::string> stringutils::split_string_at_whitespace(std::string s)
{
    std::istringstream       iss(s);
    std::vector<std::string> result(std::istream_iterator<std::string>(iss), {});
    return result;
}

char ** inputrecordutils::set_external_potential(int *ninp_p, t_inpfile **inp_p,
                                                 t_ext_pot * ep )
{

    /* variable declarations for STYPE/CTYPE macro */
    int         ninp;
    const char *tmp;
    t_inpfile  *inp;
    ninp = *ninp_p;
    inp  = *inp_p;

    CTYPE("Where to look for the external potential files");
    snew(ep->basepath, STRLEN);
    STYPE("external-potential-path", ep->basepath, "./");
    CTYPE("All registered external potential types. Seperate input for multiple files by " " : " " .");

    std::string              methodname;
    std::vector<std::string> inputfilenames;
    std::vector<std::string> outputfilenames;
    std::vector<std::string> groupnames;
    char                 *** groups_p;
    snew(groups_p, 1);
    char                  ** groups = *groups_p;
    ep->number_external_potentials = 0;
    snew(ep->inputrec_data, 1);
    ep->inputrec_data = nullptr;

    for (auto && method : Registry::methods)
    {
        methodname      = method.first;
        inputfilenames  = stringutils::split_input_at_token(&ninp, &inp, methodname +"-input", "", ':');
        outputfilenames = stringutils::split_input_at_token(&ninp, &inp, methodname +"-output", "", ':');
        groupnames      = stringutils::split_input_at_token(&ninp, &inp, methodname +"-groups", "system", ':');

        if (groupnames.size() < inputfilenames.size() )
        {
            fprintf(stderr, "Less groups than inputfiles for external potential %s. Matching as far as possible, then keeping the last value for all following. \n", methodname.c_str());
        }

        if (groupnames.size() < outputfilenames.size() )
        {
            fprintf(stderr, "Less outputfiles than inputfiles for external potential %s. Matching as far as possible, then keeping the last found value for all following. \n", methodname.c_str());
        }

        for (auto && filename : inputfilenames)
        {

            srenew(ep->inputrec_data, ep->number_external_potentials+1);
            snew(ep->inputrec_data[ep->number_external_potentials], 1);

            srenew(groups, ep->number_external_potentials+1);


            ep->inputrec_data[ep->number_external_potentials]->method         = strdup(methodname.c_str());
            ep->inputrec_data[ep->number_external_potentials]->inputfilename  = strdup(filename.c_str());
            if (outputfilenames.size() > 0)
            {
                ep->inputrec_data[ep->number_external_potentials]->outputfilename = strdup(outputfilenames.front().c_str());
            }
            else
            {
                // leaving outputfilename=nullptr here does not work because do_tpx needs at least one char* for reading/writing
                ep->inputrec_data[ep->number_external_potentials]->outputfilename = strdup("\0");
            }
            if (outputfilenames.size() > 1)
            {
                outputfilenames.erase(outputfilenames.begin());
            }
            snew(groups[ep->number_external_potentials], STRLEN);
            strcpy(groups[ep->number_external_potentials], groupnames.front().c_str());

            if (groupnames.size() > 1)
            {
                groupnames.erase(groupnames.begin());
            }

            ++(ep->number_external_potentials);

        }

        if (groupnames.size() > 1)
        {
            GMX_THROW(gmx::InconsistentInputError("More index groups given than external potentials for " + methodname + " ."));
        }
        if (outputfilenames.size() > 1)
        {
            GMX_THROW(gmx::InconsistentInputError("More ouput filenames given than external potentials for " + methodname + " ."));
        }

    }

    /* STYPE macro stuff */
    *ninp_p   = ninp;
    *inp_p    = inp;

    return groups;
};
void inputrecordutils::make_groups(t_ext_pot *ep, char **group_name, t_blocka *groups, char **all_group_names)
{
    int      group_index;
    std::vector<std::string>      curr_group_names;

    t_ext_pot_ir                 *ep_ir;


    for (int curr_ep = 0; curr_ep < ep->number_external_potentials; curr_ep++)
    {
        ep_ir                      = ep->inputrec_data[curr_ep];
        curr_group_names           = stringutils::split_string_at_whitespace(std::string(group_name[curr_ep]));
        ep_ir->number_index_groups = curr_group_names.size();

        snew(ep_ir->nat, ep_ir->number_index_groups);
        snew(ep_ir->ind, ep_ir->number_index_groups);

        for (int curr_group = 0; curr_group < ep_ir->number_index_groups; curr_group++)
        {
            group_index             = stringutils::search_string_( curr_group_names[curr_group].c_str(), groups->nr, all_group_names);
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

void stringutils::pr_str(FILE *fp, int indent, const char * title, char * s)
{
    pr_indent(fp, indent);
    fprintf(fp, "%-30s = %s\n", title, s);
};

void inputrecordutils::pr_externalpotential(FILE *fp, int indent, t_ext_pot * pot)
{

    t_ext_pot_ir * curr_ir;


    stringutils::pr_str(fp, indent, "external-potential-path", pot->basepath);

    for (int p = 0; p < pot->number_external_potentials; p++)
    {
        curr_ir = pot->inputrec_data[p];
        stringutils::pr_str(fp, indent, "external-potential-method", curr_ir->method);
        stringutils::pr_str(fp, indent, "external-potential-inputfile", curr_ir->inputfilename);
        stringutils::pr_str(fp, indent, "external-potential-outputfile", curr_ir->outputfilename);
    }

}



void Manager::do_potential(
        t_commrec      *cr,
        t_inputrec     *ir,
        matrix          box,
        rvec            x[],
        real            t,
        gmx_int64_t     step,
        gmx_wallcycle_t wcycle,
        gmx_bool        bNS)
{
    for (auto && it : potentials_)
    {
        it->do_potential(cr, ir, box, x, t, step, wcycle, bNS);
    }
    return;
};

std::vector<real> Manager::calculate_weights()
{
    std::vector<real> weights;

    // ensure numerical stability by substracting the largest potential averaging
    real V_max;  //< the largest of all potentials
    V_max = *std::max_element(V_external_.begin(), V_external_.end());
    real weight_normal;

    for (auto && V : V_external_)
    {
        weights.push_back(exp(V_max-V));
    }
    ;
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

        real                V_total = 0;
        int                 i;
        std::vector<real>   weights;
        std::vector<tensor> virials;
        i = 0;
        for (auto && it : potentials_)
        {
            V_external_[i] = it->sum_reduce_potential_virial(cr);
        }

        weights = calculate_weights();

        for (auto && it : potentials_)
        {
            it->add_forces(f, step, weights[i]);
            it->add_virial(vir, step, weights[i]);
            ++i;
        }
        return V_total;
    }
    ;
    return 0;
};

Manager::Manager(
        FILE                     *fplog,
        t_inputrec               *ir,
        gmx_mtop_t               *mtop,
        rvec                     *x,
        matrix                    box,
        t_commrec                *cr,
        const gmx_output_env_t   *oenv,
        unsigned long             Flags,
        gmx_bool                  bVerbose )
{
    t_ext_pot_ir * ir_data;

    FILE          *input_file  = nullptr;
    FILE          *output_file = nullptr;
    std::string    basepath(ir->external_potential->basepath);

    for (int i = 0; i < ir->external_potential->number_external_potentials; ++i)
    {
        ir_data = ir->external_potential->inputrec_data[i];

        if (MASTER(cr))
        {
            input_file  = gmx_ffopen((basepath + "/" + std::string(ir_data->inputfilename)).c_str(), "r");
            if (!std::string(ir_data->outputfilename).empty())
            {
                output_file = gmx_ffopen((basepath + "/" + std::string(ir_data->outputfilename)).c_str(), "w");
            }
        }

        try
        {
            potentials_.push_back(
                    Registry::methods.at(methodkey_t(ir_data->method)).factorymethod(
                            ir_data, cr, ir,  mtop, x, box, fplog, input_file, output_file,  bVerbose, oenv, Flags)
                    );
        }
        catch (std::out_of_range)
        {
            GMX_THROW(gmx::InvalidInputError("Method " ""+ std::string(ir_data->method) + "" "referenced in the .mdp file was not found registered.\n"
                                             "Most likely the gromacs version used for the mdrun is older than the one to generate the .tpr-file.\n" ));
        }

    }

    V_external_.resize(potentials_.size(), 0);

    fprintf(stderr, "\nDone initializing external potentials. Initalized %lu external potential(s).\n", potentials_.size());
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

void inputrecordutils::broadcast_inputrecord_data(const t_commrec * cr, struct ext_pot * external_potential)
{

    t_ext_pot_ir * curr_ir;
    block_bc(cr, *external_potential);

    snew_bc(cr, external_potential->inputrec_data, external_potential->number_external_potentials);

    for (int p = 0; p < external_potential->number_external_potentials; p++)
    {
        curr_ir = external_potential->inputrec_data[p];

        block_bc(cr, *curr_ir);

        snew_bc(cr, curr_ir->inputfilename, STRLEN);
        nblock_bc(cr, STRLEN, curr_ir->inputfilename);
        snew_bc(cr, curr_ir->outputfilename, STRLEN);
        nblock_bc(cr, STRLEN, curr_ir->outputfilename);

        snew_bc(cr, curr_ir->nat, curr_ir->number_index_groups);
        nblock_bc(cr, curr_ir->number_index_groups, curr_ir->nat);

        snew_bc(cr, curr_ir->ind, curr_ir->number_index_groups);
        nblock_bc(cr, curr_ir->number_index_groups, curr_ir->ind);

        for (int i = 0; i < curr_ir->number_index_groups; i++)
        {
            snew_bc(cr, curr_ir->ind[i], curr_ir->nat[i]);
            nblock_bc(cr, curr_ir->nat[i], curr_ir->ind[i] );
        }
    }
}

} //namespace externalpotential
