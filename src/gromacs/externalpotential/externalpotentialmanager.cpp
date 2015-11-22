#include "gmxpre.h"

#include "externalpotentialmanager.h"
#include "externalpotentialregistration.h"
#include "externalpotential.h"

#include "gromacs/topology/block.h"
#include "gromacs/topology/topology.h"

#include "gromacs/gmxlib/readinp.h"

#include "gromacs/timing/wallcycle.h"

#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/exceptions.h"

#include "gromacs/fileio/oenv.h"

#include "gromacs/mdlib/sim_util.h"

#include "gromacs/math/vec.h"

#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/mdtypes/inputrec.h"

#include <stdio.h>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <iterator>

#define   block_bc(cr,   d) gmx_bcast(     sizeof(d),     & (d), (cr))
/* Probably the test for (nr) > 0 in the next macro is only needed
 * on BlueGene(/L), where IBM's MPI_Bcast will segfault after
 * dereferencing a null pointer, even when no data is to be transferred. */
#define  nblock_bc(cr, nr, d) { if ((nr) > 0) {gmx_bcast((nr)*sizeof((d)[0]), (d), (cr)); }}
#define    snew_bc(cr, d, nr) { if (!MASTER(cr)) {snew((d), (nr)); }}


namespace externalpotential {

char ** Manager::set_external_potential(int *ninp_p, t_inpfile **inp_p,
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

    std::vector<std::string> inputfilenames;
    std::vector<std::string> outputfilenames;
    std::vector<std::string> groupnames;
    char                 *** groups_p;
    snew(groups_p, 1);
    char                  ** groups = *groups_p;
    ep->number_external_potentials = 0;
    snew(ep->inputrec_data,1);
    ep->inputrec_data=nullptr;

    for (size_t method_id = 0; method_id < factory.number_methods(); method_id++)
    {

        inputfilenames=split_input_at_token_(&ninp, &inp, factory.name(method_id)+"-input","",':');
        outputfilenames=split_input_at_token_(&ninp, &inp, factory.name(method_id)+"-output","",':');
        groupnames=split_input_at_token_(&ninp, &inp, factory.name(method_id)+"-groups","system",':');

        if (groupnames.size() < inputfilenames.size() )
        {
            fprintf(stderr, "Less groups than inputfiles for external potential %s. Matching as far as possible, then keeping the last value for all following. \n", factory.name(method_id).c_str());
        }

        if (groupnames.size() < outputfilenames.size() )
        {
            fprintf(stderr, "Less outputfiles than inputfiles for external potential %s. Matching as far as possible, then keeping the last found value for all following. \n", factory.name(method_id).c_str());
        }

        for (auto && filename : inputfilenames)
        {

            srenew(ep->inputrec_data, ep->number_external_potentials+1);
            snew(ep->inputrec_data[ep->number_external_potentials],1);

            srenew(groups, ep->number_external_potentials+1);


            ep->inputrec_data[ep->number_external_potentials]->method         = method_id;
            ep->inputrec_data[ep->number_external_potentials]->inputfilename  = strdup(filename.c_str());
            if (outputfilenames.size() > 0)
            {
                ep->inputrec_data[ep->number_external_potentials]->outputfilename = strdup(outputfilenames.front().c_str());
            }else{
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
            GMX_THROW(gmx::InconsistentInputError("More index groups given than external potentials for " + factory.name(method_id) + " ."));
        }
        if (outputfilenames.size() > 1)
        {
            GMX_THROW(gmx::InconsistentInputError("More ouput filenames given than external potentials for " + factory.name(method_id) + " ."));
        }

    }

    /* STYPE macro stuff */
    *ninp_p   = ninp;
    *inp_p    = inp;

    return groups;
};

std::vector<std::string> Manager::split_input_at_token_(int * ninp_p, t_inpfile ** inp_p, std::string mdp_line_lhs, std::string default_rhs,char token){
    return split_string_at_token_(get_estr(ninp_p, inp_p, mdp_line_lhs.c_str(), default_rhs.c_str()), token);
}

void Manager::make_groups(t_ext_pot *ep, char **group_name, t_blocka *groups, char **all_group_names)
{
    int      group_index;
    std::vector<std::string>      curr_group_names;

    t_ext_pot_ir                 *ep_ir;


    for (int curr_ep = 0; curr_ep < ep->number_external_potentials; curr_ep++)
    {
        ep_ir            = ep->inputrec_data[curr_ep];
        curr_group_names = split_string_at_whitespace_(std::string(group_name[curr_ep]));
        ep_ir->number_index_groups = curr_group_names.size();

        snew(ep_ir->nat, ep_ir->number_index_groups);
        snew(ep_ir->ind, ep_ir->number_index_groups);

        for (int curr_group = 0; curr_group < ep_ir->number_index_groups; curr_group++)
        {
            group_index          = search_string_( curr_group_names[curr_group].c_str(), groups->nr, all_group_names);
            ep_ir->nat[curr_group]  = groups->index[group_index+1] - groups->index[group_index];
            snew(ep_ir->ind[curr_group],ep_ir->nat[curr_group]);
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

real Manager::add_forces( rvec            f[], tensor          vir, t_commrec      *cr, gmx_int64_t     step, real            t)
{
    real V_total = 0;
    for (auto && it : potentials_)
    {
        it->add_forces(f, vir, cr, step, t, &V_total);
    }
    return V_total;
};

std::vector<std::string> Manager::split_string_at_token_(const char * names, char token)
{
    std::stringstream        iss(names);
    std::string              single_name{};
    std::vector<std::string> result;
    while (std::getline(iss, single_name, token) )
    {
        result.push_back(single_name);
    }
    return result;
};


std::vector<std::string> Manager::split_string_at_whitespace_(std::string s)
{
    std::istringstream       iss(s);
    std::vector<std::string> result(std::istream_iterator<std::string>(iss), {});
    return result;
}

void Manager::init_external_potentials(
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

    FILE *input_file  = nullptr;
    FILE *output_file = nullptr;
    std::string basepath(ir->external_potential->basepath);

    for (int i = 0; i < ir->external_potential->number_external_potentials; ++i)
    {
        ir_data = ir->external_potential->inputrec_data[i];

        fprintf(stderr, "inputfile %s\n",ir_data->inputfilename );
        fprintf(stderr, "outputfile ''%s''\n",ir_data->outputfilename );
        if (MASTER(cr))
        {
            input_file  = gmx_ffopen((basepath + "/" + std::string(ir_data->inputfilename)).c_str(), "r");
            if(!std::string(ir_data->outputfilename).empty()){
                output_file = gmx_ffopen((basepath + "/" + std::string(ir_data->outputfilename)).c_str(), "w");
            }
        }

        potentials_.push_back(
            methods_[ir_data->method](
                    ir_data, cr, ir, fplog, x, box, mtop, bVerbose, oenv, Flags)
                );
    }

    fprintf(stderr, "\nDone initializing external potentials. Initalized %lu external potential(s).\n", potentials_.size());
    return;
};

int Manager::search_string_(const char *s, int ng, char *gn[])
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


void Manager::dd_make_local_groups(
        gmx_domdec_t *dd, t_ext_pot *external_potential
        )
{
    for (auto && it : potentials_)
    {
        it->dd_make_local_groups(dd);
    }
};

void Manager::broadcast_inputrecord_data(const t_commrec * cr, struct ext_pot * external_potential)
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
        nblock_bc(cr,curr_ir->number_index_groups, curr_ir->nat);

        snew_bc(cr, curr_ir->ind, curr_ir->number_index_groups);
        nblock_bc(cr,curr_ir->number_index_groups, curr_ir->ind);

        for (int i = 0; i < curr_ir->number_index_groups; i++)
        {
            snew_bc(cr, curr_ir->ind[i], curr_ir->nat[i]);
            nblock_bc(cr, curr_ir->nat[i], curr_ir->ind[i] );
        }
    }
}

} //namespace externalpotential
