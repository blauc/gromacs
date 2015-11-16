#include "gmxpre.h"

#include "externalpotentialutil.h"
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
#include "gromacs/legacyheaders/types/inputrec.h"

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


char ** ExternalPotentialUtil::set_external_potential(int *ninp_p, t_inpfile **inp_p,
                                                      t_ext_pot * ep )
{

    /* variable declarations for STYPE/CTYPE macro */
    int         ninp;
    const char *tmp;
    t_inpfile  *inp;
    ninp = *ninp_p;
    inp  = *inp_p;

    ExternalPotentialRegistration external_potentials_registry;

    CTYPE("Where to look for the external potential files");
    snew(ep->basepath, STRLEN);
    STYPE("external-potential-path", ep->basepath, "./");
    CTYPE("All registered external potential types. Seperate input for multiple files by " " : " " .");

    char                     buf[STRLEN];

    std::vector<std::string> inputfilenames;
    std::vector<std::string> outputfilenames;
    std::vector<std::string> groupnames;
    char                 *** groups_p;
    snew(groups_p, 1);
    char                  ** groups = *groups_p;
    ep->number_external_potentials = 0;
    snew(ep->inputrec_data,1);
    ep->inputrec_data=nullptr;

    for (size_t method_id = 0; method_id < external_potentials_registry.number_methods(); method_id++)
    {

        STYPE((external_potentials_registry.name(method_id)+"-input").c_str(), buf,NULL);
        //, (external_potentials_registry.name(method_id)+"-input.dat").c_str());
        inputfilenames = split_string_at_token_(buf, ':');

        STYPE((external_potentials_registry.name(method_id)+"-output").c_str(), buf,NULL);
        //, (external_potentials_registry.name(method_id)+"-output.dat").c_str());
        outputfilenames = split_string_at_token_(buf, ':');

        STYPE((external_potentials_registry.name(method_id)+"-groups").c_str(), buf, "System");
        groupnames = split_string_at_token_(buf, ':');

        if (groupnames.size() < inputfilenames.size() )
        {
            fprintf(stderr, "Less groups than inputfiles for external potential %s. Keeping the last value for all. \n", external_potentials_registry.name(method_id).c_str());
        }

        if (groupnames.size() < outputfilenames.size() )
        {
            fprintf(stderr, "Less outputfiles than inputfiles for external potential %s. Keeping the last found value for all. \n", external_potentials_registry.name(method_id).c_str());
        }

        for (auto && filename : inputfilenames)
        {

            srenew(ep->inputrec_data, ep->number_external_potentials+1);
            snew(ep->inputrec_data[ep->number_external_potentials],1);

            srenew(groups, ep->number_external_potentials+1);


            ep->inputrec_data[ep->number_external_potentials]->method         = method_id;
            ep->inputrec_data[ep->number_external_potentials]->inputfilename  = strdup(filename.c_str());
            ep->inputrec_data[ep->number_external_potentials]->outputfilename = strdup(outputfilenames.front().c_str());
            if (outputfilenames.size() > 1)
            {
                outputfilenames.erase(outputfilenames.begin());
            }
            snew(groups[ep->number_external_potentials], STRLEN);
            strcpy(groups[ep->number_external_potentials], groupnames.front().c_str());

            if (groupnames.size() > 1)
            {
                outputfilenames.erase(outputfilenames.begin());
            }

            ++(ep->number_external_potentials);

        }

        if (groupnames.size() > 1)
        {
            GMX_THROW(gmx::InconsistentInputError("More index groups given than external potentials for " + external_potentials_registry.name(method_id) + " ."));
        }
        if (outputfilenames.size() > 1)
        {
            GMX_THROW(gmx::InconsistentInputError("More ouput filenames given than external potentials for " + external_potentials_registry.name(method_id) + " ."));
        }

    }

    /* STYPE macro stuff */
    *ninp_p   = ninp;
    *inp_p    = inp;

    return groups;
};

void ExternalPotentialUtil::make_groups(t_ext_pot *ep, char **group_name, t_blocka *groups, char **all_group_names)
{
    ExternalPotentialRegistration registry;
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
                fprintf(stderr, "Group '%s' with %d atoms subject to external potential %s.\n",
                        curr_group_names[curr_group].c_str(), ep_ir->nat[curr_group], registry.name(ep_ir->method).c_str());
                snew(ep_ir->ind[curr_group], ep_ir->nat[curr_group]);
                for (int i = 0; i < ep_ir->nat[curr_group]; i++)
                {
                    ep_ir->ind[curr_group][i] = groups->a[groups->index[group_index]+i];
                }
            }

        }
    }
};

void ExternalPotentialUtil::do_external_potentials(
        t_commrec      *cr,
        t_inputrec     *ir,
        matrix          box,
        rvec            x[],
        real            t,
        gmx_int64_t     step,
        gmx_wallcycle_t wcycle,
        gmx_bool        bNS)
{
    for (auto && it :  ir->external_potential->extpot->potentials)
    {
        it->do_potential(cr, ir, box, x, t, step, wcycle, bNS);
    }
    return;
};

real ExternalPotentialUtil::add_ext_forces( t_gmx_ext_pot * extpot, rvec            f[], tensor          vir, t_commrec      *cr, gmx_int64_t     step, real            t)
{
    real V_total = 0;
    for (auto && it : extpot->potentials)
    {
        it->add_forces(f, vir, cr, step, t, &V_total);
    }
    return V_total;
};

std::vector<std::string> ExternalPotentialUtil::split_string_at_token_(char * names, char token)
{
    std::stringstream        iss(names);
    std::string              single_name;
    std::vector<std::string> result;
    while (std::getline(iss, single_name, token) )
    {
        result.push_back(single_name);
    }
    return result;
};


std::vector<std::string> ExternalPotentialUtil::split_string_at_whitespace_(std::string s)
{
    std::istringstream       iss(s);
    std::vector<std::string> result(std::istream_iterator<std::string>(iss), {});
    return result;
}

void ExternalPotentialUtil::init_external_potentials(
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
    ExternalPotentialRegistration external_potential_factory;
    ir->external_potential->extpot = new t_gmx_ext_pot;

    t_ext_pot_ir * external_potential_ir_data;

    snew(ir->external_potential->inputrec_data,ir->external_potential->number_external_potentials);

    for (int current_potential = 0; current_potential < ir->external_potential->number_external_potentials; ++current_potential)
    {
        external_potential_ir_data = ir->external_potential->inputrec_data[current_potential];
        snew(external_potential_ir_data,1);
        ExternalPotential *external_potential = external_potential_factory.init(
                    external_potential_ir_data->method,
                    init_data_(
                            external_potential_ir_data, cr, ir, fplog, ir->external_potential->basepath,
                            x, box, mtop, bVerbose, oenv, Flags)
                    );
        ir->external_potential->extpot->potentials.push_back(external_potential);
    }

    fprintf(stderr, "\nDone initializing external potentials. Initalized %lu external potential(s).\n", ir->external_potential->extpot->potentials.size());
    return;
};

int ExternalPotentialUtil::search_string_(const char *s, int ng, char *gn[])
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
void ExternalPotentialUtil::throw_at_input_inconsistency_(t_commrec * cr, t_inputrec * ir, std::string input_file, std::string output_file, int current)
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

ExternalPotentialDataPointer ExternalPotentialUtil::init_data_(
        t_ext_pot_ir *ep_ir, t_commrec * cr, t_inputrec *ir, FILE * fplog,
        std::string basepath, rvec x[], matrix box, gmx_mtop_t *mtop,
        bool bVerbose, const gmx_output_env_t * oenv, unsigned long Flags)
{
    FILE *output_file = nullptr;
    FILE *input_file  = nullptr;

    if (MASTER(cr))
    {
        output_file = gmx_ffopen((basepath + "/" + std::string(ep_ir->outputfilename)).c_str(), "w");
        input_file  = gmx_ffopen((basepath + "/" + std::string(ep_ir->inputfilename)).c_str(), "r");
    }

    return std::make_shared<ExternalPotentialData>(ep_ir, cr, ir, mtop, x, box, input_file, output_file, fplog, bVerbose, oenv, Flags);;

};

void ExternalPotentialUtil::dd_make_local_groups(
        gmx_domdec_t *dd, t_ext_pot *external_potential
        )
{
    for (auto && it : external_potential->extpot->potentials)
    {
        it->dd_make_local_groups(dd);
    }
};

void ExternalPotentialUtil::broadcast_inputrecord_data(const t_commrec * cr, struct ext_pot * external_potential)
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
