#include "gmxpre.h"

#include "externalpotentialutil.h"
#include "externalpotentialregistration.h"
#include "externalpotential.h"

#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/gmxlib/readinp.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/mdlib/sim_util.h"
#include "gromacs/math/vec.h"

#include <iostream>
#include <string>
#include <vector>

#include <stdio.h>

char ** ExternalPotentialUtil::set_external_potential(int *ninp_p, t_inpfile **inp_p,
                                                      t_ext_pot * ep )
{
    int         ninp;
    const char *tmp; //< for STYPE macro
    t_inpfile  *inp;
    ninp = *ninp_p;
    inp  = *inp_p;
    char *** groupnames_p;
    snew(groupnames_p, 1);
    char  ** groupnames = *groupnames_p;
    ExternalPotentialRegistration external_potentials_registry;
    CTYPE("Where to look for the external potential files");
    snew(ep->basepath, STRLEN);
    STYPE("external-potential-path", ep->basepath, "./");
    fprintf(stderr, "Looking for external potential input data and writing to: %s.\n", ep->basepath);
    CTYPE("An auto-generated list of external potential types. Input files seperated by " " : " " .");
    snew(ep->filenames, external_potentials_registry.number_methods());
    snew(ep->outputfilenames, external_potentials_registry.number_methods());
    snew(groupnames, external_potentials_registry.number_methods());
    for (size_t i = 0; i < external_potentials_registry.number_methods(); i++)
    {
        snew(ep->filenames[i], STRLEN);
        snew(ep->outputfilenames[i], STRLEN);
        snew(groupnames[i], STRLEN);
        STYPE((external_potentials_registry.name(i)+"-input").c_str(), ep->filenames[i], NULL);
        STYPE((external_potentials_registry.name(i)+"-output").c_str(), ep->outputfilenames[i], NULL);
        STYPE((external_potentials_registry.name(i)+"-groups").c_str(), groupnames[i], NULL);
        fprintf(stderr, "Reading from %s.\n", ep->filenames[i]);
        fprintf(stderr, "Writing to %s %lu \n", ep->outputfilenames[i], std::string(ep->outputfilenames[i]).size());
        fprintf(stderr, "External potential group %s.\n", groupnames[i]);

        if (!matching_number_of_inputs_(std::string(ep->filenames[i]), std::string(ep->outputfilenames[i]), std::string(groupnames[i])))
        {
            GMX_THROW(gmx::InconsistentInputError("Number of external potential input file arguments does not match number of output file arguments and/or index groups."));
        }
        ;
    }
    ;
    *ninp_p   = ninp;
    *inp_p    = inp;
    return groupnames;
};

void ExternalPotentialUtil::make_groups(t_ext_pot *ext_pot, char **group_name, t_blocka *groups, char **all_group_names)
{
    ExternalPotentialRegistration registry;
    int      group_index;
    std::vector<std::string>      group_names_per_potential;
    int n_applied_potentials = 0;
    for (size_t method = 0; method < registry.number_methods(); method++)
    {
        group_names_per_potential = split_string_at_token_(group_name[method], ':');
        for (auto current_group_name : group_names_per_potential)
        {
            ++n_applied_potentials;
            group_index            = search_string_( current_group_name.c_str(), groups->nr, all_group_names);
            srenew(ext_pot->nat, n_applied_potentials);
            ext_pot->nat[n_applied_potentials]   = groups->index[group_index+1] - groups->index[group_index];
            if (ext_pot->nat[n_applied_potentials] > 0)
            {
                fprintf(stderr, "Group '%s' with %d atoms subject to external potential %s.\n",
                        current_group_name.c_str(), ext_pot->nat[n_applied_potentials], registry.name(method).c_str());
                snew(ext_pot->ind[n_applied_potentials], ext_pot->nat[n_applied_potentials]);
                for (int i = 0; i < ext_pot->nat[n_applied_potentials]; i++)
                {
                    ext_pot->ind[n_applied_potentials][i] = groups->a[groups->index[group_index]+i];
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
    ExternalPotentialRegistration external_potential_registry;
    ir->external_potential->extpot = new t_gmx_ext_pot;
    std::vector<std::string> current_output_filename;
    std::string              complete_input_filename;
    std::string              complete_output_filename;
    int current_potential = 0;
    for (size_t method = 0; method < external_potential_registry.number_methods(); ++method)
    {
        current_output_filename = split_string_at_token_(ir->external_potential->outputfilenames[method], ':');
        for (auto && inputfile : split_string_at_token_(ir->external_potential->filenames[method], ':'))
        {
            complete_input_filename = std::string(ir->external_potential->basepath)+'/'+inputfile;
            if (current_output_filename.size() > 0)
            {
                complete_output_filename = std::string(ir->external_potential->basepath)+'/'+current_output_filename.front();
                current_output_filename.erase(current_output_filename.begin(), current_output_filename.begin());
            }
            ++current_potential;
            throw_at_input_inconsistency_(cr, ir, complete_input_filename, complete_output_filename, current_potential);
            if (MASTER(cr) && bVerbose)
            {
                fprintf(stdout, "Initializing external potential %s\n", external_potential_registry.name(method).c_str());
            }
            ExternalPotential *external_potential = external_potential_registry.init(
                        method,
                        init_data_(
                                ir->external_potential->nat[current_potential],
                                ir->external_potential->ind[current_potential],
                                cr, ir, fplog,
                                complete_input_filename,
                                complete_output_filename,
                                x, box, mtop, bVerbose, oenv, Flags)
                        );
            ir->external_potential->extpot->potentials.push_back(external_potential);
        }
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
    outputfile_p = fopen(output_file.c_str(),"w");
    if (outputfile_p==NULL){
        GMX_THROW(gmx::FileIOError("Cannot open external output file " + output_file + " for writing."));
    }else{
        fclose(outputfile_p);
    };

    if (PAR(cr) && !DOMAINDECOMP(cr))
    {
        GMX_THROW(gmx::APIError("External potential modules only implemented for domain decomposition ."));
    }
};

ExternalPotentialData* ExternalPotentialUtil::init_data_(int nat, int * ind, t_commrec * cr, t_inputrec *ir, FILE * fplog, std::string inputfilename, std::string outputfilename, rvec x[], matrix box, gmx_mtop_t *mtop, bool bVerbose, const gmx_output_env_t * oenv,
                                                         unsigned long Flags)
{
    FILE *output_file;
    FILE *input_file;
    /* To be able to make the density fitting molecule whole: */
    /* the whole molecule is stored as x_assembled_old and serves as first reference on what defines a whole molecule */
    rvec *x_pbc           = NULL;
    rvec *x_assembled_old = NULL;
    if (MASTER(cr))
    {
        /* Remove pbc and prepare a whole molecule for new_t_gmx_densfit */
        snew(x_pbc, mtop->natoms); /* There ... */
        m_rveccopy(mtop->natoms, x, x_pbc);
        do_pbc_first_mtop(NULL, ir->ePBC, box, mtop, x_pbc);
        for (int i = 0; i < nat; i++)
        {
            copy_rvec(x_pbc[ind[i]], x_assembled_old[i]);
        }
        sfree(x_pbc); /* ... and back again */
    }
    if (PAR(cr))
    {
        gmx_bcast(nat * sizeof(rvec), x_assembled_old, cr);
    }
    if (MASTER(cr))
    {
        output_file = gmx_ffopen(outputfilename.c_str(), "w");
    }
    else
    {
        output_file = NULL;
    }
    if (MASTER(cr))
    {
        input_file = gmx_ffopen(inputfilename.c_str(), "r");
    }
    else
    {
        input_file = NULL;
    }
    ExternalPotentialData* data = new ExternalPotentialData(input_file, output_file, nat, ind, x_assembled_old, fplog, bVerbose, oenv, Flags);
    return data;
};

int ExternalPotentialUtil::nr_colons_in_string_(std::string s)
{
    return std::count(s.begin(), s.end(), ':');
};

bool ExternalPotentialUtil::matching_number_of_inputs_(std::string filenames, std::string output, std::string groupnames)
{
    return (nr_colons_in_string_(filenames) == nr_colons_in_string_(output)) && ( nr_colons_in_string_(filenames) == nr_colons_in_string_(groupnames));
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
