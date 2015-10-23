#ifndef _externalpotentialutil_h_
#define _externalpotentialutil_h_

#include "gromacs/externalpotential/externalpotentialregistration.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/gmxlib/readinp.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/exceptions.h"
#include <iostream>
#include <string>
#include <vector>

/*! \brief
 * Structure that keeps non-inputrecord data for external potentials.
 */
typedef struct gmx_ext_pot
{
  std::vector<ExternalPotential*> potentials;
}t_gmx_ext_pot;

class ExternalPotentialUtil
{
    public:
      /*! \brief
       * Set the input stream (file names) for the external potential in grompp.
       * TODO: generalise input possiblities
       * TODO: enforce input consistency (eg. MD5-checksums)
       */
        void set_external_potential(int *ninp_p, t_inpfile **inp_p,
                                           t_ext_pot * ep)
        {
            int ninp;
            const char *tmp;//< for STYPE macro
            t_inpfile *inp;
            ninp=*ninp_p;
            inp=*inp_p;
            ExternalPotentialRegistration external_potentials_registry;

            CTYPE("Where to look for the external potential files");
            snew(ep->basepath, STRLEN);
            STYPE("external-potential-path", ep->basepath, "./");
            fprintf(stderr,"Basepath %s ", ep->basepath);
            CTYPE("An auto-generated list of external potential types. Input files seperated by "" : "" .");

            snew(ep->filenames,external_potentials_registry.number_methods());
            snew(ep->indexfilenames,external_potentials_registry.number_methods());

            for (size_t i = 0; i < external_potentials_registry.number_methods(); i++) {
                snew(ep->filenames[i],STRLEN);
                snew(ep->indexfilenames[i],STRLEN);
                STYPE((external_potentials_registry.name(i)+"-input").c_str(),ep->filenames[i],NULL);
                STYPE((external_potentials_registry.name(i)+"-index").c_str(),ep->indexfilenames[i],NULL);
            }

            *ninp_p   = ninp;
            *inp_p    = inp;

            return;

        };

        void do_external_potentials(
          t_commrec      *cr,
          t_inputrec     *ir,
          matrix          box,
          rvec            x[],
          real            t,
          gmx_int64_t     step,
          gmx_wallcycle_t wcycle,
          gmx_bool        bNS)
        {
              for (std::vector<ExternalPotential*>::iterator it =
                  ir->external_potential->extpot->potentials.begin();
                  it != ir->external_potential->extpot->potentials.end(); it++)
              {
                  (*it)->do_potential(cr, ir, box, x, t, step, wcycle, bNS);
              }
              return;
        };


        real add_ext_forces(
          t_gmx_ext_pot * extpot,
          rvec f[],
          tensor vir,
          t_commrec *cr,
          gmx_int64_t step,
          real t)
        {

              real result=0;
              for (std::vector<ExternalPotential*>::iterator it =
                  extpot->potentials.begin();
                  it != extpot->potentials.end(); it++)
              {
                   result+= (*it)->potential();
                  (*it)->add_forces(f, vir, cr, step, t);
              }
              return result;
        };

        std::vector<std::string> split_string(char * names)
        {
            std::stringstream iss(names);
            std::string single_name;
            std::vector<std::string> result;
            while( std::getline(iss,single_name,':') ){
                result.push_back(single_name);
            }
            return result;
        }
        /*! \brief
        * Read the external potential type from an input file and throw errors if
        * the potential type is not known or the given file is not readable.
        *
        * To instantiate the corresponding external potential class,
        * gromacs learns about the type of potential in the first eight characters
        * in the external potential input files
        *
        * TODO: Parse info from input-files in xml format.
        */
        void init_external_potentials(FILE *fplog, t_inputrec *ir, int nfile,
            const t_filenm fnm[], gmx_mtop_t *mtop, rvec *x, matrix box,
            t_commrec *cr, const output_env_t oenv, unsigned long Flags,
            gmx_bool bVerbose)
        {

            ExternalPotentialRegistration external_potential_registry;
            snew(ir->external_potential->extpot,1);

            std::vector<std::string> current_index_filename;

            std::string complete_input_filename;
            std::string complete_index_filename;

            for (size_t i = 0; i < external_potential_registry.number_methods(); i++)
            {

                current_index_filename=split_string(ir->external_potential->indexfilenames[i]);


                for (auto&& inputfile : split_string(ir->external_potential->filenames[i]))
                {

                    complete_input_filename=std::string(ir->external_potential->basepath)+'/'+inputfile;
                    if (current_index_filename.size()>0){
                        complete_index_filename=std::string(ir->external_potential->basepath)+'/'+current_index_filename.front();
                        current_index_filename.erase(current_index_filename.begin(),current_index_filename.begin());
                    }

                    if (! gmx_fexist( complete_input_filename.c_str() ) )
                    {
                        GMX_THROW(gmx::FileIOError("Cannot open external input file " + complete_input_filename + " ."));
                    };

                    if (! gmx_fexist( complete_index_filename.c_str() ) )
                    {
                        GMX_THROW(gmx::FileIOError("Cannot open external index file " + complete_index_filename + " ."));
                    };

                    ir->external_potential->extpot->potentials.push_back(
                        external_potential_registry.init(
                            i, complete_input_filename.c_str() ,complete_index_filename.c_str(),
                            fplog, ir, nfile, fnm, mtop, x, box, cr,
                            oenv, Flags, bVerbose)
                        );
                }

            }

            fprintf(stderr,"\nDone initializing external potentials. Initalized %lu external potential(s).\n", ir->external_potential->extpot->potentials.size());
            return;
        };

        void dd_make_local_groups(
            gmx_domdec_t *dd, t_ext_pot *external_potential
        )
        {
            for (std::vector<ExternalPotential*>::iterator it =
                external_potential->extpot->potentials.begin();
                it != external_potential->extpot->potentials.end(); it++)
            {
                (*it)->dd_make_local_groups(dd);
            }
            return;
        }

};

#endif
