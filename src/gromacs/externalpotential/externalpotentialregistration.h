/*! \file
 * \brief
 * If you provided a new external potential module, this is the one and only
 * place to let gromacs know of your implementation.
 */
#ifndef _externalpotentialregistration_h_
#define _externalpotentialregistration_h_

#include "gromacs/externalpotential/externalpotential.h"

/* include the user-provided modules */
#include "gromacs/externalpotential/densityfitting/densityfitting.h"
#include "gromacs/externalpotential/pulling/pulling.h"


#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include <string>

/*! \brief
 * Convert string describing the external potential into a method enumerator.
 *
 */
class ExternalPotentialRegistration{
    public:

        ExternalPotential* init(
            size_t method,
            const char * inputfile,
            const char * indexfile,
            FILE *fplog,
            t_inputrec *ir,
            int nfile,
            const t_filenm fnm[],
            gmx_mtop_t *mtop,
            rvec *x,
            matrix box,
            t_commrec *cr,
            const output_env_t oenv,
            unsigned long Flags,
            gmx_bool bVerbose)
        {

            ExternalPotential* result(nullptr);

            if ( method == method_id("density-fitting") )
            {
                result = new DensityFitting(inputfile, indexfile, fplog, ir, nfile, fnm, mtop, x, box, cr, oenv, Flags, bVerbose);
            };

            if ( method ==  method_id("pulling") )
            {
                result = new Pulling(inputfile, indexfile, fplog, ir, nfile, fnm, mtop, x, box, cr, oenv, Flags, bVerbose);
            };

            if(result==nullptr)
            {
                GMX_THROW(gmx::APIError("Something went wrong with the initialisation of method #%d , reading data from %s and index from %s."));
            };

            return result;
        };

        size_t number_methods()
        {
            return method_names_.size();
        };

        size_t method_id(std::string methodstring){
            size_t pos;
            pos=0;
            for (auto it = method_names_.begin(); it != method_names_.end() && it->compare(methodstring)!=0; ++it)
            {
                ++pos;
            };
            return pos;
        };

        std::string name(size_t method_nr)
        {
            char error_msg[STRLEN];

            if ( method_nr < number_methods() )
            {
                return method_names_[method_nr];
            }
            else
            {
                sprintf(error_msg,"Requested name of registered external potential #%lu, but only %lu potentials are registered.", method_nr, number_methods() );
                GMX_THROW(gmx::InconsistentInputError(error_msg));
            }
            return "";
        };
    private:
        std::vector<std::string> method_names_={"density-fitting", "pulling"};
};
#endif
