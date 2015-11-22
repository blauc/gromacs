#include "gmxpre.h"

#include "externalpotentialregistration.h"

/* include the user-provided modules */
#include "modules/densityfitting/densityfitting.h"
#include "modules/pulling/pulling.h"

#include "gromacs/utility/exceptions.h"

ExternalPotential* ExternalPotentialRegistration::create(
        size_t                 method,
        int method, struct ext_pot_ir *ep_ir,
            t_commrec * cr, t_inputrec *ir, FILE * fplog,
            std::string basepath, rvec x[], matrix box, gmx_mtop_t *mtop,
            bool bVerbose, const gmx_output_env_t * oenv, unsigned long Flags
        )
{

    FILE *input_file  = nullptr;
    FILE *output_file = nullptr;

    fprintf(stderr, "inputfile %s\n",ep_ir->inputfilename );
    fprintf(stderr, "outputfile ''%s''\n",ep_ir->outputfilename );
    if (MASTER(cr))
    {
        input_file  = gmx_ffopen((basepath + "/" + std::string(ep_ir->inputfilename)).c_str(), "r");
        if(!std::string(ep_ir->outputfilename).empty()){
            output_file = gmx_ffopen((basepath + "/" + std::string(ep_ir->outputfilename)).c_str(), "w");
        }
    }

    return method_[method]->create(ep_ir, cr, ir, fplog,input_file, output_file, x, box, mtop, bVerbose, oenv,Flags);

};

size_t ExternalPotentialRegistration::number_methods()
{
    return methods_.size();
};

std::string ExternalPotentialRegistration::name(size_t method_nr)
{

    return std::string((methods_[method_nr])->name());
    GMX_THROW(gmx::InconsistentInputError("Cannot name method."));
};

std::vector<std::string> ExternalPotentialRegistration::names()
{
    std::vector<std::string> result;
    for (auto &&method : methods_)
    {
        result.push_back(std::string(method->name()));
    }
    return result;
}
