#include "gmxpre.h"

#include "externalpotentialregistration.h"

/* include the user-provided modules */
#include "modules/densityfitting/densityfitting.h"
#include "modules/pulling/pulling.h"


ExternalPotentialRegistration::ExternalPotentialRegistration(){
    register_<DensityFittingInfo>();
    register_<PullingInfo>();
}


ExternalPotential* ExternalPotentialRegistration::init(
        size_t                 method,
        ExternalPotentialData *data
        )
{
    return methods_[method]->create(data);
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
