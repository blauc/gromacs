#include "gmxpre.h"
#include "pulling.h"

std::string PullingInfo::name_=std::string("pulling");
std::string PullingInfo::shortDescription_=std::string("do pulling");

ExternalPotential* PullingInfo::create(ExternalPotentialData *data)
{
    return (ExternalPotential*)(new Pulling(data));
}

std::string PullingInfo::name(){return name_;};
std::string PullingInfo::shortDescription(){return shortDescription_;};
