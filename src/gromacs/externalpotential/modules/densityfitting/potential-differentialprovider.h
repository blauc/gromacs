#ifndef GMX_EXTERNALPOTENTIAL_POTENTIALDIFFERENTIALPROVIER_H
#define GMX_EXTERNALPOTENTIAL_POTENTIALDIFFERENTIALPROVIER_H
#include "potentialprovider.h"
#include "densitydifferentialprovider.h"
namespace gmx
{
namespace volumedata
{

class IDifferentialPotentialProvider : public IStructureDensityPotentialProvider,
                                       public IDensityDifferentialProvider
{
};

}
}

#endif /* end of include guard: GMX_EXTERNALPOTENTIAL_POTENTIALDIFFERENTIALPROVIER_H */
