
#include "gmxpre.h"

#include "densityfitting.h"

DensityFittingInfo::DensityFittingInfo(){
    name_=std::string("density-fitting");
    shortDescription_=std::string("do densfit");
};

ExternalPotential* DensityFittingInfo::create(ExternalPotentialDataPointer data)
{
    return (ExternalPotential*)(new DensityFitting(data));
}

DensityFitting::DensityFitting(ExternalPotentialDataPointer data)
{
    (void) data;
};

DensityFitting::~DensityFitting()
{
};

void DensityFitting::do_potential(
        t_commrec      *cr,
        t_inputrec     *ir,
        matrix          box,
        rvec            x[],
        real            t,
        gmx_int64_t     step,
        gmx_wallcycle_t wcycle,
        gmx_bool        bNS)
{
    (void) cr;
    (void) ir;
    (void) box;
    (void) x;
    (void) t;
    (void) step;
    (void) wcycle;
    (void) bNS;

};
