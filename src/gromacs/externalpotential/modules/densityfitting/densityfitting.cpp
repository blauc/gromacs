
#include "gmxpre.h"
#include "densityfitting.h"
#include "gromacs/timing/wallcycle.h"

std::string DensityFittingInfo::name_ = std::string("density-fitting");
std::string DensityFittingInfo::shortDescription_ = std::string("do density fitting");

std::string DensityFittingInfo::name(){return name_;};
std::string DensityFittingInfo::shortDescription(){return shortDescription_;};

ExternalPotential* DensityFittingInfo::create(ExternalPotentialData *data)
{
    return (ExternalPotential*)(new DensityFitting(data));
}

DensityFitting::DensityFitting(ExternalPotentialData *data)
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
    gmx_bool        bNS){
        (void) cr;
        (void) ir;
        (void) box;
        (void) x;
        (void) t;
        (void) step;
        (void) wcycle;
        (void) bNS;
};

void DensityFitting::add_forces(
    rvec f[],
    tensor vir,
    t_commrec *cr,
    gmx_int64_t step,
    real t,
    real *V_total){
        (void) f;
        (void) vir;
        (void) cr;
        (void) step;
        (void) t;
        (void) V_total;
    };
