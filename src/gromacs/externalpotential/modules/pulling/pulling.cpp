#include "gmxpre.h"
#include "pulling.h"


PullingInfo::PullingInfo(){
    name_=std::string("pulling");
    shortDescription_=std::string("do pulling");
};

ExternalPotential* PullingInfo::create(ExternalPotentialDataPointer data)
{
    return (ExternalPotential*)(new Pulling(data));
}


Pulling::Pulling(ExternalPotentialDataPointer data)
{
 (void) data;
};

Pulling::~Pulling()
{
};

void Pulling::do_potential(
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
