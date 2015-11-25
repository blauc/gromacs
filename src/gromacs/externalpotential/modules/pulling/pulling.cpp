#include "gmxpre.h"
#include "pulling.h"


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
