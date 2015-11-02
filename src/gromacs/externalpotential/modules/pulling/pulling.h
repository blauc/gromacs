#ifndef _PULLING_H
#define _PULLING_H

#include "gromacs/externalpotential/externalpotentialregistration.h"
#include <string>

class PullingInfo : public ExternalPotentialInfo {
    public:
        ExternalPotential* create(ExternalPotentialData* data);
        std::string name();
        std::string shortDescription();
    private:
        static std::string name_;
        static std::string shortDescription_;
};

class Pulling: public ExternalPotential
{
    public:
        Pulling(ExternalPotentialData *data){
            (void) data;
        }

        ~Pulling(){};

        void do_potential(
            t_commrec      *cr,
            t_inputrec     *ir,
            matrix          box,
            rvec            x[],
            real            t,
            gmx_int64_t     step,
            gmx_wallcycle_t wcycle,
            gmx_bool        bNS){
                (void)cr;
                (void)ir;
                (void)box;
                (void)x;
                (void)t;
                (void)step;
                (void)wcycle;
                (void)bNS;
            };

        void add_forces(
            rvec f[],
            tensor vir,
            t_commrec *cr,
            gmx_int64_t step,
            real t, real *V_total){
                (void) f;
                (void) vir;
                (void) cr;
                (void) step;
                (void) t;
                (void) V_total;
            } ;

        real potential(){
            return 0;
        }

};

#endif
