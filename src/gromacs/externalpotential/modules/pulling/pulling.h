#ifndef _PULLING_H
#define _PULLING_H

#include "../../externalpotential.h"
#include <string>

class PullingInfo : public ExternalPotentialInfo {
    public:
        PullingInfo();
        ExternalPotential* create(ExternalPotentialDataPointer data);
};

class Pulling: public ExternalPotential
{
    public:

        Pulling(ExternalPotentialDataPointer data);

        ~Pulling();

        void do_potential(t_commrec *cr, t_inputrec *ir, matrix box, rvec x[],
                          real t, gmx_int64_t step, gmx_wallcycle_t wcycle, gmx_bool bNS);

        void add_forces(rvec f[], tensor vir, t_commrec *cr, gmx_int64_t step,
                        real t, real *V_total);

        real potential()
        {
            return 0;
        }


};

#endif
