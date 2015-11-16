#ifndef _densityfitting_h_
#define _densityfitting_h_

#include "../../externalpotential.h"
#include <string>

class DensityFittingInfo : public ExternalPotentialInfo {
    public:
        DensityFittingInfo();
        ExternalPotential* create(ExternalPotentialDataPointer data);
};

/*! \brief
 * nothing */
class DensityFitting : public ExternalPotential
{
    public:

        DensityFitting(ExternalPotentialDataPointer data);

        ~DensityFitting();

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
