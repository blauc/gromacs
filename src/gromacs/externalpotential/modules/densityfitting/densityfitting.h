#ifndef _densityfitting_h_
#define _densityfitting_h_

#include "gromacs/externalpotential/externalpotential.h"
#include "gromacs/timing/wallcycle.h"

#include <string>



class DensityFittingInfo : public ExternalPotentialInfo {
    public:
        DensityFittingInfo(){};
        ~DensityFittingInfo(){};
        ExternalPotential* create(ExternalPotentialData* data);
        std::string name();
        std::string shortDescription();
    private:
        static std::string name_;
        static std::string shortDescription_;
};

/*! \brief
 * nothing */
class DensityFitting : public ExternalPotential
{
    public:

        DensityFitting(ExternalPotentialData *data);

        ~DensityFitting();

        void do_potential(
            t_commrec      *cr,
            t_inputrec     *ir,
            matrix          box,
            rvec            x[],
            real            t,
            gmx_int64_t     step,
            gmx_wallcycle_t wcycle,
            gmx_bool        bNS) ;

        void add_forces(
            rvec        f[],
            tensor      vir,
            t_commrec  *cr,
            gmx_int64_t step,
            real        t,
            real       *V_total) override;

        real potential()
        {
            return 0;
        }

};

#endif
