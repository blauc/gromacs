#ifndef _externalpotential_h_
#define _externalpotential_h_

#include "externalpotentialdata.h"

#include "gromacs/timing/wallcycle.h"
#include "gromacs/legacyheaders/inputrec.h"
#include "gromacs/legacyheaders/types/inputrec.h"

#include <string>

struct t_inputrec;

class ExternalPotential
{
    public:

        ExternalPotential();

        // Clean up if needed
        virtual ~ExternalPotential()=0;

        // Compute energies and forces (and possible more
        virtual void do_potential(
            t_commrec      *cr,
            t_inputrec     *ir,
            matrix          box,
            rvec            x[],
            real            t,
            gmx_int64_t     step,
            gmx_wallcycle_t wcycle,
            gmx_bool        bNS) = 0;

        virtual void add_forces(
            rvec        f[],
            tensor      vir,
            t_commrec  *cr,
            gmx_int64_t step,
            real        t,
            real       *V_total) = 0;

        virtual real potential() = 0;

        // distribute the atom indices over all nodes in the dd steps
        void dd_make_local_groups( gmx_domdec_t *dd);

        real summed_potential(t_commrec *cr);

    protected:
        ExternalPotentialData *data_;
        rvec* x_assembled(gmx_int64_t step, t_commrec *cr, rvec x[], matrix box);

};

class ExternalPotentialInfo {
    public:
        virtual ~ExternalPotentialInfo()=0;
        virtual std::string name() = 0 ;
        virtual std::string shortDescription() =0;
        virtual ExternalPotential* create(ExternalPotentialData* data) = 0;
};


#endif
