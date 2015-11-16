#ifndef _externalpotential_h_
#define _externalpotential_h_

#include "gromacs/timing/wallcycle.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"

#include <cstdio>
#include <string>
#include <memory>
#include <vector>

struct t_commrec;
struct t_inputrec;
struct ext_pot_ir;
struct gmx_output_env_t;
struct gmx_domdec_t;
struct gmx_mtop_t;

class ExternalPotential
{
    public:
        // Compute energies and forces (and possible more
        virtual void do_potential(t_commrec *cr, t_inputrec *ir, matrix box,
            rvec x[], real t, gmx_int64_t step, gmx_wallcycle_t wcycle,
            gmx_bool bNS) = 0;

        void add_forces(rvec f[], tensor vir, t_commrec *cr, gmx_int64_t step,
            real t, real *V_total, real weight);

        real potential();

        // distribute the atom indices over all nodes in the dd steps
        void dd_make_local_groups( gmx_domdec_t *dd);

    protected:
        rvec* x_assembled(gmx_int64_t step, t_commrec *cr, rvec x[], matrix box, int group_index);
    private:
        class Impl;
        std::unique_ptr<ExternalPotential::Impl> impl_;

};

class ExternalPotentialInfo
{
    public:
        std::string name();
        std::string shortDescription();
        virtual ExternalPotential* create(ExternalPotentialDataPointer data) = 0;
        std::string name_;
        std::string shortDescription_;
};



#endif
