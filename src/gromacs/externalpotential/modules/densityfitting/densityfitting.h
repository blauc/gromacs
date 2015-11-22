#ifndef _densityfitting_h_
#define _densityfitting_h_

#include "../../externalpotential.h"
#include <string>
#include <memory>
/*! \brief
 * nothing */
class DensityFitting : public ExternalPotential
{
    public:

        static std::unique_ptr<ExternalPotential> create(
            struct ext_pot_ir *ep_ir,
            t_commrec * cr,
            t_inputrec * ir,
            const gmx_mtop_t* mtop,
             rvec x[], matrix box,
             FILE *input_file,
              FILE *output_file,
             FILE *fplog,
              gmx_bool bVerbose,
             const gmx_output_env_t *oenv,
             unsigned long Flags);

        ~DensityFitting();

        void do_potential(t_commrec *cr, t_inputrec *ir, matrix box, rvec x[],
                          real t, gmx_int64_t step, gmx_wallcycle_t wcycle, gmx_bool bNS);

        void add_forces(rvec f[], tensor vir, t_commrec *cr, gmx_int64_t step,
                        real t, real *V_total);

        real potential()
        {
            return 0;
        }
private:
        DensityFitting(
            struct ext_pot_ir *ep_ir,
            t_commrec * cr,
            t_inputrec * ir,
            const gmx_mtop_t* mtop,
            rvec x[],
            matrix box,
            FILE               *input_file,
            FILE               *output_file,
            FILE               *fplog,
            gmx_bool            bVerbose,
            const gmx_output_env_t *oenv,
            unsigned long Flags);
};

#endif
