
#include "gmxpre.h"

#include "densityfitting.h"


std::unique_ptr<ExternalPotential> DensityFitting::create(struct ext_pot_ir *ep_ir, t_commrec * cr, t_inputrec * ir, const gmx_mtop_t* mtop, rvec x[], matrix box, FILE *input_file, FILE *output_file, FILE *fplog, gmx_bool bVerbose, const gmx_output_env_t *oenv, unsigned long Flags)
{
    return std::unique_ptr<ExternalPotential> (new DensityFitting(ep_ir, cr, ir, mtop, x, box, input_file, output_file, fplog, bVerbose, oenv, Flags));
}


void DensityFitting::do_potential(
        t_commrec      *cr,
        t_inputrec     *ir,
        matrix          box,
        rvec            x[],
        real            t,
        gmx_int64_t     step,
        gmx_wallcycle_t wcycle,
        gmx_bool        bNS)
{
    (void) cr;
    (void) ir;
    (void) box;
    (void) x;
    (void) t;
    (void) step;
    (void) wcycle;
    (void) bNS;

};
