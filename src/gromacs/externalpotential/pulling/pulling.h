#ifndef _PULLING_H
#define _PULLING_H

#include "gromacs/externalpotential/externalpotential.h"

class Pulling: public ExternalPotential
{
    public:
        Pulling(
            const char         *input_filename,
            const char         *index_filename,
            FILE               *fplog,
            t_inputrec         *ir,
            int                 nfile,
            const t_filenm      fnm[],
            gmx_mtop_t          *mtop,
            rvec               *x,
            matrix              box,
            t_commrec          *cr,
            const output_env_t  oenv,
            unsigned long       Flags,
            gmx_bool            bVerbose   )
        {
            (void) fplog;
            (void) input_filename;
            (void) index_filename;
            (void) ir;
            (void) nfile;
            (void) fnm;
            (void) mtop;
            (void) x;
            (void) box;
            (void) cr;
            (void) oenv;
            (void) Flags;
            (void) bVerbose;
        };

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
            real t){
                (void) f;
                (void) vir;
                (void) cr;
                (void) step;
                (void) t;
            } ;

};

#endif
