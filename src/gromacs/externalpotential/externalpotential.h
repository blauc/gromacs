#ifndef _externalpotential_h_
#define _externalpotential_h_

#include "gromacs/fileio/filenm.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/real.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/mdlib/groupcoord.h"
#include "gromacs/topology/index.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/legacyheaders/network.h"

#include <string>
class ExternalPotential
{
public:
    ExternalPotential(){};

    ExternalPotential(
        const char               *input_filename,
        const char               *index_filename,
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
        gmx_bool            bVerbose    )
    {
        (void)*input_filename;
        (void)*index_filename;
        (void)*fplog;
        (void)*ir;
        (void) nfile;
        (void) fnm;
        (void)*mtop;
        (void)*x;
        (void) box;
        (void)*cr;
        (void) oenv;
        (void) Flags;
        (void) bVerbose;
    };

    // Clean up if needed
    virtual ~ExternalPotential(){};

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
        rvec f[],
        tensor vir,
        t_commrec *cr,
        gmx_int64_t step,
        real t) = 0;

    // distribute the atom coordinates over all nodes in the dd steps
    void dd_make_local_groups( gmx_domdec_t *dd) {
        dd_make_local_group_indices(dd->ga2la, nat, ind,
                &nat_loc, &ind_loc, &nalloc_loc,c_ind);
        /* Indicate that the group's shift vectors for this structure need to be updated
         * at the next call to communicate_group_positions, since obviously we are in a NS step */
        bUpdateShifts = TRUE;
    };
    real potential(){
        return potential_;
    };

    void make_group(t_atoms *atoms, t_commrec *cr)
    {
        if (MASTER(cr))
        {
            char *grpnames;
            get_index(atoms, index_filename_, 1, &nat, &ind, &grpnames);
        };

        if (PAR(cr)) {
            gmx_bcast(sizeof(nat), &nat, cr);
            if (!MASTER(cr) )
            {
                snew(ind,nat);
            }
            gmx_bcast(nat*sizeof(atom_id), &ind, cr);
        };
    };

protected:
    char *input_filename_;
    char *index_filename_;
    real potential_;
    std::string name_;
    int        nat;                 /**< Number of atoms that are influenced by the external potential. */
    int        nat_loc;             /**< Part of the atoms that are local.           */
    atom_id   *ind;                 /**< Global indices of the atoms.                */
    atom_id   *ind_loc;             /**< Local indices of the external potential atoms.             */
    int        nalloc_loc;          /**< Allocation size for ind_loc.                */
    gmx_bool   bVerbose;            /* -v flag from command line                      */
    rvec      *x;                   /**< Positions for all IMD atoms assembled on
                                          the master node.                            */
    ivec      *x_shifts;            /**< Shifts for all IMD atoms, to make
                                          molecule(s) whole.                          */
    ivec      *x_eshifts;           /**< Extra shifts since last DD step.            */
    rvec      *x_old;               /**< Old positions for all IMD atoms on master.  */
    int       *x_ind;               /**< Position of each local atom in the
                                          collective array.                           */
    atom_id   *f_ind;               /**< Force indices.                              */
    rvec      *f_loc;               /**< The forces from external potential.                     */
    int       *c_ind;               /* at which position of the whole anrs
                                     * array is a local atom?, i.e.
                                     * c_ind[0...nr_loc-1] gives the atom index
                                     * with respect to the collective
                                     * anrs[0...nr-1] array                     */
    gmx_bool bUpdateShifts;

};

#endif
