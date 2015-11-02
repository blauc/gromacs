
#include "gmxpre.h"

#include "externalpotentialdata.h"

#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/mdlib/groupcoord.h"
#include "gromacs/legacyheaders/network.h"


ExternalPotentialData::ExternalPotentialData(
    FILE               *input_file,
    FILE               *output_file,
    size_t              nat,
    atom_id            *ind,
    rvec               *x_assembled_old,
    FILE               *fplog,
    gmx_bool            bVerbose,
    const gmx_output_env_t *oenv,
    unsigned long Flags
    )
    :
      input_file_(input_file),
      output_file_(output_file),
      logfile_(fplog),
      x_assembled_old_(x_assembled_old),
      bVerbose_(bVerbose),
      oenv_(oenv),
      Flags_(Flags),
      num_atoms_(nat),
      ind_(ind)
{
};


void ExternalPotentialData::make_local_group_indices(gmx_domdec_t *dd)
{
    dd_make_local_group_indices(dd->ga2la, num_atoms_, ind_,
                                &num_atoms_loc_, &ind_loc_, &nalloc_loc_, coll_ind_);

    /* Indicate that the group's shift vectors for this structure need to be updated
     * at the next call to communicate_group_positions, since obviously we are in a NS step */
    bUpdateShifts_ = true;
};


rvec * ExternalPotentialData::x_assembled(gmx_int64_t step, t_commrec * cr, rvec * x, matrix box)
{
    if (step == last_comm_step_)
    {
        return x_assembled_;
    }
    else
    {
        communicate_positions_all_to_all(cr, x, box);
        return x_assembled_;
    }
};


void ExternalPotentialData::communicate_positions_all_to_all(
    t_commrec  *cr,
    rvec       *x,
    matrix      box    )
{
    communicate_group_positions( cr, x_assembled_, x_shifts_, extra_shifts_, bUpdateShifts_, x, num_atoms_, num_atoms_loc_, ind_loc_, coll_ind_, x_assembled_old_, box);
    bUpdateShifts_ = false;
};

double ExternalPotentialData::summed_potential_on_ranks(t_commrec * cr)
{
    gmx_sumd(1, &potential_, cr);
    return potential_;
};
