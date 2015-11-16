#include "externalpotential.h"
#include "gromacs/topology/topology.h"
#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/sim_util.h"
#include "gromacs/mdlib/groupcoord.h"
#include "gromacs/legacyheaders/types/inputrec.h"

#include <memory>

void ExternalPotential::dd_make_local_groups( gmx_domdec_t *dd)
{
    data_->dd_make_local_groups(dd);
};



void ExternalPotentialData::dd_make_local_groups( gmx_domdec_t *dd)
{
    for ( auto && grp : atom_groups_ ) {
        grp->make_local_group_indices(dd);
    }
};

real ExternalPotential::summed_potential(t_commrec *cr)
{
    return data_->summed_potential_on_ranks(cr);
};


rvec* ExternalPotentialData::x_assembled(gmx_int64_t step, t_commrec *cr, rvec x[], matrix box, int group_index)
{
    return atom_groups_[group_index]->x_assembled(step, cr, x, box);
};


rvec* ExternalPotential::x_assembled(gmx_int64_t step, t_commrec *cr, rvec x[], matrix box, int group_index)
{
    return data_->x_assembled(step, cr, x, box, group_index);
};

std::string ExternalPotentialInfo::name(){return name_; };
std::string ExternalPotentialInfo::shortDescription(){return shortDescription_; };


ExternalPotentialGroup::ExternalPotentialGroup(
    t_inputrec * ir,
    t_commrec * cr,
    const gmx_mtop_t* mtop,
    size_t              nat,
    atom_id            *ind,
    rvec               x[],
matrix box):
     num_atoms_(nat),      ind_(ind)
{

    /* To be able to make the density fitting molecule whole: */
    /* the whole molecule is stored as x_assembled_old and serves as first reference on what defines a whole molecule */
    rvec *x_pbc            = NULL;
    rvec *x_assembled_old_ = NULL;

    if (MASTER(cr))
    {
        /* Remove pbc and prepare a whole molecule for new_t_gmx_densfit */
        snew(x_pbc, mtop->natoms); /* There ... */
        m_rveccopy(mtop->natoms, x, x_pbc);
        do_pbc_first_mtop(NULL, ir->ePBC, box, mtop, x_pbc);
        for (atom_id i = 0; i < num_atoms_; i++)
        {
            copy_rvec(x_pbc[ind[i]], x_assembled_old_[i]);
        }
        sfree(x_pbc); /* ... and back again */
    }
    if (PAR(cr))
    {
        gmx_bcast(nat * sizeof(rvec), x_assembled_old_, cr);
    }

};

ExternalPotentialData::ExternalPotentialData(
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
    unsigned long Flags
    )
    :
      input_file_(input_file),
      output_file_(output_file),
      logfile_(fplog),
      bVerbose_(bVerbose),
      Flags_(Flags),
      oenv_(oenv)
{
    for (int i = 0; i < ep_ir->number_index_groups ; i++) {
        atom_groups_.push_back(
            ExternalPotentialGroupPointer(
                new ExternalPotentialGroup(ir, cr, mtop,ep_ir->nat[i],ep_ir->ind[i],x,box)
            )
        );
    }
};


void ExternalPotentialGroup::make_local_group_indices(gmx_domdec_t *dd)
{
    dd_make_local_group_indices(dd->ga2la, num_atoms_, ind_,
                                &num_atoms_loc_, &ind_loc_, &nalloc_loc_, coll_ind_);

    /* Indicate that the group's shift vectors for this structure need to be updated
     * at the next call to communicate_group_positions, since obviously we are in a NS step */
    bUpdateShifts_ = true;
};


rvec * ExternalPotentialGroup::x_assembled(gmx_int64_t step, t_commrec * cr, rvec * x, matrix box)
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


void ExternalPotentialGroup::communicate_positions_all_to_all(
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
