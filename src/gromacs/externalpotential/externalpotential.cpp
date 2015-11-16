#include "externalpotential.h"
#include "gromacs/topology/topology.h"
#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/sim_util.h"
#include "gromacs/mdlib/groupcoord.h"
#include "gromacs/mdtypes/inputrec.h"

#include <memory>

/*! \brief
 *  Read all generic data for experimental input modules on the master node
 */
class ExternalPotential::Impl
{
    public:
        Impl(
            struct ext_pot_ir * ep_ir,
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
);

        real reduce_potential_virial(t_commrec * cr);
        void dd_make_local_groups( gmx_domdec_t *dd);
        rvec * x_assembled(gmx_int64_t step, t_commrec *cr, rvec *x, matrix box, int group_index);
        bool do_this_step(int step);
        void add_forces(rvec f[], real w);

    private:
        gmx_int64_t             nst_apply_; /**< every how many steps to apply */
        real                    potential_; /**< the (local) contribution to the potential */
        tensor                  virial_;

        real                   *mpi_inbuf_;
        real                   *mpi_outbuf_;

        FILE                   *input_file_;
        FILE                   *output_file_;
        FILE                   *logfile_;

        gmx_bool                bVerbose_; /**< -v flag from command line                      */
        unsigned long           Flags_;
        const gmx_output_env_t *oenv_;

        class Group;
        std::vector<std::unique_ptr<ExternalPotential::Impl::Group>> atom_groups_;
};

class ExternalPotential::Impl::Group
{
    public:
        Group(t_inputrec * ir,t_commrec * cr,const gmx_mtop_t * mtop,size_t              nat,
        int            *ind,
        rvec               *x, matrix box);
        void make_local_group_indices(gmx_domdec_t *dd);
        void communicate_positions_all_to_all(t_commrec *cr, rvec *x, matrix box);
        rvec * x_assembled(gmx_int64_t step, t_commrec *cr, rvec *x, matrix box);
        void add_forces(rvec f[],real w);

    private:
        rvec         *x_assembled_;     /**< the atoms of the ind_ group, made whole, using x_assembled_old_ as reference*/
        rvec         *x_assembled_old_; /**< the whole atoms from the ind group from the previous step as reference */
        ivec         *x_shifts_;        /**< helper for making the molecule whole */
        ivec         *extra_shifts_;    /**< helper variables to assemble a whole molecule,*/

        int           num_atoms_;       /**< Number of atoms that are influenced by the external potential. */
        int          *ind_;             /**< Global indices of the atoms.                */

        int           num_atoms_loc_;   /**< Part of the atoms that are local; set by make_local_group_indices. */
        int          *ind_loc_;         /**< Local indices of the external potential atoms; set by make_local_group_indices.*/
        int           nalloc_loc_;      /**< keep memory allocated; set by make_local_group_indices.*/
        rvec         *f_loc_;           /**< the forces from external potential on the local node */
        real         *weight_loc_;      /**< Weights for the local indices */
        int          *coll_ind_;        /**< the whole atom number array */
        gmx_bool      bUpdateShifts_;   /**< perfom the coordinate shifts in neighboursearching steps */
        gmx_int64_t   last_comm_step_;  /**< the last time, coordinates were communicated all to all */

};


void ExternalPotential::dd_make_local_groups( gmx_domdec_t *dd)
{
    impl_->dd_make_local_groups(dd);
};



void ExternalPotential::Impl::dd_make_local_groups( gmx_domdec_t *dd)
{
    for ( auto && grp : atom_groups_ ) {
        grp->make_local_group_indices(dd);
    }
};

rvec* ExternalPotential::Impl::x_assembled(gmx_int64_t step, t_commrec *cr, rvec x[], matrix box, int group_index)
{
    return atom_groups_[group_index]->x_assembled(step, cr, x, box);
};


rvec* ExternalPotential::x_assembled(gmx_int64_t step, t_commrec *cr, rvec x[], matrix box, int group_index)
{
    return impl_->x_assembled(step, cr, x, box, group_index);
};

std::string ExternalPotentialInfo::name(){return name_; };
std::string ExternalPotentialInfo::shortDescription(){return shortDescription_; };


ExternalPotential::Impl::Group::Group(
    t_inputrec * ir,
    t_commrec * cr,
    const gmx_mtop_t* mtop,
    size_t              nat,
    int            *ind,
    rvec               x[],
matrix box):
     num_atoms_(nat), ind_(ind)
{

    /* To be able to make the density fitting molecule whole: */
    /* the whole molecule is stored as x_assembled_old and serves as first reference on what defines a whole molecule */
    rvec *x_pbc            = NULL;
    rvec *x_assembled_old_ = NULL;
    snew(x_assembled_old_, num_atoms_);

    if (MASTER(cr))
    {
        /* Remove pbc and prepare a whole molecule for new_t_gmx_densfit */
        snew(x_pbc, mtop->natoms); /* There ... */
        m_rveccopy(mtop->natoms, x, x_pbc);
        do_pbc_first_mtop(NULL, ir->ePBC, box, mtop, x_pbc);
        for (int i = 0; i < num_atoms_; i++)
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

ExternalPotential::Impl::Impl(
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
            std::unique_ptr<ExternalPotential::Impl::Group>(
                new ExternalPotential::Impl::Group(ir, cr, mtop,ep_ir->nat[i],ep_ir->ind[i],x,box)
            )
        );
    };
};


void ExternalPotential::Impl::Group::make_local_group_indices(gmx_domdec_t *dd)
{
    dd_make_local_group_indices(dd->ga2la, num_atoms_, ind_,
                                &num_atoms_loc_, &ind_loc_, &nalloc_loc_, coll_ind_);

    /* Indicate that the group's shift vectors for this structure need to be updated
     * at the next call to communicate_group_positions, since obviously we are in a NS step */
    bUpdateShifts_ = true;
};


rvec * ExternalPotential::Impl::Group::x_assembled(gmx_int64_t step, t_commrec * cr, rvec * x, matrix box)
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

void ExternalPotential::Impl::Group::add_forces(rvec f[],real w)
{
        for (int l = 0; l < num_atoms_loc_; l++)
        {
            /* Get the right index of the local force, since typically not all local
            * atoms are subject to density fitting forces */
            /* Add to local force */
            rvec_inc(f[ind_loc_[l]], f_loc_[l]);
        }
};

void ExternalPotential::Impl::Group::communicate_positions_all_to_all(
    t_commrec  *cr,
    rvec       *x,
    matrix      box    )
{
    communicate_group_positions( cr, x_assembled_, x_shifts_, extra_shifts_, bUpdateShifts_, x, num_atoms_, num_atoms_loc_, ind_loc_, coll_ind_, x_assembled_old_, box);
    bUpdateShifts_ = false;
};

void ExternalPotential::Impl::add_forces(rvec f[], real total_weight)
{
    for ( auto && group : atom_groups_)
    {
        group->add_forces(f,total_weight);
    }
}

real ExternalPotential::Impl::reduce_potential_virial(t_commrec * cr)
{

#ifdef GMX_MPI
    MPI_Reduce(mpi_inbuf_, mpi_outbuf_, mpi_buf_size_, GMX_MPI_REAL, MPI_SUM, MASTERRANK(cr), cr->mpi_comm_mygroup);
#endif
    return potential_;
};

real ExternalPotential::potential()
{
    impl_->reduce_potential_virial(cr);

}

void ExternalPotential::add_forces( rvec f[], tensor vir, t_commrec *cr,
    gmx_int64_t step, real t, real * V_total, real weight)
{
    if (impl_->do_this_step(step))
    {
        impl_->add_forces(f, weight);
    }
};

bool ExternalPotential::Impl::do_this_step(int step)
{
    return (step % nst_apply_ == 0 );
};
