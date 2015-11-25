#include "externalpotential.h"
#include "gromacs/topology/topology.h"
#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/sim_util.h"
#include "gromacs/mdlib/groupcoord.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/utility/exceptions.h"

#include <memory>

/*! \brief
 *  Read all generic data for experimental input modules on the master node
 */
class ExternalPotential::Impl
{
    public:
        Impl(struct ext_pot_ir * ep_ir, t_commrec * cr, t_inputrec * ir,
            const gmx_mtop_t* mtop, rvec x[], matrix box, FILE *input_file,
            FILE *output_file, FILE *fplog, gmx_bool bVerbose,
            const gmx_output_env_t *oenv, unsigned long Flags);

        ~Impl();

        real reduce_potential_virial(t_commrec * cr);
        void dd_make_local_groups( gmx_domdec_t *dd);
        rvec * x_assembled(gmx_int64_t step, t_commrec *cr, rvec *x, matrix box, int group_index);
        bool do_this_step(int step);
        void add_forces(rvec f[], real w);
        void add_virial(tensor vir, real total_weight);
        real potential();

    private:

        gmx_int64_t             nst_apply_; /**< every how many steps to apply */
        real                    potential_; /**< the (local) contribution to the potential */
        tensor                  total_virial_;

        FILE                   *input_file_;
        FILE                   *output_file_;
        FILE                   *logfile_;

        gmx_bool                bVerbose_; /**< -v flag from command line                      */
        unsigned long           Flags_;
        const gmx_output_env_t *oenv_;

        class Group;
        std::vector<std::unique_ptr<ExternalPotential::Impl::Group>> atom_groups_;

        void mpi_start_(t_commrec * cr);
        void mpi_buffer_resize_(int size);
        void mpi_sum_reduce_();
        void mpi_buffer_(real value);
        void mpi_buffer_(real * vector,  int size);
        void mpi_buffer_(real matrix[3][3]);
        void mpi_finish_();

        int                   *mpi_inbuf_;
        int                   *mpi_outbuf_;
        int                    mpi_buf_size_;
        int                    mpi_current_;
        bool                   mpi_buf_write_;
        t_commrec             *cr_;

};

void ExternalPotential::add_virial(tensor vir, gmx_int64_t step,real weight){
    if (impl_->do_this_step(step))
    {
        impl_->add_virial(vir, weight);
    };
};

void ExternalPotential::Impl::add_virial(real (*vir)[3], real total_weight)
{
    for (int i = XX; i < DIM; i++) {
        for (int j = 0; j < DIM; j++) {
            vir[i][j]+=total_weight*total_virial_[i][j];
        }
    }
}

class ExternalPotential::Impl::Group
{
    public:

        Group(t_inputrec * ir,t_commrec * cr,const gmx_mtop_t * mtop,size_t nat, int *ind, rvec *x, matrix box);
        void make_local_group_indices(gmx_domdec_t *dd);
        void communicate_positions_all_to_all(t_commrec *cr, rvec *x, matrix box);
        rvec * x_assembled(gmx_int64_t step, t_commrec *cr, rvec *x, matrix box);
        void add_forces(rvec f[],real w);
        void add_virial(tensor vir);

    private:

        rvec         *x_assembled_;     /**< the atoms of the ind_ group, made whole, using x_assembled_old_ as reference*/
        rvec         *x_assembled_old_; /**< the whole atoms from the ind group from the previous step as reference */
        ivec         *x_shifts_;        /**< helper for making the molecule whole */
        ivec         *extra_shifts_;    /**< helper variables to assemble a whole molecule,*/

        int           num_atoms_;       /**< Number of atoms that are influenced by the external potential. */
        int          *ind_;             /**< Global indices of the atoms.                */
        int          *coll_ind_;        /**< the whole atom number array */

        int           num_atoms_loc_;   /**< Part of the atoms that are local; set by make_local_group_indices. */
        int          *ind_loc_;         /**< Local indices of the external potential atoms; set by make_local_group_indices.*/

        int           nalloc_loc_;      /**< keep memory allocated; set by make_local_group_indices.*/
        rvec         *f_loc_;           /**< the forces from external potential on the local node */

        tensor        virial_loc_;      /**< the local contribution to the virial */
        real         *weight_loc_;      /**< Weights for the local indices */

        gmx_bool      bUpdateShifts_;   /**< perfom the coordinate shifts in neighboursearching steps */
        gmx_int64_t   last_comm_step_;  /**< the last time, coordinates were communicated all to all */

};

ExternalPotential::Impl::~Impl()
{
    free(mpi_inbuf_);
    free(mpi_outbuf_);
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

ExternalPotential::~ExternalPotential(){};

ExternalPotential::Impl::Impl(
    struct ext_pot_ir *ep_ir, t_commrec * cr, t_inputrec * ir, const gmx_mtop_t* mtop, rvec x[], matrix box,
        FILE *input_file, FILE *output_file, FILE *fplog, gmx_bool bVerbose,
        const gmx_output_env_t *oenv, unsigned long Flags ) :
      input_file_(input_file), output_file_(output_file), logfile_(fplog),
      bVerbose_(bVerbose), Flags_(Flags), oenv_(oenv),
      mpi_inbuf_(nullptr), mpi_outbuf_(nullptr),
      mpi_buf_size_(0), mpi_current_(0), mpi_buf_write_(true)
{

    for (int i = 0; i < ep_ir->number_index_groups ; i++) {
        atom_groups_.push_back(
            std::unique_ptr<ExternalPotential::Impl::Group>(
                new ExternalPotential::Impl::Group(ir, cr, mtop,ep_ir->nat[i],ep_ir->ind[i],x,box)
            )
        );
    };

};

void ExternalPotential::Impl::mpi_buffer_resize_(int size)
{
    mpi_buf_size_=size;
    srenew(mpi_inbuf_, mpi_buf_size_);
    srenew(mpi_outbuf_, mpi_buf_size_);
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

void ExternalPotential::Impl::Group::add_virial(tensor vir)
{
    for (size_t i = XX; i < DIM; i++) {
        for (size_t j = XX; j < DIM; j++) {
            vir[i][j]+=virial_loc_[i][j];
        }
    }
};


void ExternalPotential::Impl::Group::add_forces(rvec f[],real w)
{
        for (int l = 0; l < num_atoms_loc_; l++)
        {
            /* Get the right index of the local force, since typically not all local
            * atoms are subject to density fitting forces */
            /* Add to local force */
            for (int i = XX; i < DIM; i++) {
                f_loc_[l][i]*=w;
            }
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

void ExternalPotential::Impl::mpi_buffer_(real value)
{
    if(mpi_current_ + 1< mpi_buf_size_)
    {
        if(mpi_buf_write_){
            mpi_buffer_resize_(mpi_current_ + 1);
        }else{
            GMX_THROW(gmx::InternalError("Trying to read more values from MPI output buffer than allocated in target buffer."));
        }
    }
    if(mpi_buf_write_)
    {
        mpi_inbuf_[mpi_current_++]=value;
    }else
    {
        if MASTER(cr_){
            value=mpi_outbuf_[mpi_current_++];
        }
    };
};



void ExternalPotential::Impl::mpi_buffer_( real * vector, int size)
{
    if(cr_!=nullptr)
    {
        if(mpi_current_ + size< mpi_buf_size_)
        {
            if(mpi_buf_write_){
                mpi_buffer_resize_(mpi_current_ +size  );
            }else{
                GMX_THROW(gmx::InternalError("Trying to read more values from MPI output buffer than allocated in target buffer."));
            }
        }
        if(mpi_buf_write_)
        {
            for (int i = 0; i < size; i++) {
                    mpi_inbuf_[mpi_current_++]=vector[i];
            };
        }else
        {
            if MASTER(cr_)
            {
                for (int i = 0; i < size; i++) {
                    vector[i]=mpi_outbuf_[mpi_current_++];
                };
            }
        };
};

};

void ExternalPotential::Impl::mpi_buffer_(real matrix[DIM][DIM])
{
    if (cr_!=nullptr){

        if(mpi_current_ + DIM*DIM < mpi_buf_size_)
        {
            if(mpi_buf_write_){
                mpi_buffer_resize_(mpi_current_ + DIM*DIM);
            }else{
                GMX_THROW(gmx::InternalError("Trying to read more values from MPI output buffer than allocated in target buffer."));
            }
        }
        if(mpi_buf_write_)
        {
            for (size_t i = XX; i < DIM; i++) {
                for (size_t j = XX ; j < DIM; j++) {
                    mpi_inbuf_[mpi_current_++]=matrix[i][j];
                };
            };
        }else
        {
            if MASTER(cr_){

            for (size_t i = XX; i < DIM; i++) {
                for (size_t j = XX ; j < DIM; j++) {
                    matrix[i][j]=mpi_outbuf_[mpi_current_++];
                };
            };
        };
        };
    };
};

void ExternalPotential::Impl::mpi_start_(t_commrec * cr)
{
    cr_=cr;
}

void ExternalPotential::Impl::mpi_sum_reduce_(){

#ifdef GMX_MPI
    MPI_Reduce(mpi_inbuf_, mpi_outbuf_, mpi_buf_size_, GMX_MPI_REAL, MPI_SUM, MASTERRANK(cr_), cr_->mpi_comm_mygroup);
#endif
    mpi_buf_write_=false;
    mpi_current_=0;
};

void ExternalPotential::Impl::mpi_finish_()
{
    mpi_buf_write_=true;
    cr_=nullptr;
}


real ExternalPotential::Impl::reduce_potential_virial(t_commrec * cr)
{
    clear_mat(total_virial_);

    // collect the contribution to the virial from all groups
    for ( auto && group : atom_groups_)
    {
        group->add_virial(total_virial_);
    }

    mpi_start_(cr);

    mpi_buffer_(total_virial_);
    mpi_buffer_(potential_);

    mpi_sum_reduce_();

    mpi_buffer_(total_virial_);
    mpi_buffer_(potential_);

    mpi_finish_();

    if MASTER(cr){
        return potential_;
    }

    return 0;

};

real ExternalPotential::sum_reduce_potential_virial(t_commrec * cr)
{
    return impl_->reduce_potential_virial(cr);
}

void ExternalPotential::add_forces( rvec f[], gmx_int64_t step, real weight)
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

ExternalPotential::ExternalPotential (struct ext_pot_ir *ep_ir, t_commrec * cr,
    t_inputrec * ir, const gmx_mtop_t* mtop, rvec x[], matrix box, FILE *input_file, FILE *output_file, FILE *fplog, gmx_bool bVerbose, const gmx_output_env_t *oenv, unsigned long Flags)
{
    impl_=std::unique_ptr<Impl>(new Impl(ep_ir, cr, ir, mtop, x, box, input_file, output_file, fplog, bVerbose, oenv, Flags));
};
