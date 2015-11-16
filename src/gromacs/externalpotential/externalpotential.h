#ifndef _externalpotential_h_
#define _externalpotential_h_

#include "gromacs/timing/wallcycle.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/topology/atom_id.h"
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

class ExternalPotentialGroup
{
    public:
        ExternalPotentialGroup(t_inputrec * ir,t_commrec * cr,const gmx_mtop_t * mtop,size_t              nat,
        atom_id            *ind,
        rvec               *x, matrix box);
        void make_local_group_indices(gmx_domdec_t *dd);
        void communicate_positions_all_to_all(t_commrec *cr, rvec *x, matrix box);
        rvec * x_assembled(gmx_int64_t step, t_commrec *cr, rvec *x, matrix box);

    private:
        rvec                    *x_assembled_;     /**< the atoms of the ind_ group, made whole, using x_assembled_old_ as reference*/
        rvec                    *x_assembled_old_; /**< the whole atoms from the ind group from the previous step as reference */
        ivec                    *x_shifts_;        /**< helper for making the molecule whole */
        ivec                    *extra_shifts_;    /**< helper variables to assemble a whole molecule,*/

        atom_id                  num_atoms_;       /**< Number of atoms that are influenced by the external potential. */
        atom_id                 *ind_;             /**< Global indices of the atoms.                */

        int                      num_atoms_loc_;   /**< Part of the atoms that are local; set by make_local_group_indices. */
        atom_id                 *ind_loc_;         /**< Local indices of the external potential atoms; set by make_local_group_indices.*/
        int                      nalloc_loc_;      /**< keep memory allocated; set by make_local_group_indices.*/
        rvec                    *f_loc_;           /**< the forces from external potential on the local node */
        real                    *weight_loc_;      /**< Weights for the local indices */
        int                     *coll_ind_;        /**< the whole atom number array */
        gmx_bool                 bUpdateShifts_;   /**< perfom the coordinate shifts in neighboursearching steps */
        gmx_int64_t              last_comm_step_;  /**< the last time, coordinates were communicated all to all */

};

typedef std::unique_ptr<ExternalPotentialGroup> ExternalPotentialGroupPointer;

/*! \brief
 *  Read all generic data for experimental input modules on the master node
 */
class ExternalPotentialData
{
    public:
        ExternalPotentialData(
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

        double summed_potential_on_ranks(t_commrec * cr);
        void dd_make_local_groups( gmx_domdec_t *dd);
        rvec * x_assembled(gmx_int64_t step, t_commrec *cr, rvec *x, matrix box, int group_index);

    private:

        gmx_int64_t                                nst_apply_; /**< every how many steps to apply */
        double                                     potential_; /**< the (local) contribution to the potential */

        FILE                                      *input_file_;
        FILE                                      *output_file_;
        FILE                                      *logfile_;

        gmx_bool                                   bVerbose_; /**< -v flag from command line                      */
        unsigned long                              Flags_;
        const gmx_output_env_t                    *oenv_;

        std::vector<ExternalPotentialGroupPointer> atom_groups_;
};

typedef std::shared_ptr<ExternalPotentialData> ExternalPotentialDataPointer;

class ExternalPotential
{
    public:
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
        ExternalPotentialDataPointer data_;
        rvec* x_assembled(gmx_int64_t step, t_commrec *cr, rvec x[], matrix box, int group_index);

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
