#ifndef _externalpotentialdata_h_
#define _externalpotentialdata_h_


#include <cstdio>

#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/topology/atom_id.h"

struct gmx_output_env_t;

/** !brief
    Read all generic data for experimental input modules on the master node
 */

class ExternalPotentialData
{
    public:
        ExternalPotentialData(
            FILE               *input_file,
            FILE               *output_file,
            size_t              nat,
            atom_id            *ind,
            rvec               *x_assembled_old,
            FILE               *fplog,
            gmx_bool            bVerbose,
            const gmx_output_env_t *oenv,
            unsigned long Flags
            );

        void make_local_group_indices(gmx_domdec_t *dd);

        rvec * x_assembled(gmx_int64_t step, t_commrec * cr, rvec * x, matrix box);

        void communicate_positions_all_to_all(
            t_commrec  *cr,
            rvec       *x,
            matrix      box    );

        double summed_potential_on_ranks(t_commrec * cr);

    private:

        gmx_int64_t              nst_apply_;       /**< every how many steps to apply */
        double                   potential_;       /**< the (local) contribution to the potential */

        FILE                    *input_file_;
        FILE                    *output_file_;
        FILE                    *logfile_;

        gmx_int64_t              last_comm_step_;  /**< the last time, coordinates were communicated all to all */
        rvec                    *x_assembled_;     /**< the atoms of the ind_ group, made whole, using x_assembled_old_ as reference*/
        rvec                    *x_assembled_old_; /**< the whole molecule from the previous step as reference */
        ivec                    *x_shifts_;        /**< helper for making the molecule whole */
        ivec                    *extra_shifts_;    /**< helper variables to assemble a whole molecule,*/

        gmx_bool                 bVerbose_;        /**< -v flag from command line                      */
        const gmx_output_env_t  *oenv_;
        unsigned long            Flags_;

        int                      num_atoms_;       /**< Number of atoms that are influenced by the external potential. */
        int                      num_atoms_loc_;   /**< Part of the atoms that are local; set by make_local_group_indices. */
        atom_id                 *ind_;             /**< Global indices of the atoms.                */
        atom_id                 *ind_loc_;         /**< Local indices of the external potential atoms; set by make_local_group_indices.*/
        int                      nalloc_loc_;      /**< keep memory allocated; set by make_local_group_indices.*/
        rvec                    *f_loc_;           /**< the forces from external potential on the local node */
        real                    *weight_loc_;      /**< Weights for the local indices */
        int                     *coll_ind_;        /**< the whole atom number array */
        gmx_bool                 bUpdateShifts_;   /**< perfom the coordinate shifts in neighboursearching steps */

};

#endif /* end of include guard: externalpotentialdata_h_
        */
