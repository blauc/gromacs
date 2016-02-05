/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */


namespace gmx
{
class AssembleMolecule
{
    protected:
        const std::vector<RVec> &x_assembled(const gmx_int64_t step, t_commrec *cr, const rvec x[], const matrix box, int group_index);

        /*! \brief
         * Identifies local group atoms in the assembled coordinates, such that x_assembled[collective_index(group_number)[i]] returns the ith atom in group "group_number".
         */
        const std::vector<int> &collective_index(int group_index);

}
class AssembledGroup
{
    public:
        AssembledGroup(Group* group, t_commrec * cr);
        ~AssembledGroup();
        const std::vector<RVec> &x_assembled(gmx_int64_t step, const rvec x[], const matrix box);

    private:

        void make_whole_reference(int ePBC, const gmx_mtop_t* mtop, const rvec x[], matrix  box);
        void communicate_positions_all_to_all(gmx_int64_t step, const rvec x[], const matrix box);

        std::vector<RVec> x_assembled_;     /**< all atoms in ind_ (not just the ones stored locally), made whole, using x_assembled_old_ as reference*/
        std::vector<RVec> x_assembled_old_; /**< reference positions of all atoms in ind_ from precious step, used to define "whole" */

        ivec             *x_shifts_;        /**< number of periodic boundary condition shifts, helper variable to assemble a whole molecule */
        ivec             *extra_shifts_;    /**< extra number of periodic boundary condition shifts, helper variable to assemble a whole molecule */

        gmx_int64_t       last_comm_step_;  /**< the last step when coordinates were communicated all to all, ensures all to all communication is not performed more often than necessary */
        bool              bUpdateShifts_;   /**< perfom the coordinate shifts in neighboursearching steps */
        t_commrec *       cr_;
        Group*            group_;



};
}
