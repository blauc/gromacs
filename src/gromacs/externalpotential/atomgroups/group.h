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

#include <vector>
#include <memory>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"


struct t_commrec;
struct t_inputrec;
struct gmx_ga2la_t;
struct gmx_output_env_t;
struct gmx_mtop_t;


namespace gmx
{

class AtomProperties;

struct GroupAtom{
    const real     *x;
    real           *force;
    AtomProperties *properties;
    int             i_local;
    int             i_global;
};

class Group
{
    public:

        class GroupIterator : std::iterator<std::forward_iterator_tag, GroupAtom>
        {
            public:
                friend         Group;
                GroupIterator &operator++ ();
                GroupIterator  operator++ (int);
                bool           operator== (const GroupIterator &rhs);
                bool           operator!= (const GroupIterator &rhs);
                GroupAtom &operator*();
            private:
                GroupIterator(Group * group);
                GroupIterator(Group * group, int i);
                GroupIterator(const GroupIterator * iterator);
                Group * group_;
                int     i_;
        };

        friend class GroupIterator;

        GroupIterator begin();
        GroupIterator end();
        GroupAtom &operator[] (int i);

        /*\brief
         * copy constructor used for initialisation in derived classes
         */
        Group(const Group &group);
        Group(int nat, int *ind, bool bParallel);
        ~Group();
        void set_indices(gmx_ga2la_t  *ga2la);

        void add_forces(rvec f[], real w);
        void add_virial(tensor vir);
        const std::vector<int> &collective_index();

        void parallel_loop(std::function<void(GroupAtom &, const int &)> loop_kernel_function);

        int num_atoms_loc();
        int num_atoms_global() const;

        void set_x(const rvec x[]);
        RVec local_coordinate_sum();

    protected:

        rvec                        *x_;             /**< The coordinates TODO: re-attach const attribute..*/
        GroupAtom                    atom_;
        int                          num_atoms_;     /**< Number of (global) atoms in this external potential group. */
        int                         *ind_;           /**< Global indices of the atoms in this group.*/
        std::vector<int>             coll_ind_;      /**< map local atom indices to global atom indicices of atoms of this group (i.e. global_index=ind_[coll_ind_[local_index])]*/

        int                          num_atoms_loc_; /**< number of local atoms from index group; set by set_indices. */
        std::vector<int>             ind_loc_;       /**< Local indices of the external potential atoms, used to add to local force; set by set_indices.*/

        std::vector<RVec>            f_loc_;         /**< the forces from external potential on the local node */

        std::vector<AtomProperties*> properties_;    /**< Weights for the local indices */

        bool                         bParallel_;

};



} // namespace gmx
