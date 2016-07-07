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
#ifndef GMX_ATOMSET_H
#define GMX_ATOMSET_H

#include <vector>
#include <memory>

struct gmx_ga2la_t;

namespace gmx
{

class AtomSet
{
    public:

        ~AtomSet();

        static std::unique_ptr<AtomSet> create();
        void init(const int nat, const int *ind, bool bParallel);

        /*! \brief Sets the local and collective indices from a lookup in ga2la.
         *
         * global_index_ has to be set,
         */

        void set_local_and_collective_indices(const gmx_ga2la_t  *ga2la);

        const std::vector<int> &collective_index() const;
        int num_atoms_local() const;
        int num_atoms_global() const;

        std::vector<int>::const_iterator begin_global_index_subset(int subset_number, int total_number) const;
        std::vector<int>::const_iterator end_global_index_subset(int subset_number, int total_number) const;

        std::vector<int>::const_iterator begin_collective_index_subset(int subset_number, int total_number) const;
        std::vector<int>::const_iterator end_collective_index_subset(int subset_number, int total_number) const;

        std::vector<int>::const_iterator begin_subset(int subset_number, int total_number) const;
        std::vector<int>::const_iterator end_subset(int subset_number, int total_number) const;

        std::vector<int>::const_iterator begin() const;
        std::vector<int>::const_iterator end() const;

        std::vector<int>::const_iterator begin_global_index() const;
        std::vector<int>::const_iterator end_global_index() const;

        std::vector<int>::const_iterator begin_collective_index() const;
        std::vector<int>::const_iterator end_collective_index() const;


    private:
        AtomSet();

        std::vector<int>             global_index_;      /**< Global indices of the atoms in this set, e.g., as in the index file.*/
        std::vector<int>             collective_index_;  /**< Maps indices on one node (0..num_atoms_local_) to global atom indicices. */
        std::vector<int>             local_index_;       /**< Local indices of the external potential atoms, used to add to local force;
                                                            set by set_indices. Access the coordinate of the i-th local atom of this set by x[local_index_[i]].
                                                            Constructed and updated every domain-decomposition step by set_indices. */
};



} // namespace gmx

#endif
