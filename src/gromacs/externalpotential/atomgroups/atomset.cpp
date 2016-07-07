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

#include "atomset.h"

#include <algorithm>
#include <memory>

#include "gromacs/domdec/ga2la.h"


namespace gmx
{





std::unique_ptr<AtomSet> AtomSet::create()
{
    return std::unique_ptr<AtomSet>(new AtomSet());
};


void AtomSet::init(const int nat, const int *ind, bool bParallel)
{
    num_atoms_global_ = nat;
    global_index_.assign(ind, ind+num_atoms_global_);
    if (!bParallel)
    {
        num_atoms_local_ = num_atoms_global_;
        local_index_.assign(ind, ind+num_atoms_global_);
        collective_index_.resize(num_atoms_local_);
        std::iota(collective_index_.begin(), collective_index_.end(), 0);
    }
};

AtomSet::~AtomSet()
{
}


void AtomSet::set_local_and_collective_indices(const gmx_ga2la_t *ga2la)
{
    /* Loop over all the atom indices of the group to check
     * which ones are on the local node */
    int  i_local;
    int  nalloc_loc = 0;
    int  n_local    = 0;
    local_index_.clear();
    collective_index_.clear();

    for (const auto &i_global : global_index_)
    {
        if (ga2la_get_home(ga2la, i_global, &i_local))
        {
            /* The atom with this index is a home atom ? */
            if (n_local >= nalloc_loc)  /* Check whether memory suffices */
            {
                nalloc_loc = over_alloc_dd(n_local+1);
                /* We never need more memory than the number of atoms in the group */
                nalloc_loc = std::min(nalloc_loc, n_local);

                local_index_.reserve(nalloc_loc);
                collective_index_.reserve(nalloc_loc);
            }
            /* Save the atoms index in the local atom numbers array */
            local_index_.push_back(i_local);

            /* Keep track of where this local atom belongs in the collective index array.
             * This is needed when reducing the local arrays to a collective/global array
             * in communicate_group_positions */
            collective_index_.push_back(i_global);
            n_local++;
        }
    }


};

const std::vector<int> &AtomSet::collective_index() const
{
    return collective_index_;
}

int AtomSet::num_atoms_global() const
{
    return global_index_.size();
}

int AtomSet::num_atoms_local() const
{
    return local_index_.size();
};

} // namespace gmx
