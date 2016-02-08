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

#include "assembledmolecule.h"

#include "group.h"

namespace gmx
{
AssembledGroup::AssembledGroup(Group * group, t_commrec * cr) :
    cr_(cr), group_(group)
{

}

AssembledGroup::~AssembledGroup()
{
    sfree(x_shifts_);
    sfree(extra_shifts_);
};


void AssembledGroup::communicate_positions_all_to_all(const gmx_int64_t step, const rvec x[], const matrix box)
{
    if (step != last_comm_step_)
    {
        rvec * x_to_communicate     = as_rvec_array(x_assembled_.data());
        rvec * x_old_to_communicate = as_rvec_array(x_assembled_old_.data());
        communicate_group_positions(cr_, x_to_communicate, x_shifts_, extra_shifts_, bUpdateShifts_, const_cast<rvec*>(x), group_->num_atoms_, group_->num_atoms_loc_, group_->ind_loc_.data(), group_->coll_ind_.data(), x_old_to_communicate, box);
        x_assembled_.assign(x_to_communicate, x_to_communicate + group_->num_atoms_);
        x_assembled_old_.assign(x_old_to_communicate, x_old_to_communicate + group_->num_atoms_);

        last_comm_step_ = step;
        bUpdateShifts_  = false;
    }
}

const std::vector<RVec> &AssembledGroup::x_assembled(const gmx_int64_t step, const rvec x[], const matrix box)
{
    communicate_positions_all_to_all(step, x, box);
    return x_assembled_;
};

void AssembledGroup::make_whole_reference(int ePBC, const gmx_mtop_t *mtop, const rvec x[], matrix box)
{
    /* To make the molecule whole, use the reference in x_assembled_old
     * create this reference in the first step */
    if (MASTER(cr_))
    {
        /* Remove pbc and prepare a whole molecule for new_t_gmx_densfit */
        rvec *x_pbc            = NULL;
        snew(x_pbc, mtop->natoms);       /* There ... */
        copy_rvecn(x, x_pbc, 0, mtop->natoms);

        do_pbc_first_mtop(NULL, ePBC, box, mtop, x_pbc);
        for (int i = 0; i < group_->num_atoms_; i++)
        {
            x_assembled_old_.push_back(x_pbc[group_->ind_[i]]);
        }
        sfree(x_pbc);       /* ... and back again */
    }

    if (PAR(cr_))
    {
        rvec *to_broadcast;
        if (MASTER(cr_))
        {
            to_broadcast = as_rvec_array(x_assembled_old_.data());
        }
        else
        {
            snew(to_broadcast, group_->num_atoms_*sizeof(rvec));
        }

        MPI_Bcast(to_broadcast, group_->num_atoms_ * sizeof(rvec), GMX_MPI_REAL, MASTERRANK(cr_), cr_->mpi_comm_mygroup);
        x_assembled_old_.assign(to_broadcast, to_broadcast+group_->num_atoms_);
    }

    x_assembled_.resize(x_assembled_old_.size());
    snew(x_shifts_, group_->num_atoms_loc_);
    snew(extra_shifts_, group_->num_atoms_loc_);

}


const std::vector<RVec> &ExternalPotential::Impl::x_assembled(const gmx_int64_t step, t_commrec *cr, const rvec x[], const matrix box, int group_index)
{
    return atom_groups_[group_index]->x_assembled(step, cr, x, box);
}


const std::vector<RVec> &ExternalPotential::x_assembled(const gmx_int64_t step, t_commrec *cr, const rvec x[], const matrix box, int group_index)
{
    return impl_->x_assembled(step, cr, x, box, group_index);
};

const std::vector<int> &ExternalPotential::collective_index(int group_index)
{
    return impl_->collective_index(group_index);
}



const std::vector<int> &ExternalPotential::Impl::collective_index(int group_index)
{
    return atom_groups_[group_index]->collective_index();
}

}
