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

#include "group.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/topology/topology.h"
#include "gromacs/mdlib/groupcoord.h"
#include "gromacs/mdlib/sim_util.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/math/vec.h"
#include "gromacs/domdec/domdec_struct.h"

namespace gmx
{
/*******************************************************************************
 * GroupCoordinates::Iterator
 */
GroupCoordinates::GroupIterator::GroupIterator(const GroupCoordinates * coordinates, int i)
    : coordinates_(coordinates), i_(i) {}

GroupCoordinates::GroupIterator::GroupIterator(const GroupCoordinates * coordinates)
    : coordinates_(coordinates), i_(0) {}

GroupCoordinates::GroupIterator::GroupIterator(const GroupIterator * iterator)
    : coordinates_(iterator->coordinates_), i_(iterator->i_) {}

GroupCoordinates::GroupIterator &GroupCoordinates::GroupIterator::operator++()
{
    ++i_;
    return *this;
}

GroupCoordinates::GroupIterator GroupCoordinates::GroupIterator::operator++(int)
{
    GroupIterator tmp(this);
    operator++();
    return tmp;
}

bool GroupCoordinates::GroupIterator::operator==(const GroupIterator &rhs)
{
    return (i_ == rhs.i_) && (coordinates_ == rhs.coordinates_);
}

bool GroupCoordinates::GroupIterator::operator!=(const GroupIterator &rhs)
{
    return (i_ != rhs.i_) || (coordinates_ != rhs.coordinates_);
}

const RVec GroupCoordinates::GroupIterator::operator*()
{
    return RVec(coordinates_->x_[coordinates_->index_[i_]]);
}

/*******************************************************************************
 * GroupCoordinates
 */

GroupCoordinates::GroupCoordinates(const rvec * x, int * index, int size)
    : x_(x), index_(index), size_(size) {}

GroupCoordinates::GroupIterator GroupCoordinates::begin()
{
    return GroupIterator(this);
};

GroupCoordinates::GroupIterator GroupCoordinates::end()
{
    return GroupIterator(this, size_);
};

const RVec GroupCoordinates::operator[](int i)
{
    return RVec(x_[index_[i]]);
};

/******************************************************************************
 * Group
 */

std::vector<RVec> &Group::f_local()
{
    f_loc_.resize(num_atoms_loc_, {0, 0, 0});
    return f_loc_;
}

Group::Group( int ePBC, t_commrec *cr, const gmx_mtop_t *mtop, int nat, int *ind, const rvec x[], matrix box) :
    num_atoms_(nat), ind_(ind), nalloc_loc_(1)
{
    if (!PAR(cr))
    {
        num_atoms_loc_ = num_atoms_;
        ind_loc_.assign(ind, ind+num_atoms_);
        nalloc_loc_ = num_atoms_loc_;
        coll_ind_.resize(num_atoms_loc_);
        for (int i = 0; i < num_atoms_loc_; i++)
        {
            coll_ind_[i] = i;
        }
    }
    make_whole_reference(ePBC, cr, mtop, x, box);
    x_assembled_.resize(x_assembled_old_.size());
    snew(x_shifts_, num_atoms_loc_);
    snew(extra_shifts_, num_atoms_loc_);
    return;
};

Group::~Group()
{
    sfree(x_shifts_);
    sfree(extra_shifts_);
}

void Group::make_whole_reference(int ePBC, t_commrec *cr, const gmx_mtop_t *mtop, const rvec x[], matrix box)
{
    /* To make the molecule whole, use the reference in x_assembled_old
     * create this reference in the first step */
    if (MASTER(cr))
    {
        /* Remove pbc and prepare a whole molecule for new_t_gmx_densfit */
        rvec *x_pbc            = NULL;
        snew(x_pbc, mtop->natoms);      /* There ... */
        copy_rvecn(x, x_pbc, 0, mtop->natoms);

        do_pbc_first_mtop(NULL, ePBC, box, mtop, x_pbc);
        for (int i = 0; i < num_atoms_; i++)
        {
            x_assembled_old_.push_back(x_pbc[ind_[i]]);
        }
        sfree(x_pbc);      /* ... and back again */
    }

    if (PAR(cr))
    {
        rvec *to_broadcast = as_rvec_array(x_assembled_old_.data());
        MPI_Bcast(to_broadcast, num_atoms_ * sizeof(rvec), MPI_BYTE, MASTERRANK(cr), cr->mpi_comm_mygroup);
        x_assembled_old_.assign(to_broadcast, to_broadcast+num_atoms_);
    }

}

void Group::set_indices(gmx_domdec_t *dd)
{

    int * local_ind = nullptr;
    int * coll_ind  = coll_ind_.data();
    dd_make_local_group_indices(dd->ga2la, num_atoms_, ind_, &num_atoms_loc_, &local_ind, &nalloc_loc_, coll_ind);
    ind_loc_.assign(local_ind, local_ind+num_atoms_loc_);
    sfree(local_ind);
    coll_ind_.assign(coll_ind, coll_ind+num_atoms_loc_);

    /* Indicate that the group's shift vectors for this structure need to be updated
     * at the next call to communicate_group_positions, since obviously we are in a NS step */
    bUpdateShifts_ = true;
};

const std::vector<RVec> &Group::x_assembled(const gmx_int64_t step, t_commrec * cr, const rvec x[], const matrix box)
{
    communicate_positions_all_to_all(step, cr, x, box);
    return x_assembled_;
};

void Group::add_forces(rvec f[], real w)
{
    for (int l = 0; l < num_atoms_loc_; l++)
    {
        for (int i = XX; i <= ZZ; i++)
        {
            f_loc_[l][i] *= w;
        }
        /* Get the right index of the local force, since typically not all local
         * atoms are subject to density fitting forces and add to local force */
        rvec_inc(f[ind_loc_[l]], f_loc_[l]);
    }
};

GroupCoordinates Group::x_local(const rvec x[])
{
    return GroupCoordinates(x, ind_loc_.data(), num_atoms_loc_);
}

const std::vector<int> &Group::collective_index()
{
    return coll_ind_;
}

void Group::communicate_positions_all_to_all(const gmx_int64_t step, t_commrec *cr, const rvec x[], const matrix box)
{
    if (step != last_comm_step_)
    {
        rvec * x_to_communicate     = as_rvec_array(x_assembled_.data());
        rvec * x_old_to_communicate = as_rvec_array(x_assembled_old_.data());
        communicate_group_positions( cr, x_to_communicate, x_shifts_, extra_shifts_, bUpdateShifts_, const_cast<rvec*>(x), num_atoms_, num_atoms_loc_, ind_loc_.data(), coll_ind_.data(), x_old_to_communicate, box);
        x_assembled_.assign(x_to_communicate, x_to_communicate + num_atoms_);
        x_assembled_old_.assign(x_old_to_communicate, x_old_to_communicate + num_atoms_);

        last_comm_step_ = step;
        bUpdateShifts_  = false;
    }
}
} // namespace gmx
