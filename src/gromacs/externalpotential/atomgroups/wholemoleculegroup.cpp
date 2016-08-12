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

#include <algorithm>
#include <memory>

#include "gromacs/externalpotential/mpi-helper.h"
#include "wholemoleculegroup.h"
#include "gromacs/domdec/ga2la.h"
#include "gromacs/externalpotential/forceplotter.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/groupcoord.h"
#include "gromacs/mdlib/sim_util.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/mdlib/groupcoord.h"

namespace gmx
{

WholeMoleculeGroup::WholeMoleculeGroup(const Group &base_group, std::shared_ptr<MpiHelper> mpi_helper, const matrix box, const int npbcdim) :
    Group(base_group), mpi_helper_(mpi_helper), npbcdim_(npbcdim)
{
    // in the whole molecule group, copy atom coordinates, then shift, rather then reference only
    snew(Group::x_, Group::num_atoms_);
    snew(x_reference_, Group::num_atoms_);
    snew(shifts_, Group::num_atoms_);

    set_box(box);
    if ((mpi_helper_ == nullptr) || (mpi_helper_->isMaster()))
    {
        snew(x_coll_, Group::num_atoms_);
    }
    all_group_coordinates_to_master_();
}

WholeMoleculeGroup::~WholeMoleculeGroup()
{
    sfree(Group::x_);
    sfree(x_reference_);
    sfree(shifts_);
    sfree(x_coll_);
};

void
WholeMoleculeGroup::set_x(const rvec x[])
{
    int i_local;
    int i_global;
    for (int i = 0; i < Group::num_atoms_loc_; i++)
    {
        i_local          = Group::ind_loc_[i];
        i_global         = Group::coll_ind_[i];
        Group::x_[i][XX] = x[i_local][XX]+shifts_[i_global][XX]*box_[XX][XX]+shifts_[i_global][XX]*box_[YY][XX]+shifts_[i_global][XX]*box_[ZZ][XX];
        Group::x_[i][YY] = x[i_local][YY]+shifts_[i_global][YY]*box_[YY][XX]+shifts_[i_global][YY]*box_[YY][YY]+shifts_[i_global][YY]*box_[ZZ][YY];
        Group::x_[i][ZZ] = x[i_local][ZZ]+shifts_[i_global][ZZ]*box_[XX][ZZ]+shifts_[i_global][ZZ]*box_[YY][ZZ]+shifts_[i_global][ZZ]*box_[ZZ][ZZ];
    }
}

void
WholeMoleculeGroup::set_box(const matrix box)
{
    copy_mat(box, box_);
}

/* Ensure this is called after Group::set_indices */
void
WholeMoleculeGroup::update_shifts_and_reference(const rvec x[], const matrix box)
{

    Group::set_x(x);
    set_box(box);

    all_group_coordinates_to_master_();
    if ((mpi_helper_ == nullptr) || (mpi_helper_->isMaster()) )
    {
        calculate_shifts_();
        for (int i = 0; i < Group::num_atoms_; i++)
        {
            copy_rvec((*this)[i].x, x_reference_[i]);
        }
    }
    if (mpi_helper_ != nullptr)
    {
        mpi_helper_->broadcast(shifts_, Group::num_atoms_);
    }
}

void
WholeMoleculeGroup::all_group_coordinates_to_master_()
{
    if (mpi_helper_ != nullptr)
    {
        // put all node local atoms in the right position in a global array for later sum_reduction
        for (const auto &atom : *this)
        {
            copy_rvec(atom.x, x_coll_[*(atom.i_global)]);
        }
        mpi_helper_->to_reals_buffer(*x_coll_, 3*Group::num_atoms_);
        mpi_helper_->sum_reduce();
        mpi_helper_->from_reals_buffer(*x_coll_, 3*Group::num_atoms_);
    }
}

void
WholeMoleculeGroup::make_whole_molecule_reference(const rvec x[], const gmx_mtop_t *top_global, const int ePBC)
{
    rvec * x_pbc;
    snew(x_pbc, top_global->natoms);
    copy_rvecn(x, x_pbc, 0, top_global->natoms);

    do_pbc_first_mtop(nullptr, int(ePBC), box_, top_global, x_pbc);

    for (int i = 0; i < Group::num_atoms_; i++)
    {
        copy_rvec(x_pbc[Group::ind_[i]], x_reference_[i]);
    }
    sfree(x_pbc);
    calculate_shifts_();
    if (mpi_helper_ != nullptr)
    {
        mpi_helper_->broadcast(shifts_, DIM*Group::num_atoms_);
    }

}

/** copied from groupcoords.cpp*/
void
WholeMoleculeGroup::calculate_shifts_()
{
    rvec dx;

    /* Get the shifts such that each atom is within closest
     * distance to its position at the last NS time step after shifting.
     * If we start with a whole group, and always keep track of
     * shift changes, the group will stay whole this way */
    for (int i = 0; i < Group::num_atoms_; i++)
    {
        clear_ivec(shifts_[i]);
    }

    for (int i = 0; i < Group::num_atoms_; i++)
    {
        /* The distance this atom moved since the last time step */
        /* If this is more than just a bit, it has changed its home pbc box */
        rvec_sub((*this)[i].x, x_reference_[i], dx);

        for (int m = npbcdim_-1; m >= 0; m--)
        {
            while (dx[m] < -0.5*box_[m][m])
            {
                for (int d = 0; d < DIM; d++)
                {
                    dx[d] += box_[m][d];
                }
                shifts_[i][m]++;
            }
            while (dx[m] >= 0.5*box_[m][m])
            {
                for (int d = 0; d < DIM; d++)
                {
                    dx[d] -= box_[m][d];
                }
                shifts_[i][m]--;
            }
        }
    }
}



} // namespace gmx
