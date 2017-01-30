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
#include "gromacs/math/volumedata/hilbertsort.h"
// #include "gromacs/fileio/xtcio.h"

namespace gmx
{

WholeMoleculeGroup::WholeMoleculeGroup(const Group &base_group, MpiHelper * mpi_helper, const matrix box, const int npbcdim,  const gmx_mtop_t *top_global, const int ePBC) :
    Group {base_group}, mpi_helper_ {
    mpi_helper
}, npbcdim_ {
    npbcdim
}
{
    shifts_.resize(Group::num_atoms_, {0, 0, 0});
    copy_mat(box, box_);
    /* Prepare a whole reference group on the master node */
    x_reference_.resize(Group::num_atoms_);
    if (mpi_helper_ == nullptr || mpi_helper_->isMaster())
    {
        std::vector<RVec> x_pbc(top_global->natoms);
        copy(Group::x_, Group::x_ + top_global->natoms, x_pbc.begin());

        do_pbc_first_mtop(nullptr, ePBC, box_, top_global, as_rvec_array(x_pbc.data()));

        for (int i = 0; i < Group::num_atoms_; i++)
        {
            copy_rvec(x_pbc[Group::ind_[i]], x_reference_[i]);
        }
    }
    // out_  = open_xtc("debug.xtc", "w");
    // write_xtc(out_, Group::num_atoms_, 0, 0, box, as_rvec_array(x_reference_.data()), 1000);


    // each rank keeps a global number of atoms large array for atom coordinates for easy sum reduction
    x_collective_.resize(Group::num_atoms_);
    // only on the master we keep the last known whole molecule configuration as reference to calculate the box vector shifts
    // all nodes get the shifts for all atoms
    if (mpi_helper_ != nullptr)
    {
        all_group_coordinates_to_master_(x_collective_);
    }
    else
    {
        for (const auto &atom : *this)
        {
            copy_rvec(atom.x, x_collective_[*(atom.i_global)]);
        }
    }
    // write_xtc(out_, Group::num_atoms_, 0, 0, box, as_rvec_array(x_collective_.data()), 1000);
    update_shifts_();
    if (mpi_helper_ != nullptr)
    {
        mpi_helper_->broadcast(*as_vec_array(shifts_.data()), DIM*Group::num_atoms_);
    }
    Group::xTransformed_.resize(Group::num_atoms_loc_);
    set_x(Group::x_, box);
}

WholeMoleculeGroup::~WholeMoleculeGroup()
{
    // close_xtc(out_);
};

void
WholeMoleculeGroup::reSort_(std::vector<int> &to_sort, const std::vector<int> &sortIndex)
{
    std::vector<int> sorted(sortIndex.size());
    for (std::size_t i  = 0; i < to_sort.size(); ++i)
    {
        sorted[i] = to_sort[sortIndex[i]];
    }
    to_sort = sorted;
}

void
WholeMoleculeGroup::reSortrVec_(rvec * to_sort, const std::vector<int> &sortIndex)
{
    std::vector<RVec> new_to_sort(sortIndex.size());
    for (std::size_t i  = 0; i < sortIndex.size(); ++i)
    {
        new_to_sort[i] = to_sort[sortIndex[i]];
    }
    for (std::size_t i  = 0; i < sortIndex.size(); ++i)
    {
        copy_rvec(new_to_sort[i], to_sort[i]);
    }
}

void
WholeMoleculeGroup::reSortReal_(std::vector<real> &to_sort, const std::vector<int> &sortIndex)
{
    std::vector<real> sorted(sortIndex.size());
    for (std::size_t i  = 0; i < sortIndex.size(); ++i)
    {
        sorted[i] = to_sort[sortIndex[i]];
    }
    for (std::size_t i  = 0; i < sortIndex.size(); ++i)
    {
        to_sort[i] = sorted[i];
    }
}


void
WholeMoleculeGroup::medianSort()
{
    hilbertMedianSort(Group::xTransformed_, sortIndex_);
    reSort_(ind_loc_, sortIndex_);
    reSort_(coll_ind_, sortIndex_);
    reSortrVec_(as_rvec_array(Group::xTransformed_.data()), sortIndex_);
    reSortrVec_(as_rvec_array(x_collective_.data()), sortIndex_);
    reSortrVec_(as_rvec_array(x_reference_.data()), sortIndex_);
    reSortReal_(Group::weights_, sortIndex_);
}

void
WholeMoleculeGroup::set_x(const rvec x[], const matrix box)
{
    copy_mat(box, box_);
    for (int i = 0; i < Group::num_atoms_loc_; i++)
    {
        int i_local          = Group::ind_loc_[i];
        int i_global         = Group::coll_ind_[i];
        Group::xTransformed_[i][XX] = x[i_local][XX]+shifts_[i_global][XX]*box_[XX][XX]+shifts_[i_global][XX]*box_[YY][XX]+shifts_[i_global][XX]*box_[ZZ][XX];
        Group::xTransformed_[i][YY] = x[i_local][YY]+shifts_[i_global][YY]*box_[XX][YY]+shifts_[i_global][YY]*box_[YY][YY]+shifts_[i_global][YY]*box_[ZZ][YY];
        Group::xTransformed_[i][ZZ] = x[i_local][ZZ]+shifts_[i_global][ZZ]*box_[XX][ZZ]+shifts_[i_global][ZZ]*box_[YY][ZZ]+shifts_[i_global][ZZ]*box_[ZZ][ZZ];
    }
}

const std::vector<RVec> &
WholeMoleculeGroup::xTransformed()
{
    return Group::xTransformed_;
}

const matrix *
WholeMoleculeGroup::box()
{
    return &box_;
}

RVec
WholeMoleculeGroup::local_torque_sum(RVec center)
{
    return std::accumulate(
            this->begin(), this->end(), RVec {0, 0, 0},
            [center] (RVec &torque_sum, const GroupAtom &local_atom) {
                RVec xCentered;
                RVec torque;
                rvec_sub(*local_atom.xTransformed, center, xCentered);
                cprod(*local_atom.force, xCentered, torque);
                rvec_inc(torque_sum, torque);
                return torque_sum;
            });
}

RVec
WholeMoleculeGroup::weighted_local_coordinate_sum() const
{
    return std::accumulate(
            this->begin(), this->end(), RVec {0, 0, 0},
            [](RVec &coordinate_sum, const GroupAtom &local_atom) {
                RVec weighted_coordinate;
                svmul(*local_atom.properties, *local_atom.xTransformed, weighted_coordinate);
                rvec_inc(coordinate_sum, weighted_coordinate);
                return coordinate_sum;
            });
}

RVec
WholeMoleculeGroup::weightedCenterOfMass() const
{

    RVec weightedCenterOfMass_ = weighted_local_coordinate_sum();
    real weightssum            = local_weights_sum();
    if (mpi_helper_ != nullptr)
    {
        mpi_helper_->to_reals_buffer(weightedCenterOfMass_, 3);
        mpi_helper_->to_reals_buffer(&weightssum, 1);
        mpi_helper_->sum_reduce();
        if (mpi_helper_->isMaster())
        {
            mpi_helper_->from_reals_buffer(weightedCenterOfMass_, 3);
            mpi_helper_->from_reals_buffer(&weightssum, 1);
        }
    }
    if (mpi_helper_ == nullptr || mpi_helper_->isMaster())
    {
        svmul(1.0/weightssum, weightedCenterOfMass_, weightedCenterOfMass_);
    }
    return weightedCenterOfMass_;
}
/* Ensure this is called after Group::set_indices */
void
WholeMoleculeGroup::update_shifts_and_reference(const rvec x[], const matrix box)
{
    copy_mat(box, box_);
    set_x(x, box);
    // write_xtc(out_, Group::num_atoms_, 0, 0, box, as_rvec_array(x_reference_.data()), 1000);

    if (mpi_helper_ != nullptr)
    {
        all_group_coordinates_to_master_(x_collective_);
    }
    else
    {
        for (const auto &atom : *this)
        {
            x_collective_[*(atom.i_global)] = *atom.xTransformed;
        }
    }
    // write_xtc(out_, Group::num_atoms_, 0, 0, box, as_rvec_array(x_collective_.data()), 1000);
    if ((mpi_helper_ == nullptr) || (mpi_helper_->isMaster()) )
    {
        update_shifts_();
    }
    if (mpi_helper_ != nullptr)
    {
        mpi_helper_->broadcast(shifts_.data(), Group::num_atoms_);
    }

    // re-calculated atom coordinates with the new shifts and take this as new reference
    set_x(x, box);

    if (mpi_helper_ != nullptr)
    {
        all_group_coordinates_to_master_(x_reference_);
    }
    else
    {
        for (const auto &atom : *this)
        {
            x_reference_[*(atom.i_global)] = *atom.xTransformed;
        }
    }
    // write_xtc(out_, Group::num_atoms_, 0, 0, box, as_rvec_array(x_reference_.data()), 1000);
    // write_xtc(out_, Group::num_atoms_, 0, 0, box, as_rvec_array(x_collective_.data()), 1000);
}

void
WholeMoleculeGroup::all_group_coordinates_to_master_(std::vector<RVec> &vec)
{
    // put all node local atoms in the right position in a global array for later sum_reduction
    for (const auto &atom : *this)
    {
        vec[*(atom.i_global)] = *atom.xTransformed;
    }
    mpi_helper_->to_reals_buffer(*as_rvec_array(vec.data()), 3*Group::num_atoms_);
    mpi_helper_->sum_reduce();
    mpi_helper_->from_reals_buffer(*as_rvec_array(vec.data()), 3*Group::num_atoms_);

}

/** copied from groupcoords.cpp*/
void
WholeMoleculeGroup::update_shifts_()
{
    rvec dx;

    /* Get the shifts such that each atom is within closest
     * distance to its position at the last NS time step after shifting.
     * If we start with a whole group, and always keep track of
     * shift changes, the group will stay whole this way */

    for (int i = 0; i < Group::num_atoms_; i++)
    {
        /* The distance this atom moved since the last time step */
        /* If this is more than just a bit, it has changed its home pbc box */
        rvec_sub(x_collective_[i], x_reference_[i], dx);

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
