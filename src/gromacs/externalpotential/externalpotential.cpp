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
#include "externalpotential.h"

#include <memory>
#include <vector>

#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/externalpotential/externalpotentialIO.h"
#include "gromacs/externalpotential/atomgroups/group.h"
#include "gromacs/externalpotential/atomgroups/wholemoleculegroup.h"
#include "gromacs/externalpotential/mpi-helper.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/groupcoord.h"
#include "gromacs/mdlib/sim_util.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/smalloc.h"


namespace gmx
{

namespace externalpotential
{


/******************************************************************************
 * ExternalPotential::Impl
 */

class ExternalPotential::Impl
{
    public:
        Impl();

        ~Impl();

        gmx_int64_t                          nst_apply_; /**< every how many steps to apply; \TODO for every step, since not implemented  */
        real                                 potential_; /**< the (local) contribution to the potential */
        tensor                               virial_;

        std::shared_ptr<MpiHelper>           mpi_;
        std::shared_ptr<ExternalPotentialIO> input_output_;
        std::vector < std::shared_ptr < Group>> atom_groups_;
        std::vector < std::shared_ptr < WholeMoleculeGroup>> wholemoleculegroups_;
};

ExternalPotential::Impl::Impl() :
    nst_apply_(1), potential_(0), mpi_(nullptr)
{
};

ExternalPotential::Impl::~Impl()
{
};

/*******************************************************************************
 * ExternalPotential
 */

ExternalPotential::ExternalPotential ()
    : impl_(new Impl())
{
};

ExternalPotential::~ExternalPotential(){};

void ExternalPotential::dd_make_local_groups( gmx_ga2la_t  * ga2la)
{
    for (auto && grp : impl_->atom_groups_)
    {
        grp->set_indices(ga2la);
    }
};

void ExternalPotential::add_virial(tensor vir, gmx_int64_t step, real weight)
{
    if (do_this_step(step))
    {
        for (int i = XX; i <= ZZ; i++)
        {
            for (int j = XX; j <= ZZ; j++)
            {
                vir[i][j] += weight*impl_->virial_[i][j];
            }
        }

    }
};
bool ExternalPotential::do_this_step(gmx_int64_t step)
{
    return (step % impl_->nst_apply_ == 0 );
}

real ExternalPotential::sum_reduce_potential_virial()
{
    if (impl_->mpi_ != nullptr)
    {
        impl_->mpi_->to_reals_buffer(impl_->virial_);
        impl_->mpi_->to_reals_buffer(impl_->potential_);

        impl_->mpi_->sum_reduce();

        impl_->potential_ = impl_->mpi_->from_reals_buffer();
        impl_->mpi_->from_reals_buffer(impl_->virial_);

        impl_->mpi_->finish();

    }

    if (impl_->mpi_ == nullptr || impl_->mpi_->isMaster())
    {
        return impl_->potential_;
    }

    return 0.0;
}

void ExternalPotential::add_forces( rvec f[], gmx_int64_t step, real weight)
{
    if (do_this_step(step))
    {
        for (auto && group : impl_->atom_groups_)
        {
            group->add_forces(f, weight);
        }

    }
};

void ExternalPotential::add_forces_capped( rvec f[], gmx_int64_t step, real &weight, real largest_fraction_of_f)
{
    if (do_this_step(step))
    {
        for (auto && group : impl_->atom_groups_)
        {
            real max_ff_force  = group->max_element(f);
            real max_set_force = group->max_set_f();
            if (mpi_helper() != nullptr)
            {
                max_ff_force  = mpi_helper()->max(max_ff_force);
                max_set_force = mpi_helper()->max(max_set_force);
            }
            if ( (max_ff_force != 0) && (max_set_force != 0) && max_set_force / max_ff_force  > largest_fraction_of_f)
            {
                weight *= largest_fraction_of_f * max_ff_force / max_set_force;
            }
            group->add_forces(f, weight);
        }
    }
};


void ExternalPotential::set_atom_properties(t_mdatoms * mdatoms, gmx_localtop_t * topology_loc)
{
    for (auto &group : impl_->atom_groups_)
    {
        for (auto &atom : *group)
        {
            atom.properties = single_atom_properties(mdatoms + atom.i_local, topology_loc);
        }
    }
}

void
ExternalPotential::add_wholemoleculegroup(std::shared_ptr<WholeMoleculeGroup> group)
{
    impl_->wholemoleculegroups_.push_back(group);
}


void ExternalPotential::add_group(std::shared_ptr<Group> group)
{
    impl_->atom_groups_.push_back(group);
}


void ExternalPotential::set_local_virial(tensor virial)
{
    copy_mat(virial, impl_->virial_);
}

void ExternalPotential::set_local_potential(real potential)
{
    impl_->potential_ = potential;
}

std::shared_ptr<Group> ExternalPotential::group(const rvec x[], int group_index)
{
    impl_->atom_groups_[group_index]->set_x(x);
    return impl_->atom_groups_[group_index];
}

std::shared_ptr<WholeMoleculeGroup> ExternalPotential::wholemoleculegroup(const rvec x[], const matrix box, int group_index)
{
    impl_->wholemoleculegroups_[group_index]->set_x(x);
    impl_->wholemoleculegroups_[group_index]->set_box(box);
    return impl_->wholemoleculegroups_[group_index];
}

void ExternalPotential::set_mpi_helper(std::shared_ptr<MpiHelper> mpi)
{
    impl_->mpi_ = mpi;
}

void ExternalPotential::set_input_output(std::shared_ptr<ExternalPotentialIO> &&input_output)
{
    impl_->input_output_ = std::move(input_output);
}

std::shared_ptr<ExternalPotentialIO> ExternalPotential::input_output()
{
    return impl_->input_output_;
}

std::shared_ptr<MpiHelper> ExternalPotential::mpi_helper()
{
    return impl_->mpi_;
};

} // end namespace externalpotential

} // end namespace gmx
