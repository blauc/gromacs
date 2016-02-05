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
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/groupcoord.h"
#include "gromacs/mdlib/sim_util.h"

#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/smalloc.h"
#include "group.h"
#include "mpi-helper.h"
#include "externalpotentialIO.h"

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

        real reduce_potential_virial();
        void dd_make_local_groups( gmx_ga2la_t *ga2la);
        // const std::vector<RVec> &x_assembled(const gmx_int64_t step, t_commrec *cr, const rvec x[], const matrix box, int group_index);
        GroupCoordinates x_local(const rvec x[], int group_index);
        std::vector<RVec> &f_local(int group_index);
        const std::vector<int> &collective_index(int group_index);
        bool do_this_step(int step);
        void add_forces(rvec f[], real w);
        void add_virial(tensor vir, real total_weight);
        void set_potential(real potential);
        void set_virial(tensor virial);
        void add_group(std::unique_ptr<Group> &&group);
        void set_mpi_helper(std::shared_ptr<MpiHelper> mpi);
        void set_input_output(std::shared_ptr<ExternalPotentialIO> &&input_output);
        std::shared_ptr<ExternalPotentialIO> input_output();

    private:

        gmx_int64_t                          nst_apply_; /**< every how many steps to apply; \TODO for every step, since not implemented  */
        real                                 potential_; /**< the (local) contribution to the potential */
        tensor                               virial_;

        std::shared_ptr<MpiHelper>           mpi_;
        std::shared_ptr<ExternalPotentialIO> input_output_;
        std::vector < std::unique_ptr < Group>> atom_groups_; /**< \TODO: make groups friend class, and set x[], cr, and f[] here to avoid handing down parameters all the time.?*/
};

std::shared_ptr<ExternalPotentialIO> ExternalPotential::Impl::input_output()
{
    return input_output_;
}

std::vector<RVec> &ExternalPotential::Impl::f_local(int group_index)
{
    return atom_groups_[group_index]->f_local();
}

void ExternalPotential::Impl::set_virial(tensor virial)
{
    copy_mat(virial, virial_);

}

void ExternalPotential::Impl::set_potential(real potential)
{
    potential_ = potential;
}

void ExternalPotential::Impl::dd_make_local_groups(gmx_ga2la_t *ga2la)
{
    for (auto && grp : atom_groups_)
    {
        grp->set_indices(ga2la);
    }
}

void ExternalPotential::Impl::set_mpi_helper(std::shared_ptr<MpiHelper> mpi)
{
    mpi_ = mpi;
}

void ExternalPotential::Impl::set_input_output(std::shared_ptr<ExternalPotentialIO> &&input_output)
{
    input_output_ = std::move(input_output);
}


void ExternalPotential::Impl::add_virial(real (*vir)[3], real total_weight)
{
    for (int i = XX; i <= ZZ; i++)
    {
        for (int j = XX; j <= ZZ; j++)
        {
            vir[i][j] += total_weight*virial_[i][j];
        }
    }
}

ExternalPotential::Impl::Impl() :
    nst_apply_(1), potential_(0), mpi_(nullptr)
{
};

void ExternalPotential::Impl::add_group(std::unique_ptr<Group> &&group)
{
    atom_groups_.push_back(std::move(group));
}

ExternalPotential::Impl::~Impl()
{

};

void ExternalPotential::Impl::add_forces(rvec f[], real total_weight)
{
    for (auto && group : atom_groups_)
    {
        group->add_forces(f, total_weight);
    }
}

bool ExternalPotential::Impl::do_this_step(int step)
{
    return (step % nst_apply_ == 0 );
}

GroupCoordinates ExternalPotential::Impl::x_local(const rvec x[], int group_index)
{
    return atom_groups_[group_index]->x_local(x);
}

real ExternalPotential::Impl::reduce_potential_virial()
{
    if (mpi_ != nullptr)
    {
        mpi_->buffer(virial_);
        mpi_->buffer(potential_);

        mpi_->sum_reduce();

        mpi_->buffer(virial_);
        mpi_->buffer(potential_);

        mpi_->finish();

    }

    if (mpi_ == nullptr || mpi_->isMaster())
    {
        return potential_;
    }

    return 0.0;

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
    impl_->dd_make_local_groups(ga2la);
};

void ExternalPotential::add_virial(tensor vir, gmx_int64_t step, real weight)
{
    if (impl_->do_this_step(step))
    {
        impl_->add_virial(vir, weight);
    }
};


real ExternalPotential::sum_reduce_potential_virial()
{
    return impl_->reduce_potential_virial();
}

void ExternalPotential::add_forces( rvec f[], gmx_int64_t step, real weight)
{
    if (impl_->do_this_step(step))
    {
        impl_->add_forces(f, weight);
    }
};

void ExternalPotential::add_group(std::unique_ptr<Group> &&group)
{
    impl_->add_group(std::move(group));
}


void ExternalPotential::set_local_virial(tensor virial)
{
    impl_->set_virial(virial);
}

void ExternalPotential::set_local_potential(real potential)
{
    impl_->set_potential(potential);
}

std::vector<RVec> &ExternalPotential::f_local(int group_index)
{
    return impl_->f_local(group_index);
}

GroupCoordinates ExternalPotential::x_local(const rvec x[], int group_index)
{
    return impl_->x_local(x, group_index);
}

void ExternalPotential::set_mpi_helper(std::shared_ptr<MpiHelper> mpi)
{
    impl_->set_mpi_helper(mpi);
}

void ExternalPotential::set_input_output(std::shared_ptr<ExternalPotentialIO> &&input_output)
{
    impl_->set_input_output(std::move(input_output));
}

std::shared_ptr<ExternalPotentialIO> ExternalPotential::input_output()
{
    return impl_->input_output();
}

} // end namespace externalpotential

} // end namespace gmx
