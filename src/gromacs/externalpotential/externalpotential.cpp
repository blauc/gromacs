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

namespace gmx
{

/******************************************************************************
 * ExternalPotential::Impl
 */

class ExternalPotential::Impl
{
    public:
        Impl(struct ext_pot_ir * ep_ir, t_commrec * cr, t_inputrec * ir,
             const gmx_mtop_t* mtop, const rvec x[], matrix box, FILE *input_file,
             FILE *output_file, FILE *fplog, bool bVerbose,
             const gmx_output_env_t *oenv, unsigned long Flags, int number_groups);

        ~Impl();

        real reduce_potential_virial(t_commrec * cr);
        void dd_make_local_groups( gmx_domdec_t *dd);
        const std::vector<RVec> &x_assembled(gmx_int64_t step, t_commrec *cr, const rvec x[], matrix box, int group_index);
        GroupCoordinates x_local(const rvec x[], int group_index);
        std::vector<RVec> &f_local(int group_index);
        const std::vector<int> &collective_index(int group_index);
        bool do_this_step(int step);
        void add_forces(rvec f[], real w);
        void add_virial(tensor vir, real total_weight);
        void set_potential(real potential);
        void set_virial(tensor virial);
        FILE* input_file();
        FILE* output_file();
        FILE* log_file();
        bool bVerbose();
        unsigned long Flags();
        const gmx_output_env_t *oenv();
        bool isVerbose();

    private:

        gmx_int64_t             nst_apply_; /**< every how many steps to apply; \todo for every step, since not implemented  */
        real                    potential_; /**< the (local) contribution to the potential */
        tensor                  virial_;

        FILE                   *input_file_;
        FILE                   *output_file_;
        FILE                   *log_file_;

        bool                    bVerbose_; /**< -v flag from command line                      */
        unsigned long           Flags_;
        const gmx_output_env_t *oenv_;
        MpiHelper               mpi;

        std::vector < std::unique_ptr < Group>> atom_groups_;
};

FILE*  ExternalPotential::Impl::input_file()
{
    return input_file_;
}

FILE* ExternalPotential::Impl::output_file()
{
    return output_file_;
}

FILE* ExternalPotential::Impl::log_file()
{
    return log_file_;
}

unsigned long ExternalPotential::Impl::Flags()
{
    return Flags_;
}

const gmx_output_env_t * ExternalPotential::Impl::oenv()
{
    return oenv_;
}

bool ExternalPotential::Impl::isVerbose()
{
    return bVerbose_;
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

void ExternalPotential::Impl::dd_make_local_groups( gmx_domdec_t *dd)
{
    for (auto && grp : atom_groups_)
    {
        grp->set_indices(dd);
    }
}

const std::vector<RVec> &ExternalPotential::Impl::x_assembled(gmx_int64_t step, t_commrec *cr, const rvec x[], matrix box, int group_index)
{
    return atom_groups_[group_index]->x_assembled(step, cr, x, box);
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

ExternalPotential::Impl::Impl(
        struct ext_pot_ir *ep_ir, t_commrec * cr, t_inputrec * ir, const gmx_mtop_t* mtop, const rvec x[], matrix box,
        FILE *input_file, FILE *output_file, FILE *fplog, bool bVerbose,
        const gmx_output_env_t *oenv, unsigned long Flags, int number_groups ) :
    nst_apply_(1), potential_(0), input_file_(input_file), output_file_(output_file), log_file_(fplog),
    bVerbose_(bVerbose), Flags_(Flags), oenv_(oenv), mpi(cr)
{
    if (number_groups != ep_ir->number_index_groups)
    {
        GMX_THROW(gmx::InconsistentInputError("Number of index groups found in tpr file for external potential does not match number of required index groups.\n"
                                              "This should have been caught in grompp."));
    }
    for (int i = 0; i < ep_ir->number_index_groups; i++)
    {
        atom_groups_.push_back(
                std::unique_ptr<Group>(
                        new Group(ir->ePBC, cr, mtop, ep_ir->nat[i], ep_ir->ind[i], x, box)
                        )
                );
    }
};
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

const std::vector<int> &ExternalPotential::Impl::collective_index(int group_index)
{
    return atom_groups_[group_index]->collective_index();
}

real ExternalPotential::Impl::reduce_potential_virial(t_commrec * cr)
{
    if (PAR(cr))
    {
        mpi.cr(cr);

        mpi.buffer(virial_);
        mpi.buffer(potential_);

        mpi.sum_reduce();

        mpi.buffer(virial_);
        mpi.buffer(potential_);

        mpi.finish();

    }

    if (!PAR(cr) || MASTER(cr))
    {
        return potential_;
    }

    return 0.0;

};

/*******************************************************************************
 * ExternalPotential
 */

ExternalPotential::ExternalPotential (struct ext_pot_ir *ep_ir, t_commrec * cr, t_inputrec * ir, const gmx_mtop_t* mtop, const rvec x[], matrix box, FILE *input_file_p, FILE *output_file_p, FILE *fplog, bool bVerbose, const gmx_output_env_t *oenv, unsigned long Flags, int number_groups)
{
    impl_ = std::unique_ptr<ExternalPotential::Impl>(new ExternalPotential::Impl(ep_ir, cr, ir, mtop, x, box, input_file_p, output_file_p, fplog, bVerbose, oenv, Flags, number_groups));
};

ExternalPotential::~ExternalPotential(){};

void ExternalPotential::dd_make_local_groups( gmx_domdec_t *dd)
{
    impl_->dd_make_local_groups(dd);
};

void ExternalPotential::add_virial(tensor vir, gmx_int64_t step, real weight)
{
    if (impl_->do_this_step(step))
    {
        impl_->add_virial(vir, weight);
    }
};

const std::vector<RVec> &ExternalPotential::x_assembled(gmx_int64_t step, t_commrec *cr, const rvec x[], matrix box, int group_index)
{
    return impl_->x_assembled(step, cr, x, box, group_index);
};


real ExternalPotential::sum_reduce_potential_virial(t_commrec * cr)
{
    return impl_->reduce_potential_virial(cr);
}

void ExternalPotential::add_forces( rvec f[], gmx_int64_t step, real weight)
{
    if (impl_->do_this_step(step))
    {
        impl_->add_forces(f, weight);
    }
};

FILE* ExternalPotential::input_file()
{
    return impl_->input_file();
}

FILE* ExternalPotential::output_file()
{
    return impl_->output_file();
}

FILE* ExternalPotential::log_file()
{
    return impl_->log_file();
}



unsigned long ExternalPotential::Flags()
{
    return impl_->Flags();
}

const gmx_output_env_t * ExternalPotential::oenv()
{
    return impl_->oenv();
}

bool ExternalPotential::isVerbose()
{
    return impl_->isVerbose();
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

const std::vector<int> &ExternalPotential::collective_index(int group_index)
{
    return impl_->collective_index(group_index);
}
} // end namespace gmx
