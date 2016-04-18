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
#include "gmxpre.h"

#include "densityfitting.h"

#include <cstdio>

#include <algorithm>
#include <ios>

#include "gromacs/externalpotential/forceplotter.h"
#include "gromacs/externalpotential/externalpotentialIO.h"
#include "gromacs/externalpotential/forceplotter.h"
#include "gromacs/externalpotential/group.h"
#include "gromacs/externalpotential/mpi-helper.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/math/volumedata.h"
#include "gromacs/fileio/volumedataio.h"
#include "ifgt/Ifgt.h"

namespace gmx
{

namespace externalpotential
{

std::unique_ptr<ExternalPotential> DensityFitting::create()
{
    return std::unique_ptr<ExternalPotential> (new DensityFitting());
}

gmx::AtomProperties * DensityFitting::single_atom_properties(t_mdatoms * mdatoms, gmx_localtop_t * topology_loc)
{
    (void) mdatoms;
    (void) topology_loc;
    return nullptr;
}

void DensityFitting::inv_mul(std::vector<real> &to_invert, const std::vector<real> &multiplier)
{
    auto q = to_invert.begin();
    for (auto p : multiplier)
    {
        *q = p/(*q);
        ++q;
    }
}

RVec DensityFitting::pbc_dist(const rvec x, const rvec y, const matrix box)
{

    ivec im;
    real min_dist = 1e10;
    rvec dist;
    rvec curr_im;
    rvec shift;
    rvec y_shifted;
    RVec result;
    for (im[XX] = -1; im[XX] <= 1; ++im[XX])
    {
        for (im[YY] = -1; im[YY] <= 1; ++im[YY])
        {
            for (im[ZZ] = -1; im[ZZ] <= 1; ++im[ZZ])
            {
                curr_im[XX] = im[XX];
                curr_im[YY] = im[YY];
                curr_im[ZZ] = im[ZZ];
                mvmul(box, curr_im, shift);
                rvec_add(y, shift, y_shifted);
                rvec_sub(y_shifted, x, dist);
                if (norm2(dist) < min_dist*min_dist)
                {
                    result   = dist;
                    min_dist = norm(dist);
                }
            }
        }
    }
    return result;
}

void DensityFitting::do_force_plain(const rvec x, rvec force, const matrix box)
{
    clear_rvec(force);
    RVec d;

    externalpotential::ForcePlotter plot;
    plot.start_plot_forces("plain_forces.bild");

    for (size_t i = 0; i < simulated_density_->num_gridpoints(); i++)
    {
        d = pbc_dist(simulated_density_->gridpoint_coordinate(i), x, box);
        svmul(1/difference_transform_->sigma(), d, d);
        svmul(-k_*simulated_density_->grid_cell_volume()*simulated_density_->data()[i]*exp(-norm2(d)/2.0), d, d);
        plot.plot_force(simulated_density_->gridpoint_coordinate(i), d, 0, 1e-1);
        rvec_add(force, d, force);
    }
    plot.stop_plot_forces();
}

void DensityFitting::do_potential( const matrix box, const rvec x[], const gmx_int64_t step)
{
    real pot       = 0;
    gauss_transform_->set_box(box);
    difference_transform_->set_box(box);
    volumedata::MrcFile simulated_output;
    /* compress local atom contributions to density in expansion centers*/
    for (auto atom : *group(x, 0))
    {
        gauss_transform_->compress(atom.x, 1);
    }

    gauss_transform_->sum_reduce();

    if (mpi_helper() == nullptr || mpi_helper()->isMaster())
    {
        gauss_transform_->expand(*simulated_density_);
        simulated_output.write_with_own_meta("spread.ccp4", *simulated_density_, meta_, true);
        simulated_density_->normalize();
        simulated_density_->add_offset(background_offset_);
        simulated_output.write_with_own_meta("target.ccp4", *target_density_, meta_, true);
        if (step == 0)
        {
            reference_density_ = simulated_density_->data();
        }
        pot                      = potential_reference_sum_-k_ * potential(target_density_->data(), simulated_density_->data());
        potential_reference_sum_ = pot;
        reference_density_       = simulated_density_->data();

        inv_mul(simulated_density_->data(), target_density_->data());
        simulated_output.write_with_own_meta("divided.ccp4", *simulated_density_, meta_, true);
        // difference_transform_->compress_density(*simulated_density_);
    }
    else
    {
        pot = 0;
    }

    // difference_transform_->broadcast_expansion();

    externalpotential::ForcePlotter plot;
    plot.start_plot_forces("forces.bild");
    real max_f = -1;
    for (auto atom : *group(x, 0))
    {
        // difference_transform_->do_force(atom.x, atom.force, gauss_transform_->sigma(), *simulated_density_);
        do_force_plain(atom.x, atom.force, box);
        if (norm(atom.force) > max_f)
        {
            max_f = norm(atom.force);
        }
    }
    for (auto atom : *group(x, 0))
    {
        plot.plot_force(atom.x, atom.force, 0, 1.0/max_f);
    }
    plot.stop_plot_forces();
    fprintf(stderr, "\n\nPOTENTIAL: %g \n\n", pot);
    set_local_potential(pot);
};

DensityFitting::DensityFitting() : ExternalPotential(),
                                   target_density_(std::unique_ptr<volumedata::GridReal>(new volumedata::GridReal())),
                                   simulated_density_(std::unique_ptr<volumedata::GridReal>(new volumedata::GridReal())),
                                   gauss_transform_(std::unique_ptr<Ifgt>(new Ifgt())),
                                   difference_transform_(std::unique_ptr<Ifgt>(new Ifgt()))
{
    ;
}


real DensityFitting::potential(std::vector<real> &P, std::vector<real> &Q)
{
    // for numerical stability use a reference density
    if (P.size() != Q.size())
    {
        GMX_THROW(APIError("KL-Divergence calculation requires euqally sized input vectors."));
    }
    potential_contribution_.resize(P.size());
    auto pot         = potential_contribution_.begin();
    auto q           = Q.begin();
    auto q_reference = reference_density_.begin();
    for (auto p : P)
    {
        if ((*q > 0) && (p > 0))
        {
            *pot = p * log(*q/(*q_reference));
        }
        ++q;
        ++q_reference;
        ++pot;
    }
    std::sort(potential_contribution_.begin(), potential_contribution_.end());
    real result = std::accumulate(potential_contribution_.begin(), potential_contribution_.end(), 0.0);
    return simulated_density_->grid_cell_volume()*result;
}

void DensityFitting::read_input()
{
    FILE      * inputfile       = input_output()->input_file();
    char      * line            = nullptr;
    size_t      len             = 0;
    real        sigma           = 0;
    int         expansion_order = 0;
    std::string target_density_name;
    try
    {
        getline(&line, &len, inputfile);
        k_ = strtof(line, nullptr);

        getline(&line, &len, inputfile);
        sigma = strtof(line, nullptr);

        getline(&line, &len, inputfile);
        expansion_order = strtof(line, nullptr);

        getline(&line, &len, inputfile);
        background_offset_ = strtof(line, nullptr);

        gauss_transform_->init(sigma, expansion_order);
        difference_transform_->init(sigma, expansion_order);

        getline(&line, &len, inputfile);
        target_density_name = std::string(line);
        volumedata::MrcFile         target_input_file;
        target_input_file.read_with_meta(target_density_name, *target_density_, meta_);
        target_density_->normalize();
        target_density_->add_offset(background_offset_);
        simulated_density_->copy_grid(*target_density_);
        simulated_density_->resize();
    }
    catch (const std::ios_base::failure &e)
    {
        GMX_THROW(gmx::InvalidInputError("Cannot read from externalpotential : template inputfile."));
    }

}

void DensityFitting::broadcast_internal()
{
    gauss_transform_->set_mpi_helper(mpi_helper());
    difference_transform_->set_mpi_helper(mpi_helper());
    if (mpi_helper() != nullptr)
    {
        mpi_helper()->broadcast(&k_, 1);
        gauss_transform_->broadcast_internal();
        difference_transform_->broadcast_internal();
    }
}

void DensityFitting::finish()
{

}

std::string DensityFittingInfo::name                        = "densityfitting";
std::string DensityFittingInfo::shortDescription            = "calculate forces from difference to target density";
const int   DensityFittingInfo::numberIndexGroups           = 1;
externalpotential::ModuleCreator DensityFittingInfo::create = DensityFitting::create;

} // namespace externalpotential
} // namespace gmx
