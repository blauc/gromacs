/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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

#include "densitybasedpotential.h"
#include <vector>

#include "gromacs/externalpotential/atomgroups/group.h"
#include "gromacs/externalpotential/atomgroups/wholemoleculegroup.h"
#include "gromacs/externalpotential/modules/densityfitting/densityspreader.h"
#include "gromacs/externalpotential/modules/densityfitting/forcedensity.h"

#include "gromacs/math/quaternion.h"
#include "gromacs/math/vectypes.h"

#include "gromacs/math/griddata/griddata.h"

#include "gromacs/math/griddata/operations/gridinterpolator.h"

#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/real.h"

namespace gmx
{

DensityBasedPotential::DensityBasedPotential(const DensitySpreader &spreader)
    : spreader_ {spreader}
{}

real DensityBasedPotential::potential(const GridDataReal3D       &reference) const
{
    return densityDensityPotential(spreader_.getSpreadGrid(), reference);
};

DensityBasedForce::DensityBasedForce(const DensitySpreader &spreader,
                                     real sigma_differential, int n_threads
                                     )
    : sigma_differential_ {sigma_differential},
n_threads_ {
    n_threads
}, spreader_ {
    spreader
} {};

void DensityBasedForce::force(std::vector<RVec>    &force,
                              const GridDataReal3D &reference) const
{
    densityDifferential(spreader_.getSpreadGrid(), reference);
    const auto &forceGrid =
        ForceDensity(differential_, sigma_differential_).getForce();
    // real          prefactor  = k_/(norm_simulated_*sigma_*sigma_);
    for (size_t i = 0; i != coordinates_->size(); ++i)
    {
        auto f = RVec {
            getLinearInterpolationAt(forceGrid[XX], coordinates_[i]),
            getLinearInterpolationAt(forceGrid[YY], coordinates_[i]),
            getLinearInterpolationAt(forceGrid[ZZ], coordinates_[i])
        };
        svmul(weights_[i], f, force[i]);
    }
}
//
// void DensityBasedPotentialForce::force(WholeMoleculeGroup &atoms, const
// GridDataReal3D &reference, const RVec &translation, const Quaternion
// &orientation,
//                                        const RVec &centerOfRotation )
// {
//
//     setDensityDifferential(*spread_density_, reference);
//     auto forceDensity = ForceDensity(*differential_,
//     sigma_differential_).getForce();
//     std::array<GridDataReal3D, 3> forceGrid { {
//                                                         GridDataReal3D(forceDensity[XX]),
//                                                         GridDataReal3D(forceDensity[YY]),
//                                                         GridDataReal3D(forceDensity[ZZ])
//                                                     } };
//
// #pragma omp parallel num_threads(n_threads_) shared(stderr,atoms,forceGrid,
// translation, centerOfRotation, orientation) default(none)
//     {
//         int           thread             = gmx_omp_get_thread_num();
//         // real          prefactor  = k_/(norm_simulated_*sigma_*sigma_);
//
//         auto beginThreadAtoms = atoms.begin(thread, n_threads_);
//         auto endThreadAtoms   = atoms.end(thread, n_threads_);
//
//         for (auto atom = beginThreadAtoms; atom != endThreadAtoms; ++atom)
//         {
//             auto r = orientation.shiftedAndOriented(*(*atom).xTransformed,
//             centerOfRotation, translation);
//             auto f = RVec {
//                 forceGrid[XX].getLinearInterpolationAt(r),
//                 forceGrid[YY].getLinearInterpolationAt(r),
//                 forceGrid[ZZ].getLinearInterpolationAt(r)
//             };
//             svmul(*((*atom).properties), f, *(*atom).force);
//             orientation.rotate_backwards(*(*atom).force);
//         }
//     }
// }

} /* gmx */
