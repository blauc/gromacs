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

#include "gromacs/math/volumedata/field.h"
#include "gromacs/math/volumedata/gridreal.h"
#include "gromacs/math/volumedata/volumedata.h"

#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/real.h"

namespace gmx
{
namespace volumedata
{

densityBasedPotential::densityBasedPotential(const DensitySpreader &spreader,
                                             bool                   selfSpreading)
    : spreader_ {spreader}, selfSpreading_ {
    selfSpreading
} {}

real densityBasedPotential::potential(const std::vector<RVec> &coordinates,
                                      const std::vector<real> &weights,
                                      const Field<real>       &reference,
                                      const RVec              &translation,
                                      const Quaternion        &orientation,
                                      const RVec              &centerOfRotation)
{
    if (selfSpreading_)
    {

        return densityDensityPotential(
                spreader_.spreadLocalAtoms(coordinates, weights, translation,
                                           orientation, centerOfRotation),
                reference);
    }
    else
    {
        return densityDensityPotential(spreader_.getSpreadGrid(), reference);
    }
};

densityBasedForce::densityBasedForce(const DensitySpreader &spreader,
                                     real sigma_differential, int n_threads,
                                     bool selfSpreading)
    : sigma_differential_ {sigma_differential},
n_threads_ {
    n_threads
}, spreader_ {
    spreader
}, selfSpreading_ {
    selfSpreading
} {};

std::vector<RVec> &densityBasedForce::force(
        const std::vector<RVec> &coordinates, const std::vector<real> &weights,
        const Field<real> &reference, const RVec &translation,
        const Quaternion &orientation, const RVec &centerOfRotation)
{
    if (selfSpreading_)
    {
        setDensityDifferential(spreader_.getSpreadGrid(), reference);
    }
    else
    {
        setDensityDifferential(spreader_.spreadLocalAtoms(coordinates, weights,
                                                          translation, orientation,
                                                          centerOfRotation),
                               reference);
    }
    auto forceDensity =
        ForceDensity(*differential_, sigma_differential_).getForce();
    std::array<volumedata::GridReal, 3> forceGrid {
        {
            volumedata::GridReal(forceDensity[XX]),
            volumedata::GridReal(forceDensity[YY]),
            volumedata::GridReal(forceDensity[ZZ])
        }
    };
    // real          prefactor  = k_/(norm_simulated_*sigma_*sigma_);

    for (size_t i = 0; i != coordinates.size(); ++i)
    {
        auto r = orientation.shiftedAndOriented(coordinates[i], centerOfRotation,
                                                translation);
        auto f = RVec {
            forceGrid[XX].getLinearInterpolationAt(r),
            forceGrid[YY].getLinearInterpolationAt(r),
            forceGrid[ZZ].getLinearInterpolationAt(r)
        };
        svmul(weights[i], f, force_[i]);
        orientation.rotate_backwards(force_[i]);
    }
    return force_;
}
//
// void densityBasedPotentialForce::force(WholeMoleculeGroup &atoms, const
// Field<real> &reference, const RVec &translation, const Quaternion
// &orientation,
//                                        const RVec &centerOfRotation )
// {
//
//     setDensityDifferential(*spread_density_, reference);
//     auto forceDensity = ForceDensity(*differential_,
//     sigma_differential_).getForce();
//     std::array<volumedata::GridReal, 3> forceGrid { {
//                                                         volumedata::GridReal(forceDensity[XX]),
//                                                         volumedata::GridReal(forceDensity[YY]),
//                                                         volumedata::GridReal(forceDensity[ZZ])
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

} /* volumedata */
} /* gmx */
