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
//
// #include "fouriershellcorrelation.h"
// #include <memory>
// #include <string>
//
// #include "../densityspreader.h"
// #include "gromacs/fileio/json.h"
// #include "gromacs/math/volumedata/field.h"
// #include "gromacs/math/volumedata/fouriershellcorrelation.h"
// #include "gromacs/math/volumedata/gridinterpolator.h"
// #include "gromacs/math/volumedata/gridmeasures.h"
// #include "gromacs/math/volumedata/gridreal.h"
// #include "gromacs/utility/gmxomp.h"
// #include "gromacs/utility/real.h"
//
// namespace gmx
// {
//
// const Field<real> &FourierShellCorrelationProvider::evaluateDensityDifferential(
//         const Field<real> & /*comparant*/, const Field<real> &reference)
// {
//     differential->copy_grid(reference);
//     // TODO: everything
//     return *differential;
// }
//
// real FourierShellCorrelationProvider::evaluateDensityDensityPotential(
//         const Field<real> &comparant, const Field<real> &reference,
//         const RVec &translation, const Quaternion &orientation)
// {
//     if (!comparant.sameGridInAbsTolerance(reference, 1e-10) &&
//         (norm(translation) > 1e-10) && orientation.norm() > 1e-10)
//     {
//         auto centerOfMass = GridReal(comparant).center_of_mass();
//         auto interpolated = GridInterpolator(reference).interpolateLinearly(
//                     comparant, translation, centerOfMass, orientation);
//         auto fscCurve = FourierShellCorrelation(reference).getFscCurve(
//                     *interpolated, reference);
//         return std::accumulate(std::begin(fscCurve), std::end(fscCurve), 0.) /
//                real(fscCurve.size());
//     }
//     auto fscCurve =
//         FourierShellCorrelation(reference).getFscCurve(comparant, reference);
//     return std::accumulate(std::begin(fscCurve), std::end(fscCurve), 0.) /
//            real(fscCurve.size());
//     ;
// };
//
//
// /**************************INFO CLasses****************************************/
//
//
// std::string FourierShellCorrelationProviderDensityDensityInfo::name =
//     std::string("fsc");
// std::string FourierShellCorrelationProviderDifferentialPotentialInfo::name = FourierShellCorrelationProviderDensityDensityInfo::name;
//
// std::unique_ptr<IDensityDensityPotentialProvider>
// FourierShellCorrelationProviderDensityDensityInfo::create()
// {
//     return std::unique_ptr<FourierShellCorrelationProvider>(
//             new FourierShellCorrelationProvider);
// }
//
// std::unique_ptr<IDifferentialPotentialProvider>
// FourierShellCorrelationProviderDifferentialPotentialInfo::create()
// {
//     return std::unique_ptr<FourierShellCorrelationProvider>(
//             new FourierShellCorrelationProvider);
// }
//
// } /* gmx */
