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

#include "kullbackleibler.h"
#include <memory>
#include <string>
#include "gromacs/math/volumedata/gridreal.h"
namespace gmx
{
namespace volumedata
{

const Field<real> &KullbackLeibler::evaluateDensityDifferential(
        const Field<real> &comparant, const Field<real> &reference)
{
    differential.copy_grid(reference);
    auto sumSimulatedDensity     = GridReal(comparant).properties().sum();
    auto densityGradientFunction = [sumSimulatedDensity](real densityExperiment,
                                                         real densitySimulation) {
            return (densitySimulation > 1e-15 && densityExperiment > 1e-15)
                   ? (densityExperiment / densitySimulation) *
                   (1 - densitySimulation / sumSimulatedDensity)
                   : 0;
        };

    std::transform(reference.access().begin(), reference.access().end(),
                   comparant.access().begin(), differential.access().begin(),
                   densityGradientFunction);
    return differential;
}

const std::string KullbackLeiblerInfo::name = std::string("kullback-leibler");

std::unique_ptr<IDensityDifferentialProvider>
KullbackLeiblerInfo::create()
{
    return std::unique_ptr<KullbackLeibler>(
            new KullbackLeibler);
}

} /* volumedata */
} /* gmx */
