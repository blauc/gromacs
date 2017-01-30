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

#include "crosscorrelation.h"

#include <algorithm>
#include <memory>
#include <string>

#include "gromacs/math/volumedata/field.h"
#include "gromacs/math/volumedata/gridmeasures.h"
#include "gromacs/math/volumedata/gridreal.h"
#include "gromacs/utility/real.h"

namespace gmx
{
namespace volumedata
{

CrossCorrelation::CrossCorrelation() : differential {new (Field<real>)}
{};

const Field<real> &
CrossCorrelation::evaluateDensityDifferential(const Field<real> &comparant,
                                              const Field<real> &reference)
{
    differential->copy_grid(reference);
    auto cc                      = GridMeasures(reference).correlate(comparant);
    auto normSimulation          = GridReal(comparant).properties().norm();
    auto densityGradientFunction = [normSimulation, cc](real densityExperiment,
                                                        real densitySimulation) {
            return (densityExperiment - cc * densitySimulation) / (normSimulation);
        };
    std::transform(reference.access().begin(), reference.access().end(),
                   comparant.access().begin(), differential->access().begin(),
                   densityGradientFunction);
    return *differential;
}

real
CrossCorrelation::evaluateDensityDensityPotential(
        const Field<real> &comparant, const Field<real> &reference,
        const RVec &translation,
        const Quaternion &orientation)
{
    if (norm2(translation) < 1e-10 && orientation.norm() <  1e-10)
    {
        return GridMeasures(comparant).correlate(reference, correlationThreshold_);
    }
    ;

};

real
CrossCorrelation::evaluateStructureDensityPotential(
        const std::vector<RVec> &coordinates, const std::vector<real> &weights,
        const Field<real> &reference, const RVec &translation,
        const Quaternion &orientation)
{

};

real
CrossCorrelation::evaluateGroupDensityPotential(
        const externalpotential::WholeMoleculeGroup &atoms,
        const Field<real> &reference, const RVec &translation,
        const Quaternion &orientation)
{

};

std::string CrossCorrelationDifferentialInfo::name = std::string("cross-correlation");

std::unique_ptr<IDensityDifferentialProvider> CrossCorrelationDifferentialInfo::create()
{
    return std::unique_ptr<CrossCorrelation>(new CrossCorrelation);
}

std::string CrossCorrelationDensityDensityInfo::name = std::string("cross-correlation");

std::unique_ptr<IDensityDensityPotentialProvider> CrossCorrelationDensityDensityInfo::create()
{
    return std::unique_ptr<CrossCorrelation>(new CrossCorrelation);
}

std::string CrossCorrelationStructureDensityInfo::name = std::string("cross-correlation");

std::unique_ptr<IStructureDensityPotentialProvider> CrossCorrelationStructureDensityInfo::create()
{
    return std::unique_ptr<CrossCorrelation>(new CrossCorrelation);
}

} /* volumedata */
} /* gmx */
