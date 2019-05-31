/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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

#include "densityfittingparameters.h"

namespace gmx
{

DensityFittingParameters::DensityFittingParameters(const LocalAtomSet                   &atomSet,
                                                   const TranslateAndScale              &transformationToDensityLattice,
                                                   real                                  forceConstant,
                                                   real                                  sigma,
                                                   double                                nSigma,
                                                   reference_density                     referenceDensity,
                                                   const DensityFittingAmplitudeMethod  &amplitudeMethod,
                                                   const DensitySimilarityMeasureMethod &measureMethod
                                                   ) : atomSet_ {atomSet}, forceConstant_ {
    forceConstant
},
transformationToDensityLattice_ {
    transformationToDensityLattice
},
referenceDensity_ {
    referenceDensity
},
sigma_ {
    sigma
}, nSigma_ {
    nSigma
}, amplitudeMethod_ {
    amplitudeMethod
},
measureMethod_ {
    measureMethod
}
{
}

DensitySimilarityMeasure DensityFittingParameters::makeMeasure() const
{
    return DensitySimilarityMeasure(measureMethod_, referenceDensity_);
}

const LocalAtomSet &DensityFittingParameters::atomSet() const
{
    return atomSet_;
}

real DensityFittingParameters::forceConstant() const
{
    return forceConstant_;
}

DensityFittingForce DensityFittingParameters::makeForceEvaluator() const
{
    return DensityFittingForce(makeSpreadKernel());
}

DensityFittingAmplitudeLookup DensityFittingParameters::makeAmplitudeLookup() const
{
    return DensityFittingAmplitudeLookup(amplitudeMethod_);
}

const TranslateAndScale &DensityFittingParameters::transformationToDensityLattice() const
{
    return transformationToDensityLattice_;
}

GaussianSpreadKernelParameters::Shape DensityFittingParameters::makeSpreadKernel() const
{
    RVec sigmaInLatticeCoordinates {
        sigma_, sigma_, sigma_
    };
    transformationToDensityLattice_.scaleOperationOnly() ({&sigmaInLatticeCoordinates, &sigmaInLatticeCoordinates + 1});
    return { {
                 sigmaInLatticeCoordinates[XX], sigmaInLatticeCoordinates[YY], sigmaInLatticeCoordinates[ZZ]
             }, nSigma_};
}

GaussTransform3D DensityFittingParameters::makeSpreadingTransform() const
{
    return {referenceDensity_.extents(), makeSpreadKernel()};
}

} // namespace gmx
