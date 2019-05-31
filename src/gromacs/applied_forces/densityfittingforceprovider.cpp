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
/*! \internal \file
 * \brief
 * Implements force provider for density fitting
 *
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "densityfittingforceprovider.h"

#include <algorithm>
#include <vector>

#include "gromacs/applied_forces/densityfittingamplitudelookup.h"
#include "gromacs/math/densityfit.h"
#include "gromacs/math/gausstransform.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forceoutput.h"

#include "gromacs/applied_forces/densityfittingparameters.h"

namespace gmx
{
/********************************************************************
 * DensityFittingForceProvider::Impl
 */

class DensityFittingForceProvider::Impl
{
    public:
        //! \copydoc DensityFittingForceProvider(const DensityFittingParameters &parameters)
        Impl(const DensityFittingParameters &parameters);

        void calculateForces(const ForceProviderInput &forceProviderInput, ForceProviderOutput *forceProviderOutput);

        const DensityFittingParameters            &parameters_;
        GaussTransform3D                           gaussTransform_;
        std::unique_ptr<IDensitySimilarityMeasure> measure_;
        DensityFittingForce                        densityFittingForce_;
        DensityFittingAmplitudeLookup              densityFittingAmplitudeLookup_;
        std::vector<RVec>                          transformedCoordinates_; //< the local atom coordinates transformed into the grid coordinate system
        std::vector<RVec>                          forces_;
};

DensityFittingForceProvider::Impl::Impl(const DensityFittingParameters &parameters) : parameters_(parameters),
                                                                                      gaussTransform_(parameters_.makeSpreadingTransform()),
                                                                                      measure_ {parameters_.makeMeasure()},
densityFittingForce_(parameters_.makeForceEvaluator()),
densityFittingAmplitudeLookup_(parameters_.makeAmplitudeLookup()),
transformedCoordinates_(parameters_.atomSet().numAtomsLocal())
{
}

void DensityFittingForceProvider::Impl::calculateForces(const ForceProviderInput &forceProviderInput,
                                                        ForceProviderOutput      *forceProviderOutput)
{
    // do nothing if there are no density fitting atoms on this node
    if (parameters_.atomSet().numAtomsLocal() == 0)
    {
        return;
    }
    forces_.resize(parameters_.atomSet().numAtomsLocal());

    // pick and copy atom coordinates
    std::transform(std::cbegin(parameters_.atomSet().localIndex()),
                   std::cend(parameters_.atomSet().localIndex()),
                   std::begin(transformedCoordinates_),
                   [ &x = forceProviderInput.x_](int index) { return x[index]; });
    const auto amplitudes = densityFittingAmplitudeLookup_(forceProviderInput.mdatoms_, parameters_.atomSet().localIndex());
    // transform local atom coordinates to density grid coordinates
    parameters_.transformationToDensityLattice() (transformedCoordinates_);

    // spread atoms on grid
    gaussTransform_.setZero();
    auto amplitudeIterator = amplitudes.begin();
    for (const auto &r : transformedCoordinates_)
    {
        gaussTransform_.add({r, *amplitudeIterator});
        ++amplitudeIterator;
    }

    // communicate grid
    gaussTransform_.sumReduce(forceProviderInput.cr_);

    // calculate grid derivative
    const auto &densityDerivative = measure_->gradient(gaussTransform_.view());
    // calculate forces
    std::transform(
            std::begin(transformedCoordinates_),
            std::end(transformedCoordinates_),
            std::begin(amplitudes),
            std::begin(forces_),
            [densityDerivative, this](const RVec r, real amplitude)
            {
                return densityFittingForce_.evaluateForce({r, amplitude}, densityDerivative);
            }
            );

    parameters_.transformationToDensityLattice().scaleOperationOnly().inverseIgnoringZeroScale(forces_);

    auto densfitForceIterator = forces_.cbegin();
    for (const auto localAtomIndex : parameters_.atomSet().localIndex())
    {
        forceProviderOutput->forceWithVirial_.force_[localAtomIndex] +=
            parameters_.forceConstant() * *densfitForceIterator;
        ++densfitForceIterator;
    }
    // calculate corresponding potential energy
    forceProviderOutput->enerd_.term[F_COM_PULL] +=
        parameters_.forceConstant() * measure_->similarity(gaussTransform_.view());
}

/********************************************************************
 * DensityFittingForceProvider
 */

DensityFittingForceProvider::~DensityFittingForceProvider()
{}

DensityFittingForceProvider::DensityFittingForceProvider(const DensityFittingParameters &parameters) :
    impl_(new Impl(parameters))
{}

void DensityFittingForceProvider::calculateForces(const ForceProviderInput  &forceProviderInput,
                                                  ForceProviderOutput      * forceProviderOutput)
{
    impl_->calculateForces(forceProviderInput, forceProviderOutput);
}

} // namespace gmx
