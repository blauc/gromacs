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
 * Declares parameters needed to evaluate forces and energies for density fitting
 *
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_applied_forces
 */
#ifndef GMX_APPLIED_FORCES_DENSITYFITTINGPARAMETERS_H
#define GMX_APPLIED_FORCES_DENSITYFITTINGPARAMETERS_H

#include "gromacs/applied_forces/densityfittingamplitudelookup.h"
#include "gromacs/domdec/localatomset.h"
#include "gromacs/math/coordinatetransformation.h"
#include "gromacs/math/densityfit.h"
#include "gromacs/math/densityfittingforce.h"
#include "gromacs/math/gausstransform.h"

namespace gmx
{

class DensityFittingParameters
{
    public:
        using reference_density = basic_mdspan<const float, dynamicExtents3D>;

        DensityFittingParameters(const LocalAtomSet                   &atomSet,
                                 const TranslateAndScale              &transformationToDensityLattice,
                                 real                                  forceConstant,
                                 real                                  sigma,
                                 double                                nSigma,
                                 reference_density                     referenceDensity,
                                 const DensityFittingAmplitudeMethod  &amplitudeMethod,
                                 const DensitySimilarityMeasureMethod &measureMethod,
                                 int                                   everyNSteps,
                                 bool                                  isMaster
                                 );

        DensitySimilarityMeasure makeMeasure() const;
        GaussTransform3D makeSpreadingTransform() const;
        DensityFittingForce makeForceEvaluator() const;
        const LocalAtomSet      &atomSet() const;
        real forceConstant() const;
        const TranslateAndScale &transformationToDensityLattice() const;
        GaussianSpreadKernelParameters::Shape makeSpreadKernel() const;
        DensityFittingAmplitudeLookup makeAmplitudeLookup() const;
        int everyNSteps() const;
        bool isMaster_;

    private:

        LocalAtomSet                   atomSet_;
        real                           forceConstant_;
        TranslateAndScale              transformationToDensityLattice_;
        const reference_density        referenceDensity_;
        std::vector<real>              amplitudes_;
        real                           sigma_;
        double                         nSigma_;
        DensityFittingAmplitudeMethod  amplitudeMethod_;
        DensitySimilarityMeasureMethod measureMethod_;
        int everyNSteps_;
};

}      // namespace gmx

#endif // GMX_APPLIED_FORCES_DENSITYFITTINGPARAMETERS_H
