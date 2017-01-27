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
/*!  \file
 * \brief
 * Class definition for force density.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \inpublicapi
 */

#include "forcedensity.h"
#include "gromacs/math/gmxcomplex.h"
#include "gromacs/math/volumedata/fouriertransform.h"
#include <algorithm>

/******************************************************************************
 * ForceDensity
 */
namespace gmx
{

ForceDensity::ForceDensity(const volumedata::Field<real> &grid, real sigma)
    : sigma_ {sigma}, forces_ {{
                                   grid, grid, grid
                               }}, realToComplexFT_ {
    volumedata::FourierTransformRealToComplex3D(grid)
} {
    generateFourierTransformGrids_(grid);
    generateConvolutionDensity_();
    realToComplexFT_.normalize();
};

void ForceDensity::generateFourierTransformGrids_(
        const volumedata::FiniteGrid &grid)
{
    auto fourierExtend =
        volumedata::fourierTransformGridExtendfromRealExtend(grid.extend());
    for (auto &fourierGrid : forcesFT_)
    {
        fourierGrid.set_extend(fourierExtend);
        complexToRealFTarray_.emplace_back(fourierGrid);
        complexToRealFTarray_.back().normalize();
    }
    for (auto &convolutionGrid : convolutionDensity_)
    {
        convolutionGrid.copy_grid(grid);
        convolutionGrid.convertToReciprocalSpace();
        convolutionGrid.set_extend(fourierExtend);
        convolutionGrid.resetCell();
    }
    densityGradientFT_.set_extend(fourierExtend);
}

void ForceDensity::generateConvolutionDensity_()
{
    for (int dimension = XX; dimension <= ZZ; ++dimension)
    {
        auto expPrefactor   = 2 * M_PI * M_PI * sigma_ * sigma_;
        real two_pi         = 2 * M_PI;
        auto gaussianTimesK = [dimension, expPrefactor, two_pi](t_complex &value, RVec k) {
                value.re = two_pi * k[dimension] * exp(-expPrefactor * norm2(k));
                value.im = 0;
            };
        volumedata::ApplyToUnshiftedFourierTransform(convolutionDensity_[dimension])
            .apply(gaussianTimesK);
    }
}

const std::array<volumedata::Field<real>, DIM> &
ForceDensity::getForce()
{
    realToComplexFT_.result(densityGradientFT_);
    for (int dimension = XX; dimension <= ZZ; ++dimension)
    {
        auto iWGaussianK = [](const t_complex &W, const t_complex &gaussianK) {
                return t_complex {
                           -gaussianK.re * W.im, gaussianK.re * W.re
                };
            };

        std::transform(densityGradientFT_.access().begin(),
                       densityGradientFT_.access().end(),
                       convolutionDensity_[dimension].access().begin(),
                       forcesFT_[dimension].access().begin(), iWGaussianK);
        complexToRealFTarray_[dimension].result(forces_[dimension]);

    }
    return forces_;
};
} // gmx
