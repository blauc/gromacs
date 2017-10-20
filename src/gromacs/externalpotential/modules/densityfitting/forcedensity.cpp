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
#include "gromacs/externalpotential/atomgroups/wholemoleculegroup.h"
#include "gromacs/math/gmxcomplex.h"
#include "gromacs/math/griddata/operations/fouriertransform.h"
#include "gromacs/math/griddata/field.h"
#include <algorithm>

/******************************************************************************
 * ForceDensity
 */


namespace gmx
{


ForceDensity::ForceDensity(const FieldReal3D &grid, real sigma)
    : sigma_ {sigma}, forces_ {{
                                   grid, grid, grid
                               }},
forcesFT_ {{
               convertGridToReciprocalSpace(*(grid.getGrid().duplicate())), convertGridToReciprocalSpace(*(grid.getGrid().duplicate())), convertGridToReciprocalSpace(*(grid.getGrid().duplicate()))
           }},
convolutionDensity_ {{
                         convertGridToReciprocalSpace(*(grid.getGrid().duplicate())), convertGridToReciprocalSpace(*(grid.getGrid().duplicate())), convertGridToReciprocalSpace(*(grid.getGrid().duplicate()))
                     }},
densityGradientFT_ {
    convertGridToReciprocalSpace(*(grid.getGrid().duplicate()))
},
realToComplexFT_ {
    FourierTransformRealToComplex3D(grid)
} {
    generateFourierTransformGrids_(*(grid.getGrid().duplicate()));
    generateConvolutionDensity_();
    realToComplexFT_.normalize();
};

void ForceDensity::generateFourierTransformGrids_(
        const IGrid<DIM> &realSpaceGrid)
{
    auto fourierSpaceGrid = convertGridToReciprocalSpace(realSpaceGrid);
    fourierSpaceGrid->setLatticeAndRescaleCell(fourierTransformGridExtendfromRealExtend(realSpaceGrid.lattice().getExtend()));

    for (auto &fourierField : forcesFT_)
    {
        fourierField.setGrid(fourierSpaceGrid->duplicate());
        complexToRealFTarray_.emplace_back(fourierField);
        complexToRealFTarray_.back().normalize();
    }

    for (auto &convolutionGrid : convolutionDensity_)
    {
        convolutionGrid.setGrid(fourierSpaceGrid->duplicate());
    }

    densityGradientFT_.setGrid(fourierSpaceGrid->duplicate());
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
        ApplyToUnshiftedFourierTransform(convolutionDensity_[dimension])
            .apply(gaussianTimesK);
    }
}

const std::array<FieldReal3D, DIM> &
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

        std::transform(densityGradientFT_.begin(),
                       densityGradientFT_.end(),
                       convolutionDensity_[dimension].begin(),
                       forcesFT_[dimension].begin(), iWGaussianK);
        complexToRealFTarray_[dimension].result(forces_[dimension]);

    }
    return forces_;
};

} // gmx
