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
#ifndef GMX_EXTERNALPOTENIAL_FORCEDENSITY_H
#define GMX_EXTERNALPOTENIAL_FORCEDENSITY_H

#include "gmxpre.h"
#include "gromacs/math/gmxcomplex.h"
#include "gromacs/math/griddata/field.h"
#include "gromacs/math/griddata/operations/fouriertransform.h"
#include "gromacs/utility/real.h"

#include <array>
#include <vector>


namespace gmx
{


class ForceDensity
{
    public:
        ForceDensity(const Field<real>             &grid,
                     real                           sigma);
        ~ForceDensity() = default;
        const std::array<Field<real>, DIM> &getForce();

    private:
        void generateConvolutionDensity_();
        void generateFourierTransformGrids_(const GridWithTranslation<DIM> &realSpaceGrid);
        real sigma_;
        std::array<Field<real>, DIM>                    forces_;
        std::array<Field<t_complex>, DIM>               forcesFT_;
        std::array<Field<t_complex>, DIM>               convolutionDensity_;
        Field<t_complex>                                densityGradientFT_;
        std::vector<FourierTransformComplexToReal3D>    complexToRealFTarray_;
        FourierTransformRealToComplex3D                 realToComplexFT_;
};
}
#endif /* end of include guard: GMX_EXTERNALPOTENIAL_FORCEDENSITY_H */
