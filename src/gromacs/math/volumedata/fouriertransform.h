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
 * Defines volume data containers.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \inpublicapi
 */
#ifndef GMX_MATH_FOURIERTRANSFORM_H
#define GMX_MATH_FOURIERTRANSFORM_H

#include "field.h"
#include "gromacs/math/gmxcomplex.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/real.h"
#include <fftw3.h>
#include <memory>
namespace gmx
{
namespace volumedata
{

class FourierTransform3D
{
    public:
        FourierTransform3D()  = default;
        ~FourierTransform3D() = default;
        IVec columnMajorExtendToRowMajorExtend(IVec extend);
    protected:
        fftwf_plan_s * plan_;
};

class FourierTransformRealToComplex3D : FourierTransform3D
{
    public:
        FourierTransformRealToComplex3D(const Field<real> &realInputField);
        FourierTransformRealToComplex3D &normalize();
        std::unique_ptr < Field < t_complex>> result();

    private:
        const Field<real> &realInputField_;
        std::unique_ptr < Field < t_complex>> complexTransformedField_;
};

class FourierTransformComplexToReal3D : FourierTransform3D
{
    public:
        FourierTransformComplexToReal3D(const Field<t_complex> &complexInputField);

        FourierTransformComplexToReal3D &normalize();
        std::unique_ptr < Field < real>> result();

    private:
        const Field<t_complex> &complexInputField_;
        std::unique_ptr < Field < real>> realTransformedField_;
};

class ApplyToUnshiftedFourierTransform
{
    public:
        typedef const std::function<void(t_complex &, RVec)> &FunctionOnComplexField;
        ApplyToUnshiftedFourierTransform(Field<t_complex> &field);
        const ApplyToUnshiftedFourierTransform &
        apply(FunctionOnComplexField appliedFunction);

    private:
        void applyToAllColumnsWithinSection_(
            RVec coordinateSectionBegin, int nRowsPerColumn, int nColumns,
            RVec deltakRow, RVec deltakColumn,
            std::vector<t_complex>::iterator itColumnStart,
            std::vector<t_complex>::iterator itColumnEnd,
            FunctionOnComplexField appliedFunction);

        void applyToAllRowsWithinColumn_(RVec coordinateColumnBegin, RVec deltakRow,
                                         std::vector<t_complex>::iterator itRowStart,
                                         std::vector<t_complex>::iterator itRowEnd,
                                         FunctionOnComplexField appliedFunction);
        const Field<t_complex> &field_;
};
}
}
#endif /* end of include guard: GMX_MATH_FOURIERTRANSFORM_H */
