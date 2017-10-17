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
 * Fourier transform wrapper implementation.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \inpublicapi
 */

#include "fouriertransform.h"
#include "densitypadding.h"
#include "../field.h"
#include "gromacs/fft/parallel_3dfft.h"
#include "gromacs/fileio/griddataio.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/real.h"
#include <algorithm>
#include <cassert>
#include <fftw3.h>

namespace gmx
{

std::array<int, 3> FourierTransform3D::columnMajorExtendToRowMajorExtend(const std::array<int, 3> &extend) const
{
    return { {
                 extend[ZZ], extend[YY], extend[XX]
             } };
}

void
FourierTransform3D::execute()
{
    fftwf_execute(plan_);
}

FourierTransformRealToComplex3D &
FourierTransformRealToComplex3D::normalize()
{
    bNormalize = true;
    return *this;
};

FourierTransformComplexToReal3D &
FourierTransformComplexToReal3D::normalize()
{
    bNormalize = true;
    return *this;
};

FourierTransform3D::~FourierTransform3D() { fftwf_destroy_plan(plan_); }

std::array<int, 3> fourierTransformGridExtendfromRealExtend(const std::array<int, 3> &extend)
{
    return {
               {
                   extend[XX] / 2 + 1, extend[YY], extend[ZZ]
               }
    };
};

std::array<int, 3> realGridExtendFromFourierTransfrom(const std::array<int, 3> &extend)
{
    return {
               {
                   (extend[XX] - 1) * 2, extend[YY], extend[ZZ]
               }
    };
};

std::unique_ptr < Field < t_complex>>
FourierTransformRealToComplex3D::createComplexTransformedFieldFromInput_() const {
    /*
     * set to the size of the real
     * grid first for correct
     * reciprocal basis vectors,
     * then convert to abriged FFT
     * size
     */
    auto complexTransformedGrid = realInputField_.getGrid().reciprocalGrid();
    complexTransformedGrid.setLatticeAndRescaleCell(
            fourierTransformGridExtendfromRealExtend(realInputField_.getGrid().lattice().getExtend()));
    return std::unique_ptr < Field < t_complex>>(new Field<t_complex>(complexTransformedGrid));
};

FourierTransformRealToComplex3D::FourierTransformRealToComplex3D(
        const Field<real> &realInputField)
    : realInputField_ {realInputField}
{

};
std::unique_ptr < Field < t_complex>> FourierTransformRealToComplex3D::result()
{
    auto complexTransformedField = createComplexTransformedFieldFromInput_();
    result(*complexTransformedField);
    return complexTransformedField;
}

void FourierTransformRealToComplex3D::result(
        Field<t_complex> &complexTransformedField)
{

    if (plan_ == nullptr)
    {
        plan_ = fftwf_plan_dft_r2c(
                    3, columnMajorExtendToRowMajorExtend(realInputField_.getGrid().lattice().getExtend()).data(),
                    (float *)realInputField_.data(),
                    (fftwf_complex *)complexTransformedField.data(),
                    FFTW_ESTIMATE);
    }
    fftwf_execute_dft_r2c(
            plan_, (float *)realInputField_.data(),
            (fftwf_complex *)complexTransformedField.data());
    if (bNormalize)
    {
        auto orthoscale = 1.0 / sqrt(complexTransformedField.size());

        auto normalizeComplexTransform = [orthoscale](t_complex &value) {
                value.re *= orthoscale;
                value.im *= orthoscale;
            };

        std::for_each(complexTransformedField.begin(),
                      complexTransformedField.end(),
                      normalizeComplexTransform);
    }
}

std::unique_ptr < Field < real>>
FourierTransformComplexToReal3D::createRealTransformedFieldFromInput_() const {
    auto realGrid = complexInputField_.getGrid();
    realGrid.setLatticeAndRescaleCell(realGridExtendFromFourierTransfrom(complexInputField_.getGrid().lattice().getExtend()));
    realGrid = realGrid.reciprocalGrid();
    return std::unique_ptr < Field < real>>(new Field<real>(realGrid));
}

FourierTransformComplexToReal3D::FourierTransformComplexToReal3D(
        const Field<t_complex> &complexInputField)
    : complexInputField_ {complexInputField}
{};

std::unique_ptr < Field < real>>
FourierTransformComplexToReal3D::result() {
    auto realTransformedField = createRealTransformedFieldFromInput_();
    result(*realTransformedField);
    return realTransformedField;
}

void
FourierTransformComplexToReal3D::result(
        Field<real> &realTransformedField)
{
    if (plan_ == nullptr)
    {
        plan_ = fftwf_plan_dft_c2r(
                    3, columnMajorExtendToRowMajorExtend(realTransformedField.getGrid().lattice().getExtend()).data(),
                    (fftwf_complex *)complexInputField_.data(),
                    (float *)realTransformedField.data(), FFTW_ESTIMATE);
    }
    fftwf_execute_dft_c2r(
            plan_, (fftwf_complex *)complexInputField_.data(),
            (float *)realTransformedField.data());

    if (bNormalize)
    {
        auto orthoscale         = 1 / sqrt(realTransformedField.size());
        auto normalizeTransform = [orthoscale](real &value) {
                value *= orthoscale;
            };
        std::for_each(realTransformedField.begin(),
                      realTransformedField.end(), normalizeTransform);
    }
};

ApplyToUnshiftedFourierTransform::ApplyToUnshiftedFourierTransform(
        Field<t_complex> &field)
    : field_ {field}
{};

const ApplyToUnshiftedFourierTransform &ApplyToUnshiftedFourierTransform::apply(
        FunctionOnComplexField appliedFunction)
{

    /*
     * Assume that second halves of indices denote negative frequencies
     * Use iterators to step though the grid.
     * This relies on the assumption that the grid is stored with z the slowest
     * and x the fasted changing dimension with no padding
     * (x,y,z not being linked to any coordiante system, but short-hand for
     * first, second, third dimension)
     * Loosing generality through this approach, we save substantial time when
     * we don't have to calculate the grid index.
     */

    auto fieldGrid    = field_.getGrid();
    auto extend       = fieldGrid.lattice().getExtend();

    auto k             = RVec {
        0., 0., 0.
    };
    RVec deltakRow     = {fieldGrid.unitCell().basisVectorLength(XX), 0, 0};
    RVec deltakColumn  = {0, fieldGrid.unitCell().basisVectorLength(YY), 0};
    RVec deltakSection = {0, 0, fieldGrid.unitCell().basisVectorLength(ZZ)};
    auto itMidSection  = field_.iteratorAtMultiIndex({{0, 0, extend[ZZ] / 2 + 1}});
    for (std::vector<t_complex>::iterator itValue = field_.iteratorAtMultiIndex({{0, 0, 0}}); itValue != itMidSection;
         itValue += extend[XX]*extend[YY])
    {
        applyToAllColumnsWithinSection_(
                k, extend[XX], extend[YY], deltakRow, deltakColumn,
                itValue, itValue + extend[XX]*extend[YY], appliedFunction);
        rvec_inc(k, deltakSection);
    }

    k = {0., 0., 0.};
    rvec_dec(k, deltakSection);

    for (auto itValue = field_.end(); itValue != itMidSection;
         itValue -= extend[XX]*extend[YY])
    {
        applyToAllColumnsWithinSection_(
                k, extend[XX], extend[YY], deltakRow, deltakColumn,
                itValue-extend[XX]*extend[YY], itValue, appliedFunction);
        rvec_dec(k, deltakSection);
    }
    return *this;
};

void ApplyToUnshiftedFourierTransform::applyToAllColumnsWithinSection_(
        RVec coordinateSectionBegin, int nRowsPerColumn, int nColumns,
        RVec deltakRow, RVec deltakColumn,
        std::vector<t_complex>::iterator itColumnStart,
        std::vector<t_complex>::iterator itColumnEnd,
        FunctionOnComplexField appliedFunction)
{

    auto itMiddleColumn = itColumnStart + ((nColumns + 1) / 2) * nRowsPerColumn;
    auto k              = coordinateSectionBegin;
    for (auto itValue = itColumnStart; itValue != itMiddleColumn;
         itValue += nRowsPerColumn)
    {
        applyToAllRowsWithinColumn_(k, deltakRow, itValue, itValue + nRowsPerColumn,
                                    appliedFunction);
        rvec_inc(k, deltakColumn);
    }

    k = coordinateSectionBegin;
    rvec_dec(k, deltakColumn);
    for (auto itValue = itColumnEnd; itValue != itMiddleColumn;
         itValue -= nRowsPerColumn)
    {
        applyToAllRowsWithinColumn_(k, deltakRow, itValue - nRowsPerColumn, itValue,
                                    appliedFunction);
        rvec_dec(k, deltakColumn);
    }
};

void ApplyToUnshiftedFourierTransform::applyToAllRowsWithinColumn_(
        RVec coordinateColumnBegin, RVec deltakRow,
        std::vector<t_complex>::iterator itRowStart,
        std::vector<t_complex>::iterator itRowEnd,
        FunctionOnComplexField appliedFunction)
{

    auto k = coordinateColumnBegin;
    for (auto itValue = itRowStart; itValue != itRowEnd; ++itValue)
    {
        appliedFunction(*itValue, k);
        rvec_inc(k, deltakRow);
    }
}
}
