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
#include "gridreal.h"
#include "gromacs/fft/parallel_3dfft.h"
#include "gromacs/fileio/volumedataio.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/real.h"
#include <algorithm>
#include <cassert>
#include <fftw3.h>

namespace gmx
{
namespace volumedata
{

IVec
FourierTransform3D::columnMajorExtendToRowMajorExtend(IVec extend)
{
    return {
               extend[ZZ], extend[YY], extend[XX]
    };
}

FourierTransformRealToComplex3D::FourierTransformRealToComplex3D(
        const Field<real> &realInputField)
    : realInputField_ {realInputField}, complexTransformedField_ {
    new Field<t_complex>()
} {
    /*
     * TODO: dirty trick : to get the reciprocal basis vectors correct,
     * set to the size of the real grid first
     */
    complexTransformedField_->copy_grid(realInputField_);
    complexTransformedField_->convertToReciprocalSpace();
    complexTransformedField_->set_extend({realInputField_.extend()[XX] / 2 +1,
                                          realInputField_.extend()[YY],
                                          realInputField_.extend()[ZZ] });
    complexTransformedField_->resetCell();
    plan_ = fftwf_plan_dft_r2c(3, columnMajorExtendToRowMajorExtend(realInputField_.extend()), (float *) realInputField_.access().data().data(), (fftwf_complex * )complexTransformedField_->access().data().data(), FFTW_ESTIMATE);
    fftwf_execute(plan_);
    fftwf_destroy_plan(plan_);
};

FourierTransformRealToComplex3D &FourierTransformRealToComplex3D::normalize()
{
    auto orthoscale = 1 / sqrt(complexTransformedField_->num_gridpoints());

    auto normalizeComplexTransform = [orthoscale](t_complex &value) {
            value.re *= orthoscale;
            value.im *= orthoscale;
        };

    std::for_each(complexTransformedField_->access().begin(),
                  complexTransformedField_->access().end(),
                  normalizeComplexTransform);
    return *this;
}

std::unique_ptr < Field < t_complex>> FourierTransformRealToComplex3D::result() {
    return std::move(complexTransformedField_);
};

FourierTransformComplexToReal3D::FourierTransformComplexToReal3D(
        const Field<t_complex> &complexInputField)
    : complexInputField_ {complexInputField}, realTransformedField_ {
    new Field<real>()
} {

    auto realSize = complexInputField_.extend();
    realSize[XX] = (realSize[XX] - 1) * 2;
    realTransformedField_->copy_grid(complexInputField_);
    realTransformedField_->set_extend(realSize);
    realTransformedField_->resetCell();
    realTransformedField_->convertToReciprocalSpace();

    std::vector<real> realData(realTransformedField_->num_gridpoints());

    auto              plan = fftwf_plan_dft_c2r(3, columnMajorExtendToRowMajorExtend(realTransformedField_->extend()),
                                                (fftwf_complex *)complexInputField_.access().data().data(),
                                                (float *)realTransformedField_->access().data().data(), FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
};

FourierTransformComplexToReal3D &FourierTransformComplexToReal3D::normalize()
{
    auto orthoscale         = 1 / sqrt(realTransformedField_->num_gridpoints());
    auto normalizeTransform = [orthoscale](real &value) {
            value *= orthoscale;
        };
    std::for_each(std::begin(realTransformedField_->access().data()),
                  std::end(realTransformedField_->access().data()),
                  normalizeTransform);
    return *this;
};

std::unique_ptr < Field < real>> FourierTransformComplexToReal3D::result() {
    return std::move(realTransformedField_);
}

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

    auto extend        = field_.extend();
    auto itSectionEnd  = field_.access().sectionBegin(extend[ZZ]);
    auto k             = RVec {
        0., 0., 0.
    };
    auto deltakSection = field_.unit_cell_ZZ();

    for (auto itValue = field_.access().sectionBegin(0); itValue != itSectionEnd;
         itValue = field_.access().next_slice(itValue))
    {
        applyToAllColumnsWithinSection_( k, extend[XX], extend[YY], field_.unit_cell_XX(), field_.unit_cell_YY(), itValue, field_.access().next_slice(itValue), appliedFunction);
        rvec_inc(k, deltakSection);
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

    // k = coordinateSectionBegin;
    // rvec_dec(k, deltakColumn);
    // for (auto itValue = itColumnEnd; itValue != itMiddleColumn;
    //      itValue -= nRowsPerColumn) {
    //   applyToAllRowsWithinColumn_(k, deltakRow, itValue - nRowsPerColumn, itValue,
    //                               appliedFunction);
    //   rvec_dec(k, deltakColumn);
    // }
};

void ApplyToUnshiftedFourierTransform::applyToAllRowsWithinColumn_(
        RVec coordinateColumnBegin, RVec deltakRow,
        std::vector<t_complex>::iterator itRowStart,
        std::vector<t_complex>::iterator itRowEnd,
        FunctionOnComplexField appliedFunction)
{

    auto itOnePastRowMiddle = itRowStart + (itRowEnd - itRowStart) / 2 + 1;
    auto k                  = coordinateColumnBegin;
    for (auto itValue = itRowStart; itValue != itOnePastRowMiddle; ++itValue)
    {
        appliedFunction(*itValue, k);
        rvec_inc(k, deltakRow);
    }

    // k = coordinateColumnBegin;
    // rvec_dec(k, deltakRow);
    // for (auto itValue = itRowEnd - 1; itValue != itOnePastRowMiddle - 1;
    //      --itValue) {
    //   appliedFunction(*itValue, k);
    //   rvec_dec(k, deltakRow);
    // }
}
}
}
