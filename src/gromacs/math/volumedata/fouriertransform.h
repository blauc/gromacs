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
#include <memory>

namespace gmx
{
namespace volumedata
{

class FourierTransformBackAndForth
{
    public:
        FourierTransformBackAndForth(const Field<real> &input);
        std::unique_ptr < Field < real>> transform();

    private:
        const Field<real> &input_;
};


class FourierTransformRealToComplex3D
{
    public:
        FourierTransformRealToComplex3D(const Field<real> &input);
        std::unique_ptr < Field < t_complex>> transform();

    private:
        const Field<real> &input_;
};

class FourierTransformComplexToReal3D
{
    public:
        FourierTransformComplexToReal3D(const Field<t_complex> &input);
        IVec firstNonHermitian(real tolerance = 1e-5);
        bool isHermitian(real tolerance = 1e-5);
        std::unique_ptr < Field < real>> transform();

    private:
        const Field<t_complex> &input_;
};

class ApplyToUnshiftedFourierTransform
{
    public:
        ApplyToUnshiftedFourierTransform(Field<t_complex> &field) : field_ {field}
        {};
        typedef const std::function<void(t_complex &, RVec)> &FunctionOnComplexField;
        const ApplyToUnshiftedFourierTransform &
        apply(FunctionOnComplexField appliedFunction)
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

            applyToAllSections_(field_.unit_cell_XX(), field_.unit_cell_YY(),
                                field_.unit_cell_ZZ(), field_.extend(), field_.access(),
                                appliedFunction);

            return *this;
        };

    private:
        void applyToAllSections_(RVec deltakRow, RVec deltakColumn,
                                 RVec deltakSection, IVec extend,
                                 const GridDataAccess<t_complex> &gridData,
                                 FunctionOnComplexField appliedFunction)
        {
            RVec zeroCoordinate {0., 0., 0.};
            auto itSectionMiddle = gridData.sectionBegin((extend[ZZ] + 1) / 2);
            auto k               = zeroCoordinate;

            for (auto itValue = gridData.sectionBegin(0); itValue != itSectionMiddle;
                 itValue = gridData.next_slice(itValue))
            {
                applyToAllColumnsWithinSection_(
                        k, extend[XX], extend[YY], deltakRow, deltakColumn, itValue,
                        gridData.next_slice(itValue), appliedFunction);
                rvec_inc(k, deltakSection);
            }

            k = zeroCoordinate;
            rvec_dec(k, deltakSection);
            for (auto itValue = gridData.sectionBegin(extend[ZZ]);
                 itValue != itSectionMiddle;
                 itValue = gridData.previousSection(itValue))
            {
                applyToAllColumnsWithinSection_(
                        k, extend[XX], extend[YY], deltakRow, deltakColumn,
                        gridData.previousSection(itValue), itValue, appliedFunction);
                rvec_dec(k, deltakSection);
            }
        };

        void applyToAllColumnsWithinSection_(
            RVec coordinateSectionBegin, int nRowsPerColumn, int nColumns,
            RVec deltakRow, RVec deltakColumn,
            std::vector<t_complex>::iterator itColumnStart,
            std::vector<t_complex>::iterator itColumnEnd,
            FunctionOnComplexField appliedFunction)
        {

            auto itColumnMiddle = itColumnStart + ((nColumns + 1) / 2) * nRowsPerColumn;
            auto k              = coordinateSectionBegin;
            for (auto itValue = itColumnStart; itValue != itColumnMiddle;
                 itValue += nRowsPerColumn)
            {
                applyToAllRowsWithinColumn_(k, deltakRow, itValue,
                                            itValue + nRowsPerColumn, appliedFunction);
                rvec_inc(k, deltakColumn);
            }

            k = coordinateSectionBegin;
            rvec_dec(k, deltakColumn);
            for (auto itValue = itColumnEnd; itValue != itColumnMiddle;
                 itValue -= nRowsPerColumn)
            {
                applyToAllRowsWithinColumn_(k, deltakRow, itValue - nRowsPerColumn,
                                            itValue, appliedFunction);
                rvec_dec(k, deltakColumn);
            }
        };

        void applyToAllRowsWithinColumn_(RVec coordinateColumnBegin, RVec deltakRow,
                                         std::vector<t_complex>::iterator itRowStart,
                                         std::vector<t_complex>::iterator itRowEnd,
                                         FunctionOnComplexField appliedFunction)
        {
            auto itRowMiddle = itRowStart + (itRowEnd - itRowStart) / 2 + 1;
            auto k           = coordinateColumnBegin;
            for (auto itValue = itRowStart; itValue != itRowMiddle; ++itValue)
            {
                appliedFunction(*itValue, k);
                rvec_inc(k, deltakRow);
            }

            k = coordinateColumnBegin;
            rvec_dec(k, deltakRow);
            for (auto itValue = itRowEnd - 1; itValue != itRowMiddle - 1; --itValue)
            {
                appliedFunction(*itValue, k);
                rvec_dec(k, deltakRow);
            }
        }
        const Field<t_complex> &field_;
};
}
}
#endif /* end of include guard: GMX_MATH_FOURIERTRANSFORM_H */
