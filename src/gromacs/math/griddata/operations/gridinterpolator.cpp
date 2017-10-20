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
#include "gridinterpolator.h"
#include "gromacs/math/quaternion.h"
#include "gromacs/math/griddata/field.h"
#include "gromacs/math/vec.h"

namespace gmx
{
GridInterpolator::GridInterpolator(std::unique_ptr < IGrid < DIM>> basis)
    : interpolatedGrid_ {std::unique_ptr < FieldReal3D>(new FieldReal3D(std::move(basis)))}
{
};

/*
 * for each target grid point:
 *      find rational number grid cell index in input grid
 *      use fractional part for weights
 */
std::unique_ptr < FieldReal3D>
GridInterpolator::interpolateLinearly(const FieldReal3D &other)
{
    const auto &grid   = interpolatedGrid_->getGrid();
    const auto &extend = grid.lattice().getExtend();
    for (int i_z = 0; i_z < extend[ZZ]; ++i_z)
    {
        for (int i_y = 0; i_y < extend[YY]; ++i_y)
        {
            for (int i_x = 0; i_x < extend[XX]; ++i_x)
            {
                auto r                 = grid.multiIndexToCoordinate({{i_x, i_y, i_z}});
                interpolatedGrid_->atMultiIndex({{i_x, i_y, i_z}}) = getLinearInterpolationAt(other, r);
            }
        }
    }
    return std::move(interpolatedGrid_);
};

std::unique_ptr < FieldReal3D>
GridInterpolator::interpolateLinearly(const FieldReal3D &other, const RVec &translation, const RVec &centerOfMass, const Quaternion &orientation)
{
    if (norm2(translation) < 1e-10 && orientation.norm() <  1e-10)
    {
        return interpolateLinearly(other);
    }

    const auto &grid            = interpolatedGrid_->getGrid();
    const auto &extend          = grid.lattice().getExtend();
    for (int i_z = 0; i_z < extend[ZZ]; ++i_z)
    {
        for (int i_y = 0; i_y < extend[YY]; ++i_y)
        {
            for (int i_x = 0; i_x < extend[XX]; ++i_x)
            {

                auto r                 = grid.multiIndexToCoordinate({{i_x, i_y, i_z}});
                auto shifted_r         = orientation.shiftedAndOriented(RVec(r[XX], r[YY], r[ZZ]), centerOfMass, translation);
                interpolatedGrid_->atMultiIndex({{i_x, i_y, i_z}}) = getLinearInterpolationAt(other, {{shifted_r[XX], shifted_r[YY], shifted_r[ZZ]}} );

            }
        }
    }
    return std::move(interpolatedGrid_);
}


real GridInterpolator::getLinearInterpolationAt(const FieldReal3D &field, const OrthogonalBasis<DIM>::NdVector &r) const
{
    auto iIndexInGrid = field.getGrid().coordinateToFloorMultiIndex(r);
    auto w            = field.getGrid().gridVectorFromGridPointToCoordinate(r, iIndexInGrid);

    std::array<std::array<std::array<real, 2>, 2>, 2> cube;
    const auto &lattice = field.getGrid().lattice();
    for (int ii_z = 0; ii_z <= 1; ++ii_z)
    {
        for (int ii_y = 0; ii_y <= 1; ++ii_y)
        {
            for (int ii_x = 0; ii_x <= 1; ++ii_x)
            {
                std::array<int, 3> cube_index = {{iIndexInGrid[XX]+ii_x, iIndexInGrid[YY]+ii_y, iIndexInGrid[ZZ]+ii_z}};
                if (lattice.inLattice(cube_index))
                {
                    cube[ii_x][ii_y][ii_z] = field.atMultiIndex(cube_index);
                }
                else
                {
                    cube[ii_x][ii_y][ii_z] = 0;
                }
            }
        }
    }

    std::array<std::array<real, 2>, 2> interpolated_x;
    for (int ii_z = 0; ii_z <= 1; ++ii_z)
    {
        for (int ii_y = 0; ii_y <= 1; ++ii_y)
        {
            interpolated_x[ii_y][ii_z] =
                (1 - w[XX]) * cube[0][ii_y][ii_z] + (w[XX])*cube[1][ii_y][ii_z];
        }
    }

    std::array<real, 2> interpolated_xy;
    for (int ii_z = 0; ii_z <= 1; ++ii_z)
    {
        interpolated_xy[ii_z] = (1 - w[YY]) * interpolated_x[0][ii_z] +
            (w[YY])*interpolated_x[1][ii_z];
    }

    return ((1 - w[ZZ]) * interpolated_xy[0] + w[ZZ] * interpolated_xy[1]) / 8.0;
}


}
