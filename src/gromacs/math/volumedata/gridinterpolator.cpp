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
namespace gmx
{
namespace volumedata
{
GridInterpolator::GridInterpolator(const FiniteGrid &basis)
    : interpolatedGrid_ {std::unique_ptr<GridReal>(new GridReal)}
{
    interpolatedGrid_->copy_grid(basis);
};

/*
 * for each target grid point:
 *      find rational number grid cell index in input grid
 *      use fractional part for weights
 */
std::unique_ptr<GridReal>
GridInterpolator::interpolateLinearly(const GridReal &other)
{
    auto otherAccess            = other.access();
    auto interpolatedGridAccess = interpolatedGrid_->access();

    for (int i_z = 0; i_z < interpolatedGrid_->extend()[ZZ]; ++i_z)
    {
        for (int i_y = 0; i_y < interpolatedGrid_->extend()[YY]; ++i_y)
        {
            for (int i_x = 0; i_x < interpolatedGrid_->extend()[XX]; ++i_x)
            {

                auto r                 = interpolatedGrid_->gridpoint_coordinate({i_x, i_y, i_z});
                auto rIndexInOtherGrid = other.coordinateToRealGridIndex(r);
                auto iIndexInOtherGrid = other.coordinate_to_gridindex_floor_ivec(r);

                auto w_x = rIndexInOtherGrid[XX] - (real)iIndexInOtherGrid[XX];
                auto w_y = rIndexInOtherGrid[YY] - (real)iIndexInOtherGrid[YY];
                auto w_z = rIndexInOtherGrid[ZZ] - (real)iIndexInOtherGrid[ZZ];

                std::array<std::array<std::array<real, 2>, 2>, 2> cube;

                for (int ii_z = 0; ii_z <= 1; ++ii_z)
                {
                    for (int ii_y = 0; ii_y <= 1; ++ii_y)
                    {
                        for (int ii_x = 0; ii_x <= 1; ++ii_x)
                        {
                            auto cube_index = iIndexInOtherGrid;
                            cube_index[XX] += ii_x;
                            cube_index[YY] += ii_y;
                            cube_index[ZZ] += ii_z;
                            if (other.inGrid(cube_index))
                            {
                                cube[ii_x][ii_y][ii_z] = otherAccess.at(cube_index);
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
                            (1 - w_x) * cube[0][ii_y][ii_z] + (w_x)*cube[1][ii_y][ii_z];
                    }
                }

                std::array<real, 2> interpolated_xy;
                for (int ii_z = 0; ii_z <= 1; ++ii_z)
                {
                    interpolated_xy[ii_z] = (1 - w_y) * interpolated_x[0][ii_z] +
                        (w_y)*interpolated_x[1][ii_z];
                }

                interpolatedGridAccess.at({i_x, i_y, i_z}) =
                    ((1 - w_z) * interpolated_xy[0] + w_z * interpolated_xy[1]) / 8.0;
            }
        }
    }
    return std::move(interpolatedGrid_);
};

void GridInterpolator::makeUniform() { interpolatedGrid_->makeGridUniform(); };

}
}
