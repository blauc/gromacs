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
/********************************************************************
 * GridReal
 */
#include "gridreal.h"
namespace gmx
{
namespace volumedata
{

GridReal::GridReal(const Field<real> &baseField)
    : Field<real>::Field<real>(baseField){};

GridReal::GridReal(GridReal &other)
    : Field<real>::Field<real>(other), CrystalSymmetry::CrystalSymmetry(other)
{

};

GridReal::GridReal(const GridReal &other)
    : Field<real>::Field<real>(other), CrystalSymmetry::CrystalSymmetry(other)
{

};

void GridReal::multiply(real value)
{
    std::for_each(access().data().begin(), access().data().end(),
                  [value](real &v) { v *= value; });
}

real GridReal::normalize()
{
    real integratedDensity =  properties().sum() / this->num_gridpoints();
    multiply(1/integratedDensity);
    return integratedDensity;
}

void GridReal::add_offset(real value)
{
    std::for_each(access().data().begin(), access().data().end(),
                  [value](real &datum) { datum += value; });
}

ScalarGridDataProperties<real> GridReal::properties() const
{
    return ScalarGridDataProperties<real>(access().data());
}

RVec GridReal::center_of_mass()
{
    rvec weighted_grid_coordinate;
    RVec com = {0, 0, 0};
    for (size_t i = 0; i < num_gridpoints(); i++)
    {
        svmul(access().data()[i], gridpoint_coordinate(i),
              weighted_grid_coordinate);
        rvec_inc(com, weighted_grid_coordinate);
    }
    svmul(1. / (properties().sum()), com, com);
    return com;
}

std::string GridReal::print() const
{
    std::string result;
    result += "------- real number grid -------\n";
    result += FiniteGrid::print();
    result += "  min  :" + std::to_string(properties().min()) + "\n";
    result += "  max  :" + std::to_string(properties().max()) + "\n";
    result += "  mean :" + std::to_string(properties().mean()) + "\n";
    result += "  var  :" + std::to_string(properties().var()) + "\n";
    result += "  rms  :" + std::to_string(properties().rms()) + "\n";
    result += "\n----- end real number grid -----\n\n";
    return result;
}

void GridReal::zero()
{
    std::fill(access().data().begin(), access().data().end(), 0);
};

real GridReal::getLinearInterpolationAt(RVec r) const
{
    auto rIndexInGrid = coordinateToRealGridIndex(r);
    auto iIndexInGrid = coordinate_to_gridindex_floor_ivec(r);

    auto w_x = rIndexInGrid[XX] - (real)iIndexInGrid[XX];
    auto w_y = rIndexInGrid[YY] - (real)iIndexInGrid[YY];
    auto w_z = rIndexInGrid[ZZ] - (real)iIndexInGrid[ZZ];

    std::array<std::array<std::array<real, 2>, 2>, 2> cube;
    auto data = access();

    for (int ii_z = 0; ii_z <= 1; ++ii_z)
    {
        for (int ii_y = 0; ii_y <= 1; ++ii_y)
        {
            for (int ii_x = 0; ii_x <= 1; ++ii_x)
            {
                auto cube_index = iIndexInGrid;
                cube_index[XX] += ii_x;
                cube_index[YY] += ii_y;
                cube_index[ZZ] += ii_z;
                if (inGrid(cube_index))
                {
                    cube[ii_x][ii_y][ii_z] = data.at(cube_index);
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

    return ((1 - w_z) * interpolated_xy[0] + w_z * interpolated_xy[1]) / 8.0;
}
}
}
