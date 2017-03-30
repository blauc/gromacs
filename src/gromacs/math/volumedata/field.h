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

#ifndef GMX_MATH_VOLUMEDATA_FIELD_H
#define GMX_MATH_VOLUMEDATA_FIELD_H

#include "volumedata.h"

namespace gmx
{
namespace volumedata
{

template <typename T>
class Field : public FiniteGrid
{
    public:
        Field()  = default;
        ~Field() = default;
        Field(Field<T> &other) : FiniteGrid {other}
        {
            copy_grid(other);
            data_ = other.data_;
        }
        Field(const Field<T> &other) : FiniteGrid {other}
        {
            copy_grid(other);
            data_ = other.data_;
        };

        Field(const FiniteGrid &other)
        {
            copy_grid(other);
        };

        GridDataAccess<T> access() const
        {
            return GridDataAccess<T>(extend(), data_);
        };

        GridDataAccess<T> access()
        {
            return GridDataAccess<T>(extend(), data_);
        };
        /*! \brief
         * Copy the properties from another grid to this one.
         *  \param[in] grid Pointer to the grid from which the proterties will be
         * copied
         */
        void copy_grid(const FiniteGrid &grid)
        {
            FiniteGrid::copy_grid(grid);
            data_.resize(num_gridpoints());
        };

        void multiplyGridPointNumber(const RVec factor)
        {
            FiniteGrid::multiplyGridPointNumber(factor);
            data_.resize(num_gridpoints());
        };

        void set_extend(IVec extend)
        {
            FiniteGrid::set_extend(extend);
            data_.resize(num_gridpoints());
        };

    private:
        std::vector<T> data_;
};

template <class T, class F> void ApplyToField(Field<T> &field, F function)
{
    RVec d_x = field.unit_cell_XX();
    RVec d_y = field.unit_cell_YY();
    RVec d_z = field.unit_cell_ZZ();

    /*
     * Use iterators to step though the grid.
     * This relies on the assumption that the grid is stored with z the slowest
     * and x the fasted changing dimension with no padding
     * (x,y,z not being linked to any coordiante system, but short-hand for first,
     * second, third dimension)
     * Loosing generality through this approach, we save substantial time when we
     * don't have to calculate the grid index.
     */

    auto gridCoordinate_z = field.gridpoint_coordinate({0, 0, 0});

    auto extend   = field.extend();
    auto gridData = field.access();

    for (int gridIndexZZ = 0; gridIndexZZ < extend[ZZ]; ++gridIndexZZ)
    {
        auto gridCoordinate_yz =
            gridCoordinate_z; // start at the beginning of a "y-row"
        for (int gridIndexYY = 0; gridIndexYY < extend[YY]; ++gridIndexYY)
        {
            auto gridCoordinate =
                gridCoordinate_yz; // start at the beginning of an "x-column"
            auto gridIterator = gridData.zy_column_begin(gridIndexZZ, gridIndexYY);

            for (int gridIndexXX = 0; gridIndexXX < extend[XX]; ++gridIndexXX)
            {
                function(*gridIterator, gridCoordinate);
                ++gridIterator;
                rvec_inc(gridCoordinate, d_x); // next step in grid x-direction
            }
            rvec_inc(gridCoordinate_yz, d_y);  // next step in grid y-direction
        }
        rvec_inc(gridCoordinate_z, d_z);       // next step in grid z-direction
    }
};

}      // namespace volumedata

}      // namespace gmx

#endif /* end of include guard: GMX_MATH_VOLUMEDATA_FIELD_H */
