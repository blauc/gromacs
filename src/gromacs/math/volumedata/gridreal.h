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
 #ifndef GMX_MATH_GRIDREAL_H
 #define GMX_MATH_GRIDREAL_H

#include "field.h"
namespace gmx
{
namespace volumedata
{

/*!
 * \brief
 * Real-space, real-value data on a dense grid with symmetry (symmetry can be
 * none).
 *
 * Data is stored in 1d vector following (x,y,z) convention: x fastest, y
 * medium, z slowest dimension
 */
class GridReal : public Field<real>, public CrystalSymmetry
{
    public:
        GridReal() = default;
        GridReal(Field<real> baseField);
        GridReal(GridReal &other);
        void multiply(real value);

        ScalarGridDataProperties<real> properties() const;

        /*! \brief
         * Add an offset to all values.
         * \param[in] offset value to be added
         */
        void add_offset(real value);
        /*! \brief Rescale values so their sum * grid_cell_volume is one.
         */
        real normalize();
        /*! \brief Writes all information about the grid of reals in human readable
         * form to a string.
         */
        std::string print();
        /*! \brief Set all voxel values to zero.
         */
        void zero();
        /*! \brief Center of mass as the weighted sum of gridpoint coordinates.
         */
        RVec center_of_mass();
};

}
}

 #endif /* end of include guard: GMX_MATH_GRIDREAL_H */
