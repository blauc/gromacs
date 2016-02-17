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

#ifndef GMX_MATH_VOLUMEDATA_H_
#define GMX_MATH_VOLUMEDATA_H_
#include <memory>
#include <vector>

#include "gromacs/math/vectypes.h"

namespace gmx
{

namespace volumedata
{

typedef BasicVector<int> IVec; //!< RVec equivalent for int
/*!  \brief
 * Methods for handling crystal symmetry.
 */
class CrystalSymmetry
{
    public:
        CrystalSymmetry ();
        ~CrystalSymmetry ();

        /*! \brief
         * Set space group.
         *
         * \param[in] space_group according to "International Tables for Crystallography Table 12.3.4.1 Standard space-group symbol"
         */
        void set_space_group(int space_group);

        /*! \brief
         * retreive space group
         *
         * \returns space group according to "International Tables for Crystallography Table 12.3.4.1 Standard space-group symbol"
         */
        int space_group();
    private:
        class Impl;
        std::unique_ptr<CrystalSymmetry::Impl> impl_;
};


/*!
 * \brief
 * A finite three-dimensional grid and conversion routines for different grid representations.
 *
 * Represents a three dimensional grid with integer (x,y,z) indexing.
 * Relation of the integer indices,i, to real space vector, r, are given by cell and translation:
 * i = cell . (r-translate); r = cell^{-1}.i + translate
 */
class FiniteGrid
{
    public:
        FiniteGrid();
        ~FiniteGrid();

        /*! \brief
         * The extend of the grid.
         *
         * Grid indices will alwasy run from (0,0,0) to extend = (extend[XX],extend[YY],extend[ZZ])
         */
        void   set_extend(IVec extend);
        IVec   extend();         //!< return the extend of the grid
        size_t num_gridpoints(); //!< evaluates extend[0]*extend[1]*extend[2]

        /*! \brief
         * Unique one-dimensional grid index  = x + extend[XX] * y + extend[XX] * extend[YY] * z.
         */
        int  ndx3d_to_ndx1d(IVec i_xyz);

        /*! \brief
         * Inverse for ndx3d_to_ndx1d ;
         */
        IVec ndx1d_to_ndx3d(int i);

        /*! \brief
         *
         * Set unit-cell parameters from unit-cell length and angle, aligning the first unit-cell vecor along x-axis, placing second vector in the xy-plane.
         */
        void set_cell(RVec length, RVec angle);

        /*! \brief
         *
         * Set the real-space coordinate of gridpoint (0,0,0).
         */
        void set_translation(RVec translate);
        RVec translation();  //!< return real-space coordinate of gridpoint (0,0,0).
        RVec cell_lengths(); //!< return unit-cell lengths
        RVec cell_angles();  //!< return unit-cell angles

        /*! \brief
         * Yields rotation matrix Q that will rotate grid towards canonical representation with first unit-cell vector aligned along x-axis, and second unit-cell vector aligned in xy-plane.
         * \param[in,out] Q rotation matrix
         */
        void rotation(matrix Q);

        RVec gridpoint_coordinate(int i);

        RVec gridpoint_coordinate(IVec i);

    private:
        class Impl;
        std::unique_ptr<FiniteGrid::Impl> impl_;
};


/*!
 * \brief
 * Real-space, real-value data on a dense grid with symmetry (symmetry can be none).
 *
 * Data is stored in 1d vector following (x,y,z) convention: x fastest, y medium, z slowest dimension
 */
class GridReal : public FiniteGrid, public CrystalSymmetry
{
    public:
        GridReal();
        ~GridReal ();
        /*! \brief
         * Access a section of the volume density data via start and end iterator.
         */
        std::pair<std::vector<real>::iterator, std::vector<real>::iterator> z_section(int section);
        /*! \brief
         * Access row of the volume density data through start and end iterator.
         */
        std::pair<std::vector<real>::iterator, std::vector<real>::iterator> zy_row(int section, int row);
        /*! \brief
         * Directly access an index element.
         * \throws std::out_of_range if element is out of array bounds
         */
        real &at(IVec index);
        /*! \brief The minimum grid data value. */
        real min();
        /*! \brief The maximum grid data value. */
        real max();
        /*! \brief The mean grid data value. */
        real mean();
        /*! \brief The root mean square deviation of grid data values. */
        real rms();
        /*! \brief The variance of data values.
         *
         * Two pass algorithm, where mean is calculated fist, due to large expected amount of data. (Typical data size=1,000,000)
         */
        real var();

        /*! \brief
         * Adjust data vector size to given grid.
         * \throws std::bad_alloc when unsuccesful
         */
        void resize();
        /*! \brief
         * Return the raw 1d grid data
         *
         */
        std::vector<real> &data();
        /*! \brief
         * The size of the griddata in bytes.
         */
        size_t data_size();
    private:
        class Impl;
        std::unique_ptr<GridReal::Impl> impl_;
};

}      // namespace grid_data

}      // namespace gmx

#endif // GMX_MATH_VOLUMEDATA_H_
