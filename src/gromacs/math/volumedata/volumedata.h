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
#include "gromacs/math/vec.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <memory>
#include <numeric>
#include <vector>

#include "gromacs/utility/exceptions.h"

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
        CrystalSymmetry();
        CrystalSymmetry(const CrystalSymmetry &other);
        ~CrystalSymmetry();

        /*! \brief
         * Set space group.
         *
         * \param[in] space_group according to "International Tables for
         * Crystallography Table 12.3.4.1 Standard space-group symbol"
         */
        void set_space_group(int space_group);

        /*! \brief
         * retreive space group
         *
         * \returns space group according to "International Tables for Crystallography
         * Table 12.3.4.1 Standard space-group symbol"
         */
        int space_group() const;
        /*! \brief Writes all information about the grid of reals in human readable
         * form to a string.
         */
        std::string print() const;

    private:
        class Impl;
        std::unique_ptr<CrystalSymmetry::Impl> impl_;
};

class Finite3DLatticeIndices
{
    public:
        Finite3DLatticeIndices();
        Finite3DLatticeIndices(IVec extend);
        Finite3DLatticeIndices(const Finite3DLatticeIndices &other);
        size_t num_gridpoints()
        const;                          //!< returns pre-evaluated extend[0]*extend[1]*extend[2]
        size_t numGridPointsXY() const; //!< returns pre-evaluated extend[0]*extend[1]
        IVec extend() const;            //!< return the extend of the grid
                                        /*! \brief
                                         * The extend of the grid.
                                         *
                                         * Grid indices will alwasy run from (0,0,0) to extend =
                                         * (extend[XX]-1,extend[YY]-1,extend[ZZ]-1)
                                         */
        void set_extend(IVec extend);

        /*! \brief
         * multiply grid extend by factor */
        void multiplyExtend(const RVec factor);

        /*! \brief
         * Unique one-dimensional grid index  = x + extend[XX] * y + extend[XX] *
         * extend[YY] * z.
         */
        int ndx3d_to_ndx1d(IVec i_xyz) const;

        /*! \brief
         * Inverse for ndx3d_to_ndx1d ;
         */
        IVec ndx1d_to_ndx3d(int i) const;

        bool inGrid(IVec gridIndex) const;
        bool inGrid(int gridIndex, int dimension) const;

    private:
        IVec   extend_;
        size_t numGridPointsXY_;
        size_t numGridPoints_;
};

/*!
 * \brief
 * A finite three-dimensional grid and conversion routines for different grid
 * representations.
 *
 * Represents a three dimensional grid with integer (x,y,z) indexing.
 * Relation of the integer indices,i, to real space vector, r, are given by cell
 * and translation:
 * i = cell^{-1} . (r - translate); r = cell . i + translate
 */
class FiniteGrid : public Finite3DLatticeIndices
{
    public:
        FiniteGrid();
        FiniteGrid(const FiniteGrid &other);
        FiniteGrid &operator= (const FiniteGrid &other);

        ~FiniteGrid();

        /*! \brief convert Lattice to its corresponding lattice in reciprocal space by
         * Cell inversion.
         *
         */
        void convertToReciprocalSpace();

        /*! \brief
         * Copy the properties from another grid to this one.
         *  \param[in] grid Pointer to the grid from which the proterties will be
         * copied
         */
        void copy_grid(const FiniteGrid &grid);

        /*! \brief
         *
         * Set unit-cell parameters from unit-cell length and angle [deg], aligning
         * the first unit-cell vecor along x-axis, placing second vector in the
         * xy-plane.
         */
        void set_cell(RVec length, RVec angle);

        /*! \brief
         *
         * Set the real-space coordinate of gridpoint (0,0,0).
         */
        void set_translation(RVec translate);

        /*! \brief
         * Add a vector to the translation of the map.
         */
        RVec
        translation() const;       //!< return real-space coordinate of gridpoint (0,0,0).
        RVec cell_lengths() const; //!< return unit-cell lengths
        RVec cell_angles() const;  //!< return unit-cell angles

        /*! \brief
         * Yields rotation matrix Q that will rotate grid towards canonical
         * representation with first unit-cell vector aligned along x-axis, and second
         * unit-cell vector aligned in xy-plane.
         * \param[in,out] Q rotation matrix
         */
        void rotation(matrix Q);

        void multiplyGridPointNumber(const RVec factor);

        RVec gridpoint_coordinate(int i);
        RVec coordinateToRealGridIndex(const rvec x) const;

        IVec coordinate_to_gridindex_ceil_ivec(const rvec x);
        IVec coordinate_to_gridindex_round_ivec(const rvec x);
        IVec coordinate_to_gridindex_floor_ivec(const rvec x) const;
        RVec gridpoint_coordinate(IVec i) const;

        RVec unit_cell_XX() const;
        RVec unit_cell_YY() const;
        RVec unit_cell_ZZ() const;

        real grid_cell_volume() const;

        /*! \brief Check if all cell vectors are rectangular by calling cell_angles();
         */
        bool rectangular();

        /*! \brief Check, if spacing is same in x,y,z -direction.
         */
        bool spacing_is_same_xyz();

        void makeGridUniform();

        /*! \brief The average grid spacing.
         */
        real avg_spacing();

        /*! \brief Writes all information about the grid of reals in human readable
         * form to a string.
         */
        std::string print() const;

        /*! \brief
         * set unit cell; divide cell matrix by extend in respective direction
         */
        void set_unit_cell();

        void scaleCell(RVec scale);
        /*! \brief
         * Re-evaluates cell based on unit cell and grid extend
         */
        void resetCell();

    private:
        class Impl;
        std::unique_ptr<FiniteGrid::Impl> impl_;
};

template <class T> class BasicGridDataProperties
{
    public:
        explicit BasicGridDataProperties(const std::vector<T> &data) : data_ {data}
        {}
        /*! \brief The sum of all grid values. */
        T sum() const
        {
            return std::accumulate(std::begin(data_), std::end(data_), 0.);
        }

        /*! \brief The minimum grid data value. */
        T min() const
        {
            return *std::min_element(std::begin(data_), std::end(data_));
        };
        /*! \brief The maximum grid data value. */
        T max() const
        {
            return *std::max_element(std::begin(data_), std::end(data_));
        };
        /*! \brief The mean grid data value. */
        T mean() const { return sum() / static_cast<float>(data_.size()); };
        /*! \brief
         * The size of the griddata in bytes.
         */
        size_t data_size() const { return data_.size(); };

        T normSquared() const
        {
            T data_mean        = mean();
            return std::accumulate(
                    std::begin(data_), std::end(data_), 0.,
                    [ data_mean ](const T &a, T b) { return a + gmx::square(b - data_mean); });
        }
        /*! \brief The variance of data values.
         *
         * Two pass algorithm, where mean is calculated fist, due to large expected
         * amount of data. (Typical data size=1,000,000)
         */
        T var() const
        {
            return normSquared() / static_cast<real>(data_.size());
        }

    private:
        const std::vector<T> &data_;
};

template <class T>
class ScalarGridDataProperties : public BasicGridDataProperties<T>
{
    public:
        explicit ScalarGridDataProperties(const std::vector<T> &data)
            : BasicGridDataProperties<T>{data}, data_ {data} {};
        /*! \brief The root mean square deviation of grid data values. */
        T rms() const { return sqrt(BasicGridDataProperties<T>::var()); };

    private:
        const std::vector<T> &data_;
};

template <class T> class GridDataAccess
{
    public:
        typedef typename std::vector<T>::iterator t_iterator;
        GridDataAccess(const IVec extend, const std::vector<T> &data)
            : data_ {const_cast<std::vector<T> &>(data)},
        latticeIndices3d_ {extend} {}; // TODO: this is very dirty. Think of
                                       // proper const implementation
        GridDataAccess(IVec extend, std::vector<T> &data)
            : data_ {data}, latticeIndices3d_ {extend} {};
        /*! \brief
         * Return the raw 1d grid data
         *
         */
        std::vector<T> &data() { return data_; };
        const std::vector<T> &data() const { return data_; };
        t_iterator begin(){ return data_.begin(); };
        t_iterator end(){return data_.end(); };
        /*! \brief
         * Access a section of the volume density data via start and end iterator.
         */
        t_iterator sectionBegin(int z) const { return zy_column_begin(z, 0); }
        /*! \brief
         * Access column of the volume density data through start and end iterator.
         */
        std::array<t_iterator, 2> zy_column(int z, int y) const
        {
            return {{zy_column_begin(z, y), zy_column_begin(z, y + 1)}};
        }
        /*! \brief
         * Access the start of a column through iterator.
         * No bounds checking for fast data access
         */
        t_iterator zy_column_begin(int z, int y) const
        {
            return std::begin(data_) + latticeIndices3d_.numGridPointsXY() * z +
                   latticeIndices3d_.extend()[XX] * y;
        };

        t_iterator next_column(t_iterator x) const
        {
            return x + latticeIndices3d_.extend()[XX];
        }
        t_iterator next_slice(t_iterator x) const
        {
            return x + latticeIndices3d_.numGridPointsXY();
        };
        t_iterator previousSection(t_iterator x) const
        {
            return x - latticeIndices3d_.numGridPointsXY();
        };

        /*! \brief
         * Directly access an index element.
         * \throws std::out_of_range if element is out of array bounds
         */
        T &at(IVec index)
        {
            return data_.at(latticeIndices3d_.ndx3d_to_ndx1d(index));
        };

    private:
        std::vector<T>        &data_;
        Finite3DLatticeIndices latticeIndices3d_;
};

}      // namespace volumedata

}      // namespace gmx

#endif // GMX_MATH_VOLUMEDATA_H_
