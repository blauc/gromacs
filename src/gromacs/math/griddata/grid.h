/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
 * Defines an N-dimensional grid.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \inpublicapi
 */

#ifndef GMX_MATH_FINITEGRID_H
#define GMX_MATH_FINITEGRID_H

#include "gromacs/math/vectypes.h"
#include "columnmajorlattice.h"
#include "canonicalvectorbasis.h"
#include <memory>
#include <cmath>
namespace gmx
{

/*! \brief
 * Interface to a finite N-dimensional grid that relates a finite lattice to
 * N-dimensional space via a basis.
 *
 * \tparam N The dimensionality of the grid
 */
template <int N>
class IGrid
{
    public:
        typedef typename CanonicalVectorBasis<N>::NdVector NdVector;
        typedef typename ColumnMajorLattice<N>::MultiIndex MultiIndex;
        typedef std::unique_ptr < IGrid < N>> IGridPointer;

        /*! \brief
         * Vector pointing from a gridpoint to coordinate in internal grid coordinates.
         *
         * \param[in] x untransformed input coordinate
         * \param[in] i index of lattice point
         * \returns vector from lattice point to coordinate in grid coordinates
         */
        virtual NdVector gridVectorFromGridPointToCoordinate(const NdVector &x, const MultiIndex &i) const = 0;
        /*! \brief
         * The integer part of a point in internal Grid coordinates.
         *
         * \param[in] x untransformed input coordinate
         *
         * \returns lattice point multi index
         */
        virtual MultiIndex coordinateToFloorMultiIndex(const NdVector &x) const     = 0;
        /*! \brief
         * The external coordinate of a lattice point.
         *
         * \param[in] i lattice point multi index
         *
         * \returns the coordinate cooresponding to a lattice point index by i.
         */
        virtual NdVector multiIndexToCoordinate(const MultiIndex &i) const          = 0;

        /*! \brief
         * Resets the lattice and extends or shrinks the containing grid cell correspondingly.
         * The unit cell stays the same, while grid size may change.
         * \param[in] lattice the new lattice
         */
        virtual void setLatticeAndRescaleCell(const ColumnMajorLattice<N> &lattice) = 0;

        /*! \brief
         * Lattice accessor.
         * \returns the lattice
         */
        virtual ColumnMajorLattice<N> lattice() const      = 0;

        /*! \brief
         * grid cell accessor.
         * \returns the grid cell comprising all lattice
         */
        virtual CanonicalVectorBasis<N> cell() const       = 0;

        /*! \brief
         * Unit cell accessor.
         * Unit cell times extend yields the grid cell.
         * \returns the unit cell of the grid.
         */
        virtual CanonicalVectorBasis<N> unitCell() const   = 0;
        /*! \brief
         * Return an owning pointer to the abstract base class of a copy of a specific implementation.
         *
         * Enables copying of classes holding unique_ptr's to this interface.
         *
         * \returns IGridPointer to a copy of this.
         */
        virtual IGridPointer duplicate() const = 0;

};

/*! \brief
 * A finite N-dimensional grid that is anchored at the coordinate system origin.
 *
 * Relation of the integer indices,i, to real space vector, r, are given by the unit cell.
 */
template <int N>
class Grid : public IGrid<N>
{
    public:
        typedef typename CanonicalVectorBasis<N>::NdVector NdVector;
        typedef typename ColumnMajorLattice<N>::MultiIndex MultiIndex;

        /*! \brief
         * Constructs Grid from CanonicalVectorBasis and ColumnMajorLattice.
         */
        Grid(const CanonicalVectorBasis<N> &cell, const ColumnMajorLattice<N> &lattice) : IGrid<N>(), cell_ {cell}, lattice_ {lattice}, unitCell_ {calculateUnitCellFromLatticeAndCell(lattice, cell)}
        {
        }

        /*! \brief
         * Copy constructor declared, because the default constructor is deleted.
         */
        Grid(const Grid &other) : cell_ {other.cell_}, unitCell_ {other.unitCell_}, lattice_ {other.lattice_} {}

        /*! \brief
         * Copy assignment operator declared, because the default constructor is deleted.
         */
        Grid<N> &operator=(Grid<N> other)
        {
            std::swap(cell_, other.cell_);
            std::swap(unitCell_, other.unitCell_);
            std::swap(lattice_, other.lattice_);
            return *this;
        }

        //! \copydoc IGrid::gridVectorFromGridPointToCoordinate(const NdVector &x, const MultiIndex & i)
        NdVector gridVectorFromGridPointToCoordinate(const NdVector &x, const MultiIndex &i) const override
        {
            auto     realValuedIndex = unitCell_.transformIntoBasis(x);
            NdVector result;
            std::transform(std::begin(realValuedIndex), std::end(realValuedIndex), std::begin(i), std::begin(result), [](real r, int ndx){return r-ndx; });
            return result;
        };

        //! \copydoc IGrid::coordinateToFloorMultiIndex(const NdVector &x)
        MultiIndex coordinateToFloorMultiIndex(const NdVector &x) const override
        {
            const auto &realValuedIndex = unitCell_.transformIntoBasis(x);
            MultiIndex  result;
            std::transform(std::begin(realValuedIndex), std::end(realValuedIndex), std::begin(result), [](real x){return static_cast<int>(floor(x)); });
            return result;
        }

        //! \copydoc IGrid::multiIndexToCoordinate(const MultiIndex &i)
        NdVector multiIndexToCoordinate(const MultiIndex &i) const override
        {
            return unitCell_.transformFromBasis(multiIndexToNdVector(i));
        }

        //! \copydoc IGrid::setLatticeAndRescaleCell(const ColumnMajorLattice<N> &lattice)
        void setLatticeAndRescaleCell(const ColumnMajorLattice<N> &lattice) override
        {
            lattice_ = lattice;
            cell_    = unitCell_.scaledCopy(multiIndexToNdVector(lattice_.getExtend()));
        }

        //! \copydoc IGrid::lattice()
        ColumnMajorLattice<N> lattice() const override
        {
            return lattice_;
        }

        //! \copydoc IGrid::cell()
        CanonicalVectorBasis<N> cell() const override
        {
            return cell_;
        }

        //! \copydoc IGrid::unitCell()
        CanonicalVectorBasis<N> unitCell() const override
        {
            return unitCell_;
        };

        //! \copydoc IGrid::duplicate()
        IGridPointer duplicate() const override
        {
            return IGridPointer(new Grid<N>(*this));
        }

    private:

        /*! \brief
         * convert MultiIndex to NdVector
         */
        NdVector multiIndexToNdVector(const MultiIndex &i) const
        {
            NdVector result;
            // employ implicit casting in std::copy
            std::copy(std::begin(i), std::end(i), std::begin(result));
            return result;
        }
        /*! \brief
         * Re-evaluate unit cell from cell and grid extend
         */
        CanonicalVectorBasis<N> calculateUnitCellFromLatticeAndCell(const ColumnMajorLattice<N> &lattice, const CanonicalVectorBasis<N> cell)
        {
            NdVector    cellToUnitCellScale;
            const auto &invertInteger = [](int integerExtend){return 1.0/real(integerExtend); };
            std::transform(std::begin(lattice.getExtend()), std::end(lattice.getExtend()), std::begin(cellToUnitCellScale), invertInteger);
            return cell.scaledCopy(cellToUnitCellScale);
        }
        CanonicalVectorBasis<N>    cell_;
        ColumnMajorLattice<N>      lattice_;
        CanonicalVectorBasis<N>    unitCell_;
};

/*! \brief
 * An N-dimensional grid that is shifted from the origin.
 *
 * \tparam N dimensionality of the grid
 */
template <int N>
class GridWithTranslation : public Grid<N>
{
    public:
        typedef typename CanonicalVectorBasis<N>::NdVector NdVector;
        typedef typename ColumnMajorLattice<N>::MultiIndex MultiIndex;

        GridWithTranslation(const CanonicalVectorBasis<N> &cell, const ColumnMajorLattice<N> &lattice, const NdVector &translation) : Grid<N>{cell, lattice}, translation_ {translation} { }

        GridWithTranslation(const GridWithTranslation &other) : Grid<N>{other}, translation_ {other.translation_} { }

        GridWithTranslation<N> &operator=(GridWithTranslation<N> other)
        {
            Grid<N>::operator=(other);
            std::swap(translation_, other.translation_);
            return *this;
        }

        NdVector gridVectorFromGridPointToCoordinate(const NdVector &x, const MultiIndex &i) const override
        {
            return Grid<N>::gridVectorFromGridPointToCoordinate(translateIntoGrid_(x), i);
        };

        MultiIndex coordinateToFloorMultiIndex(const NdVector &x) const override
        {
            return Grid<N>::coordinateToFloorMultiIndex(translateIntoGrid_(x));
        }

        NdVector multiIndexToCoordinate(const MultiIndex &i) const override
        {
            return translateFromGrid_(Grid<N>::multiIndexToCoordinate(i));
        }

        /*! \brief
         *
         * Set the real-space coordinate of gridpoint (0,0,0).
         */
        void setTranslation(const NdVector &translate)
        {
            translation_ = translate;
        }

        IGridPointer duplicate() const override
        {
            return IGridPointer(new GridWithTranslation<N>(*this));
        }

    private:
        NdVector translateIntoGrid_(const NdVector &x) const
        {
            CanonicalVectorBasis<DIM>::NdVector x_translated;
            std::transform(std::begin(x), std::end(x), std::begin(translation_), std::begin(x_translated), [](real x, real t){return x-t; });
            return x_translated;
        };
        NdVector translateFromGrid_(const NdVector &x) const
        {
            CanonicalVectorBasis<DIM>::NdVector x_translated;
            std::transform(std::begin(x), std::end(x), std::begin(translation_), std::begin(x_translated), [](real x, real t){return x+t; });
            return x_translated;
        };
        NdVector                     translation_;
};

}

#endif
