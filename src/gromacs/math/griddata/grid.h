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
 * Interface to a finite N-dimensional grid that relates a finite lattice to N-dimensional space via a basis.
 */
template <int N>
class IGrid
{
    public:
        typedef typename CanonicalVectorBasis<N>::NdVector NdVector;
        typedef typename ColumnMajorLattice<N>::MultiIndex MultiIndex;

        /*! \brief
         * Vector pointing from a gridpoint to coordinate in internal grid coordinates.
         */
        virtual NdVector gridVectorFromGridPointToCoordinate(const NdVector &x, const MultiIndex &i) const = 0;
        virtual MultiIndex coordinateToFloorMultiIndex(const NdVector &x) const     = 0;
        virtual NdVector multiIndexToCoordinate(const MultiIndex &i) const          = 0;
        virtual void setLatticeAndRescaleCell(const ColumnMajorLattice<N> &lattice) = 0;
        virtual ColumnMajorLattice<N> lattice() const      = 0;
        virtual CanonicalVectorBasis<N> cell() const       = 0;
        virtual CanonicalVectorBasis<N> unitCell() const   = 0;
        virtual std::unique_ptr < IGrid < N>> duplicate() const = 0;

};

/*! \brief
 * A finite N-dimensional grid relates a finite lattice to N-dimensional space via a basis.
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
        Grid(const CanonicalVectorBasis<N> &cell, const ColumnMajorLattice<N> &lattice) : IGrid<N>(), cell_ {cell}, unitCell_ {NdVector()}, lattice_ {lattice}
        {
            setUnitCell_();
        }

        /* \brief
         * Copy constructor declared, because the default constructor is deleted.
         */
        Grid(const Grid &other) : cell_ {other.cell_}, unitCell_ {other.unitCell_}, lattice_ {other.lattice_} {}

        /* \brief
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
            return unitCell_.transformFromBasis(multiIndexToNdVector_(i));
        }

        //! \copydoc IGrid::setLatticeAndRescaleCell(const ColumnMajorLattice<N> &lattice)
        void setLatticeAndRescaleCell(const ColumnMajorLattice<N> &lattice) override
        {
            lattice_ = lattice;
            cell_    = unitCell_.scaledCopy(multiIndexToNdVector_(lattice_.getExtend()));
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

        std::unique_ptr < IGrid < N>> duplicate() const override
        {
            return std::unique_ptr < IGrid < N>>(new Grid<N>(*this));
        }

    private:

        /*! \brief
         * convert MultiIndex to NdVector
         */
        NdVector multiIndexToNdVector_(const MultiIndex &i) const
        {
            NdVector result;
            std::copy(std::begin(i), std::end(i), std::begin(result));
            return result;
        }
        /*! \brief
         * Re-evaluate unit cell from cell and grid extend
         */
        void setUnitCell_()
        {
            NdVector cellToUnitCellScale;
            std::transform(std::begin(lattice_.getExtend()), std::end(lattice_.getExtend()), std::begin(cellToUnitCellScale), [](int integerExtend){return 1.0/real(integerExtend); });
            unitCell_     = cell_.scaledCopy(cellToUnitCellScale);
        }
        CanonicalVectorBasis<N>    cell_;
        CanonicalVectorBasis<N>    unitCell_;
        ColumnMajorLattice<N>      lattice_;
};

/*
 *
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

        std::unique_ptr < IGrid < N>> duplicate() const override
        {
            return std::unique_ptr < IGrid < N>>(new GridWithTranslation<N>(*this));
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
