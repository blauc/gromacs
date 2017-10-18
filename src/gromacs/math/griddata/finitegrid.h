#ifndef GMX_MATH_FINITEGRID_H
#define GMX_MATH_FINITEGRID_H

#include "gromacs/math/vectypes.h"
#include "columnmajorlattice.h"
#include "orthogonalbasis.h"
#include <memory>
#include <cmath>
namespace gmx
{

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
template <int N>
class FiniteGrid
{
    public:
        typedef typename OrthogonalBasis<N>::NdVector NdVector;
        typedef typename ColumnMajorLattice<N>::MultiIndex MultiIndex;

        FiniteGrid(const OrthogonalBasis<N> &cell, const ColumnMajorLattice<N> &lattice) : cell_ {cell}, unitCell_ {NdVector()}, lattice_ {lattice}
        {
            setUnitCell_();
        }

        FiniteGrid(const FiniteGrid &other) : cell_ {other.cell_}, unitCell_ {other.unitCell_}, lattice_ {other.lattice_} {}

        FiniteGrid<N> &operator=(FiniteGrid<N> other)
        {
            std::swap(cell_, other.cell_);
            std::swap(unitCell_, other.unitCell_);
            std::swap(lattice_, other.lattice_);
            return *this;
        }

        NdVector coordinateToRealGridIndex(const NdVector &x) const
        {
            return unitCell_.transformIntoBasis(x);
        };

        MultiIndex coordinateToFloorMultiIndex(const NdVector &x) const
        {
            auto       realValuedIndex = coordinateToRealGridIndex(x);
            MultiIndex result;
            std::transform(std::begin(realValuedIndex), std::end(realValuedIndex), std::begin(result), [](real x){return static_cast<int>(floor(x)); });
            return result;
        }

        NdVector multiIndexToCoordinate(const MultiIndex &i) const
        {
            return unitCell_.transformFromBasis(multiIndexToNdVector_(i));
        }

        void setLatticeAndRescaleCell(const ColumnMajorLattice<N> &lattice)
        {
            lattice_ = lattice;
            cell_    = unitCell_.scaledCopy(multiIndexToNdVector_(lattice_.getExtend()));
        }

        const ColumnMajorLattice<N> lattice() const
        {
            return lattice_;
        }

        OrthogonalBasis<N> cell() const
        {
            return cell_;
        }

        OrthogonalBasis<N> unitCell() const
        {
            return unitCell_;
        };

    protected:
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
    private:
        OrthogonalBasis<N>    cell_;
        OrthogonalBasis<N>    unitCell_;
        ColumnMajorLattice<N> lattice_;
};

template <int N>
class FiniteGridWithTranslation : public FiniteGrid<N>
{
    public:
        typedef typename OrthogonalBasis<N>::NdVector NdVector;
        typedef typename ColumnMajorLattice<N>::MultiIndex MultiIndex;
        const NdVector &translation() const
        {
            return translation_;
        }

        FiniteGridWithTranslation(const OrthogonalBasis<N> &cell, const ColumnMajorLattice<N> &lattice, const NdVector &translation) : FiniteGrid<N>{cell, lattice}, translation_ {translation}
        { }

        FiniteGridWithTranslation(const FiniteGridWithTranslation &other) : FiniteGrid<N>{other}, translation_ {other.translation_}
        { }

        FiniteGridWithTranslation<N> &operator=(FiniteGridWithTranslation<N> other)
        {
            FiniteGrid<N>::operator=(other);
            std::swap(translation_, other.translation_);
            return *this;
        }


        NdVector coordinateToRealGridIndex(const NdVector &x) const
        {
            OrthogonalBasis<DIM>::NdVector x_shifted;
            std::transform(std::begin(x), std::end(x), std::begin(translation_), std::begin(x_shifted), [](real x, real t){return x-t; });
            return FiniteGrid<N>::coordinateToRealGridIndex(x_shifted);
        };


        NdVector multiIndexToCoordinate(const MultiIndex &i) const
        {
            auto result = FiniteGrid<N>::multiIndexToCoordinate(i);
            std::transform(std::begin(translation_), std::end(translation_), std::begin(result), std::begin(result), [](real a, real b){return a+b; } );
            return result;
        }

        /*! \brief
         *
         * Set the real-space coordinate of gridpoint (0,0,0).
         */
        void setTranslation(const NdVector &translate)
        {
            translation_ = translate;
        }

    private:
        NdVector                     translation_;
};

}

#endif
