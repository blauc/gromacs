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
        FiniteGrid<N> &operator=(FiniteGrid<N> other)
        {
            std::swap(translation_, other.translation_);
            std::swap(cell_, other.cell_);
            std::swap(unitCell_, other.unitCell_);
            std::swap(lattice_, other.lattice_);
            return *this;
        }
        FiniteGrid(const FiniteGrid &other) : translation_ {other.translation_}, cell_ {other.cell_}, unitCell_ {other.unitCell_}, lattice_ {other.lattice_} {}

        FiniteGrid(const OrthogonalBasis<N> &cell, const ColumnMajorLattice<N> &lattice) : cell_ {cell}, unitCell_ {NdVector()}, lattice_ {lattice}
        {
            std::fill(std::begin(translation_), std::end(translation_), 0.);
            setUnitCell_();
        }

        NdVector coordinateToRealGridIndex(const NdVector &x) const
        {
            OrthogonalBasis<DIM>::NdVector x_shifted;
            std::transform(std::begin(x), std::end(x), std::begin(translation_), std::begin(x_shifted), [](real x, real t){return x-t; });
            return unitCell_.transformIntoBasis(x_shifted);
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
            auto result = unitCell_.transformFromBasis(multiIndexToNdVector_(i));
            std::transform(std::begin(translation_), std::end(translation_), std::begin(result), std::begin(result), [](real a, real b){return a+b; } );
            return result;
        }

        /*! \brief
         * Convert Lattice to its corresponding lattice in reciprocal space by
         * Cell inversion.
         */
        FiniteGrid<N> reciprocalGrid() const
        {
            return FiniteGrid(cell_.inverse().scaledCopy(multiIndexToNdVector_(lattice_.getExtend())), lattice_);
        }
        const ColumnMajorLattice<N> lattice() const
        {
            return lattice_;
        }
        const NdVector &translation() const
        {
            return translation_;
        }

        void setLatticeAndRescaleCell(const ColumnMajorLattice<N> &lattice)
        {
            lattice_ = lattice;
            cell_    = unitCell_.scaledCopy(multiIndexToNdVector_(lattice_.getExtend()));
        }

        /*! \brief
         *
         * Set the real-space coordinate of gridpoint (0,0,0).
         */
        void setTranslation(const NdVector &translate)
        {
            translation_ = translate;
        }

        void setCell(const OrthogonalBasis<N> &cell)
        {
            cell_ = cell;
            setUnitCell_();
        }

        OrthogonalBasis<N> cell() const
        {
            return cell_;
        }

        OrthogonalBasis<N> unitCell() const
        {
            return unitCell_;
        };

    private:
        /*! \brief
         * convert MultiIndex to NdVector
         */
        NdVector multiIndexToNdVector_(const MultiIndex &i) const
        {
            NdVector result;
            std::transform(std::begin(i), std::end(i), std::begin(result), [](int integerExtend){return real(integerExtend); });
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
        NdVector                     translation_;
        OrthogonalBasis<N>           cell_;
        OrthogonalBasis<N>           unitCell_;
        ColumnMajorLattice<N>        lattice_;
};

}

#endif
