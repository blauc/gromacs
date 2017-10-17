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

        FiniteGrid() : cell_ {NdVector()}, unit_cell_ {NdVector()}, lattice_ {} {
            NdVector   unity;
            MultiIndex integerUnity;
            std::fill(std::begin(unity), std::end(unity), 1.);
            std::fill(std::begin(integerUnity), std::end(integerUnity), 1);
            cell_      = OrthogonalBasis<N>(unity);
            unit_cell_ =  OrthogonalBasis<N>(unity);
            lattice_   = ColumnMajorLattice<N>(integerUnity);
            std::fill(std::begin(translation_), std::end(translation_), 0.);
        }
        FiniteGrid(const OrthogonalBasis<N> &cell, const ColumnMajorLattice<N> &lattice)
        {
            cell_    = cell;
            lattice_ = lattice;
            std::fill(std::begin(translation_), std::end(translation_), 0.);
            setUnitCell_();
        }
        NdVector coordinateToRealGridIndex(const NdVector &x) const
        {
            OrthogonalBasis<DIM>::NdVector x_shifted;
            std::transform(std::begin(x), std::end(x), std::begin(translation_), std::begin(x_shifted), [](real x, real t){return x-t; });
            return unit_cell_.transformIntoBasis(x_shifted);
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
            auto result = unit_cell_.transformFromBasis({{real(i[XX]), real(i[YY]), real(i[ZZ])}});
            std::transform(std::begin(translation_), std::end(translation_), std::begin(result), std::begin(result), [](real a, real b){return a+b; } );
            return result;
        }

        /*! \brief
         * Convert Lattice to its corresponding lattice in reciprocal space by
         * Cell inversion.
         */
        void convertToReciprocalSpace()
        {
            unit_cell_     = cell_.inverse();
            resetCell();
        }

        void scaleCell(const NdVector &scale)
        {
            cell_ = cell_.scale(scale);
            setUnitCell_();
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
            resetCell();
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
            return unit_cell_;
        };

    private:
        /*! \brief
         * Re-evaluates cell from unit cell and grid extend
         */
        void resetCell()
        {
            NdVector unitCellToCellScale;
            std::transform(std::begin(lattice_.getExtend()), std::end(lattice_.getExtend()), std::begin(unitCellToCellScale), [](int integerExtend){return real(integerExtend); });
            cell_ = unit_cell_.scale(unitCellToCellScale);
        }
        /*! \brief
         * Re-evaluate unit cell from cell and grid extend
         */
        void setUnitCell_()
        {
            NdVector cellToUnitCellScale;
            std::transform(std::begin(lattice_.getExtend()), std::end(lattice_.getExtend()), std::begin(cellToUnitCellScale), [](int integerExtend){return 1.0/real(integerExtend); });
            unit_cell_     = cell_.scale(cellToUnitCellScale);
        }
        NdVector                     translation_;
        OrthogonalBasis<N>           cell_;
        OrthogonalBasis<N>           unit_cell_;
        ColumnMajorLattice<N>        lattice_;
};

}

#endif
