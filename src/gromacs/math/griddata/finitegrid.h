#ifndef GMX_MATH_FINITEGRID_H
#define GMX_MATH_FINITEGRID_H

#include "gromacs/math/vectypes.h"
#include "columnmajorlattice.h"
#include "orthogonalbasis.h"
#include <memory>

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
class FiniteGrid
{
    public:

        OrthogonalBasis<DIM>::NdVector coordinateToRealGridIndex(const OrthogonalBasis<DIM>::NdVector &x) const;
        ColumnMajorLattice<DIM>::MultiIndex coordinate_to_gridindex_floor_ivec(const OrthogonalBasis<DIM>::NdVector &x) const;
        OrthogonalBasis<DIM>::NdVector gridpoint_coordinate(const ColumnMajorLattice<DIM>::MultiIndex &i) const;

        /*! \brief
         * Convert Lattice to its corresponding lattice in reciprocal space by
         * Cell inversion.
         */
        void convertToReciprocalSpace();

        void scaleCell(const OrthogonalBasis<DIM>::NdVector &scale);
        /*! \brief The average grid spacing.
         */
        real avg_spacing() const;

        /*! \brief Writes all information about the grid of reals in human readable
         * form to a string.
         */
        std::string print() const;

        /*! \brief
         * Compare grid spanning vectors and translation to other.
         */
        bool sameGridInAbsTolerance(const FiniteGrid &other) const;

        bool evenlySpaced() const
        {
            return unit_cell_.allVectorsSameLength(relativeTolerance_, absoluteTolerance_);
        };


        const ColumnMajorLattice<DIM> getLattice() const;
        void setLatticeAndRescaleCell(const ColumnMajorLattice<DIM> &lattice);

        /*! \brief
         *
         * Set the real-space coordinate of gridpoint (0,0,0).
         */
        void set_translation(const OrthogonalBasis<DIM>::NdVector &translate);

        void setCell(const OrthogonalBasis<DIM> &cell);
        OrthogonalBasis<DIM> getCell() const;

        OrthogonalBasis<DIM> getUnitCell() const;

    private:
        /*! \brief
         * Re-evaluates cell based on unit cell and grid extend
         */
        void resetCell();
        /*! \brief
         * set unit cell; divide cell matrix by extend in respective direction
         */
        void setUnitCell_();
        OrthogonalBasis<DIM>::NdVector translation_ = {{0., 0., 0.}};
        OrthogonalBasis<DIM>           cell_ {{{1., 1., 1.}}};
        OrthogonalBasis<DIM>           unit_cell_ {{{1., 1., 1.}}};
        ColumnMajorLattice<DIM>        lattice_ {{{1, 1, 1}}};
        real relativeTolerance_ = 1e-6;
        real absoluteTolerance_ = 1;
};

}

#endif
