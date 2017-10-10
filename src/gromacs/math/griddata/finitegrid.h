#ifndef GMX_MATH_FINITEGRID_H
#define GMX_MATH_FINITEGRID_H

#include "gromacs/math/vectypes.h"
#include "columnmajorlattice.h"
#include "gridcell.h"
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
        /*! \brief
         * Convert Lattice to its corresponding lattice in reciprocal space by
         * Cell inversion.
         */
        void convertToReciprocalSpace();

        /*! \brief
         * Compare grid spanning vectors and translation to other.
         */
        bool sameGridInAbsTolerance(const FiniteGrid &other, real tolerance = 1e-10) const;

        RVec gridpoint_coordinate(int i) const;
        RVec coordinateToRealGridIndex(const rvec x) const;

        std::array<int, 3> coordinate_to_gridindex_floor_ivec(const rvec x) const;
        RVec gridpoint_coordinate(std::array<int, 3> i) const;

        void makeGridUniform();

        /*! \brief The average grid spacing.
         */
        real avg_spacing() const;

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

        void setCell(RVec length, RVec angle);

        const ColumnMajorLattice<DIM> getLattice() const;

        void setLattice(const ColumnMajorLattice<DIM> &lattice);
        /*! \brief
         *
         * Set the real-space coordinate of gridpoint (0,0,0).
         */
        void set_translation(RVec translate);
        RVec translation() const; //!< return real-space coordinate of gridpoint (0,0,0).

        GridCell getCell() const;

        GridCell getUnitCell() const;

    private:
        GridCell                          cell_;
        GridCell                          unit_cell_;
        GridCell                          unit_cell_inv_;
        RVec                              translation_;
        ColumnMajorLattice<DIM>           lattice_ {{{1, 1, 1}}};
};

}

#endif
