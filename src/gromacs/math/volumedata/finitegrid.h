#ifndef GMX_MATH_FINITEGRID_H
#define GMX_MATH_FINITEGRID_H

#include "gromacs/math/vectypes.h"
#include "finite3dlatticeindices.h"
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
class FiniteGrid : public Finite3DLatticeIndices
{
    public:
        FiniteGrid();
        FiniteGrid(const FiniteGrid &other);
        FiniteGrid(FiniteGrid &other);
        FiniteGrid &operator= (const FiniteGrid &other);

        /*! \brief
         * Compare grid spanning vectors and translation to other.
         */

        bool sameGridInAbsTolerance(const FiniteGrid &other, real tolerance = 1e-10) const;

        ~FiniteGrid();

        /*! \brief
         * Convert Lattice to its corresponding lattice in reciprocal space by
         * Cell inversion.
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
        RVec cell_lengths() const; //!< return cell lengths
        RVec cell_angles() const;  //!< return cell angles

        /*! \brief
         * Yields rotation matrix Q that will rotate grid towards canonical
         * representation with first unit-cell vector aligned along x-axis, and second
         * unit-cell vector aligned in xy-plane.
         * \param[in,out] Q rotation matrix
         */
        void rotation(matrix Q);

        void multiplyGridPointNumber(const RVec factor);

        RVec gridpoint_coordinate(int i) const;
        RVec coordinateToRealGridIndex(const rvec x) const;

        std::vector<int> coordinate_to_gridindex_ceil_ivec(const rvec x);
        std::vector<int> coordinate_to_gridindex_round_ivec(const rvec x);
        std::vector<int> coordinate_to_gridindex_floor_ivec(const rvec x) const;
        RVec gridpoint_coordinate(std::vector<int> i) const;

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

}

#endif
