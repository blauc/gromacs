#ifndef GMX_MATH_GRIDCELL_H
#define GMX_MATH_GRIDCELL_H
#include "gromacs/math/vectypes.h"

namespace gmx
{

class GridCell
{
    public:
        /*! \brief
         *
         * Set unit-cell parameters from unit-cell length and angle [deg], aligning
         * the first unit-cell vecor along x-axis, placing second vector in the
         * xy-plane.
         */
        GridCell(RVec length, RVec angle);
        GridCell(matrix cell);
        /*! \brief
         * Convert Lattice to its corresponding lattice in reciprocal space by
         * Cell inversion.
         */
        GridCell inverse() const;

        const rvec   &operator[](int i) const;

        RVec coordinateTransformToCellSpace(RVec x)  const;
        RVec coordinateTransformFromCellSpace(RVec x)  const;


        GridCell scale(const RVec scale) const;

        RVec cell_lengths() const; //!< return cell lengths
        RVec cell_angles() const;  //!< return cell angles

        real volume() const;

        real deg_angle(const rvec a, const  rvec b) const;

        /*! \brief Check if all cell vectors are rectangular by calling cell_angles();
         */
        bool rectangular() const;

        /*! \brief Check, if spacing is same in x,y,z -direction.
         */
        bool spacing_is_same_xyz() const;
    private:
        matrix cell_;          //!< r' = cell_ . (r-translate_)
        matrix cellInverse_;


};

} // namespace gmx
#endif
