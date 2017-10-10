#include "gromacs/math/vectypes.h"

namespace gmx
{

class GridCell
{
    public:
        /*! \brief
         * Convert Lattice to its corresponding lattice in reciprocal space by
         * Cell inversion.
         */
        GridCell inverse();

        const matrix &getMatrix() const;

        const rvec   &operator[](int i) const;

        /*! \brief
         *
         * Set unit-cell parameters from unit-cell length and angle [deg], aligning
         * the first unit-cell vecor along x-axis, placing second vector in the
         * xy-plane.
         */
        void set_cell(RVec length, RVec angle);

        GridCell scale(const RVec scale) const;

        RVec cell_lengths() const; //!< return cell lengths
        RVec cell_angles() const;  //!< return cell angles
        /*! \brief
         * Yields rotation matrix Q that will rotate grid towards canonical
         * representation with first unit-cell vector aligned along x-axis, and second
         * unit-cell vector aligned in xy-plane.
         * \param[in,out] Q rotation matrix
         */
        void rotation(matrix Q) const;

        real volume() const;

        real deg_angle(const rvec a, const  rvec b) const;

        real grid_cell_volume() const;

        /*! \brief Check if all cell vectors are rectangular by calling cell_angles();
         */
        bool rectangular() const;

        /*! \brief Check, if spacing is same in x,y,z -direction.
         */
        bool spacing_is_same_xyz() const;
    private:
        matrix cell_;          //!< r' = cell_ . (r-translate_)


};

} // namespace gmx
