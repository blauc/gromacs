#ifndef GMX_MATH_FINITE3DLATTICEINDICES
#define GMX_MATH_FINITE3DLATTICEINDICES

#include "gromacs/math/vectypes.h"

#include <vector>

namespace gmx
{

class Finite3DLatticeIndices
{
    public:
        Finite3DLatticeIndices() = default;
        Finite3DLatticeIndices(const std::vector<int> &extend);
        int getNumLatticePoints() const;           //!< returns extend[0]*extend[1]*extend[2]
        const std::vector<int> &getExtend() const; //!< return the extend of the lattice

        /*! \brief
         * The extend of the lattice.
         *
         * Lattice indices will alwasy run from (0,0,0) to
         * (extend[XX]-1,extend[YY]-1,extend[ZZ]-1)
         */
        void setExtend(const std::vector<int> &extend);

        /*! \brief
         * Unique one-dimensional lattice index.
         *
         * x + extend_x * y + extend_x * extend_y * z.
         */
        std::size_t getLinearIndexFromLatticeIndex(const std::vector<int> &latticeIndex) const;

        /*! \brief
         * Inverse of getLinearIndexFromLatticeIndex
         */
        std::vector<int> getLatticeIndexFromLinearIndex(int linearIndex) const;

        /*! \brief
         * True if latticeIndex is in Lattice
         */
        bool inLattice(const std::vector<int> &latticeIndex) const;

        /*! \brief
         * True if latticeIndex is in index along given dimension
         */
        bool inLattice(int latticeIndex, int dimension) const;
        /*! \brief
         * Is true if all
         */
        bool allNonNegative(const std::vector<int> &extend) const;
    private:
        std::vector<int>  extend_;
};

} //gmx

#endif
