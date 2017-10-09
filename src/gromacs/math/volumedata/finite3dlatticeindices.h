#ifndef GMX_MATH_FINITE3DLATTICEINDICES
#define GMX_MATH_FINITE3DLATTICEINDICES

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/exceptions.h"

#include <vector>
#include <string>

namespace gmx
{

class Finite3DLatticeIndices
{
    public:
        Finite3DLatticeIndices() = default;
        Finite3DLatticeIndices(const std::vector<int> &extend){    setExtend(extend); };

        int getNumLatticePoints() const
        {
            return extend_[XX] *extend_[YY] * extend_[ZZ];
        };                                                            //!< returns extend[0]*extend[1]*extend[2]
        const std::vector<int> &getExtend() const{ return extend_; }; //!< return the extend of the lattice

        /*! \brief
         * The extend of the lattice.
         *
         * Lattice indices will alwasy run from (0,0,0) to
         * (extend[XX]-1,extend[YY]-1,extend[ZZ]-1)
         */
        void setExtend(const std::vector<int> &extend)
        {
            if (!allNonNegative(extend))
            {
                GMX_THROW(RangeError("Lattice extend must be larger or equal zero in all dimensions but is "+std::to_string(extend[XX]) + " , " + std::to_string(extend[YY])+ " , " + std::to_string(extend[ZZ]) + " ."));
            }
            extend_          = extend;
        };

        /*! \brief
         * Unique one-dimensional lattice index.
         *
         * x + extend_x * y + extend_x * extend_y * z.
         */
        std::size_t getLinearIndexFromLatticeIndex(const std::vector<int> &latticeIndex) const
        {
            auto result = latticeIndex[XX] > -1 ? latticeIndex[XX] : extend_[XX] + latticeIndex[XX];
            result += latticeIndex[YY] > -1 ? extend_[XX] * latticeIndex[YY]
                : extend_[XX] * (extend_[YY] + latticeIndex[YY]);
            result += latticeIndex[ZZ] > -1 ? extend_[XX] *extend_[YY] * latticeIndex[ZZ]
                : extend_[XX] *extend_[YY] * (extend_[ZZ] + latticeIndex[ZZ]);
            return result;
        }

        /*! \brief
         * Inverse of getLinearIndexFromLatticeIndex
         */
        std::vector<int> getLatticeIndexFromLinearIndex(int linearIndex) const
        {
            std::vector<int> result;
            result[XX] = (linearIndex % extend_[XX]) % extend_[YY];
            result[YY] = (linearIndex / extend_[XX]) % extend_[YY];
            result[ZZ] = (linearIndex / extend_[XX]) / extend_[YY];
            return result;
        }

        /*! \brief
         * True if latticeIndex is in Lattice
         */
        bool inLattice(const std::vector<int> &latticeIndex) const
        {

            if (!allNonNegative(latticeIndex))
            {
                return false;
            }
            if ((latticeIndex[XX] >= extend_[XX]) || (latticeIndex[YY] >= extend_[YY]) ||
                (latticeIndex[ZZ] >= extend_[ZZ]))
            {
                return false;
            }
            return true;
        }
        /*! \brief
         * Is true if all
         */
        bool allNonNegative(const std::vector<int> &extend) const
        {
            return (extend[XX] >= 0) && (extend[YY] >= 0) && (extend[ZZ] >= 0);
        }
    private:
        std::vector<int>  extend_;
};

} //gmx

#endif
