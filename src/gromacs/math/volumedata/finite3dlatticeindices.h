#ifndef GMX_MATH_FINITE3DLATTICEINDICES
#define GMX_MATH_FINITE3DLATTICEINDICES

#include "gromacs/utility/exceptions.h"

#include <array>
#include <algorithm>
#include <numeric>

namespace gmx
{

template <int N> class Finite3DLatticeIndices
{
    public:
        Finite3DLatticeIndices() = default;

        /*! \brief
         * Set the extend of the lattice.
         *
         * Indices span (0,...,0) to (extend[0]-1,...,extend[N]-1)
         */
        Finite3DLatticeIndices(std::array<int, N> extend)
        {
            if (!(std::all_of(std::begin(extend), std::end(extend), [](int i){return i > 0; })))
            {
                GMX_THROW(RangeError("Lattice extend must be larger zero in all dimensions."));
            }
            extend_ = extend;
        };

        /*! \brief
         * Get handle to lattice extend.
         */
        const std::array<int, N> &getExtend() const
        {
            return extend_;
        };

        int getNumLatticePoints() const
        {
            return std::accumulate(std::begin(extend_), std::end(extend_), 1, [](int lhs, int rhs){ return lhs * rhs; });
        };


        /*! \brief
         * Unique one-dimensional lattice index.
         *
         * x + extend_x * y + extend_x * extend_y * z.
         */
        int getLinearIndexFromLatticeIndex(const std::array<int, N> &latticeIndex) const
        {
            if (!inLattice(latticeIndex))
            {
                GMX_THROW(RangeError("Lattice index must be in lattice - greater zero and smaller than extend along all dimensions."));
            }

            int         result                 = 0;
            auto        currentDimensionExtend = extend_.begin();
            auto        pre_factor             = 1;
            for (auto ndx : latticeIndex)
            {
                result     += pre_factor * ndx;
                pre_factor *= *currentDimensionExtend;
                ++currentDimensionExtend;
            }

            return result;
        }

        /*! \brief
         * Inverse of getLinearIndexFromLatticeIndex
         */
        std::array<int, N> getLatticeIndexFromLinearIndex(int linearIndex) const
        {
            std::array<int, N> result {{0, 0, 0}};
            auto               currentLatticeIndex = result.rbegin();
            int                stride              = getNumLatticePoints() / extend_.back();
            for (auto currentDimensionExtend = extend_.begin(); currentDimensionExtend != extend_.end(); ++currentDimensionExtend)
            {
                *currentLatticeIndex  = linearIndex / stride;
                linearIndex          -= *currentLatticeIndex * stride;
                stride               /= *currentDimensionExtend;
                ++currentLatticeIndex;
            }
            return result;
        }

        /*! \brief
         * True if latticeIndex is in Lattice
         */
        bool inLattice(const std::array<int, N> &latticeIndex) const
        {
            auto extendInDimension = extend_.begin();
            for (const auto &indexInDimension : latticeIndex)
            {
                if (indexInDimension >= *extendInDimension)
                {
                    return false;
                }
                ++extendInDimension;
            }
            return true;
        }

    private:
        std::array<int, N>  extend_;
};

} //gmx

#endif
