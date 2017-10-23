#ifndef GMX_MATH_CanonicalVectorBasis_H
#define GMX_MATH_CanonicalVectorBasis_H
#include "gromacs/utility/compare.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/math/vectypes.h"
#include <array>
#include <numeric>
#include <algorithm>
#include <functional>

namespace gmx
{

template <int N>
class CanonicalVectorBasis
{
    public:
        typedef std::array<real, N> NdVector;
        CanonicalVectorBasis(const NdVector &cell)
        {
            cell_ = cell;
            std::transform(std::begin(cell_), std::end(cell_), std::begin(cellInverse_), [](real c){return 1/c; });
        };

        NdVector transformIntoBasis(const NdVector &x)  const
        {
            NdVector result;
            std::transform(std::begin(cellInverse_), std::end(cellInverse_), std::begin(x), std::begin(result), std::multiplies<real>());
            return result;
        };

        NdVector transformFromBasis(const NdVector &x)  const
        {
            NdVector result;
            std::transform(std::begin(cell_), std::end(cell_), std::begin(x), std::begin(result), std::multiplies<real>());
            return result;
        };

        CanonicalVectorBasis<N> inverse() const
        {
            return CanonicalVectorBasis(cellInverse_);
        };


        CanonicalVectorBasis<N> scaledCopy(const NdVector &scale) const
        {
            NdVector scaledCell;
            std::transform(std::begin(cell_), std::end(cell_), std::begin(scale), std::begin(scaledCell), std::multiplies<real>());
            return CanonicalVectorBasis(scaledCell);
        };

        real volume() const
        {
            return std::accumulate(std::begin(cell_), std::end(cell_), 1.0, std::multiplies<real>());
        }

        /*! \brief
         * True, if length is same in x,y,z -direction.
         */
        bool allVectorsSameLength(real ftol, real abstol) const
        {
            auto firstLength           = *std::begin(cell_);
            auto firstCellVectorLength = std::begin(cell_)+1;
            return std::all_of(firstCellVectorLength, std::end(cell_), [firstLength, ftol, abstol](real length){return equal_real(length, firstLength, ftol, abstol); });
        };

        real basisVectorLength(const int &dimension) const
        {
            if (dimension >= N)
            {
                GMX_THROW(RangeError("Dimension of requested length must be smaller than Basis dimension."));
            }
            ;
            return cell_[dimension];
        }

    protected:
        NdVector cell_;
        NdVector cellInverse_;
};

} // namespace gmx
#endif
