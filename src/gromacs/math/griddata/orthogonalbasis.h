#ifndef GMX_MATH_ORTHOGONALBASIS_H
#define GMX_MATH_ORTHOGONALBASIS_H
#include "gromacs/utility/compare.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/math/vectypes.h"
#include <array>
#include <numeric>
#include <algorithm>

namespace gmx
{


template <int N>
class OrthogonalBasis
{
    public:
        typedef std::array<real, N> NdVector;
        OrthogonalBasis(const NdVector &cell)
        {
            cell_ = cell;
            std::transform(std::begin(cell_), std::end(cell_), std::begin(cellInverse_), [](real c){return 1/c; });
        };

        NdVector transformIntoBasis(const NdVector &x)  const
        {
            NdVector result;
            std::transform(std::begin(cellInverse_), std::end(cellInverse_), std::begin(x), std::begin(result), [](real c, real x){return c*x; });
            return result;
        };

        NdVector transformFromBasis(const NdVector &x)  const
        {
            NdVector result;
            std::transform(std::begin(cell_), std::end(cell_), std::begin(x), std::begin(result), [](real c, real x){return c*x; });
            return result;
        };

        OrthogonalBasis<N> inverse() const
        {
            return OrthogonalBasis(cellInverse_);
        };

        bool isSameWithinTolerance(const OrthogonalBasis<N> &other, real ftol, real abstol) const
        {
            auto otherLength = std::begin(other.cell_);
            for (auto currentLength : cell_)
            {
                if (!equal_real(currentLength, *otherLength, ftol, abstol))
                {
                    return false;
                }
                ++otherLength;
            }
            return true;
        }

        OrthogonalBasis<N> scale(const NdVector &scale) const
        {
            NdVector scaledCell;
            std::transform(std::begin(cell_), std::end(cell_), std::begin(scale), std::begin(scaledCell), [](real c, real s){ return c*s; });
            return OrthogonalBasis(scaledCell);
        };

        real volume() const
        {
            return std::accumulate(std::begin(cell_), std::end(cell_), 1.0, [](real a, real b){return a*b; });
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

        real length(const int &dimension) const
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
