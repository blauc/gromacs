/*
 * ImprovedFastGaussTransform.cpp
 *
 *  Created on: Aug 25, 2015
 *      Author: cblau
 */

#include "improvedfastgausstransform.h"
#include "vec.h"
#include "volumedata.h"
#include "gausstransform.h"
#include "gromacs/simd/simd_math.h"

#include "gromacs/utility/exceptions.h"

#include <algorithm>
namespace gmx
{

namespace volumedata
{

ExpansionCenter::ExpansionCenter(ImprovedFastGaussTransform * ifgt) : ifgt_(ifgt)
{
    coefficients_.resize(ifgt_->n_I_, 0);
};

ExpansionCenter::~ExpansionCenter(){};

/*
 * Gamma((i+1)/2)
 */
real ExpansionCenter::gamma_1_2(unsigned int i)
{
    static const real result[11] = { 1.77245385090552, 1.00000000000000, 0.886226925452758, 1.00000000000000, 1.32934038817914, 2.00000000000000, 3.32335097044784, 6.00000000000000, 11.6317283965674, 24.0000000000000, 52.3427777845535 };
    return i < 11 ? result[i] : gamma_1_2(i - 2) * (i - 1) / 2;
}

bool ExpansionCenter::is_even(int i)
{
    return (i & 1) == 0;
}

void ExpansionCenter::total_density()
{
    total_density_ = 0;

    int  i_x;
    int  i_y;

    auto I_ = coefficients_.begin();

    for (int o = 0; o <= ifgt_->maximum_expansion_order_; ++o)
    {

        for (i_x = o; i_x >= 1; --i_x)
        {
            for (i_y = o - i_x; i_y != 0; --i_y)
            {
                if (is_even(i_x) && is_even(i_y) && is_even(o - i_x - i_y))
                {
                    total_density_ += *I_ * gamma_1_2(i_x) * gamma_1_2(i_y)
                        * gamma_1_2(o - i_x - i_y);
                    ++I_;
                }
            }
        }
        for (i_y = o; i_y >= 1; --i_y)
        {

            if (is_even(i_x) && is_even(i_y) && is_even(o - i_x - i_y))
            {
                total_density_ += *I_ * gamma_1_2(i_x) * gamma_1_2(i_y)
                    * gamma_1_2(o - i_x - i_y);
                ++I_;
            }
        }

        if (is_even(i_x) && is_even(i_y) && is_even(o - i_x - i_y))
        {
            total_density_ += *I_ * gamma_1_2(i_x) * gamma_1_2(i_y)
                * gamma_1_2(o - i_x - i_y);
            ++I_;
        }
    }
}

real ExpansionCenter::expand(RVec dx)
{
    auto I_     = coefficients_.begin();
    real result = *I_;
    ++I_;

    auto mn_ = ifgt_->multinomials_.begin();

    *mn_ = 1;
    auto I_x_ = mn_;
    auto I_y_ = mn_;
    auto I_z_ = mn_;

    for (int o = 1; o <= ifgt_->maximum_expansion_order_; ++o)
    {

        auto I_pre_ = I_x_;
        I_x_   = mn_ + 1;
        for (int i_x = o; i_x >= 1; --i_x)
        {
            for (int i_y = o - i_x; i_y >= 0; --i_y)
            {

                ++mn_;
                *mn_    = *I_pre_ * dx[XX];
                result += *mn_ * *I_;
                ++I_;
                ++I_pre_;
            }
        }

        I_pre_ = I_y_;
        I_y_   = mn_ + 1;
        for (int i_y = o; i_y >= 1; --i_y)
        {
            ++mn_;
            *mn_    = *I_pre_ * dx[YY];
            result += *mn_ * *I_;
            ++I_;
            ++I_pre_;
        }

        I_pre_ = I_z_;
        I_z_   = mn_ + 1;
        ++mn_;
        *mn_    = *I_pre_ * dx[ZZ];
        result += *mn_ * *I_;
        ++I_;
    }

    return result;
}

void ExpansionCenter::compress(rvec dist, real weight)
{
    real pre = weight * exp(-norm2(dist));
    int  m[3];
    auto mn_  = ifgt_->multinomials_.begin();
    *mn_ = 1;
    auto I_   = coefficients_.begin();
    *I_ += pre;
    auto I_x_ = mn_;
    auto I_y_ = mn_;
    auto I_z_ = mn_;
    for (int o = 1; o <= ifgt_->maximum_expansion_order_; ++o)
    {

        auto I_pre_ = I_x_;
        I_x_   = mn_ + 1;
        for (m[XX] = o; m[XX] >= 1; --m[XX])
        {
            for (m[YY] = o - m[XX]; m[YY] >= 0; --m[YY])
            {

                m[ZZ] = o - m[XX] - m[YY];

                ++mn_;
                ++I_;
                *mn_ = *I_pre_ * 2 * dist[XX] / (real) m[XX];
                *I_ += pre * *mn_;
                ++I_pre_;

            }
        }

        I_pre_ = I_y_;
        I_y_   = mn_ + 1;
        for (m[YY] = o; m[YY] >= 1; --m[YY])
        {

            ++mn_;
            ++I_;
            *mn_ = *I_pre_ * 2 * dist[YY] / (real) m[YY];
            *I_ += pre * *mn_;
            ++I_pre_;

        }

        I_pre_ = I_z_;
        I_z_   = mn_ + 1;
        m[ZZ]  = o;
        ++mn_;
        ++I_;
        *mn_ = *I_pre_ * 2 * dist[ZZ] / (real) m[ZZ];
        *I_ += pre * *mn_;

    }

};

uint64_t ImprovedFastGaussTransform::factorial(const int n)
{
    static const uint64_t           result[21] = { 1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600, 6227020800, 87178291200, 1307674368000, 20922789888000, 355687428096000, 6402373705728000, 121645100408832000, 2432902008176640000 };
    return n < 21 ? result[n] : n * factorial(n - 1);
}

void ImprovedFastGaussTransform::set_grid(std::unique_ptr<GridReal> grid)
{
    FiniteGrid expansionGrid;
    maximum_expansion_order_ = 3;
    h_inv_                   = 1./(::sqrt(2) * sigma_);
    float expansionCenterSpacing = sigma_ * (::sqrt(maximum_expansion_order_ / 2.0) + 1);
    n_I_                     = factorial(maximum_expansion_order_ + DIM) / (factorial(DIM) * factorial(maximum_expansion_order_));
    multinomials_.resize(n_I_);
    expansionGrid.copy_grid(*grid);
    auto gridScaleFactor = expansionGrid.avg_spacing() / (expansionCenterSpacing);
    expansionGrid.multiplyGridPointNumber({gridScaleFactor, gridScaleFactor, gridScaleFactor});
    expansionGrid.makeGridUniform();

    expansionCenterField_.copy_grid(expansionGrid);
    grid_ = std::move(grid);
}

RVec
ImprovedFastGaussTransform::distanceToExpansionCenter(IVec expansionCenterIndex, RVec x)
{
    RVec result;
    auto closestExpansionCenterCoordinate = expansionCenterField_.gridpoint_coordinate(expansionCenterIndex);
    rvec_sub(closestExpansionCenterCoordinate, x, result);
    svmul(h_inv_, result, result);
    return result;
}

void ImprovedFastGaussTransform::transform(const rvec x, real weight)
{

    auto closestExpansionCenterIndex = expansionCenterField_.coordinate_to_gridindex_round_ivec(x);
    if (expansionCenterField_.inGrid(closestExpansionCenterIndex))
    {
        auto &expansionCenter = expansionCenterField_.access().at(closestExpansionCenterIndex);
        if (expansionCenter == nullptr)
        {
            expansionCenter = std::unique_ptr<ExpansionCenter>(new ExpansionCenter(this));
        }
        expansionCenter->compress(distanceToExpansionCenter(closestExpansionCenterIndex, x), weight);
    }
}

inline real
ImprovedFastGaussTransform::exp_lookup(real x)
{
    auto result = gaussLUT_.find(int(50*x));
    if (result != gaussLUT_.end())
    {
        return result->second;
    }
    else
    {
        auto val = gmx::exp(-gmx::square(x));
        gaussLUT_[int(50*x)] = val;
        return val;
    }
};


std::unique_ptr<GridReal> &&  ImprovedFastGaussTransform::finish_and_return_grid()
{

    minimumUsedGridIndex_ = {GMX_INT32_MAX, GMX_INT32_MAX, GMX_INT32_MAX};
    maximumUsedGridIndex_ = {0, 0, 0};

    RVec     d_x       = grid_->unit_cell_XX();
    RVec     d_y       = grid_->unit_cell_YY();
    RVec     d_z       = grid_->unit_cell_ZZ();

    /* this way only cubic grids are allowed */
    RVec     dExpansion {
        h_inv_ * expansionCenterField_.unit_cell_XX()[XX], h_inv_ * expansionCenterField_.unit_cell_YY()[YY], h_inv_ * expansionCenterField_.unit_cell_ZZ()[ZZ]
    };

    /*
     * Use iterators to step though the grid.
     * This relies on the assumption that the grid is stored with z the slowest and x the fasted changing dimension with no padding
     * (x,y,z not being linked to any coordiante system, but short-hand for first, second, third dimension)
     * Loosing generality through this approach, we save substantial time when we don't have to calculate the grid index.
     */

    auto      gridCoordinate_z = grid_->gridpoint_coordinate({0, 0, 0});

    auto      extend              = grid_->extend();
    auto      gridData            = grid_->access();
    auto      expansionCenterData = expansionCenterField_.access();

    const int nNeighbours = 1;

    std::array<std::vector<int>, 3>  neighbourExpansionCenterIndex;
    std::array<std::vector<real>, 3> neighbourDistance;
    std::array<std::vector<real>, 3> neighbourExpPreFactor;
    std::array<std::vector<bool>, 3> neighbourInGrid;

    for (int dimension = XX; dimension <= ZZ; ++dimension)
    {
        neighbourExpansionCenterIndex[dimension].resize(2*nNeighbours+1);
        neighbourDistance[dimension].resize(2*nNeighbours+1);
        neighbourExpPreFactor[dimension].resize(2*nNeighbours+1);
        neighbourInGrid[dimension].resize(2*nNeighbours+1);
    }

    for (int gridIndexZZ = 0; gridIndexZZ < extend[ZZ]; ++gridIndexZZ)
    {
        auto gridCoordinate_yz = gridCoordinate_z; // start at the beginning of a "y-row"
        for (int gridIndexYY = 0; gridIndexYY <  extend[YY]; ++gridIndexYY)
        {
            auto gridCoordinate_xyz = gridCoordinate_yz;  // start at the beginning of an "x-column"
            auto gridIterator       = gridData.zy_column_begin(gridIndexZZ, gridIndexYY);

            for (int gridIndexXX = 0; gridIndexXX < extend[XX]; ++gridIndexXX)
            {
                auto expansionCenterIndex = expansionCenterField_.coordinate_to_gridindex_round_ivec(gridCoordinate_xyz);

                auto distance = distanceToExpansionCenter(expansionCenterIndex, gridCoordinate_xyz);

                // precalculate to reduce effort from (Number of Neighbours)^(Number of Dimensions) to (Number of Dimensions) * (Number of Neighbours)
                for (int dimension = XX; dimension <= ZZ; ++dimension)
                {
                    for (int iNeighbour = -nNeighbours; iNeighbour <= nNeighbours; ++iNeighbour)
                    {
                        auto neighbourIndex = iNeighbour + nNeighbours;
                        neighbourInGrid[dimension][neighbourIndex] = expansionCenterField_.inGrid(expansionCenterIndex[dimension] + iNeighbour, dimension);
                        if (neighbourInGrid[dimension][neighbourIndex])
                        {
                            neighbourExpansionCenterIndex[dimension][neighbourIndex] = expansionCenterIndex[dimension] + iNeighbour;
                            neighbourDistance[dimension][neighbourIndex]             = distance[dimension] + iNeighbour * dExpansion[dimension];
                            neighbourExpPreFactor[dimension][neighbourIndex]         = exp_lookup(neighbourDistance[dimension][neighbourIndex]);
                        }
                    }
                }
                #pragma omp simd
                for (int iNeighbourZZ = 0; iNeighbourZZ < 2*nNeighbours+1; ++iNeighbourZZ)
                {
                    for (int iNeighbourYY = 0; iNeighbourYY < 2*nNeighbours+1; ++iNeighbourYY)
                    {
                        real prefactorZZ_YY = neighbourExpPreFactor[ZZ][iNeighbourZZ] * neighbourExpPreFactor[YY][iNeighbourYY];
                        for (int iNeighbourXX = 0; iNeighbourXX < 2*nNeighbours+1; ++iNeighbourXX)
                        {
                            if (neighbourInGrid[XX][iNeighbourXX] && neighbourInGrid[YY][iNeighbourYY] && neighbourInGrid[ZZ][iNeighbourZZ])
                            {
                                auto &expansionCenter = expansionCenterData.at({neighbourExpansionCenterIndex[XX][iNeighbourXX], neighbourExpansionCenterIndex[YY][iNeighbourYY], neighbourExpansionCenterIndex[ZZ][iNeighbourZZ]});
                                if (expansionCenter != nullptr)
                                {
                                    *gridIterator += prefactorZZ_YY * neighbourExpPreFactor[XX][iNeighbourXX] *  expansionCenter->expand({neighbourDistance[XX][iNeighbourXX], neighbourDistance[YY][iNeighbourYY], neighbourDistance[ZZ][iNeighbourZZ]});
                                }
                            }
                        }
                    }
                }
                ++gridIterator;
                rvec_inc(gridCoordinate_xyz, d_x); // next step in grid x-direction
            }
            rvec_inc(gridCoordinate_yz, d_y);      // next step in grid y-direction
        }
        rvec_inc(gridCoordinate_z, d_z);           // next step in grid z-direction
    }
    return std::move(grid_);
}

} /* volumedata */

} /* gmx */
