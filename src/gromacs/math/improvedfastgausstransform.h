/*
 * ImprovedFastGaussTransform.h
 *
 *  Created on: Aug 25, 2015
 *      Author: cblau
 */
/*! \brief
 *
 * The improved fast Gauss transform with periodic boundary conditions after
 *
 * "Yang, C., Duraiswami, R., Gumerov, N., & Davis, L. (2003, October). "
 * "Improved fast gauss transform and efficient kernel density estimation. "
 * "In Computer Vision, 2003. Proceedings. "
 * "Ninth IEEE International Conference on (pp. 664-671). IEEE."
 *
 * Speedup for bandwidth > voxel-spacing achieved through expansion of the
 * density in multinomial times Gauss function at distributed centers.
 *
 * Multinomial coefficients are used for density refinement force calculation.
 *
 * ImprovedFastGaussTransform does not take care of the source coordinate, weights and box pointers
 * and their appropriate size.
 *
 */

#ifndef IFGT_H_
#define IFGT_H_

#include <vector>
#include <memory>
#include "gromacs/math/vec.h"
#include "gromacs/math/volumedata.h"
#include "gromacs/math/gausstransform.h"
#include <map>

namespace gmx
{

namespace volumedata
{

class ImprovedFastGaussTransform;

class ExpansionCenter
{
    public:
        friend class ImprovedFastGaussTransform;
        ExpansionCenter(ImprovedFastGaussTransform * ifgt);
        ~ExpansionCenter();
        void compress(rvec dist, real weight);
        real expand(RVec dx);
        void total_density();
    private:
        real gamma_1_2(unsigned int i);
        bool is_even(int i);
        std::vector<real>            coefficients_;
        ImprovedFastGaussTransform * ifgt_;
        real   total_density_ = 0;

};

class ImprovedFastGaussTransform : public GaussTransform
{
    public:

        friend class ExpansionCenter;
        void set_grid(std::unique_ptr<GridReal> grid);
        void transform(const rvec x, real weight);
        real exp_lookup(real x);

        std::unique_ptr<GridReal> && finish_and_return_grid();
        uint64_t factorial(const int n);
        real pbc_dist(ivec i, ivec n, matrix b, ivec ndx);
        RVec distanceToExpansionCenter(IVec expansionCenterIndex, RVec x);
        void integrate_densities_at_expansion_centers(volumedata::GridReal &map);

    private:
        volumedata::Field < std::unique_ptr < ExpansionCenter>> expansionCenterField_;
        int                                            n_I_;
        int                                            maximum_expansion_order_;
        std::vector<real>                              multinomials_;
        real                                           h_inv_;
        std::map<real, real> gaussLUT_;
};

}      // namesepace volumedata
}      // namespace gmx
#endif /* IFGT_H_ */
