/*
 * Ifgt.h
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
 * Ifgt does not take care of the source coordinate, weights and box pointers
 * and their appropriate size.
 *
 */

#ifndef IFGT_H_
#define IFGT_H_

#include <vector>
#include <memory>
#include "gromacs/math/vec.h"

namespace gmx
{

class MpiHelper;
class Ifgt;

class t_nb
{
    public:
        real d;
        ivec ndx;
};

namespace volumedata
{
class GaussTransform;
}

class ExpansionCenter
{
    public:
        friend class Ifgt;
        ExpansionCenter(Ifgt * ifgt);
        ~ExpansionCenter();
        void compress(rvec dist, real weight);
        real expand(rvec dx);
        void total_density();
        RVec do_force_lr_low(rvec d);
        void do_force_lr_fpoly(ivec m, rvec f, real I);
    private:
        real inner(real * a, real * b, int o);
        real gamma_1_2(unsigned int i);
        bool is_even(int i);
        std::vector<real>            coefficients_;
        std::vector<real>::iterator  I_;
        std::vector<real>::iterator  I_pre_;
        std::vector<real>::iterator  I_x_;
        std::vector<real>::iterator  I_y_;
        std::vector<real>::iterator  I_z_;
        std::vector<real>::iterator  mn_;
        Ifgt * ifgt_;
        real   voxel_density_ = 0;
        real   total_density_ = 0;

};

namespace volumedata
{
class GridReal;
}

class Ifgt
{
    public:

        friend class ExpansionCenter;

        Ifgt();
        ~Ifgt();
        void init(real sigma, int maximum_expansion_order);
        void set_box(const matrix box);
        real sigma();

        void set_mpi_helper(std::shared_ptr<MpiHelper> mpi);

        void compress(const rvec x, real weight);
        void compress_density(volumedata::GridReal &map);
        void expand(volumedata::GridReal &map);
        void matrix_inverse(const matrix m, matrix r);

        void coordinate_to_expansioncenter(const rvec r, ivec closest_e, rvec d);
        void grid_index_ivec(const rvec x, ivec i);
        void vmmul(const matrix m, const  rvec v, rvec r);
        void ivec_mmul(const matrix m, const int i[], rvec r);
        void closest_in_cell(ivec g_i, const rvec x, rvec d, ivec closest);
        void rvec_sadd(const rvec a, real b, const rvec c, rvec r);
        void rvec_splus(const real a, const rvec b, rvec r);
        void ivec_positive_modulo(ivec i, ivec n);
        int positive_modulo(int i, int n);
        int ivec2i_mod(ivec iv, ivec e_n);
        int ivec2i(ivec iv, ivec e_n);
        void i2ivec(int i, ivec e_n, ivec iv);
        uint64_t factorial(const int n);
        real pbc_dist(ivec i, ivec n, matrix b, ivec ndx);
        void broadcast_internal();
        void broadcast_expansion();
        void sum_reduce();
        void I2L_enb2v(real &voxel, ivec closest_e, int i_nb, rvec d_r);

        void do_force(const rvec x, rvec f, real map_sigma, volumedata::GridReal &map);
        void voxel_sums(volumedata::GridReal &map);

        int ivec2i_mod_const(ivec iv, ivec e_n);
        void lr_poly_init_tables(real s);
        void integrate_densities_at_expansion_centers(volumedata::GridReal &map);
    private:

        void setupExpansionCenters();
        std::vector<std::unique_ptr<ExpansionCenter> > expansion_centers_;
        std::vector<RVec>                              e_c_;
        matrix                                         box_;
        matrix                                         box_inv_;
        real                                           sigma_;
        int                                            n_I_;
        int                                            maximum_expansion_order_;
        ivec                                           e_n_;
        rvec                                           e_l_;
        int                                            e_N_;
        std::vector<real>                              multinomials_;
        real                                           h_;
        real                                           h_inv_;
        matrix                                         e_cell_;
        real                                           n_sigma_;
        std::vector<t_nb>                              e_nbs_;
        std::shared_ptr<MpiHelper>                     mpi_;
        std::vector < std::vector < real> >            c;
        std::vector < std::vector < real> >            cg;
        std::vector < std::vector < real> >            powers_;
};

}      // namespace gmx
#endif /* IFGT_H_ */
