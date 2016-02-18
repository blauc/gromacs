/*
 * Ifgt.cpp
 *
 *  Created on: Aug 25, 2015
 *      Author: cblau
 */

#include "Ifgt.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/volumedata.h"
#include "gromacs/externalpotential/mpi-helper.h"
#include <algorithm>

namespace gmx
{
ExpansionCenter::ExpansionCenter(Ifgt * ifgt) : ifgt_(ifgt)
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
    return i % 2 == 0;
}

void ExpansionCenter::total_density()
{
    real total_density_ = 0;

    int  i_x;
    int  i_y;

    I_ = coefficients_.begin();

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

real ExpansionCenter::expand(rvec dx)
{
    I_ = coefficients_.begin();
    real result = *I_;
    ++I_;

    int i_x;
    int i_y;

    mn_ = ifgt_->multinomials_.begin();

    *mn_ = 1;
    I_x_ = mn_;
    I_y_ = mn_;
    I_z_ = mn_;

    int o;

    for (o = 1; o <= ifgt_->maximum_expansion_order_; ++o)
    {

        I_pre_ = I_x_;
        I_x_   = mn_ + 1;
        for (i_x = o; i_x >= 1; --i_x)
        {
            for (i_y = o - i_x; i_y >= 0; --i_y)
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
        for (i_y = o; i_y >= 1; --i_y)
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
    mn_  = ifgt_->multinomials_.begin();
    *mn_ = 1;
    I_   = coefficients_.begin();
    I_x_ = mn_;
    I_y_ = mn_;
    I_z_ = mn_;
    for (int o = 1; o <= ifgt_->maximum_expansion_order_; ++o)
    {

        I_pre_ = I_x_;
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


Ifgt::Ifgt(){};

Ifgt::~Ifgt(){};

void Ifgt::sum_reduce()
{
    if (mpi_ != nullptr)
    {
        int i_non_empty;
        for (int i = 0; i < e_N_; i++)
        {
            if (expansion_centers_[i] != nullptr)
            {
                i_non_empty = i;
                mpi_->broadcast(&i_non_empty, 1);
            }

            if (i == i_non_empty)
            {
                if (expansion_centers_[i] == nullptr)
                {
                    expansion_centers_[i] = std::unique_ptr<ExpansionCenter>(new ExpansionCenter(this));
                }
                mpi_->to_reals_buffer(expansion_centers_[i]->coefficients_.data(), expansion_centers_[i]->coefficients_.size());
                mpi_->sum_reduce();
                mpi_->from_reals_buffer(expansion_centers_[i]->coefficients_.data(), expansion_centers_[i]->coefficients_.size());
            }
        }

    }

}

void Ifgt::broadcast_internal()
{
    if (mpi_ != nullptr)
    {
        mpi_->broadcast(&sigma_, 1);
        mpi_->broadcast(&maximum_expansion_order_, 1);
        mpi_->broadcast(&h_, 1);
        mpi_->broadcast(&h_inv_, 1);
        mpi_->broadcast(&n_sigma_, 1);
        mpi_->broadcast(&n_I_, 1);
    }
}

void Ifgt::init(real sigma, int maximum_expansion_order)
{
    sigma_                   = sigma;
    maximum_expansion_order_ = maximum_expansion_order;
    h_                       = sqrt(2) * sigma_;
    h_inv_                   = 1./h_;
    n_sigma_                 = sqrt((maximum_expansion_order_ / 2)) + 1;
    n_I_                     = factorial(maximum_expansion_order_ + DIM) / (factorial(DIM) * factorial(maximum_expansion_order_));
    multinomials_.resize(n_I_);
    powers_.resize(3);

    for (auto &power : powers_)
    {
        power.resize(maximum_expansion_order_+2);
    }
}

real Ifgt::sigma()
{
    return sigma_;
}

/** Integrated density at map over density from IFgt */
void Ifgt::integrate_densities_at_expansion_centers(volumedata::GridReal &map)
{
    ivec closest_e;
    rvec d;
    RVec v;
    int  i_e;

    for (auto &center : expansion_centers_)
    {
        if (center != nullptr)
        {
            center->total_density();
        }
    }

    for (size_t i = 0; i < map.num_gridpoints(); i++)
    {
        if (map.data()[i] > 0)
        {
            coordinate_to_expansioncenter(map.gridpoint_coordinate(i), closest_e, d);
            i_e = ivec2i_mod(closest_e, e_n_);
            if (expansion_centers_[i_e] != nullptr)
            {
                expansion_centers_[i_e]->voxel_density_ += map.data()[i];
            }
        }
    }
}


int Ifgt::ivec2i_mod_const(ivec iv, ivec e_n)
{

    ivec i_tmp = {
        positive_modulo(iv[XX], e_n[XX]),
        positive_modulo(iv[YY], e_n[YY]),
        positive_modulo(iv[ZZ], e_n[ZZ])
    };

    return ivec2i(i_tmp, e_n);
}

void Ifgt::lr_poly_init_tables(real s)
{

    cg.resize(maximum_expansion_order_ + 1);
    c.resize(maximum_expansion_order_ + 1);

    int alpha;

    alpha = 0;
    cg[alpha].resize(alpha + 2);
    c[alpha].resize(alpha + 1);

    cg[alpha][0] = 0;
    cg[alpha][1] = -sqrt(2) / (2 * s * s);

    c[alpha][0] = 1 / sqrt(2);

    ++alpha;

    cg[alpha].resize(alpha + 2);
    c[alpha].resize(alpha + 1);

    cg[alpha][0] = -h_ * sqrt(2) / (4.0 * s * s);
    cg[alpha][1] = 0;
    cg[alpha][2] = h_ * sqrt(2) / (4.0 * s * s * s * s);

    c[alpha][0] = 0;
    c[alpha][1] = -sqrt(2) * h_ / (4.0 * s * s);

    while (alpha < maximum_expansion_order_)
    {
        ++alpha;
        cg[alpha].resize(alpha + 2);
        c[alpha].resize(alpha + 1);

        cg[alpha][0] = 0.5 * (h_ * cg[alpha - 1][1] + (alpha - 1) * cg[alpha - 2][0]);
        for (int i = 1; i < alpha; ++i)
        {
            cg[alpha][i] = 0.5 * ((alpha - 1) * cg[alpha - 2][i] - h_ * cg[alpha - 1][i - 1] / (s * s) + (i + 1) * h_ * cg[alpha - 1][i + 1]);
        }
        cg[alpha][alpha]     = 0;
        cg[alpha][alpha + 1] = -h_ * cg[alpha - 1][alpha] / (2 * s * s);

        c[alpha][0] = 0.5 * (h_ * c[alpha - 1][1] + (alpha - 1) * c[alpha - 2][0]);
        for (int i = 1; i < alpha-1; ++i)
        {
            c[alpha][i] = 0.5 * ((alpha - 1) * c[alpha - 2][i] - h_ * c[alpha - 1][i - 1] / (s * s) + (i + 1) * h_ * c[alpha - 1][i + 1]);
        }
        c[alpha][alpha-1] = 0;
        c[alpha][alpha]   = -h_ * c[alpha - 1][alpha-1] / (2 * s * s);
    }
}



real ExpansionCenter::inner(real * a, real * b, int o)
{
    int  i;
    real result = 0;
    for (i = 0; i <= o; ++i)
    {
        result += a[i] * b[i];
    }
    return result;
}

void ExpansionCenter::do_force_lr_fpoly(ivec m, rvec f, real I)
{
    f[XX] += I * inner(ifgt_->cg[m[XX]].data(), ifgt_->powers_[XX].data(), m[XX] + 1) * inner(ifgt_->c[m[YY]].data(), ifgt_->powers_[YY].data(), m[YY]) * inner(ifgt_->c[m[ZZ]].data(), ifgt_->powers_[ZZ].data(), m[ZZ]);
    f[YY] += I * inner(ifgt_->c[m[XX]].data(), ifgt_->powers_[XX].data(), m[XX]) * inner(ifgt_->cg[m[YY]].data(), ifgt_->powers_[YY].data(), m[YY] + 1) * inner(ifgt_->c[m[ZZ]].data(), ifgt_->powers_[ZZ].data(), m[ZZ]);
    f[ZZ] += I * inner(ifgt_->c[m[XX]].data(), ifgt_->powers_[XX].data(), m[XX]) * inner(ifgt_->c[m[YY]].data(), ifgt_->powers_[YY].data(), m[YY]) * inner(ifgt_->cg[m[ZZ]].data(), ifgt_->powers_[ZZ].data(), m[ZZ] + 1);
}


RVec ExpansionCenter::do_force_lr_low(rvec d)
{
    RVec f {
        0, 0, 0
    };
    ivec m;
    I_ = coefficients_.begin();

    ifgt_->powers_[XX][0] = 1;
    ifgt_->powers_[YY][0] = 1;
    ifgt_->powers_[ZZ][0] = 1;

    for (int o = 1; o <= ifgt_->maximum_expansion_order_ + 1; ++o)
    {
        ifgt_->powers_[XX][o] = ifgt_->powers_[XX][o - 1] * (-d[XX]);
        ifgt_->powers_[YY][o] = ifgt_->powers_[YY][o - 1] * (-d[YY]);
        ifgt_->powers_[ZZ][o] = ifgt_->powers_[ZZ][o - 1] * (-d[ZZ]);
    }

    for (int o = 0; o <= ifgt_->maximum_expansion_order_; ++o)
    {
        for (m[XX] = o; m[XX] >= 1; --m[XX])
        {
            for (m[YY] = o - m[XX]; m[YY] >= 0; --m[YY])
            {
                m[ZZ] = o - m[XX] - m[YY];
                do_force_lr_fpoly(m, f, *I_);
                I_++;
            }
        }
        for (m[YY] = o; m[YY] >= 1; --m[YY])
        {
            m[ZZ] = o - m[XX] - m[YY];
            do_force_lr_fpoly(m, f, *I_);
            I_++;
        }
        m[ZZ] = o;
        do_force_lr_fpoly(m, f, *I_);
        I_++;
    }
    return f;
}

void Ifgt::do_force(const rvec x, rvec f, real map_sigma, volumedata::GridReal &difference_map)
{
    integrate_densities_at_expansion_centers(difference_map);

    ivec closest_e;
    real s = sqrt(sigma_ * sigma_ + map_sigma * map_sigma);
    lr_poly_init_tables(s);

    ivec e_icur;
    int  e_i;
    rvec dist;
    rvec d_curr;
    rvec c_r;
    RVec f_low;

    coordinate_to_expansioncenter(x, closest_e, dist);

    for (int i_nb = 0; i_nb < e_N_; ++i_nb)
    {

        ivec_add(closest_e, e_nbs_[i_nb].ndx, e_icur);
        e_i = ivec2i_mod_const(e_icur, e_n_);

        if (expansion_centers_[e_i] != nullptr)
        {

            ivec_mmul(e_cell_, e_nbs_[i_nb].ndx, c_r);
            rvec_sub(dist, c_r, d_curr);

            f_low = expansion_centers_[e_i]->do_force_lr_low(d_curr);

            rvec_splus(-exp(-norm2(d_curr) / (s * s)) * expansion_centers_[e_i]->voxel_density_ / expansion_centers_[e_i]->total_density_, f_low, f);

        }
    }

    svmul(det(e_cell_) / (s * s * s), f, f);
}

void Ifgt::broadcast_expansion()
{
    if (mpi_ != nullptr)
    {
        if (!mpi_->isMaster())
        {
            expansion_centers_.clear();
            expansion_centers_.resize(e_N_);
        }
        for (int i = 0; i < e_N_; i++)
        {

            if (mpi_->isMaster() && expansion_centers_[i] != nullptr)
            {
                expansion_centers_[i] = std::unique_ptr<ExpansionCenter>(new ExpansionCenter(this));
                mpi_->broadcast(expansion_centers_[i]->coefficients_.data(), expansion_centers_[i]->coefficients_.size());
            }
        }
    }
}

real Ifgt::pbc_dist(ivec i, ivec n, matrix b, ivec ndx)
{

    ivec i_curr;
    ivec im;
    real result = 1e10;
    rvec dist;
    for (im[XX] = -n[XX]; im[XX] <= n[XX]; im[XX] += n[XX])
    {
        for (im[YY] = -n[YY]; im[YY] <= n[YY]; im[YY] += n[YY])
        {
            for (im[ZZ] = -n[ZZ]; im[ZZ] <= n[ZZ]; im[ZZ] += n[ZZ])
            {
                ivec_add(i, im, i_curr);
                ivec_mmul(b, i_curr, dist);
                if (norm2(dist) < result*result)
                {
                    result = norm(dist);
                    copy_ivec(i_curr, ndx);
                }
            }
        }
    }

    return result;

}


void Ifgt::i2ivec(int i, ivec e_n, ivec iv)
{

    iv[XX] = i % e_n[XX];
    i     /= e_n[XX];
    iv[YY] = i % e_n[YY];
    i     /= e_n[YY];
    iv[ZZ] = i % e_n[ZZ];

}
void Ifgt::set_mpi_helper(std::shared_ptr<MpiHelper> mpi)
{
    mpi_ = mpi;
}

uint64_t Ifgt::factorial(const int n)
{
    static const uint64_t           result[21] = { 1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600, 6227020800, 87178291200, 1307674368000, 20922789888000, 355687428096000, 6402373705728000, 121645100408832000, 2432902008176640000 };
    return n < 21 ? result[n] : n * factorial(n - 1);
}


void Ifgt::matrix_inverse(const matrix m, matrix r)
{
    real det = m[XX][XX] * m[YY][YY] * m[ZZ][ZZ]
        + m[XX][YY] * m[YY][ZZ] * m[ZZ][XX]
        + m[XX][ZZ] * m[YY][XX] * m[ZZ][YY]
        - m[XX][XX] * m[YY][ZZ] * m[ZZ][YY]
        - m[XX][YY] * m[YY][XX] * m[ZZ][ZZ]
        - m[XX][ZZ] * m[YY][YY] * m[ZZ][XX];
    r[XX][XX] = (m[YY][YY] * m[ZZ][ZZ] - m[YY][ZZ] * m[ZZ][YY]) / det;
    r[XX][YY] = (m[XX][ZZ] * m[ZZ][YY] - m[XX][YY] * m[ZZ][ZZ]) / det;
    r[XX][ZZ] = (m[XX][YY] * m[YY][ZZ] - m[XX][ZZ] * m[YY][YY]) / det;
    r[YY][XX] = (m[YY][ZZ] * m[ZZ][XX] - m[YY][XX] * m[ZZ][ZZ]) / det;
    r[YY][YY] = (m[XX][XX] * m[ZZ][ZZ] - m[XX][ZZ] * m[ZZ][XX]) / det;
    r[YY][ZZ] = (m[XX][ZZ] * m[YY][XX] - m[XX][XX] * m[YY][ZZ]) / det;
    r[ZZ][XX] = (m[YY][XX] * m[ZZ][YY] - m[YY][YY] * m[ZZ][XX]) / det;
    r[ZZ][YY] = (m[XX][YY] * m[ZZ][XX] - m[XX][XX] * m[ZZ][YY]) / det;
    r[ZZ][ZZ] = (m[XX][XX] * m[YY][YY] - m[XX][YY] * m[YY][XX]) / det;
}

void Ifgt::vmmul(const matrix m, const rvec v, rvec r)
{
    r[XX] = m[XX][XX] * v[XX] + m[YY][XX] * v[YY] + m[ZZ][XX] * v[ZZ];
    r[YY] = m[XX][YY] * v[XX] + m[YY][YY] * v[YY] + m[ZZ][YY] * v[ZZ];
    r[ZZ] = m[XX][ZZ] * v[XX] + m[YY][ZZ] * v[YY] + m[ZZ][ZZ] * v[ZZ];
}

void Ifgt::set_box(const matrix box)
{
    copy_mat(box, box_);
    matrix_inverse(box_, box_inv_);

    e_n_[XX] = ceil((norm(box_[XX])) / (n_sigma_ * sigma_));
    e_n_[YY] = ceil((norm(box_[YY])) / (n_sigma_ * sigma_));
    e_n_[ZZ] = ceil((norm(box_[ZZ])) / (n_sigma_ * sigma_));

    e_l_[XX] = norm(box_[XX]) / ((real) e_n_[XX]);
    e_l_[YY] = norm(box_[YY]) / ((real) e_n_[YY]);
    e_l_[ZZ] = norm(box_[ZZ]) / ((real) e_n_[ZZ]);

    e_N_ = e_n_[XX] * e_n_[YY] * e_n_[ZZ];
    expansion_centers_.clear();
    expansion_centers_.resize(e_N_);

    copy_mat(box_, e_cell_);
    svmul(1./(real) e_n_[XX], e_cell_[XX], e_cell_[XX]);
    svmul(1./(real) e_n_[YY], e_cell_[YY], e_cell_[YY]);
    svmul(1./(real) e_n_[ZZ], e_cell_[ZZ], e_cell_[ZZ]);

    ivec iv;
    e_c_.resize(e_N_, {0, 0, 0});
    for (int i = 0; i < e_N_; ++i)
    {
        i2ivec(i, e_n_, iv);
        ivec_mmul(e_cell_, iv, e_c_[i]);
    }

    e_nbs_.resize(e_N_);
    ivec nb;
    std::vector<t_nb>::iterator i = e_nbs_.begin();

    for (nb[XX] = 0; nb[XX] < e_n_[XX]; ++nb[XX])
    {
        for (nb[YY] = 0; nb[YY] < e_n_[YY]; ++nb[YY])
        {
            for (nb[ZZ] = 0; nb[ZZ] < e_n_[ZZ]; ++nb[ZZ])
            {

                i->d = pbc_dist(nb, e_n_, e_cell_, i->ndx);
                ++i;
            }
        }
    }

    std::sort(e_nbs_.begin(), e_nbs_.end(), [](const t_nb &a, const t_nb &b) {return a.d > b.d; });

}

int Ifgt::ivec2i_mod(ivec iv, ivec e_n)
{
    ivec_positive_modulo(iv, e_n);
    return ivec2i(iv, e_n);
}

int Ifgt::ivec2i(ivec iv, ivec e_n)
{
    return iv[XX] + iv[YY] * e_n[XX] + iv[ZZ] * e_n[XX] * e_n[YY];
}

void Ifgt::compress(const rvec x, real weight)
{

    int    e_i;
    ivec   closest_e;
    rvec   dist;

    coordinate_to_expansioncenter(x, closest_e, dist);
    e_i = ivec2i_mod(closest_e, e_n_);

    if (expansion_centers_[e_i] == nullptr)
    {
        expansion_centers_[e_i] = std::unique_ptr<ExpansionCenter>(new ExpansionCenter(this));
    }
    expansion_centers_[e_i]->compress(dist, weight);
}

/**
 * find residence grid cell, calc distance to all corners, report closest corner
 * modulo out pbc
 */
void Ifgt::coordinate_to_expansioncenter(const rvec r, ivec closest_e, rvec d)
{
    ivec g_i; //< the front bottom left corner of the grid cell
    grid_index_ivec(r, g_i);
    closest_in_cell(g_i, r, d, closest_e);
    svmul(h_inv_, d, d);
    ivec_positive_modulo(closest_e, e_n_);
}

int Ifgt::positive_modulo(int i, int n)
{
    return (i % n + n) % n;
}

void Ifgt::ivec_positive_modulo(ivec i, ivec n)
{
    i[XX] = positive_modulo(i[XX], n[XX]);
    i[YY] = positive_modulo(i[YY], n[YY]);
    i[ZZ] = positive_modulo(i[ZZ], n[ZZ]);
}


void Ifgt::grid_index_ivec(const rvec x, ivec i)
{
    rvec g;
    vmmul(box_inv_, x, g);
    i[XX] = (int) floor(g[XX] * (real) n_g_[XX]);
    i[YY] = (int) floor(g[YY] * (real) n_g_[YY]);
    i[ZZ] = (int) floor(g[ZZ] * (real) n_g_[ZZ]);
}


void Ifgt::closest_in_cell(ivec g_i, const rvec x, rvec d, ivec closest)
{
    ivec nb;
    rvec r;
    rvec curr_r = {
        0,
        0,
        0
    };
    real dist = 1e20;
    ivec_mmul(e_cell_, g_i, r);

    for (nb[XX] = 0; nb[XX] <= 1; ++nb[XX])
    {
        for (nb[YY] = 0; nb[YY] <= 1; ++nb[YY])
        {
            for (nb[ZZ] = 0; nb[ZZ] <= 1; ++nb[ZZ])
            {

                rvec_sadd(r, (real) nb[XX], e_cell_[XX], curr_r);
                rvec_splus((real) nb[YY], e_cell_[YY], curr_r);
                rvec_splus((real) nb[ZZ], e_cell_[ZZ], curr_r);

                // FGT_FAVOUR_FIRST avoids artifacts when using
                // a single expansion center and periodic boundaries
                // but non pbc expansion
                // with FGT_FAVOURFIRST > 0, voxels/atoms right in between
                // two expansion centers are favourably assigned to the
                // bottom left center and not its (periodic) neighbours

#define FGT_FAVOUR_FIRST 1e-5

                if (distance2(curr_r, x) < dist*dist + FGT_FAVOUR_FIRST)
                {
                    dist = sqrt(distance2(curr_r, x));
                    ivec_add(g_i, nb, closest);
                    rvec_sub(x, curr_r, d);
                }
            }
        }
    }
}

void Ifgt::ivec_mmul(const matrix m, const int i[], rvec r)
{
    r[XX] = m[XX][XX] * (real) i[XX] + m[YY][XX] * (real) i[YY]
        + m[ZZ][XX] * (real) i[ZZ];
    r[YY] = m[XX][YY] * (real) i[XX] + m[YY][YY] * (real) i[YY]
        + m[ZZ][YY] * (real) i[ZZ];
    r[ZZ] = m[XX][ZZ] * (real) i[XX] + m[YY][ZZ] * (real) i[YY]
        + m[ZZ][ZZ] * (real) i[ZZ];
}


void Ifgt::rvec_sadd(const rvec a, real b, const rvec c, rvec r)
{
    r[XX] = a[XX] + b * c[XX];
    r[YY] = a[YY] + b * c[YY];
    r[ZZ] = a[ZZ] + b * c[ZZ];
}

void Ifgt::rvec_splus(const real a, const rvec b, rvec r)
{
    r[XX] += a * b[XX];
    r[YY] += a * b[YY];
    r[ZZ] += a * b[ZZ];
}


void Ifgt::I2L_enb2v(volumedata::GridReal &map, ivec closest_e, int i_nb, int i,
                     rvec d_r)
{
    ivec e_icur;
    rvec g_r;
    rvec d_curr;
    int e_i;

    ivec_add(closest_e, e_nbs_[i_nb].ndx, e_icur);
    e_i = ivec2i_mod(e_icur, e_n_);

    if (expansion_centers_[e_i] != nullptr)
    {

        ivec_mmul(e_cell_, e_nbs_[i_nb].ndx, g_r);
        rvec_sadd(d_r, -h_inv_, g_r, d_curr);

        map.data()[i] += exp(-norm2(d_curr))
            * expansion_centers_[e_i]->expand(d_curr);
    }
}

void Ifgt::compress_density(volumedata::GridReal &map)
{
    for (size_t i = 0; i < map.num_gridpoints(); ++i)
    {
        compress(map.gridpoint_coordinate(i), map.data()[i]);
    }
}

void Ifgt::expand_at_ref(volumedata::GridReal &map, volumedata::GridReal &ref, real threshold)
{
    auto refvox = ref.data();
    auto mapvox = map.data();
    ivec closest_e;
    real minimum_distance;
    rvec d_r;
    int i_nb;
    for (size_t i = 0; i < refvox.size(); i++)
    {
        if (refvox[i] > threshold)
        {

            coordinate_to_expansioncenter(map.gridpoint_coordinate(i), closest_e, d_r);
            minimum_distance = 3 * sigma_ + norm(d_r);

            // add contributinos to voxel from all expansion centers within minimum_distance
            // go beyond minmum_distance if voxel value would be zero
            for (i_nb = 0; (i_nb < e_N_) && ((e_nbs_[i_nb].d < minimum_distance) || (mapvox[i] == 0)); ++i_nb)
            {
                I2L_enb2v(map, closest_e, i_nb, i, d_r);

                // all expansion centers within the same range should also contribute
                for (++i_nb; (i_nb < e_N_) && (e_nbs_[i_nb - 1].d + 1e-10 <= e_nbs_[i_nb].d); ++i_nb)
                {
                    I2L_enb2v(map, closest_e, i_nb, i, d_r);
                }
            }
        }
    }
}

} /* gmx */
