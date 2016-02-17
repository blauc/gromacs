/*
 * fgt_vmath.c
 *
 *  Created on: Aug 4, 2015
 *      Author: cblau
 */

#include "fgt_math.h"

void vmmul(matrix m, rvec v, rvec r)
{
    r[XX] = m[XX][XX] * v[XX] + m[YY][XX] * v[YY] + m[ZZ][XX] * v[ZZ];
    r[YY] = m[XX][YY] * v[XX] + m[YY][YY] * v[YY] + m[ZZ][YY] * v[ZZ];
    r[ZZ] = m[XX][ZZ] * v[XX] + m[YY][ZZ] * v[YY] + m[ZZ][ZZ] * v[ZZ];
}

void vlmul(rvec v, matrix m, rvec r)
{
    r[XX] = m[XX][XX] * v[XX] + m[XX][YY] * v[YY] + m[XX][ZZ] * v[ZZ];
    r[YY] = m[YY][XX] * v[XX] + m[YY][YY] * v[YY] + m[YY][ZZ] * v[ZZ];
    r[ZZ] = m[ZZ][XX] * v[XX] + m[ZZ][YY] * v[YY] + m[ZZ][ZZ] * v[ZZ];
}

void print_ivec(ivec p)
{
    fprintf(stderr, "%d %d %d \n", p[XX], p[YY], p[ZZ]);
}

inline void ivec_positive_modulo(ivec i, ivec n)
{
    i[XX] = positive_modulo(i[XX], n[XX]);
    i[YY] = positive_modulo(i[YY], n[YY]);
    i[ZZ] = positive_modulo(i[ZZ], n[ZZ]);
}


real fgt_norm(rvec r)
{
    return sqrt(r[XX] * r[XX] + r[YY] * r[YY] + r[ZZ] * r[ZZ]);
}


void fgt_ivec_lmul(ivec i, matrix m, rvec r)
{
    r[XX] = m[XX][XX] * (real) i[XX] + m[XX][YY] * (real) i[YY]
        + m[XX][ZZ] * (real) i[ZZ];
    r[YY] = m[YY][XX] * (real) i[XX] + m[YY][YY] * (real) i[YY]
        + m[YY][ZZ] * (real) i[ZZ];
    r[ZZ] = m[ZZ][XX] * (real) i[XX] + m[ZZ][YY] * (real) i[YY]
        + m[ZZ][ZZ] * (real) i[ZZ];
}


void rvec_splus(real a, rvec b, rvec r)
{
    r[XX] += a * b[XX];
    r[YY] += a * b[YY];
    r[ZZ] += a * b[ZZ];
}


void fgt_rvec_add(const rvec a, const rvec b, rvec r)
{

    r[XX] = a[XX] + b[XX];
    r[YY] = a[YY] + b[YY];
    r[ZZ] = a[ZZ] + b[ZZ];
}

real fgt_dist(rvec a, rvec b)
{
    return sqrt(
            fgt_sqr(
                    a[XX] - b[XX]) + fgt_sqr(a[YY] - b[YY]) + fgt_sqr(a[ZZ] - b[ZZ]));
}
void print_rvec(rvec p)
{
    fprintf(stderr, "%f %f %f \n", p[XX], p[YY], p[ZZ]);
}

void rvec_diff(const rvec a, const rvec b, rvec r)
{

    r[XX] = a[XX] - b[XX];
    r[YY] = a[YY] - b[YY];
    r[ZZ] = a[ZZ] - b[ZZ];
}

void ivec_copy(ivec from, ivec to)
{

    to[XX] = from[XX];
    to[YY] = from[YY];
    to[ZZ] = from[ZZ];
}

void matrix_copy(matrix m, matrix r)
{

    r[XX][XX] = m[XX][XX];
    r[XX][YY] = m[XX][YY];
    r[XX][ZZ] = m[XX][ZZ];
    r[YY][XX] = m[YY][XX];
    r[YY][YY] = m[YY][YY];
    r[YY][ZZ] = m[YY][ZZ];
    r[ZZ][XX] = m[ZZ][XX];
    r[ZZ][YY] = m[ZZ][YY];
    r[ZZ][ZZ] = m[ZZ][ZZ];
}

void rvec_sdiv(rvec r, const real a)
{

    r[XX] /= a;
    r[YY] /= a;
    r[ZZ] /= a;
}

void rvec_smul(rvec r, const real a)
{

    r[XX] *= a;
    r[YY] *= a;
    r[ZZ] *= a;
}
int ivec2i(ivec iv, ivec e_n)
{

    return iv[XX] + iv[YY] * e_n[XX] + iv[ZZ] * e_n[XX] * e_n[YY];
}

void fgt_clear_rvec(rvec r)
{
    r[XX] = 0;
    r[YY] = 0;
    r[ZZ] = 0;
}

real norm_diff(rvec a, rvec b)
{
    return fgt_sqr(a[XX]-b[XX]) + fgt_sqr(a[YY]-b[YY]) + fgt_sqr(a[ZZ]-b[ZZ]);
}
void rvec_sadd(const rvec a, real b, const rvec c, rvec r)
{

    r[XX] = a[XX] + b * c[XX];
    r[YY] = a[YY] + b * c[YY];
    r[ZZ] = a[ZZ] + b * c[ZZ];
}
void fgt_ivec_add(ivec a, ivec b, ivec r)
{

    r[XX] = a[XX] + b[XX];
    r[YY] = a[YY] + b[YY];
    r[ZZ] = a[ZZ] + b[ZZ];
}

void fgt_ivec_mmul(matrix m, ivec i, rvec r)
{
    r[XX] = m[XX][XX] * (real) i[XX] + m[YY][XX] * (real) i[YY]
        + m[ZZ][XX] * (real) i[ZZ];
    r[YY] = m[XX][YY] * (real) i[XX] + m[YY][YY] * (real) i[YY]
        + m[ZZ][YY] * (real) i[ZZ];
    r[ZZ] = m[XX][ZZ] * (real) i[XX] + m[YY][ZZ] * (real) i[YY]
        + m[ZZ][ZZ] * (real) i[ZZ];
}

int ivec2i_mod(ivec iv, ivec e_n)
{
    ivec_positive_modulo(iv, e_n);
    return ivec2i(iv, e_n);
}

inline int positive_modulo(int i, int n)
{
    return (i % n + n) % n;
}

int ivec2i_mod_const(ivec iv, ivec e_n)
{

    ivec i_tmp = {
        positive_modulo(iv[XX], e_n[XX]),
        positive_modulo(iv[YY], e_n[YY]),
        positive_modulo(iv[ZZ], e_n[ZZ])
    };

    return ivec2i(i_tmp, e_n);

}


void ivec_smul(ivec a, real s, rvec r)
{

    r[XX] = (real) a[XX] * s;
    r[YY] = (real) a[YY] * s;
    r[ZZ] = (real) a[ZZ] * s;
}

uint64_t fgt_factorial(int n)
{
    static const uint64_t result[14] = {
        1,
        1,
        2,
        6,
        24,
        120,
        720,
        5040,
        40320,
        362880,
        3628800,
        39916800,
        479001600,
        6227020800
    };
    return n < 14 ? result[n] : n * fgt_factorial(n - 1);
}

void fgt_minv(matrix m, matrix r)
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



int is_even(int i)
{
    return i % 2 == 0;
}

int is_power_of_two(int x)
{
    return ((x != 0) && !(x & (x - 1)));
}


/*
 * Gamma((i+1)/2)
 */
real fgt_gamma_1_2(unsigned int i)
{
    static const real result[11] = {
        1.77245385090552,
        1.00000000000000,
        0.886226925452758,
        1.00000000000000,
        1.32934038817914,
        2.00000000000000,
        3.32335097044784,
        6.00000000000000,
        11.6317283965674,
        24.0000000000000,
        52.3427777845535
    };

    return i < 11 ? result[i] : fgt_gamma_1_2(i - 2) * (i - 1) / 2;
}



void fgt_copy_rvec(const rvec a, rvec b)
{
    b[XX] = a[XX];
    b[YY] = a[YY];
    b[ZZ] = a[ZZ];
}
