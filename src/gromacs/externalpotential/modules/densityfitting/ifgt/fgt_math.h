#include "stdint.h"
#include "stdio.h"
#include "math.h"

#ifndef XX
#define XX 0
#endif

#ifndef YY
#define YY 1
#endif

#ifndef ZZ
#define ZZ 2
#endif

#include "fgt_types.h"

void fgt_copy_rvec(const rvec a, rvec b);
void fgt_clear_rvec(rvec r);
real fgt_gamma_1_2(unsigned int i);
int is_power_of_two(int x);
real norm_diff(rvec a, rvec b);
void fgt_minv(matrix m, matrix r);
uint64_t fgt_factorial(int n);
void ivec_smul(ivec a, real s, rvec r);
void rvec_smul(rvec r, const real a);
void rvec_sdiv(rvec r, const real a);
void matrix_copy(matrix m, matrix r);
void ivec_copy(ivec from, ivec to);
void rvec_diff(const rvec a, const rvec b, rvec r);
void print_rvec(rvec p);
real fgt_dist(rvec a, rvec b);
void fgt_rvec_add(const rvec a, const rvec b, rvec r);
void rvec_sadd(const rvec a, real b, const rvec c, rvec r);
void rvec_splus(real a, rvec b, rvec r);

void fgt_ivec_add(ivec a, ivec b, ivec r);
void fgt_ivec_lmul(ivec i, matrix m, rvec r);
void fgt_ivec_mmul(matrix m, ivec i, rvec r);
int ivec2i_mod(ivec iv, ivec e_n);
int positive_modulo(int i, int n);
int ivec2i_mod_const(ivec iv, ivec e_n);
real fgt_norm(rvec r);
void ivec_positive_modulo(ivec i, ivec n);
void print_ivec(ivec p);
void vlmul(rvec v, matrix m, rvec r);
void vmmul(matrix m, rvec v, rvec r);
int is_even(int i);
