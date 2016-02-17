/*
 * fgt.c
 *
 *  Created on: Nov 21, 2014
 *      Author: cblau
 */

#ifndef _FGT_h

#include "fgt.h"
#include "fgt_math.h"

#endif

void fgt_grid_index_ivec(rvec x, ivec n_g, matrix box_inv, ivec i)
{
    rvec g;
    vmmul(box_inv, x, g);
    i[XX] = (int) floor(g[XX] * (real) n_g[XX]);
    i[YY] = (int) floor(g[YY] * (real) n_g[YY]);
    i[ZZ] = (int) floor(g[ZZ] * (real) n_g[ZZ]);

}

void fgt_closest_in_cell(ivec g_i, matrix cell, rvec x, rvec d, ivec closest)
{
    ivec nb;
    rvec r;
    rvec curr_r = {
        0,
        0,
        0
    };
    real dist = MAX_REAL;
    fgt_ivec_mmul(cell, g_i, r);

    for (nb[XX] = 0; nb[XX] <= 1; ++nb[XX])
    {
        for (nb[YY] = 0; nb[YY] <= 1; ++nb[YY])
        {
            for (nb[ZZ] = 0; nb[ZZ] <= 1; ++nb[ZZ])
            {

                rvec_sadd(r, (real) nb[XX], cell[XX], curr_r);
                rvec_splus((real) nb[YY], cell[YY], curr_r);
                rvec_splus((real) nb[ZZ], cell[ZZ], curr_r);

                // FGT_FAVOUR_FIRST avoids artifacts when using
                // a single expansion center and periodic boundaries
                // but non pbc expansion
                // with FGT_FAVOURFIRST > 0, voxels/atoms right in between
                // two expansion centers are favourably assigned to the
                // bottom left center and not its (periodic) neighbours

#define FGT_FAVOUR_FIRST 1e-10

                if (fgt_dist(curr_r, x) < dist + FGT_FAVOUR_FIRST)
                {
                    dist = fgt_dist(curr_r, x);
                    fgt_ivec_add(g_i, nb, closest);
                    rvec_diff(x, curr_r, d);
                }
            }
        }
    }
}

/**
 * find residence grid cell, calc distance to all corners, report closest corner
 * modulo out pbc
 */
void fgt_r2e(rvec r, t_fgt * ft, ivec closest_e, rvec d)
{
    ivec g_i; //< the front bottom left corner of the grid cell
    fgt_grid_index_ivec(r, ft->e_n, ft->b_inv, g_i);
    fgt_closest_in_cell(g_i, ft->e_cell, r, d, closest_e);
    rvec_sdiv(d, ft->h);
    ivec_positive_modulo(closest_e, ft->e_n);
}

void i2ivec(int i, ivec e_n, ivec iv)
{

    iv[XX] = i % e_n[XX];
    i     /= e_n[XX];
    iv[YY] = i % e_n[YY];
    i     /= e_n[YY];
    iv[ZZ] = i % e_n[ZZ];

}

static inline real norm2(rvec r)
{
    return r[XX] * r[XX] + r[YY] * r[YY] + r[ZZ] * r[ZZ];
}

int fgt_nb_compare(const void *a, const void * b)
{
    return ((t_nb*) a)->d > ((t_nb*) b)->d;
}

real pbc_dist(ivec i, ivec n, matrix b, ivec ndx)
{

    ivec i_curr;
    ivec im;
    real result = MAX_REAL;
    rvec dist;
    for (im[XX] = -n[XX]; im[XX] <= n[XX]; im[XX] += n[XX])
    {
        for (im[YY] = -n[YY]; im[YY] <= n[YY]; im[YY] += n[YY])
        {
            for (im[ZZ] = -n[ZZ]; im[ZZ] <= n[ZZ]; im[ZZ] += n[ZZ])
            {
                fgt_ivec_add(i, im, i_curr);
                fgt_ivec_mmul(b, i_curr, dist);
                if (fgt_norm(dist) < result)
                {
                    result = fgt_norm(dist);
                    ivec_copy(i_curr, ndx);
                }
            }
        }
    }

    return result;

}

int fgt_init_e_neighbours(t_fgt * ft)
{
    ft->e_nbs = malloc(ft->e_N * sizeof(t_nb));
    ivec nb;
    t_nb * i = ft->e_nbs;

    for (nb[XX] = 0; nb[XX] < ft->e_n[XX]; ++nb[XX])
    {
        for (nb[YY] = 0; nb[YY] < ft->e_n[YY]; ++nb[YY])
        {
            for (nb[ZZ] = 0; nb[ZZ] < ft->e_n[ZZ]; ++nb[ZZ])
            {

                i->d = pbc_dist(nb, ft->e_n, ft->e_cell, i->ndx);
                ++i;
            }
        }
    }

    qsort(ft->e_nbs, ft->e_N, sizeof(t_nb), &fgt_nb_compare);
    return EXIT_SUCCESS;

}

/**
 * Initalizing expansion boxes with pbc.
 * Uniformly space-filling grid by equidistant spacing along box vectors
 */
int fgt_init_e_boxes(t_fgt * ft, matrix box)
{

    int i;
    ivec iv;

    ft->e_n[XX] = ceil((fgt_norm(box[XX])) / (ft->n_sigma * ft->sigma));
    ft->e_n[YY] = ceil((fgt_norm(box[YY])) / (ft->n_sigma * ft->sigma));
    ft->e_n[ZZ] = ceil((fgt_norm(box[ZZ])) / (ft->n_sigma * ft->sigma));

    ft->e_l[XX] = fgt_norm(box[XX]) / ((real) ft->e_n[XX]);
    ft->e_l[YY] = fgt_norm(box[YY]) / ((real) ft->e_n[YY]);
    ft->e_l[ZZ] = fgt_norm(box[ZZ]) / ((real) ft->e_n[ZZ]);

    ft->e_N = ft->e_n[XX] * ft->e_n[YY] * ft->e_n[ZZ];

    ft->e_c = malloc(ft->e_N * sizeof(rvec));

    matrix_copy(ft->box, ft->e_cell);
    rvec_sdiv(ft->e_cell[XX], (real) ft->e_n[XX]);
    rvec_sdiv(ft->e_cell[YY], (real) ft->e_n[YY]);
    rvec_sdiv(ft->e_cell[ZZ], (real) ft->e_n[ZZ]);

    for (i = 0; i < ft->e_N; ++i)
    {
        i2ivec(i, ft->e_n, iv);
        fgt_ivec_mmul(ft->e_cell, iv, ft->e_c[i]);
    }

    if (ft->verbose)
    {
        fprintf(stderr, "Expanding into %dx%dx%d = %d boxes ", ft->e_n[XX],
                ft->e_n[YY], ft->e_n[ZZ], ft->e_N);
        fprintf(stderr, "with side length %5.3f,%5.3f,%5.3f and h = %5.3f.\n",
                ft->e_l[XX], ft->e_l[YY], ft->e_l[ZZ], ft->h);
        fprintf(stderr, "First box centerd at %5.3f,%5.3f,%5.3f, "
                "last box centered at %5.3f,%5.3f,%5.3f.\n", ft->e_c[0][XX],
                ft->e_c[0][YY], ft->e_c[0][ZZ], ft->e_c[ft->e_N - 1][XX],
                ft->e_c[ft->e_N - 1][YY], ft->e_c[ft->e_N - 1][ZZ]);
    }

    fgt_init_e_neighbours(ft);

    return ((ft->e_c) == NULL) ? EXIT_FAILURE : EXIT_SUCCESS;
}
void fgt_normalize(t_map * map)
{
    int i;
    real norm = 0;
    real avg;
    for (i = 0; i < map->map_dim[XX] * map->map_dim[YY] * map->map_dim[ZZ];
         ++i)
    {
        norm += map->vox[i];
    }
    avg = norm
        / (real) (map->map_dim[XX] * map->map_dim[YY] * map->map_dim[ZZ]);
    fprintf(stderr, "Map average %f ", avg);
    for (i = 0; i < map->map_dim[XX] * map->map_dim[YY] * map->map_dim[ZZ];
         ++i)
    {
        map->vox[i] -= avg;
        map->vox[i] /= norm;
        map->vox[i] /= (map->spacing * map->spacing * map->spacing);
    }
}
int fgt_alloc_I(t_fgt * ft)
{

    ft->I            = calloc(ft->e_N, sizeof(real*));
    ft->multinomials = malloc(sizeof(real) * ft->n_I);

    return (ft->I == NULL) || (ft->multinomials == NULL) ? EXIT_FAILURE :
           EXIT_SUCCESS;

}

int init_fgt(t_fgt **ft, rvec *y, real *f, int natom, matrix box, real sigma,
             int m, int verbose, int y_own, int f_own)
{

    int i;

    *ft = malloc(sizeof(t_fgt));

    (*ft)->sigma = sigma;

    (*ft)->h = sqrt(2) * sigma;
    (*ft)->n = natom;

    (*ft)->y     = y;
    (*ft)->y_own = y_own;

    (*ft)->f           = f;
    (*ft)->box[XX][XX] = box[XX][XX];
    (*ft)->box[XX][YY] = box[XX][YY];
    (*ft)->box[XX][ZZ] = box[XX][ZZ];
    (*ft)->box[YY][XX] = box[YY][XX];
    (*ft)->box[YY][YY] = box[YY][YY];
    (*ft)->box[YY][ZZ] = box[YY][ZZ];
    (*ft)->box[ZZ][XX] = box[ZZ][XX];
    (*ft)->box[ZZ][YY] = box[ZZ][YY];
    (*ft)->box[ZZ][ZZ] = box[ZZ][ZZ];

    fgt_minv((*ft)->box, (*ft)->b_inv);

    (*ft)->verbose = verbose;
    (*ft)->m       = m;
    (*ft)->n_sigma = sqrt(((*ft)->m / 2)) + 1;
    (*ft)->f_own   = f_own;
    if (f == NULL)
    {
        if ((*ft)->verbose)
        {
            fprintf(stderr,
                    "[ NOTE ] No weights given - using all weights one.\n");
        }
        (*ft)->f_own = (1 == 1);
        (*ft)->f     = malloc((*ft)->n * sizeof(real));
        for (i = 0; i < (*ft)->n; ++i)
        {
            (*ft)->f[i] = 1;
        }
    }

    if ((*ft)->verbose)
    {
        fprintf(stderr, "Initializing fgt for %d atoms, sigma of %f using ",
                (*ft)->n, (*ft)->sigma);
        fprintf(stderr, "\n");
    }

    (*ft)->n_I = fgt_factorial((*ft)->m + DIMENSIONS)
        / (fgt_factorial(DIMENSIONS) * fgt_factorial((*ft)->m));

    fgt_init_e_boxes(*ft, box);

    if (fgt_alloc_I(*ft) == EXIT_FAILURE)
    {
        fprintf(stderr,
                "ERROR: Failed to allocate memory for expansion coefficients. \n");
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

void fgt_reset_I(t_fgt * ft)
{
    int e_i;
    for (e_i = 0; e_i < ft->e_N; ++e_i)
    {
        if (ft->I[e_i] != NULL)
        {
            free(ft->I[e_i]);
            ft->I[e_i] = NULL;
        }
    }
}

int fgt_S2I(t_fgt * ft)
{
    if (ft->verbose)
    {
        fprintf(stderr, "Compressing S2I ... ");
    }

    if (ft->I != NULL)
    {
        if (ft->verbose)
        {
            fprintf(stderr,
                    " I non-empty. Clearing previous I. Proceeding ... ");
        }
        fgt_reset_I(ft);
    }

    real * I = NULL;
    int e_i;
    ivec m;

    int i_x;
    ivec closest_e;
    rvec dist;

    real pre;
    real *I_x;
    real *I_y;
    real *I_z;
    real *I_pre;
    int o;

    for (i_x = 0; i_x < ft->n; ++i_x)
    {

        fgt_r2e(ft->y[i_x], ft, closest_e, dist);

        e_i = ivec2i_mod(closest_e, ft->e_n);

        if (ft->I[e_i] == NULL)
        {
            ft->I[e_i] = calloc(ft->n_I, sizeof(real));
            if (ft->I[e_i] == NULL)
            {
                fprintf(stderr,
                        "\n Cannot allocate memory for expansion coefficients. "
                        "Aborting. \n");
                return EXIT_FAILURE;
            }
        }

        pre = ft->f[i_x] * exp(-norm2(dist));

        I = ft->I[e_i];
        real * mn = ft->multinomials;
        *mn = 1;
        *I += pre;
        I_x = mn;
        I_y = mn;
        I_z = mn;
        for (o = 1; o <= ft->m; ++o)
        {

            I_pre = I_x;
            I_x   = mn + 1;
            for (m[XX] = o; m[XX] >= 1; --m[XX])
            {
                for (m[YY] = o - m[XX]; m[YY] >= 0; --m[YY])
                {

                    m[ZZ] = o - m[XX] - m[YY];

                    ++mn;
                    ++I;
                    *mn = *I_pre * 2 * dist[XX] / (real) m[XX];
                    *I += pre * *mn;
                    ++I_pre;

                }
            }

            I_pre = I_y;
            I_y   = mn + 1;
            for (m[YY] = o; m[YY] >= 1; --m[YY])
            {

                ++mn;
                ++I;
                *mn = *I_pre * 2 * dist[YY] / (real) m[YY];
                *I += pre * *mn;
                ++I_pre;

            }

            I_pre = I_z;
            I_z   = mn + 1;
            m[ZZ] = o;
            ++mn;
            ++I;
            *mn = *I_pre * 2 * dist[ZZ] / (real) m[ZZ];
            *I += pre * *mn;

        }
    }

    if (ft->verbose)
    {
        fprintf(stderr, "done.\n");
    }

    return EXIT_SUCCESS;
}

static inline real fgt_I2L_low(t_fgt * ft, real * I, rvec dx, int m)
{

    real result = *I;
    ++I;

    int i_x;
    int i_y;

    real * mn = ft->multinomials;

    real *I_x;
    real *I_y;
    real *I_z;
    real *I_pre;
    *mn = 1;
    I_x = mn;
    I_y = mn;
    I_z = mn;

    int o;

    for (o = 1; o <= m; ++o)
    {

        I_pre = I_x;
        I_x   = mn + 1;
        for (i_x = o; i_x >= 1; --i_x)
        {
            for (i_y = o - i_x; i_y >= 0; --i_y)
            {

                ++mn;
                *mn     = *I_pre * dx[XX];
                result += *mn * *I;
                ++I;
                ++I_pre;
            }
        }

        I_pre = I_y;
        I_y   = mn + 1;
        for (i_y = o; i_y >= 1; --i_y)
        {
            ++mn;
            *mn     = *I_pre * dx[YY];
            result += *mn * *I;
            ++I;
            ++I_pre;
        }

        I_pre = I_z;
        I_z   = mn + 1;
        ++mn;
        *mn     = *I_pre * dx[ZZ];
        result += *mn * *I;
        ++I;
    }

    return result;
}

void fgt_v2r(t_map * map, int i, rvec r)
{
    ivec iv;
    i2ivec(i, map->map_dim, iv);
    fgt_ivec_add(iv, map->origin, iv);
    ivec_smul(iv, map->spacing, r);
}

void fgt_I2L_enb2v(t_fgt * ft, t_map * map, ivec closest_e, int i_nb, int i,
                   rvec d_r)
{
    ivec e_icur;
    rvec g_r;
    rvec d_curr;
    int e_i;

    fgt_ivec_add(closest_e, ft->e_nbs[i_nb].ndx, e_icur);
    e_i = ivec2i_mod(e_icur, ft->e_n);

    if (ft->I[e_i] != NULL)
    {

        fgt_ivec_mmul(ft->e_cell, ft->e_nbs[i_nb].ndx, g_r);
        rvec_sadd(d_r, -1. / (real) ft->h, g_r, d_curr);

        map->vox[i] += exp(-norm2(d_curr))
            * fgt_I2L_low(ft, ft->I[e_i], d_curr, ft->m);
    }
}
void fgt_copy_origin(t_map * map, ivec origin)
{
    origin[XX] = map->origin[XX];
    origin[YY] = map->origin[YY];
    origin[ZZ] = map->origin[ZZ];
}

int fgt_I2L_e2v(t_fgt * ft, t_map * map, int i)
{

    rvec d_r;
    real norm_d_r;
    rvec v;
    int i_nb;
    ivec closest_e;

    fgt_v2r(map, i, v);

    fgt_r2e(v, ft, closest_e, d_r);
    norm_d_r = fgt_norm(d_r);

    for (i_nb = 0;
         (i_nb < ft->e_N) && (ft->e_nbs[i_nb].d < 3 * ft->sigma + norm_d_r);
         ++i_nb)
    {
        fgt_I2L_enb2v(ft, map, closest_e, i_nb, i, d_r);
    }
    return i_nb;
}

int fgt_I2L_v2e_fill_void(t_fgt * ft, t_map * map, int i, int i_nb)
{
    rvec d_r;
    rvec v;
    ivec closest_e;

    fgt_v2r(map, i, v);

    fgt_r2e(v, ft, closest_e, d_r);

    // look for non-zero voxel contribution
    for (; i_nb < ft->e_N; ++i_nb)
    {
        fgt_I2L_enb2v(ft, map, closest_e, i_nb, i, d_r);
        ++i_nb;
        for (;
             (i_nb < ft->e_N)
             && ft->e_nbs[i_nb - 1].d + 1e-10 < ft->e_nbs[i_nb].d;
             ++i_nb)
        {
            fgt_I2L_enb2v(ft, map, closest_e, i_nb, i, d_r);
        }
        if (map->vox[i] != 0)
        {
            return i_nb;
        }
    }
    return i_nb;
}

int fgt_I2L_at_ref(t_fgt * ft, t_map * map, t_map * ref, real threshold)
{
    if (ft->verbose)
    {
        fprintf(stderr, "Expanding I2L where reference map "
                "is larger than threshold... ");
    }

    int i;
    int i_nb;

    for (i = 0; i < map->map_dim[XX] * map->map_dim[YY] * map->map_dim[ZZ];
         ++i)
    {
        if (ref->vox[i] > threshold)
        {

            i_nb = fgt_I2L_e2v(ft, map, i);

            // fill the void
            if (map->vox[i] == 0 && i_nb < ft->e_N)
            {
                fgt_I2L_v2e_fill_void(ft, map, i, i_nb);
            }
        }

    }

    if (ft->verbose)
    {
        fprintf(stderr, "done.\n");
    }
    return EXIT_SUCCESS;
}

void fgt_map_inv(t_map * map, t_map * ref, real threshold)
{
    int i_map;

    for (i_map = 0;
         i_map < map->map_dim[XX] * map->map_dim[YY] * map->map_dim[ZZ];
         ++i_map)
    {
        if (map->vox[i_map] > 0)
        {
            map->vox[i_map] = 1 / (map->vox[i_map]);
        }
    }
}

void fgt_map_mul(t_map * map, t_map * ref, real threshold)
{
    int i_map;
    for (i_map = 0;
         i_map < map->map_dim[XX] * map->map_dim[YY] * map->map_dim[ZZ];
         ++i_map)
    {
        if (map->vox[i_map] > 0)
        {
            map->vox[i_map] *= ref->vox[i_map];
        }
    }
}

void fgt_map_ratio(t_map * map, t_map * ref, real threshold)
{
    int i_map;
    for (i_map = 0;
         i_map < map->map_dim[XX] * map->map_dim[YY] * map->map_dim[ZZ];
         ++i_map)
    {
        if (map->vox[i_map] > threshold)
        {
            map->vox[i_map] = ref->vox[i_map] / (map->vox[i_map]);
        }
    }
}
/*
 * the integral over the density at expansion centers
 */
real fgt_I2int(real * I, unsigned int max_o)
{

    real result = 0;

    unsigned int i_x;
    unsigned int i_y;
    unsigned int o;

    for (o = 0; o <= max_o; ++o)
    {

        for (i_x = o; i_x >= 1; --i_x)
        {
            for (i_y = o - i_x; i_y != 0; --i_y)
            {
                if (is_even(i_x) && is_even(i_y) && is_even(o - i_x - i_y))
                {
                    result += *I * fgt_gamma_1_2(i_x) * fgt_gamma_1_2(i_y)
                        * fgt_gamma_1_2(o - i_x - i_y);
                    ++I;
                }
            }
        }
        for (i_y = o; i_y >= 1; --i_y)
        {

            if (is_even(i_x) && is_even(i_y) && is_even(o - i_x - i_y))
            {
                result += *I * fgt_gamma_1_2(i_x) * fgt_gamma_1_2(i_y)
                    * fgt_gamma_1_2(o - i_x - i_y);
                ++I;
            }
        }

        if (is_even(i_x) && is_even(i_y) && is_even(o - i_x - i_y))
        {
            result += *I * fgt_gamma_1_2(i_x) * fgt_gamma_1_2(i_y)
                * fgt_gamma_1_2(o - i_x - i_y);
            ++I;
        }
    }

    return result;

}

void fgt_map_link_voxels(t_map * map, real * v)
{
    if (map->own_vox)
    {
        free(map->vox);
    }
    map->vox     = v;
    map->own_vox = (0 == 1);

}
void fgt_map2yf(t_map * map, rvec ** y, real ** f, int * n_div, real threshold)
{

    int i;
    int allocd = 1;
    *n_div = 0;
    *y     = malloc(allocd * sizeof(rvec));
    *f     = malloc(allocd * sizeof(real));
    for (i = 0; i < map->map_dim[XX] * map->map_dim[YY] * map->map_dim[ZZ];
         ++i)
    {

        if (map->vox[i] > threshold)
        {
            ++(*n_div);
            if (*n_div >= allocd)
            {
                allocd *= 2;
                *y      = realloc(*y, allocd * sizeof(rvec));
                *f      = realloc(*f, allocd * sizeof(real));
            }
            fgt_v2r(map, i, (*y)[*n_div - 1]);
            (*f)[*n_div - 1] = map->vox[i];
        }

    }
    *y = realloc(*y, (*n_div + 1) * sizeof(rvec));
    *f = realloc(*f, (*n_div + 1) * sizeof(real));
}

// returns a mapping of expansion center indices to voxel indices
t_voxel_ndx * fgt_map2e(t_map * map, t_fgt * ft)
{
    t_voxel_ndx * result = calloc(ft->e_N, sizeof(t_voxel_ndx));
    int i;
    int i_e;

    for (i_e = 0; i_e < ft->e_N; ++i_e)
    {
        result[i_e].n   = 0;
        result[i_e].val = NULL;
        result[i_e].r   = NULL;
        result[i_e].sum = 0;
    }

    rvec v;
    ivec closest_e;
    rvec d;
    ivec g_i; //< the front bottom left corner of the grid cell

    for (i = 0; i < map->map_dim[XX] * map->map_dim[YY] * map->map_dim[ZZ];
         ++i)
    {
        if (map->vox[i] > 0)
        {

            fgt_v2r(map, i, v);
            fgt_grid_index_ivec(v, ft->e_n, ft->b_inv, g_i);
            fgt_closest_in_cell(g_i, ft->e_cell, v, d, closest_e);
            i_e = ivec2i_mod(closest_e, ft->e_n);

            result[i_e].n++;
            if (is_power_of_two(result[i_e].n))
            {
                result[i_e].val = realloc(result[i_e].val,
                                          sizeof(real) * 2 * result[i_e].n);
                result[i_e].r = realloc(result[i_e].r,
                                        sizeof(rvec) * 2 * result[i_e].n);
            }

            rvec_sadd(ft->e_c[i_e], -1, d, result[i_e].r[result[i_e].n - 1]);
            result[i_e].val[result[i_e].n - 1] = map->vox[i];
            result[i_e].sum                   += map->vox[i];
        }
    }

    return result;
}
int fgt_I2L(t_fgt * ft, t_map * map)
{
    if (ft->verbose)
    {
        fprintf(stderr,
                "Expanding I2L ... \n [\t\t\t\t\t\t\t\t\t\t\t\t\t\t]\r [");
    }

    int i;

    for (i = 0; i < map->map_dim[XX] * map->map_dim[YY] * map->map_dim[ZZ];
         ++i)
    {
        fgt_I2L_e2v(ft, map, i);
        if (ft->verbose
            && map->map_dim[XX] * map->map_dim[YY] * map->map_dim[ZZ]
            > 1000)
        {
            if (i
                % (map->map_dim[XX] * map->map_dim[YY] * map->map_dim[ZZ]
                   / 109) == 0)
            {
                fprintf(stderr, ".");
            }
        }
    }

    if (ft->verbose)
    {
        fprintf(stderr, "\ndone.\n");
    }
    return EXIT_SUCCESS;
}

real fgt_inner(real * a, real * b, int o)
{
    int i;
    real result = 0;
    for (i = 0; i <= o; ++i)
    {
        result += a[i] * b[i];
    }
    return result;
}

void fgt_do_force_lr_fpoly(ivec m, rvec f, real I, real * powers[], real ** c,
                           real ** cg)
{

    f[XX] += I * fgt_inner(cg[m[XX]], powers[XX], m[XX] + 1)
        * fgt_inner(c[m[YY]], powers[YY], m[YY])
        * fgt_inner(c[m[ZZ]], powers[ZZ], m[ZZ]);
    f[YY] += I * fgt_inner(c[m[XX]], powers[XX], m[XX])
        * fgt_inner(cg[m[YY]], powers[YY], m[YY] + 1)
        * fgt_inner(c[m[ZZ]], powers[ZZ], m[ZZ]);
    f[ZZ] += I * fgt_inner(c[m[XX]], powers[XX], m[XX])
        * fgt_inner(c[m[YY]], powers[YY], m[YY])
        * fgt_inner(cg[m[ZZ]], powers[ZZ], m[ZZ] + 1);

}

void fgt_do_force_lr_low(rvec f, real * I, int max_o, rvec d, real ** c,
                         real ** cg)
{
    int o;
    ivec m;
    fgt_clear_rvec(f);
    real * I_curr = I;

    real **powers = malloc(sizeof(real*) * 3);

    // force will have one more polynomial order than expansion polynomials
    powers[XX] = malloc(sizeof(real) * (max_o + 2));
    powers[YY] = malloc(sizeof(real) * (max_o + 2));
    powers[ZZ] = malloc(sizeof(real) * (max_o + 2));

    powers[XX][0] = 1;
    powers[YY][0] = 1;
    powers[ZZ][0] = 1;

    for (o = 1; o <= max_o + 1; ++o)
    {
        powers[XX][o] = powers[XX][o - 1] * (-d[XX]);
        powers[YY][o] = powers[YY][o - 1] * (-d[YY]);
        powers[ZZ][o] = powers[ZZ][o - 1] * (-d[ZZ]);
    }

    for (o = 0; o <= max_o; ++o)
    {
        for (m[XX] = o; m[XX] >= 1; --m[XX])
        {
            for (m[YY] = o - m[XX]; m[YY] >= 0; --m[YY])
            {
                m[ZZ] = o - m[XX] - m[YY];
                fgt_do_force_lr_fpoly(m, f, *I_curr, powers, c, cg);
                I_curr++;
            }
        }
        for (m[YY] = o; m[YY] >= 1; --m[YY])
        {
            m[ZZ] = o - m[XX] - m[YY];
            fgt_do_force_lr_fpoly(m, f, *I_curr, powers, c, cg);
            I_curr++;
        }
        m[ZZ] = o;
        fgt_do_force_lr_fpoly(m, f, *I_curr, powers, c, cg);
        I_curr++;
    }

    free(powers[XX]);
    free(powers[YY]);
    free(powers[ZZ]);
    free(powers);

}

void chimera_vector(FILE* outfile, rvec a, rvec b)
{
    fprintf(outfile, ".arrow %f %f %f %f %f %f\n", 10 * a[XX], 10 * a[YY],
            10 * a[ZZ], 10 * b[XX], 10 * b[YY], 10 * b[ZZ]);
}

void fgt_do_force_lr(rvec f, ivec closest_e, rvec d, t_fgt * div_ft,
                     t_voxel_ndx * v_ndx, real delta, real s, real ** c, real ** cg, real k)
{

    int i_nb;
    ivec e_icur;
    rvec c_r;
    rvec d_curr;
    int e_i;
    rvec f_low;
    fgt_clear_rvec(f);

    for (i_nb = 0; i_nb < div_ft->e_N; ++i_nb)
    {

        fgt_ivec_add(closest_e, div_ft->e_nbs[i_nb].ndx, e_icur);
        e_i = ivec2i_mod_const(e_icur, div_ft->e_n);

        if (div_ft->I[e_i] != NULL)
        {

            fgt_ivec_mmul(div_ft->e_cell, div_ft->e_nbs[i_nb].ndx, c_r);
            rvec_diff(d, c_r, d_curr);

            fgt_do_force_lr_low(f_low, div_ft->I[e_i], div_ft->m, d_curr, c,
                                cg);

            rvec_splus(
                    -exp(-norm2(d_curr) / (s * s)) * v_ndx[e_i].sum
                    / fgt_I2int(div_ft->I[e_i], div_ft->m), f_low, f);

        }
    }

    rvec_smul(f, k*delta * delta * delta / (s * s * s));

}

void release_vndx(t_voxel_ndx * vndx, int n)
{
    int i;
    for (i = 0; i < n; ++i)
    {
        if (vndx[i].r != NULL)
        {
            free(vndx[i].r);
        }
        if (vndx[i].val != NULL)
        {
            free(vndx[i].val);
        }
    }
    free(vndx);
}

void fgt_lr_poly_add_next(real ** arr, int alpha, real s, real h, int len)
{
    int i;
    arr[alpha][0] = 0.5
        * (h * arr[alpha - 1][1] + (alpha - 1) * arr[alpha - 2][0]);
    for (i = 1; i < len; ++i)
    {
        arr[alpha][i] = 0.5
            * ((alpha - 1) * arr[alpha - 2][i]
               - h * arr[alpha - 1][i - 1] / (s * s)
               + (i + 1) * h * arr[alpha - 1][i + 1]);
    }
    arr[alpha][len]     = 0;
    arr[alpha][len + 1] = -h * arr[alpha - 1][len] / (2 * s * s);

}

void fgt_lr_poly_init_tables(real s, real h, int m, real *** c, real *** cg)
{

    *cg = malloc((m + 1) * sizeof(real *));
    *c  = malloc((m + 1) * sizeof(real*));

    int alpha;

    alpha        = 0;
    (*cg)[alpha] = malloc((alpha + 2) * sizeof(real));
    (*c)[alpha]  = malloc((alpha + 1) * sizeof(real));

    (*cg)[alpha][0] = 0;
    (*cg)[alpha][1] = -sqrt(2) / (2 * s * s);

    (*c)[alpha][0] = 1 / sqrt(2);

    ++alpha;

    (*cg)[alpha] = malloc((alpha + 2) * sizeof(real));
    (*c)[alpha]  = malloc((alpha + 1) * sizeof(real));

    (*cg)[alpha][0] = -h * sqrt(2) / (4.0 * s * s);
    (*cg)[alpha][1] = 0;
    (*cg)[alpha][2] = h * sqrt(2) / (4.0 * s * s * s * s);

    (*c)[alpha][0] = 0;
    (*c)[alpha][1] = -sqrt(2) * h / (4.0 * s * s);

    while (alpha < m)
    {
        ++alpha;

        (*cg)[alpha] = malloc((alpha + 2) * sizeof(real));
        (*c)[alpha]  = malloc((alpha + 1) * sizeof(real));

        fgt_lr_poly_add_next(*cg, alpha, s, h, alpha);
        fgt_lr_poly_add_next(*c, alpha, s, h, alpha - 1);

    }

}

void fgt_lr_poly_release_tables(real ** c, real ** cg, int m)
{
    int alpha;
    for (alpha = 0; alpha < m + 1; ++alpha)
    {
        free(cg[alpha]);
        free(c[alpha]);
    }
    free(cg);
    free(c);
}

void fgt_do_force(t_fgt * ft, t_map * ref, real * sim, rvec ** f, real k)
{

    real threshold = 1e-10;
    t_map * map;
    int i_atoms;
    rvec d;
    rvec * y;
    real * w;
    t_voxel_ndx * v_ndx;
    int n_div      = 0;
    t_fgt * div_ft = NULL;

    init_map(&map, ref->map_dim, ref->origin, ref->spacing, 0 == 0);

    fgt_S2I(ft);
    fgt_I2L_at_ref(ft, map, ref, threshold);

    fgt_map_ratio(map, ref, threshold);

    if (*f == NULL)
    {
        *f = calloc(ft->n, sizeof(rvec));
    }
    for (i_atoms = 0; i_atoms < ft->n; ++i_atoms)
    {
        fgt_clear_rvec((*f)[i_atoms]);
    }

    fgt_map2yf(map, &y, &w, &n_div, threshold);

    if (init_fgt(&div_ft, y, w, n_div, ft->box, ft->sigma, 5, 1 == 0, 0 == 0,
                 0 == 0) != EXIT_SUCCESS)
    {
        fprintf(stderr, "ERROR: FAILED to initialise IFGT.\n");
        return;
    }

    fgt_S2I(div_ft);

//    fgt_I2L(div_ft, map);
//    fgt_copy_vox(map, sim);
    v_ndx = fgt_map2e(map, div_ft);

    ivec closest_e;
    ivec g_i; //< the front bottom left corner of the grid cell
    real s = sqrt(div_ft->sigma * div_ft->sigma + ft->sigma * ft->sigma);
    real ** c;
    real ** cg;

    fgt_lr_poly_init_tables(s, div_ft->h, div_ft->m, &c, &cg);

//    FILE * outfile = fopen("forcevec.bild", "w");
//    rvec f_to;

    for (i_atoms = 0; i_atoms < ft->n; ++i_atoms)
    {

        fgt_grid_index_ivec(ft->y[i_atoms], div_ft->e_n, div_ft->b_inv, g_i);
        fgt_closest_in_cell(g_i, div_ft->e_cell, ft->y[i_atoms], d, closest_e);
        fgt_do_force_lr((*f)[i_atoms], closest_e, d, div_ft, v_ndx,
                        map->spacing, s, c, cg, k);
//        fgt_rvec_add(ft->y[i_atoms], (*f)[i_atoms], f_to);
//        chimera_vector(outfile, ft->y[i_atoms], f_to);
    }

    fgt_lr_poly_release_tables(c, cg, div_ft->m);
    release_vndx(v_ndx, div_ft->e_N);
//    fclose(outfile);
    release_fgt(div_ft);
    release_map(map);

}

void fgt_copy_vox(t_map * map, real * vox)
{

    int i;

    for (i = 0; i < map->map_dim[XX] * map->map_dim[YY] * map->map_dim[ZZ];
         ++i)
    {

        vox[i] = map->vox[i];

    }

}

int do_fgt(t_fgt * ft, t_map * map)
{

    if (fgt_S2I(ft) != EXIT_FAILURE)
    {
        return fgt_I2L(ft, map);
    }
    return EXIT_FAILURE;

}

int release_fgt(t_fgt * ft)
{
    int i;
    for (i = 0; i < ft->e_N; ++i)
    {
        if (ft->I[i] != NULL)
        {
            free(ft->I[i]);
        }
    }
    free(ft->e_nbs);
    free(ft->I);

    free(ft->e_c);

    if (ft->f_own)
    {
        free(ft->f);
    }
    if (ft->y_own)
    {
        free(ft->y);
    }

    free(ft->multinomials);

    free(ft);

    return EXIT_SUCCESS;
}

int release_map(t_map * map)
{

    if (map->own_vox)
    {
        free(map->vox);
    }
    free(map);
    return EXIT_SUCCESS;

}

int init_map(t_map ** map, ivec nvoxels, ivec origin, real spacing, int own_vox)
{

    *map                = malloc(sizeof(t_map));
    (*map)->map_dim[XX] = nvoxels[XX];
    (*map)->map_dim[YY] = nvoxels[YY];
    (*map)->map_dim[ZZ] = nvoxels[ZZ];

    if (origin != NULL)
    {
        (*map)->origin[XX] = origin[XX];
        (*map)->origin[YY] = origin[YY];
        (*map)->origin[ZZ] = origin[ZZ];
    }
    else
    {
        (*map)->origin[XX] = 0;
        (*map)->origin[YY] = 0;
        (*map)->origin[ZZ] = 0;
    }

    (*map)->spacing = spacing;
    (*map)->own_vox = own_vox;
    if ((*map)->own_vox)
    {
        (*map)->vox = calloc(
                    (*map)->map_dim[XX] * (*map)->map_dim[YY] * (*map)->map_dim[ZZ],
                    sizeof(real));
        return (*map)->vox == NULL ? EXIT_FAILURE : EXIT_SUCCESS;
    }
    else
    {
        return EXIT_SUCCESS;
    }

}
