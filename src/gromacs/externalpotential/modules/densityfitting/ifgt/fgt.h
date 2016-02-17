#ifndef _FGT_h
#define _FGT_h

#if __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <omp.h>
#include "fgt_math.h"
#include "fgt_types.h"

int init_fgt(t_fgt ** ft, rvec * y, real * f, int natom, matrix box, real sigma,
             int m, int verbose, int y_own, int f_own);

int do_fgt(t_fgt * ft, t_map * map);
void fgt_map_link_voxels(t_map * map, real * v);
int fgt_S2I(t_fgt * ft);
void fgt_do_force(t_fgt * ft, t_map * ref, real * sim, rvec ** f, real k);

int error_map(const t_map * map, t_map * ref);
int release_map(t_map * map);
int release_fgt(t_fgt * ft);
void fgt_copy_vox(t_map * map, real * vox);
void fgt_copy_section(t_map * map, ivec map_dim, real * vox);

void fgt_normalize(t_map * map);

real max_map(t_map * map);
int init_map(t_map **map, ivec nvoxels, ivec origin, real spacing, int own_vox);
void fgt_copy_origin(t_map * map, ivec origin);


#if __cplusplus
}   // Extern C
#endif

#endif /* FGT_h */
