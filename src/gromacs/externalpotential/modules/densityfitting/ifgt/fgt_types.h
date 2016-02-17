/*
 * fgt_types.h
 *
 *  Created on: Aug 4, 2015
 *      Author: cblau
 */

#ifndef FGT_TYPES_H_
#define FGT_TYPES_H_


#define DIMENSIONS (3)
#define PI (3.14159265358979232)

#ifndef XX
#define XX 0
#endif

#ifndef YY
#define YY 1
#endif

#ifndef ZZ
#define ZZ 2
#endif

#define NO_NEIGHBOUR -1
#define MAX_REAL 1e20
#define fgt_sqr(a) ((a)*(a))

typedef float real;
typedef real rvec[DIMENSIONS];
typedef int ivec[DIMENSIONS];
typedef real matrix[DIMENSIONS][DIMENSIONS];
typedef struct complex t_cmpx;

typedef struct fgt t_fgt;
typedef struct map t_map;
typedef struct voxel_ndx t_voxel_ndx;


struct nb {
    real d;
    ivec ndx;
};

typedef struct nb t_nb;

struct fgt {
    int    n;       /* Number of sources to be expanded                   */
    matrix box;     /* (periodic) bounding box for the sources            */
    matrix b_inv;   /* the inverse of the box */
    real   sigma;   /* Width of the Gaussian for spreading                */
    real   n_sigma; /* Maximum expansion box size (sqrt(ft->m)) */
    rvec  *y;       /* Coordinates of the sources to be expanded          */
    int    y_own;   /* True if our own memory has been allocated for the coordinates*/
    real  *f;       /* Weight for spreading (scattering cross section)    */
    int    f_own;   /* True if our own memory has been allocated for the weights */
    /* determines following parameters                    */

    int    m;            /* Maximum order of Taylor expansion                  */

    real **I;            /* The improved Gauss transform multinomial coefficents */
    int    n_I;          /* The number of multinomial expansion coefficients */

    real   h;            /* sqrt(2)*sigma (compliant with delta in Greengard)  */
    ivec   e_n;          /* Number of expansion boxes in x,y,z - direction     */
    rvec   e_l;          /* Length of expansion boxes in x,y,z - direction     */
    int    e_N;          /* Total number of expansion boxes                    */

    rvec * e_c;          /* Center of expansion boxes                          */
    matrix e_cell;       /* box[XX]/e_n[XX], box[YY]/e_n[YY] etc */
    int    verbose;      /* */
    real * multinomials; /* memory for multinomial vector calculation */
    t_nb * e_nbs;        /* neighbouring e-boxes according to grid */
};

struct voxel_ndx {
    int    n;
    real * val;
    rvec * r;
    real   sum;
};

struct map {
    real *vox;     /* Value at the voxels                                */
    ivec  map_dim; /* Number of voxels in x,y,z - direction              */
    ivec  origin;  /* map origin at spacing*origin*/
    real  spacing;
    int   own_vox;
};


#endif /* FGT_TYPES_H_ */
