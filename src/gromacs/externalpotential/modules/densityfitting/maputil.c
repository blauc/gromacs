/*
 * gmx_map.c
 *
 *  Created on: Jul 25, 2011
 *      Author: ckutzne
 */
#include <stdio.h>
#include <sys/time.h>
#include "gromacs/utility/smalloc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/filenm.h"
#include "gromacs/math/vec.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/sim_util.h"
#include "gromacs/legacyheaders/types/inputrec.h"
#include "gromacs/domdec/domdec.h"
#include "assert.h"
#include "gromacs/legacyheaders/mdrun.h"
#include "gromacs/mdlib/groupcoord.h"
#include "gromacs/timing/cyclecounter.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/legacyheaders/gmx_omp_nthreads.h"
#include "gromacs/random/random.h"

/* Environment variable for setting nstfit */
static const char* NSTFIT_ENVVAR = "GMX_NSTFIT";

//#define TIME_DENSFIT

//#define GMX_DENSFIT_MKL

#ifdef GMX_DENSFIT_MKL
#include "mkl_vml.h"
#endif

#ifndef _maputil_h
#include "maputil.h"
#endif

#define WITHCCP4
#ifdef WITHCCP4
#include "cmaplib.h"
#else
typedef struct _CMMFile CMMFile;
#define FLOAT32 0;
#endif

static char* FitStrMDRUN = "Density fitting: "; /* Used if density fitting output is written from mdrun */
static char* FitStr      = "";                  /* Used in g_map output */

/* Maputil.c local density fitting data */
typedef struct gmx_densfit {
    FILE       *out_dens;      /* Output file for density fitting data           */
    const char *fn_map;        /* Output filename for simulated density maps     */
    real        spacing;       /* Map spacing in nm, derived from map data       */
    gmx_bool    bAppend;       /* Are we appending to existing output files?     */
    gmx_bool    bVerbose;      /* -v flag from command line                      */
    real        cc;            /* Correlation coefficient of the two maps        */
    rvec       *x_assembled;   /* Just the positions that are spread             */
    rvec       *f_loc;         /* Forces on the local atoms of the density fitting group*/
    real       *w_assembled;   /* Weights for the assembled positions            */
    rvec       *x_old;         /* x_assembled from last step                     */
    ivec       *x_shifts;      /* To make the molecule whole                     */
    ivec       *extra_shifts;  /* Shifts added since last NS                     */
    atom_id    *c_ind;         /* Where is this atom stored in the coll array?   */
    gmx_bool    bUpdateShifts; /* Do we have to calculate new shift vectors
                                * from x_assembled and x_old?                    */
    real       *tmpgrid;       /* Stores the erf values around a single atom.
                                * This way, we can compute a whole bunch of erf's
                                * at once, which is more efficient               */
    real  *tmpgrid2;           /* same for exp                                   */
    int    voxrange;           /* Max. number of voxels to be computed for a
                                  single atom in a single dimension x, y, or z   */
    int    vox_per_atom;       /* Total number of voxels for one atom (x,y, + z) */
    double sum_rho_exp2;       /* Term needed for the calculation of the
                                  correlation coefficient and the forces         */
} t_gmx_densfit;

/* Allocate and initialize the arrays for assembled positions, forces and
 * atomic weights. */
static t_gmx_densfit *new_t_gmx_densfit(int nat,                                      /* Size of ind and x_densfit_whole    */
                                        real sigma, real dist, rvec *x_densfit_whole, /* Can be NULL!                       */
                                        real grid_spacing, gmx_bool bVerbose, gmx_bool bAppend, gmx_bool bParallel)
{
    t_gmx_densfit *df = NULL;
    int            i, ompthreads;

    snew(df, 1);

    df->bVerbose = bVerbose;
    df->spacing  = grid_spacing;

    snew(df->x_assembled, nat);
    snew(df->f_loc, nat);
    snew(df->w_assembled, nat);
    snew(df->x_old, nat);
    snew(df->x_shifts, nat);
    snew(df->extra_shifts, nat);
    snew(df->c_ind, nat);

    if (df->bVerbose)
    {
        fprintf(stdout, "%sSetting all %d atomic weights to unity.\n", FitStr,
                nat);
    }
    for (i = 0; i < nat; i++)
    {
        df->w_assembled[i] = 1.0;
    }

    df->bUpdateShifts = TRUE;
    df->bAppend       = bAppend;

    if (!bParallel)
    {
        for (i = 0; i < nat; i++)
        {
            df->c_ind[i] = i;
        }
    }

    /* In mdruns we have to make the density fitting structure whole */
    if (NULL != x_densfit_whole)
    {
        for (i = 0; i < nat; i++)
        {
            copy_rvec(x_densfit_whole[i], df->x_old[i]);
        }
    }

    /* For erf() and exp() performance, we want all values for one atom in one
     * linear array (not in three). Therefore, we compute the maximum number of
     * voxels to be computed for a single atom. How many voxels around each atom
     * we have to spread the density on is saved in voxrange (for each dimension
     * x, y, and z */
    df->voxrange     = ceil(0.5 * sigma * dist / df->spacing);
    df->vox_per_atom = 3 * (2 * df->voxrange + 1);
    ompthreads       = max(1, gmx_omp_nthreads_get(emntDefault));
    snew(df->tmpgrid, df->vox_per_atom * ompthreads);
    snew(df->tmpgrid2, df->vox_per_atom * ompthreads);

    return df;
}

/* Make a structure containing density fitting data */
extern t_densfit *new_t_densfit(real sigma, real sigma_dist, real k, /* Spring constant */
                                int nat, atom_id *ind,               /* Which of the atoms should be used for spreading */
                                real grid_spacing, gmx_bool bVerbose)
{
    t_densfit *densfit = NULL;
    int        i;

    snew(densfit, 1); /* Global density fitting data from inputrec.h     */

    densfit->sigma = sigma;
    densfit->dist  = sigma_dist;
    densfit->k     = k;
    densfit->nat   = nat;

    if (NULL == ind)
    {
        if (bVerbose)
        {
            fprintf(stdout,
                    "%sYou did not provide an index group - will use all atoms for spreading.\n",
                    FitStr);
        }

        snew(densfit->ind, nat);
        for (i = 0; i < nat; i++)
        {
            densfit->ind[i] = i;
        }
    }
    else
    {
        densfit->ind = ind;
    }

    densfit->map_sim = NULL;
    densfit->map_ref = NULL;

    /* Initialize also the private data structure: */
    densfit->df = new_t_gmx_densfit(nat, densfit->sigma, densfit->dist, NULL,
                                    grid_spacing, bVerbose, 0, 0);

    return densfit;
}

static double gettime()
{
#ifdef HAVE_GETTIMEOFDAY
    struct timeval t;
    double         seconds;

    gettimeofday(&t, NULL);

    seconds = (double) t.tv_sec + 1e-6 * (double) t.tv_usec;

    return seconds;
#else
    double seconds;

    seconds = time(NULL);

    return seconds;
#endif
}

void dump_x(t_densfit *densfit, t_commrec *cr, gmx_int64_t step)
{
    FILE          *fpout = NULL;
    char           fn[STRLEN];
    char           buf[256];
    int            i;
    t_gmx_densfit *df;

    df = densfit->df;
    gmx_step_str(step, buf);
    sprintf(fn, "x_step%s_node%d.txt", buf, cr->nodeid);
    fpout = fopen(fn, "w");

    for (i = 0; i < densfit->nat; i++)
    {
        fprintf(fpout, "%4d %12.3e %12.3e %12.3e\n", i, df->x_assembled[i][XX],
                df->x_assembled[i][YY], df->x_assembled[i][ZZ]);
    }

    fclose(fpout);
}

extern void dump_f(const char *fn, t_densfit *densfit, t_commrec *cr)
{
    int            i;
    t_gmx_densfit *df;
    rvec          *f;
    gmx_bool       bWrite;

    FILE          *fpout = NULL;

    df = densfit->df;

    /* Assemble the forces */
    snew(f, densfit->nat);
    /* Each node copies its local forces into the collective array */
    for (i = 0; i < densfit->nat_loc; i++)
    {
        copy_rvec(df->f_loc[i], f[df->c_ind[i]]);
    }

    /* Reduce */
    if (cr && PAR(cr))
    {
        gmx_sum(3 * densfit->nat, f[0], cr);
    }

    /* Master outputs */
    bWrite = FALSE;
    if (cr)
    {
        if (MASTER(cr))
        {
            bWrite = TRUE;
        }
    }
    else
    {
        bWrite = TRUE;
    }

    if (bWrite)
    {
        fpout = gmx_fio_fopen(fn, "w");

        fprintf(stderr, "%sDumping forces to file %s\n", FitStr, fn);

        for (i = 0; i < densfit->nat; i++)
        {
            fprintf(fpout, "%4d %12.3e %12.3e %12.3e\n", densfit->ind[i] + 1,
                    f[i][XX], f[i][YY], f[i][ZZ]);
        }

        gmx_fio_fclose(fpout);
    }
    sfree(f);
}

/* Print the header information stored in file on disk */
extern void gmx_map_print_header(FILE *fpout, const char *fn)
{
    int          i;
    CMMFile     *mfile;
    char        *title;
    ivec         grid;
    float        cell[6];
    int          origin[3];
    int          axes_order[3];
    int          map_dim[3];
    int          spacegroup;
    float        min, max;
    double       mean, rms;
    unsigned int datamode;
    float        skew_mat[9], skew_trans[3];

#ifdef WITHCCP4
    mfile = ccp4_cmap_open(fn, O_RDONLY);

    /* Get the header information */
    ccp4_cmap_get_cell(mfile, cell);
    ccp4_cmap_get_grid(mfile, grid);
    ccp4_cmap_get_origin(mfile, origin);
    ccp4_cmap_get_order(mfile, axes_order);
    ccp4_cmap_get_dim(mfile, map_dim);
    ccp4_cmap_get_mapstats(mfile, &min, &max, &mean, &rms);
    ccp4_cmap_get_mask(mfile, skew_mat, skew_trans);
    title      = ccp4_cmap_get_title(mfile);
    spacegroup = ccp4_cmap_get_spacegroup(mfile);
    datamode   = ccp4_cmap_get_datamode(mfile);

    /* Print the header information: */
    fprintf(fpout, "\n");
    fprintf(fpout, "CCP4 Header:\n");
    fprintf(fpout, "Title      : %s\n", title);
    fprintf(fpout, "Datamode   : %u\n", datamode);
    fprintf(fpout, "Cell       : %g A  %g A  %g A  %g %g %g\n", cell[0],
            cell[1], cell[2], cell[3], cell[4], cell[5]);
    fprintf(fpout, "Grid       : %d x %d x %d\n", grid[XX], grid[YY], grid[ZZ]);
    fprintf(fpout, "Origin     : %d A  %d A  %d A\n", origin[XX], origin[YY],
            origin[ZZ]);
    fprintf(fpout, "Axes order : %d %d %d\n", axes_order[0], axes_order[1],
            axes_order[2]);
    fprintf(fpout, "Map_dim    : %d %d %d\n", map_dim[0], map_dim[1],
            map_dim[2]);
    fprintf(fpout, "Spacegroup : %d\n", spacegroup);
    fprintf(fpout, "Min        : %g\n", min);
    fprintf(fpout, "Max        : %g\n", max);
    fprintf(fpout, "Mean       : %g\n", mean);
    fprintf(fpout, "RMS        : %g\n", rms);
    fprintf(fpout, "Skew_trans : %g %g %g\n", skew_trans[0], skew_trans[1],
            skew_trans[2]);
    fprintf(fpout, "Skew_mat   : ");
    for (i = 0; i < 9; i++)
    {
        fprintf(fpout, "%g ", skew_mat[i]);
    }
    fprintf(fpout, "\n\n");

    ccp4_cmap_close(mfile);
#else
    gmx_fatal(FARGS, "Not compiled with CCP4 support!");
#endif
}

/* Allocate memory for the density data */
extern void allocate_density_grid(int nValues[3], float **vox)
{
    int ompthreads = max(1, gmx_omp_nthreads_get(emntDefault));

    snew(*vox, ompthreads * nValues[XX] * nValues[YY] * nValues[ZZ]);
}

/* Allocate memory for the simulated density */
extern void new_map_sim_from_ref(t_densfit *densfit, t_mapdata *map_ref)
{
    int i;

    snew(densfit->map_sim, 1);

    densfit->map_sim->datamode = FLOAT32;

    /* The simulated map copies the map_dim and grid settings from the reference map. */
    for (i = 0; i < 3; i++)
    {
        densfit->map_sim->map_dim[i]    = map_ref->map_dim[i];
        densfit->map_sim->grid[i]       = map_ref->grid[i];
        densfit->map_sim->axes_order[i] = map_ref->axes_order[i];
        densfit->map_sim->origin[i]     = map_ref->origin[i];
    }
    for (i = 0; i < 6; i++)
    {
        densfit->map_sim->cell[i] = map_ref->cell[i];
    }
    densfit->map_sim->title = strdup("Calculated map");

    allocate_density_grid(map_ref->map_dim, &densfit->map_sim->vox);
}

/* If the box turns out to be all 0's, then assign a useful box from the positions */
extern void translate_pdb(matrix box, rvec x[], int natoms, gmx_bool bTranslate)
{
    int  i;
    real x_max = -GMX_REAL_MAX, y_max = -GMX_REAL_MAX, z_max = -GMX_REAL_MAX;
    real x_min = GMX_REAL_MAX, y_min = GMX_REAL_MAX, z_min = GMX_REAL_MAX;

    if ((box[XX][YY] != 0.0) || (box[XX][ZZ] != 0.0) || (box[YY][XX] != 0.0)
        || (box[YY][ZZ] != 0.0) || (box[ZZ][XX] != 0.0)
        || (box[ZZ][YY] != 0.0))
    {
        gmx_fatal(FARGS, "All off-diagonal box elements must be zero.\n");
    }

    for (i = 0; i < natoms; i++)
    {
        x_min = MIN(x_min, x[i][XX]);
        y_min = MIN(y_min, x[i][YY]);
        z_min = MIN(z_min, x[i][ZZ]);

        x_max = MAX(x_max, x[i][XX]);
        y_max = MAX(y_max, x[i][YY]);
        z_max = MAX(z_max, x[i][ZZ]);
    }
    fprintf(stdout, "Positions were in range %g <= x <= %g\n", x_min, x_max);
    fprintf(stdout, "Positions were in range %g <= y <= %g\n", y_min, y_max);
    fprintf(stdout, "Positions were in range %g <= z <= %g\n", z_min, z_max);

    /* Optionally translate all positions to first quadrant */
    if (bTranslate)
    {
        for (i = 0; i < natoms; i++)
        {
            x[i][XX] -= x_min;
            x[i][YY] -= y_min;
            x[i][ZZ] -= z_min;
        }

        box[XX][XX] = x_max - x_min;
        box[YY][YY] = y_max - y_min;
        box[ZZ][ZZ] = z_max - z_min;
    }
    else
    {
        box[XX][XX] = x_max;
        box[YY][YY] = y_max;
        box[ZZ][ZZ] = z_max;
    }

    fprintf(stdout, "Assigning rectangular box xmax=%f ymax=%f zmax=%f nm\n",
            box[XX][XX], box[YY][YY], box[ZZ][ZZ]);
}

/* If bRead==TRUE : Write map data to fn
 *           FALSE: Read  map data from fn and store in map, allocate memory
 */
extern void gmx_do_map_ccp4(gmx_bool bRead, t_mapdata **ptr_map, const char *fn,
                            gmx_bool bOverwrite, gmx_bool bVerbose, FILE *DumpHeader)
{
#ifdef WITHCCP4
    CMMFile   *fp = NULL;
    int        n;
    t_mapdata *map;
#endif

    if (NULL == fn)
    {
        gmx_fatal(FARGS,
                  "%sgmx_do_map_ccp4: no filename given, cannot %s map!\n",
                  FitStr, bRead ? "read" : "write");
    }

    if (bRead && !gmx_fexist(fn))
    {
        gmx_fatal(FARGS,
                  "Could not find ccp4 input file '%s' with electron density map\n",
                  fn);
    }

    if (bVerbose)
    {
        fprintf(stderr, "%sOpening ccp4 cmap file '%s' for %s.\n", FitStr, fn,
                bRead ? "reading" : "writing");
    }

#ifdef WITHCCP4
    /* If mapdata is not already allocated, do it now */
    if (NULL == *ptr_map)
    {
        snew(*ptr_map, 1);
    }
    map = *ptr_map;

    if (!bRead && gmx_fexist(fn))
    {
        if (!bOverwrite)
        {
            make_backup(fn);
        }
    }

    fp = ccp4_cmap_open(fn, bRead ? O_RDONLY : O_WRONLY);

    if (FLOAT32 != ccp4_cmap_get_datamode(fp))
    {
        gmx_fatal(FARGS,
                  "The only supported ccp4 file data mode is FLOAT32.\n");
    }

    /* First read/write the header information */
    if (bRead)
    {
        map->title      = strdup(ccp4_cmap_get_title(fp));
        map->datamode   = ccp4_cmap_get_datamode(fp);
        map->spacegroup = ccp4_cmap_get_spacegroup(fp);
    }
    else
    {
        ccp4_cmap_set_title(fp, map->title);
        ccp4_cmap_set_datamode(fp, map->datamode);
        ccp4_cmap_set_spacegroup(fp, map->spacegroup);
    }

    if (bRead)
    {
        ccp4_cmap_get_cell(fp, map->cell);
        ccp4_cmap_get_grid(fp, map->grid);
        ccp4_cmap_get_origin(fp, map->origin);
        ccp4_cmap_get_order(fp, map->axes_order);
        ccp4_cmap_get_dim(fp, map->map_dim);
        ccp4_cmap_get_mask(fp, map->skew_mat, map->skew_trans);
    }
    else
    {
        ccp4_cmap_set_cell(fp, map->cell);
        ccp4_cmap_set_grid(fp, map->grid);
        ccp4_cmap_set_origin(fp, map->origin);
        ccp4_cmap_set_order(fp, map->axes_order);
        ccp4_cmap_set_dim(fp, map->map_dim);
        ccp4_cmap_set_title(fp, map->title);
        ccp4_cmap_set_mask(fp, map->skew_mat, map->skew_trans);
    }

    if (bRead)
    {
        ccp4_cmap_get_mapstats(fp, &map->min, &map->max, &map->mean, &map->rms);
    }
    else
    {
        ccp4_cmap_set_mapstats(fp, map->min, map->max, map->mean, map->rms);
    }

    /* Now read/write the actual density data */
    n = map->map_dim[XX] * map->map_dim[YY] * map->map_dim[ZZ];

    if (bRead)
    {
        allocate_density_grid(map->map_dim, &map->vox);
    }

    /* write n_items from memory to file (item determined by data mode) */
    if (bRead)
    {
        ccp4_cmap_read_data(fp, map->vox, n);
    }
    else
    {
        ccp4_cmap_write_data(fp, map->vox, n);
    }

    ccp4_cmap_close(fp);
#else
    fprintf(stderr, "Cannot write CCP4 map: Not compiled with CCP4 support!");
#endif

    if (DumpHeader)
    {
        gmx_map_print_header(DumpHeader, fn);
    }
}

extern void couple_map_spacing_to_box(matrix box, t_densfit *densfit)
{
    int        i;
    rvec       spacing;

    t_mapdata *map = densfit->map_ref;

    for (i = 0; i < DIM; i++)
    {
        map->cell[i] = NM2A * box[i][i];
        spacing[i]   = A2NM * map->cell[i] / map->map_dim[i];
    }

    /* Save the average grid spacing */
    densfit->df->spacing = (spacing[XX] + spacing[YY] + spacing[ZZ]) / 3.0;
}

/* From the map entries determine the grid spacing */
extern real get_map_spacing(t_mapdata *map, FILE *log)
{
    rvec spacing;
    real av_spacing;
    ivec points;
    int  i;

    if (map == NULL)
    {
        gmx_fatal(FARGS, "Cannot determine grid spacing without a map.\n");
    }

    /* Number of grid points in each dimension */
    for (i = 0; i < DIM; i++)
    {
        points[i]  = map->map_dim[i];
        spacing[i] = A2NM * map->cell[i] / points[i];
    }

    if ((fabs(spacing[XX] - spacing[YY]) > 0.01)
        || (fabs(spacing[XX] - spacing[ZZ]) > 0.01))
    {
        fprintf(stderr,
                "WARING: Reference grid spacing (%g %g %g) should be identical in XYZ !",
                spacing[XX], spacing[YY], spacing[ZZ]);
    }

    av_spacing = (spacing[XX] + spacing[YY] + spacing[ZZ]) / 3.0;

    if (log)
    {
        fprintf(log,
                "%sGrid spacing (nm) is %g in x-, %g in y-, and %g in z-dimension (average %g).\n",
                FitStr, spacing[XX], spacing[YY], spacing[ZZ], av_spacing);
    }

    /* Return the grid spacing */
    return av_spacing;
}

/* Get the length of the intersection of the two segments a and b */
static inline double get_intersecting_length(double x11, double x12, double x21,
                                             double x22)
{
    return max(0, min(x12, x22) - max(x11, x21));
}

/* Get the volume of the intersection of new_vox and ref_vox */
double get_intersecting_volume(dvec v1_start, dvec v1_stop, dvec v2_start,
                               dvec v2_stop)
{
    int  d;
    dvec intersection;

    for (d = 0; d < DIM; d++)
    {
        intersection[d] = get_intersecting_length(v1_start[d], v1_stop[d],
                                                  v2_start[d], v2_stop[d]);
    }

    return intersection[XX] * intersection[YY] * intersection[ZZ];
}

/* Return the average (integrated) density of the new, scaled voxel at ix, iy, iz which
 * is made of several smaller voxels of the original density map
 *
 * Variable naming convention: i for map_ref indices,
 *                             j for map_scaled indices
 */
static double get_rho_av(t_mapdata *map_ref, t_mapdata *map_new, ivec j)
{
    int    d, ix, iy, iz, ind;
    ivec   n_ref, n_new;
    dvec   mapsize;
    double rho;                               /* Average / integrated density in Angstrom^(-3) */
    double vol;
    dvec   ref_vox_size, new_vox_size;        /* Angstrom */
    dvec   new_vox_start, new_vox_stop;       /* Angstrom */
    dvec   ref_vox_start, ref_vox_stop;       /* Angstrom */
    ivec   ref_vox_startind, ref_vox_stopind; /* First and last ref voxel index per dimension
                                                 that contributes to new voxel density */
    float  refdensity;

    for (d = 0; d < DIM; d++)
    {
        /* Number of voxels per dimension of the reference and scaled map */
        n_ref[d] = map_ref->map_dim[d];
        n_new[d] = map_new->map_dim[d];

        /* Size (Angstrom) of the maps (both have the same size!) */
        mapsize[d] = map_ref->cell[d];

        /* Voxel size per dimension of reference and scaled map */
        ref_vox_size[d] = mapsize[d] / n_ref[d];
        new_vox_size[d] = mapsize[d] / n_new[d];

        /* Start and stop position (Angstrom) of the new voxel */
        new_vox_start[d] = new_vox_size[d] * j[d];
        new_vox_stop[d]  = new_vox_size[d] * (j[d] + 1);

        /* Get the first and the last voxel index of the ref map which
         * contributes to the voxel density at this point of the scaled map */
        ref_vox_startind[d] = floor(new_vox_start[d] / ref_vox_size[d]);
        ref_vox_stopind[d]  = floor(new_vox_stop[d] / ref_vox_size[d]);

        ref_vox_stopind[d] = min(ref_vox_stopind[d], n_ref[d] - 1);
    }

    /* Now loop over all small subvoxels, get the subvoxel volume plus density value,
     * and add to resulting voxel density */
    rho = 0.0;
    for (iz = ref_vox_startind[ZZ]; iz <= ref_vox_stopind[ZZ]; iz++)
    {
        for (iy = ref_vox_startind[YY]; iy <= ref_vox_stopind[YY]; iy++)
        {
            for (ix = ref_vox_startind[XX]; ix <= ref_vox_stopind[XX]; ix++)
            {
                ind        = n_ref[XX] * n_ref[YY] * iz + n_ref[XX] * iy + ix;
                refdensity = map_ref->vox[ind];

                if (refdensity != 0)
                {
                    /* Start and stop position (Angstrom) of the reference voxel with which to intersect */
                    ref_vox_start[XX] = ref_vox_size[XX] * ix;
                    ref_vox_stop[XX]  = ref_vox_size[XX] * (ix + 1);
                    ref_vox_start[YY] = ref_vox_size[YY] * iy;
                    ref_vox_stop[YY]  = ref_vox_size[YY] * (iy + 1);
                    ref_vox_start[ZZ] = ref_vox_size[ZZ] * iz;
                    ref_vox_stop[ZZ]  = ref_vox_size[ZZ] * (iz + 1);

                    /* Get the intersecting volume of the ref voxel at ix,iy,iz with
                     * the new voxel */
                    vol = get_intersecting_volume(new_vox_start, new_vox_stop,
                                                  ref_vox_start, ref_vox_stop);
                    rho += refdensity * vol;
                }
            }
        }
    }

    return rho;
}

extern void make_positive(t_mapdata *map_ref)
{
    float minimum = GMX_FLOAT_MAX;
    int   nx, ny, nz, i, count;

    nx = map_ref->map_dim[XX];
    ny = map_ref->map_dim[YY];
    nz = map_ref->map_dim[ZZ];

    /* Determine the minimum value of the voxel density */
    count = 0;
    for (i = 0; i < nx * ny * nz; i++)
    {
        minimum = min(minimum, map_ref->vox[count]);
        count++;
    }

    /* Subtract the minimum */
    fprintf(stdout, "Subtracting the minimum vale of %g from all voxels.\n",
            minimum);
    count = 0;
    for (i = 0; i < nx * ny * nz; i++)
    {
        map_ref->vox[count] -= minimum;
        count++;
    }

    map_ref->min  = 0.0;
    map_ref->max  = 0.0;
    map_ref->mean = 0.0;
    map_ref->rms  = 0.0;
}

//TODO: scale sigma by the same factor?
//TODO: allow only values > 0?
extern t_mapdata *rescale_map(t_mapdata *map_ref, real scale)
{
    int        i, ix, iy, iz;
    ivec       ind;
    t_mapdata *map_scaled = NULL;
    int        nx, ny, nz;
    int        mx, my, mz;
    float     *scaledDens = NULL;
    double     rho_av;

    snew(map_scaled, 1);
    map_scaled->datamode = FLOAT32;

    for (i = 0; i < 3; i++)
    {
        map_scaled->map_dim[i] = round(map_ref->map_dim[i] * scale);
        map_scaled->grid[i]    = round(map_ref->grid[i] * scale);
    }
    for (i = 0; i < 6; i++)
    {
        map_scaled->cell[i] = map_ref->cell[i];
    }
    map_scaled->title = strdup("Scaled map");

    allocate_density_grid(map_scaled->map_dim, &map_scaled->vox);

    nx = map_scaled->map_dim[XX];
    ny = map_scaled->map_dim[YY];
    nz = map_scaled->map_dim[ZZ];

    mx = map_ref->map_dim[XX];
    my = map_ref->map_dim[YY];
    mz = map_ref->map_dim[ZZ];

    double volratio = (double) (nx * ny * nz) / (double) (mx * my * mz);
    fprintf(stderr, "volratio = %g\n", volratio);

    scaledDens = map_scaled->vox;
    /* Determine the average density for the new voxel per dimension */
    for (iz = 0; iz < nz; iz++)
    {
        ind[ZZ] = iz;
        for (iy = 0; iy < ny; iy++)
        {
            ind[YY] = iy;
            for (ix = 0; ix < nx; ix++)
            {
                /* Now determine the average density (in Angstrom^(-3)) of the original map at the ix,iy,iz voxel */
                ind[XX] = ix;
                rho_av  = get_rho_av(map_ref, map_scaled, ind);

                /* Assign average density of original map to new voxel */
                scaledDens[nx * ny * iz + nx * iy + ix] = rho_av * volratio;
            }
        }
    }

    return map_scaled;
}

/* Like gmx_erf, but for n elements. This can be accelerated a lot when using
 * the vector math MKL library Erf functions (> factor 7) */
static void gmx_erf_vector(int n, real x[])
{
    int i;

#ifdef TIME_DENSFIT
    gmx_cycles_t c1, c2;

    c1 = gmx_cycles_read();
#endif

#ifdef GMX_DENSFIT_MKL
    /* Call the MKL vector math error function */
#ifdef GMX_DOUBLE
    vdErf(n, x, x);
#else
    vsErf(n, x, x);
#endif
#else
    /* Or the default if MKL is not present */
    for (i = 0; i < n; i++)
    {
        x[i] = gmx_erf(x[i]);
//        x[i] = erf(x[i]);
    }
#endif

#ifdef TIME_DENSFIT
    c2 = gmx_cycles_read();
    fprintf(stderr, "gmx_erf_vector %15d elements %15lld cycles  %15.3f per erf\n", n, c2-c1, (1.0*(c2-c1))/n);
#endif
}

/* As gmx_exp, using n elements. This can be accelerated by using
 * the vector math MKL library functions. */
static inline void gmx_exp_vector(int n, real x[])
{
    int i;

#ifdef TIME_DENSFIT
    gmx_cycles_t c1, c2;

    c1 = gmx_cycles_read();
#endif

#ifdef GMX_DENSFIT_MKL
#ifdef GMX_DOUBLE
    vdExp(n, x, x);
#else
    vsExp(n, x, x);
#endif
#else
    for (i = 0; i < n; i++)
    {
#ifdef GMX_DOUBLE
        x[i] = exp(x[i]);
#else
        x[i] = expf(x[i]);
#endif
    }
#endif

#ifdef TIME_DENSFIT
    c2 = gmx_cycles_read();
    fprintf(stderr, "gmx_exp_vector %15d elements %15lld cycles  %15.3f per exp\n", n, c2-c1, (1.0*(c2-c1))/n);
#endif

}

/*
 * vectorise x*exp(-OOsigma2*x*x); use mkl if possible
 */
static inline void gmx_rexp_vector(int n, real r[], real OOsigma2)
{
    int i;

#ifdef TIME_DENSFIT
    gmx_cycles_t c1, c2;

    c1 = gmx_cycles_read();
#endif

#ifdef GMX_DENSFIT_MKL
    real * tmp;
    snew(tmp, n);
#ifdef GMX_DOUBLE
    vdMul(n, r, r, tmp);
    vdMul(n, -OOsigma2, tmp, tmp);
    vdExp(n, tmp, tmp);
    vdMul(n, r, tmp, r);
#else
    vsMul(n, r, r, tmp);
    vsMul(n, OOsigma2, tmp, tmp);
    vsExp(n, tmp, tmp);
    vsMul(n, r, tmp, r);
#endif
    sfree(tmp);
#else
    for (i = 0; i < n; i++)
    {
#ifdef GMX_DOUBLE
        r[i] = r[i]*exp(-OOsigma2*r[i]*r[i]);
#else
        r[i] = r[i] * expf(-OOsigma2 * r[i] * r[i]);
#endif
    }
#endif

#ifdef TIME_DENSFIT
    c2 = gmx_cycles_read();
    fprintf(stderr, "gmx_xexp_vector %15d elements %15lld cycles  %15.3f per exp\n", n, c2-c1, (1.0*(c2-c1))/n);
#endif

}

/* For a position vector x return all direct 27 image positions given the box
 * TODO: check range when voxrange is large (segfault when istart leaves map)
 * */
static int get_images(rvec x, matrix box, rvec img[], ivec img_startindex[],
                      ivec img_stop_index[], real grid_spacing, int voxrange, ivec map_dim)
{
    int    i, j, k, m, count, icenter;
    rvec   image, shiftvecX, shiftvecY, shiftvecZ;
    ivec   istart, istop;

    double OOgrid_spacing = 1.0 / grid_spacing;

    /* Loop over all possible 3x3x3 neighboring boxes */
    /* Minimum image convention? TODO: loop over NTRICIMG = 14, not 27 */
    count = 0;
    for (i = -1; i <= 1; i++)
    {
        for (j = -1; j <= 1; j++)
        {
            for (k = -1; k <= 1; k++)
            {
                svmul(k, box[XX], shiftvecX);
                svmul(j, box[YY], shiftvecY);
                svmul(i, box[ZZ], shiftvecZ);

                copy_rvec(x, image);

                /* Move in x, y, and z directions */
                rvec_inc(image, shiftvecX);
                rvec_inc(image, shiftvecY);
                rvec_inc(image, shiftvecZ);

                /* Now check whether we will get any density from that image */
                icenter    = roundf(image[XX] * OOgrid_spacing);
                istart[XX] = icenter - voxrange;
                istop[XX]  = icenter + voxrange;

                if ((istart[XX] <= map_dim[XX]) && (istop[XX] >= 0))
                {

                    icenter    = roundf(image[YY] * OOgrid_spacing);
                    istart[YY] = icenter - voxrange;
                    istop[YY]  = icenter + voxrange;

                    if ((istart[YY] <= map_dim[YY]) && (istop[YY] >= 0))
                    {

                        icenter    = roundf(image[ZZ] * OOgrid_spacing);
                        istart[ZZ] = icenter - voxrange;
                        istop[ZZ]  = icenter + voxrange;

                        if ((istart[ZZ] <= map_dim[ZZ]) && (istop[ZZ] >= 0))
                        {
                            /* This is our new image */
                            copy_rvec(image, img[count]);

                            /* Cut away the volume where we do not have a grid */
                            for (m = 0; m < DIM; m++)
                            {
                                istart[m] = max(istart[m], 0);         /* istart >= 0       */
                                istop[m]  = min(istop[m], map_dim[m]); /* istop  <= map_dim */
                            }
                            copy_ivec(istart, img_startindex[count]);
                            copy_ivec(istop, img_stop_index[count]);
                            count++;

                        } /* ZZ */
                    }     /* YY */
                }         /* XX */
            }
        }
    }

    return count;
}

/* Make the density map, assuming that the protein is whole, and in the first
 * quadrant. Transform the positions x into a density map by replacing
 * each position by a Gaussian function of width sigma */
static void spread_atoms_low(rvec       x[],    /* Atomic positions                          */
                             int        natoms, /* Number of atomic positions                */
                             matrix     box,    /* The simulation box                        */
                             t_densfit *densfit)
{
    int            iz, iy, ix, i, j, n, th;
    int            nx, ny, nz, ngrid;
    double         V_vox;
    double         prefactor;
    double         left;
    real           ax, ay, az, ayaz;
    float         *sim_dens_1d; /* 1d array for 3d simulated density map          */
    float         *ptr;         /* Pointer in that array                          */
    ivec           elements;
    rvec           pos;
    real          *aarr[3];
    int            image, num_images;
    rvec           img[27];    /* Stores PBC images of a position x as well as x */
    ivec           istart[27]; /* Store for each image the part of the grid ...  */
    ivec           istop[27];  /* ... where stuff has to be calculated           */

    double         sigma2 = densfit->sigma * densfit->sigma;
    double         c_fac  = sqrt(1.5 / sigma2);
    t_gmx_densfit *df     = densfit->df;
    int            nth    = max(1, gmx_omp_nthreads_get(emntDefault));

//    gmx_cycles_t c1 = gmx_cycles_read();

    sim_dens_1d = densfit->map_sim->vox; /* Calculated density map */
    V_vox       = pow(df->spacing, 3);
    prefactor   = pow(M_PI * sigma2 / 6.0, 1.5);
    prefactor  /= V_vox;
    nx          = densfit->map_sim->map_dim[XX];
    ny          = densfit->map_sim->map_dim[YY];
    nz          = densfit->map_sim->map_dim[ZZ];
    ngrid       = nx * ny * nz;

    /* Zero out the grid */
#pragma omp parallel for num_threads(nth) schedule(static) \
    shared(sim_dens_1d, ngrid, nth) \
    default(none)
    for (i = 0; i < ngrid * nth; i++)
    {
        sim_dens_1d[i] = 0.0;
    }

    /* Spread the local atoms on the grid, therefore loop over all atoms of the
     * density fitting group */
#pragma omp parallel for num_threads(nth) schedule(static) \
    shared(stderr,nx, ny, nz, ngrid, natoms, x, box, df, densfit, c_fac, sim_dens_1d) \
    private(ptr, aarr, istart, istop, image, pos, num_images, img, i, ix, iy, iz, th, elements, ax, ay, az, ayaz, left) \
    default(none)
    for (n = 0; n < natoms; n++)
    {
        th = gmx_omp_get_thread_num(); /* Thread ID */

        /* Determine the grid points around the atom and its images for which we need to calculate stuff */
        num_images = get_images(x[n], box, img, istart, istop, df->spacing,
                                df->voxrange, densfit->map_sim->map_dim);

        /* Loop over all periodic images of the atom */
        for (image = 0; image < num_images; image++)
        {
            copy_rvec(img[image], pos); // periodic image of the position

            elements[XX] = istop[image][XX] + 1 - istart[image][XX];
            elements[YY] = istop[image][YY] + 1 - istart[image][YY];
            elements[ZZ] = istop[image][ZZ] + 1 - istart[image][ZZ];

            int sum = elements[XX] + elements[YY] + elements[ZZ];

            /* For erf() performance, we want all values for one atom in one linear array (not
             * int three). However, we build pointers to access it as x, y, and z arrays. */
            aarr[XX] = &df->tmpgrid[0 + th * df->vox_per_atom];
            aarr[YY] = &df->tmpgrid[0 + elements[XX] + th * df->vox_per_atom];
            aarr[ZZ] = &df->tmpgrid[0 + elements[XX] + elements[YY]
                                    + th * df->vox_per_atom];

            /* Store the x,y,z-arguments for the erf(x,y,z) evaluations at the voxel
             * boundaries in a temporary array. This way we can compute the erf's
             * for this atom all in one (vector) call. */
            i = 0;
            for (ix = istart[image][XX]; ix < istop[image][XX] + 1; ix++)
            {
                left          = (ix - 0.5) * df->spacing - pos[XX];
                aarr[XX][i++] = c_fac * left;
            }

            i = 0;
            for (iy = istart[image][YY]; iy < istop[image][YY] + 1; iy++)
            {
                left          = (iy - 0.5) * df->spacing - pos[YY];
                aarr[YY][i++] = c_fac * left;
            }

            i = 0;
            for (iz = istart[image][ZZ]; iz < istop[image][ZZ] + 1; iz++)
            {
                left          = (iz - 0.5) * df->spacing - pos[ZZ];
                aarr[ZZ][i++] = c_fac * left;
            }
            gmx_erf_vector(elements[XX] + elements[YY] + elements[ZZ],
                           aarr[XX]);

            /* Transform (in place) into a new array that directly stores the differences of
             * the error functions, i.e. the area under the Gaussian curve */
            for (ix = 0; ix < elements[XX] - 1; ix++)
            {
                aarr[XX][ix] = aarr[XX][ix + 1] - aarr[XX][ix];
            }
            for (iy = 0; iy < elements[YY] - 1; iy++)
            {
                aarr[YY][iy] = aarr[YY][iy + 1] - aarr[YY][iy];
            }
            for (iz = 0; iz < elements[ZZ] - 1; iz++)
            {
                aarr[ZZ][iz] = aarr[ZZ][iz + 1] - aarr[ZZ][iz];
            }
            /* Spread the density for this atom */
            for (iz = istart[image][ZZ]; iz < istop[image][ZZ]; iz++)
            {
                az = aarr[ZZ][iz - istart[image][ZZ]];

                for (iy = istart[image][YY]; iy < istop[image][YY]; iy++)
                {
                    ay   = aarr[YY][iy - istart[image][YY]];
                    ayaz = ay * az;
                    ptr  = &sim_dens_1d[th * ngrid + nx * (ny * iz + iy)
                                        + istart[image][XX]];
                    for (ix = istart[image][XX]; ix < istop[image][XX]; ix++)
                    {
                        ax = aarr[XX][ix - istart[image][XX]];

                        /* Assign value to voxel */
//                        sim_dens_1d[nz*ny*ix + nz*iy + iz] += ax*ay*az;
                        *ptr += ax * ayaz;

                        ptr++;
                    }
                }
            }
        } /* End of loop over images of one atom */
    }     /* End of loop over atoms */

#pragma omp parallel for num_threads(nth) schedule(static) \
    shared(nth, sim_dens_1d, ngrid, prefactor) \
    private(j) \
    default(none)
    for (i = 0; i < ngrid; i++)
    {
        for (j = 1; j < nth; j++)
        {
            sim_dens_1d[i] += sim_dens_1d[i + j * ngrid];
        }
        sim_dens_1d[i] *= prefactor;
    }

//    fprintf(stderr, "--- Spread_atoms_low: %g M cycles.\n", 0.000001*(gmx_cycles_read() - c1) );
}

/*
 * Improved fast Gauss transform
 * */

#ifdef WITHFGT
extern void spread_atoms_low_fgt(rvec       x[],    /* Atomic positions                          */
                                 int        natoms, /* Number of atomic positions                */
                                 matrix     box,    /* The simulation box                        */
                                 t_densfit *densfit)
{
    gmx_cycles_t c1 = gmx_cycles_read();
    real       * sim_dens_1d;

    t_fgt      * ft;
    t_map      * map;

    init_fgt(&ft, x, densfit->weight, natoms, box, densfit->sigma, 5, TRUE,
             FALSE, FALSE);
    init_map(&map, densfit->map_sim->map_dim, NULL, densfit->df->spacing, TRUE);

    do_fgt(ft, map);
    fgt_normalize(map);

    release_fgt(ft);
    fgt_copy_vox(map, densfit->map_sim->vox);
    fgt_copy_origin(map, densfit->map_sim->origin);
    release_map(map);

    //fprintf(stderr, "--- Spread_atoms_low_fgt: %g M cycles.\n", 0.000001*(gmx_cycles_read() - c1) );
}
#endif

/* Transform the positions x into a density map by replacing
 * each position by a Gaussian function of width sigma */
extern void spread_atoms(t_densfit *densfit, matrix box)
{
    t_gmx_densfit *df;     /* Contains maputil.c local data             */
    double         t0 = 0;
    long int       c1 = 0; /* make compiler happy */

    int           *map_dim = densfit->map_sim->map_dim;

    df = densfit->df;

    if (df->bVerbose)
    {
        fprintf(stderr,
                "%sSpreading %d atomic positions  on a %d x %d x %d grid ...\n",
                FitStr, densfit->nat, map_dim[XX], map_dim[YY], map_dim[ZZ]);
        t0 = gettime();
        c1 = gmx_cycles_read();
    }

    if (!df->spacing > 0)
    {
        gmx_fatal(FARGS, "%sSpacing must not be zero!", FitStr);
    }

    /* Make the density map */
#ifdef WITHFGT
    spread_atoms_low_fgt(df->x_assembled, densfit->nat, box, densfit);
#else
    spread_atoms_low(df->x_assembled, densfit->nat, box, densfit);
#endif

    if (df->bVerbose)
    {
        fprintf(stderr, "%sSpreading took %g seconds (%g M cycles).\n", FitStr,
                gettime() - t0, 0.000001 * (gmx_cycles_read() - c1));
    }
}

/* Multiply each voxel of the map with the corresponding voxel of the other
 * map and sum up everything, i.e. compute
 *
 *  sum_ijk [ rhoA_ijk * rhoB_ijk ]
 *
 */
static double calc_sum_rhoA_rhoB(t_mapdata *mapA, t_mapdata *mapB)
{
    int    i, count;
    int    nx, ny, nz;
    double sum_ijk;

    nx = mapA->map_dim[XX];
    ny = mapA->map_dim[YY];
    nz = mapA->map_dim[ZZ];

    sum_ijk = 0.0;
    count   = 0;
    for (i = 0; i < nx * ny * nz; i++)
    {
        sum_ijk += mapA->vox[count] * mapB->vox[count];
        count++;
    }

    return sum_ijk;
}

/* Calculate the correlation coefficient of the two maps */
extern real calc_correlation_coeff(t_densfit *densfit, FILE *log)
{
    int        i;
    double     numerator, denominator;
    real       cc;
    t_mapdata *map_exp; /* Reference map */
    t_mapdata *map_sim; /* Calculated map */

    map_exp = densfit->map_ref;
    map_sim = densfit->map_sim;

    if (log)
    {
        fprintf(log,
                "%sCalculating the correlation coefficient of the two maps.\n",
                FitStr);
    }

    for (i = 0; i < 3; i++)
    {
        if (map_exp->map_dim[i] != map_sim->map_dim[i])
        {
            gmx_fatal(FARGS, "Map dimensions must agree");
        }
    }

    numerator = calc_sum_rhoA_rhoB(map_exp, map_sim);

    denominator = sqrt(
                calc_sum_rhoA_rhoB(map_exp, map_exp)
                * calc_sum_rhoA_rhoB(map_sim, map_sim));

    cc = numerator / denominator;

    if (log)
    {
        fprintf(log, "%sThe correlation coefficient is %15.10f\n", FitStr, cc);
    }

    return cc;
}

extern void do_densfit_forces(t_densfit *densfit, matrix box)
{
    int    ix, iy, iz, i, l, ind;
    int    nx, ny, nz;
    double sum_rho_exp2;
    double sum_rho_sim2;
    double sum_rho_exp_rho_sim;
    double term1_prefactor, term2_prefactor;
    dvec   term1_sum, term2_sum;
    dvec   tmp1, tmp2;
    dvec   drho_dxl;
    dvec   force;
    rvec   pos;
    int    th;    /* ID of the thread in this team */
    int    nth;   /* Number of threads */
    ivec   elements;
    double V_vox; //TODO: with pressure coupling, the grid spacing needs to be adjusted!!!
    double prefactor;
    double left;
    double OOsigma2_fac;     /* 3 / (2*sigma*sigma) */
    double sqrtOOsigma2_fac; /* sqrt of the above */
    real   erfx, erfy, erfz; /* erf( sqrt( 3 / (2*sigma*sigma) * x)  for x */
    real   expx, expy, expz; /* exp( -3 x*x / (2*sigma*sigma) */
    real   erfy_erfz, expy_erfz, erfy_expz;

    real  *aarr[3], *garr[3];
    float *sim_dens_1d, *ref_dens_1d;
    float *ptr_ref, *ptr_sim;
    double ref_dens, sim_dens;

    int    image, num_images;
    rvec   img[27];                      /* Stores PBC images of a position x as well as x */
    ivec   istart[27];                   /* Store for each image the part of the grid ...  */
    ivec   istop[27];                    /* ... where stuff has to be calculated           */

    sim_dens_1d = densfit->map_sim->vox; /* Calculated density map */
    ref_dens_1d = densfit->map_ref->vox; /* Reference density map */
    t_gmx_densfit *df = densfit->df;

    if (densfit->nat_loc <= 0)
    {
        return;
    }

    nth = max(1, gmx_omp_nthreads_get(emntDefault));

    /* Zero out the forces array TODO */
#pragma omp parallel for num_threads(nth) schedule(static)
    for (l = 0; l < densfit->nat_loc; l++)
    {
        clear_rvec(df->f_loc[l]);
    }

    /* Calculate various prefactors */
    sum_rho_exp2        = densfit->df->sum_rho_exp2;
    sum_rho_sim2        = calc_sum_rhoA_rhoB(densfit->map_sim, densfit->map_sim);
    sum_rho_exp_rho_sim = calc_sum_rhoA_rhoB(densfit->map_ref,
                                             densfit->map_sim);

    term1_prefactor = densfit->k / sqrt(sum_rho_exp2 * sum_rho_sim2);
    term2_prefactor = -densfit->k * sum_rho_exp_rho_sim
        / (sqrt(sum_rho_exp2) * pow(sum_rho_sim2, 1.5));

    V_vox            = df->spacing * df->spacing * df->spacing;
    prefactor        = -M_PI * densfit->sigma * densfit->sigma / (6.0 * V_vox);
    OOsigma2_fac     = 1.5 / (densfit->sigma * densfit->sigma);
    sqrtOOsigma2_fac = sqrtf(OOsigma2_fac);

    nx = densfit->map_sim->map_dim[XX];
    ny = densfit->map_sim->map_dim[YY];
    nz = densfit->map_sim->map_dim[ZZ];

#pragma omp parallel for num_threads(nth) schedule(static) \
    private(term1_sum, term2_sum, pos, istart, istop, elements, aarr, garr, i, \
    ix, iy, iz, left, ind, erfx, erfy, erfz, expx, expy, expz, erfy_erfz, expy_erfz, erfy_expz, \
    ptr_ref, ptr_sim, ref_dens, sim_dens, drho_dxl, force, tmp1, tmp2, th, num_images, image, img) \
    shared(stderr,densfit, df, ref_dens_1d, sim_dens_1d, box, nx, ny, nz, \
    OOsigma2_fac, sqrtOOsigma2_fac, prefactor, term1_prefactor, term2_prefactor) \
    default(none)
    /* Loop over all atoms of the density fitting group */
    for (l = 0; l < densfit->nat_loc; l++)
    {
        th = gmx_omp_get_thread_num();

        /* Loop over all voxels and calculate term1 and term2 for the force
         * simultaneously */
        clear_dvec(term1_sum);
        clear_dvec(term2_sum);

        /* Determine the grid points around the atom and its images which
         * have a contribution to the force on this atom.
         * Get the position of this atom of the density fitting group, we must
         * use the positions from the assembled array, because only these have
         * been put in the box in do_densfit() */
        num_images = get_images(df->x_assembled[df->c_ind[l]], box, img, istart,
                                istop, df->spacing, df->voxrange, densfit->map_sim->map_dim);

        /* Loop over all periodic images of the atom */
        for (image = 0; image < num_images; image++)
        {
            copy_rvec(img[image], pos); // periodic image of the position

            /* The calculation of drho_dxl is expensive. We only calculate it, where
             * there is actually simulated density, i.e. in the given sigma-range applied
             * during spreading. */
            elements[XX] = istop[image][XX] + 1 - istart[image][XX];
            elements[YY] = istop[image][YY] + 1 - istart[image][YY];
            elements[ZZ] = istop[image][ZZ] + 1 - istart[image][ZZ];

            /* For erf() and exp() performance, we want all values for one atom in one linear array (not
             * in three). However, we build pointers to access it as x, y, and z arrays. */
            aarr[XX] = &df->tmpgrid[0 + th * df->vox_per_atom];
            aarr[YY] = &df->tmpgrid[0 + elements[XX] + th * df->vox_per_atom];
            aarr[ZZ] = &df->tmpgrid[0 + elements[XX] + elements[YY]
                                    + th * df->vox_per_atom];

            garr[XX] = &df->tmpgrid2[0 + th * df->vox_per_atom];
            garr[YY] = &df->tmpgrid2[0 + elements[XX] + th * df->vox_per_atom];
            garr[ZZ] = &df->tmpgrid2[0 + elements[XX] + elements[YY]
                                     + th * df->vox_per_atom];

            /* Store the x,y,z-arguments for the erf(x,y,z) evaluations at the voxel
             * boundaries in the tmpgrid[] array. This way we can compute the erf's
             * for this atom all in one (vector) call. */
            i = 0;
            for (ix = istart[image][XX]; ix < istop[image][XX] + 1; ix++)
            {
                left          = (ix - 0.5) * df->spacing - pos[XX];
                garr[XX][i]   = -OOsigma2_fac * left * left;
                aarr[XX][i++] = sqrtOOsigma2_fac * left;
            }

            i = 0;
            for (iy = istart[image][YY]; iy < istop[image][YY] + 1; iy++)
            {
                left          = (iy - 0.5) * df->spacing - pos[YY];
                garr[YY][i]   = -OOsigma2_fac * left * left;
                aarr[YY][i++] = sqrtOOsigma2_fac * left;
            }

            i = 0;
            for (iz = istart[image][ZZ]; iz < istop[image][ZZ] + 1; iz++)
            {
                left          = (iz - 0.5) * df->spacing - pos[ZZ];
                garr[ZZ][i]   = -OOsigma2_fac * left * left;
                aarr[ZZ][i++] = sqrtOOsigma2_fac * left;
            }

            /* Call erf() and exp() for all input values for this atom in one go */
            gmx_erf_vector(elements[XX] + elements[YY] + elements[ZZ],
                           aarr[XX]);
            gmx_exp_vector(elements[XX] + elements[YY] + elements[ZZ],
                           garr[XX]);

            /* Transform (in place) into a new array that directly stores the differences of
             * the error functions, i.e. the area under the Gaussian curve */
            for (ix = 0; ix < elements[XX] - 1; ix++)
            {
                aarr[XX][ix] = aarr[XX][ix + 1] - aarr[XX][ix];
                garr[XX][ix] = garr[XX][ix + 1] - garr[XX][ix];
            }
            for (iy = 0; iy < elements[YY] - 1; iy++)
            {
                aarr[YY][iy] = aarr[YY][iy + 1] - aarr[YY][iy];
                garr[YY][iy] = garr[YY][iy + 1] - garr[YY][iy];
            }
            for (iz = 0; iz < elements[ZZ] - 1; iz++)
            {
                aarr[ZZ][iz] = aarr[ZZ][iz + 1] - aarr[ZZ][iz];
                garr[ZZ][iz] = garr[ZZ][iz + 1] - garr[ZZ][iz];
            }

            for (iz = 0; iz < elements[ZZ] - 1; iz++)
            {
                erfz = aarr[ZZ][iz];
                expz = garr[ZZ][iz];

                for (iy = 0; iy < elements[YY] - 1; iy++)
                {
                    erfy = aarr[YY][iy];
                    expy = garr[YY][iy];

                    erfy_erfz = erfy * erfz;
                    expy_erfz = expy * erfz;
                    erfy_expz = erfy * expz;

                    ind = ny * nx * (iz + istart[image][ZZ])
                        + nx * (iy + istart[image][YY]) + istart[image][XX];
                    ptr_ref = &ref_dens_1d[ind];
                    ptr_sim = &sim_dens_1d[ind];

                    /* Calculate d/dx_l rho_sim_ijk */
                    for (ix = 0; ix < elements[XX] - 1; ix++)
                    {
                        ref_dens = (double) *ptr_ref;
                        sim_dens = (double) *ptr_sim;

                        erfx = aarr[XX][ix];
                        expx = garr[XX][ix];

                        drho_dxl[XX] = expx * erfy_erfz;
                        drho_dxl[YY] = erfx * expy_erfz;
                        drho_dxl[ZZ] = erfx * erfy_expz;

                        term1_sum[XX] += ref_dens * drho_dxl[XX];
                        term1_sum[YY] += ref_dens * drho_dxl[YY];
                        term1_sum[ZZ] += ref_dens * drho_dxl[ZZ];

                        term2_sum[XX] += sim_dens * drho_dxl[XX];
                        term2_sum[YY] += sim_dens * drho_dxl[YY];
                        term2_sum[ZZ] += sim_dens * drho_dxl[ZZ];

                        ptr_ref++;
                        ptr_sim++;
                    }
                }
            } /* end of loop over voxels of this image */
        }     /* end of loop over images */

        dsvmul(prefactor * term1_prefactor, term1_sum, tmp1);
        dsvmul(prefactor * term2_prefactor, term2_sum, tmp2);

        dvec_add(tmp1, tmp2, force);
        fprintf(stderr, "\n\n%g \n", force[XX]);
        df->f_loc[l][XX] = force[XX];
        df->f_loc[l][YY] = force[YY];
        df->f_loc[l][ZZ] = force[ZZ];

    } /* end of loop over atoms */

}

#ifdef WITHFGT
extern void do_densfit_forces_fgt(t_densfit *densfit, matrix box)
{

    t_fgt * ft;
    t_map * ref;

    if (init_fgt(&ft, densfit->df->x_assembled, densfit->weight, densfit->nat, box,
                 densfit->sigma, 5, FALSE, FALSE, FALSE) != 0)
    {
        fprintf(stderr, "ERROR: Failed IFGT initialisation.");
        return;
    }

    init_map(&ref, densfit->map_ref->map_dim, densfit->map_ref->origin, densfit->df->spacing, FALSE);
    fgt_map_link_voxels(ref, densfit->map_ref->vox);

    fgt_do_force(ft, ref, densfit->map_sim->vox, &(densfit->df->f_loc), densfit->k);

    release_fgt(ft);
    release_map(ref);

}
#endif

/* Adds 'buf' to 'str' */
static void add_to_string(char **str, char *buf)
{
    int len;

    len = strlen(*str) + strlen(buf) + 1;
    srenew(*str, len);
    strcat(*str, buf);
}

static void add_to_string_aligned(char **str, char *buf)
{
    char buf_aligned[STRLEN];

    sprintf(buf_aligned, "%12s", buf);

    add_to_string(str, buf_aligned);
}

/* Open output file and print some general information about density fitting.
 * Call on master only */
static FILE *open_densfit_out(const char *fn, t_densfit *densfit,
                              const output_env_t oenv)
{
    FILE        *fp;
    int          nsets;
    const char **setname;
    char         buf[50];
    char        *LegendStr = NULL;

    if (densfit->df->bAppend)
    {
        fp = gmx_fio_fopen(fn, "a");
    }
    else
    {
        fp = xvgropen(fn, "Density fitting correlation coefficient",
                      "Time (ps)", "correlation coefficient", oenv);
        fprintf(fp,
                "# Density fitting forces for %d atoms are updated every %d time step%s.\n",
                densfit->nat, densfit->nstfit, densfit->nstfit != 1 ? "s" : "");
        fprintf(fp,
                "#   Note that the force update frequency can be overwritten by the environment variable %s.\n",
                NSTFIT_ENVVAR);
        fprintf(fp, "# Output is written in intervals of %d time step%s.\n",
                densfit->nstout, densfit->nstout > 1 ? "s" : "");
        fprintf(fp,
                "# Fitting parameters: sigma = %g nm, sigma_dist = %g, k = %g kJ/mol\n",
                densfit->sigma, densfit->dist, densfit->k);
        if (0 != densfit->nstmapout)
        {
            fprintf(fp, "# Calculated map will be saved every %d steps\n",
                    densfit->nstmapout);
        }

        /* Print a nice legend */
        snew(LegendStr, 1);
        LegendStr[0] = '\0';
        sprintf(buf, "#     %6s", "time");
        add_to_string_aligned(&LegendStr, buf);

        nsets = 0;
        snew(setname, 2);
        setname[nsets++] = strdup("correlation coefficient");
        add_to_string_aligned(&LegendStr, "c.c.");
        setname[nsets++] = strdup("V_dens (kJ/mol)");
        add_to_string_aligned(&LegendStr, "V_dens");
        xvgr_legend(fp, nsets, setname, oenv);
        sfree(setname);

        fprintf(fp, "#\n# Legend for the following data columns:\n");
        fprintf(fp, "%s\n", LegendStr);
        sfree(LegendStr);

        fflush(fp);
    }

    return fp;

}

/* Append a number to the output file name */
extern void make_filename(const char *outf_name, int ndigit, int file_nr,
                          char newname[STRLEN])
{
    char  fmt[128], nbuf[128];
    int   nd = 0, fnr;
    char *outf_base;
    char *extpos;

    /* Strip away the extension */
    outf_base = strdup(outf_name);
    extpos    = strrchr(outf_base, '.');
    if (extpos)
    {
        *extpos = '\0';
    }

    fnr = file_nr;
    do
    {
        fnr /= 10;
        nd++;
    }
    while (fnr > 0);

    if (ndigit > 0)
    {
        snprintf(fmt, 128, "%%0%dd", ndigit);
        snprintf(nbuf, 128, fmt, file_nr);
    }
    else
    {
        snprintf(nbuf, 128, "%d", file_nr);
    }

    if (extpos)
    {
        snprintf(newname, STRLEN, "%s%s.%s", outf_base, nbuf, extpos + 1);
    }
    else
    {
        snprintf(newname, STRLEN, "%s%s", outf_base, nbuf);
    }
}

static void checkPerformance(t_inputrec *ir)
{
    int       i, j, n = 100000;
    real     *values = NULL;
    real      res;
    gmx_rng_t rng;

    fprintf(stderr,
            "\n=== Testing the performance of the erf and exp functions ===\n");
    rng = gmx_rng_init(ir->ld_seed);

    /* Fill the array with random values */
    snew(values, n);
    for (i = 0; i < n; i++)
    {
        values[i] = 10 * (gmx_rng_uniform_real(rng) - 0.5);
    }

    real        *x = NULL;
    snew(x, 2);
    gmx_cycles_t c1, c2;

    for (i = 1; i < 10; i++)
    {
        x[0] = values[0];
        x[1] = values[1];

        c1  = gmx_cycles_read();
        res = erf(x[1]) - erf(x[0]);
        c2  = gmx_cycles_read();
        fprintf(stderr, "Result: %15.10f   %lld cycles  normal erf()\n", res,
                c2 - c1);

        c1  = gmx_cycles_read();
        res = gmx_erf(x[1]) - gmx_erf(x[0]);
        c2  = gmx_cycles_read();
        fprintf(stderr, "Result: %15.10f   %lld cycles  gmx_erf()\n", res,
                c2 - c1);

        gmx_erf_vector(2, x);
        fprintf(stderr, "Result: %15.10f  vector\n\n", x[1] - x[0]);
    }

    for (j = 0; j < 2; j++)
    {
        for (i = 1; i < n; i *= 2)
        {
            gmx_erf_vector(i, values);
        }
        fprintf(stderr, "\n");

        for (i = 1; i < n; i *= 2)
        {
            gmx_exp_vector(i, values);
        }
        fprintf(stderr, "\n");
    }
    sfree(values);
}

/* If Intels vector math library is present (should be when Intel MKL is used),
 * set the precision here */
extern void gmx_set_vml_precision(FILE *fp)
{
#ifdef GMX_DENSFIT_MKL
    //  VML FUNCTION ACCURACY CONTROL
    //  VML_HA - when VML_HA is set, high accuracy VML functions are called
    //  VML_LA - when VML_LA is set, low accuracy VML functions are called
    //  VML_EP - when VML_EP is set, enhanced performance VML functions are called
    //
    //  NOTE: VML_HA, VML_LA and VML_EP must not be used in combination
    vmlSetMode( VML_EP );

    if (fp)
    {
        fprintf(fp, "# Using Intel VML acceleration for density fitting.\n");
    }
#else
    if (fp)
    {
        fprintf(fp, "# Not compiled with Intel MKL support.\n");
    }
#endif

}

/* Check whether densfit->nstfit is overwritten by environment variable */
static void get_nstfit_from_env(t_densfit *densfit, t_commrec *cr, int nstlist)
{
    char *env;
    char  buf[STRLEN];
    int   env_nstfit;

    FILE *fp = densfit->df->out_dens;

    if (MASTER(cr))
    {
        env = getenv(NSTFIT_ENVVAR);

        if (env != NULL)
        {
            sprintf(buf, "Getting nstfit from environment variable %s=%s",
                    NSTFIT_ENVVAR, env);
            fprintf(stderr, "%s%s\n", FitStr, buf);

            sscanf(env, "%d", &env_nstfit);

            if (env_nstfit >= 0)
            {
                densfit->nstfit = env_nstfit;
                if (fp)
                {
                    fprintf(fp, "# %s\n", buf);
                }
            }
            else
            {
                fprintf(stderr,
                        "WARNING: Could not get a meaningful value for environment variable %s!\n",
                        NSTFIT_ENVVAR);
            }
        }
    }

    if (PAR(cr))
    {
        gmx_bcast(sizeof(densfit->nstfit), &densfit->nstfit, cr);
    }

    /* We check nstfit versus nstlist here so we can output a warning on stderr
     * AND a note in the density fitting output file */
    if (MASTER(cr))
    {
        if (densfit->nstfit > 1 && (densfit->nstfit % nstlist != 0))
        {
            fprintf(stderr,
                    "\n%s\n"
                    "WARNING: The frequency at which the fitting forces are calculated (nstfit=%d)\n"
                    "         is not a multiple of the neighbor searching frequency (nstlist=%d).\n"
                    "         Fitting forces will be calculated when the time step is a multiple\n"
                    "         of nstfit or nstlist, that could also be irregular intervals.\n\n",
                    FitStr, densfit->nstfit, nstlist);

            fprintf(fp,
                    "# Note: will update fitting forces when the time step is a multiple of nstlist=%d or nstfit=%d.\n",
                    nstlist, densfit->nstfit);
        }

        if (densfit->nstfit < 1)
        {
            fprintf(stderr,
                    "\n%s\n"
                    "WARNING: The frequency at which the density fitting potential is calculated (nstfit=%d) is < 1!\n"
                    "         No density fitting will be performed and no fitting output will be written.\n\n",
                    FitStr, densfit->nstfit);
            fprintf(fp, "# No output follows since nstfit < 1\n");
        }
    }
}

/* Called from mdrunner */
extern void init_density_fitting(FILE *fplog, t_inputrec *ir, int nfile,
                                 const t_filenm fnm[], gmx_mtop_t *mtop, rvec *x, matrix box,
                                 t_commrec *cr, const output_env_t oenv, unsigned long Flags,
                                 gmx_bool bVerbose)
{
    int            i;
    t_densfit     *densfit;
    t_gmx_densfit *df; /* Pointer to the density fitting buffer variables */
    real           grid_spacing;

    /* To be able to make the density fitting molecule whole: */
    rvec *x_pbc         = NULL;
    rvec *x_densfit_pbc = NULL;

    FitStr = FitStrMDRUN;

    if ((PAR(cr)) && !DOMAINDECOMP(cr))
    {
        gmx_fatal(FARGS, "%sOnly implemented for domain decomposition.\n",
                  FitStr);
    }

    if (MASTER(cr) && bVerbose)
    {
        fprintf(stdout, "%sInitializing ...\n", FitStr);
    }

    densfit = ir->densfit;

    /* In serial runs, the local atom indices are equal to the global ones */
    if (!PAR(cr))
    {
        densfit->nat_loc = densfit->nat;
        densfit->ind_loc = densfit->ind;
    }

    new_map_sim_from_ref(densfit, densfit->map_ref);

    /* Make space for the density fitting private data, allocate x, w, and f arrays, ... */
    snew(x_densfit_pbc, densfit->nat); /* There ... */
    if (MASTER(cr))
    {
        /* Remove pbc and prepare a whole molecule for new_t_gmx_densfit */
        snew(x_pbc, mtop->natoms); /* There ... */
        m_rveccopy(mtop->natoms, x, x_pbc);
        do_pbc_first_mtop(NULL, ir->ePBC, box, mtop, x_pbc);
        for (i = 0; i < densfit->nat; i++)
        {
            copy_rvec(x_pbc[densfit->ind[i]], x_densfit_pbc[i]);
        }
        sfree(x_pbc); /* ... and back again */
    }
    if (PAR(cr))
    {
        gmx_bcast(densfit->nat * sizeof(rvec), x_densfit_pbc, cr);
    }

    grid_spacing = get_map_spacing(densfit->map_ref,
                                   MASTER(cr) ? stdout : NULL);
    densfit->df = new_t_gmx_densfit(densfit->nat, densfit->sigma, densfit->dist,
                                    x_densfit_pbc, grid_spacing,
                                    MASTER(cr) && bVerbose, Flags & MD_APPENDFILES, PAR(cr));

    sfree(x_densfit_pbc); /* ... and back again */
    df = densfit->df;

    if (MASTER(cr))
    {
        df->out_dens = open_densfit_out(opt2fn("-d", nfile, fnm), densfit,
                                        oenv);
        df->fn_map = opt2fn("-mo", nfile, fnm);
    }
    else
    {
        df->out_dens = NULL;
        df->fn_map   = NULL;
    }

    /* Check whether densfit->nstfit is overwritten by environment variable */
    get_nstfit_from_env(densfit, cr, ir->nstlist);

    if (MASTER(cr) && !df->bAppend)
    {
        please_cite(fplog, "Tama2008");
    }

    gmx_set_vml_precision(df->out_dens);

    /* This sum does not change, therefore we precalculate it */
    densfit->df->sum_rho_exp2 = calc_sum_rhoA_rhoB(densfit->map_ref,
                                                   densfit->map_ref);

    //checkPerformance(ir);
}

/* Make a new x array that only contains atoms to be spread */
extern void assemble_atoms_for_spread(t_densfit *densfit, rvec x[])
{
    int            i;
    t_gmx_densfit *df; /* Pointer to the density fitting buffer variables */

    df = densfit->df;

    for (i = 0; i < densfit->nat; i++)
    {
        copy_rvec(x[densfit->ind[i]], df->x_assembled[i]);
    }
}

void dump_local_indices(t_densfit *densfit, t_commrec *cr)
{
    FILE          *fpout = NULL;
    char           fn[STRLEN];
    int            i;
    t_gmx_densfit *df;

    df = densfit->df;
    sprintf(fn, "ind_node%d.txt", cr->nodeid);
    fpout = fopen(fn, "w");

    for (i = 0; i < densfit->nat; i++)
    {
        fprintf(fpout, "%4d %4d %12.3e %12.3e %12.3e\n", i, df->c_ind[i],
                df->x_assembled[i][XX], df->x_assembled[i][YY],
                df->x_assembled[i][ZZ]);
    }

    fclose(fpout);
}

extern void dd_make_local_df_indices(gmx_domdec_t *dd, t_densfit *densfit)
{

    /* Local atoms of the density fitting group (these atoms will be spread) */
    dd_make_local_group_indices(dd->ga2la, densfit->nat, densfit->ind,
                                &densfit->nat_loc, &densfit->ind_loc, &densfit->nalloc_loc,
                                densfit->df->c_ind);

    /* Indicate that the group's shift vectors for this structure need to be updated
     * at the next call to communicate_group_positions, since obviously we are in a NS step */
    densfit->df->bUpdateShifts = TRUE;
}

/* Add the density fitting forces to the MD forces and output */
extern void add_densfit_forces(t_inputrec *ir, rvec *f, t_commrec *cr,
                               gmx_int64_t step, real time)
{
    int            ii, l;
    t_densfit     *densfit = ir->densfit;
    t_gmx_densfit *df      = densfit->df;

    if (densfit->nstfit < 1)
    {
        return;
    }

    /* Diagnostic output */

    /*
       if (df->out_dens && do_per_step(step, densfit->nstout)) {
        char fn[STRLEN];
        make_filename("forces.txt", 0, step, fn);
        //dump_x(densfit, cr, step);
        dump_f(fn, densfit, cr);
       }
     */

    if (df->out_dens && do_per_step(step, densfit->nstout))
    {
        fprintf(df->out_dens, "%12.5e%12.7f%12.5e\n", time, df->cc,
                densfit->k * (1.0 - df->cc));
    }

    for (l = 0; l < densfit->nat_loc; l++)
    {
        /* Get the right index of the local force, since typically not all local
         * atoms are subject to density fitting forces */
        ii = densfit->ind_loc[l];

        /* Add to local force */
        rvec_inc(f[ii], df->f_loc[l]);
    }
}

static void sum_grid(t_mapdata *map, t_commrec *cr)
{
    int nx, ny, nz;

    nx = map->map_dim[XX];
    ny = map->map_dim[YY];
    nz = map->map_dim[ZZ];

    gmx_sum(nx * ny * nz, map->vox, cr);
}

extern void do_densfit(gmx_int64_t step, gmx_bool bOutputMap, t_inputrec *ir,
                       t_commrec *cr, rvec x[], matrix box, gmx_wallcycle_t wcycle)
{

    const char   *fn      = NULL;
    t_densfit    *densfit = NULL;
    gmx_densfit_t df      = NULL; /* Pointer to the density fitting buffer variables */
    int           istart, nspread, nnodes;
    gmx_cycles_t  cycles_comp;    /* Cycles for the density fitting computations
                                     only, does not count communication. This
                                     counter is used for load-balancing              */

    densfit = ir->densfit;
    df      = densfit->df;

    /* If the user requested to have the density fitting forces calculated only
     * every N steps, we can skip the expensive do_densfit routine. However, we
     * must always calculate the fitting forces
     * - at the first step (obviously)
     * - after domain decompositions since the local atoms have changed
     */
    if ((!do_per_step(step, densfit->nstfit) && (FALSE == df->bUpdateShifts))
        || (densfit->nstfit < 1))
    {
        return;
    }

    /**************************************************************************/
    /* ASSEMBLE THE POSITIONS OF THE ATOMS USED FOR DENSITY FITTING           */

    /* Broadcast the DF positions such that every node has all of them
     * Every node contributes its local positions x and stores it in
     * the collective df->x_assembled array. Do all the communication first!  */
    wallcycle_start(wcycle, ewcDENSFIT_COMM);

    communicate_group_positions(cr, df->x_assembled, df->x_shifts,
                                df->extra_shifts, df->bUpdateShifts, x, densfit->nat,
                                densfit->nat_loc, densfit->ind_loc, df->c_ind, df->x_old, box);

    /* If bUpdateShifts was TRUE then the shifts have just been updated in
     * communicate_group_positions. We do not need to update the shifts until
     * the next NS step */
    df->bUpdateShifts = FALSE;

    /* Put all atoms in the box (TODO: if we do that, we do not need to construct
     * a whole DF group before with communicate_group_positions!)
     */
    if (ir->ePBC != epbcNONE)
    {
        put_atoms_in_box_omp(ir->ePBC, box, densfit->nat, df->x_assembled);
    }

    wallcycle_stop(wcycle, ewcDENSFIT_COMM);
    /* Now all nodes have all of the DF positions in df->x_assembled */

    /*                                                                        */
    /**************************************************************************/

    /**************************************************************************/
    /* SPREAD THE ATOMS ON THE DENSITY GRID                                   */
    /* Produces the density map. Each node spreads an equal part of all atoms,
     * this way we do not need load balancing here                            */
    wallcycle_start(wcycle, ewcDENSFIT_SPREAD);

    if (PAR(cr))
    {
        nnodes = cr->nnodes - cr->npmenodes;

        nspread = ceil((real) densfit->nat / nnodes);
        istart  = cr->nodeid * nspread;

        if ((nnodes - 1) == cr->nodeid)
        {
            nspread = densfit->nat - nspread * (nnodes - 1);
        }
    }
    else
    {
        istart  = 0;
        nspread = densfit->nat;
    }

    assert(istart >= 0);
    assert(nspread >= 0);

/*    spread_atoms_low(&df->x_assembled[istart], nspread, box, densfit);

    wallcycle_stop(wcycle, ewcDENSFIT_SPREAD);

    if (PAR(cr)) {
        wallcycle_start(wcycle, ewcDENSFIT_SUM_GRID);
        sum_grid(densfit->map_sim, cr); // Sum the density grid across all nodes
        wallcycle_stop(wcycle, ewcDENSFIT_SUM_GRID);
    }
 */
    /*                                                                        */

    /**************************************************************************/
    /* CALCULATE THE FORCES FROM THE DENSITY FITTING POTENTIAL                */
    /* Each node calculates the forces on its home atoms                      */

    /* Start to count cycles for the dynamic load balancing now */
    cycles_comp = gmx_cycles_read();

//    long int c1 = gmx_cycles_read();
    wallcycle_start(wcycle, ewcDENSFIT_FORCES);

#ifdef WITHFGT
    do_densfit_forces_fgt(densfit, box);
#else
    do_densfit_forces(densfit, box);
#endif

    wallcycle_stop(wcycle, ewcDENSFIT_FORCES);
//    fprintf(stderr, "--- Forces took %g M cycles\n", 0.000001*(gmx_cycles_read()-c1));

    /**************************************************************************/

    if (MASTER(cr))
    {
        if (bOutputMap)
        {
            char fn_with_step[STRLEN];

            if (densfit->bKeepAndNumberMaps)
            {
                make_filename(df->fn_map, 0, step, fn_with_step);
                fn = fn_with_step;
            }
            else
            {
                fn = df->fn_map;
            }
            gmx_do_map_ccp4(FALSE, &densfit->map_sim, fn, TRUE, df->bVerbose,
                            FALSE);
            gmx_do_map_ccp4(FALSE, &densfit->map_ref, "map_ref.ccp4", TRUE,
                            df->bVerbose, FALSE);
        }

        /* Calculate the correlation coefficient versus the reference map */
        df->cc = calc_correlation_coeff(densfit, NULL);
        if (do_per_step(step, 10))
        {
            fprintf(stderr,
                    "\rDensity fitting map correlation coefficient: %g\n",
                    df->cc);
        }
    }

    /* Stop the density fitting cycle counter and add the computation-only
     * cycles to the force cycles for load balancing */
    cycles_comp = gmx_cycles_read() - cycles_comp;

    if (DOMAINDECOMP(cr) && wcycle)
    {
        dd_cycles_add(cr->dd, cycles_comp, ddCyclF);
    }
    /*                                                                        */
    /**************************************************************************/
}
