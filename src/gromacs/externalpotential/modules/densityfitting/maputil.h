#ifndef _maputil_h
#define _maputil_h

#ifdef __cplusplus
extern "C" {
#endif

#include "gromacs/legacyheaders/types/inputrec.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/topology/topology.h"

#define WITHFGT
#ifdef WITHFGT
#include "fgt.h"
#endif

extern t_densfit *new_t_densfit(
        real      sigma,
        real      sigma_dist,
        real      k,        /* Spring constant */
        int       nat,
        atom_id  *ind,      /* Which of the atoms should be used for spreading */
        real      grid_spacing,
        gmx_bool  bVerbose);

extern void gmx_map_print_header(FILE *fpout, const char *fn);

extern void dump_f(const char *fn, t_densfit *densfit, t_commrec *cr);

extern void allocate_density_grid(
        int grid[3],
        float **vox);

/* Allocate memory for the simulated density grid */
extern void new_map_sim_from_ref(
        t_densfit *densfit,
        t_mapdata *map_ref);

extern void translate_pdb(matrix box, rvec x[], int natoms, gmx_bool bTranslate);

extern void gmx_do_map_ccp4(gmx_bool bRead, t_mapdata **ptr_map, const char *fn,
                            gmx_bool bOverwrite, gmx_bool bVerbose, FILE *DumpHeader);

/* Calculate the correlation coefficient of the two maps */
extern real calc_correlation_coeff(t_densfit *densfit, FILE *log);

/* Transform the positions x into a density map by replacing
 * each position by a Gaussian function of width sigma */
extern void spread_atoms(t_densfit  *densfit, matrix box);

extern void do_densfit_forces(t_densfit *densfit, matrix box);

/* Compute the forces from fitting to an electron density map,
 * based on the Kullback-Leibler divergence/Bayesian Cryo-EM refinement
 *
 * This is the real-space version of the fitting code.
 * For much higher performance use the Fourier-transform-version.
 *
 * Uses the "fine-grid" approximation in the second summation term
 * i.e. instead of spreading the density voxel then
 * integrating over the spread density per voxel, then summing per voxel
 * voxel summation is replaced by integration over all space right away
 * This is slow in real space, but should provide a significant speed-up
 * when Fourier-transformed.
 *
 * */

/* Initialize density fitting */
extern void init_density_fitting(
        FILE *fplog, t_inputrec *ir, int nfile, const t_filenm fnm[],
        gmx_mtop_t *mtop, rvec *x, matrix box, t_commrec *cr, const output_env_t oenv,
        unsigned long Flags, gmx_bool bVerbose);

/* Allocate and initialize the atomic weights array */
extern void init_weights(gmx_bool bVerbose, t_densfit *densfit);

/* If Intels vector math library is present, set the precision here */
extern void gmx_set_vml_precision(FILE *fp);

/* Make a new x array that only contains atoms to be spread */
extern void assemble_atoms_for_spread(
        t_densfit *densfit,
        rvec       x[]);

extern void dd_make_local_df_indices(gmx_domdec_t *dd, t_densfit *densfit);
/* Make a selection of the home atoms for the density fitting group.
 * Should be called at every domain decomposition. */

/* Compute the forces from fitting to an electron density map */
extern void do_densfit(
        gmx_int64_t     step,
        gmx_bool        bOutputMap,
        t_inputrec     *ir,
        t_commrec      *cr,
        rvec            x[],
        matrix          box,
        gmx_wallcycle_t wcycle);

extern void add_densfit_forces(t_inputrec *ir, rvec *f, t_commrec *cr, gmx_int64_t step, real time);

/* From the map entries determine and return the grid spacing */
extern real get_map_spacing(t_mapdata *map, FILE *log);

extern void couple_map_spacing_to_box(matrix box, t_densfit *densfit);

/* Append a number to the output file name */
extern void make_filename(const char *outf_name, int ndigit, int file_nr, char newname[STRLEN]);

extern void make_positive(t_mapdata *map_ref);

extern t_mapdata *rescale_map(t_mapdata *map_ref, real scale);


#ifdef __cplusplus
}
#endif

#endif
