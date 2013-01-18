/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */

#ifndef _nb_free_energy_h_
#define _nb_free_energy_h_

#include "nb_kernel.h"
#include <typedefs.h>

void
gmx_nb_free_energy_kernel(t_nblist *                nlist,
                          rvec *                    x,
                          rvec *                    f,
                          t_forcerec *              fr,
                          t_mdatoms *               mdatoms,
                          nb_kernel_data_t *        kernel_data,
                          t_nrnb *                  nrnb);

real
    nb_free_energy_evaluate_single(real r2, real sc_r_power, real alpha_coul,
                                   real alpha_vdw, real tabscale, real *vftab,
                                   real qqA, real c6A, real c12A, real qqB,
                                   real c6B, real c12B, real LFC[2], real LFV[2], real DLF[2],
                                   real lfac_coul[2], real lfac_vdw[2], real dlfac_coul[2],
                                   real dlfac_vdw[2], real sigma6_def, real sigma6_min,
                                   real sigma2_def, real sigma2_min,
                                   real *velectot, real *vvdwtot, real *dvdl);

#endif
