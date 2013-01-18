/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
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

#ifndef _splitter_h
#define _splitter_h
#include "visibility.h"
#include "typedefs.h"
#include "types/inputrec.h"

#ifdef __cplusplus
extern "C" {
#endif

GMX_LIBGMX_EXPORT
void split_top(FILE *fp, int nnodes, gmx_localtop_t *top,
               t_inputrec *ir, t_block *mols,
               real *capacity, int *mulitnr_cgs, int **multinr_nre,
               int *left_range, int *right_range);
/* Split the topology (blocks and forces, based on charge groups
 * and shake blocks.
 * The capacity is releated to the capacity of each node. If all numbers are
 * equal, load will be distributed equally. If not some (the higher ones)
 * will get more than others. The sum of capacities should be 1.
 * Info is written to the file pointer fp.
 */

GMX_LIBGMX_EXPORT
void gen_sblocks(FILE *fp, int at_start, int at_end,
                 t_idef *idef, t_blocka *sblock,
                 gmx_bool bSettle);
/* Generate shake blocks from the constraint list. Set bSettle to yes for shake
 * blocks including settles. You normally do not want this.
 */

#ifdef __cplusplus
}
#endif

#endif
