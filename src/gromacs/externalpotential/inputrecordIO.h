/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
#ifndef _EXTERNALPOTENTIAL_INPUTRECORDIO_
#define _EXTERNALPOTENTIAL_INPUTRECORDIO_

#include <list>
#include <string>

#include "gromacs/fileio/readinp.h"

struct ext_pot;
struct t_blocka;
struct t_commrec;

namespace gmx
{
namespace externalpotential
{
namespace inputrecordutils
{
/*! \brief
 * Helper function, that splits a string at a token
 */
std::list<std::string> split_string_at_token(const char * names, char token);
/*! \brief
 * Broadcast inputrecord data
 */
void bc_externalpotential(const t_commrec * cr, struct ext_pot * ep);

/*! \brief
 * Set the filenames and indexgroup names for all external potentials in grompp.
 *
 * \param[in] ninp_p mdp number of input parameters for STYPE-macro magic
 * \param[in] inp_p mdp input parameters for STYPE-macro magic
 * \param[in] ep all external potential data that goes into the input-record
 * \result groupnames for all types of external potentials
 * todo: enforce input consistency (eg. by writing checksums of input files to the .tpr file / inputrecord)
 */
char ** set_external_potential(int *ninp_p, t_inpfile **inp_p, struct ext_pot * ep );

/*! \brief
 * Generate an index group per applied external potential during pre-processing in grompp.
 *
 * \param[in] ext_pot all external potential data that goes into the input-record
 * \param[in] group_name the index group names of the external potentials
 * \param[in] groups the indexgroups as read by grompp
 * \param[in] all_group_names the group names as read by grompp
 */
void make_groups(struct ext_pot *ep, char **group_name, struct t_blocka *groups, char **all_group_names);

void pr_externalpotential(FILE *fp, int indent, struct ext_pot * pot);
}      // namepsace inputrecordutils
}      // namespace externalpotential
}      // namespace gmx
#endif /* end of include guard: _EXTERNALPOTENTIAL_INPUTRECORDIO_ */
