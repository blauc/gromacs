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

#include <vector>
#include <memory>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"


struct t_commrec;
struct t_inputrec;
struct gmx_ga2la_t;
struct gmx_output_env_t;
struct gmx_mtop_t;

#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/gmxmpi.h"

namespace gmx
{
class MpiHelper;

/* \!brief
 * Adds routines to a group to ensure a whole molecule by shifting atom coordinates according to periodic boundary conditions.
 *
 * We copy the atom coordinates and add the respective box shifts to make the molecule whole are added to all atom positions.
 * An alternative approach would alter the atom iterator, and add the shifts when querying the atom's positions.
 * To the externalpotential class this class looks like a normal group that magically reports atom coordinates shifted such that molecules stay whole and atoms do not jump periodic boundaries.
 *
 * This group requires communication across all nodes in all neighborsearching steps.
 * A whole molecule is constructed using topology information in the first step,
 * then kept updated through a reference configuration of reference coordinates of all atoms from the group of the previous neighboursearching step.
 * If atoms move more than half the box size from one neighboursearching step to another, we assume they are shifted and keep track of the pbc shifts that
 * place them closest to their previous position.
 *
 * For this class to function properly
 *    update_shifts_and_reference
 *  should be called whenever atom coordinate shifts are expected to have occured.
 *
 */
class WholeMoleculeGroup : public Group
{
    public:

        WholeMoleculeGroup(const Group &base_group, std::shared_ptr<MpiHelper> mpi_helper, const matrix box, const int npbcdim);
        ~WholeMoleculeGroup();

        void update_shifts_and_reference(const rvec x[], const matrix box);
        void set_box(const matrix box);
        void set_x(const rvec x[]);
        void make_whole_molecule_reference(const rvec x[], const gmx_mtop_t *top_global, const int ePBC);

    private:
        void calculate_shifts_();
        void all_group_coordinates_to_master_();

        matrix                     box_;
        rvec                      *x_reference_;
        rvec                      *x_coll_; //< a collective array for the group coordinates used to collect coordinates on the master node
        ivec                      *shifts_;
        std::shared_ptr<MpiHelper> mpi_helper_;
        const int                  npbcdim_;
};



} // namespace gmx
