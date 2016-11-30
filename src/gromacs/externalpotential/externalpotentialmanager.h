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
#ifndef _externalpotentialmanager_h_
#define _externalpotentialmanager_h_

#include "externalpotential.h"
#include "gromacs/utility/gmxmpi.h"
#include "modules.h"

struct ext_pot;
struct t_mdatoms;
struct t_atomtypes;
struct gmx_mtop_t;

#include "gromacs/math/paddedvector.h"
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace gmx
{

class WholeMoleculeGroup;
class MpiHelper;
namespace externalpotential
{

/*! \brief
 * Manage iterations over external potentials.
 */
class Manager
{
    public:
        /*! \brief
         * Initialize the external potentials during the run (see runner.cpp).
         *
         * Check input consistency, test if files are read/writable,
         * setup coordinate communication and index groups
         *
         * \todo: Parse info from input-files in xml/json format?.
         * \todo: Check input file consistency with checksums.
         */
        explicit Manager(struct ext_pot *external_potential);

        /*! \brief
         * Trigger calculation of external potentials in the respective external
         * potential module classes.
         * \todo: this calculation could be triggered in parallel for all external
         * potentials
         *
         * \param[in] cr Communication record
         * \param[in] ir Input record
         * \result potential_ and force_ will be updated in all experimental input
         * modules, if applied at this step.
         */
        void do_potential(const matrix box, const rvec x[], const gmx_int64_t step);

        /*! \brief
         * Add the forces from the external potentials to the overall force in this
         * simulation step.
         * \todo: implement lambda-value reading for weighting the external potentials
         * \todo: implement exponential averaging for forces; make that default
         *
         * \param[in,out] f The updated forces
         * \param[in,out] vir The updated virial
         * \result contribution to the total potential from external potentials,
         * updated force and virial
         */
        real add_forces(rvec f[], tensor vir, gmx_int64_t step,
                        real forcefield_potential, real temperature);

        /*! \brief
         * Call the initialize of all external potentials that
         * sets the internal states after read_input and broadcast_internal.
         */
        void initialize(const matrix box, const PaddedRVecVector &x);

        /*! \brief
         * Keep the local indices in all applied potentials up-to-date after domain
         * decomposition.
         *
         * Is called whenever the system is partioned.
         * \param[in] dd the domain decomposition data
         * \param[in,out] the external potentials
         *
         * \result updated num_atoms_loc_, ind_loc_ and nalloc_loc_,
         */
        void dd_make_local_groups(gmx_ga2la_t *ga2la);

        void setAtomProperties(t_mdatoms *mdatom, gmx_localtop_t *types,
                               const gmx_mtop_t *mtop);

        void init_mpi(bool bMaster, MPI_Comm mpi_comm_mygroup, int masterrank);
        void init_input_output(struct ext_pot_ir **ir_data, int n_potentials,
                               const char *basepath);
        void init_groups(struct ext_pot_ir **ir_data, bool bParallel);
        void init_whole_molecule_groups(struct ext_pot_ir **ir_data,
                                        const gmx_mtop_t *top_global, const int ePBC,
                                        const matrix box, const PaddedRVecVector &x);
        void update_whole_molecule_groups(const rvec x[], const matrix box);

        void broadcast_internal();

        void finish();

    private:
        std::vector<real> calculate_weights(real V_forcefield, real temperature);
        std::vector < std::unique_ptr < ExternalPotential>> potentials_;
        std::vector<real>          V_external_;
        Modules                    modules_;
        real                       forcefield_potential_reference_;
        std::unique_ptr<MpiHelper> mpi_helper_;
        std::vector < std::shared_ptr < WholeMoleculeGroup>> whole_groups_;
};

} // namespace externalpotential
} // namespace gmx

#endif
