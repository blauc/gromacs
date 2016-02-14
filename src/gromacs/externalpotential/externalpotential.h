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

#ifndef _externalpotential_h_
#define _externalpotential_h_

#include <cstdio>

#include <memory>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/classhelpers.h"

struct t_commrec;
struct t_inputrec;
struct ext_pot_ir;
struct gmx_output_env_t;
struct gmx_ga2la_t;
struct gmx_mtop_t;
struct t_mdatoms;
struct gmx_localtop_t;

namespace gmx
{

class GroupCoordinates;
class MpiHelper;
class Group;

class AtomProperties
{

};

namespace externalpotential
{

class ExternalPotentialIO;


class ExternalPotential
{
    // TODO: make methods that should only be called by the manager also only visible to the manager
    public:

        ~ExternalPotential();
        /*! \brief
         * Compute energies, forces, and virial.
         */
        virtual void do_potential(const matrix box, const rvec x[], const gmx_int64_t step) = 0;

        /*! \brief
         * Add up potential and virial contribution from all nodes.
         */
        real sum_reduce_potential_virial();

        /*! \brief
         * Sums the external potential forces on all nodes, after they have been previously calculated in do_potential.
         */
        void add_forces(rvec f[], gmx_int64_t step, real weight);

        /*! \brief
         * Sums the contribution to the virial from this potential from all nodes
         * TODO: rename to sum
         */
        void add_virial(tensor vir, gmx_int64_t step, real weight);

        /*! \brief
         * Distribute the atom indices over all nodes in the domain decomposition steps.
         */
        void dd_make_local_groups(gmx_ga2la_t  *ga2la);

        void set_mpi_helper(std::shared_ptr<MpiHelper> mpi);

        void add_group(std::shared_ptr<Group> &&group);

        void set_input_output(std::shared_ptr<ExternalPotentialIO> &&input_output);
        virtual void read_input()         = 0;
        virtual void broadcast_internal() = 0;
        void set_atom_properties(t_mdatoms * mdatoms, gmx_localtop_t * toplogy_loc);
        bool do_this_step(gmx_int64_t step);

        virtual void finish() = 0;

    protected:

        ExternalPotential ();

        /*! \brief
         * Picks the local atoms that are in group "group_index".
         *
         */
        std::shared_ptr<Group> group(const rvec x[], int group_index);

        void set_local_potential(real potential);

        void set_local_virial(tensor virial);

        std::shared_ptr<MpiHelper> mpi_helper();

        std::shared_ptr<ExternalPotentialIO> input_output();

        /*! \brief
         * Fill provide the atom with properties that migth be needed to calculate the external potential
         */
        virtual AtomProperties* single_atom_properties(t_mdatoms * mdatoms, gmx_localtop_t * toplogy_loc) = 0;

    private:
        class Impl;
        PrivateImplPointer<Impl> impl_;
};


} // end namspace externalpotential

} // end namespace gmx

#endif
