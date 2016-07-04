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
#ifndef _PULLING_H
#define _PULLING_H

#include <memory>
#include <string>

#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/externalpotential/externalpotential.h"
#include "gromacs/externalpotential/modules.h"
#include "gromacs/mdtypes/commrec.h"

namespace gmx
{
struct GroupAtom;
namespace externalpotential
{
class Template : public ExternalPotential
{
    public:

        void do_potential(const matrix box, const rvec x[], const gmx_int64_t step);
        static std::unique_ptr<ExternalPotential> create();
        void read_input();
        void broadcast_internal();
        AtomProperties * single_atom_properties(t_mdatoms * mdatoms, gmx_localtop_t * topology_loc);
        void finish();
        bool do_this_step(gmx_int64_t step);
        void initialize(const matrix box, const rvec x[]);

    private:
        void ForceKernel_(GroupAtom &atom, const int &thread);
        Template();
        real              k_;
        RVec              com2_;
        std::vector<real> potential_;

};

class TemplateInfo
{
    public:
        static std::string name;
        static std::string shortDescription;
        static const int   numberIndexGroups;
        static const int   numberWholeMoleculeGroups;
        static externalpotential::ModuleCreator create;
};
} // namespace externalpotential
} // namespace gmx

#endif
