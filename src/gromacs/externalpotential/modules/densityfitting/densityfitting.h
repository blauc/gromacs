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
#ifndef _DENSITYFITTING_H
#define _DENSITYFITTING_H

#include <memory>
#include <string>

#include "gromacs/externalpotential/externalpotential.h"
#include "gromacs/externalpotential/modules.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/fileio/volumedataio.h"


namespace gmx
{

struct GroupAtom;
class Ifgt;

namespace volumedata
{
struct MrcMetaData;
class GridReal;
class FastGaussianGridding;
}

namespace externalpotential
{
class DensityFitting : public ExternalPotential
{
    public:

        static std::unique_ptr<ExternalPotential> create();
        void do_potential(const matrix box, const rvec x[], const gmx_int64_t step);
        void read_input();
        void initialize(const matrix box, const rvec x[]);
        void broadcast_internal();
        AtomProperties * single_atom_properties(t_mdatoms * mdatoms, gmx_localtop_t * topology_loc);
        void finish();
        real potential(std::vector<real> &P, std::vector<real> &Q);

    private:
        void do_force_plain(const rvec x, rvec force);
        void translation_removal(const matrix box, const rvec x[]);
        void translate_atoms_into_map_(const matrix box, const rvec x[]);
        void minimize_map_potential_through_translation_(const matrix box, const rvec x[]);
        RVec pbc_dist(const rvec x, const rvec y, const  matrix box);
        void inv_mul(std::vector<real> &target, const std::vector<real> & );
        void spread_density_(const rvec x[]);
        void sum_reduce_simulated_density_();
        void ForceKernel_KL(GroupAtom &atom, const int &thread);
        void plot_forces(const rvec x[]);
        DensityFitting();

        real k_;
        real sigma_;
        real n_sigma_;

        std::unique_ptr<volumedata::GridReal>             target_density_;
        std::unique_ptr<volumedata::GridReal>             simulated_density_;
        std::vector < std::unique_ptr < volumedata::GridReal>> force_density_;
        std::vector < std::unique_ptr < volumedata::FastGaussianGridding>> force_gauss_transform_;
        std::unique_ptr<volumedata::FastGaussianGridding> gauss_transform_;
        real                                              background_density_;
        volumedata::MrcMetaData                           meta_;
        std::vector<real>                                 potential_contribution_;
        RVec                                              translation_;
        std::vector<real>                                 reference_density_;
        real                                              relative_potential_ = 0;


};

class DensityFittingInfo
{
    public:
        static std::string name;
        static std::string shortDescription;
        static const int   numberIndexGroups;
        static externalpotential::ModuleCreator create;
};
} // namespace externalpotential
} // namespace gmx

#endif
