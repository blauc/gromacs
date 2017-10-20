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

#include "gmxpre.h"

#include <memory>
#include <string>

#include "gromacs/externalpotential/externalpotential.h"
#include "gromacs/externalpotential/modules.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/fileio/griddataio.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/math/quaternion.h"
#include "gromacs/math/griddata/field.h"


struct t_fileio;
struct gmx_mtop_atomlookup;

namespace gmx
{

class GroupAtom;

struct MrcMetaData;
class FastGaussianGridding;
class FastGaussianGriddingForce;
class GaussTransform;
class IStructureDensityPotentialProvider;
class PotentialEvaluator;

namespace externalpotential
{

class FitAtomProperties : public IAtomProperties
{
    public:
        explicit FitAtomProperties(real weight);
        real weight();
    private:
        real weight_; //< The scattering factor of the atom determines how much weight is has for the gaussian spreading.
};

class DensityFitting : public ExternalPotential
{
    public:

        static std::unique_ptr<ExternalPotential> create();

        void do_potential(const matrix box, const rvec x[], const gmx_int64_t step);
        void read_input();
        void initialize(const matrix box, const rvec x[]);
        void broadcast_internal();

        bool do_this_step(gmx_int64_t step);

        real single_atom_properties(GroupAtom * atom, t_mdatoms * mdatoms, gmx_localtop_t * topology_loc, const gmx_mtop_t * topology_global);
        void finish();

    private:
        DensityFitting();

        void plot_forces(WholeMoleculeGroup * plotatoms);
        std::string dumpParsedInput();
        void writeTranslatedCoordinates_(WholeMoleculeGroup * atoms, int step);
        std::string logInitialisationString_(int nAtoms, int timeStepNs);


        real        k_;
        real        sigma_;
        real        n_sigma_;
        real        norm_simulated_;
        int         every_nth_step_;
        std::string options_;

        std::unique_ptr<IStructureDensityPotentialProvider>      potentialProvider_;
        std::unique_ptr < FieldReal3D>                           target_density_;
        std::unique_ptr < FieldReal3D>                           simulated_density_;
        real                    background_density_;
        MrcMetaData             meta_;
        std::string             trajectory_name_;
        t_fileio               *out_;
        bool                    isCenterOfMassCentered_;
        RVec                    translation_;
        RVec                    centerOfRotation_;
        std::vector<real>       reference_density_;
        real                    reference_divergence_ = 0;
        real                    k_factor_;
        int                     number_of_threads_;
        real                    absolute_target_divergence_;
        real                    optimalDeltaPotentialEnergy_;
        real                    exponentialDeltaEnergyAverage_;
        real                    maximumEnergyFluctuationPerAtom_;
        std::string             target_density_name_;
        Quaternion              orientation_;
        std::string             fitMethod_;
        std::function<void(const matrix box, const rvec x[])> initialize_;
        bool                    bWriteXTC_;
        real                    totalScatteringSum_;

};

class DensityFittingInfo
{
    public:
        static std::string name;
        static std::string shortDescription;
        static const int   numberIndexGroups;
        static const int   numberWholeMoleculeGroups;
        static externalpotential::ModuleCreator create;
};

}
} // namespace gmx

#endif
