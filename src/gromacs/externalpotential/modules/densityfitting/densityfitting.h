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
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/math/quaternion.h"

struct t_fileio;
struct gmx_mtop_atomlookup;

namespace gmx
{

struct GroupAtom;

namespace volumedata
{
struct MrcMetaData;
class GridReal;
class FastGaussianGridding;
class FastGaussianGriddingForce;
class GaussTransform;
}

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

        real single_atom_properties(GroupAtom * atom, t_mdatoms * mdatoms, gmx_localtop_t * topology_loc, const gmx_mtop_t * topology_global, const gmx_mtop_atomlookup * atom_lookup);
        void finish();

        real relative_kl_divergence(std::vector<real> &P, std::vector<real> &Q, std::vector<real> &Q_reference);

    private:
        void do_force_plain(const rvec x, rvec force);
        void translate_atoms_into_map_(WholeMoleculeGroup * translationgroup);
        void minimize_map_potential_through_translation_(const matrix box, const rvec x[]);
        RVec pbc_dist(const rvec x, const rvec y, const  matrix box);
        void inv_mul(std::vector<real> &target, const std::vector<real> & );
        void spreadLocalAtoms_(WholeMoleculeGroup * spreadgroup);
        void sum_reduce_simulated_density_();
        void ForceKernel_KL(GroupAtom &atom, const int &thread);
        void plot_forces(WholeMoleculeGroup * plotatoms);
        void initialize_target_density_();
        void initialize_buffers_();
        void initialize_spreading_();
        void setCenterOfMass(WholeMoleculeGroup * atomgroup);
        void KLForceCalculation_(WholeMoleculeGroup * fitatoms);
        std::string dumpParsedInput();
        void sumReduceNormalize_();
        void invertMultiplySimulatedDensity_();
        RVec shiftedAndOriented(const RVec x);
        void writeTranslatedCoordinates_(WholeMoleculeGroup * atoms, int step);
        real KLDivergenceFromTargetOnMaster(WholeMoleculeGroup * atomgroup);
        void alignComDensityAtoms();
        bool optimizeTranslation(WholeMoleculeGroup * translationgroup, real &divergenceToCompareTo);
        bool optimizeOrientation(WholeMoleculeGroup * atomgroup, real &divergenceToCompareTo);
        void doPotentialKL_( const matrix box, const rvec x[], const gmx_int64_t step);
        void doPotentialCC_( const matrix box, const rvec x[], const gmx_int64_t step);
        void doPotentialINV_( const matrix box, const rvec x[], const gmx_int64_t step);
        void CCMethod( const matrix box, const rvec x[], const gmx_int64_t step);
        void initializeKL_(const matrix box, const rvec x[]);
        void initializeCC_(const matrix box, const rvec x[]);
        void initializeINV_(const matrix box, const rvec x[]);
        real getTotalScatteringSum_(WholeMoleculeGroup * atomgroup);

        DensityFitting();

        real k_;
        real sigma_;
        real n_sigma_;
        real norm_simulated_;
        int  every_nth_step_;


        std::unique_ptr<volumedata::GridReal>             target_density_;
        std::unique_ptr<volumedata::GridReal>             simulated_density_;
        std::vector < std::unique_ptr < volumedata::GridReal>> force_density_;
        std::vector < std::unique_ptr < volumedata::FastGaussianGriddingForce>> force_gauss_transform_;
        std::vector < std::unique_ptr < volumedata::GaussTransform>> gauss_transform_;
        std::vector < std::unique_ptr < volumedata::GridReal>>                simulated_density_buffer_;
        real                    background_density_;
        volumedata::MrcMetaData meta_;
        std::string             trajectory_name_;
        t_fileio               *out_;
        bool                    isCenterOfMassCentered_;
        RVec                    translation_;
        std::vector<real>       reference_density_;
        real                    reference_divergence_ = 0;
        real                    k_factor_;
        int                     number_of_threads_;
        real                    absolute_target_divergence_;
        real                    optimalDeltaPotentialEnergy_;
        real                    exponentialDeltaEnergyAverage_;
        real                    maximumEnergyFluctuationPerAtom_;
        std::string             target_density_name_;
        RVec                    centerOfMass_;
        Quaternion              orientation_;
        std::string             fitMethod_;
        std::function<void(const matrix box, const rvec x[])> initialize_;
        std::function<void(const matrix box, const rvec x[], const gmx_int64_t step)> doPotential_;
        std::vector < std::unique_ptr < volumedata::GridReal>> invertedDensityForces_;
        bool bWriteXTC_;
        real totalScatteringSum_;

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
} // namespace externalpotential
} // namespace gmx

#endif
