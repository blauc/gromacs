/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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

#include "common.h"

namespace gmx
{
namespace volumedata
{

void
commonDensityBased::intializeSpreaderIfNull_(const FiniteGrid &reference, const int n_threads, const int n_sigma, const real sigma)
{
    if (spreader_ == nullptr)
    {
        spreader_ = std::unique_ptr<DensitySpreader>( new DensitySpreader(reference, n_threads, n_sigma, sigma));
    }
};

real
commonDensityBased::evaluateStructureDensityPotential( const std::vector<RVec> &coordinates, const std::vector<real> &weights, const Field<real> &reference, const RVec &translation, const Quaternion &orientation)
{
    return evaluateDensityDensityPotential(spread_density_, reference);
};

real commonDensityBased::evaluateGroupDensityPotential(
        const WholeMoleculeGroup &atoms, const Field<real> &reference,
        const RVec &translation, const Quaternion &orientation)
{
    return evaluateDensityDensityPotential(spread_density_, reference);
};


void commonDensityBased::planCoordinates(const std::vector<RVec> &coordinates, const std::vector<real> &weights, const Field<real> &reference,  const RVec &translation, const Quaternion &orientation)
{
    intializeSpreaderIfNull_(reference, n_threads_, n_sigma_, sigma_);
    spread_density_ = spreader_->spreadLocalAtoms(as_rvec_array(coordinates.data()), weights, coordinates.size(), translation, orientation);
}
void commonDensityBased::planGroup(const WholeMoleculeGroup &atoms, const Field<real> &reference, const RVec &translation, const Quaternion &orientation)
{
    intializeSpreaderIfNull_(reference, n_threads_, n_sigma_, sigma_);
    spread_density_ = GridReal(spreader_->spreadLocalAtoms(atoms, translation, orientation));
    spread_density_.normalize();
}

std::vector<RVec>
commonDensityBased::evaluateCoordinateForces(const std::vector<RVec> &coordinates, const std::vector<real> &weights, const Field<real> &reference,  const RVec &translation, const Quaternion &orientation);
{
};

void commonDensityBased::evaluateGroupForces(const WholeMoleculeGroup &atoms, const Field<real> &reference, const RVec &translation, const Quaternion &orientation);
{

    differentialCalculator_(spread_density, reference);
    auto forceDensity = ForceDensity(differential, sigma_).getForce();
    std::array<volumedata::GridReal, 3> forceGrid { {
                                                        volumedata::GridReal(forceDensity[XX]), volumedata::GridReal(forceDensity[YY]), volumedata::GridReal(forceDensity[ZZ])
                                                    } };

#pragma omp parallel num_threads(number_of_threads_) shared(stderr,atoms,forceGrid) default(none)
    {
        int           thread             = gmx_omp_get_thread_num();
        // real          prefactor  = k_/(norm_simulated_*sigma_*sigma_);

        GroupIterator beginThreadAtoms = atoms->begin(thread, number_of_threads_);
        GroupIterator endThreadAtoms   = atoms->end(thread, number_of_threads_);

        for (auto atom = beginThreadAtoms; atom != endThreadAtoms; ++atom)
        {
            auto r = orientation_.shiftedAndOriented(*(*atom).xTransformed, centerOfMass_, translation_);
            auto f = RVec {
                forceGrid[XX].getLinearInterpolationAt(r), forceGrid[YY].getLinearInterpolationAt(r), forceGrid[ZZ].getLinearInterpolationAt(r)
            };
            svmul(*((*atom).properties), f, *(*atom).force);
            orientation_.rotate_backwards(*(*atom).force);
        }
    }
}

} /* volumedata */
} /* gmx */
