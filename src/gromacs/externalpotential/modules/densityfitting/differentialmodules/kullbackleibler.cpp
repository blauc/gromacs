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

#include "kullbackleibler.h"
#include "gromacs/fileio/json.h"
#include "gromacs/math/volumedata/gridinterpolator.h"
#include "gromacs/math/volumedata/gridmeasures.h"
#include "gromacs/math/volumedata/gridreal.h"
#include "gromacs/utility/gmxomp.h"
#include <memory>
#include <string>

namespace gmx
{
namespace volumedata
{

void KullbackLeiblerForce::setDensityDifferential(
        const GridReal &comparant, const Field<real> &reference)
{
    differential_->copy_grid(reference);
    auto sumSimulatedDensity     = GridReal(comparant).properties().sum();
    auto densityGradientFunction = [sumSimulatedDensity](real densityExperiment,
                                                         real densitySimulation) {
            return (densitySimulation > 1e-15 && densityExperiment > 1e-15)
                   ? (densityExperiment / densitySimulation) *
                   (1 - densitySimulation / sumSimulatedDensity)
                   : 0;
        };

    std::transform(reference.access().begin(), reference.access().end(),
                   comparant.access().begin(), differential_->access().begin(),
                   densityGradientFunction);
}
real KullbackLeiblerPotential::densityDensityPotential(
        const Field<real> &reference, const Field<real> &comparant)
{
    return GridMeasures(reference).getKLSameGrid(comparant);
}
// real
// KullbackLeibler::evaluateDensityDensityPotential(
//         const Field<real> &comparant, const Field<real> &reference,
//         const RVec &translation,
//         const Quaternion &orientation)
// {
//     if (!comparant.sameGridInAbsTolerance(reference, 1e-10) &&
//     (norm(translation) > 1e-10) && orientation.norm() > 1e-10)
//     {
//         auto centerOfMass = GridReal(comparant).center_of_mass();
//         auto interpolated =
//         GridInterpolator(reference).interpolateLinearly(comparant,
//         translation, centerOfMass, orientation);
//         return GridMeasures(reference).getKLSameGrid(*interpolated);
//     }
//     return GridMeasures(reference).getKLSameGrid(comparant);
// };

KullbackLeiblerForce::KullbackLeiblerForce(const DensitySpreader &spreader, real sigma_differential, int n_threads, bool spreadSelf) : densityBasedForce {spreader, sigma_differential, n_threads, spreadSelf}
{
}

KullbackLeiblerPotential::KullbackLeiblerPotential(const DensitySpreader &spreader, bool spreadSelf) : densityBasedPotential {spreader, spreadSelf}
{
}

std::unique_ptr<ForceEvaluator> KullbackLeiblerProvider::planForce(
        const std::vector<RVec> & /*coordinates*/, const std::vector<real> & /*weights*/,
        const Field<real> &reference, const std::string &options,
        const RVec & /*translation*/, const Quaternion & /*orientation*/,
        const RVec & /*centerOfRotation*/)
{
    parseOptions_(options);
    if (spreader_ == nullptr)
    {
        spreader_ = std::unique_ptr<DensitySpreader>(new DensitySpreader(reference, n_threads_, n_sigma_, sigma_));
        return std::unique_ptr<KullbackLeiblerForce>(new KullbackLeiblerForce(*spreader_, sigma_, n_threads_, true));
    }
    else
    {
        return std::unique_ptr<KullbackLeiblerForce>(new KullbackLeiblerForce(*spreader_, sigma_, n_threads_, false));
    }
};


std::unique_ptr<PotentialEvaluator> KullbackLeiblerProvider::planPotential(
        const std::vector<RVec> & /*coordinates*/, const std::vector<real> & /*weights*/,
        const Field<real> &reference, const std::string &options,
        const RVec & /*translation*/, const Quaternion & /*orientation*/,
        const RVec & /*centerOfRotation*/)
{
    parseOptions_(options);
    if (spreader_ != nullptr)
    {
        spreader_ = std::unique_ptr<DensitySpreader>(new DensitySpreader(reference, n_threads_, n_sigma_, sigma_));
        return std::unique_ptr<KullbackLeiblerPotential>(new KullbackLeiblerPotential(*spreader_, true));
    }
    else
    {
        return std::unique_ptr<KullbackLeiblerPotential>(new KullbackLeiblerPotential(*spreader_, false));
    }
};

void KullbackLeiblerProvider::parseOptions_(const std::string &options)
{
    json::Object parsed_json {
        options
    };
    if (parsed_json.has("sigma"))
    {
        sigma_ = std::stof(parsed_json["sigma"]);
    }
    else
    {
        fprintf(stderr,
                "\n No mobility estimate given, guessing sigma = 0.2 nm . \n");
        sigma_ = 0.2;
    }
    if (parsed_json.has("sigma_"))
    {
        sigma_ = std::stof(parsed_json["sigma"]);
    }
    else
    {
        fprintf(stderr,
                "\n No mobility estimate given, guessing sigma = 0.2 nm . \n");
    }
    if (parsed_json.has("n_threads"))
    {
        n_threads_ = std::stof(parsed_json["n_threads"]);
    }
    else
    {
        fprintf(stderr, "\n No number of threads provided, taking the maximum "
                "number of available threads.\n");
        n_threads_ = gmx_omp_get_max_threads();
    }

    if (parsed_json.has("n_sigma"))
    {
        n_sigma_ = std::stof(parsed_json["n_sigma"]);
    }
    else
    {
        fprintf(stderr,
                "\n No density spread range provided, guessing n_sigma = 5 . \n");
    }
}

/**************************INFO CLasses****************************************/
std::string KullbackLeiblerPotentialInfo::name =
    std::string("kullback-leibler");

std::unique_ptr<IStructureDensityPotentialForceProvider>
KullbackLeiblerPotentialInfo::create()
{
    return std::unique_ptr<IStructureDensityPotentialForceProvider>(
            new KullbackLeiblerProvider());
}

} /* volumedata */
} /* gmx */
