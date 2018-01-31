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
#include "gromacs/math/griddata/operations/comparefields.h"
#include "gromacs/math/griddata/operations/gridinterpolator.h"
#include "gromacs/math/griddata/operations/realfieldmeasure.h"
#include "gromacs/utility/gmxomp.h"
#include <memory>
#include <string>

namespace gmx
{
/********************************************************************
 * KullbackLeiblerPotential
 */
KullbackLeiblerPotential::KullbackLeiblerPotential(
        const DensitySpreader &spreader)
    : DensityBasedPotential {spreader}
{}

real KullbackLeiblerPotential::densityDensityPotential(
        const GridDataReal3D &reference, const GridDataReal3D &comparant) const
{
    return ComapreGridData(reference, comparant).getKLSameGrid();
}

/********************************************************************
 * KullbackLeiblerForce
 */

KullbackLeiblerForce::KullbackLeiblerForce(const DensitySpreader &spreader,
                                           real                   sigma_differential,
                                           int                    n_threads)
    : DensityBasedForce {spreader, sigma_differential, n_threads}
{}

const GridDataReal3D &KullbackLeiblerForce::densityDifferential(
        const GridDataReal3D &comparant, const GridDataReal3D &reference) const
{
    differential_.setGrid(reference.getGrid().duplicate());
    auto sumSimulatedDensity = DataVectorMeasure(comparant).sum();

    auto densityGradientFunction = [sumSimulatedDensity](real densityExperiment,
                                                         real densitySimulation) {
            return (densitySimulation > 1e-15 && densityExperiment > 1e-15)
                   ? (densityExperiment / densitySimulation) *
                   (1 - densitySimulation / sumSimulatedDensity)
                   : 0;
        };

    std::transform(reference.begin(), reference.end(), comparant.begin(),
                   differential_.begin(), densityGradientFunction);
    return differential_;
}

/********************************************************************
 * KullbackLeiblerProvider
 */

ForceEvaluatorHandle
KullbackLeiblerProvider::planForce(const GridDataReal3D &reference,
                                   const std::string    &options)
{
    parseOptions_(options);
    spreader_ = std::unique_ptr<DensitySpreader>(
                new DensitySpreader(reference.getGrid(), n_threads_, n_sigma_, sigma_));
    force_evaluator_ = std::unique_ptr<KullbackLeiblerForce>(
                new KullbackLeiblerForce(*spreader_, sigma_, n_threads_));
    return ForceEvaluatorHandle(force_evaluator_.get());
};

PotentialEvaluatorHandle KullbackLeiblerProvider::planPotential(
        const std::vector<RVec> & /*coordinates*/,
        const std::vector<real> & /*weights*/, const GridDataReal3D &reference,
        const std::string &options)
{
    parseOptions_(options);
    spreader_ = std::unique_ptr<DensitySpreader>(
                new DensitySpreader(reference.getGrid(), n_threads_, n_sigma_, sigma_));
    potential_evaluator_ = std::unique_ptr<KullbackLeiblerPotential>(
                new KullbackLeiblerPotential(*spreader_));
    return PotentialEvaluatorHandle(potential_evaluator_.get());
};

PotentialAndForceHandle
KullbackLeiblerProvider::plan(const GridDataReal3D &reference,
                              const std::string    &options)
{
    parseOptions_(options);
    spreader_ = std::unique_ptr<DensitySpreader>(
                new DensitySpreader(reference.getGrid(), n_threads_, n_sigma_, sigma_));
    force_evaluator_ = std::unique_ptr<KullbackLeiblerForce>(
                new KullbackLeiblerForce(*spreader_, sigma_, n_threads_));

    potential_evaluator_ = std::unique_ptr<KullbackLeiblerPotential>(
                new KullbackLeiblerPotential(*spreader_));

    return PotentialAndForceHandle {
               ForceEvaluatorHandle(force_evaluator_.get()),
               PotentialEvaluatorHandle(potential_evaluator_.get())
    };
}

void KullbackLeiblerProvider::setCoordinates(
        const std::vector<RVec> &coordinates, const std::vector<real> &weights)
{
    spreader_->spreadLocalAtoms(coordinates, weights);
}

void KullbackLeiblerProvider::log_(const std::string &message)
{
    fprintf(stderr, "\n%s\n", message.c_str());
}

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
        sigma_ = 0.2;
        log_("Guessing mobility = " + std::to_string(sigma_) + " nm.");
    }
    if (parsed_json.has("n_threads"))
    {
        n_threads_ = std::stof(parsed_json["n_threads"]);
    }
    else
    {
        n_threads_ = gmx_omp_get_max_threads();
        log_("Using maximum of threads: " + std::to_string(n_threads_) + ".");
    }

    if (parsed_json.has("n_sigma"))
    {
        n_sigma_ = std::stof(parsed_json["n_sigma"]);
    }
    else
    {
        n_sigma_ = 5;
        log_("Using " + std::to_string(n_sigma_) + " to spread density.");
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

} /* gmx */
