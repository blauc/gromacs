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

#include "rigidbodyfit.h"
#include "gromac/math/vectypes.h"
#include "gromacs/math/quaternion.h"
#include "gromacs/utility/real.h"
#include "densityspreader.h"

namespace gmx
{
namespace volumedata
{

std::tuple<RVec, Orientiation> RigidBodyFit::fitCoordinates(
        const Field<real> &reference, std::vector<RVec> x,
        std::vector<real> weights,
        const IDifferentialPotentialProvider &fitPotentialProvider)
{
}


std::tuple<RVec, Orientiation> RigidBodyFit::fitWholeMoleculeGroup(
        const Field<real> &reference, WholeMoleculeGroup *atoms,
        const IDifferentialPotentialProvider &fitPotentialProvider)
{

    std::vector<real> upperHalfGrid {
        1.
    };
    std::vector<real> gridPoints {
        0
    };
    Quaternion bestOrientation                  = orientation_;
    RVec       bestTranslation                  = translation_;
    real       bestDivergence                   = std::abs(KLDivergenceFromTargetOnMaster(translationgroup));
    auto       bestDivergenceCurrentOrientation = bestDivergence;
    auto       totalNumberOrientations          = upperHalfGrid.size()*gridPoints.size()*gridPoints.size()*gridPoints.size();
    decltype(totalNumberOrientations) currentOrientationNumber = 0;

    for (auto q0 : upperHalfGrid)
    {
        for (auto q1 : gridPoints)
        {
            for (auto q2 : gridPoints)
            {
                for (auto q3 : gridPoints)
                {

                    fprintf(stderr, "Optimizing Orientation: %3lu / %3lu = ", ++currentOrientationNumber, totalNumberOrientations);

                    orientation_ = Quaternion(Quaternion::QVec {{q0, q1, q2, q3}});
                    orientation_.normalize();
                    fprintf(stderr, "[ %3g %3g %3g %3g ]", orientation_[0], orientation_[1], orientation_[2], orientation_[3]);

                    bestDivergenceCurrentOrientation =  std::abs(KLDivergenceFromTargetOnMaster(translationgroup));

                    while (optimizeTranslation(translationgroup, bestDivergenceCurrentOrientation) || optimizeOrientation(translationgroup, bestDivergenceCurrentOrientation))
                    {
                        fprintf(stderr, "\r\t\t\t\t\t\t\tcurrent fit = %g ; best fit = %g \t\t", bestDivergenceCurrentOrientation, bestDivergence);
                    }

                    if (bestDivergenceCurrentOrientation < bestDivergence)
                    {
                        bestTranslation = translation_;
                        bestOrientation = orientation_;
                        bestDivergence  = bestDivergenceCurrentOrientation;
                    }
                }
            }
        }
    }
    translation_ = bestTranslation;
    orientation_ = bestOrientation;
    orientation_.normalize();

};



Quaternion RigidBodyFit::fitOrientation()
{
    const real rotationAngle               = M_PI/128.;
    bool       succeededReducingDivergence = false;

    // calculate gradient vector
    // normalize to step size
    // move along gradient vector
    // accept if new potential better than old
    // move back along gradient vector

    spreadLocalAtoms_(atomgroup);
    sumReduceNormalize_();
    invertMultiplySimulatedDensity_();
    forceCalculation(atomgroup);

    auto torque_direction = atomgroup->local_torque_sum(centerOfMass_); //TODO: check for sign of rotation

    if (mpi_helper() != nullptr)
    {
        mpi_helper()->sum_reduce_rvec(torque_direction);
    }

    if (mpi_helper() == nullptr || mpi_helper()->isMaster())
    {
        unitv(torque_direction, torque_direction);
    }

    fprintf(input_output()->output_file(), "#\tMinimization in Orientation Space \n"
            "#\tMinimization - Torque Vector:\t[ %g %g %g ]\n",
            torque_direction[XX], torque_direction[YY], torque_direction[ZZ]);

    auto orientationBeforeTrial = orientation_;

    if (mpi_helper() == nullptr || mpi_helper()->isMaster())
    {
        orientation_ *= Quaternion(torque_direction, rotationAngle);
        orientation_.normalize();
    }
    if (mpi_helper() != nullptr)
    {
        mpi_helper()->broadcast(&orientation_[0], 4);
    }

    auto new_divergence =  std::abs(KLDivergenceFromTargetOnMaster(atomgroup));

    if (mpi_helper() == nullptr || mpi_helper()->isMaster())
    {

        fprintf(input_output()->output_file(), "#\tMinimization - Divergence before rotation: %g , after rotation: %g\n", divergenceToCompareTo, new_divergence);

        if (new_divergence < divergenceToCompareTo)
        {
            succeededReducingDivergence = true;
            divergenceToCompareTo       = new_divergence;
        }
        else
        {
            orientation_ = orientationBeforeTrial;
            if (mpi_helper() != nullptr)
            {
                mpi_helper()->broadcast(&orientation_[0], 4);
            }
        }
    }

    return succeededReducingDivergence;
}

RVec RigidBodyFit::fitTranslation()
{

    const real step_size                  = 0.02;
    bool       succededReducingDivergence = false;
    // calculate gradient vector
    // normalize to step size
    // move along gradient vector
    // accept translation if new potential better
    auto translationBeforeTrial = translation_;
    auto spreader_              = std::unique_ptr<volumedata::DensitySpreader>(
                new DensitySpreader(reference, n_threads, n_sigma, sigma));
    spreadLocalAtoms_(translationgroup);
    // volumedata::MrcFile().write("test.mrc", *simulated_density_);

    forceCalculation(translationgroup);
    auto force = translationgroup->local_force_sum();
    if (mpi_helper() != nullptr)
    {
        mpi_helper()->sum_reduce_rvec(force);
    }

    unitv(force, force);
    fprintf(input_output()->output_file(),
            "#\tMinimization\n"
            "#\tMinimization - New direction:\t[ %g %g %g ]\n"
            "#\tMinimization - Step size:\t%g nm .\n",
            force[XX], force[YY], force[ZZ], step_size);
    RVec force_scaled;
    if (mpi_helper() == nullptr || mpi_helper()->isMaster())
    {
        svmul(step_size, force, force_scaled);
        rvec_inc(translation_, force_scaled);
        fprintf(input_output()->output_file(),
                "#\tMinimization - Move attemted to\t[ %g %g %g ]\n",
                force_scaled[XX], force_scaled[YY], force_scaled[ZZ]);
    }
    if (mpi_helper() != nullptr)
    {
        mpi_helper()->broadcast(&translation_, 1);
    }
    auto new_divergence =
        std::abs(KLDivergenceFromTargetOnMaster(translationgroup));
    if (mpi_helper() == nullptr || mpi_helper()->isMaster())
    {

        fprintf(input_output()->output_file(),
                "#\tMinimization - Divergence before move: %g , after move: %g\n",
                divergenceToCompareTo, new_divergence);

        if (new_divergence < divergenceToCompareTo)
        {
            succededReducingDivergence = true;
            divergenceToCompareTo      = new_divergence;
        }
        else
        {
            translation_ = translationBeforeTrial;
        }
    }
    return succededReducingDivergence;
};


} /* volumedata */
} /* gmx */
