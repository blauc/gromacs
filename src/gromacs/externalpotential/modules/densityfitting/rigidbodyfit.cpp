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

#include "gromacs/externalpotential/modules/densityfitting/potentialprovider.h"
#include "gromacs/math/griddata/griddata.h"
#include "gromacs/math/griddata/operations/gridinterpolator.h"
#include "gromacs/math/griddata/operations/realfieldmeasure.h"
#include "gromacs/math/griddata/rotatedgrid.h"
#include "gromacs/math/vec.h"
#include "gromacs/selection/centerofmass.h"

#include <algorithm>
namespace gmx
{

RigidBodyFitResult RigidBodyFit::fitCoordinates(
        const GridDataReal3D &reference, const std::vector<RVec> &coordinates,
        const std::vector<float> &weights,
        const PotentialAndForceHandle &fitPotentialProvider,
        RVec centerOfGeometry)
{
    // first guess of parameters
    auto densityCenterOfMass = RealGridDataMeasure(reference).center_of_mass();
    rvec_sub(densityCenterOfMass, centerOfGeometry, translation_);

    GridDataReal3D                 mobileReference(reference);
    GridWithTranslationOrientation mobileGrid(reference.getGrid());
    mobileGrid.setTranslation(
            {translation_[0], translation_[1], translation_[2]});
    mobileReference.setGrid(mobileGrid.duplicate());

    interpolateLinearly(reference, &mobileReference);
    real              potential;
    real              improvement = minimial_improvement_;
    std::vector<RVec> forces(coordinates.size());
    for (int step = 0;
         step < max_steps_ && (improvement >= minimial_improvement_); ++step)
    {

        // auto gradient_translation = gradientTranslation_(reference,
        // fitPotentialProvider.forceEvaluator, translation_, densityCenterOfMass,
        // orientation_);
        // auto gradient_orientation = gradientOrientation_(reference,
        // fitPotentialProvider.forceEvaluator,translation_,densityCenterOfMass,
        // orientation_);

        // update x <- x + alpha * gradient
        // calculate potential
        //   smaller: break, alpha *= 2
        // larger: alpha <- -alpha
        // alpha /=2
        //
        // RVec translation_direction ;
        // Quaternion orientation_direction;
        // svmul(alpha, gradient_translation, translation_direction);

        improvement = potential - fitPotentialProvider.potentialEvaluator.potential(mobileReference);
        potential  -= improvement;
        //   calculate potential gradient
        //   line search along potential gradient, until potential improves
    }

    return RigidBodyFitResult {
               translation_, densityCenterOfMass, orientation_,
               potential
    };
};

//
// void RigidBodyFit::fitWholeMoleculeGroup(
//         const GridDataReal3D &reference, WholeMoleculeGroup *atoms,
//         const IDifferentialPotentialProvider &fitPotentialProvider)
// {
//
//     std::vector<real> upperHalfGrid {
//         1.
//     };
//     std::vector<real> gridPoints {
//         0
//     };
//     Quaternion bestOrientation                  = orientation_;
//     RVec       bestTranslation                  = translation_;
//     real       bestDivergence                   =
//     std::abs(KLDivergenceFromTargetOnMaster(translationgroup));
//     auto       bestDivergenceCurrentOrientation = bestDivergence;
//     auto       totalNumberOrientations          =
//     upperHalfGrid.size()*gridPoints.size()*gridPoints.size()*gridPoints.size();
//     decltype(totalNumberOrientations) currentOrientationNumber = 0;
//
//     for (auto q0 : upperHalfGrid)
//     {
//         for (auto q1 : gridPoints)
//         {
//             for (auto q2 : gridPoints)
//             {
//                 for (auto q3 : gridPoints)
//                 {
//
//                     fprintf(stderr, "Optimizing orientation_: %3lu / %3lu =
//                     ",
//                     ++currentOrientationNumber, totalNumberOrientations);
//
//                     orientation_ = Quaternion(Quaternion::QVec {{q0, q1, q2,
//                     q3}});
//                     orientation_.normalize();
//                     fprintf(stderr, "[ %3g %3g %3g %3g ]", orientation_[0],
//                     orientation_[1], orientation_[2], orientation_[3]);
//
//                     bestDivergenceCurrentOrientation =
//                     std::abs(KLDivergenceFromTargetOnMaster(translationgroup));
//
//                     while (optimizeTranslation(translationgroup,
//                     bestDivergenceCurrentOrientation) ||
//                     optimizeOrientation(translationgroup,
//                     bestDivergenceCurrentOrientation))
//                     {
//                         fprintf(stderr, "\r\t\t\t\t\t\t\tcurrent fit = %g ;
//                         best fit = %g \t\t",
//                         bestDivergenceCurrentOrientation, bestDivergence);
//                     }
//
//                     if (bestDivergenceCurrentOrientation < bestDivergence)
//                     {
//                         bestTranslation = translation_;
//                         bestOrientation = orientation_;
//                         bestDivergence  = bestDivergenceCurrentOrientation;
//                     }
//                 }
//             }
//         }
//     }
//     translation_ = bestTranslation;
//     orientation_ = bestOrientation;
//     orientation_.normalize();
//
// };
//
//
//
// Quaternion RigidBodyFit::fitOrientation()
// {
//     const real rotationAngle               = M_PI/128.;
//     bool       succeededReducingDivergence = false;
//
//     // calculate gradient vector
//     // normalize to step size
//     // move along gradient vector
//     // accept if new potential better than old
//     // move back along gradient vector
//
//     spreadLocalAtoms_(atomgroup);
//     sumReduceNormalize_();
//     invertMultiplySimulatedDensity_();
//     forceCalculation(atomgroup);
//
//     auto torque_direction = atomgroup->local_torque_sum(centerOfMass_);
//     //TODO: check for sign of rotation
//
//     if (mpi_helper() != nullptr)
//     {
//         mpi_helper()->sum_reduce_rvec(torque_direction);
//     }
//
//     if (mpi_helper() == nullptr || mpi_helper()->isMaster())
//     {
//         unitv(torque_direction, torque_direction);
//     }
//
//     fprintf(input_output()->output_file(), "#\tMinimization in orientation_
//     Space \n"
//             "#\tMinimization - Torque Vector:\t[ %g %g %g ]\n",
//             torque_direction[XX], torque_direction[YY],
//             torque_direction[ZZ]);
//
//     auto orientationBeforeTrial = orientation_;
//
//     if (mpi_helper() == nullptr || mpi_helper()->isMaster())
//     {
//         orientation_ *= Quaternion(torque_direction, rotationAngle);
//         orientation_.normalize();
//     }
//     if (mpi_helper() != nullptr)
//     {
//         mpi_helper()->broadcast(&orientation_[0], 4);
//     }
//
//     auto new_divergence =
//     std::abs(KLDivergenceFromTargetOnMaster(atomgroup));
//
//     if (mpi_helper() == nullptr || mpi_helper()->isMaster())
//     {
//
//         fprintf(input_output()->output_file(), "#\tMinimization - Divergence
//         before rotation: %g , after rotation: %g\n", divergenceToCompareTo,
//         new_divergence);
//
//         if (new_divergence < divergenceToCompareTo)
//         {
//             succeededReducingDivergence = true;
//             divergenceToCompareTo       = new_divergence;
//         }
//         else
//         {
//             orientation_ = orientationBeforeTrial;
//             if (mpi_helper() != nullptr)
//             {
//                 mpi_helper()->broadcast(&orientation_[0], 4);
//             }
//         }
//     }
//
//     return succeededReducingDivergence;
// }
//
float RigidBodyFit::fitTranslation()
{

    const real step_size                  = 0.02;
    bool       succededReducingDivergence = false;
    // calculate gradient vector
    // normalize to step size
    // move along gradient vector
    auto translationBeforeTrial = translation_;

    forceCalculation(translationgroup);

    unitv(force, force);
    RVec force_scaled;
    svmul(step_size, force, force_scaled);
    rvec_inc(translation_, force_scaled);
    auto new_divergence =
        std::abs(KLDivergenceFromTargetOnMaster(translationgroup));

    if (new_divergence < divergenceToCompareTo)
    {
        succededReducingDivergence = true;
        divergenceToCompareTo      = new_divergence;
    }
    else
    {
        translation_ = translationBeforeTrial;
    }
    return succededReducingDivergence;
};

} /* gmx */
