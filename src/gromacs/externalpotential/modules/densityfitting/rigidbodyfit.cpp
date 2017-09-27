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
#include "gromacs/math/volumedata/field.h"
#include "gromacs/math/volumedata/gridreal.h"
#include "gromacs/math/volumedata/operations/realfieldmeasure.h"
#include "gromacs/selection/centerofmass.h"
#include "gromacs/math/vec.h"

#include <algorithm>
namespace gmx
{

RigidBodyFitResult RigidBodyFit::fitCoordinates(
        const Field<real> &reference, const std::vector<RVec> &x,
        const std::vector<real> &weights,
        const PotentialEvaluatorHandle &fitPotentialProvider)
{

    // first guess of parameters
    auto             center_of_mass_density = RealFieldMeasure(reference).center_of_mass();
    rvec             center_of_geometry_structure;
    std::vector<int> indices(x.size());
    std::iota(indices.begin(), indices.end(), 0);
    gmx_calc_cog(nullptr, as_rvec_array(x.data()), indices.size(), indices.data(),
                 center_of_geometry_structure);

    RVec translation;

    rvec_sub(center_of_geometry_structure, center_of_mass_density, translation);

    Quaternion orientation {{
                                0, 0, 0
                            }, 0};

    real improvement       = 0;
    int  step              = 0;
    real translation_scale = 0.1;
    real potential {
        0
    };
    do
    {

        // auto gradient_translation = gradientTranslation_(reference, x,weights,fitPotentialProvider,translation,center_of_mass_density, orientation);
        // auto gradient_orientation = gradientOrientation_(reference, x,weights,fitPotentialProvider,translation,center_of_mass_density, orientation);

        // update x <- x + alpha * gradient
        // calculate potential
        //   smaller: break, alpha *= 2
        // larger: alpha <- -alpha
        // alpha /=2
        //
        // RVec translation_direction ;
        // Quaternion orientation_direction;
        // svmul(alpha, gradient_translation, translation_direction);


        fitPotentialProvider.potential(x, weights, reference, translation,
                                       orientation, center_of_mass_density);
        orientation[0] = 1;
        //   calculate potential gradient
        //   line search along potential gradient, until potential improves
    }
    while ((improvement > minimial_improvement_) && (step < max_steps_));


    return RigidBodyFitResult(translation, center_of_mass_density, orientation,
                              potential);
};

RigidBodyFitResult::RigidBodyFitResult(const RVec       &translation,
                                       const RVec       &center_of_rotation,
                                       const Quaternion &orientation,
                                       const real        bestFitPotential)
    : translation_ {translation}, center_of_rotation_ {
    center_of_rotation
},
orientation_ {
    orientation
}, potential_ {
    bestFitPotential
} {};

real RigidBodyFitResult::potential() const { return potential_; }

Quaternion RigidBodyFitResult::orientation() const { return orientation_; }

RVec RigidBodyFitResult::translation() const { return translation_; }

//
// void RigidBodyFit::fitWholeMoleculeGroup(
//         const Field<real> &reference, WholeMoleculeGroup *atoms,
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
//                     fprintf(stderr, "Optimizing Orientation: %3lu / %3lu = ",
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
//     fprintf(input_output()->output_file(), "#\tMinimization in Orientation
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
// RVec RigidBodyFit::fitTranslation()
// {
//
//     const real step_size                  = 0.02;
//     bool       succededReducingDivergence = false;
//     // calculate gradient vector
//     // normalize to step size
//     // move along gradient vector
//     // accept translation if new potential better
//     auto translationBeforeTrial = translation_;
//     auto spreader_              =
//     std::unique_ptr<DensitySpreader>(
//                 new DensitySpreader(reference, n_threads, n_sigma, sigma));
//     spreadLocalAtoms_(translationgroup);
//     // MrcFile().write("test.mrc", *simulated_density_);
//
//     forceCalculation(translationgroup);
//     auto force = translationgroup->local_force_sum();
//     if (mpi_helper() != nullptr)
//     {
//         mpi_helper()->sum_reduce_rvec(force);
//     }
//
//     unitv(force, force);
//     fprintf(input_output()->output_file(),
//             "#\tMinimization\n"
//             "#\tMinimization - New direction:\t[ %g %g %g ]\n"
//             "#\tMinimization - Step size:\t%g nm .\n",
//             force[XX], force[YY], force[ZZ], step_size);
//     RVec force_scaled;
//     if (mpi_helper() == nullptr || mpi_helper()->isMaster())
//     {
//         svmul(step_size, force, force_scaled);
//         rvec_inc(translation_, force_scaled);
//         fprintf(input_output()->output_file(),
//                 "#\tMinimization - Move attemted to\t[ %g %g %g ]\n",
//                 force_scaled[XX], force_scaled[YY], force_scaled[ZZ]);
//     }
//     if (mpi_helper() != nullptr)
//     {
//         mpi_helper()->broadcast(&translation_, 1);
//     }
//     auto new_divergence =
//         std::abs(KLDivergenceFromTargetOnMaster(translationgroup));
//     if (mpi_helper() == nullptr || mpi_helper()->isMaster())
//     {
//
//         fprintf(input_output()->output_file(),
//                 "#\tMinimization - Divergence before move: %g , after move:
//                 %g\n",
//                 divergenceToCompareTo, new_divergence);
//
//         if (new_divergence < divergenceToCompareTo)
//         {
//             succededReducingDivergence = true;
//             divergenceToCompareTo      = new_divergence;
//         }
//         else
//         {
//             translation_ = translationBeforeTrial;
//         }
//     }
//     return succededReducingDivergence;
// };

} /* gmx */
