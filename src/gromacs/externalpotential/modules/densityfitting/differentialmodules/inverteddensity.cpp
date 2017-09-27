// /*
//  * This file is part of the GROMACS molecular simulation package.
//  *
//  * Copyright (c) 2017, by the GROMACS development team, led by
//  * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
//  * and including many others, as listed in the AUTHORS file in the
//  * top-level source directory and at http://www.gromacs.org.
//  *
//  * GROMACS is free software; you can redistribute it and/or
//  * modify it under the terms of the GNU Lesser General Public License
//  * as published by the Free Software Foundation; either version 2.1
//  * of the License, or (at your option) any later version.
//  *
//  * GROMACS is distributed in the hope that it will be useful,
//  * but WITHOUT ANY WARRANTY; without even the implied warranty of
//  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  * Lesser General Public License for more details.
//  *
//  * You should have received a copy of the GNU Lesser General Public
//  * License along with GROMACS; if not, see
//  * http://www.gnu.org/licenses, or write to the Free Software Foundation,
//  * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
//  *
//  * If you want to redistribute modifications to GROMACS, please
//  * consider that scientific software is very special. Version
//  * control is crucial - bugs must be traceable. We will be happy to
//  * consider code for inclusion in the official distribution, but
//  * derived work must not be called official GROMACS. Details are found
//  * in the README & COPYING files - if they are missing, get the
//  * official version at http://www.gromacs.org.
//  *
//  * To help us fund GROMACS development, we humbly ask that you cite
//  * the research papers on the package. Check out http://www.gromacs.org.
//  */
//
// #include "inverteddensity.h"
// #include "gromacs/externalpotential/atomgroups/wholemoleculegroup.h"
// #include "gromacs/fileio/json.h"
// #include <memory>
// #include <string>
// #include "gromacs/math/volumedata/gridmeasures.h"
//
// #include "gromacs/utility/gmxomp.h"
//
//
// namespace gmx
// {
//
// // const Field<real> &InvertedDensity::evaluateDensityDifferential(
// //         const Field<real> & /*comparant*/, const Field<real> &reference)
// // {
// //     return reference;
// // }
// //
// // real
// // InvertedDensity::evaluateStructureDensityPotential(
// //         const std::vector<RVec> &coordinates, const std::vector<real> &weights,
// //         const Field<real> &reference, const RVec &translation,
// //         const Quaternion &orientation)
// // {
// //     RVec centerOfMass {
// //         0, 0, 0
// //     };
// //     real weight_sum = 0;
// //     for (size_t atomIndex = 0; atomIndex != coordinates.size(); ++atomIndex)
// //     {
// //         RVec weighted_coordinate;
// //         svmul(weights[atomIndex], coordinates[atomIndex], weighted_coordinate);
// //         weight_sum += weights[atomIndex];
// //         rvec_inc(centerOfMass, weighted_coordinate);
// //     }
// //     svmul(1./weight_sum, centerOfMass, centerOfMass);
// //     auto referenceField<real> = Field<real>(reference);
// //     auto sumOverGrid       = [referenceField<real>, centerOfMass, translation, orientation] (const real sum, const RVec r){
// //             return sum + referenceField<real>.getLinearInterpolationAt(orientation.shiftedAndOriented(r, centerOfMass, translation));
// //         };
// //     return std::accumulate(std::begin(coordinates), std::end(coordinates), 0., sumOverGrid);
// // };
// //
// // void
// // InvertedDensity::parseOptions_(const std::string &options)
// // {
// //     json::Object parsed_json {
// //         options
// //     };
// //     if (parsed_json.has("n_threads"))
// //     {
// //         number_of_threads_            = std::stof(parsed_json["n_threads"]);
// //     }
// //     else
// //     {
// //         fprintf(stderr, "\n No number of threads provided, taking the maximum number of available threads.\n");
// //         number_of_threads_ = gmx_omp_get_max_threads();
// //     };
// // }
// //
// // real
// // InvertedDensity::evaluateGroupDensityPotential(
// //         const WholeMoleculeGroup &atoms,
// //         const Field<real> &reference, const RVec &translation,
// //         const Quaternion &orientation)
// // {
// //     real result       = 0;
// //     RVec centerOfMass = atoms.weightedCenterOfMass();
// // #pragma omp parallel num_threads(number_of_threads_) reduction(+:result)
// //     {
// //         int           thread           = gmx_omp_get_thread_num();
// //         auto          beginThreadAtoms = atoms.begin(thread, number_of_threads_);
// //         auto          endThreadAtoms   = atoms.end(thread, number_of_threads_);
// //         for (auto atom = beginThreadAtoms; atom != endThreadAtoms; ++atom)
// //         {
// //             result += Field<real>(reference).getLinearInterpolationAt(orientation.shiftedAndOriented(*(*atom).xTransformed, centerOfMass, translation));
// //         }
// //     }
// //     return result;
// // };
//
// /**************************INFO CLasses****************************************/
//
// std::string InvertedDensityDifferentialPotentialInfo::name = std::string("inverted-density");
//
// std::unique_ptr<IStructureDensityPotentialProvider>
// InvertedDensityDifferentialPotentialInfo::create()
// {
//     return std::unique_ptr<InvertedDensity>(new InvertedDensity);
// }
//
// } /* gmx */
