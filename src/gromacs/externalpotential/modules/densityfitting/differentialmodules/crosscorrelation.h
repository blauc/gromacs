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
// #ifndef GMX_EXTERNALPOTENTIAL_CROSSCORRELATION_H
// #define GMX_EXTERNALPOTENTIAL_CROSSCORRELATION_H
//
// // #include "gmxpre.h"
//
// #include "../densityspreader.h"
// #include "common.h"
// #include "gromacs/math/quaternion.h"
// #include "gromacs/utility/real.h"
// #include <memory>
// #include <string>
//
// namespace gmx
// {
// class WholeMoleculeGroup;
// {
//
// template <typename real> class Field;
// class DensitySpreader;
// class FiniteGrid;
// class CrossCorrelation : public commonDensityBased
// {
//     public:
//         CrossCorrelation() : commonDensityBased(std::bind(&CrossCorrelation::evaluateDensityDifferential, this, std::placeholders::_1, std::placeholders::_2)){};
//         ~CrossCorrelation() = default;
//         const Field<real> &evaluateDensityDifferential(const Field<real> &comparant, const Field<real> &reference);
//         real evaluateDensityDensityPotential( const Field<real> &comparant, const Field<real> &reference, const RVec &translation = {0, 0, 0}, const Quaternion &orientation = {{0, 0, 1}, 0});
//         void parseDifferentialOptionsString(const std::string &options);
//         void parseDensityDensityOptionsString(const std::string &options);
//         void parseStructureDensityOptionsString (const std::string &options);
//
//     private:
//         void parseOptions_(const std::string &options);
//         real correlationThreshold_                  = 0.;
//
// };
//
// /****************************INFO Classes**************************************/
//
// class CrossCorrelationDensityDensityInfo
// {
//     public:
//         static std::string name;
//         static std::unique_ptr<IStructureDensityPotentialProvider> create();
// };
//
// }      /* gmx */
//
// #endif /* end of include guard: GMX_EXTERNALPOTENTIAL_CROSSCORRELATION_H */
