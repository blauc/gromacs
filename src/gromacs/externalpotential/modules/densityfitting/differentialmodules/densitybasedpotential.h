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

#ifndef GMX_EXTERNALPOTENTIAL_COMMON_H
#define GMX_EXTERNALPOTENTIAL_COMMON_H
#include "gmxpre.h"
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "gromacs/math/quaternion.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/math/griddata/field.h"
#include "gromacs/utility/real.h"

#include "gromacs/externalpotential/modules/densityfitting/potentialprovider.h"

namespace gmx
{

class WholeMoleculeGroup;
class Quaternion;

class DensitySpreader;
template<int N> class GridWithTranslation;

class densityBasedPotential : public PotentialEvaluator
{
    public:
        explicit densityBasedPotential(const DensitySpreader &spreader, bool selfSpreading = true);
        ~densityBasedPotential() = default;
        real potential(const std::vector<RVec> &coordinates,
                       const std::vector<real> &weights,
                       const FieldReal3D &reference,
                       const RVec &translation = {0, 0, 0},
                       const Quaternion &orientation = {{1, 0, 0}, 0},
                       const RVec &centerOfRotation = {0, 0, 0}) const override;
        virtual real densityDensityPotential(const FieldReal3D &reference,
                                             const FieldReal3D &comparant) const = 0;
    protected:
        const DensitySpreader &spreader_;
        bool                   selfSpreading_;
};


// make a PotentialEvaluator member
class densityBasedForce : public ForceEvaluator
{
    public:
        explicit densityBasedForce(const DensitySpreader &spreader, real sigma_differential, int n_threads, bool selfSpreading = false);
        ~densityBasedForce() = default;

        void force(std::vector<RVec> &force, const std::vector<RVec> &coordinates,
                   const std::vector<real> &weights,
                   const FieldReal3D &reference,
                   const RVec &translation = {0, 0, 0},
                   const Quaternion &orientation = {{1, 0, 0}, 0},
                   const RVec &centerOfRotation = {0, 0, 0}) const override;

    protected:
        virtual void setDensityDifferential(const FieldReal3D    &reference,
                                            const FieldReal3D    &comparant) const = 0;
        std::unique_ptr < FieldReal3D> differential_;
        real                           sigma_differential_;
        int                            n_threads_;
        const DensitySpreader         &spreader_;
        bool                           selfSpreading_;
};

}      /* gmx */

#endif /* end of include guard: GMX_EXTERNALPOTENTIAL_COMMON_H */
