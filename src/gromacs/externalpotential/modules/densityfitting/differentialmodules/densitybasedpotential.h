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
#include "gromacs/math/volumedata/field.h"
#include "gromacs/utility/real.h"

#include "gromacs/externalpotential/modules/densityfitting/potentialprovider.h"

namespace gmx
{

class WholeMoleculeGroup;
class Quaternion;

namespace volumedata
{

class DensitySpreader;
class FiniteGrid;
// template<typename T> class Field;

class densityBasedPotential : public PotentialEvaluator
{
    public:
        explicit densityBasedPotential(const FiniteGrid &reference,
                                       const int n_threads = 1, const int n_sigma = 5,
                                       const real sigma = 0.2);

        real potential(const std::vector<RVec> &coordinates,
                       const std::vector<real> &weights, const Field<real> &reference,
                       const RVec &translation = {0, 0, 0},
                       const Quaternion &orientation = {{1, 0, 0}, 0},
                       const RVec &centerOfRotation = {0, 0, 0});
        virtual real densityDensityPotential(const Field<real> &reference,
                                             const Field<real> &comparant);

    protected:
        std::unique_ptr<DensitySpreader> spreader_;
        Field<real>                     *spread_density_;
};

class densityBasedPotentialForce : public densityBasedPotential,
                                   public PotentialForceEvaluator
{
    public:
        explicit densityBasedPotentialForce(const FiniteGrid &reference,
                                            const int n_threads, const int n_sigma,
                                            const real sigma,
                                            const real sigma_differential);

        std::vector<RVec> force(const std::vector<RVec> &coordinates,
                                const std::vector<real> &weights,
                                const Field<real> &reference,
                                const RVec &translation = {0, 0, 0},
                                const Quaternion &orientation = {{1, 0, 0}, 0},
                                const RVec &centerOfRotation = {0, 0, 0});

    protected:
        virtual void setDensityDifferential(const Field<real> &reference,
                                            const Field<real> &comparant) = 0;
        std::unique_ptr < Field < real>> differential_;
        real              sigma_differential_;
        int               n_threads_;
        std::vector<RVec> force_;
};

}      /* volumedata */
}      /* gmx */

#endif /* end of include guard: GMX_EXTERNALPOTENTIAL_COMMON_H */
