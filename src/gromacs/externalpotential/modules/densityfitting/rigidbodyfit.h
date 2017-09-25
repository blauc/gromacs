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
#ifndef GMX_EXTERNALPOTENTIAL_RIGIDBODYFIT_H
#define GMX_EXTERNALPOTENTIAL_RIGIDBODYFIT_H

#include "gmxpre.h"

#include <vector>

#include "gromacs/math/quaternion.h"
#include "gromacs/math/vectypes.h"

namespace gmx
{

template <class T> class Field;
class PotentialEvaluatorHandle;
class PotentialForceEvaluator;

class RigidBodyFitResult
{
    public:
        RigidBodyFitResult(const RVec &translation, const RVec &center_of_rotation,
                           const Quaternion &orientation,
                           const real bestFitPotential);
        RVec translation() const;
        Quaternion orientation() const;
        RVec centerOfRotation() const;
        real potential() const;

    private:
        const RVec       translation_;
        const RVec       center_of_rotation_;
        const Quaternion orientation_;
        const real       potential_;
};

class RigidBodyFit
{
    public:
        RigidBodyFitResult
        fitCoordinates(const Field<real> &reference, const std::vector<RVec> &x,
                       const std::vector<real> &weights,
                       const PotentialForceEvaluator &fitPotentialProvider);

        RigidBodyFitResult
        fitCoordinates(const Field<real> &reference, const std::vector<RVec> &x,
                       const std::vector<real> &weights,
                       const PotentialEvaluatorHandle &fitPotentialProvider);

    private:
        const real minimial_improvement_ = 1e-10;
        const int  max_steps_            = 1e5;
        RVec
        gradientTranslation_(const Field<real> &reference, const std::vector<RVec> &x,
                             const std::vector<real> &weights,
                             const PotentialEvaluatorHandle &fitPotentialProvider,
                             const RVec &translation,
                             const RVec &center_of_rotation,
                             const Quaternion &orientation);
        Quaternion
        gradientOrientation_(const Field<real> &reference, const std::vector<RVec> &x,
                             const std::vector<real> &weights,
                             const PotentialEvaluatorHandle &fitPotentialProvider,
                             const RVec &translation,
                             const RVec &center_of_rotation,
                             const Quaternion &orientation);
};
}      // namespace gmx
#endif /* end of include guard: GMX_EXTERNALPOTENTIAL_RIGIDBODYFIT_H */
