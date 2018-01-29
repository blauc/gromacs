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
#ifndef GMX_EXTERNALPOTENTIAL_POTENTIALPROVIDER_H
#define GMX_EXTERNALPOTENTIAL_POTENTIALPROVIDER_H

#include "gromacs/math/griddata/griddata.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

#include <memory>
#include <string>
#include <vector>

namespace gmx
{
class WholeMoleculeGroup;

class IPotentialEvaluator
{
    public:
        virtual ~IPotentialEvaluator() = default;
        virtual real potential(const GridDataReal3D &reference) const = 0;
};

class PotentialEvaluatorHandle
{
    public:
        PotentialEvaluatorHandle() = default;
        PotentialEvaluatorHandle(const IPotentialEvaluator *evaluator)
            : evaluator_ {evaluator}
        {}
        real potential(const GridDataReal3D &reference) const
        {
            return evaluator_->potential(reference);
        };

    private:
        const IPotentialEvaluator *evaluator_;
};

class IStructureDensityPotentialProvider
{
    public:
        virtual ~IStructureDensityPotentialProvider() = default;
        virtual void setCoordinates(const std::vector<RVec> &coordinates,
                                    const std::vector<real> &weights) = 0;
        virtual PotentialEvaluatorHandle planPotential(
            const std::vector<RVec> &coordinates, const std::vector<real> &weights,
            const GridDataReal3D &reference, const std::string &options) = 0;
};

class IForceEvaluator
{
    public:
        virtual ~IForceEvaluator() = default;
        virtual void force(std::vector<RVec>        &force,
                           const std::vector<RVec>  &coordinates_,
                           const std::vector<float> &weights_,
                           const GridDataReal3D     &reference) const = 0;
};

class ForceEvaluatorHandle
{
    public:
        ForceEvaluatorHandle() = default;
        ForceEvaluatorHandle(const IForceEvaluator *evaluator)
            : evaluator_ {evaluator}
        {};
        void force(std::vector<RVec> &force, const std::vector<RVec> &coordinates,
                   const std::vector<float> &weights,
                   const GridDataReal3D &reference) const
        {
            evaluator_->force(force, coordinates, weights, reference);
        };

    private:
        const IForceEvaluator *evaluator_;
};

class IStructureDensityForceProvider
{
    public:
        virtual ~IStructureDensityForceProvider() = default;
        virtual void setCoordinates(const std::vector<RVec> &coordinates,
                                    const std::vector<real> &weights) = 0;
        virtual ForceEvaluatorHandle planForce(const GridDataReal3D &reference,
                                               const std::string    &options) = 0;
};

class IStructureDensityPotentialForceProvider
    : public IStructureDensityForceProvider,
      public IStructureDensityPotentialProvider
{
    public:
        virtual void setCoordinates(const std::vector<RVec> &coordinates,
                                    const std::vector<real> &weights) override;

        virtual ~IStructureDensityPotentialForceProvider() = default;
        virtual std::pair<PotentialEvaluatorHandle, ForceEvaluatorHandle> plan(const GridDataReal3D &reference,
                                                                               const std::string    &options) = 0;
};

}      /* gmx */
#endif /* end of include guard: GMX_EXTERNALPOTENTIAL_POTENTIALPROVIDER_H */
