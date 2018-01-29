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
#ifndef GMX_EXTERNALPOTENTIAL_KULLBACKLEIBLER_H
#define GMX_EXTERNALPOTENTIAL_KULLBACKLEIBLER_H

#include "gmxpre.h"

#include "../densityspreader.h"
#include "../potentialprovider.h"
#include "densitybasedpotential.h"
#include "gromacs/math/griddata/griddata.h"

#include "gromacs/math/quaternion.h"

#include <functional>
#include <memory>
#include <string>
namespace gmx
{
class WholeMoleculeGroup;

class KullbackLeiblerPotential : public DensityBasedPotential
{
    public:
        KullbackLeiblerPotential(const DensitySpreader &spreader);
        ~KullbackLeiblerPotential() = default;

        real densityDensityPotential(const GridDataReal3D &reference,
                                     const GridDataReal3D &comparant) const override;
};

class KullbackLeiblerForce : public DensityBasedForce
{
    public:
        KullbackLeiblerForce(const DensitySpreader &spreader, real sigma_differential, int n_threads);
        ~KullbackLeiblerForce() = default;

        const GridDataReal3D &densityDifferential(const GridDataReal3D    &reference, const GridDataReal3D    &comparant) const override;
};

class KullbackLeiblerProvider : public IStructureDensityPotentialForceProvider
{
    public:
        ~KullbackLeiblerProvider() = default;
        ForceEvaluatorHandle
        planForce(const GridDataReal3D &reference, const std::string &options) override;
        PotentialEvaluatorHandle
        planPotential(const std::vector<RVec> &coordinates,
                      const std::vector<real> &weights, const GridDataReal3D &reference,
                      const std::string &options) override;
        void setCoordinates(const std::vector<RVec> &coordinates,
                            const std::vector<real> &weights) override;

    private:
        void log_(const std::string &message);
        void parseOptions_(const std::string &options);
        real sigma_;
        int  n_threads_;
        int  n_sigma_;
        std::unique_ptr<IForceEvaluator>     force_evaluator_;
        std::unique_ptr<IPotentialEvaluator> potential_evaluator_;
        std::unique_ptr<DensitySpreader>     spreader_;

};
/****************************INFO Classes**************************************/

class KullbackLeiblerPotentialInfo
{
    public:
        static std::string name;
        static std::unique_ptr<IStructureDensityPotentialForceProvider> create();
};

}      /* gmx */

#endif /* end of include guard: GMX_EXTERNALPOTENTIAL_KULLBACKLEIBLER_H */
