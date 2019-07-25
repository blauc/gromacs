/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Tests for functionality of the density fitting module.
 *
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "gromacs/applied_forces/densityfittingforceprovider.h"

#include <numeric>

#include <gtest/gtest.h>

#include "gromacs/applied_forces/densityfittingparameters.h"
#include "gromacs/domdec/localatomsetmanager.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/coordinatetransformation.h"
#include "gromacs/math/gausstransform.h"
#include "gromacs/math/multidimarray.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/mdatom.h"

namespace gmx
{

namespace
{

class DensityFittingForceProviderTest : public ::testing::Test
{
    public:
        DensityFittingForceProviderTest()
        {
            referenceDensityGenerator.add({{5, 5, 5}, 1});
        }

        ~DensityFittingForceProviderTest()
        {
            done_commrec(cr_);
        }

        ForceProviderInput generateForceProviderInput()
        {
            double       t = 0;
            t_mdatoms    mdatoms;
            const matrix box {};
            return {coordinates_, mdatoms, t, box, *cr_};
        }

        ForceWithVirial generateForceWithVirial()
        {
            std::vector<RVec> force = std::vector<RVec>(coordinates_.size(), RVec(0, 0, 0));
            return {force, false};
        }

        DensityFittingParameters generateParameters(std::vector<int> atomIndices)
        {
            const auto        densityFittingAtoms  = atomSets.add(atomIndices);
            std::vector<real> amplitudes(coordinates_.size(), 1.);
            TranslateAndScale intoLatticeTransform = {scale, translation};
            return {
                       densityFittingAtoms, intoLatticeTransform, forceConstant_,
                       sigma, nSigma, referenceDensityGenerator.view(),
                       DensityFittingAmplitudeMethod::Unity,
                       DensitySimilarityMeasureMethod::innerProduct,
                       everyNSteps_, true, 10
            };
        }

        gmx_enerdata_t generateEnergdata()
        {return {F_NRE, 0}; }
    protected:
        t_commrec *cr_ = init_commrec();
        LocalAtomSetManager                   atomSets;
        real                                  sigma                     = 1;
        float                                 nSigma                    = 5;
        dynamicExtents3D                      ex                        = {11, 11, 11};
        GaussianSpreadKernelParameters::Shape gaussShape                = {{sigma, sigma, sigma}, nSigma};
        GaussTransform3D                      referenceDensityGenerator = { ex, gaussShape};
        RVec                                  scale                     = {1, 1, 1};
        RVec                                  translation               = {0, 0, 0};
        std::vector<RVec>                     coordinates_              = {{0, 0, 0}, {5, 5, 5}, {4, 4, 4}, {6, 4, 4}, {4, 6, 4}, {4, 4, 6}, {6, 6, 6}, {10, 10, 10}};
        real                                  forceConstant_            = 1e8;
        int everyNSteps_ = 1;
};

TEST_F(DensityFittingForceProviderTest, CanConstruct)
{
    std::vector<int>            atomIndicesToFit = {0};
    auto                        parameters       = generateParameters(atomIndicesToFit);
    DensityFittingForceProvider densfitForces(parameters);
}

TEST_F(DensityFittingForceProviderTest, Forces)
{
    std::vector<int>            atomIndicesToFit(coordinates_.size());
    std::iota(std::begin(atomIndicesToFit), std::end(atomIndicesToFit), 0);
    auto                        parameters       = generateParameters(atomIndicesToFit);
    DensityFittingForceProvider densfitForces(parameters);

    auto                        forceProviderInput = generateForceProviderInput();
    auto                        forceWithVirial    = generateForceWithVirial();
    auto                        enerd              = generateEnergdata();
    ForceProviderOutput         output             = { &forceWithVirial, &enerd};
    densfitForces.calculateForces(forceProviderInput, &output);

    for (const auto fitAtomIndex : atomIndicesToFit)
    {
        fprintf(stderr, "\n %d  = [%10.10g %10.10g %10.10g] : [%10.10g %10.10g %10.10g] \n",
                fitAtomIndex,
                coordinates_[fitAtomIndex][XX],
                coordinates_[fitAtomIndex][YY],
                coordinates_[fitAtomIndex][ZZ],
                output.forceWithVirial_.force_[fitAtomIndex][XX],
                output.forceWithVirial_.force_[fitAtomIndex][YY],
                output.forceWithVirial_.force_[fitAtomIndex][ZZ]);
    }
}

} // namespace

} // namespace gmx
