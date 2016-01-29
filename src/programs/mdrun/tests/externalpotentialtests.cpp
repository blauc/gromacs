/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
 * Test for external potential modularity.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "config.h"

#include <gtest/gtest.h>

#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textwriter.h"
#include "energyreader.h"
#include "moduletest.h"

namespace
{

//! Test fixture for exernal potential modules
class ExternalPotentialTemplateTest :
    public gmx::test::MdrunTestFixture,
    public ::testing::WithParamInterface<const char *>
{
    public:
        //! The file name of the MDP file
        std::string                     theMdpFile;

        std::string                     ExternalPotentialInputData;
        std::string                     theExternalPotentialInputFile  = fileManager_.getTemporaryFilePath("argon5832-input.dat");
        std::string                     theExternalPotentialOutputFile = fileManager_.getTemporaryFilePath("argon5832-output.dat");
        //! Execute the external potential module test
        gmx::test::EnergyFrameReaderPtr energy;
        void runTest()
        {
            runner_.useStringAsMdpFile(theMdpFile);
            runner_.useTopGroAndNdxFromDatabase("argon5832");
            EXPECT_EQ(0, runner_.callGrompp());
            gmx::TextWriter::writeFileFromString(theExternalPotentialInputFile, ExternalPotentialInputData);
            runner_.edrFileName_ = fileManager_.getTemporaryFilePath("argon5832-ener.edr");
            ASSERT_EQ(0, runner_.callMdrun());
            // energy = gmx::test::openEnergyFileToReadFields(runner_.edrFileName_, {{"Potential"}});
        }
};

//! Helper typedef for naming test cases like sentences
typedef ExternalPotentialTemplateTest ExternalPotential;

/* Ensure grompp and mdrun run */
TEST_F(ExternalPotential, CanRun)
{
    theMdpFile = gmx::formatString("integrator              = steep\n"
                                   "nsteps                  = 10\n"
                                   "external-potential      = yes\n"
                                   "external-potential-path = \n"
                                   "template-input          = ") + theExternalPotentialInputFile +
        gmx::formatString("\ntemplate-output         = ") + theExternalPotentialOutputFile +
        gmx::formatString("\ntemplate-groups         = first_half second_half\n");

    ExternalPotentialInputData = std::string("100");

    runTest();
}


/* This test ensures mdrun can write various quantities at various frequencies */
/*TEST_F(ExternalPotential, PullsAtomsToReferencePoints)
   {
    theMdpFile = gmx::formatString("integrator              = steep\n"
                                   "nsteps                  = 100\n"
                                   "external-potential      = yes\n"
                                   "external-potential-path = \n"
                                   "template-input          = ") + theExternalPotentialInputFile +
        gmx::formatString("\ntemplate-output         = ") + theExternalPotentialOutputFile +
        gmx::formatString("\ntemplate-groups         = first_half second_half\n");

    ExternalPotentialInputData = std::string("100");

    runTest();
   }*/

#ifdef __INTEL_COMPILER
#pragma warning( disable : 177 )
#endif

// TEST_F(ExternalPotential,NoFatalErrorFrom);

// TEST_F(ExternalPotential,EnergyReproducedFrom);

} // namespace
