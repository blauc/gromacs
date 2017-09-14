/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2011,2012,2013,2014,2015, by the GROMACS development team, led
 * by
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
 * Implements gmx::analysismodules::map.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "densitymorph.h"

#include <string>

#include "gromacs/externalpotential/modules/densityfitting/emscatteringfactors.h"
#include "gromacs/externalpotential/modules/densityfitting/potentiallibrary.h"
#include "gromacs/externalpotential/modules/densityfitting/potentialprovider.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/fileio/volumedataio.h"
#include "gromacs/math/volumedata/gridreal.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"


namespace gmx
{

namespace analysismodules
{

namespace
{

class DensityMorph : public TrajectoryAnalysisModule
{
    public:
        DensityMorph()  = default;
        ~DensityMorph() = default;

        virtual void initOptions(IOptionsContainer          *options,
                                 TrajectoryAnalysisSettings *settings) override;
        virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
                                  const TopologyInformation        &top) override;
        virtual void initAfterFirstFrame(const TrajectoryAnalysisSettings &settings,
                                         const t_trxframe                 &fr) override;
        virtual void optionsFinished(TrajectoryAnalysisSettings *settings) override;

        virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                                  TrajectoryAnalysisModuleData *pdata) override;

        virtual void finishAnalysis(int nframes) override;
        virtual void writeOutput() override;

    private:

        std::string          fnmobile_;
        std::string          fntarget_;
        std::string          fnforcedensity_;
        std::string          fnmorphtraj_;

        int                  every_ = 10;

        volumedata::GridReal mobile_;
        volumedata::GridReal target_;
        volumedata::GridReal morph_;
        int                  nFr_           = 10;
        std::string          potentialType_;
        std::unique_ptr<volumedata::IStructureDensityPotentialProvider>
        potentialProvider_;
        volumedata::PotentialEvaluatorHandle potentialEvaluator;
};

void DensityMorph::initOptions(IOptionsContainer          *options,
                               TrajectoryAnalysisSettings *settings)
{
    auto        potentialNames = volumedata::PotentialLibrary().available();
    const char *c_potentialTypes[1]; // TODO: this fixed size array is required
                                     // for StringOption enumValue
    for (size_t i = 0; i < potentialNames.size(); i++)
    {
        /* code */
        c_potentialTypes[i] = potentialNames[i].c_str();
    }
    potentialType_ = potentialNames[0];
    static const char *const desc[] = {
        "[THISMODULE] Morphs mobile density into target density."
    };

    settings->setHelpText(desc);

    options->addOption(FileNameOption("from")
                           .filetype(eftVolumeData)
                           .inputFile()
                           .store(&fnmobile_)
                           .defaultBasename("mobile")
                           .description("CCP4 density map to morph."));

    options->addOption(IntegerOption("every")
                           .store(&every_)
                           .description("Write -every intermediate steps."));

    options->addOption(IntegerOption("n-frames")
                           .store(&nFr_)
                           .description("Number of morphing iterations."));

    options->addOption(FileNameOption("to")
                           .filetype(eftVolumeData)
                           .inputFile()
                           .store(&fntarget_)
                           .defaultBasename("mobile")
                           .description("CCP4 density map to morph into."));

    options->addOption(FileNameOption("forces")
                           .filetype(eftVolumeData)
                           .outputFile()
                           .store(&fnforcedensity_)
                           .defaultBasename("forces")
                           .description("CCP4 density maps with forces."));

    options->addOption(FileNameOption("morph")
                           .filetype(eftVolumeData)
                           .outputFile()
                           .store(&fnmorphtraj_)
                           .defaultBasename("morph")
                           .description("CCP4 density maps during morphing."));

    options->addOption(StringOption("potential-type")
                           .store(&potentialType_)
                           .description("potential type.")
                           .enumValue(c_potentialTypes));

}

void DensityMorph::initAnalysis( const TrajectoryAnalysisSettings & /*settings*/, const TopologyInformation        & /*top*/) {}
void DensityMorph::optionsFinished(
        TrajectoryAnalysisSettings * /*settings*/)
{

    volumedata::MrcFile ccp4inputfile;
    ccp4inputfile.read(fnmobile_, mobile_);

    potentialProvider_ = volumedata::PotentialLibrary().create(potentialType_)();


}
void DensityMorph::initAfterFirstFrame( const TrajectoryAnalysisSettings & /*settings*/, const t_trxframe & /*fr*/) { }
void DensityMorph::analyzeFrame(int /*frnr*/, const t_trxframe & /*fr*/, t_pbc * /*pbc*/, TrajectoryAnalysisModuleData * /*pdata*/) {}
void DensityMorph::finishAnalysis(int /*nframes*/) {}
void DensityMorph::writeOutput() {}

}       // namespace

const char DensityMorphInfo::name[]             = "densitymorph";
const char DensityMorphInfo::shortDescription[] =
    "Morph mobile density into target.";

TrajectoryAnalysisModulePointer DensityMorphInfo::create()
{
    return TrajectoryAnalysisModulePointer(new DensityMorph);
}

} // namespace analysismodules

} // namespace gmx
