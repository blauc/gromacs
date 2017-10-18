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

#include "rigidbodyfit.h"

#include <string>

#include "gromacs/externalpotential/modules/densityfitting/emscatteringfactors.h"
#include "gromacs/externalpotential/modules/densityfitting/potentiallibrary.h"
#include "gromacs/externalpotential/modules/densityfitting/potentialprovider.h"
#include "gromacs/externalpotential/modules/densityfitting/rigidbodyfit.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/fileio/griddataio.h"

#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"

#include <algorithm>

namespace gmx
{

namespace analysismodules
{

namespace
{

class RigidBodyFit : public TrajectoryAnalysisModule
{
    public:
        RigidBodyFit()  = default;
        ~RigidBodyFit() = default;

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
        void frameToForceDensity_(const t_trxframe &fr);

        std::string             fnmapinput_;
        std::string             fnpotential_ = std::string("potential.dat");
        std::string             forcedensity_;
        std::string             fnoptions_;
        std::string             fntrajoutput_ = std::string("fitted.xtc");
        std::string             optionsstring_;

        std::unique_ptr < Field < real>> inputdensity_;
        bool                    bRigidBodyFit_ = true;
        std::vector<float>      weight_;
        int                     every_ = 1;
        FILE                   *potentialFile_;
        int                     nFr_ = 0;
        std::string             potentialType_;
        std::unique_ptr<IStructureDensityPotentialProvider>
        potentialProvider_;
        PotentialEvaluatorHandle potentialEvaluator;
};

void RigidBodyFit::initOptions(IOptionsContainer          *options,
                               TrajectoryAnalysisSettings *settings)
{
    auto        potentialNames = PotentialLibrary().available();
    const char *c_potentialTypes[1]; // TODO: this fixed size array is required
                                     // for StringOption enumValue
    for (size_t i = 0; i < potentialNames.size(); i++)
    {
        /* code */
        c_potentialTypes[i] = potentialNames[i].c_str();
    }
    potentialType_ = potentialNames[0];
    static const char *const desc[] = {
        "[THISMODULE] Fits trajectory to a three-dimensional density."
    };

    settings->setHelpText(desc);
    settings->setFlag(TrajectoryAnalysisSettings::efUseTopX, true);

    options->addOption(FileNameOption("fo").filetype(eftTrajectory).outputFile().store(&fntrajoutput_).description("The trajectory fitted to the input density."));
    options->addOption(FileNameOption("mi")
                           .filetype(eftgriddata)
                           .inputFile()
                           .store(&fnmapinput_)
                           .defaultBasename("ccp4in")
                           .description("CCP4 density map to compare to"));
    options->addOption(IntegerOption("every").store(&every_).description(
                               "Analyse only -every frame."));
    options->addOption(FileNameOption("options")
                           .filetype(eftGenericData)
                           .inputFile()
                           .store(&fnoptions_)
                           .description("Options for the potential."));
    options->addOption(BooleanOption("rigidBodyFit")
                           .store(&bRigidBodyFit_)
                           .description("Use rigid body fitting in all steps."));
    options->addOption(StringOption("potential-type")
                           .store(&potentialType_)
                           .description("potential type.")
                           .enumValue(c_potentialTypes));
}

void RigidBodyFit::initAnalysis(
        const TrajectoryAnalysisSettings & /*settings*/,
        const TopologyInformation       &top)
{

    if (top.topology())
    {
        auto atomprop = gmx_atomprop_init();
        get_pdb_atomnumber(&(top.topology()->atoms), atomprop);
        for (int i_atom = 0; i_atom < top.topology()->atoms.nr; i_atom++)
        {
            weight_.push_back(atomicNumber2EmScatteringFactor(
                                      top.topology()->atoms.atom[i_atom].atomnumber));
        }
        gmx_atomprop_destroy(atomprop);
    }
}

void RigidBodyFit::optionsFinished(
        TrajectoryAnalysisSettings * /*settings*/)
{

    inputdensity_ = std::unique_ptr < Field < real>>(new Field<real>(MrcFile().read(fnmapinput_)));

    // set negative values to zero
    std::for_each(std::begin(*inputdensity_), std::end(*inputdensity_), [](real &value) { value = std::max(value, (real)0.); });
    // potentialProvider_ = PotentialLibrary().create(potentialType_)();
    // fprintf(potentialFile_, "time        %s", potentialType_.c_str());
    //
    // auto optionsFile = fopen(fnoptions_.c_str(), "r");
    // if (!fnoptions_.empty() && optionsFile)
    // {
    //
    //     fseek(optionsFile, 0, SEEK_END);
    //     auto  length = ftell(optionsFile);
    //     fseek(optionsFile, 0, SEEK_SET);
    //     char *buffer = (char *)malloc(length + 1);
    //     if (buffer)
    //     {
    //         fread(buffer, 1, length, optionsFile);
    //     }
    //     buffer[length] = '\0';
    //     optionsstring_ = buffer;
    // }
}

void RigidBodyFit::initAfterFirstFrame(
        const TrajectoryAnalysisSettings & /*settings*/, const t_trxframe &fr)
{
    if (weight_.size() < std::size_t(fr.natoms))
    {
        weight_.assign(fr.natoms, 1);
        fprintf(stderr, "\nDid not find topology information to derive atom "
                "scattering weights. Setting atom weights to unity.\n");
    }

    std::vector<RVec> rVecCoordinates(fr.x, fr.x + fr.natoms);
    potentialEvaluator = potentialProvider_->planPotential(
                rVecCoordinates, weight_, *inputdensity_, optionsstring_);
}

void RigidBodyFit::analyzeFrame(int /*frnr*/, const t_trxframe &fr,
                                t_pbc * /*pbc*/,
                                TrajectoryAnalysisModuleData * /*pdata*/)
{
    std::vector<RVec>        rVecCoordinates(fr.x, fr.x + fr.natoms);

    // RigidBodyFit             rigidbodyfit;
    // rigidbodyfit.fitCoordinates(inputdensity_, rVecCoordinates, weight_,
    // potentialEvaluator);
}

void RigidBodyFit::finishAnalysis(int /*nframes*/) {}

void RigidBodyFit::writeOutput()
{
}

}       // namespace

const char RigidBodyFitInfo::name[]             = "rigidbodyfit";
const char RigidBodyFitInfo::shortDescription[] =
    "Fit a structure or trajectory to an input density.";

TrajectoryAnalysisModulePointer RigidBodyFitInfo::create()
{
    return TrajectoryAnalysisModulePointer(new RigidBodyFit);
}

} // namespace analysismodules

} // namespace gmx
