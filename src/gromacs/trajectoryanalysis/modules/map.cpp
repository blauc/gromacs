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

#include "map.h"

#include <algorithm>
#include <numeric>
#include <string>
#include <vector>
#include <iterator>

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/modules/average.h"
#include "gromacs/analysisdata/modules/histogram.h"
#include "gromacs/analysisdata/modules/plot.h"
#include "gromacs/externalpotential/forceplotter.h"
#include "gromacs/externalpotential/modules/densityfitting/emscatteringfactors.h"
#include "gromacs/externalpotential/modules/densityfitting/forcedensity.h"
#include "gromacs/externalpotential/modules/densityfitting/provider.h"
#include "gromacs/externalpotential/modules/densityfitting/densityspreader.h"
#include "gromacs/externalpotential/modules/densityfitting/rigidbodyfit.h"
#include "gromacs/externalpotential/modules/densityfitting/potentialprovider.h"

#include "gromacs/math/quaternion.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/fileio/volumedataio.h"
#include "gromacs/math/do_fit.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/math/volumedata/convolution.h"
#include "gromacs/math/volumedata/densitypadding.h"
#include "gromacs/math/volumedata/gausstransform.h"
#include "gromacs/math/volumedata/gridinterpolator.h"
#include "gromacs/math/volumedata/gridmeasures.h"
#include "gromacs/math/volumedata/improvedfastgausstransform.h"
#include "gromacs/math/volumedata/fouriershellcorrelation.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/timing/cyclecounter.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/math/vectypes.h"

namespace gmx
{

namespace analysismodules
{

namespace
{

class Map : public TrajectoryAnalysisModule
{
    public:
        Map() = default;

        virtual void initOptions(IOptionsContainer          *options,
                                 TrajectoryAnalysisSettings *settings);
        virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
                                  const TopologyInformation        &top);
        virtual void optionsFinished(TrajectoryAnalysisSettings *settings);

        virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                                  TrajectoryAnalysisModuleData *pdata);

        virtual void finishAnalysis(int nframes);
        void set_finitegrid_from_box(matrix box, rvec translation);
        void set_box_from_frame(const t_trxframe &fr, matrix box, rvec translation);
        virtual void writeOutput();

    private:
        void frameToDensity_(const t_trxframe &fr, int nFr);
        void frameToForceDensity_(const t_trxframe &fr);
        void frameToPotentials_(const t_trxframe &fr, int nFr);
        void openPotentialFileAndPrintHeader_();
        void openFscFileAndPrintHeader_();

        std::string                                                     fnmapinput_;
        std::string                                                     fnmapoutput_;
        std::string                                                     fnfrmapoutput_;
        std::string                                                     fnpotential_;
        std::string                                                     forcedensity_;

        float                                                           sigma_   = 0.4;
        float                                                           n_sigma_ = 5;
        AnalysisData                                                    mapdata_;
        volumedata::GridReal                                            inputdensity_;
        volumedata::GridReal                                            outputdensity_;
        volumedata::GridReal                                            outputDensityBuffer_;
        real                                                            spacing_        = 0.2;
        bool                                                            bPrint_         = false;
        bool                                                            bRigidBodyFit_  = true;
        std::vector<float>                                              weight_;
        int                                                             every_                = 1;
        real                                                            expAverage_           = 1.;
        real                                                            correlationThreshold_ = 0.;
        FILE                                                           *potentialFile_;
        int                                                             nFr_                          = 0;
        bool                                                            bFitFramesToTopStructure_     = false;
        int                                                             nAtomsReference_;
        std::vector<RVec>                                               referenceX_;
        matrix                                                          referenceBox;
        std::vector<real>                                               fitWeights;
        bool                                                            bUseBox_       = false;
        std::string                                                     potentialType_ = "kullbackleibler";
        volumedata::FourierShellCorrelation                             fsc_;
        std::unique_ptr<volumedata::DensitySpreader>                    spreader_;
        FILE                                                           *fscFile_;
        std::unique_ptr<volumedata::IStructureDensityPotentialProvider> potentialProvider_;
        std::unique_ptr<volumedata::PotentialForceEvaluator>            potentialForceEvaluator_;
};

void Map::initOptions(IOptionsContainer          *options,
                      TrajectoryAnalysisSettings *settings)
{
    auto          potentialNames = volumedata::PotentialLibrary<volumedata::IStructureDensityPotentialProvider>().available();
    const char *  c_potentialTypes[4]; // TODO: this fixed size array is required for StringOption enumValue
    for (size_t i = 0; i < potentialNames.size(); i++)
    {
        /* code */
        c_potentialTypes[i] = potentialNames[i].c_str();
    }
    static const char *const desc[] = {
        "[THISMODULE] is a tool to read in and write out (electron) density "
        "maps.",
        "With this tool you can produce a density map from an input ",
        "coordinate file to be embedded in a .tpr file with grompp for",
        "example.[PAR]", "Possible uses are:[PAR]",
        "* Provide an input structure with [TT]-f[tt] and output a density map "
        "based on",
        "the grid settings. Note that you can also input a whole trajectory, in "
        "that",
        "case, a series of map files will be output. Use and index file to "
        "spread just",
        "a part of the atoms[BR]", "* Provide a map with [TT]-mi[tt] and and "
        "output some map characteristics[BR]",
        "* Provide a [TT].tpr[tt] file with [TT]-s[tt] and output the embedded "
        "map with [TT]-mo[tt]"
    };

    settings->setHelpText(desc);
    settings->setFlag(TrajectoryAnalysisSettings::efUseTopX, true);

    options->addOption(FileNameOption("mi")
                           .filetype(eftVolumeData)
                           .inputFile()
                           .store(&fnmapinput_)
                           .defaultBasename("ccp4in")
                           .description("CCP4 density map input file"));
    options->addOption(FileNameOption("mo")
                           .filetype(eftVolumeData)
                           .outputFile()
                           .store(&fnmapoutput_)
                           .defaultBasename("ccp4out")
                           .description("CCP4 density map output file"));
    options->addOption(FileNameOption("mofr")
                           .filetype(eftVolumeData)
                           .outputFile()
                           .store(&fnfrmapoutput_)
                           .defaultBasename("ccp4out")
                           .description("CCP4 density map output file per frame"));
    options->addOption(FloatOption("sigma").store(&sigma_).description(
                               "Create a simulated density by replacing the atoms by Gaussian functions "
                               "of width sigma (nm)"));
    options->addOption(FloatOption("N_sigma").store(&n_sigma_).description(
                               "How many Gaussian width shall be used for spreading?"));
    options->addOption(FloatOption("spacing").store(&spacing_).description(
                               "Spacing of the density grid (nm)"));
    options->addOption(IntegerOption("every").store(&every_).description(
                               "Analyse only -every frame."));
    options->addOption(FileNameOption("potential")
                           .filetype(eftGenericData)
                           .outputFile()
                           .store(&fnpotential_)
                           .description("Calculate potential."));
    options->addOption(BooleanOption("print").store(&bPrint_).description(
                               "Output density information to terminal, then exit."));
    options->addOption(BooleanOption("rigidBodyFit").store(&bRigidBodyFit_).description(
                               "Use rigid body fitting in all steps."));
    options->addOption(BooleanOption("useBox").store(&bUseBox_).description(
                               "Use the box information in the structure file to setup the grid."));
    options->addOption(
            BooleanOption("fit").store(&bFitFramesToTopStructure_).description("Fit to structure."));
    options->addOption(
            FloatOption("expaverage")
                .store(&expAverage_)
                .description(
                    "Factor for exponential averaging new_map = f*curr+(1-f)*old."));
    options->addOption(
            FloatOption("corrthreshold")
                .store(&correlationThreshold_)
                .description("Density threshold when calculating correlations."));
    options->addOption(FileNameOption("forcedensity")
                           .filetype(eftVolumeData)
                           .outputFile()
                           .store(&forcedensity_)
                           .description("Force density."));
    options->addOption(StringOption("potential-type").store(&potentialType_)
                           .description("potential type.").enumValue(c_potentialTypes));
}

void Map::initAnalysis(const TrajectoryAnalysisSettings & /*settings*/,
                       const TopologyInformation       &top)
{
    if (top.topology())
    {
        auto atomprop = gmx_atomprop_init();
        get_pdb_atomnumber(&(top.topology()->atoms), atomprop);
        for (int i_atom = 0; i_atom < top.topology()->atoms.nr; i_atom++)
        {
            weight_.push_back(externalpotential::atomicNumber2EmScatteringFactor(
                                      top.topology()->atoms.atom[i_atom].atomnumber));
        }
        gmx_atomprop_destroy(atomprop);
    }
    if (bFitFramesToTopStructure_)
    {
        if (top.topology())
        {

            nAtomsReference_ = top.topology()->atoms.nr;
            fitWeights.resize(nAtomsReference_, 1);

            rvec *x = as_rvec_array(referenceX_.data());
            top.getTopologyConf(&x, referenceBox);

            reset_x(nAtomsReference_, nullptr, nAtomsReference_, nullptr, x,
                    fitWeights.data());
            referenceX_.assign(x, x + nAtomsReference_);
        }
        else
        {
            GMX_THROW(InvalidInputError("Please provide a topology to fit to."));
        }
    }
}

void Map::optionsFinished(TrajectoryAnalysisSettings * /*settings*/)
{

    if (!fnmapinput_.empty())
    {
        volumedata::MrcFile ccp4inputfile;
        ccp4inputfile.read(fnmapinput_, inputdensity_);
        if (bPrint_)
        {
            fprintf(stderr, "\n%s\n", ccp4inputfile.print_to_string().c_str());
            fprintf(stderr, "\n%s\n", inputdensity_.print().c_str());
        }
    }
    if (!fnmapoutput_.empty() || !fnpotential_.empty() || !forcedensity_.empty())
    {
        spreader_ =
            std::unique_ptr<volumedata::DensitySpreader>(new volumedata::DensitySpreader(inputdensity_, 1, n_sigma_, sigma_));
    }
    if (!fnpotential_.empty() || !forcedensity_.empty())
    {
        if (fnmapinput_.empty())
        {
            fprintf(stderr, "Please provide map to correlate to with -mi.\n");
            std::exit(0);
        }
    }
    if (!fnpotential_.empty())
    {
        // set negative values to zero
        std::for_each(std::begin(inputdensity_.access().data()), std::end(inputdensity_.access().data()), [](real &value) { value = std::max(value, (real)0.); });
        inputdensity_.normalize();
        potentialProvider_ = volumedata::PotentialLibrary<volumedata::IStructureDensityPotentialProvider>().create(potentialType_)();
        std::string optionsString        = "{\"sigma\":"+std::to_string(sigma_)+",\"n_sigma\":"+std::to_string(n_sigma_)+"}\n";

        openPotentialFileAndPrintHeader_();
        if (potentialType_.compare("fsc") == 0)
        {
            openFscFileAndPrintHeader_();
        }

    }
    if (!forcedensity_.empty())
    {
        inputdensity_.normalize();
    }
}

void Map::set_box_from_frame(const t_trxframe &fr, matrix box,
                             rvec translation)
{
    if ((det(fr.box) > 1e-6) && bUseBox_)
    {
        copy_mat(fr.box, box);
        clear_rvec(translation); // TODO: more user options for setting translation
        return;
    }

    fprintf(stderr, "Did not find suitable box for atom in structure file, "
            "guessing from structure extend.\n");
    clear_mat(box);
    const std::vector<RVec> x_RVec(fr.x, fr.x + fr.natoms);
    for (int i = XX; i <= ZZ; i++)
    {
        auto compareIthComponent = [i](RVec a, RVec b) {
                return a[i] < b[i];
            };
        auto minMaxX             =
            minmax_element(x_RVec.begin(), x_RVec.end(), compareIthComponent);
        box[i][i] =
            2 * n_sigma_ * sigma_ + (*minMaxX.second)[i] - (*minMaxX.first)[i];
        translation[i] = -n_sigma_ * sigma_ + (*minMaxX.first)[i];
    }
}

void Map::set_finitegrid_from_box(matrix box, rvec translation)
{
    gmx::volumedata::IVec extend({(int)ceil(box[XX][XX] / spacing_),
                                  (int)ceil(box[YY][YY] / spacing_),
                                  (int)ceil(box[ZZ][ZZ] / spacing_)});
    outputdensity_.set_extend(extend);
    outputDensityBuffer_.set_extend(extend);
    outputdensity_.set_cell(
            {extend[XX] * spacing_, extend[YY] * spacing_, extend[ZZ] * spacing_},
            {90, 90, 90});
    outputDensityBuffer_.set_cell(
            {extend[XX] * spacing_, extend[YY] * spacing_, extend[ZZ] * spacing_},
            {90, 90, 90});
    outputdensity_.set_translation(
            {roundf(translation[XX] / spacing_) * spacing_,
             roundf(translation[YY] / spacing_) * spacing_,
             roundf(translation[ZZ] / spacing_) * spacing_});
    outputDensityBuffer_.set_translation(
            {roundf(translation[XX] / spacing_) * spacing_,
             roundf(translation[YY] / spacing_) * spacing_,
             roundf(translation[ZZ] / spacing_) * spacing_});
}

void Map::frameToDensity_(const t_trxframe &fr, int nFr)
{
    if (nFr == 1)
    {
        if (fnmapinput_.empty())
        {
            // Guess the extend of the map if no box is given in structure file
            // and no other input density is given.
            matrix box;
            rvec   translation;
            set_box_from_frame(fr, box, translation);
            set_finitegrid_from_box(box, translation);
        }
        else
        {
            // If a reference input density is given, copy the grid properties from
            // here
            outputdensity_.copy_grid(inputdensity_);
            outputDensityBuffer_.copy_grid(inputdensity_);
        }

        outputDensityBuffer_.zero();
    }

    std::vector<RVec> coordinates(fr.x, fr.x+fr.natoms);

    outputdensity_ = *spreader_->spreadLocalAtoms(coordinates, weight_);

    if (nFr_ == 1)
    {
        std::copy(std::begin(outputdensity_.access().data()),
                  std::end(outputdensity_.access().data()),
                  std::begin(outputDensityBuffer_.access().data()));
    }
    else
    {
        // auto linearAveraging = [this](const real & current, const real &
        // average){return (average * (nFr_ - 1) + current) / nFr_;};
        auto exponentialAveraging = [this](const real &current,
                                           const real &average) {
                return expAverage_ * current + (1 - expAverage_) * average;
            };
        std::transform(std::begin(outputdensity_.access().data()),
                       std::end(outputdensity_.access().data()),
                       std::begin(outputdensity_.access().data()),
                       std::begin(outputDensityBuffer_.access().data()),
                       exponentialAveraging);
    }
}

void Map::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc * /*pbc*/,
                       TrajectoryAnalysisModuleData * /*pdata*/)
{
    if (frnr % every_ == 0)
    {
        ++nFr_;
        if (bFitFramesToTopStructure_)
        {
            reset_x(nAtomsReference_, nullptr, nAtomsReference_, nullptr, fr.x,
                    fitWeights.data());
            do_fit(nAtomsReference_, fitWeights.data(),
                   as_rvec_array(referenceX_.data()), fr.x);
        }

        if (nFr_ == 1)
        {
            if (weight_.size() < std::size_t(fr.natoms))
            {
                weight_.assign(fr.natoms, 1);
                fprintf(stderr, "\nSetting atom weights to unity.\n");
            }
        }

        bool requiresDensity = !fnmapoutput_.empty() || !fnpotential_.empty()  || !forcedensity_.empty() || !fnfrmapoutput_.empty();
        if (requiresDensity)
        {
            frameToDensity_(fr, nFr_);
        }
        if (!fnpotential_.empty())
        {
            frameToPotentials_(fr, nFr_);
        }
        if (!fnmapoutput_.empty())
        {
            volumedata::MrcFile().write(fnmapoutput_.substr(0, fnmapoutput_.size() - 5) + ".ccp4", outputdensity_);
        }

        if (!fnfrmapoutput_.empty())
        {
            volumedata::MrcFile().write(fnfrmapoutput_.substr(0, fnmapoutput_.size() - 5) + std::to_string(frnr) + ".ccp4", outputDensityBuffer_);
        }
        if (!forcedensity_.empty())
        {
            frameToForceDensity_(fr);
        }
    }
}

void Map::openPotentialFileAndPrintHeader_()
{
    potentialFile_   = fopen(fnpotential_.c_str(), "w");

    fprintf(potentialFile_, "time        %s", potentialType_.c_str());
}

void Map::openFscFileAndPrintHeader_()
{
    fsc_                   = volumedata::FourierShellCorrelation(inputdensity_);
    fscFile_               = fopen("fsccurve.dat", "w");
    fprintf(fscFile_, "%s", (std::accumulate(std::begin(fsc_.getBinEdges()), std::end(fsc_.getBinEdges()), std::string(""),  [](std::string accumulant, const real &value){
                                                 return accumulant + std::string(" ") + std::to_string(value);
                                             })+ "\n").c_str());
}

void Map::frameToPotentials_(const t_trxframe &fr, int /*nFr*/)
{
    using namespace volumedata;
    real potential;

    fprintf(potentialFile_, "\n%8g ", fr.time);

    std::vector<RVec> RVecCoordinates;
    RVecCoordinates.assign(fr.x, fr.x+fr.natoms);

    if (bRigidBodyFit_)
    {
        RigidBodyFit rigidbodyfit;
        rigidbodyfit.fitCoordinates(inputdensity_, RVecCoordinates, weight_,  *potentialForceEvaluator_);
        potential = rigidbodyfit.bestFitPotential();
    }
    else
    {
        potential = potentialForceEvaluator_->potential(RVecCoordinates, weight_, inputdensity_);
    }

    fprintf(potentialFile_, "%8g ", potential);

    if (potentialType_.compare("fsc") == 0)
    {
        auto fscCurve = fsc_.getFscCurve(outputDensityBuffer_, inputdensity_);
        fprintf(fscFile_, "%s", (std::accumulate(std::begin(fscCurve), std::end(fscCurve), std::string(""),  [](std::string accumulant, const real &value){ return accumulant + std::string(" ") + std::to_string(value); })+ "\n").c_str());
    }

}

void Map::frameToForceDensity_(const t_trxframe &fr)
{
    const std::vector<RVec> coordinates(fr.x, fr.x+fr.natoms);

    auto                    forcePlotterForces = potentialForceEvaluator_->force(coordinates, weight_, inputdensity_);

    auto                    plotter = externalpotential::ForcePlotter();
    plotter.start_plot_forces("forces.bild");
    plotter.plot_forces(fr.x, as_rvec_array(forcePlotterForces.data()),
                        fr.natoms, 1);
    plotter.stop_plot_forces();

}

void Map::finishAnalysis(int /*nframes*/) {}

void Map::writeOutput()
{
    if (!fnpotential_.empty())
    {
        fclose(potentialFile_);
        fclose(fscFile_);
    }
}

}       // namespace

const char MapInfo::name[]             = "map";
const char MapInfo::shortDescription[] = "Spread atoms on grid";

TrajectoryAnalysisModulePointer MapInfo::create()
{
    return TrajectoryAnalysisModulePointer(new Map);
}

} // namespace analysismodules

} // namespace gmx
