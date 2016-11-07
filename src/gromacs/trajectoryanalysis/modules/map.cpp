/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2011,2012,2013,2014,2015, by the GROMACS development team, led by
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

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/modules/average.h"
#include "gromacs/analysisdata/modules/histogram.h"
#include "gromacs/analysisdata/modules/plot.h"
#include "gromacs/fileio/volumedataio.h"
#include "gromacs/externalpotential/modules/densityfitting/emscatteringfactors.h"
// #include "gromacs/externalpotential/modules/densityfitting/ifgt/Ifgt.h"
#include "gromacs/math/gausstransform.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/topology.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/math/do_fit.h"
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
        Map();

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
        std::string                           fnmapinput_;
        std::string                           fnmapoutput_;
        std::string                           fncorrelation_;
        std::string                           fnkldivergence_;
        std::string                           fnmapundersample_;
        std::string                           fouriertransform_;

        float                                 sigma_;
        float                                 n_sigma_;
        AnalysisData                          mapdata_;
        volumedata::FastGaussianGridding      gt_;
        volumedata::GridReal                  inputdensity_;
        std::unique_ptr<volumedata::GridReal> outputdensity_;
        volumedata::GridReal                  outputDensityBuffer_;
        real                                  spacing_;
        bool                                  bPrint_;
        std::vector<float>                    weight_;
        int                                   every_;
        real                                  expAverage_;
        real                                  correlationThreshold_;
        FILE                                * correlationFile_;
        FILE                                * kldivergenceFile_;
        real                                  klOffset_;
        int                                   nFr_;
        bool                                  bFit_;
        int                                   nAtomsReference_;
        std::vector<RVec>                     referenceX;
        matrix                                referenceBox;
        std::vector<real>                     fitWeights;
        real                                  inputGridEntropy;
        bool                                  bUseBox_;
};

Map::Map()
    : sigma_(0.4), n_sigma_(5), spacing_(0.2), bPrint_(false), every_(1), expAverage_ {1.}, correlationThreshold_ {
    0.
}, klOffset_ {
    0.
}, nFr_ {
    0
}, bFit_ {
    false
}, bUseBox_ {
    false
}
{
}


void
Map::initOptions(IOptionsContainer *options, TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] = {
        "[THISMODULE] is a tool to read in and write out (electron) density maps.",
        "With this tool you can produce a density map from an input ",
        "coordinate file to be embedded in a .tpr file with grompp for",
        "example.[PAR]",
        "Possible uses are:[PAR]",
        "* Provide an input structure with [TT]-f[tt] and output a density map based on",
        "the grid settings. Note that you can also input a whole trajectory, in that",
        "case, a series of map files will be output. Use and index file to spread just",
        "a part of the atoms[BR]",
        "* Provide a map with [TT]-mi[tt] and and output some map characteristics[BR]",
        "* Provide a [TT].tpr[tt] file with [TT]-s[tt] and output the embedded map with [TT]-mo[tt]"
    };

    settings->setHelpText(desc);
    settings->setFlag(TrajectoryAnalysisSettings::efUseTopX, true);

    options->addOption(FileNameOption("mi").filetype(eftVolumeData).inputFile()
                           .store(&fnmapinput_).defaultBasename("ccp4in")
                           .description("CCP4 density map input file"));
    options->addOption(FileNameOption("moversample").filetype(eftVolumeData).outputFile()
                           .store(&fnmapundersample_).defaultBasename("undersampledccp4")
                           .description("CCP4 density map interpolated file"));
    options->addOption(FileNameOption("mo").filetype(eftVolumeData).outputFile()
                           .store(&fnmapoutput_).defaultBasename("ccp4out")
                           .description("CCP4 density map output file"));
    options->addOption(FloatOption("sigma").store(&sigma_)
                           .description("Create a simulated density by replacing the atoms by Gaussian functions of width sigma (nm)"));
    options->addOption(FloatOption("N_sigma").store(&n_sigma_)
                           .description("How many Gaussian width shall be used for spreading?"));
    options->addOption(FloatOption("spacing").store(&spacing_)
                           .description("Spacing of the density grid (nm)"));
    options->addOption(IntegerOption("every").store(&every_)
                           .description("Analyse only -every frame."));
    options->addOption(FileNameOption("corr").filetype(eftGenericData).outputFile().store(&fncorrelation_)
                           .description("Correlation."));
    options->addOption(FileNameOption("kl").filetype(eftGenericData).outputFile().store(&fnkldivergence_)
                           .description("KL-divergence."));
    options->addOption(BooleanOption("print").store(&bPrint_)
                           .description("Output density information to terminal, then exit."));
    options->addOption(BooleanOption("useBox").store(&bUseBox_)
                           .description("Use the box information in the structure file to setup the grid."));
    options->addOption(BooleanOption("fit").store(&bFit_)
                           .description("Fit to structure."));
    options->addOption(FloatOption("expaverage").store(&expAverage_)
                           .description("Factor for exponential averaging new_map = f*curr+(1-f)*old."));
    options->addOption(FloatOption("corrthreshold").store(&correlationThreshold_)
                           .description("Density threshold when calculating correlations."));
    options->addOption(FloatOption("klOffset").store(&klOffset_)
                           .description("Threshold for KL-Divergence."));
    options->addOption(FileNameOption("fourier").filetype(eftVolumeData).outputFile()
                           .store(&fouriertransform_).description("Fourier Transform."));

}

void
Map::initAnalysis(const TrajectoryAnalysisSettings & /*settings*/, const TopologyInformation &top)
{
    if (top.topology())
    {
        get_pdb_atomnumber(&(top.topology()->atoms), gmx_atomprop_init() );
        for (int i_atom = 0; i_atom < top.topology()->atoms.nr; i_atom++)
        {
            weight_.push_back(externalpotential::atomicNumber2EmScatteringFactor(top.topology()->atoms.atom[i_atom].atomnumber));
        }
    }
    if (top.topology() && bFit_)
    {
        nAtomsReference_ = top.topology()->atoms.nr;
        fitWeights.resize(nAtomsReference_, 1);

        rvec * x = as_rvec_array(referenceX.data());
        top.getTopologyConf(&x, referenceBox);

        reset_x(nAtomsReference_, nullptr, nAtomsReference_, nullptr, x, fitWeights.data());
        referenceX.assign(x, x+nAtomsReference_);
    }
}

void
Map::optionsFinished(TrajectoryAnalysisSettings * /*settings*/)
{
    gt_.set_sigma(sigma_);
    gt_.set_n_sigma(n_sigma_);
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
    if (!fnmapoutput_.empty() || !fncorrelation_.empty() || !fnkldivergence_.empty())
    {
        outputdensity_ = std::unique_ptr<volumedata::GridReal>(new volumedata::GridReal());
    }
    if (!fncorrelation_.empty() || !fnkldivergence_.empty())
    {
        if (fnmapinput_.empty())
        {
            fprintf(stderr, "Please provide map to correlate to with -mi.\n");
            std::exit(0);
        }
    }
    if (!fncorrelation_.empty())
    {
        // set negative values to zero
        inputdensity_.add_offset(klOffset_);
        std::for_each(std::begin(inputdensity_.access().data()), std::end(inputdensity_.access().data()), [](real &value){value = std::max(value, (real)0.); });
        correlationFile_ = fopen(fncorrelation_.c_str(), "w");
    }

    if (!fnkldivergence_.empty())
    {
        inputdensity_.add_offset(klOffset_);
        inputdensity_.normalize();
        kldivergenceFile_ = fopen(fnkldivergence_.c_str(), "w");
        inputGridEntropy  = volumedata::GridMeasures(inputdensity_).entropy();
    }
    if (!fouriertransform_.empty())
    {
        auto result {
            volumedata::FourierTransformRealToComplex3D(inputdensity_).transform()
        };
        fprintf(stderr, " \n");
        // fprintf(stderr,"%d ", result->numGridPointsXY();
        for (int i = 0; i < 1; ++i)
        {

            for (auto x = result->access().zy_column(0, 0)[0]; x != result->access().zy_column(0, 0)[1]; ++x)
            {
                fprintf(stderr, "%g\n", std::sqrt((x)->re*(x)->re + (x)->im*(x)->im));
            }
            fprintf(stderr, "\n\n");
        }
    }
}

void
Map::set_box_from_frame(const t_trxframe &fr, matrix box, rvec translation)
{
    if ((det(fr.box) > 1e-6) && bUseBox_)
    {
        copy_mat(fr.box, box);
        clear_rvec(translation); // TODO: more user options for setting translation
        return;
    }

    fprintf(stderr, "Did not find suitable box for atom in structure file, guessing from structure extend.\n");
    clear_mat(box);
    std::vector<RVec> x_RVec(fr.natoms);
    x_RVec.assign(fr.x, fr.x+fr.natoms);
    for (int i = XX; i <= ZZ; i++)
    {
        auto compareIthComponent =  [ i ](RVec a, RVec b){
                return a[i] < b[i];
            };
        auto minMaxX             = minmax_element(x_RVec.begin(), x_RVec.end(), compareIthComponent);
        box[i][i]      = 2*n_sigma_*sigma_ + (*minMaxX.second)[i] - (*minMaxX.first)[i];
        translation[i] = -n_sigma_*sigma_ + (*minMaxX.first)[i];
    }

}

void
Map::set_finitegrid_from_box(matrix box, rvec translation)
{
    gmx::volumedata::IVec extend({(int)ceil(box[XX][XX]/spacing_), (int)ceil(box[YY][YY]/spacing_), (int)ceil(box[ZZ][ZZ]/spacing_)});
    outputdensity_->set_extend(extend);
    outputDensityBuffer_.set_extend(extend);
    outputdensity_->set_cell({extend[XX]*spacing_, extend[YY]*spacing_, extend[ZZ]*spacing_}, {90, 90, 90});
    outputDensityBuffer_.set_cell({extend[XX]*spacing_, extend[YY]*spacing_, extend[ZZ]*spacing_}, {90, 90, 90});
    outputdensity_->set_translation({roundf(translation[XX]/spacing_)*spacing_, roundf(translation[YY]/spacing_)*spacing_, roundf(translation[ZZ]/spacing_)*spacing_});
    outputDensityBuffer_.set_translation({roundf(translation[XX]/spacing_)*spacing_, roundf(translation[YY]/spacing_)*spacing_, roundf(translation[ZZ]/spacing_)*spacing_});
}

void
Map::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc * /*pbc*/,
                  TrajectoryAnalysisModuleData * /*pdata*/)
{
    if (frnr % every_ == 0)
    {
        ++nFr_;
        if (bFit_)
        {
            reset_x(nAtomsReference_, nullptr, nAtomsReference_, nullptr, fr.x, fitWeights.data());
            do_fit(nAtomsReference_, fitWeights.data(), as_rvec_array(referenceX.data()), fr.x);
        }

        if (weight_.size() < std::size_t(fr.natoms))
        {
            weight_.assign(fr.natoms, 1);
            fprintf(stderr, "\nSetting atom weights to unity.\n");
        }
        if (!fnmapoutput_.empty() || !fncorrelation_.empty() || !fnkldivergence_.empty())
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
                // If a reference input density is given, copy the grid properties from here
                outputdensity_->copy_grid(inputdensity_);
                outputDensityBuffer_.copy_grid(inputdensity_);
            }
            if (frnr == 0)
            {
                outputDensityBuffer_.zero();
            }
            outputdensity_->zero();

            gt_.set_grid(std::move(outputdensity_));
            for (int i = 0; i < fr.natoms; ++i)
            {
                gt_.transform(fr.x[i], weight_[i]);
            }
            outputdensity_ = std::move(gt_.finish_and_return_grid());
            if (frnr == 0)
            {
                std::copy(std::begin(outputdensity_->access().data()), std::end(outputdensity_->access().data()), std::begin(outputDensityBuffer_.access().data()));
            }
            else
            {
                // auto linearAveraging = [this](const real & current, const real & average){return (average * (nFr_ - 1) + current) / nFr_;};
                auto exponentialAveraging = [this](const real &current, const real &average){
                        return expAverage_ * current + (1-expAverage_) * average;
                    };
                std::transform(std::begin(outputdensity_->access().data()), std::end(outputdensity_->access().data()),
                               std::begin(outputdensity_->access().data()), std::begin(outputDensityBuffer_.access().data()), exponentialAveraging);
            }
        }
        if (!fncorrelation_.empty())
        {
            fprintf(correlationFile_, "%g\n", volumedata::GridMeasures(inputdensity_).correlate(outputDensityBuffer_, correlationThreshold_));
        }
        if (!fnkldivergence_.empty())
        {
            outputdensity_->add_offset(this->klOffset_);
            outputdensity_->normalize();
            fprintf(kldivergenceFile_, "%g\n", volumedata::GridMeasures(inputdensity_).getKLCrossTermSameGrid(*outputdensity_)-inputGridEntropy);
        }
        if (!fnmapoutput_.empty())
        {
            volumedata::MrcFile ccp4outputfile;
            ccp4outputfile.write(fnmapoutput_.substr(0, fnmapoutput_.size()-5) + std::to_string(frnr) + ".ccp4", outputDensityBuffer_);
            if (frnr == 0)
            {
                ccp4outputfile.write(fnmapoutput_.c_str(), outputDensityBuffer_);
            }
        }
        if (!fnmapundersample_.empty())
        {
            volumedata::MrcFile    ccp4undersample;
            volumedata::FiniteGrid undersampleGrid;
            undersampleGrid.copy_grid(outputDensityBuffer_);
            undersampleGrid.set_extend({2*outputDensityBuffer_.extend()[XX], 2*outputDensityBuffer_.extend()[YY], 2*outputDensityBuffer_.extend()[ZZ]});
            undersampleGrid.set_cell({2*outputDensityBuffer_.cell_lengths()[XX], 2*outputDensityBuffer_.cell_lengths()[YY], 2*outputDensityBuffer_.cell_lengths()[ZZ]}, {90, 90, 90});
            volumedata::GridInterpolator interpolator(undersampleGrid);
            ccp4undersample.write(fnmapundersample_, interpolator.interpolateLinearly(outputDensityBuffer_));
        }
    }
}


void
Map::finishAnalysis(int /*nframes*/)
{
}


void
Map::writeOutput()
{
    if (!fncorrelation_.empty())
    {
        fclose(correlationFile_);
    }
}

}       // namespace

const char MapInfo::name[]             = "map";
const char MapInfo::shortDescription[] =
    "Spread atoms on grid";

TrajectoryAnalysisModulePointer MapInfo::create()
{
    return TrajectoryAnalysisModulePointer(new Map);
}

} // namespace analysismodules

} // namespace gmx
