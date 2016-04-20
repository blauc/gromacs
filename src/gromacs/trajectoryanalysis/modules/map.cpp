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
 * Implements gmx::analysismodules::Angle.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "map.h"

#include <algorithm>
#include <string>
#include <vector>

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/modules/average.h"
#include "gromacs/analysisdata/modules/histogram.h"
#include "gromacs/analysisdata/modules/plot.h"
#include "gromacs/fileio/volumedataio.h"
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
        float                                 sigma_;
        float                                 n_sigma_;
        AnalysisData                          mapdata_;
        volumedata::FastGaussianGridding      gt_;
        volumedata::GridReal                  inputdensity_;
        std::unique_ptr<volumedata::GridReal> outputdensity_;
        real                                  spacing_;
        bool                                  bPrint_;


};

Map::Map()
    : sigma_(0.4), n_sigma_(5), spacing_(0.2), bPrint_(false)
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

    options->addOption(FileNameOption("mi").filetype(eftVolumeData).inputFile()
                           .store(&fnmapinput_).defaultBasename("ccp4in")
                           .description("CCP4 density map input file"));
    options->addOption(FileNameOption("mo").filetype(eftVolumeData).outputFile()
                           .store(&fnmapoutput_).defaultBasename("ccp4out")
                           .description("CCP4 density map output file"));
    options->addOption(FloatOption("sigma").store(&sigma_)
                           .description("Create a simulated density by replacing the atoms by Gaussian functions of width sigma (nm)"));
    options->addOption(FloatOption("N_sigma").store(&n_sigma_)
                           .description("How many Gaussian width shall be used for spreading?"));
    options->addOption(FloatOption("spacing").store(&spacing_)
                           .description("Spacing of the density grid (nm)"));
    options->addOption(BooleanOption("print").store(&bPrint_)
                           .description("Output density information to terminal, then exit."));


}

void
Map::initAnalysis(const TrajectoryAnalysisSettings & /*settings*/, const TopologyInformation & /*top*/)
{

}

void
Map::optionsFinished(TrajectoryAnalysisSettings * /*settings*/)
{
    gt_.set_sigma(sigma_);
    gt_.set_n_sigma(5);
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
    if (!fnmapoutput_.empty())
    {
        outputdensity_ = std::unique_ptr<volumedata::GridReal>(new volumedata::GridReal());
    }
}

void
Map::set_box_from_frame(const t_trxframe &fr, matrix box, rvec translation)
{
    if (det(fr.box) > 1e-6)
    {
        copy_mat(fr.box, box);
        return;
    }

    fprintf(stderr, "Did not find suitable box for atom in structure file, guessing from structure extend.\n");
    clear_mat(box);
    std::vector<RVec> x_RVec(fr.natoms);
    x_RVec.assign(fr.x, fr.x+fr.natoms);
    for (int i = XX; i <= ZZ; i++)
    {
        box[i][i]      = 2*n_sigma_*sigma_;
        box[i][i]     += (*max_element(x_RVec.begin(), x_RVec.end(), [ = ](RVec a, RVec b){ return a[i] < b[i]; }))[i];
        box[i][i]     -= (*min_element(x_RVec.begin(), x_RVec.end(), [ = ](RVec a, RVec b){ return a[i] < b[i]; }))[i];
        translation[i] = -n_sigma_*sigma_+(*min_element(x_RVec.begin(), x_RVec.end(), [ = ](RVec a, RVec b){ return a[i] < b[i]; }))[i];
    }

}

void
Map::set_finitegrid_from_box(matrix box, rvec translation)
{
    gmx::volumedata::IVec extend({(int)ceil(box[XX][XX]/spacing_), (int)ceil(box[YY][YY]/spacing_), (int)ceil(box[ZZ][ZZ]/spacing_)});
    outputdensity_->set_extend(extend);
    outputdensity_->set_cell({extend[XX]*spacing_, extend[YY]*spacing_, extend[ZZ]*spacing_}, {90, 90, 90});
    outputdensity_->set_translation({roundf(translation[XX]/spacing_)*spacing_, roundf(translation[YY]/spacing_)*spacing_, roundf(translation[ZZ]/spacing_)*spacing_});
}

void
Map::analyzeFrame(int /*frnr*/, const t_trxframe &fr, t_pbc * /*pbc*/,
                  TrajectoryAnalysisModuleData * /*pdata*/)
{
    if (!fnmapoutput_.empty())
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
        }
        outputdensity_->zero();

        gt_.set_grid(std::move(outputdensity_));
        for (int i = 0; i < fr.natoms; ++i)
        {
            gt_.transform(fr.x[i], 1);
        }
        outputdensity_ = std::move(gt_.finish_and_return_grid());

        volumedata::MrcFile ccp4outputfile;
        ccp4outputfile.write(fnmapoutput_, *outputdensity_);
    }
}


void
Map::finishAnalysis(int /*nframes*/)
{
}


void
Map::writeOutput()
{

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
