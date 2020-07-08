/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
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
 *
 * \brief Implements routines in mapcompare.h .
 *
 * \author Christian Blau <blau@kth.se>
 */

#include "gmxpre.h"

#include "mapcompare.h"

#include <algorithm>
#include <numeric>
#include <stdexcept>

#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/fileio/mrcdensitymap.h"
#include "gromacs/math/fsc.h"
#include "gromacs/math/densityfit.h"
#include "gromacs/math/coordinatetransformation.h"

namespace gmx
{

namespace
{
class MapCompare : public ICommandLineOptionsModule
{
public:
    MapCompare() {}

    // From ICommandLineOptionsModule
    void init(CommandLineModuleSettings* /*settings*/) override {}
    void initOptions(IOptionsContainer* options, ICommandLineOptionsModuleSettings* settings) override;
    void optionsFinished() override;
    int  run() override;

private:
    std::string fnReferenceMap_  = "ref.mrc";
    std::string fnComparisonMap_ = "comp.mrc";
    std::string fnResults_       = "out.xvg";
    int         numFscShells_    = 61;
};

void MapCompare::initOptions(IOptionsContainer* options, ICommandLineOptionsModuleSettings* settings)
{
    const char* const desc[] = {
        "[THISMODULE] compares two density maps to another with different methods."
    };

    settings->setHelpText(desc);

    options->addOption(
            StringOption("refmap").required().store(&fnReferenceMap_).description("Reference to compare to."));

    options->addOption(
            StringOption("compmap").required().store(&fnComparisonMap_).description("Map to be compared."));

    options->addOption(FileNameOption("o")
                               .filetype(eftPlot)
                               .outputFile()
                               .required()
                               .store(&fnResults_)
                               .defaultBasename("similarity")
                               .description("Similarity between maps."));
}

void MapCompare::optionsFinished() {}

void normalizeMap(gmx::basic_mdspan<float, gmx::dynamicExtents3D, gmx::layout_right> data)
{
    const auto norm = std::abs(std::accumulate(begin(data), end(data), 0.));
    std::transform(begin(data), end(data), begin(data), [norm](float value) { return value / norm; });
}

int MapCompare::run()
{
    MrcDensityMapOfFloatFromFileReader readReference(fnReferenceMap_);
    auto                               refData = readReference.densityDataCopy();
    MrcDensityMapOfFloatFromFileReader readComparison(fnComparisonMap_);
    auto                               compdata = readComparison.densityDataCopy();

    normalizeMap(refData);
    normalizeMap(compdata);

    std::vector<DensitySimilarityMeasure> measure;
    for (const auto& method : EnumerationWrapper<DensitySimilarityMeasureMethod>{})
    {
        measure.emplace_back(method, refData);
    }

    RVec unitVector = { 1, 0, 0 };
    readReference.transformationToDensityLattice().scaleOperationOnly().inverseIgnoringZeroScale(
            { &unitVector, &unitVector + 1 });

    FourierShellCorrelation fsc(refData.asConstView(), unitVector[0], numFscShells_);

    const FourierShellCorrelationCurve& curve = fsc.fscCurve(compdata.asConstView());

    const size_t fscAvgIndexTwoTimesPixelsize =
            static_cast<size_t>(std::floor(1.0 / (2.0 * unitVector[0] * fsc.spacing())) - 1);
    fprintf(stderr, " pixelsize: %12.5g, spacing: %12.5g, index: %lu", unitVector[0], fsc.spacing(),
            fscAvgIndexTwoTimesPixelsize);
    const auto fscAverageCurve = fscAverage(curve);
    if (fscAvgIndexTwoTimesPixelsize > fscAverageCurve.size())
    {
        throw std::range_error(
                "FSC curve not long enough to evaluate fscavg score at 2 * pixelsize.");
    }

    auto resultsFile = fopen(fnResults_.c_str(), "w");
    fprintf(resultsFile, "inner-product\trelative-entropy\tcross-correlation\tjensen-shannon\tfscavg-two-pixelsize\tfscavg-four-pixelsize\n");
    fprintf(resultsFile, "%15.5g\t%15.5g\t%15.5g\t%15.5g", measure[0].similarity(compdata),
            measure[1].similarity(compdata), measure[2].similarity(compdata),
            measure[3].similarity(compdata));
    fprintf(resultsFile, "%15.5g\t%15.5g\n", fscAverageCurve[fscAvgIndexTwoTimesPixelsize],
            fscAverageCurve[fscAvgIndexTwoTimesPixelsize / 2]);
    fclose(resultsFile);
    return EXIT_SUCCESS;
}

} // namespace

const char                       MapCompareInfo::name[]             = "mapcompare";
const char                       MapCompareInfo::shortDescription[] = "Compare two density maps.";
ICommandLineOptionsModulePointer MapCompareInfo::create()
{
    return ICommandLineOptionsModulePointer(std::make_unique<MapCompare>());
}


} // namespace gmx