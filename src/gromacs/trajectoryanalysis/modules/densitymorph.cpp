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

#include "gromacs/externalpotential/modules/densityfitting/potentiallibrary.h"
#include "gromacs/externalpotential/modules/densityfitting/potentialprovider.h"

#include "gromacs/externalpotential/modules/densityfitting/forcedensity.h"

#include "gromacs/fileio/volumedataio.h"

#include "gromacs/math/volumedata/operations/densitypadding.h"
#include "gromacs/math/volumedata/operations/gridmeasures.h"
#include "gromacs/math/volumedata/operations/realfieldmeasure.h"
#include "gromacs/math/volumedata/operations/modifygriddata.h"
#include "gromacs/math/volumedata/gridreal.h"
#include "gromacs/utility/exceptions.h"

#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
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

        void evaluateDensityDifferential_(const GridReal &morph, GridReal &differential);
        void evaluateFlow_(const GridReal &differential, std::array<GridReal, DIM> &densityflow);
        void applyFlowOnVoxel_(RVec f_vec, IVec gridIndex, GridDataAccess<real> &d_new, const GridDataAccess<real> &d_old);
        void scaleFlow_(std::array<GridReal, DIM> &densityflow, real scale);
        std::string          fnmobile_       = "from.ccp4";
        std::string          fntarget_       = "to.ccp4";
        std::string          fnforcedensity_ = "forcedensity.ccp4";
        std::string          fnmorphtraj_    = "morph.ccp4";
        GridReal             mobile_;
        GridReal             target_;
        int                  every_          = 10;
        int                  morphSteps_     = 10;
        real                 morphstepscale_ = 1;
        real                 sigma_          = 1;
        std::string          potentialType_;
        std::unique_ptr<IStructureDensityPotentialProvider>
        potentialProvider_;
        PotentialEvaluatorHandle potentialEvaluator;
};

void DensityMorph::initOptions(IOptionsContainer          *options,
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
        "[THISMODULE] Morphs mobile density into target density."
    };

    settings->setHelpText(desc);

    options->addOption(FileNameOption("from")
                           .filetype(eftVolumeData)
                           .inputFile()
                           .store(&fnmobile_)
                           .defaultBasename("mobile")
                           .description("CCP4 density map to oldMorph."));

    options->addOption(IntegerOption("every").store(&every_).description(
                               "Write every -every step."));

    options->addOption(IntegerOption("steps")
                           .store(&morphSteps_)
                           .description("Number of morphing iterations."));

    options->addOption(FloatOption("step-scale")
                           .store(&morphstepscale_)
                           .description("Scale for density flow."));

    options->addOption(FloatOption("sigma").store(&sigma_).description(
                               "Feature size for morph."));

    options->addOption(FileNameOption("to")
                           .filetype(eftVolumeData)
                           .inputFile()
                           .defaultBasename("mobile")
                           .store(&fntarget_)
                           .description("CCP4 density map to oldMorph into."));

    options->addOption(FileNameOption("forces")
                           .filetype(eftVolumeData)
                           .outputFile()
                           .defaultBasename("forces")
                           .store(&fnforcedensity_)
                           .description("CCP4 density maps with forces."));

    options->addOption(FileNameOption("morph")
                           .filetype(eftVolumeData)
                           .outputFile()
                           .store(&fnmorphtraj_)
                           .description("CCP4 density maps during morphing."));

    options->addOption(StringOption("potential-type")
                           .store(&potentialType_)
                           .description("potential type.")
                           .enumValue(c_potentialTypes));
}

void DensityMorph::initAnalysis(const TrajectoryAnalysisSettings & /*settings*/,
                                const TopologyInformation        & /*top*/) {}
void DensityMorph::optionsFinished(TrajectoryAnalysisSettings * /*settings*/) {}
void DensityMorph::initAfterFirstFrame(
        const TrajectoryAnalysisSettings & /*settings*/,
        const t_trxframe                 & /*fr*/) {}
void DensityMorph::analyzeFrame(int /*frnr*/, const t_trxframe & /*fr*/,
                                t_pbc * /*pbc*/,
                                TrajectoryAnalysisModuleData * /*pdata*/) {}

void DensityMorph::evaluateDensityDifferential_(const GridReal &morph, GridReal &differential)
{

    // define the derivative of the density
    auto sumSimulatedDensity =
        RealFieldMeasure(morph).sum();
    //
    // auto cc = GridMeasures(morph).correlate(target_);
    // auto normSimulation = GridReal(morph).properties().norm();
    // auto densityGradientFunction = [normSimulation, cc](
    //     real densityExperiment, real densitySimulation) {
    //   return (densityExperiment - cc * densitySimulation) / (normSimulation);
    // };
    auto densityGradientFunction =
        [sumSimulatedDensity](real densityExperiment, real densitySimulation)
        {
            return ((densityExperiment) / (densitySimulation)) *
                   (1 - densitySimulation / sumSimulatedDensity);
        };

    // evaluate density derivative
    std::transform(morph.access().begin(), morph.access().end(),
                   target_.access().begin(), differential.access().begin(),
                   densityGradientFunction);

}

void DensityMorph::evaluateFlow_(const GridReal &differential, std::array<GridReal, DIM> &densityflow)
{
    // calculate force direction on grid from density derivative
    auto paddedDifferential =
        DensityPadding(differential).pad({2.0, 2.0, 2.0});
    auto paddedDensityFlow =
        ForceDensity(*paddedDifferential, sigma_).getForce();
    for (size_t dimension = XX; dimension <= ZZ; dimension++)
    {
        densityflow[dimension] =
            *DensityPadding(paddedDensityFlow[dimension])
                .unpad(differential.extend());
    }

    // find the scale for the flow
    std::vector<real> flow_rms;
    for (auto fd : densityflow)
    {
        auto measure = RealFieldMeasure(fd);
        flow_rms.push_back(std::max(fabs(measure.min()), fabs(measure.max())));
    }

    scaleFlow_(densityflow, 1 / (*std::max_element(std::begin(flow_rms),
                                                   std::end(flow_rms))));

}

void DensityMorph::scaleFlow_(std::array<GridReal, DIM> &densityflow, real scale)
{
    for (size_t dimension = XX; dimension <= ZZ; dimension++)
    {
        ModifyGridData(densityflow[dimension]).multiply(scale);
    }
}

void DensityMorph::applyFlowOnVoxel_(RVec f_vec, IVec gridIndex, GridDataAccess<real> &d_new, const GridDataAccess<real> &d_old)
{
    auto f = norm(f_vec);

// flow into voxels
    std::array<int, 3> direction;
    for (size_t dimension = 0; dimension < DIM; dimension++)
    {
        if (f_vec[dimension] > 0)
        {
            direction[dimension] = 1;
        }
        else
        {
            direction[dimension] = -1;
        }
    }
    svmul(1. / f, f_vec, f_vec);
    IVec df {
        0, 0, 0
    };
    real flow_sum = 0;
    for (df[ZZ] = 0; df[ZZ] <= 1; df[ZZ]++)
    {
        for (df[YY] = 0; df[YY] <= 1; df[YY]++)
        {
            for (df[XX] = 0; df[XX] <= 1; df[XX]++)
            {

                auto receivingVoxel {
                    gridIndex
                };
                for (size_t dimension = 0; dimension < DIM; dimension++)
                {
                    receivingVoxel[dimension] -= df[dimension] * direction[dimension];
                }

                // NOTE: checking boundaries here,
                // but not for the flow away from voxels violates flow conversation
                if (d_new.indices().inGrid(receivingVoxel))
                {
                    RVec share;
                    for (size_t dimension = 0; dimension < DIM; dimension++)
                    {
                        share[dimension] = df[dimension] == 0 ? (1-fabs(f_vec[dimension])) : fabs(f_vec[dimension]);
                    }
                    auto voxel_voxel_flow =  share[XX]*share[YY]* share[ZZ] * fabs(f) * d_old.get(gridIndex);
                    flow_sum                 += voxel_voxel_flow;
                    d_new.at(receivingVoxel) += voxel_voxel_flow;
                }
            }
        }
    }

    // flow away from voxels
    d_new.at(gridIndex) += d_old.get(gridIndex) - flow_sum;

}

void DensityMorph::finishAnalysis(int /*nframes*/)
{
    if (every_ < 1)
    {
        GMX_THROW(
                gmx::InconsistentInputError("Every must be integer bigger than one."));
    }

    MrcFile().read(fnmobile_, mobile_);
    MrcFile().read(fntarget_, target_);

    ModifyGridData(mobile_).add_offset(0.001);
    ModifyGridData(target_).add_offset(0.001);

    GridReal                              oldMorph(mobile_);

    auto                                  common_grid = (FiniteGrid)mobile_;
    GridReal                              differential(common_grid);
    GridReal                              newMorph(common_grid);
    std::array<GridReal, DIM>             densityflow;
    fprintf(stderr, "\n");
    fflush(stderr);
    auto basesum = RealFieldMeasure(mobile_).sum();
    for (int iMorphIterations = 0; iMorphIterations < morphSteps_;
         iMorphIterations++)
    {
        auto densityDistance = GridMeasures(target_).getKLSameGrid(oldMorph);
        fprintf(stderr, "\r Iteration [%7d/%7d] : d = %7g ", iMorphIterations+1, morphSteps_, densityDistance);

        evaluateDensityDifferential_(oldMorph, differential);
        evaluateFlow_(differential, densityflow);
        auto newDensityDistance = GMX_FLOAT_MAX;
        while (densityDistance < newDensityDistance)
        {
            morphstepscale_ /= 2.;
            if (morphstepscale_ < 1e2*GMX_FLOAT_EPS)
            {
                GMX_THROW(ToleranceError("Cannot move closer to target density before maximum number of steps reached."));
            }
            scaleFlow_(densityflow, morphstepscale_);
            /*
             * new, morphed state is:
             *          density(r) <- (1 - s * |flow(r)|) * density(r)
             *          density(r+|flow_direction|_max) <- s * |flow(r)| * density(r)
             */
            // go through grid and use densityflow to shift density

            auto d_new = newMorph.access();
            auto d_old = oldMorph.access();

            auto extend = common_grid.extend();

            std::array<GridDataAccess<real>, 3> flow {
                densityflow[XX].access(), densityflow[YY].access(), densityflow[ZZ].access()
            };
            // MrcFile().write("flowxx.ccp4", densityflow[XX]);
            // MrcFile().write("flowyy.ccp4", densityflow[YY]);
            // MrcFile().write("flowzz.ccp4", densityflow[ZZ]);
            //
            // GridReal flowintensity(densityflow[XX]);
            // std::transform(oldMorph.access().begin(), oldMorph.access().end(), densityflow[XX].access().begin(), flowintensity.access().begin(),[](real a,real b){return a*b;});
            // MrcFile().write("intensityxx.ccp4", flowintensity);


            // create a new, morphed density
            ModifyGridData(newMorph).zero();

            for (int iz = 0; iz < extend[ZZ]; iz++)
            {
                for (int iy = 0; iy < extend[YY]; iy++)
                {
                    for (int ix = 0; ix < extend[XX]; ix++)
                    {

                        IVec gridIndex {
                            ix, iy, iz
                        };
                        RVec f_vec {
                            flow[XX].at(gridIndex), flow[YY].at(gridIndex),
                            flow[ZZ].at(gridIndex)
                        };
                        applyFlowOnVoxel_(f_vec, gridIndex, d_new, d_old);

                    }
                }
            }
            newDensityDistance = GridMeasures(target_).getKLSameGrid(newMorph);
            fprintf(stderr, "\r\t\t\t\t\t\td = %7g , delta = %7g sum = %13g ", newDensityDistance, morphstepscale_, RealFieldMeasure(newMorph).sum()-basesum);

        }

        if (morphstepscale_ < 0.5)
        {
            morphstepscale_ *= 4.0;
        }
        newMorph.swapData(oldMorph);
        if (iMorphIterations % every_ == 0)
        {
            std::string s {
                fnmorphtraj_
            };
            auto fnmorph = s.insert(s.size() - 5, std::to_string(iMorphIterations));
            MrcFile().write(fnmorph, oldMorph);
            auto f = fopen((fnmorph+".dat").c_str(), "w+");
            fprintf(f, "%s\n", RealFieldMeasure(oldMorph).to_string().c_str());
            fclose(f);
        }
    }
    fprintf(stderr, "\n");
}
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
