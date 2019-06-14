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
 * Declares data structure and utilities for density fitting
 *
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "densityfitting.h"

#include "gromacs/domdec/localatomsetmanager.h"
#include "gromacs/fileio/gmxfio_xdr.h"
#include "gromacs/fileio/mrcdensitymap.h"
#include "gromacs/mdlib/broadcaststructs.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/mdtypes/imdoutputprovider.h"
#include "gromacs/mdtypes/imdpoptionprovider.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/optionsection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/selection/selectioncollection.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreetransform.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/pbcutil/pbc.h"

#include "gromacs/fileio/mrcdensitymapheader.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/trxio.h"

#include "densityfittingforceprovider.h"
#include "densityfittingparameters.h"
#include "densityfittingamplitudelookup.h"
namespace gmx
{

namespace
{

t_trxframe makeTrxFrameFromCoordinates(ArrayRef<const RVec> coordinates)
{
    t_trxframe frame;
    clear_trxframe(&frame, true);
    frame.bX     = true;
    frame.natoms = coordinates.ssize();
    frame.x      = const_cast<rvec *>(as_rvec_array(coordinates.data()));
    return frame;
}

/*! \internal
 * \brief Input data storage for density fitting
 */
class DensityFittingOptions : public IMdpOptionProvider
{
    public:
        DensityFittingOptions()
        { }

        //! From IMdpOptionProvider
        void initMdpTransform(IKeyValueTreeTransformRules * /*transform*/ ) override
        {}

        /*! \brief
         * Build mdp parameters for density fitting to be output after pre-processing.
         * \param[in, out] builder the builder for the mdp options output KV-tree.
         * \note This should be symmetrical to option initialization without
         *       employing manual prefixing with the section name string once
         *       the legacy code blocking this design is removed.
         */
        void buildMdpOutput(KeyValueTreeObjectBuilder *builder) const override
        {
            builder->addValue<bool>(inputSectionName_ + "-" +
                                    c_activeTag_, active_);
            builder->addValue<std::string>(inputSectionName_ + "-" +
                                           c_similarityMeasureTag_, similarityMeasure_);
            builder->addValue<std::string>(inputSectionName_ + "-" +
                                           c_amplitudeMethodTag_, c_densityFittingAmplitudeMethodNames[static_cast<int>(amplitudeMethod_)]);
            builder->addValue<std::string>(inputSectionName_ + "-" +
                                           c_referenceDensityFileNameTag_, referenceDensityFileName);
            builder->addValue<std::string>(inputSectionName_ + "-" +
                                           c_fittingGroupTag_, fitGroupString_);
            builder->addValue<float>(inputSectionName_ + "-" +
                                     c_forceConstantTag_, forceConstant_);
            builder->addValue<float>(inputSectionName_ + "-" +
                                     c_sigmaTag_, sigma_);
        }

        /*! \brief
         * Connect option name and data.
         */
        void initMdpOptions(IOptionsContainerWithSections *options) override
        {
            auto section = options->addSection(OptionSection(inputSectionName_.c_str()));
            section.addOption(BooleanOption(c_activeTag_.c_str()).store(&active_));
            section.addOption(StringOption(c_similarityMeasureTag_.c_str()).enumValue({"inner-product"}).store(&similarityMeasure_));
            section.addOption(EnumOption<DensityFittingAmplitudeMethod>(c_amplitudeMethodTag_.c_str()).enumValue(c_densityFittingAmplitudeMethodNames).store(&amplitudeMethod_));
            section.addOption(StringOption(c_referenceDensityFileNameTag_.c_str()).store(&referenceDensityFileName));
            section.addOption(StringOption(c_fittingGroupTag_.c_str()).store(&fitGroupString_));
            section.addOption(FloatOption(c_forceConstantTag_.c_str()).store(&forceConstant_));
            section.addOption(FloatOption(c_sigmaTag_.c_str()).store(&sigma_));
        }

        //! Report if this set of options is active
        bool active() const
        {
            return active_;
        }

        //! Process input options to parameters, including input file reading.
        DensityFittingParameters buildParameters()
        {
            GMX_ASSERT(atomSet != nullptr, "Atom set needs to be set before initializing force provider");
            t_fileio                                   *mrcFile = gmx_fio_open(referenceDensityFileName.c_str(), "r");
            FileIOXdrSerializer                         serializer(mrcFile);
            mapReader_ = std::make_unique<MrcDensityMapOfFloatReader>(&serializer);
            TranslateAndScale                           transformationToDensityLattice = getCoordinateTransformationToLattice(mapReader_->header());
            dynamicExtents3D                            ext = getDynamicExtents3D(mapReader_->header());
            basic_mdspan<const float, dynamicExtents3D> referenceDensity(mapReader_->data().data(), ext);
            return {*atomSet, transformationToDensityLattice, forceConstant_, sigma_, nSigma_, referenceDensity, amplitudeMethod_, similarityMeasure_};
        }

        void buildAtomSets(LocalAtomSetManager *atomSets)
        {
            GMX_ASSERT(fitGroup_ != nullptr, "Fit group selection needs to be built before atom sets.");
            GMX_ASSERT(fitGroup_->size() == 1, "Density guided simulations may only have one selection group. Use the merge functionality in the selection syntax.");
            *atomSet = atomSets->add(fitGroup_->front().atomIndices());
        }

        void buildSelection(SelectionCollection * selectionCollection)
        {
            GMX_ASSERT(pbc_ != nullptr, "Periodic boundaries need to be set before selection.");
            GMX_ASSERT(box_ != nullptr, "Box needs to be set before selection.");
            GMX_ASSERT(inputCoordinates_ != nullptr, "Needs to receive input coordinates before building selection.");

            *fitGroup_ = selectionCollection->parseFromString(fitGroupString_);
            selectionCollection->compile();

            t_trxframe frame = makeTrxFrameFromCoordinates(*inputCoordinates_);
            t_pbc      pbcOptions;
            set_pbc(&pbcOptions, *pbc_, *box_);
            selectionCollection->evaluate(&frame, &pbcOptions);
        }

        void setCoordinates(ArrayRef<const RVec> x)
        {
            GMX_ASSERT(commrec_ != nullptr, "Communication record needs to be set before setting coordinates.");

            int nAtoms = ssize(x);
            inputCoordinates_->resize(x.size());
            std::copy(x.begin(), x.end(), std::begin(*inputCoordinates_));
            gmx_bcast(sizeof(nAtoms), &nAtoms, commrec_.get());
            nblock_abc(commrec_.get(), nAtoms, inputCoordinates_.get());

        }

        void setPbc(int pbc)
        {
            *pbc_ = pbc;
        }

        void setBox(const matrix &box)
        {
            copy_mat(box, *box_);
        }

        void setCommrec(const t_commrec &commrec)
        {
            *commrec_ = commrec;
        }

    private:

        std::unique_ptr<t_commrec>                  commrec_;
        std::unique_ptr < std::vector < RVec>> inputCoordinates_;
        std::unique_ptr<int>                        pbc_;
        std::unique_ptr<matrix>                     box_;

        std::unique_ptr<LocalAtomSet>               atomSet;
        std::unique_ptr<MrcDensityMapOfFloatReader> mapReader_;
        //! The name of the density-fitting module
        const std::string inputSectionName_ = "density-guided-simulation";

        //! Forces are only calculated if density guided simulation is active
        const std::string c_activeTag_      = "active";
        bool              active_           = false;

        //! The type of the fitting potential
        const std::string              c_similarityMeasureTag_ = "similarity-measure";
        std::string                    similarityMeasure_;

        const std::string              c_amplitudeMethodTag_ = "amplitude-weight";
        DensityFittingAmplitudeMethod  amplitudeMethod_;

        std::unique_ptr<SelectionList> fitGroup_;
        const std::string              c_fittingGroupTag_ = "group";
        std::string                    fitGroupString_;

        const std::string              c_referenceDensityFileNameTag_ = "reference";
        std::string                    referenceDensityFileName       = "reference.mrc";

        const std::string              c_forceConstantTag_ = "force-constant";
        float                          forceConstant_;

        const std::string              c_sigmaTag_ = "sigma";
        float sigma_;

        const real                    nSigma_ = 5.;

};

/*! \brief Class that handles the output to files for density guided simulations.
 */
class DensityFittingOutputProvider final : public IMDOutputProvider
{
    public:
        //! Initialize output
        void initOutput(FILE * /*fplog*/, int /*nfile*/, const t_filenm /*fnm*/[],
                        bool /*bAppendFiles*/, const gmx_output_env_t * /*oenv*/) override
        {}
        //! Finalizes output from a simulation run.
        void finishOutput() override {}
};

/*! \internal
 * \brief Density fitting
 *
 * Class that implements the density fitting forces and potential
 * \note the virial calculation is not yet implemented
 */
class DensityFitting final : public IMDModule
{
    public:
        DensityFitting() = default;

        //! return the class that provides the options for this modules
        IMdpOptionProvider *mdpOptionProvider() override { return &densityFittingOptions_; }

        //! Add this module to the force providers if active
        void initForceProviders(ForceProviders *forceProviders) override
        {
            if (densityFittingOptions_.active())
            {
                const auto &parameters = densityFittingOptions_.buildParameters();
                forceProvider_ = std::make_unique<DensityFittingForceProvider>(parameters);
                forceProviders->addForceProvider(forceProvider_.get());
            }
        }

        void callback(LocalAtomSetManager * atomSets)
        {
            densityFittingOptions_.buildAtomSets(atomSets);
        }

        void callback(GlobalCoordinatesProvidedOnMaster coordinatesOnMaster)
        {
            densityFittingOptions_.setCoordinates(coordinatesOnMaster.coordinates);
        }

        void callback(SelectionCollection * selectionCollection)
        {
            densityFittingOptions_.buildSelection(selectionCollection);
        }

        void callback(PeriodicBoundaryConditionOptionIsSetup pbc)
        {
            densityFittingOptions_.setPbc(pbc.ePBC_);
        }

        void callback(CommunicationIsSetup communicationSetup)
        {
            densityFittingOptions_.setCommrec(communicationSetup.communicationRecord_);
        }

        void callback(BoxIsSetup boxSetup)
        {
            densityFittingOptions_.setBox(boxSetup.box_);
        }
        //! This MDModule provides its own output
        IMDOutputProvider *outputProvider() override { return &densityFittingOutputProvider_; }

    private:
        //! The output provider
        DensityFittingOutputProvider                 densityFittingOutputProvider_;
        //! The options provided for density fitting
        DensityFittingOptions                        densityFittingOptions_;
        //! Object that evaluates the forces
        std::unique_ptr<DensityFittingForceProvider> forceProvider_;
};

}   // namespace

std::unique_ptr<IMDModule> createDensityFittingModule(MdModuleCallBack * mdModuleMessageTriggers)
{
    auto densityFittingModule = std::make_unique<DensityFitting>();
    mdModuleMessageTriggers->subscribe<LocalAtomSetManager *>(densityFittingModule.get());
    mdModuleMessageTriggers->subscribe<GlobalCoordinatesProvidedOnMaster>(densityFittingModule.get());
    mdModuleMessageTriggers->subscribe<SelectionCollection *>(densityFittingModule.get());
    mdModuleMessageTriggers->subscribe<PeriodicBoundaryConditionOptionIsSetup>(densityFittingModule.get());
    mdModuleMessageTriggers->subscribe<CommunicationIsSetup>(densityFittingModule.get());
    mdModuleMessageTriggers->subscribe<BoxIsSetup>(densityFittingModule.get());
    return std::move(densityFittingModule);
}

} // namespace gmx
