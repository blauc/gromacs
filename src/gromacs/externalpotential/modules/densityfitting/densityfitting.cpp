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
#include "gmxpre.h"

#include "densityfitting.h"
#include "forcedensity.h"

#include <cstdio>


#include <algorithm>
#include <numeric>
#include <ios>

#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/externalpotential/forceplotter.h"
#include "gromacs/externalpotential/externalpotentialIO.h"
#include "gromacs/externalpotential/forceplotter.h"
#include "gromacs/math/quaternion.h"
#include "potentiallibrary.h"
#include "rigidbodyfit.h"

#include "gromacs/externalpotential/atomgroups/wholemoleculegroup.h"
#include "gromacs/externalpotential/mpi-helper.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/fileio/volumedataio.h"
#include "gromacs/fileio/json.h"
#include "gromacs/math/volumedata/operations/gausstransform.h"
#include "gromacs/fileio/xtcio.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/mtop_lookup.h"
#include "gromacs/topology/atoms.h"
#include "emscatteringfactors.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/math/volumedata/operations/modifygriddata.h"

#include "gromacs/externalpotential/modules/densityfitting/potentialprovider.h"


namespace gmx
{

namespace externalpotential
{

std::unique_ptr<ExternalPotential> DensityFitting::create()
{
    return std::unique_ptr<ExternalPotential> (new DensityFitting());
}

real DensityFitting::single_atom_properties(GroupAtom * atom, t_mdatoms * /*mdatoms*/, gmx_localtop_t * /*topology_loc*/, const gmx_mtop_t * topology_global)
{
    int           molblock        = 0;
    const t_atom &atomPropeties   = mtopGetAtomParameters(topology_global, *(atom->i_global), &molblock);

    return atomicNumber2EmScatteringFactor(atomPropeties.atomnumber);
}

std::string
DensityFitting::logInitialisationString_(int nAtoms, int timeStepNs)
{
    std::string result;
    result  = "#---Step-size estimate in DensityFitting potential space---.\n";
    result += "#   Maximum energy fluctuation per atom: " + std::to_string( maximumEnergyFluctuationPerAtom_) + "\n";
    result += "#   Number of refined atoms:  " + std::to_string(nAtoms) + "\n";
    result += "#   Fitting every " + std::to_string(every_nth_step_) + "th step\n";
    result += "#   With time step " + std::to_string(timeStepNs) + "ns \n";
    result += "#   Thus estimating optimal change in refinement potental energy: " + std::to_string(optimalDeltaPotentialEnergy_) + "\n";
    return result;
}

void
DensityFitting::initialize(const matrix box, const rvec x[])
{
    potentialProvider_ = PotentialLibrary().create(fitMethod_)();

    if (bWriteXTC_)
    {
        out_  = open_xtc(trajectory_name_.c_str(), "w");
    }
    auto fitatoms = wholemoleculegroup(x, box, 0);
    ModifyGridData(*target_density_).normalize();

    // potentialEvaluator_ = potentialProvider_->planPotential(fitatoms->xTransformed(), fitatoms->weights(), *target_density_, options_);
    // if (isCenterOfMassCentered_)
    // {
    // RigidBodyFit().fitCoordinates(*target_density_, fitatoms->xTransformed(), fitatoms->weights(), *potentialEvaluator_);
    // }

    const real timeStepNs = 0.004; // TODO: pass this from simulation input
    optimalDeltaPotentialEnergy_ = fitatoms->num_atoms_global()*std::sqrt((float)every_nth_step_) * timeStepNs * maximumEnergyFluctuationPerAtom_;

    fprintf(input_output()->output_file(), "%s", logInitialisationString_(fitatoms->num_atoms_global(), timeStepNs).c_str());
    // fitatoms->medianSort();

}

void DensityFitting::plot_forces(WholeMoleculeGroup * plotatoms)
{

    externalpotential::ForcePlotter plot;
    static int i = 0;
    plot.start_plot_forces("forces" + std::to_string(++i) + ".bild");
    real       max_f    = -1;
    for (auto &atom : *plotatoms)
    {
        if (norm(*atom.force) > max_f)
        {
            max_f = norm(*atom.force);
        }
    }
    fprintf(stderr, "max force:  %g \n", max_f);
    for (auto &atom : *plotatoms)
    {
        plot.plot_force(orientation_.shiftedAndOriented(*atom.xTransformed, centerOfRotation_, translation_), *atom.force, 0, 1.0/max_f);
    }
    plot.stop_plot_forces();

}

void
DensityFitting::writeTranslatedCoordinates_(WholeMoleculeGroup * atoms, int step)
{
    if (bWriteXTC_)
    {
        std::vector<RVec> x_out(atoms->num_atoms_global());
        int               i = 0;
        for (auto atom : *atoms)
        {
            x_out[i++] = orientation_.shiftedAndOriented(*atom.xTransformed, centerOfRotation_, translation_);
        }
        write_xtc(out_, atoms->num_atoms_global(), step, 0, *(atoms->box()), as_vec_array(x_out.data()), 1000);
    }
}

void DensityFitting::do_potential( const matrix /*box*/, const rvec /*x*/[], const gmx_int64_t /*step*/)
{
    // auto fitatoms = wholemoleculegroup(x, box, 0);
    // // setCenterOfMass(fitatoms);
    //
    // if ( (step%(every_nth_step_) == 0) && (bWriteXTC_))
    // {
    //     writeTranslatedCoordinates_(fitatoms, step);
    // }
    //
    // if (mpi_helper() == nullptr || mpi_helper()->isMaster())
    // {
    //
    //
    //     set_local_potential(k_*reference_divergence_);
    //
    //
    //     exponentialDeltaEnergyAverage_ = 0.3* (deltaPotential) + 0.7 * exponentialDeltaEnergyAverage_;
    //     fprintf(input_output()->output_file(), "%8g\t%8g\t%8g\t", reference_divergence_, k_, exponentialDeltaEnergyAverage_);
    //
    //     if (exponentialDeltaEnergyAverage_ < optimalDeltaPotentialEnergy_)
    //     {
    //         k_ *= k_factor_;
    //     }
    //     else
    //     {
    //         k_ /= k_factor_;
    //     }
    // }
    // else
    // {
    //     set_local_potential(0);
    // }
    //
    // forceCalculation(fitatoms);
    //
    // // MrcFile simulated_output;
    // if (step % (every_nth_step_) == 0)
    // {
    //     // simulated_output.write("simulated.mrc", *simulated_density_);
    //     plot_forces(fitatoms);
    // }
    // // fitatoms->parallel_loop(std::bind( &DensityFitting::ForceKernel_KL, this, std::placeholders::_1, std::placeholders::_2));

}

DensityFitting::DensityFitting() : ExternalPotential(),
                                   target_density_(std::unique_ptr < Field < real>>(new Field<real>())),
                                   simulated_density_(std::unique_ptr < Field < real>>(new Field<real>())),
                                   translation_({0, 0, 0}
                                                ), number_of_threads_ {
    std::max(1, gmx_omp_nthreads_get(emntDefault))
},
exponentialDeltaEnergyAverage_ {
    0.0
}, orientation_ {
    Quaternion::QVec({{1., 0., 0., 0.}})
}
{
}

bool DensityFitting::do_this_step(gmx_int64_t step)
{
    return (step % every_nth_step_ == 0 );
}


void DensityFitting::read_input()
{
    FILE      * inputfile       = input_output()->input_file();
    if (inputfile == nullptr)
    {
        GMX_THROW(InconsistentInputError("Please provide an external potential input file."));
    }

    // read whole inputfile as a single string
    fseek(inputfile, 0, SEEK_END);
    long fsize = ftell(inputfile);
    fseek(inputfile, 0, SEEK_SET);  //same as rewind(f);
    auto line = (char*) malloc(fsize + 1);;
    fread(line, fsize, 1, inputfile);

    // TODO: implement JSON scheme for checking input consistency
    options_ = std::string(line);
    json::Object parsed_json {
        options_
    };

    if (parsed_json.has("k"))
    {
        k_                  = std::stof(parsed_json["k"]);
    }
    else
    {
        GMX_THROW(InconsistentInputError("Force constant needed for densityfitting. Please provide a real number " "k" " in JSON input file."));
    };


    if (parsed_json.has("method"))
    {
        fitMethod_ = parsed_json["method"];
    }
    else
    {
        fitMethod_ = "kullback-leibler";
    }

    if (parsed_json.has("energy_step"))
    {
        maximumEnergyFluctuationPerAtom_ = std::stof(parsed_json.at("energy_step"));
    }
    else
    {
        maximumEnergyFluctuationPerAtom_ = 1e-4;
    }

    if (parsed_json.has("k_factor"))
    {
        k_factor_ = std::stof(parsed_json.at("k_factor"));
    }
    else
    {
        k_factor_ = 1;
    }
    if (parsed_json.has("center_to_density"))
    {
        isCenterOfMassCentered_ = parsed_json.at("center_to_density") == "true" || parsed_json.at("center_to_density") == "yes";
    }
    else
    {
        isCenterOfMassCentered_ = false;
    }

    if (parsed_json.has("write"))
    {
        bWriteXTC_ = (parsed_json.at("write") == "true") || (parsed_json.at("write") == "yes");
    }
    else
    {
        bWriteXTC_ = false;
    }

    if (parsed_json.has("every_nth_step"))
    {
        every_nth_step_ = std::stoi(parsed_json.at("every_nth_step"));
        k_             *= every_nth_step_;
    }
    else
    {
        every_nth_step_ = 1;
    }

    if (parsed_json.has("target_density"))
    {
        MrcFile  target_input_file;
        target_density_name_ = parsed_json["target_density"];
        target_input_file.read(target_density_name_, *target_density_);
        if (parsed_json.has("target_density_used"))
        {
            MrcFile target_output_file;
            target_output_file.write(parsed_json["target_density_used"], *target_density_);
        }
        else
        {
            MrcFile target_output_file;
            target_output_file.write("target.ccp4", *target_density_);
        }
    }
    else
    {
        GMX_THROW(InconsistentInputError("No target_density name given in input file. A target em-density map is required for refining atom coordinates against it."));
    }


    if (parsed_json.has("trajectory"))
    {
        trajectory_name_ = parsed_json["trajectory"];
    }
    else
    {
        trajectory_name_ = "out.xtc";
    }

    fprintf(input_output()->output_file(), "%s", dumpParsedInput().c_str());
    fprintf(input_output()->output_file(), "KL-div\tk\t%%match\tF_max\tweight\n");


}

std::string
DensityFitting::dumpParsedInput()
{
    std::string result("#input: \n");
    result += "#input:\tk\t\t = " + std::to_string(k_) + "\n";
    result += "#input:\tsigma\t\t = " + std::to_string(sigma_) + "\n";
    result += "#input:\tn_sigma\t\t = " + std::to_string(n_sigma_) + "\n";
    result += "#input:\tbackground_density\t\t = " + std::to_string(background_density_) + "\n";
    result += "#input:\tenergy_step\t\t = " + std::to_string(maximumEnergyFluctuationPerAtom_) + "\n";
    result += "#input:\tk_factor\t\t = " + std::to_string(k_factor_) + "\n";
    result += "#input:\tcenter_to_density\t\t = " + std::to_string(isCenterOfMassCentered_) + "\n";
    result += "#input:\tevery_nth_step\t\t = " + std::to_string(every_nth_step_) + "\n";
    result += "#input:\tmethod\t\t = " + fitMethod_ + "\n";
    result += "#input:\ttarget_density\t\t = " + target_density_name_ + "\n";

    return result;
}

void DensityFitting::broadcast_internal()
{
    if (mpi_helper() != nullptr)
    {
        mpi_helper()->broadcast(&k_, 1);
        mpi_helper()->broadcast(&k_factor_, 1);
        mpi_helper()->broadcast(&every_nth_step_, 1);
        mpi_helper()->broadcast(&isCenterOfMassCentered_, 1);
        mpi_helper()->broadcast(&sigma_, 1);
        mpi_helper()->broadcast(&n_sigma_, 1);
    }
}

void DensityFitting::finish()
{
    if (bWriteXTC_)
    {
        close_xtc(out_);
    }
}

std::string DensityFittingInfo::name                        = "densityfitting";
std::string DensityFittingInfo::shortDescription            = "calculate forces from difference to target density";
const int   DensityFittingInfo::numberIndexGroups           = 1;
const int   DensityFittingInfo::numberWholeMoleculeGroups   = 1;
externalpotential::ModuleCreator DensityFittingInfo::create = DensityFitting::create;

} // namespace externalpotential
} // namespace gmx
