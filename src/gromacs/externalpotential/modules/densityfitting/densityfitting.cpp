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
#include "densitydifferential.h"
#include "densitydifferentialprovider.h"

#include "gromacs/externalpotential/atomgroups/wholemoleculegroup.h"
#include "gromacs/externalpotential/mpi-helper.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/fileio/volumedataio.h"
#include "gromacs/fileio/json.h"
#include "gromacs/math/volumedata/gausstransform.h"
#include "gromacs/math/volumedata/gridreal.h"
#include "gromacs/math/volumedata/gridmeasures.h"
#include "gromacs/fileio/xtcio.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/mtop_lookup.h"
#include "gromacs/topology/atoms.h"
#include "emscatteringfactors.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/utility/gmxomp.h"



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

void DensityFitting::inv_mul(std::vector<real> &to_invert, const std::vector<real> &multiplier)
{
#pragma omp parallel for num_threads(number_of_threads_) schedule(static, to_invert.size()/number_of_threads_ + 1)
    for (std::size_t i = 0; i < to_invert.size(); ++i)
    {
        to_invert[i] = multiplier[i]/(to_invert[i]);
    }
}

RVec
DensityFitting::shiftedAndOriented(const RVec x)
{
    RVec x_translated;

    rvec_sub(x, centerOfMass_, x_translated);
    orientation_.rotate(x_translated);
    rvec_inc(x_translated, centerOfMass_);

    rvec_inc(x_translated, translation_);
    return x_translated;
}

real DensityFitting::getTotalScatteringSum_(WholeMoleculeGroup * atomgroup)
{
    real weightssum = atomgroup->local_weights_sum();
    if (mpi_helper() != nullptr)
    {
        mpi_helper()->to_reals_buffer(&weightssum, 1);
        mpi_helper()->sum_allReduce();
    }
    return weightssum;
}

void DensityFitting::setCenterOfMass(WholeMoleculeGroup * atomgroup)
{
    centerOfMass_ = atomgroup->weighted_local_coordinate_sum();
    real weightssum = atomgroup->local_weights_sum();
    if (mpi_helper() != nullptr)
    {
        mpi_helper()->to_reals_buffer(centerOfMass_, 3);
        mpi_helper()->to_reals_buffer(&weightssum, 1);
        mpi_helper()->sum_reduce();
        if (mpi_helper()->isMaster())
        {
            mpi_helper()->from_reals_buffer(centerOfMass_, 3);
            mpi_helper()->from_reals_buffer(&weightssum, 1);
        }
    }
    if (mpi_helper() == nullptr || mpi_helper()->isMaster())
    {
        svmul(1.0/weightssum, centerOfMass_, centerOfMass_);
    }

}

void DensityFitting::sumReduceNormalize_()
{
    if (mpi_helper() != nullptr)
    {
        sum_reduce_simulated_density_();
    }
    if (mpi_helper() == nullptr || mpi_helper()->isMaster())
    {
        simulated_density_->multiply(simulated_density_->grid_cell_volume() * totalScatteringSum_);
    }
}

real DensityFitting::KLDivergenceFromTargetOnMaster(WholeMoleculeGroup * atomgroup)
{
    spreadLocalAtoms_(atomgroup);
    sumReduceNormalize_();
    if (mpi_helper() == nullptr || mpi_helper()->isMaster())
    {
        return volumedata::GridMeasures(*target_density_).getKLSameGrid(*simulated_density_);
    }
    else
    {
        return 0;
    }
}

void
DensityFitting::alignComDensityAtoms()
{
    if (mpi_helper() == nullptr || mpi_helper()->isMaster())
    {
        rvec_sub(target_density_->center_of_mass(), centerOfMass_, translation_);
    }

    if (mpi_helper() != nullptr)
    {
        mpi_helper()->broadcast(translation_, 3);
    }
}

bool
DensityFitting::optimizeTranslation(WholeMoleculeGroup * translationgroup, real &divergenceToCompareTo)
{
    const real step_size                  = 0.02;
    bool       succededReducingDivergence = false;
    // calculate gradient vector
    // normalize to step size
    // move along gradient vector
    // accept translation if new potential better
    auto translationBeforeTrial = translation_;
    spreadLocalAtoms_(translationgroup);
    volumedata::MrcFile().write("test.mrc", *simulated_density_);
    sumReduceNormalize_();

    invertMultiplySimulatedDensity_();
    forceCalculation(translationgroup);
    auto force = translationgroup->local_force_sum();
    if (mpi_helper() != nullptr)
    {
        mpi_helper()->sum_reduce_rvec(force);
    }

    if (mpi_helper() == nullptr || mpi_helper()->isMaster())
    {
        unitv(force, force);
    }

    fprintf(input_output()->output_file(), "#\tMinimization\n"
            "#\tMinimization - New direction:\t[ %g %g %g ]\n"
            "#\tMinimization - Step size:\t%g nm .\n", force[XX], force[YY], force[ZZ], step_size);
    RVec force_scaled;
    if (mpi_helper() == nullptr || mpi_helper()->isMaster())
    {
        svmul(step_size, force, force_scaled);
        rvec_inc(translation_, force_scaled);
        fprintf(input_output()->output_file(), "#\tMinimization - Move attemted to\t[ %g %g %g ]\n", force_scaled[XX], force_scaled[YY], force_scaled[ZZ]);
    }
    if (mpi_helper() != nullptr)
    {
        mpi_helper()->broadcast(&translation_, 1);
    }
    auto new_divergence =  std::abs(KLDivergenceFromTargetOnMaster(translationgroup));
    if (mpi_helper() == nullptr || mpi_helper()->isMaster())
    {

        fprintf(input_output()->output_file(), "#\tMinimization - Divergence before move: %g , after move: %g\n", divergenceToCompareTo, new_divergence);

        if (new_divergence < divergenceToCompareTo)
        {
            succededReducingDivergence   = true;
            divergenceToCompareTo        = new_divergence;
        }
        else
        {
            translation_ = translationBeforeTrial;
        }
    }
    return succededReducingDivergence;
}

bool
DensityFitting::optimizeOrientation(WholeMoleculeGroup * atomgroup, real &divergenceToCompareTo)
{
    const real rotationAngle               = M_PI/128.;
    bool       succeededReducingDivergence = false;

    // calculate gradient vector
    // normalize to step size
    // move along gradient vector
    // accept if new potential better than old
    // move back along gradient vector

    spreadLocalAtoms_(atomgroup);
    sumReduceNormalize_();
    invertMultiplySimulatedDensity_();
    forceCalculation(atomgroup);

    auto torque_direction = atomgroup->local_torque_sum(centerOfMass_); //TODO: check for sign of rotation

    if (mpi_helper() != nullptr)
    {
        mpi_helper()->sum_reduce_rvec(torque_direction);
    }

    if (mpi_helper() == nullptr || mpi_helper()->isMaster())
    {
        unitv(torque_direction, torque_direction);
    }

    fprintf(input_output()->output_file(), "#\tMinimization in Orientation Space \n"
            "#\tMinimization - Torque Vector:\t[ %g %g %g ]\n",
            torque_direction[XX], torque_direction[YY], torque_direction[ZZ]);

    auto orientationBeforeTrial = orientation_;

    if (mpi_helper() == nullptr || mpi_helper()->isMaster())
    {
        orientation_ *= Quaternion(torque_direction, rotationAngle);
        orientation_.normalize();
    }
    if (mpi_helper() != nullptr)
    {
        mpi_helper()->broadcast(&orientation_[0], 4);
    }

    auto new_divergence =  std::abs(KLDivergenceFromTargetOnMaster(atomgroup));

    if (mpi_helper() == nullptr || mpi_helper()->isMaster())
    {

        fprintf(input_output()->output_file(), "#\tMinimization - Divergence before rotation: %g , after rotation: %g\n", divergenceToCompareTo, new_divergence);

        if (new_divergence < divergenceToCompareTo)
        {
            succeededReducingDivergence = true;
            divergenceToCompareTo       = new_divergence;
        }
        else
        {
            orientation_ = orientationBeforeTrial;
            if (mpi_helper() != nullptr)
            {
                mpi_helper()->broadcast(&orientation_[0], 4);
            }
        }
    }

    return succeededReducingDivergence;
}


void DensityFitting::translate_atoms_into_map_(WholeMoleculeGroup * translationgroup)
{

    std::vector<real> upperHalfGrid {
        1.
    };
    std::vector<real> gridPoints {
        0
    };
    Quaternion bestOrientation                  = orientation_;
    RVec       bestTranslation                  = translation_;
    real       bestDivergence                   = std::abs(KLDivergenceFromTargetOnMaster(translationgroup));
    auto       bestDivergenceCurrentOrientation = bestDivergence;
    auto       totalNumberOrientations          = upperHalfGrid.size()*gridPoints.size()*gridPoints.size()*gridPoints.size();
    decltype(totalNumberOrientations) currentOrientationNumber = 0;

    for (auto q0 : upperHalfGrid)
    {
        for (auto q1 : gridPoints)
        {
            for (auto q2 : gridPoints)
            {
                for (auto q3 : gridPoints)
                {

                    fprintf(stderr, "Optimizing Orientation: %3lu / %3lu = ", ++currentOrientationNumber, totalNumberOrientations);

                    orientation_ = Quaternion(Quaternion::QVec {{q0, q1, q2, q3}});
                    orientation_.normalize();
                    fprintf(stderr, "[ %3g %3g %3g %3g ]", orientation_[0], orientation_[1], orientation_[2], orientation_[3]);

                    bestDivergenceCurrentOrientation =  std::abs(KLDivergenceFromTargetOnMaster(translationgroup));

                    while (optimizeTranslation(translationgroup, bestDivergenceCurrentOrientation) || optimizeOrientation(translationgroup, bestDivergenceCurrentOrientation))
                    {
                        fprintf(stderr, "\r\t\t\t\t\t\t\tcurrent fit = %g ; best fit = %g \t\t", bestDivergenceCurrentOrientation, bestDivergence);
                    }

                    if (bestDivergenceCurrentOrientation < bestDivergence)
                    {
                        bestTranslation = translation_;
                        bestOrientation = orientation_;
                        bestDivergence  = bestDivergenceCurrentOrientation;
                    }
                }
            }
        }
    }
    translation_ = bestTranslation;
    orientation_ = bestOrientation;
    orientation_.normalize();
}

void
DensityFitting::initialize_target_density_()
{
    target_density_->normalize();
}

void
DensityFitting::initializeThreadLocalBuffers_()
{
#pragma omp parallel for ordered num_threads(number_of_threads_)
    for (int thread = 0; thread < number_of_threads_; ++thread)
    {
        simulated_density_buffer_.emplace_back(new volumedata::GridReal(*target_density_));
        gauss_transform_.emplace_back(new volumedata::FastGaussianGridding());
    }
}

void
DensityFitting::initialize_spreading_()
{
    simulated_density_->copy_grid(*target_density_);
    for (auto &transform : gauss_transform_)
    {
        transform->set_sigma(sigma_);
        transform->set_n_sigma(n_sigma_);
    }
}

void
DensityFitting::initialize(const matrix box, const rvec x[])
{
    if (bWriteXTC_)
    {
        out_  = open_xtc(trajectory_name_.c_str(), "w");
    }
    auto fitatoms = wholemoleculegroup(x, box, 0);
    setCenterOfMass(fitatoms);
    alignComDensityAtoms();

    initialize_target_density_();
    initializeThreadLocalBuffers_();

    totalScatteringSum_ = getTotalScatteringSum_(fitatoms);
    initialize_spreading_();

    spreadLocalAtoms_(fitatoms);
    sumReduceNormalize_();

    reference_density_ = simulated_density_->access().data();

    if (isCenterOfMassCentered_)
    {
        translate_atoms_into_map_(fitatoms);
    }

    absolute_target_divergence_ = volumedata::GridMeasures(*target_density_).getKLSameGrid(*simulated_density_);
    const real timeStepNs = 0.004; // TODO: pass this from simulation input
    optimalDeltaPotentialEnergy_ = fitatoms->num_atoms_global()*std::sqrt((float)every_nth_step_) * timeStepNs * maximumEnergyFluctuationPerAtom_;
    // fitatoms->medianSort();

    fprintf(input_output()->output_file(), "#KL-divergence(target|simulated): %g.\n", absolute_target_divergence_);
    fprintf(input_output()->output_file(), "#---Step-size estimate in DensityFitting potential space---.\n");
    fprintf(input_output()->output_file(), "#   Maximum energy fluctuation per atom: %g.\n", maximumEnergyFluctuationPerAtom_);
    fprintf(input_output()->output_file(), "#   Number of refined atoms: %d .\n", fitatoms->num_atoms_global());
    fprintf(input_output()->output_file(), "#   Fitting every %dth step.\n", every_nth_step_);
    fprintf(input_output()->output_file(), "#   With time step %g.\n", timeStepNs);
    fprintf(input_output()->output_file(), "#   Thus estimating optimal change in refinement potental energy: %g .\n", optimalDeltaPotentialEnergy_);
}


void DensityFitting::sum_reduce_simulated_density_()
{
    mpi_helper()->to_reals_buffer(simulated_density_->access().data().data(), simulated_density_->access().data().size() );
    mpi_helper()->sum_reduce();
    if (mpi_helper()->isMaster())
    {
        mpi_helper()->from_reals_buffer(simulated_density_->access().data().data(), simulated_density_->access().data().size());
    }
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
        plot.plot_force(shiftedAndOriented(*atom.xTransformed), *atom.force, 0, 1.0/max_f);
    }
    plot.stop_plot_forces();

}

void DensityFitting::forceCalculation(WholeMoleculeGroup * fitatoms)
{
    auto densityDifferential = volumedata::PotentialLibrary<volumedata::IDensityDifferentialProvider>().create(fitMethod_)();
    auto densityGradient     = densityDifferential->evaluateDensityDifferential(*simulated_density_, *target_density_);
    auto forceDensity        = ForceDensity(densityGradient, sigma_).getForce();

    std::array<volumedata::GridReal, 3> forceGrid { {
                                                        volumedata::GridReal(forceDensity[XX]), volumedata::GridReal(forceDensity[YY]), volumedata::GridReal(forceDensity[ZZ])
                                                    } };

#pragma omp parallel num_threads(number_of_threads_) shared(stderr,fitatoms,forceGrid) default(none)
    {
        int           thread             = gmx_omp_get_thread_num();
        // real          prefactor  = k_/(norm_simulated_*sigma_*sigma_);

        GroupIterator beginThreadAtoms = fitatoms->begin(thread, number_of_threads_);
        GroupIterator endThreadAtoms   = fitatoms->end(thread, number_of_threads_);

        for (auto atom = beginThreadAtoms; atom != endThreadAtoms; ++atom)
        {
            auto r = shiftedAndOriented(*(*atom).xTransformed);
            auto f = RVec {
                forceGrid[XX].getLinearInterpolationAt(r), forceGrid[YY].getLinearInterpolationAt(r), forceGrid[ZZ].getLinearInterpolationAt(r)
            };
            svmul(*((*atom).properties), f, *(*atom).force);
            orientation_.rotate_backwards(*(*atom).force);
        }
    }
}

void DensityFitting::invertMultiplySimulatedDensity_()
{
    if (mpi_helper() == nullptr || mpi_helper()->isMaster())
    {
        inv_mul(simulated_density_->access().data(), target_density_->access().data());
    }

    if (mpi_helper() != nullptr)
    {
        mpi_helper()->broadcast(simulated_density_->access().data().data(), simulated_density_->access().data().size());
        mpi_helper()->broadcast(&norm_simulated_, 1);
    }

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
            x_out[i++] = shiftedAndOriented(*atom.xTransformed);
        }
        write_xtc(out_, atoms->num_atoms_global(), step, 0, *(atoms->box()), as_vec_array(x_out.data()), 1000);
    }
}

void DensityFitting::do_potential( const matrix box, const rvec x[], const gmx_int64_t step)
{
    auto fitatoms = wholemoleculegroup(x, box, 0);
    setCenterOfMass(fitatoms);

    if ( (step%(every_nth_step_) == 0) && (bWriteXTC_))
    {
        writeTranslatedCoordinates_(fitatoms, step);
    }

    /*
     * Spread all local atoms on a grid
     */
    spreadLocalAtoms_(fitatoms);
    /*
     * Gather the local contributions to the overall spread density on the master node.
     */
    sumReduceNormalize_();
    if (mpi_helper() == nullptr || mpi_helper()->isMaster())
    {
        auto delta_divergence  = volumedata::GridMeasures(*target_density_).getRelativeKLCrossTermSameGrid(*simulated_density_, reference_density_);
        reference_divergence_ -= delta_divergence;
        set_local_potential(k_*reference_divergence_);
        reference_density_             = simulated_density_->access().data();
        exponentialDeltaEnergyAverage_ = 0.3* (k_ * delta_divergence) + 0.7 * exponentialDeltaEnergyAverage_;
        fprintf(input_output()->output_file(), "%8g\t%8g\t%8g\t", reference_divergence_, k_, exponentialDeltaEnergyAverage_);

        if (exponentialDeltaEnergyAverage_ < optimalDeltaPotentialEnergy_)
        {
            k_ *= k_factor_;
        }
        else
        {
            k_ /= k_factor_;
        }
    }
    else
    {
        set_local_potential(0);
    }

    forceCalculation(fitatoms);

    // volumedata::MrcFile simulated_output;
    if (step % (every_nth_step_) == 0)
    {
        // simulated_output.write("simulated.mrc", *simulated_density_);
        plot_forces(fitatoms);
    }
    // fitatoms->parallel_loop(std::bind( &DensityFitting::ForceKernel_KL, this, std::placeholders::_1, std::placeholders::_2));

}

DensityFitting::DensityFitting() : ExternalPotential(),
                                   target_density_(std::unique_ptr<volumedata::GridReal>(new volumedata::GridReal())),
                                   simulated_density_(std::unique_ptr<volumedata::GridReal>(new volumedata::GridReal())),
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
    fseek(inputfile, 0, SEEK_END);
    long fsize = ftell(inputfile);
    fseek(inputfile, 0, SEEK_SET);  //same as rewind(f);

    auto line = (char*) malloc(fsize + 1);;
    fsize = fread(line, fsize, 1, inputfile);
    // TODO: implement JSON scheme for checking input consistency
    json::Object parsed_json {
        std::string(line)
    };

    if (parsed_json.has("k"))
    {
        k_                  = std::stof(parsed_json["k"]);
    }
    else
    {
        GMX_THROW(InconsistentInputError("Force constant needed for densityfitting. Please provide a real number " "k" " in JSON input file."));
    };

    if (parsed_json.has("sigma"))
    {
        sigma_              = std::stof(parsed_json["sigma"]);
    }
    else
    {
        fprintf(stderr, "\n No mobility estimate given, guessing sigma = 0.2 nm . \n");
        sigma_ = 0.2;
    }

    if (parsed_json.has("n_sigma"))
    {
        n_sigma_            = std::stof(parsed_json["n_sigma"]);
    }
    else
    {
        fprintf(stderr, "\n No density spread range provided, guessing n_sigma = 5 . \n");
        n_sigma_ = 5;
    }
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
        volumedata::MrcFile  target_input_file;
        target_density_name_ = parsed_json["target_density"];
        target_input_file.read(target_density_name_, *target_density_);
        if (parsed_json.has("target_density_used"))
        {
            volumedata::MrcFile target_output_file;
            target_output_file.write(parsed_json["target_density_used"], *target_density_);
        }
        else
        {
            volumedata::MrcFile target_output_file;
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
