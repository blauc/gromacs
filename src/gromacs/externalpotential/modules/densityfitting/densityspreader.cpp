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
#include "densityspreader.h"
#include "gromacs/externalpotential/atomgroups/wholemoleculegroup.h"
#include "gromacs/math/volumedata/field.h"
#include "gromacs/math/volumedata/volumedata.h"
#include "gromacs/math/volumedata/gridreal.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/math/volumedata/gausstransform.h"
#include "gromacs/math/quaternion.h"


namespace gmx
{
namespace volumedata
{

DensitySpreader::DensitySpreader(const FiniteGrid &grid, int numberOfThreads, int n_sigma, int sigma) :
    simulated_density_ {new GridReal(grid)},
number_of_threads_ {
    numberOfThreads
}
{
#pragma omp parallel for ordered num_threads(number_of_threads_)
    for (int thread = 0; thread < number_of_threads_; ++thread)
    {
        simulated_density_buffer_.emplace_back(new volumedata::GridReal(grid));
        gauss_transform_.emplace_back(new volumedata::FastGaussianGridding());
    }
    for (auto &transform : gauss_transform_)
    {
        transform->set_sigma(sigma);
        transform->set_n_sigma(n_sigma);
    }
}

//TODO use wholemoleculegroup concept throughout
const Field<real> &
DensitySpreader::spreadLocalAtoms(const rvec * x, const std::vector<real> &weights, const int nAtoms,  const RVec &translation, const Quaternion &orientation)
{
    std::vector<volumedata::IVec>            minimumUsedGridIndex(number_of_threads_);
    std::vector<volumedata::IVec>            maximumUsedGridIndex(number_of_threads_);
    RVec centerOfMass {
        0, 0, 0
    };
    real weight_sum = 0;
    for (int atomIndex = 0; atomIndex != nAtoms; ++atomIndex)
    {
        RVec weighted_coordinate;
        svmul(weights[atomIndex], x[atomIndex], weighted_coordinate);
        weight_sum += weights[atomIndex];
        rvec_inc(centerOfMass, weighted_coordinate);
    }
    svmul(1./weight_sum, centerOfMass, centerOfMass);
#pragma omp parallel num_threads(number_of_threads_)
    {
        int           thread     = gmx_omp_get_thread_num();
        simulated_density_buffer_[thread]->zero();
        gauss_transform_[thread]->set_grid(std::move(simulated_density_buffer_[thread]));
        auto beginThreadAtoms = thread * (nAtoms / number_of_threads_ );
        auto endThreadAtoms   = (thread == number_of_threads_-1) ?  nAtoms : (thread+1) * ( nAtoms / number_of_threads_);

        for (auto atomIndex = beginThreadAtoms; atomIndex != endThreadAtoms; ++atomIndex)
        {
            gauss_transform_[thread]->transform(orientation.shiftedAndOriented(x[atomIndex], centerOfMass, translation), weights[atomIndex]);
        }
        simulated_density_buffer_[thread] = gauss_transform_[thread]->finish_and_return_grid(); //TODO:std:move ?
        minimumUsedGridIndex[thread]      = gauss_transform_[thread]->getMinimumUsedGridIndex();
        maximumUsedGridIndex[thread]      = gauss_transform_[thread]->getMaximumUsedGridIndex();
    }
    return sumThreadLocalGrids_(minimumUsedGridIndex, maximumUsedGridIndex);
}

const Field<real> &
DensitySpreader::sumThreadLocalGrids_(const std::vector<IVec> &minimumUsedGridIndex, const std::vector<IVec> &maximumUsedGridIndex)
{
    std::vector<std::vector<real>::iterator> contributingThreadLocalVoxelIterators(number_of_threads_);
    std::vector<real>::iterator              simulatedDensityVoxelIterator;

    volumedata::IVec gridStart;
    volumedata::IVec gridEnd;
    for (int i = XX; i <= ZZ; ++i)
    {
        auto ithElementLarger = [i](volumedata::IVec a, volumedata::IVec b){
                return a[i] < b[i];
            };
        gridStart[i] = (*std::min_element(std::begin(minimumUsedGridIndex), std::end(minimumUsedGridIndex), ithElementLarger))[i];
        gridEnd[i]   = (*std::max_element(std::begin(maximumUsedGridIndex), std::end(maximumUsedGridIndex), ithElementLarger))[i];
    }
    // add together all thread-local grids, using only density values, where there is
    // actual spread density
    int  nGridPointsXX        = gridEnd[XX]-gridStart[XX];
    auto simulatedDensityData = simulated_density_->access();

    // store the data accessors for threadlocal grid data in a vector
    std::vector<volumedata::GridDataAccess<real> > simulatedDensityThreadGridData;
    for (int thread = 0; thread < number_of_threads_; ++thread)
    {
        simulatedDensityThreadGridData.emplace_back(simulated_density_buffer_[thread]->extend(), simulated_density_buffer_[thread]->access().data());
    }

    //
    for (int gridIndexZZ = gridStart[ZZ]; gridIndexZZ < gridEnd[ZZ]; ++gridIndexZZ)
    {
        for (int gridIndexYY = gridStart[YY]; gridIndexYY < gridEnd[YY]; ++gridIndexYY)
        {
            contributingThreadLocalVoxelIterators.resize(0);
            // access rows in thread local grids if they contribute to the density
            for (int thread = 0; thread < number_of_threads_; ++thread)
            {
                if ((minimumUsedGridIndex[thread][ZZ] <= gridIndexZZ) && ( gridIndexZZ <= maximumUsedGridIndex[thread][ZZ] ) && (minimumUsedGridIndex[thread][YY] <= gridIndexYY) && (gridIndexYY <= maximumUsedGridIndex[thread][YY]))
                {
                    contributingThreadLocalVoxelIterators.push_back(simulatedDensityThreadGridData[thread].zy_column_begin(gridIndexZZ, gridIndexYY)+gridStart[XX]);
                }
            }
            simulatedDensityVoxelIterator = simulatedDensityData.zy_column_begin(gridIndexZZ, gridIndexYY) + gridStart[XX];
            // step though grid row by row
            for (int gridIndexXX = 0; gridIndexXX < nGridPointsXX; ++gridIndexXX)
            {
                // loop over the local threads, collect voxel contributions and advance all threadlocal iterators to next voxel in row
                for (auto &threadLocalVoxel : contributingThreadLocalVoxelIterators)
                {
                    *simulatedDensityVoxelIterator += *threadLocalVoxel;
                    ++threadLocalVoxel;
                }
                ++simulatedDensityVoxelIterator;
            }
        }
    }
    return *simulated_density_;
};

const Field<real> &
DensitySpreader::spreadLocalAtoms(const WholeMoleculeGroup * spreadgroup, const RVec &translation, const Quaternion &orientation)
{
    simulated_density_->zero();
    std::vector<volumedata::IVec>            minimumUsedGridIndex(number_of_threads_);
    std::vector<volumedata::IVec>            maximumUsedGridIndex(number_of_threads_);
    RVec centerOfMass = spreadgroup->weightedCenterOfMass();

#pragma omp parallel num_threads(number_of_threads_)
    {
        int           thread     = gmx_omp_get_thread_num();
        simulated_density_buffer_[thread]->zero();
        gauss_transform_[thread]->set_grid(std::move(simulated_density_buffer_[thread]));
        auto beginThreadAtoms = spreadgroup->begin(thread, number_of_threads_);
        auto endThreadAtoms   = spreadgroup->end(thread, number_of_threads_);
        for (auto atom = beginThreadAtoms; atom != endThreadAtoms; ++atom)
        {
            gauss_transform_[thread]->transform(orientation.shiftedAndOriented(*(*atom).xTransformed, centerOfMass, translation), *(*atom).properties);
        }
        simulated_density_buffer_[thread] = gauss_transform_[thread]->finish_and_return_grid(); //TODO:std:move ?
        minimumUsedGridIndex[thread]      = gauss_transform_[thread]->getMinimumUsedGridIndex();
        maximumUsedGridIndex[thread]      = gauss_transform_[thread]->getMaximumUsedGridIndex();
    }

    return sumThreadLocalGrids_(minimumUsedGridIndex, maximumUsedGridIndex);

}

} /* volumedata */
} /* gmx */
