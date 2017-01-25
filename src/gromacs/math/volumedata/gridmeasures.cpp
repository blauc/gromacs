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

#include "gridmeasures.h"
#include "fouriertransform.h"
#include "gridreal.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "volumedata.h"
#include <map>
#include <iterator>
#include <set>

namespace gmx
{
namespace volumedata
{

GridMeasures::GridMeasures(const GridReal &reference)
    : reference_ {reference}
{};

real GridMeasures::correlate_(const std::vector<real> &a,
                              const std::vector<real> &b) const
{
    auto              aMean           = ScalarGridDataProperties<real>(a).mean();
    auto              aSSE            = ScalarGridDataProperties<real>(a).sumOfSquareDeviation();
    auto              bMean           = ScalarGridDataProperties<real>(b).mean();
    auto              bSSE            = ScalarGridDataProperties<real>(b).sumOfSquareDeviation();
    std::vector<real> mulArray;
    std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(mulArray), [aMean, bMean](real aValue, real bValue) { return (aValue - aMean) * (bValue - bMean); });
    return std::accumulate(mulArray.begin(), mulArray.end(), 0.) / sqrt(aSSE*bSSE);
}

real GridMeasures::correlate(const GridReal &other, real threshold) const
{
    std::vector<real> referenceAboveThreshold;
    std::vector<real> otherWhereReferenceAboveThreshold;

    auto              otherDatum = other.access().begin();
    for (auto &referenceDatum : reference_.access())
    {
        if (referenceDatum > threshold)
        {
            referenceAboveThreshold.push_back(referenceDatum);
            otherWhereReferenceAboveThreshold.push_back(*otherDatum);
        }
        ++otherDatum;
    }


    return correlate_(referenceAboveThreshold, otherWhereReferenceAboveThreshold);
};

real GridMeasures::gridSumAtCoordiantes(const std::vector<RVec> &coordinates)
{
    return std::accumulate(std::begin(coordinates), std::end(coordinates), 0.,
                           [this] (const real sum, const RVec r){
                               return sum + reference_.getLinearInterpolationAt(r);
                           });
};

real GridMeasures::getRelativeKLCrossTermSameGrid(
        const GridReal &other, const std::vector<real> &other_reference) const
{
    auto P = reference_.access().data();
    auto Q = other.access().data();
    // for numerical stability use a reference density
    if (P.size() != Q.size())
    {
        GMX_THROW(APIError(
                          "KL-Divergence calculation requires euqally sized input vectors."));
    }
    auto p           = P.begin();
    auto q           = Q.begin();
    auto q_reference = other_reference.begin();
    int  size        = Q.size();
    real sum         = 0;
#pragma omp parallel for num_threads(std::max( \
    1, gmx_omp_nthreads_get(emntDefault))) reduction(+ : sum) schedule( \
    static, size / std::max(1, gmx_omp_nthreads_get(emntDefault)) + 1)
    for (int i = 0; i < size; ++i)
    {
        if ((p[i] > 0) && (q[i] > 0) && (q_reference[i] > 0))
        {
            sum += p[i] * log(q[i] / (q_reference[i]));
        }
    }
    return sum;
}

real GridMeasures::getKLSameGrid(const GridReal &other) const
{
    auto P = reference_.access().data();
    auto Q = other.access().data();
    if (P.size() != Q.size())
    {
        GMX_THROW(APIError(
                          "KL-CrossTerm calculation requires euqally sized input vectors."));
    }
    auto p    = P.begin();
    auto q    = Q.begin();
    int  size = Q.size();
    real sum  = 0;
#pragma omp parallel for num_threads(std::max( \
    1, gmx_omp_nthreads_get(emntDefault))) reduction(+ : sum) schedule( \
    static, size / std::max(1, gmx_omp_nthreads_get(emntDefault)) + 1)
    for (int i = 0; i < size; ++i)
    {
        if ((p[i] > 0) && (q[i] > 0))
        {
            sum += p[i] * log(q[i] / p[i]);
        }
    }
    return sum;
};

real GridMeasures::getKLCrossTermSameGrid(const GridReal &other) const
{
    auto P = reference_.access().data();
    auto Q = other.access().data();
    if (P.size() != Q.size())
    {
        GMX_THROW(APIError(
                          "KL-CrossTerm calculation requires euqally sized input vectors."));
    }
    auto p    = P.begin();
    auto q    = Q.begin();
    int  size = Q.size();
    real sum  = 0;
#pragma omp parallel for num_threads(std::max( \
    1, gmx_omp_nthreads_get(emntDefault))) reduction(+ : sum) schedule( \
    static, size / std::max(1, gmx_omp_nthreads_get(emntDefault)) + 1)
    for (int i = 0; i < size; ++i)
    {
        if ((p[i] > 0) && (q[i] > 0))
        {
            sum += p[i] * log(q[i]);
        }
    }
    return sum;
};

real GridMeasures::entropy() const
{
    auto P    = reference_.access().data();
    auto p    = P.begin();
    int  size = P.size();
    real sum  = 0;
#pragma omp parallel for num_threads(std::max( \
    1, gmx_omp_nthreads_get(emntDefault))) reduction(+ : sum) schedule( \
    static, size / std::max(1, gmx_omp_nthreads_get(emntDefault)) + 1)
    for (int i = 0; i < size; ++i)
    {
        if (p[i] > 0)
        {
            sum += p[i] * log(p[i]);
        }
    }
    return sum;
};
}
}
