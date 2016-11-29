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
#include "gridreal.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"

namespace gmx
{
namespace volumedata
{


GridMeasures::GridMeasures(const GridReal &reference)
    : reference_ {reference}
{};

real GridMeasures::correlate(const GridReal &other, real threshold) const
{
    real   this_mean        = 0;
    real   this_var         = 0;
    real   other_mean       = 0;
    real   other_var        = 0;
    auto   current_value    = reference_.access().data().begin();
    auto   other_curr_value = other.access().data().begin();
    size_t count            = 0;
    real   result           = 0;
    for (size_t i = 0; i < reference_.access().data().size(); i++)
    {
        if ((*other_curr_value > threshold) || (*current_value > threshold))
        {
            this_mean  += *current_value;
            other_mean += *other_curr_value;
            ++count;
        }
        ++current_value;
        ++other_curr_value;
    }
    this_mean  /= count;
    other_mean /= count;

    current_value    = reference_.access().data().begin();
    other_curr_value = other.access().data().begin();
    for (size_t i = 0; i < reference_.access().data().size(); i++)
    {
        if ((*other_curr_value > threshold) || (*current_value > threshold))
        {
            this_var  += (*current_value) * (*current_value);
            other_var += (*other_curr_value) * (*other_curr_value);
        }
        ++current_value;
        ++other_curr_value;
    }
    this_var  = sqrt(this_var);
    other_var = sqrt(other_var);

    current_value    = reference_.access().data().begin();
    other_curr_value = other.access().data().begin();
    for (size_t i = 0; i < reference_.access().data().size(); i++)
    {
        if ((*other_curr_value > threshold) || (*current_value > threshold))
        {
            result += (*current_value) * (*other_curr_value);
        }
        ++current_value;
        ++other_curr_value;
    }
    result /= (this_var * other_var);
    return result;
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
