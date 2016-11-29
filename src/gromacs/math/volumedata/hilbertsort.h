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
/* After CGAL */
#ifndef HILBERT_SORT_H
#define HILBERT_SORT_H

#include <functional>
#include <algorithm>
#include "gromacs/math/vec.h"
#include <vector>
#include <numeric>

namespace gmx
{

std::vector<int>::iterator
median_split (std::vector<int>::iterator begin, std::vector<int>::iterator end, const std::vector<RVec> &coordinates,
              int dimension)
{
    if (begin >= end)
    {
        return begin;
    }
    auto middle = begin + (end - begin) / 2;
    std::nth_element(begin, middle, end, [dimension, coordinates](const int i, const int j){return (coordinates[i][dimension] < coordinates[j][dimension]); });
    return middle;
}

void hilbertMedianSortRecursion3d(std::vector<int>::iterator begin, std::vector<int>::iterator end, const std::vector<RVec> &coordinates, int x = 0)
{
    const int y = (x + 1) % 3;
    const int z = (x + 2) % 3;
    if (end - begin <= 1)
    {
        return;
    }

    auto m0 = begin;
    auto m8 = end;

    auto m4 = median_split(m0, m8, coordinates, x);
    auto m2 = median_split(m0, m4, coordinates, y);
    auto m1 = median_split(m0, m2, coordinates, z);
    auto m3 = median_split(m2, m4, coordinates, z);
    auto m6 = median_split(m4, m8, coordinates, y);
    auto m5 = median_split(m4, m6, coordinates, z);
    auto m7 = median_split(m6, m8, coordinates, z);

    hilbertMedianSortRecursion3d(m0, m1, coordinates, z);
    hilbertMedianSortRecursion3d(m1, m2, coordinates, y);
    hilbertMedianSortRecursion3d(m2, m3, coordinates, y);
    hilbertMedianSortRecursion3d(m3, m4, coordinates, x);
    hilbertMedianSortRecursion3d(m4, m5, coordinates, x);
    hilbertMedianSortRecursion3d(m5, m6, coordinates, y);
    hilbertMedianSortRecursion3d(m6, m7, coordinates, y);
    hilbertMedianSortRecursion3d(m7, m8, coordinates, z);
}

void hilbertMedianSortRecursion2d(std::vector<int>::iterator begin, std::vector<int>::iterator end, const std::vector<RVec> &coordinates, int x = 0, int y = 1)
{
    std::swap(x, y);
    if (end - begin <= 1)
    {
        return;
    }

    auto m0 = begin;
    auto m4 = end;

    auto m2 = median_split(m0, m4, coordinates, x);
    auto m1 = median_split(m0, m2, coordinates, y);
    auto m3 = median_split(m2, m4, coordinates, y);

    hilbertMedianSortRecursion2d(m0, m1, coordinates, y);
    hilbertMedianSortRecursion2d(m1, m2, coordinates, x);
    hilbertMedianSortRecursion2d(m2, m3, coordinates, x);
    hilbertMedianSortRecursion2d(m3, m4, coordinates, y);
}

void
hilbertMedianSort(const std::vector<RVec> &coordinates,  std::vector<int> &result)
{
    result.resize(coordinates.size());
    std::iota(std::begin(result), std::end(result), 0);
    hilbertMedianSortRecursion2d(std::begin(result), std::end(result), coordinates, YY, ZZ);
}

}      // namespace gmx

#endif
