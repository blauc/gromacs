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
#include "vec.h"

namespace gmx
{
template <class RandomAccessIterator, class Cmp>
RandomAccessIterator
hilbert_split (RandomAccessIterator begin, RandomAccessIterator end,
               Cmp cmp = Cmp ())
{
    if (begin >= end)
    {
        return begin;
    }
    RandomAccessIterator middle = begin + (end - begin) / 2;
    std::nth_element(begin, middle, end, cmp);
    return middle;
}

template <int x> struct Hilbert_cmp_3;

template <>
struct Hilbert_cmp_3<XX>
{
    bool operator() (const rvec &p, const rvec &q) const
    {
        return (p[XX] < q[XX]);
    }
};

template <>
struct Hilbert_cmp_3<YY>
{
    bool operator() (const rvec &p, const rvec &q) const
    {
        return (p[YY] < q[YY]);
    }
};

template <>
struct Hilbert_cmp_3<ZZ>
{
    bool operator() (const rvec &p, const rvec &q) const
    {
        return (p[ZZ] < q[ZZ]);
    }
};

class Hilbert_sort_median_3
{

    private:
        std::ptrdiff_t _limit;

        template <int x> struct Cmp : public Hilbert_cmp_3<x>
        {
            Cmp () : Hilbert_cmp_3<x> (){}
        };

    public:

        template <int x, class RandomAccessIterator>
        void sort (RandomAccessIterator begin, RandomAccessIterator end) const
        {
            const int y = (x + 1) % 3, z = (x + 2) % 3;
            if (end - begin <= _limit) { return; }

            RandomAccessIterator m0 = begin, m8 = end;

            RandomAccessIterator m4 = hilbert_split (m0, m8, Cmp<x>());
            RandomAccessIterator m2 = hilbert_split (m0, m4, Cmp<y>());
            RandomAccessIterator m1 = hilbert_split (m0, m2, Cmp<z>());
            RandomAccessIterator m3 = hilbert_split (m2, m4, Cmp<z>());
            RandomAccessIterator m6 = hilbert_split (m4, m8, Cmp<y>());
            RandomAccessIterator m5 = hilbert_split (m4, m6, Cmp<z>());
            RandomAccessIterator m7 = hilbert_split (m6, m8, Cmp<z>());

            sort<z>(m0, m1);
            sort<y>(m1, m2);
            sort<y>(m2, m3);
            sort<x>(m3, m4);
            sort<x>(m4, m5);
            sort<y>(m5, m6);
            sort<y>(m6, m7);
            sort<z>(m7, m8);
        }

        template <class RandomAccessIterator>
        void operator() (RandomAccessIterator begin, RandomAccessIterator end) const
        {
            sort<0>(begin, end);
        }
};

}      // namespace gmx

#endif //CGAL_HILBERT_SORT_MEDIAN_3_H
