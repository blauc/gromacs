/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
/*! \file
 * \brief
 * Implements gmx::ColumnMajorLattice template class.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \inpublicapi
 * \ingroup module_griddata
 */
#ifndef GMX_MATH_COLUMNMAJORLATTICE
#define GMX_MATH_COLUMNMAJORLATTICE

#include "gromacs/utility/exceptions.h"

#include <array>
#include <algorithm>
#include <numeric>

namespace gmx
{

template <int N> class ColumnMajorLattice
{
    public:
        /*! \brief
         * Set the extend of the lattice.
         *
         * Indices span (0,...,0) to (extend[0]-1,...,extend[N]-1)
         */
        ColumnMajorLattice(std::array<int, N> extend)
        {
            if (!(std::all_of(std::begin(extend), std::end(extend), [](int i){return i > 0; })))
            {
                GMX_THROW(RangeError("Lattice extend must be larger zero in all dimensions."));
            }
            extend_ = extend;
        };

        /*! \brief
         * Get handle to lattice extend.
         */
        const std::array<int, N> &getExtend() const
        {
            return extend_;
        };

        int getNumLatticePoints() const
        {
            return std::accumulate(std::begin(extend_), std::end(extend_), 1, [](int lhs, int rhs){ return lhs * rhs; });
        };


        /*! \brief
         * Unique one-dimensional lattice index.
         *
         * x + extend_x * y + extend_x * extend_y * z.
         */
        int lineariseVectorIndex(const std::array<int, N> &latticeIndex) const
        {
            if (!inLattice(latticeIndex))
            {
                GMX_THROW(RangeError("Lattice index must be in lattice - greater zero and smaller than extend along all dimensions."));
            }

            int         result                 = 0;
            auto        currentDimensionExtend = extend_.begin();
            auto        pre_factor             = 1;
            for (auto ndx : latticeIndex)
            {
                result     += pre_factor * ndx;
                pre_factor *= *currentDimensionExtend;
                ++currentDimensionExtend;
            }

            return result;
        }

        /*! \brief
         * Inverse of lineariseVectorIndex
         */
        std::array<int, N> vectoriseLinearIndex(int linearIndex) const
        {
            std::array<int, N> result {{0, 0, 0}};
            auto               currentLatticeIndex = result.rbegin();
            int                stride              = getNumLatticePoints() / extend_.back();
            for (auto currentDimensionExtend = extend_.begin(); currentDimensionExtend != extend_.end(); ++currentDimensionExtend)
            {
                *currentLatticeIndex  = linearIndex / stride;
                linearIndex          -= *currentLatticeIndex * stride;
                stride               /= *currentDimensionExtend;
                ++currentLatticeIndex;
            }
            return result;
        }

        /*! \brief
         * True if latticeIndex is in Lattice
         */
        bool inLattice(const std::array<int, N> &latticeIndex) const
        {
            auto extendInCurrentDimension = extend_.begin();
            for (const auto &indexInCurrentDimension : latticeIndex)
            {
                if ((indexInCurrentDimension < 0 ) || (indexInCurrentDimension >= *extendInCurrentDimension))
                {
                    return false;
                }
                ++extendInCurrentDimension;
            }
            return true;
        }

    private:
        std::array<int, N>  extend_;
};

} //gmx

#endif
