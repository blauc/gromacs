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


#ifndef GMX_EXTERNALPOTENTIAL_DENSITYSPREADER_H
#define GMX_EXTERNALPOTENTIAL_DENSITYSPREADER_H

#include "gromacs/utility/real.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/math/quaternion.h"


#include <memory>
#include <vector>

#include "gromacs/math/volumedata/gridreal.h"

namespace gmx
{

class WholeMoleculeGroup;

template<typename real> class Field;
class GridReal;
class FiniteGrid;
class GaussTransform;

class DensitySpreader
{
    public:
        explicit DensitySpreader(const FiniteGrid &grid, int numberOfThreads, int n_sigma, real sigma);
        ~DensitySpreader();
        const GridReal &spreadLocalAtoms(const std::vector<RVec> &x, const std::vector<real> &weights, const RVec &translation = {0, 0, 0}, const Quaternion &orientation = {{1, 0, 0}, 0}, const RVec &centerOfRotation = {0, 0, 0}) const;
        void zero() const;
        const GridReal &getSpreadGrid() const;

    private:
        std::unique_ptr < std::vector < std::unique_ptr<GaussTransform>>> gauss_transform_;
        std::unique_ptr < std::vector < std::unique_ptr< GridReal>>> simulated_density_buffer_;
        std::unique_ptr< GridReal > simulated_density_;
        int number_of_threads_;
        const GridReal             &sumThreadLocalGrids_(const std::vector<IVec> &minimumUsedGridIndex, const std::vector<IVec> &maximumUsedGridIndex) const;
};

} /* gmx */


 #endif /* end of include guard: GMX_EXTERNALPOTENTIAL_DENSITYSPREADER_H */
