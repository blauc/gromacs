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
/*!  \file
 * \brief
 * Defines volume data containers.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \inpublicapi
 */

#ifndef GMX_MATH_griddata_FIELD_H
#define GMX_MATH_griddata_FIELD_H

#include "grid.h"
#include "gromacs/math/gmxcomplex.h"
#include "gromacs/utility/real.h"

#include <vector>

namespace gmx
{

template <typename T, int N>
class Field : public std::vector<T>
{
    public:
        Field(std::unique_ptr < IGrid < N>> grid) : grid_ {std::move(grid)}
        {
            this->resize(grid_->lattice().getNumLatticePoints());
        }

        Field(const Field &other) : std::vector<T>(other), grid_ {other.grid_->duplicate()}
        {
            this->resize(grid_->lattice().getNumLatticePoints());
        }

        const IGrid<N> &getGrid() const
        {
            return *grid_;
        }

        void setGrid(std::unique_ptr < IGrid < N>> &&grid)
        {
            grid_ = std::move(grid);
            this->resize(grid_->lattice().getNumLatticePoints());
        }

        /*! \brief
         * Access an element via multi index.
         * \throws std::out_of_range if element is out of array bounds
         */
        T &atMultiIndex(const typename IGrid<N>::MultiIndex &index)
        {
            return this->at(grid_->lattice().lineariseVectorIndex(index));
        };

        /*! \brief
         * Const version to access an element via multi index.
         * \throws std::out_of_range if element is out of array bounds
         */
        const T &atMultiIndex(const typename IGrid<N>::MultiIndex &index) const
        {
            return this->at(grid_->lattice().lineariseVectorIndex(index));
        };

        typename std::vector<T>::iterator iteratorAtMultiIndex(const typename IGrid<N>::MultiIndex &index)
        {
            return this->begin() + grid_->lattice().lineariseVectorIndex(index);
        }

    private:
        std::unique_ptr < IGrid < N>> grid_;
};

typedef Field<real, 3> FieldReal3D;
typedef Field<t_complex, 3> FieldComplex3D;

}      // namespace gmx

#endif /* end of include guard: GMX_MATH_griddata_FIELD_H */
