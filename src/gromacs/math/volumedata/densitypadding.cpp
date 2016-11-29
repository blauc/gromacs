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

#include "densitypadding.h"
#include "field.h"

namespace gmx
{
namespace volumedata
{

std::unique_ptr < Field < real>>
DensityPadding::pad(RVec paddingFactor) {
    std::unique_ptr < Field < real>> padded(new Field<real>);
    padded->copy_grid(toPad_);
    padded->multiplyGridPointNumber(paddingFactor);
    padded->scaleCell(paddingFactor);
    std::fill(std::begin(padded->access().data()),
              std::end(padded->access().data()), 0.);
    IVec extend = toPad_.extend();
    for (int iZZ = 0; iZZ < extend[ZZ]; ++iZZ)
    {
        for (int iYY = 0; iYY < extend[YY]; ++iYY)
        {
            for (int iXX = 0; iXX < extend[XX]; ++iXX)
            {
                padded->access().at({iXX, iYY, iZZ}) =
                    toPad_.access().at({iXX, iYY, iZZ});
            }
        }
    }
    return padded;
}
DensityPadding::DensityPadding(const Field<real> &toPad) : toPad_ {toPad}
{}

std::unique_ptr < Field < real>>
DensityPadding::unpad(IVec unPadExtend) {
    std::unique_ptr < Field < real>> unpadded(new Field<real>);
    unpadded->copy_grid(toPad_);
    unpadded->set_extend(unPadExtend);
    unpadded->resetCell();
    for (int iZZ = 0; iZZ < unPadExtend[ZZ]; ++iZZ)
    {
        for (int iYY = 0; iYY < unPadExtend[YY]; ++iYY)
        {
            for (int iXX = 0; iXX < unPadExtend[XX]; ++iXX)
            {
                unpadded->access().at({iXX, iYY, iZZ}) =
                    toPad_.access().at({iXX, iYY, iZZ});
            }
        }
    }

    return unpadded;
}
}
}
