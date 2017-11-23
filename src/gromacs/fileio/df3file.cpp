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
/*! \internal \file
 * \brief
 * Implements methods from voxels.h
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_fileio
 * \ingroup module_griddata
 */
#include "gmxpre.h"

#include "griddataio.h"

#include "config.h"

#include <cstdio>

#include <algorithm>
#include <string>
#include <vector>
#include <set>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/futil.h"
#include "gromacs/math/units.h"
#include "gromacs/math/griddata/operations/realfieldmeasure.h"

namespace gmx
{

Df3File::SuccessfulDf3Write
Df3File::write(std::string filename, const gmx::GridDataReal3D &grid_data)
{
    const auto &grid  = grid_data.getGrid();
    auto        file_ = gmx_fio_fopen(filename.c_str(), "w");
    int16_t     xExtendShort {
        int16_t(grid.lattice().extend()[XX])
    };
    int16_t yExtendShort {
        int16_t(grid.lattice().extend()[YY])
    };
    int16_t zExtendShort {
        int16_t(grid.lattice().extend()[ZZ])
    };
    fputc(xExtendShort >> 8, file_);
    fputc(xExtendShort & 0xff, file_);
    fputc(yExtendShort >> 8, file_);
    fputc(yExtendShort & 0xff, file_);
    fputc(zExtendShort >> 8, file_);
    fputc(zExtendShort & 0xff, file_);
    auto measure = DataVectorMeasure(grid_data);
    auto gridMax = measure.max();
    auto gridMin = measure.min();
    for (auto voxel : grid_data)
    {
        auto scaledValue = float(voxel - gridMin)/float(gridMax-gridMin);
        char datum       = static_cast<char>(255*scaledValue);
        fputc(datum, file_);
    }
    gmx_fio_fclose(file_);
    return Df3File::SuccessfulDf3Write(filename, grid_data);
}

Df3File::SuccessfulDf3Write::SuccessfulDf3Write(std::string filename, const gmx::GridDataReal3D &grid_data) : filename_ {filename}, gridData_ {
    grid_data
} {};

void Df3File::SuccessfulDf3Write::writePovray()
{
    const auto &cell         = gridData_.getGrid().cell();
    auto        file         = gmx_fio_fopen((filename_+".pov").c_str(), "w");
    auto        translation  = cell.transformFromBasis({{0., 0., 0.}});
    std::string povRayString = "#declare DD = <" + std::to_string(NM2A * cell.basisVectorLength(XX)) + "," + std::to_string(NM2A *cell.basisVectorLength(YY)) + "," + std::to_string(NM2A *cell.basisVectorLength(ZZ)) + ">;\n";
    povRayString += std::string("#declare theinterior = interior {\n")
        +"\tmedia {\n"
        +"\t\temission <1,1,1> / 10\n"
        +"\t\tabsorption <1,1,1> / 3\n"
        +"\t\tscattering {1 0.5}\n"
        +"\t\tdensity {\n"
        +"\t\t\tdensity_file df3 \"" + filename_+"\"\n"
        +"\t\t\tinterpolate 1\n"
        +"\t\t\tcolor_map {\n"
        +"\t\t\t\t[0 rgb <0,0,0>]\n"
        +"\t\t\t\t[1 rgb <1,1,1>]\n"
        +"\t\t\t}\n"
        +"\t\t}\n"
        +"\t}\n"
        +"}\n"
        +"\n"
        + "box {\n"
        +"\t\t<0,0,0>, <1,1,1>\n"
        +"\t\tpigment { rgbf 1 }\n"
        +"\t\tinterior { theinterior }\n"
        +"\t\thollow\n"
        +"\t\tscale DD\n"
        +"\t\ttranslate <" + std::to_string(NM2A *translation[XX]) + ","+std::to_string(NM2A *translation[YY])+ "," + std::to_string(NM2A *translation[ZZ]) + ">\n"
        +"\t}\n";
    fprintf(file, "%s", povRayString.c_str());
    gmx_fio_fclose(file);

}

} // namespace gmx
