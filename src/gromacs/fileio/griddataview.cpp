/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
 * Implement file data to GridDataFloat3D and back conversions.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_fileio
 */
#include "griddataview.h"

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/griddataio.h"
#include "gromacs/fileio/mrcmetadata.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/fileio/xplor.h"
#include "gromacs/math/container/containermeasure.h"
#include "gromacs/math/griddata/griddata.h"

namespace gmx {

GridDataFloat3D gridDataFromXplor(const XplorData &xplor) {
  RVec resolution;
  Grid<DIM>::NdVector translation;
  MdFloatVector<DIM> cellLength;
  bounds<DIM> extend;
  for (size_t i = 0; i < DIM; i++) {
    extend[i] = xplor.crs_end[i] - xplor.crs_start[i] + 1;
    resolution[i] = xplor.cell_length[i] / xplor.extend[i];
    cellLength[i] = (extend[i]) * resolution[i];
    // origin of the lattice
    translation[i] = xplor.crs_start[i] * resolution[i];
  }

  CanonicalVectorBasis<DIM> cell(cellLength);

  GridDataFloat3D result({cell, extend, translation});
  std::copy(std::begin(xplor.data), std::end(xplor.data), std::begin(result));

  return result;
}

Grid<DIM> cellFromXplor(const XplorData &xplor) {
  CanonicalVectorBasis<DIM> cell(MdFloatVector<DIM>(
      {xplor.cell_length[XX], xplor.cell_length[YY], xplor.cell_length[ZZ]}));
  return Grid<DIM>(cell,
                   {xplor.extend[XX], xplor.extend[YY], xplor.extend[ZZ]});
}

XplorData xplorFromGridData(const GridDataFloat3D &map,
                            std::array<int, DIM> cellExtend) {
  XplorData xplor;
  bool cellExtendIsValid =
      std::all_of(std::begin(cellExtend), std::end(cellExtend),
                  [](int i) { return i != 0; });
  xplor.cell_angles = {90., 90., 90};

  /* With valid cell extend, set up the unit cell from grid resolution and unit
   * cell lengths
   * The cell and the griddata from gromacs may have different extends
   */
  if (cellExtendIsValid) {
    for (size_t dim = XX; dim < DIM; ++dim) {
      xplor.cell_length[dim] = map.getGrid().unitCell()[dim] * cellExtend[dim];
    }
    xplor.extend = cellExtend;
  } else {
    for (size_t dim = XX; dim < DIM; ++dim) {
      xplor.cell_length[dim] = map.getGrid().cell()[dim];
      xplor.extend[dim] = static_cast<int>(map.getGrid().lattice()[XX]);
    }
  }

  auto originCoordinate = map.getGrid().multiIndexToCoordinate({{0, 0, 0}});
  for (size_t dim = XX; dim < DIM; ++dim) {
    xplor.crs_start[dim] =
        originCoordinate[dim] / map.getGrid().unitCell()[dim];
  }
  for (int dimension = XX; dimension <= ZZ; ++dimension) {
    xplor.crs_end[dimension] =
        xplor.crs_start[dimension] + map.getGrid().lattice()[dimension] - 1;
  }

  xplor.mean_value = containermeasure::mean(map);
  xplor.rms_value = containermeasure::rms(map);

  std::copy(std::begin(map), std::end(map), std::back_inserter(xplor.data));
  return xplor;
}

namespace {

void pdbFromGridData(const GridDataFloat3D &grid, FILE *outputfile) {
  const char *space = " ";
  for (int x = 0; x != grid.getGrid().lattice()[XX]; ++x) {
    for (int y = 0; y != grid.getGrid().lattice()[YY]; ++y) {
      for (int z = 0; z != grid.getGrid().lattice()[ZZ]; ++z) {
        const auto coordinate =
            grid.getGrid().multiIndexToCoordinate({x, y, z});
        gmx_fprintf_pdb_atomline(outputfile, epdbATOM, 0, "V", ' ', "VOX", ' ',
                                 0, ' ', coordinate[XX], coordinate[YY],
                                 coordinate[ZZ], 1, grid[{x, y, z}], space);
      }
    }
  }
}

void dumpFromGridData(const GridDataFloat3D &grid, FILE *outputfile) {
  for (int x = 0; x != grid.getGrid().lattice()[XX]; ++x) {
    for (int y = 0; y != grid.getGrid().lattice()[YY]; ++y) {
      for (int z = 0; z != grid.getGrid().lattice()[ZZ]; ++z) {
        const auto coordinate =
            grid.getGrid().multiIndexToCoordinate({x, y, z});
        fprintf(outputfile, "%14.7g %14.7g %14.7g %14.7g\n", coordinate[XX],
                coordinate[YY], coordinate[ZZ], grid[{x, y, z}]);
      }
    }
  }
}

} // namespace

MapConverter::MapConverter(std::string filename) {
  const auto &extension = filename.substr(filename.find_last_of(".") + 1);
  if (extension == "ccp4" || extension == "map" || extension == "mrc") {
    grid_ = MrcFile().read(filename);
  }
  if (extension == "xplor") {
    auto outputfile = gmx_fio_fopen(filename.c_str(), "r");
    grid_ = gridDataFromXplor(XplorData(outputfile));
    gmx_fio_fclose(outputfile);
  }
}

void MapConverter::MapConverter::to(std::string filename) {
  const auto &extension = filename.substr(filename.find_last_of(".") + 1);
  if (extension == "ccp4" || extension == "map" || extension == "mrc") {
    MrcFile().write(filename, grid_);
    return;
  }

  auto outputfile = gmx_fio_fopen(filename.c_str(), "w");

  if (extension == "xplor") {
    xplorFromGridData(grid_).write(outputfile);
  }
  if (extension == "dat") {
    const auto textDump = MrcMetaData()
                              .fromGrid(grid_.getGrid())
                              .set_grid_stats(grid_)
                              .to_string();
    fprintf(outputfile, "%s", textDump.c_str());
  }
  if (extension == "pdb") {
    pdbFromGridData(grid_, outputfile);
  }
  if (extension == "dump") {
    dumpFromGridData(grid_, outputfile);
  }

  gmx_fio_fclose(outputfile);
}

MapConverter::MapConverter(const GridDataFloat3D &grid) { grid_ = grid; }

MapConverter::MapConverter() = default;

const GridDataFloat3D &MapConverter::map() { return grid_; }
} // namespace gmx
