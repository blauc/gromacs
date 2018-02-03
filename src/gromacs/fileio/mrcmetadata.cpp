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
/*! \file
 * \brief
 * Implements mrc header.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \inpublicapi
 * \ingroup module_fileio
 */

#include "mrcmetadata.h"

namespace gmx
{

void MrcMetaData::setEMDBDefaults()
{
    {
        swap_bytes               = false;
        space_group              = 1;
        mrc_data_mode            = 2;
        num_bytes_extened_header = 0;
        has_skew_matrix          = false;
        crs_start                = {{0, 0, 0}};
        crs_to_xyz               = {{0, 1, 2}};
        xyz_to_crs               = {{0, 1, 2}};
        skew_matrix              = {{0, 0, 0, 0, 0, 0, 0, 0, 0}};
        skew_translation         = {0, 0, 0};
        is_crystallographic      = false;
        extra                    = {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
        extraskew                = {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, }};
        format_identifier        = "MAP ";

        // check endianess
        #ifdef GMX_INTEGER_BIG_ENDIAN
        machine_stamp            = 1145110528;
        #else
        machine_stamp            = 4369;
        #endif

        labels                   = {
            {"::::EMDataBank.org::::EMD-xxxx::::Own Data Following EMDB convention::::::::::::"},
            {"                                                                                "},
            {"                                                                                "},
            {"                                                                                "},
            {"                                                                                "},
            {"                                                                                "},
            {"                                                                                "},
            {"                                                                                "},
            {"                                                                                "},
            {"                                                                                "}
        };
        num_labels               = labels.size();
        extended_header          = {};
    }
}

std::string MrcMetaData::to_string() const
{
    std::string result;

    result += "swap_bytes        : " + (swap_bytes == true ? std::string("true") : std::string("false"))    + " | swap bytes upon reading/writing (applied, when endianess is different between file and machine architecture\n";
    result += "mrc_data_mode     : " + std::to_string(mrc_data_mode) + " | data mode, currently only mode 2 is supported (32-bit float real values)\n";
    result += "machine_stamp     : " + std::to_string(machine_stamp) + " | endianess of map writing architecture (big endian " + std::to_string(0x44411111) + " little endian "+ std::to_string(0x11111444) + "\n";
    result += "format_identifier : '"+ format_identifier;
    result += "' | for all density formats: four 1-byte chars reading MAP \n";
    result += "num_labels        : " + std::to_string(num_labels)    + " number of used crystallographic labels, 0 for imagestacks, 1 for em data\n";
    result += "labels            : '" + std::string(labels[0].c_str()) + " crystallographic labels or ::::EMDataBank.org::::EMD-1234:::: for EMDB entries\n";
    result += "crs_to_xyz        : " + std::to_string(crs_to_xyz[0]) + " " + std::to_string(crs_to_xyz[1]) + " " + std::to_string(crs_to_xyz[2]) + " Axis order\n";
    result += "xyz_to_crs        : " + std::to_string(xyz_to_crs[0]) + " " + std::to_string(xyz_to_crs[1]) + " " + std::to_string(xyz_to_crs[2]) + " reversed Axis order\n";
    result += "num_crs           : " + std::to_string(num_crs[0])    + " " + std::to_string(num_crs[1]) + " " + std::to_string(num_crs[2]) + " extend in column row section redundand entry, we use the grid extend (NX,NY,NZ) from header words 8-10\n";
    result += "extend            : " + std::to_string(extend[0])    + " " + std::to_string(extend[1]) + " " + std::to_string(extend[2]) + " grid extend in x, y, z \n";
    result += "crs_start         : " + std::to_string(crs_start[0])  + " " + std::to_string(crs_start[0]) + " " + std::to_string(crs_start[0])+" Start of values in grid, typically 0,0,0\n";
    result += "min_value         : " + std::to_string(min_value)     + " minimum voxel value may be used to scale values in currently unsupported compressed data mode (mrc_data_mode=0)\n";
    result += "max_value         : " + std::to_string(max_value)     + " maximum voxel value may be used to scale values in currently unsupported compressed data mode (mrc_data_mode=0)\n";
    result += "mean_value        : " + std::to_string(mean_value)    + " mean voxel value   (not always reported,as evident from density)\n";
    result += "rms_value         : " + std::to_string(rms_value)     + " rms of the density (not always reported,as evident from density)\n";
    result += "is_crystallographic: " + std::to_string(is_crystallographic) +" true if crystallographic data is to be read\n";
    result += "has_skew_matrix   : " + (has_skew_matrix == true ? std::string("true") : std::string("false"))  + " only crystallographic data: true if skew matrix is stored\n";
    // result+= "skew_matrix: " + std::to_string(skew_matrix)              +" only crystallographic data: skew matrix or, if skew flag is zero, data in place \n";
    // result+= "skew_translation: " + std::to_string(skew_translation)         +" only crystallographic data: skew translatation or, if skew flag is zero, data in place of skew translation\n";
    result += "num_bytes_extened_header : " + std::to_string(num_bytes_extened_header) +" only crystallographic data: the size of the symbol table in bytes\n";
    // result+= "extended_header: " + std::to_string(extended_header)          +" only crystallographic data: extended header, usually symbol tables\n";
    // result+= "extraskew: " + std::to_string(extraskew)                +" fields unused in EMDB standard, but used for skew matrix and translation in crystallogrphic data (skew flag, skew matrix and skew translation)\n";
    // result+= "extra: " + std::to_string(extra)                    +" extra data in header, currently unused\n";
    return result;
}

}      // namespace gmx
