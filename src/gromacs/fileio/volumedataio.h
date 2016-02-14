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
 * Reading and writing routines for volume data formats ccp4, mrc and imod.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \inpublicapi
 * \ingroup module_fileio
 */

#ifndef GMX_FILEIO_VOLUMEDATAIO_H_
#define GMX_FILEIO_VOLUMEDATAIO_H_
#include <array>
#include <memory>
#include <string>

#include "gromacs/math/volumedata.h"

namespace gmx
{

/*! \brief Classes for writing, reading, and storing volumetric data.
 */
namespace volumedata
{

/*! \brief
 * A container for the metadata in mrc file formats (compatible with ccp4 and map and mostly imod).
 *
 * For a detailed decription see
 * "EMDB Map Distribution Format Description Version 1.01 (c) emdatabank.org 2014"
 */
struct MrcMetaData{
    bool                       swap_bytes;               //!< swap bytes upon reading/writing (applied, when endianess is different between file and machine architecture)
    int                        mrc_data_mode;            //!< data mode, currently only mode 2 is supported (32-bit float real values)
    int                        machine_stamp;            //!< endianess of map writing architecture (big endian; 0x44410000 , little endian: 0x11110000)
    std::string                format_identifier;        //!< for all density formats: four 1-byte chars reading "MAP " (a little pointless, I know)

    int                        num_labels;               //!< number of used crystallographic labels, 0 for imagestacks, 1 for emdb data
    std::vector < std::string> labels;                   //!< crystallographic labels or \:\:\:\:EMDataBank.org\:\:\:\:EMD-1234\:\:\:\: for EMDB entries

    IVec                       crs_to_xyz;               //!< Axis order
    IVec                       xyz_to_crs;               //!< reversed Axis order
    IVec                       num_crs;                  //!< redundand entry, we use the grid extend (NX,NY,NZ) from header words 8-10
    IVec                       crs_start;                //!< Start of values in grid, typically 0,0,0

    float                      min_value;                //!< minimum voxel value; may be used to scale values in currently unsupported compressed data mode (mrc_data_mode=0)
    float                      max_value;                //!< maximum voxel value; may be used to scale values in currently unsupported compressed data mode (mrc_data_mode=0)
    float                      mean_value;               //!< mean voxel value   (not always reported,as evident from density)
    float                      rms_value;                //!< rms of the density (not always reported,as evident from density)

    bool                       is_crystallographic;      //!< true if crystallographic data is to be read
    bool                       has_skew_matrix;          //!< only crystallographic data: true if skew matrix is stored
    std::array<float, 9>       skew_matrix;              //!< only crystallographic data: skew matrix or, if skew flag is zero, data in place of skew matrix
    RVec                       skew_translation;         //!< only crystallographic data: skew translatation or, if skew flag is zero, data in place of skew translation
    int                        num_bytes_extened_header; //!< only crystallographic data: the size of the symbol table in bytes
    std::vector<char>          extended_header;          //!< only crystallographic data: extended header, usually symbol tables

    std::array<float, 13>      extraskew;                //!< fields unused in EMDB standard, but used for skew matrix and translation in crystallogrphic data (skew flag, skew matrix and skew translation)
    std::array<float, 15>      extra;                    //!< extra data in header, currently unused
};

/*! \brief
 * Read and write real-space real-valued volume data files
 * according to the electron microscopy data bank (EMDB) standard.
 *
 * The formatting guraranties compliance with 3D EM maps described in
 * "EMDB Map Distribution Format Description Version 1.01 (c) emdatabank.org 2014"
 * However, other ccp4, mrc, imod and map formats might be compatible.
 *
 * Future implementations for reading crystallographic data or image stacks might
 * want to split MrcFile into an abstract volumedata base class and respective
 * child implementatons, if demand exists.
 */
class MrcFile
{
    public:
        MrcFile();
        ~MrcFile();

        /*! \brief Write real-spaced, real-valued griddata to file
         * with default metadata for 3D electron microscopy data.
         *
         * \param[in] filename name of the file to write the griddata to, typically *.cpp4, *.mrc or *.map
         * \param[in] grid_data real-valued, real-space data on a grid
         */
        void write(std::string filename, GridReal *grid_data);

        /*! \brief Write real-spaced, real-valued griddata to file with user-defined metadata.
         *
         * \param[in] filename name of the file to write the griddata to, typically *.cpp4, *.mrc or *.map
         * \param[in] grid_data real-valued, real-space data on a grid
         * \param[in] meta struct with own metadata
         * \param[in] bOwnGridStats calculate min, max, mean and rms self if true, otherwise copy from metadata
         */
        void write_with_own_meta(std::string filename, GridReal *grid_data, MrcMetaData *meta, bool bOwnGridStats);

        /*! \brief Reads real-spaced, real-valued griddata from file.
         *
         * \param[in] filename name of the file from which to read the griddata, typically *.cpp4, *.mrc or *.map
         * \param[in] grid_data will be filled with real-valued, real-space data on a grid upon succesful reading; previous content will be discarded
         */
        void read(std::string filename, GridReal *grid_data);
        /*! \brief Reads real-spaced, real-valued griddata into voxel data.
         *
         * \param[in] filename name of the file from which to read the griddata, typically *.cpp4, *.mrc or *.map
         * \param[in] grid_data will be filled with real-valued, real-space data on a grid upon succesful reading; previous content will be discarded
         * \param[in] meta returns the metadata from reading; previous content will be overwritten
         */
        void read_meta(std::string filename, GridReal *grid_data, MrcMetaData *meta);

    private:
        class Impl;
        std::unique_ptr<Impl> impl_;
};

}      // namespace volumedata

}      // namespace gmx
#endif /* end of include guard: GMX_FILEIO_VOLUMEDATAIO_H_ */
