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

/*******************************************************************************
 * MrcFile::Impl
 */
class MrcFile::Impl
{
    public:
        Impl();
        ~Impl();

        /*! \brief Opens file stream in file_.
         * \param[in] filename typically *.mrc, *.map, or *.ccp4
         * \param[in] bRead Open the file for reading if true, otherwise for writing
         */
        void open_file(std::string filename, bool bRead);
        void close_file();
        std::string filetype(std::string filename);
        bool known_extension(std::string filename);
        void read_file_size();

        std::string print_to_string();

        void do_mrc_data_(Field<real> &grid_data, bool bRead);
        void do_mrc_header_(bool bRead);
        FiniteGridWithTranslation<DIM> setFiniteGridWithTranslationFromMrcMeta();
        void setMrcMetaFromFiniteGridWithTranslation(const FiniteGridWithTranslation<DIM> &grid );

        /*! \brief Set the mrc metadata to default values for 3d cryo-EM.
         */
        void set_metadata_mrc_default();

        std::array<int, 3> xyz_to_crs(const std::array<int, 3> &order);
        std::array<int, 3> to_xyz_order(const std::array<int, 3> &i_crs);
        std::array<int, 3> to_crs_order(const std::array<int, 3> &xyz_order);

        /*! \brief Guess, whether endianess differs between input file and reading architecture .
         *
         * If the number of columns in the density file is negative or larger than 65534,
         * assume endianess missmatch between input file and reading machine architecture.*/
        void check_swap_bytes();

        bool has_skew_matrix();
        void set_skew_matrix(matrix skew);


        template <typename T> void do_(T * result, bool bRead)
        {
            if (bRead)
            {
                if (fread(result, sizeof(T), 1, file_) != 1)
                {
                    GMX_THROW(gmx::FileIOError("Cannot read from volume data at " + std::to_string(ftell(file_)) + "."));
                }
            }
            // swap bytes for correct endianness
            if (meta_.swap_bytes)
            {
                if (sizeof(T) == 4)
                {
                    gmx_int32_t int_tmp = gmx_int32_t(*result);
                    *result = (int_tmp & 0xFF000000) >> 24 | (int_tmp & 0x00FF0000) >> 8 | (int_tmp & 0x0000FF00) << 8 | (int_tmp & 0x000000FF) << 24;
                }
                if (sizeof(T) == 2)
                {
                    int16_t int_tmp = int16_t(*result);
                    *result = (int_tmp & 0xFF00) >> 8 | (int_tmp & 0x00FF) << 8;
                }
            }

            if (!bRead)
            {
                fwrite(result, sizeof(T), 1, file_);
            }
        }

        void do_float32_rvec_(RVec * result, bool bRead);
        void do_int32_ivec_(std::array<int, 3> * result, bool bRead);

        void set_grid_stats(const Field<real> &grid_data);

        bool colummn_row_section_order_valid_(std::array<int, 3> crs_to_xyz);

        void swap_int32_(gmx_int32_t *result);
        void swap_float32_(float *result);

        FILE     *file_;
        long      file_size_;
        const int size_extrarecord;
        const int number_labels;
        const int label_size;
        const int header_bytes;
        const std::vector<std::string> filetypes;

        MrcMetaData                    meta_;
};

std::string MrcFile::Impl::print_to_string()
{
    meta_.labels[0].erase(std::remove_if(meta_.labels[0].begin(), meta_.labels[0].end(), [](int i){return i == '\000'; }), meta_.labels[0].end());
    std::string result;
    result += "\n---MrcFile Info--file size:" + std::to_string(file_size_) + "---\n";

    result += "swap_bytes        : " + (meta_.swap_bytes == true ? std::string("true") : std::string("false"))    + " | swap bytes upon reading/writing (applied, when endianess is different between file and machine architecture\n";
    result += "mrc_data_mode     : " + std::to_string(meta_.mrc_data_mode) + " | data mode, currently only mode 2 is supported (32-bit float real values)\n";
    result += "machine_stamp     : " + std::to_string(meta_.machine_stamp) + " | endianess of map writing architecture (big endian " + std::to_string(0x44411111) + " little endian "+ std::to_string(0x11111444) + "\n";
    result += "format_identifier : '"+ meta_.format_identifier;
    result += "' | for all density formats: four 1-byte chars reading MAP \n";
    result += "num_labels        : " + std::to_string(meta_.num_labels)    + " number of used crystallographic labels, 0 for imagestacks, 1 for em data\n";
    result += "labels            : '" + std::string(meta_.labels[0].c_str()) + " crystallographic labels or ::::EMDataBank.org::::EMD-1234:::: for EMDB entries\n";
    result += "crs_to_xyz        : " + std::to_string(meta_.crs_to_xyz[0]) + " " + std::to_string(meta_.crs_to_xyz[1]) + " " + std::to_string(meta_.crs_to_xyz[2]) + " Axis order\n";
    result += "xyz_to_crs        : " + std::to_string(meta_.xyz_to_crs[0]) + " " + std::to_string(meta_.xyz_to_crs[1]) + " " + std::to_string(meta_.xyz_to_crs[2]) + " reversed Axis order\n";
    result += "num_crs           : " + std::to_string(meta_.num_crs[0])    + " " + std::to_string(meta_.num_crs[1]) + " " + std::to_string(meta_.num_crs[2]) + " extend in column row section redundand entry, we use the grid extend (NX,NY,NZ) from header words 8-10\n";
    result += "extend            : " + std::to_string(meta_.extend[0])    + " " + std::to_string(meta_.extend[1]) + " " + std::to_string(meta_.extend[2]) + " grid extend in x, y, z \n";
    result += "crs_start         : " + std::to_string(meta_.crs_start[0])  + " " + std::to_string(meta_.crs_start[0]) + " " + std::to_string(meta_.crs_start[0])+" Start of values in grid, typically 0,0,0\n";
    result += "min_value         : " + std::to_string(meta_.min_value)     + " minimum voxel value may be used to scale values in currently unsupported compressed data mode (mrc_data_mode=0)\n";
    result += "max_value         : " + std::to_string(meta_.max_value)     + " maximum voxel value may be used to scale values in currently unsupported compressed data mode (mrc_data_mode=0)\n";
    result += "mean_value        : " + std::to_string(meta_.mean_value)    + " mean voxel value   (not always reported,as evident from density)\n";
    result += "rms_value         : " + std::to_string(meta_.rms_value)     + " rms of the density (not always reported,as evident from density)\n";
    result += "is_crystallographic: " + std::to_string(meta_.is_crystallographic) +" true if crystallographic data is to be read\n";
    result += "has_skew_matrix   : " + (meta_.has_skew_matrix == true ? std::string("true") : std::string("false"))  + " only crystallographic data: true if skew matrix is stored\n";
    // result+= "skew_matrix: " + std::to_string(meta_.skew_matrix)              +" only crystallographic data: skew matrix or, if skew flag is zero, data in place \n";
    // result+= "skew_translation: " + std::to_string(meta_.skew_translation)         +" only crystallographic data: skew translatation or, if skew flag is zero, data in place of skew translation\n";
    result += "num_bytes_extened_header : " + std::to_string(meta_.num_bytes_extened_header) +" only crystallographic data: the size of the symbol table in bytes\n";
    // result+= "extended_header: " + std::to_string(meta_.extended_header)          +" only crystallographic data: extended header, usually symbol tables\n";
    // result+= "extraskew: " + std::to_string(meta_.extraskew)                +" fields unused in EMDB standard, but used for skew matrix and translation in crystallogrphic data (skew flag, skew matrix and skew translation)\n";
    // result+= "extra: " + std::to_string(meta_.extra)                    +" extra data in header, currently unused\n";
    return result;
}

FiniteGridWithTranslation<DIM> MrcFile::Impl::setFiniteGridWithTranslationFromMrcMeta()
{
    RVec            cell_length;
    svmul(A2NM, meta_.cell_length, cell_length);

    FiniteGrid<DIM>                tmpgrid({{{cell_length[XX], cell_length[YY], cell_length[ZZ]}}}, meta_.extend);
    auto                           translation =  tmpgrid.multiIndexToCoordinate( {{ meta_.crs_start[meta_.xyz_to_crs[XX]], meta_.crs_start[meta_.xyz_to_crs[YY]], meta_.crs_start[meta_.xyz_to_crs[ZZ]] }});
    FiniteGridWithTranslation<DIM> result({{{cell_length[XX], cell_length[YY], cell_length[ZZ]}}}, meta_.extend, translation);
    /* If this the map origin is shifted, because the grid indexing starts at other values than zero,
     * values read here are ignored.
     * Convention is not clear at this point, whether the translation here should be treated as extra shift,
     * in observed cases this was not the case.
     *
     * Silently ignore if the map translation due to grid-index start shift
     * does not match the shift reported here.
     */
    if (!meta_.is_crystallographic && meta_.crs_start[XX] == 0 && meta_.crs_start[YY] == 0 && meta_.crs_start[ZZ] == 0)
    {
        result.setTranslation({{(float) A2NM * meta_.extra[size_extrarecord-3], (float) A2NM * meta_.extra[size_extrarecord-2], (float) A2NM * meta_.extra[size_extrarecord-1]}});
    }

    return result;
}

void MrcFile::Impl::setMrcMetaFromFiniteGridWithTranslation(const FiniteGridWithTranslation<DIM> &grid )
{
    auto index_of_origin = grid.coordinateToFloorMultiIndex({{1e-6, 1e-6, 1e-6}});
    meta_.crs_start   = {{-index_of_origin[XX], -index_of_origin[YY], -index_of_origin[ZZ]}};
    meta_.num_crs     = to_crs_order(grid.lattice().getExtend());
    meta_.extend      = grid.lattice().getExtend();
    meta_.cell_angles = { 90, 90, 90 };

    for (int dimension = 0; dimension <= ZZ; ++dimension)
    {
        meta_.cell_length[dimension] = NM2A * grid.cell().basisVectorLength(dimension);
    }
    meta_.num_labels = meta_.labels.size();
}

std::array<int, 3> MrcFile::Impl::to_xyz_order(const std::array<int, 3> &i_crs)
{
    std::array<int, 3> i_xyz;

    i_xyz[meta_.crs_to_xyz[XX]] = i_crs[XX];
    i_xyz[meta_.crs_to_xyz[YY]] = i_crs[YY];
    i_xyz[meta_.crs_to_xyz[ZZ]] = i_crs[ZZ];

    return i_xyz;
}


void MrcFile::Impl::do_int32_ivec_(std::array<int, 3> * result, bool bRead)
{
    do_<gmx_int32_t>(&(*result)[XX], bRead);
    do_<gmx_int32_t>(&(*result)[YY], bRead);
    do_<gmx_int32_t>(&(*result)[ZZ], bRead);
}

void MrcFile::Impl::do_float32_rvec_(RVec * result, bool bRead)
{
    do_<float>(&((*result)[XX]), bRead);
    do_<float>(&((*result)[YY]), bRead);
    do_<float>(&((*result)[ZZ]), bRead);
}


std::string MrcFile::Impl::filetype(std::string filename)
{
    std::string result = "";
    size_t      pos    = filename.rfind(".");
    if (pos == 0 || pos == std::string::npos)
    {
        return result;
    }
    else
    {
        return filename.substr(pos+1);
    }
}

std::array<int, 3> MrcFile::Impl::xyz_to_crs(const std::array<int, 3> &order)
{
    std::array<int, 3> result;
    result[order[XX]]  = XX;
    result[order[YY]]  = YY;
    result[order[ZZ]]  = ZZ;
    return result;
}

std::array<int, 3> MrcFile::Impl::to_crs_order(const std::array<int, 3> &xyz_order)
{
    std::array<int, 3> result;
    result[meta_.xyz_to_crs[XX]] = xyz_order[XX];
    result[meta_.xyz_to_crs[YY]] = xyz_order[YY];
    result[meta_.xyz_to_crs[ZZ]] = xyz_order[ZZ];
    return result;
}

bool MrcFile::Impl::colummn_row_section_order_valid_(std::array<int, 3> crs_to_xyz)
{
    const std::set<int> valid_crs_set {
        0, 1, 2
    };
    std::set<int> crs_set {
        crs_to_xyz[XX], crs_to_xyz[YY], crs_to_xyz[ZZ]
    };
    return valid_crs_set == crs_set;
};

void MrcFile::Impl::do_mrc_header_(bool bRead)
{
    if (bRead)
    {
        check_swap_bytes();
        read_file_size();
        if (file_size_ < header_bytes)
        {
            GMX_THROW(gmx::FileIOError("Density format file is smaller than expected header."));
        }
    }

    /* Supports reading according to
       ftp://ftp.wwpdb.org/pub/emdb/doc/Map-format/current/EMDB_map_format.pdf
       note, that
       http://www.ccpem.ac.uk/mrc_format/mrc2014.php
       differs slightly in definition */

    /* 1-3 | NC, NR, NS | signed int >0
     * # of columns (fastest changing),rows, sections (slowest changing)
     * emdb convention: NC=NR=NS                     */

    do_int32_ivec_(&meta_.num_crs, bRead);

    /* 4   | MODE | signed int | 0,1,2,3,4
     * voxel datatype
     * emdb convention: 2       */
    do_<gmx_int32_t>(&meta_.mrc_data_mode, bRead);

    /* MODE = 0: 8 bits, density stored as a signed byte (range -128 to 127, ISO/IEC 10967)
     * MODE = 1: 16 bits, density stored as a signed integer (range -32768 to 32767, ISO/IEC 10967)
     * MODE = 2: 32 bits, density stored as a floating point number (IEEE 754)
     * MODE = 3: 32 bits, Fourier transform stored as complex signed integers (ISO/IEC 10967)
     * MODE = 4: 64 bits, Fourier transform stored as complex floating point numbers (IEEE 754)     */
    if (meta_.mrc_data_mode < 0 || meta_.mrc_data_mode > 4)
    {
        GMX_THROW(gmx::FileIOError("Read invalid mrc/cpp4/imod data mode. Mode " + std::to_string(meta_.mrc_data_mode) + " not in [0..4] ."));
    }
    if (meta_.mrc_data_mode != 2)
    {
        GMX_THROW(gmx::NotImplementedError("Other mrc/ccp4/imod data modes than 32 bit float not currently implemented. Upgrading to the newest gromacs verison might possibly help."));
    }

    /* 5-7 | NCSTART, NRSTART, NSSTART | signed int
     * position of first column, first row, and first section (voxel grid units)
     *
     * The position of the first voxel is defined in grid units by NCSTART, NRSTART, and NSSTART.
     * The center of the voxel with grid position (0,0,0) corresponds to the Cartesian coordinate origin.*/
    do_int32_ivec_(&meta_.crs_start, bRead);

    /* 8-10 | NX, NY, NZ | signed int >0 |
     * intervals per unit cell repeat along X,Y Z
     * intervals per map length along X,Y,Z;
     * emdb convention: same as NC, NR, NS
     *
     * Lengths in Aangstroms for a single voxel are as follows:
     * Xvoxel = X_LENGTH/NX Yvoxel = Y_LENGTH/NY Zvoxel = Z_LENGTH/NZ */

    do_int32_ivec_(&meta_.extend, bRead);

    /* 11-13 | X_LENGTH, Y_LENGTH, Z_LENGTH | floating pt >0
     * Unit Cell repeats along X, Y, Z In Aangstrom
     * emdb Map lengths along X,Y,Z in Aangstrom
     *
     * Lengths in Aangstroms for a single voxel are as follows:
     * Xvoxel = X_LENGTH/NX Yvoxel = Y_LENGTH/NY Zvoxel = Z_LENGTH/NZ */
    do_float32_rvec_(&meta_.cell_length, bRead);

    /* 14-16 | ALPHA,BETA,GAMMA | floating pt >0, <180
     * Unit Cell angles (degrees)
     * emdb convention: 90, 90, 90
     *
     * By convention, cell angles (ALPHA, BETA, GAMMA)
     * are 90 degrees for single particle or tomogram EM maps;
     * they follow IUCr space group conventions for crystals.*/
    do_float32_rvec_(&meta_.cell_angles, bRead);
    // By convention, unset cell angles (all 0) are interpreted as 90 deg.
    if (meta_.cell_angles[XX]*meta_.cell_angles[YY]*meta_.cell_angles[ZZ] < 1e-5)
    {
        meta_.cell_angles = {90, 90, 90};
    }

    /* 17-19 | MAPC, MAPR, MAPS | signed int | 1 (=X) 2 (=Y) 3 (=Z)
     * relationship of X,Y,Z axes to columns, rows, sections
     * emdb convention: 1, 2, 3 */
    std::array<int, 3> crs_to_xyz {{
                                       meta_.crs_to_xyz[XX]+1, meta_.crs_to_xyz[YY]+1, meta_.crs_to_xyz[ZZ]+1
                                   }};
    do_int32_ivec_(&crs_to_xyz, bRead);

    if (bRead)
    {
        meta_.crs_to_xyz[XX] = crs_to_xyz[XX]-1;
        meta_.crs_to_xyz[YY] = crs_to_xyz[YY]-1;
        meta_.crs_to_xyz[ZZ] = crs_to_xyz[ZZ]-1;
        if (!colummn_row_section_order_valid_(meta_.crs_to_xyz))
        {
            meta_.crs_to_xyz = {{0, 1, 2}};
        }
        meta_.xyz_to_crs = xyz_to_crs(meta_.crs_to_xyz);
    }


    /* 20-22 | AMIN, AMAX, AMEAN | floating pt
     * Minimum, maximum, average density */
    do_<float>(&(meta_.min_value ), bRead);
    do_<float>(&(meta_.max_value ), bRead);
    do_<float>(&(meta_.mean_value), bRead);

    /* 23 | ISPG | signed int 1-230 |
     * space group #
     * emdb convention 1
     *
     * Space Group Numbers are defined by IUCr conventions
     * (Table 12.3.4.1 Standard space-group symbols”, pages 824-831,
     * International Tables for Crystallography, Volume A, fifth edition).
     *
     * For 3D volumes of single particle or tomogram entries, ISPG=1 and NSYMBT=0.
     * For image stacks ISPG = 0 */
    do_<gmx_int32_t>(&meta_.space_group, bRead);

    /* 24 | NSYMBT | signed int | 80n
     * # of bytes in symmetry table (multiple of 80)
     * emdb convention 0 */
    do_<gmx_int32_t>(&meta_.num_bytes_extened_header, bRead);
    if (meta_.num_bytes_extened_header%80 != 0)
    {
        GMX_THROW(gmx::FileIOError("Read invalid number of bytes in symbol table from mrc/cpp4/imod file. Should be 80, but is " + std::to_string(meta_.num_bytes_extened_header) + "instead."));
    }

    if (meta_.is_crystallographic)
    {
        /* 25 | LSKFLG | signed int | 0,1
         * flag for skew matrix
         * emdb convention 0 */
        gmx_int32_t hasSkewMatrix = meta_.has_skew_matrix ? 1 : 0;
        do_<gmx_int32_t>(&hasSkewMatrix, bRead);
        if (bRead)
        {
            if (!(hasSkewMatrix == 0 || hasSkewMatrix == 1))
            {
                GMX_THROW(gmx::FileIOError("Skew matrix flag set to invalid value in mrc/cpp4/imod file. Should be 0 or 1 but is " + std::to_string(hasSkewMatrix) + "instead."));
            }
            meta_.has_skew_matrix = (hasSkewMatrix == 1) ? true : false;
        }

        if (meta_.has_skew_matrix)
        {
            /* TODO: A2NM conversion for skew matrix if necessary */
            /* 26-34 | SKWMAT | floating pt
             * skew matrix-S11, S12, S13, S21, S22, S23, S31, S32, S33
             * emdb convention: not set
             *
             * SKWMAT, SKWTRN, and EXTRA fields are not currently used by EMDB. */

            /* 35-37 | SKWTRN | floating pt
             * skew translation-T1, T2, T3
             * emdb convention: not set
             *
             * SKWMAT, SKWTRN, and EXTRA fields are not currently used by EMDB. */
            for (auto && i : meta_.skew_matrix)
            {
                do_<float>(&i, bRead);
            }
            do_float32_rvec_(&meta_.skew_translation, bRead);
        }
    }
    else
    {
        /* 25-37 not used in EMDB */
        for (auto && i : meta_.extraskew)
        {
            do_<float>(&i, bRead);
        }
    }

    /* 38-52 | EXTRA | 32 bit binary
     * user-defined metadata
     *
     * SKWMAT, SKWTRN, and EXTRA fields are not currently used by EMDB.
     * EMDB might use fields 50,51 and 52 for setting the coordinate system origin */
    for (auto && i : meta_.extra)
    {
        do_<float>(&i, bRead);
    }

    /* 53 | MAP | ASCII char
     * "MAP "
     * MRC/CCP4 MAP format identifier */
    meta_.format_identifier.resize(4);
    for (size_t i = 0; i < 4; i++)
    {
        do_<char>(&meta_.format_identifier[i], bRead);
    }
    if (!(meta_.format_identifier.compare("MAP ") == 0))
    {
        fprintf(stderr, "\nWARNING: Expected " "MAP " " as format identifier.\n");
    }
    /* 54 | MACHST | 32 bit
     * binary machine stamp
     *
     * MACHST is (written/read as 4 hex byte sequence)
     * 0x44,0x41,0x00,0x00  for little endian machines
     * 0x11,0x11,0x00,0x00  for big endian machines
     */
    do_<gmx_int32_t>(&(meta_.machine_stamp), bRead);

    /* 55 | RMS | floating pt
     * Density root-mean-square deviation */
    do_<float>(&(meta_.rms_value), bRead);

    /* 56 | NLABL | signed int | 0-10
     * # of labels
     *
     * Following the 2010 remediation, maps distributed by EMDB
     * now have a single label of form “::::EMDataBank.org::::EMD-1234::::”.  */
    do_<gmx_int32_t>(&meta_.num_labels, bRead);

    /* 57-256 | LABEL_N | ASCII char
     * Up to 10 user-defined labels */
    if (bRead)
    {
        std::string mrc_label(label_size, ' ');
        meta_.labels.clear();
        for (int i = 0; i < number_labels; i++)
        {
            int read_size = fread(&mrc_label[0], 1, label_size, file_);
            if (read_size != label_size)
            {
                GMX_THROW(gmx::FileIOError("Could not read label from file."));
            }
            meta_.labels.push_back(std::string(mrc_label));
        }
    }
    else
    {
        for (auto label : meta_.labels)
        {
            fprintf(file_, "%.*s", label_size, label.c_str());
        }
    }

    /* 257-257+NSYMBT | anything
     */
    if (bRead)
    {
        if (file_size_ < meta_.num_bytes_extened_header + header_bytes)
        {
            GMX_THROW(gmx::FileIOError("Density format file is smaller than expected extended header."));
        }
        meta_.extended_header.resize(meta_.num_bytes_extened_header);
        if (fgets(&((meta_.extended_header)[0]), meta_.num_bytes_extened_header, file_) != nullptr)
        {
            GMX_THROW(gmx::FileIOError("Cannot read extended header from file."));
        }
    }
    else
    {
        fwrite(&((meta_.extended_header)[0]), sizeof(char), meta_.extended_header.size(), file_);
    }

};

void MrcFile::Impl::swap_int32_(gmx_int32_t *result)
{

    gmx_int32_t src = *result;
    *result = 0;

    *result |= (src & 0xFF000000) >> 24;
    *result |= (src & 0x00FF0000) >> 8;
    *result |= (src & 0x0000FF00) << 8;
    *result |= (src & 0x000000FF) << 24;
}

void MrcFile::Impl::do_mrc_data_(Field<real> &grid_data, bool bRead)
{
    const auto &lattice = grid_data.getGrid().lattice();
    if (bRead)
    {
        if (file_size_ < header_bytes + meta_.num_bytes_extened_header + long(lattice.getNumLatticePoints()*sizeof(float)) )
        {
            GMX_THROW(gmx::FileIOError("Density format file is smaller than indicated in its header."));
        }
        if (file_size_ > header_bytes + meta_.num_bytes_extened_header +  long(lattice.getNumLatticePoints()*sizeof(float)) )
        {
            fprintf(stderr, "WARNING : Density format file size is %ld, however header (%d) + symbol table (%d) + data  (%ld) is larger than indicated in its header. Reading anyway.. \n", file_size_, header_bytes, meta_.num_bytes_extened_header, lattice.getNumLatticePoints()*sizeof(float));
        }
    }

    auto num_crs = to_crs_order(lattice.getExtend());
    for (int section = 0; section  < num_crs[ZZ]; section++)
    {
        for (int row  = 0; row  < num_crs[YY]; row++)
        {
            for (int column  = 0; column  < num_crs[XX]; column++)
            {
                do_<float>(&(grid_data.atMultiIndex(to_xyz_order({{column, row, section}}))), bRead);
            }
        }
    }
}

void MrcFile::Impl::read_file_size()
{
    fpos_t current;
    fgetpos(file_, &current);
    fseek(file_, 0, SEEK_END);
    file_size_ = ftell(file_);
    fsetpos(file_, &current);
}

void MrcFile::Impl::check_swap_bytes()
{

    fpos_t      current;
    meta_.swap_bytes = false;
    gmx_int32_t number_columns;
    fgetpos(file_, &current);

    fseek(file_, 0, SEEK_SET);
    do_<gmx_int32_t>(&number_columns, true);
    if (number_columns <= 0 || number_columns >= 65536)
    {
        meta_.swap_bytes = true;
    }

    // rewind the file
    fsetpos(file_, &current);

}

bool MrcFile::Impl::known_extension(std::string filename)
{
    return std::find(filetypes.begin(), filetypes.end(), filetype(filename)) != filetypes.end();
}

void MrcFile::Impl::open_file(std::string filename, bool bRead)
{

    if (filename.empty())
    {
        GMX_THROW(gmx::FileIOError("Filename empty."));
    }

    filename.erase(std::remove(filename.begin(), filename.end(), '\n'), filename.end());

    if (bRead)
    {
        if (!known_extension(filename) )
        {
            GMX_THROW(gmx::FileIOError("Cannot read filetype " + filetype(filename) + "."));
        }
        file_ = gmx_fio_fopen(filename.c_str(), "r");
    }
    else
    {
        file_ = gmx_fio_fopen(filename.c_str(), "w");
    }
    if (file_ == nullptr)
    {
        GMX_THROW(gmx::FileIOError("Cannot open " + filename + ". "));
    }

}

void MrcFile::Impl::set_metadata_mrc_default()
{
    meta_.swap_bytes               = false;
    meta_.space_group              = 1;
    meta_.mrc_data_mode            = 2;
    meta_.num_bytes_extened_header = 0;
    meta_.has_skew_matrix          = false;
    meta_.crs_start                = {{0, 0, 0}};
    meta_.crs_to_xyz               = {{0, 1, 2}};
    meta_.xyz_to_crs               = {{0, 1, 2}};
    meta_.skew_matrix              = {{0, 0, 0, 0, 0, 0, 0, 0, 0}};
    meta_.skew_translation         = {0, 0, 0};
    meta_.is_crystallographic      = false;
    meta_.extra                    = {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
    meta_.extraskew                = {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, }};
    meta_.format_identifier        = "MAP ";

    // check endianess
#if GMX_INTEGER_BIG_ENDIAN
    meta_.machine_stamp            = 1145110528;
#else
    meta_.machine_stamp            = 4369;
#endif

    meta_.labels                   = {
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
    meta_.num_labels               = meta_.labels.size();
    meta_.extended_header          = {};
}

void MrcFile::Impl::close_file()
{
    gmx_fio_fclose(file_);
}

MrcFile::Impl::Impl() : file_(nullptr), file_size_(0), size_extrarecord(15),
                        number_labels(10), label_size(80), header_bytes(1024),
                        filetypes({"mrc", "ccp4", "imod", "map"}
                                  )
{
    MrcFile::Impl::set_metadata_mrc_default();
};

MrcFile::Impl::~Impl()
{
    if (file_ != nullptr)
    {
        close_file();
    }
};

void
MrcFile::Impl::set_grid_stats(const Field<real> &grid_data)
{
    auto properties = RealFieldMeasure(grid_data);
    meta_.min_value  = properties.min();
    meta_.max_value  = properties.max();
    meta_.mean_value = properties.mean();
    meta_.rms_value  = properties.rms();
}


/*******************************************************************************
 * MrcFile
 */
MrcFile::MrcFile() : impl_(new MrcFile::Impl)
{
}

MrcFile::~MrcFile()
{
    impl_->close_file();
}

std::string MrcFile::print_to_string()
{
    return impl_->print_to_string();
}

void MrcFile::write_with_own_meta(std::string filename, Field<real> &grid_data, MrcMetaData &meta, bool bOwnGridStats)
{
    bool bRead = false;

    impl_->meta_ = meta;

    if (bOwnGridStats)
    {
        impl_->set_grid_stats(grid_data);
    }

    impl_->open_file(filename, bRead);
    impl_->do_mrc_header_(bRead);
    impl_->do_mrc_data_(*(const_cast<Field<real>*>(&grid_data)), bRead);
    impl_->close_file();
}

void MrcFile::write(std::string filename, const Field<real> &grid_data)
{
    bool bRead = false;
    impl_->set_grid_stats(grid_data);
    impl_->open_file(filename, bRead);
    impl_->setMrcMetaFromFiniteGridWithTranslation(grid_data.getGrid());
    impl_->do_mrc_header_(bRead);
    impl_->do_mrc_data_(*(const_cast<Field<real>*>(&grid_data)), bRead);
    impl_->close_file();
}


void MrcFile::read_meta(std::string filename, MrcMetaData &meta)
{
    bool        bRead = true;
    impl_->open_file(filename, bRead);
    impl_->do_mrc_header_(bRead);
    impl_->close_file();
    meta = impl_->meta_;
}

Field<real> MrcFile::read_with_meta(std::string filename, MrcMetaData &meta)
{
    auto result = read(filename);
    meta = impl_->meta_;
    return result;
}

Field<real> MrcFile::read(std::string filename)
{
    bool bRead = true;
    impl_->open_file(filename, bRead);
    impl_->do_mrc_header_(bRead);
    auto grid   = impl_->setFiniteGridWithTranslationFromMrcMeta();
    auto result = Field<real>(grid);
    impl_->do_mrc_data_(result, bRead);
    impl_->close_file();
    return result;
}

Df3File::SuccessfulDf3Write
Df3File::write(std::string filename, const gmx::Field<real> &grid_data)
{
    const auto &grid  = grid_data.getGrid();
    auto        file_ = gmx_fio_fopen(filename.c_str(), "w");
    int16_t     xExtendShort {
        int16_t(grid.lattice().getExtend()[XX])
    };
    int16_t yExtendShort {
        int16_t(grid.lattice().getExtend()[YY])
    };
    int16_t zExtendShort {
        int16_t(grid.lattice().getExtend()[ZZ])
    };
    fputc(xExtendShort >> 8, file_);
    fputc(xExtendShort & 0xff, file_);
    fputc(yExtendShort >> 8, file_);
    fputc(yExtendShort & 0xff, file_);
    fputc(zExtendShort >> 8, file_);
    fputc(zExtendShort & 0xff, file_);
    auto measure = RealFieldMeasure(grid_data);
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

Df3File::SuccessfulDf3Write::SuccessfulDf3Write(std::string filename, const gmx::Field<real> &grid_data) : filename_ {filename}, gridData_ {
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
