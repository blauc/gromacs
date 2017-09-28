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
 */
#include "gmxpre.h"

#include "volumedataio.h"

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
#include "gromacs/math/volumedata/operations/realfieldmeasure.h"

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
        void read_format_identifier();

        std::string print_to_string();

        void do_mrc(const Field<real> &grid_data, bool bRead);
        void do_mrc_data_(Field<real> &grid_data, bool bRead);
        void do_mrc_header_(Field<real> &grid_data, bool bRead);

        void set_mrc_data_mode(int mode);

        /*! \brief Set the mrc metadata to default values for 3d cryo-EM.
         */
        void set_metadata_mrc_default();
        void set_num_bytes_extened_header(int n);

        IVec xyz_to_crs(IVec order);
        IVec to_xyz_order(IVec i_crs);
        IVec to_crs_order(IVec xyz_order);
        void set_crs_to_xyz(IVec order);

        /*! \brief Guess, whether endianess differs between input file and reading architecture .
         *
         * If the number of columns in the density file is negative or larger than 65534,
         * assume endianess missmatch between input file and reading machine architecture.*/
        void check_swap_bytes();

        void set_labels();
        void write_labels();
        void set_extended_header();
        void write_extended_header();

        bool is_crystallographic();
        bool has_skew_matrix();
        void set_has_skew_matrix(int flag);
        void set_skew_matrix(matrix skew);

        char read_char_();
        gmx_int32_t read_int32_();
        float read_float32_();
        RVec read_float32_rvec_();
        IVec read_int32_ivec_();
        void write_int32_(int data);
        void write_int32_ivec_(const IVec &i);
        void write_float32_(real data);
        void write_float32_rvec_(const RVec &i);

        void set_meta(const Field<real> &grid_data);
        void set_grid_stats(const Field<real> &grid_data);


        bool colummn_row_section_order_valid_(IVec crs_to_xyz);

        void swap_int32_(gmx_int32_t &result);
        void swap_float32_(float &result);

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
    result += "format_identifier : '";
    result.push_back(meta_.format_identifier[0]);
    result.push_back(meta_.format_identifier[1]);
    result.push_back(meta_.format_identifier[2]);
    result.push_back(meta_.format_identifier[3]);
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

IVec MrcFile::Impl::xyz_to_crs(IVec order)
{
    IVec result;
    result[order[XX]]  = XX;
    result[order[YY]]  = YY;
    result[order[ZZ]]  = ZZ;
    return result;
}

void MrcFile::Impl::set_crs_to_xyz(IVec order)
{
    meta_.crs_to_xyz = order;
    meta_.xyz_to_crs = xyz_to_crs(order);
}

IVec MrcFile::Impl::to_xyz_order(IVec i_crs)
{
    IVec i_xyz;

    i_xyz[meta_.crs_to_xyz[XX]] = i_crs[XX];
    i_xyz[meta_.crs_to_xyz[YY]] = i_crs[YY];
    i_xyz[meta_.crs_to_xyz[ZZ]] = i_crs[ZZ];

    return i_xyz;
}


IVec MrcFile::Impl::read_int32_ivec_()
{
    IVec result;
    result[0] = read_int32_();
    result[1] = read_int32_();
    result[2] = read_int32_();
    return result;
}

RVec MrcFile::Impl::read_float32_rvec_()
{
    RVec result;
    result[0] = read_float32_();
    result[1] = read_float32_();
    result[2] = read_float32_();
    return result;
}


bool MrcFile::Impl::is_crystallographic()
{
    return meta_.is_crystallographic;
}

void MrcFile::Impl::write_float32_rvec_(const RVec &i)
{
    write_float32_(static_cast<float>(i[XX]));
    write_float32_(static_cast<float>(i[YY]));
    write_float32_(static_cast<float>(i[ZZ]));
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

void MrcFile::Impl::write_extended_header()
{
    fwrite(&((meta_.extended_header)[0]), sizeof(char), meta_.extended_header.size(), file_);
}

void MrcFile::Impl::set_labels()
{
    std::string mrc_label(label_size, ' ');
    meta_.labels.clear();
    int         read_size;
    for (int i = 0; i < number_labels; i++)
    {
        read_size = fread(&mrc_label[0], 1, label_size, file_);
        if (read_size != label_size)
        {
            GMX_THROW(gmx::FileIOError("Could not read label from file."));
        }
        meta_.labels.push_back(std::string(mrc_label));
    }
}

void MrcFile::Impl::write_labels()
{
    for (auto label : meta_.labels)
    {
        fprintf(file_, "%.*s", label_size, label.c_str());
    }
}

void MrcFile::Impl::set_extended_header()
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
    ;
}

bool MrcFile::Impl::has_skew_matrix()
{
    return meta_.has_skew_matrix;
}

void MrcFile::Impl::set_has_skew_matrix(int flag)
{
    if (!(flag == 0 || flag == 1))
    {
        GMX_THROW(gmx::FileIOError("Skew matrix flag set to invalid value in mrc/cpp4/imod file. Should be 0 or 1 but is " + std::to_string(flag) + "instead."));
    }
    if (flag == 1)
    {
        meta_.has_skew_matrix = true;
    }
    if (flag == 0)
    {
        meta_.has_skew_matrix = false;
    }
}

void MrcFile::Impl::set_num_bytes_extened_header(int n)
{
    if (n%80 != 0)
    {
        GMX_THROW(gmx::FileIOError("Read invalid number of bytes in symbol table from mrc/cpp4/imod file. Should be 80, but is " + std::to_string(n) + "instead."));
    }
    meta_.num_bytes_extened_header = n;

};

void MrcFile::Impl::set_mrc_data_mode(int mode)
{
/* MODE = 0: 8 bits, density stored as a signed byte (range -128 to 127, ISO/IEC 10967)
 * MODE = 1: 16 bits, density stored as a signed integer (range -32768 to 32767, ISO/IEC 10967)
 * MODE = 2: 32 bits, density stored as a floating point number (IEEE 754)
 * MODE = 3: 32 bits, Fourier transform stored as complex signed integers (ISO/IEC 10967)
 * MODE = 4: 64 bits, Fourier transform stored as complex floating point numbers (IEEE 754)     */
    if (mode < 0 || mode > 4)
    {
        GMX_THROW(gmx::FileIOError("Read invalid mrc/cpp4/imod data mode. Mode " + std::to_string(mode) + " not in [0..4] ."));
    }
    if (mode != 2)
    {
        GMX_THROW(gmx::NotImplementedError("Other mrc/ccp4/imod data modes than 32 bit float not currently implemented. Upgrading to the newest gromacs verison might possibly help."));
    }
    meta_.mrc_data_mode = mode;
}

void MrcFile::Impl::read_format_identifier()
{
    meta_.format_identifier.clear();
    meta_.format_identifier.push_back(read_char_());
    meta_.format_identifier.push_back(read_char_());
    meta_.format_identifier.push_back(read_char_());
    meta_.format_identifier.push_back(read_char_());
    if (!(meta_.format_identifier.compare("MAP ") == 0))
    {
        fprintf(stderr, "\nWARNING: Expected " "MAP " " as format identifier.\n");
    }

}

IVec MrcFile::Impl::to_crs_order(IVec xyz_order)
{
    IVec result;
    result[meta_.xyz_to_crs[XX]] = xyz_order[XX];
    result[meta_.xyz_to_crs[YY]] = xyz_order[YY];
    result[meta_.xyz_to_crs[ZZ]] = xyz_order[ZZ];
    return result;
}

void MrcFile::Impl::write_int32_ivec_(const IVec &i)
{
    write_int32_(i[XX]);
    write_int32_(i[YY]);
    write_int32_(i[ZZ]);
}

bool MrcFile::Impl::colummn_row_section_order_valid_(IVec crs_to_xyz)
{
    const std::set<int> valid_crs_set {
        0, 1, 2
    };
    std::set<int> crs_set {
        crs_to_xyz[XX], crs_to_xyz[YY], crs_to_xyz[ZZ]
    };
    return valid_crs_set == crs_set;
};

void MrcFile::Impl::do_mrc_header_(Field<real> &grid_data, bool bRead)
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
    if (bRead)
    {
        meta_.num_crs = read_int32_ivec_();
    }
    else
    {
        write_int32_ivec_(to_crs_order(grid_data.getExtend()));
    }

    /* 4   | MODE | signed int | 0,1,2,3,4
     * voxel datatype
     * emdb convention: 2       */
    if (bRead)
    {
        set_mrc_data_mode(read_int32_());
    }
    else
    {
        write_int32_(meta_.mrc_data_mode);
    }

    /* 5-7 | NCSTART, NRSTART, NSSTART | signed int
     * position of first column, first row, and first section (voxel grid units)
     *
     * The position of the first voxel is defined in grid units by NCSTART, NRSTART, and NSSTART.
     * The center of the voxel with grid position (0,0,0) corresponds to the Cartesian coordinate origin.*/
    if (bRead)
    {
        meta_.crs_start = read_int32_ivec_();
    }
    else
    {
        write_int32_ivec_(meta_.crs_start);
    }
    /* 8-10 | NX, NY, NZ | signed int >0 |
     * intervals per unit cell repeat along X,Y Z
     * intervals per map length along X,Y,Z;
     * emdb convention: same as NC, NR, NS
     *
     * Lengths in Aangstroms for a single voxel are as follows:
     * Xvoxel = X_LENGTH/NX Yvoxel = Y_LENGTH/NY Zvoxel = Z_LENGTH/NZ */

    if (bRead)
    {
        meta_.extend = read_int32_ivec_();
        grid_data.set_extend(meta_.extend);
    }
    else
    {
        write_int32_ivec_(grid_data.getExtend());
    }

    /* 11-13 | X_LENGTH, Y_LENGTH, Z_LENGTH | floating pt >0
     * Unit Cell repeats along X, Y, Z In Aangstrom
     * emdb Map lengths along X,Y,Z in Aangstrom
     *
     * Lengths in Aangstroms for a single voxel are as follows:
     * Xvoxel = X_LENGTH/NX Yvoxel = Y_LENGTH/NY Zvoxel = Z_LENGTH/NZ */

    /* 14-16 | ALPHA,BETA,GAMMA | floating pt >0, <180
     * Unit Cell angles (degrees)
     * emdb convention: 90, 90, 90
     *
     * By convention, cell angles (ALPHA, BETA, GAMMA)
     * are 90 degrees for single particle or tomogram EM maps;
     * they follow IUCr space group conventions for crystals.*/
    if (bRead)
    {
        RVec cell_length = read_float32_rvec_();
        svmul(A2NM, cell_length, cell_length);
        RVec cell_angle  = read_float32_rvec_();
        // By convention, unset cell angles (all 0) are interpreted as 90 deg.
        if (cell_angle[XX]*cell_angle[YY]*cell_angle[ZZ] < 1e-5)
        {
            cell_angle = {90, 90, 90};
        }
        grid_data.set_cell(cell_length, cell_angle);
    }
    else
    {
        RVec cell_length_AA;
        svmul(NM2A, grid_data.cell_lengths(), cell_length_AA);
        write_float32_rvec_(cell_length_AA);
        write_float32_rvec_(grid_data.cell_angles());
    }

    /* 17-19 | MAPC, MAPR, MAPS | signed int | 1 (=X) 2 (=Y) 3 (=Z)
     * relationship of X,Y,Z axes to columns, rows, sections
     * emdb convention: 1, 2, 3 */
    if (bRead)
    {
        IVec crs_to_xyz = read_int32_ivec_();

        crs_to_xyz[XX] -= 1;
        crs_to_xyz[YY] -= 1;
        crs_to_xyz[ZZ] -= 1;
        if (colummn_row_section_order_valid_(crs_to_xyz))
        {
            set_crs_to_xyz(crs_to_xyz);
        }
        else
        {
            crs_to_xyz[XX] = 0;
            crs_to_xyz[YY] = 1;
            crs_to_xyz[ZZ] = 2;
            set_crs_to_xyz(crs_to_xyz);
        }
        grid_data.set_translation(grid_data.gridpoint_coordinate(
                                          IVec {
                                              meta_.crs_start[meta_.xyz_to_crs[XX]],
                                              meta_.crs_start[meta_.xyz_to_crs[YY]],
                                              meta_.crs_start[meta_.xyz_to_crs[ZZ]]
                                          }));
    }
    else
    {
        write_int32_ivec_({meta_.crs_to_xyz[XX]+1, meta_.crs_to_xyz[YY]+1, meta_.crs_to_xyz[ZZ]+1});
    }


    /* 20-22 | AMIN, AMAX, AMEAN | floating pt
     * Minimum, maximum, average density */
    if (bRead)
    {
        meta_.min_value  = read_float32_();
        meta_.max_value  = read_float32_();
        meta_.mean_value = read_float32_();
    }
    else
    {
        write_float32_(meta_.min_value);
        write_float32_(meta_.max_value);
        write_float32_(meta_.mean_value);
    }
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
    if (bRead)
    {
        meta_.space_group = read_int32_();
    }
    else
    {
        write_int32_(meta_.space_group);
    }

    /* 24 | NSYMBT | signed int | 80n
     * # of bytes in symmetry table (multiple of 80)
     * emdb convention 0 */
    if (bRead)
    {
        set_num_bytes_extened_header(read_int32_());
    }
    else
    {
        write_int32_(meta_.num_bytes_extened_header);
    }

    if (is_crystallographic())
    {
        /* 25 | LSKFLG | signed int | 0,1
         * flag for skew matrix
         * emdb convention 0 */
        if (bRead)
        {
            set_has_skew_matrix(read_int32_());
        }
        else
        {
            write_int32_(has_skew_matrix());
        }

        if (has_skew_matrix())
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
            if (bRead)
            {
                for (auto && i : meta_.skew_matrix)
                {
                    i = read_float32_();
                }
                meta_.skew_translation = read_float32_rvec_();
            }
            else
            {
                for (auto i : meta_.skew_matrix)
                {
                    write_float32_(i);
                }
                write_float32_rvec_(meta_.skew_translation);
            }
        }
    }
    else
    {
        /* 25-37 not used in EMDB */
        if (bRead)
        {
            for (auto && i : meta_.extraskew)
            {
                i = read_int32_();
            }
        }
        else
        {
            for (auto i : meta_.extraskew)
            {
                write_float32_(i);
            }
        }
    }

    /* 38-52 | EXTRA | 32 bit binary
     * user-defined metadata
     *
     * SKWMAT, SKWTRN, and EXTRA fields are not currently used by EMDB.
     * EMDB might use fields 50,51 and 52 for setting the coordinate system origin */

    if (bRead)
    {
        for (auto && i : meta_.extra)
        {
            i = read_float32_();
        }
    }
    else
    {
        for (auto i : meta_.extra)
        {
            write_float32_(i);
        }
    }

    if (!is_crystallographic() && bRead)
    {
        /* If this the map origin is shifted, because the grid indexing starts at other values than zero,
         * values read here are ignored.
         * Convention is not clear at this point, whether the translation here should be treated as extra shift,
         * in observed cases this was not the case.
         *
         * Silently ignore if the map translation due to grid-index start shift
         * does not match the shift reported here.
         */
        if (meta_.crs_start[XX] == 0 && meta_.crs_start[YY] == 0 && meta_.crs_start[ZZ] == 0)
        {
            grid_data.set_translation({(float) A2NM * meta_.extra[size_extrarecord-3], (float) A2NM * meta_.extra[size_extrarecord-2], (float) A2NM * meta_.extra[size_extrarecord-1]});
        }
    }


    /* 53 | MAP | ASCII char
     * "MAP "
     * MRC/CCP4 MAP format identifier */
    if (bRead)
    {
        read_format_identifier();
    }
    else
    {
        for (auto i : meta_.format_identifier)
        {
            fwrite(&i, 1, 1, file_);
        }
    }
    /* 54 | MACHST | 32 bit
     * binary machine stamp
     *
     * MACHST is (written/read as 4 hex byte sequence)
     * 0x44,0x41,0x00,0x00  for little endian machines
     * 0x11,0x11,0x00,0x00  for big endian machines
     */
    if (bRead)
    {
        meta_.machine_stamp = read_int32_();
    }
    else
    {
        write_int32_(meta_.machine_stamp);
    }
    /* 55 | RMS | floating pt
     * Density root-mean-square deviation */
    if (bRead)
    {
        meta_.rms_value = read_float32_();
    }
    else
    {
        write_float32_(meta_.rms_value);
    }
    /* 56 | NLABL | signed int | 0-10
     * # of labels
     *
     * Following the 2010 remediation, maps distributed by EMDB
     * now have a single label of form “::::EMDataBank.org::::EMD-1234::::”.  */
    if (bRead)
    {
        meta_.num_labels = read_int32_();
    }
    else
    {
        write_int32_(meta_.labels.size());
    }

    /* 57-256 | LABEL_N | ASCII char
     * Up to 10 user-defined labels */
    if (bRead)
    {
        set_labels();
    }
    else
    {
        write_labels();
    }

    /* 257-257+NSYMBT | anything
     */
    if (bRead)
    {
        set_extended_header();
    }
    else
    {
        write_extended_header();
    }

};

void MrcFile::Impl::swap_int32_(gmx_int32_t &result)
{

    gmx_int32_t src = result;
    result = 0;

    result |= (src & 0xFF000000) >> 24;
    result |= (src & 0x00FF0000) >> 8;
    result |= (src & 0x0000FF00) << 8;
    result |= (src & 0x000000FF) << 24;
}


int MrcFile::Impl::read_int32_()
{
    gmx_int32_t result;
    if (fread(&result, sizeof(result), 1, file_) != 1)
    {
        GMX_THROW(gmx::FileIOError("Cannot read integer from volume data at " + std::to_string(ftell(file_)) + "."));
    }
    ;

    /* byte swap in case of different endianness */
    if (meta_.swap_bytes)
    {
        swap_int32_(result);
    }
    return result;
}

char MrcFile::Impl::read_char_()
{
    char result;
    if (fread(&result, 1, 1, file_) != 1)
    {
        GMX_THROW(gmx::FileIOError("Cannot read char from volume data at " + std::to_string(ftell(file_)) + "."));
    }
    return result;
}

void MrcFile::Impl::swap_float32_(float &result)
{
    swap_int32_(reinterpret_cast<gmx_int32_t &>(result));
}


float MrcFile::Impl::read_float32_()
{
    float result;
    if (fread(&result, sizeof(result), 1, file_) != 1)
    {
        GMX_THROW(gmx::FileIOError("Cannot read float from volume data at " + std::to_string(ftell(file_)) + "."));
    }
    ;

    if (meta_.swap_bytes)
    {
        swap_float32_(result);
    }
    return result;
}

void MrcFile::Impl::write_int32_(int data)
{
    if (meta_.swap_bytes)
    {
        swap_int32_(data);
    }

    fwrite(&data, sizeof(data), 1, file_);

}


void MrcFile::Impl::write_float32_(real data)
{
    float to_write = static_cast<float> (data);
    if (meta_.swap_bytes)
    {
        swap_float32_(to_write);
    }

    fwrite(&to_write, sizeof(to_write), 1, file_);
}

void MrcFile::Impl::do_mrc_data_(Field<real> &grid_data, bool bRead)
{
    if (bRead)
    {
        if (file_size_ < header_bytes + meta_.num_bytes_extened_header + long(grid_data.getNumLatticePoints()*sizeof(float)) )
        {
            GMX_THROW(gmx::FileIOError("Density format file is smaller than indicated in its header."));
        }
        else
        {
            if (file_size_ > header_bytes + meta_.num_bytes_extened_header +  long(grid_data.getNumLatticePoints()*sizeof(float)) )
            {
                fprintf(stderr, "WARNING : Density format file size is %ld, however header (%d) + symbol table (%d) + data  (%ld) is larger than indicated in its header. Reading anyway.. \n", file_size_, header_bytes, meta_.num_bytes_extened_header, grid_data.getNumLatticePoints()*sizeof(float));
            }
        }
    }

    IVec num_crs = to_crs_order(grid_data.getExtend());
    if (bRead)
    {
        auto gridDataAccess = grid_data.access();
        for (int section = 0; section  < num_crs[ZZ]; section++)
        {
            for (int row  = 0; row  < num_crs[YY]; row++)
            {
                for (int column  = 0; column  < num_crs[XX]; column++)
                {
                    gridDataAccess.at(to_xyz_order({column, row, section})) = read_float32_();
                }
            }
        }

    }
    else
    {
        auto gridDataAccess = grid_data.access();
        for (int section = 0; section  < num_crs[ZZ]; section++)
        {
            for (int row  = 0; row  < num_crs[YY]; row++)
            {
                for (int column  = 0; column  < num_crs[XX]; column++)
                {
                    write_float32_(gridDataAccess.at(to_xyz_order({column, row, section})));
                }
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
    number_columns = read_int32_();
    if (number_columns <= 0 || number_columns >= 65536)
    {
        meta_.swap_bytes = true;
    }

    // rewind the file
    fsetpos(file_, &current);

}

void MrcFile::Impl::do_mrc(const Field<real> &grid_data, bool bRead)
{
    //TODO: avoid const cast
    do_mrc_header_(*(const_cast<Field<real>*>(&grid_data)), bRead);
    do_mrc_data_(*(const_cast<Field<real>*>(&grid_data)), bRead);
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
    meta_.crs_start                = {0, 0, 0};
    meta_.crs_to_xyz               = {0, 1, 2};
    meta_.xyz_to_crs               = {0, 1, 2};
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

void MrcFile::Impl::set_meta(const Field<real> &grid_data)
{
    set_grid_stats(grid_data);
    IVec index_of_origin = grid_data.coordinate_to_gridindex_floor_ivec(RVec {1e-6, 1e-6, 1e-6});
    meta_.crs_start = {-index_of_origin[XX], -index_of_origin[YY], -index_of_origin[ZZ]};
}

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
    impl_->do_mrc(grid_data, bRead);
    impl_->close_file();
}

void MrcFile::write(std::string filename, const Field<real> &grid_data)
{
    bool bRead = false;
    impl_->set_meta(grid_data);
    impl_->open_file(filename, bRead);
    impl_->do_mrc(grid_data, bRead);
    impl_->close_file();
}


void MrcFile::read_meta(std::string filename, MrcMetaData &meta)
{
    Field<real> griddata;
    bool        bRead = true;
    impl_->open_file(filename, bRead);
    impl_->do_mrc_header_(griddata, bRead);
    impl_->close_file();
    meta = impl_->meta_;
}

void MrcFile::read_with_meta(std::string filename, Field<real> &grid_data, MrcMetaData &meta)
{
    read(filename, grid_data);
    meta = impl_->meta_;
}

void MrcFile::read(std::string filename, Field<real> &grid_data)
{
    bool bRead = true;
    impl_->open_file(filename, bRead);
    impl_->do_mrc(grid_data, bRead);
    impl_->close_file();
}

Df3File::SuccessfulDf3Write
Df3File::write(std::string filename, const gmx::Field<real> &grid_data)
{
    auto    file_ = gmx_fio_fopen(filename.c_str(), "w");
    int16_t xExtendShort {
        int16_t(grid_data.getExtend()[XX])
    };
    int16_t yExtendShort {
        int16_t(grid_data.getExtend()[YY])
    };
    int16_t zExtendShort {
        int16_t(grid_data.getExtend()[ZZ])
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
    for (auto voxel : grid_data.access().data())
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
    auto        file         = gmx_fio_fopen((filename_+".pov").c_str(), "w");
    std::string povRayString = "#declare DD = <" + std::to_string(NM2A * gridData_.cell_lengths()[XX]) + "," + std::to_string(NM2A *gridData_.cell_lengths()[YY]) + "," + std::to_string(NM2A *gridData_.cell_lengths()[ZZ]) + ">;\n";
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
        +"\t\ttranslate <" + std::to_string(NM2A *gridData_.translation()[XX]) + ","+std::to_string(NM2A *gridData_.translation()[YY])+ "," + std::to_string(NM2A *gridData_.translation()[ZZ]) + ">\n"
        +"\t}\n";
    fprintf(file, "%s", povRayString.c_str());
    gmx_fio_fclose(file);

}

} // namespace gmx
