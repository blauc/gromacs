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

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/futil.h"


namespace gmx
{
namespace volumedata
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

        void do_mrc(volumedata::GridReal * grid_data, bool bRead);
        void do_mrc_data_(volumedata::GridReal * grid_data, bool bRead);
        void do_mrc_header_(volumedata::GridReal * grid_data, bool bRead);

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
        void write_int32_ivec_(IVec i);
        void write_float32_(real data);
        void write_float32_rvec_(RVec i);

        void swap_int32_(gmx_int32_t &result);
        void swap_float32_(float &result);

        FILE     *file_;
        long      file_size_;
        const int size_extrarecord;
        const int number_labels;
        const int label_size;
        const int header_bytes;
        const std::vector<std::string> filetypes;


        MrcMetaData meta_;
};


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

    i_xyz[meta_.crs_to_xyz[XX]] = i_crs[XX] - meta_.crs_start[XX];
    i_xyz[meta_.crs_to_xyz[YY]] = i_crs[YY] - meta_.crs_start[YY];
    i_xyz[meta_.crs_to_xyz[ZZ]] = i_crs[ZZ] - meta_.crs_start[ZZ];

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

void MrcFile::Impl::write_float32_rvec_(RVec i)
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

void MrcFile::Impl::write_int32_ivec_(IVec i)
{
    write_int32_(i[XX]);
    write_int32_(i[YY]);
    write_int32_(i[ZZ]);
}

void MrcFile::Impl::do_mrc_header_(volumedata::GridReal * grid_data, bool bRead)
{

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
        write_int32_ivec_(to_crs_order(grid_data->extend()));
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
        grid_data->set_extend(read_int32_ivec_());
    }
    else
    {
        write_int32_ivec_(grid_data->extend());
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
        RVec cell_angle  = read_float32_rvec_();
        grid_data->set_cell(cell_length, cell_angle);
    }
    else
    {
        write_float32_rvec_(grid_data->cell_lengths());
        write_float32_rvec_(grid_data->cell_angles());
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
        set_crs_to_xyz(crs_to_xyz);
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
        grid_data->set_space_group(read_int32_());
    }
    else
    {
        write_int32_(grid_data->space_group());
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

    if (!is_crystallographic())
    {
        grid_data->set_translation({meta_.extra[size_extrarecord-3], meta_.extra[size_extrarecord-2], meta_.extra[size_extrarecord-1]});
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
        write_int32_(meta_.num_labels);
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

void MrcFile::Impl::do_mrc_data_(volumedata::GridReal * grid_data, bool bRead)
{

    IVec num_crs = to_crs_order(grid_data->extend());
    if (bRead)
    {
        grid_data->resize();
        for (int section = 0; section  < num_crs[ZZ]; section++)
        {
            for (int row  = 0; row  < num_crs[YY]; row++)
            {
                for (int column  = 0; column  < num_crs[XX]; column++)
                {
                    grid_data->at(to_xyz_order({column, row, section})) = read_float32_();
                }
            }
        }

    }
    else
    {
        for (int section = 0; section  < num_crs[ZZ]; section++)
        {
            for (int row  = 0; row  < num_crs[YY]; row++)
            {
                for (int column  = 0; column  < num_crs[XX]; column++)
                {
                    write_float32_(grid_data->at(to_xyz_order({column, row, section})));
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

void MrcFile::Impl::do_mrc(volumedata::GridReal * grid_data, bool bRead)
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

    do_mrc_header_(grid_data, bRead);

    if (bRead)
    {
        if (file_size_ < header_bytes + meta_.num_bytes_extened_header + long(grid_data->num_gridpoints()*sizeof(float)) )
        {
            GMX_THROW(gmx::FileIOError("Density format file is smaller than indicated in its header."));
        }
        else
        {
            if (file_size_ > header_bytes + meta_.num_bytes_extened_header +  long(grid_data->num_gridpoints()*sizeof(float)) )
            {
                fprintf(stderr, "WARNING : Density format file size is %ld, however header (%d) + symbol table (%d) + data  (%ld) is larger than indicated in its header. Reading anyway.. \n", file_size_, header_bytes, meta_.num_bytes_extened_header, grid_data->num_gridpoints()*sizeof(float));
            }
        }
    }

    do_mrc_data_(grid_data, bRead);
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
    meta_.mrc_data_mode            = 2;
    meta_.num_bytes_extened_header = 0;
    meta_.has_skew_matrix          = false;
    meta_.crs_start                = {0, 0, 0};
    meta_.crs_to_xyz               = {0, 1, 2};
    meta_.skew_matrix              = {{0, 0, 0, 0, 0, 0, 0, 0, 0}};
    meta_.skew_translation         = {0, 0, 0};
    meta_.is_crystallographic      = false;
    meta_.extra                    = {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
    meta_.extraskew                = {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, }};
    meta_.format_identifier        = "MAP ";
#if GMX_INTEGER_BIG_ENDIAN
    meta_.machine_stamp            = 0x11110000;
#else
    meta_.machine_stamp            = 0x44410000;
#endif

    meta_.num_labels               = 1;
    meta_.labels                   = {
        {"::::EMDataBank.org::::EMD-xxxx::::Own Data Following EMDB convention:::::::::::"},
        {"                                                                               "},
        {"                                                                               "},
        {"                                                                               "},
        {"                                                                               "},
        {"                                                                               "},
        {"                                                                               "},
        {"                                                                               "},
        {"                                                                               "},
        {"                                                                               "}
    };
    meta_.extended_header = {};
}

void MrcFile::Impl::close_file()
{
    gmx_fio_fclose(file_);
    file_ = nullptr;
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

/*******************************************************************************
 * MrcFile
 */
MrcFile::MrcFile() : impl_(new MrcFile::Impl)
{
}

MrcFile::~MrcFile()
{
}

void MrcFile::write_with_own_meta(std::string filename, GridReal *grid_data, MrcMetaData *meta, bool bOwnGridStats)
{
    bool bRead = false;

    if (meta == nullptr)
    {
        GMX_THROW(gmx::InternalError("Cannot write MrcFile with own metadata, MrcMetaData pointer is null."));
    }
    if (grid_data == nullptr)
    {
        GMX_THROW(gmx::InternalError("Cannot write MrcFile, GridReal pointer to grid_data is null."));
    }

    impl_->meta_ = *meta;

    if (bOwnGridStats)
    {
        impl_->meta_.min_value  = grid_data->min();
        impl_->meta_.max_value  = grid_data->max();
        impl_->meta_.mean_value = grid_data->mean();
        impl_->meta_.rms_value  = grid_data->rms();
    }

    impl_->open_file(filename, bRead);
    impl_->do_mrc(grid_data, bRead);
    impl_->close_file();
}


void MrcFile::write(std::string filename, volumedata::GridReal *grid_data)
{
    if (grid_data == nullptr)
    {
        GMX_THROW(gmx::InternalError("Cannot write MrcFile, GridReal pointer to grid_data is null."));
    }
    bool bRead = false;
    impl_->open_file(filename, bRead);
    impl_->do_mrc(grid_data, bRead);
    impl_->close_file();
}

void MrcFile::read_meta(std::string filename, GridReal *grid_data, MrcMetaData *meta)
{
    if (meta == nullptr)
    {
        GMX_THROW(gmx::InternalError("Cannot read MrcFile with own metadata, MrcMetaData pointer is null."));
    }
    read(filename, grid_data);
    *meta = impl_->meta_;
}

void MrcFile::read(std::string filename, GridReal *grid_data)
{
    if (grid_data == nullptr)
    {
        GMX_THROW(gmx::InternalError("Cannot read MrcFile with own metadata, GridReal pointer to grid_data is null."));
    }
    bool bRead = true;
    impl_->open_file(filename, bRead);
    impl_->do_mrc(grid_data, bRead);
    impl_->close_file();
}

} // namespace volumedata

} // namespace gmx
