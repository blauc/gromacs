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
 * Tests for file I/O routines
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_fileio
 */
#include "gmxpre.h"

#include <string>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/volumedataio.h"
#include "gromacs/math/volumedata/volumedata.h"
#include "gromacs/utility/path.h"

#include "testutils/testfilemanager.h"

namespace
{

class MapTest : public ::testing::Test
{
    public:
        MapTest()
        {
        }
        gmx::test::TestFileManager      fileManager_;
};

TEST_F(MapTest, CanReadMapFile)
{
    gmx::volumedata::GridReal grid_data;
    gmx::volumedata::MrcFile  map_file;
    std::string               filename = fileManager_.getInputFilePath("EMD-2578.map");

    map_file.read(filename, grid_data);

}

TEST_F(MapTest, IORoundTripMapFile)
{
    gmx::volumedata::GridReal    grid_data;
    gmx::volumedata::MrcFile     map_file;
    gmx::volumedata::MrcMetaData metadata;
    const int                    header_byte_size = 1024;
    int                          fread_size;
    bool                         bOwnGridStats    = false;

    std::string                  filename = fileManager_.getInputFilePath("EMD-2578.map");
    std::string                  write_filename {
        filename + ".test"
    };

    map_file.read_with_meta(filename, grid_data, metadata);
    map_file.write_with_own_meta(write_filename, grid_data, metadata, bOwnGridStats);

    // Files should be of same size and have at least 1024 bytes to contain the header

    t_fileio *read_map    = gmx_fio_open(filename.c_str(), "r");
    t_fileio *written_map = gmx_fio_open(write_filename.c_str(), "r");
    fseek(gmx_fio_getfp(read_map), 0, SEEK_END);
    fseek(gmx_fio_getfp(written_map), 0, SEEK_END);

    gmx_off_t read_size    = gmx_fio_ftell(read_map);
    gmx_off_t written_size = gmx_fio_ftell(written_map);

    ASSERT_EQ(read_size, written_size);

    ASSERT_GE(read_size, header_byte_size);
    ASSERT_GE(written_size, header_byte_size);

    // Byte by byte comparison of header, since gmx MD5 - checksum implementation does not work on files < 1MB

    gmx_fio_rewind(read_map);
    gmx_fio_rewind(written_map);

    unsigned char read_header[header_byte_size];
    unsigned char written_header[header_byte_size];

    fread_size = fread(read_header, 1, header_byte_size, gmx_fio_getfp(read_map));
    ASSERT_EQ(fread_size, header_byte_size);
    fread_size = fread(written_header, 1, header_byte_size, gmx_fio_getfp(written_map));
    ASSERT_EQ(fread_size, header_byte_size);

    ASSERT_THAT(read_header, testing::ContainerEq(written_header));

    std::vector<unsigned char> read_data(read_size-header_byte_size);
    std::vector<unsigned char> written_data(written_size-header_byte_size);

    fread_size = fread(&read_data[0], 1, read_size-header_byte_size, gmx_fio_getfp(read_map));
    ASSERT_EQ(fread_size, read_size-header_byte_size);
    fread_size = fread(&written_data[0], 1, written_size-header_byte_size, gmx_fio_getfp(written_map));
    ASSERT_EQ(fread_size, written_size-header_byte_size);

    ASSERT_THAT(read_data, testing::ContainerEq(written_data));

    gmx_fio_close(read_map);
    gmx_fio_close(written_map);

    remove(write_filename.c_str());
}

} // namespace
