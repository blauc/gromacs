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

#include "externalpotentialIO.h"
#include <stdio.h>
#include <string>

#include "gromacs/utility/futil.h"

namespace gmx
{

namespace externalpotential
{
FILE*  ExternalPotentialIO::input_file()
{
    return input_file_;
}

FILE* ExternalPotentialIO::output_file()
{
    return output_file_;
}

FILE* ExternalPotentialIO::log_file()
{
    return log_file_;
}

unsigned long ExternalPotentialIO::Flags()
{
    return Flags_;
}

const gmx_output_env_t * ExternalPotentialIO::oenv()
{
    return oenv_;
}

bool ExternalPotentialIO::isVerbose()
{
    return bVerbose_;
}

void ExternalPotentialIO::set_input_file(std::string basepath, std::string filename)
{
    input_file_ = open_(basepath, filename, "r");
}

void ExternalPotentialIO::set_output_file(std::string basepath, std::string filename)
{
    output_file_ = open_(basepath, filename, "w");
}

void ExternalPotentialIO::set_log_file(std::string basepath, std::string filename)
{
    log_file_ = open_(basepath, filename, "w");
}

FILE * ExternalPotentialIO::open_(std::string basename, std::string filename, const char * mode)
{
    if (gmx_fexist((basename + "/" + filename).c_str()))
    {
        return gmx_ffopen((basename + "/" + filename).c_str(), mode);
    }
    else
    {
        return nullptr;
    }
}

}       // namespace externalpotential

}       // namespace gmx
