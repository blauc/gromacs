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
#include "mpi-helper.h"
#include "gromacs/utility/real.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/exceptions.h"

namespace gmx
{

MpiHelper::MpiHelper(MPI_Comm mpi_comm_mygroup, int masterrank, bool bMaster) :
    buf_write_(true), mpi_comm_mygroup_(mpi_comm_mygroup), masterrank_(masterrank),  bMaster_(bMaster){}

void MpiHelper::buffer(real value)
{
    if (buf_write_)
    {
        inbuf_.push_back(value);
    }
    else
    {
        if (bMaster_)
        {
            value = outbuf_.back();
        }
        outbuf_.pop_back();
    }
}

void MpiHelper::buffer( real * vector, int size)
{
    for (int i = 0; i < size; i++)
    {
        buffer(vector[i]);
    }
};

void MpiHelper::buffer(real matrix[DIM][DIM])
{
    for (int i = 0; i < DIM; i++)
    {
        buffer(matrix[i], DIM);
    }
};

void MpiHelper::sum_reduce()
{
    if (inbuf_.empty())
    {
        return;
    }
    outbuf_.resize(inbuf_.size());
#ifdef GMX_MPI
    MPI_Reduce(&inbuf_[0], &outbuf_[0], inbuf_.size(), GMX_MPI_REAL, MPI_SUM, masterrank_, mpi_comm_mygroup_);
#endif
    inbuf_.clear();
    buf_write_ = false;
};

void MpiHelper::finish()
{
    if (outbuf_.size() != 0)
    {
        GMX_THROW(APIError("MpiHelper : Buffer not read completely before switching to writing."));
    }
    else
    {
        buf_write_ = true;
    }
}

bool MpiHelper::isMaster()
{
    return bMaster_;
}

} // namespace gmx
