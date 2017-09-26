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
#ifndef GMX_MPIHELPER_H
#define GMX_MPIHELPER_H


#include <vector>
#include "gromacs/utility/real.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/gmxmpi.h"

namespace gmx
{
class MpiHelper
{
    public:
        explicit MpiHelper(MPI_Comm mpi_comm_mygroup, int masterrank, bool bMaster);
        void cr(t_commrec *cr);
        void sum_reduce();
        void sum_allReduce();
        bool isMaster();
        void broadcast(void * ptr, size_t size);
        real from_reals_buffer();
        void from_reals_buffer(matrix result);
        void from_reals_buffer(real * vector, int size);
        void to_reals_buffer(real value);
        void to_reals_buffer(const real * vector, int size);
        void to_reals_buffer(real matrix[3][3]);
        real max(real value);
        void finish();
        void sum_reduce_rvec(rvec vector);
    private:
        std::vector<real> inbuf_;
        std::vector<real> outbuf_;
        bool              buf_write_;
        MPI_Comm          mpi_comm_mygroup_;
        int               masterrank_;
        bool              bMaster_;
};

} // namespace gmx
#endif
