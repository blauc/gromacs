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
#include "gmxpre.h"
#include <algorithm>
#include "template.h"
#include <cstdio>

#include "gromacs/math/vec.h"
#include "gromacs/externalpotential/group.h"

namespace gmx
{

std::unique_ptr<ExternalPotential> Template::create(struct ext_pot_ir *ep_ir, t_commrec * cr, t_inputrec * ir, const gmx_mtop_t* mtop, const rvec x[], matrix box, FILE *input_file, FILE *output_file, FILE *fplog, bool bVerbose, const gmx_output_env_t *oenv, unsigned long Flags, int number_groups)
{
    return std::unique_ptr<ExternalPotential> (new Template(ep_ir, cr, ir, mtop, x, box, input_file, output_file, fplog, bVerbose, oenv, Flags, number_groups));
}

void Template::do_potential( t_commrec *cr, const matrix box, const rvec x[], const gmx_int64_t step)
{
    GroupCoordinates   r1_local = x_local(x, 0);
    std::vector<RVec>  r1       = x_assembled(step, cr, x, box, 0);
    RVec               com1;
    // = std::accumulate( r1.begin(), r1.end(), RVec(0, 0, 0), [](const RVec &a, const RVec &b) {return RVec(a[XX]+b[XX], a[YY]+b[YY], a[ZZ]+b[ZZ]); } );
    GroupCoordinates   r2_local = x_local(x, 1);
    // std::vector<RVec> r2       = x_assembled(step, cr, x, box, 1);
    RVec               com2;
    //   = std::accumulate( r2.begin(), r2.end(), RVec(0, 0, 0), [](const RVec &a, const RVec &b) {return RVec(a[XX]+b[XX], a[YY]+b[YY], a[ZZ]+b[ZZ]); } );
    real               potential = 0;

    com1[XX] = box[XX][XX]*0.25;
    com1[YY] = box[YY][YY]*0.25;
    com1[ZZ] = box[ZZ][ZZ]*0.25;

    com2[XX] = 3*com1[XX];
    com2[YY] = 3*com1[YY];
    com2[ZZ] = 3*com1[ZZ];

    std::vector<RVec> &f1_local = f_local(0);

    for (size_t i = 0; i < f1_local.size(); i++)
    {
        rvec_sub(r1_local[i], com1, f1_local[i]);
        svmul(-k_, f1_local[i], f1_local[i]);

        potential      += k_*distance2(r1_local[i], com1)/2.0;
    }

    std::vector<RVec> &f2_local = f_local(1);
    for (size_t i = 0; i < f2_local.size(); i++)
    {
        rvec_sub(r2_local[i], com2, f2_local[i]);
        svmul(-k_, f2_local[i], f2_local[i]);
        potential      += k_*distance2(r2_local[i], com2)/2.0;
    }

    set_local_potential(potential);

    (void) cr;
    (void) box;
    (void) x;
    (void) step;
};

Template::Template(struct ext_pot_ir *ep_ir, t_commrec * cr, t_inputrec * ir, const gmx_mtop_t* mtop, const rvec x[], matrix box, FILE *input_file_p, FILE *output_file_p, FILE *fplog, bool bVerbose, const gmx_output_env_t *oenv, unsigned long Flags, int number_groups) :
    ExternalPotential(ep_ir, cr, ir, mtop, x, box, input_file_p, output_file_p, fplog, bVerbose, oenv, Flags, number_groups)
{
    if (MASTER(cr) || !PAR(cr))
    {
        char * line = nullptr;
        size_t len  = 0;
        if (input_file() != nullptr)
        {
            getline(&line, &len, input_file());
            k_ = strtof(line, nullptr);
        }
        else
        {
            fprintf(stderr, "\n\n\nNo inputfile found, setting k to default value..\n\n\n");
            k_ = 0;
        }
    }

    if (PAR(cr))
    {
        MPI_Bcast(&k_, 1, GMX_MPI_REAL, MASTERRANK(cr), cr->mpi_comm_mygroup);
    }
}

std::string TemplateInfo::name                        = "Template";
std::string TemplateInfo::shortDescription            = "basic external potential example";
const int   TemplateInfo::numberIndexGroups           = 2;
externalpotential::ModuleCreator TemplateInfo::create = Template::create;


} // namespace gmx
