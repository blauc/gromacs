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
/*!  \file
 * \brief
 * Fourier transform wrapper implementation.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \inpublicapi
 */

#include "fouriertransform.h"
#include "gromacs/fft/parallel_3dfft.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/utility/real.h"
#include "densitypadding.h"
#include "gridreal.h"
#include "gromacs/fileio/volumedataio.h"
#include <cassert>
#include <algorithm>

namespace gmx
{
namespace volumedata
{

FourierTransformBackAndForth::FourierTransformBackAndForth(
        const Field<real> &input)
    : input_ {input}
{};

std::unique_ptr < Field < real>> FourierTransformBackAndForth::transform() {
    gmx_parallel_3dfft_t fft = nullptr;

    MPI_Comm             mpiCommunicatorForFFT[] = {MPI_COMM_NULL, MPI_COMM_NULL};
    auto                 paddedInput             = *volumedata::DensityPadding(input_).pad({3, 3, 3});
    real                *rdata                   = (real*)(malloc(sizeof(real)*paddedInput.num_gridpoints()));
    t_complex           *cdata                   = (t_complex*)malloc(sizeof(t_complex)*paddedInput.num_gridpoints());
    gmx_parallel_3dfft_init(&fft, paddedInput.extend(), &rdata, &cdata,
                            mpiCommunicatorForFFT, true, 1);
    fprintf(stderr, "\nsize t_complex: %lu size real: %lu \n", sizeof(t_complex), sizeof(real));
    IVec ndata;
    IVec nsize;
    gmx_parallel_3dfft_real_limits(fft, ndata, IVec(), nsize);
    fprintf(stderr, "SIZE  %d %d \n", ndata[XX], nsize[XX]);
    fprintf(stderr, "SIZE  %d %d \n", ndata[YY], nsize[YY]);
    fprintf(stderr, "SIZE  %d %d \n", ndata[ZZ], nsize[ZZ]);

    assert(ndata[ZZ] * ndata[YY] * ndata[XX] ==
           int(std::distance(paddedInput.access().begin(), paddedInput.access().end())));
    std::copy(paddedInput.access().begin(), paddedInput.access().end(), rdata);

    gmx_parallel_3dfft_execute(fft, GMX_FFT_REAL_TO_COMPLEX, 1, nullptr);

    std::unique_ptr < Field < t_complex>> output(new Field<t_complex>());
    output->copy_grid(paddedInput);
    output->convertToReciprocalSpace();

    std::copy(cdata, cdata + output->num_gridpoints()/2 + output->extend()[XX] *output->extend()[YY] +1, output->access().begin());

    auto orthoscale                = 1 / sqrt(input_.num_gridpoints());
    auto normalizeComplexTransform = [orthoscale](t_complex &value) {
            value.re *= orthoscale;
            value.im *= orthoscale;
        };
    std::for_each(output->access().begin(), output->access().end(), normalizeComplexTransform);

    std::copy(output->access().begin(), output->access().begin()+output->num_gridpoints()/2+ output->extend()[XX] *output->extend()[YY] + 1, cdata);

    volumedata::GridReal plotComplex;
    plotComplex.copy_grid(*output);
    std::transform(output->access().begin(), output->access().begin()+output->num_gridpoints()/2+ output->extend()[XX] *output->extend()[YY] + 1, plotComplex.access().begin(), [](const t_complex &o){return o.re; });
    MrcFile().write("RealValues.ccp4", plotComplex);

    std::transform(output->access().begin(), output->access().begin()+output->num_gridpoints()/2+ output->extend()[XX] *output->extend()[YY] + 1, plotComplex.access().begin(), [](const t_complex &o){return o.im; });
    MrcFile().write("ImValues.ccp4", plotComplex);

    std::transform(output->access().begin(), output->access().begin()+output->num_gridpoints()/2+ output->extend()[XX] *output->extend()[YY] + 1, plotComplex.access().begin(), [](const t_complex &o){return o.re; });

    gmx_parallel_3dfft_execute(fft, GMX_FFT_COMPLEX_TO_REAL, 1, nullptr);

    std::unique_ptr < Field < real>> result(new Field<real>());
    result->copy_grid(paddedInput);
    for (auto &x  : result->access())
    {
        x = 0;
    }
    std::copy(rdata, rdata + result->num_gridpoints(),
              result->access().begin());

    orthoscale = 1 / sqrt(input_.num_gridpoints());
    auto normalizeTransform = [orthoscale](real &value) {
            value *= orthoscale;
        };
    std::for_each(std::begin(result->access().data()),
                  std::end(result->access().data()), normalizeTransform);

    gmx_parallel_3dfft_destroy(fft);

    return result;
};

FourierTransformRealToComplex3D::FourierTransformRealToComplex3D(
        const Field<real> &input)
    : input_ {input}
{};

std::unique_ptr < Field < t_complex>> FourierTransformRealToComplex3D::transform() {
    gmx_parallel_3dfft_t fft = nullptr;
    real                *rdata;
    t_complex           *cdata;
    MPI_Comm             mpiCommunicatorForFFT[] = {MPI_COMM_NULL, MPI_COMM_NULL};

    gmx_parallel_3dfft_init(&fft, input_.extend(), &rdata, &cdata,
                            mpiCommunicatorForFFT, false,
                            std::max(1, gmx_omp_nthreads_get(emntDefault)));

    IVec ndata;
    gmx_parallel_3dfft_real_limits(fft, ndata, IVec(), IVec());
    assert(ndata[ZZ] * ndata[YY] * ndata[XX] ==
           int(std::distance(input_.access().begin(), input_.access().end())));
    std::copy(input_.access().begin(), input_.access().end(), rdata);

    gmx_parallel_3dfft_execute(fft, GMX_FFT_REAL_TO_COMPLEX, 0, nullptr);

    IVec size;
    gmx_parallel_3dfft_complex_limits(fft, IVec(), ndata, IVec(), size);


    std::unique_ptr < Field < t_complex>> output(new Field<t_complex>());
    output->set_extend(size);
    output->convertToReciprocalSpace();

    std::copy(cdata, cdata + output->num_gridpoints()+1, output->access().begin());
    auto orthoscale = 1 / sqrt(input_.num_gridpoints());

    auto normalizeComplexTransform = [orthoscale](t_complex &value) {
            value.re *= orthoscale;
            value.im *= orthoscale;
        };

    std::for_each(output->access().begin(), output->access().end(),
                  normalizeComplexTransform);

    gmx_parallel_3dfft_destroy(fft);

    return output;
};

FourierTransformComplexToReal3D::FourierTransformComplexToReal3D(
        const Field<t_complex> &input)
    : input_ {input}
{};

IVec FourierTransformComplexToReal3D::firstNonHermitian(real tolerance)
{
    auto extend = input_.extend();
    for (int ix = 0; ix < extend[XX] / 2; ++ix)
    {
        for (int iy = 0; iy < extend[YY]; ++iy)
        {

            for (int iz = 0; iz < extend[ZZ]; ++iz)
            {
                if (ix != 0 || iy != 0 || iz != 0)
                {
                    auto a = input_.access().at({ix, iy, iz});
                    auto b = input_.access().at({-ix, -iy, -iz});
                    if (abs(a.re - b.re) < tolerance)
                    {
                        return {
                                   ix, iy, iz
                        };
                    }
                    ;
                    if (abs(a.im + b.im) < tolerance)
                    {
                        return {
                                   ix, iy, iz
                        };
                    }
                }
            }
        }
    }
    return {
               extend[XX], extend[YY], extend[ZZ]
    };
}

bool FourierTransformComplexToReal3D::isHermitian(real tolerance)
{
    return firstNonHermitian(tolerance) ==
           IVec({input_.extend()[XX], input_.extend()[YY], input_.extend()[ZZ]});
}

std::unique_ptr < Field < real>> FourierTransformComplexToReal3D::transform() {
    gmx_parallel_3dfft_t fft = nullptr;
    real                *rdata;
    t_complex           *cdata;
    MPI_Comm             mpiCommunicatorForFFT[] = {MPI_COMM_NULL, MPI_COMM_NULL};
    auto                 realSize                = input_.extend();
    realSize[ZZ] = (realSize[ZZ]-1)*2;
    gmx_parallel_3dfft_init(&fft, realSize, &rdata, &cdata,
                            mpiCommunicatorForFFT, false,
                            std::max(1, gmx_omp_nthreads_get(emntDefault)));

    // IVec complex_order, ndata, offset, size;
    IVec size;
    gmx_parallel_3dfft_complex_limits(fft, IVec(), IVec(), IVec(), size);
    assert(size[ZZ] * size[YY] * size[XX] ==
           int(std::distance(input_.access().begin(), input_.access().end())));
    std::copy(input_.access().begin(), input_.access().end(), cdata);

    gmx_parallel_3dfft_execute(fft, GMX_FFT_COMPLEX_TO_REAL, 0, nullptr);

    std::unique_ptr < Field < real>> output(new Field<real>());
    output->copy_grid(input_);
    output->set_extend(realSize);
    output->resetCell();
    output->convertToReciprocalSpace();

    std::copy(rdata, rdata + output->num_gridpoints(),
              output->access().data().begin());

    auto orthoscale = 1 / sqrt(output->num_gridpoints());

    auto normalizeTransform = [orthoscale](real &value) {
            value *= orthoscale;
        };
    std::for_each(std::begin(output->access().data()),
                  std::end(output->access().data()), normalizeTransform);

    gmx_parallel_3dfft_destroy(fft);

    return output;
};
}
}
