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
 * Defines volume data containers.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \inpublicapi
 */
#include "volumedata.h"
#include "gausstransform.h"
#include "gromacs/utility/exceptions.h"

namespace gmx
{
namespace volumedata
{

void
GaussTransform::set_grid(std::unique_ptr<GridReal> grid)
{
    grid_ = std::move((grid));
};

void
GaussTransform::set_sigma(real sigma)
{
    sigma_ = sigma;
};

void
GaussTransform::set_n_sigma(real n_sigma)
{
    n_sigma_ = n_sigma;
};

void
FastGaussianGridding::set_grid(std::unique_ptr<GridReal> grid)
{
    if (grid->rectangular())
    {
        GaussTransform::set_grid(std::move(grid));
    }
    else
    {
        GMX_THROW(gmx::InconsistentInputError("Grid needs to be rectangular to use the current implementation of fast gaussian gridding."));
    }
    if (grid_->spacing_is_same_xyz())
    {
        nu_       = grid_->avg_spacing()/sigma_;
        m_spread_ = int(ceil(n_sigma_/nu_)); // number of grid cells for spreading
        E3_.resize(2*m_spread_+1);

        for (int i = -m_spread_; i <= m_spread_; ++i)
        {
            E3_[i+m_spread_] = exp( -i * i * nu_ * nu_ / 2. );
        }
    }
    else
    {
        GMX_THROW(gmx::InconsistentInputError("Grid needs to be evently spaced to use the current implementation of fast gaussian gridding."));
    }
    grid_->resize();
};

std::array<std::vector<real>, 3>
FastGaussianGridding::spread_1d_xyz_(real weight, int m_spread, rvec dx, real nu, const std::vector<real> &E3)
{
    std::array<std::vector<real>, 3> result;

    for (size_t i = XX; i <= ZZ; i++)
    {

        result[i].resize(2*m_spread+1);

        real E1         = weight*exp(-dx[i]*dx[i]/2.0); //< exp(-dx_*dx_/2) , following the naming convention of Greengard et al., ;
        real E2         = exp(dx[i]*nu);                //< exp(dx_*nu_) , following the naming convention of Greengard et al., ;
        real E2_power_l = E2;

        result[i][m_spread] = E1;
        for (int l = 1; l < m_spread; l++)
        {
            result[i][m_spread-l] = (E1 / E2_power_l) * E3[m_spread-l];
            result[i][m_spread+l] = (E1 * E2_power_l) * E3[m_spread+l];

            E2_power_l *= E2;
        }
    }
    return result;
};

void
FastGaussianGridding::tensor_product_2d_()
{
    spread_2d_.resize(2*m_spread_+1);
    real spread_z;
    for (int l_z = 0; l_z < 2 * m_spread_ + 1; ++l_z)
    {
        spread_2d_[l_z].resize(2*m_spread_+1);
        spread_z = spread_1d_[ZZ][l_z];
        for (int l_y = 0; l_y < 2 * m_spread_ + 1; ++l_y)
        {
            spread_2d_[l_z][l_y] = spread_z*spread_1d_[YY][l_y];
        }
    }
}

void
FastGaussianGridding::tensor_product_()
{

    ivec l_min;
    ivec l_max;
    for (size_t i = XX; i <= ZZ; ++i)
    {
        l_min[i] = grid_index_[i]-m_spread_ < 0 ? 0 : grid_index_[i]-m_spread_;
        l_max[i] = grid_index_[i]+m_spread_ >= grid_->extend()[i] ? grid_->extend()[i]-1 : grid_index_[i]+m_spread_;
    }
    real   spread_zy;
    std::vector<real>::iterator voxel;
    int    l_x_min = l_min[XX]-grid_index_[XX]+m_spread_;
    int    n_l_x   = l_max[XX]-l_min[XX];
    real * spread_1d_XX;

    for (int l_grid_z = l_min[ZZ]; l_grid_z <= l_max[ZZ]; ++l_grid_z)
    {
        int l_z = l_grid_z - grid_index_[ZZ]+m_spread_;
        for (int l_grid_y = l_min[YY]; l_grid_y <= l_max[YY]; ++l_grid_y)
        {
            int l_y          = l_grid_y - grid_index_[YY]+m_spread_;
            spread_zy    = spread_2d_[l_z][l_y];
            voxel        = grid_->zy_column_begin(l_grid_z, l_grid_y)+l_min[XX];
            spread_1d_XX = &(spread_1d_[XX][l_x_min]);

            for (int l_x = 0; l_x <= n_l_x; ++l_x)
            {
                *voxel += spread_zy * (*spread_1d_XX);
                ++spread_1d_XX;
                ++voxel;
            }
        }
    }
};

void FastGaussianGridding::prepare_2d_grid(const rvec x, const real weight)
{
    grid_index_ = grid_->coordinate_to_gridindex_floor_ivec(x);
    RVec dx; // (x-nearest voxel)/sigma
    for (size_t i = XX; i <= ZZ; ++i)
    {
        dx[i]  = (x[i]-grid_->gridpoint_coordinate(grid_index_)[i])/sigma_;
    }
    spread_1d_ = spread_1d_xyz_(weight, m_spread_, dx, nu_, E3_);
    tensor_product_2d_();
}

void
FastGaussianGridding::transform(const rvec x, real weight)
{
    prepare_2d_grid(x, weight);
    tensor_product_();
}

std::unique_ptr<GridReal> &&
FastGaussianGridding::finish_and_return_grid()
{
    return std::move(grid_);
}

}    // namespace volumedata

}    // namespace gmx
