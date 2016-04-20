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
GaussTransform::set_n_sigma(int n_sigma)
{
    n_sigma_ = n_sigma;
};

std::unique_ptr<GridReal> grid_;
real                      sigma_;
int                       n_sigma_;


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
        nu_ = grid_->avg_spacing()/sigma_;
    }
    else
    {
        GMX_THROW(gmx::InconsistentInputError("Grid needs to be evently spaced to use the current implementation of fast gaussian gridding."));
    }
    grid_->resize();
    m_spread_ = int(ceil(n_sigma_/nu_));
    for (int i = -m_spread_; i <= m_spread_; ++i)
    {
        E3_.push_back(exp( -i * i * nu_ * nu_ / 2. ));
    }
};

void
FastGaussianGridding::gauss_1d_()
{
    real E2_power_l;
    for (size_t i = XX; i <= ZZ; i++)
    {

        spread_1d_[i].resize(2*m_spread_+1);

        E1_        = exp(-dx_[i]*dx_[i]/2.0);
        E2_        = exp(dx_[i]*nu_);
        E2_power_l = E2_;

        spread_1d_[i][m_spread_] = E1_;
        for (int l = 1; l < m_spread_; l++)
        {

            spread_1d_[i][m_spread_-l] = (E1_ / E2_power_l) * E3_[m_spread_-l];
            spread_1d_[i][m_spread_+l] = (E1_ * E2_power_l) * E3_[m_spread_+l];

            E2_power_l *= E2_;
        }
    }

};

void
FastGaussianGridding::tensor_product_1d_(real weight)
{
    spread_2d_.resize(2*m_spread_+1);
    real spread_z;
    for (int l_z = 0; l_z < 2 * m_spread_ + 1; ++l_z)
    {
        spread_2d_[l_z].resize(2*m_spread_+1);
        spread_z = weight*spread_1d_[ZZ][l_z];
        for (int l_y = 0; l_y < 2 * m_spread_ + 1; ++l_y)
        {
            spread_2d_[l_z][l_y] = spread_z*spread_1d_[YY][l_y];
        }
    }

    ivec l_min;
    ivec l_max;
    for (size_t i = 0; i <= ZZ; ++i)
    {
        l_min[i] = grid_index_[i]-m_spread_ < 0 ? 0 : grid_index_[i]-m_spread_;
        l_max[i] = grid_index_[i]+m_spread_ >= grid_->extend()[i] ? grid_->extend()[i]-1 : grid_index_[i]+m_spread_;
    }
    real   spread_zy;
    int    l_y, l_z;
    std::vector<real>::iterator voxel;
    int    l_x_min = l_min[XX]-grid_index_[XX]+m_spread_;
    int    n_l_x   = l_max[XX]-l_min[XX];
    real * spread_1d_XX;

    for (int l_grid_z = l_min[ZZ]; l_grid_z <= l_max[ZZ]; ++l_grid_z)
    {
        l_z = l_grid_z - grid_index_[ZZ]+m_spread_;
        for (int l_grid_y = l_min[YY]; l_grid_y <= l_max[YY]; ++l_grid_y)
        {
            l_y          = l_grid_y - grid_index_[YY]+m_spread_;
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

void
FastGaussianGridding::set_grid_index_and_dx_(real * x)
{
    grid_index_ = grid_->coordinate_to_gridindex_floor_ivec(x);
    RVec nearest_grid_point_coordinate = grid_->gridpoint_coordinate(grid_index_);
    for (size_t i = XX; i <= ZZ; ++i)
    {
        dx_[i]         = (x[i]-nearest_grid_point_coordinate[i])/sigma_;
    }
    ;
};

void
FastGaussianGridding::transform(real *x, real weight)
{
    set_grid_index_and_dx_(x);
    gauss_1d_();
    tensor_product_1d_(weight);
}

std::unique_ptr<GridReal> &&
FastGaussianGridding::finish_and_return_grid()
{
    return std::move(grid_);
}

}    // namespace volumedata

}    // namespace volumedata
