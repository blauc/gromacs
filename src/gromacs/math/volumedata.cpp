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
 * Implements methods from volumedata.h
 *
 * \author Christian Blau <cblau@gwdg.de>
 */
#include "gmxpre.h"

#include "volumedata.h"

#include <algorithm>
#include <numeric>
#include <vector>

#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"

namespace gmx
{

namespace volumedata
{

/********************************************************************
 * CrystalSymmetry::Impl
 */

/*! \internal \brief
 * Private implementation class for CrystalSymmetry.
 *
 */
class CrystalSymmetry::Impl
{
    public:

        Impl();
        ~Impl();

        int         space_group_; //!< space group as defined by IUCr conventions (Table 12.3.4.1 Standard space-group symbols”, pages 824-831, International Tables for Crystallography, Volume A, fifth edition)

};

CrystalSymmetry::Impl::Impl() : space_group_(1)
{
}

CrystalSymmetry::Impl::~Impl()
{
}

/********************************************************************
 * CrystalSymmetry
 */

void CrystalSymmetry::set_space_group(int space_group)
{
    impl_->space_group_ = space_group;
}
int CrystalSymmetry::space_group()
{
    return impl_->space_group_;
}

CrystalSymmetry::CrystalSymmetry() : impl_(new CrystalSymmetry::Impl())
{

};

CrystalSymmetry::~CrystalSymmetry()
{

};

/********************************************************************
 * FiniteGrid::Impl
 */

/*! \internal \brief
 * Private implementation class for FiniteGrid.
 *
 */

class FiniteGrid::Impl
{
    public:
        Impl();
        ~Impl();

        IVec             extend_;    //!< The grid-size.
        matrix           cell_;      //!< r' = cell_ . (r-translate_)
        matrix           unit_cell_; //!< cell_/extend_
        RVec             translate_; //!< r' = cell_ . (r-translate_)


        /*! \brief
         * Angle between two vectors in degree.
         */
        real deg_angle(rvec a, rvec b);

        /*! \brief
         * Gram-Schmidt orthonormalisation.
         *
         * The numerical shortcommings of Gram-Schmidt
         * should not matter for three dimensions and cell matrices.
         */
        void QRDecomposition(const matrix A, matrix Q, matrix R);

        /*! \brief
         * Projects vector a onto vector b, stores this projection in c. Conveninence function for Gram-Schmidt method.
         */
        void project(const rvec a, const rvec b, rvec c);

        /*! \brief
         * set unit cell; divide cell matrix by extend in respective direction
         */
        void set_unit_cell();
};

real FiniteGrid::Impl::deg_angle(rvec a, rvec b)
{
    return gmx_angle(a, b)*180.0/M_PI;
}

void FiniteGrid::Impl::project(const rvec a, const rvec b, rvec c)
{
    svmul(iprod(a, b)/iprod(a, a), a, c);
}

void FiniteGrid::Impl::QRDecomposition(const matrix A, matrix Q, matrix R)
{
    rvec u;
    copy_rvec(A[XX], Q[XX]);

    copy_rvec(A[YY], Q[YY]);
    project(Q[YY], Q[XX], u);
    rvec_dec(Q[YY], u);

    copy_rvec(A[ZZ], Q[ZZ]);
    project(Q[ZZ], Q[XX], u);
    rvec_dec(Q[ZZ], u);
    project(Q[ZZ], Q[YY], u);
    rvec_dec(Q[ZZ], u);

    unitv(Q[XX], Q[XX]);
    unitv(Q[YY], Q[YY]);
    unitv(Q[ZZ], Q[ZZ]);

    tmmul(A, Q, R);

}

void FiniteGrid::Impl::set_unit_cell()
{
    svmul(1./extend_[XX], cell_[XX], unit_cell_[XX]);
    svmul(1./extend_[YY], cell_[YY], unit_cell_[YY]);
    svmul(1./extend_[ZZ], cell_[ZZ], unit_cell_[ZZ]);
}

FiniteGrid::Impl::Impl()
{
    // to avoid MSBuild Error C2536, use this intialization instead of braced list
    for (size_t i = XX; i <= ZZ; i++)
    {
        for (size_t j = XX; j <= ZZ; j++)
        {
            cell_[i][j] = 0;
        }
    }
    cell_[XX][XX] = 1;
    cell_[YY][YY] = 1;
    cell_[ZZ][ZZ] = 1;
};

FiniteGrid::Impl::~Impl()
{
};

/********************************************************************
 * FiniteGrid
 */
FiniteGrid::FiniteGrid() : impl_(new FiniteGrid::Impl())
{
};


FiniteGrid::~FiniteGrid()
{
};

void FiniteGrid::set_extend(IVec extend)
{
    impl_->extend_ = extend;
}
IVec FiniteGrid::extend()
{
    return impl_->extend_;
}

real FiniteGrid::grid_cell_volume()
{
    return det(impl_->cell_)/(real) num_gridpoints();
}

int FiniteGrid::ndx3d_to_ndx1d(IVec i_xyz)
{
    if (i_xyz[XX] >= impl_->extend_[XX] || i_xyz[YY] >= impl_->extend_[YY] || i_xyz[ZZ] >= impl_->extend_[ZZ])
    {
        return -1;
    }
    return i_xyz[XX] + impl_->extend_[XX] * i_xyz[YY]  + impl_->extend_[XX] * impl_->extend_[YY] * i_xyz[ZZ];
}

IVec FiniteGrid::ndx1d_to_ndx3d(int i)
{
    IVec result;
    result[XX] = (i % impl_->extend_[XX]) % impl_->extend_[YY];
    result[YY] = (i / impl_->extend_[XX]) % impl_->extend_[YY];
    result[ZZ] = (i / impl_->extend_[XX]) / impl_->extend_[YY];
    return result;
}

void FiniteGrid::set_translation(RVec translate)
{
    impl_->translate_ = translate;
}
RVec FiniteGrid::translation()
{
    return impl_->translate_;
}

size_t FiniteGrid::num_gridpoints()
{
    return impl_->extend_[XX]*impl_->extend_[YY]*impl_->extend_[ZZ];
}

RVec FiniteGrid::cell_lengths()
{
    return {
               norm(impl_->cell_[XX]), norm(impl_->cell_[YY]), norm(impl_->cell_[ZZ])
    };
}

RVec FiniteGrid::cell_angles()
{
    return {
               impl_->deg_angle(impl_->cell_[XX], impl_->cell_[ZZ]),
               impl_->deg_angle(impl_->cell_[XX], impl_->cell_[YY]),
               impl_->deg_angle(impl_->cell_[YY], impl_->cell_[ZZ])
    };
}


void FiniteGrid::set_cell(RVec length, RVec angle)
{


    real cos_beta  = cos(M_PI*angle[YY]/180.);
    real cos_gamma = cos(M_PI*angle[ZZ]/180.);
    real sin_gamma = sin(M_PI*angle[ZZ]/180.);

    impl_->cell_[XX][XX] = length[XX];
    impl_->cell_[XX][YY] = 0;
    impl_->cell_[XX][ZZ] = 0;

    impl_->cell_[YY][XX] = length[YY] * cos_gamma;
    impl_->cell_[YY][YY] = length[YY] * sin_gamma;
    impl_->cell_[YY][ZZ] = 0;

    impl_->cell_[ZZ][XX] = cos_beta;
    impl_->cell_[ZZ][YY] = ( cos(M_PI*angle[XX]/180.) - cos_beta * cos_gamma ) / sin_gamma;
    impl_->cell_[ZZ][ZZ] = sqrt(1-impl_->cell_[ZZ][XX]*impl_->cell_[ZZ][XX]-impl_->cell_[ZZ][YY]*impl_->cell_[ZZ][YY]);

    impl_->cell_[ZZ][XX] *= length[ZZ];
    impl_->cell_[ZZ][YY] *= length[ZZ];
    impl_->cell_[ZZ][ZZ] *= length[ZZ];

    if ((impl_->extend_[XX] > 0) && (impl_->extend_[YY] > 0) && (impl_->extend_[ZZ] > 0) )
    {
        impl_->set_unit_cell();
    }

};

RVec FiniteGrid::gridpoint_coordinate(IVec i)
{
    RVec result;
    mvmul(impl_->unit_cell_, RVec(i[XX], i[YY], i[ZZ]), result);
    return result;
};

RVec FiniteGrid::gridpoint_coordinate(int i)
{
    return gridpoint_coordinate(ndx1d_to_ndx3d(i));
}

void FiniteGrid::rotation(matrix Q)
{
    matrix R;
    impl_->QRDecomposition(impl_->cell_, Q, R);
}

void FiniteGrid::copy_grid(FiniteGrid &grid)
{
    copy_mat(grid.impl_->cell_, this->impl_->cell_);
    copy_ivec(grid.impl_->extend_, this->impl_->extend_);
    copy_rvec(grid.impl_->translate_, this->impl_->translate_);
    copy_mat(grid.impl_->unit_cell_, this->impl_->unit_cell_);
}

/********************************************************************
 * GridReal::Impl
 */

/*! \internal \brief
 * Private implementation class for GridReal.
 *
 */
class GridReal::Impl
{
    public:
        Impl();
        ~Impl();

        std::vector<real> data_; //!< The data on the grid, represented as one-dimensional vector
};

GridReal::Impl::Impl()
{

};

GridReal::Impl::~Impl()
{

};


/********************************************************************
 * GridReal
 */
real GridReal::min()
{
    return *std::min(std::begin(impl_->data_), std::end(impl_->data_));
};

real GridReal::max()
{
    return *std::max(std::begin(impl_->data_), std::end(impl_->data_));
};

real GridReal::mean()
{
    return std::accumulate(std::begin(impl_->data_), std::end(impl_->data_), 0.)/ static_cast<float>(impl_->data_.size());
};

real GridReal::var()
{
    real data_mean = mean();
    return std::accumulate(std::begin(impl_->data_), std::end(impl_->data_), 0.,
                           [ = ](const real &a, real b){return a + (b - data_mean)*(b - data_mean); } );
}

real GridReal::rms()
{
    return sqrt(var());
};

size_t GridReal::data_size()
{
    return impl_->data_.size();
};

std::pair<std::vector<real>::iterator, std::vector<real>::iterator> GridReal::z_section(int section)
{
    std::vector<real>::iterator section_begin = impl_->data_.begin()+extend()[XX]*extend()[YY]*section;
    std::vector<real>::iterator section_end   = impl_->data_.begin()+extend()[XX]*extend()[YY]*(section+1);
    return std::pair<std::vector<real>::iterator, std::vector<real>::iterator> (section_begin, section_end);
}
std::pair<std::vector<real>::iterator, std::vector<real>::iterator> GridReal::zy_row(int section, int row)
{
    std::vector<real>::iterator section_begin = impl_->data_.begin()+extend()[XX]*extend()[YY]*section + extend()[XX]*row;
    std::vector<real>::iterator section_end   =   impl_->data_.begin()+extend()[XX]*extend()[YY]*section + extend()[XX]*(row+1);
    return std::pair<std::vector<real>::iterator, std::vector<real>::iterator> (section_begin, section_end);
};

void GridReal::resize()
{
    impl_->data_.resize(num_gridpoints());
}
void GridReal::normalize()
{
    real scale =  this->grid_cell_volume() / std::accumulate(impl_->data_.begin(), impl_->data_.end(), 0.);
    std::for_each(impl_->data_.begin(), impl_->data_.end(), [ = ] (real &v) {v *= scale; });
}

real &GridReal::at(IVec index)
{
    return impl_->data_.at(ndx3d_to_ndx1d(index));
};

std::vector<real> &GridReal::data()
{
    return impl_->data_;
}

void GridReal::add_offset(real value)
{
    std::for_each(impl_->data_.begin(), impl_->data_.end(), [ = ](real &datum){ datum += value; });
}

GridReal::GridReal() : impl_(new GridReal::Impl())
{
}

GridReal::~GridReal()
{
}

void GridReal::copy_grid(FiniteGrid &grid)
{
    FiniteGrid::copy_grid(grid);
    resize();
}

} //namespace grid_data

} //namespace gmx
