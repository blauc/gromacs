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
#include "invertmatrix.h"

#include <algorithm>
#include <numeric>
#include <vector>
#include <string>
#include <cmath>

#include "gromacs/fft/parallel_3dfft.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"

#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/exceptions.h"


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

std::string
CrystalSymmetry::print()
{
    return "---crystal symmetry---\nspace group : "+std::to_string(space_group()) + "---\n";
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

        matrix           cell_;          //!< r' = cell_ . (r-translate_)
        matrix           unit_cell_;     //!< cell_/extend_
        matrix           unit_cell_inv_; //!< cell_/extend_
        RVec             translate_;     //!< r' = cell_ . (r-translate_)


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

void
FiniteGrid::multiplyGridPointNumber(const RVec factor)
{
    Finite3DLatticeIndices::multiplyExtend(factor);
    set_unit_cell();
};

void FiniteGrid::set_unit_cell()
{
    svmul(1./extend()[XX], impl_->cell_[XX], impl_->unit_cell_[XX]);
    svmul(1./extend()[YY], impl_->cell_[YY], impl_->unit_cell_[YY]);
    svmul(1./extend()[ZZ], impl_->cell_[ZZ], impl_->unit_cell_[ZZ]);
    invertMatrix(impl_->unit_cell_, impl_->unit_cell_inv_);
}

FiniteGrid::Impl::Impl() : unit_cell_
{
    {
        1, 0, 0
    }, {
        0, 1, 0
    }, {
        0, 0, 1
    }
}, unit_cell_inv_ {{
                       1, 0, 0
                   }, {
                       0, 1, 0
                   }, {
                       0, 0, 1
                   }}, translate_ {
    0, 0, 0
}
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

Finite3DLatticeIndices::Finite3DLatticeIndices()
{
}


Finite3DLatticeIndices::Finite3DLatticeIndices(IVec extend)
{
    set_extend(extend);
}


void Finite3DLatticeIndices::set_extend(IVec extend)
{
    extend_          = extend;
    numGridPointsXY_ = extend[XX] * extend[YY];
    numGridPoints_   = numGridPointsXY_ * extend[ZZ];
}

IVec
Finite3DLatticeIndices::extend() const
{
    return extend_;
}

void
Finite3DLatticeIndices::multiplyExtend(const RVec factor)
{
    for (size_t i = 0; i <= ZZ; i++)
    {
        extend_[i] = std::ceil(extend_[i] * factor[i]);
    }
}

real FiniteGrid::grid_cell_volume()
{
    return det(impl_->cell_)/(real) num_gridpoints();
}

int Finite3DLatticeIndices::ndx3d_to_ndx1d(IVec ndx3D) const
{
    return ndx3D[XX] + extend_[XX] * ndx3D[YY]  + numGridPointsXY_ * ndx3D[ZZ];
}

IVec Finite3DLatticeIndices::ndx1d_to_ndx3d(int i) const
{
    IVec result;
    result[XX] = (i % extend_[XX]) % extend_[YY];
    result[YY] = (i / extend_[XX]) % extend_[YY];
    result[ZZ] = (i / extend_[XX]) / extend_[YY];
    return result;
}

void FiniteGrid::set_translation(RVec translate)
{
    impl_->translate_ = translate;
}

RVec FiniteGrid::translation() const
{
    return impl_->translate_;
}

size_t Finite3DLatticeIndices::numGridPointsXY() const
{
    return numGridPointsXY_;
}


size_t Finite3DLatticeIndices::num_gridpoints() const
{
    return numGridPoints_;
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

void
FiniteGrid::extendCellByUnitCellInXX()
{
    rvec_inc(impl_->cell_[XX], impl_->unit_cell_[XX]);
    set_unit_cell();
}


void
FiniteGrid::extendCellByUnitCellInYY()
{
    rvec_inc(impl_->cell_[YY], impl_->unit_cell_[YY]);
    set_unit_cell();
}

void
FiniteGrid::extendCellByUnitCellInZZ()
{
    rvec_inc(impl_->cell_[ZZ], impl_->unit_cell_[ZZ]);
    set_unit_cell();
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

    if ((extend()[XX] > 0) && (extend()[YY] > 0) && (extend()[ZZ] > 0) )
    {
        set_unit_cell();
    }

};

bool FiniteGrid::rectangular()
{
    RVec angles = cell_angles();
    for (int i = XX; i <= ZZ; ++i)
    {
        if ((angles[i] < 89.999) || (angles[i] > 90.001))
        {
            return false;
        }
    }
    return true;
};

RVec FiniteGrid::unit_cell_XX() const
{
    return impl_->unit_cell_[XX];
}

RVec FiniteGrid::unit_cell_YY() const
{
    return impl_->unit_cell_[YY];
}

RVec FiniteGrid::unit_cell_ZZ() const
{
    return impl_->unit_cell_[ZZ];
}

bool FiniteGrid::spacing_is_same_xyz()
{
    return (std::abs(norm2(impl_->unit_cell_[XX]) - norm2(impl_->unit_cell_[YY])) < 1e-5) && (std::abs(norm2(impl_->unit_cell_[XX])-norm2(impl_->unit_cell_[ZZ])) < 1e-5);
};

IVec FiniteGrid::coordinate_to_gridindex_round_ivec(const rvec x)
{
    RVec result = coordinateToRealGridIndex(x);
    return {
               (int)round(result[XX]), (int)round(result[YY]), (int)round(result[ZZ])
    };
}


IVec FiniteGrid::coordinate_to_gridindex_ceil_ivec(const rvec x)
{
    RVec result = coordinateToRealGridIndex(x);
    return {
               (int)ceil(result[XX]), (int)ceil(result[YY]), (int)ceil(result[ZZ])
    };
}

IVec FiniteGrid::coordinate_to_gridindex_floor_ivec(const rvec x) const
{
    RVec result = coordinateToRealGridIndex(x);
    return {
               (int)floor(result[XX]), (int)floor(result[YY]), (int)floor(result[ZZ])
    };
}

RVec FiniteGrid::coordinateToRealGridIndex(const rvec x) const
{
    RVec result;
    rvec x_shifted;
    rvec_sub(x, impl_->translate_, x_shifted);
    mvmul(impl_->unit_cell_inv_, x_shifted, result);
    return result;
}

real
FiniteGrid::avg_spacing()
{
    return (impl_->unit_cell_[XX][XX] + impl_->unit_cell_[YY][YY]+impl_->unit_cell_[ZZ][ZZ])/3;
}

bool
Finite3DLatticeIndices::inGrid(IVec gridIndex) const
{
    for (size_t dimension = XX; dimension <= ZZ; dimension++)
    {
        if ((gridIndex[dimension] >= extend()[dimension]) || gridIndex[dimension] < 0)
        {
            return false;
        }
    }
    return true;
}

RVec FiniteGrid::gridpoint_coordinate(IVec i) const
{
    RVec result;
    mvmul(impl_->unit_cell_, RVec(i[XX], i[YY], i[ZZ]), result);
    rvec_inc(result, impl_->translate_);
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

void FiniteGrid::copy_grid(const FiniteGrid &grid)
{
    copy_mat(grid.impl_->cell_, this->impl_->cell_);
    Finite3DLatticeIndices::set_extend(grid.extend());
    set_translation(grid.translation());
    set_unit_cell();
}

void
FiniteGrid::makeGridUniform()
{
    if (!spacing_is_same_xyz())
    {
        real l_XX = norm(impl_->unit_cell_[XX]);
        real l_YY = norm(impl_->unit_cell_[YY]);
        real l_ZZ = norm(impl_->unit_cell_[ZZ]);
        svmul(l_XX/l_YY, impl_->cell_[YY], impl_->cell_[YY]);
        svmul(l_XX/l_ZZ, impl_->cell_[ZZ], impl_->cell_[ZZ]);
        set_unit_cell();
    }
}

std::string
FiniteGrid::print()
{
    std::string result("\n  ------- finite grid -------\n");
    result += "    extend       : " + std::to_string(extend()[0]) + " " + std::to_string(extend()[1]) + " " + std::to_string(extend()[2]) + "\n";
    result += "    ngridpoints  : " + std::to_string(num_gridpoints()) + "\n";
    result += "    translation  : " + std::to_string(translation()[0]) + " " + std::to_string(translation()[1]) + " " + std::to_string(translation()[2]) + "\n";
    result += "    cell_lengths : " + std::to_string(cell_lengths()[0]) + " " + std::to_string(cell_lengths()[1]) + " " + std::to_string(cell_lengths()[2]) + "\n";
    result += "    cell_angles  : " + std::to_string(cell_angles()[0]) + " " + std::to_string(cell_angles()[1]) + " " + std::to_string(cell_angles()[2]) + "\n";
    result += "    V_cell       : " + std::to_string(grid_cell_volume()) + "\n";
    return result+"  ----- end finite grid -----\n\n";
}


GridInterpolator::GridInterpolator(const FiniteGrid &basis) : interpolatedGrid_ {std::unique_ptr<GridReal>(new GridReal)}
{
    interpolatedGrid_->copy_grid(basis);
};

/*
 * for each target grid point:
 *      find rational number grid cell index in input grid
 *      use fractional part for weights
 */
std::unique_ptr<GridReal>
GridInterpolator::interpolateLinearly(const GridReal &other)
{
    auto otherAccess            = other.access();
    auto interpolatedGridAccess = interpolatedGrid_->access();

    for (int i_z = 0; i_z < interpolatedGrid_->extend()[ZZ]; ++i_z)
    {
        for (int i_y = 0; i_y < interpolatedGrid_->extend()[YY]; ++i_y)
        {
            for (int i_x = 0; i_x < interpolatedGrid_->extend()[XX]; ++i_x)
            {

                auto r                 = interpolatedGrid_->gridpoint_coordinate({i_x, i_y, i_z});
                auto rIndexInOtherGrid = other.coordinateToRealGridIndex(r);
                auto iIndexInOtherGrid = other.coordinate_to_gridindex_floor_ivec(r);

                auto w_x = rIndexInOtherGrid[XX] - (real)iIndexInOtherGrid[XX];
                auto w_y = rIndexInOtherGrid[YY] - (real)iIndexInOtherGrid[YY];
                auto w_z = rIndexInOtherGrid[ZZ] - (real)iIndexInOtherGrid[ZZ];

                std::array<std::array<std::array<real, 2>, 2>, 2> cube;

                for (int ii_z = 0; ii_z <= 1; ++ii_z)
                {
                    for (int ii_y = 0; ii_y <= 1; ++ii_y)
                    {
                        for (int ii_x = 0; ii_x <= 1; ++ii_x)
                        {
                            auto cube_index = iIndexInOtherGrid;
                            cube_index[XX] += ii_x;
                            cube_index[YY] += ii_y;
                            cube_index[ZZ] += ii_z;
                            if (other.inGrid(cube_index))
                            {
                                cube[ii_x][ii_y][ii_z] = otherAccess.at(cube_index);
                            }
                            else
                            {
                                cube[ii_x][ii_y][ii_z] = 0;
                            }
                        }
                    }
                }

                std::array<std::array<real, 2>, 2> interpolated_x;
                for (int ii_z = 0; ii_z <= 1; ++ii_z)
                {
                    for (int ii_y = 0; ii_y <= 1; ++ii_y)
                    {
                        interpolated_x[ii_y][ii_z] = (1-w_x)*cube[0][ii_y][ii_z] + (w_x)*cube[1][ii_y][ii_z];
                    }
                }

                std::array<real, 2> interpolated_xy;
                for (int ii_z = 0; ii_z <= 1; ++ii_z)
                {
                    interpolated_xy[ii_z] = (1-w_y)*interpolated_x[0][ii_z] + (w_y) * interpolated_x[1][ii_z];
                }

                interpolatedGridAccess.at({i_x, i_y, i_z}) = ((1-w_z) * interpolated_xy[0] + w_z * interpolated_xy[1])/8.0;

            }
        }
    }
    return std::move(interpolatedGrid_);
};

void
GridInterpolator::makeUniform()
{
    interpolatedGrid_->makeGridUniform();
};

GridMeasures::GridMeasures(const GridReal &reference) : reference_ {reference}
{
};

real
GridMeasures::correlate(const GridReal &other, real threshold) const
{
    real   this_mean        = 0;
    real   this_var         = 0;
    real   other_mean       = 0;
    real   other_var        = 0;
    auto   current_value    = reference_.access().data().begin();
    auto   other_curr_value = other.access().data().begin();
    size_t count            = 0;
    real   result           = 0;
    for (size_t i = 0; i < reference_.access().data().size(); i++)
    {
        if ((*other_curr_value > threshold) || (*current_value > threshold))
        {
            this_mean  +=  *current_value;
            other_mean += *other_curr_value;
            ++count;
        }
        ++current_value;
        ++other_curr_value;
    }
    this_mean  /= count;
    other_mean /= count;

    current_value    = reference_.access().data().begin();
    other_curr_value = other.access().data().begin();
    for (size_t i = 0; i < reference_.access().data().size(); i++)
    {
        if ((*other_curr_value > threshold) || (*current_value > threshold))
        {
            this_var  +=  (*current_value )*(*current_value );
            other_var += (*other_curr_value )*(*other_curr_value );
        }
        ++current_value;
        ++other_curr_value;
    }
    this_var  = sqrt(this_var);
    other_var = sqrt(other_var);

    current_value    = reference_.access().data().begin();
    other_curr_value = other.access().data().begin();
    for (size_t i = 0; i < reference_.access().data().size(); i++)
    {
        if ((*other_curr_value > threshold) || (*current_value > threshold))
        {
            result += (*current_value ) * (*other_curr_value );
        }
        ++current_value;
        ++other_curr_value;
    }
    result /= (this_var * other_var );
    return result;
};

real
GridMeasures::getRelativeKLCrossTermSameGrid(const GridReal &other, const std::vector<real> &other_reference) const
{
    auto P = reference_.access().data();
    auto Q = other.access().data();
    // for numerical stability use a reference density
    if (P.size() != Q.size())
    {
        GMX_THROW(APIError("KL-Divergence calculation requires euqally sized input vectors."));
    }
    auto p           = P.begin();
    auto q           = Q.begin();
    auto q_reference = other_reference.begin();
    int  size        = Q.size();
    real sum         = 0;
#pragma omp parallel for num_threads(std::max(1, gmx_omp_nthreads_get(emntDefault))) reduction(+:sum) schedule(static, size/std::max(1, gmx_omp_nthreads_get(emntDefault)) + 1)
    for (int i = 0; i < size; ++i)
    {
        if ((p[i] > 0) && (q[i] > 0 ) && ( q_reference[i] > 0 ))
        {
            sum += p[i] * log(q[i]/(q_reference[i]));
        }
    }
    return sum;
}


real
GridMeasures::getKLSameGrid(const GridReal &other) const
{
    auto P = reference_.access().data();
    auto Q = other.access().data();
    if (P.size() != Q.size())
    {
        GMX_THROW(APIError("KL-CrossTerm calculation requires euqally sized input vectors."));
    }
    auto p           = P.begin();
    auto q           = Q.begin();
    int  size        = Q.size();
    real sum         = 0;
#pragma omp parallel for num_threads(std::max(1, gmx_omp_nthreads_get(emntDefault))) reduction(+:sum) schedule(static, size/std::max(1, gmx_omp_nthreads_get(emntDefault)) + 1)
    for (int i = 0; i < size; ++i)
    {
        if ((p[i] > 0) && (q[i] > 0 ))
        {
            sum += p[i] * log(q[i]/p[i]);
        }
    }
    return sum;
};

real
GridMeasures::getKLCrossTermSameGrid(const GridReal &other) const
{
    auto P = reference_.access().data();
    auto Q = other.access().data();
    if (P.size() != Q.size())
    {
        GMX_THROW(APIError("KL-CrossTerm calculation requires euqally sized input vectors."));
    }
    auto p           = P.begin();
    auto q           = Q.begin();
    int  size        = Q.size();
    real sum         = 0;
#pragma omp parallel for num_threads(std::max(1, gmx_omp_nthreads_get(emntDefault))) reduction(+:sum) schedule(static, size/std::max(1, gmx_omp_nthreads_get(emntDefault)) + 1)
    for (int i = 0; i < size; ++i)
    {
        if ((p[i] > 0) && (q[i] > 0 ))
        {
            sum += p[i] * log(q[i]);
        }
    }
    return sum;
};

real
GridMeasures::entropy() const
{
    auto P           = reference_.access().data();
    auto p           = P.begin();
    int  size        = P.size();
    real sum         = 0;
#pragma omp parallel for num_threads(std::max(1, gmx_omp_nthreads_get(emntDefault))) reduction(+:sum) schedule(static, size/std::max(1, gmx_omp_nthreads_get(emntDefault)) + 1)
    for (int i = 0; i < size; ++i)
    {
        if (p[i] > 0)
        {
            sum += p[i] * log(p[i]);
        }
    }
    return sum;
};

/********************************************************************
 * GridReal
 */

void
GridReal::multiply(real value)
{
    std::for_each(access().data().begin(), access().data().end(), [ value ] (real &v) {v *= value; });
}

real GridReal::normalize()
{
    real integratedDensity =  this->grid_cell_volume() / properties().sum();
    multiply(1/integratedDensity);
    return integratedDensity;
}

void GridReal::add_offset(real value)
{
    std::for_each(access().data().begin(), access().data().end(), [ value ](real &datum){ datum += value; });
}

ScalarGridDataProperties<real>
GridReal::properties()
{
    return ScalarGridDataProperties<real>(access().data());
}

RVec
GridReal::center_of_mass()
{
    rvec weighted_grid_coordinate;
    RVec com = {0, 0, 0};
    for (size_t i = 0; i < num_gridpoints(); i++)
    {
        svmul(access().data()[i], gridpoint_coordinate(i), weighted_grid_coordinate);
        rvec_inc(com, weighted_grid_coordinate);
    }
    svmul(1./(properties().sum()), com, com);
    return com;
}

std::string
GridReal::print()
{
    std::string result;
    result += "------- real number grid -------\n";
    result += FiniteGrid::print();
    result += "  min  :" + std::to_string(properties().min()) + "\n";
    result += "  max  :" + std::to_string(properties().max()) + "\n";
    result += "  mean :" + std::to_string(properties().mean()) + "\n";
    result += "  var  :" + std::to_string(properties().var()) + "\n";
    result += "  rms  :" + std::to_string(properties().rms()) + "\n";
    result += "\n----- end real number grid -----\n\n";
    return result;
}


void GridReal::zero()
{
    std::fill(access().data().begin(), access().data().end(), 0);
};

FourierTransformRealToComplex3D::FourierTransformRealToComplex3D(const Field<real> &input) : input_ {input}
{};

std::unique_ptr < Field < t_complex>>
FourierTransformRealToComplex3D::transform()
{
    gmx_parallel_3dfft_t fft = nullptr;
    real               * rdata;
    t_complex          * cdata;
    MPI_Comm             mpiCommunicatorForFFT[]  = {MPI_COMM_NULL, MPI_COMM_NULL};

    gmx_parallel_3dfft_init(&fft, input_.extend(), &rdata, &cdata, mpiCommunicatorForFFT, false, 1);

    std::copy(input_.access().data().begin(), input_.access().data().end(), rdata);
    gmx_parallel_3dfft_execute(fft, GMX_FFT_REAL_TO_COMPLEX, 0, nullptr);


    std::unique_ptr < Field < t_complex>> output(new Field<t_complex>());
    output->copy_grid(input_);

    std::copy(cdata, cdata+input_.num_gridpoints(), output->access().data().begin());

    gmx_parallel_3dfft_destroy(fft);

    return output;
};

} //namespace grid_data

} //namespace gmx
