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

#include "gromacs/math/invertmatrix.h"
#include "volumedata.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>
#include <string>
#include <vector>

#include "gromacs/fft/parallel_3dfft.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"

#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxomp.h"

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
        Impl()  = default;
        ~Impl() = default;
        Impl(const Impl &other);
        int space_group_ = 1; //!< space group as defined by IUCr conventions (Table
                              //! 12.3.4.1 Standard space-group symbols”, pages
        //! 824-831, International Tables for Crystallography,
        //! Volume A, fifth edition)
};

CrystalSymmetry::Impl::Impl(const Impl &other)
{
    space_group_ = other.space_group_;
};

/********************************************************************
 * CrystalSymmetry
 */

void CrystalSymmetry::set_space_group(int space_group)
{
    impl_->space_group_ = space_group;
}
int CrystalSymmetry::space_group() { return impl_->space_group_; }

CrystalSymmetry::CrystalSymmetry() : impl_(new CrystalSymmetry::Impl()){};

std::string CrystalSymmetry::print()
{
    return "---crystal symmetry---\nspace group : " +
           std::to_string(space_group()) + "---\n";
};

CrystalSymmetry::CrystalSymmetry(const CrystalSymmetry &other)
    : impl_(new CrystalSymmetry::Impl(*other.impl_)){};

CrystalSymmetry::~CrystalSymmetry() = default;

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

        matrix cell_;          //!< r' = cell_ . (r-translate_)
        matrix unit_cell_;     //!< cell_/extend_
        matrix unit_cell_inv_; //!< cell_/extend_
        RVec   translate_;     //!< r' = cell_ . (r-translate_)

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
         * Projects vector a onto vector b, stores this projection in c. Conveninence
         * function for Gram-Schmidt method.
         */
        void project(const rvec a, const rvec b, rvec c);
};

real FiniteGrid::Impl::deg_angle(rvec a, rvec b)
{
    return gmx_angle(a, b) * 180.0 / M_PI;
}

void FiniteGrid::Impl::project(const rvec a, const rvec b, rvec c)
{
    svmul(iprod(a, b) / iprod(a, a), a, c);
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

void FiniteGrid::multiplyGridPointNumber(const RVec factor)
{
    Finite3DLatticeIndices::multiplyExtend(factor);
    set_unit_cell();
};

void FiniteGrid::set_unit_cell()
{
    svmul(1. / extend()[XX], impl_->cell_[XX], impl_->unit_cell_[XX]);
    svmul(1. / extend()[YY], impl_->cell_[YY], impl_->unit_cell_[YY]);
    svmul(1. / extend()[ZZ], impl_->cell_[ZZ], impl_->unit_cell_[ZZ]);
    invertMatrix(impl_->unit_cell_, impl_->unit_cell_inv_);
}

void FiniteGrid::scaleCell(RVec scale)
{
    svmul(scale[XX], impl_->cell_[XX], impl_->cell_[XX]);
    svmul(scale[YY], impl_->cell_[YY], impl_->cell_[YY]);
    svmul(scale[ZZ], impl_->cell_[ZZ], impl_->cell_[ZZ]);
    set_unit_cell();
}

void FiniteGrid::resetCell()
{
    svmul(extend()[XX], impl_->unit_cell_[XX], impl_->cell_[XX]);
    svmul(extend()[YY], impl_->unit_cell_[YY], impl_->cell_[YY]);
    svmul(extend()[ZZ], impl_->unit_cell_[ZZ], impl_->cell_[ZZ]);
};

FiniteGrid::Impl::Impl()
    : unit_cell_
{
    {
        1, 0, 0
    }, {
        0, 1, 0
    }, {
        0, 0, 1
    }
},
unit_cell_inv_ {{
                    1, 0, 0
                }, {
                    0, 1, 0
                }, {
                    0, 0, 1
                }}, translate_ {
    0, 0, 0
} {
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

FiniteGrid::Impl::~Impl() = default;

/********************************************************************
 * FiniteGrid
 */
FiniteGrid::FiniteGrid() : impl_(new FiniteGrid::Impl()){};

FiniteGrid::~FiniteGrid() = default;

void FiniteGrid::convertToReciprocalSpace()
{
    copy_mat(impl_->cell_, impl_->unit_cell_inv_);
    invertMatrix(impl_->cell_, impl_->unit_cell_);
    svmul(extend()[XX], impl_->unit_cell_[XX], impl_->cell_[XX]);
    svmul(extend()[YY], impl_->unit_cell_[YY], impl_->cell_[YY]);
    svmul(extend()[ZZ], impl_->unit_cell_[ZZ], impl_->cell_[ZZ]);
}

Finite3DLatticeIndices::Finite3DLatticeIndices() = default;

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

IVec Finite3DLatticeIndices::extend() const { return extend_; }

void Finite3DLatticeIndices::multiplyExtend(const RVec factor)
{
    for (size_t i = 0; i <= ZZ; i++)
    {
        extend_[i] = std::ceil(extend_[i] * factor[i]);
    }
    numGridPointsXY_ = extend_[XX] * extend_[YY];
    numGridPoints_   = numGridPointsXY_ * extend_[ZZ];
}

real FiniteGrid::grid_cell_volume()
{
    return det(impl_->cell_) / (real)num_gridpoints();
}

int Finite3DLatticeIndices::ndx3d_to_ndx1d(IVec ndx3D) const
{
    auto result = ndx3D[XX] > -1 ? ndx3D[XX] : extend_[XX] + ndx3D[XX];
    result += ndx3D[YY] > -1 ? extend_[XX] * ndx3D[YY]
        : extend_[XX] * (extend_[YY] + ndx3D[YY]);
    result += ndx3D[ZZ] > -1 ? numGridPointsXY_ * ndx3D[ZZ]
        : numGridPointsXY_ * (extend_[ZZ] + ndx3D[ZZ]);
    return result;
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

RVec FiniteGrid::translation() const { return impl_->translate_; }

size_t Finite3DLatticeIndices::numGridPointsXY() const
{
    return numGridPointsXY_;
}

size_t Finite3DLatticeIndices::num_gridpoints() const { return numGridPoints_; }

RVec FiniteGrid::cell_lengths()
{
    return {
               norm(impl_->cell_[XX]), norm(impl_->cell_[YY]),
               norm(impl_->cell_[ZZ])
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

void FiniteGrid::extendCellByUnitCellInXX()
{
    rvec_inc(impl_->cell_[XX], impl_->unit_cell_[XX]);
    set_unit_cell();
}

void FiniteGrid::extendCellByUnitCellInYY()
{
    rvec_inc(impl_->cell_[YY], impl_->unit_cell_[YY]);
    set_unit_cell();
}

void FiniteGrid::extendCellByUnitCellInZZ()
{
    rvec_inc(impl_->cell_[ZZ], impl_->unit_cell_[ZZ]);
    set_unit_cell();
}

void FiniteGrid::set_cell(RVec length, RVec angle)
{

    real cos_beta  = cos(M_PI * angle[YY] / 180.);
    real cos_gamma = cos(M_PI * angle[ZZ] / 180.);
    real sin_gamma = sin(M_PI * angle[ZZ] / 180.);

    impl_->cell_[XX][XX] = length[XX];
    impl_->cell_[XX][YY] = 0;
    impl_->cell_[XX][ZZ] = 0;

    impl_->cell_[YY][XX] = length[YY] * cos_gamma;
    impl_->cell_[YY][YY] = length[YY] * sin_gamma;
    impl_->cell_[YY][ZZ] = 0;

    impl_->cell_[ZZ][XX] = cos_beta;
    impl_->cell_[ZZ][YY] =
        (cos(M_PI * angle[XX] / 180.) - cos_beta * cos_gamma) / sin_gamma;
    impl_->cell_[ZZ][ZZ] = sqrt(1 - impl_->cell_[ZZ][XX] * impl_->cell_[ZZ][XX] -
                                impl_->cell_[ZZ][YY] * impl_->cell_[ZZ][YY]);

    impl_->cell_[ZZ][XX] *= length[ZZ];
    impl_->cell_[ZZ][YY] *= length[ZZ];
    impl_->cell_[ZZ][ZZ] *= length[ZZ];

    if ((extend()[XX] > 0) && (extend()[YY] > 0) && (extend()[ZZ] > 0))
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

RVec FiniteGrid::unit_cell_XX() const { return impl_->unit_cell_[XX]; }

RVec FiniteGrid::unit_cell_YY() const { return impl_->unit_cell_[YY]; }

RVec FiniteGrid::unit_cell_ZZ() const { return impl_->unit_cell_[ZZ]; }

bool FiniteGrid::spacing_is_same_xyz()
{
    return (std::abs(norm2(impl_->unit_cell_[XX]) -
                     norm2(impl_->unit_cell_[YY])) < 1e-5) &&
           (std::abs(norm2(impl_->unit_cell_[XX]) -
                     norm2(impl_->unit_cell_[ZZ])) < 1e-5);
};

IVec FiniteGrid::coordinate_to_gridindex_round_ivec(const rvec x)
{
    RVec result = coordinateToRealGridIndex(x);
    return {
               (int)round(result[XX]), (int)round(result[YY]),
               (int)round(result[ZZ])
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
               (int)floor(result[XX]), (int)floor(result[YY]),
               (int)floor(result[ZZ])
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

real FiniteGrid::avg_spacing()
{
    return (impl_->unit_cell_[XX][XX] + impl_->unit_cell_[YY][YY] +
            impl_->unit_cell_[ZZ][ZZ]) /
           3;
}

bool Finite3DLatticeIndices::inGrid(int gridIndex, int dimension) const
{
    return ((gridIndex >= 0) && (gridIndex < extend_[dimension]));
}

bool Finite3DLatticeIndices::inGrid(IVec gridIndex) const
{
    if ((gridIndex[XX] < 0) || (gridIndex[YY] < 0) || (gridIndex[ZZ] < 0))
    {
        return false;
    }
    if ((gridIndex[XX] >= extend_[XX]) || (gridIndex[YY] >= extend_[YY]) ||
        (gridIndex[ZZ] >= extend_[ZZ]))
    {
        return false;
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

void FiniteGrid::makeGridUniform()
{
    if (!spacing_is_same_xyz())
    {
        real l_XX = norm(impl_->unit_cell_[XX]);
        real l_YY = norm(impl_->unit_cell_[YY]);
        real l_ZZ = norm(impl_->unit_cell_[ZZ]);
        svmul(l_XX / l_YY, impl_->cell_[YY], impl_->cell_[YY]);
        svmul(l_XX / l_ZZ, impl_->cell_[ZZ], impl_->cell_[ZZ]);
        set_unit_cell();
    }
}

std::string FiniteGrid::print()
{
    std::string result("\n  ------- finite grid -------\n");
    result += "    extend       : " + std::to_string(extend()[0]) + " " +
        std::to_string(extend()[1]) + " " + std::to_string(extend()[2]) +
        "\n";
    result += "    ngridpoints  : " + std::to_string(num_gridpoints()) + "\n";
    result += "    translation  : " + std::to_string(translation()[0]) + " " +
        std::to_string(translation()[1]) + " " +
        std::to_string(translation()[2]) + "\n";
    result += "    cell_lengths : " + std::to_string(cell_lengths()[0]) + " " +
        std::to_string(cell_lengths()[1]) + " " +
        std::to_string(cell_lengths()[2]) + "\n";
    result += "    cell_angles  : " + std::to_string(cell_angles()[0]) + " " +
        std::to_string(cell_angles()[1]) + " " +
        std::to_string(cell_angles()[2]) + "\n";
    result += "    V_cell       : " + std::to_string(grid_cell_volume()) + "\n";
    return result + "  ----- end finite grid -----\n\n";
}

GridInterpolator::GridInterpolator(const FiniteGrid &basis)
    : interpolatedGrid_ {std::unique_ptr<GridReal>(new GridReal)}
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
                        interpolated_x[ii_y][ii_z] =
                            (1 - w_x) * cube[0][ii_y][ii_z] + (w_x)*cube[1][ii_y][ii_z];
                    }
                }

                std::array<real, 2> interpolated_xy;
                for (int ii_z = 0; ii_z <= 1; ++ii_z)
                {
                    interpolated_xy[ii_z] = (1 - w_y) * interpolated_x[0][ii_z] +
                        (w_y)*interpolated_x[1][ii_z];
                }

                interpolatedGridAccess.at({i_x, i_y, i_z}) =
                    ((1 - w_z) * interpolated_xy[0] + w_z * interpolated_xy[1]) / 8.0;
            }
        }
    }
    return std::move(interpolatedGrid_);
};

void GridInterpolator::makeUniform() { interpolatedGrid_->makeGridUniform(); };

GridMeasures::GridMeasures(const GridReal &reference)
    : reference_ {reference}
{};

real GridMeasures::correlate(const GridReal &other, real threshold) const
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
            this_mean  += *current_value;
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
            this_var  += (*current_value) * (*current_value);
            other_var += (*other_curr_value) * (*other_curr_value);
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
            result += (*current_value) * (*other_curr_value);
        }
        ++current_value;
        ++other_curr_value;
    }
    result /= (this_var * other_var);
    return result;
};

real GridMeasures::getRelativeKLCrossTermSameGrid(
        const GridReal &other, const std::vector<real> &other_reference) const
{
    auto P = reference_.access().data();
    auto Q = other.access().data();
    // for numerical stability use a reference density
    if (P.size() != Q.size())
    {
        GMX_THROW(APIError(
                          "KL-Divergence calculation requires euqally sized input vectors."));
    }
    auto p           = P.begin();
    auto q           = Q.begin();
    auto q_reference = other_reference.begin();
    int  size        = Q.size();
    real sum         = 0;
#pragma omp parallel for num_threads(std::max( \
    1, gmx_omp_nthreads_get(emntDefault))) reduction(+ : sum) schedule( \
    static, size / std::max(1, gmx_omp_nthreads_get(emntDefault)) + 1)
    for (int i = 0; i < size; ++i)
    {
        if ((p[i] > 0) && (q[i] > 0) && (q_reference[i] > 0))
        {
            sum += p[i] * log(q[i] / (q_reference[i]));
        }
    }
    return sum;
}

real GridMeasures::getKLSameGrid(const GridReal &other) const
{
    auto P = reference_.access().data();
    auto Q = other.access().data();
    if (P.size() != Q.size())
    {
        GMX_THROW(APIError(
                          "KL-CrossTerm calculation requires euqally sized input vectors."));
    }
    auto p    = P.begin();
    auto q    = Q.begin();
    int  size = Q.size();
    real sum  = 0;
#pragma omp parallel for num_threads(std::max( \
    1, gmx_omp_nthreads_get(emntDefault))) reduction(+ : sum) schedule( \
    static, size / std::max(1, gmx_omp_nthreads_get(emntDefault)) + 1)
    for (int i = 0; i < size; ++i)
    {
        if ((p[i] > 0) && (q[i] > 0))
        {
            sum += p[i] * log(q[i] / p[i]);
        }
    }
    return sum;
};

real GridMeasures::getKLCrossTermSameGrid(const GridReal &other) const
{
    auto P = reference_.access().data();
    auto Q = other.access().data();
    if (P.size() != Q.size())
    {
        GMX_THROW(APIError(
                          "KL-CrossTerm calculation requires euqally sized input vectors."));
    }
    auto p    = P.begin();
    auto q    = Q.begin();
    int  size = Q.size();
    real sum  = 0;
#pragma omp parallel for num_threads(std::max( \
    1, gmx_omp_nthreads_get(emntDefault))) reduction(+ : sum) schedule( \
    static, size / std::max(1, gmx_omp_nthreads_get(emntDefault)) + 1)
    for (int i = 0; i < size; ++i)
    {
        if ((p[i] > 0) && (q[i] > 0))
        {
            sum += p[i] * log(q[i]);
        }
    }
    return sum;
};

real GridMeasures::entropy() const
{
    auto P    = reference_.access().data();
    auto p    = P.begin();
    int  size = P.size();
    real sum  = 0;
#pragma omp parallel for num_threads(std::max( \
    1, gmx_omp_nthreads_get(emntDefault))) reduction(+ : sum) schedule( \
    static, size / std::max(1, gmx_omp_nthreads_get(emntDefault)) + 1)
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

GridReal::GridReal(Field<real> baseField)
    : Field<real>::Field<real>(baseField){};

GridReal::GridReal(GridReal &other)
    : Field<real>::Field<real>(other), CrystalSymmetry::CrystalSymmetry(other)
{

};

void GridReal::multiply(real value)
{
    std::for_each(access().data().begin(), access().data().end(),
                  [value](real &v) { v *= value; });
}

real GridReal::normalize()
{
    real integratedDensity = this->grid_cell_volume() / properties().sum();
    multiply(1 / integratedDensity);
    return integratedDensity;
}

void GridReal::add_offset(real value)
{
    std::for_each(access().data().begin(), access().data().end(),
                  [value](real &datum) { datum += value; });
}

ScalarGridDataProperties<real> GridReal::properties() const
{
    return ScalarGridDataProperties<real>(access().data());
}

RVec GridReal::center_of_mass()
{
    rvec weighted_grid_coordinate;
    RVec com = {0, 0, 0};
    for (size_t i = 0; i < num_gridpoints(); i++)
    {
        svmul(access().data()[i], gridpoint_coordinate(i),
              weighted_grid_coordinate);
        rvec_inc(com, weighted_grid_coordinate);
    }
    svmul(1. / (properties().sum()), com, com);
    return com;
}

std::string GridReal::print()
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

    std::unique_ptr < Field < t_complex>> output(new Field<t_complex>());
    output->copy_grid(input_);
    output->convertToReciprocalSpace();

    std::copy(cdata, cdata + input_.num_gridpoints(), output->access().begin());
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

IVec
FourierTransformComplexToReal3D::firstNonHermitian(real tolerance)
{
    auto extend = input_.extend();
    for (int ix = 0; ix < extend[XX]/2; ++ix)
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

bool
FourierTransformComplexToReal3D::isHermitian(real tolerance)
{
    return firstNonHermitian(tolerance) == IVec({input_.extend()[XX], input_.extend()[YY], input_.extend()[ZZ]});
}

std::unique_ptr < Field < real>> FourierTransformComplexToReal3D::transform() {
    gmx_parallel_3dfft_t fft = nullptr;
    real                *rdata;
    t_complex           *cdata;
    MPI_Comm             mpiCommunicatorForFFT[] = {MPI_COMM_NULL, MPI_COMM_NULL};

    gmx_parallel_3dfft_init(&fft, input_.extend(), &rdata, &cdata,
                            mpiCommunicatorForFFT, false,
                            std::max(1, gmx_omp_nthreads_get(emntDefault)));

    // IVec complex_order, ndata, offset, size;
    IVec size;
    gmx_parallel_3dfft_complex_limits(fft, IVec(), IVec(), IVec(), size);
    auto halfDataEnd =
        input_.access().sectionBegin((input_.extend()[ZZ] + 1) / 2 + 1);
    assert(size[ZZ] * size[YY] * size[XX] ==
           int(std::distance(input_.access().begin(), halfDataEnd)));
    std::copy(input_.access().begin(), halfDataEnd, cdata);

    gmx_parallel_3dfft_execute(fft, GMX_FFT_COMPLEX_TO_REAL, 0, nullptr);

    std::unique_ptr < Field < real>> output(new Field<real>());
    output->copy_grid(input_);
    output->convertToReciprocalSpace();

    std::copy(rdata, rdata + input_.num_gridpoints(),
              output->access().data().begin());

    auto orthoscale = 1 / sqrt(input_.num_gridpoints());

    auto normalizeTransform = [orthoscale](real &value) {
            value *= orthoscale;
        };
    std::for_each(std::begin(output->access().data()),
                  std::end(output->access().data()), normalizeTransform);

    gmx_parallel_3dfft_destroy(fft);

    return output;
};

GaussConvolution::GaussConvolution(const Field<real> &input)
    : input_(input), padded_input_ {nullptr}
{
    extendBeforePadding_ = input_.extend();
};

GaussConvolution &GaussConvolution::pad(RVec paddingFactor)
{
    padded_input_ = DensityPadding(input_).pad(paddingFactor);
    return *this;
}

std::unique_ptr < Field < real>> GaussConvolution::convolute(real sigma) {
    if (padded_input_ != nullptr)
    {
        fourierTransform_ =
            FourierTransformRealToComplex3D(*padded_input_).transform();
    }
    else
    {
        fourierTransform_ = FourierTransformRealToComplex3D(input_).transform();
    }
    auto sigma2                = gmx::square(sigma);
    auto convoluteWithGaussian = [sigma2](t_complex &value, RVec k) {
            auto prefactor = exp(-2 * sigma2 * norm2(k));
            value.re *= prefactor;
            value.im *= prefactor;
        };

    volumedata::ApplyToUnshiftedFourierTransform(*fourierTransform_)
        .apply(convoluteWithGaussian);
    auto result = FourierTransformComplexToReal3D(*fourierTransform_).transform();

    if (padded_input_ != nullptr)
    {
        result = DensityPadding(*result).unpad(extendBeforePadding_);
    }
    return result;
};

DensityPadding::DensityPadding(const Field<real> &toPad) : toPad_ {toPad}
{}

std::unique_ptr < Field < real>> DensityPadding::unpad(IVec unPadExtend) {
    std::unique_ptr < Field < real>> unpadded(new Field<real>);
    unpadded->copy_grid(toPad_);
    unpadded->set_extend(unPadExtend);
    unpadded->resetCell();
    for (int iZZ = 0; iZZ < unPadExtend[ZZ]; ++iZZ)
    {
        for (int iYY = 0; iYY < unPadExtend[YY]; ++iYY)
        {
            for (int iXX = 0; iXX < unPadExtend[XX]; ++iXX)
            {
                unpadded->access().at({iXX, iYY, iZZ}) =
                    toPad_.access().at({iXX, iYY, iZZ});
            }
        }
    }

    return unpadded;
}

std::unique_ptr < Field < real>> DensityPadding::pad(RVec paddingFactor) {
    std::unique_ptr < Field < real>> padded(new Field<real>);
    padded->copy_grid(toPad_);
    padded->multiplyGridPointNumber(paddingFactor);
    padded->scaleCell(paddingFactor);
    std::fill(std::begin(padded->access().data()),
              std::end(padded->access().data()), 0.);
    IVec extend = toPad_.extend();
    for (int iZZ = 0; iZZ < extend[ZZ]; ++iZZ)
    {
        for (int iYY = 0; iYY < extend[YY]; ++iYY)
        {
            for (int iXX = 0; iXX < extend[XX]; ++iXX)
            {
                padded->access().at({iXX, iYY, iZZ}) =
                    toPad_.access().at({iXX, iYY, iZZ});
            }
        }
    }
    return padded;
}

} // namespace volumedata

} // namespace gmx
