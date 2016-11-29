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

#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"

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


} // namespace volumedata

} // namespace gmx
