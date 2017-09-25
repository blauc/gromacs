#include "finitegrid.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/invertmatrix.h"
#include <string>

namespace gmx
{

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
        Impl(const Impl &other);
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

bool FiniteGrid::sameGridInAbsTolerance(const FiniteGrid &other, real tolerance) const
{
    rvec translationDifference;
    rvec_sub(translation(), other.translation(), translationDifference);
    if (norm(translationDifference) > tolerance)
    {
        return false;
    }
    for (int dim = 0; dim <= ZZ; dim++)
    {
        rvec cellDifference;
        rvec_sub(impl_->cell_[dim], other.impl_->cell_[dim], cellDifference);
        if (norm(cellDifference) > tolerance)
        {
            return false;
        }
    }

    return true;
}
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
FiniteGrid::Impl::Impl(const Impl &other)
{
    copy_mat(other.cell_, cell_);
    copy_mat(other.unit_cell_, unit_cell_);
    copy_mat(other.unit_cell_inv_, unit_cell_inv_);
    translate_ = other.translate_;
}

/********************************************************************
 * FiniteGrid
 */
FiniteGrid::FiniteGrid() : impl_(new FiniteGrid::Impl()){};
FiniteGrid::FiniteGrid(const FiniteGrid &other) :  Finite3DLatticeIndices(other), impl_(new FiniteGrid::Impl(*(other.impl_)))
{
};

FiniteGrid::FiniteGrid(FiniteGrid &other) :
    Finite3DLatticeIndices(other), impl_(new FiniteGrid::Impl(*(other.impl_)))
{
};

FiniteGrid &FiniteGrid::operator= (const FiniteGrid &other)
{
    set_extend(other.extend());
    copy_mat(other.impl_->cell_, impl_->cell_);
    copy_mat(other.impl_->unit_cell_, impl_->unit_cell_);
    copy_mat(other.impl_->unit_cell_inv_, impl_->unit_cell_inv_);
    impl_->translate_ = other.impl_->translate_;
    return *this;
}


FiniteGrid::~FiniteGrid() = default;

void FiniteGrid::convertToReciprocalSpace()
{
    invertMatrix(impl_->cell_, impl_->unit_cell_);
    invertMatrix(impl_->unit_cell_, impl_->unit_cell_inv_);
    resetCell();
}


real FiniteGrid::grid_cell_volume() const
{
    return det(impl_->cell_) / (real)num_gridpoints();
}

void FiniteGrid::set_translation(RVec translate)
{
    impl_->translate_ = translate;
}

RVec FiniteGrid::translation() const { return impl_->translate_; }

RVec FiniteGrid::cell_lengths() const
{
    return {
               norm(impl_->cell_[XX]), norm(impl_->cell_[YY]),
               norm(impl_->cell_[ZZ])
    };
}

RVec FiniteGrid::cell_angles() const
{
    return {
               impl_->deg_angle(impl_->cell_[XX], impl_->cell_[ZZ]),
               impl_->deg_angle(impl_->cell_[XX], impl_->cell_[YY]),
               impl_->deg_angle(impl_->cell_[YY], impl_->cell_[ZZ])
    };
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

std::string FiniteGrid::print() const
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

} // gmx
