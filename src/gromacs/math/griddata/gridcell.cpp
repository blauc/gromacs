#ifndef GMX_MATH_GRIDCELL
#define GMX_MATH_GRIDCELL

#include "gridcell.h"

#include "gromacs/math/vec.h"
#include "gromacs/math/invertmatrix.h"
namespace gmx
{

namespace
{


/*! \brief
 * Projects vector a onto vector b, stores this projection in c. Conveninence
 * function for Gram-Schmidt method.
 */
void project(const rvec a, const rvec b, rvec c)
{
    svmul(iprod(a, b) / iprod(a, a), a, c);
}


/*! \brief
 * Gram-Schmidt orthonormalisation.
 *
 * The numerical shortcommings of Gram-Schmidt
 * should not matter for three dimensions and cell matrices.
 */
void QRDecomposition(const matrix A, matrix Q, matrix R)
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

}

GridCell GridCell::scale(const RVec scale) const
{
    GridCell result;
    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            result.cell_[i][j] = scale[i] * cell_[i][j];
        }
    }
    return result;
}

const matrix &GridCell::getMatrix() const
{
    return cell_;
}

void GridCell::rotation(matrix Q) const
{
    matrix R;
    QRDecomposition(cell_, Q, R);
}

bool GridCell::spacing_is_same_xyz() const
{
    return (std::abs(norm2(cell_[XX]) -
                     norm2(cell_[YY])) < 1e-5) &&
           (std::abs(norm2(cell_[XX]) -
                     norm2(cell_[ZZ])) < 1e-5);
};

bool GridCell::rectangular() const
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

RVec GridCell::cell_lengths() const
{
    return {
               norm(cell_[XX]), norm(cell_[YY]),
               norm(cell_[ZZ])
    };
}

RVec GridCell::cell_angles() const
{
    return {
               deg_angle(cell_[XX], cell_[ZZ]),
               deg_angle(cell_[XX], cell_[YY]),
               deg_angle(cell_[YY], cell_[ZZ])
    };
}


/*! \brief
 * Angle between two vectors in degree.
 */
real GridCell::deg_angle(const rvec a, const rvec b) const
{
    return gmx_angle(a, b) * 180.0 / M_PI;
}

real GridCell::volume() const
{
    return det(cell_);
};

const rvec &GridCell::operator[](int i) const
{
    return cell_[i];
}


GridCell GridCell::inverse()
{
    matrix   inv;
    invertMatrix(cell_, inv);
    GridCell result;
    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            result.cell_[i][j] = inv[i][j];
        }
    }
    return result;
};

void GridCell::set_cell(RVec length, RVec angle)
{

    real cos_beta  = cos(M_PI * angle[YY] / 180.);
    real cos_gamma = cos(M_PI * angle[ZZ] / 180.);
    real sin_gamma = sin(M_PI * angle[ZZ] / 180.);

    cell_[XX][XX] = length[XX];
    cell_[XX][YY] = 0;
    cell_[XX][ZZ] = 0;

    cell_[YY][XX] = length[YY] * cos_gamma;
    cell_[YY][YY] = length[YY] * sin_gamma;
    cell_[YY][ZZ] = 0;

    cell_[ZZ][XX] = cos_beta;
    cell_[ZZ][YY] =
        (cos(M_PI * angle[XX] / 180.) - cos_beta * cos_gamma) / sin_gamma;
    cell_[ZZ][ZZ] = sqrt(1 - cell_[ZZ][XX] * cell_[ZZ][XX] -
                         cell_[ZZ][YY] * cell_[ZZ][YY]);

    cell_[ZZ][XX] *= length[ZZ];
    cell_[ZZ][YY] *= length[ZZ];
    cell_[ZZ][ZZ] *= length[ZZ];


};


}

#endif
