#include "finite3dlatticeindices.h"

#include <cmath>


namespace gmx
{

Finite3DLatticeIndices::Finite3DLatticeIndices() = default;

Finite3DLatticeIndices::Finite3DLatticeIndices(const Finite3DLatticeIndices &other)
{
    extend_          = other.extend_;
    numGridPointsXY_ = other.numGridPointsXY_;
    numGridPoints_   = other.numGridPoints_;
};

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
    for (int i = 0; i <= ZZ; i++)
    {
        extend_[i] = std::ceil(extend_[i] * factor[i]);
    }
    numGridPointsXY_ = extend_[XX] * extend_[YY];
    numGridPoints_   = numGridPointsXY_ * extend_[ZZ];
}

IVec Finite3DLatticeIndices::ndx1d_to_ndx3d(int i) const
{
    IVec result;
    result[XX] = (i % extend_[XX]) % extend_[YY];
    result[YY] = (i / extend_[XX]) % extend_[YY];
    result[ZZ] = (i / extend_[XX]) / extend_[YY];
    return result;
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


int Finite3DLatticeIndices::numGridPointsXY() const
{
    return numGridPointsXY_;
}

int Finite3DLatticeIndices::num_gridpoints() const { return numGridPoints_; }

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

}
