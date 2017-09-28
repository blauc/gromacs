#include "finite3dlatticeindices.h"

#include <cmath>

namespace gmx
{

Finite3DLatticeIndices::Finite3DLatticeIndices(const Finite3DLatticeIndices &other)
{
    extend_          = other.extend_;
};

Finite3DLatticeIndices::Finite3DLatticeIndices(const IVec &extend)
{
    setExtend(extend);
}

void Finite3DLatticeIndices::setExtend(const IVec extend)
{
    extend_          = extend;
}

IVec Finite3DLatticeIndices::getExtend() const { return extend_; }

void Finite3DLatticeIndices::multiplyExtend(const RVec factor)
{
    for (int i = XX; i <= ZZ; i++)
    {
        extend_[i] = std::ceil(extend_[i] * factor[i]);
    }
}

IVec Finite3DLatticeIndices::getLatticeIndexFromLinearIndex(int linearIndex) const
{
    IVec result;
    result[XX] = (linearIndex % extend_[XX]) % extend_[YY];
    result[YY] = (linearIndex / extend_[XX]) % extend_[YY];
    result[ZZ] = (linearIndex / extend_[XX]) / extend_[YY];
    return result;
}


int Finite3DLatticeIndices::getLinearIndexFromLatticeIndex(const IVec &latticeIndex) const
{
    auto result = latticeIndex[XX] > -1 ? latticeIndex[XX] : extend_[XX] + latticeIndex[XX];
    result += latticeIndex[YY] > -1 ? extend_[XX] * latticeIndex[YY]
        : extend_[XX] * (extend_[YY] + latticeIndex[YY]);
    result += latticeIndex[ZZ] > -1 ? extend_[XX] *extend_[YY] * latticeIndex[ZZ]
        : extend_[XX] *extend_[YY] * (extend_[ZZ] + latticeIndex[ZZ]);
    return result;
}


int Finite3DLatticeIndices::getNumLatticePointsXY() const
{
    return extend_[XX] *extend_[YY];
}

int Finite3DLatticeIndices::getNumLatticePoints() const { return extend_[XX] *extend_[YY] * extend_[ZZ]; }

bool Finite3DLatticeIndices::inLattice(int latticeIndex, int dimension) const
{
    return ((latticeIndex >= 0) && (latticeIndex < extend_[dimension]));
}

bool Finite3DLatticeIndices::inLattice(const IVec &latticeIndex) const
{
    if ((latticeIndex[XX] < 0) || (latticeIndex[YY] < 0) || (latticeIndex[ZZ] < 0))
    {
        return false;
    }
    if ((latticeIndex[XX] >= extend_[XX]) || (latticeIndex[YY] >= extend_[YY]) ||
        (latticeIndex[ZZ] >= extend_[ZZ]))
    {
        return false;
    }
    return true;
}

}
