#include "finite3dlatticeindices.h"

#include <cmath>
#include <string>
#include "gromacs/utility/exceptions.h"
namespace gmx
{

Finite3DLatticeIndices::Finite3DLatticeIndices(const Finite3DLatticeIndices &other)
{
    extend_          = other.extend_;
};

Finite3DLatticeIndices::Finite3DLatticeIndices(const std::vector<int> &extend)
{
    setExtend(extend);
}

void Finite3DLatticeIndices::setExtend(const std::vector<int> &extend)
{
    if (!allNonNegative(extend))
    {
        GMX_THROW(RangeError("Lattice extend must be larger or equal zero in all dimensions but is "+std::to_string(extend[XX]) + " , " + std::to_string(extend[YY])+ " , " + std::to_string(extend[ZZ]) + " ."));
    }
    extend_          = extend;
}

const std::vector<int> &Finite3DLatticeIndices::getExtend() const { return extend_; }

void Finite3DLatticeIndices::multiplyExtend(const RVec factor)
{
    for (int i = XX; i <= ZZ; i++)
    {
        extend_[i] = std::ceil(extend_[i] * factor[i]);
    }
}

std::vector<int> Finite3DLatticeIndices::getLatticeIndexFromLinearIndex(int linearIndex) const
{
    std::vector<int> result;
    result[XX] = (linearIndex % extend_[XX]) % extend_[YY];
    result[YY] = (linearIndex / extend_[XX]) % extend_[YY];
    result[ZZ] = (linearIndex / extend_[XX]) / extend_[YY];
    return result;
}

std::size_t Finite3DLatticeIndices::getLinearIndexFromLatticeIndex(const std::vector<int> &latticeIndex) const
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
    if ((dimension < 0) || (dimension > ZZ))
    {
        GMX_THROW(RangeError("Lattice dimension must be between "+std::to_string(XX) + " and " + std::to_string(YY)+ " , but is " + std::to_string(dimension) + " instead ."));
    }

    return ((latticeIndex >= 0) && (latticeIndex < extend_[dimension]));
};

bool Finite3DLatticeIndices::inLattice(const std::vector<int> &latticeIndex) const
{
    if (!allNonNegative(latticeIndex))
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

bool Finite3DLatticeIndices::allNonNegative(const std::vector<int> &latticeIndex) const
{
    return (latticeIndex[XX] >= 0) && (latticeIndex[YY] >= 0) && (latticeIndex[ZZ] >= 0);
}

}
