#include "finitegrid.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/invertmatrix.h"
#include "gromacs/utility/compare.h"
#include <string>

#include <vector>

namespace gmx
{

bool FiniteGrid::sameGridInAbsTolerance(const FiniteGrid &other) const
{
    auto otherTranslationIterator = std::begin(other.translation_);
    for (auto t : translation_)
    {
        equal_real(*otherTranslationIterator, t, relativeTolerance_, absoluteTolerance_);
        ++otherTranslationIterator;
    }
    return cell_.isSameWithinTolerance(other.getCell(), relativeTolerance_, absoluteTolerance_);
}

void FiniteGrid::setUnitCell_()
{
    unit_cell_     = cell_.scale({{real(1. / lattice_.getExtend()[XX]), real(1. / lattice_.getExtend()[YY]), real(1. / lattice_.getExtend()[ZZ])}} );
}

void FiniteGrid::scaleCell(const OrthogonalBasis<DIM>::NdVector &scale)
{
    cell_ = cell_.scale(scale);
    setUnitCell_();
}

void FiniteGrid::resetCell()
{
    cell_ = unit_cell_.scale({{real(lattice_.getExtend()[XX]), real(lattice_.getExtend()[YY]), real(lattice_.getExtend()[ZZ])}});
};

void FiniteGrid::convertToReciprocalSpace()
{
    unit_cell_     = cell_.inverse();
    resetCell();
}

void FiniteGrid::setTranslation(const OrthogonalBasis<DIM>::NdVector &translate)
{
    translation_ = translate;
}

ColumnMajorLattice<DIM>::MultiIndex FiniteGrid::coordinateToFloorMultiIndex(const OrthogonalBasis<DIM>::NdVector &x) const
{
    auto realValuedIndex = coordinateToRealGridIndex(x);
    ColumnMajorLattice<DIM>::MultiIndex result;
    std::transform(std::begin(realValuedIndex), std::end(realValuedIndex), std::begin(result), [](real x){return static_cast<int>(floor(x)); });
    return result;
}

OrthogonalBasis<DIM>::NdVector FiniteGrid::coordinateToRealGridIndex(const OrthogonalBasis<DIM>::NdVector &x) const
{
    OrthogonalBasis<DIM>::NdVector x_shifted;
    std::transform(std::begin(x), std::end(x), std::begin(translation_), std::begin(x_shifted), [](real x, real t){return x-t; });
    return unit_cell_.transformIntoBasis(x_shifted);
}

OrthogonalBasis<DIM>::NdVector FiniteGrid::multiIndexToCoordinate(const ColumnMajorLattice<DIM>::MultiIndex &i) const
{
    auto result = unit_cell_.transformFromBasis({{real(i[XX]), real(i[YY]), real(i[ZZ])}});
    std::transform(std::begin(translation_), std::end(translation_), std::begin(result), std::begin(result), [](real a, real b){return a+b; } );
    return result;
};

void FiniteGrid::setCell(const OrthogonalBasis<DIM> &cell)
{
    cell_ = cell;
    setUnitCell_();
}


const ColumnMajorLattice<DIM> FiniteGrid::getLattice() const
{
    return lattice_;
}

void FiniteGrid::setLatticeAndRescaleCell(const ColumnMajorLattice<DIM> &lattice)
{
    lattice_ = lattice;
    resetCell();
}

OrthogonalBasis<DIM> FiniteGrid::getCell() const
{
    return cell_;
}

OrthogonalBasis<DIM> FiniteGrid::getUnitCell() const
{
    return unit_cell_;
};


std::string FiniteGrid::to_string() const
{
    std::string result("\n  ------- finite grid -------\n");
    result += "    getExtend       : " + std::to_string(lattice_.getExtend()[0]) + " " +
        std::to_string(lattice_.getExtend()[1]) + " " + std::to_string(lattice_.getExtend()[2]) +
        "\n";
    result += "    ngridpoints  : " + std::to_string(lattice_.getNumLatticePoints()) + "\n";
    result += "    translation  : " + std::to_string(translation_[0]) + " " +
        std::to_string(translation_[1]) + " " +
        std::to_string(translation_[2]) + "\n";
    result += "    cell_lengths : " + std::to_string(cell_.basisVectorLength(0)) + " " +
        std::to_string(cell_.basisVectorLength(1)) + " " +
        std::to_string(cell_.basisVectorLength(2)) + "\n";
    result += "    V_cell       : " + std::to_string(unit_cell_.volume()) + "\n";
    return result + "  ----- end finite grid -----\n\n";
}

} // gmx
