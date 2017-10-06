#include "finitegrid.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/invertmatrix.h"
#include <string>

#include <vector>

namespace gmx
{

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
        rvec_sub(cell_[dim], other.cell_[dim], cellDifference);
        if (norm(cellDifference) > tolerance)
        {
            return false;
        }
    }

    return true;
}
void FiniteGrid::set_unit_cell()
{

    RVec scale {
        real(1. / lattice_.getExtend()[XX]), real(1. / lattice_.getExtend()[YY]), real(1. / lattice_.getExtend()[ZZ])
    };
    unit_cell_     = cell_.scale(scale);
    unit_cell_inv_ = unit_cell_.inverse();
}

void FiniteGrid::scaleCell(RVec scale)
{
    cell_ = cell_.scale(scale);
    set_unit_cell();
}

void FiniteGrid::resetCell()
{
    RVec scale {
        real(lattice_.getExtend()[XX]), real(lattice_.getExtend()[YY]), real(lattice_.getExtend()[ZZ])
    };
    cell_ = unit_cell_.scale(scale);
};

void FiniteGrid::convertToReciprocalSpace()
{
    unit_cell_     = cell_.inverse();
    unit_cell_inv_ = unit_cell_.inverse();
    resetCell();
}


real FiniteGrid::grid_cell_volume() const
{
    return cell_.volume() / (real)lattice_.getNumLatticePoints();
}

void FiniteGrid::set_translation(RVec translate)
{
    translation_ = translate;
}

RVec FiniteGrid::translation() const { return translation_; }

std::vector<int> FiniteGrid::coordinate_to_gridindex_round_ivec(const rvec x)
{
    RVec result = coordinateToRealGridIndex(x);
    return {
               (int)round(result[XX]), (int)round(result[YY]),
               (int)round(result[ZZ])
    };
}

std::vector<int> FiniteGrid::coordinate_to_gridindex_ceil_ivec(const rvec x)
{
    RVec result = coordinateToRealGridIndex(x);
    return {
               (int)ceil(result[XX]), (int)ceil(result[YY]), (int)ceil(result[ZZ])
    };
}

std::vector<int> FiniteGrid::coordinate_to_gridindex_floor_ivec(const rvec x) const
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
    rvec_sub(x, translation_, x_shifted);
    mvmul(unit_cell_inv_.getMatrix(), x_shifted, result);
    return result;
}

real FiniteGrid::avg_spacing() const
{
    return (unit_cell_[XX][XX] + unit_cell_[YY][YY] +
            unit_cell_[ZZ][ZZ]) /
           3;
}

RVec FiniteGrid::gridpoint_coordinate(std::vector<int> i) const
{
    RVec result;
    mvmul(unit_cell_.getMatrix(), RVec(i[XX], i[YY], i[ZZ]), result);
    rvec_inc(result, translation_);
    return result;
};

RVec FiniteGrid::gridpoint_coordinate(int linearIndex) const
{
    return gridpoint_coordinate(lattice_.getLatticeIndexFromLinearIndex(linearIndex));
}

void FiniteGrid::makeGridUniform()
{
    if (!cell_.spacing_is_same_xyz())
    {
        RVec scale {
            1, norm(unit_cell_[XX])/norm(unit_cell_[YY]), norm(unit_cell_[XX])/norm(unit_cell_[ZZ])
        };
        cell_ = cell_.scale(scale);
        set_unit_cell();
    }
}

void FiniteGrid::setCell(RVec length, RVec angle)
{
    cell_.set_cell(length, angle);
    if ((lattice_.getExtend()[XX] > 0) && (lattice_.getExtend()[YY] > 0) && (lattice_.getExtend()[ZZ] > 0))
    {
        set_unit_cell();
    }
}


const Finite3DLatticeIndices FiniteGrid::getLattice() const
{
    return lattice_;
}

void FiniteGrid::setLattice(const Finite3DLatticeIndices &lattice)
{
    lattice_ = lattice;
}

GridCell FiniteGrid::getCell() const
{
    return cell_;
}

GridCell FiniteGrid::getUnitCell() const
{
    return unit_cell_;
};


std::string FiniteGrid::print() const
{
    std::string result("\n  ------- finite grid -------\n");
    result += "    getExtend       : " + std::to_string(lattice_.getExtend()[0]) + " " +
        std::to_string(lattice_.getExtend()[1]) + " " + std::to_string(lattice_.getExtend()[2]) +
        "\n";
    result += "    ngridpoints  : " + std::to_string(lattice_.getNumLatticePoints()) + "\n";
    result += "    translation  : " + std::to_string(translation()[0]) + " " +
        std::to_string(translation()[1]) + " " +
        std::to_string(translation()[2]) + "\n";
    result += "    cell_lengths : " + std::to_string(cell_.cell_lengths()[0]) + " " +
        std::to_string(cell_.cell_lengths()[1]) + " " +
        std::to_string(cell_.cell_lengths()[2]) + "\n";
    result += "    cell_angles  : " + std::to_string(cell_.cell_angles()[0]) + " " +
        std::to_string(cell_.cell_angles()[1]) + " " +
        std::to_string(cell_.cell_angles()[2]) + "\n";
    result += "    V_cell       : " + std::to_string(grid_cell_volume()) + "\n";
    return result + "  ----- end finite grid -----\n\n";
}

} // gmx
