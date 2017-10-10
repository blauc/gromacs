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
    rvec_sub(translation_, other.translation_, translationDifference);
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
void FiniteGrid::setUnitCell_()
{

    RVec scale {
        real(1. / lattice_.getExtend()[XX]), real(1. / lattice_.getExtend()[YY]), real(1. / lattice_.getExtend()[ZZ])
    };
    unit_cell_     = cell_.scale(scale);
}

void FiniteGrid::scaleCell(RVec scale)
{
    cell_ = cell_.scale(scale);
    setUnitCell_();
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
    resetCell();
}

void FiniteGrid::set_translation(RVec translate)
{
    translation_ = translate;
}

std::array<int, 3> FiniteGrid::coordinate_to_gridindex_floor_ivec(const rvec x) const
{
    RVec result = coordinateToRealGridIndex(x);
    return { {
                 static_cast<int>(floor(result[XX])), static_cast<int>(floor(result[YY])), static_cast<int>(floor(result[ZZ]))
             } };
}

RVec FiniteGrid::coordinateToRealGridIndex(const rvec x) const
{
    RVec x_shifted;
    rvec_sub(x, translation_, x_shifted);
    return unit_cell_.coordinateTransformToCellSpace(x_shifted);
}

real FiniteGrid::avg_spacing() const
{
    return (unit_cell_[XX][XX] + unit_cell_[YY][YY] + unit_cell_[ZZ][ZZ]) / 3;
}

RVec FiniteGrid::gridpoint_coordinate(std::array<int, 3> i) const
{
    RVec result;
    result = unit_cell_.coordinateTransformFromCellSpace(RVec(i[XX], i[YY], i[ZZ]));
    rvec_inc(result, translation_);
    return result;
};

RVec FiniteGrid::gridpoint_coordinate(int linearIndex) const
{
    return gridpoint_coordinate(lattice_.vectoriseLinearIndex(linearIndex));
}

void FiniteGrid::makeGridUniform()
{
    if (!cell_.spacing_is_same_xyz())
    {
        RVec scale {
            1, norm(unit_cell_[XX])/norm(unit_cell_[YY]), norm(unit_cell_[XX])/norm(unit_cell_[ZZ])
        };
        cell_ = cell_.scale(scale);
        setUnitCell_();
    }
}

void FiniteGrid::setCell(RVec length, RVec angle)
{
    cell_ = GridCell(length, angle);
    if ((lattice_.getExtend()[XX] > 0) && (lattice_.getExtend()[YY] > 0) && (lattice_.getExtend()[ZZ] > 0))
    {
        setUnitCell_();
    }
}


const ColumnMajorLattice<3> FiniteGrid::getLattice() const
{
    return lattice_;
}

void FiniteGrid::setLattice(const ColumnMajorLattice<3> &lattice)
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
    result += "    translation  : " + std::to_string(translation_[0]) + " " +
        std::to_string(translation_[1]) + " " +
        std::to_string(translation_[2]) + "\n";
    result += "    cell_lengths : " + std::to_string(cell_.cell_lengths()[0]) + " " +
        std::to_string(cell_.cell_lengths()[1]) + " " +
        std::to_string(cell_.cell_lengths()[2]) + "\n";
    result += "    cell_angles  : " + std::to_string(cell_.cell_angles()[0]) + " " +
        std::to_string(cell_.cell_angles()[1]) + " " +
        std::to_string(cell_.cell_angles()[2]) + "\n";
    result += "    V_cell       : " + std::to_string(unit_cell_.volume()) + "\n";
    return result + "  ----- end finite grid -----\n\n";
}

} // gmx
