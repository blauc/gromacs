#include "finitegrid.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/invertmatrix.h"
#include "gromacs/utility/compare.h"
#include <string>

#include <vector>

namespace gmx
{

// std::string FiniteGrid::to_string() const
// {
//     std::string result = "\tgetExtend     \t: " +
//
//     std::to_string(lattice_.getExtend()[0]) + " " +
//         std::to_string(lattice_.getExtend()[1]) + " " + std::to_string(lattice_.getExtend()[2]) +
//         "\n";
//     result += "\tngridpoints\t: " + std::to_string(lattice_.getNumLatticePoints()) + "\n";
//     result += "\ttranslation\t: " + std::to_string(translation_[0]) + " " +
//         std::to_string(translation_[1]) + " " +
//         std::to_string(translation_[2]) + "\n";
//     result += "\tcell_lengths\t: " + std::to_string(cell_.basisVectorLength(0)) + " " +
//         std::to_string(cell_.basisVectorLength(1)) + " " +
//         std::to_string(cell_.basisVectorLength(2)) + "\n";
//     result += "\tV_cell     \t: " + std::to_string(unit_cell_.volume()) + "\n";
// }

} // gmx
