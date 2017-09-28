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

#include "realfieldmeasure.h"
namespace gmx
{

RealFieldMeasure::RealFieldMeasure(const Field<real> &realfield)
    : realfield_ {realfield}
{};

real RealFieldMeasure::sum() const
{
    return std::accumulate(realfield_.cbegin(), realfield_.cend(), 0.);
}

/*! \brief The minimum grid data value. */
real RealFieldMeasure::min() const
{
    return *std::min_element(realfield_.cbegin(), realfield_.cend());
};
/*! \brief The maximum grid data value. */
real RealFieldMeasure::max() const
{
    return *std::max_element(realfield_.cbegin(), realfield_.cend());
};
/*! \brief The mean grid data value. */
real RealFieldMeasure::mean() const
{
    return sum() / static_cast<float>(realfield_.getNumLatticePoints());
};
/*! \brief
 * The size of the griddata in bytes.
 */
size_t RealFieldMeasure::data_size() const
{
    return realfield_.access().data().size();
};

real RealFieldMeasure::sumOfSquareDeviation() const
{
    real data_mean = mean();
    return std::accumulate(realfield_.cbegin(), realfield_.cend(), 0.,
                           [data_mean](const real &a, real b) {
                               return a + gmx::square(b - data_mean);
                           });
}

real RealFieldMeasure::normSquared() const
{
    return std::accumulate(
            realfield_.cbegin(), realfield_.cend(), 0.,
            [](const real &a, const real &b) { return a + square(b); });
}

real RealFieldMeasure::norm() const { return sqrt(normSquared()); }
/*! \brief The variance of data values.
 *
 * Two pass algorithm, where mean is calculated fist, due to large expected
 * amount of data. (Typical data size=1,000,000)
 */
real RealFieldMeasure::var() const
{
    return sumOfSquareDeviation() / realfield_.getNumLatticePoints();
}
/*! \brief The root mean square deviation of grid data values. */
real RealFieldMeasure::rms() const { return sqrt(var()); };

RVec RealFieldMeasure::center_of_mass() const
{
    rvec weighted_grid_coordinate;
    RVec com = {0, 0, 0};
    for (int i = 0; i < realfield_.getNumLatticePoints(); i++)
    {
        svmul(realfield_.access().data()[i], realfield_.gridpoint_coordinate(i),
              weighted_grid_coordinate);
        rvec_inc(com, weighted_grid_coordinate);
    }
    svmul(1. / (RealFieldMeasure(*this).sum()), com, com);
    return com;
}

std::string RealFieldMeasure::to_string() const
{
    std::string result;
    result += "------- real number grid -------\n";
    result += realfield_.print();
    result += "  min  :" + std::to_string(min()) + "\n";
    result += "  max  :" + std::to_string(max()) + "\n";
    result += "  mean :" + std::to_string(mean()) + "\n";
    result += "  var  :" + std::to_string(var()) + "\n";
    result += "  rms  :" + std::to_string(rms()) + "\n";
    result += "\n----- end real number grid -----\n\n";
    return result;
}
real RealFieldMeasure::entropy() const
{
    auto p    = realfield_.cbegin();
    int  size = realfield_.getNumLatticePoints();
    real sum  = 0;
#pragma omp parallel for num_threads(std::max( \
    1, gmx_omp_get_max_threads())) reduction(+ : sum) schedule( \
    static, size / std::max(1, gmx_omp_get_max_threads()) + 1)
    for (int i = 0; i < size; ++i)
    {
        if (p[i] > 0)
        {
            sum += p[i] * log(p[i]);
        }
    }
    return -1 * realfield_.grid_cell_volume() * sum;
};

} // namespace gmx
