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
#ifndef GMX_MATH_GRIDDATAPROPERTIES_H
#define GMX_MATH_GRIDDATAPROPERTIES_H

#include <vector>
#include <numeric>
#include <algorithm>
#include "gromacs/math/functions.h"

namespace gmx
{

template <class T> class BasicGridDataProperties
{
    public:
        explicit BasicGridDataProperties(const std::vector<T> &data) : data_ {data}
        {}
        /*! \brief The sum of all grid values. */
        T sum() const
        {
            return std::accumulate(std::begin(data_), std::end(data_), 0.);
        }

        /*! \brief The minimum grid data value. */
        T min() const
        {
            return *std::min_element(std::begin(data_), std::end(data_));
        };
        /*! \brief The maximum grid data value. */
        T max() const
        {
            return *std::max_element(std::begin(data_), std::end(data_));
        };
        /*! \brief The mean grid data value. */
        T mean() const { return sum() / static_cast<float>(data_.size()); };
        /*! \brief
         * The size of the griddata in bytes.
         */
        size_t data_size() const { return data_.size(); };

        T sumOfSquareDeviation() const
        {
            T data_mean        = mean();
            return std::accumulate(
                    std::begin(data_), std::end(data_), 0.,
                    [ data_mean ](const T &a, T b) { return a + gmx::square(b - data_mean); });
        }

        T normSquared() const
        {
            return std::accumulate( std::begin(data_), std::end(data_), 0., [](const T &a, const T &b) { return a + square(b); });
        }

        T norm() const
        {
            return sqrt(normSquared());
        }
        /*! \brief The variance of data values.
         *
         * Two pass algorithm, where mean is calculated fist, due to large expected
         * amount of data. (Typical data size=1,000,000)
         */
        T var() const
        {
            return sumOfSquareDeviation() / static_cast<real>(data_.size());
        }

    private:
        const std::vector<T> &data_;
};

template <class T>
class ScalarGridDataProperties : public BasicGridDataProperties<T>
{
    public:
        explicit ScalarGridDataProperties(const std::vector<T> &data)
            : BasicGridDataProperties<T>{data}, data_ {data} {};
        /*! \brief The root mean square deviation of grid data values. */
        T rms() const { return sqrt(BasicGridDataProperties<T>::var()); };

    private:
        const std::vector<T> &data_;
};

}      // namespace gmx
#endif
