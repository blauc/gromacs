/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
/*! \libinternal \file
 * \brief
 * implements quaterions
 *
 * \author Christian Blau <cblau@gwdg.de>
 */

#include "quaternion.h"

#include <cmath>
#include <numeric>
#include <array>

namespace gmx
{

Quaternion::Quaternion(RVec &x) : q_
{
    {
        0, x[0], x[1], x[2]
    }
}
{}

Quaternion::Quaternion(QVec q) : q_ {q}
{}


Quaternion::Quaternion(RVec direction, real phi)
{
    real sin_phi_half = std::sin(phi/2.0);
    real cos_phi_half = std::cos(phi/2.0); //< avoid type narrowing warning in initializer list
    q_ = QVec {{
                   cos_phi_half,
                   sin_phi_half * direction[0],
                   sin_phi_half * direction[1],
                   sin_phi_half * direction[2]
               }};
    normalize();
}

Quaternion
Quaternion::operator=(const Quaternion &other)
{
    q_ = other.q_;
    return *this;
}

real &
Quaternion::operator[](int i )
{
    return q_[i];
}

Quaternion
Quaternion::operator*(Quaternion other) const
{
    return {{{
                 q_[0] * other[0] - q_[1] * other[1] - q_[2] * other[2] - q_[3] * other[3],
                 q_[0] * other[1] + q_[1] * other[0] + q_[2] * other[3] - q_[3] * other[2],
                 q_[0] * other[2] - q_[1] * other[3] + q_[2] * other[0] + q_[3] * other[1],
                 q_[0] * other[3] + q_[1] * other[2] - q_[2] * other[1] + q_[3] * other[0]
             }}};
}

Quaternion &
Quaternion::operator*=(Quaternion other)
{
    QVec result = {{
                       q_[0] * other[0] - q_[1] * other[1] - q_[2] * other[2] - q_[3] * other[3],
                       q_[0] * other[1] + q_[1] * other[0] + q_[2] * other[3] - q_[3] * other[2],
                       q_[0] * other[2] - q_[1] * other[3] + q_[2] * other[0] + q_[3] * other[1],
                       q_[0] * other[3] + q_[1] * other[2] - q_[2] * other[1] + q_[3] * other[0]
                   }};
    std::copy(std::begin(result), std::end(result), std::begin(q_));
    return *this;
}


void
Quaternion::rotate(RVec &x)
{
    auto result = this->inverse()*(Quaternion {x} **this);
    x = {result[1], result[2], result[3]};
}

void
Quaternion::rotate_backwards(RVec &x)
{
    auto result = *this * (Quaternion {x} *this->inverse());
    x = {result[1], result[2], result[3]};
}

Quaternion
Quaternion::inverse()
{
    return {{{
                 q_[0], -q_[1], -q_[2], -q_[3]
             }}};
}

void
Quaternion::invert()
{
    q_ = {{q_[0], -q_[1], -q_[2], -q_[3]}};
}

void
Quaternion::normalize()
{
    auto norm = this->norm();
    for (auto &component : q_)
    {
        component /= norm;
    }
}

real
Quaternion::norm()
{
    return std::sqrt(std::accumulate(begin(q_), end(q_), 0., [](const real squaresum, real component){return squaresum + component*component; }));
}


} /* gmx */
