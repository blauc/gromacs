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
 * Declares quaternions
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_math
 * \inlibraryapi
 */
#ifndef GMX_MATH_QUATERNION_H_
#define GMX_MATH_QUATERNION_H_

#include "vectypes.h"
#include <array>
#include "gromacs/utility/classhelpers.h"

namespace gmx
{

class Quaternion
{
    public:
        typedef std::array<real, 4> QVec;
        /*! Construct a QVec from orientation direction and rotation angle.
         * \param[in] direction Vector pointing into the orientation direction.
         * \param[in] phi Rotation around direction in rad.
         */
        Quaternion(RVec direction, real phi);
        Quaternion(const RVec &x);
        Quaternion(QVec q);
        real &operator[](int i);
        /*! \brief Invert this QVec.
         * q^(-1) = q*
         */
        Quaternion operator*(Quaternion other) const;
        Quaternion &operator*=(Quaternion other);
        Quaternion operator=(const Quaternion &other);
        void invert();
        Quaternion inverse() const;
        void normalize();
        real norm() const;
        /*! \brief Multiply q from right.
         */
        void rotate(RVec &x) const;
        void rotate_backwards(RVec &x) const;
        RVec shiftedAndOriented(const RVec &x, const RVec &centerOfMass, const RVec &shift) const;
    private:
        QVec q_;
};


}

#endif /* end of include guard: GMX_MATH_QUATERNIONS_H_ */
