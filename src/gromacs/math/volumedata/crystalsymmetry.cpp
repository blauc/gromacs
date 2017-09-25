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
/*! \internal \file
 * \brief
 * Implements methods from crystalsymmetry.h
 *
 * \author Christian Blau <cblau@gwdg.de>
 */
#include "gmxpre.h"

#include "crystalsymmetry.h"
/********************************************************************
 * CrystalSymmetry::Impl
 */

/*! \internal \brief
 * Private implementation class for CrystalSymmetry.
 *
 */
class CrystalSymmetry::Impl
{

    public:
        Impl()  = default;
        ~Impl() = default;
        Impl(const Impl &other);
        Impl &operator=(const Impl &other);
        int space_group_ = 1; //!< space group as defined by IUCr conventions (Table
                              //! 12.3.4.1 Standard space-group symbols”, pages
        //! 824-831, International Tables for Crystallography,
        //! Volume A, fifth edition)
};

CrystalSymmetry::Impl::Impl(const Impl &other)
{
    space_group_ = other.space_group_;
};

CrystalSymmetry::Impl &
CrystalSymmetry::Impl::operator=(const Impl &other)
{
    space_group_ = other.space_group_;
    return *this;
};


/********************************************************************
 * CrystalSymmetry
 */

void CrystalSymmetry::set_space_group(int space_group)
{
    impl_->space_group_ = space_group;
}
int CrystalSymmetry::space_group() const { return impl_->space_group_; }

CrystalSymmetry::CrystalSymmetry() : impl_(new CrystalSymmetry::Impl()){};

std::string CrystalSymmetry::print() const
{
    return "---crystal symmetry---\nspace group : " +
           std::to_string(space_group()) + "---\n";
};

CrystalSymmetry::CrystalSymmetry(const CrystalSymmetry &other)
    : impl_ {new CrystalSymmetry::Impl(*other.impl_)}
{};

CrystalSymmetry::~CrystalSymmetry() = default;

CrystalSymmetry &
CrystalSymmetry::operator=(const CrystalSymmetry &other)
{
    *impl_ = *other.impl_;
    return *this;
}
