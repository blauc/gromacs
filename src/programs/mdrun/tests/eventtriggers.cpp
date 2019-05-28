/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
 * Tests MdModuleResourceProvider
 *
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include <gmock/gmock.h>

#include "gromacs/mdrun/eventtriggers.h"

namespace gmx
{

namespace
{

class MockCommrecConsumer
{
    public:
        void callback(CommunicationIsSetup /*cr*/)
        {
            notifiedCommrec_ = true;
        };

        bool notifiedCommrec() { return notifiedCommrec_; }

    private:
        bool notifiedCommrec_ = false;
};

class MockLocalAtomSetManagerConsumer final
{
    public:
        void callback(LocalAtomSetManager * /* manager */)
        {
            notifiedlocalAtomSetManager_ = true;
        };

        bool notifiedLocalAtomSetManager() { return notifiedlocalAtomSetManager_; }

    private:
        bool notifiedlocalAtomSetManager_ = false;
};

class MockCombinedConsumer final
{
    public:
        void callback(LocalAtomSetManager * /* manager */)
        {
            notifiedlocalAtomSetManager_ = true;
        };

        void callback(CommunicationIsSetup /*cr*/)
        {
            notifiedCommrec_ = true;
        };

        bool notifiedLocalAtomSetManager() { return notifiedlocalAtomSetManager_; }
        bool notifiedCommrec() { return notifiedCommrec_; }

    private:
        bool notifiedlocalAtomSetManager_ = false;
        bool notifiedCommrec_             = false;
};

TEST(MdDModuleCallBackregisterCallBackTest, addCommrecConsumer)
{
    gmx::detail::registerCallBack<CommunicationIsSetup>::type events;
    MockCommrecConsumer commrecConsumer;

    EXPECT_FALSE(commrecConsumer.notifiedCommrec());

    events.subscribe<CommunicationIsSetup>(&commrecConsumer);
    t_commrec * cr = nullptr;
    events({*cr});

    EXPECT_TRUE(commrecConsumer.notifiedCommrec());
}

TEST(MdDModuleCallBackregisterCallBackTest, addLocalAtomSetManagerConsumer)
{
    gmx::detail::registerCallBack<LocalAtomSetManager *>::type events;
    MockLocalAtomSetManagerConsumer localAtomSetManagerConsumer;

    EXPECT_FALSE(localAtomSetManagerConsumer.notifiedLocalAtomSetManager());

    events.subscribe<LocalAtomSetManager *>(&localAtomSetManagerConsumer);
    LocalAtomSetManager * atomSets = nullptr;
    events(atomSets);

    EXPECT_TRUE(localAtomSetManagerConsumer.notifiedLocalAtomSetManager());
}

TEST(MdDModuleCallBackregisterCallBackTest, addTwoDifferentConsumers)
{
    detail::registerCallBack<CommunicationIsSetup, LocalAtomSetManager *>::type events;
    MockLocalAtomSetManagerConsumer localAtomSetManagerConsumer;
    MockCommrecConsumer             commrecConsumer;

    EXPECT_FALSE(localAtomSetManagerConsumer.notifiedLocalAtomSetManager());
    EXPECT_FALSE(commrecConsumer.notifiedCommrec());

    events.subscribe<LocalAtomSetManager *>(&localAtomSetManagerConsumer);
    events.subscribe<CommunicationIsSetup>(&commrecConsumer);

    LocalAtomSetManager * atomSets = nullptr;
    events(atomSets);

    t_commrec *cr = nullptr;
    events(CommunicationIsSetup {*cr});

    EXPECT_TRUE(localAtomSetManagerConsumer.notifiedLocalAtomSetManager());
    EXPECT_TRUE(commrecConsumer.notifiedCommrec());
}

TEST(MdDModuleCallBackregisterCallBackTest, consumerOfTwoResources)
{
    MdModuleCallBack         events;

    MockCombinedConsumer     consumer;

    EXPECT_FALSE(consumer.notifiedLocalAtomSetManager());
    EXPECT_FALSE(consumer.notifiedCommrec());

    // requires a template parameter here, because call is ambiguous otherwise
    events.subscribe<CommunicationIsSetup>(&consumer);
    events.subscribe<LocalAtomSetManager *>(&consumer);

    LocalAtomSetManager * atomSets = nullptr;
    events(atomSets);

    t_commrec *cr = nullptr;
    events(CommunicationIsSetup {*cr});

    EXPECT_TRUE(consumer.notifiedLocalAtomSetManager());
    EXPECT_TRUE(consumer.notifiedCommrec());
}


} // namespace

} // namespace gmx
