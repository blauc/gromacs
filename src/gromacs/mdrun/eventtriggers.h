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
/*! \libinternal \file
 * \brief
 * Declares gmx::CallBack.
 *
 * \author Christian Blau <blau@kth.se>
 * \inlibraryapi
 * \ingroup module_mdrun
 */

#ifndef GMX_MDRUN_MESSAGEPUBLISHER_H
#define GMX_MDRUN_MESSAGEPUBLISHER_H

#include <functional>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"

struct t_commrec;

namespace gmx
{

class LocalAtomSetManager;
class SelectionCollection;

/*! \libinternal \brief
 * Signals that the communication record is set up and provides self record.
 */
struct CommunicationIsSetup
{
    //! The communication record that is set up.
    const t_commrec &communicationRecord_;
};

/*!\libinternal \brief
 * Signals that periodic boundary option is chosen and provides it.
 */
struct PeriodicBoundaryConditionOptionIsSetup
{
    int ePBC_; // the periodic boundary options
};

struct BoxIsSetup
{
    BoxIsSetup(const matrix box)
    {
        copy_mat(box, box_);
    }
    matrix box_; // the box
};

/*!\libinternal \brief Provid coordinates of input structure.*/
struct GlobalCoordinatesProvidedOnMaster
{
    ArrayRef<const RVec> coordinates;
};
/* \libinternal \brief
 * Extents base event notify with new event notify and routine to add event listeners.
 *
 * \tparam Message the event to be triggered
 * \tparam BasePublisher a event notify base class to be extended with Message
 *
 * This class triggers call-backs in classes that registersered to be listening to
 * an event by adding a reference to subscribe.
 *
 * \msc
 *   MOD [label = "MdModules"],
 *   CallBack [label="CallBack"],
 *   MessageListener [label="IMessageSubscriber<Message>"],
 *   Message [label = "Message"];
 *   MOD box MessageListener [label = "MdModules owns CallBack and IMessageListeners"];
 *   MOD => CallBack [ label="subscribe(IMessageConsumer<Message>)"];
 *   ...;
 *   --- [label = "Message is triggered"];
 *   ...;
 *   CallBack =>> CA [label="notify()"];
 * \endmsc
 *
 * \note All added IMessageListeners are required to outlive the CallBack.
 */
template <class Message, class BasePublisher>
class CallBack : public BasePublisher
{
    public:
        //! Make base class notifications, subscription and unsubscription available to this class
        using BasePublisher::operator();
        using BasePublisher::subscribe;

        /*! \brief Publisher the event, calling the callback in the MessageListeners.
         * \param[in] event to be triggered
         */
        void operator()(Message message) const
        {
            for (auto &callback : callbackFunction_)
            {
                callback(message);
            }
        }

        /*! \brief Add callback function.
         *
         * Callback functions are distinguished by their call-signature.
         * To allow for different overloaded subscribe functiions, add
         * the unused Message to the call signature.
         * Templating for Msg allows to pick the correct overloaded subscribe function.
         *
         * \tparam Msg the call signature to be subscribed to
         * \tparam Subscriber the class to be subscribed
         *
         * \param[in] subscriber to be called back from this class
         * \param[in, unused] msg disambiguates overloaded functions
         */
        template <typename Msg, class Subscriber>
        typename std::enable_if_t<std::is_same<Msg, Message>::value>
        subscribe(Subscriber *subscriber, Message gmx_unused *msg = nullptr)
        {
            callbackFunction_.emplace_back([subscriber](Message msg) { subscriber->callback(msg); });
        }

    private:
        std::vector < std::function < void(Message)>> callbackFunction_;
};

namespace detail
{

/*! \internal
 * \brief Helper struct to avoid spelling out a nested CallBack
 *        type definition.
 *
 * Instead of
 * CallBack<MessageA, CallBack<MessageB, etc ... >>
 * this allows to write registerCallBack<MessageA, MessageB, ...>::type
 *
 * \tparam Messages all the event types to be registered
 */
template <class ... Messages> struct registerCallBack;

/*! \internal \brief Template specialization to end parameter unpacking recursion.
 */
template <>
struct registerCallBack<>
{
    /*! \internal \brief Do nothing but be base class of CallBack.
     * Required so that using BasePublisher::notify and BasePublisher::subscribe are
     * valid if no events are registered.
     */
    class NoMessage
    {
        public:
            //! Do nothing but provide provideConsumersWithMessage() to CallBack
            void operator()() {}
            //! Do nothing but provide subscribe() to CallBack
            void subscribe() {}
    };
    /*! \brief Define a type even if no events are manged, so that code works
     * with MdModuleMessageManagement that does not actually manage any event.
     */
    using type = NoMessage;
};

/*! \internal \brief Template specialization that assembles the CallBack type
 *                   by taking off the front of the Messages parameter pack.
 *
 * \tparam CurrentMessage front of the template parameter pack
 * \tparam Messages rest of the event types
 */
template <class CurrentMessage, class ... Messages>
struct registerCallBack<CurrentMessage, Messages...>
{
    //! Typing aid that determine the following type
    using next_type = typename registerCallBack<Messages...>::type;
    //! The type of the CallBack
    using type      = CallBack<CurrentMessage, next_type>;
};
}   // namespace detail

//! Registers the event triggers that might occur for MDModules
using MdModuleCallBack = detail::registerCallBack<
            CommunicationIsSetup,
            LocalAtomSetManager *,
            SelectionCollection *,
            GlobalCoordinatesProvidedOnMaster,
            PeriodicBoundaryConditionOptionIsSetup,
            BoxIsSetup
            >::type;

} // namespace gmx

#endif
