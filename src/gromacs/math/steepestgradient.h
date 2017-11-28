/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
#ifndef GMX_STEEPESTGRADIENT_H_
#define GMX_STEEPESTGRADIENT_H_

#include <limits>
#include <string>
#include <vector>

class ISteepestGradientBase
{
    protected:
        virtual void stash()                        = 0;
        virtual void applyGradientToStash()         = 0;
        virtual void steepestGradientMoveAccepted() = 0;
        virtual void steepestGradientMoveRejected() = 0;
        virtual bool steepestGradientConverged()    = 0;
};

template <class DataType>
class ISteepestGradient : public ISteepestGradientBase
{
    public:
        using ParameterList = std::vector<std::string>;
        virtual float logLikelihood(const DataType &data) const = 0;

        void steepestGradientRun(const DataType &data, ParameterList names,
                                 int maxSteps = 100000)
        {

            steepestGradientIntialize(data, names);

            float newLogLikelihood = std::numeric_limits<float>::min();
            float oldLogLikelihood = logLikelihood(data);

            for (int step = 0; (step < maxSteps) && (!steepestGradientConverged());
                 ++step)
            {

                stash();
                evaluateGradient(data);
                applyGradientToStash();
                newLogLikelihood = logLikelihood(data);

                // backtrack if move along gradient was not successful
                while (newLogLikelihood < oldLogLikelihood)
                {

                    // backtracking also counts as step
                    ++step;

                    // make a smaller step on the stored parameters with the already
                    // calculated gradient
                    steepestGradientMoveRejected();
                    applyGradientToStash();

                    // evaluate likelihood of new position
                    newLogLikelihood = logLikelihood(data);
                }

                // when move succeed, new step sizes may become bolder
                steepestGradientMoveAccepted();

                // after successful descend, update the likelihood history
                oldLogLikelihood = newLogLikelihood;
            }
        }

    private:
        virtual void steepestGradientIntialize(const DataType &data,
                                               ParameterList   parametersToMove) = 0;
        virtual void evaluateGradient(const DataType &data)                      = 0;
};

#endif /* end of include guard: GMX_STEEPESTGRADIENT_H_ */
