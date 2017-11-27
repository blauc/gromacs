#ifndef RFT_STEEPESTGRADIENT
#define RFT_STEEPESTGRADIENT

#include <cmath>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <vector>

namespace gradient
{

using ParameterList = std::vector<std::string>;
template <class DataType> class ISteepestGradient
{
    public:
        virtual void intializeSteepestGradient(const DataType &data,
                                               ParameterList   parametersToMove) = 0;
        virtual float logLikelihood(const DataType &data) = 0;
        virtual void stash() = 0;
        virtual void evaluteGradient(const DataType &data)                = 0;
        virtual void applyGradientToStash()                               = 0;
        virtual void adjustSteepestGradientStepSize(bool moveWasAccepted) = 0;
        virtual bool convergedSteepestGradient()                          = 0;
};

template <class DataType>
void runSteepestGradient(const DataType &data,
                         ISteepestGradient<DataType> *parameters,
                         ParameterList names, int maxSteps)
{

    parameters->intializeSteepestGradient(data, names);

    float newLogLikelihood = std::numeric_limits<float>::min();
    float oldLogLikelihood = parameters->logLikelihood(data);

    for (int step = 0;
         (step < maxSteps) && (!parameters->convergedSteepestGradient());
         ++step)
    {

        parameters->stash();
        parameters->evaluteGradient(data);
        parameters->applyGradientToStash();
        newLogLikelihood = parameters->logLikelihood(data);

        // backtrack if move along gradient was not successful
        while (newLogLikelihood < oldLogLikelihood)
        {

            // backtracking also counts as step
            ++step;

            // make a smaller step on the stored parameters with the already
            // calculated gradient
            parameters->adjustSteepestGradientStepSize(false);
            parameters->applyGradientToStash();

            // evaluate likelihood of new position
            newLogLikelihood = parameters->logLikelihood(data);
        }

        // when move succeed, new step sizes may become bolder
        parameters->adjustSteepestGradientStepSize(true);

        // after successful descend, update the likelihood history
        oldLogLikelihood = newLogLikelihood;
    }
};

}      // namespace gradient
#endif /* end of include guard: RFT_STEEPESTGRADIENT */
