#pragma once

#include "MantidDataObjects/EventListBaseFunctionsTemplate.h"

namespace Mantid {
namespace DataObjects {
template <typename T, typename SELF>
class EventListErrorFunctionsTemplate
{
public:
// using BASE = typename EventListErrorFunctionsTemplate<T, SELF>::EventListBaseFunctionsTemplate<T, SELF>;

EventListErrorFunctionsTemplate(std::shared_ptr<std::vector<T>> events): 
EventListBaseFunctionsTemplate<T, SELF>(events){}
  

protected:

    // --------------------------------------------------------------------------
    /** Get the weight error member of all events in a list
     *
     * @param events :: source vector of events
     * @param weightErrors :: vector to fill
     */
    void getWeightErrorsHelper(const std::vector<T> &events, std::vector<double> &weightErrors) {
    weightErrors.clear();
    weightErrors.reserve(events.size());
    std::transform(events.cbegin(), events.cend(), std::back_inserter(weightErrors),
                    [](const auto &event) { return event.error(); });
    }

    friend T;
    friend SELF;
    EventListErrorFunctionsTemplate() = default;

    inline SELF & as_underlying()
    {
        return static_cast<SELF&>(*this);
    }
};
}
}