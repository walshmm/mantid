#pragma once

#include "MantidDataObjects/EventListBaseFunctionsTemplate.h"

namespace Mantid {
namespace DataObjects {
template <typename T, typename SELF>
class EventListWeightFunctionsTemplate : public EventListBaseFunctionsTemplate<T, SELF>
{
  public:


protected:

    // --------------------------------------------------------------------------
/** Get the weight member of all events in a list
 *
 * @param events :: source vector of events
 * @param weights :: vector to fill
 */
void getWeightsHelper(const std::vector<T> &events, std::vector<double> &weights) {
  weights.clear();
  weights.reserve(events.size());
  std::transform(events.cbegin(), events.cend(), std::back_inserter(weights),
                 [](const auto &event) { return event.weight(); });
}

    friend T;
    friend SELF;

    inline SELF & as_underlying()
    {
        return static_cast<SELF&>(*this);
    }
};
}
}