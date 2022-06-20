
#include "MantidDataObjects/EventListBaseFunctionsTemplate.h"

namespace Mantid {
namespace DataObjects {
template <typename T, typename SELF>
class EventListErrorFunctionsTemplate : public EventListBaseFunctionsTemplate<T, SELF>
{
  

private:

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
    EventListErrorFunctionsTemplate() = default;

    inline T & as_underlying()
    {
        return static_cast<T&>(*this);
    }
};
}
}