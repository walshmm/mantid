#include "MantidDataObjects/EventListWeightedEventNoTime.h"


namespace Mantid::DataObjects {
using Types::Core::DateAndTime;
using Types::Event::TofEvent;
using namespace Mantid::API;

bool EventListWeightedEventNoTime::equals(const EventListBase &rhs, const double tolTof, const double tolWeight,
                       const int64_t tolPulse) const {
    if (this->getNumberEvents() != rhs.getNumberEvents())
        return false;
    if (this->eventType != rhs.eventType)
        return false;

    // loop over the events
    size_t numEvents = this->getNumberEvents();
    for (size_t i = 0; i < numEvents; ++i) {
    if (!((WeightedEventNoTime)this->events[i]).equals(rhs.events[i], tolTof, tolWeight))
        return false;
    }                         
}

WeightedEvent EventListWeightedEventNoTime::getEvent(size_t event_number) {
    return WeightedEvent(events[event_number].tof(), 0, 0, 0);
}

void EventListWeightedEventNoTime::sortTimeAtSample(const double &tofFactor, const double &tofShift, bool forceResort) const {
    // Check pre-cached sort flag.
    if (this->order == TIMEATSAMPLE_SORT && !forceResort)
        return;

    // Avoid sorting from multiple threads
    std::lock_guard<std::mutex> _lock(m_sortMutex);
    // If the list was sorted while waiting for the lock, return.
    if (this->order == TIMEATSAMPLE_SORT && !forceResort)
        return;

    // Perform sort.

    CompareTimeAtSample<TofEventNoTime> comparitor(tofFactor, tofShift);
    tbb::parallel_sort(events.begin(), events.end(), comparitor);
    // Save the order to avoid unnecessary re-sorting.
    this->order = TIMEATSAMPLE_SORT;    
}



}
