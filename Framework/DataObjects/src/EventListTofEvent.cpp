#include "MantidDataObjects/EventListTofEvent.h"


namespace Mantid::DataObjects {
using Types::Core::DateAndTime;
using Types::Event::TofEvent;
using namespace Mantid::API;

bool EventListTofEvent::equals(const EventListBase &rhs, const double tolTof, const double tolWeight,
                       const int64_t tolPulse) const {
    if (this->getNumberEvents() != rhs.getNumberEvents())
        return false;
    if (this->eventType != rhs.eventType)
        return false;

    // loop over the events
    size_t numEvents = this->getNumberEvents();
    for (size_t i = 0; i < numEvents; ++i) {
    if (!((TofEvent)*(this->events)[i]).equals(*(rhs.events)[i], tolTof, tolPulse))
        return false;
    }                   
}

WeightedEvent EventListTofEvent::getEvent(size_t event_number) {
    return WeightedEvent(events[event_number]);
}

void EventListTofEvent::sortTimeAtSample(const double &tofFactor, const double &tofShift, bool forceResort) const {
    // Check pre-cached sort flag.
    if (this->order == TIMEATSAMPLE_SORT && !forceResort)
        return;

    // Avoid sorting from multiple threads
    std::lock_guard<std::mutex> _lock(m_sortMutex);
    // If the list was sorted while waiting for the lock, return.
    if (this->order == TIMEATSAMPLE_SORT && !forceResort)
        return;

    // Perform sort.

    CompareTimeAtSample<TofEvent> comparitor(tofFactor, tofShift);
    tbb::parallel_sort(events->begin(), events->end(), comparitor);
    // Save the order to avoid unnecessary re-sorting.
    this->order = TIMEATSAMPLE_SORT;    
}

size_t EventListTofEvent::getMemorySize() const {
    return this->events->capacity() * sizeof(TofEvent) + sizeof(EventListTofEvent);
}

void EventListTofEvent::generateHistogramTimeAtSample(const MantidVec &X, MantidVec &Y, MantidVec &E, const double &tofFactor,
                                              const double &tofOffset, bool skipError) {
// All types of weights need to be sorted by time at sample
    this->sortTimeAtSample(tofFactor, tofOffset);
    // Make the single ones
    this->generateCountsHistogramTimeAtSample(X, Y, tofFactor, tofOffset);
    if (!skipError)
      this->generateErrorsHistogram(Y, E);
}

void EventListTofEvent::generateHistogramPulseTime(const MantidVec &X, MantidVec &Y, MantidVec &E, bool skipError) const {
    this->sortPulseTime();
    this->generateCountsHistogramPulseTime(X, Y);
    if (!skipError)
      this->generateErrorsHistogram(Y, E);
}

}