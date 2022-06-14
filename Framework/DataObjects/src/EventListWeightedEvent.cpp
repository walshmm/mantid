#include "MantidDataObjects/EventListWeightedEvent.h"


namespace Mantid::DataObjects {
using Types::Core::DateAndTime;
using Types::Event::TofEvent;
using namespace Mantid::API;

bool EventListWeightedEvent::equals(const EventListBase &rhs, const double tolTof, const double tolWeight,
                       const int64_t tolPulse) const {
    if (this->getNumberEvents() != rhs.getNumberEvents())
        return false;
    if (this->eventType != rhs.eventType)
        return false;

    // loop over the events
    size_t numEvents = this->getNumberEvents();
    for (size_t i = 0; i < numEvents; ++i) {
    if (!((WeightedEvent)this->events[i]).equals(rhs.events[i], tolTof, tolWeight, tolPulse))
        return false;
    }                         
}

WeightedEvent EventListWeightedEvent::getEvent(size_t event_number) {
    return events[event_number];
}

void EventListWeightedEvent::sortTimeAtSample(const double &tofFactor, const double &tofShift, bool forceResort) const {
    // Check pre-cached sort flag.
    if (this->order == TIMEATSAMPLE_SORT && !forceResort)
        return;

    // Avoid sorting from multiple threads
    std::lock_guard<std::mutex> _lock(m_sortMutex);
    // If the list was sorted while waiting for the lock, return.
    if (this->order == TIMEATSAMPLE_SORT && !forceResort)
        return;

    // Perform sort.

    CompareTimeAtSample<WeightedEvent> comparitor(tofFactor, tofShift);
    tbb::parallel_sort(events->begin(), events->end(), comparitor);
    // Save the order to avoid unnecessary re-sorting.
    this->order = TIMEATSAMPLE_SORT;    
}

size_t EventListWeightedEvent::getMemorySize() const {
    return this->events->capacity() * sizeof(WeightedEvent) + sizeof(EventListWeightedEvent);
}

void EventListWeightedEvent::generateHistogramTimeAtSample(const MantidVec &X, MantidVec &Y, MantidVec &E, const double &tofFactor,
                                              const double &tofOffset, bool skipError) {
    throw std::runtime_error("Cannot histogram by time at sample on Weighted "
                             "Events currently"); // This could be supported.
}

void EventListWeightedEvent::generateHistogramPulseTime(const MantidVec &X, MantidVec &Y, MantidVec &E, bool skipError) const {
    throw std::runtime_error("Cannot histogram by pulse time on Weighted "
                             "Events currently"); // This could be supported.
}

void EventListWeightedEvent::generateHistogram(const MantidVec &X, MantidVec &Y, MantidVec &E, bool skipError) const {
    this->sortTof();
    histogramForWeightsHelper(*(this->events.get()), X, Y, E);
}

void EventListWeightedEvent::getWeights(std::vector<double> &weights) const {
    weights.reserve(this->getNumberEvents());
    this->getWeightsHelper(*(this->events.get()), weights);
}

void EventListWeightedEvent::getWeightErrors(std::vector<double> &weightErrors) const {
    // Set the capacity of the vector to avoid multiple resizes
    weightErrors.reserve(this->getNumberEvents());
    this->getWeightErrorsHelper(*(this->events.get()), weightErrors);
}

}