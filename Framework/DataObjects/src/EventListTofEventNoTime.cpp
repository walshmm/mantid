#include "MantidDataObjects/EventListTofEventNoTime.h"


namespace Mantid::DataObjects {
using Types::Core::DateAndTime;
using Types::Event::TofEvent;
using namespace Mantid::API;

bool EventListTofEventNoTime::equals(const EventListBase &rhs, const double tolTof, const double tolWeight,
                       const int64_t tolPulse) const {
    if (this->getNumberEvents() != rhs.getNumberEvents())
        return false;
    if (this->eventType != rhs.eventType)
        return false;

    // loop over the events
    size_t numEvents = this->getNumberEvents();
    for (size_t i = 0; i < numEvents; ++i) {
    if (!((TofEventNoTime)this->events[i]).equals(rhs.events[i], tolTof))
        return false;
    }                         
}

WeightedEvent EventListTofEventNoTime::getEvent(size_t event_number) {
    return WeightedEvent(events[event_number].tof(), 0, events[event_number].weight(),
                         events[event_number].errorSquared());
}

void EventListTofEventNoTime::sortTimeAtSample(const double &tofFactor, const double &tofShift, bool forceResort) const {
    // Check pre-cached sort flag.
    if (this->order == TIMEATSAMPLE_SORT && !forceResort)
        return;

    // Avoid sorting from multiple threads
    std::lock_guard<std::mutex> _lock(m_sortMutex);
    // If the list was sorted while waiting for the lock, return.
    if (this->order == TIMEATSAMPLE_SORT && !forceResort)
        return;

    // Perform sort.

    CompareTimeAtSample<WeightedEventNoTime> comparitor(tofFactor, tofShift);
    tbb::parallel_sort(events.begin(), events.end(), comparitor);
    // Save the order to avoid unnecessary re-sorting.
    this->order = TIMEATSAMPLE_SORT;    
}

void EventListTofEventNoTime::sortPulseTime() const {
    // do nothing, no time to sort;
}

void EventListTofEventNoTime::sortPulseTimeTOF() const {
    // do nothing, no time to sort;
}

void EventListTofEventNoTime::sortPulseTimeTOFDelta(const Types::Core::DateAndTime &start, const double seconds) const override {
    // do nothing, no time to sort;
}

size_t EventListTofEventNoTime::getMemorySize() const {
    return this->events.capacity() * sizeof(WeightedEventNoTime) + sizeof(EventListTofEventNoTime);
}

void EventListTofEventNoTime::compressFatEvents(const double tolerance, const Mantid::Types::Core::DateAndTime &timeStart,
                                  const double seconds, EventListBase *destination) {
    throw std::invalid_argument("Cannot compress events that do not have pulsetime");
}

void EventListTofEventNoTime::generateHistogramTimeAtSample(const MantidVec &X, MantidVec &Y, MantidVec &E, const double &tofFactor,
                                              const double &tofOffset, bool skipError) {
    throw std::runtime_error("TODO: Determine if I can histogram with only pulse time.");
}

void EventListTofEventNoTime::generateHistogramPulseTime(const MantidVec &X, MantidVec &Y, MantidVec &E, bool skipError) const {
    throw std::runtime_error("TODO: Determine if I can histogram with only pulse time.")
}

void EventListTofEventNoTime::filterByPulseTime(DateAndTime start, DateAndTime stop, EventListBase &output) const {
    throw std::runtime_error("EventListBase::filterByTimeAtSample() called on an "
                             "EventListBase that no longer has full time "
                             "information.");
}

void EventListTofEventNoTime::filterByTimeAtSample(Types::Core::DateAndTime start, Types::Core::DateAndTime stop, double tofFactor,
                                     double tofOffset, EventListBase &output) const {
    throw std::runtime_error("EventListBase::filterByTimeAtSample() called on an "
                             "EventListBase that no longer has full time "
                             "information.");
}

void EventListTofEventNoTime::filterInPlace(Kernel::TimeSplitterType &splitter) {
    throw std::runtime_error("EventListBase::filterInPlace() called on an "
                             "EventListBase that no longer has time information.");
}

void EventListTofEventNoTime::splitByTime(Kernel::TimeSplitterType &splitter, std::vector<EventListBase *> outputs) const {
    throw std::runtime_error("EventListBase::splitByTime() called on an EventListBase "
                             "that no longer has time information.");
}

void EventListTofEventNoTime::splitByFullTime(Kernel::TimeSplitterType &splitter, std::map<int, EventListBase *> outputs,
                                bool docorrection, double toffactor, double tofshift) const {
    throw std::runtime_error("EventListBase::splitByTime() called on an EventListBase "
                             "that no longer has time information.");
}

std::string EventListTofEventNoTime::splitByFullTimeMatrixSplitter(const std::vector<int64_t> &vec_splitters_time,
                                                     const std::vector<int> &vecgroups,
                                                     std::map<int, EventListBase *> vec_outputEventList, bool docorrection,
                                                     double toffactor, double tofshift) const {
    throw std::runtime_error("EventListBase::splitByTime() called on an EventListBase "
                             "that no longer has time information.");                                                     
}

void EventListTofEventNoTime::splitByPulseTime(Kernel::TimeSplitterType &splitter, std::map<int, EventListBase *> outputs) const {
    throw std::runtime_error("EventListBase::splitByTime() called on an EventListBase "
                             "that no longer has time information.");
}

void EventListTofEventNoTime::splitByPulseTimeWithMatrix(const std::vector<int64_t> &vec_times, const std::vector<int> &vec_target,
                                           std::map<int, EventListBase *> outputs) const {
    throw std::runtime_error("EventListBase::splitByTime() called on an EventListBase "
                             "that no longer has time information.");
}

}