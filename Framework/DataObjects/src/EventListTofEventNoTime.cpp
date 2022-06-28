#include "MantidDataObjects/EventListTofEventNoTime.h"


namespace Mantid::DataObjects {
using Types::Core::DateAndTime;
using Types::Event::TofEvent;
using namespace Mantid::API;


EventListTofEventNoTime::EventListTofEventNoTime() {
}

/** Constructor, taking a vector of events.
 * @param events :: Vector of TofEventNoTime's */
EventListTofEventNoTime::EventListTofEventNoTime(const std::vector<TofEventNoTime> &events) : EventListBase() {
     this->events.assign(events.begin(), events.end());
}

/** Constructor copying from an existing event list
 * @param rhs :: EventListBase object to copy*/
EventListTofEventNoTime::EventListTofEventNoTime(const EventList &rhs) :EventListBase() {
  // Note that operator= also assigns m_histogram, but the above use of the copy
  // constructor avoid a memory allocation and is thus faster.
  this->operator=(rhs);
}

/// Destructor
EventListTofEventNoTime::~EventListTofEventNoTime() {
  // Note: These two lines do not seem to have an effect on releasing memory
  //  at least on Linux. (Memory usage seems to increase event after deleting
  //  EventWorkspaces.
  //  Therefore, for performance, they are kept commented:
  
  
  clear();

  // this->events->clear();
  // std::vector<TofEvent>().swap(events); //Trick to release the vector memory.
}


bool EventListTofEventNoTime::equals(const EventList &rhs, const double tolTof, const double tolWeight,
                       const int64_t tolPulse) const {
    if (this->events.size() != rhs.getNumberEvents())
        return false;
    if (this->eventType != rhs.getEventType())
        return false;

    // loop over the events
    size_t numEvents = this->events.size();
    for (size_t i = 0; i < numEvents; ++i) {
    if (!((TofEventNoTime)this->events[i]).equals(rhs.getEventsNoTime()[i], tolTof))
        return false;
    }                         
}

WeightedEvent EventListTofEventNoTime::getEvent(size_t event_number) {
    throw std::runtime_error("EventList: Unimplemented function");
    // return WeightedEvent(events[event_number].tof(), 0, events[event_number].weight(),
    //                      events[event_number].errorSquared());
}

void EventListTofEventNoTime::sortTimeAtSample(const double &tofFactor, const double &tofShift, bool forceResort) const {
    throw std::invalid_argument("Cannot sort events that do not have pulsetime");   
}

void EventListTofEventNoTime::sortPulseTime() const {
    // do nothing, no time to sort;
}

void EventListTofEventNoTime::sortPulseTimeTOF() const {
    // do nothing, no time to sort;
}

void EventListTofEventNoTime::sortPulseTimeTOFDelta(const Types::Core::DateAndTime &start, const double seconds) const {
    // do nothing, no time to sort;
}

size_t EventListTofEventNoTime::getMemorySize() const {
    return this->events.capacity() * sizeof(WeightedEventNoTime) + sizeof(EventListTofEventNoTime);
}

void EventListTofEventNoTime::compressFatEvents(const double tolerance, const Mantid::Types::Core::DateAndTime &timeStart,
                                  const double seconds, EventList *destination) {
    throw std::invalid_argument("Cannot compress events that do not have pulsetime");
}

void EventListTofEventNoTime::generateHistogramTimeAtSample(const MantidVec &X, MantidVec &Y, MantidVec &E, const double &tofFactor,
                                              const double &tofOffset, bool skipError) {
    throw std::runtime_error("TODO: Determine if I can histogram with only pulse time.");
}

void EventListTofEventNoTime::generateHistogramPulseTime(const MantidVec &X, MantidVec &Y, MantidVec &E, bool skipError) const {
    throw std::runtime_error("TODO: Determine if I can histogram with only pulse time.");
}

void EventListTofEventNoTime::filterByPulseTime(Types::Core::DateAndTime start, Types::Core::DateAndTime stop, EventList &output) const {
    throw std::runtime_error("EventListTofEventNoTime::filterByTimeAtSample() called on an "
                             "EventListBase that no longer has full time "
                             "information.");
}

void EventListTofEventNoTime::filterByTimeAtSample(Types::Core::DateAndTime start, Types::Core::DateAndTime stop, double tofFactor,
                                     double tofOffset, EventList &output) const {
    throw std::runtime_error("EventListTofEventNoTime::filterByTimeAtSample() called on an "
                             "EventListBase that no longer has full time "
                             "information.");
}

void EventListTofEventNoTime::filterInPlace(Kernel::TimeSplitterType &splitter) {
    throw std::runtime_error("EventListTofEventNoTime::filterInPlace() called on an "
                             "EventListBase that no longer has time information.");
}

void EventListTofEventNoTime::splitByTime(Kernel::TimeSplitterType &splitter, std::vector<EventList *> outputs) const {
    throw std::runtime_error("EventListTofEventNoTime::splitByTime() called on an EventListBase "
                             "that no longer has time information.");
}

void EventListTofEventNoTime::splitByFullTime(Kernel::TimeSplitterType &splitter, std::map<int, EventList *> outputs,
                                bool docorrection, double toffactor, double tofshift) const {
    throw std::runtime_error("EventListTofEventNoTime::splitByTime() called on an EventListBase "
                             "that no longer has time information.");
}

std::string EventListTofEventNoTime::splitByFullTimeMatrixSplitter(const std::vector<int64_t> &vec_splitters_time,
                                                     const std::vector<int> &vecgroups,
                                                     std::map<int, EventList *> vec_outputEventList, bool docorrection,
                                                     double toffactor, double tofshift) const {
    throw std::runtime_error("EventListTofEventNoTime::splitByTime() called on an EventListBase "
                             "that no longer has time information.");                                                     
}

void EventListTofEventNoTime::splitByPulseTime(Kernel::TimeSplitterType &splitter, std::map<int, EventList *> outputs) const {
    throw std::runtime_error("EventListTofEventNoTime::splitByTime() called on an EventListBase "
                             "that no longer has time information.");
}

void EventListTofEventNoTime::splitByPulseTimeWithMatrix(const std::vector<int64_t> &vec_times, const std::vector<int> &vec_target,
                                           std::map<int, EventList *> outputs) const {
    throw std::runtime_error("EventListTofEventNoTime::splitByTime() called on an EventListBase "
                             "that no longer has time information.");
}

}