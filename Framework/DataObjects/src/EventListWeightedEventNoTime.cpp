#include "MantidDataObjects/EventListWeightedEventNoTime.h"


namespace Mantid::DataObjects {
using Types::Core::DateAndTime;
using Types::Event::TofEvent;
using namespace Mantid::API;

EventListWeightedEventNoTime::EventListWeightedEventNoTime() : EventListWeightErrorTofFunctionsTemplate(std::make_shared<std::vector<WeightedEventNoTime>>()){
    events = EventListWeightErrorTofFunctionsTemplate::events;
}

/** Constructor copying from an existing event list
 * @param rhs :: EventListBase object to copy*/
EventListWeightedEventNoTime::EventListBase(const EventListBase &rhs) :EventListWeightedEventNoTime(), IEventList(rhs), m_histogram(rhs.m_histogram), mru{nullptr} {
  // Note that operator= also assigns m_histogram, but the above use of the copy
  // constructor avoid a memory allocation and is thus faster.
  this->operator=(rhs);
}

/** Constructor, taking a vector of events->
 * @param events :: Vector of WeightedEventNoTime's */
EventListWeightedEventNoTime::EventListBase(const std::vector<WeightedEventNoTime> &events)
    :EventListWeightedEventNoTime(), m_histogram(HistogramData::Histogram::XMode::BinEdges, HistogramData::Histogram::YMode::Counts), mru(nullptr) {
  this->events->assign(events.begin(), events.end());
  this->eventType = WEIGHTED_NOTIME;
  this->order = UNSORTED;
}

/// Destructor
EventListWeightedEventNoTime::~EventListBase() {
  // Note: These two lines do not seem to have an effect on releasing memory
  //  at least on Linux. (Memory usage seems to increase event after deleting
  //  EventWorkspaces.
  //  Therefore, for performance, they are kept commented:
  
  
  clear();

  // this->events->clear();
  // std::vector<TofEvent>().swap(events); //Trick to release the vector memory.
}


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
    tbb::parallel_sort(events->begin(), events->end(), comparitor);
    // Save the order to avoid unnecessary re-sorting.
    this->order = TIMEATSAMPLE_SORT;    
}

size_t EventListWeightedEventNoTime::getMemorySize() const {
    return this->events->capacity() * sizeof(TofEventNoTime) + sizeof(EventListWeightedEventNoTime);
}

void EventListWeightedEventNoTime::compressFatEvents(const double tolerance, const Mantid::Types::Core::DateAndTime &timeStart,
                                  const double seconds, EventListBase *destination) {
    throw std::invalid_argument("Cannot compress events that do not have pulsetime");
}

void EventListWeightedEventNoTime::generateHistogramTimeAtSample(const MantidVec &X, MantidVec &Y, MantidVec &E, const double &tofFactor,
                                              const double &tofOffset, bool skipError) {
    throw std::runtime_error("Cannot histogram by time at sample on Weighted Events NoTime");
}

void EventListWeightedEventNoTime::generateHistogramPulseTime(const MantidVec &X, MantidVec &Y, MantidVec &E, bool skipError) const {
    throw std::runtime_error("Cannot histogram by pulse time on Weighted Events NoTime");
}

void EventListWeightedEventNoTime::generateHistogram(const MantidVec &X, MantidVec &Y, MantidVec &E, bool skipError) const {
    this->sortTof();
    histogramForWeightsHelper((*this->events), X, Y, E);
}

void EventListWeightedEventNoTime::addPulsetime(const double seconds) {
    throw std::runtime_error("EventListWeightedEventNoTime::addPulsetime() called on an event "
                             "list with no pulse times. You must call this "
                             "algorithm BEFORE CompressEvents.");
}

void EventListWeightedEventNoTime::addPulsetimes(const std::vector<double> &seconds) {
    throw std::runtime_error("EventListWeightedEventNoTime::addPulsetime() called on an event "
                             "list with no pulse times. You must call this "
                             "algorithm BEFORE CompressEvents.");
}

void EventListWeightedEventNoTime::getWeights(std::vector<double> &weights) const {
    weights.reserve(this->getNumberEvents());
    this->getWeightsHelper((*this->events), weights);
}

void EventListWeightedEventNoTime::getWeightErrors(std::vector<double> &weightErrors) const {
    // Set the capacity of the vector to avoid multiple resizes
    weightErrors.reserve(this->getNumberEvents());
    this->getWeightErrorsHelper((*this->events), weightErrors);
}

void EventListWeightedEventNoTime::filterByPulseTime(Types::Core::DateAndTime start, Types::Core::DateAndTime stop, EventListBase &output) const {
    throw std::runtime_error("EventListWeightedEventNoTime::filterByTimeAtSample() called on an "
                             "EventListBase that no longer has full time "
                             "information.");
}

void EventListWeightedEventNoTime::filterByTimeAtSample(Types::Core::DateAndTime start, Types::Core::DateAndTime stop, double tofFactor,
                                     double tofOffset, EventListBase &output) const {
    throw std::runtime_error("EventListWeightedEventNoTime::filterByTimeAtSample() called on an "
                             "EventListBase that no longer has full time "
                             "information.");
}

void EventListWeightedEventNoTime::filterInPlace(Kernel::TimeSplitterType &splitter) {
    throw std::runtime_error("EventListWeightedEventNoTime::filterInPlace() called on an "
                             "EventListBase that no longer has time information.");
}

void EventListWeightedEventNoTime::splitByTime(Kernel::TimeSplitterType &splitter, std::vector<EventListBase *> outputs) const {
    throw std::runtime_error("EventListWeightedEventNoTime::splitByTime() called on an EventListBase "
                             "that no longer has time information.");
}

void EventListWeightedEventNoTime::splitByFullTime(Kernel::TimeSplitterType &splitter, std::map<int, EventListBase *> outputs,
                                bool docorrection, double toffactor, double tofshift) const {
    throw std::runtime_error("EventListWeightedEventNoTime::splitByTime() called on an EventListBase "
                             "that no longer has time information.");
}

std::string EventListWeightedEventNoTime::splitByFullTimeMatrixSplitter(const std::vector<int64_t> &vec_splitters_time,
                                                     const std::vector<int> &vecgroups,
                                                     std::map<int, EventListBase *> vec_outputEventList, bool docorrection,
                                                     double toffactor, double tofshift) const {
    throw std::runtime_error("EventListWeightedEventNoTime::splitByTime() called on an EventListBase "
                             "that no longer has time information.");                                                     
}

void EventListWeightedEventNoTime::splitByPulseTime(Kernel::TimeSplitterType &splitter, std::map<int, EventListBase *> outputs) const {
    throw std::runtime_error("EventListWeightedEventNoTime::splitByTime() called on an EventListBase "
                             "that no longer has time information.");
}

void EventListWeightedEventNoTime::splitByPulseTimeWithMatrix(const std::vector<int64_t> &vec_times, const std::vector<int> &vec_target,
                                           std::map<int, EventListBase *> outputs) const {
    throw std::runtime_error("EventListWeightedEventNoTime::splitByTime() called on an EventListBase "
                             "that no longer has time information.");
}

}
