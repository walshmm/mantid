// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#include "MantidDataObjects/EventList.h"
#include "MantidDataObjects/EventListTofEvent.h"
#include "MantidDataObjects/EventListWeightedEvent.h"
#include "MantidDataObjects/EventListWeightedEventNoTime.h"
#include "MantidDataObjects/EventListTofEventNoTime.h"
#include "MantidAPI/MatrixWorkspace.h"
#include "MantidDataObjects/EventWorkspaceMRU.h"
#include "MantidDataObjects/Histogram1D.h"
#include "MantidKernel/DateAndTime.h"
#include "MantidKernel/DateAndTimeHelpers.h"
#include "MantidKernel/Exception.h"
#include "MantidKernel/Logger.h"
#include "MantidKernel/Unit.h"

#ifdef _MSC_VER
// qualifier applied to function type has no meaning; ignored
#pragma warning(disable : 4180)
#endif
#include "tbb/parallel_sort.h"
#ifdef _MSC_VER
#pragma warning(default : 4180)
#endif

#include <cfloat>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>

using std::ostream;
using std::runtime_error;
using std::size_t;
using std::vector;

namespace Mantid::DataObjects {
using Types::Core::DateAndTime;
using Types::Event::TofEvent;
using namespace Mantid::API;

/// Constructor (empty)
// EventWorkspace is always histogram data and so is thus EventList
EventList::EventList() {
  this->eventList = EventListTofEvent();
}

/** Constructor with a MRU list
 * @param mru :: pointer to the MRU of the parent EventWorkspace
 * @param specNo :: the spectrum number for the event list
 */
EventList::EventList(EventWorkspaceMRU *mru, specnum_t specNo) {
  this->eventList = EventListTofEvent(mru, specNo);
}

/** Constructor copying from an existing event list
 * @param rhs :: EventList object to copy*/
EventList::EventList(const EventList &rhs) : IEventList(rhs), m_histogram(rhs.m_histogram), mru{nullptr} {
  // Note that operator= also assigns m_histogram, but the above use of the copy
  // constructor avoid a memory allocation and is thus faster.
  this->operator=(rhs);
}

/** Constructor, taking a vector of events.
 * @param events :: Vector of TofEvent's */
EventList::EventList(const std::vector<TofEvent> &events){
  this->eventList = EventListTofEvent(events);
}

/** Constructor, taking a vector of events.
 * @param events :: Vector of WeightedEvent's */
EventList::EventList(const std::vector<WeightedEvent> &events){
  this->eventList = EventListWeightedEvent(events);
}

/** Constructor, taking a vector of events.
 * @param events :: Vector of WeightedEventNoTime's */
EventList::EventList(const std::vector<WeightedEventNoTime> &events) {
  this->eventList = EventListWeightedEventNoTime(events);
}

/** Constructor, taking a vector of events.
 * @param events :: Vector of TofEventNoTime's */
EventList::EventList(const std::vector<TofEventNoTime> &events) {
  this->eventList = EventListTofEventNoTime(events);
}

/// Destructor
EventList::~EventList() {
  // Note: These two lines do not seem to have an effect on releasing memory
  //  at least on Linux. (Memory usage seems to increase event after deleting
  //  EventWorkspaces.
  //  Therefore, for performance, they are kept commented:
  clear();

  // this->events.clear();
  // std::vector<TofEvent>().swap(events); //Trick to release the vector memory.
}

/// Copy data from another EventList, via ISpectrum reference.
void EventList::copyDataFrom(const ISpectrum &source) { this->eventList->copyDataInto(source); }

/// Used by copyDataFrom for dynamic dispatch for its `source`.
void EventList::copyDataInto(EventList &sink) const {
  this->eventList->copyDataInto(sink);
}

/// Used by Histogram1D::copyDataFrom for dynamic dispatch for `other`.
void EventList::copyDataInto(Histogram1D &sink) const { this->eventList->copyDataInto(sink); }

// --------------------------------------------------------------------------
/** Create an EventList from a histogram. This converts bins to weighted
 * events.
 * Any existing events are cleared.
 *
 * @param inSpec :: ISpectrum ptr to histogram data.
 * @param GenerateZeros :: if true, generate event(s) for empty bins
 * @param GenerateMultipleEvents :: if true, create several evenly-spaced fake
 *events inside the bin
 * @param MaxEventsPerBin :: max number of events to generate in one bin, if
 *GenerateMultipleEvents
 */
void EventList::createFromHistogram(const ISpectrum *inSpec, bool GenerateZeros, bool GenerateMultipleEvents,
                                    int MaxEventsPerBin) {
  this->eventList->createFromHistogram(inSpec, GenerateZeros, GenerateMultipleEvents, MaxEventsPerBin);
}

// --------------------------------------------------------------------------
// --- Operators
// -------------------------------------------------------------------

/** Copy into this event list from another
 * @param rhs :: We will copy all the events from that into this object.
 * @return reference to this
 * */
EventList &EventList::operator=(const EventList &rhs) {
  return this->eventList = rhs;
}

// --------------------------------------------------------------------------
/** Append an event to the histogram.
 * @param event :: TofEvent to add at the end of the list.
 * @return reference to this
 * */
EventList &EventList::operator+=(const TofEvent &event) {
  this->eventList+=event;
  return *this;
}

// --------------------------------------------------------------------------
/** Append a list of events to the histogram.
 * The internal event list will switch to the required type.
 *
 * @param more_events :: A vector of events to append.
 * @return reference to this
 * */
EventList &EventList::operator+=(const std::vector<TofEvent> &more_events) {
  this->eventList+=more_events;
  return *this;
}

// --------------------------------------------------------------------------
/** Append a WeightedEvent to the histogram.
 * Note: The whole list will switch to weights (a possibly lengthy operation)
 *  if it did not have weights before.
 *
 * @param event :: WeightedEvent to add at the end of the list.
 * @return reference to this
 * */
EventList &EventList::operator+=(const WeightedEvent &event) {
  this->eventList+=event;
  return *this;
}

// --------------------------------------------------------------------------
/** Append a list of events to the histogram.
 * Note: The whole list will switch to weights (a possibly lengthy operation)
 *  if it did not have weights before.
 *
 * @param more_events :: A vector of events to append.
 * @return reference to this
 * */
EventList &EventList::operator+=(const std::vector<WeightedEvent> &more_events) {
  this->eventList=more_events;
  return *this;
}

// --------------------------------------------------------------------------
/** Append a list of events to the histogram.
 * Note: The whole list will switch to weights (a possibly lengthy operation)
 *  if it did not have weights before.
 *
 * @param more_events :: A vector of events to append.
 * @return reference to this
 * */
EventList &EventList::operator+=(const std::vector<WeightedEventNoTime> &more_events) {
  this->eventList+= more_events;
  return *this;
}

// --------------------------------------------------------------------------
/** Append another EventList to this event list.
 * The event lists are concatenated, and a union of the sets of detector ID's is
 *done.
 * Switching of event types may occur if the two are different.
 *
 * @param more_events :: Another EventList.
 * @return reference to this
 * */
EventList &EventList::operator+=(const EventList &more_events) {
  this->eventList+=more_events;
  return *this;
}

// --------------------------------------------------------------------------
/** SUBTRACT another EventList from this event list.
 * The event lists are concatenated, but the weights of the incoming
 *    list are multiplied by -1.0.
 *
 * @param more_events :: Another EventList.
 * @return reference to this
 * */
EventList &EventList::operator-=(const EventList &more_events) {
  this->eventList-=more_events;
  // NOTE: What to do about detector ID's?
  return *this;
}

// --------------------------------------------------------------------------
/** Equality operator between EventList's
 * @param rhs :: other EventList to compare
 * @return :: true if equal.
 */
bool EventList::operator==(const EventList &rhs) const {
  return this->eventList==rhs;
}

/** Inequality comparator
 * @param rhs :: other EventList to compare
 * @return :: true if not equal.
 */
bool EventList::operator!=(const EventList &rhs) const { return this->eventList!=rhs; }

bool EventList::equals(const EventList &rhs, const double tolTof, const double tolWeight,
                       const int64_t tolPulse) const {
  return this->eventList->equals(rhs, tolTof, tolWeight, tolPulse);
}

// -----------------------------------------------------------------------------------------------
/** Return the type of Event vector contained within.
 * @return :: a EventType value.
 */
EventType EventList::getEventType() const { return this->eventList->getEventType(); }

// -----------------------------------------------------------------------------------------------
/** Switch the EventList to use the given EventType (TOF, WEIGHTED, or
 * WEIGHTED_NOTIME)
 */
void EventList::switchTo(EventType newType) {
  //TODO: may be a nessessary evil, will think of a way to remove this another time
  switch (newType) {
  case TOF:
    if (eventType != TOF)
      throw std::runtime_error("EventListBase::switchTo() called on an EventListBase "
                               "with weights to go down to TofEvent's. This "
                               "would remove weight information and therefore "
                               "is not possible.");
    break;

  case WEIGHTED:
    switchToWeightedEvents();
    break;

  case WEIGHTED_NOTIME:
    switchToWeightedEventsNoTime();
    break;

  case UNWEIGHTED:
    switchToUnweightedEvents();
    break;
  }
}


// ==============================================================================================
// --- Testing functions (mostly)
// ---------------------------------------------------------------
// ==============================================================================================

/** Return the given event in the list.
 * Handles the different types of events by converting to WeightedEvent (the
 * most general type).
 * @param event_number :: the index of the event to retrieve
 * @return a WeightedEvent
 */
WeightedEvent EventList::getEvent(size_t event_number) {
  return this->eventList->getEvent(event_number);
}

// ==============================================================================================
// --- Handling the event list
// -------------------------------------------------------------------
// ==============================================================================================

/** Return the const list of TofEvents contained.
 * NOTE! This should be used for testing purposes only, as much as possible. The
 *EventList
 * may contain weighted events, requiring use of getWeightedEvents() instead.
 *
 * @return a const reference to the list of non-weighted events
 * */
const std::vector<TofEvent> &EventList::getEvents() const {
  return this->eventList->getEvents();
}

/** Return the list of TofEvents contained.
 * NOTE! This should be used for testing purposes only, as much as possible. The
 *EventList
 * may contain weighted events, requiring use of getWeightedEvents() instead.
 *
 * @return a reference to the list of non-weighted events
 * */
std::vector<TofEvent> &EventList::getEvents() {
  return this->eventList->getEvents();
}

// std::vector<Event> &EventList::getEventsTyped() {
//   //TODO probably throw, else possible inf loop?;
//   return this->eventList->getEventsTyped();
// }

// const std::vector<Event> &EventList::getEventsTyped() const {
//   return this->eventList->getEventsTyped();
// }

/** Return the list of WeightedEvent contained.
 * NOTE! This should be used for testing purposes only, as much as possible. The
 *EventList
 * may contain un-weighted events, requiring use of getEvents() instead.
 *
 * @return a reference to the list of weighted events
 * */
std::vector<WeightedEvent> &EventList::getWeightedEvents() {
  return this->eventList->getWeightedEvents();
}

/** Return the list of WeightedEvent contained.
 * NOTE! This should be used for testing purposes only, as much as possible. The
 *EventList
 * may contain un-weighted events, requiring use of getEvents() instead.
 *
 * @return a const reference to the list of weighted events
 * */
const std::vector<WeightedEvent> &EventList::getWeightedEvents() const {
  return this->eventList->getWeightedEvents();
}

/** Return the list of WeightedEvent contained.
 * NOTE! This should be used for testing purposes only, as much as possible.
 *
 * @return a reference to the list of weighted events
 * */
std::vector<WeightedEventNoTime> &EventList::getWeightedEventsNoTime() {
  return this->eventList->getWeightedEventsNoTime();
}

// /** Return the list of WeightedEvent contained.
//  * NOTE! This should be used for testing purposes only, as much as possible.
//  *
//  * @return a reference to the list of weighted events
//  * */
// std::vector<WeightedEventNoTime> &EventList::getUnweightedEvents() {
//   return this->eventList->getUnweightedEvents();
// }

/** Return the list of WeightedEventNoTime contained.
 * NOTE! This should be used for testing purposes only, as much as possible.
 *
 * @return a const reference to the list of weighted events
 * */
const std::vector<WeightedEventNoTime> &EventList::getWeightedEventsNoTime() const {
  return this->eventList->getWeightedEventsNoTime();
}

/** Clear the list of events and any
 * associated detector ID's.
 * */
void EventList::clear(const bool removeDetIDs) {
  this->eventList->clear(removeDetIDs);
}

/** Clear any unused event lists (the ones that do not
 * match the currently used type).
 * Memory is freed.
 * */
void EventList::clearUnused() {
  this->eventList->clearUnused();
}

/// Mask the spectrum to this value. Removes all events.
void EventList::clearData() { this->eventList->clearData(); }

/** Sets the MRU list for this event list
 *
 * @param newMRU :: new MRU for the workspace containing this EventList
 */
void EventList::setMRU(EventWorkspaceMRU *newMRU) { this->eventList->setMRU(newMRU);}

/** Reserve a certain number of entries in event list of the specified eventType
 *
 * Calls std::vector<>::reserve() in order to pre-allocate the length of the
 *event list vector.
 *
 * @param num :: number of events that will be in this EventList
 */
void EventList::reserve(size_t num) {
  this->eventList->reserve(num);
}

// ==============================================================================================
// --- Sorting functions -----------------------------------------------------
// ==============================================================================================

// --------------------------------------------------------------------------
/** Sort events by TOF or Frame
 * @param order :: Order by which to sort.
 * */
void EventList::sort(const EventSortType order) const {
  this->eventList->sort(order);
}

// --------------------------------------------------------------------------
/** Manually set the event list sort order value. No actual sorting takes place.
 * SHOULD ONLY BE USED IN TESTS or if you know what you are doing.
 * @param order :: sort order to set.
 */
void EventList::setSortOrder(const EventSortType order) const { this->eventList->setSortOrder(order); }

// --------------------------------------------------------------------------
/** Sort events by TOF in one thread */
void EventList::sortTof() const {
  this->eventList->sortTof();
}

// --------------------------------------------------------------------------
/**
 * Sort events by time at sample
 * @param tofFactor : For the elastic case, L1 / (L1 + L2)
 * @param tofShift : Tof offset in Seconds
 * @param forceResort : If the tofFactor, or tofShift are different
 *   from a previous run of the same sort type, you need to trigger a full
 * resort using forceResort = true. False by default.
 */
void EventList::sortTimeAtSample(const double &tofFactor, const double &tofShift, bool forceResort) const {
  this->eventList->sortTimeAtSample(tofFactor, tofShift, forceResort);
}

// --------------------------------------------------------------------------
/** Sort events by Frame */
void EventList::sortPulseTime() const {
  this->eventList->sortPulseTime();
}

/*
 * Sort events by pulse time + TOF
 * (the absolute time)
 */
void EventList::sortPulseTimeTOF() const {
  this->eventList->sortPulseTimeTOF();
}

// --------------------------------------------------------------------------
/** Return true if the event list is sorted by TOF */
bool EventList::isSortedByTof() const { return this->eventList->isSortedByTof(); }

// --------------------------------------------------------------------------
/** Return the type of sorting used in this event list */
EventSortType EventList::getSortType() const { return this->eventList->getSortType(); }

// --------------------------------------------------------------------------
/** Reverse the histogram boundaries and the associated events if they are
 * sorted
 * by time-of-flight.
 * Does nothing if sorted otherwise or unsorted.
 * */
void EventList::reverse() {
  this->eventList->reverse();
}

// --------------------------------------------------------------------------
/** Return the number of events in the list.
 * NOTE: If the events have weights, this returns the NUMBER of WeightedEvent's
 *in the
 * list, and NOT the sum of their weights (which may be two different numbers).
 *
 * @return the number of events in the list.
 *  */
size_t EventList::getNumberEvents() const {
  return this->eventList->getNumberEvents();
}

/**
 * Much like stl containers, returns true if there is nothing in the event list.
 */
bool EventList::empty() const {
  return this->eventList->empty();
}

// --------------------------------------------------------------------------
/** Memory used by this event list. Note: It reports the CAPACITY of the
 * vectors, rather than their size, since that is a more accurate
 * representation of the size used.
 *
 * @return :: the memory used by the EventList, in bytes.
 * */
size_t EventList::getMemorySize() const {
  return this->eventList->getMemorySize();
}

// --------------------------------------------------------------------------
/** Return the size of the histogram data.
 * @return the size of the histogram representation of the data (size of Y) **/
size_t EventList::histogram_size() const {
  return this->eventList->histogram_size();
}

// ==============================================================================================
// --- Setting the Histogram X axis, without recalculating the histogram
// -----------------------
// ==============================================================================================

/** Deprecated, use setSharedX() instead. Set the x-component for the histogram
 * view. This will NOT cause the histogram to be calculated.
 * @param X :: The vector of doubles to set as the histogram limits.
 */
void EventList::setX(const Kernel::cow_ptr<HistogramData::HistogramX> &X) {
  this->eventList->setX(X);
}

/** Deprecated, use mutableX() instead. Returns a reference to the x data.
 *  @return a reference to the X (bin) vector.
 */
MantidVec &EventList::dataX() {
  return this->eventList->dataX();
}

/** Deprecated, use x() instead. Returns a const reference to the x data.
 *  @return a reference to the X (bin) vector. */
const MantidVec &EventList::dataX() const { return this->eventList->dataX(); }

/// Deprecated, use x() instead. Returns the x data const
const MantidVec &EventList::readX() const { return this->eventList->readX(); }

/// Deprecated, use sharedX() instead. Returns a pointer to the x data
Kernel::cow_ptr<HistogramData::HistogramX> EventList::ptrX() const { return this->eventList->ptrX(); }

/// Deprecated, use mutableDx() instead.
MantidVec &EventList::dataDx() { return this->eventList->dataDx(); }
/// Deprecated, use dx() instead.
const MantidVec &EventList::dataDx() const { return this->eventList->dataDx(); }
/// Deprecated, use dx() instead.
const MantidVec &EventList::readDx() const { return this->eventList->readDx(); }

// ==============================================================================================
// --- Return Data Vectors --------------------------------------------------
// ==============================================================================================

/** Calculates and returns a pointer to the Y histogrammed data.
 * Remember to delete your pointer after use!
 *
 * @return a pointer to a MantidVec
 */
MantidVec *EventList::makeDataY() const {
  return this->eventList->makeDataY();
}

/** Calculates and returns a pointer to the E histogrammed data.
 * Remember to delete your pointer after use!
 *
 * @return a pointer to a MantidVec
 */
MantidVec *EventList::makeDataE() const {
  return this->eventList->makeDataE();
}

HistogramData::Histogram EventList::histogram() const {
  return this->eventList->histogram();
}

HistogramData::Counts EventList::counts() const { return this->eventList->counts(); }

HistogramData::CountVariances EventList::countVariances() const { return this->eventList->countVariances(); }

HistogramData::CountStandardDeviations EventList::countStandardDeviations() const {
  return this->eventList->countStandardDeviations();
}

HistogramData::Frequencies EventList::frequencies() const { return this->eventList->frequencies(); }

HistogramData::FrequencyVariances EventList::frequencyVariances() const { return this->eventList->frequencyVariances(); }

HistogramData::FrequencyStandardDeviations EventList::frequencyStandardDeviations() const {
  return this->eventList->frequencyStandardDeviations();
}

const HistogramData::HistogramY &EventList::y() const {
  return this->eventList->y();
}
const HistogramData::HistogramE &EventList::e() const {
  return this->eventList->e();
}
Kernel::cow_ptr<HistogramData::HistogramY> EventList::sharedY() const {
  return this->eventList->sharedY();
}
Kernel::cow_ptr<HistogramData::HistogramE> EventList::sharedE() const {
  return this->eventList->sharedE();
}
/** Look in the MRU to see if the Y histogram has been generated before.
 * If so, return that. If not, calculate, cache and return it.
 *
 * @return reference to the Y vector.
 */
const MantidVec &EventList::dataY() const {
  this->eventList->dataY();
}

/** Look in the MRU to see if the E histogram has been generated before.
 * If so, return that. If not, calculate, cache and return it.
 *
 * @return reference to the E vector.
 */
const MantidVec &EventList::dataE() const {
  return this->eventList->dataE();
}
// --------------------------------------------------------------------------
/** Compress the event list by grouping events with the same
 * TOF (within a given tolerance). PulseTime is ignored.
 * The event list will be switched to WeightedEventNoTime.
 *
 * @param tolerance :: how close do two event's TOF have to be to be considered
 *the same.
 * @param destination :: EventList that will receive the compressed events. Can
 *be == this.
 */
void EventList::compressEvents(double tolerance, EventList *destination) {
  this->eventList->compressEvents(tolerance, destination);
}

void EventList::compressFatEvents(const double tolerance, const Mantid::Types::Core::DateAndTime &timeStart,
                                  const double seconds, EventList *destination) {
  this->eventList->compressFatEvents(tolerance, timeStart, seconds, destination);
}

// --------------------------------------------------------------------------
/** Utility function:
 * Returns the iterator into events of the first TofEvent with
 * tof() > seek_tof
 * Will return events.end() if nothing is found!
 *
 * @param events :: event vector in which to look.
 * @param seek_tof :: tof to find (typically the first bin X[0])
 * @return iterator where the first event matching it is.
 */
template <class T>
typename std::vector<T>::const_iterator static findFirstEvent(const std::vector<T> &events, T seek_tof) {
  return std::find_if_not(events.cbegin(), events.cend(), [seek_tof](const T &x) { return x < seek_tof; });
}

// --------------------------------------------------------------------------
/** Utility function:
 * Returns the iterator into events of the first TofEvent with
 * tof() > seek_tof
 * Will return events.end() if nothing is found!
 *
 * @param events :: event vector in which to look.
 * @param seek_tof :: tof to find (typically the first bin X[0])
 * @return iterator where the first event matching it is.
 */
template <class T> typename std::vector<T>::iterator static findFirstEvent(std::vector<T> &events, T seek_tof) {
  return std::find_if_not(events.begin(), events.end(), [seek_tof](const T &x) { return x < seek_tof; });
}

// --------------------------------------------------------------------------
/** Generates both the Y and E (error) histograms
 * for an EventListBase with Weightedevents.
 *
 * @param events: vector of events (with weights)
 * @param X: X-bins supplied
 * @param Y: counts returned
 * @param E: errors returned
 * @throw runtime_error if the EventListBase does not have weighted events
 */
template <class T>
void EventListBase::histogramForWeightsHelper(const std::vector<T> &events, const MantidVec &X, MantidVec &Y,
                                          MantidVec &E) {
  // For slight speed=up.
  size_t x_size = X.size();

  if (x_size <= 1) {
    // X was not set. Return an empty array.
    Y.resize(0, 0);
    return;
  }

  // If the sizes are the same, then the "resize" command will NOT clear the
  // original values.
  bool mustFill = (Y.size() == x_size - 1);
  // Clear the Y data, assign all to 0.
  Y.resize(x_size - 1, 0.0);
  // Clear the Error data, assign all to 0.
  // Note: Errors will be squared until the last step.
  E.resize(x_size - 1, 0.0);

  if (mustFill) {
    // We must make sure the starting point is 0.0
    std::fill(Y.begin(), Y.end(), 0.0);
    std::fill(E.begin(), E.end(), 0.0);
  }

  //---------------------- Histogram without weights
  //---------------------------------

  // Do we even have any events to do?
  if (!events.empty()) {
    // Iterate through all events (sorted by tof)
    auto itev = findFirstEvent(events, T(X[0]));
    auto itev_end = events.cend();
    // The above can still take you to end() if no events above X[0], so check
    // again.
    if (itev == itev_end)
      return;

    // Find the first bin
    size_t bin = 0;
    // The tof is greater the first bin boundary, so we need to find the first
    // bin
    double tof = itev->tof();
    while (bin < x_size - 1) {
      // Within range?
      if ((tof >= X[bin]) && (tof < X[bin + 1])) {
        // Add up the weight (convert to double before adding, to preserve
        // precision)
        Y[bin] += double(itev->m_weight);
        E[bin] += double(itev->m_errorSquared); // square of error
        break;
      }
      ++bin;
    }
    // Go to the next event, we've already binned this first one.
    ++itev;

    // Keep going through all the events
    while ((itev != itev_end) && (bin < x_size - 1)) {
      tof = itev->tof();
      while (bin < x_size - 1) {
        // Within range? Since both events and X are sorted, they are going to
        // have
        // tof >= X[bin] because the previous event was.
        if (tof < X[bin + 1]) {
          // Add up the weight (convert to double before adding, to preserve
          // precision)
          Y[bin] += double(itev->m_weight);
          E[bin] += double(itev->m_errorSquared); // square of error
          break;
        }
        ++bin;
      }
      ++itev;
    }
  } // end if (there are any events to histogram)

  // Now do the sqrt of all errors
  std::transform(E.begin(), E.end(), E.begin(), static_cast<double (*)(double)>(sqrt));
}
// --------------------------------------------------------------------------
/** Generates both the Y and E (error) histograms w.r.t Pulse Time
 * for an EventList with or without WeightedEvents.
 *
 * @param X: x-bins supplied
 * @param Y: counts returned
 * @param E: errors returned
 * @param skipError: skip calculating the error. This has no effect for weighted
 *        events; you can just ignore the returned E vector.
 */
void EventList::generateHistogramPulseTime(const MantidVec &X, MantidVec &Y, MantidVec &E, bool skipError) const {
  this->eventList->generateHistogramPulseTime(X, Y, E, skipError);
}

/**
 * Generates both the Y and E (error) histograms w.r.t Time at sample position.
 * @param X: X - axis supplied as reference
 * @param Y: counts to fill
 * @param E: errors to fill
 * @param tofFactor : Time of flight factor. Usually L1/(L1 + L2)
 * @param tofOffset : Time of flight offset.
 * @param skipError : skip calculating the error. This has no effect for
 * weighted
 *          events; you can just ignore the returned E vector.
 */
void EventList::generateHistogramTimeAtSample(const MantidVec &X, MantidVec &Y, MantidVec &E, const double &tofFactor,
                                              const double &tofOffset, bool skipError) const {
  this->eventList->generateHistogramTimeAtSample(X, Y, E, tofFactor, tofOffset, skipError);
}

// --------------------------------------------------------------------------
/** Generates both the Y and E (error) histograms w.r.t TOF
 * for an EventList with or without WeightedEvents.
 *
 * @param X: x-bins supplied
 * @param Y: counts returned
 * @param E: errors returned
 * @param skipError: skip calculating the error. This has no effect for weighted
 *        events; you can just ignore the returned E vector.
 */
void EventList::generateHistogram(const MantidVec &X, MantidVec &Y, MantidVec &E, bool skipError) const {
  this->eventList->generateHistogram(X, Y, E, skipError);
}


/** With respect to PulseTime fill a histogram given equal histogram
 *   bins.
 * Number of bins is equal to number of elements in vector Y.
 * Appends values to existing Y values.
 *
 * @param xMin :: Minimal Pulse time (in nanoseconds,
 *                i.e. DateTime->totalNanoseconds()) value to include
 *                in binning.
 * @param xMax :: Maximal Pulse time value to constrain binning by (include the
 *                times smaller than right boundary, excluding equal)
 * @param Y :: The generated counts histogram
 * @param TOF_min -- min TOF to include in histogram.
 * @param TOF_max -- max TOF to constrain values included in histogram.
 */
void EventList::generateCountsHistogramPulseTime(const double &xMin, const double &xMax, MantidVec &Y,
                                                 const double TOF_min, const double TOF_max) const {

  this->eventList->generateCountsHistogramPulseTime(xMin, xMax, Y, TOF_min, TOF_max);
}


// --------------------------------------------------------------------------
/** Integrate the events between a range of X values, or all events.
 *
 * @param minX :: minimum X bin to use in integrating.
 * @param maxX :: maximum X bin to use in integrating.
 * @param entireRange :: set to true to use the entire range. minX and maxX are
 *then ignored!
 * @return the integrated number of events.
 */
double EventList::integrate(const double minX, const double maxX, const bool entireRange) const {
  return this->eventList->integrate(minX, maxX, entireRange);
}

/** Integrate the events between a range of X values, or all events.
 *
 * @param minX :: minimum X bin to use in integrating.
 * @param maxX :: maximum X bin to use in integrating.
 * @param entireRange :: set to true to use the entire range. minX and maxX are
 *then ignored!
 * @param sum :: place holder for the resulting sum
 * @param error :: place holder for the resulting sum of errors
 * @return the integrated number of events.
 */
void EventList::integrate(const double minX, const double maxX, const bool entireRange, double &sum,
                          double &error) const {
  this->eventList->integrate(minX, maxX, entireRange, sum, error);
}

// ==============================================================================================
// ----------- Conversion Functions (changing tof values)
// ---------------------------------------
// ==============================================================================================

/**
 * @param func Function to do the conversion.
 * @param sorting How the events are sorted after the operation. 0 = unsorted
 * (default),
 * positive = unchanged, negative = reverse.
 */
void EventList::convertTof(std::function<double(double)> func, const int sorting) {
  this->eventList->convertTof(func, sorting);
}

// --------------------------------------------------------------------------
/**
 * Convert the time of flight by tof'=tof*factor+offset
 * @param factor :: The value to scale the time-of-flight by
 * @param offset :: The value to shift the time-of-flight by
 */
void EventList::convertTof(const double factor, const double offset) {
  this->eventList->convertTof(factor, offset);
}


// --------------------------------------------------------------------------
/**
 * Convert the units in the TofEvent's m_tof field to
 *  some other value, by scaling by a multiplier.
 * @param factor :: conversion factor (e.g. multiply TOF by this to get
 * d-spacing)
 */
void EventList::scaleTof(const double factor) { this->eventList->scaleTof(factor); }

// --------------------------------------------------------------------------
/** Add an offset to the TOF of each event in the list.
 *
 * @param offset :: The value to shift the time-of-flight by
 */
void EventList::addTof(const double offset) { this->eventList->addTof(offset); }

// --------------------------------------------------------------------------
/** Add an offset to the pulsetime (wall-clock time) of each event in the list.
 *
 * @param seconds :: The value to shift the pulsetime by, in seconds
 */
void EventList::addPulsetime(const double seconds) {
  this->eventList->addPulsetime(seconds);
}

// --------------------------------------------------------------------------
/** Add an offset to the pulsetime (wall-clock time) of each event in the list.
 *
 * @param seconds :: A set of values to shift the pulsetime by, in seconds
 */
void EventList::addPulsetimes(const std::vector<double> &seconds) {
  this->eventList->addPulsetimes(seconds);
}

// --------------------------------------------------------------------------
/**
 * Mask out events that have a tof between tofMin and tofMax (inclusively).
 * Events are removed from the list.
 * @param tofMin :: lower bound of TOF to filter out
 * @param tofMax :: upper bound of TOF to filter out
 */
void EventList::maskTof(const double tofMin, const double tofMax) {
  this->eventList->maskTof(tofMin, tofMax);
}
// --------------------------------------------------------------------------
/**
 * Mask out events by the condition vector.
 * Events are removed from the list.
 * @param mask :: condition vector
 */
void EventList::maskCondition(const std::vector<bool> &mask) {
  this->eventList->maskCondition(mask);
}

/** Fill a vector with the list of TOFs
 *  @param tofs :: A reference to the vector to be filled
 */
void EventList::getTofs(std::vector<double> &tofs) const {
  this->eventList->getTofs(tofs);
}

/** Get the times-of-flight of each event in this EventList.
 *
 * @return by copy a vector of doubles of the tof() value
 */
std::vector<double> EventList::getTofs() const {
  return this->eventList->getTofs();
}

/** Fill a vector with the list of Weights
 *  @param weights :: A reference to the vector to be filled
 */
void EventList::getWeights(std::vector<double> &weights) const {
  this->eventList->getWeights(weights);
}

/** Get the weight of each event in this EventList.
 *
 * @return by copy a vector of doubles of the weight() value
 */
std::vector<double> EventList::getWeights() const {
  return this->eventList->getWeights();
}

/** Fill a vector with the list of Weight Errors
 *  @param weightErrors :: A reference to the vector to be filled
 */
void EventList::getWeightErrors(std::vector<double> &weightErrors) const {
  this->eventList->getWeightErrors(weightErrors);
}

/** Get the weight error of each event in this EventList.
 *
 * @return by copy a vector of doubles of the weight() value
 */
std::vector<double> EventList::getWeightErrors() const {
  return this->eventList->getWeightErrors();
}

/** Get the pulse times of each event in this EventList.
 *
 * @return by copy a vector of DateAndTime times
 */
std::vector<Mantid::Types::Core::DateAndTime> EventList::getPulseTimes() const {
  return this->eventList->getPulseTimes();
}

// --------------------------------------------------------------------------
/**
 * @return The minimum tof value for the list of the events.
 */
double EventList::getTofMin() const {
  return this->eventList->getTofMin();
}

/**
 * @return The maximum tof value for the list of events.
 */
double EventList::getTofMax() const {
  return this->eventList->getTofMax();
}

// --------------------------------------------------------------------------
/**
 * @return The minimum tof value for the list of the events.
 */
DateAndTime EventList::getPulseTimeMin() const {
  return this->eventList->getPulseTimeMin();
}

/**
 * @return The maximum tof value for the list of events.
 */
DateAndTime EventList::getPulseTimeMax() const {
  return this->eventList->getPulseTimeMax();
}

void EventList::getPulseTimeMinMax(Mantid::Types::Core::DateAndTime &tMin,
                                   Mantid::Types::Core::DateAndTime &tMax) const {
  return this->eventList->getPulseTimeMinMax(tMin, tMax);
}

DateAndTime EventList::getTimeAtSampleMax(const double &tofFactor, const double &tofOffset) const {
  return this->eventList->getTimeAtSampleMax(tofFactor, tofOffset);
}

DateAndTime EventList::getTimeAtSampleMin(const double &tofFactor, const double &tofOffset) const {
  return this->eventList->getTimeAtSampleMin(tofFactor, tofOffset);
}


// --------------------------------------------------------------------------
/**
 * Set a list of TOFs to the current event list. Modify the units if necessary.
 *
 * @param tofs :: The vector of doubles to set the tofs to.
 */
void EventList::setTofs(const MantidVec &tofs) {
  this->eventList->setTofs(tofs);
}

//------------------------------------------------------------------------------------------------
/** Operator to multiply the weights in this EventList by an error-less scalar.
 * Use multiply(value,error) if you wish to multiply by a real variable with an
 *error!
 *
 * The event list switches to WeightedEvent's if needed.
 * Note that if the multiplier is exactly 1.0, the list is NOT switched to
 *WeightedEvents - nothing happens.
 *
 * @param value :: multiply by this
 * @return reference to this
 */
EventList &EventList::operator*=(const double value) {
  *(this->eventList)*=value;
  return *this;
}

//------------------------------------------------------------------------------------------------
/** Multiply the weights in this event list by a scalar variable with an error;
 * though the error can be 0.0
 *
 * The event list switches to WeightedEvent's if needed.
 * Note that if the multiplier is exactly 1.0 and the error is exactly 0.0, the
 *list is NOT switched to WeightedEvents - nothing happens.
 *
 * Given:
 *  - A is the weight, variance \f$\sigma_A \f$
 *  - B is the scalar multiplier, variance \f$\sigma_B \f$
 *
 * The error propagation formula used is:
 *
 * \f[ \left(\frac{\sigma_f}{f}\right)^2 = \left(\frac{\sigma_A}{A}\right)^2 +
 *\left(\frac{\sigma_B}{B}\right)^2 + 2\frac{\sigma_A\sigma_B}{AB}\rho_{AB} \f]
 *
 * \f$ \rho_{AB} \f$ is the covariance between A and B, which we take to be 0
 *(uncorrelated variables).
 * Therefore, this reduces to:
 * \f[ \sigma_{AB}^2 = B^2 \sigma_A^2 + A^2 \sigma_B ^ 2  \f]
 *
 * In the case of no error:
 *  - The weight is simply \f$ aA \f$
 *  - The error \f$ \sigma_A \f$ becomes \f$ \sigma_{aA} = a \sigma_{A} \f$
 *
 * @param value: multiply all weights by this amount.
 * @param error: error on 'value'. Can be 0.
 */
void EventList::multiply(const double value, const double error) {
  this->eventList->multiply(value, error);
}

//------------------------------------------------------------------------------------------------
/** Multiply the weights in this event list by a histogram.
 * The event list switches to WeightedEvent's if needed.
 * NOTE: no unit checks are made (or possible to make) to compare the units of X
 *and tof() in the EventList.
 *
 * The formula used for calculating the error on the neutron weight is:
 * \f[ \sigma_{f}^2 = B^2 \sigma_A^2 + A^2 \sigma_B ^ 2  \f]
 *
 * where:
 *  * A is the weight of the event
 *  * B is the weight of the BIN that the event falls in
 *  * \f$\sigma_A\f$ is the error (not squared) of the weight of the event
 *  * \f$\sigma_B\f$ is the error (not squared) of the bin B
 *  * f is the resulting weight of the multiplied event
 *
 * @param X: bins of the multiplying histogram.
 * @param Y: value to multiply the weights.
 * @param E: error on the value to multiply.
 * @throw invalid_argument if the sizes of X, Y, E are not consistent.
 */
void EventList::multiply(const MantidVec &X, const MantidVec &Y, const MantidVec &E) {
  this->eventList->multiply(X, Y, E);
}

//------------------------------------------------------------------------------------------------
/** Divide the weights in this event list by a histogram.
 * The event list switches to WeightedEvent's if needed.
 * NOTE: no unit checks are made (or possible to make) to compare the units of X
 *and tof() in the EventList.
 *
 * The formula used for calculating the error on the neutron weight is:
 * \f[ \sigma_{f}^2 = (A / B)^2 * (\sigma_A^2 / A^2 + \sigma_B^2 / B^2) \f]
 *
 * where:
 *  * A is the weight of the event
 *  * B is the weight of the BIN that the event falls in
 *  * \f$\sigma_A\f$ is the error (not squared) of the weight of the event
 *  * \f$\sigma_B\f$ is the error (not squared) of the bin B
 *  * f is the resulting weight of the divided event
 *
 *
 * @param X: bins of the multiplying histogram.
 * @param Y: value to multiply the weights.
 * @param E: error on the value to multiply.
 * @throw invalid_argument if the sizes of X, Y, E are not consistent.
 */
void EventList::divide(const MantidVec &X, const MantidVec &Y, const MantidVec &E) {
  this->eventList->divide(X, Y, E);
}

//------------------------------------------------------------------------------------------------
/** Operator to divide the weights in this EventList by an error-less scalar.
 * Use divide(value,error) if your scalar has an error!
 * This simply calls the equivalent function: multiply(1.0/value).
 *
 * @param value :: divide by this
 * @return reference to this
 * @throw std::invalid_argument if value == 0; cannot divide by zero.
 */
EventList &EventList::operator/=(const double value) {
  *(this->eventList)/=value;
  return *this;
}

//------------------------------------------------------------------------------------------------
/** Divide the weights in this event list by a scalar with an (optional) error.
 * The event list switches to WeightedEvent's if needed.
 * This simply calls the equivalent function: multiply(1.0/value,
 *error/(value*value)).
 *
 * @param value: divide all weights by this amount.
 * @param error: error on 'value'. Can be 0.
 * @throw std::invalid_argument if value == 0; cannot divide by zero.
 */
void EventList::divide(const double value, const double error) {
  this->eventList->divide(value, error);
}

//------------------------------------------------------------------------------------------------
/** Filter this EventList into an output EventList, using
 * keeping only events within the >= start and < end pulse times.
 * Detector IDs and the X axis are copied as well.
 *
 * @param start :: start time (absolute)
 * @param stop :: end time (absolute)
 * @param output :: reference to an event list that will be output.
 * @throws std::invalid_argument If output is a reference to this EventList
 */
void EventList::filterByPulseTime(DateAndTime start, DateAndTime stop, EventList &output) const {
  this->eventList->filterByPulseTime(start, stop, output);
}

void EventList::filterByTimeAtSample(Types::Core::DateAndTime start, Types::Core::DateAndTime stop, double tofFactor,
                                     double tofOffset, EventList &output) const {
  this->eventList->filterByTimeAtSample(start, stop, tofFactor, tofOffset, output);
}

//------------------------------------------------------------------------------------------------
/** Use a TimeSplitterType to filter the event list in place.
 *
 * @param splitter :: a TimeSplitterType where all the entries (start/end time)
 *indicate events
 *     that will be kept. Any other events will be deleted.
 */
void EventList::filterInPlace(Kernel::TimeSplitterType &splitter) {
  this->eventList->filterInPlace(splitter);
}

//------------------------------------------------------------------------------------------------
/** Split the event list into n outputs
 *
 * @param splitter :: a TimeSplitterType giving where to split
 * @param outputs :: a vector of where the split events will end up. The # of
 *entries in there should
 *        be big enough to accommodate the indices.
 */
void EventList::splitByTime(Kernel::TimeSplitterType &splitter, std::vector<EventList *> outputs) const {
  this->eventList->splitByTime(splitter, outputs);
}


//------------------------------------------------------------------------------------------------
/** Split the event list into n outputs by event's full time (tof + pulse time)
 *
 * @param splitter :: a TimeSplitterType giving where to split
 * @param outputs :: a map of where the split events will end up. The # of
 *entries in there should
 *        be big enough to accommodate the indices.
 * @param docorrection :: a boolean to indiciate whether it is need to do
 *correction
 * @param toffactor:  a correction factor for each TOF to multiply with
 * @param tofshift:  a correction shift for each TOF to add with
 */
void EventList::splitByFullTime(Kernel::TimeSplitterType &splitter, std::map<int, EventList *> outputs,
                                bool docorrection, double toffactor, double tofshift) const {
  this->eventList->splitByFullTime(splitter, outputs, docorrection, toffactor, tofshift);
}


//----------------------------------------------------------------------------------------------
/**
 * @brief EventList::splitByFullTimeMatrixSplitter
 * @param vec_splitters_time  :: vector of splitting times
 * @param vecgroups :: vector of index group for splitters
 * @param vec_outputEventList :: vector of groups of splitted events
 * @param docorrection :: flag to do TOF correction from detector to sample
 * @param toffactor :: factor multiplied to TOF for correction
 * @param tofshift :: shift to TOF in unit of SECOND for correction
 * @return
 */
// TODO/FIXME/NOW - Consider to use vector to replace vec_outputEventList and
// have an option to ignore the un-filtered events!
std::string EventList::splitByFullTimeMatrixSplitter(const std::vector<int64_t> &vec_splitters_time,
                                                     const std::vector<int> &vecgroups,
                                                     std::map<int, EventList *> vec_outputEventList, bool docorrection,
                                                     double toffactor, double tofshift) const {
  return this->eventList->splitByFullTimeMatrixSplitter(vec_splitters_time, vecgroups, vec_outputEventList, docorrection, toffactor, tofshift);
  }

//----------------------------------------------------------------------------------------------
/** Split the event list by pulse time
 */
void EventList::splitByPulseTime(Kernel::TimeSplitterType &splitter, std::map<int, EventList *> outputs) const {
  this->eventList->splitByPulseTime(splitter, outputs);
}

//----------------------------------------------------------------------------------------------
/** Split the event list by pulse time
 */
// TODO/NOW - TEST
void EventList::splitByPulseTimeWithMatrix(const std::vector<int64_t> &vec_times, const std::vector<int> &vec_target,
                                           std::map<int, EventList *> outputs) const {
  this->eventList->splitByPulseTimeWithMatrix(vec_times, vec_target, outputs);
}

//--------------------------------------------------------------------------
/** Converts the X units in each event by going through TOF.
 * Note: if the unit conversion reverses the order, use "reverse()" to flip it
 *back.
 *
 * @param fromUnit :: the Unit describing the input unit. Must be initialized.
 * @param toUnit :: the Unit describing the output unit. Must be initialized.
 */
void EventList::convertUnitsViaTof(Mantid::Kernel::Unit *fromUnit, Mantid::Kernel::Unit *toUnit) {
  this->eventList->convertUnitsViaTof(fromUnit, toUnit);
}

//--------------------------------------------------------------------------
/** Convert the event's TOF (x) value according to a simple output = a *
 * (input^b) relationship
 *  @param factor :: the conversion factor a to apply
 *  @param power :: the Power b to apply to the conversion
 */
void EventList::convertUnitsQuickly(const double &factor, const double &power) {
  this->eventList->convertUnitsQuickly(factor, power);
}


void EventList::checkAndSanitizeHistogram(HistogramData::Histogram &histogram) {
  this->eventList->checkAndSanitizeHistogram(histogram);
}

void EventList::checkWorksWithPoints() const {
  this->eventList->checkWorksWithPoints();
}

void EventList::checkIsYAndEWritable() const {
  this->eventList->checkIsYAndEWritable();
}

//--------------------------------------------------------------------------
/** Get the vector of events contained in an EventListBase;
 * this is overloaded by event type.
 *
 * @param el :: The EventListBase to retrieve
 * @param[out] events :: reference to a pointer to a vector of this type of
 *event.
 *             The pointer will be set to point to the vector.
 * @throw runtime_error if you call this on the wrong type of EventListBase.
 */

// TODO: Fix before pr, disable for now to concentrate on main issues


void getEventsFrom(EventList &el, std::vector<TofEvent> *&events) { 
  // events = std::dynamic_pointer_cast<EventListTofEvent>(&el.eventList)->events;
}
void getEventsFrom(const EventList &el, std::vector<TofEvent> const *&events) {  
  
  // events = std::dynamic_pointer_cast<EventListTofEvent>(&el.eventList)->events;
  }

// void getEventsFrom(EventListBase &el, std::vector<Event> *&events) { events = &el.getEventsTyped(); }
// void getEventsFrom(const EventListBase &el, std::vector<Event> const *&events) { events = &el.getEventsTyped(); }


//--------------------------------------------------------------------------
/** Get the vector of events contained in an EventListBase;
 * this is overloaded by event type.
 *
 * @param el :: The EventListBase to retrieve
 * @param[out] events :: reference to a pointer to a vector of this type of
 *event.
 *             The pointer will be set to point to the vector.
 * @throw runtime_error if you call this on the wrong type of EventListBase.
 */
void getEventsFrom(EventList &el, std::vector<WeightedEvent> *&events) { 
  // events = std::dynamic_pointer_cast<EventListWeightedEvent>(&el.eventList)->events; 
  }
void getEventsFrom(const EventList &el, std::vector<WeightedEvent> const *&events) { 
  // events = std::dynamic_pointer_cast<EventListWeightedEvent>(&el.eventList)->events; 
  }

//--------------------------------------------------------------------------
/** Get the vector of events contained in an EventListBase;
 * this is overloaded by event type.
 *
 * @param el :: The EventListBase to retrieve
 * @param[out] events :: reference to a pointer to a vector of this type of
 *event.
 *             The pointer will be set to point to the vector.
 * @throw runtime_error if you call this on the wrong type of EventListBase.
 */
void getEventsFrom(EventList &el, std::vector<WeightedEventNoTime> *&events) { 
  // events = std::dynamic_pointer_cast<EventListWeightedEventNoTime>(&el.eventList)->events; 
  }
void getEventsFrom(const EventList &el, std::vector<WeightedEventNoTime> const *&events) {
  // events = std::dynamic_pointer_cast<EventListWeightedEventNoTime>(&el.eventList)->events;
}

} // namespace Mantid::DataObjects
