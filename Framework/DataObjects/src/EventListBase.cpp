// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#include "MantidDataObjects/EventListBase.h"
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

//==========================================================================
/// --------------------- TofEvent Comparators
/// ----------------------------------
//==========================================================================
/** Compare two events' FRAME id, return true if e1 should be before e2.
 * @param e1 :: first event
 * @param e2 :: second event
 *  */
bool compareEventPulseTime(const TofEvent &e1, const TofEvent &e2) { return (e1.pulseTime() < e2.pulseTime()); }

/** Compare two events' FRAME id, return true if e1 should be before e2.
 *  Assuming that if e1's pulse time is earlier than e2's, then e1 must be
 * earlier regardless TOF value
 * @param e1 :: first event
 * @param e2 :: second event
 *  */
bool compareEventPulseTimeTOF(const TofEvent &e1, const TofEvent &e2) {

  if (e1.pulseTime() < e2.pulseTime()) {
    return true;
  } else if ((e1.pulseTime() == e2.pulseTime()) && (e1.tof() < e2.tof())) {
    return true;
  }

  return false;
}

// comparator for pulse time with tolerance
struct comparePulseTimeTOFDelta {
  explicit comparePulseTimeTOFDelta(const Types::Core::DateAndTime &start, const double seconds)
      : startNano(start.totalNanoseconds()), deltaNano(static_cast<int64_t>(seconds * SEC_TO_NANO)) {}

  bool operator()(const TofEvent &e1, const TofEvent &e2) {
    // get the pulse times converted into bin number from start time
    const int64_t e1Pulse = (e1.pulseTime().totalNanoseconds() - startNano) / deltaNano;
    const int64_t e2Pulse = (e2.pulseTime().totalNanoseconds() - startNano) / deltaNano;

    // compare with the calculated bin information
    if (e1Pulse < e2Pulse) {
      return true;
    } else if ((e1Pulse == e2Pulse) && (e1.tof() < e2.tof())) {
      return true;
    }

    return false;
  }

  int64_t startNano;
  int64_t deltaNano;
};

/// Constructor (empty)
// EventWorkspace is always histogram data and so is thus EventListBase
EventListBase::EventListBase()
    : m_histogram(HistogramData::Histogram::XMode::BinEdges, HistogramData::Histogram::YMode::Counts), eventType(TOF),
      order(UNSORTED), mru(nullptr) {}

/** Constructor with a MRU list
 * @param mru :: pointer to the MRU of the parent EventWorkspace
 * @param specNo :: the spectrum number for the event list
 */
EventListBase::EventListBase(EventWorkspaceMRU *mru, specnum_t specNo)
    : IEventList(specNo),
      m_histogram(HistogramData::Histogram::XMode::BinEdges, HistogramData::Histogram::YMode::Counts), eventType(TOF),
      order(UNSORTED), mru(mru) {}

/** Constructor copying from an existing event list
 * @param rhs :: EventListBase object to copy*/
EventListBase::EventListBase(const EventListBase &rhs) : IEventList(rhs), m_histogram(rhs.m_histogram), mru{nullptr} {
  // Note that operator= also assigns m_histogram, but the above use of the copy
  // constructor avoid a memory allocation and is thus faster.
  this->operator=(rhs);
}

/** Constructor, taking a vector of events->
 * @param events :: Vector of TofEvent's */
EventListBase::EventListBase(const std::vector<TofEvent> &events)
    : m_histogram(HistogramData::Histogram::XMode::BinEdges, HistogramData::Histogram::YMode::Counts), eventType(TOF),
      mru(nullptr) {
  this->events->assign(events.begin(), events.end());
  this->eventType = TOF;
  this->order = UNSORTED;
}

/** Constructor, taking a vector of events->
 * @param events :: Vector of WeightedEvent's */
EventListBase::EventListBase(const std::vector<WeightedEvent> &events)
    : m_histogram(HistogramData::Histogram::XMode::BinEdges, HistogramData::Histogram::YMode::Counts), mru(nullptr) {
  this->events->assign(events.begin(), events.end());
  this->eventType = WEIGHTED;
  this->order = UNSORTED;
}

/** Constructor, taking a vector of events->
 * @param events :: Vector of WeightedEventNoTime's */
EventListBase::EventListBase(const std::vector<WeightedEventNoTime> &events)
    : m_histogram(HistogramData::Histogram::XMode::BinEdges, HistogramData::Histogram::YMode::Counts), mru(nullptr) {
  this->events->assign(events.begin(), events.end());
  this->eventType = WEIGHTED_NOTIME;
  this->order = UNSORTED;
}

/** Constructor, taking a vector of events->
 * @param events :: Vector of TofEventNoTime's */
EventListBase::EventListBase(const std::vector<TofEventNoTime> &events)
    : m_histogram(HistogramData::Histogram::XMode::BinEdges, HistogramData::Histogram::YMode::Counts), mru(nullptr) {
  this->events->assign(events.begin(), events.end());
  this->eventType = UNWEIGHTED;
  this->order = UNSORTED;
}

/// Destructor
EventListBase::~EventListBase() {
  // Note: These two lines do not seem to have an effect on releasing memory
  //  at least on Linux. (Memory usage seems to increase event after deleting
  //  EventWorkspaces.
  //  Therefore, for performance, they are kept commented:
  clear();

  // this->events->clear();
  // std::vector<TofEvent>().swap(events); //Trick to release the vector memory.
}

/// Copy data from another EventListBase, via ISpectrum reference.
void EventListBase::copyDataFrom(const ISpectrum &source) { source.copyDataInto(*this); }

/// Used by copyDataFrom for dynamic dispatch for its `source`.
void EventListBase::copyDataInto(EventListBase &sink) const {
  sink.m_histogram = m_histogram;
  sink.events = events;
  sink.eventType = eventType;
  sink.order = order;
}

/// Used by Histogram1D::copyDataFrom for dynamic dispatch for `other`.
void EventListBase::copyDataInto(Histogram1D &sink) const { sink.setHistogram(histogram()); }

// --------------------------------------------------------------------------
/** Create an EventListBase from a histogram. This converts bins to weighted
 * events->
 * Any existing events are cleared.
 *
 * @param inSpec :: ISpectrum ptr to histogram data.
 * @param GenerateZeros :: if true, generate event(s) for empty bins
 * @param GenerateMultipleEvents :: if true, create several evenly-spaced fake
 *events inside the bin
 * @param MaxEventsPerBin :: max number of events to generate in one bin, if
 *GenerateMultipleEvents
 */
void EventListBase::createFromHistogram(const ISpectrum *inSpec, bool GenerateZeros, bool GenerateMultipleEvents,
                                    int MaxEventsPerBin) {
  // Fresh start
  this->clear(true);
  // Get the input histogram
  const MantidVec &X = inSpec->readX();
  const MantidVec &Y = inSpec->readY();
  const MantidVec &E = inSpec->readE();
  if (Y.size() + 1 != X.size())
    throw std::runtime_error("Expected a histogram (X vector should be 1 longer than the Y vector)");

  // Copy detector IDs and spectra
  this->copyInfoFrom(*inSpec);
  // We need weights but have no way to set the time. So use weighted, no time
  this->switchTo(WEIGHTED_NOTIME);
  if (GenerateZeros)
    this->events->reserve(Y.size());

  for (size_t i = 0; i < X.size() - 1; i++) {
    double weight = Y[i];
    if ((weight != 0.0 || GenerateZeros) && std::isfinite(weight)) {
      double error = E[i];
      // Also check that the error is not a bad number
      if (std::isfinite(error)) {
        if (GenerateMultipleEvents) {
          // --------- Multiple events per bin ----------
          double errorSquared = error * error;
          // Find how many events to fake
          double val = weight / E[i];
          val *= val;
          // Convert to int with slight rounding up. This is to avoid rounding
          // errors
          auto numEvents = int(val + 0.2);
          if (numEvents < 1)
            numEvents = 1;
          if (numEvents > MaxEventsPerBin)
            numEvents = MaxEventsPerBin;
          // Scale the weight and error for each
          weight /= numEvents;
          errorSquared /= numEvents;

          // Spread the TOF. e.g. 2 events = 0.25, 0.75.
          double tofStep = (X[i + 1] - X[i]) / (numEvents);
          for (size_t j = 0; j < size_t(numEvents); j++) {
            double tof = X[i] + tofStep * (0.5 + double(j));
            // Create and add the event
            // TODO: try emplace_back() here.
            events->emplace_back(tof, weight, errorSquared);
          }
        } else {
          // --------- Single event per bin ----------
          // TOF = midpoint of the bin
          double tof = (X[i] + X[i + 1]) / 2.0;
          // Error squared is carried in the event
          double errorSquared = E[i];
          errorSquared *= errorSquared;
          // Create and add the event
          events->emplace_back(tof, weight, errorSquared);
        }
      } // error is nont NAN or infinite
    }   // weight is non-zero, not NAN, and non-infinite
  }     // (each bin)

  // Set the X binning parameters
  this->setX(inSpec->ptrX());

  // Manually set that this is sorted by TOF, since it is. This will make it
  // "threadSafe" in other algos.
  this->setSortOrder(TOF_SORT);
}

// --------------------------------------------------------------------------
// --- Operators
// -------------------------------------------------------------------

/** Copy into this event list from another
 * @param rhs :: We will copy all the events from that into this object.
 * @return reference to this
 * */
EventListBase &EventListBase::operator=(const EventListBase &rhs) {
  // Note that we are NOT copying the MRU pointer
  // the EventWorkspace that posseses the EventListBase has already configured the mru
  IEventList::operator=(rhs);
  m_histogram = rhs.m_histogram;
  events = rhs.events;
  eventType = rhs.eventType;
  order = rhs.order;
  return *this;
}

// --------------------------------------------------------------------------
/** Append an event to the histogram.
 * @param event :: TofEvent to add at the end of the list.
 * @return reference to this
 * */
EventListBase &EventListBase::operator+=(const TofEvent &event) {
  this->events->emplace_back(event);
  this->order = UNSORTED;
  return *this;
}

// --------------------------------------------------------------------------
/** Append a list of events to the histogram.
 * The internal event list will switch to the required type.
 *
 * @param more_events :: A vector of events to append.
 * @return reference to this
 * */
EventListBase &EventListBase::operator+=(const std::vector<TofEvent> &more_events) {

  // case TOF:
  //   // Simply push the events
  //   this->events->insert(this->events->end(), more_events->begin(), more_events->end());
  //   break;

  this->events->reserve(this->events->size() + more_events->size());
  for (const auto &event : more_events) {
    this->events->emplace_back(event);
  }

  this->order = UNSORTED;
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
EventListBase &EventListBase::operator+=(const WeightedEvent &event) {
  this->switchTo(WEIGHTED);
  this->events->emplace_back(event);
  this->order = UNSORTED;
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
EventListBase &EventListBase::operator+=(const std::vector<WeightedEvent> &more_events) {
  
  // TOF seems to have this thing where if it recieves a weighted event it just switches to weighted mode


  // case TOF:
  //   // Need to switch to weighted
  //   this->switchTo(WEIGHTED);
  //   // Fall through to the insertion!

  // case WEIGHTED:
  //   // Append the two lists
  //   this->weightedevents->insert(weightedevents->end(), more_events->begin(), more_events->end());
  //   break;
  this->events->reserve(this->events->size() + more_events->size());
  for (const auto &event : more_events) {
    this->events->emplace_back(event);
  }


  this->order = UNSORTED;
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
EventListBase &EventListBase::operator+=(const std::vector<WeightedEventNoTime> &more_events) {
  // TOF and Weighted switch to weighted no time when given events from such, perhaps unweighted with time should follow suit

  // case TOF:
  // case WEIGHTED:
  //   // Need to switch to weighted with no time
  //   this->switchTo(WEIGHTED_NOTIME);
  //   // Fall through to the insertion!

  this->events->insert(events->end(), more_events->begin(), more_events->end());
  this->order = UNSORTED;
  return *this;
}

// --------------------------------------------------------------------------
/** Append another EventListBase to this event list.
 * The event lists are concatenated, and a union of the sets of detector ID's is
 *done.
 * Switching of event types may occur if the two are different.
 *
 * @param more_events :: Another EventListBase.
 * @return reference to this
 * */
EventListBase &EventListBase::operator+=(const EventListBase &more_events) {
  // We'll let the += operator for the given vector of event lists handle it
  this->operator+=(more_events->events);


  // No guaranteed order
  this->order = UNSORTED;
  // Do a union between the detector IDs of both lists
  addDetectorIDs(more_events->getDetectorIDs());

  return *this;
}

// --------------------------------------------------------------------------
/** SUBTRACT another EventListBase from this event list.
 * The event lists are concatenated, but the weights of the incoming
 *    list are multiplied by -1.0.
 *
 * @param more_events :: Another EventListBase.
 * @return reference to this
 * */
EventListBase &EventListBase::operator-=(const EventListBase &more_events) {
  if (this == &more_events) {
    // Special case, ticket #3844 part 2.
    // When doing this = this - this,
    // simply clear the input event list. Saves memory!
    this->clearData();
    return *this;
  }

  // We'll let the -= operator for the given vector of event lists handle it
  // case TOF:
  //   this->switchTo(WEIGHTED);
  //   // Fall through


    // case WEIGHTED_NOTIME:
    //   // TODO: Should this throw?
    //   minusHelper(this->weightedEvents, more_events->weightedEventsNoTime);
    //   break;
    // case UNWEIGHTED:
    //   // TODO: Should this throw?
    //   minusHelper(this->weightedEvents, more_events->unweightedEvents);
    //   break;

  minusHelper(*(this->events), more_events->events);
  // No guaranteed order
  this->order = UNSORTED;

  // NOTE: What to do about detector ID's?
  return *this;
}

// --------------------------------------------------------------------------
/** Equality operator between EventListBase's
 * @param rhs :: other EventListBase to compare
 * @return :: true if equal.
 */
bool EventListBase::operator==(const EventListBase &rhs) const {
  return this->getNumberEvents() == rhs.getNumberEvents()
  && this->eventType == rhs.eventType
  && *events == *(rhs.events);
}

/** Inequality comparator
 * @param rhs :: other EventListBase to compare
 * @return :: true if not equal.
 */
bool EventListBase::operator!=(const EventListBase &rhs) const { return (!this->operator==(rhs)); }

void throwUnimplementedError() {
  throw std::logic_error("Function not implemented, please use a derived class and not the base class.");
}

bool EventListBase::equals(const EventListBase &rhs, const double tolTof, const double tolWeight,
                       const int64_t tolPulse) const {
  
  // This is overriden in implemented classes, should not be used
  throwUnimplementedError();
}

// -----------------------------------------------------------------------------------------------
/** Return the type of Event vector contained within.
 * @return :: a EventType value.
 */
EventType EventListBase::getEventType() const { return eventType; }

// -----------------------------------------------------------------------------------------------
/** Switch the EventListBase to use the given EventType (TOF, WEIGHTED, or
 * WEIGHTED_NOTIME)
 */


// TODO: Extract out to wrapper
void EventListBase::switchTo(EventType newType) {
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
  // Make sure to free memory
  this->clearUnused();
}

// -----------------------------------------------------------------------------------------------
/** Switch the EventListBase to use WeightedEvents instead
 * of TofEvent.
 */
void EventListBase::switchToWeightedEvents() {
  switch (eventType) {
  case WEIGHTED:
    // Do nothing; it already is weighted
    return;

  case WEIGHTED_NOTIME:
    throw std::runtime_error("EventListBase::switchToWeightedEvents() called on an "
                             "EventListBase with WeightedEventNoTime's. It has "
                             "lost the pulse time information and can't go "
                             "back to WeightedEvent's.");
    break;

  case UNWEIGHTED:
    throw std::runtime_error("TestTestTestTest");
    break;

  case TOF:
    weightedEventsNoTime.clear();
    // Convert and copy all TofEvents to the weightedEvents list.
    this->weightedevents->assign(events->cbegin(), events->cend());
    // Get rid of the old events
    events->clear();
    eventType = WEIGHTED;
    break;
  }
}

// -----------------------------------------------------------------------------------------------
/** Switch the EventListBase to use WeightedEventNoTime's instead
 * of TofEvent.
 */
void EventListBase::switchToWeightedEventsNoTime() {
  switch (eventType) {
  case WEIGHTED_NOTIME:
    // Do nothing if already there
    return;

  case TOF: {
    // Convert and copy all TofEvents to the weightedEvents list.
    this->weightedEventsNoTime.assign(events->cbegin(), events->cend());
    // Get rid of the old events
    events->clear();
    weightedevents->clear();
    eventType = WEIGHTED_NOTIME;
  } break;

  case WEIGHTED: {
    // Convert and copy all TofEvents to the weightedEvents list.
    this->weightedEventsNoTime.assign(weightedevents->cbegin(), weightedevents->cend());
    // Get rid of the old events
    events->clear();
    weightedevents->clear();
    eventType = WEIGHTED_NOTIME;
  } break;

  case UNWEIGHTED: {
    // Convert and copy all TofEvents to the unweightedEvents list.
    this->unweightedevents->assign(unweightedevents->cbegin(), unweightedevents->cend());
    // Get rid of the old events
    events->clear();
    unweightedevents->clear();
    eventType = UNWEIGHTED;
  } break;
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
WeightedEvent EventListBase::getEvent(size_t event_number) {
  throwUnimplementedError();
}

// ==============================================================================================
// --- Handling the event list
// -------------------------------------------------------------------
// ==============================================================================================

/** Return the const list of TofEvents contained.
 * NOTE! This should be used for testing purposes only, as much as possible. The
 *EventListBase
 * may contain weighted events, requiring use of getWeightedEvents() instead.
 *
 * @return a const reference to the list of non-weighted events
 * */
const std::vector<TofEvent> &EventListBase::getEvents() const {
  if (eventType != TOF)
    throw std::runtime_error("EventListBase::getEvents() called for an EventListBase "
                             "that has weights. Use getWeightedEvents() or "
                             "getWeightedEventsNoTime().");
  return this->events;
}

/** Return the list of TofEvents contained.
 * NOTE! This should be used for testing purposes only, as much as possible. The
 *EventListBase
 * may contain weighted events, requiring use of getWeightedEvents() instead.
 *
 * @return a reference to the list of non-weighted events
 * */
std::vector<TofEvent> &EventListBase::getEvents() {
  if (eventType != TOF)
    throw std::runtime_error("EventListBase::getEvents() called for an EventListBase "
                             "that has weights. Use getWeightedEvents() or "
                             "getWeightedEventsNoTime().");
  return this->events;
}

// std::vector<Event> &EventListBase::getEventsTyped() {
//   return this->events;
// }

// const std::vector<Event> &EventListBase::getEventsTyped() const {
//   return this->events;
// }

/** Return the list of WeightedEvent contained.
 * NOTE! This should be used for testing purposes only, as much as possible. The
 *EventListBase
 * may contain un-weighted events, requiring use of getEvents() instead.
 *
 * @return a reference to the list of weighted events
 * */
std::vector<WeightedEvent> &EventListBase::getWeightedEvents() {
  if (eventType != WEIGHTED)
    throw std::runtime_error("EventListBase::getWeightedEvents() called for an "
                             "EventListBase not of type WeightedEvent. Use "
                             "getEvents() or getWeightedEventsNoTime().");
  return this->events;
}

/** Return the list of WeightedEvent contained.
 * NOTE! This should be used for testing purposes only, as much as possible. The
 *EventListBase
 * may contain un-weighted events, requiring use of getEvents() instead.
 *
 * @return a const reference to the list of weighted events
 * */
const std::vector<WeightedEvent> &EventListBase::getWeightedEvents() const {
  if (eventType != WEIGHTED)
    throw std::runtime_error("EventListBase::getWeightedEvents() called for an "
                             "EventListBase not of type WeightedEvent. Use "
                             "getEvents() or getWeightedEventsNoTime().");
  return this->events;
}

/** Return the list of WeightedEvent contained.
 * NOTE! This should be used for testing purposes only, as much as possible.
 *
 * @return a reference to the list of weighted events
 * */
std::vector<WeightedEventNoTime> &EventListBase::getWeightedEventsNoTime() {
  if (eventType != WEIGHTED_NOTIME)
    throw std::runtime_error("EventListBase::getWeightedEvents() called for an "
                             "EventListBase not of type WeightedEventNoTime. Use "
                             "getEvents() or getWeightedEvents().");
  return this->events;
}

// /** Return the list of WeightedEvent contained.
//  * NOTE! This should be used for testing purposes only, as much as possible.
//  *
//  * @return a reference to the list of weighted events
//  * */
// std::vector<WeightedEventNoTime> &EventListBase::getUnweightedEvents() {
//   if (eventType != UNWEIGHTED)
//     throw std::runtime_error("EventListBase::getWeightedEvents() called for an "
//                              "EventListBase not of type WeightedEventNoTime. Use "
//                              "getEvents() or getWeightedEvents().");
//   return this->events;
// }

/** Return the list of WeightedEventNoTime contained.
 * NOTE! This should be used for testing purposes only, as much as possible.
 *
 * @return a const reference to the list of weighted events
 * */
const std::vector<WeightedEventNoTime> &EventListBase::getWeightedEventsNoTime() const {
  if (eventType != WEIGHTED_NOTIME)
    throw std::runtime_error("EventListBase::getWeightedEventsNoTime() called for "
                             "an EventListBase not of type WeightedEventNoTime. "
                             "Use getEvents() or getWeightedEvents().");
  return this->events;
}

/** Clear the list of events and any
 * associated detector ID's.
 * */
void EventListBase::clear(const bool removeDetIDs) {
  if (mru)
    mru->deleteIndex(this);
  this->events->clear();
  std::vector<Event>().swap(this->events); // STL Trick to release memory
  if (removeDetIDs)
    this->clearDetectorIDs();
}

/** Clear any unused event lists (the ones that do not
 * match the currently used type).
 * Memory is freed.
 * */
void EventListBase::clearUnused() {
  // should be unnecessary now that the different event types have been split
  // else move this up to the wrapper
}

/// Mask the spectrum to this value. Removes all events->
void EventListBase::clearData() { this->clear(false); }

/** Sets the MRU list for this event list
 *
 * @param newMRU :: new MRU for the workspace containing this EventListBase
 */
void EventListBase::setMRU(EventWorkspaceMRU *newMRU) { mru = newMRU; }

/** Reserve a certain number of entries in event list of the specified eventType
 *
 * Calls std::vector<>::reserve() in order to pre-allocate the length of the
 *event list vector.
 *
 * @param num :: number of events that will be in this EventListBase
 */
void EventListBase::reserve(size_t num) {
  this->events->reserve(num);
}

// ==============================================================================================
// --- Sorting functions -----------------------------------------------------
// ==============================================================================================

// --------------------------------------------------------------------------
/** Sort events by TOF or Frame
 * @param order :: Order by which to sort.
 * */
void EventListBase::sort(const EventSortType order) const {
  if (order == UNSORTED) {
    return; // don't bother doing anything. Why did you ask to unsort?
  } else if (order == TOF_SORT) {
    this->sortTof();
  } else if (order == PULSETIME_SORT) {
    this->sortPulseTime();
  } else if (order == PULSETIMETOF_SORT) {
    this->sortPulseTimeTOF();
  } else if (order == PULSETIMETOF_DELTA_SORT) {
    throw std::invalid_argument("sorting by pulse time with delta requires "
                                "extra parameters. Use sortPulseTimeTOFDelta "
                                "instead.");
  } else if (order == TIMEATSAMPLE_SORT) {
    throw std::invalid_argument("sorting by time at sample requires extra "
                                "parameters. Use sortTimeAtSample instead.");
  } else {
    throw runtime_error("Invalid sort type in EventListBase::sort(EventSortType)");
  }
}

// --------------------------------------------------------------------------
/** Manually set the event list sort order value. No actual sorting takes place.
 * SHOULD ONLY BE USED IN TESTS or if you know what you are doing.
 * @param order :: sort order to set.
 */
void EventListBase::setSortOrder(const EventSortType order) const { this->order = order; }

// --------------------------------------------------------------------------
/** Sort events by TOF in one thread */
void EventListBase::sortTof() const {
  // nothing to do
  if (this->order == TOF_SORT)
    return;

  // Avoid sorting from multiple threads
  std::lock_guard<std::mutex> _lock(m_sortMutex);
  // If the list was sorted while waiting for the lock, return.
  if (this->order == TOF_SORT)
    return;

  // TODO determine how these are setup to compare

  tbb::parallel_sort(events->begin(), events->end());
  // Save the order to avoid unnecessary re-sorting.
  this->order = TOF_SORT;
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
void EventListBase::sortTimeAtSample(const double &tofFactor, const double &tofShift, bool forceResort) const {
  throwUnimplementedError();
}

// --------------------------------------------------------------------------
/** Sort events by Frame */
void EventListBase::sortPulseTime() const {
  if (this->order == PULSETIME_SORT)
    return; // nothing to do

  // Avoid sorting from multiple threads
  std::lock_guard<std::mutex> _lock(m_sortMutex);
  // If the list was sorted while waiting for the lock, return.
  if (this->order == PULSETIME_SORT)
    return;

  tbb::parallel_sort(events->begin(), events->end(), compareEventPulseTime);
 
  // Save the order to avoid unnecessary re-sorting.
  this->order = PULSETIME_SORT;
}

/*
 * Sort events by pulse time + TOF
 * (the absolute time)
 */
void EventListBase::sortPulseTimeTOF() const {
  if (this->order == PULSETIMETOF_SORT)
    return; // already ordered

  // Avoid sorting from multiple threads
  std::lock_guard<std::mutex> _lock(m_sortMutex);
  // If the list was sorted while waiting for the lock, return.
  if (this->order == PULSETIMETOF_SORT)
    return;

  tbb::parallel_sort(events->begin(), events->end(), compareEventPulseTimeTOF);

  // Save
  this->order = PULSETIMETOF_SORT;
}

/**
 * Sort by the pulse time with a tolerance. The pulsetime to compare is a
 * constant binning of seconds from start. This will set the sort order to
 * UNSORTED upon completion rather than storing the call parameters.
 * @param start The absolute start time
 * @param seconds The tolerance of pulse time in seconds.
 */
void EventListBase::sortPulseTimeTOFDelta(const Types::Core::DateAndTime &start, const double seconds) const {
  // Avoid sorting from multiple threads
  std::lock_guard<std::mutex> _lock(m_sortMutex);

  std::function<bool(const TofEvent &, const TofEvent &)> comparator = comparePulseTimeTOFDelta(start, seconds);

  tbb::parallel_sort(events->begin(), events->end(), comparator);
  
  this->order = UNSORTED; // so the function always re-runs
}

// --------------------------------------------------------------------------
/** Return true if the event list is sorted by TOF */
bool EventListBase::isSortedByTof() const { return (this->order == TOF_SORT); }

// --------------------------------------------------------------------------
/** Return the type of sorting used in this event list */
EventSortType EventListBase::getSortType() const { return this->order; }

// --------------------------------------------------------------------------
/** Reverse the histogram boundaries and the associated events if they are
 * sorted
 * by time-of-flight.
 * Does nothing if sorted otherwise or unsorted.
 * */
void EventListBase::reverse() {
  // reverse the histogram bin parameters
  MantidVec &x = dataX();
  std::reverse(x.begin(), x.end());

  // flip the events if they are tof sorted
  if (this->isSortedByTof()) {
    std::reverse(this->events->begin(), this->events->end());
    // And we are still sorted! :)
  }
  // Otherwise, do nothing. If it was sorted by pulse time, then it still is
}

// --------------------------------------------------------------------------
/** Return the number of events in the list.
 * NOTE: If the events have weights, this returns the NUMBER of WeightedEvent's
 *in the
 * list, and NOT the sum of their weights (which may be two different numbers).
 *
 * @return the number of events in the list.
 *  */
size_t EventListBase::getNumberEvents() const {
    return this->events->size();
}

/**
 * Much like stl containers, returns true if there is nothing in the event list.
 */
bool EventListBase::empty() const {
  return this->events->empty();
}

// --------------------------------------------------------------------------
/** Memory used by this event list. Note: It reports the CAPACITY of the
 * vectors, rather than their size, since that is a more accurate
 * representation of the size used.
 *
 * @return :: the memory used by the EventListBase, in bytes.
 * */
size_t EventListBase::getMemorySize() const {
  throwUnimplementedError();
}

// --------------------------------------------------------------------------
/** Return the size of the histogram data.
 * @return the size of the histogram representation of the data (size of Y) **/
size_t EventListBase::histogram_size() const {
  size_t x_size = readX().size();
  if (x_size > 1)
    return x_size - 1;
  else
    return 0;
}

// ==============================================================================================
// --- Setting the Histogram X axis, without recalculating the histogram
// -----------------------
// ==============================================================================================

/** Deprecated, use setSharedX() instead. Set the x-component for the histogram
 * view. This will NOT cause the histogram to be calculated.
 * @param X :: The vector of doubles to set as the histogram limits.
 */
void EventListBase::setX(const Kernel::cow_ptr<HistogramData::HistogramX> &X) {
  m_histogram.setX(X);
  if (mru)
    mru->deleteIndex(this);
}

/** Deprecated, use mutableX() instead. Returns a reference to the x data.
 *  @return a reference to the X (bin) vector.
 */
MantidVec &EventListBase::dataX() {
  if (mru)
    mru->deleteIndex(this);
  return m_histogram.dataX();
}

/** Deprecated, use x() instead. Returns a const reference to the x data.
 *  @return a reference to the X (bin) vector. */
const MantidVec &EventListBase::dataX() const { return m_histogram.dataX(); }

/// Deprecated, use x() instead. Returns the x data const
const MantidVec &EventListBase::readX() const { return m_histogram.readX(); }

/// Deprecated, use sharedX() instead. Returns a pointer to the x data
Kernel::cow_ptr<HistogramData::HistogramX> EventListBase::ptrX() const { return m_histogram.ptrX(); }

/// Deprecated, use mutableDx() instead.
MantidVec &EventListBase::dataDx() { return m_histogram.dataDx(); }
/// Deprecated, use dx() instead.
const MantidVec &EventListBase::dataDx() const { return m_histogram.dataDx(); }
/// Deprecated, use dx() instead.
const MantidVec &EventListBase::readDx() const { return m_histogram.readDx(); }

// ==============================================================================================
// --- Return Data Vectors --------------------------------------------------
// ==============================================================================================

/** Calculates and returns a pointer to the Y histogrammed data.
 * Remember to delete your pointer after use!
 *
 * @return a pointer to a MantidVec
 */
MantidVec *EventListBase::makeDataY() const {
  auto Y = new MantidVec();
  MantidVec E;
  // Generate the Y histogram while skipping the E if possible.
  generateHistogram(readX(), *Y, E, true);
  return Y;
}

/** Calculates and returns a pointer to the E histogrammed data.
 * Remember to delete your pointer after use!
 *
 * @return a pointer to a MantidVec
 */
MantidVec *EventListBase::makeDataE() const {
  MantidVec Y;
  auto E = new MantidVec();
  generateHistogram(readX(), Y, *E);
  // Y is unused.
  return E;
}

HistogramData::Histogram EventListBase::histogram() const {
  HistogramData::Histogram ret(m_histogram);
  ret.setSharedY(sharedY());
  ret.setSharedE(sharedE());
  return ret;
}

HistogramData::Counts EventListBase::counts() const { return histogram().counts(); }

HistogramData::CountVariances EventListBase::countVariances() const { return histogram().countVariances(); }

HistogramData::CountStandardDeviations EventListBase::countStandardDeviations() const {
  return histogram().countStandardDeviations();
}

HistogramData::Frequencies EventListBase::frequencies() const { return histogram().frequencies(); }

HistogramData::FrequencyVariances EventListBase::frequencyVariances() const { return histogram().frequencyVariances(); }

HistogramData::FrequencyStandardDeviations EventListBase::frequencyStandardDeviations() const {
  return histogram().frequencyStandardDeviations();
}

const HistogramData::HistogramY &EventListBase::y() const {
  if (!mru)
    throw std::runtime_error("'EventListBase::y()' called with no MRU set. This is not allowed.");

  return *sharedY();
}
const HistogramData::HistogramE &EventListBase::e() const {
  if (!mru)
    throw std::runtime_error("'EventListBase::e()' called with no MRU set. This is not allowed.");

  return *sharedE();
}
Kernel::cow_ptr<HistogramData::HistogramY> EventListBase::sharedY() const {
  // This is the thread number from which this function was called.
  int thread = PARALLEL_THREAD_NUMBER;

  Kernel::cow_ptr<HistogramData::HistogramY> yData(nullptr);

  // Is the data in the mrulist?
  if (mru) {
    mru->ensureEnoughBuffersY(thread);
    yData = mru->findY(thread, this);
  }

  if (!yData) {
    MantidVec Y;
    MantidVec E;
    this->generateHistogram(readX(), Y, E);

    // Create the MRU object
    yData = Kernel::make_cow<HistogramData::HistogramY>(std::move(Y));

    // Lets save it in the MRU
    if (mru) {
      mru->insertY(thread, yData, this);
      auto eData = Kernel::make_cow<HistogramData::HistogramE>(std::move(E));
      mru->ensureEnoughBuffersE(thread);
      mru->insertE(thread, eData, this);
    }
  }
  return yData;
}
Kernel::cow_ptr<HistogramData::HistogramE> EventListBase::sharedE() const {
  // This is the thread number from which this function was called.
  int thread = PARALLEL_THREAD_NUMBER;

  Kernel::cow_ptr<HistogramData::HistogramE> eData(nullptr);

  // Is the data in the mrulist?
  if (mru) {
    mru->ensureEnoughBuffersE(thread);
    eData = mru->findE(thread, this);
  }

  if (!eData) {
    // Now use that to get E -- Y values are generated from another function
    MantidVec Y_ignored;
    MantidVec E;
    this->generateHistogram(readX(), Y_ignored, E);
    eData = Kernel::make_cow<HistogramData::HistogramE>(std::move(E));

    // Lets save it in the MRU
    if (mru)
      mru->insertE(thread, eData, this);
  }
  return eData;
}
/** Look in the MRU to see if the Y histogram has been generated before.
 * If so, return that. If not, calculate, cache and return it.
 *
 * @return reference to the Y vector.
 */
const MantidVec &EventListBase::dataY() const {
  if (!mru)
    throw std::runtime_error("'EventListBase::dataY()' called with no MRU set. This is not allowed.");

  // WARNING: The Y data of sharedY() is stored in MRU, returning reference fine
  // as long as it stays there.
  return sharedY()->rawData();
}

/** Look in the MRU to see if the E histogram has been generated before.
 * If so, return that. If not, calculate, cache and return it.
 *
 * @return reference to the E vector.
 */
const MantidVec &EventListBase::dataE() const {
  if (!mru)
    throw std::runtime_error("'EventListBase::dataE()' called with no MRU set. This is not allowed.");

  // WARNING: The E data of sharedE() is stored in MRU, returning reference fine
  // as long as it stays there.
  return sharedE()->rawData();
}

namespace {
inline double calcNorm(const double errorSquared) {
  if (errorSquared == 0.)
    return 0;
  else if (errorSquared == 1.)
    return 1.;
  else
    return 1. / std::sqrt(errorSquared);
}
} // namespace




// --------------------------------------------------------------------------
/** Compress the event list by grouping events with the same
 * TOF (within a given tolerance). PulseTime is ignored.
 * The event list will be switched to WeightedEventNoTime.
 *
 * @param tolerance :: how close do two event's TOF have to be to be considered
 *the same.
 * @param destination :: EventListBase that will receive the compressed events-> Can
 *be == this.
 */
void EventListBase::compressEvents(double tolerance, EventListBase *destination) {
  if (!this->empty()) {
    this->sortTof();
    if (destination == this) {
            // Put results in a temp output
            std::vector<Event> out;
            compressEventsHelper(*(this->events), out, tolerance);
            // Put it back
            this->events->swap(out);
        } else {
            compressEventsHelper(*(this->events), destination->events, tolerance);
    }
    
  }
  // In all cases, you end up WEIGHTED_NOTIME.
  destination->eventType = WEIGHTED_NOTIME;
  // The sort is still valid!
  destination->order = TOF_SORT;
  // Empty out storage for vectors that are now unused.
  destination->clearUnused();
}

void EventListBase::compressFatEvents(const double tolerance, const Mantid::Types::Core::DateAndTime &timeStart,
                                  const double seconds, EventListBase *destination) {

  // only worry about non-empty EventLists
  if (!this->empty()) {
    if (destination == this) {
      // Put results in a temp output
      std::vector<Event> out;
      compressFatEventsHelper(*(this->events), out, tolerance, timeStart, seconds);
      // Put it back
      this->weightedEveventsents.swap(out);
    } else {
      compressFatEventsHelper(*(this->events), destination->events, tolerance, timeStart, seconds);
    }
  }
  // In all cases, you end up WEIGHTED_NOTIME.
  destination->eventType = WEIGHTED;
  // The sort order is pulsetimetof as we've compressed out the tolerance
  destination->order = PULSETIMETOF_SORT;
  // Empty out storage for vectors that are now unused.
  destination->clearUnused();
}


// --------------------------------------------------------------------------
/** Utility function:
 * Returns the iterator into events of the first TofEvent with
 * time at sample > seek_time
 * Will return events.end() if nothing is found!
 *
 * @param events :: event vector in which to look.
 * @param seek_time :: seek time to find (typically the first bin X[0]). Seek
 *time in nanoseconds.
 * @param tofFactor :: Time of flight factor
 * @param tofOffset :: Time of flight offset
 * @return iterator where the first event matching it is.
 */
template <class T>
typename std::vector<T>::const_iterator
EventListBase::findFirstTimeAtSampleEvent(const std::vector<T> &events, const double seek_time, const double &tofFactor,
                                      const double &tofOffset) const {
  auto itev = events.cbegin();
  auto itev_end = events.cend(); // cache for speed

  // if tof < X[0], that means that you need to skip some events
  while ((itev != itev_end) &&
         (static_cast<double>(calculateCorrectedFullTime(*itev, tofFactor, tofOffset)) < seek_time))
    itev++;
  // Better fix would be to use a binary search instead of the linear one used
  // here.
  return itev;
}



// --------------------------------------------------------------------------
/** Generates both the Y and E (error) histograms w.r.t Pulse Time
 * for an EventListBase with or without Weightedevents->
 *
 * @param X: x-bins supplied
 * @param Y: counts returned
 * @param E: errors returned
 * @param skipError: skip calculating the error. This has no effect for weighted
 *        events; you can just ignore the returned E vector.
 */
void EventListBase::generateHistogramPulseTime(const MantidVec &X, MantidVec &Y, MantidVec &E, bool skipError) const {
  throwUnimplementedError();
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
void EventListBase::generateHistogramTimeAtSample(const MantidVec &X, MantidVec &Y, MantidVec &E, const double &tofFactor,
                                              const double &tofOffset, bool skipError) const {
  throwUnimplementedError();
}

// --------------------------------------------------------------------------
/** Generates both the Y and E (error) histograms w.r.t TOF
 * for an EventListBase with or without Weightedevents->
 *
 * @param X: x-bins supplied
 * @param Y: counts returned
 * @param E: errors returned
 * @param skipError: skip calculating the error. This has no effect for weighted
 *        events; you can just ignore the returned E vector.
 */
void EventListBase::generateHistogram(const MantidVec &X, MantidVec &Y, MantidVec &E, bool skipError) const {
  // All types of weights need to be sorted by TOF
  this->sortTof();
  this->generateCountsHistogram(X, Y);
  if (!skipError)
    this->generateErrorsHistogram(Y, E);
}

// --------------------------------------------------------------------------
/** With respect to PulseTime Fill a histogram given specified histogram bounds.
 * Does not modify
 * the EventListBase (const method).
 * @param X :: The x bins
 * @param Y :: The generated counts histogram
 */
void EventListBase::generateCountsHistogramPulseTime(const MantidVec &X, MantidVec &Y) const {
  // For slight speed=up.
  size_t x_size = X.size();

  if (x_size <= 1) {
    // X was not set. Return an empty array.
    Y.resize(0, 0);
    return;
  }

  // Sort the events by pulsetime
  this->sortPulseTime();
  // Clear the Y data, assign all to 0.
  Y.resize(x_size - 1, 0);

  //---------------------- Histogram without weights
  //---------------------------------

  if (!this->events->empty()) {
    // Iterate through all events (sorted by pulse time)
    auto itev = findFirstPulseEvent(this->events, X[0]);
    auto itev_end = events->cend(); // cache for speed
    // The above can still take you to end() if no events above X[0], so check
    // again.
    if (itev == itev_end)
      return;

    // Find the first bin
    size_t bin = 0;

    // The tof is greater the first bin boundary, so we need to find the first
    // bin
    double pulsetime = static_cast<double>(itev->pulseTime().totalNanoseconds());
    while (bin < x_size - 1) {
      // Within range?
      if ((pulsetime >= X[bin]) && (pulsetime < X[bin + 1])) {
        Y[bin]++;
        break;
      }
      ++bin;
    }
    // Go to the next event, we've already binned this first one.
    ++itev;

    // Keep going through all the events
    while ((itev != itev_end) && (bin < x_size - 1)) {
      pulsetime = static_cast<double>(itev->pulseTime().totalNanoseconds());
      while (bin < x_size - 1) {
        // Within range?
        if ((pulsetime >= X[bin]) && (pulsetime < X[bin + 1])) {
          Y[bin]++;
          break;
        }
        ++bin;
      }
      ++itev;
    }
  } // end if (there are any events to histogram)
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
void EventListBase::generateCountsHistogramPulseTime(const double &xMin, const double &xMax, MantidVec &Y,
                                                 const double TOF_min, const double TOF_max) const {

  if (this->events->empty())
    return;

  size_t nBins = Y.size();

  if (nBins == 0)
    return;

  double step = (xMax - xMin) / static_cast<double>(nBins);

  for (const TofEvent &ev : this->events) {
    double pulsetime = static_cast<double>(ev.pulseTime().totalNanoseconds());
    if (pulsetime < xMin || pulsetime >= xMax)
      continue;
    if (ev.tof() < TOF_min || ev.tof() >= TOF_max)
      continue;

    auto n_bin = static_cast<size_t>((pulsetime - xMin) / step);
    Y[n_bin]++;
  }
}

// --------------------------------------------------------------------------
/** With respect to Time at Sample, fill a histogram given specified histogram
 * bounds. Does not modify
 * the EventListBase (const method).
 * @param X :: The x bins
 * @param Y :: The generated counts histogram
 * @param tofFactor :: time of flight factor
 * @param tofOffset :: time of flight offset
 */
void EventListBase::generateCountsHistogramTimeAtSample(const MantidVec &X, MantidVec &Y, const double &tofFactor,
                                                    const double &tofOffset) const {
  // For slight speed=up.
  const size_t x_size = X.size();

  if (x_size <= 1) {
    // X was not set. Return an empty array.
    Y.resize(0, 0);
    return;
  }

  // Sort the events by pulsetime
  this->sortTimeAtSample(tofFactor, tofOffset);
  // Clear the Y data, assign all to 0.
  Y.resize(x_size - 1, 0);

  //---------------------- Histogram without weights
  //---------------------------------

  if (!this->events->empty()) {
    // Iterate through all events (sorted by pulse time)
    auto itev = findFirstTimeAtSampleEvent(this->events, X[0], tofFactor, tofOffset);
    std::vector<TofEvent>::const_iterator itev_end = events->end(); // cache for speed
    // The above can still take you to end() if no events above X[0], so check
    // again.
    if (itev == itev_end)
      return;

    // Find the first bin
    size_t bin = 0;

    auto tAtSample = static_cast<double>(calculateCorrectedFullTime(*itev, tofFactor, tofOffset));
    while (bin < x_size - 1) {
      // Within range?
      if ((tAtSample >= X[bin]) && (tAtSample < X[bin + 1])) {
        Y[bin]++;
        break;
      }
      ++bin;
    }
    // Go to the next event, we've already binned this first one.
    ++itev;

    // Keep going through all the events
    while ((itev != itev_end) && (bin < x_size - 1)) {
      tAtSample = static_cast<double>(calculateCorrectedFullTime(*itev, tofFactor, tofOffset));
      while (bin < x_size - 1) {
        // Within range?
        if ((tAtSample >= X[bin]) && (tAtSample < X[bin + 1])) {
          Y[bin]++;
          break;
        }
        ++bin;
      }
      ++itev;
    }
  } // end if (there are any events to histogram)
}

// --------------------------------------------------------------------------
/** Fill a histogram given specified histogram bounds. Does not modify
 * the EventListBase (const method).
 * @param X :: The x bins
 * @param Y :: The generated counts histogram
 */
void EventListBase::generateCountsHistogram(const MantidVec &X, MantidVec &Y) const {
  // For slight speed=up.
  size_t x_size = X.size();

  if (x_size <= 1) {
    // X was not set. Return an empty array.
    Y.resize(0, 0);
    return;
  }

  // Sort the events by tof
  this->sortTof();
  // Clear the Y data, assign all to 0.
  Y.resize(x_size - 1, 0);

  //---------------------- Histogram without weights
  //---------------------------------

  // Do we even have any events to do?
  if (!this->events->empty()) {
    // Iterate through all events (sorted by tof) placing them in the correct
    // bin.
    auto itev = findFirstEvent(this->events, TofEvent(X[0]));
    // Go through all the events,
    for (auto itx = X.cbegin(); itev != events->end(); ++itev) {
      double tof = itev->tof();
      itx = std::find_if(itx, X.cend(), [tof](const double x) { return tof < x; });
      if (itx == X.cend()) {
        break;
      }
      auto bin = std::max(std::distance(X.cbegin(), itx) - 1, std::ptrdiff_t{0});
      ++Y[bin];
    }
  } // end if (there are any events to histogram)
}

// --------------------------------------------------------------------------
/**
 * Generate the Error histogram for the provided counts histogram.
 * It simply returns the sqrt of the number of counts for each bin.
 *
 * @param Y :: The counts histogram
 * @param E :: The generated error histogram
 */
void EventListBase::generateErrorsHistogram(const MantidVec &Y, MantidVec &E) const {
  // Fill the vector for the errors, containing sqrt(count)
  E.resize(Y.size(), 0);

  // windows can get confused about std::sqrt
  std::transform(Y.begin(), Y.end(), E.begin(), static_cast<double (*)(double)>(sqrt));

} //----------------------------------------------------------------------------------

// --------------------------------------------------------------------------
/** Integrate the events between a range of X values, or all events->
 *
 * @param minX :: minimum X bin to use in integrating.
 * @param maxX :: maximum X bin to use in integrating.
 * @param entireRange :: set to true to use the entire range. minX and maxX are
 *then ignored!
 * @return the integrated number of events->
 */
double EventListBase::integrate(const double minX, const double maxX, const bool entireRange) const {
  double sum(0), error(0);
  integrate(minX, maxX, entireRange, sum, error);
  return sum;
}

/** Integrate the events between a range of X values, or all events->
 *
 * @param minX :: minimum X bin to use in integrating.
 * @param maxX :: maximum X bin to use in integrating.
 * @param entireRange :: set to true to use the entire range. minX and maxX are
 *then ignored!
 * @param sum :: place holder for the resulting sum
 * @param error :: place holder for the resulting sum of errors
 * @return the integrated number of events->
 */
void EventListBase::integrate(const double minX, const double maxX, const bool entireRange, double &sum,
                          double &error) const {
  sum = 0;
  error = 0;
  if (!entireRange) {
    // The event list must be sorted by TOF!
    this->sortTof();
  }

  integrateHelper(*(this->events), minX, maxX, entireRange, sum, error);
    
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
void EventListBase::convertTof(std::function<double(double)> func, const int sorting) {
  // fix the histogram parameter
  MantidVec &x = dataX();
  transform(x.begin(), x.end(), x.begin(), func);

  // do nothing if sorting > 0
  if (sorting == 0) {
    this->setSortOrder(UNSORTED);
  } else if ((sorting < 0) && (this->getSortType() == TOF_SORT)) {
    this->reverse();
  }

  if (this->getNumberEvents() <= 0)
    return;

  this->convertTofHelper(*(this->events), func);

}


// --------------------------------------------------------------------------
/**
 * Convert the time of flight by tof'=tof*factor+offset
 * @param factor :: The value to scale the time-of-flight by
 * @param offset :: The value to shift the time-of-flight by
 */
void EventListBase::convertTof(const double factor, const double offset) {
  // fix the histogram parameter
  auto &x = mutableX();
  x *= factor;
  x += offset;

  if ((factor < 0.) && (this->getSortType() == TOF_SORT))
    this->reverse();

  if (this->getNumberEvents() <= 0)
    return;
  
  this->convertTofHelper(*(this->events), factor, offset);

}

// --------------------------------------------------------------------------
/**
 * Convert the units in the TofEvent's m_tof field to
 *  some other value, by scaling by a multiplier.
 * @param factor :: conversion factor (e.g. multiply TOF by this to get
 * d-spacing)
 */
void EventListBase::scaleTof(const double factor) { this->convertTof(factor, 0.0); }

// --------------------------------------------------------------------------
/** Add an offset to the TOF of each event in the list.
 *
 * @param offset :: The value to shift the time-of-flight by
 */
void EventListBase::addTof(const double offset) { this->convertTof(1.0, offset); }




// --------------------------------------------------------------------------
/** Add an offset to the pulsetime (wall-clock time) of each event in the list.
 *
 * @param seconds :: The value to shift the pulsetime by, in seconds
 */
void EventListBase::addPulsetime(const double seconds) {
  if (this->getNumberEvents() <= 0)
    return;

  // Convert the list
  this->addPulsetimeHelper(*(this->events), seconds);
}

// --------------------------------------------------------------------------
/** Add an offset to the pulsetime (wall-clock time) of each event in the list.
 *
 * @param seconds :: A set of values to shift the pulsetime by, in seconds
 */
void EventListBase::addPulsetimes(const std::vector<double> &seconds) {
  if (this->getNumberEvents() <= 0)
    return;
  if (this->getNumberEvents() != seconds.size()) {
    throw std::runtime_error("");
  }
  this->addPulsetimesHelper(*(this->events), seconds);
}


// --------------------------------------------------------------------------
/**
 * Mask out events that have a tof between tofMin and tofMax (inclusively).
 * Events are removed from the list.
 * @param tofMin :: lower bound of TOF to filter out
 * @param tofMax :: upper bound of TOF to filter out
 */
void EventListBase::maskTof(const double tofMin, const double tofMax) {
  if (tofMax <= tofMin)
    throw std::runtime_error("EventListBase::maskTof: tofMax must be > tofMin");

  // don't do anything with an emply list
  if (this->getNumberEvents() == 0)
    return;

  // Start by sorting by tof
  this->sortTof();

  // Convert the list
  size_t numOrig = 0;
  size_t numDel = 0;

  numOrig = this->events->size();
  numDel = this->maskTofHelper(*(this->events), tofMin, tofMax);
   

  if (numDel >= numOrig)
    this->clear(false);
}



// --------------------------------------------------------------------------
/**
 * Mask out events by the condition vector.
 * Events are removed from the list.
 * @param mask :: condition vector
 */
void EventListBase::maskCondition(const std::vector<bool> &mask) {

  // mask size must match the number of events
  if (this->getNumberEvents() != mask.size())
    throw std::runtime_error("EventListBase::maskTof: tofMax must be > tofMin");

  // don't do anything with an emply list
  if (this->getNumberEvents() == 0)
    return;

  // Convert the list
  size_t numOrig = 0;
  size_t numDel = 0;
  
  numOrig = this->events->size();
  numDel = this->maskConditionHelper(*(this->events), mask);
   

  if (numDel >= numOrig)
    this->clear(false);
}


/** Fill a vector with the list of TOFs
 *  @param tofs :: A reference to the vector to be filled
 */
void EventListBase::getTofs(std::vector<double> &tofs) const {
  // Set the capacity of the vector to avoid multiple resizes
  tofs.reserve(this->getNumberEvents());
  this->getTofsHelper(*(this->events), tofs);
}

/** Get the times-of-flight of each event in this EventListBase.
 *
 * @return by copy a vector of doubles of the tof() value
 */
std::vector<double> EventListBase::getTofs() const {
  std::vector<double> tofs;
  this->getTofs(tofs);
  return tofs;
}

/** Fill a vector with the list of Weights
 *  @param weights :: A reference to the vector to be filled
 */
void EventListBase::getWeights(std::vector<double> &weights) const {
  // Set the capacity of the vector to avoid multiple resizes
  weights.reserve(this->getNumberEvents());
    // not a weighted event type, return 1.0 for all.
  weights.assign(this->getNumberEvents(), 1.0);
}

/** Get the weight of each event in this EventListBase.
 *
 * @return by copy a vector of doubles of the weight() value
 */
std::vector<double> EventListBase::getWeights() const {
  std::vector<double> weights;
  this->getWeights(weights);
  return weights;
}

/** Fill a vector with the list of Weight Errors
 *  @param weightErrors :: A reference to the vector to be filled
 */
void EventListBase::getWeightErrors(std::vector<double> &weightErrors) const {
  // Set the capacity of the vector to avoid multiple resizes
  weightErrors.reserve(this->getNumberEvents());
  weightErrors.assign(this->getNumberEvents(), 1.0);
}

/** Get the weight error of each event in this EventListBase.
 *
 * @return by copy a vector of doubles of the weight() value
 */
std::vector<double> EventListBase::getWeightErrors() const {
  std::vector<double> weightErrors;
  this->getWeightErrors(weightErrors);
  return weightErrors;
}



/** Get the pulse times of each event in this EventListBase.
 *
 * @return by copy a vector of DateAndTime times
 */
std::vector<Mantid::Types::Core::DateAndTime> EventListBase::getPulseTimes() const {
  std::vector<Mantid::Types::Core::DateAndTime> times;
  // Set the capacity of the vector to avoid multiple resizes
  times.reserve(this->getNumberEvents());

  // Convert the list
  this->getPulseTimesHelper(*(this->events), times);

  return times;
}

// --------------------------------------------------------------------------
/**
 * @return The minimum tof value for the list of the events->
 */
double EventListBase::getTofMin() const {
  // set up as the maximum available double
  double tMin = std::numeric_limits<double>::max();

  // no events is a soft error
  if (this->empty())
    return tMin;

  // when events are ordered by tof just need the first value
  if (this->order == TOF_SORT) {
    return this->events->begin()->tof();
  }

  // now we are stuck with a linear search
  double temp = tMin; // start with the largest possible value
  size_t numEvents = this->getNumberEvents();
  for (size_t i = 0; i < numEvents; i++) {
    temp = this->events[i].tof();
    if (temp < tMin)
      tMin = temp;
  }
  return tMin;
}

/**
 * @return The maximum tof value for the list of events->
 */
double EventListBase::getTofMax() const {
  // set up as the minimum available double
  double tMax = std::numeric_limits<double>::lowest();

  // no events is a soft error
  if (this->empty())
    return tMax;

  // when events are ordered by tof just need the first value
  if (this->order == TOF_SORT) {
    return this->events->rbegin()->tof();
  }

  // now we are stuck with a linear search
  size_t numEvents = this->getNumberEvents();
  double temp = tMax; // start with the smallest possible value
  for (size_t i = 0; i < numEvents; i++) {
    temp = this->events[i].tof();
    if (temp > tMax)
      tMax = temp;
  }
  return tMax;
}

// --------------------------------------------------------------------------
/**
 * @return The minimum tof value for the list of the events->
 */
DateAndTime EventListBase::getPulseTimeMin() const {
  // set up as the maximum available date time.
  DateAndTime tMin = DateAndTime::maximum();

  // no events is a soft error
  if (this->empty())
    return tMin;

  // when events are ordered by pulse time just need the first value
  if (this->order == PULSETIME_SORT) {
    return this->events->begin()->pulseTime();
  }

  // now we are stuck with a linear search
  DateAndTime temp = tMin; // start with the largest possible value
  size_t numEvents = this->getNumberEvents();
  for (size_t i = 0; i < numEvents; i++) {
    temp = this->events[i].pulseTime();
    if (temp < tMin)
      tMin = temp;
  }
  return tMin;
}

/**
 * @return The maximum tof value for the list of events->
 */
DateAndTime EventListBase::getPulseTimeMax() const {
  // set up as the minimum available date time.
  DateAndTime tMax = DateAndTime::minimum();

  // no events is a soft error
  if (this->empty())
    return tMax;

  // when events are ordered by pulse time just need the first value
  if (this->order == PULSETIME_SORT) {
    return this->events->rbegin()->pulseTime();
  }

  // now we are stuck with a linear search
  size_t numEvents = this->getNumberEvents();
  DateAndTime temp = tMax; // start with the smallest possible value
  for (size_t i = 0; i < numEvents; i++) {
    temp = this->events[i].pulseTime();
    if (temp > tMax)
      tMax = temp;
  }
  return tMax;
}

void EventListBase::getPulseTimeMinMax(Mantid::Types::Core::DateAndTime &tMin,
                                   Mantid::Types::Core::DateAndTime &tMax) const {
  // set up as the minimum available date time.
  tMax = DateAndTime::minimum();
  tMin = DateAndTime::maximum();

  // no events is a soft error
  if (this->empty())
    return;

  // when events are ordered by pulse time just need the first/last values
  if (this->order == PULSETIME_SORT) {
    Min = this->events->begin()->pulseTime();
    tMax = this->events->rbegin()->pulseTime();
    return;
  }

  // now we are stuck with a linear search
  size_t numEvents = this->getNumberEvents();
  DateAndTime temp = tMax; // start with the smallest possible value
  for (size_t i = 0; i < numEvents; i++) {
    temp = this->events[i].pulseTime();
    
    if (temp > tMax)
      tMax = temp;
    if (temp < tMin)
      tMin = temp;
  }
}

DateAndTime EventListBase::getTimeAtSampleMax(const double &tofFactor, const double &tofOffset) const {
  // set up as the minimum available date time.
  DateAndTime tMax = DateAndTime::minimum();

  // no events is a soft error
  if (this->empty())
    return tMax;

  // when events are ordered by time at sample just need the first value
  if (this->order == TIMEATSAMPLE_SORT) {
    return calculateCorrectedFullTime(*(this->events->rbegin()), tofFactor, tofOffset);
  }

  // now we are stuck with a linear search
  size_t numEvents = this->getNumberEvents();
  DateAndTime temp = tMax; // start with the smallest possible value
  for (size_t i = 0; i < numEvents; i++) {
    temp = calculateCorrectedFullTime(this->events[i], tofFactor, tofOffset);
      
    if (temp > tMax)
      tMax = temp;
  }
  return tMax;
}

DateAndTime EventListBase::getTimeAtSampleMin(const double &tofFactor, const double &tofOffset) const {
  // set up as the minimum available date time.
  DateAndTime tMin = DateAndTime::maximum();

  // no events is a soft error
  if (this->empty())
    return tMin;

  // when events are ordered by time at sample just need the first value
  if (this->order == TIMEATSAMPLE_SORT) {
    return calculateCorrectedFullTime(*(this->events->begin()), tofFactor, tofOffset);
  }

  // now we are stuck with a linear search
  size_t numEvents = this->getNumberEvents();
  DateAndTime temp = tMin; // start with the smallest possible value
  for (size_t i = 0; i < numEvents; i++) {
    temp = calculateCorrectedFullTime(this->events[i], tofFactor, tofOffset);
     
    if (temp < tMin)
      tMin = temp;
  }
  return tMin;
}



// --------------------------------------------------------------------------
/**
 * Set a list of TOFs to the current event list. Modify the units if necessary.
 *
 * @param tofs :: The vector of doubles to set the tofs to.
 */
void EventListBase::setTofs(const MantidVec &tofs) {
  this->order = UNSORTED;

  // Convert the list
  this->setTofsHelper(*(this->events), tofs);
    
}

// ==============================================================================================
// ----------- MULTIPLY AND DIVIDE ---------------------------------------
// ==============================================================================================


//------------------------------------------------------------------------------------------------
/** Operator to multiply the weights in this EventListBase by an error-less scalar.
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
EventListBase &EventListBase::operator*=(const double value) {
  this->multiply(value);
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
void EventListBase::multiply(const double value, const double error) {
  // Do nothing if multiplying by exactly one and there is no error
  if ((value == 1.0) && (error == 0.0))
    return;

  // TODO abstract out to wrapper?
  this->switchTo(WEIGHTED);

  multiplyHelper(*(this->events), value, error);

}


//------------------------------------------------------------------------------------------------
/** Multiply the weights in this event list by a histogram.
 * The event list switches to WeightedEvent's if needed.
 * NOTE: no unit checks are made (or possible to make) to compare the units of X
 *and tof() in the EventListBase.
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
void EventListBase::multiply(const MantidVec &X, const MantidVec &Y, const MantidVec &E) {

    // Switch to weights if needed.
    // TODO abstract out to wrapper?
    this->switchTo(WEIGHTED);
    // Fall through

    this->sortTof();
    multiplyHistogramHelper(*(this->events), X, Y, E);

}

//------------------------------------------------------------------------------------------------
/** Divide the weights in this event list by a histogram.
 * The event list switches to WeightedEvent's if needed.
 * NOTE: no unit checks are made (or possible to make) to compare the units of X
 *and tof() in the EventListBase.
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
void EventListBase::divide(const MantidVec &X, const MantidVec &Y, const MantidVec &E) {
    // Switch to weights if needed.
    // TODO abstract out to wrapper?
    this->switchTo(WEIGHTED);
    // Fall through

  this->sortTof();
  divideHistogramHelper(*(this->events), X, Y, E);
    
}

//------------------------------------------------------------------------------------------------
/** Operator to divide the weights in this EventListBase by an error-less scalar.
 * Use divide(value,error) if your scalar has an error!
 * This simply calls the equivalent function: multiply(1.0/value).
 *
 * @param value :: divide by this
 * @return reference to this
 * @throw std::invalid_argument if value == 0; cannot divide by zero.
 */
EventListBase &EventListBase::operator/=(const double value) {
  if (value == 0.0)
    throw std::invalid_argument("EventListBase::divide() called with value of 0.0. Cannot divide by zero.");
  this->multiply(1.0 / value, 0.0);
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
void EventListBase::divide(const double value, const double error) {
  if (value == 0.0)
    throw std::invalid_argument("EventListBase::divide() called with value of 0.0. Cannot divide by zero.");
  // Do nothing if dividing by exactly 1.0, no error
  else if (value == 1.0 && error == 0.0)
    return;

  // We'll multiply by 1/value
  double invValue = 1.0 / value;
  // Relative error remains the same
  double invError = (error / value) * invValue;

  this->multiply(invValue, invError);
}

// ==============================================================================================
// ----------- SPLITTING AND FILTERING ---------------------------------------
// ==============================================================================================


/** Filter a vector of events into another based on time at sample.
 * TODO: Make this more efficient using STL-fu.
 * @param events :: input events
 * @param start :: start time (absolute)
 * @param stop :: end time (absolute)
 * @param tofFactor :: scaling factor for tof
 * @param tofOffset :: offset for tof
 * @param output :: reference to an event list that will be output.
 */
template <class T>
void EventListBase::filterByTimeAtSampleHelper(std::vector<T> &events, DateAndTime start, DateAndTime stop,
                                           double tofFactor, double tofOffset, std::vector<T> &output) {
  auto itev = events.begin();
  auto itev_end = events.end();
  // Find the first event with m_pulsetime >= start
  while ((itev != itev_end) && (calculateCorrectedFullTime(*itev, tofFactor, tofOffset) < start.totalNanoseconds()))
    itev++;

  while ((itev != itev_end) && (calculateCorrectedFullTime(*itev, tofFactor, tofOffset) < stop.totalNanoseconds())) {
    // Add the copy to the output
    output.emplace_back(*itev);
    ++itev;
  }
}

//------------------------------------------------------------------------------------------------
/** Filter this EventListBase into an output EventListBase, using
 * keeping only events within the >= start and < end pulse times.
 * Detector IDs and the X axis are copied as well.
 *
 * @param start :: start time (absolute)
 * @param stop :: end time (absolute)
 * @param output :: reference to an event list that will be output.
 * @throws std::invalid_argument If output is a reference to this EventListBase
 */
void EventListBase::filterByPulseTime(DateAndTime start, DateAndTime stop, EventListBase &output) const {
  if (this == &output) {
    throw std::invalid_argument("In-place filtering is not allowed");
  }

  // Start by sorting the event list by pulse time.
  this->sortPulseTime();
  // Clear the output
  output.clear();
  // Has to match the given type
  output.switchTo(eventType);
  output.setDetectorIDs(this->getDetectorIDs());
  output.setHistogram(m_histogram);
  output.setSortOrder(this->order);

  // Iterate through all events (sorted by pulse time)
  filterByPulseTimeHelper(*(this->events), start, stop, output.events);    break;
  
}

void EventListBase::filterByTimeAtSample(Types::Core::DateAndTime start, Types::Core::DateAndTime stop, double tofFactor,
                                     double tofOffset, EventListBase &output) const {
  if (this == &output) {
    throw std::invalid_argument("In-place filtering is not allowed");
  }

  // Start by sorting
  this->sortTimeAtSample(tofFactor, tofOffset);
  // Clear the output
  output.clear();
  // Has to match the given type
  output.switchTo(eventType);
  output.setDetectorIDs(this->getDetectorIDs());
  output.setHistogram(m_histogram);
  output.setSortOrder(this->order);

  filterByTimeAtSampleHelper(*(this->events), start, stop, tofFactor, tofOffset, output.events);
  
}


//------------------------------------------------------------------------------------------------
/** Use a TimeSplitterType to filter the event list in place.
 *
 * @param splitter :: a TimeSplitterType where all the entries (start/end time)
 *indicate events
 *     that will be kept. Any other events will be deleted.
 */
void EventListBase::filterInPlace(Kernel::TimeSplitterType &splitter) {
  // Start by sorting the event list by pulse time.
  this->sortPulseTime();

  // Iterate through all events (sorted by pulse time)
  filterInPlaceHelper(splitter, this->events);
   
}



//------------------------------------------------------------------------------------------------
/** Split the event list into n outputs
 *
 * @param splitter :: a TimeSplitterType giving where to split
 * @param outputs :: a vector of where the split events will end up. The # of
 *entries in there should
 *        be big enough to accommodate the indices.
 */
void EventListBase::splitByTime(Kernel::TimeSplitterType &splitter, std::vector<EventListBase *> outputs) const {


  // Start by sorting the event list by pulse time.
  this->sortPulseTime();

  // Initialize all the outputs
  size_t numOutputs = outputs.size();
  for (size_t i = 0; i < numOutputs; i++) {
    outputs[i]->clear();
    outputs[i]->setDetectorIDs(this->getDetectorIDs());
    outputs[i]->setHistogram(m_histogram);
    // Match the output event type.
    outputs[i]->switchTo(eventType);
  }

  // Do nothing if there are no entries
  if (splitter.empty())
    return;

  splitByTimeHelper(splitter, outputs, this->events);

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
void EventListBase::splitByFullTime(Kernel::TimeSplitterType &splitter, std::map<int, EventListBase *> outputs,
                                bool docorrection, double toffactor, double tofshift) const {
  // 1. Start by sorting the event list by pulse time.
  this->sortPulseTimeTOF();

  // 2. Initialize all the outputs
  std::map<int, EventListBase *>::iterator outiter;
  for (outiter = outputs.begin(); outiter != outputs.end(); ++outiter) {
    EventListBase *opeventlist = outiter->second;
    opeventlist->clear();
    opeventlist->setDetectorIDs(this->getDetectorIDs());
    opeventlist->setHistogram(m_histogram);
    // Match the output event type.
    opeventlist->switchTo(eventType);
  }

  // Do nothing if there are no entries
  if (splitter.empty()) {
    // 3A. Copy all events to group workspace = -1
    (*outputs[-1]) = (*this);
    // this->duplicate(outputs[-1]);
  } else {
    // 3B. Split
    splitByFullTimeHelper(splitter, outputs, this->events, docorrection, toffactor, tofshift);
    
  }
}


//----------------------------------------------------------------------------------------------
/**
 * @brief EventListBase::splitByFullTimeMatrixSplitter
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
std::string EventListBase::splitByFullTimeMatrixSplitter(const std::vector<int64_t> &vec_splitters_time,
                                                     const std::vector<int> &vecgroups,
                                                     std::map<int, EventListBase *> vec_outputEventList, bool docorrection,
                                                     double toffactor, double tofshift) const {

  // Start by sorting the event list by pulse time, if its flag is not set up
  // right
  sortPulseTimeTOF();

  // Initialize all the output event list
  std::map<int, EventListBase *>::iterator outiter;
  for (outiter = vec_outputEventList.begin(); outiter != vec_outputEventList.end(); ++outiter) {
    EventListBase *opeventlist = outiter->second;
    opeventlist->clear();
    opeventlist->setDetectorIDs(this->getDetectorIDs());
    opeventlist->setHistogram(m_histogram);
    // Match the output event type.
    opeventlist->switchTo(eventType);
  }

  std::string debugmessage;

  // Do nothing if there are no entries
  if (vecgroups.empty()) {
    // Copy all events to group workspace = -1
    (*vec_outputEventList[-1]) = (*this);
    // this->duplicate(outputs[-1]);
  } else {
    // Split

    // Try to find out which filtering algorithm to use by comparing number of
    // splitters and number of events
    bool sparse_splitter = vec_splitters_time.size() < this->getNumberEvents();
    if (sparse_splitter)
      debugmessage = splitByFullTimeSparseVectorSplitterHelper(vec_splitters_time, vecgroups, vec_outputEventList,
                                                                this->events, docorrection, toffactor, tofshift);
    else
      debugmessage = splitByFullTimeVectorSplitterHelper(vec_splitters_time, vecgroups, vec_outputEventList,
                                                          this->events, docorrection, toffactor, tofshift);
      
  }

  return debugmessage;
}



//----------------------------------------------------------------------------------------------
/** Split the event list by pulse time
 */
void EventListBase::splitByPulseTime(Kernel::TimeSplitterType &splitter, std::map<int, EventListBase *> outputs) const {
  // Start by sorting the event list by pulse time.
  this->sortPulseTimeTOF();

  // Initialize all the output event lists
  std::map<int, EventListBase *>::iterator outiter;
  for (outiter = outputs.begin(); outiter != outputs.end(); ++outiter) {
    EventListBase *opeventlist = outiter->second;
    opeventlist->clear();
    opeventlist->setDetectorIDs(this->getDetectorIDs());
    opeventlist->setHistogram(m_histogram);
    // Match the output event type.
    opeventlist->switchTo(eventType);
  }

  // Split
  if (splitter.empty()) {
    // No splitter: copy all events to group workspace = -1
    (*outputs[-1]) = (*this);
  } else {
    splitByPulseTimeHelper(splitter, outputs, this->events);    
  }
}

//----------------------------------------------------------------------------------------------
/** Split the event list by pulse time
 */
// TODO/NOW - TEST
void EventListBase::splitByPulseTimeWithMatrix(const std::vector<int64_t> &vec_times, const std::vector<int> &vec_target,
                                           std::map<int, EventListBase *> outputs) const {
  // Start by sorting the event list by pulse time.
  this->sortPulseTimeTOF();

  // Initialize all the output event lists
  std::map<int, EventListBase *>::iterator outiter;
  for (outiter = outputs.begin(); outiter != outputs.end(); ++outiter) {
    EventListBase *opeventlist = outiter->second;
    opeventlist->clear();
    opeventlist->setDetectorIDs(this->getDetectorIDs());
    opeventlist->setHistogram(m_histogram);
    // Match the output event type.
    opeventlist->switchTo(eventType);
  }

  // Split
  if (vec_target.empty()) {
    // No splitter: copy all events to group workspace = -1
    (*outputs[-1]) = (*this);
  } else {
    // Split
    splitByPulseTimeWithMatrixHelper(vec_times, vec_target, outputs, this->events);
  }
}


//--------------------------------------------------------------------------
/** Converts the X units in each event by going through TOF.
 * Note: if the unit conversion reverses the order, use "reverse()" to flip it
 *back.
 *
 * @param fromUnit :: the Unit describing the input unit. Must be initialized.
 * @param toUnit :: the Unit describing the output unit. Must be initialized.
 */
void EventListBase::convertUnitsViaTof(Mantid::Kernel::Unit *fromUnit, Mantid::Kernel::Unit *toUnit) {
  // Check for initialized
  if (!fromUnit || !toUnit)
    throw std::runtime_error("EventListBase::convertUnitsViaTof(): one of the units is NULL!");
  if (!fromUnit->isInitialized())
    throw std::runtime_error("EventListBase::convertUnitsViaTof(): fromUnit is not initialized!");
  if (!toUnit->isInitialized())
    throw std::runtime_error("EventListBase::convertUnitsViaTof(): toUnit is not initialized!");


  convertUnitsViaTofHelper(*(this->events), fromUnit, toUnit);
  
}

//--------------------------------------------------------------------------
/** Convert the event's TOF (x) value according to a simple output = a *
 * (input^b) relationship
 *  @param factor :: the conversion factor a to apply
 *  @param power :: the Power b to apply to the conversion
 */
void EventListBase::convertUnitsQuickly(const double &factor, const double &power) {
  auto convertEvents = *(this->events);
  convertUnitsQuicklyHelper(convertEvents, factor, power);  
}

HistogramData::Histogram &EventListBase::mutableHistogramRef() {
  if (mru)
    mru->deleteIndex(this);
  return m_histogram;
}

void EventListBase::checkAndSanitizeHistogram(HistogramData::Histogram &histogram) {
  if (histogram.xMode() != HistogramData::Histogram::XMode::BinEdges)
    throw std::runtime_error("EventListBase: setting histogram with storage mode "
                             "other than BinEdges is not possible");
  if (histogram.sharedY() || histogram.sharedE())
    throw std::runtime_error("EventListBase: setting histogram data with non-null "
                             "Y or E data is not possible");
  // Avoid flushing of YMode: we only change X but YMode depends on events.
  if (histogram.yMode() == HistogramData::Histogram::YMode::Uninitialized)
    histogram.setYMode(m_histogram.yMode());
  if (histogram.yMode() != m_histogram.yMode())
    throw std::runtime_error("EventListBase: setting histogram data with different "
                             "YMode is not possible");
}

void EventListBase::checkWorksWithPoints() const {
  throw std::runtime_error("EventListBase: setting Points as X data is not "
                           "possible, only BinEdges are supported");
}

void EventListBase::checkIsYAndEWritable() const {
  throw std::runtime_error("EventListBase: Cannot set Y or E data, these data are "
                           "generated automatically based on the events");
}

} // namespace Mantid::DataObjects
