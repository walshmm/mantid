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




/// Constructor (empty)
// EventWorkspace is always histogram data and so is thus EventListBase
EventListBase::EventListBase(): m_histogram(HistogramData::Histogram::XMode::BinEdges, HistogramData::Histogram::YMode::Counts), eventType(TOF),
      order(UNSORTED), mru(nullptr)
      {}

/** Constructor with a MRU list
 * @param mru :: pointer to the MRU of the parent EventWorkspace
 * @param specNo :: the spectrum number for the event list
 */
EventListBase::EventListBase(EventWorkspaceMRU *mru, specnum_t specNo): 
      m_histogram(HistogramData::Histogram::XMode::BinEdges, HistogramData::Histogram::YMode::Counts), eventType(TOF),
      order(UNSORTED), mru(mru) {}

/** Constructor copying from an existing event list
 * @param rhs :: EventListBase object to copy*/
EventListBase::EventListBase(const EventListBase &rhs):  m_histogram(rhs.m_histogram), mru{nullptr} {

}

/** Constructor, taking a vector of events->
 * @param events :: Vector of TofEvent's */
EventListBase::EventListBase(const std::vector<TofEvent> &events) 
: m_histogram(HistogramData::Histogram::XMode::BinEdges, HistogramData::Histogram::YMode::Counts), eventType(TOF),
      mru(nullptr) {
}

/** Constructor, taking a vector of events->
 * @param events :: Vector of WeightedEvent's */
EventListBase::EventListBase(const std::vector<WeightedEvent> &events)
: m_histogram(HistogramData::Histogram::XMode::BinEdges, HistogramData::Histogram::YMode::Counts), mru(nullptr){
}

/** Constructor, taking a vector of events->
 * @param events :: Vector of WeightedEventNoTime's */
EventListBase::EventListBase(const std::vector<WeightedEventNoTime> &events)
: m_histogram(HistogramData::Histogram::XMode::BinEdges, HistogramData::Histogram::YMode::Counts), mru(nullptr) {
}

/** Constructor, taking a vector of events->
 * @param events :: Vector of TofEventNoTime's */
EventListBase::EventListBase(const std::vector<TofEventNoTime> &events) 
: m_histogram(HistogramData::Histogram::XMode::BinEdges, HistogramData::Histogram::YMode::Counts), mru(nullptr){
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

void throwUnimplementedError() {
  throw std::logic_error("Function not implemented, please use a derived class and not the base class.");
}


/// Copy data from another EventListBase, via ISpectrum reference.
void EventListBase::copyDataFrom(const ISpectrum &source) { throwUnimplementedError(); }

/// Used by copyDataFrom for dynamic dispatch for its `source`.
void EventListBase::copyDataInto(EventListBase &sink) const {
  throwUnimplementedError();
}

/// Used by Histogram1D::copyDataFrom for dynamic dispatch for `other`.
void EventListBase::copyDataInto(Histogram1D &sink) const { throwUnimplementedError(); }

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
  throwUnimplementedError();
}

// --------------------------------------------------------------------------
// --- Operators
// -------------------------------------------------------------------

/** Copy into this event list from another
 * @param rhs :: We will copy all the events from that into this object.
 * @return reference to this
 * */
EventListBase &EventListBase::operator=(const EventListBase &rhs) {
 throwUnimplementedError();
}

// --------------------------------------------------------------------------
/** Append an event to the histogram.
 * @param event :: TofEvent to add at the end of the list.
 * @return reference to this
 * */
EventListBase &EventListBase::operator+=(const TofEvent &event) {
  throwUnimplementedError();
}

// --------------------------------------------------------------------------
/** Append a list of events to the histogram.
 * The internal event list will switch to the required type.
 *
 * @param more_events :: A vector of events to append.
 * @return reference to this
 * */
EventListBase &EventListBase::operator+=(const std::vector<TofEvent> &more_events) {
  throwUnimplementedError();
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
  throwUnimplementedError();
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
  throwUnimplementedError();
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

  throwUnimplementedError();
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
  throwUnimplementedError();
}

// --------------------------------------------------------------------------
/** Equality operator between EventListBase's
 * @param rhs :: other EventListBase to compare
 * @return :: true if equal.
 */
bool EventListBase::operator==(const EventListBase &rhs) const {
  //should go to the more specific == opertator if of the same type
  return false;
}

/** Inequality comparator
 * @param rhs :: other EventListBase to compare
 * @return :: true if not equal.
 */
bool EventListBase::operator!=(const EventListBase &rhs) const { return (!this->operator==(rhs)); }

bool EventListBase::equals(const EventListBase &rhs, const double tolTof, const double tolWeight,
                       const int64_t tolPulse) const {
  
  // This is overriden in implemented classes, should not be used
  throwUnimplementedError();
}

// -----------------------------------------------------------------------------------------------
/** Return the type of Event vector contained within.
 * @return :: a EventType value.
 */
EventType EventListBase::getEventType() const { throwUnimplementedError();}

// -----------------------------------------------------------------------------------------------
/** Switch the EventListBase to use the given EventType (TOF, WEIGHTED, or
 * WEIGHTED_NOTIME)
 */


// TODO: Extract out to wrapper
void EventListBase::switchTo(EventType newType) {
  throwUnimplementedError();
}

// -----------------------------------------------------------------------------------------------
/** Switch the EventListBase to use WeightedEvents instead
 * of TofEvent.
 */
void EventListBase::switchToWeightedEvents() {
  throwUnimplementedError();
}

// -----------------------------------------------------------------------------------------------
/** Switch the EventListBase to use WeightedEventNoTime's instead
 * of TofEvent.
 */
void EventListBase::switchToWeightedEventsNoTime() {
  throwUnimplementedError();
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
  throwUnimplementedError();
}

/** Return the list of TofEvents contained.
 * NOTE! This should be used for testing purposes only, as much as possible. The
 *EventListBase
 * may contain weighted events, requiring use of getWeightedEvents() instead.
 *
 * @return a reference to the list of non-weighted events
 * */
std::vector<TofEvent> &EventListBase::getEvents() {
  throwUnimplementedError();
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
  throwUnimplementedError();
}

/** Return the list of WeightedEvent contained.
 * NOTE! This should be used for testing purposes only, as much as possible. The
 *EventListBase
 * may contain un-weighted events, requiring use of getEvents() instead.
 *
 * @return a const reference to the list of weighted events
 * */
const std::vector<WeightedEvent> &EventListBase::getWeightedEvents() const {
  throwUnimplementedError();
}

/** Return the list of WeightedEvent contained.
 * NOTE! This should be used for testing purposes only, as much as possible.
 *
 * @return a reference to the list of weighted events
 * */
std::vector<WeightedEventNoTime> &EventListBase::getWeightedEventsNoTime() {
  throwUnimplementedError();
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
  throwUnimplementedError();
}

/** Clear the list of events and any
 * associated detector ID's.
 * */
void EventListBase::clear(const bool removeDetIDs) {
  throwUnimplementedError();
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
void EventListBase::clearData() { throwUnimplementedError(); }

/** Sets the MRU list for this event list
 *
 * @param newMRU :: new MRU for the workspace containing this EventListBase
 */
void EventListBase::setMRU(EventWorkspaceMRU *newMRU) {throwUnimplementedError(); }

/** Reserve a certain number of entries in event list of the specified eventType
 *
 * Calls std::vector<>::reserve() in order to pre-allocate the length of the
 *event list vector.
 *
 * @param num :: number of events that will be in this EventListBase
 */
void EventListBase::reserve(size_t num) {
   throwUnimplementedError();
}

// ==============================================================================================
// --- Sorting functions -----------------------------------------------------
// ==============================================================================================

// --------------------------------------------------------------------------
/** Sort events by TOF or Frame
 * @param order :: Order by which to sort.
 * */
void EventListBase::sort(const EventSortType order) const {
  throwUnimplementedError();
}

// --------------------------------------------------------------------------
/** Manually set the event list sort order value. No actual sorting takes place.
 * SHOULD ONLY BE USED IN TESTS or if you know what you are doing.
 * @param order :: sort order to set.
 */
void EventListBase::setSortOrder(const EventSortType order) const {  throwUnimplementedError(); }

// --------------------------------------------------------------------------
/** Sort events by TOF in one thread */
void EventListBase::sortTof() const {
  throwUnimplementedError();
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
/** Return true if the event list is sorted by TOF */
bool EventListBase::isSortedByTof() const { throwUnimplementedError(); }

// --------------------------------------------------------------------------
/** Return the type of sorting used in this event list */
EventSortType EventListBase::getSortType() const { throwUnimplementedError(); }

// --------------------------------------------------------------------------
/** Reverse the histogram boundaries and the associated events if they are
 * sorted
 * by time-of-flight.
 * Does nothing if sorted otherwise or unsorted.
 * */
void EventListBase::reverse() {
  throwUnimplementedError();
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
    throwUnimplementedError();
}

/**
 * Much like stl containers, returns true if there is nothing in the event list.
 */
bool EventListBase::empty() const {
  throwUnimplementedError();
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
  throwUnimplementedError();
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
  throwUnimplementedError();;
}

/** Deprecated, use mutableX() instead. Returns a reference to the x data.
 *  @return a reference to the X (bin) vector.
 */
MantidVec &EventListBase::dataX() {
  throwUnimplementedError();
}

/** Deprecated, use x() instead. Returns a const reference to the x data.
 *  @return a reference to the X (bin) vector. */
const MantidVec &EventListBase::dataX() const { throwUnimplementedError(); }

/// Deprecated, use x() instead. Returns the x data const
const MantidVec &EventListBase::readX() const { throwUnimplementedError(); }

/// Deprecated, use sharedX() instead. Returns a pointer to the x data
Kernel::cow_ptr<HistogramData::HistogramX> EventListBase::ptrX() const { throwUnimplementedError(); }

/// Deprecated, use mutableDx() instead.
MantidVec &EventListBase::dataDx() { throwUnimplementedError(); }
/// Deprecated, use dx() instead.
const MantidVec &EventListBase::dataDx() const { throwUnimplementedError(); }
/// Deprecated, use dx() instead.
const MantidVec &EventListBase::readDx() const { throwUnimplementedError(); }

// ==============================================================================================
// --- Return Data Vectors --------------------------------------------------
// ==============================================================================================

/** Calculates and returns a pointer to the Y histogrammed data.
 * Remember to delete your pointer after use!
 *
 * @return a pointer to a MantidVec
 */
MantidVec *EventListBase::makeDataY() const {
  throwUnimplementedError();
}

/** Calculates and returns a pointer to the E histogrammed data.
 * Remember to delete your pointer after use!
 *
 * @return a pointer to a MantidVec
 */
MantidVec *EventListBase::makeDataE() const {
  throwUnimplementedError();
}

HistogramData::Histogram EventListBase::histogram() const {
 throwUnimplementedError();
}

HistogramData::Counts EventListBase::counts() const { throwUnimplementedError(); }

HistogramData::CountVariances EventListBase::countVariances() const { throwUnimplementedError(); }

HistogramData::CountStandardDeviations EventListBase::countStandardDeviations() const {
  throwUnimplementedError();
}

HistogramData::Frequencies EventListBase::frequencies() const { throwUnimplementedError(); }

HistogramData::FrequencyVariances EventListBase::frequencyVariances() const { throwUnimplementedError(); }

HistogramData::FrequencyStandardDeviations EventListBase::frequencyStandardDeviations() const {
  throwUnimplementedError();
}

const HistogramData::HistogramY &EventListBase::y() const {
throwUnimplementedError();
}
const HistogramData::HistogramE &EventListBase::e() const {
 throwUnimplementedError();
}
Kernel::cow_ptr<HistogramData::HistogramY> EventListBase::sharedY() const {
throwUnimplementedError();
}
Kernel::cow_ptr<HistogramData::HistogramE> EventListBase::sharedE() const {
throwUnimplementedError();
}
/** Look in the MRU to see if the Y histogram has been generated before.
 * If so, return that. If not, calculate, cache and return it.
 *
 * @return reference to the Y vector.
 */
const MantidVec &EventListBase::dataY() const {
  throwUnimplementedError();
}

/** Look in the MRU to see if the E histogram has been generated before.
 * If so, return that. If not, calculate, cache and return it.
 *
 * @return reference to the E vector.
 */
const MantidVec &EventListBase::dataE() const {
  throwUnimplementedError();
}

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

//NOTE: Dont like that this is operating with a parent type as a parameter
void EventListBase::compressEvents(double tolerance, EventList *destination) {
  throwUnimplementedError();
}

void EventListBase::compressFatEvents(const double tolerance, const Mantid::Types::Core::DateAndTime &timeStart,
                                  const double seconds, EventList *destination) {

  throwUnimplementedError();
}


// --------------------------------------------------------------------------




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
  throwUnimplementedError();
}

// --------------------------------------------------------------------------
/** With respect to PulseTime Fill a histogram given specified histogram bounds.
 * Does not modify
 * the EventListBase (const method).
 * @param X :: The x bins
 * @param Y :: The generated counts histogram
 */
void EventListBase::generateCountsHistogramPulseTime(const MantidVec &X, MantidVec &Y) const {
  throwUnimplementedError();
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
  throwUnimplementedError();
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
  throwUnimplementedError();
}

// --------------------------------------------------------------------------
/** Fill a histogram given specified histogram bounds. Does not modify
 * the EventListBase (const method).
 * @param X :: The x bins
 * @param Y :: The generated counts histogram
 */
void EventListBase::generateCountsHistogram(const MantidVec &X, MantidVec &Y) const {
  throwUnimplementedError();
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
  throwUnimplementedError();
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
  throwUnimplementedError();
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
  throwUnimplementedError();
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
  throwUnimplementedError();
}


// --------------------------------------------------------------------------
/**
 * Convert the time of flight by tof'=tof*factor+offset
 * @param factor :: The value to scale the time-of-flight by
 * @param offset :: The value to shift the time-of-flight by
 */
void convertTof(const double factor, const double offset) {
  throwUnimplementedError();
}

// --------------------------------------------------------------------------
/**
 * Convert the units in the TofEvent's m_tof field to
 *  some other value, by scaling by a multiplier.
 * @param factor :: conversion factor (e.g. multiply TOF by this to get
 * d-spacing)
 */
void EventListBase::scaleTof(const double factor) { throwUnimplementedError(); }

// --------------------------------------------------------------------------
/** Add an offset to the TOF of each event in the list.
 *
 * @param offset :: The value to shift the time-of-flight by
 */
void EventListBase::addTof(const double offset) { throwUnimplementedError(); }




// --------------------------------------------------------------------------
/** Add an offset to the pulsetime (wall-clock time) of each event in the list.
 *
 * @param seconds :: The value to shift the pulsetime by, in seconds
 */
void EventListBase::addPulsetime(const double seconds) {
  throwUnimplementedError();
}

// --------------------------------------------------------------------------
/** Add an offset to the pulsetime (wall-clock time) of each event in the list.
 *
 * @param seconds :: A set of values to shift the pulsetime by, in seconds
 */
void EventListBase::addPulsetimes(const std::vector<double> &seconds) {
  throwUnimplementedError();
}


// --------------------------------------------------------------------------
/**
 * Mask out events that have a tof between tofMin and tofMax (inclusively).
 * Events are removed from the list.
 * @param tofMin :: lower bound of TOF to filter out
 * @param tofMax :: upper bound of TOF to filter out
 */
void EventListBase::maskTof(const double tofMin, const double tofMax) {
  throwUnimplementedError();
}



// --------------------------------------------------------------------------
/**
 * Mask out events by the condition vector.
 * Events are removed from the list.
 * @param mask :: condition vector
 */
void EventListBase::maskCondition(const std::vector<bool> &mask) {
  throwUnimplementedError();
}


/** Fill a vector with the list of TOFs
 *  @param tofs :: A reference to the vector to be filled
 */
void EventListBase::getTofs(std::vector<double> &tofs) const {
  throwUnimplementedError();
}

/** Get the times-of-flight of each event in this EventListBase.
 *
 * @return by copy a vector of doubles of the tof() value
 */
std::vector<double> EventListBase::getTofs() const {
  throwUnimplementedError();
}

/** Fill a vector with the list of Weights
 *  @param weights :: A reference to the vector to be filled
 */
void EventListBase::getWeights(std::vector<double> &weights) const {
  throwUnimplementedError();
}

/** Get the weight of each event in this EventListBase.
 *
 * @return by copy a vector of doubles of the weight() value
 */
std::vector<double> EventListBase::getWeights() const {
  throwUnimplementedError();
}

/** Fill a vector with the list of Weight Errors
 *  @param weightErrors :: A reference to the vector to be filled
 */
void EventListBase::getWeightErrors(std::vector<double> &weightErrors) const {
  throwUnimplementedError();
}

/** Get the weight error of each event in this EventListBase.
 *
 * @return by copy a vector of doubles of the weight() value
 */
std::vector<double> EventListBase::getWeightErrors() const {
  throwUnimplementedError();
}



/** Get the pulse times of each event in this EventListBase.
 *
 * @return by copy a vector of DateAndTime times
 */
std::vector<Mantid::Types::Core::DateAndTime> EventListBase::getPulseTimes() const {
  throwUnimplementedError();
}

// --------------------------------------------------------------------------
/**
 * @return The minimum tof value for the list of the events->
 */
double EventListBase::getTofMin() const {
  throwUnimplementedError();
}

/**
 * @return The maximum tof value for the list of events->
 */
double EventListBase::getTofMax() const {
  throwUnimplementedError();
}

// --------------------------------------------------------------------------
/**
 * @return The minimum tof value for the list of the events->
 */
DateAndTime EventListBase::getPulseTimeMin() const {
  throwUnimplementedError();
}

/**
 * @return The maximum tof value for the list of events->
 */
DateAndTime EventListBase::getPulseTimeMax() const {
  throwUnimplementedError();
}

void EventListBase::getPulseTimeMinMax(Mantid::Types::Core::DateAndTime &tMin,
                                   Mantid::Types::Core::DateAndTime &tMax) const {
  throwUnimplementedError();
}

DateAndTime EventListBase::getTimeAtSampleMax(const double &tofFactor, const double &tofOffset) const {
  throwUnimplementedError();
}

DateAndTime EventListBase::getTimeAtSampleMin(const double &tofFactor, const double &tofOffset) const {
  throwUnimplementedError();
}



// --------------------------------------------------------------------------
/**
 * Set a list of TOFs to the current event list. Modify the units if necessary.
 *
 * @param tofs :: The vector of doubles to set the tofs to.
 */
void EventListBase::setTofs(const MantidVec &tofs) {
  throwUnimplementedError();   
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
  throwUnimplementedError();  
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
  throwUnimplementedError();  
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
  throwUnimplementedError();
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
   throwUnimplementedError();  
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
  throwUnimplementedError();
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
  throwUnimplementedError();
}

// ==============================================================================================
// ----------- SPLITTING AND FILTERING ---------------------------------------
// ==============================================================================================




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
void EventListBase::filterByPulseTime(DateAndTime start, DateAndTime stop, EventList &output) const {
  throwUnimplementedError();
}

void EventListBase::filterByTimeAtSample(Types::Core::DateAndTime start, Types::Core::DateAndTime stop, double tofFactor,
                                     double tofOffset, EventList &output) const {
  throwUnimplementedError();
}


//------------------------------------------------------------------------------------------------
/** Use a TimeSplitterType to filter the event list in place.
 *
 * @param splitter :: a TimeSplitterType where all the entries (start/end time)
 *indicate events
 *     that will be kept. Any other events will be deleted.
 */
void EventListBase::filterInPlace(Kernel::TimeSplitterType &splitter) {
  throwUnimplementedError();
}



//------------------------------------------------------------------------------------------------
/** Split the event list into n outputs
 *
 * @param splitter :: a TimeSplitterType giving where to split
 * @param outputs :: a vector of where the split events will end up. The # of
 *entries in there should
 *        be big enough to accommodate the indices.
 */
void EventListBase::splitByTime(Kernel::TimeSplitterType &splitter, std::vector<EventList *> outputs) const {
  throwUnimplementedError();
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
void EventListBase::splitByFullTime(Kernel::TimeSplitterType &splitter, std::map<int, EventList *> outputs,
                                bool docorrection, double toffactor, double tofshift) const {
  throwUnimplementedError();
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
                                                     std::map<int, EventList *> vec_outputEventList, bool docorrection,
                                                     double toffactor, double tofshift) const {
  throwUnimplementedError();
}



//----------------------------------------------------------------------------------------------
/** Split the event list by pulse time
 */
void EventListBase::splitByPulseTime(Kernel::TimeSplitterType &splitter, std::map<int, EventList *> outputs) const {
  throwUnimplementedError();
}

//----------------------------------------------------------------------------------------------
/** Split the event list by pulse time
 */
// TODO/NOW - TEST
void EventListBase::splitByPulseTimeWithMatrix(const std::vector<int64_t> &vec_times, const std::vector<int> &vec_target,
                                           std::map<int, EventList *> outputs) const {
  throwUnimplementedError();
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
  throwUnimplementedError();
}

//--------------------------------------------------------------------------
/** Convert the event's TOF (x) value according to a simple output = a *
 * (input^b) relationship
 *  @param factor :: the conversion factor a to apply
 *  @param power :: the Power b to apply to the conversion
 */
void EventListBase::convertUnitsQuickly(const double &factor, const double &power) {
  throwUnimplementedError();
}

HistogramData::Histogram &EventListBase::mutableHistogramRef() {
  throwUnimplementedError();
}

void EventListBase::checkAndSanitizeHistogram(HistogramData::Histogram &histogram) {
  throwUnimplementedError();
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
