
#pragma once

#include "MantidDataObjects/EventList.h"
#include "MantidAPI/IEventList.h"
#include "MantidDataObjects/Events.h"
#include "MantidKernel/MultiThreaded.h"
#include "MantidKernel/System.h"
#include "MantidKernel/cow_ptr.h"
#include "MantidDataObjects/CompareTimeAtSample.h"
#include <iosfwd>
#include <memory>
#include <vector>

namespace Mantid {
namespace DataObjects {
template <typename T, typename SELF>
class EventListBaseFunctionsTemplate
{
  public:

  EventListBaseFunctionsTemplate(std::shared_ptr<std::vector<T>> events) {
    this->events = events;
  }

  void checkAndSanitizeHistogram(HistogramData::Histogram &histogram) {
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

HistogramData::Histogram &mutableHistogramRef() {
  if (mru)
    mru->deleteIndex(this);
  return m_histogram;
}

  // --------------------------------------------------------------------------
/**
 * Convert the units in the TofEvent's m_tof field to
 *  some other value, by scaling by a multiplier.
 * @param factor :: conversion factor (e.g. multiply TOF by this to get
 * d-spacing)
 */
void scaleTof(const double factor) { this->convertTof(factor, 0.0); }

// --------------------------------------------------------------------------
/** Add an offset to the TOF of each event in the list.
 *
 * @param offset :: The value to shift the time-of-flight by
 */
void addTof(const double offset) { this->convertTof(1.0, offset); }


// --------------------------------------------------------------------------
/** With respect to PulseTime Fill a histogram given specified histogram bounds.
 * Does not modify
 * the EventListBase (const method).
 * @param X :: The x bins
 * @param Y :: The generated counts histogram
 */
void  generateCountsHistogramPulseTime(const MantidVec &X, MantidVec &Y) const {
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

  //NOTE:  In the original implementation this only really did stuff for TofEvents? Is that right?

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
void  generateCountsHistogramPulseTime(const double &xMin, const double &xMax, MantidVec &Y,
                                                 const double TOF_min, const double TOF_max) const {

  if (this->events->empty())
    return;

  size_t nBins = Y.size();

  if (nBins == 0)
    return;

  double step = (xMax - xMin) / static_cast<double>(nBins);

    //NOTE:  In the original implementation this only really did stuff for TofEvents? Is that right?
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
void  generateCountsHistogramTimeAtSample(const MantidVec &X, MantidVec &Y, const double &tofFactor,
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

  //NOTE:  In the original implementation this only really did stuff for TofEvents? Is that right?
}

// --------------------------------------------------------------------------
/** Fill a histogram given specified histogram bounds. Does not modify
 * the EventListBase (const method).
 * @param X :: The x bins
 * @param Y :: The generated counts histogram
 */
void  generateCountsHistogram(const MantidVec &X, MantidVec &Y) const {
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

  //NOTE:  In the original implementation this only really did stuff for TofEvents? Is that right?

}

// --------------------------------------------------------------------------
/**
 * Generate the Error histogram for the provided counts histogram.
 * It simply returns the sqrt of the number of counts for each bin.
 *
 * @param Y :: The counts histogram
 * @param E :: The generated error histogram
 */
void  generateErrorsHistogram(const MantidVec &Y, MantidVec &E) const {
  // Fill the vector for the errors, containing sqrt(count)
  E.resize(Y.size(), 0);

  // windows can get confused about std::sqrt
  std::transform(Y.begin(), Y.end(), E.begin(), static_cast<double (*)(double)>(sqrt));

} //----------------------------------------------------------------------------------


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
void generateHistogram(const MantidVec &X, MantidVec &Y, MantidVec &E, bool skipError) const {
  // All types of weights need to be sorted by TOF
  this->sortTof();
  this->generateCountsHistogram(X, Y);
  if (!skipError)
    this->generateErrorsHistogram(Y, E);
}

// ==============================================================================================
// --- Setting the Histogram X axis, without recalculating the histogram
// -----------------------
// ==============================================================================================

/** Deprecated, use setSharedX() instead. Set the x-component for the histogram
 * view. This will NOT cause the histogram to be calculated.
 * @param X :: The vector of doubles to set as the histogram limits.
 */
void  setX(const Kernel::cow_ptr<HistogramData::HistogramX> &X) {
  m_histogram.setX(X);
  if (mru)
    mru->deleteIndex(this);
}

/** Deprecated, use mutableX() instead. Returns a reference to the x data.
 *  @return a reference to the X (bin) vector.
 */
MantidVec & dataX() {
  if (mru)
    mru->deleteIndex(this);
  return m_histogram.dataX();
}

/** Deprecated, use x() instead. Returns a const reference to the x data.
 *  @return a reference to the X (bin) vector. */
const MantidVec & dataX() const { return m_histogram.dataX(); }

/// Deprecated, use x() instead. Returns the x data const
const MantidVec & readX() const { return m_histogram.readX(); }

/// Deprecated, use sharedX() instead. Returns a pointer to the x data
Kernel::cow_ptr<HistogramData::HistogramX>  ptrX() const { return m_histogram.ptrX(); }

/// Deprecated, use mutableDx() instead.
MantidVec & dataDx() { return m_histogram.dataDx(); }
/// Deprecated, use dx() instead.
const MantidVec & dataDx() const { return m_histogram.dataDx(); }
/// Deprecated, use dx() instead.
const MantidVec & readDx() const { return m_histogram.readDx(); }

// ==============================================================================================
// --- Return Data Vectors --------------------------------------------------
// ==============================================================================================

/** Calculates and returns a pointer to the Y histogrammed data.
 * Remember to delete your pointer after use!
 *
 * @return a pointer to a MantidVec
 */
MantidVec * makeDataY() const {
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
MantidVec * makeDataE() const {
  MantidVec Y;
  auto E = new MantidVec();
  generateHistogram(readX(), Y, *E);
  // Y is unused.
  return E;
}

HistogramData::Histogram  histogram() const {
  HistogramData::Histogram ret(m_histogram);
  ret.setSharedY(sharedY());
  ret.setSharedE(sharedE());
  return ret;
}

HistogramData::Counts  counts() const { return histogram().counts(); }

HistogramData::CountVariances  countVariances() const { return histogram().countVariances(); }

HistogramData::CountStandardDeviations  countStandardDeviations() const {
  return histogram().countStandardDeviations();
}

HistogramData::Frequencies  frequencies() const { return histogram().frequencies(); }

HistogramData::FrequencyVariances  frequencyVariances() const { return histogram().frequencyVariances(); }

HistogramData::FrequencyStandardDeviations  frequencyStandardDeviations() const {
  return histogram().frequencyStandardDeviations();
}

const HistogramData::HistogramY & y() const {
  if (!mru)
    throw std::runtime_error("' y()' called with no MRU set. This is not allowed.");

  return *sharedY();
}
const HistogramData::HistogramE & e() const {
  if (!mru)
    throw std::runtime_error("' e()' called with no MRU set. This is not allowed.");

  return *sharedE();
}
Kernel::cow_ptr<HistogramData::HistogramY>  sharedY() const {
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
Kernel::cow_ptr<HistogramData::HistogramE>  sharedE() const {
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
const MantidVec &dataY() const {
  if (!mru)
    throw std::runtime_error("' dataY()' called with no MRU set. This is not allowed.");

  // WARNING: The Y data of sharedY() is stored in MRU, returning reference fine
  // as long as it stays there.
  return sharedY()->rawData();
}

/** Look in the MRU to see if the E histogram has been generated before.
 * If so, return that. If not, calculate, cache and return it.
 *
 * @return reference to the E vector.
 */
const MantidVec &dataE() const {
  if (!mru)
    throw std::runtime_error("' dataE()' called with no MRU set. This is not allowed.");

  // WARNING: The E data of sharedE() is stored in MRU, returning reference fine
  // as long as it stays there.
  return sharedE()->rawData();
}


// --------------------------------------------------------------------------
/** Return the size of the histogram data.
 * @return the size of the histogram representation of the data (size of Y) **/
size_t histogram_size() const {
  size_t x_size = readX().size();
  if (x_size > 1)
    return x_size - 1;
  else
    return 0;
}

// --------------------------------------------------------------------------
/** Return true if the event list is sorted by TOF */
bool isSortedByTof() const { return (this->order == TOF_SORT); }

// --------------------------------------------------------------------------
/** Return the type of sorting used in this event list */
EventSortType getSortType() const { return this->order; }

// --------------------------------------------------------------------------
/** Reverse the histogram boundaries and the associated events if they are
 * sorted
 * by time-of-flight.
 * Does nothing if sorted otherwise or unsorted.
 * */
void reverse() {
  // reverse the histogram bin parameters
  MantidVec &x = dataX();
  std::reverse(x.begin(), x.end());

  // flip the events if they are tof sorted
  if (this->isSortedByTof()) {
    std::reverse(this->events.begin(), this->events.end());
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
size_t getNumberEvents() const {
    return this->events.size();
}

/**
 * Much like stl containers, returns true if there is nothing in the event list.
 */
bool empty() const {
  return this->events.empty();
}

// --------------------------------------------------------------------------
/** Manually set the event list sort order value. No actual sorting takes place.
 * SHOULD ONLY BE USED IN TESTS or if you know what you are doing.
 * @param order :: sort order to set.
 */
void setSortOrder(const EventSortType order) const { this->order = order; }

/** Sets the MRU list for this event list
 *
 * @param newMRU :: new MRU for the workspace containing this EventListBase
 */
void setMRU(EventWorkspaceMRU *newMRU) { mru = newMRU; }


//TODO: These series of get event type should probably be in their respective classes and not in this template

/** Return the list of WeightedEvent contained.
 * NOTE! This should be used for testing purposes only, as much as possible.
 *
 * @return a reference to the list of weighted events
 * */
std::vector<T> & getWeightedEventsNoTime() {
  if (eventType != WEIGHTED_NOTIME)
    throw std::runtime_error(" getWeightedEvents() called for an "
                             "EventListBase not of type WeightedEventNoTime. Use "
                             "getEvents() or getWeightedEvents().");
  return this->events;
}

/** Return the list of WeightedEventNoTime contained.
 * NOTE! This should be used for testing purposes only, as much as possible.
 *
 * @return a const reference to the list of weighted events
 * */
const std::vector<T> & getWeightedEventsNoTime() const {
  if (eventType != WEIGHTED_NOTIME)
    throw std::runtime_error(" getWeightedEventsNoTime() called for "
                             "an EventListBase not of type WeightedEventNoTime. "
                             "Use getEvents() or getWeightedEvents().");
  return this->events;
}


/** Return the list of WeightedEvent contained.
 * NOTE! This should be used for testing purposes only, as much as possible. The
 *EventListBase
 * may contain un-weighted events, requiring use of getEvents() instead.
 *
 * @return a reference to the list of weighted events
 * */
std::vector<T> & getWeightedEvents() {
  if (eventType != WEIGHTED)
    throw std::runtime_error(" getWeightedEvents() called for an "
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
const std::vector<T> & getWeightedEvents() const {
  if (eventType != WEIGHTED)
    throw std::runtime_error(" getWeightedEvents() called for an "
                             "EventListBase not of type WeightedEvent. Use "
                             "getEvents() or getWeightedEventsNoTime().");
  return this->events;
}

/** Return the list of TofEvents contained.
 * NOTE! This should be used for testing purposes only, as much as possible. The
 *EventListBase
 * may contain weighted events, requiring use of getWeightedEvents() instead.
 *
 * @return a reference to the list of non-weighted events
 * */
std::vector<T> & getEvents() {
  if (eventType != TOF)
    throw std::runtime_error(" getEvents() called for an EventListBase "
                             "that has weights. Use getWeightedEvents() or "
                             "getWeightedEventsNoTime().");
  return this->events;
}

/** Return the const list of TofEvents contained.
 * NOTE! This should be used for testing purposes only, as much as possible. The
 *EventListBase
 * may contain weighted events, requiring use of getWeightedEvents() instead.
 *
 * @return a const reference to the list of non-weighted events
 * */
const std::vector<T> & getEvents() const {
  if (eventType != TOF)
    throw std::runtime_error(" getEvents() called for an EventListBase "
                             "that has weights. Use getWeightedEvents() or "
                             "getWeightedEventsNoTime().");
  return this->events;
}

EventType getEventType() const { return eventType; }
  /** Inequality comparator
 * @param rhs :: other EventListBase to compare
 * @return :: true if not equal.
 */
bool operator!=(const SELF &rhs) const { return (!this->operator==(rhs)); }

  // --------------------------------------------------------------------------
/** Equality operator between EventListBase's
 * @param rhs :: other EventListBase to compare
 * @return :: true if equal.
 */
bool operator==(const SELF &rhs) const {
  return this->getNumberEvents() == rhs.getNumberEvents()
  && *events == *(rhs.events);
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
SELF &operator+=(const SELF &more_events) {
  // We'll let the += operator for the given vector of event lists handle it
  this->operator+=(more_events->events);


  // No guaranteed order
  this->order = UNSORTED;
  // Do a union between the detector IDs of both lists
  addDetectorIDs(more_events->getDetectorIDs());

  return *this;
}

// --------------------------------------------------------------------------
/** Append a list of events to the histogram.
 * The internal event list will switch to the required type.
 *
 * @param more_events :: A vector of events to append.
 * @return reference to this
 * */
SELF &operator+=(const std::vector<T> &more_events) {

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
/** Append an event to the histogram.
 * @param event :: TofEvent to add at the end of the list.
 * @return reference to this
 * */
SELF &operator+=(const T &event) {
  this->events->emplace_back(event);
  this->order = UNSORTED;
  return *this;
}

  /** Copy into this event list from another
 * @param rhs :: We will copy all the events from that into this object.
 * @return reference to this
 * */
SELF &operator=(const T &rhs) {
  // Note that we are NOT copying the MRU pointer
  // the EventWorkspace that posseses the EventListBase has already configured the mru
  
  //have wrapper perform type swap before attempting copy into
  //TODO: Check if actual issue later, why would they lazily ovewrite a typed event list with a non corresponding type?
  //this.swapTo(rhs.eventType) // will do nothing if same type
  IEventList::operator=(rhs);
  m_histogram = rhs.m_histogram;
  events = rhs.events;
  eventType = rhs.eventType;
  order = rhs.order;
  return *this;
}

/// Copy data from another EventListBase, via ISpectrum reference.
void copyDataFrom(const ISpectrum &source) { source.copyDataInto(*this); }

/// Mask the spectrum to this value. Removes all events->
void clearData() { this->clear(false); }
  /** Clear the list of events and any
 * associated detector ID's.
 * */
void clear(const bool removeDetIDs) {
  if (mru)
    mru->deleteIndex(this);
  this->events.clear();
  std::vector<T>().swap(this->events); // STL Trick to release memory
  if (removeDetIDs)
    this->clearDetectorIDs();
}


/** Reserve a certain number of entries in event list of the specified eventType
 *
 * Calls std::vector<>::reserve() in order to pre-allocate the length of the
 *event list vector.
 *
 * @param num :: number of events that will be in this EventListBase
 */
void reserve(size_t num) {
  this->events.reserve(num);
}

  // --------------------------------------------------------------------------
  /** Append an event to the histogram, without clearing the cache, to make it
   *faster.
   * NOTE: Only call this on a un-weighted event list!
   *
   * @param event :: TofEvent to add at the end of the list.
   * */
  inline void addEventQuickly(const T &event) {
    events->emplace_back(event);
    this->order = UNSORTED;
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
void compressEvents(double tolerance, EventList *destination) {
  if (!this->empty()) {
    this->sortTof();

    //in the wrapper do this before calling compress Events
    destination->clear()
    destination->switchTo(WEIGHTED_NOTIME)


    //should work because the wrapper's equal operator refers to the wrapped EventList for equality
    if (destination == this) {
            // Put results in a temp output
            std::vector<WeightedEventNoTime> out;
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

void compressFatEvents(const double tolerance, const Mantid::Types::Core::DateAndTime &timeStart,
                                  const double seconds, EventListBase *destination) {
                                      
  // only worry about non-empty EventLists
  if (!this->empty()) {
    if (destination == this) {
      // Put results in a temp output
      std::vector<WeightedEvent> out;
      compressFatEventsHelper(*(this->events), out, tolerance, timeStart, seconds);
      // Put it back
      *(this->events).swap(out);
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
/** Sort events by TOF in one thread */
void sortTof() const {
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

protected:

    //Should be a shared pointer that gets passed down 
    //to parents so they all point to the same list
    std::shared_ptr<std::vector<T>> events;


    /// Histogram object holding the histogram data. Currently only X.
    HistogramData::Histogram m_histogram;

    /// What type of event is in our list.
    Mantid::API::EventType eventType;

    /// Last sorting order
    mutable EventSortType order;

    /// MRU lists of the parent EventWorkspace
    mutable EventWorkspaceMRU *mru;

    /// Mutex that is locked while sorting an event list
    mutable std::mutex m_sortMutex;



// private:

const Histogram &histogramRef() const override { return m_histogram; }

/// Used by Histogram1D::copyDataFrom for dynamic dispatch for `other`.
void copyDataInto(Histogram1D &sink) const { sink.setHistogram(histogram()); }

/// Used by copyDataFrom for dynamic dispatch for its `source`.
void copyDataInto(T &sink) const {
  //sink.clear()
  //sink.switchTo(eventType)
  //have wrapper perform type swap before attempting copy into
  //TODO: Check if actual issue later, why would they lazily ovewrite a typed event list with a non corresponding type?
  sink.m_histogram = m_histogram;
  sink.events = events;
  sink.eventType = eventType;
  sink.order = order;
}

/** Get the weight error of each event in this EventListBase.
 *
 * @return by copy a vector of doubles of the weight() value
 */
std::vector<double> getWeightErrors() const {
  std::vector<double> weightErrors;
  this->getWeightErrors(weightErrors);
  return weightErrors;
}

/** Fill a vector with the list of Weight Errors
 *  @param weightErrors :: A reference to the vector to be filled
 */
void getWeightErrors(std::vector<double> &weightErrors) const {
  // Set the capacity of the vector to avoid multiple resizes
  weightErrors.reserve(this->events.size());
  weightErrors.assign(this->events.size(), 1.0);
}

/** Get the weight of each event in this EventListBase.
 *
 * @return by copy a vector of doubles of the weight() value
 */
std::vector<double> getWeights() const {
  std::vector<double> weights;
  this->getWeights(weights);
  return weights;
}
/** Fill a vector with the list of Weights
 *  @param weights :: A reference to the vector to be filled
 */
void getWeights(std::vector<double> &weights) const {
  // Set the capacity of the vector to avoid multiple resizes
  weights.reserve(this->events.size());
    // not a weighted event type, return 1.0 for all.
  weights.assign(this->events.size(), 1.0);
}

  // --------------------------------------------------------------------------
/** Mask out events by the condition vector.
 * Events are removed from the list.
 * @param events :: reference to a vector of events to change.
 * @param mask :: condition vector
 * @returns The number of events deleted.
 */
std::size_t maskConditionHelper(std::vector<T> &events, const std::vector<bool> &mask) {

  // runs through the two synchronized vectors and delete elements
  // for condition false
  auto itm = std::find(mask.begin(), mask.end(), false);
  auto first = events.begin() + (itm - mask.begin());

  if (itm != mask.end()) {
    for (auto ite = first; ++ite != events.end() && ++itm != mask.end();) {
      if (*itm != false) {
        *first++ = std::move(*ite);
      }
    }
  }

  auto n = events.end() - first;
  if (n != 0)
    events.erase(first, events.end());

  return n;
}

    friend T;
    friend SELF;



    // EventListBaseFunctionsTemplate() = default;

    inline T & as_underlying()
    {
        return static_cast<T&>(*this);
    }
};
}
}
