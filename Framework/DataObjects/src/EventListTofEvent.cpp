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

// --------------------------------------------------------------------------
/** SUBTRACT another EventListBase from this event list.
 * The event lists are concatenated, but the weights of the incoming
 *    list are multiplied by -1.0.
 *
 * @param more_events :: Another EventListBase.
 * @return reference to this
 * */
EventListTofEvent &EventListTofEvent::operator-=(const EventListBase &more_events) {
  throw std::logic_error("EventListTofEvent's must be converted to EventListWeightedEvent, how did you get here?");
}


// NOTE/TODO: May need to move this stuff to a WeightlessTemplate so that more than just TofEvent Benefit from it


// --------------------------------------------------------------------------
/** With respect to PulseTime Fill a histogram given specified histogram bounds.
 * Does not modify
 * the EventListBase (const method).
 * @param X :: The x bins
 * @param Y :: The generated counts histogram
 */
void EventListTofEvent::generateCountsHistogramPulseTime(const MantidVec &X, MantidVec &Y) const {
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
void EventListTofEvent::generateCountsHistogramPulseTime(const double &xMin, const double &xMax, MantidVec &Y,
                                                 const double TOF_min, const double TOF_max) const {

  if (this->events->empty())
    return;

  size_t nBins = Y.size();

  if (nBins == 0)
    return;

  double step = (xMax - xMin) / static_cast<double>(nBins);

    //NOTE:  In the original implementation this only really did stuff for TofEvents? Is that right?


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
void EventListTofEvent::generateCountsHistogramTimeAtSample(const MantidVec &X, MantidVec &Y, const double &tofFactor,
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
void EventListTofEvent::generateCountsHistogram(const MantidVec &X, MantidVec &Y) const {
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



}