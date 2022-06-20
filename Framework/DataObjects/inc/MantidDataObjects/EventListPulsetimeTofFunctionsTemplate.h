#include "MantidDataObjects/EventListPulsetimeFunctionsTemplate.h"
#include "MantidDataObjects/EventListTofFunctionsTemplate.h"

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

namespace Mantid {
namespace DataObjects {
template <typename T, typename SELF>
class EventListPulsetimeTofFunctionsTemplate : public EventListPulsetimeFunctionsTemplate<T, SELF>, public EventListTofFunctionsTemplate<T, SELF>
{
public:
// --------------------------------------------------------------------------
/** Sort events by TOF or Frame
 * @param order :: Order by which to sort.
 * */
void sort(const EventSortType order) const {
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

private:

/**
 * Sort by the pulse time with a tolerance. The pulsetime to compare is a
 * constant binning of seconds from start. This will set the sort order to
 * UNSORTED upon completion rather than storing the call parameters.
 * @param start The absolute start time
 * @param seconds The tolerance of pulse time in seconds.
 */
void sortPulseTimeTOFDelta(const Types::Core::DateAndTime &start, const double seconds) const {
  // Avoid sorting from multiple threads
  std::lock_guard<std::mutex> _lock(m_sortMutex);

  std::function<bool(const TofEvent &, const TofEvent &)> comparator = comparePulseTimeTOFDelta(start, seconds);

  tbb::parallel_sort(events->begin(), events->end(), comparator);
  
  this->order = UNSORTED; // so the function always re-runs
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
void filterByPulseTime(DateAndTime start, DateAndTime stop, EventListBase &output) const {
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

void filterByTimeAtSample(Types::Core::DateAndTime start, Types::Core::DateAndTime stop, double tofFactor,
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


/** Filter a vector of events into another based on time at sample.
 * TODO: Make this more efficient using STL-fu.
 * @param events :: input events
 * @param start :: start time (absolute)
 * @param stop :: end time (absolute)
 * @param tofFactor :: scaling factor for tof
 * @param tofOffset :: offset for tof
 * @param output :: reference to an event list that will be output.
 */
void filterByTimeAtSampleHelper(std::vector<T> &events, DateAndTime start, DateAndTime stop,
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
typename std::vector<T>::const_iterator
findFirstTimeAtSampleEvent(const std::vector<T> &events, const double seek_time, const double &tofFactor,
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

DateAndTime getTimeAtSampleMin(const double &tofFactor, const double &tofOffset) const {
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
  size_t numEvents = this->events.size();
  DateAndTime temp = tMin; // start with the smallest possible value
  for (size_t i = 0; i < numEvents; i++) {
    temp = calculateCorrectedFullTime(this->events[i], tofFactor, tofOffset);
     
    if (temp < tMin)
      tMin = temp;
  }
  return tMin;
}


DateAndTime getTimeAtSampleMax(const double &tofFactor, const double &tofOffset) const {
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
  size_t numEvents = this->events.size();
  DateAndTime temp = tMax; // start with the smallest possible value
  for (size_t i = 0; i < numEvents; i++) {
    temp = calculateCorrectedFullTime(this->events[i], tofFactor, tofOffset);
      
    if (temp > tMax)
      tMax = temp;
  }
  return tMax;
}


        /**
     * Calculate the corrected full time in nanoseconds
     * @param event : The event with pulse time and time-of-flight
     * @param tofFactor : Time of flight coefficient factor
     * @param tofShift : Tof shift in seconds
     * @return Corrected full time at sample in Nanoseconds.
     */
    int64_t calculateCorrectedFullTime(const T &event, const double tofFactor, const double tofShift) {
    return event.pulseTime().totalNanoseconds() +
            static_cast<int64_t>(tofFactor * (event.tof() * 1.0E3) + (tofShift * 1.0E9));
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
void splitByFullTime(Kernel::TimeSplitterType &splitter, std::map<int, EventListBase *> outputs,
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


//------------------------------------------------------------------------------------------------
/** Split the event list into n outputs, operating on a vector of either
 *TofEvent's or WeightedEvent's
 *  The comparison between neutron event and splitter is based on neutron
 *event's pulse time plus
 *
 * @param splitter :: a TimeSplitterType giving where to split
 * @param outputs :: a vector of where the split events will end up. The # of
 *entries in there should
 *        be big enough to accommodate the indices.
 * @param events :: either this->events or this->weightedevents.
 * @param docorrection :: flag to determine whether or not to apply correction
 * @param toffactor :: factor to correct TOF in formula toffactor*tof+tofshift
 * @param tofshift :: amount to shift (in SECOND) to correct TOF in formula:
 *toffactor*tof+tofshift
 */
void splitByFullTimeHelper(Kernel::TimeSplitterType &splitter, std::map<int, EventListBase *> outputs,
                                      typename std::vector<T> &events, bool docorrection, double toffactor,
                                      double tofshift) const {
  // 1. Prepare to Iterate through the splitter at the same time
  auto itspl = splitter.begin();
  auto itspl_end = splitter.end();

  // 2. Prepare to Iterate through all events (sorted by tof)
  auto itev = events.begin();
  auto itev_end = events.end();

  // 3. This is the time of the first section. Anything before is thrown out.
  while (itspl != itspl_end) {
    // Get the splitting interval times and destination
    int64_t start = itspl->start().totalNanoseconds();
    int64_t stop = itspl->stop().totalNanoseconds();
    const int index = itspl->index();

    // a) Skip the events before the start of the time
    // TODO This step can be
    EventListBase *myOutput = outputs[-1];
    while (itev != itev_end) {
      int64_t fulltime;
      if (docorrection)
        fulltime = calculateCorrectedFullTime(*itev, toffactor, tofshift);
      else
        fulltime = itev->m_pulsetime.totalNanoseconds() + static_cast<int64_t>(itev->m_tof * 1000);
      if (fulltime < start) {
        // a1) Record to index = -1 space
        const T eventCopy(*itev);
        myOutput->addEventQuickly(eventCopy);
        itev++;
      } else {
        break;
      }
    }

    // b) Go through all the events that are in the interval (if any)
    while (itev != itev_end) {
      int64_t fulltime;
      if (docorrection)
        fulltime = itev->m_pulsetime.totalNanoseconds() +
                   static_cast<int64_t>(toffactor * itev->m_tof * 1000 + tofshift * 1.0E9);
      else
        fulltime = itev->m_pulsetime.totalNanoseconds() + static_cast<int64_t>(itev->m_tof * 1000);
      if (fulltime < stop) {
        // b1) Add a copy to the output
        outputs[index]->addEventQuickly(*itev);
        ++itev;
      } else {
        break;
      }
    }

    // Go to the next interval
    ++itspl;
    // But if we reached the end, then we are done.
    if (itspl == itspl_end)
      break;

    // No need to keep looping through the filter if we are out of events
    if (itev == itev_end)
      break;
  } // END-WHILE Splitter
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
std::string splitByFullTimeMatrixSplitter(const std::vector<int64_t> &vec_splitters_time,
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

//------------------------------------------------------------------------------------------------
/** Split the event list into n outputs, operating on a vector of either
 *TofEvent's or WeightedEvent's
 *  The comparison between neutron event and splitter is based on neutron
 *event's pulse time plus
 *
 * @param vectimes :: a vector of absolute time in nanoseconds serving as
 *boundaries of splitters
 * @param vecgroups :: a vector of integer serving as the target workspace group
 *for splitters
 * @param outputs :: a vector of where the split events will end up. The # of
 *entries in there should
 *        be big enough to accommodate the indices.
 * @param vecEvents :: either this->events or this->weightedevents.
 * @param docorrection :: flag to determine whether or not to apply correction
 * @param toffactor :: factor multiplied to TOF for correcting event time from
 *detector to sample
 * @param tofshift :: shift in SECOND to TOF for correcting event time from
 *detector to sample
 */
std::string splitByFullTimeVectorSplitterHelper(const std::vector<int64_t> &vectimes, const std::vector<int> &vecgroups,
                                               std::map<int, EventListBase *> outputs, typename std::vector<T> &vecEvents,
                                               bool docorrection, double toffactor, double tofshift) const {
  // Define variables for events
  // size_t numevents = events.size();
  typename std::vector<T>::iterator eviter;
  std::stringstream msgss;

  // Loop through events
  for (eviter = vecevents.begin(); eviter != vecevents.end(); ++eviter) {
    // Obtain time of event
    int64_t evabstimens;
    if (docorrection)
      evabstimens = eviter->m_pulsetime.totalNanoseconds() +
                    static_cast<int64_t>(toffactor * eviter->m_tof * 1000 + tofshift * 1.0E9);
    else
      evabstimens = eviter->m_pulsetime.totalNanoseconds() + static_cast<int64_t>(eviter->m_tof * 1000);

    // Search in vector
    int index = static_cast<int>(lower_bound(vectimes.begin(), vectimes.end(), evabstimens) - vectimes.begin());
    int group;
    // FIXME - whether lower_bound() equal to vectimes.size()-1 should be
    // filtered out?
    if (index == 0 || index > static_cast<int>(vectimes.size() - 1)) {
      // Event is before first splitter or after last splitter.  Put to -1
      group = -1;
    } else {
      group = vecgroups[index - 1];
    }

    // Copy event to the proper group
    EventListBase *myOutput = outputs[group];
    if (!myOutput) {
      std::stringstream errss;
      errss << "Group " << group << " has a NULL output EventListBase. "
            << "\n";
      msgss << errss.str();
    } else {
      const T eventCopy(*eviter);
      myOutput->addEventQuickly(eventCopy);
    }
  }

  return (msgss.str());
}

//------------------------------------------------------------------------------------------------
/** Split the event list into n outputs, operating on a vector of either
 *TofEvent's or WeightedEvent's
 *  The comparison between neutron event and splitter is based on neutron
 *event's pulse time plus
 *
 * @param vectimes :: a vector of absolute time in nanoseconds serving as
 *boundaries of splitters
 * @param vecgroups :: a vector of integer serving as the target workspace group
 *for splitters
 * @param outputs :: a vector of where the split events will end up. The # of
 *entries in there should
 *        be big enough to accommodate the indices.
 * @param vecEvents :: either this->events or this->weightedevents.
 * @param docorrection :: flag to determine whether or not to apply correction
 * @param toffactor :: factor multiplied to TOF for correcting event time from
 *detector to sample
 * @param tofshift :: shift in SECOND to TOF for correcting event time from
 *detector to sample
 */
std::string splitByFullTimeSparseVectorSplitterHelper(const std::vector<int64_t> &vectimes,
                                                                 const std::vector<int> &vecgroups,
                                                                 std::map<int, EventListBase *> outputs,
                                                                 typename std::vector<T> &vecEvents, bool docorrection,
                                                                 double toffactor, double tofshift) const {
  // Define variables for events
  // size_t numevents = events.size();
  // typename std::vector<T>::iterator eviter;
  std::stringstream msgss;

  size_t num_splitters = vecgroups.size();
  // prepare to Iterate through all events (sorted by tof)
  auto iter_events = vecevents.begin();
  auto iter_events_end = vecevents.end();

  // std::stringstream debug_ss;
  // debug_ss << "\nFilter events...:\n";

  for (size_t i = 0; i < num_splitters; ++i) {
    // get one splitter
    int64_t start_i64 = vectimes[i];
    int64_t stop_i64 = vectimes[i + 1];
    int group = vecgroups[i];
    // debug_ss << "working on splitter: " << i << " from " << start_i64 << " to
    // " << stop_i64 << "\n";

    // go over events
    while (iter_events != iter_events_end) {
      int64_t absolute_time;
      if (docorrection)
        absolute_time = iter_events->m_pulsetime.totalNanoseconds() +
                        static_cast<int64_t>(toffactor * iter_events->m_tof * 1000 + tofshift * 1.0E9);
      else
        absolute_time = iter_events->m_pulsetime.totalNanoseconds() + static_cast<int64_t>(iter_events->m_tof * 1000);

      // debug_ss << "  event " << iter_events - vecevents->begin() << " abs.time
      // = " << absolute_time << "\n";

      if (absolute_time < start_i64) {
        // event occurs before the splitter. only can happen with first
        // splitter. Then ignore and move to next
        ++iter_events;
        continue;
      }

      if (absolute_time < stop_i64) {
        // in the splitter, then copy the event into another
        const T eventCopy(*iter_events);
        // Copy event to the proper group
        EventListBase *myOutput = outputs[group];
        if (!myOutput) {
          // there is no such group defined. quit for this group
          std::stringstream errss;
          errss << "Group " << group << " has a NULL output EventListBase. "
                << "\n";
          msgss << errss.str();
          throw std::runtime_error(errss.str());
        }
        // Add the copy to the output
        myOutput->addEventQuickly(eventCopy);
        ++iter_events;
      } else {
        // event occurs after the stop time, it should belonged to the next
        // splitter
        break;
      }
    } // while

    // quit the loop if there is no more event left
    if (iter_events == iter_events_end)
      break;
  } // for splitter

  // std::cout << debug_ss.str();

  return (msgss.str());
}

    friend T;
    EventListPulsetimeTofFunctionsTemplate() = default;

    inline T & as_underlying()
    {
        return static_cast<T&>(*this);
    }
};
}
}