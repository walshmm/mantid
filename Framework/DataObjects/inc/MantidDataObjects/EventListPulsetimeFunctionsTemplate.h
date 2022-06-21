#pragma once

#include "MantidDataObjects/EventListBaseFunctionsTemplate.h"

namespace Mantid {
namespace DataObjects {
template <typename T, typename SELF>
class EventListPulsetimeFunctionsTemplate : public EventListBaseFunctionsTemplate<T, SELF>
{
public:

  EventListPulsetimeFunctionsTemplate(std::shared_ptr<std::vector<T>> events): 
  EventListBaseFunctionsTemplate<T, SELF>(events){}


void sort(const EventSortType order) const {
  if (order == UNSORTED) {
    return; // don't bother doing anything. Why did you ask to unsort?
  } else if (order == PULSETIME_SORT) {
    this->sortPulseTime();
  } else {
    throw runtime_error("Invalid sort type in EventListPulsetimeFunctionsTemplate::sort(EventSortType)");
  }
}
// --------------------------------------------------------------------------
/** Sort events by Frame */
void sortPulseTime() const {
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

protected:
/** Compare two events' FRAME id, return true if e1 should be before e2.
 * @param e1 :: first event
 * @param e2 :: second event
 *  */
bool compareEventPulseTime(const TofEvent &e1, const TofEvent &e2) { return (e1.pulseTime() < e2.pulseTime()); }


    // --------------------------------------------------------------------------
/** Utility function:
 * Returns the iterator into events of the first TofEvent with
 * pulsetime() > seek_pulsetime
 * Will return events.end() if nothing is found!
 *
 * @param events :: event vector in which to look.
 * @param seek_pulsetime :: pulse time to find (typically the first bin X[0])
 * @return iterator where the first event matching it is.
 */
typename std::vector<T>::const_iterator findFirstPulseEvent(const std::vector<T> &events,
                                                                       const double seek_pulsetime) {
  auto itev = events.begin();
  auto itev_end = events.end(); // cache for speed

  // if tof < X[0], that means that you need to skip some events
  while ((itev != itev_end) && (static_cast<double>(itev->pulseTime().totalNanoseconds()) < seek_pulsetime))
    itev++;
  // Better fix would be to use a binary search instead of the linear one used
  // here.
  return itev;
}

// --------------------------------------------------------------------------
/** Add an offset to the pulsetime (wall-clock time) of each event in the list.
 *
 * @param seconds :: A set of values to shift the pulsetime by, in seconds
 */
void addPulsetimes(const std::vector<double> &seconds) {
  if (this->getNumberEvents() <= 0)
    return;
  if (this->getNumberEvents() != seconds.size()) {
    throw std::runtime_error("");
  }
  this->addPulsetimesHelper(*(this->events), seconds);
}

// --------------------------------------------------------------------------
/** Add an offset to the pulsetime (wall-clock time) of each event in the list.
 *
 * @param seconds :: The value to shift the pulsetime by, in seconds
 */
void addPulsetime(const double seconds) {
  if (this->getNumberEvents() <= 0)
    return;

  // Convert the list
  this->addPulsetimeHelper(*(this->events), seconds);
}

// --------------------------------------------------------------------------
/** Add an offset to the pulsetime (wall-clock time) of each event in the list.
 *
 * @param events :: reference to a vector of events to change.
 * @param seconds :: The value to shift the pulsetime by, in seconds
 */
void addPulsetimeHelper(std::vector<T> &events, const double seconds) {
  // iterate through all events
  for (auto &event : events) {
    event.m_pulsetime += seconds;
  }
}


/** Add an offset per event to the pulsetime (wall-clock time) of each event in
 * the list. It is assumed that the vector sizes match.
 *
 * @param events :: reference to a vector of events to change.
 * @param seconds :: The set of values to shift the pulsetime by, in seconds
 */
void addPulsetimesHelper(std::vector<T> &events, const std::vector<double> &seconds) {
  auto eventIterEnd{events->end()};
  auto secondsIter{seconds.cbegin()};
  for (auto eventIter = events->begin(); eventIter < eventIterEnd; ++eventIter, ++secondsIter) {
    eventIter->m_pulsetime += *secondsIter;
  }
}

void getPulseTimeMinMax(Mantid::Types::Core::DateAndTime &tMin,
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
  size_t numEvents = this->events.size();
  DateAndTime temp = tMax; // start with the smallest possible value
  for (size_t i = 0; i < numEvents; i++) {
    temp = this->events[i].pulseTime();
    
    if (temp > tMax)
      tMax = temp;
    if (temp < tMin)
      tMin = temp;
  }
}


/**
 * @return The maximum tof value for the list of events->
 */
DateAndTime getPulseTimeMax() const {
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
  size_t numEvents = this->events.size();
  DateAndTime temp = tMax; // start with the smallest possible value
  for (size_t i = 0; i < numEvents; i++) {
    temp = this->events[i].pulseTime();
    if (temp > tMax)
      tMax = temp;
  }
  return tMax;
}

// --------------------------------------------------------------------------
/**
 * @return The minimum tof value for the list of the events->
 */
DateAndTime getPulseTimeMin() const {
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
  size_t numEvents = this->events.size()();
  for (size_t i = 0; i < numEvents; i++) {
    temp = this->events[i].pulseTime();
    if (temp < tMin)
      tMin = temp;
  }
  return tMin;
}

/** Get the pulse times of each event in this EventListBase.
 *
 * @return by copy a vector of DateAndTime times
 */
std::vector<Mantid::Types::Core::DateAndTime> getPulseTimes() const {
  std::vector<Mantid::Types::Core::DateAndTime> times;
  // Set the capacity of the vector to avoid multiple resizes
  times.reserve(this->events.size());

  // Convert the list
  this->getPulseTimesHelper(*(this->events), times);

  return times;
}

// --------------------------------------------------------------------------
/** Get the pulsetimes of all events in a list
 *
 * @param events :: source vector of events
 * @param times :: vector to fill
 */
void getPulseTimesHelper(const std::vector<T> &events,
                                    std::vector<Mantid::Types::Core::DateAndTime> &times) {
  times.clear();
  times.reserve(events.size());
  std::transform(events.cbegin(), events.cend(), std::back_inserter(times),
                 [](const auto &event) { return event.pulseTime(); });
}

/** Filter a vector of events into another based on pulse time.
 * TODO: Make this more efficient using STL-fu.
 * @param events :: input events
 * @param start :: start time (absolute)
 * @param stop :: end time (absolute)
 * @param output :: reference to an event list that will be output.
 */
void filterByPulseTimeHelper(std::vector<T> &events, DateAndTime start, DateAndTime stop,
                                        std::vector<T> &output) {
  auto itev = events.begin();
  auto itev_end = events.end();
  // Find the first event with m_pulsetime >= start
  while ((itev != itev_end) && (itev->m_pulsetime < start))
    itev++;

  while ((itev != itev_end) && (itev->m_pulsetime < stop)) {
    // Add the copy to the output
    output.emplace_back(*itev);
    ++itev;
  }
}

//------------------------------------------------------------------------------------------------
/** Use a TimeSplitterType to filter the event list in place.
 *
 * @param splitter :: a TimeSplitterType where all the entries (start/end time)
 *indicate events
 *     that will be kept. Any other events will be deleted.
 */
void filterInPlace(Kernel::TimeSplitterType &splitter) {
  // Start by sorting the event list by pulse time.
  this->sortPulseTime();

  // Iterate through all events (sorted by pulse time)
  filterInPlaceHelper(splitter, this->events);
   
}

//------------------------------------------------------------------------------------------------
/** Perform an in-place filtering on a vector of either TofEvent's or
 *WeightedEvent's
 *
 * @param splitter :: a TimeSplitterType where all the entries (start/end time)
 *indicate events
 *     that will be kept. Any other events will be deleted.
 * @param events :: either this->events or this->weightedevents.
 */
void filterInPlaceHelper(Kernel::TimeSplitterType &splitter, typename std::vector<T> &events) {
  // Iterate through the splitter at the same time
  auto itspl = splitter.begin();
  auto itspl_end = splitter.end();
  DateAndTime start, stop;

  // Iterate for the input
  auto itev = events.begin();
  auto itev_end = events.end();

  // Iterator for the outputted list; will follow the input except when events
  // are dropped.
  auto itOut = events.begin();

  // This is the time of the first section. Anything before is thrown out.
  while (itspl != itspl_end) {
    // Get the splitting interval times and destination
    start = itspl->start();
    stop = itspl->stop();
    const int index = itspl->index();

    // Skip the events before the start of the time
    while ((itev != itev_end) && (itev->m_pulsetime < start))
      itev++;

    // Are we aligned in the input vs output?
    bool copyingInPlace = (itOut == itev);
    if (copyingInPlace) {
      while ((itev != itev_end) && (itev->m_pulsetime < stop))
        ++itev;
      // Make sure the iterators still match
      itOut = itev;
    } else {
      // Go through all the events that are in the interval (if any)
      while ((itev != itev_end) && (itev->m_pulsetime < stop)) {
        if (index >= 0) {
          // Copy the input Event to the output iterator position.
          // Strictly speaking, this is not necessary if itOut == itev; but the
          // extra check would likely
          //  slow down the filtering in the 99% of cases where itOut != itev.
          *itOut = *itev;
          // And go up a spot in the output iterator.
          ++itOut;
        }
        ++itev;
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

  } // Looping through entries in the splitter vector

  // Ok, now resize the event list to reflect the fact that it (probably) shrank
  events.resize((itOut - events.begin()));
}

//------------------------------------------------------------------------------------------------
/** Split the event list into n outputs
 *
 * @param splitter :: a TimeSplitterType giving where to split
 * @param outputs :: a vector of where the split events will end up. The # of
 *entries in there should
 *        be big enough to accommodate the indices.
 */
void splitByTime(Kernel::TimeSplitterType &splitter, std::vector<EventListBase *> outputs) const {


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
/** Split the event list into n outputs, operating on a vector of either
 *TofEvent's or WeightedEvent's
 *  Only event's pulse time is used to compare with splitters.
 *  It is a faster and simple version of splitByFullTimeHelper
 *
 * @param splitter :: a TimeSplitterType giving where to split
 * @param outputs :: a vector of where the split events will end up. The # of
 *entries in there should
 *        be big enough to accommodate the indices.
 * @param events :: either this->events or this->weightedevents.
 */
void splitByTimeHelper(Kernel::TimeSplitterType &splitter, std::vector<EventListBase *> outputs,
                                  typename std::vector<T> &events) const {
  size_t numOutputs = outputs.size();

  // Iterate through the splitter at the same time
  auto itspl = splitter.begin();
  auto itspl_end = splitter.end();
  DateAndTime start, stop;

  // Iterate through all events (sorted by tof)
  auto itev = events.begin();
  auto itev_end = events.end();

  // This is the time of the first section. Anything before is thrown out.
  while (itspl != itspl_end) {
    // Get the splitting interval times and destination
    start = itspl->start();
    stop = itspl->stop();
    const size_t index = itspl->index();

    // Skip the events before the start of the time
    while ((itev != itev_end) && (itev->m_pulsetime < start))
      itev++;

    // Go through all the events that are in the interval (if any)
    while ((itev != itev_end) && (itev->m_pulsetime < stop)) {
      // Copy the event into another
      const T eventCopy(*itev);
      if (index < numOutputs) {
        EventListBase *myOutput = outputs[index];
        // Add the copy to the output
        myOutput->addEventQuickly(eventCopy);
      }
      ++itev;
    }

    // Go to the next interval
    ++itspl;
    // But if we reached the end, then we are done.
    if (itspl == itspl_end)
      break;

    // No need to keep looping through the filter if we are out of events
    if (itev == itev_end)
      break;
  }
  // Done!
}

//----------------------------------------------------------------------------------------------
/** Split the event list by pulse time
 */
void splitByPulseTime(Kernel::TimeSplitterType &splitter, std::map<int, EventListBase *> outputs) const {
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

//-------------------------------------------
//--------------------------------------------------
/** Split the event list into n outputs by each event's pulse time only
 */
void splitByPulseTimeHelper(Kernel::TimeSplitterType &splitter, std::map<int, EventListBase *> outputs,
                                       typename std::vector<T> &events) const {
  // Prepare to TimeSplitter Iterate through the splitter at the same time
  auto itspl = splitter.begin();
  auto itspl_end = splitter.end();
  Types::Core::DateAndTime start, stop;

  // Prepare to Events Iterate through all events (sorted by tof)
  auto itev = events.begin();
  auto itev_end = events.end();

  // Iterate (loop) on all splitters
  while (itspl != itspl_end) {
    // Get the splitting interval times and destination group
    start = itspl->start().totalNanoseconds();
    stop = itspl->stop().totalNanoseconds();
    const int index = itspl->index();

    // Skip the events before the start of the time and put to 'unfiltered'
    // EventListBase
    EventListBase *myOutput = outputs[-1];
    while (itev != itev_end) {
      if (itev->m_pulsetime < start) {
        // Record to index = -1 space
        const T eventCopy(*itev);
        myOutput->addEventQuickly(eventCopy);
        ++itev;
      } else {
        // Event within a splitter interval
        break;
      }
    }

    // Go through all the events that are in the interval (if any)
    while (itev != itev_end) {

      if (itev->m_pulsetime < stop) {
        outputs[index]->addEventQuickly(*itev);
        ++itev;
      } else {
        // Out of interval
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
/** Split the event list by pulse time
 */
// TODO/NOW - TEST
void splitByPulseTimeWithMatrix(const std::vector<int64_t> &vec_times, const std::vector<int> &vec_target,
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


void splitByPulseTimeWithMatrixHelper(const std::vector<int64_t> &vec_split_times,
                                                 const std::vector<int> &vec_split_target,
                                                 std::map<int, EventListBase *> outputs,
                                                 typename std::vector<T> &events) const {
  // Prepare to TimeSplitter Iterate through the splitter at the same time
  if (vec_split_times.size() != vec_split_target.size() + 1)
    throw std::runtime_error("Splitter time vector size and splitter target "
                             "vector size are not correct.");

  // Prepare to Events Iterate through all events (sorted by tof)
  auto itev = events.begin();
  auto itev_end = events.end();

  // Iterate (loop) on all splitters
  for (size_t i_target = 0; i_target < vec_split_target.size(); ++i_target) {
    // Get the splitting interval times and destination group
    int64_t start = vec_split_times[i_target];
    int64_t stop = vec_split_times[i_target + 1];
    const int index = vec_split_target[i_target];

    // Skip the events before the start of the time and put to 'unfiltered'
    // EventListBase
    EventListBase *myOutput = outputs[-1];
    while (itev != itev_end) {
      if (itev->m_pulsetime < start) {
        // Record to index = -1 space
        const T eventCopy(*itev);
        myOutput->addEventQuickly(eventCopy);
        ++itev;
      } else {
        // Event within a splitter interval
        break;
      }
    }

    // Go through all the events that are in the interval (if any)
    while (itev != itev_end) {

      if (itev->m_pulsetime < stop) {
        outputs[index]->addEventQuickly(*itev);
        ++itev;
      } else {
        // Out of interval
        break;
      }
    }

    // No need to keep looping through the filter if we are out of events
    if (itev == itev_end)
      break;
  } // END-WHILE Splitter
}

    friend T;
    // EventListPulsetimeFunctionsTemplate() = default;

    inline T & as_underlying()
    {
        return static_cast<T&>(*this);
    }
};
}
}