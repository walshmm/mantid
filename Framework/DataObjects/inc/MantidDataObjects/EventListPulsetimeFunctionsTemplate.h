
namespace Mantid {
namespace DataObjects {
template <typename T>
class EventListPulsetimeFunctionsTemplate
{
  

private:
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
template <class T>
typename std::vector<T>::const_iterator EventListBase::findFirstPulseEvent(const std::vector<T> &events,
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
 * @param events :: reference to a vector of events to change.
 * @param seconds :: The value to shift the pulsetime by, in seconds
 */
template <class T> void EventListBase::addPulsetimeHelper(std::vector<T> &events, const double seconds) {
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
template <class T> void EventListBase::addPulsetimesHelper(std::vector<T> &events, const std::vector<double> &seconds) {
  auto eventIterEnd{events->end()};
  auto secondsIter{seconds.cbegin()};
  for (auto eventIter = events->begin(); eventIter < eventIterEnd; ++eventIter, ++secondsIter) {
    eventIter->m_pulsetime += *secondsIter;
  }
}

// --------------------------------------------------------------------------
/** Get the pulsetimes of all events in a list
 *
 * @param events :: source vector of events
 * @param times :: vector to fill
 */
template <class T>
void EventListBase::getPulseTimesHelper(const std::vector<T> &events,
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
template <class T>
void EventListBase::filterByPulseTimeHelper(std::vector<T> &events, DateAndTime start, DateAndTime stop,
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
/** Perform an in-place filtering on a vector of either TofEvent's or
 *WeightedEvent's
 *
 * @param splitter :: a TimeSplitterType where all the entries (start/end time)
 *indicate events
 *     that will be kept. Any other events will be deleted.
 * @param events :: either this->events or this->weightedevents.
 */
template <class T>
void EventListBase::filterInPlaceHelper(Kernel::TimeSplitterType &splitter, typename std::vector<T> &events) {
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
template <class T>
void EventListBase::splitByTimeHelper(Kernel::TimeSplitterType &splitter, std::vector<EventListBase *> outputs,
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

//-------------------------------------------
//--------------------------------------------------
/** Split the event list into n outputs by each event's pulse time only
 */
template <class T>
void EventListBase::splitByPulseTimeHelper(Kernel::TimeSplitterType &splitter, std::map<int, EventListBase *> outputs,
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

template <class T>
void EventListBase::splitByPulseTimeWithMatrixHelper(const std::vector<int64_t> &vec_split_times,
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
    EventListTemplate() = default;

    inline T & as_underlying()
    {
        return static_cast<T&>(*this);
    }
};
}
}