#include "MantidDataObjects/EventListPulsetimeFunctionsTemplate.h"
#include "MantidDataObjects/EventListTofFunctionsTemplate.h"


namespace Mantid {
namespace DataObjects {
template <typename T>
class EventListPulsetimeTofFunctionsTemplate : public EventListPulsetimeFunctionsTemplate<T>, public EventListTofFunctionsTemplate<T>
{
private:
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