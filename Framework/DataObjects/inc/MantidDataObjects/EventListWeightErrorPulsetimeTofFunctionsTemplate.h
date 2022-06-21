#pragma once

#include "MantidDataObjects/EventListWeightErrorTofFunctionsTemplate.h"
#include "MantidDataObjects/EventListPulsetimeTofFunctionsTemplate.h"
#include "MantidKernel/DateAndTimeHelpers.h"

namespace Mantid {
namespace DataObjects {
using Types::Core::DateAndTime;

template <typename T, typename SELF>
class EventListWeightErrorPulsetimeTofFunctionsTemplate : public EventListWeightErrorTofFunctionsTemplate<T, SELF>, public EventListPulsetimeTofFunctionsTemplate<T, SELF>
{
    using EventListWeightErrorTofFunctionsTemplate<T, SELF>
    ::EventListWeightErrorFunctionsTemplate<T, SELF>
    ::EventListWeightFunctionsTemplate<T, SELF>
    ::EventListBaseFunctionsTemplate<T, SELF>
    ::events;
  
  public:
     EventListWeightErrorPulsetimeTofFunctionsTemplate(std::shared_ptr<std::vector<T>> events): 
  EventListWeightErrorTofFunctionsTemplate<T, SELF>(events), 
  EventListPulsetimeTofFunctionsTemplate<T, SELF>(events){}

protected:
    // --------------------------------------------------------------------------
    /** SUBTRACT another EventListBase from this event list.
     * The event lists are concatenated, but the weights of the incoming
     *    list are multiplied by -1.0.
     *
     * @tparam T1, T2 :: TofEvent, WeightedEvent or WeightedEventNoTime
     * @param events :: The event vector being changed.
     * @param more_events :: Another event vector being subtracted from this.
     * @return reference to this
     * */
    template <class T1, class T2> void minusHelper(std::vector<T1> &events, const std::vector<T2> &more_events) {
    // Make the end vector big enough in one go (avoids repeated re-allocations).
    events.reserve(events.size() + more_events.size());
    /* In the event of subtracting in place, calling the end() vector would make
    * it point at the wrong place
    * Using it caused a segault, Ticket #2306.
    * So we cache the end (this speeds up too).
    */
    for (const auto &ev : more_events) {
        // We call the constructor for T1. In the case of WeightedEventNoTime, the
        // pulse time will just be ignored.
        events.emplace_back(ev.tof(), ev.pulseTime(), ev.weight() * (-1.0), ev.errorSquared());
    }
    }

   
inline void compressFatEventsHelper(const std::vector<T> &events, std::vector<WeightedEvent> &out,
                                               const double tolerance, const Types::Core::DateAndTime &timeStart,
                                               const double seconds) {
  // Clear the output. We can't know ahead of time how much space to reserve :(
  out.clear();
  // We will make a starting guess of 1/20th of the number of input events.
  out.reserve(events.size() / 20);

  // The last TOF to which we are comparing.
  double lastTof = std::numeric_limits<double>::lowest();
  // For getting an accurate average TOF
  double totalTof = 0;

  // pulsetime bin information - stored as int nanoseconds because it
  // is the implementation type for DateAndTime object
  const int64_t pulsetimeStart = timeStart.totalNanoseconds();
  const auto pulsetimeDelta = static_cast<int64_t>(seconds * SEC_TO_NANO);

  // pulsetime information
  std::vector<DateAndTime> pulsetimes; // all the times for new event
  std::vector<double> pulsetimeWeights;

  // Carrying weight and error
  double weight = 0.;
  double errorSquared = 0.;
  double tofNormalization = 0.;

  // Move up to first event that has a large enough pulsetime. This is just in
  // case someone starts from after the starttime of the run. It is expected
  // that users will normally use the default which means this will only check
  // the first event.
  auto it = events.cbegin();
  for (; it != events.cend(); ++it) {
    if (it->m_pulsetime >= timeStart)
      break;
  }

  // bin if the pulses are histogrammed
  int64_t lastPulseBin = (it->m_pulsetime.totalNanoseconds() - pulsetimeStart) / pulsetimeDelta;
  // loop through events and accumulate weight
  for (; it != events.cend(); ++it) {
    const int64_t eventPulseBin = (it->m_pulsetime.totalNanoseconds() - pulsetimeStart) / pulsetimeDelta;
    if ((eventPulseBin <= lastPulseBin) && (std::fabs(it->m_tof - lastTof) <= tolerance)) {
      // Carry the error and weight
      weight += it->weight();
      errorSquared += it->errorSquared();
      double norm = calcNorm(it->errorSquared());
      tofNormalization += norm;
      // Track the average tof
      totalTof += it->m_tof * norm;
      // Accumulate the pulse times
      pulsetimes.emplace_back(it->m_pulsetime);
      pulsetimeWeights.emplace_back(norm);
    } else {
      // We exceeded the tolerance
      if (!pulsetimes.empty()) {
        // Create a new event with the average TOF and summed weights and
        // squared errors. 1 event used doesn't need to average
        if (pulsetimes.size() == 1) {
          out.emplace_back(lastTof, pulsetimes.front(), weight, errorSquared);
        } else {
          out.emplace_back(totalTof / tofNormalization,
                           Kernel::DateAndTimeHelpers::averageSorted(pulsetimes, pulsetimeWeights), weight,
                           errorSquared);
        }
      }
      // Start a new combined object
      double norm = calcNorm(it->errorSquared());
      totalTof = it->m_tof * norm;
      weight = it->weight();
      errorSquared = it->errorSquared();
      tofNormalization = norm;
      lastTof = it->m_tof;
      lastPulseBin = eventPulseBin;
      pulsetimes.clear();
      pulsetimes.emplace_back(it->m_pulsetime);
      pulsetimeWeights.clear();
      pulsetimeWeights.emplace_back(norm);
    }
  }

  // Put the last event in there too.
  if (!pulsetimes.empty()) {
    // Create a new event with the average TOF and summed weights and
    // squared errors. 1 event used doesn't need to average
    if (pulsetimes.size() == 1) {
      out.emplace_back(lastTof, pulsetimes.front(), weight, errorSquared);
    } else {
      out.emplace_back(totalTof / tofNormalization,
                       Kernel::DateAndTimeHelpers::averageSorted(pulsetimes, pulsetimeWeights), weight, errorSquared);
    }
  }

  // If you have over-allocated by more than 5%, reduce the size.
  size_t excess_limit = out.size() / 20;
  if ((out.capacity() - out.size()) > excess_limit) {
    out.shrink_to_fit();
  }
}

    friend T;
    // EventListWeightErrorPulsetimeTofFunctionsTemplate() = default;

    inline T & as_underlying()
    {
        return static_cast<T&>(*this);
    }
};
}
}