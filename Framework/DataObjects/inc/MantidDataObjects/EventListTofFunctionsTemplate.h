
#include "MantidDataObjects/EventListBaseFunctionsTemplate.h"

namespace Mantid {
namespace DataObjects {
template <typename T>
class EventListTofFunctionsTemplate : public EventListBaseFunctionsTemplate<T>
{
  

private:

    /**
     * @param events
     * @param func
     */
    void convertTofHelper(std::vector<T> &events, const std::function<double(double)> &func) {
    // iterate through all events
    for (auto &ev : events)
        ev.m_tof = func(ev.m_tof);
    }

    // --------------------------------------------------------------------------
    /** Function to do the conversion factor work
     * on either the TofEvent list or the WeightedEvent list.
     * Does NOT reverse the event list if the factor < 0
     *
     * @param events :: reference to a vector of events to change.
     * @param factor :: multiply by this
     * @param offset :: add this
     */
    void convertTofHelper(std::vector<T> &events, const double factor, const double offset) {
    // iterate through all events
    for (auto &event : events) {
        event.m_tof = event.m_tof * factor + offset;
    }
    }


// --------------------------------------------------------------------------
/** Mask out events that have a tof between tofMin and tofMax (inclusively).
 * Events are removed from the list.
 * @param events :: reference to a vector of events to change.
 * @param tofMin :: lower bound of TOF to filter out
 * @param tofMax :: upper bound of TOF to filter out
 * @returns The number of events deleted.
 */
std::size_t maskTofHelper(std::vector<T> &events, const double tofMin, const double tofMax) {
  // quick checks to make sure that the masking range is even in the data
  if (tofMin > events.rbegin()->tof())
    return 0;
  if (tofMax < events.begin()->tof())
    return 0;

  // Find the index of the first tofMin
  auto it_first = std::lower_bound(events.begin(), events.end(), tofMin);
  if ((it_first != events.end()) && (it_first->tof() < tofMax)) {
    // Something was found
    // Look for the first one > tofMax
    auto it_last = std::upper_bound(it_first, events.end(), T(tofMax));

    if (it_first >= it_last) {
      throw std::runtime_error("Event filter is all messed up"); // TODO
    }

    size_t tmp = (it_last - it_first);
    // it_last will either be at the end (if not found) or before it.
    // Erase this range from the vector
    events.erase(it_first, it_last);

    // Done! Sorting is still valid, no need to redo.
    return tmp; //(it_last - it_first); the iterators get invalid after erase
                //(on my machine)
  }
  return 0;
}

// --------------------------------------------------------------------------
/** Get the m_tof member of all events in a list
 *
 * @param events :: source vector of events
 * @param tofs :: vector to fill
 */
 void getTofsHelper(const std::vector<T> &events, std::vector<double> &tofs) {
  tofs.clear();
  for (auto itev = events.cbegin(); itev != events.cend(); ++itev)
    tofs.emplace_back(itev->m_tof);
}

// --------------------------------------------------------------------------
/** Set a list of TOFs to the current event list.
 *
 * @param events :: source vector of events
 * @param tofs :: The vector of doubles to set the tofs to.
 */
void setTofsHelper(std::vector<T> &events, const std::vector<double> &tofs) {
  if (tofs.empty())
    return;

  size_t x_size = tofs.size();
  if (events.size() != x_size)
    return; // should this throw an exception?

  for (size_t i = 0; i < x_size; ++i)
    events[i].m_tof = tofs[i];
}

//--------------------------------------------------------------------------
/** Helper function for the conversion to TOF. This handles the different
 *  event types.
 *
 * @param events the list of events
 * @param fromUnit the unit to convert from
 * @param toUnit the unit to convert to
 */
void convertUnitsViaTofHelper(typename std::vector<T> &events, Mantid::Kernel::Unit *fromUnit,
                                         Mantid::Kernel::Unit *toUnit) {
  auto itev = events.begin();
  auto itev_end = events.end();
  for (; itev != itev_end; itev++) {
    // Conver to TOF
    double tof = fromUnit->singleToTOF(itev->m_tof);
    // And back from TOF to whatever
    itev->m_tof = toUnit->singleFromTOF(tof);
  }
}

//--------------------------------------------------------------------------
/** Convert the event's TOF (x) value according to a simple output = a *
 * (input^b) relationship
 *  @param events :: templated class for the list of events
 *  @param factor :: the conversion factor a to apply
 *  @param power :: the Power b to apply to the conversion
 */
void convertUnitsQuicklyHelper(typename std::vector<T> &events, const double &factor, const double &power) {
  for (auto &event : events) {
    // Output unit = factor * (input) ^ power
    event.m_tof = factor * std::pow(event.m_tof, power);
  }
}

    friend T;
    EventListTofFunctionsTemplate() = default;

    inline T & as_underlying()
    {
        return static_cast<T&>(*this);
    }
};
}
}