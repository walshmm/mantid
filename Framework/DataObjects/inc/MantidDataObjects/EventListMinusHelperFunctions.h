
#include "MantidDataObjects/Events.h"

namespace Mantid {
namespace DataObjects {
template <class T1, class T2>
class EventListMinusHelperFunctions
{
public:

// --------------------------------------------------------------------------
/** SUBTRACT another EventListBase from this event list.
 * The event lists are concatenated, but the weights of the incoming
 *    list are multiplied by -1.0.
 *
 * @param more_events :: Another EventListBase.
 * @return reference to this
 * */
T1 &operator-=(const T2 &more_events) {
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


  minusHelper(*(this->events), more_events->events);
  // No guaranteed order
  this->order = UNSORTED;

  // NOTE: What to do about detector ID's?
  return *this;
}

private:
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

    friend T1;
    friend T2;
    EventListMinusHelperFunctions() = default;
};



//TofEvent doesnt get one as it gets converted to a WeightedEvent
template <class T>
class EventListPermutationsMinusHelperFunctions :
public  EventListMinusHelperFunctions<T, Types::Event::TofEvent>,
public  EventListMinusHelperFunctions<T, WeightedEvent>,
public  EventListMinusHelperFunctions<T, WeightedEventNoTime>,
public  EventListMinusHelperFunctions<T, TofEventNoTime>
{

};

}
}