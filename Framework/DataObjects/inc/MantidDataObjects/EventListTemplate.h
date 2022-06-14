#include "MantidDataObjects/Events.h"

namespace Mantid {
namespace DataObjects {
template <typename T>
class EventListTemplate
{
    // T &get()
    // {
    //     return instance;
    // }
    // void set(T newInstance){
    //     this->instance = newInstance
    // }

      // --------------------------------------------------------------------------
  /** Append an event to the histogram, without clearing the cache, to make it
   *faster.
   * NOTE: Only call this on a un-weighted event list!
   *
   * @param event :: TofEvent to add at the end of the list.
   * */
  inline void addEventQuickly(const Types::Event::TofEvent &event) {
    this->eventList.get().addEventQuickly(event);
  }

  // --------------------------------------------------------------------------
  /** Append an event to the histogram, without clearing the cache, to make it
   * faster.
   * @param event :: WeightedEvent to add at the end of the list.
   * */
  inline void addEventQuickly(const WeightedEvent &event) {
    this->eventList.get().addEventQuickly(event);
  }

  // --------------------------------------------------------------------------
  /** Append an event to the histogram, without clearing the cache, to make it
   * faster.
   * @param event :: WeightedEventNoTime to add at the end of the list.
   * */
  inline void addEventQuickly(const WeightedEventNoTime &event) {
    this->eventList.addEventQuickly(event);
  }

private:
    friend T;
    EventListTemplate() = default;

    T eventList;

    inline T & as_underlying()
    {
        return static_cast<T&>(*this);
    }
};
}
}