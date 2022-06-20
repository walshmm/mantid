
#include "MantidDataObjects/EventListBaseFunctionsTemplate.h"

namespace Mantid {
namespace DataObjects {
template <typename T, typename SELF>
class EventListTofFunctionsTemplate : public EventListBaseFunctionsTemplate<T, SELF>
{
  public:
  
  void sort(const EventSortType order) const {
    if (order == UNSORTED) {
      return; // don't bother doing anything. Why did you ask to unsort?
    } else if (order == TOF_SORT) {
      this->sortTof();
    }  else {
      throw runtime_error("Invalid sort type in EventListTofFunctionsTemplate::sort(EventSortType)");
    }
  }

private:

/**
 * @return The maximum tof value for the list of events->
 */
double getTofMax() const {
  // set up as the minimum available double
  double tMax = std::numeric_limits<double>::lowest();

  // no events is a soft error
  if (this->empty())
    return tMax;

  // when events are ordered by tof just need the first value
  if (this->order == TOF_SORT) {
    return this->events->rbegin()->tof();
  }

  // now we are stuck with a linear search
  size_t numEvents = this->events.size();
  double temp = tMax; // start with the smallest possible value
  for (size_t i = 0; i < numEvents; i++) {
    temp = this->events[i].tof();
    if (temp > tMax)
      tMax = temp;
  }
  return tMax;
}

// --------------------------------------------------------------------------
/**
 * @return The minimum tof value for the list of the events->
 */
double getTofMin() const {
  // set up as the maximum available double
  double tMin = std::numeric_limits<double>::max();

  // no events is a soft error
  if (this->empty())
    return tMin;

  // when events are ordered by tof just need the first value
  if (this->order == TOF_SORT) {
    return this->events->begin()->tof();
  }

  // now we are stuck with a linear search
  double temp = tMin; // start with the largest possible value
  size_t numEvents = this->events.size();
  for (size_t i = 0; i < numEvents; i++) {
    temp = this->events[i].tof();
    if (temp < tMin)
      tMin = temp;
  }
  return tMin;
}

// --------------------------------------------------------------------------
/**
 * Convert the time of flight by tof'=tof*factor+offset
 * @param factor :: The value to scale the time-of-flight by
 * @param offset :: The value to shift the time-of-flight by
 */
void convertTof(const double factor, const double offset) {
  // fix the histogram parameter
  auto &x = mutableX();
  x *= factor;
  x += offset;

  if ((factor < 0.) && (this->getSortType() == TOF_SORT))
    this->reverse();

  if (this->getNumberEvents() <= 0)
    return;
  
  this->convertTofHelper(*(this->events), factor, offset);

}

/**
 * @param func Function to do the conversion.
 * @param sorting How the events are sorted after the operation. 0 = unsorted
 * (default),
 * positive = unchanged, negative = reverse.
 */
void convertTof(std::function<double(double)> func, const int sorting) {
  // fix the histogram parameter
  MantidVec &x = dataX();
  transform(x.begin(), x.end(), x.begin(), func);

  // do nothing if sorting > 0
  if (sorting == 0) {
    this->setSortOrder(UNSORTED);
  } else if ((sorting < 0) && (this->getSortType() == TOF_SORT)) {
    this->reverse();
  }

  if (this->getNumberEvents() <= 0)
    return;

  this->convertTofHelper(*(this->events), func);

}

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


/**
 * Mask out events by the condition vector.
 * Events are removed from the list.
 * @param mask :: condition vector
 */
void maskCondition(const std::vector<bool> &mask) {

  // mask size must match the number of events
  if (this->getNumberEvents() != mask.size())
    throw std::runtime_error("EventListBase::maskTof: tofMax must be > tofMin");

  // don't do anything with an emply list
  if (this->getNumberEvents() == 0)
    return;

  // Convert the list
  size_t numOrig = 0;
  size_t numDel = 0;
  
  numOrig = this->events->size();
  numDel = this->maskConditionHelper(*(this->events), mask);
   

  if (numDel >= numOrig)
    this->clear(false);
}


/**
 * Mask out events that have a tof between tofMin and tofMax (inclusively).
 * Events are removed from the list.
 * @param tofMin :: lower bound of TOF to filter out
 * @param tofMax :: upper bound of TOF to filter out
 */
void maskTof(const double tofMin, const double tofMax) {
  if (tofMax <= tofMin)
    throw std::runtime_error("EventListBase::maskTof: tofMax must be > tofMin");

  // don't do anything with an emply list
  if (this->getNumberEvents() == 0)
    return;

  // Start by sorting by tof
  this->sortTof();

  // Convert the list
  size_t numOrig = 0;
  size_t numDel = 0;

  numOrig = this->events->size();
  numDel = this->maskTofHelper(*(this->events), tofMin, tofMax);
   

  if (numDel >= numOrig)
    this->clear(false);
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


/** Get the times-of-flight of each event in this EventListBase.
 *
 * @return by copy a vector of doubles of the tof() value
 */
std::vector<double> getTofs() const {
  std::vector<double> tofs;
  this->getTofs(tofs);
  return tofs;
}

/** Fill a vector with the list of TOFs
 *  @param tofs :: A reference to the vector to be filled
 */
void getTofs(std::vector<double> &tofs) const {
  // Set the capacity of the vector to avoid multiple resizes
  tofs.reserve(this->getNumberEvents());
  this->getTofsHelper(*(this->events), tofs);
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
/**
 * Set a list of TOFs to the current event list. Modify the units if necessary.
 *
 * @param tofs :: The vector of doubles to set the tofs to.
 */
void setTofs(const MantidVec &tofs) {
  this->order = UNSORTED;

  // Convert the list
  this->setTofsHelper(*(this->events), tofs);
    
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
/** Converts the X units in each event by going through TOF.
 * Note: if the unit conversion reverses the order, use "reverse()" to flip it
 *back.
 *
 * @param fromUnit :: the Unit describing the input unit. Must be initialized.
 * @param toUnit :: the Unit describing the output unit. Must be initialized.
 */
void convertUnitsViaTof(Mantid::Kernel::Unit *fromUnit, Mantid::Kernel::Unit *toUnit) {
  // Check for initialized
  if (!fromUnit || !toUnit)
    throw std::runtime_error("EventListBase::convertUnitsViaTof(): one of the units is NULL!");
  if (!fromUnit->isInitialized())
    throw std::runtime_error("EventListBase::convertUnitsViaTof(): fromUnit is not initialized!");
  if (!toUnit->isInitialized())
    throw std::runtime_error("EventListBase::convertUnitsViaTof(): toUnit is not initialized!");


  convertUnitsViaTofHelper(*(this->events), fromUnit, toUnit);
  
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
 *  @param factor :: the conversion factor a to apply
 *  @param power :: the Power b to apply to the conversion
 */
void convertUnitsQuickly(const double &factor, const double &power) {
  auto convertEvents = *(this->events);
  convertUnitsQuicklyHelper(convertEvents, factor, power);  
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