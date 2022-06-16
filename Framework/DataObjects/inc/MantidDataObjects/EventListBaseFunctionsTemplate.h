
namespace Mantid {
namespace DataObjects {
template <typename T>
class EventListBaseFunctionsTemplate
{
  public:

/// Mask the spectrum to this value. Removes all events->
void clearData() { this->clear(false); }
  /** Clear the list of events and any
 * associated detector ID's.
 * */
void clear(const bool removeDetIDs) {
  if (mru)
    mru->deleteIndex(this);
  this->events.clear();
  std::vector<T>().swap(this->events); // STL Trick to release memory
  if (removeDetIDs)
    this->clearDetectorIDs();
}


/** Reserve a certain number of entries in event list of the specified eventType
 *
 * Calls std::vector<>::reserve() in order to pre-allocate the length of the
 *event list vector.
 *
 * @param num :: number of events that will be in this EventListBase
 */
void reserve(size_t num) {
  this->events.reserve(num);
}

  // --------------------------------------------------------------------------
  /** Append an event to the histogram, without clearing the cache, to make it
   *faster.
   * NOTE: Only call this on a un-weighted event list!
   *
   * @param event :: TofEvent to add at the end of the list.
   * */
  inline void addEventQuickly(const T &event) {
    events->emplace_back(event);
    this->order = UNSORTED;
  }

  // --------------------------------------------------------------------------
/** Compress the event list by grouping events with the same
 * TOF (within a given tolerance). PulseTime is ignored.
 * The event list will be switched to WeightedEventNoTime.
 *
 * @param tolerance :: how close do two event's TOF have to be to be considered
 *the same.
 * @param destination :: EventListBase that will receive the compressed events-> Can
 *be == this.
 */

//NOTE: Dont like that this is operating with a parent type as a parameter
void compressEvents(double tolerance, EventList *destination) {
  if (!this->empty()) {
    this->sortTof();

    //in the wrapper do this before calling compress Events
    destination->clear()
    destination->switchTo(WEIGHTED_NOTIME)


    //should work because the wrapper's equal operator refers to the wrapped EventList for equality
    if (destination == this) {
            // Put results in a temp output
            std::vector<WeightedEventNoTime> out;
            compressEventsHelper(*(this->events), out, tolerance);
            // Put it back
            this->events->swap(out);
        } else {
            compressEventsHelper(*(this->events), destination->events, tolerance);
    }
    
  }
  // In all cases, you end up WEIGHTED_NOTIME.
  destination->eventType = WEIGHTED_NOTIME;
  // The sort is still valid!
  destination->order = TOF_SORT;
  // Empty out storage for vectors that are now unused.
  destination->clearUnused();
}

void compressFatEvents(const double tolerance, const Mantid::Types::Core::DateAndTime &timeStart,
                                  const double seconds, EventListBase *destination) {
                                      
  // only worry about non-empty EventLists
  if (!this->empty()) {
    if (destination == this) {
      // Put results in a temp output
      std::vector<WeightedEvent> out;
      compressFatEventsHelper(*(this->events), out, tolerance, timeStart, seconds);
      // Put it back
      *(this->events).swap(out);
    } else {
      compressFatEventsHelper(*(this->events), destination->events, tolerance, timeStart, seconds);
    }
  }
  // In all cases, you end up WEIGHTED_NOTIME.
  destination->eventType = WEIGHTED;
  // The sort order is pulsetimetof as we've compressed out the tolerance
  destination->order = PULSETIMETOF_SORT;
  // Empty out storage for vectors that are now unused.
  destination->clearUnused();
}

// --------------------------------------------------------------------------
/** Sort events by TOF in one thread */
void sortTof() const {
  // nothing to do
  if (this->order == TOF_SORT)
    return;

  // Avoid sorting from multiple threads
  std::lock_guard<std::mutex> _lock(m_sortMutex);
  // If the list was sorted while waiting for the lock, return.
  if (this->order == TOF_SORT)
    return;

  // TODO determine how these are setup to compare

  tbb::parallel_sort(events->begin(), events->end());
  // Save the order to avoid unnecessary re-sorting.
  this->order = TOF_SORT;
}

private:

/** Get the weight error of each event in this EventListBase.
 *
 * @return by copy a vector of doubles of the weight() value
 */
std::vector<double> getWeightErrors() const {
  std::vector<double> weightErrors;
  this->getWeightErrors(weightErrors);
  return weightErrors;
}

/** Fill a vector with the list of Weight Errors
 *  @param weightErrors :: A reference to the vector to be filled
 */
void getWeightErrors(std::vector<double> &weightErrors) const {
  // Set the capacity of the vector to avoid multiple resizes
  weightErrors.reserve(this->events.size());
  weightErrors.assign(this->events.size(), 1.0);
}

/** Get the weight of each event in this EventListBase.
 *
 * @return by copy a vector of doubles of the weight() value
 */
std::vector<double> getWeights() const {
  std::vector<double> weights;
  this->getWeights(weights);
  return weights;
}
/** Fill a vector with the list of Weights
 *  @param weights :: A reference to the vector to be filled
 */
void getWeights(std::vector<double> &weights) const {
  // Set the capacity of the vector to avoid multiple resizes
  weights.reserve(this->events.size());
    // not a weighted event type, return 1.0 for all.
  weights.assign(this->events.size(), 1.0);
}

  // --------------------------------------------------------------------------
/** Mask out events by the condition vector.
 * Events are removed from the list.
 * @param events :: reference to a vector of events to change.
 * @param mask :: condition vector
 * @returns The number of events deleted.
 */
std::size_t maskConditionHelper(std::vector<T> &events, const std::vector<bool> &mask) {

  // runs through the two synchronized vectors and delete elements
  // for condition false
  auto itm = std::find(mask.begin(), mask.end(), false);
  auto first = events.begin() + (itm - mask.begin());

  if (itm != mask.end()) {
    for (auto ite = first; ++ite != events.end() && ++itm != mask.end();) {
      if (*itm != false) {
        *first++ = std::move(*ite);
      }
    }
  }

  auto n = events.end() - first;
  if (n != 0)
    events.erase(first, events.end());

  return n;
}

    friend T;
    std::vector<T> events;
    EventListBaseFunctionsTemplate() = default;

    inline T & as_underlying()
    {
        return static_cast<T&>(*this);
    }
};
}
}