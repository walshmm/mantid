
namespace Mantid {
namespace DataObjects {
template <typename T>
class EventListBaseFunctionsTemplate
{
  

private:

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
    EventListBaseFunctionsTemplate() = default;

    inline T & as_underlying()
    {
        return static_cast<T&>(*this);
    }
};
}
}