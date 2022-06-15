
namespace Mantid {
namespace DataObjects {
template <typename T>
class EventListWeightFunctionsTemplate
{
  

private:

    // --------------------------------------------------------------------------
/** Get the weight member of all events in a list
 *
 * @param events :: source vector of events
 * @param weights :: vector to fill
 */
template <class T> void EventListBase::getWeightsHelper(const std::vector<T> &events, std::vector<double> &weights) {
  weights.clear();
  weights.reserve(events.size());
  std::transform(events.cbegin(), events.cend(), std::back_inserter(weights),
                 [](const auto &event) { return event.weight(); });
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