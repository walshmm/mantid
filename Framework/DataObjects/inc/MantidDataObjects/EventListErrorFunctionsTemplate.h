
namespace Mantid {
namespace DataObjects {
template <typename T>
class EventListErrorFunctionsTemplate
{
  

private:

    // --------------------------------------------------------------------------
    /** Get the weight error member of all events in a list
     *
     * @param events :: source vector of events
     * @param weightErrors :: vector to fill
     */
    template <class T>
    void EventListBase::getWeightErrorsHelper(const std::vector<T> &events, std::vector<double> &weightErrors) {
    weightErrors.clear();
    weightErrors.reserve(events.size());
    std::transform(events.cbegin(), events.cend(), std::back_inserter(weightErrors),
                    [](const auto &event) { return event.error(); });
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