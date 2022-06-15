#include "MantidDataObjects/EventListWeightFunctionsTemplate.h"
#include "MantidDataObjects/EventListErrorFunctionsTemplate.h"


namespace Mantid {
namespace DataObjects {
template <typename T>
class EventListWeightErrorFunctionsTemplate : public EventListWeightFunctionsTemplate<T>, public EventListErrorFunctionsTemplate<T>
{
  

private:

/** Integrate the events between a range of X values, or all events->
 *
 * @param events :: reference to a vector of events to change.
 * @param minX :: minimum X bin to use in integrating.
 * @param maxX :: maximum X bin to use in integrating.
 * @param entireRange :: set to true to use the entire range. minX and maxX are
 *then ignored!
 * @return the integrated number of events->
 */
double integrateHelper(std::vector<T> &events, const double minX, const double maxX,
                                  const bool entireRange) {
  double sum(0), error(0);
  integrateHelper(events, minX, maxX, entireRange, sum, error);
  return sum;
}


    /** Integrate the events between a range of X values, or all events->
 *
 * @param events :: reference to a vector of events to change.
 * @param minX :: minimum X bin to use in integrating.
 * @param maxX :: maximum X bin to use in integrating.
 * @param entireRange :: set to true to use the entire range. minX and maxX are
 *then ignored!
 * @param sum :: reference to a double to put the sum in.
 * @param error :: reference to a double to put the error in.
 */
void integrateHelper(std::vector<T> &events, const double minX, const double maxX, const bool entireRange,
                                double &sum, double &error) {
  sum = 0;
  error = 0;
  // Nothing in the list?
  if (events.empty())
    return;

  // Iterators for limits - whole range by default
  typename std::vector<T>::iterator lowit, highit;
  lowit = events.begin();
  highit = events.end();

  // But maybe we don't want the entire range?
  if (!entireRange) {
    // If a silly range was given, return 0.
    if (maxX < minX)
      return;

    // If the first element is lower that the xmin then search for new lowit
    if (lowit->tof() < minX)
      lowit = std::lower_bound(events.begin(), events.end(), minX);
    // If the last element is higher that the xmax then search for new lowit
    if ((highit - 1)->tof() > maxX) {
      highit = std::upper_bound(lowit, events.end(), T(maxX));
    }
  }

  // Sum up all the weights
  for (auto it = lowit; it != highit; ++it) {
    sum += it->weight();
    error += it->errorSquared();
  }
  error = std::sqrt(error);
}

//------------------------------------------------------------------------------------------------
/** Helper method for multiplying an event list by a scalar value with/without
 *error
 *
 * @param events: vector of events (with weights)
 * @param value: multiply all weights by this amount.
 * @param error: error on 'value'. Can be 0.
 * */
void multiplyHelper(std::vector<T> &events, const double value, const double error) {
  // Square of the value's error
  double errorSquared = error * error;
  double valueSquared = value * value;

  auto itev_end = events.end();

  if (error == 0) {
    // Error-less calculation
    for (auto itev = events.begin(); itev != itev_end; itev++) {
      itev->m_errorSquared = static_cast<float>(itev->m_errorSquared * valueSquared);
      itev->m_weight *= static_cast<float>(value);
    }
  } else {
    // Carry the scalar error
    for (auto itev = events.begin(); itev != itev_end; itev++) {
      itev->m_errorSquared =
          static_cast<float>(itev->m_errorSquared * valueSquared + errorSquared * itev->m_weight * itev->m_weight);
      itev->m_weight *= static_cast<float>(value);
    }
  }
}



    friend T;
    EventListWeightErrorFunctionsTemplate() = default;

    inline T & as_underlying()
    {
        return static_cast<T&>(*this);
    }
};
}
}