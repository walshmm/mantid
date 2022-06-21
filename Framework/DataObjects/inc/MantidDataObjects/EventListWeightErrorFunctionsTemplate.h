#pragma once

#include "MantidDataObjects/EventListWeightFunctionsTemplate.h"
#include "MantidDataObjects/EventListErrorFunctionsTemplate.h"


namespace Mantid {
namespace DataObjects {
template <typename T, typename SELF>
class EventListWeightErrorFunctionsTemplate : public EventListWeightFunctionsTemplate<T, SELF>, public EventListErrorFunctionsTemplate<T, SELF>
{
    using EventListWeightFunctionsTemplate
  ::EventListBaseFunctionsTemplate
  ::events;

  public:
    EventListWeightErrorFunctionsTemplate(std::shared_ptr<std::vector<T>> events): 
  EventListWeightFunctionsTemplate<T, SELF>(events), 
  EventListErrorFunctionsTemplate<T, SELF>(events){}

protected:

// --------------------------------------------------------------------------
/** Integrate the events between a range of X values, or all events->
 *
 * @param minX :: minimum X bin to use in integrating.
 * @param maxX :: maximum X bin to use in integrating.
 * @param entireRange :: set to true to use the entire range. minX and maxX are
 *then ignored!
 * @return the integrated number of events->
 */
double integrate(const double minX, const double maxX, const bool entireRange) const {
  double sum(0), error(0);
  integrate(minX, maxX, entireRange, sum, error);
  return sum;
}

/** Integrate the events between a range of X values, or all events->
 *
 * @param minX :: minimum X bin to use in integrating.
 * @param maxX :: maximum X bin to use in integrating.
 * @param entireRange :: set to true to use the entire range. minX and maxX are
 *then ignored!
 * @param sum :: place holder for the resulting sum
 * @param error :: place holder for the resulting sum of errors
 * @return the integrated number of events->
 */
void integrate(const double minX, const double maxX, const bool entireRange, double &sum,
                          double &error) const {
  sum = 0;
  error = 0;
  if (!entireRange) {
    // The event list must be sorted by TOF!
    this->sortTof();
  }

  integrateHelper(*(this->events), minX, maxX, entireRange, sum, error);
    
}

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
/** Operator to multiply the weights in this EventListBase by an error-less scalar.
 * Use multiply(value,error) if you wish to multiply by a real variable with an
 *error!
 *
 * The event list switches to WeightedEvent's if needed.
 * Note that if the multiplier is exactly 1.0, the list is NOT switched to
 *WeightedEvents - nothing happens.
 *
 * @param value :: multiply by this
 * @return reference to this
 */
EventListBase &operator*=(const double value) {
  this->multiply(value);
  return *this;
}

//------------------------------------------------------------------------------------------------
/** Operator to divide the weights in this EventListBase by an error-less scalar.
 * Use divide(value,error) if your scalar has an error!
 * This simply calls the equivalent function: multiply(1.0/value).
 *
 * @param value :: divide by this
 * @return reference to this
 * @throw std::invalid_argument if value == 0; cannot divide by zero.
 */
EventListBase &operator/=(const double value) {
  if (value == 0.0)
    throw std::invalid_argument("EventListBase::divide() called with value of 0.0. Cannot divide by zero.");
  this->multiply(1.0 / value, 0.0);
  return *this;
}


//------------------------------------------------------------------------------------------------
/** Multiply the weights in this event list by a scalar variable with an error;
 * though the error can be 0.0
 *
 * The event list switches to WeightedEvent's if needed.
 * Note that if the multiplier is exactly 1.0 and the error is exactly 0.0, the
 *list is NOT switched to WeightedEvents - nothing happens.
 *
 * Given:
 *  - A is the weight, variance \f$\sigma_A \f$
 *  - B is the scalar multiplier, variance \f$\sigma_B \f$
 *
 * The error propagation formula used is:
 *
 * \f[ \left(\frac{\sigma_f}{f}\right)^2 = \left(\frac{\sigma_A}{A}\right)^2 +
 *\left(\frac{\sigma_B}{B}\right)^2 + 2\frac{\sigma_A\sigma_B}{AB}\rho_{AB} \f]
 *
 * \f$ \rho_{AB} \f$ is the covariance between A and B, which we take to be 0
 *(uncorrelated variables).
 * Therefore, this reduces to:
 * \f[ \sigma_{AB}^2 = B^2 \sigma_A^2 + A^2 \sigma_B ^ 2  \f]
 *
 * In the case of no error:
 *  - The weight is simply \f$ aA \f$
 *  - The error \f$ \sigma_A \f$ becomes \f$ \sigma_{aA} = a \sigma_{A} \f$
 *
 * @param value: multiply all weights by this amount.
 * @param error: error on 'value'. Can be 0.
 */
void multiply(const double value, const double error) {
  // Do nothing if multiplying by exactly one and there is no error
  if ((value == 1.0) && (error == 0.0))
    return;

  // TODO abstract out to wrapper?
  this->switchTo(WEIGHTED);

  multiplyHelper(*(this->events), value, error);

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
    // EventListWeightErrorFunctionsTemplate() = default;

    inline T & as_underlying()
    {
        return static_cast<T&>(*this);
    }
};
}
}