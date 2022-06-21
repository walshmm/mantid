#pragma once

#include "MantidDataObjects/EventListWeightErrorFunctionsTemplate.h"
#include "MantidDataObjects/EventListTofFunctionsTemplate.h"


namespace Mantid {
namespace DataObjects {
template <typename T, typename SELF>
class EventListWeightErrorTofFunctionsTemplate : public EventListWeightErrorFunctionsTemplate<T, SELF>, public EventListTofFunctionsTemplate<T, SELF>
{
  using EventListWeightErrorFunctionsTemplate
  ::EventListWeightFunctionsTemplate
  ::EventListBaseFunctionsTemplate
  ::events;

  public:

  EventListWeightErrorTofFunctionsTemplate(std::shared_ptr<std::vector<T>> events): 
  EventListWeightErrorFunctionsTemplate<T, SELF>(events), 
  EventListTofFunctionsTemplate<T, SELF>(events) 
  {}

protected:

// --------------------------------------------------------------------------
/** Create an EventListBase from a histogram. This converts bins to weighted
 * events->
 * Any existing events are cleared.
 *
 * @param inSpec :: ISpectrum ptr to histogram data.
 * @param GenerateZeros :: if true, generate event(s) for empty bins
 * @param GenerateMultipleEvents :: if true, create several evenly-spaced fake
 *events inside the bin
 * @param MaxEventsPerBin :: max number of events to generate in one bin, if
 *GenerateMultipleEvents
 */
void createFromHistogram(const ISpectrum *inSpec, bool GenerateZeros, bool GenerateMultipleEvents,
                                    int MaxEventsPerBin) {
  // Fresh start
  this->clear(true);
  // Get the input histogram
  const MantidVec &X = inSpec->readX();
  const MantidVec &Y = inSpec->readY();
  const MantidVec &E = inSpec->readE();
  if (Y.size() + 1 != X.size())
    throw std::runtime_error("Expected a histogram (X vector should be 1 longer than the Y vector)");

  // Copy detector IDs and spectra
  this->copyInfoFrom(*inSpec);
  // We need weights but have no way to set the time. So use weighted, no time
  this->switchTo(WEIGHTED_NOTIME);
  if (GenerateZeros)
    this->events->reserve(Y.size());

  for (size_t i = 0; i < X.size() - 1; i++) {
    double weight = Y[i];
    if ((weight != 0.0 || GenerateZeros) && std::isfinite(weight)) {
      double error = E[i];
      // Also check that the error is not a bad number
      if (std::isfinite(error)) {
        if (GenerateMultipleEvents) {
          // --------- Multiple events per bin ----------
          double errorSquared = error * error;
          // Find how many events to fake
          double val = weight / E[i];
          val *= val;
          // Convert to int with slight rounding up. This is to avoid rounding
          // errors
          auto numEvents = int(val + 0.2);
          if (numEvents < 1)
            numEvents = 1;
          if (numEvents > MaxEventsPerBin)
            numEvents = MaxEventsPerBin;
          // Scale the weight and error for each
          weight /= numEvents;
          errorSquared /= numEvents;

          // Spread the TOF. e.g. 2 events = 0.25, 0.75.
          double tofStep = (X[i + 1] - X[i]) / (numEvents);
          for (size_t j = 0; j < size_t(numEvents); j++) {
            double tof = X[i] + tofStep * (0.5 + double(j));
            // Create and add the event
            // TODO: try emplace_back() here.
            events->emplace_back(tof, weight, errorSquared);
          }
        } else {
          // --------- Single event per bin ----------
          // TOF = midpoint of the bin
          double tof = (X[i] + X[i + 1]) / 2.0;
          // Error squared is carried in the event
          double errorSquared = E[i];
          errorSquared *= errorSquared;
          // Create and add the event
          events->emplace_back(tof, weight, errorSquared);
        }
      } // error is nont NAN or infinite
    }   // weight is non-zero, not NAN, and non-infinite
  }     // (each bin)

  // Set the X binning parameters
  this->setX(inSpec->ptrX());

  // Manually set that this is sorted by TOF, since it is. This will make it
  // "threadSafe" in other algos.
  this->setSortOrder(TOF_SORT);
}

inline double calcNorm(const double errorSquared) {
  if (errorSquared == 0.)
    return 0;
  else if (errorSquared == 1.)
    return 1.;
  else
    return 1. / std::sqrt(errorSquared);
}

     // --------------------------------------------------------------------------
/** Compress the event list by grouping events with the same TOF.
 *
 * @param events :: input event list.
 * @param out :: output WeightedEventNoTime vector.
 * @param tolerance :: how close do two event's TOF have to be to be considered
 *the same.
 */

template <class T>
inline void compressEventsHelper(const std::vector<T> &events, std::vector<WeightedEventNoTime> &out,
                                            double tolerance) {
  // Clear the output. We can't know ahead of time how much space to reserve :(
  out.clear();
  // We will make a starting guess of 1/20th of the number of input events->
  out.reserve(events.size() / 20);

  // The last TOF to which we are comparing.
  double lastTof = std::numeric_limits<double>::lowest();
  // For getting an accurate average TOF
  double totalTof = 0;
  int num = 0;
  // Carrying weight, error, and normalization
  double weight = 0;
  double errorSquared = 0;
  double normalization = 0.;

  for (auto it = events.cbegin(); it != events.cend(); it++) {
    if ((it->m_tof - lastTof) <= tolerance) {
      // Carry the error and weight
      weight += it->weight();
      errorSquared += it->errorSquared();
      // Track the average tof
      num++;
      const double norm = calcNorm(it->errorSquared());
      normalization += norm;
      totalTof += it->m_tof * norm;
    } else {
      // We exceeded the tolerance
      // Create a new event with the average TOF and summed weights and
      // squared errors.
      if (num == 1) {
        // last time-of-flight is the only one contributing
        out.emplace_back(lastTof, weight, errorSquared);
      } else if (num > 1) {
        out.emplace_back(totalTof / normalization, weight, errorSquared);
      }
      // Start a new combined object
      num = 1;
      const double norm = calcNorm(it->errorSquared());
      normalization = norm;
      totalTof = it->m_tof * norm;
      weight = it->weight();
      errorSquared = it->errorSquared();
      lastTof = it->m_tof;
    }
  }

  // Put the last event in there too with the average TOF and summed weights and
  // squared errors.
  if (num == 1) {
    // last time-of-flight is the only one contributing
    out.emplace_back(lastTof, weight, errorSquared);
  } else if (num > 1) {
    out.emplace_back(totalTof / normalization, weight, errorSquared);
  }

  // If you have over-allocated by more than 5%, reduce the size.
  size_t excess_limit = out.size() / 20;
  if ((out.capacity() - out.size()) > excess_limit) {
    out.shrink_to_fit();
  }
}

// --------------------------------------------------------------------------
/** Compress the event list by grouping events with the same TOF.
 * Performs the compression in parallel.
 *
 * @param events :: input event list.
 * @param out :: output WeightedEventNoTime vector.
 * @param tolerance :: how close do two event's TOF have to be to be considered
 *the same.
 */

template <class T>
void compressEventsParallelHelper(const std::vector<T> &events, std::vector<WeightedEventNoTime> &out,
                                             double tolerance) {
  // Create a local output vector for each thread
  int numThreads = PARALLEL_GET_MAX_THREADS;
  std::vector<std::vector<WeightedEventNoTime>> outputs(numThreads);
  // This is how many events to process in each thread.
  size_t numPerBlock = events.size() / numThreads;

  // Do each block in parallel
  PARALLEL_FOR_NO_WSP_CHECK()
  for (int thread = 0; thread < numThreads; thread++) {
    // The local output vector
    std::vector<WeightedEventNoTime> &localOut = outputs[thread];
    // Reserve a bit of space to avoid excess copying
    localOut.clear();
    localOut.reserve(numPerBlock / 20);

    // The last TOF to which we are comparing.
    double lastTof = std::numeric_limits<double>::lowest();
    // For getting an accurate average TOF
    double totalTof = 0;
    int num = 0;
    // Carrying weight, error, and normalization
    double weight = 0;
    double errorSquared = 0;
    double normalization = 0.;

    // Separate the
    typename std::vector<T>::const_iterator it = events.begin() + thread * numPerBlock;
    typename std::vector<T>::const_iterator it_end = events.begin() + (thread + 1) * numPerBlock; // cache for speed
    if (thread == numThreads - 1)
      it_end = events.end();
    for (; it != it_end; ++it) {
      if ((it->m_tof - lastTof) <= tolerance) {
        // Carry the error and weight
        weight += it->weight();
        errorSquared += it->errorSquared();
        // Track the average tof
        num++;
        const double norm = calcNorm(it->errorSquared());
        normalization += norm;
        totalTof += it->m_tof * norm;
      } else {
        // We exceeded the tolerance
        if (num > 0) {
          // Create a new event with the average TOF and summed weights and
          // squared errors.
          localOut.emplace_back(totalTof / normalization, weight, errorSquared);
        }
        // Start a new combined object
        num = 1;
        const double norm = calcNorm(it->errorSquared());
        normalization = norm;
        totalTof = it->m_tof * norm;
        weight = it->weight();
        errorSquared = it->errorSquared();
        lastTof = it->m_tof;
      }
    }

    // Put the last event in there too.
    if (num > 0) {
      // Create a new event with the average TOF and summed weights and squared
      // errors.
      localOut.emplace_back(totalTof / normalization, weight, errorSquared);
    }
  }

  // Clear the output. Reserve the required size
  out.clear();
  size_t numEvents = 0;
  for (int thread = 0; thread < numThreads; thread++)
    numEvents += outputs[thread].size();
  out.reserve(numEvents);

  // Re-join all the outputs
  for (int thread = 0; thread < numThreads; thread++)
    out.insert(out.end(), outputs[thread].begin(), outputs[thread].end());
}

//------------------------------------------------------------------------------------------------
/** Multiply the weights in this event list by a histogram.
 * The event list switches to WeightedEvent's if needed.
 * NOTE: no unit checks are made (or possible to make) to compare the units of X
 *and tof() in the EventListBase.
 *
 * The formula used for calculating the error on the neutron weight is:
 * \f[ \sigma_{f}^2 = B^2 \sigma_A^2 + A^2 \sigma_B ^ 2  \f]
 *
 * where:
 *  * A is the weight of the event
 *  * B is the weight of the BIN that the event falls in
 *  * \f$\sigma_A\f$ is the error (not squared) of the weight of the event
 *  * \f$\sigma_B\f$ is the error (not squared) of the bin B
 *  * f is the resulting weight of the multiplied event
 *
 * @param X: bins of the multiplying histogram.
 * @param Y: value to multiply the weights.
 * @param E: error on the value to multiply.
 * @throw invalid_argument if the sizes of X, Y, E are not consistent.
 */
void multiply(const MantidVec &X, const MantidVec &Y, const MantidVec &E) {

    // Switch to weights if needed.
    // TODO abstract out to wrapper?
    this->switchTo(WEIGHTED);
    // Fall through

    this->sortTof();
    multiplyHistogramHelper(*(this->events), X, Y, E);

}

//------------------------------------------------------------------------------------------------
/** Helper method for multiplying an event list by a histogram with error
 *
 * @param events: vector of events (with weights)
 * @param X: bins of the multiplying histogram.
 * @param Y: value to multiply the weights.
 * @param E: error on the value to multiply.
 * @throw invalid_argument if the sizes of X, Y, E are not consistent.
 * */
void multiplyHistogramHelper(std::vector<T> &events, const MantidVec &X, const MantidVec &Y,
                                        const MantidVec &E) {
  // Validate inputs
  if ((X.size() < 2) || (Y.size() != E.size()) || (X.size() != 1 + Y.size())) {
    std::stringstream msg;
    msg << "EventListBase::multiply() was given invalid size or "
           "inconsistent histogram arrays: X["
        << X.size() << "] "
        << "Y[" << Y.size() << " E[" << E.size() << "]";
    throw std::invalid_argument(msg.str());
  }

  size_t x_size = X.size();

  // Iterate through all events (sorted by tof)
  auto itev = findFirstEvent(events, T(X[0]));
  auto itev_end = events.end();
  // The above can still take you to end() if no events above X[0], so check
  // again.
  if (itev == itev_end)
    return;

  // Find the first bin
  size_t bin = 0;

  // Multiplier values
  double value;
  double error;
  double valueSquared;
  double errorSquared;

  // If the tof is greater the first bin boundary, so we need to find the first
  // bin
  double tof = itev->tof();
  while (bin < x_size - 1) {
    // Within range?
    if ((tof >= X[bin]) && (tof < X[bin + 1]))
      break; // Stop increasing bin
    ++bin;
  }

  // New bin! Find what you are multiplying!
  value = Y[bin];
  error = E[bin];
  valueSquared = value * value;
  errorSquared = error * error;

  // Keep going through all the events
  while ((itev != itev_end) && (bin < x_size - 1)) {
    tof = itev->tof();
    while (bin < x_size - 1) {
      // Event is Within range?
      if ((tof >= X[bin]) && (tof < X[bin + 1])) {
        // Process this event. Multiply and calculate error.
        itev->m_errorSquared =
            static_cast<float>(itev->m_errorSquared * valueSquared + errorSquared * itev->m_weight * itev->m_weight);
        itev->m_weight *= static_cast<float>(value);
        break; // out of the bin-searching-while-loop
      }
      ++bin;
      if (bin >= x_size - 1)
        break;

      // New bin! Find what you are multiplying!
      value = Y[bin];
      error = E[bin];
      valueSquared = value * value;
      errorSquared = error * error;
    }
    ++itev;
  }
}

//------------------------------------------------------------------------------------------------
/** Divide the weights in this event list by a scalar with an (optional) error.
 * The event list switches to WeightedEvent's if needed.
 * This simply calls the equivalent function: multiply(1.0/value,
 *error/(value*value)).
 *
 * @param value: divide all weights by this amount.
 * @param error: error on 'value'. Can be 0.
 * @throw std::invalid_argument if value == 0; cannot divide by zero.
 */
void divide(const double value, const double error) {
  if (value == 0.0)
    throw std::invalid_argument("EventListBase::divide() called with value of 0.0. Cannot divide by zero.");
  // Do nothing if dividing by exactly 1.0, no error
  else if (value == 1.0 && error == 0.0)
    return;

  // We'll multiply by 1/value
  double invValue = 1.0 / value;
  // Relative error remains the same
  double invError = (error / value) * invValue;

  this->multiply(invValue, invError);
}


//------------------------------------------------------------------------------------------------
/** Divide the weights in this event list by a histogram.
 * The event list switches to WeightedEvent's if needed.
 * NOTE: no unit checks are made (or possible to make) to compare the units of X
 *and tof() in the EventListBase.
 *
 * The formula used for calculating the error on the neutron weight is:
 * \f[ \sigma_{f}^2 = (A / B)^2 * (\sigma_A^2 / A^2 + \sigma_B^2 / B^2) \f]
 *
 * where:
 *  * A is the weight of the event
 *  * B is the weight of the BIN that the event falls in
 *  * \f$\sigma_A\f$ is the error (not squared) of the weight of the event
 *  * \f$\sigma_B\f$ is the error (not squared) of the bin B
 *  * f is the resulting weight of the divided event
 *
 *
 * @param X: bins of the multiplying histogram.
 * @param Y: value to multiply the weights.
 * @param E: error on the value to multiply.
 * @throw invalid_argument if the sizes of X, Y, E are not consistent.
 */
void divide(const MantidVec &X, const MantidVec &Y, const MantidVec &E) {
    // Switch to weights if needed.
    // TODO abstract out to wrapper?
    this->switchTo(WEIGHTED);
    // Fall through

  this->sortTof();
  divideHistogramHelper(*(this->events), X, Y, E);
    
}

//------------------------------------------------------------------------------------------------
/** Helper method for dividing an event list by a histogram with error
 *
 * @param events: vector of events (with weights)
 * @param X: bins of the dividing histogram.
 * @param Y: value to dividing the weights.
 * @param E: error on the value to dividing.
 * @throw invalid_argument if the sizes of X, Y, E are not consistent.
 * */
void divideHistogramHelper(std::vector<T> &events, const MantidVec &X, const MantidVec &Y,
                                      const MantidVec &E) {
  // Validate inputs
  if ((X.size() < 2) || (Y.size() != E.size()) || (X.size() != 1 + Y.size())) {
    std::stringstream msg;
    msg << "EventListBase::divide() was given invalid size or "
           "inconsistent histogram arrays: X["
        << X.size() << "] "
        << "Y[" << Y.size() << " E[" << E.size() << "]";
    throw std::invalid_argument(msg.str());
  }

  size_t x_size = X.size();

  // Iterate through all events (sorted by tof)
  auto itev = findFirstEvent(events, T(X[0]));
  auto itev_end = events.end();
  // The above can still take you to end() if no events above X[0], so check
  // again.
  if (itev == itev_end)
    return;

  // Find the first bin
  size_t bin = 0;

  // Multiplier values
  double value;
  double error;
  double valError_over_value_squared;

  // If the tof is greater the first bin boundary, so we need to find the first
  // bin
  double tof = itev->tof();
  while (bin < x_size - 1) {
    // Within range?
    if ((tof >= X[bin]) && (tof < X[bin + 1]))
      break; // Stop increasing bin
    ++bin;
  }

  // New bin! Find what you are multiplying!
  value = Y[bin];
  error = E[bin];

  // --- Division case ---
  if (value == 0) {
    value = std::numeric_limits<float>::quiet_NaN(); // Avoid divide by zero
    valError_over_value_squared = 0;
  } else
    valError_over_value_squared = error * error / (value * value);

  // Keep going through all the events
  while ((itev != events.end()) && (bin < x_size - 1)) {
    tof = itev->tof();
    while (bin < x_size - 1) {
      // Event is Within range?
      if ((tof >= X[bin]) && (tof < X[bin + 1])) {
        // Process this event. Divide and calculate error.
        double newWeight = itev->m_weight / value;
        itev->m_errorSquared = static_cast<float>(
            newWeight * newWeight *
            ((itev->m_errorSquared / (itev->m_weight * itev->m_weight)) + valError_over_value_squared));
        itev->m_weight = static_cast<float>(newWeight);
        break; // out of the bin-searching-while-loop
      }
      ++bin;
      if (bin >= x_size - 1)
        break;

      // New bin! Find what you are multiplying!
      value = Y[bin];
      error = E[bin];

      // --- Division case ---
      if (value == 0) {
        value = std::numeric_limits<float>::quiet_NaN(); // Avoid divide by zero
        valError_over_value_squared = 0;
      } else
        valError_over_value_squared = error * error / (value * value);
    }
    ++itev;
  }
}


    friend T;
    // EventListWeightErrorTofFunctionsTemplate() = default;

    inline T & as_underlying()
    {
        return static_cast<T&>(*this);
    }
};
}
}