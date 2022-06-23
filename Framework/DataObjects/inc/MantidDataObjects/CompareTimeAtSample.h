#pragma once

const double SEC_TO_NANO = 1.e9;

/**
 * Type for comparing events in terms of time at sample
 */
template <typename EventType> class CompareTimeAtSample {
private:
  const double m_tofFactor;
  const double m_tofShift;

public:
  CompareTimeAtSample(const double tofFactor, const double tofShift) : m_tofFactor(tofFactor), m_tofShift(tofShift) {}

  /**
   * Compare two events based on the time they arrived at the sample.
   * Coefficient is used to provide scaling.
   * For elastic scattering coefficient is L1 / (L1 + L2)
   * @param e1 :: first event to compare
   * @param e2 :: second event to compare
   * @param coefficient :: scaling coefficient
   * @return True if first event evaluates to be < second event, otherwise false
   */
  bool operator()(const EventType &e1, const EventType &e2) const {
    const auto tAtSample1 = calculateCorrectedFullTime(e1, m_tofFactor, m_tofShift);
    const auto tAtSample2 = calculateCorrectedFullTime(e2, m_tofFactor, m_tofShift);
    return (tAtSample1 < tAtSample2);
  }

     /**
     * Calculate the corrected full time in nanoseconds
     * @param event : The event with pulse time and time-of-flight
     * @param tofFactor : Time of flight coefficient factor
     * @param tofShift : Tof shift in seconds
     * @return Corrected full time at sample in Nanoseconds.
     */
    template <typename T>
    int64_t calculateCorrectedFullTime(const T &event, const double tofFactor, const double tofShift) {
    return event.pulseTime().totalNanoseconds() +
            static_cast<int64_t>(tofFactor * (event.tof() * 1.0E3) + (tofShift * 1.0E9));
    }
};