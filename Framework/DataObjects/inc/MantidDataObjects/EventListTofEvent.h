// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2010 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "MantidDataObjects/EventListBase.h"

namespace Mantid {
namespace DataObjects {
// Share the same definition of wrapper/interface
class DLLExport EventListTofEvent : public EventListBase {
    bool equals(const EventListBase &rhs, const double tolTof, const double tolWeight,
                       const int64_t tolPulse) const override;
    WeightedEvent getEvent(size_t event_number) override;
    void sortTimeAtSample(const double &tofFactor, const double &tofShift, bool forceResort) const override;
    size_t getMemorySize() const override;
    void generateHistogramTimeAtSample(const MantidVec &X, MantidVec &Y, MantidVec &E, const double &tofFactor,
                                              const double &tofOffset, bool skipError) override;
    void generateHistogramPulseTime(const MantidVec &X, MantidVec &Y, MantidVec &E, bool skipError) const override;
};

} // namespace DataObjects
} // namespace Mantid