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
class DLLExport EventListWeightedEventNoTime : public EventListBase {
    bool equals(const EventListBase &rhs, const double tolTof, const double tolWeight,
                       const int64_t tolPulse) const override;
    WeightedEvent getEvent(size_t event_number) override;
    void sortTimeAtSample(const double &tofFactor, const double &tofShift, bool forceResort) const override;
    void sortPulseTime() const override;
    void sortPulseTimeTOF() const override;
    void sortPulseTimeTOFDelta(const Types::Core::DateAndTime &start, const double seconds) const override;
    size_t getMemorySize() const override;
    void compressEvents(double tolerance, EventListBase *destination) override;
    void compressFatEvents(const double tolerance, const Mantid::Types::Core::DateAndTime &timeStart,
                                  const double seconds, EventListBase *destination) override;
    void generateHistogramTimeAtSample(const MantidVec &X, MantidVec &Y, MantidVec &E, const double &tofFactor,
                                              const double &tofOffset, bool skipError) override;                             
    void generateHistogramPulseTime(const MantidVec &X, MantidVec &Y, MantidVec &E, bool skipError) const override;
    void generateHistogram(const MantidVec &X, MantidVec &Y, MantidVec &E, bool skipError) const override;
    void addPulsetime(const double seconds) override;
    void addPulsetimes(const std::vector<double> &seconds) override;
    void getWeights(std::vector<double> &weights) const override;
    void getWeightErrors(std::vector<double> &weightErrors) const override;
    void filterByPulseTime(DateAndTime start, DateAndTime stop, EventListBase &output) const override;
    void filterByTimeAtSample(Types::Core::DateAndTime start, Types::Core::DateAndTime stop, double tofFactor,
                                     double tofOffset, EventListBase &output) const override;
    void filterInPlace(Kernel::TimeSplitterType &splitter) override;
    void splitByTime(Kernel::TimeSplitterType &splitter, std::vector<EventListBase *> outputs) const override;
    void splitByFullTime(Kernel::TimeSplitterType &splitter, std::map<int, EventListBase *> outputs,
                                bool docorrection, double toffactor, double tofshift) const override;
    std::string splitByFullTimeMatrixSplitter(const std::vector<int64_t> &vec_splitters_time,
                                                     const std::vector<int> &vecgroups,
                                                     std::map<int, EventListBase *> vec_outputEventList, bool docorrection,
                                                     double toffactor, double tofshift) const override;
};

} // namespace DataObjects
} // namespace Mantid