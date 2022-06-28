// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2010 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "MantidDataObjects/EventListBase.h"
#include "MantidDataObjects/EventListTofFunctionsTemplate.h"
#include "MantidDataObjects/EventListMinusHelperFunctions.h"

namespace Mantid {
namespace DataObjects {
// Share the same definition of wrapper/interface
class DLLExport EventListTofEventNoTime : 
public EventListBase,
// extract this out to another class?
public EventListTofFunctionsTemplate<TofEventNoTime, EventListTofEventNoTime>,
public EventListBaseFunctionsTemplate<TofEventNoTime, EventListTofEventNoTime>,
// END: extract this out to another class?
public EventListPermutationsMinusHelperFunctions<TofEventNoTime> {
    friend class EventListBaseFunctionsTemplate<TofEventNoTime, EventListTofEventNoTime>;
    public:
    using EventListBaseFunctionsTemplate<TofEventNoTime, EventListTofEventNoTime>::clear;

    EventListTofEventNoTime();
    EventListTofEventNoTime(const std::vector<TofEventNoTime> &events);
    EventListTofEventNoTime(const EventList &rhs);
    ~EventListTofEventNoTime();
    bool equals(const EventList &rhs, const double tolTof, const double tolWeight,
                       const int64_t tolPulse) const  ;
    WeightedEvent getEvent(size_t event_number)  ;
    void sortTimeAtSample(const double &tofFactor, const double &tofShift, bool forceResort) const  ;
    size_t getMemorySize() const  ;
    void compressFatEvents(const double tolerance, const Mantid::Types::Core::DateAndTime &timeStart,
                                  const double seconds, EventList *destination)  ;
    void generateHistogramTimeAtSample(const MantidVec &X, MantidVec &Y, MantidVec &E, const double &tofFactor,
                                              const double &tofOffset, bool skipError)  ;
    void generateHistogramPulseTime(const MantidVec &X, MantidVec &Y, MantidVec &E, bool skipError) const  ;
    void filterByPulseTime(Types::Core::DateAndTime start, Types::Core::DateAndTime stop, EventList &output) const;
    void filterByTimeAtSample(Types::Core::DateAndTime start, Types::Core::DateAndTime stop, double tofFactor,
                                     double tofOffset, EventList &output) const  ;
    void filterInPlace(Kernel::TimeSplitterType &splitter)  ;
    void splitByTime(Kernel::TimeSplitterType &splitter, std::vector<EventList *> outputs) const  ;
    void splitByFullTime(Kernel::TimeSplitterType &splitter, std::map<int, EventList *> outputs,
                                bool docorrection, double toffactor, double tofshift) const  ;
    std::string splitByFullTimeMatrixSplitter(const std::vector<int64_t> &vec_splitters_time,
                                                     const std::vector<int> &vecgroups,
                                                     std::map<int, EventList *> vec_outputEventList, bool docorrection,
                                                     double toffactor, double tofshift) const  ;
    void splitByPulseTime(Kernel::TimeSplitterType &splitter, std::map<int, EventList *> outputs) const  ;
    void splitByPulseTimeWithMatrix(const std::vector<int64_t> &vec_times, const std::vector<int> &vec_target,
                                           std::map<int, EventList *> outputs) const  ;
    EventListTofEventNoTime &operator-=(const EventList &more_events);
    void sortPulseTimeTOFDelta(const Types::Core::DateAndTime &start, const double seconds) const;
    void sortPulseTime() const ;
    void sortPulseTimeTOF() const ;


    protected:
     /// List of Events
    mutable std::vector<TofEventNoTime> events;
};

} // namespace DataObjects
} // namespace Mantid