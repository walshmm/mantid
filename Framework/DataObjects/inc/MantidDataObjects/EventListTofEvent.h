// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2010 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "MantidDataObjects/EventListBase.h"
#include "MantidDataObjects/EventListWeightErrorPulsetimeTofFunctionsTemplate.h"

namespace Mantid {
namespace DataObjects {
// Share the same definition of wrapper/interface
class DLLExport EventListTofEvent : 
public EventListBase,
// extract this out to another class?
public EventListWeightErrorPulsetimeTofFunctionsTemplate<Types::Event::TofEvent, EventListTofEvent>,
public EventListWeightErrorTofFunctionsTemplate<Types::Event::TofEvent, EventListTofEvent>,
public EventListWeightErrorFunctionsTemplate<Types::Event::TofEvent, EventListTofEvent>,
public EventListPulsetimeTofFunctionsTemplate<Types::Event::TofEvent, EventListTofEvent>,
public EventListWeightFunctionsTemplate<Types::Event::TofEvent, EventListTofEvent>,
public EventListErrorFunctionsTemplate<Types::Event::TofEvent, EventListTofEvent>,
public EventListTofFunctionsTemplate<Types::Event::TofEvent, EventListTofEvent>,
public EventListPulsetimeFunctionsTemplate<Types::Event::TofEvent, EventListTofEvent>,
public EventListBaseFunctionsTemplate<Types::Event::TofEvent, EventListTofEvent>
// END: extract this out to another class?
{
    friend class EventListBaseFunctionsTemplate<Types::Event::TofEvent, EventListTofEvent>;
    public:
    using EventListBaseFunctionsTemplate<Types::Event::TofEvent, EventListTofEvent>::clear;

    EventListTofEvent(const std::vector<TofEvent> &events);
    EventListTofEvent(EventWorkspaceMRU *mru, specnum_t specNo);
    EventListTofEvent();
    EventListTofEvent(const EventList &rhs) ;
    ~EventListTofEvent() ;
    bool equals(const EventList &rhs, const double tolTof, const double tolWeight,
                       const int64_t tolPulse) const  ;
    WeightedEvent getEvent(size_t event_number)  ;
    void sortTimeAtSample(const double &tofFactor, const double &tofShift, bool forceResort) const  ;
    size_t getMemorySize() const  ;
    void generateHistogramTimeAtSample(const MantidVec &X, MantidVec &Y, MantidVec &E, const double &tofFactor,
                                              const double &tofOffset, bool skipError)  ;
    void generateHistogramPulseTime(const MantidVec &X, MantidVec &Y, MantidVec &E, bool skipError) const  ;
    EventListTofEvent &operator-=(const EventList &more_events);
    void generateCountsHistogramPulseTime(const MantidVec &X, MantidVec &Y) const; 
    void generateCountsHistogramPulseTime(const double &xMin, const double &xMax, MantidVec &Y,
                                                 const double TOF_min, const double TOF_max) const;
    void generateCountsHistogramTimeAtSample(const MantidVec &X, MantidVec &Y, const double &tofFactor,
                                                    const double &tofOffset) const;
    void generateCountsHistogram(const MantidVec &X, MantidVec &Y) const;
    protected:
    /// List of Events
    mutable std::vector<Types::Event::TofEvent> events;
    // Declared in EventListBaseFunctionsTemplate
};

} // namespace DataObjects
} // namespace Mantid