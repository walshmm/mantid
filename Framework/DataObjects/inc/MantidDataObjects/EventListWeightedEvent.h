// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2010 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "MantidDataObjects/EventListBase.h"
#include "MantidDataObjects/EventListWeightErrorPulsetimeTofFunctionsTemplate.h"
#include "MantidDataObjects/EventListMinusHelperFunctions.h"

namespace Mantid {
namespace DataObjects {
// Share the same definition of wrapper/interface
class DLLExport EventListWeightedEvent : 
public EventListBase,
// extract this out to another class?
public EventListWeightErrorPulsetimeTofFunctionsTemplate<WeightedEvent, EventListWeightedEvent>,
public EventListWeightErrorTofFunctionsTemplate<WeightedEvent, EventListWeightedEvent>,
public EventListWeightErrorFunctionsTemplate<WeightedEvent, EventListWeightedEvent>,
public EventListPulsetimeTofFunctionsTemplate<WeightedEvent, EventListWeightedEvent>,
public EventListWeightFunctionsTemplate<WeightedEvent, EventListWeightedEvent>,
public EventListErrorFunctionsTemplate<WeightedEvent, EventListWeightedEvent>,
public EventListTofFunctionsTemplate<WeightedEvent, EventListWeightedEvent>,
public EventListPulsetimeFunctionsTemplate<WeightedEvent, EventListWeightedEvent>,
public EventListBaseFunctionsTemplate<WeightedEvent, EventListWeightedEvent>,
// END: extract this out to another class?
public EventListPermutationsMinusHelperFunctions<WeightedEvent> {

    friend class EventListBaseFunctionsTemplate<WeightedEvent, EventListWeightedEvent>;
    public:
    using EventListBaseFunctionsTemplate<WeightedEvent, EventListWeightedEvent>::clear;
    using EventListBaseFunctionsTemplate<WeightedEvent, EventListWeightedEvent>::getNumberEvents;
    using EventListBaseFunctionsTemplate<WeightedEvent, EventListWeightedEvent>::sortTof;
    using EventListWeightFunctionsTemplate<WeightedEvent, EventListWeightedEvent>::getWeightsHelper;

    EventListWeightedEvent(const std::vector<WeightedEvent> &events);
    EventListWeightedEvent(const EventList &rhs);
    EventListWeightedEvent();
    ~EventListWeightedEvent();
    bool equals(const EventList &rhs, const double tolTof, const double tolWeight,
                       const int64_t tolPulse) const  ;
    WeightedEvent getEvent(size_t event_number)  ;
    void sortTimeAtSample(const double &tofFactor, const double &tofShift, bool forceResort) const  ;
    size_t getMemorySize() const  ;
    void generateHistogramTimeAtSample(const MantidVec &X, MantidVec &Y, MantidVec &E, const double &tofFactor,
                                              const double &tofOffset, bool skipError)  ;
    void generateHistogramPulseTime(const MantidVec &X, MantidVec &Y, MantidVec &E, bool skipError) const  ;
    void generateHistogram(const MantidVec &X, MantidVec &Y, MantidVec &E, bool skipError) const  ;
    void getWeights(std::vector<double> &weights) const  ;
    void getWeightErrors(std::vector<double> &weightErrors) const  ;
    EventListWeightedEvent &operator-=(const EventList &more_events);

    protected:
     /// List of Events
    mutable std::vector<WeightedEvent> events;
};

} // namespace DataObjects
} // namespace Mantid