// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2010 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "MantidDataObjects/EventList.h"
#include "MantidAPI/IEventList.h"
#include "MantidDataObjects/Events.h"
#include "MantidKernel/MultiThreaded.h"
#include "MantidKernel/System.h"
#include "MantidKernel/cow_ptr.h"
#include "MantidDataObjects/CompareTimeAtSample.h"
#include <iosfwd>
#include <memory>
#include <vector>


// Children:
// * EventListWeightedEvent
// * EventListTofEvent
// * EventListTofEventNoTime
// * EventListWeightedEventNoTime

namespace Mantid {
namespace DataObjects {

class EventList;    
class EventListTofEvent;
class EventListWeightedEvent;
class EventListWeightedEventNoTime;
class EventListTofEventNoTime;
template<typename T, typename SELF>
class EventListBaseFunctionsTemplate;
// Share the same definition of wrapper/interface
class DLLExport EventListBase : public EventList {

friend class EventList;
friend class EventListBaseFunctionsTemplate<Types::Event::TofEvent,EventListTofEvent>;
friend class EventListBaseFunctionsTemplate<WeightedEvent,EventListWeightedEvent>;
friend class EventListBaseFunctionsTemplate<WeightedEventNoTime,EventListWeightedEventNoTime>;
friend class EventListBaseFunctionsTemplate<TofEventNoTime,EventListTofEventNoTime>;
friend class EventListTofEvent;
friend class EventListWeightedEvent;
friend class EventListWeightedEventNoTime;
friend class EventListTofEventNoTime;

public:
  EventListBase();

  EventListBase(EventWorkspaceMRU *mru, specnum_t specNo);

  EventListBase(const EventList &rhs);

  EventListBase(const std::vector<Types::Event::TofEvent> &events);

  EventListBase(const std::vector<WeightedEvent> &events);

  EventListBase(const std::vector<WeightedEventNoTime> &events);

  EventListBase(const std::vector<TofEventNoTime> &events);

  ~EventListBase() override;

  void copyDataFrom(const ISpectrum &source) override;

  void createFromHistogram(const ISpectrum *inSpec, bool GenerateZeros, bool GenerateMultipleEvents,
                           int MaxEventsPerBin);

  EventListBase &operator=(const EventListBase &);

  EventListBase &operator+=(const Types::Event::TofEvent &event);

  EventListBase &operator+=(const std::vector<Types::Event::TofEvent> &more_events);

  EventListBase &operator+=(const WeightedEvent &event);

  EventListBase &operator+=(const std::vector<WeightedEvent> &more_events);

  EventListBase &operator+=(const std::vector<WeightedEventNoTime> &more_events);

  EventListBase &operator+=(const EventListBase &more_events);

  EventListBase &operator-=(const EventListBase &more_events);

  EventListBase &operator*=(const double value);

  EventListBase &operator/=(const double value);

  bool operator==(const EventListBase &rhs) const;
  bool operator!=(const EventListBase &rhs) const;
  bool equals(const EventListBase &rhs, const double tolTof, const double tolWeight, const int64_t tolPulse) const;

  virtual Mantid::API::EventType getEventType() const override;

  virtual void switchTo(Mantid::API::EventType newType) override;

  virtual  WeightedEvent getEvent(size_t event_number);

  virtual std::vector<Types::Event::TofEvent> &getEvents();
  virtual const std::vector<Types::Event::TofEvent> &getEvents() const;

//   std::vector<Types::Event::Event> &getEventsTyped();
//   const std::vector<Types::Event::Event> &getEventsTyped() const;

  virtual std::vector<WeightedEvent> &getWeightedEvents();
  virtual const std::vector<WeightedEvent> &getWeightedEvents() const;

  virtual std::vector<WeightedEventNoTime> &getWeightedEventsNoTime();
  virtual const std::vector<WeightedEventNoTime> &getWeightedEventsNoTime() const;

  virtual std::vector<TofEventNoTime> &getEventsNoTime();
  virtual const std::vector<TofEventNoTime> &getEventsNoTime() const;

  virtual void clear(const bool removeDetIDs = true) override;
  virtual void clearUnused();

  virtual void setMRU(EventWorkspaceMRU *newMRU);

  virtual void clearData() override;

  virtual void reserve(size_t num) override;

  virtual void sort(const EventSortType order) const;

  virtual void setSortOrder(const EventSortType order) const;

  virtual void sortTof() const;

  virtual void sortPulseTime() const;
  virtual void sortPulseTimeTOF() const;
  virtual void sortTimeAtSample(const double &tofFactor, const double &tofShift, bool forceResort = false) const;

  virtual bool isSortedByTof() const override;

  virtual EventSortType getSortType() const;

  // X-vector accessors. These reset the MRU for this spectrum
  virtual void setX(const Kernel::cow_ptr<HistogramData::HistogramX> &X) override;
  virtual MantidVec &dataX() override;
  virtual const MantidVec &dataX() const override;
  virtual const MantidVec &readX() const override;
  virtual Kernel::cow_ptr<HistogramData::HistogramX> ptrX() const override;

  virtual MantidVec &dataDx() override;
  virtual const MantidVec &dataDx() const override;
  virtual const MantidVec &readDx() const override;

  /// Deprecated, use mutableY() instead. Disallowed data accessors - can't
  /// modify Y/E on a EventList
  virtual MantidVec &dataY() override { throw std::runtime_error("EventList: non-const access to Y data is not possible."); }
  /// Deprecated, use mutableE() instead. Disallowed data accessors - can't
  /// modify Y/E on a EventList
  virtual MantidVec &dataE() override { throw std::runtime_error("EventList: non-const access to E data is not possible."); }

  // Allowed data accessors - read-only Y/E histogram VIEWS of an event list
  /// Deprecated, use y() instead. Return a read-only Y histogram view of an
  /// event list
  virtual const MantidVec &dataY() const override;
  /// Deprecated, use e() instead. Return a read-only E histogram view of an
  /// event list
  virtual const MantidVec &dataE() const override;

  virtual MantidVec *makeDataY() const;
  virtual MantidVec *makeDataE() const;

  virtual std::size_t getNumberEvents() const override;
  virtual bool empty() const;

  virtual size_t getMemorySize() const override;

  virtual size_t histogram_size() const;

  virtual void compressEvents(double tolerance, EventList *destination);
  virtual void compressFatEvents(const double tolerance, const Types::Core::DateAndTime &timeStart, const double seconds,
                         EventList *destination);
  // get EventType declaration
  virtual void generateHistogram(const MantidVec &X, MantidVec &Y, MantidVec &E, bool skipError = false) const override;
  virtual void generateHistogramPulseTime(const MantidVec &X, MantidVec &Y, MantidVec &E,
                                  bool skipError = false) const override;

  virtual void generateHistogramTimeAtSample(const MantidVec &X, MantidVec &Y, MantidVec &E, const double &tofFactor,
                                     const double &tofOffset, bool skipError = false) const override;

  virtual void integrate(const double minX, const double maxX, const bool entireRange, double &sum, double &error) const;

  double integrate(const double minX, const double maxX, const bool entireRange) const override;

  virtual void convertTof(std::function<double(double)> func, const int sorting = 0) override;

  virtual void convertTof(const double factor, const double offset = 0.) override;

  virtual void scaleTof(const double factor) override;

  virtual void addTof(const double offset) override;

  virtual void addPulsetime(const double seconds) override;

  virtual void addPulsetimes(const std::vector<double> &seconds) override;

  virtual void maskTof(const double tofMin, const double tofMax) override;
  virtual void maskCondition(const std::vector<bool> &mask) override;

  virtual void getTofs(std::vector<double> &tofs) const override;
  virtual double getTofMin() const override;
  virtual double getTofMax() const override;
  virtual Mantid::Types::Core::DateAndTime getPulseTimeMax() const override;
  virtual Mantid::Types::Core::DateAndTime getPulseTimeMin() const override;
  virtual void getPulseTimeMinMax(Mantid::Types::Core::DateAndTime &tMin, Mantid::Types::Core::DateAndTime &tM) const;
  virtual Mantid::Types::Core::DateAndTime getTimeAtSampleMax(const double &tofFactor, const double &tofOffset) const override;
  virtual Mantid::Types::Core::DateAndTime getTimeAtSampleMin(const double &tofFactor, const double &tofOffset) const override;

  virtual std::vector<double> getTofs() const override;

  /// Return the list of event weight  values
  virtual std::vector<double> getWeights() const override;
  /// Return the list of event weight values
  virtual void getWeights(std::vector<double> &weights) const override;

  /// Return the list of event weight  error values
  virtual std::vector<double> getWeightErrors() const override;
  /// Return the list of event weight error values
  virtual void getWeightErrors(std::vector<double> &weightErrors) const override;

  virtual std::vector<Mantid::Types::Core::DateAndTime> getPulseTimes() const override;

  virtual void setTofs(const MantidVec &tofs) override;

  virtual void reverse();

  virtual void filterByPulseTime(Types::Core::DateAndTime start, Types::Core::DateAndTime stop, EventList &output) const;

  virtual void filterByTimeAtSample(Types::Core::DateAndTime start, Types::Core::DateAndTime stop, double tofFactor,
                            double tofOffset, EventList &output) const;

  virtual void filterInPlace(Kernel::TimeSplitterType &splitter);

  virtual void splitByTime(Kernel::TimeSplitterType &splitter, std::vector<EventList *> outputs) const;

  virtual void splitByFullTime(Kernel::TimeSplitterType &splitter, std::map<int, EventList *> outputs, bool docorrection,
                       double toffactor, double tofshift) const;

  /// Split ...
  virtual std::string splitByFullTimeMatrixSplitter(const std::vector<int64_t> &vec_splitters_time,
                                            const std::vector<int> &vecgroups,
                                            std::map<int, EventList *> vec_outputEventList, bool docorrection,
                                            double toffactor, double tofshift) const;

  /// Split events by pulse time
  virtual void splitByPulseTime(Kernel::TimeSplitterType &splitter, std::map<int, EventList *> outputs) const;

  /// Split events by pulse time with Matrix splitters
  virtual void splitByPulseTimeWithMatrix(const std::vector<int64_t> &vec_times, const std::vector<int> &vec_target,
                                  std::map<int, EventList *> outputs) const;

  virtual void multiply(const double value, const double error = 0.0) override;

  virtual void multiply(const MantidVec &X, const MantidVec &Y, const MantidVec &E) override;

  virtual void divide(const double value, const double error = 0.0) override;

  virtual void divide(const MantidVec &X, const MantidVec &Y, const MantidVec &E) override;

  virtual void convertUnitsViaTof(Mantid::Kernel::Unit *fromUnit, Mantid::Kernel::Unit *toUnit);
  virtual void convertUnitsQuickly(const double &factor, const double &power);

  /// Returns the Histogram associated with this spectrum. Y and E data is
  /// computed from the event list.
  virtual HistogramData::Histogram histogram() const override;
  virtual HistogramData::Counts counts() const override;
  virtual HistogramData::CountVariances countVariances() const override;
  virtual HistogramData::CountStandardDeviations countStandardDeviations() const override;
  virtual HistogramData::Frequencies frequencies() const override;
  virtual HistogramData::FrequencyVariances frequencyVariances() const override;
  virtual HistogramData::FrequencyStandardDeviations frequencyStandardDeviations() const override;
  virtual const HistogramData::HistogramY &y() const override;
  virtual const HistogramData::HistogramE &e() const override;
  virtual Kernel::cow_ptr<HistogramData::HistogramY> sharedY() const override;
  virtual Kernel::cow_ptr<HistogramData::HistogramE> sharedE() const override;

  virtual void generateCountsHistogramPulseTime(const double &xMin, const double &xMax, MantidVec &Y,
                                        const double TofMin = std::numeric_limits<double>::lowest(),
                                        const double TofMax = std::numeric_limits<double>::max()) const;

protected:
  virtual void checkAndSanitizeHistogram(HistogramData::Histogram &histogram) override;
  virtual void checkWorksWithPoints() const override;
  virtual void checkIsYAndEWritable() const override;



    /// Histogram object holding the histogram data. Currently only X.
    HistogramData::Histogram m_histogram;

    /// What type of event is in our list.
    Mantid::API::EventType eventType;

    /// Last sorting order
    mutable EventSortType order;

    /// MRU lists of the parent EventWorkspace
    mutable EventWorkspaceMRU *mru;

    /// Mutex that is locked while sorting an event list
    mutable std::mutex m_sortMutex;
    
private:

    //Should be a shared pointer that gets passed down 
    //to parents so they all point to the same list
    // std::vector<T> events;


  using ISpectrum::copyDataInto;
  void copyDataInto(EventListBase &sink) const;
  void copyDataInto(Histogram1D &sink) const;

  const HistogramData::Histogram &histogramRef() const override { throw std::runtime_error("unimplemented, use a derived class"); }
  HistogramData::Histogram &mutableHistogramRef() override;



  void generateCountsHistogram(const MantidVec &X, MantidVec &Y) const;

  void generateCountsHistogramPulseTime(const MantidVec &X, MantidVec &Y) const;

  void generateCountsHistogramTimeAtSample(const MantidVec &X, MantidVec &Y, const double &tofFactor,
                                           const double &tofOffset) const;

  void generateErrorsHistogram(const MantidVec &Y, MantidVec &E) const;

  void switchToWeightedEvents();
  void switchToWeightedEventsNoTime();
  // should not be called externally
  void sortPulseTimeTOFDelta(const Types::Core::DateAndTime &start, const double seconds) const;

};
// Methods overloaded to get event vectors.
} // namespace DataObjects
} // namespace Mantid