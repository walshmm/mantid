// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2010 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "MantidAPI/IEventList.h"
#include "MantidDataObjects/Events.h"
#include "MantidKernel/MultiThreaded.h"
#include "MantidKernel/System.h"
#include "MantidKernel/cow_ptr.h"
#include <iosfwd>
#include <vector>

namespace Mantid {
namespace Types {
namespace Core {
class DateAndTime;
}
} // namespace Types
namespace Kernel {
class SplittingInterval;
using TimeSplitterType = std::vector<SplittingInterval>;
class Unit;
} // namespace Kernel
namespace DataObjects {
class EventWorkspaceMRU;
class EventListBase;
/// How the event list is sorted.
enum EventSortType {
  UNSORTED,
  TOF_SORT,
  PULSETIME_SORT,
  PULSETIMETOF_SORT,
  PULSETIMETOF_DELTA_SORT,
  TIMEATSAMPLE_SORT
};

//==========================================================================================
/** @class Mantid::DataObjects::EventList

    A class for holding :
      - a list of neutron detection events (TofEvent or WeightedEvent).
      - a list of associated detector ID's.

    This class can switch from holding regular TofEvent's (implied weight of
   1.0)
    or WeightedEvent (where each neutron can have a non-1 weight).
    This is done transparently.

    @author Janik Zikovsky, SNS ORNL
    @date 4/02/2010
*/

class DLLExport EventList : public Mantid::API::IEventList {
public:
  EventList();

  EventList(EventWorkspaceMRU *mru, specnum_t specNo);

  EventList(const EventList &rhs);

  EventList(const std::vector<Types::Event::TofEvent> &events);

  EventList(const std::vector<WeightedEvent> &events);

  EventList(const std::vector<WeightedEventNoTime> &events);

  EventList(const std::vector<TofEventNoTime> &events);

  ~EventList() override;

  void copyDataFrom(const ISpectrum &source) override;

  void createFromHistogram(const ISpectrum *inSpec, bool GenerateZeros, bool GenerateMultipleEvents,
                           int MaxEventsPerBin);

  EventList &operator=(const EventList &);

  EventList &operator+=(const Types::Event::TofEvent &event);

  EventList &operator+=(const std::vector<Types::Event::TofEvent> &more_events);

  EventList &operator+=(const WeightedEvent &event);

  EventList &operator+=(const std::vector<WeightedEvent> &more_events);

  EventList &operator+=(const std::vector<WeightedEventNoTime> &more_events);

  EventList &operator+=(const EventList &more_events);

  EventList &operator-=(const EventList &more_events);

  bool operator==(const EventList &rhs) const;
  bool operator!=(const EventList &rhs) const;
  bool equals(const EventList &rhs, const double tolTof, const double tolWeight, const int64_t tolPulse) const;

  // --------------------------------------------------------------------------
  /** Append an event to the histogram, without clearing the cache, to make it
   *faster.
   * NOTE: Only call this on a un-weighted event list!
   *
   * @param event :: TofEvent to add at the end of the list.
   * */
  inline void addEventQuickly(const Types::Event::TofEvent &event) {
    this->eventList->addEventQuickly(event);
  }

  // --------------------------------------------------------------------------
  /** Append an event to the histogram, without clearing the cache, to make it
   * faster.
   * @param event :: WeightedEvent to add at the end of the list.
   * */
  inline void addEventQuickly(const WeightedEvent &event) {
    this->eventList->addEventQuickly(event);
  }

  // --------------------------------------------------------------------------
  /** Append an event to the histogram, without clearing the cache, to make it
   * faster.
   * @param event :: WeightedEventNoTime to add at the end of the list.
   * */
  inline void addEventQuickly(const WeightedEventNoTime &event) {
    this->eventList->addEventQuickly(event);
  }

  Mantid::API::EventType getEventType() const override;

  void switchTo(Mantid::API::EventType newType) override;

  WeightedEvent getEvent(size_t event_number);

  std::vector<Types::Event::TofEvent> &getEvents();
  const std::vector<Types::Event::TofEvent> &getEvents() const;

  std::vector<Types::Event::Event> &getEventsTyped();
  const std::vector<Types::Event::Event> &getEventsTyped() const;

  std::vector<WeightedEvent> &getWeightedEvents();
  const std::vector<WeightedEvent> &getWeightedEvents() const;

  std::vector<WeightedEventNoTime> &getWeightedEventsNoTime();
  const std::vector<WeightedEventNoTime> &getWeightedEventsNoTime() const;

  void clear(const bool removeDetIDs = true) override;
  void clearUnused();

  void setMRU(EventWorkspaceMRU *newMRU);

  void clearData() override;

  void reserve(size_t num) override;

  void sort(const EventSortType order) const;

  void setSortOrder(const EventSortType order) const;

  void sortTof() const;

  void sortPulseTime() const;
  void sortPulseTimeTOF() const;
  void sortTimeAtSample(const double &tofFactor, const double &tofShift, bool forceResort = false) const;

  bool isSortedByTof() const override;

  EventSortType getSortType() const;

  // X-vector accessors. These reset the MRU for this spectrum
  void setX(const Kernel::cow_ptr<HistogramData::HistogramX> &X) override;
  MantidVec &dataX() override;
  const MantidVec &dataX() const override;
  const MantidVec &readX() const override;
  Kernel::cow_ptr<HistogramData::HistogramX> ptrX() const override;

  MantidVec &dataDx() override;
  const MantidVec &dataDx() const override;
  const MantidVec &readDx() const override;

  /// Deprecated, use mutableY() instead. Disallowed data accessors - can't
  /// modify Y/E on a EventList
  MantidVec &dataY() override { throw std::runtime_error("EventList: non-const access to Y data is not possible."); }
  /// Deprecated, use mutableE() instead. Disallowed data accessors - can't
  /// modify Y/E on a EventList
  MantidVec &dataE() override { throw std::runtime_error("EventList: non-const access to E data is not possible."); }

  // Allowed data accessors - read-only Y/E histogram VIEWS of an event list
  /// Deprecated, use y() instead. Return a read-only Y histogram view of an
  /// event list
  const MantidVec &dataY() const override;
  /// Deprecated, use e() instead. Return a read-only E histogram view of an
  /// event list
  const MantidVec &dataE() const override;

  MantidVec *makeDataY() const;
  MantidVec *makeDataE() const;

  std::size_t getNumberEvents() const override;
  bool empty() const;

  size_t getMemorySize() const override;

  virtual size_t histogram_size() const;

  void compressEvents(double tolerance, EventList *destination);
  void compressFatEvents(const double tolerance, const Types::Core::DateAndTime &timeStart, const double seconds,
                         EventList *destination);
  // get EventType declaration
  void generateHistogram(const MantidVec &X, MantidVec &Y, MantidVec &E, bool skipError = false) const override;
  void generateHistogramPulseTime(const MantidVec &X, MantidVec &Y, MantidVec &E,
                                  bool skipError = false) const override;

  void generateHistogramTimeAtSample(const MantidVec &X, MantidVec &Y, MantidVec &E, const double &tofFactor,
                                     const double &tofOffset, bool skipError = false) const override;

  void integrate(const double minX, const double maxX, const bool entireRange, double &sum, double &error) const;

  double integrate(const double minX, const double maxX, const bool entireRange) const override;

  void convertTof(std::function<double(double)> func, const int sorting = 0) override;

  void convertTof(const double factor, const double offset = 0.) override;

  void scaleTof(const double factor) override;

  void addTof(const double offset) override;

  void addPulsetime(const double seconds) override;

  void addPulsetimes(const std::vector<double> &seconds) override;

  void maskTof(const double tofMin, const double tofMax) override;
  void maskCondition(const std::vector<bool> &mask) override;

  void getTofs(std::vector<double> &tofs) const override;
  double getTofMin() const override;
  double getTofMax() const override;
  Mantid::Types::Core::DateAndTime getPulseTimeMax() const override;
  Mantid::Types::Core::DateAndTime getPulseTimeMin() const override;
  void getPulseTimeMinMax(Mantid::Types::Core::DateAndTime &tMin, Mantid::Types::Core::DateAndTime &tM) const;
  Mantid::Types::Core::DateAndTime getTimeAtSampleMax(const double &tofFactor, const double &tofOffset) const override;
  Mantid::Types::Core::DateAndTime getTimeAtSampleMin(const double &tofFactor, const double &tofOffset) const override;

  std::vector<double> getTofs() const override;

  /// Return the list of event weight  values
  std::vector<double> getWeights() const override;
  /// Return the list of event weight values
  void getWeights(std::vector<double> &weights) const override;

  /// Return the list of event weight  error values
  std::vector<double> getWeightErrors() const override;
  /// Return the list of event weight error values
  void getWeightErrors(std::vector<double> &weightErrors) const override;

  std::vector<Mantid::Types::Core::DateAndTime> getPulseTimes() const override;

  void setTofs(const MantidVec &tofs) override;

  void reverse();

  void filterByPulseTime(Types::Core::DateAndTime start, Types::Core::DateAndTime stop, EventList &output) const;

  void filterByTimeAtSample(Types::Core::DateAndTime start, Types::Core::DateAndTime stop, double tofFactor,
                            double tofOffset, EventList &output) const;

  void filterInPlace(Kernel::TimeSplitterType &splitter);

  void splitByTime(Kernel::TimeSplitterType &splitter, std::vector<EventList *> outputs) const;

  void splitByFullTime(Kernel::TimeSplitterType &splitter, std::map<int, EventList *> outputs, bool docorrection,
                       double toffactor, double tofshift) const;

  /// Split ...
  std::string splitByFullTimeMatrixSplitter(const std::vector<int64_t> &vec_splitters_time,
                                            const std::vector<int> &vecgroups,
                                            std::map<int, EventList *> vec_outputEventList, bool docorrection,
                                            double toffactor, double tofshift) const;

  /// Split events by pulse time
  void splitByPulseTime(Kernel::TimeSplitterType &splitter, std::map<int, EventList *> outputs) const;

  /// Split events by pulse time with Matrix splitters
  void splitByPulseTimeWithMatrix(const std::vector<int64_t> &vec_times, const std::vector<int> &vec_target,
                                  std::map<int, EventList *> outputs) const;

  void multiply(const double value, const double error = 0.0) override;
  EventList &operator*=(const double value);

  void multiply(const MantidVec &X, const MantidVec &Y, const MantidVec &E) override;

  void divide(const double value, const double error = 0.0) override;
  EventList &operator/=(const double value);

  void divide(const MantidVec &X, const MantidVec &Y, const MantidVec &E) override;

  void convertUnitsViaTof(Mantid::Kernel::Unit *fromUnit, Mantid::Kernel::Unit *toUnit);
  void convertUnitsQuickly(const double &factor, const double &power);

  /// Returns the Histogram associated with this spectrum. Y and E data is
  /// computed from the event list.
  HistogramData::Histogram histogram() const override;
  HistogramData::Counts counts() const override;
  HistogramData::CountVariances countVariances() const override;
  HistogramData::CountStandardDeviations countStandardDeviations() const override;
  HistogramData::Frequencies frequencies() const override;
  HistogramData::FrequencyVariances frequencyVariances() const override;
  HistogramData::FrequencyStandardDeviations frequencyStandardDeviations() const override;
  const HistogramData::HistogramY &y() const override;
  const HistogramData::HistogramE &e() const override;
  Kernel::cow_ptr<HistogramData::HistogramY> sharedY() const override;
  Kernel::cow_ptr<HistogramData::HistogramE> sharedE() const override;

  void generateCountsHistogramPulseTime(const double &xMin, const double &xMax, MantidVec &Y,
                                        const double TofMin = std::numeric_limits<double>::lowest(),
                                        const double TofMax = std::numeric_limits<double>::max()) const;
mutable std::shared_ptr<EventList> eventList;
protected:
  void checkAndSanitizeHistogram(HistogramData::Histogram &histogram) override;
  void checkWorksWithPoints() const override;
  void checkIsYAndEWritable() const override;

  

private:
  const HistogramData::Histogram &histogramRef() const override { throw std::runtime_error("Not Implemented"); };
  HistogramData::Histogram &mutableHistogramRef() override { throw std::runtime_error("Not Implemented"); };
 
  void switchToUnweightedEvents();

  void switchToWeightedEvents();

  void switchToWeightedEventsNoTime();

};

// Methods overloaded to get event vectors.
DLLExport void getEventsFrom(EventList &el, std::vector<Types::Event::TofEvent> *&events);
DLLExport void getEventsFrom(const EventList &el, std::vector<Types::Event::TofEvent> const *&events);
DLLExport void getEventsFrom(EventList &el, std::vector<WeightedEvent> *&events);
DLLExport void getEventsFrom(const EventList &el, std::vector<WeightedEvent> const *&events);
DLLExport void getEventsFrom(EventList &el, std::vector<WeightedEventNoTime> *&events);
DLLExport void getEventsFrom(const EventList &el, std::vector<WeightedEventNoTime> const *&events);

} // namespace DataObjects
} // namespace Mantid