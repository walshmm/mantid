// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2010 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "MantidDataObjects/EventList.h"


// Children:
// * EventListWeightedEvent
// * EventListTofEvent
// * EventListTofEventNoTime
// * EventListWeightedEventNoTime

namespace Mantid {
namespace DataObjects {
// Share the same definition of wrapper/interface
class DLLExport EventListBase : public EventList {
private:
  using ISpectrum::copyDataInto;
  void copyDataInto(EventList &sink) const override;
  void copyDataInto(Histogram1D &sink) const override;

  const HistogramData::Histogram &histogramRef() const override { return m_histogram; }
  HistogramData::Histogram &mutableHistogramRef() override;

  /// Histogram object holding the histogram data. Currently only X.
  HistogramData::Histogram m_histogram;

  /// List of Events
  mutable std::vector<Types::Event::Event> events;

  /// What type of event is in our list.
  Mantid::API::EventType eventType;

  /// Last sorting order
  mutable EventSortType order;

  /// MRU lists of the parent EventWorkspace
  mutable EventWorkspaceMRU *mru;

  /// Mutex that is locked while sorting an event list
  mutable std::mutex m_sortMutex;

  template <class T>
  static typename std::vector<T>::const_iterator findFirstPulseEvent(const std::vector<T> &events,
                                                                     const double seek_pulsetime);

  template <class T>
  typename std::vector<T>::const_iterator findFirstTimeAtSampleEvent(const std::vector<T> &events,
                                                                     const double seek_time, const double &tofFactor,
                                                                     const double &tofOffset) const;

  void generateCountsHistogram(const MantidVec &X, MantidVec &Y) const;

  void generateCountsHistogramPulseTime(const MantidVec &X, MantidVec &Y) const;

  void generateCountsHistogramTimeAtSample(const MantidVec &X, MantidVec &Y, const double &tofFactor,
                                           const double &tofOffset) const;

  void generateErrorsHistogram(const MantidVec &Y, MantidVec &E) const;

  void switchToWeightedEvents();
  void switchToWeightedEventsNoTime();
  // should not be called externally
  void sortPulseTimeTOFDelta(const Types::Core::DateAndTime &start, const double seconds) const;

  // helper functions are all internal to simplify the code
  template <class T1, class T2> static void minusHelper(std::vector<T1> &events, const std::vector<T2> &more_events);
  template <class T>
  static void compressEventsHelper(const std::vector<T> &events, std::vector<WeightedEventNoTime> &out,
                                   double tolerance);
  template <class T>
  void compressEventsParallelHelper(const std::vector<T> &events, std::vector<WeightedEventNoTime> &out,
                                    double tolerance);
  template <class T>
  static void compressFatEventsHelper(const std::vector<T> &events, std::vector<WeightedEvent> &out,
                                      const double tolerance, const Mantid::Types::Core::DateAndTime &timeStart,
                                      const double seconds);

  template <class T>
  static void histogramForWeightsHelper(const std::vector<T> &events, const MantidVec &X, MantidVec &Y, MantidVec &E);
  template <class T>
  static void integrateHelper(std::vector<T> &events, const double minX, const double maxX, const bool entireRange,
                              double &sum, double &error);
  template <class T>
  static double integrateHelper(std::vector<T> &events, const double minX, const double maxX, const bool entireRange);
  template <class T> void convertTofHelper(std::vector<T> &events, const std::function<double(double)> &func);

  template <class T> void convertTofHelper(std::vector<T> &events, const double factor, const double offset);
  template <class T> void addPulsetimeHelper(std::vector<T> &events, const double seconds);
  template <class T> void addPulsetimesHelper(std::vector<T> &events, const std::vector<double> &seconds);
  template <class T> static std::size_t maskTofHelper(std::vector<T> &events, const double tofMin, const double tofMax);
  template <class T> static std::size_t maskConditionHelper(std::vector<T> &events, const std::vector<bool> &mask);

  template <class T> static void getTofsHelper(const std::vector<T> &events, std::vector<double> &tofs);
  template <class T> static void getWeightsHelper(const std::vector<T> &events, std::vector<double> &weights);
  template <class T> static void getWeightErrorsHelper(const std::vector<T> &events, std::vector<double> &weightErrors);
  template <class T>
  static void getPulseTimesHelper(const std::vector<T> &events, std::vector<Mantid::Types::Core::DateAndTime> &times);
  template <class T> static void setTofsHelper(std::vector<T> &events, const std::vector<double> &tofs);
  template <class T>
  static void filterByPulseTimeHelper(std::vector<T> &events, Types::Core::DateAndTime start,
                                      Types::Core::DateAndTime stop, std::vector<T> &output);
  template <class T>
  static void filterByTimeAtSampleHelper(std::vector<T> &events, Types::Core::DateAndTime start,
                                         Types::Core::DateAndTime stop, double tofFactor, double tofOffset,
                                         std::vector<T> &output);
  template <class T> void filterInPlaceHelper(Kernel::TimeSplitterType &splitter, typename std::vector<T> &events);
  template <class T>
  void splitByTimeHelper(Kernel::TimeSplitterType &splitter, std::vector<EventList *> outputs,
                         typename std::vector<T> &events) const;
  template <class T>
  void splitByFullTimeHelper(Kernel::TimeSplitterType &splitter, std::map<int, EventList *> outputs,
                             typename std::vector<T> &events, bool docorrection, double toffactor,
                             double tofshift) const;
  /// Split events by pulse time
  template <class T>
  void splitByPulseTimeHelper(Kernel::TimeSplitterType &splitter, std::map<int, EventList *> outputs,
                              typename std::vector<T> &events) const;

  /// Split events (template) by pulse time with matrix splitters
  template <class T>
  void splitByPulseTimeWithMatrixHelper(const std::vector<int64_t> &vec_split_times,
                                        const std::vector<int> &vec_split_target, std::map<int, EventList *> outputs,
                                        typename std::vector<T> &events) const;

  template <class T>
  std::string splitByFullTimeVectorSplitterHelper(const std::vector<int64_t> &vectimes,
                                                  const std::vector<int> &vecgroups, std::map<int, EventList *> outputs,
                                                  typename std::vector<T> &vecEvents, bool docorrection,
                                                  double toffactor, double tofshift) const;

  template <class T>
  std::string
  splitByFullTimeSparseVectorSplitterHelper(const std::vector<int64_t> &vectimes, const std::vector<int> &vecgroups,
                                            std::map<int, EventList *> outputs, typename std::vector<T> &vecEvents,
                                            bool docorrection, double toffactor, double tofshift) const;

  template <class T> static void multiplyHelper(std::vector<T> &events, const double value, const double error = 0.0);
  template <class T>
  static void multiplyHistogramHelper(std::vector<T> &events, const MantidVec &X, const MantidVec &Y,
                                      const MantidVec &E);
  template <class T>
  static void divideHistogramHelper(std::vector<T> &events, const MantidVec &X, const MantidVec &Y, const MantidVec &E);
  template <class T>
  void convertUnitsViaTofHelper(typename std::vector<T> &events, Mantid::Kernel::Unit *fromUnit,
                                Mantid::Kernel::Unit *toUnit);
  template <class T>
  void convertUnitsQuicklyHelper(typename std::vector<T> &events, const double &factor, const double &power);
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