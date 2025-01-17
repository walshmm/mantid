// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include <cxxtest/TestSuite.h>

#include "MantidAPI/AnalysisDataService.h"
#include "MantidAPI/Axis.h"
#include "MantidAPI/Run.h"
#include "MantidAlgorithms/CorelliCrossCorrelate.h"
#include "MantidDataObjects/EventWorkspace.h"
#include "MantidKernel/DateAndTime.h"
#include "MantidKernel/TimeSeriesProperty.h"

using Mantid::Algorithms::CorelliCrossCorrelate;
using namespace Mantid::API;
using namespace Mantid::DataObjects;
using namespace Mantid::Kernel;
using Mantid::Types::Core::DateAndTime;
using Mantid::Types::Event::TofEvent;

class CorelliCrossCorrelateTest : public CxxTest::TestSuite {
public:
  // This pair of boilerplate methods prevent the suite being created statically
  // This means the constructor isn't called when running other tests
  static CorelliCrossCorrelateTest *createSuite() { return new CorelliCrossCorrelateTest(); }
  static void destroySuite(CorelliCrossCorrelateTest *suite) { delete suite; }

  void test_Init() {
    CorelliCrossCorrelate alg;
    TS_ASSERT_THROWS_NOTHING(alg.initialize())
    TS_ASSERT(alg.isInitialized())
  }

  void test_exec() {
    // Name of the output workspace.
    std::string outWSName("CorelliCrossCorrelateTest_OutputWS");

    IAlgorithm_sptr lei = AlgorithmFactory::Instance().create("LoadEmptyInstrument", 1);
    lei->initialize();
    lei->setPropertyValue("Filename", "CORELLI_Definition.xml");
    lei->setPropertyValue("OutputWorkspace", "CorelliCrossCorrelateTest_OutputWS");
    lei->setPropertyValue("MakeEventWorkspace", "1");
    lei->execute();

    EventWorkspace_sptr ws;
    ws = AnalysisDataService::Instance().retrieveWS<EventWorkspace>("CorelliCrossCorrelateTest_OutputWS");

    DateAndTime startTime("2007-11-30T16:17:00");
    auto &evlist = ws->getSpectrum(0);

    // Add some events to the workspace.
    evlist.addEventQuickly(TofEvent(10.0, startTime + 0.007));
    evlist.addEventQuickly(TofEvent(100.0, startTime + 0.012));
    evlist.addEventQuickly(TofEvent(1000.0, startTime + 0.012));
    evlist.addEventQuickly(TofEvent(10000.0, startTime + 0.012));
    evlist.addEventQuickly(TofEvent(1222.0, startTime + 0.03));

    ws->getAxis(0)->setUnit("TOF");

    ws->sortAll(PULSETIME_SORT, nullptr);

    // Add some chopper TDCs to the workspace.
    double period = 1 / 293.383;
    auto tdc = new TimeSeriesProperty<int>("chopper4_TDC");
    for (int i = 0; i < 10; i++) {
      double tdcTime = i * period;
      tdc->addValue(startTime + tdcTime, 1);
    }
    ws->mutableRun().addLogData(tdc);

    // Add motorSpeed to the workspace
    auto motorSpeed = new TimeSeriesProperty<double>("BL9:Chop:Skf4:MotorSpeed");
    motorSpeed->addValue(startTime, 293.383);
    ws->mutableRun().addLogData(motorSpeed);

    CorelliCrossCorrelate alg;
    TS_ASSERT_THROWS_NOTHING(alg.initialize())
    TS_ASSERT(alg.isInitialized())
    TS_ASSERT_THROWS_NOTHING(alg.setPropertyValue("InputWorkspace", "CorelliCrossCorrelateTest_OutputWS"));
    TS_ASSERT_THROWS_NOTHING(alg.setPropertyValue("OutputWorkspace", "CorelliCrossCorrelateTest_OutputWS"));
    TS_ASSERT_THROWS_NOTHING(alg.setPropertyValue("TimingOffset", "20000"));
    TS_ASSERT_THROWS_NOTHING(alg.execute(););
    TS_ASSERT(alg.isExecuted());

    // Retrieve the workspace from data service.
    TS_ASSERT_THROWS_NOTHING(
        ws = AnalysisDataService::Instance().retrieveWS<EventWorkspace>("CorelliCrossCorrelateTest_OutputWS"));
    TS_ASSERT(ws);
    if (!ws)
      return;

    std::vector<WeightedEvent> &events = evlist.getWeightedEvents();

    TS_ASSERT_DELTA(events[0].weight(), -0.993919, 0.00001)
    TS_ASSERT_DELTA(events[1].weight(), -0.993919, 0.00001)
    TS_ASSERT_DELTA(events[2].weight(), 1.0, 0.00001)
    TS_ASSERT_DELTA(events[3].weight(), -0.993919, 0.00001)
    TS_ASSERT_DELTA(events[4].weight(), 1.0, 0.00001)

    // Remove workspace from the data service.
    AnalysisDataService::Instance().remove(outWSName);
  }
};
