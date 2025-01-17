// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include <cxxtest/TestSuite.h>

#include "MantidAPI/AlgorithmManager.h"
#include "MantidKernel/ConfigService.h"
#include "MantidKernel/FacilityInfo.h"
#include "MantidRemoteAlgorithms/DownloadRemoteFile.h"

using namespace Mantid::RemoteAlgorithms;

class DownloadRemoteFileTest : public CxxTest::TestSuite {
public:
  // This pair of boilerplate methods prevent the suite being created statically
  // This means the constructor isn't called when running other tests
  static DownloadRemoteFileTest *createSuite() { return new DownloadRemoteFileTest(); }
  static void destroySuite(DownloadRemoteFileTest *suite) { delete suite; }

  void test_algorithm() {
    testAlg = Mantid::API::AlgorithmManager::Instance().create("DownloadRemoteFile", 1);
    TS_ASSERT(testAlg);
    TS_ASSERT_EQUALS(testAlg->name(), "DownloadRemoteFile");
    TS_ASSERT_EQUALS(testAlg->version(), 1);
  }

  void test_castAlgorithm() {
    // can create
    std::shared_ptr<DownloadRemoteFile> a;
    TS_ASSERT(a = std::make_shared<DownloadRemoteFile>());

    // can cast to inherited interfaces and base classes
    TS_ASSERT(dynamic_cast<Mantid::RemoteAlgorithms::DownloadRemoteFile *>(a.get()));
    TS_ASSERT(dynamic_cast<Mantid::API::Algorithm *>(a.get()));
    TS_ASSERT(dynamic_cast<Mantid::Kernel::PropertyManagerOwner *>(a.get()));
    TS_ASSERT(dynamic_cast<Mantid::API::IAlgorithm *>(a.get()));
    TS_ASSERT(dynamic_cast<Mantid::Kernel::IPropertyManager *>(a.get()));
  }

  void test_init() {
    if (!testAlg->isInitialized())
      TS_ASSERT_THROWS_NOTHING(testAlg->initialize());

    TS_ASSERT(testAlg->isInitialized());

    DownloadRemoteFile dl;
    TS_ASSERT_THROWS_NOTHING(dl.initialize());
  }

  // TODO: when we have a RemoteJobManager capable of creating
  // algorithms for different types of compute resources (example:
  // Fermi@SNS), create different algorithms for them
  void test_propertiesMissing() {
    DownloadRemoteFile alg1;
    TS_ASSERT_THROWS_NOTHING(alg1.initialize());
    // Transaction id missing
    TS_ASSERT_THROWS(alg1.setPropertyValue("ComputeResource", "missing!"), const std::invalid_argument &);
    TS_ASSERT_THROWS_NOTHING(alg1.setPropertyValue("RemoteFileName", "file name"));

    TS_ASSERT_THROWS(alg1.execute(), const std::runtime_error &);
    TS_ASSERT(!alg1.isExecuted());

    DownloadRemoteFile alg2;
    TS_ASSERT_THROWS_NOTHING(alg2.initialize());
    // file name missing
    TS_ASSERT_THROWS(alg2.setPropertyValue("ComputeResource", "missing!"), const std::invalid_argument &);
    TS_ASSERT_THROWS_NOTHING(alg2.setPropertyValue("TransactionID", "id001"));

    TS_ASSERT_THROWS(alg2.execute(), const std::runtime_error &);
    TS_ASSERT(!alg2.isExecuted());

    DownloadRemoteFile alg3;
    TS_ASSERT_THROWS_NOTHING(alg3.initialize());
    // compute resource missing
    TS_ASSERT_THROWS_NOTHING(alg3.setPropertyValue("RemoteFileName", "file name"));
    TS_ASSERT_THROWS_NOTHING(alg3.setPropertyValue("TransactionID", "id001"));

    TS_ASSERT_THROWS(alg3.execute(), const std::runtime_error &);
    TS_ASSERT(!alg3.isExecuted());
  }

  void test_wrongProperty() {
    DownloadRemoteFile dl;
    TS_ASSERT_THROWS_NOTHING(dl.initialize();)
    TS_ASSERT_THROWS(dl.setPropertyValue("Compute", "anything"), const std::runtime_error &);
    TS_ASSERT_THROWS(dl.setPropertyValue("TransID", "anything"), const std::runtime_error &);
    TS_ASSERT_THROWS(dl.setPropertyValue("FileName", "anything"), const std::runtime_error &);
  }

  void test_propertiesOK() {
    testFacilities.emplace_back("SNS", "Fermi");

    const Mantid::Kernel::FacilityInfo &prevFac = Mantid::Kernel::ConfigService::Instance().getFacility();
    for (auto &testFacility : testFacilities) {
      const std::string facName = testFacility.first;
      const std::string compName = testFacility.second;

      Mantid::Kernel::ConfigService::Instance().setFacility(facName);
      DownloadRemoteFile dl;
      TS_ASSERT_THROWS_NOTHING(dl.initialize());
      TS_ASSERT_THROWS_NOTHING(dl.setPropertyValue("ComputeResource", compName));
      TS_ASSERT_THROWS_NOTHING(dl.setPropertyValue("TransactionID", "anything"));
      TS_ASSERT_THROWS_NOTHING(dl.setPropertyValue("RemoteFileName", "anything"));
      // TODO: this would run the algorithm and do a remote
      // connection. uncomment only when/if we have a mock up for this
      // TS_ASSERT_THROWS(dl.execute(), std::exception);
      TS_ASSERT(!dl.isExecuted());
    }
    Mantid::Kernel::ConfigService::Instance().setFacility(prevFac.name());
  }

  // TODO: void test_runOK() - with a mock when we can add it.
  // ideally, with different compute resources to check the remote job
  // manager factory, etc.

private:
  Mantid::API::IAlgorithm_sptr testAlg;
  std::vector<std::pair<std::string, std::string>> testFacilities;
};
