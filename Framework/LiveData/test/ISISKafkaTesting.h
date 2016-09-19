#ifndef MANTID_LIVEDATA_ISISKAFKATESTING_H_
#define MANTID_LIVEDATA_ISISKAFKATESTING_H_

#include "MantidLiveData/Kafka/IKafkaBroker.h"
#include "MantidKernel/WarningSuppressions.h"
#include <gmock/gmock.h>

GCC_DIAG_OFF(conversion)
#include "ISIS/private/Kafka/Schema/det_spec_mapping_schema_generated.h"
#include "ISIS/private/Kafka/Schema/event_schema_generated.h"
#include "ISIS/private/Kafka/Schema/run_info_schema_generated.h"
GCC_DIAG_ON(conversion)

#include <ctime>

namespace ISISKafkaTesting {

// -----------------------------------------------------------------------------
// Mock broker to inject fake subscribers
// -----------------------------------------------------------------------------
class MockKafkaBroker : public Mantid::LiveData::IKafkaBroker {
public:
  using IKafkaStreamSubscriber_uptr =
      std::unique_ptr<Mantid::LiveData::IKafkaStreamSubscriber>;
  using IKafkaStreamSubscriber_ptr = Mantid::LiveData::IKafkaStreamSubscriber *;

  GCC_DIAG_OFF_SUGGEST_OVERRIDE
  // GMock cannot mock non-copyable return types so we resort to a small
  // adapter method. Users have to use EXPECT_CALL(subscribe_) instead
  MOCK_CONST_METHOD1(subscribe_,
                     IKafkaStreamSubscriber_ptr(const std::string &));
  IKafkaStreamSubscriber_uptr subscribe(const std::string &s) const override {
    return std::unique_ptr<Mantid::LiveData::IKafkaStreamSubscriber>(
        this->subscribe_(s));
  }
  GCC_DIAG_ON_SUGGEST_OVERRIDE
};

// -----------------------------------------------------------------------------
// Fake stream to raise error to tests
// -----------------------------------------------------------------------------
class FakeExceptionThrowingStreamSubscriber
    : public Mantid::LiveData::IKafkaStreamSubscriber {
public:
  void subscribe() override {}
  void consumeMessage(std::string *buffer) override {
    buffer->clear();
    throw std::runtime_error("FakeExceptionThrowingStreamSubscriber");
  }
};

// -----------------------------------------------------------------------------
// Fake stream to provide empty stream to client
// -----------------------------------------------------------------------------
class FakeEmptyStreamSubscriber
    : public Mantid::LiveData::IKafkaStreamSubscriber {
public:
  void subscribe() override {}
  void consumeMessage(std::string *buffer) override { buffer->clear(); }
};

// -----------------------------------------------------------------------------
// Fake ISIS event stream to provide event data
// -----------------------------------------------------------------------------
class FakeISISSinglePeriodEventSubscriber
    : public Mantid::LiveData::IKafkaStreamSubscriber {
public:
  void subscribe() override {}
  void consumeMessage(std::string *buffer) override {
    assert(buffer);

    flatbuffers::FlatBufferBuilder builder;
    std::vector<int32_t> spec = {5, 4, 3, 2, 1, 2};
    std::vector<float> tof = {11000, 10000, 9000, 8000, 7000, 6000};
    auto messageNEvents = ISISStream::CreateNEvents(
        builder, builder.CreateVector(tof), builder.CreateVector(spec));

    int32_t frameNumber(2), period(1);
    float frameTime(1.f), protonCharge(0.5f);
    bool endOfFrame(false), endOfRun(false);
    // No SE events
    auto messageFramePart = ISISStream::CreateFramePart(
        builder, frameNumber, frameTime, ISISStream::RunState_RUNNING,
        protonCharge, period, endOfFrame, endOfRun, messageNEvents);
    auto messageFlatbuf = ISISStream::CreateEventMessage(
        builder, ISISStream::MessageTypes_FramePart, messageFramePart.Union());
    builder.Finish(messageFlatbuf);

    // Copy to provided buffer
    buffer->assign(reinterpret_cast<const char *>(builder.GetBufferPointer()),
                   builder.GetSize());
  }
};

// -----------------------------------------------------------------------------
// Fake ISIS event stream to provide event data
// -----------------------------------------------------------------------------
class FakeISISMultiplePeriodEventSubscriber
    : public Mantid::LiveData::IKafkaStreamSubscriber {
public:
  void subscribe() override {}
  void consumeMessage(std::string *buffer) override { assert(buffer); }
};

// -----------------------------------------------------------------------------
// Fake ISIS run data stream
// -----------------------------------------------------------------------------
class FakeISISRunInfoStreamSubscriber
    : public Mantid::LiveData::IKafkaStreamSubscriber {
public:
  FakeISISRunInfoStreamSubscriber(int32_t nperiods) : m_nperiods(nperiods) {}
  void subscribe() override {}
  void consumeMessage(std::string *buffer) override {
    assert(buffer);

    // Convert date to time_t
    std::tm tmb;
    strptime(m_startTime.c_str(), "%Y-%m-%dT%H:%M:%S", &tmb);
    uint64_t startTime = static_cast<uint64_t>(std::mktime(&tmb));

    // Serialize data with flatbuffers
    flatbuffers::FlatBufferBuilder builder;
    auto runInfo = ISISStream::CreateRunInfo(builder, startTime, m_runNumber,
                                             builder.CreateString(m_instName),
                                             m_streamOffset, m_nperiods);
    builder.Finish(runInfo);
    // Copy to provided buffer
    buffer->assign(reinterpret_cast<const char *>(builder.GetBufferPointer()),
                   builder.GetSize());
  }

private:
  std::string m_startTime = "2016-08-31T12:07:42";
  int32_t m_runNumber = 1000;
  std::string m_instName = "HRPDTEST";
  int64_t m_streamOffset = 0;
  int32_t m_nperiods = 1;
};

// -----------------------------------------------------------------------------
// Fake ISIS spectra-detector stream
// -----------------------------------------------------------------------------
class FakeISISSpDetStreamSubscriber
    : public Mantid::LiveData::IKafkaStreamSubscriber {
public:
  void subscribe() override {}
  void consumeMessage(std::string *buffer) override {
    assert(buffer);

    // Serialize data with flatbuffers
    flatbuffers::FlatBufferBuilder builder;
    auto specVector = builder.CreateVector(m_spec);
    auto detIdsVector = builder.CreateVector(m_detid);
    auto spdet = ISISStream::CreateSpectraDetectorMapping(
        builder, specVector, detIdsVector, static_cast<int32_t>(m_spec.size()));
    builder.Finish(spdet);
    // Copy to provided buffer
    buffer->assign(reinterpret_cast<const char *>(builder.GetBufferPointer()),
                   builder.GetSize());
  }

private:
  std::vector<int32_t> m_spec = {1, 2, 3, 4, 5};
  // These match the detector numbers in HRPDTEST_Definition.xml
  std::vector<int32_t> m_detid = {1001, 1002, 1100, 901000, 10100};
};
}

#endif // MANTID_LIVEDATA_ISISKAFKAEVENTSTREAMDECODERTESTMOCKS_H_
