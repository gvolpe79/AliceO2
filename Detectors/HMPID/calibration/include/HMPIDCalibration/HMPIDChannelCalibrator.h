// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef HMPID_CHANNEL_CALIBRATOR_H_
#define HMPID_CHANNEL_CALIBRATOR_H_

#include "DetectorsCalibration/TimeSlotCalibration.h"
#include "DetectorsCalibration/TimeSlot.h"
#include "DataFormatsHMPID/CalibInfoHMPID.h"
#include "DataFormatsHMPID/CalibLHCphaseHMPID.h"
#include "HMPIDBase/Geo.h"
#include "CCDB/CcdbObjectInfo.h"
#include "HMPIDCalibration/CalibHMPIDapi.h"

#include <array>
#include <boost/histogram.hpp>

namespace o2
{
namespace hmpid
{

class HMPIDChannelData
{

  using Slot = o2::calibration::TimeSlot<o2::hmpid::HMPIDChannelData>;
  using CalibHMPIDapi = o2::hmpid::CalibHMPIDapi;
  using boostHisto = boost::histogram::histogram<std::tuple<boost::histogram::axis::regular<double, boost::use_default, boost::use_default, boost::use_default>, boost::histogram::axis::integer<>>, boost::histogram::unlimited_storage<std::allocator<char>>>;

 public:
  static constexpr int NCHANNELSXSECTOR = o2::hmpid::Geo::NCHANNELS / o2::hmpid::Geo::NSECTORS;

  HMPIDChannelData()
  {
    LOG(INFO) << "Default c-tor, not to be used";
  }

  HMPIDChannelData(int nb, float r, CalibHMPIDapi* cta) : mNBins(nb), mRange(r), mCalibHMPIDapi(cta)
  {
    if (r <= 0. || nb < 1) {
      throw std::runtime_error("Wrong initialization of the histogram");
    }
    mV2Bin = mNBins / (2 * mRange);
    for (int isect = 0; isect < 18; isect++) {
      mHisto[isect] = boost::histogram::make_histogram(boost::histogram::axis::regular<>(mNBins, -mRange, mRange, "t-texp"),
                                                       boost::histogram::axis::integer<>(0, o2::hmpid::Geo::NPADSXSECTOR, "channel index in sector" + std::to_string(isect))); // bin is defined as [low, high[
    }
    mEntries.resize(o2::hmpid::Geo::NCHANNELS, 0);
  }

  ~HMPIDChannelData() = default;

  void print() const;
  void print(int isect) const;
  void printEntries() const;
  void fill(const gsl::span<const o2::dataformats::CalibInfoHMPID> data);
  void merge(const HMPIDChannelData* prev);
  int findBin(float v) const;
  float integral(int chmin, int chmax, float binmin, float binmax) const;
  float integral(int chmin, int chmax, int binxmin, int binxmax) const;
  float integral(int ch, float binmin, float binmax) const;
  float integral(int ch, int binxmin, int binxmax) const;
  float integral(int ch) const;
  bool hasEnoughData(int minEntries) const;

  float getRange() const { return mRange; }
  void setRange(float r) { mRange = r; }

  int getNbins() const { return mNBins; }
  void setNbins(int nb) { mNBins = nb; }

  boostHisto& getHisto(int isect) { return mHisto[isect]; }
  const boostHisto& getHisto(int isect) const { return mHisto[isect]; }
  //const boostHisto getHisto() const { return &mHisto[0]; }
  // boostHisto* getHisto(int isect) const { return &mHisto[isect]; }

  std::vector<int> getEntriesPerChannel() const { return mEntries; }

 private:
  float mRange = o2::hmpid::Geo::BC_TIME_INPS * 0.5;
  int mNBins = 1000;
  float mV2Bin;
  std::array<boostHisto, 18> mHisto;
  std::vector<int> mEntries; // vector containing number of entries per channel

  CalibHMPIDapi* mCalibHMPIDapi = nullptr; // calibHMPIDapi to correct the t-text

  ClassDefNV(HMPIDChannelData, 1);
};

class HMPIDChannelCalibrator final : public o2::calibration::TimeSlotCalibration<o2::dataformats::CalibInfoHMPID, o2::hmpid::HMPIDChannelData>
{
  using TFType = uint64_t;
  using Slot = o2::calibration::TimeSlot<o2::hmpid::HMPIDChannelData>;
  using CalibHMPIDapi = o2::hmpid::CalibHMPIDapi;
  using TimeSlewing = o2::dataformats::CalibTimeSlewingParamHMPID;
  using CcdbObjectInfo = o2::ccdb::CcdbObjectInfo;
  using CcdbObjectInfoVector = std::vector<CcdbObjectInfo>;
  using TimeSlewingVector = std::vector<TimeSlewing>;

 public:
  static const int NCHANNELSXSECTOR = o2::hmpid::Geo::NCHANNELS / o2::hmpid::Geo::NSECTORS;
  HMPIDChannelCalibrator(int minEnt = 500, int nb = 1000, float r = 24400) : mMinEntries(minEnt), mNBins(nb), mRange(r){};

  ~HMPIDChannelCalibrator() final = default;

  bool hasEnoughData(const Slot& slot) const final;
  void initOutput() final;
  void finalizeSlot(Slot& slot) final;
  Slot& emplaceNewSlot(bool front, TFType tstart, TFType tend) final;

  const TimeSlewingVector& getTimeSlewingVector() const { return mTimeSlewingVector; }
  const CcdbObjectInfoVector& getTimeSlewingInfoVector() const { return mInfoVector; }
  CcdbObjectInfoVector& getTimeSlewingInfoVector() { return mInfoVector; }

  void setIsTest(bool isTest) { mTest = isTest; }
  bool isTest() const { return mTest; }

  void setCalibHMPIDapi(CalibHMPIDapi* api) { mCalibHMPIDapi = api; }
  CalibHMPIDapi* getCalibHMPIDapi() const { return mCalibHMPIDapi; }

  void setRange(float r) { mRange = r; }
  float getRange() const { return mRange; }

 private:
  int mMinEntries = 0; // min number of entries to calibrate the TimeSlot
  int mNBins = 0;      // bins of the histogram with the t-text per channel
  float mRange = 0.;   // range of the histogram with the t-text per channel
  bool mTest = false;  // flag to be used when running in test mode: it simplify the processing (e.g. does not go through all channels)

  CalibHMPIDapi* mCalibHMPIDapi = nullptr; // CalibHMPIDapi needed to get the previous calibrations read from CCDB (do we need that it is a pointer?)

  // output
  CcdbObjectInfoVector mInfoVector;     // vector of CCDB Infos , each element is filled with the CCDB description of the accompanying TimeSlewing object
  TimeSlewingVector mTimeSlewingVector; // vector of TimeSlewing, each element is filled in "process"
                                        // when we finalize one slot (multiple can be finalized
                                        // during the same "process", which is why we have a vector).
                                        // Each element is to be considered the output of the device,
                                        // and will go to the CCDB. Note that for the channel offset
                                        // we still fill the TimeSlewing object

  ClassDefOverride(HMPIDChannelCalibrator, 1);
};

} // end namespace hmpid
} // end namespace o2

#endif /* HMPID_CHANNEL_CALIBRATOR_H_ */
