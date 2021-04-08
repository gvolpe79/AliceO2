
// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef HMPID_CALIB_COLLECTOR_H_
#define HMPID_CALIB_COLLECTOR_H_

#include "DetectorsCalibration/TimeSlotCalibration.h"
#include "DetectorsCalibration/TimeSlot.h"
#include "DataFormatsHMP/CalibInfoHMPID.h"
#include "HMPIDBase/Geo.h"
#include "DataFormatsHMP/CalibInfoHMPIDshort.h"

#include <array>

namespace o2
{
namespace hmpid
{

class HMPIDCalibInfoSlot
{

  using Slot = o2::calibration::TimeSlot<o2::hmpid::HMPIDCalibInfoSlot>;
  using Geo = o2::hmpid::Geo;

 public:
  static constexpr int NCHANNELSXSECTOR = o2::hmpid::Geo::NCHANNELS / o2::hmpid::Geo::NSECTORS;

  HMPIDCalibInfoSlot()
  {
    for (int ch = 0; ch < Geo::NCHANNELS; ch++) {
      mEntriesSlot[ch] = 0;
    }
  }

  ~HMPIDCalibInfoSlot() = default;

  void print() const;
  void printEntries() const;
  void fill(const gsl::span<const o2::dataformats::CalibInfoHMPID> data);
  void merge(const HMPIDCalibInfoSlot* prev);

  auto& getEntriesPerChannel() const { return mEntriesSlot; }
  auto& getEntriesPerChannel() { return mEntriesSlot; }
  auto& getCollectedCalibInfoSlot() { return mHMPIDCollectedCalibInfoSlot; }
  auto& getCollectedCalibInfoSlot() const { return mHMPIDCollectedCalibInfoSlot; }

 private:
  std::array<int, Geo::NCHANNELS> mEntriesSlot;                               // vector containing number of entries per channel
  std::vector<o2::dataformats::CalibInfoHMPIDshort> mHMPIDCollectedCalibInfoSlot; ///< output HMPID calibration info

  ClassDefNV(HMPIDCalibInfoSlot, 1);
};

class HMPIDCalibCollector final : public o2::calibration::TimeSlotCalibration<o2::dataformats::CalibInfoHMPID, o2::hmpid::HMPIDCalibInfoSlot>
{
  using TFType = uint64_t;
  using Slot = o2::calibration::TimeSlot<o2::hmpid::HMPIDCalibInfoSlot>;

 public:
  HMPIDCalibCollector(bool TFsendingPolicy, int maxNumOfHits, bool test = false) : mTFsendingPolicy(TFsendingPolicy), mMaxNumOfHits(maxNumOfHits), mTest(test){};

  ~HMPIDCalibCollector() final = default;

  bool hasEnoughData(const Slot& slot) const final;
  void initOutput() final;
  void finalizeSlot(Slot& slot) final;
  Slot& emplaceNewSlot(bool front, TFType tstart, TFType tend) final;
  void setIsTest(bool istest) { mTest = istest; }
  auto& getCollectedCalibInfo() const { return mHMPIDCollectedCalibInfo; }
  auto& getEntriesPerChannel() const { return mEntries; }
  void setIsMaxNumberOfHitsAbsolute(bool absNumber) { mAbsMaxNumOfHits = absNumber; }

 private:
  bool mTFsendingPolicy = false;                                          // whether we will send information at every TF or only when we have a certain statistics
  int mMaxNumOfHits = 500;                                                // maximum number of hits for one single channel to trigger the sending of the information (if mTFsendingPolicy = false)
  bool mTest = false;                                                     // flag to say whether we are in test mode or not
  bool mAbsMaxNumOfHits = true;                                           // to decide if the mMaxNumOfHits should be multiplied by the number of HMPID channels
  std::array<int, Geo::NCHANNELS> mEntries;                               // vector containing number of entries per channel
  std::vector<o2::dataformats::CalibInfoHMPIDshort> mHMPIDCollectedCalibInfo; ///< output HMPID calibration info

  ClassDefOverride(HMPIDCalibCollector, 1);
};

} // end namespace hmpid
} // end namespace o2

#endif /* HMPID_CHANNEL_CALIBRATOR_H_ */

