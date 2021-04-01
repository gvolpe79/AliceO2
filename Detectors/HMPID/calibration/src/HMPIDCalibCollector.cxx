
// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "HMPIDCalibration/HMPIDCalibCollector.h"
#include "Framework/Logger.h"
#include <cassert>
#include <iostream>
#include <sstream>
#include <TStopwatch.h>

namespace o2
{
namespace hmpid
{

using Slot = o2::calibration::TimeSlot<o2::hmpid::HMPIDCalibInfoSlot>;

//_____________________________________________
void HMPIDCalibInfoSlot::fill(const gsl::span<const o2::dataformats::CalibInfoHMPID> data)
{
  // fill container
  // we do not apply any calibration at this stage, it will be applied when we
  // process the data before filling the CCDB in the separate process

  // we first order the data that arrived, to improve speed when filling
  int nd = data.size();
  LOG(DEBUG) << "entries in incoming data = " << nd;
  std::vector<int> ord(nd);
  std::iota(ord.begin(), ord.end(), 0);
  std::sort(ord.begin(), ord.end(), [&data](int i, int j) { return data[i].getHMPIDChIndex() < data[j].getHMPIDChIndex(); });
  int chPrev = 0, offsPrev = 0;
  for (int i = 0; i < nd; i++) {
    const auto& dti = data[ord[i]];
    auto ch = dti.getHMPIDChIndex();
    auto offset = offsPrev;
    if (ch > chPrev) {
      offset += std::accumulate(mEntriesSlot.begin() + chPrev, mEntriesSlot.begin() + ch, 0);
    }
    offsPrev = offset;
    chPrev = ch;
    mHMPIDCollectedCalibInfoSlot.emplace(mHMPIDCollectedCalibInfoSlot.begin() + offset, data[ord[i]].getTimestamp(), data[ord[i]].getDeltaTimePi(), data[ord[i]].getTot(), data[ord[i]].getFlags());
    mEntriesSlot[ch]++;
  }
}
//_____________________________________________
void HMPIDCalibInfoSlot::merge(const HMPIDCalibInfoSlot* prev)
{
  // merge data of 2 slots

  LOG(DEBUG) << "Merging two slots with entries: current slot -> " << mHMPIDCollectedCalibInfoSlot.size() << " , previous slot -> " << prev->mHMPIDCollectedCalibInfoSlot.size();

  int offset = 0, offsetPrev = 0;
  std::vector<o2::dataformats::CalibInfoHMPIDshort> tmpVector;
  for (int ch = 0; ch < Geo::NCHANNELS; ch++) {
    if (mEntriesSlot[ch] != 0) {
      for (int i = offset; i < offset + mEntriesSlot[ch]; i++) {
        tmpVector.emplace_back(mHMPIDCollectedCalibInfoSlot[i]);
      }
      offset += mEntriesSlot[ch];
    }
    if (prev->mEntriesSlot[ch] != 0) {
      for (int i = offsetPrev; i < offsetPrev + prev->mEntriesSlot[ch]; i++) {
        tmpVector.emplace_back(prev->mHMPIDCollectedCalibInfoSlot[i]);
      }
      offsetPrev += prev->mEntriesSlot[ch];
      mEntriesSlot[ch] += prev->mEntriesSlot[ch];
    }
  }
  mHMPIDCollectedCalibInfoSlot.swap(tmpVector);
  LOG(DEBUG) << "After merging the size is " << mHMPIDCollectedCalibInfoSlot.size();
  return;
}
//_____________________________________________
void HMPIDCalibInfoSlot::print() const
{
  // to print number of entries in the tree and the channel with the max number of entries

  LOG(INFO) << "Total number of entries " << mHMPIDCollectedCalibInfoSlot.size();
  auto maxElementIndex = std::max_element(mEntriesSlot.begin(), mEntriesSlot.end());
  auto channelIndex = std::distance(mEntriesSlot.begin(), maxElementIndex);
  LOG(INFO) << "The maximum number of entries per channel in the current mHMPIDCollectedCalibInfo is " << *maxElementIndex << " for channel " << channelIndex;
  return;
}

//_____________________________________________
void HMPIDCalibInfoSlot::printEntries() const
{
  // to print number of entries in the tree and per channel

  LOG(INFO) << "Total number of entries " << mHMPIDCollectedCalibInfoSlot.size();
  for (int i = 0; i < mEntriesSlot.size(); ++i) {
    if (mEntriesSlot[i] != 0) {
      LOG(INFO) << "channel " << i << " has " << mEntriesSlot[i] << " entries";
    }
  }
  return;
}

//===================================================================

//_____________________________________________
void HMPIDCalibCollector::initOutput()
{
  // emptying the vectors

  mHMPIDCollectedCalibInfo.clear();
  for (int ch = 0; ch < Geo::NCHANNELS; ch++) {
    mEntries[ch] = 0;
  }

  return;
}

//_____________________________________________
bool HMPIDCalibCollector::hasEnoughData(const Slot& slot) const
{

  // We define that we have enough data if the tree is big enough.
  // each CalibInfoHMPIDShort is composed of one int, two floats, one unsigned char --> 13 bytes
  // E.g. supposing that we have 256 entries per channel (which is an upper limit ) --> ~523 MB
  // we can check if we have at least 1 GB of data --> 500*o2::hmpid::Geo::NCHANNELS entries in the vector
  // (see header file for the fact that mMaxNumOfHits = 500)
  // The case in which mScaleMaxNumOfHits = false allows for a fast check

  if (mTest) {
    return true;
  }
  const o2::hmpid::HMPIDCalibInfoSlot* c = slot.getContainer();
  LOG(INFO) << "we have " << c->getCollectedCalibInfoSlot().size() << " entries";
  int maxNumberOfHits = mAbsMaxNumOfHits ? mMaxNumOfHits : mMaxNumOfHits * o2::hmpid::Geo::NCHANNELS;
  if (mTFsendingPolicy || c->getCollectedCalibInfoSlot().size() > maxNumberOfHits) {
    return true;
  }
  return false;
}

//_____________________________________________
void HMPIDCalibCollector::finalizeSlot(Slot& slot)
{
  // here we fill the tree with the remaining stuff that was not filled before

  o2::hmpid::HMPIDCalibInfoSlot* c = slot.getContainer();
  mHMPIDCollectedCalibInfo = c->getCollectedCalibInfoSlot();
  LOG(DEBUG) << "vector of CalibHMPIDInfoShort received with size = " << mHMPIDCollectedCalibInfo.size();
  mEntries = c->getEntriesPerChannel();
  return;
}

//_____________________________________________
Slot& HMPIDCalibCollector::emplaceNewSlot(bool front, TFType tstart, TFType tend)
{

  auto& cont = getSlots();
  auto& slot = front ? cont.emplace_front(tstart, tend) : cont.emplace_back(tstart, tend);
  slot.setContainer(std::make_unique<HMPIDCalibInfoSlot>());
  return slot;
}

} // end namespace hmpid
} // end namespace o2

