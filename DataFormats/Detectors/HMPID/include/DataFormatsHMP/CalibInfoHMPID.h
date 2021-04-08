// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file CalibInfoHMPID.h
/// \brief Class to store the output of the matching to HMPID for calibration

#ifndef ALICEO2_CALIBINFOHMPID_H
#define ALICEO2_CALIBINFOHMPID_H

#include "Rtypes.h"

namespace o2
{
namespace dataformats
{
class CalibInfoHMPID
{
 public:
  CalibInfoHMPID(int indexHMPIDCh, int timestamp, float DeltaTimePi, float tot, int flags = 0) : mHMPIDChIndex(indexHMPIDCh), mTimestamp(timestamp), mDeltaTimePi(DeltaTimePi), mTot(tot), mFlags(flags){};
  CalibInfoHMPID() = default;
  ~CalibInfoHMPID() = default;

  void setHMPIDChIndex(int index) { mHMPIDChIndex = index; }
  int getHMPIDChIndex() const { return mHMPIDChIndex; }

  void setTimestamp(int ts) { mTimestamp = ts; }
  int getTimestamp() const { return mTimestamp; }

  void setDeltaTimePi(float time) { mDeltaTimePi = time; }
  float getDeltaTimePi() const { return mDeltaTimePi; }

  void setTot(float tot) { mTot = tot; }
  float getTot() const { return mTot; }

  void setFlags(int flags) { mFlags = flags; }
  float getFlags() const { return mFlags; }

 private:
  int mHMPIDChIndex;      // index of the HMPID channel
  int mTimestamp;       // timestamp in seconds
  float mDeltaTimePi;   // raw HMPID time - expected time for pi hypotesis
  float mTot;           // time-over-threshold
  unsigned char mFlags; // bit mask with quality flags (to be defined)

  ClassDefNV(CalibInfoHMPID, 1);
};
} // namespace dataformats
} // namespace o2
#endif
