// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

///
/// \file    DatDecoderSpec.h
/// \author  Andrea Ferrero
///
/// \brief Definition of a data processor to run the raw decoding
///

#ifndef DETECTORS_HMPID_WORKFLOW_INCLUDE_HMPIDWORKFLOW_DATADECODERSPEC_H_
#define DETECTORS_HMPID_WORKFLOW_INCLUDE_HMPIDWORKFLOW_DATADECODERSPEC_H_

#include "Framework/DataProcessorSpec.h"
#include "Framework/Task.h"

#include "HMPIDReconstruction/HmpidDecodeRawMem.h"

namespace o2
{
namespace hmpid
{

  class DataDecoderTask : public framework::Task
  {
    public:
      DataDecoderTask() = default;
      ~DataDecoderTask() override = default;
      void init(framework::InitContext& ic) final;
      void run(framework::ProcessingContext& pc) final;
      void decodeTF(framework::ProcessingContext& pc);
      void decodeReadout(framework::ProcessingContext& pc);
      void decodeRawFile(framework::ProcessingContext& pc);

    private:
      HmpidDecodeRawDigit *mDeco;

  };

o2::framework::DataProcessorSpec getDecodingSpec(std::string inputSpec = "TF:HMP/RAWDATA");
//o2::framework::DataProcessorSpec getDecodingSpec();
} // end namespace hmpid
} // end namespace o2

#endif