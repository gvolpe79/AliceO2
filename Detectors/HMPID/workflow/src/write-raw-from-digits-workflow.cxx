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
/// \file    raw-to-digits-workflow.cxx
/// \author  Andrea Ferrero
///
/// \brief This is an executable that runs the decoding via DPL.
///
/// This is an executable that takes a raw buffer from the Data Processing Layer, runs the decoding and sends the digits via the Data Processing Layer.
/// The decoder expects an input buffer in the format generated by o2-raw-file-reader-workflow
///

#include "Framework/WorkflowSpec.h"
#include "Framework/DataSpecUtils.h"
#include "Framework/CallbackService.h"
#include "Framework/ControlService.h"
#include "Framework/Task.h"
#include "Framework/runDataProcessing.h"

#include "HMPIDWorkflow/WriteRawFromDigitsSpec.h"

using namespace o2;
using namespace o2::framework;

WorkflowSpec defineDataProcessing(const ConfigContext&)
{
  WorkflowSpec specs;

  DataProcessorSpec consumer = o2::hmpid::getWriteRawFromDigitsSpec();
//  DataProcessorSpec consumer = o2::hmpid::getDecodingSpec();
  specs.push_back(consumer);
//  specs.push_back(consumer);

  return specs;
}
