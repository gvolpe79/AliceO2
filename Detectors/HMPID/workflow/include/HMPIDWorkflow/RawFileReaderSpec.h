/*
 * RawFilereader.h
 *
 *  Created on: 18 gen 2021
 *      Author: fap
 */

#ifndef DETECTORS_HMPID_WORKFLOW_INCLUDE_HMPIDWORKFLOW_RAWFILEREADERSPEC_H_
#define DETECTORS_HMPID_WORKFLOW_INCLUDE_HMPIDWORKFLOW_RAWFILEREADERSPEC_H_

#include <random>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include "Framework/CallbackService.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/ControlService.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/Lifetime.h"
#include "Framework/Output.h"
#include "Framework/Task.h"
#include "Framework/WorkflowSpec.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/runDataProcessing.h"

namespace o2
{
namespace hmpid
{

class RawFileReaderTask : public framework::Task
{
  public:
    RawFileReaderTask() = default;
    ~RawFileReaderTask() override = default;
    void init(framework::InitContext& ic) final;
    void run(framework::ProcessingContext& pc) final;

  private:
    std::ifstream mInputFile{}; ///< input file
};

o2::framework::DataProcessorSpec getRawFileReaderSpec(std::string inputSpec = "HMP/rawfile");
} // end namespace hmpid
} // end namespace o2



#endif /* DETECTORS_HMPID_WORKFLOW_INCLUDE_HMPIDWORKFLOW_RAWFILEREADERSPEC_H_ */
