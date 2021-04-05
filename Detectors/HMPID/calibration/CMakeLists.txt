o2_add_library(HMPIDCalibration
               SOURCES  src/HMPIDCalibCollector.cxx
			src/HMPIDDCSProcessor.cxx
			src/HMPIDChannelCalibrator.cxx
               PUBLIC_LINK_LIBRARIES O2::DataFormatsHMP O2::HMPIDBase
				     O2::CCDB
				     O2::DetectorsCalibration
                                     O2::DetectorsDCS)

o2_add_executable(hmpid-dcs-sim-workflow
                  COMPONENT_NAME calibration
                  SOURCES testWorkflow/hmpid-dcs-sim-workflow.cxx
                  PUBLIC_LINK_LIBRARIES O2::DCStestWorkflow)
