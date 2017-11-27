#ifndef O2_MC_APPLICATION_H
#define O2_MC_APPLICATION_H

#include "FairMCApplication.h"
#include "Rtypes.h"                     // for Int_t, Bool_t, Double_t, etc
#include <iostream>

namespace o2 {
namespace steer {

// O2 specific changes/overrides to FairMCApplication
// (like for device based processing in which we
//  forward the data instead of using FairRootManager::Fill())
class O2MCApplication : public FairMCApplication
{
  public:
    using FairMCApplication::FairMCApplication;
    ~O2MCApplication() override = default;
    
    /** Define actions at the end of event */
    void  FinishEvent() override {
      // insead of filling the tree we can just send the data
      // for instance
      std::cout << "FINISHING EVENT " << std::endl;
    }
    
    /** Define actions at the end of run */
    void  FinishRun();

    /** Generate primary particles */
    void  GeneratePrimaries() override {
      std::cout << "Generating primaries" << std::endl;
      FairMCApplication::GeneratePrimaries();
    }

    ///** Initialize MC engine */
    //void                  InitMC(const char* setup,  const char* cuts);
    
    ClassDefOverride(O2MCApplication,4)  //Interface to MonteCarlo application
};

} // end namespace steer
} // end namespace o2

#endif
