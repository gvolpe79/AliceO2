// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef O2_DEVICES_PRIMSERVDEVICE_H_
#define O2_DEVICES_PRIMSERVDEVICE_H_

#include <FairMQDevice.h>
#include <FairLogger.h>
#include <FairGenerator.h>
#include <FairBoxGenerator.h>
#include <FairPrimaryGenerator.h>
#include <FairMQMessage.h>
#include <SimulationDataFormat/Stack.h>
#include <FairMCEventHeader.h>
#include <TMessage.h>
#include <TClass.h>
#include <SimulationDataFormat/PrimaryChunk.h>
#include <Generators/GeneratorFromFile.h>
#include <Configuration/SimConfig.h>
#include <typeinfo>

namespace o2 {
namespace devices {

class O2PrimaryServerDevice : public FairMQDevice
{
  public:
    /// Default constructor
    O2PrimaryServerDevice() {
      OnData("primary-get", &O2PrimaryServerDevice::HandleRequest);
    }

    /// Default destructor
    ~O2PrimaryServerDevice() final = default;

  protected:
    void InitTask() final {
      LOG(INFO) << "Init Server device " << FairLogger::endl;

      auto boxGen = new FairBoxGenerator(211, 1000); /*protons*/
      boxGen->SetEtaRange(-0.9, 0.9);
      boxGen->SetPRange(0.1, 5);
      boxGen->SetPhiRange(0., 360.);
      boxGen->SetDebug(kFALSE);
      mPrimGen.AddGenerator(boxGen);
      mPrimGen.SetEvent(&mEventHeader);

      //auto extGen =  new o2::eventgen::GeneratorFromFile(confref.getExtKinematicsFileName().c_str());
      //auto extGen =  new o2::eventgen::GeneratorFromFile("~/Downloads/Kinematics_HijingCent1_N100.root");
      //extGen->SetStartEvent(0);
      //mPrimGen.AddGenerator(extGen);
      //mPrimGen.SetEvent(&mEventHeader);
    }

    /// Overloads the ConditionalRun() method of FairMQDevice
    bool HandleRequest(FairMQMessagePtr& request, int /*index*/)
    {
      static int counter = 0;
      if(counter >= mMaxEvents && mNeedNewEvent) return false; // max 3 events for the moment

      LOG(INFO) << "Received request for work " << FairLogger::endl;

      if (mNeedNewEvent) {
        mStack.Reset();
        mPrimGen.GenerateEvent(&mStack);
        mNeedNewEvent = false;
	mPartCounter = 0;
        counter++;
      }



      auto& prims = mStack.getPrimaries();
      auto numberofparts = (int)std::ceil( prims.size() / (1.*mChunkGranularity));

      o2::Data::PrimaryChunk m;
      o2::Data::SubEventInfo i;
      i.eventID = counter;
      i.maxEvents = mMaxEvents;
      i.part = mPartCounter + 1;
      i.nparts = numberofparts;
      i.seed = counter + 10;
      i.index = m.mParticles.size();
      m.mEventIDs.emplace_back(i);

      auto startindex = mPartCounter*mChunkGranularity;
      auto endindex = startindex + mChunkGranularity;
      auto startiter = prims.begin() + mPartCounter*mChunkGranularity;
      auto enditer = endindex < prims.size() ? startiter + mChunkGranularity : prims.end();
      std::copy(startiter, enditer, std::back_inserter(m.mParticles));

      LOG(WARNING) << "Sending " << m.mParticles.size() << " particles\n";
      LOG(WARNING) << "treating ev " << counter << " part " << i.part << " out of " << i.nparts << "\n";

      mPartCounter++;
      if (mPartCounter == numberofparts) {
        mNeedNewEvent = true;
      }

      TMessage* tmsg = new TMessage(kMESS_OBJECT);
      tmsg->WriteObjectAny((void*)&m, TClass::GetClass("o2::Data::PrimaryChunk"));

      auto free_tmessage = [](void* data, void* hint) { delete static_cast<TMessage*>(hint); };

      std::unique_ptr<FairMQMessage> message(
        fTransportFactory->CreateMessage(tmsg->Buffer(), tmsg->BufferSize(), free_tmessage, tmsg));

      // send answer
      if (Send(message, "primary-get") > 0) {
        LOG(INFO) << "reply send " << FairLogger::endl;
        return true;
      }
      return true;
    }

   private:
    std::string mOutChannelName = "";
    FairPrimaryGenerator mPrimGen;
    FairMCEventHeader mEventHeader;
    o2::Data::Stack mStack; // the stack which is filled
    int mChunkGranularity = 500; // how many primaries to send to a worker
    int mLastPosition = 0; // last position in stack vector
    int mPartCounter = 0;
    bool mNeedNewEvent = true;
    int mMaxEvents = 2;
};

} // namespace devices
} // namespace o2

#endif
