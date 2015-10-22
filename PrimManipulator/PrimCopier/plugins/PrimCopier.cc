// -*- C++ -*-
//
// Package:    PrimManipulator/PrimCopier
// Class:      PrimCopier
// 
/**\class PrimCopier PrimCopier.cc PrimManipulator/PrimCopier/plugins/PrimCopier.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
 */
//
// Original Author:  Zachary Lesko
//         Created:  Thu, 22 Oct 2015 22:11:49 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HcalDigi/interface/HcalTriggerPrimitiveDigi.h"
#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "Geometry/HcalTowerAlgo/interface/HcalTrigTowerGeometry.h"
#include "DataFormats/HcalDetId/interface/HcalTrigTowerDetId.h"

#include "CalibCalorimetry/CaloTPG/src/CaloTPGTranscoderULUT.h"
#include "CalibFormats/HcalObjects/interface/HcalTPGRecord.h"
#include "CalibFormats/HcalObjects/interface/HcalTPGCoder.h"
#include "CalibFormats/CaloTPG/interface/CaloTPGRecord.h"
//
#include <iostream>
using namespace std;

// class declaration
//

class PrimCopier : public edm::EDProducer {
public:
  explicit PrimCopier(const edm::ParameterSet&);
  ~PrimCopier();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  virtual void produce(edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//

PrimCopier::PrimCopier(const edm::ParameterSet& iConfig) {
  //register your products
  /* Examples
     produces<ExampleData2>();

     //if do put with a label
     produces<ExampleData2>("label");
 
     //if you want to put into the Run
     produces<ExampleData2,InRun>();
   */
   produces<HcalTrigPrimDigiCollection>();
  //now do what ever other initialization is needed

}

PrimCopier::~PrimCopier() {

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------

void
PrimCopier::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  cout<<"test 1"<<endl;
  
   edm::Handle<HcalTrigPrimDigiCollection> hfpr_digi;

  std::auto_ptr<HcalTrigPrimDigiCollection> CopiedTP(new HcalTrigPrimDigiCollection());
 CopiedTP->reserve(56*72+18*8);
 cout<<"test 2"<<endl;
  iEvent.getByLabel("simHcalTriggerPrimitiveDigis", hfpr_digi);
  for (HcalTrigPrimDigiCollection::const_iterator tp = hfpr_digi->begin(); tp != hfpr_digi->end(); ++tp) {
  CopiedTP->push_back(*tp);
  }
  cout<<"test 3"<<endl;
  iEvent.put(CopiedTP);
  cout<<"test 4"<<endl;
  
  /* This is an event example
     //Read 'ExampleData' from the Event
     Handle<ExampleData> pIn;
     iEvent.getByLabel("example",pIn);

     //Use the ExampleData to create an ExampleData2 which 
     // is put into the Event
     std::auto_ptr<ExampleData2> pOut(new ExampleData2(*pIn));
     iEvent.put(pOut);
   */

  /* this is an EventSetup example
     //Read SetupData from the SetupRecord in the EventSetup
     ESHandle<SetupData> pSetup;
     iSetup.get<SetupRecord>().get(pSetup);
   */

}

// ------------ method called once each job just before starting event loop  ------------

void
PrimCopier::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------

void
PrimCopier::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
PrimCopier::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
 */

// ------------ method called when ending the processing of a run  ------------
/*
void
PrimCopier::endRun(edm::Run const&, edm::EventSetup const&)
{
}
 */

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
PrimCopier::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
 */

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
PrimCopier::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
 */

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------

void
PrimCopier::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PrimCopier);
