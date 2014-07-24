// -*- C++ -*-
//
// Package:    TestEDA/JetAlgorithm
// Class:      JetAlgorithm
// 
/**\class JetAlgorithm JetAlgorithm.cc TestEDA/DataAnalizer/src/JetAlgorithm.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
 */
//
// Original Author:  Andrew Evans with lots of borrowing from Zachary Lesko's DataAnalizer.cc
//         Created:  Thu, 12 Jun 2014 15:30:54 GMT
//
//
#ifndef JetAlgorithm_included
#define JetAlgorithm_included 1


// system include files
#include <memory>
#include <utility>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TFile.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/HcalDigi/interface/HcalTriggerPrimitiveDigi.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
//Lumi stuff
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/NoDataException.h"
#include "FWCore/Framework/interface/NoRecordException.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/EventSetupRecord.h"
#include "RecoLuminosity/LumiProducer/interface/LumiCorrectionParam.h"
#include "RecoLuminosity/LumiProducer/interface/LumiCorrectionParamRcd.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"




//lumi stuff^
#include "CalibCalorimetry/CaloTPG/src/CaloTPGTranscoderULUT.h"
#include "CalibFormats/HcalObjects/interface/HcalTPGRecord.h"
#include "CalibFormats/HcalObjects/interface/HcalTPGCoder.h"
#include "CalibFormats/CaloTPG/interface/CaloTPGRecord.h"
#include <iostream>
using namespace std;





//
// class declaration
//

class JetAlgorithm : public edm::EDAnalyzer
{
public:
    explicit JetAlgorithm(const edm::ParameterSet&);
    ~JetAlgorithm();



private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;


    //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

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

JetAlgorithm::JetAlgorithm(const edm::ParameterSet& iConfig)
{
    edm::Service<TFileService> fs;
    //now do what ever initialization is needed
}

JetAlgorithm::~JetAlgorithm() {

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------

void
JetAlgorithm::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
//     cout << "beginning";
     using namespace edm;

//     cout << "test \n";              //this cout statement does display

     double EtArrayPos [40][72];
     double EtArrayNeg [40][72];

     edm::ESHandle<CaloTPGTranscoder> outTranscoder;
     iSetup.get<CaloTPGRecord>().get(outTranscoder);
     outTranscoder->setup(iSetup, CaloTPGTranscoder::HcalTPG);
//     cout << "finished outTranscoder";
     // This initialization and loop is really all I care about.  I just need to save the MaxEt in an array. by eta and phi.
     // Also, I need to convert from compressed Et to real Et: I think this is done by the outTranscoder.  This sequence may need modification pending talking to Zach
     edm::Handle<HcalTrigPrimDigiCollection> hfpr_digi;
     iEvent.getByLabel("simHcalTriggerPrimitiveDigis", hfpr_digi);
//     cout << "Just before loop \n";
     if(hfpr_digi->begin()==hfpr_digi->end())
     {
         cout << "empty";
     }

     
     for(HcalTrigPrimDigiCollection::const_iterator tp = hfpr_digi->begin(); tp != hfpr_digi->end(); ++tp)
     {
//         cout << "HI";                    //this cout statement does not display
         if(abs(tp->id().ieta()) < 30 || abs(tp->id().ieta()) > 39)
         {
//             cout << "useful part of HF \n";
             continue;    //Grabs the useful part of HF
         }
         //We need the sample of interest (SOI)
         if(tp->id().ieta() < 0)
	 {
             EtArrayNeg [abs(tp->id().ieta())][tp->id().iphi()] = outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt());
	 }
	 else
	 {
	     EtArrayPos [tp->id().ieta()][tp->id().iphi()] = outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt());
         }
	 cout << tp->id().ieta() << " and " << tp->id().iphi() <<endl;
	 cout << EtArrayPos [1][1]*EtArrayNeg [1][1]*0;          //this cout statement does not display
     }
//     cout << "\nJust after loop \n";
}
void
JetAlgorithm::beginJob() { }

// ------------ method called once each job just after ending the event loop  ------------

void
JetAlgorithm::endJob() { }

DEFINE_FWK_MODULE(JetAlgorithm);
#endif
