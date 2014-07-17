// -*- C++ -*-
//
// Package:    TestEDA/LSEn
// Class:      LSEn
// 
/**\class LSEn LSEn.cc TestEDA/LSEn/plugins/LSEn.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
 */
//
// Original Author:  Zachary Lesko
//         Created:  Fri, 13 Jun 2014 19:49:06 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HcalDigi/interface/HcalUpgradeQIESample.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
//
// class declaration
//

#include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "FWCore/Framework/interface/EventSetup.h"


#include "DataFormats/HcalDigi/interface/HFDataFrame.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "CalibFormats/HcalObjects/interface/HcalCoderDb.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TH2.h"
#include <iostream>
using namespace std;

class LSEn : public edm::EDAnalyzer
{
public:
    explicit LSEn(const edm::ParameterSet&);
    ~LSEn();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;
    TH2 *LgStEDep;

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

LSEn::LSEn(const edm::ParameterSet& iConfig) {
    edm::Service<TFileService> fs;
    LgStEDep=fs->make<TH2F>("Long_vrs_Short_Energy","Long_vrs_Short_Energy;LongEnergy;ShortEnergy",150,50,1550,150,0,1500);
    //now do what ever initialization is needed
}

LSEn::~LSEn() {

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------

void
LSEn::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    edm::ESHandle<HcalDbService> conditions;

    iSetup.get<HcalDbRecord>().get(conditions);

    
  edm::Handle<HFDigiCollection>   hf_digi;

  if(!iEvent.getByLabel("hcalDigis",hf_digi)) return;
  
    for(unsigned ihit = 0; ihit < hf_digi->size(); ++ihit)
    {
        const HFDataFrame& longfib = (*hf_digi)[ihit]; //should be long fiber short depth
        
        
        
        
        if(longfib.id().depth() != 1||abs(longfib.id().ieta())<30||abs(longfib.id().ieta())>39)continue; //so this makes it so that if we are not at depth zero, that it will continue, so we don't double
        //int ieta = longfib.id().ieta();
        //int iphi = longfib.id().iphi();
        
        HcalDetId shortFibId(longfib.id().subdet(), longfib.id().ieta(), longfib.id().iphi(), 2); // this should be short fiber id
        
        HFDigiCollection::const_iterator shortFib = hf_digi->find(shortFibId);
        //cout<<"test 6"<<endl;
        //if(shortFib->id().ieta()==0)continue;
        if(shortFib == hf_digi->end())continue;
         //cout<<"test 5"<<endl;
         //cout<<"depth"<<longfib.id().depth()<<endl;
         
         if(shortFib->id().ieta()!=longfib.id().ieta()||shortFib->id().iphi()!=longfib.id().iphi()||shortFib->id().depth()!=2)continue;
        //should do the HCalQIESample for the short fiber
        //inputs data into digi.
        //cout<<"iphi :"<<shortFib->id().iphi()<<"  and ieta :"<<longfib.id().ieta()<<endl;
        const HcalCalibrations& calibrations = conditions->getHcalCalibrations(shortFib->id());
        
        // so this next part is a guess not 100 percent sure how to make this work
        const HcalQIECoder* channelCoderS = conditions->getHcalCoder(shortFib->id());
        
        if(channelCoderS==0)
        {
            cout<<"channelCoderS "<<endl;
            return;
        }
        const HcalQIEShape* shapeS = conditions->getHcalShape(channelCoderS);

        if(shapeS==0)
        {cout<<" shapeS"<<endl;
            return;
        }
        HcalCoderDb coders(*channelCoderS, *shapeS); //magic piece of code

  
        CaloSamples tools;
        coders.adc2fC((*shortFib), tools); // this step takes short fiber digi and tools and fills tools with the linadc
        


        //seperation between long and short



        //all this should be long, hopefully gets the long energy

        //so I think I need this part but I am unsure, depends what it does,
        const HcalCalibrations& calibrationL = conditions->getHcalCalibrations(longfib.id());

        CaloSamples tool_l;

        const HcalQIECoder* channelCoderL = conditions->getHcalCoder(longfib.id());
        if(channelCoderL==0)
        {
            cout<<" channelCoderL "<<endl;
            return;
        }
        const HcalQIEShape* shapeL = conditions->getHcalShape(channelCoderL);
        if(shapeL==0)
        {
            cout<<" shapeL "<<endl;
            return;
        }

        HcalCoderDb coderL(*channelCoderL, *shapeL);
        coderL.adc2fC(longfib, tool_l); // this fills tool_l[0] with linearized adc
        int MaxLE = 0;
        int MaxSE = 0;
        for(int i = 0; i < longfib.size() && i < shortFib->size(); i++)
        {
            float LongE = (tool_l[i] - calibrationL.pedestal(longfib[i].capid())) * calibrationL.respcorrgain(longfib[i].capid());
            float ShortE = (tools[i] - calibrations.pedestal((*shortFib)[i].capid())) * calibrations.respcorrgain((*shortFib)[i].capid());
            //if(ShortE!=0)cout<<"so this is the tools[i] :"<<ShortE<<endl;
            if(MaxLE<LongE)MaxLE=LongE;
            if(MaxSE<ShortE)MaxSE=ShortE;
        }
        LgStEDep->Fill(MaxLE,MaxSE);

    }
}


// ------------ method called once each job just before starting event loop  ------------

void
LSEn::beginJob() { }

// ------------ method called once each job just after ending the event loop  ------------

void
LSEn::endJob() { }

// ------------ method called when starting to processes a run  ------------
/*
void 
LSEn::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
 */

// ------------ method called when ending the processing of a run  ------------
/*
void 
LSEn::endRun(edm::Run const&, edm::EventSetup const&)
{
}
 */

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
LSEn::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
 */

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
LSEn::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
 */

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------

void
LSEn::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(LSEn);
