// -*- C++ -*-
//
// Package:    TestEDA/tester
// Class:      tester
// 
/**\class tester tester.cc TestEDA/tester/plugins/tester.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
 */
//
// Original Author:  Zachary Lesko
//         Created:  Mon, 12 May 2014 15:26:07 GMT
//
//
#ifndef EdAnalizer_Tester
#define EdAnalizer_Tester

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

#include "CondFormats/HcalObjects/interface/HcalQIECoder.h"
#include "CalibFormats/HcalObjects/interface/HcalCoderDb.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
//and this is where everything is tested. should have differnt
#include "DataFormats/HcalDigi/interface/HcalTriggerPrimitiveDigi.h"

#include <memory>
#include <vector>
#include <map>
#include <utility>
#define _USE_MATH_DEFINES 
#include <math.h>
#include <TDirectory.h> 
//obviously trying to get directories to work so I could look at the electrons in the correct regions
#include <TSystem.h>
//maybe
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "TestEDA/CutAssessment/interface/Genconverter.h"

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include <iostream>
using namespace std;
//
// class declaration
//

class tester : public edm::EDAnalyzer
{
public:
    explicit tester(const edm::ParameterSet&);
    ~tester();

    //static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;
    //gSystem->mkdir("remember");
    Genconverter* Real2I;
    TFile* makingtest;
    std::string fname;
    TH2 *LSDeposited;
    bool ElecIniEvent;
    int TotalReal;

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

tester::tester(const edm::ParameterSet& iConfig)
{
    //TSystem      *gSystem   = 0;
    //string ofile = "progresses";
    //gSystem->mkdir((ofile).c_str());

    //makingtest = new TFile("Thow.root", "RECREATE");
    edm::Service<TFileService> fs;
    TFileDirectory passed = fs->mkdir("testing");
    LSDeposited = passed.make<TH2F>("Long_fiber_vrs_ShortEDeposit", "long fiber vrs short; long fiber EDep; short fiber EDep", 200, 0, 100, 200, 0, 100);
    //LSDeposited = new TH2F("short_fiber_vrs_long_EDeposit", "short fiber vrs long; long fiber EDep; short fiber EDep", 200, 0, 50, 200, 0, 50);
}
//now do what ever initialization is needed

tester::~tester() {

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------

void
tester::analyze(const edm::Event& iEvent, const edm::EventSetup& eventsetup)
{
    ElecIniEvent = false;
    using namespace edm;
    int numElect = 0;
    set<pair<int, int> > gparts; //makes a set that contains locations of all generated electrons
    edm::Handle<reco::GenParticleCollection> genInfo;
    iEvent.getByLabel("genParticles", genInfo);
    Real2I = new Genconverter();




    for(unsigned ighit = 0; ighit < genInfo->size(); ++ighit)// so this function goes and figures out the location of all electron location to allow matching4
    {
        const reco::GenParticle& genin = (* genInfo)[ighit];
        int giphi;

        giphi = Real2I->Phi2Iphi(genin.phi(), genin.eta()); //so this converts into degrees since everything is nicer there. followed by conversion it iphi




        int gieta = Real2I->Eta2IEta(genin.eta()); //the function that gets ieta from eta(genin.eta());
        int pdgid = abs(genin.pdgId());






        if(pdgid == 11 && 1 == genin.status() && genin.pt() > 10 && abs(gieta) > 28 && abs(gieta) < 42)
        {


            numElect++;
            ElecIniEvent = true;
            TotalReal++; // counts total number of real electrons
            if(gieta > 0) gparts.insert(make_pair(gieta - 30, giphi)); // should make say whether a electron hit the sensor
            else gparts.insert(make_pair(gieta + 30, giphi));
        }

        //end generated stuff
        using namespace edm;


        edm::ESHandle < HcalDbService > pSetup;
        eventsetup.get<HcalDbRecord> ().get(pSetup);

        edm::Handle<HFDigiCollection> hf_digi; //
        edm::Handle<HcalTrigPrimDigiCollection> hfpr_digi;
        iEvent.getByLabel("simHcalUnsuppressedDigis", hf_digi);
        iEvent.getByLabel("simHcalTriggerPrimitiveDigis", hfpr_digi);


        //int primdigi = 0;
        //int realdigi = 0;
        bool primedigi[100][72][5];
        bool normdigi[100][72][5];
        for(HcalTrigPrimDigiCollection::const_iterator tp = hfpr_digi->begin(); tp != hfpr_digi->end(); ++tp)
        {






            for(int i = 0; i < 4; i++)
            {
                primedigi[tp->id().ieta() + 50][tp->id().iphi()][i] = (*tp)[i].fineGrain();
            }
        }




        //probably going to get rid of this.
        for(unsigned ihit = 0; ihit < hf_digi->size(); ++ihit)
        {
            const HFDataFrame& LongDigi = (*hf_digi)[ihit];
            int ieta = LongDigi.id().ieta();
            int iphi = LongDigi.id().iphi();
            //int depth = LongDigi.id().depth();


            if((LongDigi.id().depth() != 1) || (abs(ieta) < 30 || abs(ieta) > 39))continue;

            HcalDetId bob(LongDigi.id().subdet(), ieta, iphi, 2);
            // HFDigiCollection::const_iterator theit = hf_digi->find(bob);


            HcalDetId shortdetid(LongDigi.id().subdet(), LongDigi.id().ieta(), LongDigi.id().iphi(), 2);
            HFDigiCollection::const_iterator ShortDigi = hf_digi->find(bob);
            if(LongDigi.id().ieta()!=(*ShortDigi).id().ieta()||LongDigi.id().iphi()!=(*ShortDigi).id().iphi()||(*ShortDigi).id().depth()!=2)continue;

            const HcalCalibrations& calibrationl = (*pSetup).getHcalCalibrations(LongDigi.id());





            double ShortEE;
            double LongEE;
            double ShortLadc;
            double LongLadc;
            for(int w = 0; w < 5; w++)
            {


                //mess for getting the bloody stupid lin adc
                CaloSamples tool_l;

                const HcalQIECoder* channelCoderL = pSetup->getHcalCoder(LongDigi.id());

                const HcalQIEShape* shapeL = (*pSetup).getHcalShape(channelCoderL);

                HcalCoderDb coderL(*channelCoderL, *shapeL);

                coderL.adc2fC(LongDigi, tool_l);


                const HcalCalibrations& calibrations = (*pSetup).getHcalCalibrations(shortdetid);

                const HcalQIECoder* channelCoderS = pSetup->getHcalCoder(shortdetid);

                const HcalQIEShape* shapeS = (*pSetup).getHcalShape(channelCoderS);

                HcalCoderDb coders(*channelCoderS, *shapeS);

                CaloSamples tools;
                coders.adc2fC((*ShortDigi), tools);
                //


                double ShortE = (tools[w] - calibrations.pedestal((*ShortDigi)[w].capid())) * calibrations.respcorrgain((*ShortDigi)[w].capid());

                double LongE = (tool_l[w] - calibrationl.pedestal(LongDigi[w].capid())) * calibrationl.respcorrgain(LongDigi[w].capid());

                normdigi[ieta + 50][iphi][w] = ((ShortE < (LongE - 10.1) * 10.2) && (LongE > 10) && (ShortE > 11));
                ShortLadc = tools[2];
                LongLadc = tool_l[2];

                ShortEE = (tools[2] - calibrations.pedestal((*ShortDigi)[2].capid())) * calibrations.respcorrgain((*ShortDigi)[2].capid());

                LongEE = (tool_l[2] - calibrationl.pedestal(LongDigi[2].capid())) * calibrationl.respcorrgain(LongDigi[2].capid());


            }
            LSDeposited->Fill(LongEE, ShortEE);

            for(int i = 0; i < 4; i++)
            {
                if(normdigi[ieta + 50][iphi][i] != primedigi[ieta + 50][iphi][i]&&false)
                {
                    cout << "okay testing, short  adc " << (*ShortDigi)[0].adc() << endl;
                    cout << " and our long adc is " << LongDigi[0].adc() << endl;
                    cout << "well it broke.  well then ieta is " << ieta << " our phi is " << iphi << endl;
                    cout << " so which is broke I guess becomes the question. norm digi is" << normdigi[ieta + 50][iphi][i] << endl;
                    cout << " and our primedigi is " << primedigi[ieta + 50][iphi][i] << endl;
                    cout << "and our energy is long " << LongEE << endl;
                    cout << "energy short is " << ShortEE << endl;
                    cout << "and last but not least Longfc " << LongLadc << endl;
                    cout << "and our shortfc is " << ShortLadc << endl;
                    int j;
                    cin >> j;

                }
            }
        }
    }
}


/*
#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}
 */

// ------------ method called once each job just before starting event loop  ------------

void
tester::beginJob() { }

// ------------ method called once each job just after ending the event loop  ------------

void
tester::endJob() { }

// ------------ method called when starting to processes a run  ------------
/*
void 
tester::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
 */

// ------------ method called when ending the processing of a run  ------------
/*
void 
tester::endRun(edm::Run const&, edm::EventSetup const&)
{
}
 */

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
tester::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
 */

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
tester::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
 */

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
/*void
tester::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
 */
//define this as a plug-in
DEFINE_FWK_MODULE(tester);
#endif