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
//ROOT stuff
#include "TH2D.h"


//lumi stuff^
#include "CalibCalorimetry/CaloTPG/src/CaloTPGTranscoderULUT.h"
#include "CalibFormats/HcalObjects/interface/HcalTPGRecord.h"
#include "CalibFormats/HcalObjects/interface/HcalTPGCoder.h"
#include "CalibFormats/CaloTPG/interface/CaloTPGRecord.h"
#include <iostream>
#include <vector>

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
    TH2D* HFEtPos;
    TH2D* HFEtNeg;
    TH1D* HFEt;
    TH1D* HFPosJets;
    TH1D* HFNegJets;
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
//Generic Constructor: (const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup) I'm choosing x to be ieta and y to iphi 
    HFEtPos = fs->make<TH2D>("HFEtPos", "Forward HF Et", 40, 0.0, 40, 72, 0.0, 72);
    HFEtNeg = fs->make<TH2D>("HFEtNeg", "Backward HF Et", 40, 0.0, 40, 72, 0.0, 72);
    HFEt = fs->make<TH1D>("HFEt", "Et in HF", 100, 0.0, 100);
    HFPosJets =fs->make<TH1D>("HFPosJets", "Positive HF 3x3 Jet Et", 100, 0.0, 100);
    HFNegJets = fs->make<TH1D>("HFNegJets", "Negative HF 3x3 Jet Et", 100, 0.0, 100);
//    PosHFHits = fs->make<vector<TH2D>>();
//    NegHFHits = fs->make<vector<TH2D>>();   This is how it is done...
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
    double PosJets [40][72];
    double NegJets [40][72];
    bool nonZero = false;

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

    // Make Histograms
    edm::Service<TFileService> fs;
    // TODO-ANDREW: You will probably want to rename these so they don't all have the same name!

    for(HcalTrigPrimDigiCollection::const_iterator tp = hfpr_digi->begin(); tp != hfpr_digi->end(); ++tp)
    {
        //         cout << "HI";                    //this cout statement does not display
        if(abs(tp->id().ieta()) < 30 || abs(tp->id().ieta()) > 39)
        {
            //             cout << "useful part of HF \n";
            continue;    //Grabs the useful part of HF
        }
        //We need the sample of interest (SOI)

        HFEt->Fill(outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt()));

        if(tp->id().ieta() < 0)
        {
            EtArrayNeg [abs(tp->id().ieta())][tp->id().iphi()] = outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt());
            EtArrayNeg [abs(tp->id().ieta())][tp->id().iphi()+1] = outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt());
            HFEtNeg->SetBinContent(abs(tp->id().ieta()),tp->id().iphi(), outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt()));
            HFEtNeg->SetBinContent(abs(tp->id().ieta()),tp->id().iphi()+1, outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt()));
        }
        else
        {
            EtArrayPos [tp->id().ieta()][tp->id().iphi()] = outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt());
            EtArrayPos [tp->id().ieta()][tp->id().iphi()+1] = outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt());
            HFEtPos->SetBinContent(tp->id().ieta(),tp->id().iphi(), outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt()));
            HFEtPos->SetBinContent(tp->id().ieta(),tp->id().iphi()+1, outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt()));
        }
    /*    cout << tp->id().ieta() << " and " << tp->id().iphi() <<endl;
        cout << EtArrayPos [tp->id().ieta()][tp->id().iphi()] <<endl;
        cout << EtArrayNeg [abs(tp->id().ieta())][tp->id().iphi()] << endl << endl; */

        if(EtArrayPos [tp->id().ieta()][tp->id().iphi()] + EtArrayNeg [abs(tp->id().ieta())][tp->id().iphi()] < 0) //to make compiler happy
        {
            cout << "uh oh! negative number! :(" << endl;
        }

        if(outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt()) > .01)
        {
            nonZero = true;
        }
    }
    //     cout << "\nJust after loop \n";
    if(nonZero)
    {
        HFEtNeg = fs->make<TH2D>("HFEtNeg", "Backward HF Et", 40, 0.0, 40, 72, 0.0, 72);
        HFEtPos = fs->make<TH2D>("HFEtPos", "Forward HF Et", 40, 0.0, 40, 72, 0.0, 72);
    }

    //Now that the DIGIs are read in and graphed, this part of the code actually tries to find jets.  It makes a 3x6 ietaxiphi sliding window of jets
    if(true)
    {    
        
        for(int i = 31; i < 39; i++) //loops over ieta
        {
            for(int j = 3; j < 70; j+=2) //loops over iphi
            {
                NegJets [i][j] = EtArrayNeg [i-1][j-2] + EtArrayNeg [i-1][j] + EtArrayNeg [i-1][j+2] + 
                                 EtArrayNeg [i][j-2]   + EtArrayNeg [i][j]   + EtArrayNeg [i][j+2]   +
                                 EtArrayNeg [i+1][j-2] + EtArrayNeg [i+1][j] + EtArrayNeg [i+1][j+2];
                
                HFNegJets->Fill(NegJets [i][j]);

                PosJets [i][j] = EtArrayPos [i-1][j-2] + EtArrayPos [i-1][j] + EtArrayPos [i-1][j+2] +
                                 EtArrayPos [i][j-2]   + EtArrayPos [i][j]   + EtArrayPos [i][j+2]   +
                                 EtArrayPos [i+1][j-2] + EtArrayPos [i+1][j] + EtArrayPos [i+1][j+2];

                HFPosJets->Fill(PosJets [i][j]);
            }
         }
    }

    

    if(PosJets [1][1] + NegJets [1][1] < 0)   //random usage of arrays to make compiler happy (designed to never print)
    {
        cout << "Happy Birthday!" << endl;
    }
                                
    
    









}
void
JetAlgorithm::beginJob() { }

// ------------ method called once each job just after ending the event loop  ------------

void
JetAlgorithm::endJob() { }

DEFINE_FWK_MODULE(JetAlgorithm);
#endif
