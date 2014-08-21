// -*- C++ -*-
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

//GenJets
#include "DataFormats/JetReco/interface/GenJet.h"


using namespace std;





//
// class declaration
//

class JetAlgorithm : public edm::EDAnalyzer
{
public:
    explicit JetAlgorithm(const edm::ParameterSet&);
    ~JetAlgorithm();

    struct Jet
    {
        float seedEt;     //Energy of center cell
        float jetEt;      //Energy of 3x3 area centerd on cell
        bool pass;        //true means its a good jet seed cell false means its not
        bool match;
        float matchEta;
        float matchPhi;
    };
    struct genJet
    {
        float Et;
        bool matchPass;
        float matchiEta;
        float matchiPhi;
    };

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
    TH1D* HFNegRatio;
    TH1D* HFPosRatio;
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
//Generic Constructor For TH2D: (const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup) I'm choosing x to be ieta and y to iphi
    HFEtPos = fs->make<TH2D>("HFEtPos", "Forward HF Et", 10, 30, 40, 72, 0.0, 72);
    HFEtNeg = fs->make<TH2D>("HFEtNeg", "Backward HF Et", 10, 30, 40, 72, 0.0, 72);
    HFEt = fs->make<TH1D>("HFEt", "Et in HF", 100, 0.0, 300);
    HFPosJets =fs->make<TH1D>("HFPosJets", "Positive HF 3x3 Jet Et", 100, 0.0, 300);
    HFNegJets = fs->make<TH1D>("HFNegJets", "Negative HF 3x3 Jet Et", 100, 0.0, 300);
    HFNegRatio = fs->make<TH1D>("HFNegRatio", "Neg HF 3x3 vs Seed Ratio", 50, 0.0, 1.1);
    HFPosRatio = fs->make<TH1D>("HFPosRatio", "Pos HF 3x3 vs Seed Ratio", 50, 0.0, 1.1);
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
    using namespace reco;

    //     cout << "test \n";              //this cout statement does display

    double EtArrayPos [40][72];
    double EtArrayNeg [40][72];
    double PosJets [40][72];
    double NegJets [40][72];
    bool nonZero = false;
    bool isPosJet = false;
    bool isNegJet = false;
    bool isData = iEvent.isRealData();
    Jet myJet = Jet();
    myJet.pass = false;
    Jet HFarrayPos [40][72];
    Jet HFarrayNeg [40][72];
    genJet myGenJet = genJet();
    myGenJet.matchPass = false;

    cout << "is data: ";
    cout << isData << endl;

    Handle<std::vector<reco::GenJet>> pIn;
    iEvent.getByLabel("ak5GenJets", pIn); //gets truth jets
    edm::Service<TFileService> fs;


    edm::ESHandle<CaloTPGTranscoder> outTranscoder;
    iSetup.get<CaloTPGRecord>().get(outTranscoder);
    outTranscoder->setup(iSetup, CaloTPGTranscoder::HcalTPG);
    //     cout << "finished outTranscoder";
    // This initialization and loop is really all I care about.  I just need to save the MaxEt in an array. by eta and phi.
    // Also, I need to convert from compressed Et to real Et: I think this is done by the outTranscoder.  This sequence may need modification pending talking to Zach
    edm::Handle<HcalTrigPrimDigiCollection> hfpr_digi;
    iEvent.getByLabel("simHcalTriggerPrimitiveDigis", hfpr_digi);

    for(unsigned int k = 0; k < pIn->size(); ++k)  //checks to see if any truth jets are in HF
    {
  //      cout << "Eta Values: ";
  //      cout << pIn->at(k).eta() << endl;
        if(pIn->at(k).eta() > 2.964 && pIn->at(k).eta() < 4.716)
        {
            isPosJet = true;
            cout << "PosJet with Eta: ";
            cout << pIn->at(k).eta() << endl;
        }
        else if(pIn->at(k).eta() < -2.964 && pIn->at(k).eta() > -4.716)
        {
            isNegJet = true;
            cout << "NegJet with Eta: ";
            cout << pIn->at(k).eta() << endl;
        }

    }
    if(isPosJet)
    {
        HFEtPos = fs->make<TH2D>("HFEtPos", "Forward HF Et", 10, 30, 40, 72, 0.0, 72);
        cout << "Writing new Histograms" << endl;
    }
    if(isNegJet)
    {
        HFEtNeg = fs->make<TH2D>("HFEtNeg", "Backward HF Et", 10, 30, 40, 72, 0.0, 72);
        cout << "Writing new Histograms" << endl;
    }

    //     cout << "Just before loop \n";
    if(hfpr_digi->begin()==hfpr_digi->end())
    {
        cout << "empty" << endl;
    }

    // Make Histograms


    if(isPosJet || isNegJet) //checks that there is an HF jet
    {
        for(HcalTrigPrimDigiCollection::const_iterator tp = hfpr_digi->begin(); tp != hfpr_digi->end(); ++tp)
        {

            if(abs(tp->id().ieta()) < 30 || abs(tp->id().ieta()) > 39)  //Doesn't use the non useful part of HF or the other parts of the detector
            {

                continue;
            }

            //We need the sample of interest (SOI)
            //Fills out the Et distribution
            HFEt->Fill(outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt()));

            if(tp->id().ieta() < 0 && isNegJet)
            {
                EtArrayNeg [abs(tp->id().ieta())][tp->id().iphi()] = outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt());
                EtArrayNeg [abs(tp->id().ieta())][tp->id().iphi()+1] = outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt());
                HFEtNeg->SetBinContent(abs(tp->id().ieta())-30,tp->id().iphi(), outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt()));
                HFEtNeg->SetBinContent(abs(tp->id().ieta())-30,tp->id().iphi()+1, outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt()));
                if(outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt()) > .01)
                {
                    cout << "Setting Bin Content for HFEtNeg" << endl;
                    cout << outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt()) << endl;
                }
            }
            else if (tp->id().ieta() > 0 && isPosJet)
            {
                EtArrayPos [tp->id().ieta()][tp->id().iphi()] = outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt());
                EtArrayPos [tp->id().ieta()][tp->id().iphi()+1] = outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt());
                HFEtPos->SetBinContent(tp->id().ieta()-30,tp->id().iphi(), outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt()));
                HFEtPos->SetBinContent(tp->id().ieta()-30,tp->id().iphi()+1, outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt()));
                if(outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt()) > .01)
                {
                    cout << "Setting Bin Content for HFEtPos" << endl;
                    cout << outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt()) << endl;
                }
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
                cout << nonZero << endl;
            }
        }
    }
    //     cout << "\nJust after loop \n";


    //Now that the DIGIs are read in and graphed, this part of the code actually tries to find jets.
    if(true)
    {

        for(int i = 31; i < 39; i++) //loops over ieta
        {
            for(int j = 3; j < 70; j+=2) //loops over iphi
            {
                if(EtArrayNeg [i][j] > 10 && isNegJet)  //10GeV seed cell threshhold
                {
                    NegJets [i][j] = EtArrayNeg [i-1][j-2] + EtArrayNeg [i-1][j] + EtArrayNeg [i-1][j+2] +
                                     EtArrayNeg [i][j-2]   + EtArrayNeg [i][j]   + EtArrayNeg [i][j+2]   +
                                     EtArrayNeg [i+1][j-2] + EtArrayNeg [i+1][j] + EtArrayNeg [i+1][j+2];

                    //checks that the seed is greater than the periphery
                    if(EtArrayNeg [i][j] > EtArrayNeg [i-1][j-2] && EtArrayNeg [i][j] > EtArrayNeg [i-1][j]   && EtArrayNeg [i][j] > EtArrayNeg [i-1][j+2] &&
                       EtArrayNeg [i][j] > EtArrayNeg [i][j-2]                                                && EtArrayNeg [i][j] > EtArrayNeg [i][j+2]   &&
                       EtArrayNeg [i][j] > EtArrayNeg [i+1][j-2] && EtArrayNeg [i][j] > EtArrayNeg [i+1][j]   && EtArrayNeg [i][j] > EtArrayNeg [i+1][j+2])
                    {
                        HFarrayNeg[i][j].jetEt = NegJets [i][j];
                        HFarrayNeg[i][j].seedEt = EtArrayNeg [i][j];
                        HFarrayNeg[i][j].pass = true;                       //initial true setting for whether it is a good jet

                    }
                    HFNegJets->Fill(NegJets [i][j]);
                    HFNegRatio->Fill(EtArrayNeg [i][j]/NegJets [i][j]);
                }
                if(EtArrayPos [i][j] > 10 && isPosJet)   //10GeV seed cell threshhold
                {
                    PosJets [i][j] = EtArrayPos [i-1][j-2] + EtArrayPos [i-1][j] + EtArrayPos [i-1][j+2] +
                                     EtArrayPos [i][j-2]   + EtArrayPos [i][j]   + EtArrayPos [i][j+2]   +
                                     EtArrayPos [i+1][j-2] + EtArrayPos [i+1][j] + EtArrayPos [i+1][j+2];

                    //checks that the seed is greater than the periphery
                    if(EtArrayPos [i][j] > EtArrayPos [i-1][j-2] && EtArrayPos [i][j] > EtArrayPos [i-1][j]   && EtArrayPos [i][j] > EtArrayPos [i-1][j+2] &&
                       EtArrayNeg [i][j] > EtArrayPos [i][j-2]                                                && EtArrayNeg [i][j] > EtArrayPos [i][j+2]   &&
                       EtArrayNeg [i][j] > EtArrayPos [i+1][j-2] && EtArrayPos [i][j] > EtArrayPos [i+1][j]   && EtArrayNeg [i][j] > EtArrayPos [i+1][j+2])
                    {
                        HFarrayPos[i][j].jetEt = PosJets [i][j];
                        HFarrayPos[i][j].seedEt = EtArrayPos [i][j];
                        HFarrayPos[i][j].pass = true;                       //initial true setting for whether it is a good jet

                    }
                    HFPosJets->Fill(PosJets [i][j]);

                    HFPosRatio->Fill(EtArrayPos [i][j]/PosJets [i][j]);
                }
            }
        }
    }
    for(int i = 31; i < 39; i++) //loops over ieta
    {
        for(int j = 3; j < 70; j+=2) //loops over iphi
        {
            if(HFarrayPos[i][j].pass)
            {
               // ->Fill(HFarrayPos[i][j].jetEt);
            }
        }
    }
    for(unsigned int k = 0; k < pIn->size(); ++k)  //checks to see if any truth jets are in HF
    {
             
      
  
    }













}
void
JetAlgorithm::beginJob() { }

// ------------ method called once each job just after ending the event loop  ------------

void
JetAlgorithm::endJob() { }

DEFINE_FWK_MODULE(JetAlgorithm);
#endif
