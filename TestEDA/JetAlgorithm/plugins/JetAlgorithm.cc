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
//deltaR matching
#include "DataFormats/Math/interface/deltaR.h"
//C++
#include <math.h>
//Zach's Converter
#include "TestEDA/CutAssessment/interface/Genconverter.h"



using namespace std;





//
// class declaration
//

class JetAlgorithm : public edm::EDAnalyzer
{
    public:
        explicit JetAlgorithm(const edm::ParameterSet&);
        ~JetAlgorithm();

        int genJetsNegHF = 0;
        int genJetsPosHF = 0;
        int trigNegJets = 0;
        int trigPosJets = 0;

        struct Jet
        {
            float seedEt;     //Energy of center cell
            float jetEt;      //Energy of 3x3 area centerd on cell
            bool pass;        //true means its a good jet seed cell false means its not
            bool match;
            float matchEta;
            float matchPhi;
            int genId;
        };
        struct genJet
        {
            double Et;
            bool matchPass;
            float matchiEta;
            float matchiPhi;
            float eta;
            float phi;
            int id;
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

        TH1D* Et;
        TH1D* Eta;
        TH1D* Phi;
        TH1D* IPhi;
        TH1D* absEta;
        TH1D* EtMatched;
        TH1D* EtaMatched;
        TH1D* PhiMatched;
        TH1D* IPhiMatched;
        TH1D* absEtaMatched;
        TH1D* trigJetEt;
        TH1D* trigJetEta;
        TH1D* trigJetPhi;
        TH1D* trigJetIPhi;

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
    HFEt = fs->make<TH1D>("HFEt", "GenJet Et", 100, 0.0, 100);
    HFPosJets = fs->make<TH1D>("HFPosJets", "Positive HF 3x3 Jet Et", 100, 0.0, 300);
    HFNegJets = fs->make<TH1D>("HFNegJets", "Negative HF 3x3 Jet Et", 100, 0.0, 300);
    HFNegRatio = fs->make<TH1D>("HFNegRatio", "Neg HF 3x3 vs Seed Ratio", 50, 0.0, 1.1);
    HFPosRatio = fs->make<TH1D>("HFPosRatio", "Pos HF 3x3 vs Seed Ratio", 50, 0.0, 1.1);

    Et = fs->make<TH1D>("Et", "genJet Et", 100, 0.0, 100.0);
    Eta = fs->make<TH1D>("Eta", "GenJet Eta", 50, -5.0, 5.0);
    Phi = fs->make<TH1D>("Phi", "GenJet Phi", 100, -5.0, 5.0);
    IPhi = fs->make<TH1D>("IPhi", "GenJet IPhi", 76, -1.0, 74);
    absEta = fs->make<TH1D>("absEta", "GenJet |eta|", 100, 0, 5.0);
    EtMatched = fs->make<TH1D>("EtMatched", "Matched genJet Et", 100, 0.0, 100.0);
    EtaMatched = fs->make<TH1D>("EtaMatched", "Matched genJet Eta", 50, -5.0, 5.0);
    PhiMatched = fs->make<TH1D>("PhiMatched", "Matched genJet Phi", 100, -5.0, 5.0);
    IPhiMatched = fs->make<TH1D>("IPhiMatched", "Matched genJet IPhi", 76, -1.0, 74);
    absEtaMatched = fs->make<TH1D>("absEtaMatched", "Matched genJet |eta|", 100, 0, 5.0);
    trigJetEt = fs->make<TH1D>("trigJetEt", "trigger Jet Et", 100, 0.0, 100.0);
    trigJetEta = fs->make<TH1D>("trigJetEta", "trigger Jet Eta", 50, -5.0, 5.0);
    trigJetPhi = fs->make<TH1D>("trigJetPhi", "trigger Jet Phi", 100, -5.0, 5.0);
    trigJetIPhi = fs->make<TH1D>("trigJetIPhi", "trigger Jet IPhi", 76, -1.0, 74);


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
    using namespace edm;
    using namespace reco;

    double EtArrayPos [40][72];
    double EtArrayNeg [40][72];
    double PosJets [40][72];
    double NegJets [40][72];
    double dR;

    int passed;
    int left;
    int right;
    //    bool nonZero = false;
    bool isPosJet = false;
    bool isNegJet = false;
    bool isData = iEvent.isRealData();
    vector<genJet> negGenJets;
    vector<genJet> posGenJets;

    Jet myJet = Jet();
    myJet.pass = false;
    Jet HFarrayPos [40][72];
    Jet HFarrayNeg [40][72];
    genJet myGenJet = genJet();
    myGenJet.matchPass = false;

    if(myGenJet.matchPass)
    {
        cout << "hi";
    }
    if(myJet.pass)
    {
        cout << "hi";
    }

    if(isData)
    {
        cout << "Running on Data\n";
    }
    Handle<std::vector<reco::GenJet>> pIn;
    if(!isData)
    {
        iEvent.getByLabel("ak5GenJets", pIn); //gets truth ak5 jets
    }
    edm::Service<TFileService> fs;

    edm::ESHandle<CaloTPGTranscoder> outTranscoder;
    iSetup.get<CaloTPGRecord>().get(outTranscoder);
    outTranscoder->setup(iSetup, CaloTPGTranscoder::HcalTPG);
    //     cout << "finished outTranscoder";
    // This initialization and loop is really all I care about.  I just need to save the MaxEt in an array. by eta and phi.
    // Also, I need to convert from compressed Et to real Et: I think this is done by the outTranscoder.  This sequence may need modification pending talking to Zach
    edm::Handle<HcalTrigPrimDigiCollection> hfpr_digi;
    iEvent.getByLabel("simHcalTriggerPrimitiveDigis", hfpr_digi);

    if(!isData)
    {
        for(unsigned int k = 0; k < pIn->size(); ++k)  //loops over all ak5 jest and checks if they are in HF
        {
            //      cout << "Eta Values: ";
            //      cout << pIn->at(k).eta() << endl;
            if(pIn->at(k).et() > 10)                   //makes sure that each jet has at least 10GeV of transverse energy
            {    
                if(pIn->at(k).eta() > 2.964 && pIn->at(k).eta() < 4.716)
                {
                    isPosJet = true; 
                    genJetsPosHF++;

                    if(pIn->at(k).eta() < 4.4 && pIn->at(k).eta() > 3.4)
                    {
                        Et->Fill(pIn->at(k).et());
                        Eta->Fill(pIn->at(k).eta());
                        Phi->Fill(pIn->at(k).phi());
                        absEta->Fill(fabs(pIn->at(k).eta()));
                    }
                    IPhi->Fill(Genconverter().Phi2Iphi(pIn->at(k).phi(), pIn->at(k).eta()));

                    //          cout << "PosJet with Eta: ";
                    //        cout << pIn->at(k).eta() << endl;
                }
                if(pIn->at(k).eta() < -2.964 && pIn->at(k).eta() > -4.716)
                {
                    isNegJet = true;
                    genJetsNegHF++;

                    if(pIn->at(k).eta() > -4.4 && pIn->at(k).eta() < -3.4)
                    {
                        Et->Fill(pIn->at(k).et());
                        Eta->Fill(pIn->at(k).eta());
                        Phi->Fill(pIn->at(k).phi());
                        absEta->Fill(fabs(pIn->at(k).eta()));
                    }
                    IPhi->Fill(Genconverter().Phi2Iphi(pIn->at(k).phi(), pIn->at(k).eta()));

                    //          cout << "NegJet with Eta: ";
                    //        cout << pIn->at(k).eta() << endl;
                }
            }
        }
    }
    /*
       if(isPosJet)   //checks that somewhere in the current event, there is an ak5 jet in postive HF
       {
       HFEtPos = fs->make<TH2D>("HFEtPos", "Forward HF Et", 10, 30, 40, 72, 0.0, 72);
    //        cout << "Writing new Histograms" << endl;
    }
    if(isNegJet)   //checks the same for negative
    {
    HFEtNeg = fs->make<TH2D>("HFEtNeg", "Backward HF Et", 10, 30, 40, 72, 0.0, 72);
    //        cout << "Writing new Histograms" << endl;
    }
    */
    //     cout << "Just before loop \n";
    if(hfpr_digi->begin()==hfpr_digi->end()) //checks if the event is empty
    {
        cout << "empty" << endl;
    }


    //Now the code loops over cell by cell in the event to begin to reconstruct jets.  It only looks at events that have an ak5 genJet in HF.

    if(isPosJet || isNegJet || isData) //checks that there is an HF ak5 truth jet
    {
        for(HcalTrigPrimDigiCollection::const_iterator tp = hfpr_digi->begin(); tp != hfpr_digi->end(); ++tp)
        {

            if(abs(tp->id().ieta()) < 30 || abs(tp->id().ieta()) > 39)  //Doesn't use the part of HF which has a different cell size or the other parts of the detector
            {

                continue;
            }

            //We need the sample of interest (SOI)
            //Fills out the Et distribution
            HFEt->Fill(outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt()));

            if(tp->id().ieta() < 0 && (isNegJet || isData))
            {
                EtArrayNeg [abs(tp->id().ieta())][tp->id().iphi()] = outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt());
                EtArrayNeg [abs(tp->id().ieta())][tp->id().iphi()+1] = outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt());
                HFEtNeg->SetBinContent(abs(tp->id().ieta())-30,tp->id().iphi(), outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt()));
                HFEtNeg->SetBinContent(abs(tp->id().ieta())-30,tp->id().iphi()+1, outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt()));
                if(outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt()) > .01 && tp->id().iphi() == 71)
                {
                //    cout << "Setting Bin Content for HFEtNeg" << endl;
                //    cout << outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt()) << endl;
                }
            }
            else if (tp->id().ieta() > 0 && (isPosJet || isData))
            {
                EtArrayPos [tp->id().ieta()][tp->id().iphi()] = outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt());
                EtArrayPos [tp->id().ieta()][tp->id().iphi()+1] = outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt());
                HFEtPos->SetBinContent(tp->id().ieta()-30,tp->id().iphi(), outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt()));
                HFEtPos->SetBinContent(tp->id().ieta()-30,tp->id().iphi()+1, outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt()));
                if(outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt()) > .01)
                {
                //    cout << "Setting Bin Content for HFEtPos" << endl;
                //    cout << outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt()) << endl;
                }
            }
            /*    cout << tp->id().ieta() << " and " << tp->id().iphi() <<endl;
                  cout << EtArrayPos [tp->id().ieta()][tp->id().iphi()] <<endl;
                  cout << EtArrayNeg [abs(tp->id().ieta())][tp->id().iphi()] << endl << endl; */

            if(EtArrayPos [tp->id().ieta()][tp->id().iphi()] + EtArrayNeg [abs(tp->id().ieta())][tp->id().iphi()] < 0) //to make compiler happy
            {
                //              cout << EtArrayPos [tp->id().ieta()][tp->id().iphi()] << endl;
                //              cout << EtArrayNeg [abs(tp->id().ieta())][tp->id().iphi()] << endl;
                //              cout << "uh oh! negative number! :(" << endl;

            }

            if(outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt()) > .01)
            {
                //              nonZero = true;
                //      cout << nonZero << endl;
            }
        }
    }
    //     cout << "\nJust after loop \n";


    //Now that the DIGIs are read in and graphed, this part of the code actually tries to find jets.
    if(true)
    {

        for(int i = 31; i < 39; i++) //loops over ieta
        {
            for(int j = 1; j < 72; j+=2) //loops over iphi
            {
                left = j-2;
                right = j+2;
                if(j == 1)
                {
                    left = 71;
                    //    cout << left << endl;
                    /*     if(isNegJet && EtArrayNeg[i][j] > 0)
                           {
                           cout << EtArrayNeg [i][left] << endl;
                           }*/
                }
                else if(j == 71)
                {
                    right = 1;
                    //     cout << right << endl;
                    /*    if(isNegJet && EtArrayNeg[i][j] > 0)
                          {    
                          cout << EtArrayNeg [i][right] << endl;
                          } */   

                }

                if(EtArrayNeg [i][j] > 10 && (isNegJet || isData))  //10GeV seed cell threshhold
                {
                    NegJets [i][j] = EtArrayNeg [i-1][left] + EtArrayNeg [i-1][j] + EtArrayNeg [i-1][right] +
                        EtArrayNeg [i][left]   + EtArrayNeg [i][j]   + EtArrayNeg [i][right]   +
                        EtArrayNeg [i+1][left] + EtArrayNeg [i+1][j] + EtArrayNeg [i+1][right];

                    //checks that the seed is greater than the periphery
                    if(EtArrayNeg [i][j] > EtArrayNeg [i-1][left] && EtArrayNeg [i][j] > EtArrayNeg [i-1][j]   && EtArrayNeg [i][j] > EtArrayNeg [i-1][right] &&
                            EtArrayNeg [i][j] > EtArrayNeg [i][left]                                                && EtArrayNeg [i][j] > EtArrayNeg [i][right]   &&
                            EtArrayNeg [i][j] > EtArrayNeg [i+1][left] && EtArrayNeg [i][j] > EtArrayNeg [i+1][j]   && EtArrayNeg [i][j] > EtArrayNeg [i+1][right] 
                            && NegJets [i][j] > 10)
                    {
                        cout << "Made Trigger Jet in NegHF" << endl;
                        trigNegJets++;
                        HFarrayNeg[i][j].jetEt = NegJets [i][j];
                        HFarrayNeg[i][j].seedEt = EtArrayNeg [i][j];
                        HFarrayNeg[i][j].pass = true; //initial true setting for whether it is a good jet
                        trigJetEt->Fill(NegJets [i][j]);
                        trigJetEta->Fill(Genconverter().IEta2Eta(i));
                        trigJetPhi->Fill(Genconverter().IPhi2Phi(j));
                        trigJetIPhi->Fill(j);
                        if(HFarrayNeg[i][j].jetEt < 0)
                        {
                            cout << "hi";
                        }
                    }
                    else
                    {
                        HFarrayNeg[i][j].pass = false;
                    }
                    HFNegJets->Fill(NegJets [i][j]);
                    HFNegRatio->Fill(EtArrayNeg [i][j]/NegJets [i][j]);
                }
                else
                {
                    HFarrayNeg[i][j].pass = false;
                }
                if(EtArrayPos [i][j] > 10 && (isPosJet || isData))   //10GeV seed cell threshhold
                {
                    PosJets [i][j] = EtArrayPos [i-1][left] + EtArrayPos [i-1][j] + EtArrayPos [i-1][right] +
                        EtArrayPos [i][left]   + EtArrayPos [i][j]   + EtArrayPos [i][right]   +
                        EtArrayPos [i+1][left] + EtArrayPos [i+1][j] + EtArrayPos [i+1][right];

                    //checks that the seed is greater than the periphery
                    if(EtArrayPos [i][j] > EtArrayPos [i-1][left] && EtArrayPos [i][j] > EtArrayPos [i-1][j]   && EtArrayPos [i][j] > EtArrayPos [i-1][right] &&
                            EtArrayPos [i][j] > EtArrayPos [i][left]                                                && EtArrayPos [i][j] > EtArrayPos [i][right]   &&
                            EtArrayPos [i][j] > EtArrayPos [i+1][left] && EtArrayPos [i][j] > EtArrayPos [i+1][j]   && EtArrayPos [i][j] > EtArrayPos [i+1][right]
                            && PosJets[i][j] > 10)
                    {
                        cout << "Made Trigger Jet in PosHF" << endl;
                        trigPosJets++;
                        HFarrayPos[i][j].jetEt = PosJets [i][j];
                        HFarrayPos[i][j].seedEt = EtArrayPos [i][j];
                        HFarrayPos[i][j].pass = true;                       //initial true setting for whether it is a good jet
                        trigJetEt->Fill(PosJets [i][j]);
                        trigJetEta->Fill(Genconverter().IEta2Eta(-i));
                        trigJetPhi->Fill(Genconverter().IPhi2Phi(j));
                        trigJetIPhi->Fill(j);


                    }
                    else
                    {
                        HFarrayPos[i][j].pass = false;
                    }
                    HFPosJets->Fill(PosJets [i][j]);

                    HFPosRatio->Fill(EtArrayPos [i][j]/PosJets [i][j]);
                }
                else
                {
                    HFarrayPos[i][j].pass = false;
                }

            }
        }
    }
    if(!isData)
    {
        for(unsigned int n = 0; n < pIn->size(); ++n)  //checks to see if any truth jets are in HF and fills genJets vectors
        { 
            if(pIn->at(n).eta() < -2.964 && pIn->at(n).eta() > -4.716 && pIn->at(n).et() > 10)
            {
                genJet aGenJet = genJet();
                aGenJet.Et = pIn->at(n).et();
                aGenJet.eta = pIn->at(n).eta();
                aGenJet.phi = pIn->at(n).phi();
                aGenJet.matchPass = false;
                aGenJet.id = n;
                negGenJets.push_back(aGenJet);
            }
            if(pIn->at(n).eta() > 2.964 && pIn->at(n).eta() < 4.716 && pIn->at(n).et() > 10)
            {
                genJet aGenJet = genJet();
                aGenJet.Et = pIn->at(n).et();
                aGenJet.eta = pIn->at(n).eta();
                aGenJet.phi = pIn->at(n).phi();
                aGenJet.matchPass = false;
                aGenJet.id = n;
                posGenJets.push_back(aGenJet);
            }

        }



        for(int l = 31; l < 39; l++) //loops over ieta
        {
            for(int m = 1; m < 72; m+=2) //loops over iphi
            {
                if(HFarrayPos[l][m].pass)  //checks for a good reco jet which passes criteria enforced above
                {

                    passed = 0; 
                    dR = 100;
                    for(unsigned int n = 0; n < posGenJets.size(); ++n)
                    {

                        if(deltaR(posGenJets.at(n).eta, posGenJets.at(n).phi, Genconverter().IEta2Eta(l), Genconverter().IPhi2Phi(m)) < dR)
                        {
                            dR = deltaR(posGenJets.at(n).eta, posGenJets.at(n).phi, Genconverter().IEta2Eta(l), Genconverter().IPhi2Phi(m));
                            passed = n;
                        }
                    }
                    //    cout << "Best Match Pos Jet:";
                    //    cout << dR << endl;
                    if(dR < 1)
                    {
                        posGenJets.at(passed).matchPass = true;   //says that the gen jet at 'passed' is close enough to a reco jet
                    }

                }
                if(HFarrayNeg[l][m].pass)
                {



                    dR = 100;
                    passed = 0;
                    for(unsigned int n = 0; n < negGenJets.size(); ++n)
                    {

                        if(deltaR(negGenJets.at(n).eta, negGenJets.at(n).phi, Genconverter().IEta2Eta(-l), Genconverter().IPhi2Phi(m)) < dR)
                        {
                            dR = deltaR(negGenJets.at(n).eta, negGenJets.at(n).phi, Genconverter().IEta2Eta(-l), Genconverter().IPhi2Phi(m));
                            passed = n;
                        }
                    }
                    cout << "Best Match Neg Jet:";
                    cout << dR << endl;
                    if(dR < .5)
                    {
                        negGenJets.at(passed).matchPass = true;
                    }
                }
            }
        }



        for(unsigned int i = 0; i < negGenJets.size(); ++i)  //loops over gen jets again
        {
            if(negGenJets.at(i).matchPass)  //checks for a gen jet which was a close delta R match
            {
                if(negGenJets.at(i).eta < -3.4 && negGenJets.at(i).eta > -4.4) //Apply some cuts on eta, eta, phi, etc before making filling the histo.  For now this justs keeps things deep within acceptance in eta. (to turn it off just set it all to true)
                {
                    EtMatched->Fill(negGenJets.at(i).Et);           
                    EtaMatched->Fill(negGenJets.at(i).eta);
                    PhiMatched->Fill(negGenJets.at(i).phi);
                    absEtaMatched->Fill(fabs(negGenJets.at(i).eta));
                    IPhiMatched->Fill(Genconverter().Phi2Iphi(negGenJets.at(i).phi, negGenJets.at(i).eta));
                    //      cout << "Hi! EtaMatched should be getting:";
                    //      cout << negGenJets.at(i).eta << endl;
                } 
            }
        }    
        for(unsigned int i = 0; i < posGenJets.size(); ++i)
        {
            if(posGenJets.at(i).matchPass)
            {
                if(posGenJets.at(i).eta < 4.4 && posGenJets.at(i).eta > 3.4)
                {
                    EtMatched->Fill(posGenJets.at(i).Et);
                    EtaMatched->Fill(posGenJets.at(i).eta);
                    PhiMatched->Fill(posGenJets.at(i).phi);
                    absEtaMatched->Fill(fabs(posGenJets.at(i).eta));
                    IPhiMatched->Fill(Genconverter().Phi2Iphi(posGenJets.at(i).phi, posGenJets.at(i).eta));
                    //        cout << "Hi!" << endl;
                }
            }
        }
    }
    cout << "Number of negative and positive trigger jets so far:" << endl;
    cout << trigNegJets << endl;
    cout << trigPosJets << endl;
    //         cout << genJetsNegHF << endl;
    //         cout << genJetsPosHF << endl;

    //A silly bit of code to check the deltaR algorithm
    /*  int cand1 [16] = {1, 1, 1, 1, 2, 2, 2, 2, 72, 72, 72, 72, 71, 71, 71, 71};
        int cand2 [16] = {1, 2, 72, 71, 1, 2, 72, 71, 1, 2, 72, 71, 1, 2, 72, 71};
        double delR [16];
        double etaVal = 3.6;
        for(unsigned int cand = 0; cand < 16; ++cand)
        {
        delR [cand] = deltaR(etaVal, Genconverter().IPhi2Phi(cand1[cand]), etaVal, Genconverter().IPhi2Phi(cand2[cand]));
        std::cout << delR[cand] << std::endl;
        }
        */
}
void
JetAlgorithm::beginJob() { }

// ------------ method called once each job just after ending the event loop  ------------

void
JetAlgorithm::endJob() { }

DEFINE_FWK_MODULE(JetAlgorithm);
#endif
