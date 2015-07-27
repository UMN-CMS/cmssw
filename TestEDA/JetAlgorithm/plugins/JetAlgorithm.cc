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
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
//deltaR matching
#include "DataFormats/Math/interface/deltaR.h"
//C++
#include <math.h>
#include <stdlib.h>
//Zach's Converter
#include "TestEDA/JetAlgorithm/interface/Genconverter.h"
//L1 Jets
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"

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
            float seedTotalR;
            std::vector<unsigned int> oldGenIds;   //not implemented
            bool isLoser;                          //not implemented 
            bool isFlat;                           //whether or not there are any equal height peaks inside the jet size
            bool coPeaks [3][3];                   //where the co-peaks are (coPeak[1][1] is the jet itself)
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
            vector<const reco::GenParticle*> genConstituents;
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
        TH2D* HFEtPoswGen;
        TH2D* HFEtNegwGen;
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
        TH1D* EtMatchedwHEcoin;
        TH1D* EtaMatchedwHEcoin;
        TH1D* PhiMatchedwHEcoin;
        TH1D* IPhiMatchedwHEcoin;
        TH1D* absEtaMatchedwHEcoin;
        TH1D* EtMatchedwHFcoin;
        TH1D* EtaMatchedwHFcoin;
        TH1D* PhiMatchedwHFcoin;
        TH1D* IPhiMatchedwHFcoin;
        TH1D* absEtaMatchedwHFcoin;
        TH1D* EtMatched;
        TH1D* EtaMatched;
        TH1D* PhiMatched;
        TH1D* IPhiMatched;
        TH1D* absEtaMatched;
        TH1D* trigJetEt;
        TH1D* trigJetEta;
        TH1D* trigJetPhi;
        TH1D* trigJetEtwHEcoin;
        TH1D* trigJetEtawHEcoin;
        TH1D* trigJetPhiwHEcoin;
        TH1D* trigJetEtwHFcoin;
        TH1D* trigJetEtawHFcoin;
        TH1D* trigJetPhiwHFcoin;
        TH1D* trigJetIPhi;
        TH1D* numTrigJets;
        TH1D* EtPUHist;
        TH1D* L1ForJetPt;
        TH1D* trigJetSeedTotalR;
        TH1D* matchedTrigJetEt;
        TH2D* matchedGenTrigEtR;
        TH1D* trigJetSeedTotalRwHEcoin;
        TH1D* matchedTrigJetEtwHEcoin;
        TH2D* matchedGenTrigEtRwHEcoin;
        TH1D* trigJetSeedTotalRwHFcoin;
        TH1D* matchedTrigJetEtwHFcoin;
        TH2D* matchedGenTrigEtRwHFcoin;
        TH2D* genHadJetEnergy;
        TH2D* genEtDist;
        TH2D* genEtDistIntegral;
        //    TH1D* EtPUHistwJets;
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
    EtMatchedwHEcoin = fs->make<TH1D>("EtMatchedwHEcoin", "Matched genJet Et with opposite eta sign HE coincidence", 100, 0.0, 100.0);
    EtaMatchedwHEcoin = fs->make<TH1D>("EtaMatchedwHEcoin", "Matched genJet Eta with opposite eta sign HE coincidence", 50, -5.0, 5.0);
    PhiMatchedwHEcoin = fs->make<TH1D>("PhiMatchedwHEcoin", "Matched genJet Phi with opposite eta sign HE coincidence", 100, -5.0, 5.0);
    IPhiMatchedwHEcoin = fs->make<TH1D>("IPhiMatchedwHEcoin", "Matched genJet IPhi with opposite eta sign HE coincidence", 76, -1.0, 74);
    absEtaMatchedwHEcoin = fs->make<TH1D>("absEtaMatchedwHEcoin", "Matched genJet |eta| with opposite eta sign HE coincidence", 100, 0, 5.0);
    EtMatchedwHFcoin = fs->make<TH1D>("EtMatchedwHEcoin", "Matched genJet Et with opposite eta sign HE coincidence", 100, 0.0, 100.0);
    EtaMatchedwHFcoin = fs->make<TH1D>("EtaMatchedwHEcoin", "Matched genJet Eta with opposite eta sign HE coincidence", 50, -5.0, 5.0);
    PhiMatchedwHFcoin = fs->make<TH1D>("PhiMatchedwHEcoin", "Matched genJet Phi with opposite eta sign HE coincidence", 100, -5.0, 5.0);
    IPhiMatchedwHFcoin = fs->make<TH1D>("IPhiMatchedwHEcoin", "Matched genJet IPhi with opposite eta sign HE coincidence", 76, -1.0, 74);
    absEtaMatchedwHFcoin = fs->make<TH1D>("absEtaMatchedwHEcoin", "Matched genJet |eta| with opposite eta sign HE coincidence", 100, 0, 5.0);
    trigJetEt = fs->make<TH1D>("trigJetEt", "trigger Jet Et", 100, 0.0, 100.0);
    trigJetEta = fs->make<TH1D>("trigJetEta", "trigger Jet Eta", 50, -5.0, 5.0);
    trigJetPhi = fs->make<TH1D>("trigJetPhi", "trigger Jet Phi", 100, -5.0, 5.0);
    trigJetIPhi = fs->make<TH1D>("trigJetIPhi", "trigger Jet IPhi", 76, -1.0, 74);
    trigJetEtwHEcoin = fs->make<TH1D>("trigJetEtwHEcoin", "trigger Jet Et with opposite eta sign HE coincidence", 100, 0.0, 100.0);
    trigJetEtawHEcoin= fs->make<TH1D>("trigJetEtawHEcoin", "trigger Jet Eta with opposite eta sign HE coincidence", 50, -5.0, 5.0);
    trigJetPhiwHEcoin = fs->make<TH1D>("trigJetPhiwHEcoin", "trigger Jet Phi with opposite eta sign HE coincidence", 100, -5.0, 5.0);
    trigJetEtwHFcoin = fs->make<TH1D>("trigJetEtwHEcoin", "trigger Jet Et with opposite eta sign HE coincidence", 100, 0.0, 100.0);
    trigJetEtawHFcoin= fs->make<TH1D>("trigJetEtawHEcoin", "trigger Jet Eta with opposite eta sign HE coincidence", 50, -5.0, 5.0);
    trigJetPhiwHFcoin = fs->make<TH1D>("trigJetPhiwHEcoin", "trigger Jet Phi with opposite eta sign HE coincidence", 100, -5.0, 5.0);
    numTrigJets = fs->make<TH1D>("numTrigJets", "Number of Trigger Jets In Event", 10, 0, 10);
    EtPUHist = fs->make<TH1D>("EtPU", "Average Et Outside of Trigger Jets", 80, 0, 20);
    L1ForJetPt = fs->make<TH1D>("L1ForJetPt","Forward L1 Jet Pt", 20, 0.0, 100.0);
    trigJetSeedTotalR = fs->make<TH1D>("trigJetSeedTotalR","trigger Jet Seed to Total Ratio",10, 0.0, 1.0);
    matchedTrigJetEt = fs->make<TH1D>("matchedTrigJetEt","Et of matched trigger jets",100,0.0,100.0);
    matchedGenTrigEtR = fs->make<TH2D>("matchedGenTrigEtR","Ratio of matched Gen and Trig Jet Et",100,0.0,10.0,100,0.0,100.0);
    trigJetSeedTotalRwHEcoin = fs->make<TH1D>("trigJetSeedTotalRwHEcoin","trigger Jet Seed to Total Ratio with opposite eta sign HE coincidence",10, 0.0, 1.0);
    matchedTrigJetEtwHEcoin = fs->make<TH1D>("matchedTrigJetEtwHEcoin","Et of matched trigger jets with opposite eta sign HE coincidence",100,0.0,100.0);
    matchedGenTrigEtRwHEcoin = fs->make<TH2D>("matchedGenTrigEtRwHEcoin","Ratio of matched Gen and Trig Jet Et with opposite eta sign HE coincidence",100,0.0,10.0,100,0.0,100.0);
    trigJetSeedTotalRwHFcoin = fs->make<TH1D>("trigJetSeedTotalRwHEcoin","trigger Jet Seed to Total Ratio with opposite eta sign HE coincidence",10, 0.0, 1.0);
    matchedTrigJetEtwHFcoin = fs->make<TH1D>("matchedTrigJetEtwHEcoin","Et of matched trigger jets with opposite eta sign HE coincidence",100,0.0,100.0);
    matchedGenTrigEtRwHFcoin = fs->make<TH2D>("matchedGenTrigEtRwHEcoin","Ratio of matched Gen and Trig Jet Et with opposite eta sign HE coincidence",100,0.0,10.0,100,0.0,100.0);
    genHadJetEnergy = fs->make<TH2D>("genEMHadJetEnergy","GenJet EM and Hadronic Energy",200,0.0,2.0,100,0.0,100.0);
    genEtDist = fs->make<TH2D>("genEtDist","Distribution of gen jet constituent energy ratio",10,0.0,1.0,20,0.0,2.0);
    genEtDistIntegral = fs->make<TH2D>("genEtDistIntegral","Integral Distribution of gen jet constituent energy ratio",10,0.0,1.0,10,0.0,1.0);
    //  EtPUHistwJets = fs->make<TH1D>("EtPUwJets", "Average Et Outside of Trigger Jets in Events with a Trigger Jet", 80, 0, 20);
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
    // Parameters to set or unset
    double seedThreshold = 5;
    double jetThreshold = 10;
    double genJetThreshold = 10;   
    double cenJetThreshold = 16;   //endcap jet minimum et for coincidence
    double cenGenJetThresh = 10;   //endcap jet minimum em+had et for coincidence in MC
    bool pileUpSubtracting = true; //true means pile up Et subtraction will be done 
    bool splitPU = true;           //true means pile up is averaged in negative and postive HF separately (still goes in the same histo)
    bool coPeaks = true;          //true means the algorithm includes coPeaks and attempts to resolve degeneracy by highest jet Et
    bool useCoin = true;          //records trig jets in coincidence with opposite eta sign HE l1 jet
    double endCapEtaMin = 0.00;   //eta ranges for the cen jet coincidence (roughly corresponds to endcap)
    double endCapEtaMax = 2.44;   //

    double totalPosHFEt = 0;
    double totalNegHFEt = 0;
    double totalPosTrigJetEt = 0;
    double totalNegTrigJetEt = 0;
    double EtPUNeg = 0;
    double EtPUPos = 0;

    int passed;
    int left;
    int right;
    int eventPosTrigJets = 0;
    int eventNegTrigJets = 0;

    bool isPosJet = false;
    bool isNegJet = false;
    bool isData = iEvent.isRealData();
    bool posHEjet = false;
    bool negHEjet = false;
    bool negHFjet = false;
    bool posHFjet = false;
    if(posHFjet && negHFjet) cout << "hi";
    vector<genJet> negGenJets;
    vector<genJet> posGenJets;

    Jet myJet = Jet();
    myJet.pass = false;
    Jet HFarrayPos [40][72];
    Jet HFarrayNeg [40][72];
    genJet myGenJet = genJet();
    myGenJet.matchPass = false;
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
    edm::Handle<HcalTrigPrimDigiCollection> hfpr_digi;
    iEvent.getByLabel("simHcalTriggerPrimitiveDigis", hfpr_digi);
//////////////////////////////////////////LOOPS OVER OLD L1 JETS//////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(isData)
    {
        cout << "Getting l1jets!" << endl;
        Handle<std::vector<l1extra::L1JetParticle>> l1Forjets;
        iEvent.getByLabel(edm::InputTag("l1extraParticles","Forward"), l1Forjets);
        cout << "Looping over them" << endl;
        for(unsigned int iJet = 0; iJet < l1Forjets->size(); ++iJet)
        {
            if(fabs(l1Forjets->at(iJet).eta()) >= 3.314 && fabs(l1Forjets->at(iJet).eta()) <= 4.539)
            {
                cout << "Filling Jet: ";
                cout << iJet << endl;
                L1ForJetPt->Fill(l1Forjets->at(iJet).pt());
                if(useCoin)
                {
                    if(l1Forjets->at(iJet).et() >= cenJetThreshold && l1Forjets->at(iJet).eta()/fabs(l1Forjets->at(iJet).eta()) < 0) negHFjet = true;
                    if(l1Forjets->at(iJet).et() >= cenJetThreshold && l1Forjets->at(iJet).eta()/fabs(l1Forjets->at(iJet).eta()) > 0) posHFjet = true;
                }
            }

        }
        if(useCoin)
        {
            Handle<std::vector<l1extra::L1JetParticle>> l1Cenjets;
            iEvent.getByLabel(edm::InputTag("l1extraParticles","Central"), l1Cenjets);
            for(unsigned int aJet = 0; aJet < l1Cenjets->size(); ++aJet)
            {
                if(l1Cenjets->at(aJet).eta() <= -1*endCapEtaMin && l1Cenjets->at(aJet).eta() >= -1*endCapEtaMax && l1Cenjets->at(aJet).et() >= cenJetThreshold) negHEjet = true;
                if(l1Cenjets->at(aJet).eta() <= endCapEtaMax && l1Cenjets->at(aJet).eta() >= endCapEtaMin && l1Cenjets->at(aJet).et() >= cenJetThreshold) posHEjet = true;
            }

        }
    }
////////////////////////////////////////////////LOOPS OVER GEN JETS AND MAKES PLOTS/////////////////////////////////////////////////////////////////////////////////////    
    if(!isData)
    {
        for(unsigned int k = 0; k < pIn->size(); ++k)  //loops over all ak5 jest and checks if they are in HF
        {
            if(useCoin && 
                    ((pIn->at(k).hadEnergy() + pIn->at(k).emEnergy())/pIn->at(k).energy())*pIn->at(k).et() > cenGenJetThresh &&
                    pIn->at(k).eta() < -1*endCapEtaMin && 
                    pIn->at(k).eta() > -1*endCapEtaMax)
                negHEjet = true;
            if(useCoin && 
                    ((pIn->at(k).hadEnergy() + pIn->at(k).emEnergy())/pIn->at(k).energy())*pIn->at(k).et() > cenGenJetThresh &&
                    pIn->at(k).eta() < endCapEtaMax && 
                    pIn->at(k).eta() > endCapEtaMin)
                posHEjet = true;
            if(((pIn->at(k).hadEnergy() + pIn->at(k).emEnergy())/pIn->at(k).energy())*pIn->at(k).et() > genJetThreshold)
                //makes sure that each jet has at least roughly 10GeV of transverse hadronic and em energy
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
                        genHadJetEnergy->Fill((pIn->at(k).hadEnergy() + pIn->at(k).emEnergy())/pIn->at(k).energy(),pIn->at(k).et());
                    }
                    IPhi->Fill(Genconverter().Phi2Iphi(pIn->at(k).phi(), pIn->at(k).eta()));
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
                        genHadJetEnergy->Fill((pIn->at(k).hadEnergy() + pIn->at(k).emEnergy())/pIn->at(k).energy(),pIn->at(k).et());
                    }
                    IPhi->Fill(Genconverter().Phi2Iphi(pIn->at(k).phi(), pIn->at(k).eta()));
                }
            }
        }
    }
//////////////////////////////////////////////////////////////ET ARRAY OF EVENTS MADE//////////////////////////////////////////////////////////////////////////////////////////
    //Now the code loops over cell by cell in the event to begin to reconstruct jets.  It only looks at events that have an ak5 genJet in HF.
    if(isPosJet || isNegJet || isData) //checks that there is an HF ak5 truth jet
    {
        for(HcalTrigPrimDigiCollection::const_iterator tp = hfpr_digi->begin(); tp != hfpr_digi->end(); ++tp)
        {
            if(abs(tp->id().ieta()) < 30 || abs(tp->id().ieta()) > 39) continue; //Doesn't use the part of HF which has a different cell size or the other parts of the detector
            //We need the sample of interest (SOI)
            //Fills out the Et distribution

            if(tp->id().ieta() < 0 && (isNegJet || isData))
            {
                EtArrayNeg [abs(tp->id().ieta())][tp->id().iphi()] = outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt());
                EtArrayNeg [abs(tp->id().ieta())][tp->id().iphi()+1] = outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt());
            }
            else if (tp->id().ieta() > 0 && (isPosJet || isData))
            {
                EtArrayPos [tp->id().ieta()][tp->id().iphi()] = outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt());
                EtArrayPos [tp->id().ieta()][tp->id().iphi()+1] = outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt());
            }
        }
    }

//////////////////////////////////////////////MAKES TRIGGER JETS////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    for(int i = 31; i < 40; i++) //loops over ieta
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
            totalNegHFEt += EtArrayNeg [i][j];
            totalPosHFEt += EtArrayPos [i][j];
            if(EtArrayNeg [i][j] > seedThreshold && (isNegJet || isData))  //Seed cell threshold
            {
                NegJets [i][j] = EtArrayNeg [i-1][left] + EtArrayNeg [i-1][j] + EtArrayNeg [i-1][right] +
                    EtArrayNeg [i][left]   + EtArrayNeg [i][j]   + EtArrayNeg [i][right]   +
                    EtArrayNeg [i+1][left] + EtArrayNeg [i+1][j] + EtArrayNeg [i+1][right];

                //checks that the seed is greater than the periphery
                if(EtArrayNeg [i][j] >= EtArrayNeg [i-1][left] && EtArrayNeg [i][j] >= EtArrayNeg [i-1][j]   && EtArrayNeg [i][j] >= EtArrayNeg [i-1][right] &&
                        EtArrayNeg [i][j] >= EtArrayNeg [i][left]                                                && EtArrayNeg [i][j] >= EtArrayNeg [i][right]   &&
                        EtArrayNeg [i][j] >= EtArrayNeg [i+1][left] && EtArrayNeg [i][j] >= EtArrayNeg [i+1][j]   && EtArrayNeg [i][j] >= EtArrayNeg [i+1][right] 
                        && NegJets [i][j] > jetThreshold) //total energy threshold
                {
                    cout << "Made Trigger Jet in NegHF" << endl;
                    trigNegJets++;
                    eventNegTrigJets++;
                    HFarrayNeg[i][j].jetEt = NegJets [i][j];
                    totalNegTrigJetEt += NegJets[i][j];
                    HFarrayNeg[i][j].seedEt = EtArrayNeg [i][j];
                    HFarrayNeg[i][j].seedTotalR = EtArrayNeg[i][j]/NegJets[i][j];
                    HFarrayNeg[i][j].pass = true; //initial true setting for whether it is a good jet
                    HFarrayNeg[i][j].isFlat = false;  //initial assumption
                    if(HFarrayNeg[i][j].jetEt < 0)
                    {
                        cout << "hi";
                    }
                    if(EtArrayNeg [i][j] == EtArrayNeg [i][left])
                    {
                        HFarrayNeg[i][j].coPeaks[1][0] = true;
                        HFarrayNeg[i][j].isFlat = true;
                    }
                    if(EtArrayNeg [i][j] == EtArrayNeg [i][right])
                    {
                        HFarrayNeg[i][j].coPeaks[1][2] = true;
                        HFarrayNeg[i][j].isFlat = true;
                    }
                    if(EtArrayNeg [i][j] == EtArrayNeg [i-1][left])
                    {
                        HFarrayNeg[i][j].coPeaks[0][0] = true;
                        HFarrayNeg[i][j].isFlat = true;
                    }
                    if(EtArrayNeg [i][j] == EtArrayNeg [i-1][j])
                    {
                        HFarrayNeg[i][j].coPeaks[0][1] = true;
                        HFarrayNeg[i][j].isFlat = true;
                    }
                    if(EtArrayNeg [i][j] == EtArrayNeg [i-1][right])
                    {
                        HFarrayNeg[i][j].coPeaks[0][2] = true;
                        HFarrayNeg[i][j].isFlat = true;
                    }
                    if(EtArrayNeg [i][j] == EtArrayNeg [i+1][left])
                    {
                        HFarrayNeg[i][j].coPeaks[2][0] = true;
                        HFarrayNeg[i][j].isFlat = true;
                    }
                    if(EtArrayNeg [i][j] == EtArrayNeg [i+1][j])
                    {
                        HFarrayNeg[i][j].coPeaks[2][1] = true;
                        HFarrayNeg[i][j].isFlat = true;
                    }
                    if(EtArrayNeg [i][j] == EtArrayNeg [i+1][right])
                    {
                        HFarrayNeg[i][j].coPeaks[2][2] = true;
                        HFarrayNeg[i][j].isFlat = true;
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
            if(EtArrayPos [i][j] > seedThreshold && (isPosJet || isData))   //Seed cell threshold
            {
                PosJets [i][j] = EtArrayPos [i-1][left] + EtArrayPos [i-1][j] + EtArrayPos [i-1][right] +
                    EtArrayPos [i][left]   + EtArrayPos [i][j]   + EtArrayPos [i][right]   +
                    EtArrayPos [i+1][left] + EtArrayPos [i+1][j] + EtArrayPos [i+1][right];

                //checks that the seed is greater than the periphery
                if(EtArrayPos [i][j] >= EtArrayPos [i-1][left] && EtArrayPos [i][j] >= EtArrayPos [i-1][j]   && EtArrayPos [i][j] >= EtArrayPos [i-1][right] &&
                        EtArrayPos [i][j] >= EtArrayPos [i][left]                                                && EtArrayPos [i][j] >= EtArrayPos [i][right]   &&
                        EtArrayPos [i][j] >= EtArrayPos [i+1][left] && EtArrayPos [i][j] >= EtArrayPos [i+1][j]   && EtArrayPos [i][j] >= EtArrayPos [i+1][right]
                        && PosJets[i][j] > jetThreshold) //total energy threshold
                {
                    cout << "Made Trigger Jet in PosHF" << endl;
                    trigPosJets++;
                    eventPosTrigJets++;
                    HFarrayPos[i][j].jetEt = PosJets [i][j];
                    totalPosTrigJetEt += EtArrayPos [i][j];
                    HFarrayPos[i][j].seedEt = EtArrayPos [i][j];
                    HFarrayPos[i][j].seedTotalR = EtArrayPos[i][j]/PosJets[i][j];
                    HFarrayPos[i][j].pass = true;                       //initial true setting for whether it is a good jet
                    HFarrayPos[i][j].isFlat = false;
                    if(EtArrayPos [i][j] == EtArrayPos [i][left])
                    {
                        HFarrayPos[i][j].coPeaks[1][0] = true;
                        HFarrayPos[i][j].isFlat = true;
                    }
                    if(EtArrayPos [i][j] == EtArrayPos [i][right])
                    {
                        HFarrayPos[i][j].coPeaks[1][2] = true;
                        HFarrayPos[i][j].isFlat = true;
                    }
                    if(EtArrayPos [i][j] == EtArrayPos [i-1][left])
                    {
                        HFarrayPos[i][j].coPeaks[0][0] = true;
                        HFarrayPos[i][j].isFlat = true;
                    }
                    if(EtArrayPos [i][j] == EtArrayPos [i-1][j])
                    {
                        HFarrayPos[i][j].coPeaks[0][1] = true;
                        HFarrayPos[i][j].isFlat = true;
                    }
                    if(EtArrayPos [i][j] == EtArrayPos [i-1][right])
                    {
                        HFarrayPos[i][j].coPeaks[0][2] = true;
                        HFarrayPos[i][j].isFlat = true;
                    }
                    if(EtArrayPos [i][j] == EtArrayPos [i+1][left])
                    {
                        HFarrayPos[i][j].coPeaks[2][0] = true;
                        HFarrayPos[i][j].isFlat = true;
                    }
                    if(EtArrayPos [i][j] == EtArrayPos [i+1][j])
                    {
                        HFarrayPos[i][j].coPeaks[2][1] = true;
                        HFarrayPos[i][j].isFlat = true;
                    }
                    if(EtArrayPos [i][j] == EtArrayPos [i+1][right])
                    {
                        HFarrayPos[i][j].coPeaks[2][2] = true;
                        HFarrayPos[i][j].isFlat = true;
                    }

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
////////////////////////////////////////////////////////RESOLVE FLAT TOPPEDNESS////////////////////////////////////////////////////////////////////////////////
    //Now it is time to resolve overlapping and or flat topped jets
    for(int iphi = 1; iphi < 72; iphi+=2)
    {
        for(int ieta = 31; ieta < 40; ieta++)
        {
            //check and resolve flat tops
            if(HFarrayPos[ieta][iphi].isFlat && HFarrayPos[ieta][iphi].pass)
            {
                for(int relphi = 0; relphi < 3; relphi++)
                {
                    for(int releta = 0; releta < 3; releta++)
                    {
                        if(HFarrayPos[ieta][iphi].coPeaks[releta][relphi]  && coPeaks)
                        {
                            int coPhi = iphi+2*relphi-2;
                            if(coPhi > 71)
                            {
                                coPhi = 1;
                            }
                            if(coPhi < 1)
                            {
                                coPhi = 71;
                            }
                            if(HFarrayPos[ieta][iphi].jetEt < HFarrayPos[ieta+releta-1][coPhi].jetEt)
                            {
                                HFarrayPos[ieta][iphi].pass = false;
                                //   HFarrayPos[ieta+releta-1][coPhi].pass = true; //unclear whether this is a good idea
                            }
                        }
                        if(HFarrayPos[ieta][iphi].isFlat  && !coPeaks)
                        {
                            HFarrayPos[ieta][iphi].pass = false;
                        }
                    }
                }
            }
            if(HFarrayNeg[ieta][iphi].isFlat && HFarrayNeg[ieta][iphi].pass)
            {
                for(int relphi = 0; relphi < 3; relphi++)
                {
                    for(int releta = 0; releta < 3; releta++)
                    {
                        if(HFarrayNeg[ieta][iphi].coPeaks[releta][relphi] && coPeaks)
                        {
                            int coPhi = iphi+2*relphi-2;
                            if(coPhi > 71)
                            {
                                coPhi = 1;
                            }
                            if(coPhi < 1)
                            {
                                coPhi = 71;
                            }
                            if(HFarrayNeg[ieta][iphi].jetEt < HFarrayNeg[ieta+releta-1][coPhi].jetEt)
                            {
                                HFarrayNeg[ieta][iphi].pass = false;
                                //  HFarrayPos[ieta+releta-1][coPhi].pass = true;  //unclear whether this is a good idea
                            }
                        }
                        if(HFarrayNeg[ieta][iphi].isFlat  && !coPeaks)
                        {
                            HFarrayNeg[ieta][iphi].pass = false;
                        }

                    }
                }
            }

        }
    }
/////////////////////////////////////////////PILE-UP SUBTRACTION//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
////Now, it loops over every cell to find the total energy in HF, then subtracts the newly made jet energies.  Then averages over the number of cells in HF (minus the ones taken up by a jet).  It then subtracts that average energy (multiplied by 9) from each trigger jet and evaluates whether or not the trigger jet is still valid (passes seed threshold and total threshold) if it no longer passes it throws it in the average and loops back through (and adjusts the number of jets in the event accordingly)
    while(pileUpSubtracting)
    {
        //Initial average Et per bin
        if(splitPU)
        {
            EtPUNeg = (totalNegHFEt - totalNegTrigJetEt)/(162-(eventNegTrigJets)*9);
            EtPUPos = (totalPosHFEt - totalPosTrigJetEt)/(162-(eventPosTrigJets)*9);
        }
        else
        {
            EtPUNeg = (totalNegHFEt + totalPosHFEt - totalPosTrigJetEt - totalNegTrigJetEt)/(324-(eventPosTrigJets + eventNegTrigJets)*9);
            EtPUPos = EtPUNeg;
        }
        pileUpSubtracting = false;
        for(int j = 31; j < 40; j++)
        {
            for(int k = 1; k < 72; k+=2)
            {
                if(HFarrayPos[j][k].pass)
                {
                    if(((HFarrayPos[j][k].seedEt - EtPUPos) < seedThreshold) || ((HFarrayPos[j][k].jetEt - 9*EtPUPos) < jetThreshold))
                    {
                        HFarrayPos[j][k].pass = false;
                        eventPosTrigJets--;
                        trigPosJets--;
                        totalPosTrigJetEt -= HFarrayPos[j][k].jetEt;
                        totalPosHFEt += HFarrayPos[j][k].jetEt;
                        //                  pileUpSubtracting = true;
                    }
                    else
                    {
                        HFarrayPos[j][k].jetEt -= 9*EtPUPos;
                        HFarrayPos[j][k].seedEt -= EtPUPos;
                    }

                }
                if(HFarrayNeg[j][k].pass)
                {
                    if(((HFarrayNeg[j][k].seedEt - EtPUNeg) < seedThreshold) || ((HFarrayNeg[j][k].jetEt - 9*EtPUNeg) < jetThreshold))
                    {
                        HFarrayNeg[j][k].pass = false;
                        eventNegTrigJets--;
                        trigNegJets--;
                        totalNegTrigJetEt -= HFarrayNeg[j][k].jetEt;
                        totalNegHFEt += HFarrayNeg[j][k].jetEt;
                        //                    pileUpSubtracting = true;  //Setting pileUpSubtracting to true will make the algorithm recalculate after a jet is removed.
                    }
                    else
                    {
                        HFarrayNeg[j][k].jetEt -= 9*EtPUNeg;
                        HFarrayNeg[j][k].seedEt -= EtPUNeg;
                    }

                }

            }
        }
    }
    EtPUHist->Fill(EtPUNeg);
    EtPUHist->Fill(EtPUPos);
    numTrigJets->Fill(eventPosTrigJets+eventNegTrigJets);
/////////////////////////////////////////MAKES MY OWN CUSTOM GEN JET VECTORS/////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
                aGenJet.genConstituents = pIn->at(n).getGenConstituents();
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
                aGenJet.genConstituents = pIn->at(n).getGenConstituents();
                posGenJets.push_back(aGenJet);
            }
        }

////////////////////////////////////////////////////////A VERY LONG WAY TO MAKE A SINGLE PLOT///////////////////////////////////////////////////////////////////////////////////////////////
        double energyWithinOne;
        double energyWithinTwo;
        double energyWithinThree;
        double energyWithinFour;
        double energyWithinFive;
        double energyWithinSix;
        double energyWithinSeven;
        double energyWithinEight;
        double energyWithinNine;
        double energyWithinTen;
        for(unsigned int ijet = 0; ijet < posGenJets.size(); ++ijet)
        {
            energyWithinOne = 0;
            energyWithinTwo = 0;
            energyWithinThree = 0;
            energyWithinFour = 0;
            energyWithinFive = 0;
            energyWithinSix = 0;
            energyWithinSeven = 0;
            energyWithinEight = 0;
            energyWithinNine = 0;
            energyWithinTen = 0;
            for(unsigned int iconst = 0; iconst < posGenJets.at(ijet).genConstituents.size(); iconst++)
            {
                if(posGenJets.at(ijet).genConstituents.at(iconst)->status()==1){
                    if(deltaR(posGenJets.at(ijet).genConstituents.at(iconst)->eta(),posGenJets.at(ijet).genConstituents.at(iconst)->phi(),posGenJets.at(ijet).eta,posGenJets.at(ijet).phi) < .1)
                        energyWithinOne += posGenJets.at(ijet).genConstituents.at(iconst)->et();
                    if(deltaR(posGenJets.at(ijet).genConstituents.at(iconst)->eta(),posGenJets.at(ijet).genConstituents.at(iconst)->phi(),posGenJets.at(ijet).eta,posGenJets.at(ijet).phi)< .2)
                        energyWithinTwo += posGenJets.at(ijet).genConstituents.at(iconst)->et();
                    if(deltaR(posGenJets.at(ijet).genConstituents.at(iconst)->eta(),posGenJets.at(ijet).genConstituents.at(iconst)->phi(),posGenJets.at(ijet).eta,posGenJets.at(ijet).phi) < .3)
                        energyWithinThree += posGenJets.at(ijet).genConstituents.at(iconst)->et();
                    if(deltaR(posGenJets.at(ijet).genConstituents.at(iconst)->eta(),posGenJets.at(ijet).genConstituents.at(iconst)->phi(),posGenJets.at(ijet).eta,posGenJets.at(ijet).phi) < .4)
                        energyWithinFour += posGenJets.at(ijet).genConstituents.at(iconst)->et();
                    if(deltaR(posGenJets.at(ijet).genConstituents.at(iconst)->eta(),posGenJets.at(ijet).genConstituents.at(iconst)->phi(),posGenJets.at(ijet).eta,posGenJets.at(ijet).phi) < .5)
                        energyWithinFive += posGenJets.at(ijet).genConstituents.at(iconst)->et();
                    if(deltaR(posGenJets.at(ijet).genConstituents.at(iconst)->eta(),posGenJets.at(ijet).genConstituents.at(iconst)->phi(),posGenJets.at(ijet).eta,posGenJets.at(ijet).phi) < .6)
                        energyWithinSix += posGenJets.at(ijet).genConstituents.at(iconst)->et();
                    if(deltaR(posGenJets.at(ijet).genConstituents.at(iconst)->eta(),posGenJets.at(ijet).genConstituents.at(iconst)->phi(),posGenJets.at(ijet).eta,posGenJets.at(ijet).phi) < .7)
                        energyWithinSeven += posGenJets.at(ijet).genConstituents.at(iconst)->et();
                    if(deltaR(posGenJets.at(ijet).genConstituents.at(iconst)->eta(),posGenJets.at(ijet).genConstituents.at(iconst)->phi(),posGenJets.at(ijet).eta,posGenJets.at(ijet).phi) < .8)
                        energyWithinEight += posGenJets.at(ijet).genConstituents.at(iconst)->et();
                    if(deltaR(posGenJets.at(ijet).genConstituents.at(iconst)->eta(),posGenJets.at(ijet).genConstituents.at(iconst)->phi(),posGenJets.at(ijet).eta,posGenJets.at(ijet).phi) < .9)
                        energyWithinNine += posGenJets.at(ijet).genConstituents.at(iconst)->et();
                    if(deltaR(posGenJets.at(ijet).genConstituents.at(iconst)->eta(),posGenJets.at(ijet).genConstituents.at(iconst)->phi(),posGenJets.at(ijet).eta,posGenJets.at(ijet).phi) < 1.)
                        energyWithinTen += posGenJets.at(ijet).genConstituents.at(iconst)->et();
                }

            }
            genEtDist->Fill(.09,energyWithinOne/posGenJets.at(ijet).Et);
            genEtDist->Fill(.19,energyWithinTwo/posGenJets.at(ijet).Et);
            genEtDist->Fill(.29,energyWithinThree/posGenJets.at(ijet).Et);
            genEtDist->Fill(.39,energyWithinFour/posGenJets.at(ijet).Et);
            genEtDist->Fill(.49,energyWithinFive/posGenJets.at(ijet).Et);
            genEtDist->Fill(.59,energyWithinSix/posGenJets.at(ijet).Et);
            genEtDist->Fill(.69,energyWithinSeven/posGenJets.at(ijet).Et);
            genEtDist->Fill(.79,energyWithinEight/posGenJets.at(ijet).Et);
            genEtDist->Fill(.89,energyWithinNine/posGenJets.at(ijet).Et);
            genEtDist->Fill(.99,energyWithinTen/posGenJets.at(ijet).Et);
        }
        for(unsigned int ijet = 0; ijet < negGenJets.size(); ++ijet)
        {
            energyWithinOne = 0;
            energyWithinTwo = 0;
            energyWithinThree = 0;
            energyWithinFour = 0;
            energyWithinFive = 0;
            energyWithinSix = 0;
            energyWithinSeven = 0;
            energyWithinEight = 0;
            energyWithinNine = 0;
            energyWithinTen = 0;
            for(unsigned int iconst = 0; iconst < negGenJets.at(ijet).genConstituents.size(); iconst++)
            {
                if(negGenJets.at(ijet).genConstituents.at(iconst)->status()==1){
                    if(deltaR(negGenJets.at(ijet).genConstituents.at(iconst)->eta(),negGenJets.at(ijet).genConstituents.at(iconst)->phi(),negGenJets.at(ijet).eta,negGenJets.at(ijet).phi) < .1)
                        energyWithinOne += negGenJets.at(ijet).genConstituents.at(iconst)->et();
                    if(deltaR(negGenJets.at(ijet).genConstituents.at(iconst)->eta(),negGenJets.at(ijet).genConstituents.at(iconst)->phi(),negGenJets.at(ijet).eta,negGenJets.at(ijet).phi) < .2)
                        energyWithinTwo += negGenJets.at(ijet).genConstituents.at(iconst)->et();
                    if(deltaR(negGenJets.at(ijet).genConstituents.at(iconst)->eta(),negGenJets.at(ijet).genConstituents.at(iconst)->phi(),negGenJets.at(ijet).eta,negGenJets.at(ijet).phi) < .3)
                        energyWithinThree += negGenJets.at(ijet).genConstituents.at(iconst)->et();
                    if(deltaR(negGenJets.at(ijet).genConstituents.at(iconst)->eta(),negGenJets.at(ijet).genConstituents.at(iconst)->phi(),negGenJets.at(ijet).eta,negGenJets.at(ijet).phi) < .4)
                        energyWithinFour += negGenJets.at(ijet).genConstituents.at(iconst)->et();
                    if(deltaR(negGenJets.at(ijet).genConstituents.at(iconst)->eta(),negGenJets.at(ijet).genConstituents.at(iconst)->phi(),negGenJets.at(ijet).eta,negGenJets.at(ijet).phi) < .5)
                        energyWithinFive += negGenJets.at(ijet).genConstituents.at(iconst)->et();
                    if(deltaR(negGenJets.at(ijet).genConstituents.at(iconst)->eta(),negGenJets.at(ijet).genConstituents.at(iconst)->phi(),negGenJets.at(ijet).eta,negGenJets.at(ijet).phi) < .6)
                        energyWithinSix += negGenJets.at(ijet).genConstituents.at(iconst)->et();
                    if(deltaR(negGenJets.at(ijet).genConstituents.at(iconst)->eta(),negGenJets.at(ijet).genConstituents.at(iconst)->phi(),negGenJets.at(ijet).eta,negGenJets.at(ijet).phi) < .7)
                        energyWithinSeven += negGenJets.at(ijet).genConstituents.at(iconst)->et();
                    if(deltaR(negGenJets.at(ijet).genConstituents.at(iconst)->eta(),negGenJets.at(ijet).genConstituents.at(iconst)->phi(),negGenJets.at(ijet).eta,negGenJets.at(ijet).phi) < .8)
                        energyWithinEight += negGenJets.at(ijet).genConstituents.at(iconst)->et();
                    if(deltaR(negGenJets.at(ijet).genConstituents.at(iconst)->eta(),negGenJets.at(ijet).genConstituents.at(iconst)->phi(),negGenJets.at(ijet).eta,negGenJets.at(ijet).phi) < .9)
                        energyWithinNine += negGenJets.at(ijet).genConstituents.at(iconst)->et();
                    if(deltaR(negGenJets.at(ijet).genConstituents.at(iconst)->eta(),negGenJets.at(ijet).genConstituents.at(iconst)->phi(),negGenJets.at(ijet).eta,negGenJets.at(ijet).phi) < 1.)
                        energyWithinTen += negGenJets.at(ijet).genConstituents.at(iconst)->et();

                }
            }

            genEtDist->Fill(.09,energyWithinOne/negGenJets.at(ijet).Et);
            genEtDist->Fill(.19,energyWithinTwo/negGenJets.at(ijet).Et);
            genEtDist->Fill(.29,energyWithinThree/negGenJets.at(ijet).Et);
            genEtDist->Fill(.39,energyWithinFour/negGenJets.at(ijet).Et);
            genEtDist->Fill(.49,energyWithinFive/negGenJets.at(ijet).Et);
            genEtDist->Fill(.59,energyWithinSix/negGenJets.at(ijet).Et);
            genEtDist->Fill(.69,energyWithinSeven/negGenJets.at(ijet).Et);
            genEtDist->Fill(.79,energyWithinEight/negGenJets.at(ijet).Et);
            genEtDist->Fill(.89,energyWithinNine/negGenJets.at(ijet).Et);
            genEtDist->Fill(.99,energyWithinTen/negGenJets.at(ijet).Et);
        }

//////////////////////////////////////MATCHING ALGORITHM///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        for(int l = 31; l < 40; l++) //loops over ieta
        {
            for(int m = 1; m < 72; m+=2) //loops over iphi
            {
                if(HFarrayPos[l][m].pass)  //checks for a good trig jet which passes criteria enforced above
                {
                    passed = 0; 
                    dR = 100;
                    HFarrayPos[l][m].match = false;
                    for(unsigned int n = 0; n < posGenJets.size(); ++n)
                    {
                        if(deltaR(posGenJets.at(n).eta, posGenJets.at(n).phi, Genconverter().IEta2Eta(l), Genconverter().IPhi2Phi(m)) < dR)
                        {
                            dR = deltaR(posGenJets.at(n).eta, posGenJets.at(n).phi, Genconverter().IEta2Eta(l), Genconverter().IPhi2Phi(m));
                            passed = n;
                        }
                    }
                    if(dR < .5)
                    {
                        posGenJets.at(passed).matchPass = true;   //says that the gen jet at 'passed' is close enough to a trig jet
                        HFarrayPos[l][m].match = true;
                        HFarrayPos[l][m].genId = passed;
                    }

                }
                if(HFarrayNeg[l][m].pass)
                {
                    dR = 100;
                    passed = 0;
                    HFarrayNeg[l][m].match = false;
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
                        HFarrayNeg[l][m].match = true;
                        HFarrayNeg[l][m].genId = passed;
                    }
                }
            }
        }
///////////////////////////////////////////CODE NOW MAKES LOTS OF PLOTS////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        
        for(unsigned int i = 0; i < negGenJets.size(); ++i)  //loops over gen jets again
        {
            if(negGenJets.at(i).matchPass)  //checks for a gen jet which was a close delta R match
            {
                if(negGenJets.at(i).eta < -3.4 && negGenJets.at(i).eta > -4.4) //Apply some cuts on eta, eta, phi, etc before making filling the histo.  For now this justs keeps things deep within acceptance in eta. (to turn it off just set it all to true)
                {
                    if(posHEjet)
                    {
                        EtMatchedwHEcoin->Fill(negGenJets.at(i).Et);           
                        EtaMatchedwHEcoin->Fill(negGenJets.at(i).eta);
                        PhiMatchedwHEcoin->Fill(negGenJets.at(i).phi);
                        absEtaMatchedwHEcoin->Fill(fabs(negGenJets.at(i).eta));
                        IPhiMatchedwHEcoin->Fill(Genconverter().Phi2Iphi(negGenJets.at(i).phi, negGenJets.at(i).eta));
                    }
                    EtMatched->Fill(negGenJets.at(i).Et);           
                    EtaMatched->Fill(negGenJets.at(i).eta);
                    PhiMatched->Fill(negGenJets.at(i).phi);
                    absEtaMatched->Fill(fabs(negGenJets.at(i).eta));
                    IPhiMatched->Fill(Genconverter().Phi2Iphi(negGenJets.at(i).phi, negGenJets.at(i).eta));
                } 
            }
        }    
        for(unsigned int i = 0; i < posGenJets.size(); ++i)
        {
            if(posGenJets.at(i).matchPass)
            {
                if(posGenJets.at(i).eta < 4.4 && posGenJets.at(i).eta > 3.4)
                {
                    if(negHEjet)
                    {
                        EtMatchedwHEcoin->Fill(posGenJets.at(i).Et);
                        EtaMatchedwHEcoin->Fill(posGenJets.at(i).eta);
                        PhiMatchedwHEcoin->Fill(posGenJets.at(i).phi);
                        absEtaMatchedwHEcoin->Fill(fabs(posGenJets.at(i).eta));
                        IPhiMatchedwHEcoin->Fill(Genconverter().Phi2Iphi(posGenJets.at(i).phi, posGenJets.at(i).eta));
                    }
                    EtMatched->Fill(posGenJets.at(i).Et);
                    EtaMatched->Fill(posGenJets.at(i).eta);
                    PhiMatched->Fill(posGenJets.at(i).phi);
                    absEtaMatched->Fill(fabs(posGenJets.at(i).eta));
                    IPhiMatched->Fill(Genconverter().Phi2Iphi(posGenJets.at(i).phi, posGenJets.at(i).eta));
                }
            }
        }
    }
    for(int q = 31; q < 40; q++)
    {
        for(int p = 1; p < 72; p+=2)
        {
            if(HFarrayNeg[q][p].pass)
            {
                trigJetEt->Fill(HFarrayNeg[q][p].jetEt);
                trigJetEta->Fill(Genconverter().IEta2Eta(-q));
                trigJetPhi->Fill(Genconverter().IPhi2Phi(p));
                trigJetIPhi->Fill(p);
                if(isData)
                {
                    trigJetSeedTotalR->Fill(HFarrayNeg[q][p].seedTotalR);
                }
                if(!isData && HFarrayNeg[q][p].match)
                {
                    trigJetSeedTotalR->Fill(HFarrayNeg[q][p].seedTotalR);
                    matchedTrigJetEt->Fill(HFarrayNeg[q][p].jetEt);
                    matchedGenTrigEtR->Fill(HFarrayNeg[q][p].jetEt/negGenJets.at(HFarrayNeg[q][p].genId).Et,negGenJets.at(HFarrayNeg[q][p].genId).Et);
                }
            }
            if(HFarrayNeg[q][p].pass && posHEjet)
            {
                trigJetEtwHEcoin->Fill(HFarrayNeg[q][p].jetEt);
                trigJetEtawHEcoin->Fill(Genconverter().IEta2Eta(-q));
                trigJetPhiwHEcoin->Fill(Genconverter().IPhi2Phi(p));
                if(isData)
                {
                    trigJetSeedTotalRwHEcoin->Fill(HFarrayNeg[q][p].seedTotalR);
                }
                if(!isData && HFarrayNeg[q][p].match)
                {
                    trigJetSeedTotalRwHEcoin->Fill(HFarrayNeg[q][p].seedTotalR);
                    matchedTrigJetEtwHEcoin->Fill(HFarrayNeg[q][p].jetEt);
                    matchedGenTrigEtRwHEcoin->Fill(HFarrayNeg[q][p].jetEt/negGenJets.at(HFarrayNeg[q][p].genId).Et,negGenJets.at(HFarrayNeg[q][p].genId).Et);
                }
            }
            if(HFarrayPos[q][p].pass && negHEjet)
            {
                trigJetEtwHEcoin->Fill(HFarrayPos[q][p].jetEt);
                trigJetEtawHEcoin->Fill(Genconverter().IEta2Eta(-q));
                trigJetPhiwHEcoin->Fill(Genconverter().IPhi2Phi(p));
                if(isData)
                {
                    trigJetSeedTotalRwHEcoin->Fill(HFarrayPos[q][p].seedTotalR);
                }
                if(!isData && HFarrayPos[q][p].match)
                {
                    trigJetSeedTotalRwHEcoin->Fill(HFarrayPos[q][p].seedTotalR);
                    matchedTrigJetEtwHEcoin->Fill(HFarrayPos[q][p].jetEt);
                    matchedGenTrigEtRwHEcoin->Fill(HFarrayPos[q][p].jetEt/posGenJets.at(HFarrayPos[q][p].genId).Et,posGenJets.at(HFarrayPos[q][p].genId).Et);
                }
            }
            if(HFarrayPos[q][p].pass)
            {
                trigJetEt->Fill(HFarrayPos[q][p].jetEt);
                trigJetEta->Fill(Genconverter().IEta2Eta(q));
                trigJetPhi->Fill(Genconverter().IPhi2Phi(p));
                trigJetIPhi->Fill(p);
                if(isData)
                {
                    trigJetSeedTotalR->Fill(HFarrayPos[q][p].seedTotalR);
                }
                if(!isData && HFarrayPos[q][p].match)
                {
                    trigJetSeedTotalR->Fill(HFarrayPos[q][p].seedTotalR);
                    matchedTrigJetEt->Fill(HFarrayPos[q][p].jetEt);
                    matchedGenTrigEtR->Fill(HFarrayPos[q][p].jetEt/posGenJets.at(HFarrayPos[q][p].genId).Et,posGenJets.at(HFarrayPos[q][p].genId).Et);
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
