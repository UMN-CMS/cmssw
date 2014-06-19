// -*- C++ -*-
//
// Package:    TestEDA/CutAssessment
// Class:      CutAssessment
// 
/**\class CutAssessment CutAssessment.cc TestEDA/CutAssessment/plugins/CutAssessment.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
 */
//
// Original Author:  Zachary Lesko
//         Created:  Wed, 14 May 2014 21:41:06 GMT
//
//


// system include files
#ifndef CutAssessmentClass_included
#define CutAssessmentClass_included 1


#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"




#include "TestEDA/CutAssessment/interface/Genconverter.h"

//
// class declaration
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"

#include <TDirectory.h> 
//obviously trying to get directories to work so I could look at the electrons in the correct regions
#include <TSystem.h>
//maybe
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HcalDigi/interface/HcalTriggerPrimitiveDigi.h"
#include "DataFormats/Common/interface/Handle.h"



#include <memory>

// user include files
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "Geometry/HcalTowerAlgo/interface/HcalTrigTowerGeometry.h"
#include "DataFormats/HcalDetId/interface/HcalTrigTowerDetId.h"

#include "CalibCalorimetry/CaloTPG/src/CaloTPGTranscoderULUT.h"
#include "CalibFormats/HcalObjects/interface/HcalTPGRecord.h"
#include "CalibFormats/HcalObjects/interface/HcalTPGCoder.h"
#include "CalibFormats/CaloTPG/interface/CaloTPGRecord.h"

//debug stuff
#include <iostream>
using namespace std;

class CutAssessment : public edm::EDAnalyzer
{
public:
    explicit CutAssessment(const edm::ParameterSet&);
    ~CutAssessment();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:

    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;
    TFile* PrimDigiCut;
    TH1 *EtPassedElectron, *EtPassedNonElectron, *EtRatioPassedElectron, *EtRatioPassedNonElectron;
    TH1 *EtAllPassedElectron, *EtAllRatioPassedElectron;
    int Passed, Failed, Backround, TotalReal;
    int PassedAll, BackroundAll;
    Genconverter* Real2I;
    bool ElecIniEvent;

    struct PrimDigiUse
    {
        bool LocMaxPos[50][73];
        double MaxEtPos[50][73];
        bool LocMaxNeg[50][73];
        double MaxEtNeg[50][73];
    };

    PrimDigiUse PrimInfo(const edm::Event& iEvent, const edm::EventSetup&, bool UseFine);

    //struct PrimDigiUse PrimInfo(const edm::Event& iEvent);


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

CutAssessment::CutAssessment(const edm::ParameterSet& iConfig)
{
    //now do what ever initialization is needed
    //PrimDigiCut = new TFile("PrimCut.root", "RECREATE");
    //PrimDigiCut->cd();
    edm::Service<TFileService> fs;
    TFileDirectory fine = fs->mkdir("finebit");
    TFileDirectory nofine = fs->mkdir("no_fine");
    EtPassedElectron = fine.make<TH1F>("Et_of_Electron", "Et_of_Electron", 120, 0, 120);
    EtPassedNonElectron = fine.make<TH1F>("Et_of_NonElectron", "Et_of_NonElectron", 120, 0, 120);
    EtAllPassedElectron = nofine.make<TH1F>("Et_of_Electron", "Et_of_Electron", 120, 0, 120);

    EtRatioPassedElectron = fine.make<TH1F>("Et_Ratio_of_Electron", "Et_Ratio_of_Electron", 100, 0, 1.01);
    EtRatioPassedNonElectron = fine.make<TH1F>("Et_Ratio_of_NonElectron", "Et_Ratio_of_NonElectron", 100, 0, 1.01);
    EtAllRatioPassedElectron = nofine.make<TH1F>("Et_Ratio_of_Electron", "Et_Ratio_of_Electron", 100, 0, 1.01);
    Passed = 0;
    Failed = 0;
    Backround = 0;
    TotalReal = 0;
    cout << "test 2" << endl;
}

CutAssessment::~CutAssessment()
{
    cout << " total number of electrons is " << TotalReal << endl;
    cout << "so the total number of passed electrons are" << Passed << endl;
    cout << "total that passed without featurebit" << PassedAll<<endl;
    cout << " and the total number of backround is " << Backround << endl;
    cout << " and the total number of backround without feature is " << BackroundAll << endl;
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------

void
CutAssessment::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    ElecIniEvent = false;
    bool fine = false;
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





    }
    for(int q = 0; q < 2; q++)
    {
        PrimDigiUse PCuts = PrimInfo(iEvent, iSetup, fine);

        int numEPass = 0;

        for(int ieta = 0; ieta < 10; ieta++)
        {
            for(int iphi = 1; iphi < 72; iphi += 2)
            {




                if(PCuts.LocMaxPos[ieta][iphi] && gparts.find(make_pair(ieta, iphi)) != gparts.end())//says if an electron hit in the isolated region.
                {
                    numEPass++;
                    if(fine)Passed++;
                    else PassedAll++;



                    if(numEPass > numElect)
                    {
                        cout << "our ieta is" << ieta;
                        cout << "our total number of numElect is " << numElect << endl;
                        int k;
                        cin >> k;
                    }

                    //cout << "and the center is " << PCuts.MaxEtPos[ieta][iphi] << endl;
                    double ratio = 0;
                    double bottom = 0;
                    if(iphi != 1 && iphi != 71)
                    {
                        bottom = double (PCuts.MaxEtPos[ieta][iphi - 2] + PCuts.MaxEtPos[ieta][iphi] + PCuts.MaxEtPos[ieta][iphi + 2]);
                        if(ieta != 9) bottom += double (PCuts.MaxEtPos[ieta + 1][iphi - 2] + PCuts.MaxEtPos[ieta + 1][iphi] + PCuts.MaxEtPos[ieta + 1][iphi + 2]);
                        if(ieta != 0)bottom += double (PCuts.MaxEtPos[ieta - 1][iphi - 2] + PCuts.MaxEtPos[ieta - 1][iphi] + PCuts.MaxEtPos[ieta - 1][iphi + 2]);
                    } else if(iphi == 1)
                    {
                        bottom = double (PCuts.MaxEtPos[ieta][71] + PCuts.MaxEtPos[ieta][1] + PCuts.MaxEtPos[ieta][3]);
                        if(ieta != 9) bottom += double (PCuts.MaxEtPos[ieta + 1][71] + PCuts.MaxEtPos[ieta + 1][1] + PCuts.MaxEtPos[ieta + 1][3]);
                        if(ieta != 0)bottom += double (PCuts.MaxEtPos[ieta - 1][71] + PCuts.MaxEtPos[ieta - 1][1] + PCuts.MaxEtPos[ieta - 1][3]);

                    } else if(iphi == 71)
                    {
                        bottom = double (PCuts.MaxEtPos[ieta][69] + PCuts.MaxEtPos[ieta][71] + PCuts.MaxEtPos[ieta][1]);
                        if(ieta != 9) bottom += double (PCuts.MaxEtPos[ieta + 1][69] + PCuts.MaxEtPos[ieta + 1][71] + PCuts.MaxEtPos[ieta + 1][1]);
                        if(ieta != 0)bottom += double (PCuts.MaxEtPos[ieta - 1][69] + PCuts.MaxEtPos[ieta - 1][71] + PCuts.MaxEtPos[ieta - 1][1]);
                    }
                    ratio = double(PCuts.MaxEtPos[ieta][iphi]) / (bottom);
                    if(fine)EtPassedElectron->Fill(PCuts.MaxEtPos[ieta][iphi]);
                    else EtAllPassedElectron->Fill(PCuts.MaxEtPos[ieta][iphi]);
                    if(fine)EtRatioPassedElectron->Fill(ratio);
                    else EtAllRatioPassedElectron->Fill(ratio);
                    if(ratio > 1 || ratio < 0)
                    {
                        cout << " Positive ieta " << ieta << endl;
                        cout << " and our phi is " << iphi << endl;
                    }

                } else if(PCuts.LocMaxPos[ieta][iphi]&&!ElecIniEvent)
                {
                    if(fine)Backround++;
                    else BackroundAll++;
                    EtPassedNonElectron->Fill(PCuts.MaxEtPos[ieta][iphi]);

                    double bottom = 0;
                    double ratio = 0;
                    if(iphi != 1 && iphi != 71)
                    {
                        bottom = double (PCuts.MaxEtPos[ieta][iphi - 2] + PCuts.MaxEtPos[ieta][iphi] + PCuts.MaxEtPos[ieta][iphi + 2]);
                        if(ieta != 9) bottom += double (PCuts.MaxEtPos[ieta + 1][iphi - 2] + PCuts.MaxEtPos[ieta + 1][iphi] + PCuts.MaxEtPos[ieta + 1][iphi + 2]);
                        if(ieta != 0)bottom += double (PCuts.MaxEtPos[ieta - 1][iphi - 2] + PCuts.MaxEtPos[ieta - 1][iphi] + PCuts.MaxEtPos[ieta - 1][iphi + 2]);
                    } else if(iphi == 1)
                    {
                        bottom = double (PCuts.MaxEtPos[ieta][71] + PCuts.MaxEtPos[ieta][1] + PCuts.MaxEtPos[ieta][3]);
                        if(ieta != 9) bottom += double (PCuts.MaxEtPos[ieta + 1][71] + PCuts.MaxEtPos[ieta + 1][1] + PCuts.MaxEtPos[ieta + 1][3]);
                        if(ieta != 0)bottom += double (PCuts.MaxEtPos[ieta - 1][71] + PCuts.MaxEtPos[ieta - 1][1] + PCuts.MaxEtPos[ieta - 1][3]);

                    } else if(iphi == 71)
                    {
                        bottom = double (PCuts.MaxEtPos[ieta][69] + PCuts.MaxEtPos[ieta][71] + PCuts.MaxEtPos[ieta][1]);
                        if(ieta != 9) bottom += double (PCuts.MaxEtPos[ieta + 1][69] + PCuts.MaxEtPos[ieta + 1][71] + PCuts.MaxEtPos[ieta + 1][1]);
                        if(ieta != 0)bottom += double (PCuts.MaxEtPos[ieta - 1][69] + PCuts.MaxEtPos[ieta - 1][71] + PCuts.MaxEtPos[ieta - 1][1]);
                    }
                    ratio = double(PCuts.MaxEtPos[ieta][iphi]) / (bottom);
                    EtRatioPassedNonElectron->Fill(ratio);
                    if(ratio < .1)
                    {
                        cout << "so our ratio is " << ratio << endl;
                        cout << "bottom is " << bottom << endl;
                        cout << "and our ieta is " << ieta << endl;
                        cout << " our iphi is " << iphi << endl << endl;
                    }
                }

                if(PCuts.LocMaxNeg[ieta][iphi] && gparts.find(make_pair(-ieta, iphi)) != gparts.end())//says if an electron hit in the isolated region.
                {
                    numEPass++;

                    if(numEPass > numElect)
                    {
                        cout << "our ieta is -" << ieta;
                        cout << "our total number of numElect is " << numElect << endl;
                        int k;
                        cin >> k;
                    }
                    if(fine)Passed++;
                    else PassedAll++;


                    double ratio = 0;
                    double bottom = 0;
                    if(iphi != 1 && iphi != 71)
                    {
                        bottom = double (PCuts.MaxEtNeg[ieta][iphi - 2] + PCuts.MaxEtNeg[ieta][iphi] + PCuts.MaxEtNeg[ieta][iphi + 2]);
                        if(ieta != 9) bottom += double (PCuts.MaxEtNeg[ieta + 1][iphi - 2] + PCuts.MaxEtNeg[ieta + 1][iphi] + PCuts.MaxEtNeg[ieta + 1][iphi + 2]);
                        if(ieta != 0)bottom += double (PCuts.MaxEtNeg[ieta - 1][iphi - 2] + PCuts.MaxEtNeg[ieta - 1][iphi] + PCuts.MaxEtNeg[ieta - 1][iphi + 2]);
                    } else if(iphi == 1)
                    {
                        bottom = double (PCuts.MaxEtNeg[ieta][71] + PCuts.MaxEtNeg[ieta][1] + PCuts.MaxEtNeg[ieta][3]);
                        if(ieta != 9) bottom += double (PCuts.MaxEtNeg[ieta + 1][71] + PCuts.MaxEtNeg[ieta + 1][1] + PCuts.MaxEtNeg[ieta + 1][3]);
                        if(ieta != 0)bottom += double (PCuts.MaxEtNeg[ieta - 1][71] + PCuts.MaxEtNeg[ieta - 1][1] + PCuts.MaxEtNeg[ieta - 1][3]);

                    } else if(iphi == 71)
                    {
                        bottom = double (PCuts.MaxEtNeg[ieta][69] + PCuts.MaxEtNeg[ieta][71] + PCuts.MaxEtNeg[ieta][1]);
                        if(ieta != 9) bottom += double (PCuts.MaxEtNeg[ieta + 1][69] + PCuts.MaxEtNeg[ieta + 1][71] + PCuts.MaxEtNeg[ieta + 1][1]);
                        if(ieta != 0)bottom += double (PCuts.MaxEtNeg[ieta - 1][69] + PCuts.MaxEtNeg[ieta - 1][71] + PCuts.MaxEtNeg[ieta - 1][1]);
                    }
                    ratio = double(PCuts.MaxEtNeg[ieta][iphi]) / (bottom);
                    if(fine)EtPassedElectron->Fill(PCuts.MaxEtNeg[ieta][iphi]);
                    else EtAllPassedElectron->Fill(PCuts.MaxEtNeg[ieta][iphi]);
                    if(fine) EtRatioPassedElectron->Fill(ratio);
                    else EtAllRatioPassedElectron->Fill(ratio);
                    if(ratio > 1 || ratio < 0)
                    {
                        cout << " Positive -ieta " << ieta << endl;
                        cout << " and our phi is " << iphi << endl;
                    }



                } else if(PCuts.LocMaxNeg[ieta][iphi]&&!ElecIniEvent)
                {
                    if(fine)Backround++;
                    else BackroundAll++;
                    EtPassedNonElectron->Fill(PCuts.MaxEtNeg[ieta][iphi]);
                    double ratio = 0;
                    double bottom = 0;
                    if(iphi != 1 && iphi != 71)
                    {
                        bottom = double (PCuts.MaxEtNeg[ieta][iphi - 2] + PCuts.MaxEtNeg[ieta][iphi] + PCuts.MaxEtNeg[ieta][iphi + 2]);
                        if(ieta != 9) bottom += double (PCuts.MaxEtNeg[ieta + 1][iphi - 2] + PCuts.MaxEtNeg[ieta + 1][iphi] + PCuts.MaxEtNeg[ieta + 1][iphi + 2]);
                        if(ieta != 0)bottom += double (PCuts.MaxEtNeg[ieta - 1][iphi - 2] + PCuts.MaxEtNeg[ieta - 1][iphi] + PCuts.MaxEtNeg[ieta - 1][iphi + 2]);
                    } else if(iphi == 1)
                    {
                        bottom = double (PCuts.MaxEtNeg[ieta][71] + PCuts.MaxEtNeg[ieta][1] + PCuts.MaxEtNeg[ieta][3]);
                        if(ieta != 9) bottom += double (PCuts.MaxEtNeg[ieta + 1][71] + PCuts.MaxEtNeg[ieta + 1][1] + PCuts.MaxEtNeg[ieta + 1][3]);
                        if(ieta != 0)bottom += double (PCuts.MaxEtNeg[ieta - 1][71] + PCuts.MaxEtNeg[ieta - 1][1] + PCuts.MaxEtNeg[ieta - 1][3]);

                    } else if(iphi == 71)
                    {
                        bottom = double (PCuts.MaxEtNeg[ieta][69] + PCuts.MaxEtNeg[ieta][71] + PCuts.MaxEtNeg[ieta][1]);
                        if(ieta != 9) bottom += double (PCuts.MaxEtNeg[ieta + 1][69] + PCuts.MaxEtNeg[ieta + 1][71] + PCuts.MaxEtNeg[ieta + 1][1]);
                        if(ieta != 0)bottom += double (PCuts.MaxEtNeg[ieta - 1][69] + PCuts.MaxEtNeg[ieta - 1][71] + PCuts.MaxEtNeg[ieta - 1][1]);
                    }
                    ratio = double(PCuts.MaxEtNeg[ieta][iphi]) / (bottom);
                    EtRatioPassedNonElectron->Fill(ratio);

                    if(ratio < .1)
                    {
                        cout << "so our ratio is " << ratio << endl;
                        cout << "bottom is " << bottom << endl;
                        cout << "and our -ieta is " << ieta << endl;
                        cout << " our iphi is " << iphi << endl << endl;
                    }
                }
            }
        }
        fine = true;
    }
}

// ------------ method called once each job just before starting event loop  ------------

void
CutAssessment::beginJob() { }

// ------------ method called once each job just after ending the event loop  ------------

void
CutAssessment::endJob() { }

CutAssessment::PrimDigiUse CutAssessment::PrimInfo(const edm::Event& iEvent, const edm::EventSetup& eventSetup, bool UseFine)
{
    edm::ESHandle<CaloTPGTranscoder> outTranscoder;
    eventSetup.get<CaloTPGRecord>().get(outTranscoder);
    outTranscoder->setup(eventSetup, CaloTPGTranscoder::HcalTPG);

    PrimDigiUse CutInfo;
    for(int ietac = 0; ietac < 10; ietac++)
    {
        for(int iphi = 1; iphi < 72; iphi += 2)
        {
            //cout<<" okay last here is iphi"<< iphi<<endl;
            CutInfo.MaxEtPos[ietac][iphi] = 0;
            CutInfo.LocMaxPos[ietac][iphi] = false;
            CutInfo.MaxEtNeg[ietac][iphi] = 0;
            CutInfo.LocMaxNeg[ietac][iphi] = false;
        }
    }

    edm::Handle<HcalTrigPrimDigiCollection> hfpr_digi;
    iEvent.getByLabel("simHcalTriggerPrimitiveDigis", hfpr_digi);
    for(HcalTrigPrimDigiCollection::const_iterator tp = hfpr_digi->begin(); tp != hfpr_digi->end(); ++tp)
    {
        if(abs(tp->id().ieta()) < 30 || abs(tp->id().ieta()) > 39) continue;

        double MaxEt = 0;

        for(int i = 0; i < tp->size(); i++)
        {
            if(MaxEt < outTranscoder->hcaletValue(tp->id(), (*tp)[i])&&((*tp)[i].fineGrain() || !UseFine))
            {
                MaxEt = outTranscoder->hcaletValue(tp->id(), (*tp)[i]);

            }
        }
        if(tp->id().ieta() > 0)//delete
        {


            CutInfo.MaxEtPos[(tp->id().ieta() - 30)][(tp->id().iphi())] = MaxEt;
        } else
        {
            CutInfo.MaxEtNeg[(abs(tp->id().ieta()) - 30)][(tp->id().iphi())] = MaxEt;

        }
    }

    for(int ietac = 0; ietac < 10; ietac++)//this handles the positive eta's
    {
        for(int iphi = 3; iphi < 71; iphi += 2)
        {
            bool upleft = (CutInfo.MaxEtPos[ietac][iphi] > CutInfo.MaxEtPos[ietac + 1][iphi - 2]) || ietac == 9;
            bool upcent = CutInfo.MaxEtPos[ietac][iphi] > CutInfo.MaxEtPos[ietac + 1][iphi] || ietac == 9;
            bool upright = CutInfo.MaxEtPos[ietac][iphi] > CutInfo.MaxEtPos[ietac + 1][iphi + 2] || ietac == 9;
            bool left = CutInfo.MaxEtPos[ietac][iphi] > CutInfo.MaxEtPos[ietac][iphi - 2];
            bool right = CutInfo.MaxEtPos[ietac][iphi] > CutInfo.MaxEtPos[ietac][iphi + 2];
            bool botleft = CutInfo.MaxEtPos[ietac][iphi] > CutInfo.MaxEtPos[ietac - 1][iphi - 2] || ietac == 0;
            bool botcent = CutInfo.MaxEtPos[ietac][iphi] > CutInfo.MaxEtPos[ietac - 1][iphi] || ietac == 0;
            bool botright = CutInfo.MaxEtPos[ietac][iphi] > CutInfo.MaxEtPos[ietac - 1][iphi + 2] || ietac == 0;

            CutInfo.LocMaxPos[ietac][iphi] = (upleft && upcent && upright && right && left && botleft && botcent && botright && CutInfo.MaxEtPos[ietac][iphi] > 5);
            //if(CutInfo.LocMax[ietac][iphi])cout<<"we have a ietac of "<<ietac<<endl;

        }
        bool upleft = CutInfo.MaxEtPos[ietac][1] > CutInfo.MaxEtPos[ietac + 1][71] || ietac == 9;
        bool upcent = CutInfo.MaxEtPos[ietac][1] > CutInfo.MaxEtPos[ietac + 1][1] || ietac == 9;
        bool upright = CutInfo.MaxEtPos[ietac][1] > CutInfo.MaxEtPos[ietac + 1][3] || ietac == 9;
        bool left = CutInfo.MaxEtPos[ietac][1] > CutInfo.MaxEtPos[ietac][71];
        bool right = CutInfo.MaxEtPos[ietac][1] > CutInfo.MaxEtPos[ietac][3];
        bool botleft = CutInfo.MaxEtPos[ietac][1] > CutInfo.MaxEtPos[ietac - 1][71] || ietac == 0;
        bool botcent = CutInfo.MaxEtPos[ietac][1] > CutInfo.MaxEtPos[ietac - 1][1] || ietac == 0;
        bool botright = CutInfo.MaxEtPos[ietac][1] > CutInfo.MaxEtPos[ietac - 1][3] || ietac == 0;
        CutInfo.LocMaxPos[ietac][1] = (upleft && upcent && upright && right && left && botleft && botcent && botright);

        upleft = CutInfo.MaxEtPos[ietac][71] > CutInfo.MaxEtPos[ietac + 1][69] || ietac == 9;
        upcent = CutInfo.MaxEtPos[ietac][71] > CutInfo.MaxEtPos[ietac + 1][71] || ietac == 9;
        upright = CutInfo.MaxEtPos[ietac][71] > CutInfo.MaxEtPos[ietac + 1][1] || ietac == 9;
        left = CutInfo.MaxEtPos[ietac][71] > CutInfo.MaxEtPos[ietac][69];
        right = CutInfo.MaxEtPos[ietac][71] > CutInfo.MaxEtPos[ietac][1];
        botleft = CutInfo.MaxEtPos[ietac][71] > CutInfo.MaxEtPos[ietac - 1][69] || ietac == 0;
        botcent = CutInfo.MaxEtPos[ietac][71] > CutInfo.MaxEtPos[ietac - 1][71] || ietac == 0;
        botright = CutInfo.MaxEtPos[ietac][71] > CutInfo.MaxEtPos[ietac - 1][1] || ietac == 0;
        CutInfo.LocMaxPos[ietac][71] = (upleft && upcent && upright && left && botleft && botcent && botright && CutInfo.MaxEtPos[ietac][71] > 5);

    }

    for(int ietac = 0; ietac < 10; ietac++)//this handles the negative eta's
    {
        for(int iphi = 3; iphi < 71; iphi += 2)
        {
            bool upleft = CutInfo.MaxEtNeg[ietac][iphi] > CutInfo.MaxEtNeg[ietac + 1][iphi - 2] || ietac == 9;
            bool upcent = CutInfo.MaxEtNeg[ietac][iphi] > CutInfo.MaxEtNeg[ietac + 1][iphi] || ietac == 9;
            bool upright = CutInfo.MaxEtNeg[ietac][iphi] > CutInfo.MaxEtNeg[ietac + 1][iphi + 2] || ietac == 9;
            bool left = CutInfo.MaxEtNeg[ietac][iphi] > CutInfo.MaxEtNeg[ietac][iphi - 2];
            bool right = CutInfo.MaxEtNeg[ietac][iphi] > CutInfo.MaxEtNeg[ietac][iphi + 2];
            bool botleft = CutInfo.MaxEtNeg[ietac][iphi] > CutInfo.MaxEtNeg[ietac - 1][iphi - 2] || ietac == 0 || ietac == 19;
            bool botcent = CutInfo.MaxEtNeg[ietac][iphi] > CutInfo.MaxEtNeg[ietac - 1][iphi] || ietac == 0 || ietac == 19;
            bool botright = CutInfo.MaxEtNeg[ietac][iphi] > CutInfo.MaxEtNeg[ietac - 1][iphi + 2] || ietac == 0 || ietac == 19;



            CutInfo.LocMaxNeg[ietac][iphi] = (upleft && upcent && upright && right && left && botleft && botcent && botright && CutInfo.MaxEtNeg[ietac][iphi] > 5);

        }
        bool upleft = CutInfo.MaxEtNeg[ietac][1] > CutInfo.MaxEtNeg[ietac + 1][71] || ietac == 10 || ietac == 9;
        bool upcent = CutInfo.MaxEtNeg[ietac][1] > CutInfo.MaxEtNeg[ietac + 1][1] || ietac == 10 || ietac == 9;
        bool upright = CutInfo.MaxEtNeg[ietac][1] > CutInfo.MaxEtNeg[ietac + 1][3] || ietac == 10 || ietac == 9;
        bool left = CutInfo.MaxEtNeg[ietac][1] > CutInfo.MaxEtNeg[ietac][71];
        bool right = CutInfo.MaxEtNeg[ietac][1] > CutInfo.MaxEtNeg[ietac][3];
        bool botleft = CutInfo.MaxEtNeg[ietac][1] > CutInfo.MaxEtNeg[ietac - 1][71] || ietac == 0;
        bool botcent = CutInfo.MaxEtNeg[ietac][1] > CutInfo.MaxEtNeg[ietac - 1][1] || ietac == 0;
        bool botright = CutInfo.MaxEtNeg[ietac][1] > CutInfo.MaxEtNeg[ietac - 1][3] || ietac == 0;
        CutInfo.LocMaxNeg[ietac][1] = (upleft && upcent && upright && right && left && botleft && botcent && botright && CutInfo.MaxEtNeg[ietac][1] > 5);

        upleft = CutInfo.MaxEtNeg[ietac][71] > CutInfo.MaxEtNeg[ietac + 1][69] || ietac == 9;
        upcent = CutInfo.MaxEtNeg[ietac][71] > CutInfo.MaxEtNeg[ietac + 1][71] || ietac == 9;
        upright = CutInfo.MaxEtNeg[ietac][71] > CutInfo.MaxEtNeg[ietac + 1][1] || ietac == 9;
        left = CutInfo.MaxEtNeg[ietac][71] > CutInfo.MaxEtNeg[ietac][69];
        right = CutInfo.MaxEtNeg[ietac][71] > CutInfo.MaxEtNeg[ietac][1];
        botleft = CutInfo.MaxEtNeg[ietac][71] > CutInfo.MaxEtNeg[ietac - 1][69] || ietac == 0;
        botcent = CutInfo.MaxEtNeg[ietac][71] > CutInfo.MaxEtNeg[ietac - 1][71] || ietac == 0;
        botright = CutInfo.MaxEtNeg[ietac][71] > CutInfo.MaxEtNeg[ietac - 1][1] || ietac == 0;
        CutInfo.LocMaxNeg[ietac][71] = (upleft && upcent && upright && left && botleft && botcent && botright && CutInfo.MaxEtNeg[ietac][71] > 5);

    }

    return (CutInfo);
}

void
CutAssessment::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(CutAssessment);

#endif
