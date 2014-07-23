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
    TH1 *Seedoverbackround;
    TH2 *Etvrs3X3Fine, *Etvrs3X3NoFine;
    TH2 *SeedPrimVrs3X3Fine, *SeedPrimVrs3X3NoFine;
    TH2 *ElectonEtCompEt;
    TH2 *comp3X3VS5X5Fine, *comp3X3VS5X5NoFine;
    double GenEtPos[20][72];
    double GenEtNeg[20][72];
    double CompEtPos[20][72];
    double CompEtNeg[20][72];
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
    Etvrs3X3Fine = fine.make<TH2F>("centervrsdenominator", "centervrsdenominator;et;denominator", 120, 0, 120, 120, 0, 120);
    Etvrs3X3NoFine = nofine.make<TH2F>("centervrsdenominator", "centervrsdenominator;et;denominator", 120, 0, 120, 120, 0, 120);
    ElectonEtCompEt = fs->make<TH2F>("compressed_VS_Generated", "compressed_VS_Generated; Comp; Gen(Et)", 120, 0, 120, 150, 0, 150);
    //ElectronETunC = fs->make<TH2F>("ETuncomp_VS_Generated", "ETuncomp_VS_Generated; Uncompressed; Gen(Et)", 120, 0, 120, 150, 0, 150);
    SeedPrimVrs3X3Fine = fine.make<TH2F>("Seed_Plus_EdgeVs3X3", "Seed_Plus_EdgeVs3X3; Seed Energy; 3X3", 100, 0, 100, 100, 0, 100);
    SeedPrimVrs3X3NoFine = nofine.make<TH2F>("Seed_Plus_EdgeVs3X3", "Seed_Plus_EdgeVs3X3; Seed Energy; 3X3", 100, 0, 100, 100, 0, 100);
    comp3X3VS5X5Fine = fine.make<TH2F>("3X3VS5X5F", "3X3VS5X5F; 3X3; 5X5", 100, 0, 100, 100, 0, 100);
    comp3X3VS5X5NoFine = nofine.make<TH2F>("3X3VS5X5F", "3X3VS5X5F; 3X3; 5X5", 100, 0, 100, 100, 0, 100);
    Passed = 0;
    Failed = 0;
    Backround = 0;
    TotalReal = 0;
    PassedAll = 0;
    BackroundAll = 0;
    cout << "test 2" << endl;
}

CutAssessment::~CutAssessment()
{
    cout << " total number of electrons is " << TotalReal << endl;
    cout << "so the total number of passed electrons are" << Passed << endl;
    cout << "total that passed without featurebit" << PassedAll << endl;
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

    if(TotalReal < PassedAll)
    {
        cout << "more passed then is there was" << endl;
        int k;
        cin >> k;
    }
    ElecIniEvent = false;
    bool fine = false;
    using namespace edm;
    int numElect = 0;
    set<pair<int, int> > gparts; //makes a set that contains locations of all generated electrons
    edm::Handle<reco::GenParticleCollection> genInfo;
    iEvent.getByLabel("genParticles", genInfo);
    Real2I = new Genconverter();



    for(int i = 0; i < 10; i++)
    {
        for(int k = 0; k < 72; k += 2)
        {
            GenEtPos[i][k] = 0;
            GenEtNeg[i][k] = 0;
            CompEtPos[i][k] = 0;
            CompEtNeg[i][k] = 0;
        }
    }



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
            if(gieta > 0)
            {
                GenEtPos[gieta - 30][giphi] = genin.energy() / cosh(genin.eta());
                gparts.insert(make_pair(gieta, giphi));
            }// should make say whether a electron hit the sensor
            else
            {
                GenEtNeg[abs(gieta) - 30][giphi] = genin.energy() / cosh(genin.eta());
                gparts.insert(make_pair(gieta, giphi));
            }
        }

        if(TotalReal < PassedAll)
        {
            cout << "test 4" << endl;
            int k;
            cin >> k;
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




                if(PCuts.LocMaxPos[ieta][iphi] && gparts.find(make_pair(ieta + 30, iphi)) != gparts.end())//says if an electron hit in the isolated region.
                {
                    numEPass++;
                    if(fine)Passed++;
                    else PassedAll++;
                    if(PCuts.MaxEtPos[ieta][iphi] < 10)
                        cout << "well thats odd" << endl;


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

                    double T5X5 = 0;
                    double seed = 0;
                    for(int i = -2; i < 3; i++)
                    {
                        for(int w = -2; w < 3; w++)
                        {
                            int iphip = iphi + i * 2;
                            int ietap = ieta + w;
                            if(iphip == 73)iphip = 1;
                            if(iphip == 75)iphip = 3;
                            if(iphip == -1)iphip = 71;
                            if(iphip == -3)iphip = 69;
                            if(ietap > 9 || ietap < 0)
                            {
                                ietap = 0;
                                iphip = 0;
                            }
                            if(PCuts.MaxEtPos[ietap][iphip] > seed && abs(i) < 2 && abs(w) < 2)
                            {
                              if(i!=0&&w!=0) seed = PCuts.MaxEtPos[ietap][iphip];
                                bottom += PCuts.MaxEtPos[ietap][iphip];
                            }
                            T5X5 += PCuts.MaxEtPos[ietap][iphip];
                        }
                    }

                    if(seed*2<PCuts.MaxEtPos[ieta][iphi])seed=0;
                    ratio = double(PCuts.MaxEtPos[ieta][iphi]) / (bottom);
                    if(fine)
                    {
                        EtPassedElectron->Fill(PCuts.MaxEtPos[ieta][iphi]);
                        EtRatioPassedElectron->Fill(ratio);
                        comp3X3VS5X5Fine->Fill(bottom, T5X5);
                    } else
                    {
                        EtAllPassedElectron->Fill(PCuts.MaxEtPos[ieta][iphi]);

                        EtAllRatioPassedElectron->Fill(ratio);


                        ElectonEtCompEt->Fill(PCuts.MaxEtPos[ieta][iphi], GenEtPos[ieta][iphi]);
                        comp3X3VS5X5NoFine->Fill(bottom, T5X5);

                    }
                    if(ratio > 1 || ratio < 0)
                    {
                        cout << " Positive ieta " << ieta << endl;
                        cout << " and our phi is " << iphi << endl;
                    }

                }
                
                
                
                
                
                
                
                
                ///////////////
                if(PCuts.LocMaxNeg[ieta][iphi] && gparts.find(make_pair(-ieta - 30, iphi)) != gparts.end())//says if an electron hit in the isolated region.
                {
                    numEPass++;

                    if(PCuts.MaxEtNeg[ieta][iphi] < 10)
                        cout << "so the negatives are the problem" << endl;

                    if(numEPass > numElect)
                    {
                        cout << "our ieta is -" << ieta << " and our iphi is :" << iphi << endl;
                        cout << "our total number of numElect is " << numElect << endl << endl;
                        //int k;
                        //cin >> k;
                    }
                    if(fine)Passed++;
                    else PassedAll++;


                    double ratio = 0;
                    double bottom = 0;
                    double T5X5 = 0;
                    double seed = 0;
                    
                    for(int i = -2; i < 3; i++)
                    {
                        for(int w = -2; w < 3; w++)
                        {
                            int iphip = iphi + i * 2;
                            int ietap = ieta + w;
                            if(iphip == 73)iphip = 1;
                            if(iphip == 75)iphip = 3;
                            if(iphip == -1)iphip = 71;
                            if(iphip == -3)iphip = 69;
                            if(ietap > 9 || ietap < 0)
                            {
                                ietap = 0;
                                iphip = 0;
                            }
                            if(PCuts.MaxEtNeg[ietap][iphip] > seed && abs(i) < 2 && abs(w) < 2)
                            {
                              if(i!=0&&w!=0) seed = PCuts.MaxEtNeg[ietap][iphip];
                                bottom += PCuts.MaxEtNeg[ietap][iphip];
                            }
                            T5X5 += PCuts.MaxEtNeg[ietap][iphip];
                        }
                    }
                    if(seed*2<PCuts.MaxEtNeg[ieta][iphi])seed=0;
                    ratio = double(PCuts.MaxEtNeg[ieta][iphi]) / (bottom);
                    if(fine)EtPassedElectron->Fill(PCuts.MaxEtNeg[ieta][iphi]);
                    else EtAllPassedElectron->Fill(PCuts.MaxEtNeg[ieta][iphi]);



                    if(fine)
                    {
                        EtRatioPassedElectron->Fill(ratio);
                        Etvrs3X3Fine->Fill(PCuts.MaxEtNeg[ieta][iphi], bottom);
                        comp3X3VS5X5Fine->Fill(bottom, T5X5);
                    } else
                    {
                        EtAllRatioPassedElectron->Fill(ratio);
                        comp3X3VS5X5NoFine->Fill(bottom, T5X5);
                        Etvrs3X3NoFine->Fill(PCuts.MaxEtNeg[ieta][iphi], bottom);
                        if(PCuts.MaxEtNeg[ieta][iphi] > 10)
                        {
                            ElectonEtCompEt->Fill(PCuts.MaxEtNeg[ieta][iphi], GenEtNeg[ieta][iphi]);
                            // ElectronETunC->Fill(PCuts.MaxEtNeg[ieta][iphi], GenEtNeg[ieta][iphi]);
                        }

                    }
                    if(ratio > 1 || ratio < 0)
                    {
                        cout << " Positive -ieta " << ieta << endl;
                        cout << " and our phi is " << iphi << endl;
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
    double MinSeed = 10;
    bool featurebitPos[11][74];
    bool featurebitNeg[11][74];
    //double MinEt =0;
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
            featurebitPos[ietac][iphi] = false;
            featurebitNeg[ietac][iphi] = false;
        }
    }
    CutInfo.MaxEtPos[0][0] = 0; //so I am putting one of these to zero so that I could just change to a specific ieta iphi if I want zero.
    CutInfo.MaxEtNeg[0][0] = 0;
    edm::Handle<HcalTrigPrimDigiCollection> hfpr_digi;
    iEvent.getByLabel("simHcalTriggerPrimitiveDigis", hfpr_digi);
    //cout<<" test 2"<<endl;
    for(HcalTrigPrimDigiCollection::const_iterator tp = hfpr_digi->begin(); tp != hfpr_digi->end(); ++tp)
    {
        if(abs(tp->id().ieta()) < 30 || abs(tp->id().ieta()) > 39) continue;

        double MaxEt = 0;
        //int MaxCEt = 0;

        //if(!tp->SOI_fineGrain() && UseFine)continue;

        MaxEt = outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt());
        if(MaxEt * 2 != tp->SOI_compressedEt())cout << " the correction factor is :" << MaxEt / tp->SOI_compressedEt() << endl;
        if(tp->id().ieta() > 0)//delete
        {

            featurebitPos[(tp->id().ieta() - 30)][(tp->id().iphi())] = (tp->SOI_fineGrain() || !UseFine);
            CutInfo.MaxEtPos[(tp->id().ieta() - 30)][(tp->id().iphi())] = MaxEt;
        } else
        {
            featurebitNeg[(abs(tp->id().ieta()) - 30)][(tp->id().iphi())] = (tp->SOI_fineGrain() || !UseFine);
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

            CutInfo.LocMaxPos[ietac][iphi] = (upleft && upcent && upright && right && left && botleft && botcent && botright && CutInfo.MaxEtPos[ietac][iphi] > MinSeed && featurebitPos[ietac][iphi]);
            //if(featurebitPos[ietac][iphi])cout<<"well we got here"<<endl;

        }
        bool upleft = CutInfo.MaxEtPos[ietac][1] > CutInfo.MaxEtPos[ietac + 1][71] || ietac == 9;
        bool upcent = CutInfo.MaxEtPos[ietac][1] > CutInfo.MaxEtPos[ietac + 1][1] || ietac == 9;
        bool upright = CutInfo.MaxEtPos[ietac][1] > CutInfo.MaxEtPos[ietac + 1][3] || ietac == 9;
        bool left = CutInfo.MaxEtPos[ietac][1] > CutInfo.MaxEtPos[ietac][71];
        bool right = CutInfo.MaxEtPos[ietac][1] > CutInfo.MaxEtPos[ietac][3];
        bool botleft = CutInfo.MaxEtPos[ietac][1] > CutInfo.MaxEtPos[ietac - 1][71] || ietac == 0;
        bool botcent = CutInfo.MaxEtPos[ietac][1] > CutInfo.MaxEtPos[ietac - 1][1] || ietac == 0;
        bool botright = CutInfo.MaxEtPos[ietac][1] > CutInfo.MaxEtPos[ietac - 1][3] || ietac == 0;
        CutInfo.LocMaxPos[ietac][1] = (upleft && upcent && upright && right && left && botleft && botcent && botright && CutInfo.MaxEtPos[ietac][1] > MinSeed && featurebitPos[ietac][1]);

        upleft = CutInfo.MaxEtPos[ietac][71] > CutInfo.MaxEtPos[ietac + 1][69] || ietac == 9;
        upcent = CutInfo.MaxEtPos[ietac][71] > CutInfo.MaxEtPos[ietac + 1][71] || ietac == 9;
        upright = CutInfo.MaxEtPos[ietac][71] > CutInfo.MaxEtPos[ietac + 1][1] || ietac == 9;
        left = CutInfo.MaxEtPos[ietac][71] > CutInfo.MaxEtPos[ietac][69];
        right = CutInfo.MaxEtPos[ietac][71] > CutInfo.MaxEtPos[ietac][1];
        botleft = CutInfo.MaxEtPos[ietac][71] > CutInfo.MaxEtPos[ietac - 1][69] || ietac == 0;
        botcent = CutInfo.MaxEtPos[ietac][71] > CutInfo.MaxEtPos[ietac - 1][71] || ietac == 0;
        botright = CutInfo.MaxEtPos[ietac][71] > CutInfo.MaxEtPos[ietac - 1][1] || ietac == 0;
        CutInfo.LocMaxPos[ietac][71] = (upleft && upcent && upright && left && botleft && botcent && botright && CutInfo.MaxEtPos[ietac][71] > MinSeed && featurebitPos[ietac][71]);

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
            bool botleft = CutInfo.MaxEtNeg[ietac][iphi] > CutInfo.MaxEtNeg[ietac - 1][iphi - 2] || ietac == 0;
            bool botcent = CutInfo.MaxEtNeg[ietac][iphi] > CutInfo.MaxEtNeg[ietac - 1][iphi] || ietac == 0;
            bool botright = CutInfo.MaxEtNeg[ietac][iphi] > CutInfo.MaxEtNeg[ietac - 1][iphi + 2] || ietac == 0;



            CutInfo.LocMaxNeg[ietac][iphi] = (upleft && upcent && upright && right && left && botleft && botcent && botright && CutInfo.MaxEtNeg[ietac][iphi] > MinSeed && featurebitNeg[ietac][iphi]);

        }
        bool upleft = CutInfo.MaxEtNeg[ietac][1] > CutInfo.MaxEtNeg[ietac + 1][71] || ietac == 9;
        bool upcent = CutInfo.MaxEtNeg[ietac][1] > CutInfo.MaxEtNeg[ietac + 1][1] || ietac == 9;
        bool upright = CutInfo.MaxEtNeg[ietac][1] > CutInfo.MaxEtNeg[ietac + 1][3] || ietac == 9;
        bool left = CutInfo.MaxEtNeg[ietac][1] > CutInfo.MaxEtNeg[ietac][71];
        bool right = CutInfo.MaxEtNeg[ietac][1] > CutInfo.MaxEtNeg[ietac][3];
        bool botleft = CutInfo.MaxEtNeg[ietac][1] > CutInfo.MaxEtNeg[ietac - 1][71] || ietac == 0;
        bool botcent = CutInfo.MaxEtNeg[ietac][1] > CutInfo.MaxEtNeg[ietac - 1][1] || ietac == 0;
        bool botright = CutInfo.MaxEtNeg[ietac][1] > CutInfo.MaxEtNeg[ietac - 1][3] || ietac == 0;
        CutInfo.LocMaxNeg[ietac][1] = (upleft && upcent && upright && right && left && botleft && botcent && botright && CutInfo.MaxEtNeg[ietac][1] > MinSeed && featurebitNeg[ietac][1]);

        upleft = CutInfo.MaxEtNeg[ietac][71] > CutInfo.MaxEtNeg[ietac + 1][69] || ietac == 9;
        upcent = CutInfo.MaxEtNeg[ietac][71] > CutInfo.MaxEtNeg[ietac + 1][71] || ietac == 9;
        upright = CutInfo.MaxEtNeg[ietac][71] > CutInfo.MaxEtNeg[ietac + 1][1] || ietac == 9;
        left = CutInfo.MaxEtNeg[ietac][71] > CutInfo.MaxEtNeg[ietac][69];
        right = CutInfo.MaxEtNeg[ietac][71] > CutInfo.MaxEtNeg[ietac][1];
        botleft = CutInfo.MaxEtNeg[ietac][71] > CutInfo.MaxEtNeg[ietac - 1][69] || ietac == 0;
        botcent = CutInfo.MaxEtNeg[ietac][71] > CutInfo.MaxEtNeg[ietac - 1][71] || ietac == 0;
        botright = CutInfo.MaxEtNeg[ietac][71] > CutInfo.MaxEtNeg[ietac - 1][1] || ietac == 0;
        CutInfo.LocMaxNeg[ietac][71] = (upleft && upcent && upright && left && botleft && botcent && botright && CutInfo.MaxEtNeg[ietac][71] > MinSeed && featurebitNeg[ietac][71]);

    }

    return (CutInfo);
}

void
CutAssessment::fillDescriptions(edm::ConfigurationDescriptions & descriptions)
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
