// -*- C++ -*-
//
// Package:    TestEDA/DataAnalizer
// Class:      DataAnalizer
// 
/**\class DataAnalizer DataAnalizer.cc TestEDA/DataAnalizer/plugins/DataAnalizer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
 */
//
// Original Author:  Zachary Lesko
//         Created:  Thu, 12 Jun 2014 15:30:54 GMT
//
//
#ifndef DataAnalizer_included
#define DataAnalizer_included 1


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

////////
#include "DataFormats/HLTReco/interface/TriggerEvent.h" // trigger::TriggerEvent
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h" // reco::PhotonCollection
#include "DataFormats/EgammaReco/interface/HFEMClusterShape.h" // reco::HFEMClusterShape
#include "DataFormats/EgammaReco/interface/HFEMClusterShapeAssociation.h" // reco::HFEMClusterShapeAssociationCollection
#include "DataFormats/EgammaReco/interface/HFEMClusterShapeFwd.h" // reco::HFEMClusterShapeRef,
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h" // reco::SuperClusterCollection, reco::SuperClusterRef
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h" // reco::RecoEcalCandidateCollection

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h" // reco::GsfElectron
#include "DataFormats/EgammaCandidates/interface/Photon.h" // reco::Photon
#include "DataFormats/HLTReco/interface/TriggerObject.h" // trigger::TriggerObject
#include "DataFormats/HepMCCandidate/interface/GenParticle.h" // reco::GenParticle
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h" // reco::RecoEcalCandidate

//////////reco stuff?
//lumi stuff^
#include "CalibCalorimetry/CaloTPG/src/CaloTPGTranscoderULUT.h"
#include "CalibFormats/HcalObjects/interface/HcalTPGRecord.h"
#include "CalibFormats/HcalObjects/interface/HcalTPGCoder.h"
#include "CalibFormats/CaloTPG/interface/CaloTPGRecord.h"
#include "TH1.h"
#include "TH2.h"
#include <iostream>
using namespace std;





//
// class declaration
//

class DataAnalizer : public edm::EDAnalyzer
{
public:
    explicit DataAnalizer(const edm::ParameterSet&);
    ~DataAnalizer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;


    int passed;
    int failed;
    int lumi;
    TH1 *PrimPassEt, *PrimPassRatioEt;
    TH1 *HitsEventfine, *HitsEventNo;
    TH1 *PrimAllPassEt, *PrimAllPassRatioEt;
    TH2 *NumofVetVrsTrigCanPrim, *NumofVetVrsTrigCanNoPrim;
    TH2 *PrimFineTrigCanLum, *PrimNoFineTrigCanLum;
    TH2 *SeedPrimVrs3X3Fine, *SeedPrimVrs3X3NoFine;
    TH2 *comp3X3VS5X5Fine, *comp3X3VS5X5NoFine;
    TH2 *CentVs3X3Fine, *CentVs3X3NoFine;
    struct PrimDigiUse
    {
        bool LocMaxPos[50][73];
        double MaxEtPos[50][73];
        bool LocMaxNeg[50][73];
        double MaxEtNeg[50][73];
    };

    double GenEtPos[20][72];
    double GenEtNeg[20][72];
    PrimDigiUse PrimInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup, bool UseFine);

    //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
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
TH2 *PrimPassEtvrsiEta, *PrimPassEtvrsiPhi;

DataAnalizer::DataAnalizer(const edm::ParameterSet& iConfig)
{
    edm::Service<TFileService> fs;
    TFileDirectory fine = fs->mkdir("finebit");
    TFileDirectory nofine = fs->mkdir("no_fine");
    PrimPassEt = fine.make<TH1F>("Et", "Et;Et;Count", 200, 0, 100);
    PrimPassRatioEt = fine.make<TH1F>("EtRatio", "EtRatio;Et/3X3", 200, 0, 1.001);
    PrimFineTrigCanLum = fine.make<TH2F>("Luminosity_Vrs_Trig_Candidates", "Luminosity_Vrs_Trig_Candidates; Luminosity ;Trigger_Candidates", 30, 1900, 7900, 11, 0, 11.);
    NumofVetVrsTrigCanPrim = fine.make<TH2F>("Number_Of_vertices_Vs_Trig_Candidates", "Number_Of_vertices_Vs_Trig_Candidates; Number of Vertices ;Trigger_Candidates", 50, 0, 50, 11, 0, 11.);
    NumofVetVrsTrigCanNoPrim = nofine.make<TH2F>("Number_Of_vertices_Vs_Trig_Candidates", "Number_Of_vertices_Vs_Trig_Candidates; Number of Vertices ;Trigger_Candidates", 50, 0, 50, 11, 0, 11.);
    HitsEventfine = fine.make<TH1I>("Hits", "Hits", 10.1, 0, 10);

    PrimAllPassEt = nofine.make<TH1F>("Et", "Et;Et;Count", 200, 0, 100);
    PrimAllPassRatioEt = nofine.make<TH1F>("EtRatio", "EtRatio;Et/3X3", 200, 0, 1.001);

    HitsEventNo = nofine.make<TH1I>("Hits", "Hits", 10.1, 0, 10);
    PrimNoFineTrigCanLum = nofine.make<TH2F>("Luminosity_Vrs_Trig_Candidates", "Luminosity_Vrs_Trig_Candidates; Luminosity ;Trigger_Candidates", 30, 1900, 7900, 11, 0, 11.);
    SeedPrimVrs3X3Fine = fine.make<TH2F>("Seed_Plus_EdgeVs3X3", "Seed_Plus_EdgeVs3X3; Seed Energy; 3X3", 100, 0, 100, 100, 0, 100);
    SeedPrimVrs3X3NoFine = nofine.make<TH2F>("Seed_Plus_EdgeVs3X3", "Seed_Plus_EdgeVs3X3; Seed Energy; 3X3", 100, 0, 100, 100, 0, 100);
    comp3X3VS5X5Fine = fine.make<TH2F>("3X3VS5X5F", "3X3VS5X5F; 3X3; 5X5", 100, 0, 100, 100, 0, 100);
    comp3X3VS5X5NoFine = nofine.make<TH2F>("3X3VS5X5F", "3X3VS5X5F; 3X3; 5X5", 100, 0, 100, 100, 0, 100);
    CentVs3X3Fine = fine.make<TH2F>("CentVs3X3", "CentVs3X3; cent; 3X3", 100, 0, 100, 100, 0, 100);
    CentVs3X3NoFine = nofine.make<TH2F>("CentVs3X3", "CentVs3X3; cent; 3X3", 100, 0, 100, 100, 0, 100);
    passed = 0;
    failed = 0;
    lumi = 0;
    //now do what ever initialization is needed
}

DataAnalizer::~DataAnalizer() {

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------

void
DataAnalizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    bool fine = false;
    //int test=0;
    //cout<<" ls number :"<<iEvent.luminosityBlock()<<endl;
    edm::Handle<reco::VertexCollection> reco_vertices;
    iEvent.getByLabel("offlinePrimaryVertices", reco_vertices);
    int NumbVert = reco_vertices->size();

    for(int q = 0; q < 2; q++)
    {
        int hit = 0;
        PrimDigiUse PCuts = PrimInfo(iEvent, iSetup, fine);
        double bottom = 0;
        double ratio = 0;
        double seed = 0;
        for(int ieta = 0; ieta < 10; ieta++)
        {
            for(int iphi = 1; iphi < 72; iphi += 2)
            {
                //if(!fine)cout<<"energy :"<<PCuts.MaxEtPos[ieta][iphi]<<endl;


                double T5X5 = 0;
                if(PCuts.LocMaxPos[ieta][iphi])
                {
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
                            if(abs(i) < 2 && abs(w) < 2)
                            {
                                if(PCuts.MaxEtPos[ietap][iphip] > seed && i != 0 && w != 0) seed = PCuts.MaxEtPos[ietap][iphip];
                                bottom += PCuts.MaxEtPos[ietap][iphip];
                            }
                            T5X5 += PCuts.MaxEtPos[ietap][iphip];
                        }
                    }

                    if(seed * 2 < PCuts.MaxEtPos[ieta][iphi])seed = 0;
                    ratio= PCuts.MaxEtPos[ieta][iphi]/bottom;
                    if(fine)
                    {
                        hit++;
                        //cout<<" hit number :"<<hit<<endl;
                        cout << "OKAY WE GOT THIS :" << PCuts.MaxEtPos[ieta][iphi] + seed << endl;
                        CentVs3X3Fine->Fill(PCuts.MaxEtPos[ieta][iphi], bottom);
                        SeedPrimVrs3X3Fine->Fill(PCuts.MaxEtPos[ieta][iphi] + seed, bottom);
                        PrimPassEt->Fill(PCuts.MaxEtPos[ieta][iphi]);
                        PrimPassRatioEt->Fill(ratio);
                        comp3X3VS5X5Fine->Fill(bottom, T5X5);

                    } else
                    {
                        hit++;
                        CentVs3X3NoFine->Fill(PCuts.MaxEtPos[ieta][iphi], bottom);
                        SeedPrimVrs3X3NoFine->Fill(PCuts.MaxEtPos[ieta][iphi] + seed, bottom);
                        PrimAllPassEt->Fill(PCuts.MaxEtPos[ieta][iphi]);
                        PrimAllPassRatioEt->Fill(ratio);
                        comp3X3VS5X5NoFine->Fill(bottom, T5X5);
                    }


                }






                ////////POS to Neg


                seed=0;
                if(PCuts.LocMaxNeg[ieta][iphi])
                {
                    cout << "Test energy :" << PCuts.MaxEtNeg[ieta][iphi] << endl;
                    double ratio = 0;
                    double bottom = 0;
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
                            if(abs(i) < 2 && abs(w) < 2)
                            {
                                if(PCuts.MaxEtNeg[ietap][iphip] > seed && i != 0 && w != 0)seed = PCuts.MaxEtNeg[ietap][iphip];
                                bottom += PCuts.MaxEtNeg[ietap][iphip];
                            }
                            T5X5 += PCuts.MaxEtNeg[ietap][iphip];
                            
                        }
                    }
                    if(seed * 2 < PCuts.MaxEtNeg[ieta][iphi])seed = 0;
                    ratio = double(PCuts.MaxEtNeg[ieta][iphi]) / (bottom);
                    if(fine)
                    {
                        hit++;

                        
                        PrimPassEt->Fill(PCuts.MaxEtNeg[ieta][iphi]);
                        PrimPassRatioEt->Fill(ratio);
                        comp3X3VS5X5Fine->Fill(bottom, T5X5);
                        CentVs3X3Fine->Fill(PCuts.MaxEtNeg[ieta][iphi], bottom);
                        SeedPrimVrs3X3Fine->Fill(seed+PCuts.MaxEtNeg[ieta][iphi], bottom);
                    } else
                    {
                        hit++;

                        CentVs3X3NoFine->Fill(PCuts.MaxEtNeg[ieta][iphi], bottom);
                        PrimAllPassEt->Fill(PCuts.MaxEtNeg[ieta][iphi]);
                        SeedPrimVrs3X3NoFine->Fill(PCuts.MaxEtNeg[ieta][iphi]+seed, bottom);
                        PrimAllPassRatioEt->Fill(ratio);
                        comp3X3VS5X5NoFine->Fill(bottom, T5X5);

                    }

                    ratio = double(PCuts.MaxEtNeg[ieta][iphi]) / (bottom);
                    if(PCuts.MaxEtNeg[ieta][iphi] == 0)
                    {
                        cout << "this should not happen" << endl;
                        cout << " our ieta is :-" << ieta + 30 << " and our iphi is :" << iphi << endl;
                        int k;
                        cin >> k;
                    }
                    if(seed * 2 < PCuts.MaxEtNeg[ieta][iphi])seed = 0;
                }


            }
        }

        if(fine)
        {
            NumofVetVrsTrigCanPrim->Fill(NumbVert, hit);
            HitsEventfine->Fill(hit);
            PrimFineTrigCanLum->Fill(lumi, hit);

        } else
        {
            HitsEventNo->Fill(hit);
            PrimNoFineTrigCanLum->Fill(lumi, hit);
            NumofVetVrsTrigCanNoPrim->Fill(NumbVert, hit);

        }
        hit = 0;
        fine = true;
    }
}


// ------------ method called once each job just before starting event loop  ------------

void
DataAnalizer::beginJob() { }

// ------------ method called once each job just after ending the event loop  ------------

void
DataAnalizer::endJob() { }

// ------------ method called when starting to processes a run  ------------
/*
void 
DataAnalizer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
 */

// ------------ method called when ending the processing of a run  ------------
/*
void 
DataAnalizer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
 */

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
DataAnalizer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
 */

// ------------ method called when ending the processing of a luminosity block  ------------

void
DataAnalizer::beginLuminosityBlock(edm::LuminosityBlock const& lumiBlock, edm::EventSetup const& es)
{
    edm::Handle<LumiSummary> lumisummary; //get the raw lumi data from edm::LumiSummary
    lumiBlock.getByLabel("lumiProducer", lumisummary);
    float instlumi = lumisummary->avgInsDelLumi();
    float correctedinstlumi = instlumi;
    //float recinstlumi = lumisummary->avgInsRecLumi();
    float corrfac = 1.;
    edm::ESHandle<LumiCorrectionParam> datahandle; //get LumiCorrectionParam object from event setup
    es.getData(datahandle);
    if(datahandle.isValid())
    {
        const LumiCorrectionParam* mydata = datahandle.product();
        //std::cout << "correctionparams " << *mydata << std::endl;
        corrfac = mydata->getCorrection(instlumi); //get lumi correction factor. Note: the corrfac is dependent of the average inst lumi value of the LS
    } else
    {
        std::cout << "no valid record found" << std::endl;
    }
    correctedinstlumi = instlumi*corrfac; //final lumi value=correction*raw lumi value
    //float correctedinstRecLumi = recinstlumi*corrfac; //Note: the correction factor is a function of raw inst delivered lumi. 
    //cout << "okay so this is the recinstlumi assuming no correction" << recinstlumi << endl;
    //cout << " and the instlumi is :" << instlumi << endl;
    //cout<<"correctedinstlumi :"<<correctedinstlumi;
    //cout << "correctedinstRecLumi :" << correctedinstRecLumi << endl;
    lumi = correctedinstlumi;
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------

void
DataAnalizer::fillDescriptions(edm::ConfigurationDescriptions & descriptions)
{
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

DataAnalizer::PrimDigiUse DataAnalizer::PrimInfo(const edm::Event& iEvent, const edm::EventSetup& eventSetup, bool UseFine)
{

    edm::ESHandle<HcalTPGCoder> inputCoder;
    eventSetup.get<HcalTPGRecord>().get(inputCoder);

    edm::ESHandle<CaloTPGTranscoder> outTranscoder;
    eventSetup.get<CaloTPGRecord>().get(outTranscoder);
    outTranscoder->setup(eventSetup, CaloTPGTranscoder::HcalTPG);
    double MinSeed = 10;
    bool featurebitPos[10][72];
    bool featurebitNeg[10][72];
    //double MinEt =0;
    struct PrimDigiUse CutInfo;
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
        if(abs(tp->id().ieta()) < 30 || abs(tp->id().ieta()) > 39)continue;

        double MaxEt = 0;

        for(int i = 0; i < tp->size(); i++)
        {
            //if((MaxEt < outTranscoder->hcaletValue(tp->id(), (*tp)[i])&&MinEt < outTranscoder->hcaletValue(tp->id(), (*tp)[i]))&&((*tp)[i].fineGrain() || !UseFine))
            //{
            //cout<<"okay hopfully this works we have a compression bit of :"<<(*tp)[i].compressedEt()<<endl;
            //cout<<"and our real Et is "<<outTranscoder->hcaletValue(tp->id(), (*tp)[i])<<endl;
            //  MaxEt = outTranscoder->hcaletValue(tp->id(), (*tp)[i]);
            // cout << "okay this should be non zero ptest"<<MaxEt << endl;



            //}
        }


        MaxEt = outTranscoder->hcaletValue(tp->id(), tp->SOI_compressedEt());
        //if((!tp->SOI_fineGrain())&&UseFine)MaxEt=0;
        if(tp->id().ieta() > 0)
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

            CutInfo.LocMaxPos[ietac][iphi] = (upleft && upcent && upright && right && left && botleft && botcent && botright && CutInfo.MaxEtPos[ietac][iphi] > MinSeed && featurebitPos[ietac][71]);
            //if(ietac==0&&CutInfo.LocMaxPos[ietac][iphi])cout<<" huh this is odd."<<endl;

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
            bool upleft = CutInfo.MaxEtNeg[ietac][iphi] > CutInfo.MaxEtNeg[ietac + 1][iphi - 2] || ietac == 9 || ietac == 10;
            bool upcent = CutInfo.MaxEtNeg[ietac][iphi] > CutInfo.MaxEtNeg[ietac + 1][iphi] || ietac == 10 || ietac == 9;
            bool upright = CutInfo.MaxEtNeg[ietac][iphi] > CutInfo.MaxEtNeg[ietac + 1][iphi + 2] || ietac == 10 || ietac == 9;
            bool left = CutInfo.MaxEtNeg[ietac][iphi] > CutInfo.MaxEtNeg[ietac][iphi - 2];
            bool right = CutInfo.MaxEtNeg[ietac][iphi] > CutInfo.MaxEtNeg[ietac][iphi + 2];
            bool botleft = CutInfo.MaxEtNeg[ietac][iphi] > CutInfo.MaxEtNeg[ietac - 1][iphi - 2] || ietac == 0 || ietac == 19;
            bool botcent = CutInfo.MaxEtNeg[ietac][iphi] > CutInfo.MaxEtNeg[ietac - 1][iphi] || ietac == 0 || ietac == 19;
            bool botright = CutInfo.MaxEtNeg[ietac][iphi] > CutInfo.MaxEtNeg[ietac - 1][iphi + 2] || ietac == 0 || ietac == 19;



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


//define this as a plug-in
DEFINE_FWK_MODULE(DataAnalizer);
#endif
