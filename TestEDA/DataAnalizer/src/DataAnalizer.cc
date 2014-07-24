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
    TH1 *PrimPassEt, *PrimPassRatioEt;
    TH1 *HitsEventfine, *HitsEventNo;
    TH1 *PrimAllPassEt, *PrimAllPassRatioEt;
    TH2 *PrimPassEtvrsiEta, *PrimPassEtvrsiPhi;

    struct PrimDigiUse
    {
        bool LocMaxPos[50][73];
        double MaxEtPos[50][73];
        bool LocMaxNeg[50][73];
        double MaxEtNeg[50][73];
    };

    PrimDigiUse PrimInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup, bool UseFine);

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

DataAnalizer::DataAnalizer(const edm::ParameterSet& iConfig)
{
    edm::Service<TFileService> fs;
    cout<<"test 1"<<endl;
    TFileDirectory fine = fs->mkdir("finebit");
    TFileDirectory nofine = fs->mkdir("no_fine");
    PrimPassEt = fine.make<TH1F>("Et", "Et;Et;Count", 200, 0, 100);
    PrimPassRatioEt = fine.make<TH1F>("EtRatio", "EtRatio;Et/3X3", 200, 0, 1.001);
    PrimPassEtvrsiEta = fine.make<TH2F>("EtaVrsEt", "EtaVrsEt;Eta;Et", 10, 29.9, 39.9, 120, 0, 60);
    PrimPassEtvrsiPhi = fine.make<TH2F>("PhiVrsEt", "PhiVrsEt;Phi;Et", 36, 0, 72, 120, 0, 60);
    HitsEventfine = fine.make<TH1I>("Hits", "Hits", 100, 0, 100);

    PrimAllPassEt = nofine.make<TH1F>("Et", "Et;Et;Count", 200, 0, 100);
    PrimAllPassRatioEt = nofine.make<TH1F>("EtRatio", "EtRatio;Et/3X3", 200, 0, 1.001);
    HitsEventNo = nofine.make<TH1I>("Hits", "Hits", 100, 0, 100);
    passed = 0;
    failed = 0;
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

    for(int q = 0; q < 2; q++)
    {
        int hit = 0;
        PrimDigiUse PCuts = PrimInfo(iEvent, iSetup, fine);
        double bottom = 0;
        double ratio = 0;
        for(int ieta = 0; ieta < 10; ieta++)
        {
            for(int iphi = 1; iphi < 72; iphi += 2)
            {
                if(PCuts.LocMaxPos[ieta][iphi])
                {
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

                    if(fine)
                    {

                        PrimPassEtvrsiEta->Fill(30 + ieta, PCuts.MaxEtPos[ieta][iphi]);
                        PrimPassEtvrsiPhi->Fill(iphi, PCuts.MaxEtPos[ieta][iphi]);
                        PrimPassEtvrsiEta->Fill(30 + ieta, PCuts.MaxEtNeg[ieta][iphi]);
                        PrimPassEtvrsiPhi->Fill(iphi, PCuts.MaxEtNeg[ieta][iphi]);
                    } else
                    {
                        PrimAllPassEt->Fill(PCuts.MaxEtPos[ieta][iphi]);
                        PrimAllPassRatioEt->Fill(ratio);
                    }
                    hit++;
                    if(PCuts.MaxEtPos[ieta][iphi] == 0)
                    {
                        cout << "this should not happen" << endl;
                        cout << " our ieta is :" << ieta + 30 << " and our iphi is :" << iphi << endl;
                        int k;
                        cin >> k;
                    }
                }
                if(PCuts.LocMaxNeg[ieta][iphi])
                {
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
                    if(PCuts.MaxEtNeg[ieta][iphi] == 0)
                    {
                        cout << "this should not happen" << endl;
                        cout << " our ieta is :-" << ieta + 30 << " and our iphi is :" << iphi << endl;
                        int k;
                        cin >> k;
                    }
                    if(fine)
                    {
                        PrimPassEt->Fill(PCuts.MaxEtNeg[ieta][iphi]);
                        PrimPassRatioEt->Fill(ratio);
                        PrimPassEtvrsiEta->Fill(30 + ieta, PCuts.MaxEtNeg[ieta][iphi]);
                        PrimPassEtvrsiPhi->Fill(iphi, PCuts.MaxEtNeg[ieta][iphi]);
                    } else
                    {
                        PrimAllPassEt->Fill(PCuts.MaxEtNeg[ieta][iphi]);
                        PrimAllPassRatioEt->Fill(ratio);
                    }
                    hit++;
                }
            }
        }
        if(fine)HitsEventfine->Fill(hit);
        else HitsEventNo->Fill(hit);
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
DataAnalizer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
 */

// ------------ method called when ending the processing of a luminosity block  ------------

/*void DataAnalizer::endLuminosityBlock(edm::LuminosityBlock const& lumiBlock, edm::EventSetup const& es)
{
    cout << " I AM IN RUN NUMBER " << lumiBlock.run() << " LS NUMBER " << lumiBlock.luminosityBlock() << std::endl;
    edm::eventsetup::EventSetupRecordKey recordKey(edm::eventsetup::EventSetupRecordKey::TypeTag::findType("LumiCorrectionParamRcd"));
    if(recordKey.type() == edm::eventsetup::EventSetupRecordKey::TypeTag())
    {
        std::cout << "Record \"LumiCorrectionParamRcd" << "\" does not exist " << std::endl;
    }
    try
    {
        cout<<"test 0"<<endl;
        edm::Handle<LumiSummary> lumisummary;
        cout<<"test 1"<<endl;
        lumiBlock.getByLabel("lumiProducer", lumisummary);
        float instlumi = lumisummary->avgInsDelLumi();
        std::cout << "raw data tag " << lumisummary->lumiVersion() << std::endl;
        ;
        cout<<"test 2"<<endl;
        float correctedinstlumi = instlumi;
        float recinstlumi = lumisummary->avgInsRecLumi();
        cout<<"test 3"<<endl;
        float corrfac = 1.;
        edm::ESHandle<LumiCorrectionParam> datahandle;
        es.getData(datahandle);
        cout<<"test 4"<<endl;
        if(datahandle.isValid())
        {
            const LumiCorrectionParam* mydata = datahandle.product();
            std::cout << "correctionparams " << *mydata << std::endl;
            corrfac = mydata->getCorrection(instlumi);
        } else
        {
            std::cout << "no valid record found" << std::endl;
        }
        correctedinstlumi = instlumi*corrfac;
        std::cout << "correctedinstlumi " << correctedinstlumi << std::endl;
        float correctedinstRecLumi = recinstlumi*corrfac;
        std::cout << "corrected rec instlumi " << correctedinstRecLumi << std::endl;
    } catch(const edm::eventsetup::NoRecordException<LumiCorrectionParamRcd>& er)
    {
        std::cout << "no data found" << std::endl;
    } catch(const cms::Exception& ee)
    {
        std::cout << ee.what() << std::endl;
    }
}
*/

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
            if(MaxEt < outTranscoder->hcaletValue(tp->id(), (*tp)[i])&&((*tp)[i].fineGrain() || !UseFine))
            {
                //cout<<"okay hopfully this works we have a compression bit of :"<<(*tp)[i].compressedEt()<<endl;
                //cout<<"and our real Et is "<<outTranscoder->hcaletValue(tp->id(), (*tp)[i])<<endl;
                MaxEt = outTranscoder->hcaletValue(tp->id(), (*tp)[i]);
                // cout << "okay this should be non zero ptest"<<MaxEt << endl;



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
        CutInfo.LocMaxPos[ietac][1] = (upleft && upcent && upright && right && left && botleft && botcent && botright && CutInfo.MaxEtPos[ietac][1] > 5);

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
            bool upleft = CutInfo.MaxEtNeg[ietac][iphi] > CutInfo.MaxEtNeg[ietac + 1][iphi - 2] || ietac == 9 || ietac == 10;
            bool upcent = CutInfo.MaxEtNeg[ietac][iphi] > CutInfo.MaxEtNeg[ietac + 1][iphi] || ietac == 10 || ietac == 9;
            bool upright = CutInfo.MaxEtNeg[ietac][iphi] > CutInfo.MaxEtNeg[ietac + 1][iphi + 2] || ietac == 10 || ietac == 9;
            bool left = CutInfo.MaxEtNeg[ietac][iphi] > CutInfo.MaxEtNeg[ietac][iphi - 2];
            bool right = CutInfo.MaxEtNeg[ietac][iphi] > CutInfo.MaxEtNeg[ietac][iphi + 2];
            bool botleft = CutInfo.MaxEtNeg[ietac][iphi] > CutInfo.MaxEtNeg[ietac - 1][iphi - 2] || ietac == 0 || ietac == 19;
            bool botcent = CutInfo.MaxEtNeg[ietac][iphi] > CutInfo.MaxEtNeg[ietac - 1][iphi] || ietac == 0 || ietac == 19;
            bool botright = CutInfo.MaxEtNeg[ietac][iphi] > CutInfo.MaxEtNeg[ietac - 1][iphi + 2] || ietac == 0 || ietac == 19;



            CutInfo.LocMaxNeg[ietac][iphi] = (upleft && upcent && upright && right && left && botleft && botcent && botright && CutInfo.MaxEtNeg[ietac][iphi] > 5);

        }
        bool upleft = CutInfo.MaxEtNeg[ietac][1] > CutInfo.MaxEtNeg[ietac + 1][71] || ietac == 9;
        bool upcent = CutInfo.MaxEtNeg[ietac][1] > CutInfo.MaxEtNeg[ietac + 1][1] || ietac == 9;
        bool upright = CutInfo.MaxEtNeg[ietac][1] > CutInfo.MaxEtNeg[ietac + 1][3] || ietac == 9;
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


//define this as a plug-in
DEFINE_FWK_MODULE(DataAnalizer);
#endif