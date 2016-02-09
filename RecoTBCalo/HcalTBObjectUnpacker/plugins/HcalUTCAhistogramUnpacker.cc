using namespace std;

#include "DataFormats/HcalDigi/interface/HcalHistogramDigi.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "RecoTBCalo/HcalTBObjectUnpacker/plugins/HcalUTCAhistogramUnpacker.h"
#include <iostream>


  HcalUTCAhistogramUnpacker::HcalUTCAhistogramUnpacker(edm::ParameterSet const& conf)
  {

    tok_raw_ = consumes<FEDRawDataCollection>(conf.getParameter<edm::InputTag>("fedRawDataCollectionTag"));
    produces<HcalHistogramDigiCollection>();
  }

  // Virtual destructor needed.
  HcalUTCAhistogramUnpacker::~HcalUTCAhistogramUnpacker() { }  

  // Functions that gets called by framework every event
  void HcalUTCAhistogramUnpacker::produce(edm::Event& e, const edm::EventSetup& es)
  {
    edm::Handle<FEDRawDataCollection> rawraw;  
    edm::ESHandle<HcalElectronicsMap>   item;
    edm::ESHandle<HcalDbService> pSetup;

    e.getByToken(tok_raw_, rawraw);          
    es.get<HcalDbRecord>().get(pSetup);
    es.get<HcalElectronicsMapRcd>().get(item);

    const HcalElectronicsMap* readoutMap = item.product();
    std::auto_ptr<HcalHistogramDigiCollection> hd(new HcalHistogramDigiCollection);

    const FEDRawData& fed = rawraw->FEDData(18);

    histoUnpacker_.unpack(fed, *readoutMap, hd);

    e.put(hd);
  }

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(HcalUTCAhistogramUnpacker);
