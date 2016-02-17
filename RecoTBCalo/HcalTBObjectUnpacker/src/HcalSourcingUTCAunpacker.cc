#include "RecoTBCalo/HcalTBObjectUnpacker/interface/HcalSourcingUTCAunpacker.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include <iostream>
#include <string>


/// Per Event Header Structure
struct eventHeader {
  uint32_t cdf0;
  uint32_t cdf1;
  uint32_t cdf2;
  uint32_t cdf3;
  uint32_t h0;
  uint32_t h1;
  uint32_t h2;
  uint32_t h3;
};

  
void HcalSourcingUTCAunpacker::unpack(const FEDRawData&  raw, const HcalElectronicsMap emap, std::auto_ptr<HcalUHTRhistogramDigiCollection>& histoDigiCollection, bool DEBUG) const {
  if(DEBUG) std::cout << "Unpacker Time!" << std::endl; 

  if (raw.size()<4*38) {
//    throw cms::Exception("Missing Data") << "Less than 1 histogram in event";
    std::cout << "Warning: raw size in bytes: " << raw.size() << std::endl;
  }
  const uint32_t* pData = (const uint32_t*) raw.data(); 
  if(DEBUG) {
  int nwords=raw.size()/4;
  for (int iw=0; iw<nwords; iw++)
    printf("%04d %04x\n",iw,pData[iw]);
  }
  const struct eventHeader* eh =
    (const struct eventHeader*)(raw.data());
  
//  if (raw.size()<sizeof(xdaqSourcePositionDataFormat)) {
//    throw cms::Exception("DataFormatError","Fragment too small");
//  }
//Read event header
  uint32_t numHistos  = ((eh->h3)>>16)&0xFFFF;
  if(DEBUG) std::cout << "Number of Histograms: " << numHistos << std::endl;
  uint32_t numBins    = ((eh->h3)>>1)&0x0000FFFE; //includes overflow and header word
  if(DEBUG) std::cout << "Bins per Histogram: " << numBins << std::endl;
  bool sepCapIds = eh->h3&0x00000001;
  if(DEBUG) std::cout << "Separate CapIds: " << sepCapIds << std::endl;

  histoDigiCollection.reset(new HcalUHTRhistogramDigiCollection(numBins+1, sepCapIds));
//Set histogram word pointer to first histogram    
  uint32_t crate   = 0;
  uint32_t slot    = 0;
  uint32_t fiber   = 0;
  uint32_t channel = 0;
  uint32_t cap     = 0;
//Loop over data
  pData+=8;
  for (unsigned int iHist = 0; iHist<numHistos; iHist++) {
    if(DEBUG) std::cout << "Histogram " << iHist <<" header: "<< *pData << std::endl;
    crate   = ((*pData)>>16)&0x00FF;
    if(DEBUG) std::cout << "Crate: " << crate << std::endl;
    slot    = ((*pData)>>12)&0x0000F;
    if(DEBUG) std::cout << "Slot: " << slot << std::endl;
    fiber   = (*pData>>7)&0x1F;
    if(DEBUG) std::cout << "Fiber: " << fiber << " "<< (*pData&0x00000F80)<<std::endl;
    channel = (*pData>>2)&0x1F;
    if(DEBUG) std::cout << "Channel: " << channel << std::endl;
    cap     = *pData&0x00000003;
    if(DEBUG) std::cout << "CapId: " << cap << " from: " << 0x3 << std::endl;
    HcalElectronicsId eid(crate, slot, fiber, channel, false);
  //  eid.setHTR(htr_cr,htr_slot,htr_tb);
    DetId did=emap.lookup(eid);
    if (did.null() || did.det()!=DetId::Hcal || did.subdetId()==0) {
      if (unknownIds_.find(eid)==unknownIds_.end()) {
        edm::LogWarning("HCAL") << "HcalHistogramUnpacker: No match found for electronics id :" << eid;
      }
      if(iHist<(numHistos-1)) {
        pData+=(numBins+2);
      }
      continue;
    }
    if(DEBUG) std::cout << "Det Id: " << ((HcalDetId)did) << std::endl;
    HcalUHTRhistogramDigiMutable digi = histoDigiCollection->addHistogram( (HcalDetId)did );
    for(unsigned int iBin = 0; iBin<numBins+1; iBin++) {
      digi.fillBin(cap, iBin, pData[iBin+1]);
      if(DEBUG) std::cout << "CapId: " << cap << "Bin: " << iBin << "Val: " << pData[iBin+1] << std::endl;
        
    }
    if(iHist<(numHistos-1)) {
      pData+=(numBins+2);
    }
  }
  if(DEBUG) std::cout << "DONE" << std::endl;
}
