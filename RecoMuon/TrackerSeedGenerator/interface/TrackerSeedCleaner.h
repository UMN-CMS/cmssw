#ifndef RecoMuon_TrackerSeedCleaner_H
#define RecoMuon_TrackerSeedCleaner_H

/** \class TrackerSeedCleaner
 *  Seeds Cleaner based on direction
 *  $Date: 2008/03/03 15:36:12 $
 *  $Revision: 1.1 $
    \author A. Grelli -  Purdue University, Pavia University
 */ 

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TkTransientTrackingRecHitBuilder.h"
#include "RecoTracker/TkTrackingRegions/interface/RectangularEtaPhiTrackingRegion.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

class MuonServiceProxy;
class TSGFromL2Muon;
class MuonTrackingRegionBuilder;

namespace edm {class ParameterSet; class Event; class EventSetup;}

//              ---------------------
//              -- Class Interface --
//              ---------------------

class TrackerSeedCleaner{

public:

  typedef std::vector<TrajectorySeed> tkSeeds;
  /// constructor  
  TrackerSeedCleaner(const edm::ParameterSet& pset) : theProxyService(0),theEvent(0) {
                   builderName_ = pset.getParameter<std::string>("TTRHBuilder");
                   theBeamSpotTag = pset.getParameter<edm::InputTag>("beamSpot");
                   usePt_Cleaner = pset.getParameter<bool>("ptCleaner");

  }

  ///intizialization
  virtual void init(const MuonServiceProxy *service);
  
  /// destructor
  virtual ~TrackerSeedCleaner() {}
  /// the cleaner
  virtual void  clean(const reco::TrackRef& , const RectangularEtaPhiTrackingRegion& region, tkSeeds&);  
  /// setEvent
  virtual void setEvent(const edm::Event&);

private:

  const MuonServiceProxy * theProxyService;
  const edm::Event * theEvent;

  edm::InputTag theBeamSpotTag; //beam spot
  edm::Handle<reco::BeamSpot> bsHandle_;

  std::string builderName_;
  edm::ESHandle<TransientTrackingRecHitBuilder> theTTRHBuilder;
  bool usePt_Cleaner; 
};

#endif

