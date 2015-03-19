#ifndef HcalTrigTowerGeometry_h
#define HcalTrigTowerGeometry_h

#include "Geometry/CaloTopology/interface/HcalTopology.h"
#include "DataFormats/DetId/interface/DetId.h"
#include <vector>
class HcalTrigTowerDetId;
class HcalDetId;

class HcalTrigTowerGeometry {
private:
  // The versions
  static const int RCT_VERSION = 0;
  static const int HF_1X1_VERSION = 1;

public:

  HcalTrigTowerGeometry( const HcalTopology* topology );

  /// the mapping to and from DetIds
  std::vector<HcalTrigTowerDetId> towerIds(const HcalDetId & cellId) const;
  std::vector<HcalDetId> detIds(const HcalTrigTowerDetId &) const;

  void setupHFTowers(bool enableRCT, bool enable1x1) {
    useRCT_=enableRCT;
    use1x1_=enable1x1;
  }

  int firstHFTower(int version) const {return (version==RCT_VERSION)?(29):(30);}

  // where this tower begins and ends in eta
  void towerEtaBounds(int ieta, int version, double & eta1, double & eta2) const;

  // number of towers (version dependent)
  int nTowers(int version) const {return (version==RCT_VERSION)?(32):(41);}

  // get the topology pointer
  const HcalTopology& topology() const { return *theTopology; }

  // Get the useRCT and use1x1 values
  bool useRCT() const { return useRCT_; }
  bool use1x1() const { return use1x1_; }

  /// the number of phi bins in this eta ring
  int nPhiBins(int ieta, int version) const {
    int nPhiBinsHF = ( 18 );
    return (abs(ieta) < firstHFTower(version)) ? 72 : nPhiBinsHF;
  }

  /// the number of HF eta rings in this trigger tower
  /// ieta starts at firstHFTower()
  int hfTowerEtaSize(int ieta) const;

  /// since the towers are irregular in eta in HF
  int firstHFRingInTower(int ietaTower) const;

private:

  const HcalTopology* theTopology;
  bool useRCT_;
  bool use1x1_;

};

#endif

