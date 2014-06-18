/* 
 * File:   FindSquare.h
 * Author: lesko
 *
 * Created on May 19, 2014, 9:42 AM
 */
#include "DataFormats/HcalDigi/interface/HcalTriggerPrimitiveDigi.h"
#include "Geometry/HcalTowerAlgo/interface/HcalTrigTowerGeometry.h"
#include "DataFormats/HcalDetId/interface/HcalTrigTowerDetId.h"


#ifndef FINDSQUARE_H
#define	FINDSQUARE_H
//note this only works for HF at this point in time. may turn this into a child class for hf in future. 
class FindSquare {
public:
    FindSquare(const edm::Event& iEvent);
   // FindSquare(const FindSquare& orig);
    virtual ~FindSquare();
    HcalTrigTowerDetId *UpperLeft();
    HcalTrigTowerDetId *UpperCenter();
    HcalTrigTowerDetId *UpperRight();
    HcalTrigTowerDetId *MidLeft();
    HcalTrigTowerDetId *MidRight();
    HcalTrigTowerDetId *BottomLeft();
    HcalTrigTowerDetId *BottomCenter();
    HcalTrigTowerDetId *BottomRight();
    
private:
const edm::Event& iEvent_;
bool basic_;
struct PrimDigiUse{
    bool LocMax[18][72];
    int MaxEt[18][72];
};
};

#endif	/* FINDSQUARE_H */

