#ifndef StRPEvent_h
#define StRPEvent_h

//_____________________________________________________________________________
//    Class for making picoDst RP data 2019
//    Author: Truhlar Tomas
//_____________________________________________________________________________

class StUPCRpsCluster;
class StUPCRpsTrack;
class StUPCRpsTrackPoint;

class TClonesArray;
class TIterator;
class TParticle;

#include "TArrayI.h"

class StRPEvent
{

public:

  StRPEvent();
  StRPEvent(const StRPEvent &evt);
  virtual ~StRPEvent(); 

  void clearEvent(); 

  //setters
  StUPCRpsCluster *addCluster();
  StUPCRpsTrack *addTrack();
  StUPCRpsTrackPoint *addTrackPoint();


  //getters
  UInt_t getNumberOfClusters() const;
  StUPCRpsCluster *getCluster(Int_t iCluster) const;  
  TClonesArray *getClusters() const { return mClusters; }

  UInt_t getNumberOfTracks() const;
  StUPCRpsTrack *getTrack(Int_t iTrack) const;
  TClonesArray *getTracks() const { return mTracks; }

  UInt_t getNumberOfTrackPoints() const;
  StUPCRpsTrackPoint *getTrackPoint(Int_t iTrackPoint) const;
  TClonesArray *getTrackPoints() const { return mTrackPoints; }

// Setters, romanPotID == RomanPot ID, planeID == Plane Id, val == value of current setter
void setSiliconBunch(UChar_t val) 
{
  mSiliconBunch = val;
}
void setStatus(UInt_t romanPotID, UChar_t val) 
{
  if( romanPotID > mNumberOfRomanPots -1) return;
  mStatusRomanPot[romanPotID] = val;
}
void setNumberPlanes(UInt_t romanPotID, UShort_t val) 
{
  if( romanPotID > mNumberOfRomanPots -1) return;
  mNumberPlanes[romanPotID] = val;
}

void setNumberPlanesWithCluster(UInt_t romanPotID, UInt_t val) 
{
  if( romanPotID > mNumberOfRomanPots -1) return;
  mNumberPlanesWithCluster[romanPotID] = val;
}

void setAdc(UInt_t romanPotID, UInt_t adc0, UShort_t adc1) 
{
  if( romanPotID > mNumberOfRomanPots -1) return;

  mADC[romanPotID][0] = adc0;
  mADC[romanPotID][1] = adc1;
}
void setTac(UInt_t romanPotID, UInt_t tac0, UShort_t tac1)
{ 
  if( romanPotID > mNumberOfRomanPots -1) return;
  mTAC[romanPotID][0] = tac0;
  mTAC[romanPotID][1] = tac1;
}   
void setOffset(UInt_t romanPotID, UInt_t planeID, Float_t val) 
{
  if( romanPotID > mNumberOfRomanPots -1) return;
  if( planeID > mNumberOfPlanesInRP -1) return;

  mOffsetPlane[romanPotID][planeID] = val;
}
void setZ(UInt_t romanPotID, UInt_t planeID, Float_t val) 
{
  if( romanPotID > mNumberOfRomanPots -1) return;
  if( planeID > mNumberOfPlanesInRP -1) return;
  mzPlane[romanPotID][planeID] = val;
}
void setAngle(UInt_t romanPotID, UInt_t planeID, Float_t val) 
{
  if( romanPotID > mNumberOfRomanPots -1) return;
  if( planeID > mNumberOfPlanesInRP -1) return;
  mAnglePlane[romanPotID][planeID] = val;
}
void setOrientation(UInt_t romanPotID, UInt_t planeID, Short_t val) 
{
  if( romanPotID > mNumberOfRomanPots -1) return;
  if( planeID > mNumberOfPlanesInRP -1) return;
  mOrientationPlane[romanPotID][planeID] = val;
}
void setStatus(UInt_t romanPotID, UInt_t planeID, UChar_t val) 
{
  if( romanPotID > mNumberOfRomanPots -1) return;
  if( planeID > mNumberOfPlanesInRP -1) return;
  mStatusPlane[romanPotID][planeID] = val;
}

// Getters, romanPotID == RomanPot ID, planeID == Plane Id
UChar_t siliconBunch() const 
{ 
  return mSiliconBunch;
}

UShort_t numberOfPlanes(UInt_t romanPotID) const 
{ 
  if( romanPotID > mNumberOfRomanPots -1) return 0;
  return mNumberPlanes[romanPotID];
}

UShort_t numberOfPlanesWithClusters(UInt_t romanPotID) const 
{ 
  if( romanPotID > mNumberOfRomanPots -1) return 0;
  return mNumberPlanesWithCluster[romanPotID];
}

UChar_t status(UInt_t romanPotID) const 
{ 
  if( romanPotID > mNumberOfRomanPots -1) return 0;
  return mStatusRomanPot[romanPotID];
}
UShort_t adc(UInt_t romanPotID, UInt_t pmtID) const // return adc value from first or second pmt
{ 
  if( romanPotID > mNumberOfRomanPots -1) return 0;
  if( pmtID > mNumberOfPmtsInRp -1) return 0;
  return mADC[romanPotID][pmtID];
}
UShort_t tac(UInt_t romanPotID, UInt_t pmtID) const  // return tac value from first or second pmt
{ 
  if( romanPotID > mNumberOfRomanPots -1) return 0;
  if( pmtID > mNumberOfPmtsInRp -1) return 0;
  return mTAC[romanPotID][pmtID];
}
Float_t offset(UInt_t romanPotID, UInt_t planeID) const 
{ 
  if( romanPotID > mNumberOfRomanPots -1) return 0;
  if( planeID > mNumberOfPlanesInRP -1) return 0;
  return mOffsetPlane[romanPotID][planeID];
}
Float_t z(UInt_t romanPotID, UInt_t planeID) const 
{ 
  if( romanPotID > mNumberOfRomanPots -1) return 0;
  if( planeID > mNumberOfPlanesInRP -1) return 0;
  return mzPlane[romanPotID][planeID];
}
Float_t angle(UInt_t romanPotID, UInt_t planeID) const 
{ 
  if( romanPotID > mNumberOfRomanPots -1) return 0;
  if( planeID > mNumberOfPlanesInRP -1) return 0;
  return mAnglePlane[romanPotID][planeID];
}
Short_t  orientation(UInt_t romanPotID, UInt_t planeID) const 
{ 
  if( romanPotID > mNumberOfRomanPots -1) return 0;
  if( planeID > mNumberOfPlanesInRP -1) return 0;
  return mOrientationPlane[romanPotID][planeID];
}

UChar_t status(UInt_t romanPotID, UInt_t planeID) const 
{ 
  if( romanPotID > mNumberOfRomanPots -1) return 0;
  if( planeID > mNumberOfPlanesInRP -1) return 0;
  return mStatusPlane[romanPotID][planeID];
}


private:

  StRPEvent &operator=(const StRPEvent &o); // not implemented

  static const Int_t mNumberOfRomanPots = 8;
  static const Int_t mNumberOfPlanesInRP = 4;
  static const Int_t mNumberOfPmtsInRp = 2;

  Float_t  mOffsetPlane[mNumberOfRomanPots][mNumberOfPlanesInRP];
  Float_t  mzPlane[mNumberOfRomanPots][mNumberOfPlanesInRP];
  Float_t  mAnglePlane[mNumberOfRomanPots][mNumberOfPlanesInRP];

  Short_t  mOrientationPlane[mNumberOfRomanPots][mNumberOfPlanesInRP];
  UShort_t mNumberPlanesWithCluster[mNumberOfRomanPots]; 
  UShort_t mADC[mNumberOfRomanPots][mNumberOfPmtsInRp]; // ADC in 2 pmt
  UShort_t mTAC[mNumberOfRomanPots][mNumberOfPmtsInRp]; // TAC in 2 pmt
  UShort_t mNumberPlanes[mNumberOfRomanPots]; 

  UChar_t  mStatusPlane[mNumberOfRomanPots][mNumberOfPlanesInRP];  
  UChar_t  mSiliconBunch;
  UChar_t  mStatusRomanPot[mNumberOfRomanPots]; 

  static TClonesArray *mgClusters; // array of upc tracks
  TClonesArray *mClusters; //-> array of upc tracks
  UInt_t mNClusters; //! number of upc tracks in event, local use when filling

  static TClonesArray *mgTracks; // array of upc tracks
  TClonesArray *mTracks; //-> array of upc tracks
  UInt_t mNTracks; //! number of upc tracks in event, local use when filling

  static TClonesArray *mgTrackPoints; // array of upc tracks
  TClonesArray *mTrackPoints; //-> array of upc tracks
  UInt_t mNTrackPoints; //! number of upc tracks in event, local use when filling

  ClassDef(StRPEvent, 2);
};

#endif













