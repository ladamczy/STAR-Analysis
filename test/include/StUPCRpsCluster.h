#ifndef StUPCRpsCluster_hh
#define StUPCRpsCluster_hh

//_____________________________________________________________________________
//    Class for making picoDst RP data 2019
//    Author: Truhlar Tomas
//_____________________________________________________________________________


#include <iostream> 
#include "TObject.h" 


class StUPCRpsCluster : public TObject {
public:
    StUPCRpsCluster();
    ~StUPCRpsCluster() { /* noop */ };

Float_t position() const { return mPosition; }

Float_t positionRMS() const { return mPositionRMS; }

Short_t  length() const { return mLength; }

Float_t energy() const { return mEnergy; }

Float_t xy() const { return mXY; }

UChar_t quality() const { return mQuality; }

unsigned int romanPotId() const { return mRomanPotId; }

unsigned int planeId() const { return mPlaneId; }

void setPosition(Float_t val) { mPosition = val; }

void setPositionRMS(Float_t val) { mPositionRMS = val; }

void setLength(Short_t val) { mLength = val; }

void setEnergy(Float_t val) { mEnergy = val; }

void setXY(Float_t val) { mXY = val; }

void setQuality(UChar_t val) { mQuality = val; }

void setPlaneId(UChar_t val) { mPlaneId = val; }

void setRomanPotId(UChar_t val) { mRomanPotId = val; };


 private:
    Float_t       mPosition;
    Float_t       mPositionRMS;    
    Float_t       mEnergy;
    Float_t       mXY;
    Short_t        mLength;
    UChar_t        mQuality;
    UChar_t        mPlaneId;
    UChar_t        mRomanPotId;

    ClassDef(StUPCRpsCluster,2)
};

#endif
