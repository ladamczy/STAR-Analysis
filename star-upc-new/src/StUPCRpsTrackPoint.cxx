//_____________________________________________________________________________
//    Class for making picoDst RP data 2019
//    Author: Truhlar Tomas
//_____________________________________________________________________________

#include "StUPCRpsTrackPoint.h"

ClassImp(StUPCRpsTrackPoint)

StUPCRpsTrackPoint::StUPCRpsTrackPoint()
{
    mRpId = -1;
    for (unsigned int i=0; i<mNumberOfPlanesInRp; ++i) mClusterId[i] = -1;
    for (unsigned int i=0; i<mNumberOfPmtsInRp; ++i) mTime[i] = -1;
    mQuality = rpsNotSet;
}

StUPCRpsTrackPoint::~StUPCRpsTrackPoint() { /* no op */ }

void StUPCRpsTrackPoint::Clear(Option_t *)
{
  //clear rp trackPOINT object, overridden from TObject
  //for use in TClonesArray

    mRpId = -1;
    for (unsigned int i=0; i<mNumberOfPlanesInRp; ++i) mClusterId[i] = -1;
    for (unsigned int i=0; i<mNumberOfPmtsInRp; ++i) mTime[i] = -1;
    mQuality = rpsNotSet;

}//Clear

StUPCRpsTrackPoint& StUPCRpsTrackPoint::operator=(const StUPCRpsTrackPoint& trackPoint)
{
    if (this != &trackPoint) {
        mPosition = trackPoint.positionVec();
        mRpId = trackPoint.rpId();
        for (unsigned int i=0; i<mNumberOfPlanesInRp; ++i ) mClusterId[i] = trackPoint.clusterId(i);
        for (unsigned int i=0; i<mNumberOfPmtsInRp; ++i ) mTime[i] = trackPoint.time(i);
        mQuality = trackPoint.quality();
    }
    return *this;
}

Int_t StUPCRpsTrackPoint::planesUsed() const
{
    unsigned int nPlanesUsed = 0;
    for(unsigned int i=0; i<mNumberOfPlanesInRp; ++i)
        if (mClusterId[i]>-1) ++nPlanesUsed;
    return nPlanesUsed;
}