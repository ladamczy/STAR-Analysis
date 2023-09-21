//_____________________________________________________________________________
//    Class for making picoDst RP data 2019
//    Author: Truhlar Tomas
//_____________________________________________________________________________

#include "StUPCRpsTrackPoint.h"
#include "StUPCRpsTrack.h"
#include "StRPEvent.h"
#include <cmath>

ClassImp(StUPCRpsTrack)

StUPCRpsTrack::StUPCRpsTrack(): TObject(),
  mEvt(0x0)
{
    mFirstTrackPointId = -1;
    mSecondTrackPointId = -1;
    mBranch = -1;
    mType = rpsUndefined;
}

StUPCRpsTrack::~StUPCRpsTrack() { /* no op */ }

void StUPCRpsTrack::Clear(Option_t *)
{
  //clear rp track object, overridden from TObject
  //for use in TClonesArray

    mFirstTrackPointId = -1;
    mSecondTrackPointId = -1;
    mBranch = -1;
    mType = rpsUndefined;

}//Clear

StUPCRpsTrack& StUPCRpsTrack::operator=(const StUPCRpsTrack& track) {
    if (this != &track) {
        mP = track.pVec();
        mType = track.type();
    }
    return *this;
}

StUPCRpsTrackPoint* StUPCRpsTrack::getTrackPoint(unsigned int station) const {
    if(!mEvt) return 0x0;

    if(station==0){
        return mFirstTrackPointId == -1 ? nullptr : mEvt->getTrackPoint(mFirstTrackPointId);
    }else if(station==1){
        return mSecondTrackPointId == -1 ? nullptr : mEvt->getTrackPoint(mSecondTrackPointId);
    }else{
        return 0x0;
    }
}

unsigned int StUPCRpsTrack::planesUsed() const {
    unsigned int nPlanes = 0;
    for(unsigned int i=0; i<mNumberOfStationsInBranch; ++i)
        nPlanes += getTrackPoint(i) ? getTrackPoint(i)->planesUsed() : 0;
    return nPlanes;
}

double StUPCRpsTrack::thetaRp(unsigned int coordinate) const {
    if(coordinate>rpsAngleTheta) return 0.0;
    if(mType==rpsLocal) return theta(coordinate);
    TVector3 deltaVector = getTrackPoint(1)->positionVec() - getTrackPoint(0)->positionVec();
    return atan((coordinate<rpsAngleTheta ? deltaVector[coordinate] : deltaVector.Perp())/abs(deltaVector.z()));
}

double StUPCRpsTrack::phiRp() const{
    if(mType==rpsLocal) return phi();
    TVector3 deltaVector = getTrackPoint(1)->positionVec() - getTrackPoint(0)->positionVec();
    return deltaVector.Phi();
}

double StUPCRpsTrack::time() const{
    double sumTime=0.0;
    unsigned int numberOfPmtsWithSignal=0;
    for(unsigned int i=0; i<mNumberOfStationsInBranch; ++i){
        if(getTrackPoint(i))
            for(int j=0; j< getTrackPoint(i)->mNumberOfPmtsInRp; ++j){
                if(getTrackPoint(i)->time(j)>0){
                    sumTime += getTrackPoint(i)->time(j);
                    ++numberOfPmtsWithSignal;
                }
            }
    }
    return numberOfPmtsWithSignal>0 ? sumTime/numberOfPmtsWithSignal : -1;
}

double StUPCRpsTrack::theta(unsigned int coordinate) const
{
    return coordinate < mNumberOfAngleTypes ? atan((coordinate<rpsAngleTheta ? mP[coordinate] : mP.Perp())/abs(mP.z())) : 0.0;
}
