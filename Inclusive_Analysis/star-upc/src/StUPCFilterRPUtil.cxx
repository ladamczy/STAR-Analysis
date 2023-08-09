//_____________________________________________________________________________
//    Utility class for making picoDst RP data 2019
//    Author: Truhlar Tomas
//_____________________________________________________________________________


//c++ headers
#include "string.h"
#include <vector>
#include <iostream>

//root headers
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TH1I.h"

//StRoot headers
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuRpsCollection.h"
#include "StMuDSTMaker/COMMON/StMuRpsTrackPoint.h"
#include "StMuDSTMaker/COMMON/StMuRpsTrack.h"


//local headers
#include "StRPEvent.h"
#include "StUPCRpsCluster.h"
#include "StUPCRpsTrackPoint.h"
#include "StUPCRpsTrack.h"
#include "StUPCFilterRPUtil.h"
#include "StUPCFilterMaker.h"


//_____________________________________________________________________________
StUPCFilterRPUtil::StUPCFilterRPUtil() {
  //constructor

}//StUPCFilterRPUtil


void StUPCFilterRPUtil::clear() {

  //clear hits and clusters vectors
  mMuTrackPoints.clear();

}//clear

//_____________________________________________________________________________
void StUPCFilterRPUtil::processEvent(StRPEvent *rpEvt, StMuDst *mMuDst) {


  StMuRpsCollection *collection = mMuDst->RpsCollection();
  if( !collection ){
    cout<< "StUPCFilterRPUtil::processEvent() no muDst RpsCollection input" << endl;
    return;
  }
  rpEvt->setSiliconBunch(collection->siliconBunch());

  for(UInt_t iRomanPotId = 0; iRomanPotId < collection->numberOfRomanPots(); ++iRomanPotId){

    rpEvt->setStatus(iRomanPotId, collection->status(iRomanPotId));
    rpEvt->setAdc(iRomanPotId, collection->adc(iRomanPotId, 0), collection->adc(iRomanPotId, 1));
    rpEvt->setTac(iRomanPotId, collection->tac(iRomanPotId, 0), collection->tac(iRomanPotId, 1)); 
    rpEvt->setNumberPlanes(iRomanPotId, collection->numberOfPlanes());
    rpEvt->setNumberPlanesWithCluster(iRomanPotId, collection->numberOfPlanesWithClusters(iRomanPotId));
    for(UInt_t iPlaneId=0; iPlaneId < collection->numberOfPlanes(); ++iPlaneId){
      rpEvt->setOffset(iRomanPotId, iPlaneId, collection->offsetPlane(iRomanPotId, iPlaneId));
      rpEvt->setZ(iRomanPotId, iPlaneId, collection->zPlane(iRomanPotId, iPlaneId));  
      rpEvt->setAngle(iRomanPotId, iPlaneId, collection->anglePlane(iRomanPotId, iPlaneId));  
      rpEvt->setOrientation(iRomanPotId, iPlaneId, collection->orientationPlane(iRomanPotId, iPlaneId));  
      rpEvt->setStatus(iRomanPotId, iPlaneId, collection->statusPlane(iRomanPotId, iPlaneId));
      for(Int_t iCluster=0; iCluster < collection->numberOfClusters(iRomanPotId, iPlaneId); ++iCluster){
        StUPCRpsCluster *rpCluster = rpEvt->addCluster();
        rpCluster->setPosition(collection->positionCluster(iRomanPotId, iPlaneId, iCluster));
        rpCluster->setPositionRMS(collection->positionRMSCluster(iRomanPotId, iPlaneId, iCluster));
        rpCluster->setLength(collection->lengthCluster(iRomanPotId, iPlaneId, iCluster)); 
        rpCluster->setEnergy(collection->energyCluster(iRomanPotId, iPlaneId, iCluster)); 
        rpCluster->setXY(collection->xyCluster(iRomanPotId, iPlaneId, iCluster)); 
        rpCluster->setQuality(collection->qualityCluster(iRomanPotId, iPlaneId, iCluster));
        rpCluster->setPlaneId(iPlaneId);
        rpCluster->setRomanPotId(iRomanPotId);
      }
    }
  }

  for(Int_t iTrackPoint=0; iTrackPoint < collection->numberOfTrackPoints(); ++iTrackPoint){ 
    const StMuRpsTrackPoint *trackPoint = collection->trackPoint(iTrackPoint);
    StUPCRpsTrackPoint *rpTrackPoint = rpEvt->addTrackPoint();
    rpTrackPoint->setPosition(trackPoint->positionVec());
    rpTrackPoint->setRpId(trackPoint->rpId());
    for(UInt_t iPlane=0; iPlane < collection->numberOfPlanes(); ++iPlane)
      rpTrackPoint->setClusterId(trackPoint->clusterId(iPlane), iPlane);
    rpTrackPoint->setTime(trackPoint->time(0), 0); // pmtId = 0
    rpTrackPoint->setTime(trackPoint->time(1), 1); // pmtId = 1
    switch(trackPoint->quality()){
      case StMuRpsTrackPoint::rpsNormal: rpTrackPoint->setQuality(StUPCRpsTrackPoint::rpsNormal);
      break;      
      case StMuRpsTrackPoint::rpsGolden: rpTrackPoint->setQuality(StUPCRpsTrackPoint::rpsGolden);
      break; 
      default: rpTrackPoint->setQuality(StUPCRpsTrackPoint::rpsNotSet);
      break; 
    }
    mMuTrackPoints.push_back(trackPoint);
  }

  for(Int_t iTrack=0; iTrack < collection->numberOfTracks(); ++iTrack){
    StMuRpsTrack *track = collection->track(iTrack);
    StUPCRpsTrack *rpTrack = rpEvt->addTrack();
    for(Int_t iStation = 0; iStation < 2; ++iStation){ // Number Of Stations In Branch = 2, one trackpoint in each station
      const StMuRpsTrackPoint *trackPoint = track->trackPoint(iStation);
      if(!trackPoint)
        continue;
      for(Int_t i = 0; i < collection->numberOfTrackPoints(); ++i){
			if(trackPoint == mMuTrackPoints[i]){
				if(iStation==0){
					rpTrack->setFirstTrackPointId(i);
				}else{
					rpTrack->setSecondTrackPointId(i);
				}
				break;
			} 
      }
    }     
    rpTrack->setP(track->pVec()); 
    rpTrack->setBranch(track->branch());
    switch(track->type()){
      case StMuRpsTrack::rpsLocal: rpTrack->setType(StUPCRpsTrack::rpsLocal);
      break;      
      case StMuRpsTrack::rpsGlobal: rpTrack->setType(StUPCRpsTrack::rpsGlobal);
      break; 
      default: rpTrack->setType(StUPCRpsTrack::rpsUndefined);
      break; 
    } 
  } 
}//processEvent

























