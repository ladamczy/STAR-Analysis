//_____________________________________________________________________________
//    Class for making picoDst RP data 2019
//    Author: Truhlar Tomas
//_____________________________________________________________________________

//root headers
#include "TClonesArray.h"
#include "TIterator.h"
#include "TParticle.h"

//local headers
#include "StUPCRpsCluster.h"
#include "StUPCRpsTrack.h"
#include "StUPCRpsTrackPoint.h"
#include "StRPEvent.h"

ClassImp(StRPEvent);


TClonesArray *StRPEvent::mgTrackPoints = 0;
TClonesArray *StRPEvent::mgTracks = 0;
TClonesArray *StRPEvent::mgClusters = 0;


//_____________________________________________________________________________
StRPEvent::StRPEvent():
   mSiliconBunch(0),
   mClusters(0x0), mNClusters(0),
   mTracks(0x0), mNTracks(0),
   mTrackPoints(0x0), mNTrackPoints(0)
{
  //default constructor

  for (UInt_t iRP = 0; iRP < mNumberOfRomanPots; ++iRP){
    mNumberPlanes[iRP] = 0;
    mNumberPlanesWithCluster[iRP] = 0;
    mStatusRomanPot[iRP] = 0;
    for (UInt_t iPMT = 0; iPMT < mNumberOfPmtsInRp; ++iPMT){
      mADC[iRP][iPMT] = 0;
      mTAC[iRP][iPMT] = 0;
    }
    for (UInt_t iPlane = 0; iPlane < mNumberOfPlanesInRP; ++iPlane){
      mOffsetPlane[iRP][iPlane] = 0;
      mzPlane[iRP][iPlane] = 0;
      mAnglePlane[iRP][iPlane] = 0;
      mOrientationPlane[iRP][iPlane] = 0;
      mStatusPlane[iRP][iPlane] = 0;  
    }
  }
  if(!mgClusters) {
    mgClusters = new TClonesArray("StUPCRpsCluster");
    mClusters = mgClusters;
    mClusters->SetOwner(kTRUE);
  }

  if(!mgTrackPoints) {
    mgTrackPoints = new TClonesArray("StUPCRpsTrackPoint");
    mTrackPoints = mgTrackPoints;
    mTrackPoints->SetOwner(kTRUE);
  }

  if(!mgTracks) {
    mgTracks = new TClonesArray("StUPCRpsTrack");
    mTracks = mgTracks;
    mTracks->SetOwner(kTRUE);
  }

}//StRPEvent

//_____________________________________________________________________________
StRPEvent::StRPEvent(const StRPEvent &evt) {

  //copy constructor

  //one-to-one copy of data values

  mSiliconBunch = evt.mSiliconBunch;

  //Roman Pot loop
  for (UInt_t iRP = 0; iRP < mNumberOfRomanPots; ++iRP){
    mNumberPlanes[iRP] = evt.mNumberPlanes[iRP];
    mNumberPlanesWithCluster[iRP] = evt.mNumberPlanesWithCluster[iRP];
    mStatusRomanPot[iRP] = evt.mStatusRomanPot[iRP];
    //PMT loop
    for (UInt_t iPMT = 0; iPMT < mNumberOfPmtsInRp; ++iPMT){
      mADC[iRP][iPMT] = evt.mADC[iRP][iPMT];
      mTAC[iRP][iPMT] = evt.mTAC[iRP][iPMT];
    }//PMT loop
    //planes loop
    for (UInt_t iPlane = 0; iPlane < mNumberOfPlanesInRP; ++iPlane){
      mOffsetPlane[iRP][iPlane] = evt.mOffsetPlane[iRP][iPlane];
      mzPlane[iRP][iPlane] = evt.mzPlane[iRP][iPlane];
      mAnglePlane[iRP][iPlane] = evt.mAnglePlane[iRP][iPlane];
      mOrientationPlane[iRP][iPlane] = evt.mOrientationPlane[iRP][iPlane];
      mStatusPlane[iRP][iPlane] = evt.mStatusPlane[iRP][iPlane];
    }//planes loop
  }//Roman Pot loop

  //make fresh clones arrays
  mClusters = new TClonesArray("StUPCRpsCluster");
  mClusters->SetOwner(kTRUE);

  mTrackPoints = new TClonesArray("StUPCRpsTrackPoint");
  mTrackPoints->SetOwner(kTRUE);

  mTracks = new TClonesArray("StUPCRpsTrack");
  mTracks->SetOwner(kTRUE);

}//StRPEvent

//_____________________________________________________________________________
StRPEvent::~StRPEvent()
{
  //destructor

  if(mTrackPoints) {delete mTrackPoints; mTrackPoints = 0x0;}
  if(mTracks) {delete mTracks; mTracks = 0x0;}
  if(mClusters) {delete mClusters; mClusters = 0x0;}


}//~StRPEvent

//_____________________________________________________________________________
void StRPEvent::clearEvent()
{
  // clear event variables

  mTrackPoints->Clear("C"); 
  mNTrackPoints = 0;
  mTracks->Clear("C"); 
  mNTracks = 0;
  mClusters->Clear("C"); 
  mNClusters = 0;

}//clearEvent

//_____________________________________________________________________________
StUPCRpsCluster *StRPEvent::addCluster()
{
  // construct new  RP cluster

  return dynamic_cast<StUPCRpsCluster*>( mClusters->ConstructedAt(mNClusters++) );

}//addCluster

//_____________________________________________________________________________
StUPCRpsTrack *StRPEvent::addTrack()
{
  // construct new  RP track

  return dynamic_cast<StUPCRpsTrack*>( mTracks->ConstructedAt(mNTracks++) );

}//addTrack


//_____________________________________________________________________________
StUPCRpsTrackPoint *StRPEvent::addTrackPoint()
{
  // construct new  RP trackPoint

  return dynamic_cast<StUPCRpsTrackPoint*>( mTrackPoints->ConstructedAt(mNTrackPoints++) );

}//addTrackPoint

//_____________________________________________________________________________
UInt_t StRPEvent::getNumberOfClusters() const {

  //number of RP clusters in event

  if( !mClusters ) return 0;

  return mClusters->GetEntriesFast();

}//getNumberOfClusters

//_____________________________________________________________________________
StUPCRpsCluster *StRPEvent::getCluster(Int_t iCluster) const
{
  // get RP cluster

  return dynamic_cast<StUPCRpsCluster*>( mClusters->At(iCluster) );

}//getCluster

//_____________________________________________________________________________
UInt_t StRPEvent::getNumberOfTracks() const {

  //number of RP tracks in event

  if( !mTracks ) return 0;

  return mTracks->GetEntriesFast();

}//getNumberOfTracks

//_____________________________________________________________________________
StUPCRpsTrack *StRPEvent::getTrack(Int_t iTrack) const
{
  // get RP track

  return dynamic_cast<StUPCRpsTrack*>( mTracks->At(iTrack) );

}//getTrack

//_____________________________________________________________________________
UInt_t StRPEvent::getNumberOfTrackPoints() const {

  //number of RP trackPoints in event

  if( !mTrackPoints ) return 0;

  return mTrackPoints->GetEntriesFast();

}//getNumberOfTrackPoints

//_____________________________________________________________________________
StUPCRpsTrackPoint *StRPEvent::getTrackPoint(Int_t iTrackPoint) const
{
  // get RP trackPoint

  return dynamic_cast<StUPCRpsTrackPoint*>( mTrackPoints->At(iTrackPoint) );

}//getTrackPoint


