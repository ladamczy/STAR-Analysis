
//_____________________________________________________________________________
//    Class for UPC data
//    Author: Jaroslav Adam
//
//    UPC event
//_____________________________________________________________________________

//root headers
#include "TClonesArray.h"
#include "TIterator.h"
#include "TParticle.h"

//local headers
#include "StUPCTrack.h"
#include "StUPCBemcCluster.h"
#include "StUPCTofHit.h"
#include "StUPCVertex.h"

#include "StUPCEvent.h"

ClassImp(StUPCEvent);

TClonesArray *StUPCEvent::mgUPCTracks = 0;
TClonesArray *StUPCEvent::mgUPCBemcClusters = 0;
TClonesArray *StUPCEvent::mgUPCTOFHits = 0;
TClonesArray *StUPCEvent::mgUPCVertices = 0;
TClonesArray *StUPCEvent::mgMCParticles = 0;

//_____________________________________________________________________________
StUPCEvent::StUPCEvent():
  mRunNum(0), mEvtNum(0), mFillNum(0), mbCrossId(0), mbCrossId7bit(0),
  mMagField(0), mZdcEastRate(0), mZdcWestRate(0), mZdcCoincRate(0),
  mLastDSM0(0), mLastDSM1(0), mLastDSM3(0),
  mZdcEastUA(0), mZdcWestUA(0), mZdcEastTDC(0), mZdcWestTDC(0),
  mZdcTimeDiff(0), mZdcVertexZ(0),
  mBBCSmallEast(0), mBBCSmallWest(0), mBBCLargeEast(0), mBBCLargeWest(0),
  mVPDSumEast(0), mVPDSumWest(0), mVPDTimeDiff(0),
  mTofMult(0), mBemcMult(0),
  mNGlobTracks(0), mNPrimTracks(0),mNPrimVertices(0),
  mUPCTracks(0x0), mNtracks(0),
  mUPCBemcClusters(0x0), mNclusters(0),
  mUPCTOFHits(0x0), mNhits(0),
  mUPCVertices(0x0), mNvertices(0),
  mMCParticles(0x0), mNmc(0)
{
  //default constructor

  mTrgIDs.Set(0);

  for (Int_t ipmt=0; ipmt<mNZdcPmt; ipmt++) {
    mZdcEastADC[ipmt] = 0;
    mZdcWestADC[ipmt] = 0;
  }
  for(Int_t i=0; i<mNZdcSmd; i++) {
    mZdcSmdEast[i] = 0;
    mZdcSmdWest[i] = 0;
  }


  if(!mgUPCTracks) {
    mgUPCTracks = new TClonesArray("StUPCTrack");
    mUPCTracks = mgUPCTracks;
    mUPCTracks->SetOwner(kTRUE);
  }

  if(!mgUPCBemcClusters) {
    mgUPCBemcClusters = new TClonesArray("StUPCBemcCluster");
    mUPCBemcClusters = mgUPCBemcClusters;
    mUPCBemcClusters->SetOwner(kTRUE);
  }

  if(!mgUPCTOFHits) {
    mgUPCTOFHits = new TClonesArray("StUPCTofHit");
    mUPCTOFHits = mgUPCTOFHits;
    mUPCTOFHits->SetOwner(kTRUE);
  }

  if(!mgUPCVertices) {
    mgUPCVertices = new TClonesArray("StUPCVertex");
    mUPCVertices = mgUPCVertices;
    mgUPCVertices->SetOwner(kTRUE);
  }

}//StUPCEvent

//_____________________________________________________________________________
StUPCEvent::~StUPCEvent()
{
  //destructor

  if(mUPCTracks) {delete mUPCTracks; mUPCTracks = 0x0;}
  if(mUPCBemcClusters) {delete mUPCBemcClusters; mUPCBemcClusters = 0x0;}
  if(mUPCTOFHits) {delete mUPCTOFHits; mUPCTOFHits = 0x0;}
  if(mUPCVertices) {delete mUPCVertices; mUPCVertices = 0x0;}
  if(mMCParticles) {delete mMCParticles; mMCParticles = 0x0;}

}//~StUPCEvent

//_____________________________________________________________________________
void StUPCEvent::clearEvent()
{
  // clear event variables

  mTrgIDs.Set(0);

  mBemcMult = 0;

  mUPCTracks->Clear("C");
  mNtracks = 0;

  mUPCBemcClusters->Clear("C");
  mNclusters = 0;

  mUPCTOFHits->Clear("C");
  mNhits = 0;

  mUPCVertices->Clear("C");
  mNvertices = 0;

  if(mMCParticles) {
    mMCParticles->Clear("C");
    mNmc = 0;
  }

}//clearEvent

//_____________________________________________________________________________
void StUPCEvent::addTriggerId(Int_t id) {
  //add fired trigger ID

  //position to put the ID
  Int_t pos = mTrgIDs.GetSize();

  //extend the array
  mTrgIDs.Set(pos+1);

  //put the ID
  mTrgIDs.AddAt(id, pos);

}//addTriggerId

//_____________________________________________________________________________
void StUPCEvent::setZDCEastADC(UShort_t signal, Int_t pmt) {

  //set ZDC PMT signal, check for bounds, east side
  Int_t idx = pmt - 1;
  if( idx < 0 or idx > mNZdcPmt-1 ) return;

  mZdcEastADC[idx] = signal;

}//setZDCEastADC

//_____________________________________________________________________________
void StUPCEvent::setZDCWestADC(UShort_t signal, Int_t pmt) {

  //set ZDC PMT signal, check for bounds, west side
  Int_t idx = pmt - 1;
  if( idx < 0 or idx > mNZdcPmt-1 ) return;

  mZdcWestADC[idx] = signal;

}//setZDCWestADC

//_____________________________________________________________________________
void StUPCEvent::setZdcSMDEast(Int_t verthori, Int_t strip, UShort_t smd) {

  //ZDC SMD data, indexing for verthori is 0=vertical, 1=Horizontal and strips go as strip=1, strip<9

  if( verthori < 0 or verthori > 1 ) return;
  if( strip < 1 or strip > 8 ) return;

  mZdcSmdEast[strip-1 + verthori*8] = smd;

}//setZdcSMDEast

//_____________________________________________________________________________
void StUPCEvent::setZdcSMDWest(Int_t verthori, Int_t strip, UShort_t smd) {

  //ZDC SMD data, indexing for verthori is 0=vertical, 1=Horizontal and strips go as strip=1, strip<9

  if( verthori < 0 or verthori > 1 ) return;
  if( strip < 1 or strip > 8 ) return;

  mZdcSmdWest[strip-1 + verthori*8] = smd;

}//setZdcSMDWest

//_____________________________________________________________________________
StUPCTrack *StUPCEvent::addTrack()
{
  // construct new upc track

  return dynamic_cast<StUPCTrack*>( mUPCTracks->ConstructedAt(mNtracks++) );

}//addTrack

//_____________________________________________________________________________
StUPCBemcCluster *StUPCEvent::addCluster()
{
  // construct new BEMC cluster

  return dynamic_cast<StUPCBemcCluster*>( mUPCBemcClusters->ConstructedAt(mNclusters++) );

}//addCluster

//_____________________________________________________________________________
StUPCTofHit *StUPCEvent::addHit()
{
  // construct new TOF hit

  return dynamic_cast<StUPCTofHit*>( mUPCTOFHits->ConstructedAt(mNhits++) );

}//addHit

//_____________________________________________________________________________
StUPCVertex *StUPCEvent::addVertex()
{
  //construct new UPC vertex

  return dynamic_cast<StUPCVertex*>( mUPCVertices->ConstructedAt(mNvertices++) );

}//addVertex

//_____________________________________________________________________________
void StUPCEvent::setIsMC(Bool_t mc) {

  // set the event as MC, initialize mc array

  if(!mc) return;

  if(!mgMCParticles) {
    mgMCParticles = new TClonesArray("TParticle");
    mMCParticles = mgMCParticles;
    mMCParticles->SetOwner(kTRUE);
  }

}//setIsMC

//_____________________________________________________________________________
TParticle *StUPCEvent::addMCParticle()
{
  // construct new mc TParticle

  if(!mMCParticles) return 0x0;

  return dynamic_cast<TParticle*>( mMCParticles->ConstructedAt(mNmc++) );

}//addMCParticle

//_____________________________________________________________________________
Bool_t StUPCEvent::isTrigger(Int_t id) const
{
  //get fired trigger ID
  for(Int_t i=0; i<mTrgIDs.GetSize(); i++) {
    if( mTrgIDs.At(i) == id ) return kTRUE;
  }

  return kFALSE;

}//setTrigger

//_____________________________________________________________________________
Int_t StUPCEvent::getNumberOfTracks() const {

  //number of tracks in event

  if( !mUPCTracks ) return 0;

  return mUPCTracks->GetEntriesFast();

}//getNumberOfTracks

//_____________________________________________________________________________
StUPCTrack *StUPCEvent::getTrack(Int_t iTrack) const
{
  // get upc track

  StUPCTrack *track = dynamic_cast<StUPCTrack*>( mUPCTracks->At(iTrack) );
  if(track) track->setEvent( const_cast<StUPCEvent*>(this) );

  return track;

}//getTrack

//_____________________________________________________________________________
Int_t StUPCEvent::getNumberOfClusters() const {

  //number of BEMC clusters in event

  if( !mUPCBemcClusters ) return 0;

  return mUPCBemcClusters->GetEntriesFast();

}//getNumberOfClusters

//_____________________________________________________________________________
StUPCBemcCluster *StUPCEvent::getCluster(Int_t iCls) const
{
  // get BEMC cluster

  return dynamic_cast<StUPCBemcCluster*>( mUPCBemcClusters->At(iCls) );

}//getCluster

//_____________________________________________________________________________
StUPCBemcCluster *StUPCEvent::getClusterId(UInt_t clsId) const
{
  // get BEMC cluster according to origial ID

  //clusters loop
  StUPCBemcCluster *cls=0x0;
  TIterator *clsIter = makeClustersIter();
  while( (cls = dynamic_cast<StUPCBemcCluster*>( clsIter->Next() )) != NULL ) {

    //get cluster with a given ID
    if( cls->getId() == clsId ) {delete clsIter; return cls;}
  }//clusters loop

  delete clsIter;
  return 0x0;

}//getCluster

//_____________________________________________________________________________
TIterator *StUPCEvent::makeClustersIter() const
{
  // create iterator for BEMC clusters clones array

  return mUPCBemcClusters->MakeIterator();

}//getClustersIter

//_____________________________________________________________________________
Int_t StUPCEvent::getNumberOfHits() const {

  //number of TOF hits in event

  if( !mUPCTOFHits ) return 0;

  return mUPCTOFHits->GetEntriesFast();

}//getNumberOfHits

//_____________________________________________________________________________
StUPCTofHit *StUPCEvent::getHit(Int_t iHit) const
{
  // get TOF hit

  return dynamic_cast<StUPCTofHit*>( mUPCTOFHits->At(iHit) );

}//getHit

//_____________________________________________________________________________
Int_t StUPCEvent::getNumberOfVertices() const {

  //number of primary vertices in event

  if( !mUPCVertices ) return 0;

  return mUPCVertices->GetEntriesFast();

}//getNumberOfVertices

//_____________________________________________________________________________
StUPCVertex *StUPCEvent::getVertex(Int_t iVtx) const
{
  //get UPC vertex from clones array

  return dynamic_cast<StUPCVertex*>( mUPCVertices->At(iVtx) );

}//getVertex

//_____________________________________________________________________________
StUPCVertex *StUPCEvent::getVertexId(UInt_t vtxId) const
{
  // get vertex according to original ID

  //vertex loop
  StUPCVertex *vtx = 0x0;
  TIterator *vtxIter = makeVerticesIter();
  while( (vtx = dynamic_cast<StUPCVertex*>( vtxIter->Next() )) != NULL ) {

    //get vertex with a given ID
    if( vtx->getId() == vtxId ) {delete vtxIter; return vtx;}
  }//vertex loop

  delete vtxIter;
  return 0x0;

}//getVertexId

//_____________________________________________________________________________
TIterator *StUPCEvent::makeVerticesIter() const
{
  // create iterator for vertices clones array

  return mUPCVertices->MakeIterator();

}//makeVerticesIter

//_____________________________________________________________________________
Int_t StUPCEvent::getNumberOfMCParticles() const {

  //number of MC particles in event

  if( !mMCParticles ) return 0;

  return mMCParticles->GetEntriesFast();

}//getNumberOfMCParticles

//_____________________________________________________________________________
TParticle *StUPCEvent::getMCParticle(Int_t iMC) const
{
  // get mc TParticle

  if(!mMCParticles) return 0x0;

  return dynamic_cast<TParticle*>( mMCParticles->At(iMC) );

}//getMCParticle

















