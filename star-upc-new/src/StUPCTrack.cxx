
//_____________________________________________________________________________
//    Class for UPC data
//    Author: Jaroslav Adam
//
//    contains parameters of the central track relevant for UPC analysis
//_____________________________________________________________________________

//root headers
#include "TLorentzVector.h"

//local headers
#include "StUPCTrack.h"
#include "StUPCEvent.h"
#include "StUPCBemcCluster.h"
#include "StUPCVertex.h"

ClassImp(StUPCTrack);

//_____________________________________________________________________________
StUPCTrack::StUPCTrack(): TObject(),
  mFlags(0), mPt(0), mEta(0), mPhi(0),
  mCurvature(0), mDipAngle(), mPhase(0),
  mDcaXY(0), mDcaZ(0), mCharge(0), mNhits(0), mNhitsFit(0), mChi2(0),
  mNhitsDEdx(0), mDEdxSignal(0),
  mBemcPt(-999), mBemcEta(-999), mBemcPhi(-999), mBemcClsId(0), mBemcHitE(-999),
  mTofTime(0), mTofPathLength(0), mVtxId(0), mEvt(0x0), mIdTruth(0), mQATruth(0)
{
  //default constructor

  for(Int_t i=0; i<mNpart; i++) mNSigmasTPC[i] = -999;

}//StUPCEvent

//_____________________________________________________________________________
StUPCTrack::~StUPCTrack() {
  //destructor

}//~StUPCTrack

//_____________________________________________________________________________
void StUPCTrack::Clear(Option_t *)
{
  //clear track object, overridden from TObject
  //for use in TClonesArray

  mFlags = 0;

  //nullify BEMC and TOF members because not every
  //track is matched to them
  mBemcPt = -999.;
  mBemcEta = -999.;
  mBemcPhi = -999.;
  mBemcClsId = 0;
  mBemcHitE = -999.;

  mTofTime = 0;
  mTofPathLength = 0;

}//Clear

//_____________________________________________________________________________
void StUPCTrack::setFlag(Flag flg)
{
  //set track flag

  mFlags |= (1 << flg);

}//set flag

//_____________________________________________________________________________
Bool_t StUPCTrack::getFlag(Flag flg) const
{
  //get track flag

  if( mFlags & (1 << flg) ) return kTRUE;
  return kFALSE;

}//getFlag

//_____________________________________________________________________________
void StUPCTrack::getLorentzVector(TLorentzVector &lvec, Double_t mass) const
{
  // get track 4-momentum at DCA to primary vertex

  lvec.SetPtEtaPhiM(mPt, mEta, mPhi, mass);

}//getLorentzVector

//_____________________________________________________________________________
void StUPCTrack::getMomentum(TVector3 &vec) const
{
  // get track 3-momentum at DCA to primary vertex

  vec.SetPtEtaPhi(mPt, mEta, mPhi);

}//getMomentum

//_____________________________________________________________________________
Double_t StUPCTrack::getBemcPmag() const
{
  // get total momentum at BEMC

  TVector3 vec;
  vec.SetPtEtaPhi(mBemcPt, mBemcEta, mBemcPhi);

  return vec.Mag();

}//getBemcPmag

//_____________________________________________________________________________
StUPCBemcCluster *StUPCTrack::getBemcCluster() const
{
  //get BEMC cluster to which the track is matched to

  if(!mEvt) return 0x0;

  return mEvt->getClusterId(mBemcClsId);

}//getBemcCluster

//_____________________________________________________________________________
void StUPCTrack::getBemcLorentzVector(TLorentzVector &blvec, Double_t mass) const
{
  // get track 4-momentum at BEMC position

  blvec.SetPtEtaPhiM(mBemcPt, mBemcEta, mBemcPhi, mass);

}//getBemcLorentzVector

Double_t StUPCTrack::getT0(Double_t mass) const{
  //guards against crash without tof data:
  //positive ToF time
  //positive ToF length
  //track matched with ToF hit
  if(getTofTime()<=0||getTofPathLength()<=0||!getFlag(StUPCTrack::kTof)){
    return -999.;
  }
  //calculations
  TVector3 mom;
  mom.SetPtEtaPhi(mPt, mEta, mPhi);
  //t1-t0 (flight time) is calculated from momentum, known mass and known length of flight path
  //t0 = t1-(t1-t0), in ns
  return getTofTime()-getTofPathLength()/100.0/0.299792458*sqrt(1+pow(mass/mom.Mag(), 2));
}//getT0

Double_t StUPCTrack::getMass(Double_t T0) const{
  //guards against crash without tof data:
  //positive ToF time
  //positive ToF length
  //track matched with ToF hit
  if(getTofTime()<=0||getTofPathLength()<=0||!getFlag(StUPCTrack::kTof)){
    return -999.;
  }
  //calculations
  TVector3 mom;
  mom.SetPtEtaPhi(mPt, mEta, mPhi);
  //t1-t0, in ns
  Double_t time = getTofTime()-T0;
  //distance in m (conversion cm->m)
  Double_t dist = getTofPathLength()/100.0;
  //velocity, in m/ns
  Double_t vel = dist/time;
  Double_t c = 0.299792458;
  //beta & gamma factors
  Double_t beta = vel/c;
  Double_t gamma = 1/sqrt(1-beta*beta);
  //return p/(beta*gamma)
  return mom.Mag()/beta/gamma;
}//getMass

//_____________________________________________________________________________
StUPCVertex *StUPCTrack::getVertex() const
{
  // get primary vertex associated with track

  if(!mEvt) return 0x0;

  return mEvt->getVertexId(mVtxId);

}//getVertex
































