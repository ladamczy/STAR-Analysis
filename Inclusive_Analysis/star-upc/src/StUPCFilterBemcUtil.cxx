
//_____________________________________________________________________________
//
//    Utility class for BEMC in UPC filter maker
//    Author: Jaroslav Adam
//
//_____________________________________________________________________________

//c++ headers
#include <vector>

//StRoot headers
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StEmcUtil/projection/StEmcPosition.h"
#include "StEmcUtil/geometry/StEmcGeom.h"
#include "StEvent/StEmcCollection.h"
#include "StEvent/StEmcDetector.h"
#include "StEvent/StEmcClusterCollection.h"
#include "StEvent/StEmcCluster.h"
#include "StEvent/StEmcRawHit.h"

//local headers
#include "StUPCEvent.h"
#include "StUPCBemcCluster.h"

#include "StUPCFilterBemcUtil.h"

//_____________________________________________________________________________
StUPCFilterBemcUtil::StUPCFilterBemcUtil():
  mMagField(0), mEmcPos(0x0), mEmcGeom(0x0),
  mEmcHit(), mVecEmcHits(0x0), mEmcCluster(), mVecEmcCluster(0x0) {
  //constructor

  //initialize BEMC position and geometry
  mEmcPos = new StEmcPosition();
  mEmcGeom = StEmcGeom::instance("bemc");

  //vector of hits and clusters
  mVecEmcHits = new vector<emcHits>;
  mVecEmcCluster = new vector<emcCluster>;

}//StUPCFilterBemcUtil

//_____________________________________________________________________________
StUPCFilterBemcUtil::~StUPCFilterBemcUtil() {
  //destructor

  delete mVecEmcHits; mVecEmcHits=0;
  delete mVecEmcCluster; mVecEmcCluster=0;
  delete mEmcPos; mEmcPos=0;

}//~StUPCFilterBemcUtil

//_____________________________________________________________________________
void StUPCFilterBemcUtil::clear() {

  //clear hits and clusters vectors
  mVecEmcHits->clear();
  mVecEmcCluster->clear();

}//clear

//_____________________________________________________________________________
Bool_t StUPCFilterBemcUtil::processEvent(StMuDst *muDst, StUPCEvent *upcEvt) {

  //fills BEMC clusters vector, sets BEMC multiplicity in UPC event
  //and fills vector of raw hits for track matching

  //BEMC clusters
  StEmcCollection *emcCollection = muDst->emcCollection();
  if( !emcCollection ) return kFALSE;

  //take BEMC from collection
  StEmcDetector *barrel = emcCollection->detector(kBarrelEmcTowerId);
  if( !barrel ) return kFALSE;
  if( !barrel->cluster() ) return kFALSE;

  //BEMC clusters
  const StSPtrVecEmcCluster &clusters = barrel->cluster()->clusters();

  //set BEMC multiplicity in output UPC event as number of all clusters
  upcEvt->setBEMCMultiplicity( clusters.size() );

  UInt_t nCls=0; // number of clusters written to clusters vector

  //clusters loop
  for(UInt_t icl=0; icl<clusters.size(); icl++) {
    StEmcCluster *cls = clusters[icl];
    if( !cls ) continue;

    //hits in cluster
    const StPtrVecEmcRawHit &hits = cls->hit();

    Float_t htEn = -9999.; // high tower hit energy
    Int_t htID = -9999; // high tower softID

    //hits loop
    for(UInt_t ihit=0; ihit<hits.size(); ihit++) {
      StEmcRawHit *rawhit = hits[ihit];
      if( !rawhit ) continue;

      //get soft id of raw hit
      Int_t emcModule = rawhit->module();
      Int_t emcEtabin = rawhit->eta();
      Int_t emcSub = rawhit->sub();
      Int_t emcSoftId=0;
      mEmcGeom->getId(emcModule, emcEtabin, emcSub, emcSoftId);

      //fill vector of emc hits
      mEmcHit.hitSoftId = emcSoftId;
      mEmcHit.clsId = nCls; //cluster ID in hit as position of cluster in clusters vector
      //hit energy
      Float_t en = rawhit->energy();
      mEmcHit.hitE = en;
      //high tower energy and softID
      if( en > htEn ) {
        htEn = en;
        htID = emcSoftId;
      }
      mVecEmcHits->push_back( mEmcHit );

    }//hits loop

    //fill vector of clusters
    mEmcCluster.clsEta = cls->eta();
    mEmcCluster.clsPhi = cls->phi();
    mEmcCluster.clsSigmaEta = cls->sigmaEta();
    mEmcCluster.clsSigmaPhi = cls->sigmaPhi();
    mEmcCluster.clsE = cls->energy();
    mEmcCluster.clsHT = htEn;
    mEmcCluster.clsHTsoftID = htID;
    mEmcCluster.isMatched = kFALSE;
    mVecEmcCluster->push_back( mEmcCluster ); //write cluster to vector
    nCls++; //increment number of clusters writtent to vector

  }//clusters loop

  return kTRUE;

}//processEvent

//_____________________________________________________________________________
Short_t StUPCFilterBemcUtil::matchBEMC(StMuTrack *track,
                                    Double_t &emcPhi, Double_t &emcEta, Double_t &emcPt, Bool_t &proj, UInt_t &clsId, Float_t &hitE) {
  //function to match track to BEMC cluster, return number of matched raw hits

  //projection to BEMC
  StThreeVectorD emcPos, emcMom;
  if( !mEmcPos->trackOnEmc(&emcPos, &emcMom, track, mMagField) ) return 0;

  //projection successfull
  proj = kTRUE;
  emcPhi = emcPos.phi();
  emcEta = emcPos.pseudoRapidity();
  emcPt = emcMom.perp();

  if( mVecEmcHits->size() == 0 ) return 0; //no BEMC hits in event

  //get BEMC soft id
  Int_t emcSoftId=0;
  mEmcGeom->getId(emcPos.phi(), emcPos.pseudoRapidity(), emcSoftId);

  //find projected soft id in vector of raw BEMC hits
  Short_t nhits=0;
  //hits vector loop
  for(vector<emcHits>::const_iterator hit = mVecEmcHits->cbegin(); hit != mVecEmcHits->cend(); ++hit) {
    if( hit->hitSoftId != emcSoftId ) continue;

    clsId = hit->clsId;
    mVecEmcCluster->at(clsId).isMatched = kTRUE; // mark cluster as matched to track
    hitE = hit->hitE;

    nhits++;

  }//hits vector loop

  return nhits;

}//matchBEMC

//_____________________________________________________________________________
void StUPCFilterBemcUtil::writeBEMC(StUPCEvent *upcEvt) {
  //function to write BEMC clusters to output UPC event,
  //only clusters matched to a track are written

  //clusters vector loop
  for(vector<emcCluster>::const_iterator cluster = mVecEmcCluster->cbegin(); cluster != mVecEmcCluster->cend(); ++cluster) {

    //take only clusters to which any track is matching
    if( !cluster->isMatched ) continue;

    //put cluster in clones array in UPC event
    StUPCBemcCluster *upcCls = upcEvt->addCluster();
    upcCls->setEta( cluster->clsEta );
    upcCls->setPhi( cluster->clsPhi );
    upcCls->setSigmaEta( cluster->clsSigmaEta );
    upcCls->setSigmaPhi( cluster->clsSigmaPhi );
    upcCls->setEnergy( cluster->clsE );
    upcCls->setHTEnergy( cluster->clsHT );
    upcCls->setHTsoftID( cluster->clsHTsoftID );
    //cluster ID as position of cluster in clusters vector
    upcCls->setId( (UInt_t) distance(mVecEmcCluster->cbegin(), cluster) );

  }//clusters vector loop

}//writeBEMC
















































