
//_____________________________________________________________________________
//    Class for UPC filter
//    Author: Jaroslav Adam
//
//    Fills structure of StUPCEvent
//_____________________________________________________________________________

//c++ headers
#include "string.h"
#include <vector>
#include <map>

//root headers
#include "TObjArray.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1I.h"
#include "TClonesArray.h"
#include "TDatabasePDG.h"
#include "TParticle.h"

//StRoot headers
#include "StMessMgr.h"
#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StMuDSTMaker/COMMON/StMuPrimaryVertex.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuMcTrack.h"
#include "StMuDSTMaker/COMMON/StMuMcVertex.h"
#include "StMuDSTMaker/COMMON/StMuBTofHit.h"

#include "StEvent/StTriggerData.h"
#include "StEvent/StRunInfo.h"
#include "StEvent/StEventSummary.h"

//local headers
#include "StUPCEvent.h"
#include "StUPCTrack.h"
#include "StUPCBemcCluster.h"
#include "StUPCVertex.h"
#include "StUPCTofHit.h"

#include "StRPEvent.h"

#include "StUPCFilterTrgUtil.h"
#include "StUPCFilterBemcUtil.h"
#include "StUPCFilterRPUtil.h"
#include "StUPCFilterMaker.h"

ClassImp(StUPCFilterMaker);

//_____________________________________________________________________________
StUPCFilterMaker::StUPCFilterMaker(StMuDstMaker *maker, string outnam) : StMaker("StReadMuDstMaker"),
  mMaker(maker), mMuDst(0x0), mIsMC(0), mOutName(outnam), mOutFile(0x0),
  mHistList(0x0), mCounter(0x0), mErrCounter(0x0),
  mUPCEvent(0x0), mUPCTree(0x0),  mTrgUtil(0x0), mBemcUtil(0x0),
  mMakeRP(0), mRPUtil(0x0), mRPEvent(0x0)
{
  //constructor

  LOG_INFO << "StUPCFilterMaker::StUPCFilterMaker() called" << endm;

}//StUPCFilterMaker

//_____________________________________________________________________________
StUPCFilterMaker::~StUPCFilterMaker()
{
  //destructor

  LOG_INFO << "StUPCFilterMaker::~StUPCFilterMaker() destructor called" << endm;

  delete mTrgUtil; mTrgUtil=0;
  delete mBemcUtil; mBemcUtil=0;
  delete mRPUtil; mRPUtil=0;
  delete mRPEvent; mRPEvent=0;
  delete mHistList; mHistList=0;
  delete mCounter; mCounter=0;
  delete mErrCounter; mErrCounter=0;
  delete mUPCTree; mUPCTree=0;
  delete mUPCEvent; mUPCEvent=0;
  delete mOutFile; mOutFile=0;

}//~StUPCFilterMaker

//_____________________________________________________________________________
void StUPCFilterMaker::addTriggerId(UInt_t id, Int_t rmin, Int_t rmax) {

  //add triggger ID to the table of IDs, set run range optionally

  mTrgIDs.push_back(id);

  mTrgRanLo.push_back(rmin);
  mTrgRanHi.push_back(rmax);

}//StUPCFilterMaker::addTriggerId

//_____________________________________________________________________________
Int_t StUPCFilterMaker::Init() {
  //called at the beginning

  //show active trigger IDs
  LOG_INFO << "StUPCFilterMaker::Init(), using the following trigger IDs:" << endm;
  for(UInt_t i=0; i<mTrgIDs.size(); i++) {
    LOG_INFO << "  ID: " << mTrgIDs[i] << ", run range: " << mTrgRanLo[i] << " - " << mTrgRanHi[i] << endm;
  }

  //utility for trigger data, BBC and ZDC
  mTrgUtil = new StUPCFilterTrgUtil();

  //initialize BEMC matching utility
  mBemcUtil = new StUPCFilterBemcUtil();

  //create the output file
  mOutFile = new TFile(mOutName.c_str(), "recreate");

  //output UPC event
  mUPCEvent = new StUPCEvent();
  //configure the UPC event
  if( mIsMC > 0 ) mUPCEvent->setIsMC( kTRUE );

  //configure for Roman Pot event
  if( mMakeRP ) {
    mRPUtil = new StUPCFilterRPUtil();
    mRPEvent = new StRPEvent();
  }

  //create the tree
  mUPCTree = new TTree("mUPCTree", "mUPCTree");
  //add branch with event object
  mUPCTree->Branch("mUPCEvent", &mUPCEvent);

  //add RP branch
  if( mMakeRP ) {
    mUPCTree->Branch("mRPEvent", &mRPEvent);
  }

  //output histograms
  mHistList = new TList();
  mHistList->SetOwner();

  //counter of processed events
  mCounter = new TH1I("mCounter", "mCounter", kMaxCnt-1, 1, kMaxCnt);
  mHistList->Add(mCounter);

  //counter for errors encountered during the analysis
  mErrCounter = new TH1I("mErrCounter", "mErrCounter", kMaxErrCnt-1, 1, kMaxErrCnt);
  mHistList->Add(mErrCounter);

  return kStOk;

}//Init

//_____________________________________________________________________________
Int_t StUPCFilterMaker::Make()
{
  //called for each event

  mUPCEvent->clearEvent(); //clear the output UPC event
  mBemcUtil->clear(); //clear data structures in BEMC util
  if( mRPUtil ) mRPUtil->clear(); // clear data structures in RP util
  if( mRPEvent ) mRPEvent->clearEvent(); // clear output RP event if present

  //input muDst data
  mMuDst = mMaker->muDst();
  if( !mMuDst ) {
    LOG_INFO << "StUPCFilterMaker::Make() no muDst input" << endm;
    mErrCounter->Fill( kErrNoEvt ); // no muDst input, same err flag as event
    return kStErr;
  }
  //input mu event
  StMuEvent *evt = mMuDst->event();
  if( !evt ) {
    LOG_INFO << "StUPCFilterMaker::Make() no input event" << endm;
    mErrCounter->Fill( kErrNoEvt ); // no input event
    return kStErr;
  }

  mCounter->Fill( kAna ); // analyzed events

  //mc
  if( mIsMC > 0 ) {
    if( !runMC() ) {
      LOG_INFO << "StUPCFilterMaker::Make() no MC" << endm;
      mErrCounter->Fill( kNoMC );
    }
  }

  //trigger
  const StTriggerId trgId = evt->triggerIdCollection().nominal();
  Int_t runnum = evt->runNumber();

  Bool_t isTrg = kFALSE; //determine whether at least one of trigger IDs was fired
  for(UInt_t i=0; i<mTrgIDs.size(); i++) {
    // run range for a given trigger ID
    if( mTrgRanLo[i] != 0 && runnum < mTrgRanLo[i] ) continue;
    if( mTrgRanHi[i] != 0 && runnum > mTrgRanHi[i] ) continue;
    //test trigger ID at 'i'
    if( !trgId.isTrigger( mTrgIDs[i] ) ) continue;

    //trigger ID was fired
    mUPCEvent->addTriggerId( mTrgIDs[i] );
    isTrg = kTRUE;

  }
  if( mIsMC > 0 ) isTrg = kTRUE; //override for MC
  if( !isTrg ) return kStOk;
  //event passed the trigger

  mCounter->Fill( kTrg ); // events after trigger

  //run number
  mUPCEvent->setRunNumber( runnum );
  //event number
  mUPCEvent->setEventNumber( evt->eventNumber() );

  //beam fill number
  const StRunInfo &runInfo = evt->runInfo();
  Int_t fillY = (Int_t) runInfo.beamFillNumber(yellow);
  Int_t fillB = (Int_t) runInfo.beamFillNumber(blue);
  if( fillY != fillB ) {
    LOG_INFO << "StUPCFilterMaker::Make() fill number mismatch" << endm;
    //fill number mismatch
    mErrCounter->Fill( kErrFillMsc );
  }
  mUPCEvent->setFillNumber( fillY );

  //bunch crossing ID
  const StL0Trigger &l0trig = evt->l0Trigger();
  mUPCEvent->setBunchCrossId( l0trig.bunchCrossingId() );
  mUPCEvent->setBunchCrossId7bit( l0trig.bunchCrossingId7bit(runnum) );

  //magnetic field in UPC event
  const StEventSummary &evtSummary = evt->eventSummary();
  mUPCEvent->setMagneticField( evtSummary.magneticField() );

  //ZDC rates
  mUPCEvent->setZDCEastRate( runInfo.zdcEastRate() );
  mUPCEvent->setZDCWestRate( runInfo.zdcWestRate() );
  mUPCEvent->setZDCCoincRate( runInfo.zdcCoincidenceRate() );

  //trigger data for DSM, ZDC, BBC and TOF
  const StTriggerData *trgdat = evt->triggerData();
  if( !trgdat && mIsMC==0 ) {
    mErrCounter->Fill( kNoTrgDat ); // no trigger data
    LOG_INFO << "StUPCFilterMaker::Make() no trigger data" << endm;
    return kStErr;
  }
  if( trgdat ) {
    mTrgUtil->processEvent(trgdat, mUPCEvent);
    mCounter->Fill( kTrgDat ); // events having trigger data
  }

  //number of global, primary tracks and vertices
  mUPCEvent->setNGlobTracks( evtSummary.numberOfGoodTracks() );
  mUPCEvent->setNPrimTracks( evtSummary.numberOfGoodPrimaryTracks() );
  mUPCEvent->setNPrimVertices( mMuDst->numberOfPrimaryVertices() );

  //BEMC util
  //magnetic field for track projection to BEMC
  mBemcUtil->setMagField( evtSummary.magneticField()/10. ); //conversion ok, checked by track->pt() and emcPt
  //fill structures with clusters and hits
  if( mBemcUtil->processEvent(mMuDst, mUPCEvent) ) mCounter->Fill( kBemc ); // events having BEMC clusters

  //mark primary vertices with at least one TOF or BEMC matched track
  map<UInt_t, Bool_t> vtxMap;
  //vertex map loop
  for(UInt_t ivtx=0; ivtx<mMuDst->numberOfPrimaryVertices(); ivtx++) {
    StMuDst::setVertexIndex(ivtx);

    //initialize the map
    vtxMap[ivtx] = kFALSE;

    //run over tracks for this vertex
    TObjArray *trkArray = mMuDst->primaryTracks();
    for(Int_t itrk=0; itrk<trkArray->GetEntriesFast(); itrk++) {
      StMuTrack *track = dynamic_cast<StMuTrack*>( trkArray->At(itrk) );
      //evaluate the matching
      //BEMC
      UInt_t clsId=0;
      Double_t emcPhi=-999., emcEta=-999., emcPt=-999.;
      Bool_t emcProj=kFALSE;
      Float_t hitE=-999.;
      Short_t nhitsBemc = mBemcUtil->matchBEMC(track, emcPhi, emcEta, emcPt, emcProj, clsId, hitE);
      Bool_t matchBemc = nhitsBemc > 0 ? kTRUE : kFALSE;
      //TOF
      const StMuBTofPidTraits &tofPid = track->btofPidTraits();
      Bool_t matchTof = tofPid.matchFlag() != 0 ? kTRUE : kFALSE;

      //mark current primary vertex as having matched track
      if( matchBemc == kTRUE or matchTof == kTRUE ) {
        vtxMap[ivtx] = kTRUE;
        break;
      }
    }//tracks for this vertex
  }//vertex map loop

  //vertex loop
  for(UInt_t ivtx=0; ivtx<mMuDst->numberOfPrimaryVertices(); ivtx++) {
    //select only vertices with at least one matched track
    //only in data or embedding MC
    if( mIsMC != 1 && vtxMap[ivtx] == kFALSE ) continue;

    //static call to set current primary vertex
    StMuDst::setVertexIndex(ivtx);

    //get array of primary tracks
    TObjArray *trkArray = mMuDst->primaryTracks();
    if( !trkArray ) continue;

    Int_t nSelTracks = 0; //number of tracks selected to write to output UPC event

    //tracks loop
    for(Int_t itrk=0; itrk<trkArray->GetEntriesFast(); itrk++) {
      StMuTrack *track = dynamic_cast<StMuTrack*>( trkArray->At(itrk) );
      if( !track ) continue;

      //matching to BEMC cluster
      UInt_t clsId=0;
      Double_t emcPhi=-999., emcEta=-999., emcPt=-999.;
      Bool_t emcProj=kFALSE;
      Float_t hitE=-999.;
      Short_t nhitsBemc = mBemcUtil->matchBEMC(track, emcPhi, emcEta, emcPt, emcProj, clsId, hitE);
      Bool_t matchBemc = nhitsBemc > 0 ? kTRUE : kFALSE;

      //TOF matching
      const StMuBTofPidTraits &tofPid = track->btofPidTraits();
      Bool_t matchTof = tofPid.matchFlag() != 0 ? kTRUE : kFALSE;

      //require at least one match, only in data or embedding MC
      //if( mIsMC!=1 && !matchBemc && !matchTof ) continue;

      //track matched to BEMC or TOF and selected to write to output UPC event
      nSelTracks++;

      //UPC track
      StUPCTrack *upcTrack = mUPCEvent->addTrack();
      upcTrack->setPtEtaPhi(track->pt(), track->eta(), track->phi());
      TVector3 origin(0.0,0.0,0.0);
      if( track->globalTrack()){
        origin.SetX(track->globalTrack()->helix().origin().x());
        origin.SetY(track->globalTrack()->helix().origin().y());
        origin.SetZ(track->globalTrack()->helix().origin().z());
        upcTrack->setCurvatureDipAnglePhase(track->globalTrack()->helix().curvature(), track->globalTrack()->helix().dipAngle(), track->globalTrack()->helix().phase() );
      }
      upcTrack->setOrigin(origin);
      upcTrack->setDcaXY( track->dcaGlobal().perp() );
      upcTrack->setDcaZ( track->dcaGlobal().z() );
      upcTrack->setCharge( track->charge() );
      upcTrack->setNhits( track->nHits() );
      upcTrack->setNhitsFit( track->nHitsFit() );
      upcTrack->setChi2( track->chi2() );
      upcTrack->setNhitsDEdx( track->nHitsDedx() );
      upcTrack->setDEdxSignal( track->dEdx() );
      upcTrack->setNSigmasTPC( StUPCTrack::kElectron, track->nSigmaElectron() );
      upcTrack->setNSigmasTPC( StUPCTrack::kPion, track->nSigmaPion() );
      upcTrack->setNSigmasTPC( StUPCTrack::kKaon, track->nSigmaKaon() );
      upcTrack->setNSigmasTPC( StUPCTrack::kProton, track->nSigmaProton() );
      upcTrack->setVertexId( ivtx );
      upcTrack->setFlag( StUPCTrack::kPrimary );
      if( emcProj ) {
        upcTrack->setFlag( StUPCTrack::kBemcProj );
        upcTrack->setBemcPtEtaPhi(emcPt, emcEta, emcPhi);
      }
      if( matchBemc ) {
        upcTrack->setFlag( StUPCTrack::kBemc );
        upcTrack->setBemcPtEtaPhi(emcPt, emcEta, emcPhi);
        upcTrack->setBemcClusterId(clsId);
        upcTrack->setBemcHitE(hitE);
      }
      if( matchTof ) {
        upcTrack->setFlag( StUPCTrack::kTof );
        upcTrack->setTofTime( tofPid.timeOfFlight() );
        upcTrack->setTofPathLength( tofPid.pathLength() );
      }

    }//tracks loop

    //if( nSelTracks <= 0 ) continue; //no selected tracks for this vertex

    //position of current primary vertex
    StThreeVectorF vtxPos = evt->primaryVertexPosition();
    StThreeVectorF vtxPosErr = evt->primaryVertexErrors();

    //write vertex to output UPC event
    StUPCVertex *upcVtx = mUPCEvent->addVertex();
    upcVtx->setPosX( vtxPos.x() );
    upcVtx->setPosY( vtxPos.y() );
    upcVtx->setPosZ( vtxPos.z() );
    upcVtx->setErrX( vtxPosErr.x() );
    upcVtx->setErrY( vtxPosErr.y() );
    upcVtx->setErrZ( vtxPosErr.z() );
    upcVtx->setNPrimaryTracks( trkArray->GetEntriesFast() );
    upcVtx->setNTracksUsed( mMuDst->primaryVertex()->nTracksUsed() );
    upcVtx->setId( ivtx );

  }//vertex loop


  // TOF hit loop
  for(UInt_t ihit=0; ihit < mMuDst->numberOfBTofHit(); ihit++) {
    // get hit from muDst
    StMuBTofHit  *tofhit = mMuDst->btofHit(ihit);
    if( !tofhit ) continue;
    //put hit in clones array in UPC event
    StUPCTofHit *upcHit = mUPCEvent->addHit();
    upcHit->setTray(tofhit->tray());
    upcHit->setModule(tofhit->module());
    upcHit->setCell(tofhit->cell());
    upcHit->setLeadingEdgeTime(tofhit->leadingEdgeTime());
    upcHit->setTrailingEdgeTime(tofhit->trailingEdgeTime());

  }


  //Roman Pot data
  if( mRPUtil ) mRPUtil->processEvent(mRPEvent, mMuDst);

  //write BEMC clusters
  mBemcUtil->writeBEMC(mUPCEvent);

  //event accepted to write to output tree
  mCounter->Fill( kWritten ); // events written to output tree
  mUPCTree->Fill();

  return kStOk;

}//Make

//_____________________________________________________________________________
Bool_t StUPCFilterMaker::runMC() {

  TClonesArray *muMcVtx = mMuDst->mcArray(0);
  TClonesArray *muMcTracks = mMuDst->mcArray(1);
  if( !muMcTracks ) return kFALSE;
  if( !muMcVtx ) return kFALSE;

  TDatabasePDG *pdgdat = TDatabasePDG::Instance();

  //vertex map from id to position in clones array
  map<Int_t, Int_t> vmap;
  //MC vertex loop
  for(Int_t ivtx=0; ivtx<muMcVtx->GetEntries(); ivtx++) {
    StMuMcVertex *vtx = dynamic_cast<StMuMcVertex*>(muMcVtx->At(ivtx));
    vmap[vtx->Id()] = ivtx;
  }//MC vertex loop

  //mc tracks loop
  for(Int_t i=0; i<muMcTracks->GetEntries(); i++) {
    StMuMcTrack *mcTrk = dynamic_cast<StMuMcTrack*>( muMcTracks->At(i) );
    if(!mcTrk) continue;

    const StThreeVectorF pxyz = mcTrk->Pxyz();
    Float_t px = pxyz.x();
    Float_t py = pxyz.y();
    Float_t pz = pxyz.z();
    Float_t energy = mcTrk->E();

    TParticle *part = mUPCEvent->addMCParticle();
    if(!part) return kFALSE;

    part->SetMomentum(px, py, pz, energy);
    part->SetPdgCode( pdgdat->ConvertGeant3ToPdg( mcTrk->GePid() ) );

    //MC vertex
    StMuMcVertex *vtx = dynamic_cast<StMuMcVertex*>(muMcVtx->At( vmap[mcTrk->IdVx()] ));
    const StThreeVectorF vxyz = vtx->XyzV();

    part->SetProductionVertex(vxyz.x(), vxyz.y(), vxyz.z(), 0.);
    //set original vertex id
    part->SetFirstMother(mcTrk->IdVx());

  }//mc tracks loop

  return kTRUE;

}//runMC

//_____________________________________________________________________________
Int_t StUPCFilterMaker::Finish()
{
  //called at the end

  //write the output file
  mOutFile->cd();

  mUPCTree->Write();
  mHistList->Write("HistList", TObject::kSingleKey);

  mOutFile->Close();

  return kStOk;

}//Finish
















































