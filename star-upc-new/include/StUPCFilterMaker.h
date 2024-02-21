#ifndef StUPCFilterMaker_h
#define StUPCFilterMaker_h

#include "StMaker.h"

class StMuDstMaker;
class StMuDst;
class StUPCEvent;
class StTriggerData;
class StEmcPosition;
class StEmcGeom;
class StMuTrack;
class StUPCFilterTrgUtil;
class StUPCFilterBemcUtil;
class StUPCFilterRPUtil;
class StRPEvent;

class StUPCFilterMaker: public StMaker {

public:

  StUPCFilterMaker(StMuDstMaker *maker, string outnam="StUPC.root");
  ~StUPCFilterMaker();

  void setIsMC(Int_t mc) { mIsMC = mc; }
  void setMakeRPEvent(Bool_t makerp) { mMakeRP = makerp; }

  void addTriggerId(UInt_t id, Int_t rmin=-1e8, Int_t rmax=1e8);

  Int_t Init();   //called at the beginning of the analysis
  Int_t Make();   //called for each event
  Int_t Finish(); //called at the end

private:

  StUPCFilterMaker(const StUPCFilterMaker &o); // not implemented
  StUPCFilterMaker &operator=(const StUPCFilterMaker &o); // not implemented

  Bool_t runMC();

  StMuDstMaker *mMaker;  //StMuDstMaker provided to the constructor
  StMuDst *mMuDst; // input muDst data

  Int_t mIsMC; // MC or data

  string mOutName;  // name of the output file
  TFile *mOutFile;  // output file

  TList *mHistList; // list of output histograms

  TH1I *mCounter; // analysis counter
  enum EvtCount{ kAna=1, kTrg, kTrgDat, kBemc, kWritten, kMaxCnt };
  TH1I *mErrCounter; // error counter
  enum ErrCount{ kErrNoEvt=1, kErrFillMsc, kNoTrgDat, kNoMC, kMaxErrCnt };

  StUPCEvent *mUPCEvent; // output UPC event
  TTree *mUPCTree; // output tree

  //table of trigger IDs
  vector<UInt_t> mTrgIDs; // trigger IDs
  vector<Int_t> mTrgRanLo; // lower limit for run range for a given trigger ID
  vector<Int_t> mTrgRanHi; // upper limit for run range for a given trigger ID

  StUPCFilterTrgUtil *mTrgUtil;
  StUPCFilterBemcUtil *mBemcUtil; //utility class for BEMC matching

  Bool_t mMakeRP; // flag to configure for RP event
  StUPCFilterRPUtil *mRPUtil; // utility class for Roman Pots
  StRPEvent *mRPEvent; // output RP event

  ClassDef(StUPCFilterMaker, 3);

};






#endif

















