#ifndef StUPCFilterBemcUtil_h
#define StUPCFilterBemcUtil_h

class StEmcPosition;
class StEmcGeom;
class StMuDst;
class StMuTrack;
class StUPCEvent;

class StUPCFilterBemcUtil {

public:

  StUPCFilterBemcUtil();
  virtual ~StUPCFilterBemcUtil();

  void clear();
  void setMagField(Double_t mfield) {mMagField = mfield;}
  Bool_t processEvent(StMuDst *muDst, StUPCEvent *upcEvt);
  Short_t matchBEMC(StMuTrack *track, Double_t &emcPhi, Double_t &emcEta, Double_t &emcPt, Bool_t &proj, UInt_t &clsId, Float_t &hitE);
  void writeBEMC(StUPCEvent *upcEvt);

private:

  StUPCFilterBemcUtil(const StUPCFilterBemcUtil &o); // not implemented
  StUPCFilterBemcUtil &operator=(const StUPCFilterBemcUtil &o); // not implemented

  Double_t mMagField; // magnetic field
  StEmcPosition *mEmcPos; // BEMC position for track projection
  StEmcGeom *mEmcGeom; // BEMC geometry for track projection

  //structure and vector for BEMC raw hits, used for track-BEMC matching
  struct emcHits {
    emcHits(): hitSoftId(0), clsId(0), hitE(0) {}
    Int_t hitSoftId;
    UInt_t clsId;
    Float_t hitE;
  };
  emcHits mEmcHit; // object of hits structure

  vector<emcHits> *mVecEmcHits; // vector of all BEMC raw hits for track-BEMC matching

  //BEMC clusters, structure and vector for writing clusters to UPC trees,
  //it allows to write only clusters matched to a track
  struct emcCluster {
    emcCluster(): clsEta(0), clsPhi(0), clsSigmaEta(0), clsSigmaPhi(0), clsE(0), isMatched(0),
      clsHT(0), clsHTsoftID(0) {}
    Float_t clsEta;
    Float_t clsPhi;
    Float_t clsSigmaEta;
    Float_t clsSigmaPhi;
    Float_t clsE;
    Bool_t isMatched;
    Float_t clsHT;
    Int_t clsHTsoftID;
  };
  emcCluster mEmcCluster; // object of cluster structure

  vector<emcCluster> *mVecEmcCluster; // vector of all BEMC clusters in event

};

#endif

