#ifndef StUPCTrack_h
#define StUPCTrack_h

//_____________________________________________________________________________
//    Class for UPC data
//    Author: Jaroslav Adam
//_____________________________________________________________________________

class TLorentzVector;
class StUPCEvent;
class StUPCBemcCluster;
class StUPCVertex;

#include "TObject.h"
#include "TVector3.h"

class StUPCTrack: public TObject
{

public:

  StUPCTrack();
  virtual ~StUPCTrack();

  //public types
  //flags for mFlags, up to 8 (uchar)
  enum Flag{kBemc=0, kTof, kBemcProj, kPrimary};
  //particles for dE/dx, 4 particle species
  enum Part{kElectron=0, kPion, kKaon, kProton};

  void Clear(Option_t * /*option*/ ="");

  //setters
  void setFlag(Flag flg);

  void setPtEtaPhi(Double_t pt, Double_t eta, Double_t phi) {mPt = pt; mEta = eta; mPhi = phi;}

  void setCurvatureDipAnglePhase(Double_t curvature, Double_t dipAngle, Double_t phase) {mCurvature = curvature; mDipAngle = dipAngle; mPhase = phase;}
  void setOrigin(const TVector3& origin) { mOrigin = origin; }

  void setDcaXY(Float_t dca) { mDcaXY = dca; }
  void setDcaZ(Float_t dca) { mDcaZ = dca; }

  void setCharge(Short_t ch) { mCharge = ch; }

  void setNhits(UShort_t nh) { mNhits = nh; }
  void setNhitsFit(UShort_t nf) { mNhitsFit = nf; }

  void setChi2(Double_t chi2) { mChi2 = chi2; }

  void setNhitsDEdx(UShort_t nh) { mNhitsDEdx = nh; }
  void setDEdxSignal(Double_t de) { mDEdxSignal = de; }
  void setNSigmasTPC(Part pt, Float_t nsig) { mNSigmasTPC[pt] = nsig; }

  void setBemcPtEtaPhi(Double_t pt, Double_t eta, Double_t phi) {mBemcPt = pt; mBemcEta = eta; mBemcPhi = phi;}
  void setBemcClusterId(UInt_t id) { mBemcClsId = id; }
  void setBemcHitE(Float_t he) { mBemcHitE = he; }

  void setTofTime(Float_t time) { mTofTime = time; }
  void setTofPathLength(Float_t len) { mTofPathLength = len; }

  void setVertexId(UInt_t id) { mVtxId = id; }

  void setEvent(StUPCEvent *evt) { mEvt = evt; }

  //getters
  Bool_t getFlag(Flag flg) const;

  void getPtEtaPhi(Double_t &pt, Double_t &eta, Double_t &phi) const {pt = mPt; eta = mEta; phi = mPhi;}
  Double_t getPt() const { return mPt; }
  Double_t getEta() const { return mEta; }
  Double_t getPhi() const { return mPhi; }
  void getLorentzVector(TLorentzVector &lvec, Double_t mass) const;
  void getMomentum(TVector3 &vec) const;

  Double_t getCurvature() const { return mCurvature; }
  Double_t getDipAngle() const { return mDipAngle; }
  Double_t getPhase() const { return mPhase; }
  TVector3 getOrigin() const { return mOrigin; }

  Float_t getDcaXY() const { return mDcaXY; }
  Float_t getDcaZ() const { return mDcaZ; }

  Short_t getCharge() const { return mCharge; }

  UShort_t getNhits() const { return mNhits; }
  UShort_t getNhitsFit() const { return mNhitsFit; }

  Double_t getChi2() const { return mChi2; }

  UShort_t getNhitsDEdx() const { return mNhitsDEdx; }
  Double_t getDEdxSignal() const { return mDEdxSignal; }
  Float_t getNSigmasTPC(Part pt) const { return mNSigmasTPC[pt]; }
  Float_t getNSigmasTPCElectron() const { return mNSigmasTPC[kElectron]; }
  Float_t getNSigmasTPCPion() const { return mNSigmasTPC[kPion]; }
  Float_t getNSigmasTPCKaon() const { return mNSigmasTPC[kKaon]; }
  Float_t getNSigmasTPCProton() const { return mNSigmasTPC[kProton]; }

  void getBemcPtEtaPhi(Double_t &pt, Double_t &eta, Double_t &phi) const {pt = mBemcPt; eta = mBemcEta; phi = mBemcPhi;}
  Double_t getBemcPt() const { return mBemcPt; }
  Double_t getBemcEta() const { return mBemcEta; }
  Double_t getBemcPhi() const { return mBemcPhi; }
  Double_t getBemcPmag() const;
  UInt_t getBemcClusterId() const { return mBemcClsId; }
  Float_t getBemcHitE() const { return mBemcHitE; }
  StUPCBemcCluster *getBemcCluster() const;
  void getBemcLorentzVector(TLorentzVector &blvec, Double_t mass) const;

  Float_t getTofTime() const { return mTofTime; }
  Float_t getTofPathLength() const { return mTofPathLength; }

  UInt_t getVertexId() const { return mVtxId; }
  StUPCVertex *getVertex() const;

  StUPCEvent *getEvent() const { return mEvt; }

private:

  StUPCTrack(const StUPCTrack &o); //not implemented
  StUPCTrack &operator=(const StUPCTrack &o); //not implemented

  UChar_t mFlags; // track flags

  Double32_t mPt; // pT at point of dca to primary vertex
  Double32_t mEta; // pseudorapidity at point of dca to primary vertex
  Double32_t mPhi; // phi at point of dca to primary vertex

  Double32_t mCurvature; // curvature of track from StHelix
  Double32_t mDipAngle; // dip angle of track from StHelix
  Double32_t mPhase; // phase of track from StHelix
  TVector3 mOrigin; // origin of track from StHelix

  Float_t mDcaXY; // perpendicular dca to primary vertex of associated global track
  Float_t mDcaZ; // longitudinal dca to primary vertex of associated global track

  Short_t mCharge; // track electrical charge

  UShort_t mNhits; // total number of hits on track
  UShort_t mNhitsFit; // number of hits used in fit

  Double32_t mChi2; // chi2 of fit

  UShort_t mNhitsDEdx; // number of hits used for dE/dx measurement
  Double32_t mDEdxSignal; // measured dE/dx value
  static const Int_t mNpart = 4; // number of particle species with dE/dx identification
  Float16_t mNSigmasTPC[mNpart]; // dE/dx n sigmas for particle species

  Double32_t mBemcPt; // pT at BEMC radius
  Double32_t mBemcEta; // pseudorapidity at BEMC radius
  Double32_t mBemcPhi; // phi at BEMC radius
  UInt_t mBemcClsId; // id of BEMC cluster which track is matched to
  Float_t mBemcHitE; // energy of matched BEMC cluster

  Float_t mTofTime; // time of flight by TOF
  Float_t mTofPathLength; // path lenght from TOF

  UInt_t mVtxId; // ID of primary vertex associated with track

  StUPCEvent *mEvt; //! pointer to current event, local use only

  ClassDef(StUPCTrack, 2)

};

#endif


















