#ifndef StUPCV0_hh
#define StUPCV0_hh

/* **************************************************
 *  Generic class calculating and storing pairs in StUPCV0 analysis
 *  It is based on StHFPair code ported to StUPCDst format and traying
 *  to do the analysis when there is no good condidate for primary vertex  
 *  Allows to combine:
 *  - two particles, using
 *      StUPCV0(StUPCTrack const * particle1, StUPCTrack const * particle2, ...
 ***************************************************
 *  Initial Authors: 
 *            Leszek Adamczyk (leszek.adamczyk@agh.edu.pl)
 *
 *  ** Code Maintainer 
 *
 * **************************************************
 */

#include "TObject.h"
#include "TClonesArray.h"

#include "TVector3.h"
#include "TLorentzVector.h"

class StUPCTrack;
class StPicoPhysicalHelix;

class StUPCV0 : public TObject
{
 public:
  StUPCV0();
  StUPCV0(StUPCV0 const *);

  StUPCV0(StUPCTrack const * particle1, StUPCTrack const * particle2, 
	   float p1MassHypo, float p2MassHypo,
	   unsigned short p1Idx, unsigned short p2Idx,
	   TVector3 const & vtx, double * beamLine, float bField, bool useStraightLine = true);

  ~StUPCV0() {;}
  

  TLorentzVector const & lorentzVector() const;
  TVector3 const & decayVertex() const;    // V0 decay vertex
  TVector3 const & prodVertexHypo() const; // V0 production vertex 
  TVector3 const & prodPlane() const;
  float m()    const;
  float pt()   const;
  float eta()  const;
  float phi()  const;
  float pointingAngle() const;       // Pointing angle using input vertex (vtx)
  float pointingAngleHypo() const;   // Pointing angle using reconstructed V0 production vertex 
  float pointingAngle(TVector3 const & vtx2) const;
  float decayLength() const; // decay Length using input vertex (vtx)
  float decayLengthHypo() const; // decay Length using reconstructed V0 production vertex
  float decayLength(TVector3 const & vtx2) const;
  float DCABeamLine() const; // DCA between V0 straight-line and beam-line
  float particle1Dca() const;
  float particle1Dca(StPicoPhysicalHelix  p1Helix, TVector3 const & vtx2, float bField) const;
  float particle2Dca() const;
  float particle2Dca(StPicoPhysicalHelix  p1Helix, TVector3 const & vtx2, float bField) const;
  unsigned short particle1Idx() const;
  unsigned short particle2Idx() const;
  float dcaDaughters() const;
  float cosThetaStar() const;
  float thetaProdPlane() const;
  float v0x() const;
  float v0y() const;
  float v0z() const;
  float px() const;
  float py() const;
  float pz() const;
  float DcaToPrimaryVertex() const;
  float alphaAP() const;
  float ptAP() const;

 private:
  StUPCV0(StUPCV0 const &);
  StUPCV0& operator=(StUPCV0 const &);
  TLorentzVector mLorentzVector; 
  TVector3   mDecayVertex;

  TVector3 mProdPlane; //production plane vector: unity vector perpendicular to production plane

  float mDCABeamLine;
  TVector3 mProdVertexHypo;
  float mPointingAngleHypo;
  float mDecayLengthHypo;
  float mPointingAngle;
  float mDecayLength;
  float mParticle1Dca;
  float mParticle2Dca;
  float mDcaToPrimaryVertex;

  unsigned short  mParticle1Idx; // index of track in StUPCDstEvent 
  unsigned short  mParticle2Idx; // index of track in StUPCDstEvent for particle, idx in tertiary vertex array for pair 

  float mDcaDaughters;
  float mCosThetaStar;
  float mThetaProdPlane; //anlgle between momentum of particle1 and mProdPlane vector in mother rest frame
  float mAlphaAP;
  float mPtAP;
  ClassDef(StUPCV0,2)
};

inline TLorentzVector const & StUPCV0::lorentzVector() const { return mLorentzVector;}
inline float StUPCV0::m()    const { return mLorentzVector.M();}
inline float StUPCV0::pt()   const { return mLorentzVector.Perp();}
inline float StUPCV0::eta()  const { return mLorentzVector.PseudoRapidity();}
inline float StUPCV0::phi()  const { return mLorentzVector.Phi();}
inline float StUPCV0::px()   const { return mLorentzVector.Px();}
inline float StUPCV0::py()   const { return mLorentzVector.Py();}
inline float StUPCV0::pz()   const { return mLorentzVector.Pz();}
inline float StUPCV0::pointingAngle() const { return mPointingAngle;}
inline float StUPCV0::decayLength()   const { return mDecayLength;}
inline TVector3 const &  StUPCV0::prodVertexHypo() const { return mProdVertexHypo;}
inline float StUPCV0::pointingAngleHypo() const { return mPointingAngleHypo;}
inline float StUPCV0::decayLengthHypo()   const { return mDecayLengthHypo;}
inline float StUPCV0::DCABeamLine()   const { return mDCABeamLine;}
inline float StUPCV0::particle1Dca()  const { return mParticle1Dca;}
inline float StUPCV0::particle2Dca()  const { return mParticle2Dca;}
inline unsigned short StUPCV0::particle1Idx() const { return mParticle1Idx;}
inline unsigned short StUPCV0::particle2Idx() const { return mParticle2Idx;}
inline float StUPCV0::dcaDaughters() const { return mDcaDaughters;}
inline float StUPCV0::cosThetaStar() const { return mCosThetaStar;}
inline float StUPCV0::thetaProdPlane() const { return mThetaProdPlane;}
inline TVector3 const & StUPCV0::prodPlane() const { return mProdPlane;}
inline TVector3 const & StUPCV0::decayVertex() const { return mDecayVertex;}
inline float StUPCV0::v0x() const { return mDecayVertex.x();}
inline float StUPCV0::v0y() const { return mDecayVertex.y();}
inline float StUPCV0::v0z() const { return mDecayVertex.z();}
inline float StUPCV0::DcaToPrimaryVertex() const { return mDcaToPrimaryVertex; }
inline float StUPCV0::alphaAP() const { return mAlphaAP; }
inline float StUPCV0::ptAP() const { return mPtAP; }
#endif

