#include <limits>
#include <cmath>
#include <iostream>

#include "StUPCV0.h"

#include "StPicoPhysicalHelix.h"
#include "SystemOfUnits.h"
#include "StUPCTrack.h"

using namespace std;

ClassImp(StUPCV0)

// _________________________________________________________
StUPCV0::StUPCV0(): mLorentzVector(TLorentzVector()), mDecayVertex(TVector3()),
  mPointingAngle(std::numeric_limits<float>::quiet_NaN()), mDecayLength(std::numeric_limits<float>::quiet_NaN()),
  mParticle1Dca(std::numeric_limits<float>::quiet_NaN()), mParticle2Dca(std::numeric_limits<float>::quiet_NaN()),
  mParticle1Idx(std::numeric_limits<unsigned short>::max()), mParticle2Idx(std::numeric_limits<unsigned short>::max()),
  mDcaDaughters(std::numeric_limits<float>::max()), mCosThetaStar(std::numeric_limits<float>::quiet_NaN()),
  mThetaProdPlane(std::numeric_limits<float>::quiet_NaN()), mProdPlane(TVector3()) {
}

// _________________________________________________________
StUPCV0::StUPCV0(StUPCV0 const * t) : mLorentzVector(t->mLorentzVector), mDecayVertex(t->mDecayVertex),
   mPointingAngle(t->mPointingAngle), mDecayLength(t->mDecayLength),
   mParticle1Dca(t->mParticle1Dca), mParticle2Dca(t->mParticle2Dca),
   mParticle1Idx(t->mParticle1Idx), mParticle2Idx(t->mParticle2Idx),
   mDcaDaughters(t->mDcaDaughters), mCosThetaStar(t->mCosThetaStar),
   mThetaProdPlane(t->mThetaProdPlane), mProdPlane(t->mProdPlane) {
}

// _________________________________________________________
StUPCV0::StUPCV0(StUPCTrack const * const particle1, StUPCTrack const * const particle2,
		   float p1MassHypo, float p2MassHypo, unsigned short const p1Idx, unsigned short const p2Idx,
		   TVector3 const & vtx, double * beamLine, float const bField, bool const useStraightLine) : 
  mLorentzVector(TLorentzVector()), mDecayVertex(TVector3()),
  mPointingAngle(std::numeric_limits<float>::quiet_NaN()), mDecayLength(std::numeric_limits<float>::quiet_NaN()),
  mParticle1Dca(std::numeric_limits<float>::quiet_NaN()), mParticle2Dca(std::numeric_limits<float>::quiet_NaN()),
  mParticle1Idx(p1Idx), mParticle2Idx(p2Idx),
  mDcaDaughters(std::numeric_limits<float>::max()), mCosThetaStar(std::numeric_limits<float>::quiet_NaN()) {
  // -- Create pair out of 2 tracks
  //     prefixes code:
  //      p1 means particle 1
  //      p2 means particle 2
  //      pair means particle1-particle2  pair

// StUPCTrack do not have id() method :: to do  
/*
  if ((!particle1 || !particle2) || (particle1->id() == particle2->id())) {
    mParticle1Idx = std::numeric_limits<unsigned short>::max();
    mParticle2Idx = std::numeric_limits<unsigned short>::max();
    return;
  }
*/

// StUPCTrack do not have helix() method :: to do
/*
  StPicoPhysicalHelix p1Helix = particle1->helix(bField); 
  StPicoPhysicalHelix p2Helix = particle2->helix(bField); //bFiled not in kilogauss - is properly computed inside helix(double B) function in StUPCTrack.h
*/
 
  StPicoPhysicalHelix p1Helix = StPicoPhysicalHelix (particle1->getCurvature(), 
                                                     particle1->getDipAngle(), 
                                                     particle1->getPhase(), 
                                                     particle1->getOrigin(),
                                                     particle1->getCharge());

  StPicoPhysicalHelix p2Helix =	StPicoPhysicalHelix (particle2->getCurvature(), 
                                                     particle2->getDipAngle(), 
                                                     particle2->getPhase(), 
                                                     particle2->getOrigin(),
                                                     particle2->getCharge());
  p1Helix.moveOrigin(p1Helix.pathLength(vtx,false));
  p2Helix.moveOrigin(p2Helix.pathLength(vtx,false));

// -- use straight lines approximation to get point of DCA of particle1-particle2 pair
// bField is in kilogauss

  TVector3 const p1Mom = p1Helix.momentum(bField * kilogauss);
  TVector3 const p2Mom = p2Helix.momentum(bField * kilogauss);

  StPicoPhysicalHelix const p1StraightLine(p1Mom, p1Helix.origin(), 0, particle1->getCharge());
  StPicoPhysicalHelix const p2StraightLine(p2Mom, p2Helix.origin(), 0, particle2->getCharge());

  pair<double, double> const ss = (useStraightLine) ? p1StraightLine.pathLengths(p2StraightLine) : p1Helix.pathLengths(p2Helix);
  TVector3 const p1AtDcaToP2 = (useStraightLine) ? p1StraightLine.at(ss.first) : p1Helix.at(ss.first);
  TVector3 const p2AtDcaToP1 = (useStraightLine) ? p2StraightLine.at(ss.second) : p2Helix.at(ss.second);

  // -- calculate DCA of particle1 to particle2 at their DCA
  mDcaDaughters = (p1AtDcaToP2 - p2AtDcaToP1).Mag();

  // -- calculate Lorentz vector of particle1-particle2 pair

  TVector3 const p1MomAtDca = p1Helix.momentumAt(ss.first,  bField *kilogauss);
  TVector3 const p2MomAtDca = p2Helix.momentumAt(ss.second, bField *kilogauss );

  TLorentzVector const p1FourMom(p1MomAtDca, sqrt(p1MomAtDca.Mag2() + p1MassHypo*p1MassHypo));
  TLorentzVector const p2FourMom(p2MomAtDca, sqrt(p2MomAtDca.Mag2() + p2MassHypo*p2MassHypo));

  mLorentzVector = p1FourMom + p2FourMom; // K0 momentum

  // -- calculate cosThetaStar
  TLorentzVector const pairFourMomReverse(-mLorentzVector.Px(), -mLorentzVector.Py(), -mLorentzVector.Pz(), mLorentzVector.E());
  TLorentzVector FourMomStar;
  if ( particle1->getCharge() == 1 ) {
    FourMomStar = p1FourMom; }
  else { 
    FourMomStar = p2FourMom; }
  FourMomStar.Boost(pairFourMomReverse.BoostVector());  
  mCosThetaStar = std::cos(FourMomStar.Vect().Angle(mLorentzVector.Vect()));

//  TVector3 beamVector(0.,0.,1.); //unity vector along the beam axis
  TVector3 beamVector(beamLine[2],beamLine[3],1.); //unity vector along the beam axis with beamLine3D

  TVector3 mProdPlane_work = beamVector.Cross(mLorentzVector.Vect());
  mProdPlane = ( mProdPlane_work )*(1./mProdPlane_work.Mag() ); //unity normal vector to production plane

  mThetaProdPlane = mProdPlane.Angle(FourMomStar.Vect());
  
  // -- calculate decay vertex (secondary or tertiary) 
  mDecayVertex = (p1AtDcaToP2 + p2AtDcaToP1) * 0.5 ;
  mDCABeamLine = abs(mProdPlane*mDecayVertex);

  float t = mDecayVertex.Cross(beamVector).Mag()/mProdPlane_work.Mag();
  mProdVertexHypo = mDecayVertex - t*mLorentzVector.Vect();
  mProdVertexHypo.SetX(beamLine[0]+beamLine[2]* mProdVertexHypo.Z());
  mProdVertexHypo.SetY(beamLine[1]+beamLine[3]* mProdVertexHypo.Z());

  // -- calculate pointing angle and decay length with respect to primary vertex 
  //    if decay vertex is a tertiary vertex
  //    -> only rough estimate -> needs to be updated after secondary vertex is found
  TVector3 const vtxToV0 = mDecayVertex - vtx;
  TVector3 const ProdvtxToV0 = mDecayVertex - mProdVertexHypo;
  mPointingAngle = vtxToV0.Angle(mLorentzVector.Vect());
  mPointingAngleHypo = cos(ProdvtxToV0.Angle(mLorentzVector.Vect()));
  mDecayLength = vtxToV0.Mag();
  mDecayLengthHypo = ProdvtxToV0.Mag();
  mDcaToPrimaryVertex = mDecayLength*sin(mPointingAngle); // sine law: DcaToPrimaryVertex/sin(pointingAngle) = decayLength/sin(90°)

  // -- calculate DCA of tracks to primary vertex
  //    if decay vertex is a tertiary vertex
  //    -> only rough estimate -> needs to be updated after secondary vertex is found
  mParticle1Dca = (p1Helix.origin() - vtx).Mag();
  mParticle2Dca = (p2Helix.origin() - vtx).Mag();
}
// _________________________________________________________
float StUPCV0::pointingAngle(TVector3 const & vtx2) const{
  // -- Overloaded function recalculates pointing angle given secondary vertex
  TVector3 const vtx2ToTertiary(mDecayVertex - vtx2);
  float const nPointingAngle = vtx2ToTertiary.Angle(mLorentzVector.Vect());
  return nPointingAngle;
}
// _________________________________________________________
float StUPCV0::decayLength(TVector3 const & vtx2) const{
  // -- Overloaded function recalculates decayLength given secondary vertex
  TVector3 const vtx2ToTertiary(mDecayVertex - vtx2); 
  float const nDecayLength = vtx2ToTertiary.Mag();  
  return nDecayLength;
}
// _________________________________________________________
float StUPCV0::particle1Dca(StPicoPhysicalHelix  p1Helix, TVector3 const & vtx2, float bField) const{
  // -- Overloaded function recalculates daughter dca 2 updated vertex
  // StPicoPhysicalHelix p1Helix = p1track->helix(bField);
  // -- move origins of helices to the primary vertex origin
  p1Helix.moveOrigin(p1Helix.pathLength(vtx2));

  float const nParticle1Dca = (p1Helix.origin() - vtx2).Mag();
  return nParticle1Dca;
}
// _________________________________________________________
float StUPCV0::particle2Dca(StPicoPhysicalHelix  p2Helix, TVector3 const & vtx2, float bField) const{
  // -- Overloaded function recalculates daughter dca 2 updated vertex
  // StPicoPhysicalHelix p2Helix = p2track->helix(bField);
  // -- move origins of helices to the primary vertex origin
  p2Helix.moveOrigin(p2Helix.pathLength(vtx2));

  float const nParticle2Dca = (p2Helix.origin()  - vtx2).Mag();
  return nParticle2Dca;
}



