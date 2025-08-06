#ifndef Efficiencies_hh
#define Efficiencies_hh

//#define RUN_OVER_EMBEDDED_MC
//#define TEST_GENEX_PION_EFFICIENCIES

#include <limits>
#include <iostream>
#include <vector>
#include <cmath>
#include <TF1.h>
#include <TH3.h>
#include "Util.hh"
#include "Parameters.hh"

class TVector3;
class TEfficiency;
class TGraphAsymmErrors;
class TH2Poly;

class Efficiencies{
  public:
    static Efficiencies* instance( const Parameters * = nullptr, TString = TString(""), Int_t = 0 );
    
    TVector3 correctedMomentumVec(UInt_t, UInt_t, Double_t, const TVector3 &) const;
    Double_t tpcRecoEff(UInt_t, UInt_t, Double_t) const;
    Double_t tofEff(UInt_t, UInt_t, Double_t) const;
    Double_t tpcRecoEff3D(UInt_t, UInt_t, UInt_t, Double_t, Double_t, Double_t, Bool_t = false /*with or without reco cuts*/) const;
    Double_t tpcRecoEff3D_phi(UInt_t, UInt_t, UInt_t, Double_t, Double_t, Double_t) const;
    Double_t tofEff3D(UInt_t, UInt_t, Double_t, Double_t, Double_t, Bool_t = false /*with or without reco cuts*/) const;
    Double_t tofEff3D_phi(UInt_t, UInt_t, Double_t, Double_t, Double_t) const;
    Double_t tofEff2DCorrection( UInt_t, Double_t, Double_t) const;
    Double_t tofEff2DExtraCorrection( UInt_t, Double_t, Double_t) const;
    Double_t tofEff2DSystematicError( UInt_t, Double_t, Double_t) const;
    Double_t tofEffDifferenceHftMC3D( UInt_t, Double_t, Double_t, Double_t) const;
    inline TH3F* tpcRecoEff3DHist(UInt_t a, UInt_t b, UInt_t rr = Util::RUN_RANGE_2, Bool_t c = false /*with or without reco cuts*/) const { return c ? mhTpcRecoEff3D_withRecoCuts[a][b] : mhTpcRecoEff3D[a][b][rr]; };
    inline TH3F* tofEff3DHist(UInt_t a, UInt_t b, Bool_t c = false /*with or without reco cuts*/) const { return c ? mhTofEff3D_withRecoCuts[a][b] : mhTofEff3D[a][b]; };
    inline TH2F* tofEffCorrection2DHist(UInt_t a) const { return mhTofEffCorrection2D_ptVsEta[a]; };
    Double_t deadMatEffectTpcRecoEff3D(UInt_t, UInt_t, Double_t, Double_t, Double_t) const;
    Double_t getLowHighPileUpSystematics(Double_t, /*Double_t, */UInt_t, UInt_t) const;
    
    Double_t tpcEffMultiplicativeCorrectionFactor(UInt_t, Double_t) const;
    
    Double_t purityPid(UInt_t, Double_t) const;
    inline TEfficiency* pidEff(UInt_t a, UInt_t b) const{ return (a<Util::nDefinedParticles && b<Util::nDefinedParticles) ? mPidEfficiency[a][b] : nullptr; };
    inline TEfficiency* pidEff_Pt(UInt_t a, UInt_t b) const{ return (a<Util::nDefinedParticles && b<Util::nDefinedParticles) ? mPidEfficiency_Pt[a][b] : nullptr; };
    Double_t pidEff(UInt_t, Double_t, Double_t) const;
    
    inline TEfficiency* ptMissCutEff() const{ return mPtMissCutEff; };
    
    Double_t trackEffZVtxPxPy(UInt_t, Double_t, Double_t, Double_t, Bool_t) const;
    inline TH2D *trackEffPxPyHist(UInt_t br) const { return br<Util::nBranches ? mRpTrackRecoEffiIntegratedOverZVtxPerBranch[br] : nullptr; };
    inline TH3F *RpTrackRecoEffiPerBranch(UInt_t br) const {  return br<Util::nBranches ? mRpTrackRecoEffiPerBranch[br] : nullptr; };
    inline TH3F *RpTrackRecoEffiPerBranch_FullRecoEfficiency(UInt_t br) const {  return br<Util::nBranches ? mRpTrackRecoEffiPerBranch_FullRecoEfficiency[br] : nullptr; };
    inline TH3F *RpTrackRecoEffiPerBranch_DeadMaterialVetoEff(UInt_t br) const {  return br<Util::nBranches ? mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff[br] : nullptr; };
    inline TH1F *SimultaneousEastWestRpTrackLossProbability() const {  return mhSimultaneousEastWestRpTrackLossProbability; };
    Double_t rpTrackEffPxPySystematicError(UInt_t, Double_t, Double_t) const;
    Double_t rpDeadMatVetoEffPxPySystematicError(UInt_t, Double_t, Double_t) const;
    Double_t trackEffZVtxPxPy_FullRecoEfficiency(UInt_t, Double_t, Double_t, Double_t, Bool_t) const;
    Double_t trackEffZVtxPxPy_DeadMaterialVetoEff(UInt_t, Double_t, Double_t, Double_t, Bool_t) const;
    Double_t trackEffZVtxPxPy_DeadMaterialVetoEff_ProtonTrackRequired(UInt_t, Double_t, Double_t, Double_t, Bool_t) const;
    Double_t trackEffTypeZVtxPxPy(UInt_t, UInt_t, Double_t, Double_t, Double_t) const;
    Double_t zVertexCutEff(UInt_t, Double_t, Double_t, Double_t = 0.0, Double_t = 0.0) const;
    Double_t zVertexCutEff(Double_t, Double_t, Double_t = 0.0, Double_t = 50.0) const;
    Double_t onlineAndOfflineVetoEff(UInt_t, Double_t) const;
    Double_t onlineAndOfflineVetoEff(UInt_t, UInt_t) const;
    Double_t vertexFinderEff_InvMass(UInt_t, Double_t) const;
    Double_t vertexingEff(Double_t, Bool_t = true) const;
    inline Double_t tofClusterLimitEff() const{ return mTofClusterLimitEff; };
    Double_t goodRecoProtonProbability(UInt_t) const;
    Double_t simultaneousBothSidesProtonLossProbability(UInt_t) const;
    Double_t rpTrackRecoEffEastWestCorrelation(UInt_t) const;
    Double_t rpDeadMatVetoEffEastWestCorrelation(UInt_t) const;
    Double_t rpTotalEffEastWestCorrelation(UInt_t) const;
//     Double_t rpTotalEffEastWestCorrelation(UInt_t, Double_t) const;
    Double_t rpTotalEffEastWestCorrelation(UInt_t, Double_t, Double_t) const;
    Double_t tpcTotalEffPlusMinusCorrelation(Double_t) const;
    Double_t tpcRecoEffPlusMinusCorrelation(Double_t) const;
    Double_t rpTpcTotalEfficiencyCorrelation(Double_t) const;
    Double_t simultaneousBothSidesProtonLossProbability(UInt_t, Double_t) const;
    Double_t correctedDEdxMC(UInt_t, Double_t, Double_t) const;
    Double_t correctedNSigmaMC(UInt_t, UInt_t, Double_t, Double_t, Double_t, Double_t) const;
    Double_t meanDEdx(UInt_t, UInt_t, Double_t) const;
    Double_t sigmaDEdx(UInt_t, UInt_t, Double_t) const;
    Double_t migrationsAndFakesCorrectionCentralTracks(Double_t, Double_t) const;
    Double_t migrationsAndFakesCorrectionForwardProtons(Double_t, Double_t) const;
    Double_t migrationsAndFakesCorrectionForwardProtons(Int_t, Double_t, Double_t) const;
    Double_t migrationsAndFakesCorrectionForwardProtons_LimitedMandelstamT(Double_t, Double_t) const;
    Double_t migrationsAndFakesCorrectionForwardProtons_LimitedMandelstamT(Double_t) const;
    Double_t migrationsAndFakesCorrectionForwardProtons_LimitedMandelstamT(Int_t, Double_t) const;
    Double_t migrationsCorrectionWithinFiducialRegionForwardProtons(UInt_t, Double_t, Double_t) const;
    inline Double_t vertexFinderEff() const{ return mVertexFinderEff; };
    inline Double_t zVertexCutEff() const{ return mZVertexCutEff; };
    inline Double_t missingPtCutEff() const{ return mMissingPtCutEff; };
    Double_t missingPtCutEff(Double_t, Double_t, Bool_t = kTRUE) const;
    Double_t geomAcc(UInt_t, UInt_t, Double_t) const;
    Double_t deltaPhiCutEff(UInt_t, Double_t, Double_t) const;
    Double_t deltaPhiCutEff_GenEx(UInt_t, Double_t, Double_t) const;
    Double_t deltaPhiCutEff_DiMe(UInt_t, Double_t, Double_t) const;
    inline TH2F* deltaPhiCutEff(UInt_t bin) const{ return bin<Util::nDeltaPhiRanges ? mAcceptanceForDeltaPhiCuts_Vs_t1t2_TH2F[bin] : nullptr; }
    inline TH2F* deltaPhiCutEff_GenEx(UInt_t bin) const{ return bin<Util::nDeltaPhiRanges ? mAcceptanceForDeltaPhiCuts_Vs_t1t2_GenEx_TH2F[bin] : nullptr; }
    inline TH2F* deltaPhiCutEff_DiMe(UInt_t bin) const{ return bin<Util::nDeltaPhiRanges ? mAcceptanceForDeltaPhiCuts_Vs_t1t2_DiMe_TH2F[bin] : nullptr; }
    inline Double_t tofTrigEff() const{ return mTofTrigEff; };
    inline Double_t deltaZ0CutEff() const{ return mDeltaZ0CutEff; };
    Double_t deltaZ0CutEff(Double_t, Bool_t) const;
    inline TGraphAsymmErrors* TagAndProbeEfficiency_1TofTrack_Graph(UInt_t triggerId) const{ return triggerId < Util::nTriggers ? mTagAndProbeEfficiency_1TofTrack[triggerId] : nullptr; };
    inline TGraphAsymmErrors* TagAndProbeEfficiency_2TofTracks_Graph(UInt_t triggerId) const{ return triggerId < Util::nTriggers ? mTagAndProbeEfficiency_2TofTracks[triggerId] : nullptr; };
    inline TGraphAsymmErrors* TagAndProbe_VertexingEffVsD0_1TofTrack(UInt_t triggerId) const{ return triggerId < Util::nTriggers ? mTagAndProbe_VertexingEffVsD0_1TofTrack[triggerId] : nullptr; };
    inline TH2F* TagAndProbe_AttachingToVertexingEffVsDcaRVsDcaZ_1TofTrack(UInt_t triggerId) const{ return triggerId < Util::nTriggers ? mTagAndProbe_AttachingToVertexingEffVsDcaRVsDcaZ_1TofTrack[triggerId] : nullptr; };
    Double_t probability_TofL0MultAtLeast1VsClusterSize( unsigned int ) const;
    Double_t probability_TofL0MultAtLeast2VsClusterSize( unsigned int ) const;
    Double_t probability_TofL0MultAtLeast1VsNRecoTofHits_1TofVtx_1TofTrack( unsigned int ) const;
    Double_t probability_TofL0MultAtLeast1VsNRecoTofHits_1TofVtx_2TofTracks( unsigned int ) const;
    Double_t probability_TofL0MultAtLeast2VsNRecoTofHits_1TofVtx_1TofTrack( unsigned int ) const;
    Double_t probability_TofL0MultAtLeast2VsNRecoTofHits_1TofVtx_2TofTracks( unsigned int ) const;
    Double_t probability_TofL0MultAtLeast1VsTrackEta_1TofVtx_1TofTrack( double ) const;
    Double_t probability_TofL0MultAtLeast2VsTrackEta_1TofVtx_1TofTrack( double ) const;
    const TF1* tpcEffTF1(UInt_t, UInt_t) const;
    const TF1* tofEffTF1(UInt_t, UInt_t) const;
    TF1* tpcEff3DTF1(UInt_t, UInt_t, Double_t, Double_t) const;
    TF1* tofEff3DTF1(UInt_t, UInt_t, Double_t, Double_t) const;
    void multiplyTofEffPlateu(UInt_t, UInt_t, Double_t);
    void shiftTofEffPlateu(UInt_t, UInt_t, Double_t);
    void multiplyTpcRecoEffPlateu(UInt_t, UInt_t, Double_t);
    void shiftTpcEffPlateu(UInt_t, UInt_t, Double_t);
    void shiftTofEffInflation(UInt_t, UInt_t, Double_t);
    void shiftTpcRecoEffInflation(UInt_t, UInt_t, Double_t);
    void resetTofEffToDefault(UInt_t, UInt_t);
    void resetTpcRecoEffToDefault(UInt_t, UInt_t);

    void shiftTofEff3DPlateu(UInt_t, UInt_t, Double_t);
    void multiplyTofEff3DPlateu(UInt_t, UInt_t, Double_t);
    void shiftTpcRecoEff3DPlateu(UInt_t, UInt_t, Double_t);
    void multiplyTpcRecoEff3DPlateu(UInt_t, UInt_t, Double_t);
    void shiftTofEff3DInflation(UInt_t, UInt_t, Double_t);
    void shiftTpcRecoEff3DInflation(UInt_t, UInt_t, Double_t);
    void resetTofEff3DToDefault(UInt_t, UInt_t);
    void resetTpcRecoEff3DToDefault(UInt_t, UInt_t);
    
    inline void interpolate(Bool_t option){ mInterpolate = option; };
  
    enum GEOM_ACCEPTANCE_CORRECTIONS { S0m_PHASE_SPACE, D0m_PHASE_SPACE, S0mD0mMix_PHASE_SPACE, S0m_GENEX, D0m_GENEX, GENEX, S0m_DIME, D0m_DIME, DIME, nGeomAcceptances };
    
  protected:
    
    inline void setParametersObject( const Parameters *par ){ mParams = par; };
    
  private:
      
    enum ETA_BINS { PLUS_07_03_ETA, PLUS_03_0_ETA, MINUS_03_0_ETA, MINUS_07_03_ETA, nEtaBins };
    
    void loadGeometricalAcceptance( TString = TString(""));
    void loadEfficiencies3D( TString = TString(""), Int_t = 0 );
    Double_t tofEff3DForTF1(Double_t*, Double_t*) const;
    Double_t tpcEff3DForTF1(Double_t*, Double_t*) const;

    TF1* mTpcReconstructionEfficiency[Util::nSigns][Util::nDefinedParticles];
    TF1* mTofMatchingEfficiency[Util::nSigns][Util::nDefinedParticles];
    TF1* mFullVetosEfficiency[Util::nBranchesConfigurations];
    
    TF1* mTpcEffMultiplicativeCorrectionFactor;
    
    TH3F* mhTofEff3D[Util::nSigns][Util::nDefinedParticles];
    TH3F* mhTpcRecoEff3D[Util::nSigns][Util::nDefinedParticles][Util::nRunRanges];
    TH3F* mhTofEff3D_withRecoCuts[Util::nSigns][Util::nDefinedParticles];
    TH3F* mhTpcRecoEff3D_withRecoCuts[Util::nSigns][Util::nDefinedParticles];
    TH2F* mhTofEffCorrection2D_ptVsEta[Util::nDefinedParticles];
    TH3F* mhTofEffDifferenceHftMC3D_ptVsEtaVsZVtx[Util::nDefinedParticles];
    TH2D *mhTofEffExtraCorr2D_ptVsEta[Util::nDefinedParticles];
    TH2D *mhTofEffSystErr2D_ptVsEta[Util::nDefinedParticles];
    
    TH3F* mhTpcRecoEff3D_Phi[Util::nSigns][Util::nDefinedParticles][Util::nRunRanges];
    TH3F* mhTofEff3D_Phi[Util::nSigns][Util::nDefinedParticles];
    
    TH3F* mEnergyLossCorrection3D[Util::nSigns][Util::nDefinedParticles][20];
    TH3F* mhDeadMaterialEffectTpcRecoEff3D[Util::nSigns][Util::nDefinedParticles];
    TH1D* mhSystematicsPileUpTPC[Util::nRunRanges][Util::nSigns][2];
    
    TH1D* mhZVtxMeanVsFill;
    TH1D* mhZVtxSigmaVsFill;
    
    TH2Poly* mMigrationsAndFakesCorrection_PtEta;
    TH2F* mMigrationsAndFakesCorrection_PxPy;
    TH2F* mMigrationsAndFakesCorrection_PxPy_LimitedMandelstamT;
    TH1F* mMigrationsAndFakesCorrection_MandelstamT_LimitedMandelstamT;
    TH2F* mMigrationsWithinFiducialRegionCorrection_PxPy[Util::nSides];

    TH2F* mMigrationsAndFakesCorrection_PxPy_PerSide[Util::nSides];
    TH1F* mMigrationsAndFakesCorrection_MandelstamT_LimitedMandelstamT_PerBranch[Util::nBranches];
    
    TH3F* mRpTrackRecoEffiPerBranch[Util::nBranches];
    TH2D* mRpTrackRecoEffiIntegratedOverZVtxPerBranch[Util::nBranches];
    TH3F* mRpTrackRecoEffiPerBranch_FullRecoEfficiency[Util::nBranches];
    TH3F* mRpTrackRecoEffiPerBranchPerTrackType[Util::nBranches][2];
    TH3F *mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff[Util::nBranches];
    TH3F *mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff_ProtonTrackRequired[Util::nBranches];
    TH2F* mRpTrackRecoEffiCorrectionPerBranch[Util::nBranches];
    TH2F* mRpTrackRecoEffiSystematicsPerBranch[Util::nBranches];
    
    Double_t mEnergyLossParameters[Util::nSigns][Util::nDefinedParticles][3];
    Double_t mTpcRecoEffParameters[Util::nSigns][Util::nDefinedParticles][3];
    Double_t mTofMatchingEffParameters[Util::nSigns][Util::nDefinedParticles][3];
    
    Double_t mTofEff3D_PlateuShift[Util::nSigns][Util::nDefinedParticles];
    Double_t mTofEff3D_MultiplicationFactor[Util::nSigns][Util::nDefinedParticles];
    Double_t mTofEff3D_InflationShift[Util::nSigns][Util::nDefinedParticles];
    Double_t mTpcRecoEff3D_PlateuShift[Util::nSigns][Util::nDefinedParticles];
    Double_t mTpcRecoEff3D_MultiplicationFactor[Util::nSigns][Util::nDefinedParticles];
    Double_t mTpcRecoEff3D_InflationShift[Util::nSigns][Util::nDefinedParticles];
    
    Bool_t mTofEff3D_isUnchanged[Util::nSigns][Util::nDefinedParticles];
    Bool_t mTpcRecoEff3D_isUnchanged[Util::nSigns][Util::nDefinedParticles];
    
    Double_t mVertexFinderEff;
    Double_t mZVertexCutEff;
    Double_t mMissingPtCutEff; // estimate 0.95
    Double_t mTofTrigEff;
    Double_t mDeltaZ0CutEff;
    Double_t mTofClusterLimitEff;
    
    TH1F* mhRpCoincidence;
    TH1F* mhRpCoincidence_ETITVeto;
    TH1F* mhSimultaneousEastWestRpTrackLossProbability;
    TH2F* mhSimultaneousEastWestRpTrackLossProbability_zVertex;
    TH1F* mhEastWestRpTrackEffCorrelation;
    TH1F* mhEastWestVetoEffCorrelation;
    TH1F* mhEastWestTotalEffCorrelation;
    TH1F* mhEastWestTotalEffCorrelationVsMandelstamTSum[Util::nBranchesConfigurations];
    TH2F *mhEastWestTotalEffCorrelationVsMandelstamT[Util::nBranchesConfigurations];
    
    TH1F *mhPlusMinusTpcTrackRecoEfficiencyCorrelationVsAveragePt;
    TH1F *mhPlusMinusTpcTrackTotalEfficiencyCorrelationVsAveragePt;
    
    TH1F *mhRpTpcTotalEfficiencyCorrelationVsInvMass;
    TH2F* mhRpRecoEffDifference[Util::nBranches];
  
    TEfficiency* mGeometricalAcceptanceVsMass[nGeomAcceptances][2];
    
    TEfficiency *mAcceptanceForDeltaPhiCuts_Vs_t1t2[Util::nDeltaPhiRanges];
    TEfficiency *mAcceptanceForDeltaPhiCuts_Vs_t1t2_GenEx[Util::nDeltaPhiRanges];
    TEfficiency *mAcceptanceForDeltaPhiCuts_Vs_t1t2_DiMe[Util::nDeltaPhiRanges];
    TH2F *mAcceptanceForDeltaPhiCuts_Vs_t1t2_TH2F[Util::nDeltaPhiRanges];
    TH2F *mAcceptanceForDeltaPhiCuts_Vs_t1t2_GenEx_TH2F[Util::nDeltaPhiRanges];
    TH2F *mAcceptanceForDeltaPhiCuts_Vs_t1t2_DiMe_TH2F[Util::nDeltaPhiRanges];
    
    TH3F *mhFullVetoEfficiencyParameters[2];
    TEfficiency* mFullVetoEff[Util::nBranchesConfigurations];
    
    TH2F *mhDeadMatSystematics[Util::nBranches];
    TH2F *mhDeadMatCorrectionMultiplicativeFactor[Util::nBranches];
    
    TH1F *mhSquaredMassPurity[Util::nDefinedParticles];
    
    TEfficiency* mPidEfficiency[Util::nDefinedParticles][Util::nDefinedParticles];
    TEfficiency* mPidEfficiency_Pt[Util::nDefinedParticles][Util::nDefinedParticles];
    
    TEfficiency* mPtMissCutEff;
    TEfficiency* mPtMissCutEff_MC;

    
    Double_t mDEdxMeanOffsetParametersData[Util::nDefinedParticlesExtended][6];
    Double_t mDEdxMeanOffsetParametersMC[Util::nDefinedParticlesExtended][6];
    Double_t mDEdxWidthParametersData[Util::nDefinedParticlesExtended][6];
    Double_t mDEdxWidthParametersMC[Util::nDefinedParticlesExtended][6];
    
    TEfficiency *mVertexingEfficiencyVsDeltaZ0;
    TEfficiency *mVertexingEfficiencyVsDeltaZ0_MC;
    
    TEfficiency *mDeltaZ0CutEfficiencyVsLowerPt;
    TEfficiency *mDeltaZ0CutEfficiencyVsLowerPt_MC;
    
    TGraphAsymmErrors *mTagAndProbeEfficiency_1TofTrack[Util::nTriggers];
    TGraphAsymmErrors *mTagAndProbeEfficiency_2TofTracks[Util::nTriggers];
    TGraphAsymmErrors *mTagAndProbe_VertexingEffVsD0_1TofTrack[Util::nTriggers];
    TH2F *mTagAndProbe_AttachingToVertexingEffVsDcaRVsDcaZ_1TofTrack[Util::nTriggers];
    TEfficiency *mProbability_TofL0MultAtLeast1VsClusterSize;
    TEfficiency *mProbability_TofL0MultAtLeast2VsClusterSize;
    TEfficiency *mProbability_TofL0MultAtLeast1VsNRecoTofHits_1TofVtx_1TofTrack;
    TEfficiency *mProbability_TofL0MultAtLeast1VsNRecoTofHits_1TofVtx_2TofTracks;
    TEfficiency *mProbability_TofL0MultAtLeast2VsNRecoTofHits_1TofVtx_1TofTrack;
    TEfficiency *mProbability_TofL0MultAtLeast2VsNRecoTofHits_1TofVtx_2TofTracks;
    TEfficiency *mProbability_TofL0MultAtLeast1VsTrackEta_1TofVtx_1TofTrack;
    TEfficiency *mProbability_TofL0MultAtLeast2VsTrackEta_1TofVtx_1TofTrack;
    
    TH1F *mhEmbeddingVetoCorrectionFactor;
    
    TF1* mVertexFinderEff_InvMass[Util::nDefinedParticles];
    Bool_t mInterpolate;
    
    const Parameters *mParams;

    Efficiencies( const Parameters * = nullptr, TString = TString(""), Int_t = 0 );
    ~Efficiencies();

    static Efficiencies* mInst;

};

inline Double_t Efficiencies::tpcRecoEff(UInt_t sign, UInt_t pid, Double_t pt) const{
  if(sign<Util::nSigns && pid<Util::nDefinedParticles && pt>0) return mTpcReconstructionEfficiency[sign][pid]->Eval( pt );
  else{ std::cerr << "ERROR in Efficiencies::tpcRecoEff(UInt_t, UInt_t, Double_t)" << std::endl; return 9e99; }
}

inline Double_t Efficiencies::tofEff(UInt_t sign, UInt_t pid, Double_t pt) const{
  if(sign<Util::nSigns && pid<Util::nDefinedParticles && pt>0) return mTofMatchingEfficiency[sign][pid]->Eval( pt );
  else{ std::cerr << "ERROR in Efficiencies::tofEff(UInt_t, UInt_t, Double_t)" << std::endl; return 9e99; }
}

inline Double_t Efficiencies::onlineAndOfflineVetoEff(UInt_t branchesConf, Double_t instLumi /* in (ub*s)^-1 !!! */) const{
  if(instLumi>0 && branchesConf<Util::nBranchesConfigurations) return mFullVetosEfficiency[branchesConf]->Eval( instLumi );
  else{ std::cerr << "ERROR in Efficiencies::onlineAndOfflineVetoEff(UInt_t, Double_t)" << std::endl; return 9e99; }
}

inline Double_t Efficiencies::vertexFinderEff_InvMass(UInt_t pid, Double_t mass) const{
  if(pid<Util::nDefinedParticles && /*mass>0 &&*/ mVertexFinderEff_InvMass[pid])
    return mVertexFinderEff_InvMass[pid]->Eval( mass<2 ? mass : 2 );
  else{ std::cerr << "ERROR in Efficiencies::vertexFinderEff_InvMass(UInt_t, Double_t)" << std::endl; return 9e99; }
}

inline const TF1* Efficiencies::tpcEffTF1(UInt_t sign, UInt_t pid) const{
  if(sign<Util::nSigns && pid<Util::nDefinedParticles) return mTpcReconstructionEfficiency[sign][pid];
  else{ std::cerr << "ERROR in Efficiencies::tpcEffTF1(UInt_t, UInt_t)" << std::endl; return nullptr; }
}

inline const TF1* Efficiencies::tofEffTF1(UInt_t sign, UInt_t pid) const{
  if(sign<Util::nSigns && pid<Util::nDefinedParticles) return mTofMatchingEfficiency[sign][pid];
  else{ std::cerr << "ERROR in Efficiencies::tofEffTF1(UInt_t, UInt_t)" << std::endl; return nullptr; }
}

Efficiencies* Efficiencies::mInst = nullptr;

#endif
