#define Efficiencies_cxx

#include "Efficiencies.hh"
#include "Util.hh"
#include "TVector3.h"
#include "TFile.h"
#include "TArrayD.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TH2Poly.h"

Efficiencies::Efficiencies( const Parameters *par, TString codePath, int opt ){
  
  this->setParametersObject( par );
  mInterpolate = true;
  
  #ifdef RUN_OVER_EMBEDDED_MC
  cout << "\nCAUTION: Code is compiled using RUN_OVER_EMBEDDED_MC preprocessor directive!\n" << endl;
  #endif

//---- below parameters from embedded MC (20.01.17)
  // setting parameters of TPC track reconstruction efficiency functions
  mTpcRecoEffParameters[Util::MINUS][Util::PION][0] = 0.7525;
  mTpcRecoEffParameters[Util::MINUS][Util::PION][1] = 0.1532;
  mTpcRecoEffParameters[Util::MINUS][Util::PION][2] = 5.589;
  mTpcRecoEffParameters[Util::MINUS][Util::KAON][0] = 0.6131;
  mTpcRecoEffParameters[Util::MINUS][Util::KAON][1] = 0.2611;
  mTpcRecoEffParameters[Util::MINUS][Util::KAON][2] = 2.934;
  mTpcRecoEffParameters[Util::MINUS][Util::PROTON][0] = 0.766158;
  mTpcRecoEffParameters[Util::MINUS][Util::PROTON][1] = 0.29129;
  mTpcRecoEffParameters[Util::MINUS][Util::PROTON][2] = 7.00598;
  mTpcRecoEffParameters[Util::PLUS][Util::PION][0] = 0.7551;
  mTpcRecoEffParameters[Util::PLUS][Util::PION][1] = 0.1533;
  mTpcRecoEffParameters[Util::PLUS][Util::PION][2] = 4.894;
  mTpcRecoEffParameters[Util::PLUS][Util::KAON][0] = 0.619;
  mTpcRecoEffParameters[Util::PLUS][Util::KAON][1] = 0.2598;
  mTpcRecoEffParameters[Util::PLUS][Util::KAON][2] = 2.823;
  mTpcRecoEffParameters[Util::PLUS][Util::PROTON][0] = 0.855837;
  mTpcRecoEffParameters[Util::PLUS][Util::PROTON][1] = 0.288003;
  mTpcRecoEffParameters[Util::PLUS][Util::PROTON][2] = 7.64921;
//----
  
  // new TOF parameters 21.01.17
  mTofMatchingEffParameters[Util::MINUS][Util::PION][0] = 0.6439;
  mTofMatchingEffParameters[Util::MINUS][Util::PION][1] = 0.1312;
  mTofMatchingEffParameters[Util::MINUS][Util::PION][2] = 3.808;
  mTofMatchingEffParameters[Util::MINUS][Util::KAON][0] = 0.6084;
  mTofMatchingEffParameters[Util::MINUS][Util::KAON][1] = 0.2294;
  mTofMatchingEffParameters[Util::MINUS][Util::KAON][2] = 3.498;
  mTofMatchingEffParameters[Util::MINUS][Util::PROTON][0] = 0.6464;
  mTofMatchingEffParameters[Util::MINUS][Util::PROTON][1] = 0.3153;
  mTofMatchingEffParameters[Util::MINUS][Util::PROTON][2] = 4.588;
  mTofMatchingEffParameters[Util::PLUS][Util::PION][0] = 0.6499;
  mTofMatchingEffParameters[Util::PLUS][Util::PION][1] = 0.1313;
  mTofMatchingEffParameters[Util::PLUS][Util::PION][2] = 3.783;
  mTofMatchingEffParameters[Util::PLUS][Util::KAON][0] = 0.6171;
  mTofMatchingEffParameters[Util::PLUS][Util::KAON][1] = 0.2278;
  mTofMatchingEffParameters[Util::PLUS][Util::KAON][2] = 4.093;
  mTofMatchingEffParameters[Util::PLUS][Util::PROTON][0] = 0.6655;
  mTofMatchingEffParameters[Util::PLUS][Util::PROTON][1] = 0.3524;
  mTofMatchingEffParameters[Util::PLUS][Util::PROTON][2] = 8.802;
  //---
  
  for(int i=0; i<Util::nSigns; ++i)
    for(int j=0; j<Util::nDefinedParticles; ++j){
      mTpcRecoEff3D_MultiplicationFactor[i][j] = 1.0;
      mTpcRecoEff3D_InflationShift[i][j] = 0.0;
      mTofEff3D_MultiplicationFactor[i][j] = 1.0;
      mTofEff3D_InflationShift[i][j] = 0.0;
      mTofEff3D_isUnchanged[i][j] = true;
      mTpcRecoEff3D_isUnchanged[i][j] = true;
    }

    
  for(int i=0; i<Util::nSigns; ++i)
    for(int j=0; j<Util::nDefinedParticles; ++j){
      TString idStr; idStr.Form("%d%d", i, j);
      mTpcReconstructionEfficiency[i][j] = new TF1("mTpcReconstructionEfficiency_"+idStr, "[0]*TMath::Exp(-TMath::Power([1]/x, [2]))", 0.001, 100);
      for(int par=0; par<3; ++par)
	mTpcReconstructionEfficiency[i][j]->FixParameter(par, mTpcRecoEffParameters[i][j][par]);
    }

  for(int i=0; i<Util::nSigns; ++i)
    for(int j=0; j<Util::nDefinedParticles; ++j){
      TString idStr; idStr.Form("%d%d", i, j);
      mTofMatchingEfficiency[i][j] = new TF1("mTofMatchingEfficiency_"+idStr, "[0]*TMath::Exp(-TMath::Power([1]/x, [2]))", 0.001, 100);
      for(int par=0; par<3; ++par)
	mTofMatchingEfficiency[i][j]->FixParameter(par, mTofMatchingEffParameters[i][j][par]);
    }
  
  mVertexFinderEff = 0.35;
  mZVertexCutEff = 0.9545; //2sigma
  mMissingPtCutEff = 0.968; //average
  mTofTrigEff = 0.987236;
  
  mDeltaZ0CutEff = 0.975;
  
  mTofClusterLimitEff = 0.957;
  
  TFile *file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) +"MC/PidAnalysisHistograms.root");
  if( file->IsOpen() ){
    for(int i=0; i<Util::nDefinedParticles; ++i){
      for(int j=0; j<Util::nDefinedParticles; ++j){
        mPidEfficiency[i][j] = dynamic_cast<TEfficiency*>( file->Get( "pidEfficiency_true_"+TString(i==0?"pion":(i==1?"kaon":"proton"))+"_reco_"+TString(j==0?"pion":(j==1?"kaon":"proton"))) );
        mPidEfficiency[i][j]->SetDirectory(0);
        
        mPidEfficiency_Pt[i][j] = dynamic_cast<TEfficiency*>( file->Get( "pidEfficiency_pT_true_"+TString(i==0?"pion":(i==1?"kaon":"proton"))+"_reco_"+TString(j==0?"pion":(j==1?"kaon":"proton"))) );
        mPidEfficiency_Pt[i][j]->SetDirectory(0);
      }
    }
    file->Close();
  } else
    cout << "ERROR in Efficiencies when loading the PID efficiency file!!!" << endl;
  
  
  
  
  
  file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) +"MC/ptMissCutEff.root");
  if( file->IsOpen() ){
    mPtMissCutEff = dynamic_cast<TEfficiency*>( file->Get( "PmaxPMin_PtMissAll_clone" ) );
    mPtMissCutEff->SetDirectory(0);
    file->Close();
  } else
    cout << "ERROR in Efficiencies when loading the PtMiss cut efficiency file!!!" << endl;
  
  
  
  file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) +"MC/ptMissCutEffMC.root");
  if( file->IsOpen() ){
    mPtMissCutEff_MC = dynamic_cast<TEfficiency*>( file->Get( "PmaxPMin_PtMissAll_clone" ) );
    mPtMissCutEff_MC->SetName("PmaxPMin_PtMissAll_clone_MC");
    mPtMissCutEff_MC->SetDirectory(0);
    file->Close();
  } else
    cout << "ERROR in Efficiencies when loading the PtMiss cut efficiency file (for MC)!!!" << endl;
  


  
  file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) +"vertexingEff/VertexingEfficiency.root");
  if( file->IsOpen() ){
    mVertexingEfficiencyVsDeltaZ0 = dynamic_cast<TEfficiency*>( file->Get("VertexingEfficiencyVsDeltaZ0") );
    mVertexingEfficiencyVsDeltaZ0->SetDirectory(0);
    file->Close();
  } else
    cout << "ERROR in Efficiencies when loading the vertexing efficiency file!!!" << endl;
  
  file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) +"vertexingEff/VertexingEfficiency_MC.root");
  if( file->IsOpen() ){
    mVertexingEfficiencyVsDeltaZ0_MC = dynamic_cast<TEfficiency*>( file->Get("DeltaZ0_All_CutSet6_clone;1") );
    mVertexingEfficiencyVsDeltaZ0_MC ->SetDirectory(0);
    file->Close();
  } else
    cout << "ERROR in Efficiencies when loading the vertexing efficiency file (MC)!!!" << endl;
  
  
  
  
  
  file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) +"vertexingEff/DeltaZ0Eff.root");
  if( file->IsOpen() ){
    mDeltaZ0CutEfficiencyVsLowerPt = dynamic_cast<TEfficiency*>( file->Get("deltaZ0Eff_bkgdSub") );
    mDeltaZ0CutEfficiencyVsLowerPt->SetDirectory(0);
    file->Close();
  } else
    cout << "ERROR in Efficiencies when loading the deltaZ0 efficiency file!!!" << endl;
  
  file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) +"vertexingEff/DeltaZ0EfficiencyMC.root");
  if( file->IsOpen() ){
    mDeltaZ0CutEfficiencyVsLowerPt_MC = dynamic_cast<TEfficiency*>( file->Get("effDeltaZ0_vs_lowerPt") );
    mDeltaZ0CutEfficiencyVsLowerPt_MC->SetDirectory(0);
    file->Close();
  } else
    cout << "ERROR in Efficiencies when loading the deltaZ0 efficiency file (MC)!!!" << endl;
  
  
  
  
  file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) +"vertexingEff/zVtxVsFill.root");
  if( file->IsOpen() ){
    mhZVtxMeanVsFill = dynamic_cast<TH1D*>( file->Get("fillMean") );
    mhZVtxMeanVsFill->SetDirectory(0);
    mhZVtxSigmaVsFill = dynamic_cast<TH1D*>( file->Get("fillWidth") );
    mhZVtxSigmaVsFill->SetDirectory(0);
    file->Close();
  } else
    cout << "ERROR in Efficiencies when loading the zVtx vs. fill number file!!!" << endl;
  
  
  // ------------------- for TAG and PROBE ----------------------
  
//   file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) +"vertexingEff/FastTofEfficiencySelector_PiPiEmbedding.CPT.root");
  file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) +"vertexingEff/FastTofEfficiencySelector_PiPiEmbedding.root");
  if( file->IsOpen() ){
    mTagAndProbeEfficiency_1TofTrack[Util::CPT] = new TGraphAsymmErrors( *dynamic_cast<TGraphAsymmErrors*>( file->Get("gvertexingEffVsDeltaZ0_1TofTrack") ) );
    mTagAndProbeEfficiency_2TofTracks[Util::CPT] = new TGraphAsymmErrors( *dynamic_cast<TGraphAsymmErrors*>( file->Get("gvertexingEffVsDeltaZ0_2TofTracks") ) );
    mTagAndProbe_VertexingEffVsD0_1TofTrack[Util::CPT] = new TGraphAsymmErrors( *dynamic_cast<TGraphAsymmErrors*>( file->Get("gvertexingEffVsD0_1TofTrack") ) );
    mTagAndProbe_AttachingToVertexingEffVsDcaRVsDcaZ_1TofTrack[Util::CPT] = dynamic_cast<TH2F*>( file->Get("attachingToVertexingEffVsDcaRVsDcaZ_1TofTrack") );
    mTagAndProbe_AttachingToVertexingEffVsDcaRVsDcaZ_1TofTrack[Util::CPT]->SetDirectory(0);
    file->Close();
  } else
//     cout << "ERROR in Efficiencies when loading the FastTofEfficiencySelector_PiPiEmbedding.CPT.root file!!!" << endl;
    cout << "ERROR in Efficiencies when loading the FastTofEfficiencySelector_PiPiEmbedding.root file!!!" << endl;
//   file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) +"vertexingEff/FastTofEfficiencySelector_PiPiEmbedding.CPT2.root");
  file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) +"vertexingEff/FastTofEfficiencySelector_PiPiEmbedding.root");
  if( file->IsOpen() ){
    mTagAndProbeEfficiency_1TofTrack[Util::CPT2] = new TGraphAsymmErrors( *dynamic_cast<TGraphAsymmErrors*>( file->Get("gvertexingEffVsDeltaZ0_1TofTrack") ) );
    mTagAndProbeEfficiency_2TofTracks[Util::CPT2] = new TGraphAsymmErrors( *dynamic_cast<TGraphAsymmErrors*>( file->Get("gvertexingEffVsDeltaZ0_2TofTracks") ) );
    mTagAndProbe_VertexingEffVsD0_1TofTrack[Util::CPT2] = new TGraphAsymmErrors( *dynamic_cast<TGraphAsymmErrors*>( file->Get("gvertexingEffVsD0_1TofTrack") ) );
    mTagAndProbe_AttachingToVertexingEffVsDcaRVsDcaZ_1TofTrack[Util::CPT2] = dynamic_cast<TH2F*>( file->Get("attachingToVertexingEffVsDcaRVsDcaZ_1TofTrack") );
    mTagAndProbe_AttachingToVertexingEffVsDcaRVsDcaZ_1TofTrack[Util::CPT2]->SetDirectory(0);
    file->Close();
  } else
//     cout << "ERROR in Efficiencies when loading the FastTofEfficiencySelector_PiPiEmbedding.CPT2.root file!!!" << endl;
    cout << "ERROR in Efficiencies when loading the FastTofEfficiencySelector_PiPiEmbedding.root file!!!" << endl;
  
  file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) +"zeroBiasCorrections/ProbabilitiesForTagAndProbe.root");
  if( file->IsOpen() ){
    mProbability_TofL0MultAtLeast1VsClusterSize = dynamic_cast<TEfficiency*>( file->Get("probability_TofL0MultAtLeast1VsClusterSize") );
    mProbability_TofL0MultAtLeast1VsClusterSize->SetDirectory(0);
    mProbability_TofL0MultAtLeast2VsClusterSize = dynamic_cast<TEfficiency*>( file->Get("probability_TofL0MultAtLeast2VsClusterSize") );
    mProbability_TofL0MultAtLeast2VsClusterSize->SetDirectory(0);
    
    mProbability_TofL0MultAtLeast1VsNRecoTofHits_1TofVtx_1TofTrack = dynamic_cast<TEfficiency*>( file->Get("probability_TofL0MultAtLeast1VsNRecoTofHits_1TofVtx_1TofTrack") );
    mProbability_TofL0MultAtLeast1VsNRecoTofHits_1TofVtx_1TofTrack->SetDirectory(0);
    mProbability_TofL0MultAtLeast1VsNRecoTofHits_1TofVtx_2TofTracks = dynamic_cast<TEfficiency*>( file->Get("probability_TofL0MultAtLeast1VsNRecoTofHits_1TofVtx_2TofTracks") );
    mProbability_TofL0MultAtLeast1VsNRecoTofHits_1TofVtx_2TofTracks->SetDirectory(0);
    mProbability_TofL0MultAtLeast2VsNRecoTofHits_1TofVtx_1TofTrack = dynamic_cast<TEfficiency*>( file->Get("probability_TofL0MultAtLeast2VsNRecoTofHits_1TofVtx_1TofTrack") );
    mProbability_TofL0MultAtLeast2VsNRecoTofHits_1TofVtx_1TofTrack->SetDirectory(0);
    mProbability_TofL0MultAtLeast2VsNRecoTofHits_1TofVtx_2TofTracks = dynamic_cast<TEfficiency*>( file->Get("probability_TofL0MultAtLeast2VsNRecoTofHits_1TofVtx_2TofTracks") );
    mProbability_TofL0MultAtLeast2VsNRecoTofHits_1TofVtx_2TofTracks->SetDirectory(0);
    
    mProbability_TofL0MultAtLeast1VsTrackEta_1TofVtx_1TofTrack = dynamic_cast<TEfficiency*>( file->Get("probability_TofL0MultAtLeast1VsTrackEta_1TofVtx_1TofTrack") );
    mProbability_TofL0MultAtLeast1VsTrackEta_1TofVtx_1TofTrack->SetDirectory(0);
    mProbability_TofL0MultAtLeast2VsTrackEta_1TofVtx_1TofTrack = dynamic_cast<TEfficiency*>( file->Get("probability_TofL0MultAtLeast2VsTrackEta_1TofVtx_1TofTrack") );
    mProbability_TofL0MultAtLeast2VsTrackEta_1TofVtx_1TofTrack->SetDirectory(0);
    
    file->Close();
  } else
    cout << "ERROR in Efficiencies when loading the zeroBiasCorrections/ProbabilitiesForTagAndProbe.root file!!!" << endl;
  
  // -------------------------------- -------------------------
  
  
  
  file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) +"tpcEff_tofEff/TpcTofEffCorrelation.root");
  if( file->IsOpen() ){
    mhPlusMinusTpcTrackRecoEfficiencyCorrelationVsAveragePt = dynamic_cast<TH1F*>( file->Get( "hPlusMinusTpcEffCorrelationVsAveragePt" ) );
    mhPlusMinusTpcTrackRecoEfficiencyCorrelationVsAveragePt->SetDirectory(0);
    mhPlusMinusTpcTrackRecoEfficiencyCorrelationVsAveragePt->SetBinContent( mhPlusMinusTpcTrackRecoEfficiencyCorrelationVsAveragePt->GetNbinsX()+1, mhPlusMinusTpcTrackRecoEfficiencyCorrelationVsAveragePt->GetBinContent( mhPlusMinusTpcTrackRecoEfficiencyCorrelationVsAveragePt->GetNbinsX() ) );
    mhPlusMinusTpcTrackTotalEfficiencyCorrelationVsAveragePt = dynamic_cast<TH1F*>( file->Get( "hPlusMinusTotalEffCorrelationVsAveragePt" ) );
    mhPlusMinusTpcTrackTotalEfficiencyCorrelationVsAveragePt->SetDirectory(0);
    mhPlusMinusTpcTrackTotalEfficiencyCorrelationVsAveragePt->SetBinContent( mhPlusMinusTpcTrackTotalEfficiencyCorrelationVsAveragePt->GetNbinsX()+1, mhPlusMinusTpcTrackTotalEfficiencyCorrelationVsAveragePt->GetBinContent( mhPlusMinusTpcTrackTotalEfficiencyCorrelationVsAveragePt->GetNbinsX() ) );
    file->Close();
  } else
    cout << "ERROR in Efficiencies when loading the tpcEff_tofEff/TpcTofEffCorrelation.root file!!!" << endl;
  
  
  #ifdef RUN_OVER_EMBEDDED_MC
  file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) +/*"rpAcceptance/FastRpTrackRecoEfficiencySelector.root.3"*/"rpAcceptance/FastRpTrackRecoEfficiencySelector.root.5");
  #else
  file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) +/*"rpAcceptance/FastRpTrackRecoEfficiencySelector.root.2"*/"rpAcceptance/FastRpTrackRecoEfficiencySelector.root.4");
  #endif
  if( file->IsOpen() ){
    for(int i=0; i<Util::nBranches; ++i){
      mRpTrackRecoEffiPerBranch[i] = dynamic_cast<TH3F*>( file->Get( Form("rpTrackRecoEffiPerBranch_AfterDivergence_%d",i) ) );
      mRpTrackRecoEffiPerBranch[i]->SetDirectory(0);
      mRpTrackRecoEffiIntegratedOverZVtxPerBranch[i] = dynamic_cast<TH2D*>( file->Get( Form("recoEff_AfterDivergence_IntegratedOverZVtx_%d",i) ) ); //ALERT this is not in the *root.backup file
      mRpTrackRecoEffiIntegratedOverZVtxPerBranch[i]->SetDirectory(0);
      mRpTrackRecoEffiPerBranch_FullRecoEfficiency[i] = dynamic_cast<TH3F*>( file->Get( Form("rpTrackRecoEffiPerBranch_AfterDivergence_FullRecoEfficiency%d",i) ) );
      mRpTrackRecoEffiPerBranch_FullRecoEfficiency[i]->SetDirectory(0);
      mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff[i] = dynamic_cast<TH3F*>( file->Get( Form("rpTrackRecoEffiPerBranch_AfterDivergence_DeadMaterialVetoEff%d",i) ) );
      mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff[i]->SetDirectory(0);
      for(int j=0; j<2; ++j){
        mRpTrackRecoEffiPerBranchPerTrackType[i][j] = dynamic_cast<TH3F*>( file->Get( Form("rpTrackRecoEffiPerBranchPerType_AfterDivergence_%d%d",i,j) ) );
        mRpTrackRecoEffiPerBranchPerTrackType[i][j]->SetDirectory(0);
      }
      mhSimultaneousEastWestRpTrackLossProbability = dynamic_cast<TH1F*>( file->Get( "simultaneousEastWestRpTrackLossProbability" ) );
      mhSimultaneousEastWestRpTrackLossProbability->SetDirectory(0);
      mhSimultaneousEastWestRpTrackLossProbability_zVertex = dynamic_cast<TH2F*>( file->Get( "simultaneousEastWestRpTrackLossProbability_zVertex" ) );
      mhSimultaneousEastWestRpTrackLossProbability_zVertex->SetDirectory(0);
      mhEastWestRpTrackEffCorrelation = dynamic_cast<TH1F*>( file->Get( "eastWestRpTrackEffCorrelation" ) );
      mhEastWestRpTrackEffCorrelation->SetDirectory(0);
      mhEastWestVetoEffCorrelation = dynamic_cast<TH1F*>( file->Get( "eastWestVetoEffCorrelation" ) );
      mhEastWestVetoEffCorrelation->SetDirectory(0);
      mhEastWestTotalEffCorrelation = dynamic_cast<TH1F*>( file->Get( "eastWestTotalEffCorrelation" ) );
      mhEastWestTotalEffCorrelation->SetDirectory(0);
    }
    for(int i=0; i<Util::nBranchesConfigurations; ++i){
      mhEastWestTotalEffCorrelationVsMandelstamTSum[i] = dynamic_cast<TH1F*>( file->Get( Form("hEastWestTotalEffCorrelationVsMandelstamTSum_%d", i) ) );
      mhEastWestTotalEffCorrelationVsMandelstamTSum[i]->SetDirectory(0);
      mhEastWestTotalEffCorrelationVsMandelstamTSum[i]->SetBinContent(1, mhEastWestTotalEffCorrelationVsMandelstamTSum[i]->GetBinContent(2) );
      mhEastWestTotalEffCorrelationVsMandelstamTSum[i]->SetBinContent(mhEastWestTotalEffCorrelationVsMandelstamTSum[i]->GetNbinsX(), mhEastWestTotalEffCorrelationVsMandelstamTSum[i]->GetBinContent(mhEastWestTotalEffCorrelationVsMandelstamTSum[i]->GetNbinsX()-1) );
    }
    for(int i=0; i<Util::nBranchesConfigurations; ++i){
      mhEastWestTotalEffCorrelationVsMandelstamT[i] = dynamic_cast<TH2F*>( file->Get( Form("hEastWestTotalEffCorrelationVsMandelstamT_%d", i) ) );
      if( mhEastWestTotalEffCorrelationVsMandelstamT[i] ){
        mhEastWestTotalEffCorrelationVsMandelstamT[i]->SetDirectory(0);
      }
    }
    
    file->Close();
  } else
    #ifdef RUN_OVER_EMBEDDED_MC
    cout << "ERROR in Efficiencies when loading the rpAcceptance/FastRpTrackRecoEfficiencySelector.root.3 file!!!" << endl;
    #else
    cout << "ERROR in Efficiencies when loading the rpAcceptance/FastRpTrackRecoEfficiencySelector.root.2 file!!!" << endl;
    #endif
  

  file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) +"MC/PidPurity.root");
  if( file->IsOpen() ){
    Util *util = Util::instance();
    for(int i=0; i<Util::nDefinedParticles; ++i){
      mhSquaredMassPurity[i] = dynamic_cast<TH1F*>( file->Get( "SquaredMassPurity_"+util->particleName(i) ) );
      mhSquaredMassPurity[i]->SetDirectory(0);
    }
    file->Close();
  }
  
  
  
  file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) +"tpcEff_tofEff/TotalCorrelation.root");
  if( file->IsOpen() ){
    mhRpTpcTotalEfficiencyCorrelationVsInvMass = dynamic_cast<TH1F*>( file->Get( "hRpTpcTotalEfficiencyCorrelationVsInvMass" ) );
    mhRpTpcTotalEfficiencyCorrelationVsInvMass->SetDirectory(0);
    mhRpTpcTotalEfficiencyCorrelationVsInvMass->SetBinContent(1, mhRpTpcTotalEfficiencyCorrelationVsInvMass->GetBinContent(2) );
    mhRpTpcTotalEfficiencyCorrelationVsInvMass->SetBinContent(mhRpTpcTotalEfficiencyCorrelationVsInvMass->GetNbinsX()+1, mhRpTpcTotalEfficiencyCorrelationVsInvMass->GetBinContent(mhRpTpcTotalEfficiencyCorrelationVsInvMass->GetNbinsX()) );
    file->Close();
  }
  
  
  
  file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) +"rpAcceptance/FastRpTrackRecoEfficiencySelector.root.test");
  if( file->IsOpen() ){
    for(int i=0; i<Util::nBranches; ++i){
      mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff_ProtonTrackRequired[i] = dynamic_cast<TH3F*>( file->Get( Form("rpTrackRecoEffiPerBranch_AfterDivergence_GoodProtonRequired_DeadMaterialVetoEff%d",i) ) );
      mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff_ProtonTrackRequired[i]->SetDirectory(0);
    }
    file->Close();
  } else
    cout << "ERROR in Efficiencies when loading the rpAcceptance/FastRpTrackRecoEfficiencySelector.root.test file!!!" << endl;
  
  
  
  file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) +"rpAcceptance/RpEffSystematics.root");
  if( file->IsOpen() ){
    for(int i=0; i<Util::nBranches; ++i){
      mRpTrackRecoEffiSystematicsPerBranch[i] = dynamic_cast<TH2F*>( file->Get( Form("SystematicUncertaintyRpEff_%d", i) ) );
      mRpTrackRecoEffiSystematicsPerBranch[i]->SetDirectory(0);
      mRpTrackRecoEffiCorrectionPerBranch[i] = dynamic_cast<TH2F*>( file->Get( Form("RpEffCorrection_%d", i) ) );
      mRpTrackRecoEffiCorrectionPerBranch[i]->SetDirectory(0);
    }
    file->Close();
  } else
    cout << "ERROR in Efficiencies when loading the rpAcceptance/RpEffSystematics.root file!!!" << endl;
  
  
  file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) +"rpAcceptance/DeadMatSystematics.root");
  if( file->IsOpen() ){
    Util *util = Util::instance();
    for(int i=0; i<Util::nBranches; ++i){
      mhDeadMatSystematics[i] = dynamic_cast<TH2F*>( file->Get( "DeadMatSystematics_"+util->branchName(i) ) );
      mhDeadMatSystematics[i]->SetDirectory(0);
      mhDeadMatCorrectionMultiplicativeFactor[i] = dynamic_cast<TH2F*>( file->Get( "DeadMatCorrectionMultiplicativeFactor_"+util->branchName(i) ) );
      mhDeadMatCorrectionMultiplicativeFactor[i]->SetDirectory(0);
    }
    file->Close();
  } else
    cout << "ERROR in Efficiencies when loading the rpAcceptance/DeadMatSystematics.root file!!!" << endl;

  
  file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) +"zeroBiasCorrections/ZeroBiasCorrectionsSelector.root");
  if( file->IsOpen() ){
    mhRpCoincidence = dynamic_cast<TH1F*>( file->Get( "RpCoincidence" ) );
    mhRpCoincidence->SetDirectory(0);
    mhRpCoincidence_ETITVeto = dynamic_cast<TH1F*>( file->Get( "RpCoincidence_ETITVeto" ) );
    mhRpCoincidence_ETITVeto->SetDirectory(0);
    file->Close();
  } else
    cout << "ERROR in Efficiencies when loading the ZeroBiasCorrectionsSelectorfile!!!" << endl;
  
  
  
  file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) +"migrationCorrections/migrationsCorrection.root");
  if( file->IsOpen() ){
    Util *util = Util::instance();
    mMigrationsAndFakesCorrection_PtEta = dynamic_cast<TH2Poly*>( file->Get( "PtVsEta_MigrationsAndFakesTotalCorrection" ) );
    mMigrationsAndFakesCorrection_PtEta->SetDirectory(0);
    mMigrationsAndFakesCorrection_PxPy = dynamic_cast<TH2F*>( file->Get( "PyVsPx_MigrationsAndFakesTotalCorrection" ) );
    mMigrationsAndFakesCorrection_PxPy->SetDirectory(0);
    mMigrationsAndFakesCorrection_PxPy_LimitedMandelstamT = dynamic_cast<TH2F*>( file->Get( "PyVsPx_MigrationsAndFakesTotalCorrection_LimitedMandelstamT" ) );
    mMigrationsAndFakesCorrection_PxPy_LimitedMandelstamT->SetDirectory(0);
    mMigrationsAndFakesCorrection_MandelstamT_LimitedMandelstamT = dynamic_cast<TH1F*>( file->Get( "MandelstamT_MigrationsAndFakesTotalCorrection_LimitedMandelstamT" ) );
    mMigrationsAndFakesCorrection_MandelstamT_LimitedMandelstamT->SetDirectory(0);
    for(int i=0; i<Util::nSides; ++i){
      mMigrationsWithinFiducialRegionCorrection_PxPy[i] = dynamic_cast<TH2F*>( file->Get( "PyVsPx_MigrationsWithinFiducialRegion_"+util->sideName(i) ) );
      mMigrationsWithinFiducialRegionCorrection_PxPy[i]->SetDirectory(0);
    }
    file->Close();
  } else
    cout << "ERROR in Efficiencies when loading the MigrationsCorrection.root file!!!" << endl;
  
  #ifdef RUN_OVER_EMBEDDED_MC
  file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) +"migrationCorrections/FastRpTrackRecoEfficiencySelector_MigrationsCorrection.root");
  #else
  file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) +"migrationCorrections/FastRpTrackRecoEfficiencySelector_MigrationsCorrection_PP2PPEmbedding.root");
  #endif
  if( file->IsOpen() ){
    for(int i=0; i<Util::nSides; ++i){
      mMigrationsAndFakesCorrection_PxPy_PerSide[i] = dynamic_cast<TH2F*>( file->Get( Form("PyVsPx_MigrationsAndFakesTotalCorrection_PerSide_%d", i) ) );
      mMigrationsAndFakesCorrection_PxPy_PerSide[i]->SetDirectory(0);
    }
    for(int i=0; i<Util::nBranches; ++i){
      mMigrationsAndFakesCorrection_MandelstamT_LimitedMandelstamT_PerBranch[i] = dynamic_cast<TH1F*>( file->Get( Form("MandelstamT_MigrationsAndFakesTotalCorrection_LimitedMandelstamT_PerBranch_%d", i) ) );
      mMigrationsAndFakesCorrection_MandelstamT_LimitedMandelstamT_PerBranch[i]->SetDirectory(0);
    }
    file->Close();
  } else
    #ifdef RUN_OVER_EMBEDDED_MC
    cout << "ERROR in Efficiencies when loading the FastRpTrackRecoEfficiencySelector_MigrationsCorrection.root file!!!" << endl;
    #else
    cout << "ERROR in Efficiencies when loading the FastRpTrackRecoEfficiencySelector_MigrationsCorrection_PP2PPEmbedding.root file!!!" << endl;
    #endif
  
  
  file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) +"zeroBiasCorrections/ZeroBiasVetoCorrections.root");
  if( file->IsOpen() ){
    Util *util = Util::instance();
    for(int i=0; i<2; ++i){
      mhFullVetoEfficiencyParameters[i] = dynamic_cast<TH3F*>( file->Get( Form("FullVetoEfficiencyParameters_%d", i) ) );
      mhFullVetoEfficiencyParameters[i]->SetDirectory(0);
    }
    for(int bconf=0; bconf<Util::nBranchesConfigurations; ++bconf){
      mFullVetoEff[bconf] = dynamic_cast<TEfficiency*>( file->Get( Form("fullVetoEff_MaxNTofClusters%d_BbcLargeAdcThr%d_", mParams->maxNExtraTofClusters(), static_cast<int>(mParams->bbcLargeNoiseRateThreshold()*1e4)) + util->branchesConfigurationName(bconf) ) );
      mFullVetoEff[bconf]->SetDirectory(0);
    }
    file->Close();
  } else
    cout << "ERROR in Efficiencies when loading the ZeroBiasVetoCorrections file!!!" << endl;
  
  
  for(int bconf=0; bconf<Util::nBranchesConfigurations; ++bconf){
    mFullVetosEfficiency[bconf] = new TF1(Form("mFullVetosEfficiency_%d",bconf), "[0]*TMath::Exp(-[1]*x)", 0, 100);
    for(int i=0; i<2; ++i){
      int binNumber = mhFullVetoEfficiencyParameters[i]->FindBin( mParams->maxNExtraTofClusters(), mParams->bbcLargeNoiseRateThreshold()*1e4, bconf );
      double val = 0;
      if( binNumber > 0 && binNumber <= (mhFullVetoEfficiencyParameters[i]->GetNbinsX()+2)*(mhFullVetoEfficiencyParameters[i]->GetNbinsY()+2)*(mhFullVetoEfficiencyParameters[i]->GetNbinsZ()+2) )
        val = mhFullVetoEfficiencyParameters[i]->GetBinContent( binNumber );
      mFullVetosEfficiency[bconf]->FixParameter( i, val );
//       cout << "bconf = " << bconf << " NTofClMax = " << mParams->maxNExtraTofClusters() << " bbcLargeNoiseRateThreshold = " << mParams->bbcLargeNoiseRateThreshold()*1e4 << "        " << val << endl;
    }
  }
  
  
  file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) +"zeroBiasCorrections/ZeroBiasCorrRatio.root");
  if( file->IsOpen() ){
    mhEmbeddingVetoCorrectionFactor = dynamic_cast<TH1F*>( file->Get( "DataToMC" ) );
    mhEmbeddingVetoCorrectionFactor->SetDirectory(0);
    file->Close();
  } else
    cout << "ERROR in Efficiencies when loading the ZeroBiasCorrRatio.root file!!!" << endl;
  
  
  
  
  
  mVertexFinderEff_InvMass[Util::PION] = new TF1("mVertexFinderEff_InvMass_pion", "pol4", 0, 100);
  mVertexFinderEff_InvMass[Util::PION]->FixParameter(0, 0.578089);
  mVertexFinderEff_InvMass[Util::PION]->FixParameter(1, -1.80336);
  mVertexFinderEff_InvMass[Util::PION]->FixParameter(2, 2.83227);
  mVertexFinderEff_InvMass[Util::PION]->FixParameter(3, -1.43645);
  mVertexFinderEff_InvMass[Util::PION]->FixParameter(4, 0.240502);
  mVertexFinderEff_InvMass[Util::KAON] = nullptr;
  mVertexFinderEff_InvMass[Util::PROTON] = nullptr;
  
  
  mTpcEffMultiplicativeCorrectionFactor = new TF1("mTpcEffMultiplicativeCorrectionFactor", "[0]*TMath::ATan([1]*x)+[2]", -1, 1); // x = track eta
  mTpcEffMultiplicativeCorrectionFactor->SetParameters( 4.34579e-02, 3.108, 9.45834e-01 );  
//   mTpcEffMultiplicativeCorrectionFactor = new TF1("mTpcEffMultiplicativeCorrectionFactor2", "[0]*TMath::Erf(x/[1]+[2]) + [3]", -1, 1); // x = track eta
//   mTpcEffMultiplicativeCorrectionFactor->SetParameters( 4.533e-02, 2.16546e-1, 5.1617e-1, 0.9411 );
  
 
  mDEdxMeanOffsetParametersMC[Util::PIOn][0] = 7.183e-8;
  mDEdxMeanOffsetParametersMC[Util::PIOn][1] = -1.647e-4;
  mDEdxMeanOffsetParametersMC[Util::PIOn][2] = 41.678;
  mDEdxMeanOffsetParametersMC[Util::PIOn][3] = 0;
  mDEdxMeanOffsetParametersMC[Util::PIOn][4] = 0;
  mDEdxMeanOffsetParametersMC[Util::PIOn][5] = 0;
  mDEdxMeanOffsetParametersMC[Util::KAOn][0] = 4.35888e-8;
  mDEdxMeanOffsetParametersMC[Util::KAOn][1] = -9.285e-6;
  mDEdxMeanOffsetParametersMC[Util::KAOn][2] = 7.6973;
  mDEdxMeanOffsetParametersMC[Util::KAOn][3] = 0;
  mDEdxMeanOffsetParametersMC[Util::KAOn][4] = 0;
  mDEdxMeanOffsetParametersMC[Util::KAOn][5] = 0;
  mDEdxMeanOffsetParametersMC[Util::PROTOn][0] = 3.556e-8;
  mDEdxMeanOffsetParametersMC[Util::PROTOn][1] = -8.621e-6;
  mDEdxMeanOffsetParametersMC[Util::PROTOn][2] = 3.97979;
  mDEdxMeanOffsetParametersMC[Util::PROTOn][3] = 0;
  mDEdxMeanOffsetParametersMC[Util::PROTOn][4] = 0;
  mDEdxMeanOffsetParametersMC[Util::PROTOn][5] = 0;
  mDEdxMeanOffsetParametersMC[Util::ELECTROn][0] = -6.21854e-8;
  mDEdxMeanOffsetParametersMC[Util::ELECTROn][1] = 2.0648e-7;
  mDEdxMeanOffsetParametersMC[Util::ELECTROn][2] = 3.2414;
  mDEdxMeanOffsetParametersMC[Util::ELECTROn][3] = 0;
  mDEdxMeanOffsetParametersMC[Util::ELECTROn][4] = 0;
  mDEdxMeanOffsetParametersMC[Util::ELECTROn][5] = 0;
  mDEdxMeanOffsetParametersMC[Util::DEUTEROn][0] = -1.30462e-6;
  mDEdxMeanOffsetParametersMC[Util::DEUTEROn][1] = -5.26838e-6;
  mDEdxMeanOffsetParametersMC[Util::DEUTEROn][2] = 3.48682;
  mDEdxMeanOffsetParametersMC[Util::DEUTEROn][3] = 0;
  mDEdxMeanOffsetParametersMC[Util::DEUTEROn][4] = 0;
  mDEdxMeanOffsetParametersMC[Util::DEUTEROn][5] = 0; 

  mDEdxWidthParametersMC[Util::PIOn][0] = 0.070463;
  mDEdxWidthParametersMC[Util::PIOn][1] = 0;
  mDEdxWidthParametersMC[Util::PIOn][2] = 0;
  mDEdxWidthParametersMC[Util::PIOn][3] = -0.0014197;
  mDEdxWidthParametersMC[Util::PIOn][4] = 9.8603;
  mDEdxWidthParametersMC[Util::PIOn][5] = 0.95092;
  mDEdxWidthParametersMC[Util::KAOn][0] = 0.05108;
  mDEdxWidthParametersMC[Util::KAOn][1] = 0.033894;
  mDEdxWidthParametersMC[Util::KAOn][2] = 1.6751;
  mDEdxWidthParametersMC[Util::KAOn][3] = 0.01014;
  mDEdxWidthParametersMC[Util::KAOn][4] = 4.9341;
  mDEdxWidthParametersMC[Util::KAOn][5] = 0.5284;
  mDEdxWidthParametersMC[Util::PROTOn][0] = 0.06304;
  mDEdxWidthParametersMC[Util::PROTOn][1] = -7.72567;
  mDEdxWidthParametersMC[Util::PROTOn][2] = 27.1704;
  mDEdxWidthParametersMC[Util::PROTOn][3] = 0.00337034;
  mDEdxWidthParametersMC[Util::PROTOn][4] = 5.24451;
  mDEdxWidthParametersMC[Util::PROTOn][5] = 0.670272;
  mDEdxWidthParametersMC[Util::ELECTROn][0] = 0.03542;
  mDEdxWidthParametersMC[Util::ELECTROn][1] = 0.982;
  mDEdxWidthParametersMC[Util::ELECTROn][2] = 26.5747;
  mDEdxWidthParametersMC[Util::ELECTROn][3] = 0.0179093;
  mDEdxWidthParametersMC[Util::ELECTROn][4] = 41.5149;
  mDEdxWidthParametersMC[Util::ELECTROn][5] = 0.0950489;
  mDEdxWidthParametersMC[Util::DEUTEROn][0] = 0.0967217;
  mDEdxWidthParametersMC[Util::DEUTEROn][1] = -1525.93;
  mDEdxWidthParametersMC[Util::DEUTEROn][2] = 18.7519;
  mDEdxWidthParametersMC[Util::DEUTEROn][3] = 0;
  mDEdxWidthParametersMC[Util::DEUTEROn][4] = 0;
  mDEdxWidthParametersMC[Util::DEUTEROn][5] = 0;
  
  
  mDEdxMeanOffsetParametersData[Util::PIOn][0] = -1.39946e-8;
  mDEdxMeanOffsetParametersData[Util::PIOn][1] = 2.01233e-7;
  mDEdxMeanOffsetParametersData[Util::PIOn][2] = 10.3873;
  mDEdxMeanOffsetParametersData[Util::PIOn][3] = 0;
  mDEdxMeanOffsetParametersData[Util::PIOn][4] = 0;
  mDEdxMeanOffsetParametersData[Util::PIOn][5] = 0;
  mDEdxMeanOffsetParametersData[Util::KAOn][0] = 2.32521e-9;
  mDEdxMeanOffsetParametersData[Util::KAOn][1] = -3.68831e-6;
  mDEdxMeanOffsetParametersData[Util::KAOn][2] = 8.71201;
  mDEdxMeanOffsetParametersData[Util::KAOn][3] = 0;
  mDEdxMeanOffsetParametersData[Util::KAOn][4] = 0;
  mDEdxMeanOffsetParametersData[Util::KAOn][5] = 0;
  mDEdxMeanOffsetParametersData[Util::PROTOn][0] = -1.45811e-7;
  mDEdxMeanOffsetParametersData[Util::PROTOn][1] = 0.665524;
  mDEdxMeanOffsetParametersData[Util::PROTOn][2] = 59.0551;
  mDEdxMeanOffsetParametersData[Util::PROTOn][3] = 1.17146e-7;
  mDEdxMeanOffsetParametersData[Util::PROTOn][4] = 4.66009;
  mDEdxMeanOffsetParametersData[Util::PROTOn][5] = 0.643503;
  mDEdxMeanOffsetParametersData[Util::ELECTROn][0] = 9.0049e-8;
  mDEdxMeanOffsetParametersData[Util::ELECTROn][1] = 2.4937e-7;
  mDEdxMeanOffsetParametersData[Util::ELECTROn][2] = 8.83355;
  mDEdxMeanOffsetParametersData[Util::ELECTROn][3] = 0;
  mDEdxMeanOffsetParametersData[Util::ELECTROn][4] = 0;
  mDEdxMeanOffsetParametersData[Util::ELECTROn][5] = 0;
  mDEdxMeanOffsetParametersData[Util::DEUTEROn][0] = -1.91064e-7;
  mDEdxMeanOffsetParametersData[Util::DEUTEROn][1] = 0.0056367;
  mDEdxMeanOffsetParametersData[Util::DEUTEROn][2] = 14.4817;
  mDEdxMeanOffsetParametersData[Util::DEUTEROn][3] = 0;
  mDEdxMeanOffsetParametersData[Util::DEUTEROn][4] = 0;
  mDEdxMeanOffsetParametersData[Util::DEUTEROn][5] = 0;
  
    
  mDEdxWidthParametersData[Util::PIOn][0] = 0.0734176;
  mDEdxWidthParametersData[Util::PIOn][1] = 1.90731;
  mDEdxWidthParametersData[Util::PIOn][2] = 31.8635;
  mDEdxWidthParametersData[Util::PIOn][3] = -0.000819566;
  mDEdxWidthParametersData[Util::PIOn][4] = 22.7876;
  mDEdxWidthParametersData[Util::PIOn][5] = -0.653612;
  mDEdxWidthParametersData[Util::KAOn][0] = 0.0808044;
  mDEdxWidthParametersData[Util::KAOn][1] = -0.0400346;
  mDEdxWidthParametersData[Util::KAOn][2] = 7.95099;
  mDEdxWidthParametersData[Util::KAOn][3] = 0.00561903;
  mDEdxWidthParametersData[Util::KAOn][4] = -17.0792;
  mDEdxWidthParametersData[Util::KAOn][5] = -0.269016;
  mDEdxWidthParametersData[Util::PROTOn][0] = 0.07954;
  mDEdxWidthParametersData[Util::PROTOn][1] = 0.181;
  mDEdxWidthParametersData[Util::PROTOn][2] = 12.1176;
  mDEdxWidthParametersData[Util::PROTOn][3] = 0;
  mDEdxWidthParametersData[Util::PROTOn][4] = 0;
  mDEdxWidthParametersData[Util::PROTOn][5] = 0;
  mDEdxWidthParametersData[Util::ELECTROn][0] = 0.0679527;
  mDEdxWidthParametersData[Util::ELECTROn][1] = -0.000881926;
  mDEdxWidthParametersData[Util::ELECTROn][2] = 1.54877;
  mDEdxWidthParametersData[Util::ELECTROn][3] = 0;
  mDEdxWidthParametersData[Util::ELECTROn][4] = 0;
  mDEdxWidthParametersData[Util::ELECTROn][5] = 0;
  mDEdxWidthParametersData[Util::DEUTEROn][0] = 0.116086;
  mDEdxWidthParametersData[Util::DEUTEROn][1] = -0.147449;
  mDEdxWidthParametersData[Util::DEUTEROn][2] = 2.8896;
  mDEdxWidthParametersData[Util::DEUTEROn][3] = 0;
  mDEdxWidthParametersData[Util::DEUTEROn][4] = 0;
  mDEdxWidthParametersData[Util::DEUTEROn][5] = 0;
 
  loadEfficiencies3D( codePath, opt );
  loadGeometricalAcceptance( codePath );
}



Double_t Efficiencies::tpcEffMultiplicativeCorrectionFactor(UInt_t rr, Double_t eta) const{
  return (rr == Util::RUN_RANGE_1 ) ? mTpcEffMultiplicativeCorrectionFactor->Eval(eta) : 1.0 ;
}


Double_t Efficiencies::geomAcc(UInt_t accType, UInt_t deltaPhiBin, Double_t m) const{
  if(accType < nGeomAcceptances && deltaPhiBin<Util::nDeltaPhiRanges){
    TEfficiency *acc = mGeometricalAcceptanceVsMass[accType][deltaPhiBin];
    if( m > 2.0 ) m = 1.9999; //ALERT
    return acc->GetEfficiency( acc->FindFixBin(m) );
  }
  else{ std::cerr << "ERROR in Efficiencies::geomAcc(UInt_t, UInt_t, Double_t)" << std::endl; return 0.00000000001; }
}


void Efficiencies::loadGeometricalAcceptance( TString codePath ){
  for(int i=0; i<nGeomAcceptances; ++i)
    for(int j=0; j<2; ++j)
      mGeometricalAcceptanceVsMass[i][j] = nullptr;
    
  TFile *file;
  file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) + "tpcEff_tofEff/massAccForFit.root", "READ" );
  for(int i=0; i<2; ++i){
    mGeometricalAcceptanceVsMass[S0m_PHASE_SPACE][i] = dynamic_cast<TEfficiency*>( file->Get( Form("mMassAcceptance_LimitedMandelstamT_DeltaPhiNarrowerBin%d_S0m", i+1) ) );
    mGeometricalAcceptanceVsMass[S0m_PHASE_SPACE][i]->SetDirectory(0);
    mGeometricalAcceptanceVsMass[D0m_PHASE_SPACE][i] = dynamic_cast<TEfficiency*>( file->Get( Form("mMassAcceptance_LimitedMandelstamT_DeltaPhiNarrowerBin%d_D0m", i+1) ) );
    mGeometricalAcceptanceVsMass[D0m_PHASE_SPACE][i]->SetDirectory(0);
  }
  file->Close();
  
  file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) + "tpcEff_tofEff/massAccForFit_GenEx.root", "READ" );
  for(int i=0; i<2; ++i){
    mGeometricalAcceptanceVsMass[S0m_GENEX][i] = dynamic_cast<TEfficiency*>( file->Get( Form("mMassAcceptance_LimitedMandelstamT_DeltaPhiNarrowerBin%d_S0m", i+1) ) );
    mGeometricalAcceptanceVsMass[S0m_GENEX][i]->SetDirectory(0);
    mGeometricalAcceptanceVsMass[D0m_GENEX][i] = dynamic_cast<TEfficiency*>( file->Get( Form("mMassAcceptance_LimitedMandelstamT_DeltaPhiNarrowerBin%d_D0m", i+1) ) );
    mGeometricalAcceptanceVsMass[D0m_GENEX][i]->SetDirectory(0);
    mGeometricalAcceptanceVsMass[GENEX][i] = dynamic_cast<TEfficiency*>( file->Get( Form("mMassAcceptance_LimitedMandelstamT_DeltaPhiNarrowerBin%d_GenEx", i+1) ) );
    mGeometricalAcceptanceVsMass[GENEX][i]->SetDirectory(0);
  }
  file->Close();
  
  file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) + "tpcEff_tofEff/massAccForFit_DiMe.root", "READ" );
  for(int i=0; i<2; ++i){
    mGeometricalAcceptanceVsMass[S0m_DIME][i] = dynamic_cast<TEfficiency*>( file->Get( Form("mMassAcceptance_LimitedMandelstamT_DeltaPhiNarrowerBin%d_S0m", i+1) ) );
    mGeometricalAcceptanceVsMass[S0m_DIME][i]->SetDirectory(0);
    mGeometricalAcceptanceVsMass[D0m_DIME][i] = dynamic_cast<TEfficiency*>( file->Get( Form("mMassAcceptance_LimitedMandelstamT_DeltaPhiNarrowerBin%d_D0m", i+1) ) );
    mGeometricalAcceptanceVsMass[D0m_DIME][i]->SetDirectory(0);
    mGeometricalAcceptanceVsMass[DIME][i] = dynamic_cast<TEfficiency*>( file->Get( Form("mMassAcceptance_LimitedMandelstamT_DeltaPhiNarrowerBin%d_DiMe", i+1) ) );
    mGeometricalAcceptanceVsMass[DIME][i]->SetDirectory(0);
  }
  file->Close();
  
  
  for(int i=0; i<2; ++i){
    const int nBins = 200;
    const double massMin = 0.0;
    const double massMax = 2.0;
    mGeometricalAcceptanceVsMass[S0mD0mMix_PHASE_SPACE][i] = new TEfficiency( Form("mMassAcceptance_LimitedMandelstamT_DeltaPhiNarrowerBin%d_S0mD0mMix", i+1),Form("mMassAcceptance_LimitedMandelstamT_DeltaPhiNarrowerBin%d_S0mD0mMix", i+1), nBins, massMin, massMax);
    const int nEvts=1e5;
    const double binWidth = (massMax - massMin)/nBins;
    for(int bin=1; bin<=nBins; ++bin){
      const double m = binWidth * (bin-0.5);
      double accS0m = mGeometricalAcceptanceVsMass[S0m_PHASE_SPACE][i]->GetEfficiency( mGeometricalAcceptanceVsMass[S0m_PHASE_SPACE][i]->FindFixBin(m) );
      double accD0m = mGeometricalAcceptanceVsMass[D0m_PHASE_SPACE][i]->GetEfficiency( mGeometricalAcceptanceVsMass[D0m_PHASE_SPACE][i]->FindFixBin(m) );
      double errorFuncVal = 0.5*(1.0 + TMath::Erf(15*(m-1.1)));
      double acceptance = errorFuncVal * accD0m + (1.0-errorFuncVal) * accS0m;
      for(int jj=0; jj<nEvts; ++jj)
        mGeometricalAcceptanceVsMass[S0mD0mMix_PHASE_SPACE][i]->Fill( jj<acceptance*nEvts, m );
    }
  }
  
}




Efficiencies::~Efficiencies(){
  for(int i=0; i<Util::nSigns; ++i)
    for(int j=0; j<Util::nDefinedParticles; ++j)
      if(mTpcReconstructionEfficiency[i][j]) delete mTpcReconstructionEfficiency[i][j];
  for(int i=0; i<Util::nSigns; ++i)
    for(int j=0; j<Util::nDefinedParticles; ++j)
      if(mTofMatchingEfficiency[i][j]) delete mTofMatchingEfficiency[i][j];
  for(int i=0; i<Util::nBranchesConfigurations; ++i)
    if(mFullVetosEfficiency[i]) delete mFullVetosEfficiency[i];
  for(int j=0; j<Util::nDefinedParticles; ++j)
    if(mVertexFinderEff_InvMass[j]) delete mVertexFinderEff_InvMass[j];
//   if(mTagAndProbeEfficiency_1TofTrack) delete mTagAndProbeEfficiency_1TofTrack;
//   if(mTagAndProbeEfficiency_2TofTracks) delete mTagAndProbeEfficiency_2TofTracks;
}

Efficiencies* Efficiencies::instance( const Parameters *par, TString codePath, int opt ){
  if(!mInst) mInst = new Efficiencies(par, codePath, opt);
  return mInst;
}

void Efficiencies::loadEfficiencies3D( TString codePath, int opt ){
  TString histNames[Util::nSigns][Util::nDefinedParticles];
  histNames[Util::MINUS][Util::PION] = TString("0");
  histNames[Util::MINUS][Util::KAON] = TString("1");
  histNames[Util::MINUS][Util::PROTON] = TString("2");
  histNames[Util::PLUS][Util::PION] = TString("3");
  histNames[Util::PLUS][Util::KAON] = TString("4");
  histNames[Util::PLUS][Util::PROTON] = TString("5");

  
//   Przygotowałem wydajność 2D TPC-TOF (eta, pT) w binach Vz, o których Pan Doktor pisał. Dodatkowo zrobiłem to dla różnych cięć d0 oraz dla dwóch wariantów:
// 1. Bez cięć na pT^reco, eta^reco
// 2. Z cięciami na pT^reco, eta^reco
// 
// 
// Nazewnictwo plików pdf to tofEffi_[proces]_d0_[ciecie_d0]_etapt_[cięcie na reco/lub nie].pdf
// 
// W plikach root są histogramy 3D (Vz, pT, eta):
// hTPCEffi[nr czastki 0-6][ciecie 1/2][d0 0-4]
// Cząstki pi-,K-,p_bar, pi+,K+,p, naladowane
  
  #ifdef TIGHT_D0_CUT
  TFile *file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) + "tpcEff_tofEff/tof_d0Down.root" );
  #elif defined LOOSE_D0_CUT
  TFile *file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) + "tpcEff_tofEff/tof_d0Up.root" );
  #elif defined LOOSE_NHITS_CUT
  TFile *file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) + "tpcEff_tofEff/tof_nHitsDown.root" );
  #elif defined TIGHT_NHITS_CUT
  TFile *file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) + "tpcEff_tofEff/tof_nHitsUp.root" );
  #else
  TFile *file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) + (opt==0 ? /*"tpcEff_tofEff/interpolated_FilledEmptyBins_BBC.root"*/"tpcEff_tofEff/etaPhiEfficiency_4_1_19_delta015.root" : "tpcEff_tofEff/effWithBinningForSystematics.root" ) );
  #endif
  
  if( file->IsOpen() ){
    cout << "TOF efficiency file name = " << file->GetName() << endl;
    
    for(int i=0; i<Util::nSigns; ++i){
      for(int j=0; j<Util::nDefinedParticles; ++j){
        mhTofEff3D[i][j] = dynamic_cast<TH3F* >( file->Get("hTOFEffiCD"+histNames[i][j]+"12") );
        mhTofEff3D[i][j]->SetDirectory(0);
        mhTofEff3D_withRecoCuts[i][j] = nullptr;
        // 	mhTofEff3D_withRecoCuts[i][j] = dynamic_cast<TH3F* >( file->Get("hTOFEffiCD"+histNames[i][j]+"22") );
        // 	mhTofEff3D_withRecoCuts[i][j]->SetDirectory(0);
      }
    }
    file->Close();
  }
  
  
  #ifdef TIGHT_D0_CUT
  file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) + "tpcEff_tofEff/tpc_d0Down.root" );
  #elif defined LOOSE_D0_CUT
  file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) + "tpcEff_tofEff/tpc_d0Up.root" );
  #elif defined LOOSE_NHITS_CUT
  file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) + "tpcEff_tofEff/tpc_nHitsDown.root" );
  #elif defined TIGHT_NHITS_CUT
  file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) + "tpcEff_tofEff/tpc_nHitsUp.root" );
  #else
  file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) + "tpcEff_tofEff/etaPhiEfficiency_16_01_19_delta015_twoRuns.root" );
  #endif
  
  if( file->IsOpen() ){
    cout << "TPC efficiency file name = " << file->GetName() << endl;
    
    for(int i=0; i<Util::nSigns; ++i){
      for(int j=0; j<Util::nDefinedParticles; ++j){
        for(int k=0; k<Util::nRunRanges; ++k){
          mhTpcRecoEff3D[i][j][k] = dynamic_cast<TH3F* >( file->Get("hTPCEffiCD"+histNames[i][j]+Form("12%d", k)) );
          mhTpcRecoEff3D[i][j][k]->SetDirectory(0);
        }
        mhTpcRecoEff3D_withRecoCuts[i][j] = nullptr;
        // 	mhTpcRecoEff3D_withRecoCuts[i][j] = dynamic_cast<TH3F* >( file->Get("hTPCEffiCD"+histNames[i][j]+"22") );
        // 	mhTpcRecoEff3D_withRecoCuts[i][j]->SetDirectory(0);
      }
    }
    file->Close();
  }
  
  
  
  #ifdef TEST_GENEX_PION_EFFICIENCIES
//   file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) + "tpcEff_tofEff/genex_tof.root" );
  file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) + "tpcEff_tofEff/tof_genex_1.root" );
  if( file->IsOpen() ){
    cout << "TOF efficiency file name = " << file->GetName() << endl;
    
    for(int i=0; i<Util::nSigns; ++i){
        const int j=Util::PION;
        mhTofEff3D[i][j] = dynamic_cast<TH3F* >( file->Get("hTOFEffiCD"+histNames[i][j]+"12") );
        mhTofEff3D[i][j]->SetDirectory(0);
        mhTofEff3D_withRecoCuts[i][j] = nullptr;      
    }
    file->Close();
  }
  
//   file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) + "tpcEff_tofEff/genex_tpc.root" );
  file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) + "tpcEff_tofEff/tpc_genex_2.root" );
  if( file->IsOpen() ){
    cout << "TPC efficiency file name = " << file->GetName() << endl;
    
    for(int i=0; i<Util::nSigns; ++i){
        const int j=Util::PION;
        for(int k=0; k<Util::nRunRanges; ++k){
          mhTpcRecoEff3D[i][j][k] = dynamic_cast<TH3F* >( file->Get("hTPCEffiCD"+histNames[i][j]+Form("12%d", k)) );
          mhTpcRecoEff3D[i][j][k]->SetDirectory(0);
        }
        mhTpcRecoEff3D_withRecoCuts[i][j] = nullptr;
    }
    file->Close();
  }
  #endif
  
  
  
  
    file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) + Form("tpcEff_tofEff/tpc_phi.root") );  
    if( file->IsOpen() ){
      cout << "TPC efficiency file name #2 (phi-dependent) = " << file->GetName() << endl;
      
      for(int i=0; i<Util::nSigns; ++i){
        for(int j=0; j<Util::nDefinedParticles; ++j){
          for(int k=0; k<Util::nRunRanges; ++k){
            mhTpcRecoEff3D_Phi[i][j][k] = dynamic_cast<TH3F* >( file->Get("hTPCEffiCD"+histNames[i][j]+Form("12%d", k)) );
            mhTpcRecoEff3D_Phi[i][j][k]->SetDirectory(0);
          }
        }
      }
      file->Close();
    }
  
  
    file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) + Form("tpcEff_tofEff/tof_phi.root") );  
    if( file->IsOpen() ){
      cout << "TOF efficiency file name #2 (phi-dependent) = " << file->GetName() << endl;
      
      for(int i=0; i<Util::nSigns; ++i){
        for(int j=0; j<Util::nDefinedParticles; ++j){
          mhTofEff3D_Phi[i][j] = dynamic_cast<TH3F* >( file->Get("hTOFEffiCD"+histNames[i][j]+"12") );
          mhTofEff3D_Phi[i][j]->SetDirectory(0);
        }
      }
      file->Close();
    }
  
  
  
  file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) +"tpcEff_tofEff/tofEffCorrectionParameters.root");
  if( file->IsOpen() ){
      
      TH2F *tofEffCorrectionParameters = dynamic_cast<TH2F* >( file->Get("tofEffCorrectionParameters") );
      tofEffCorrectionParameters->SetDirectory(0);
      
      auto effDiff = [](double* x, double* p) {
          return (TMath::Erf((x[0] - p[0]) / p[1]) + 1) / 2. * p[2] - (TMath::Erf((x[0] - p[0+3]) / p[1+3]) + 1) / 2. * p[2+3];
      };
      
      TF1* effDataMinusMc[nEtaBins];
      for(int i=0; i<nEtaBins; ++i){
        effDataMinusMc[i] = new TF1(Form("effDataMinusMc_%d", i), effDiff, 0, 10.0, 6);
        for(int par=0; par<6; ++par)
            effDataMinusMc[i]->SetParameter( par, tofEffCorrectionParameters->GetBinContent(i+1, par+1) );
      }
      
      const int nBinsEta = 80;
      const int nBinsPt = 200;
      
      for(int j=0; j<Util::nDefinedParticles; ++j){
          mhTofEffCorrection2D_ptVsEta[j] = new TH2F( "tofEffCorrection2D_ptVsEta_"+TString(j==0?"pion":(j==1?"kaon":"proton")), "tofEffCorrection2D_ptVsEta_"+TString(j==0?"pion":(j==1?"kaon":"proton")), nBinsEta, -0.8, 0.8, nBinsPt, 0, 2);
          mhTofEffCorrection2D_ptVsEta[j]->SetDirectory(0);
            for(int binEta=1; binEta <= nBinsEta; ++binEta){
                const double etaVal = mhTofEffCorrection2D_ptVsEta[j]->GetXaxis()->GetBinCenter(binEta);
                const int etaBinCorrection = etaVal < -0.3 ? MINUS_07_03_ETA : ( etaVal < 0.0 ? MINUS_03_0_ETA : ( etaVal < 0.3 ? PLUS_03_0_ETA : PLUS_07_03_ETA ) );
                for(int binPt=1; binPt <= nBinsPt; ++binPt){
                    const double ptVal = mhTofEffCorrection2D_ptVsEta[j]->GetYaxis()->GetBinCenter(binPt);
                    const double correction = effDataMinusMc[etaBinCorrection]->Eval( ptVal );
                    mhTofEffCorrection2D_ptVsEta[j]->SetBinContent(binEta, binPt, correction);
                }
            }
      }
      file->Close();
      
  } else
      cout << "ERROR in Efficiencies when loading the 3D tof eff. corrections file!!!" << endl;
  
  
  
  file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) +"tpcEff_tofEff/TofEffDiffHft.root");
  if( file->IsOpen() ){
    for(int j=0; j<Util::nDefinedParticles; ++j){
      mhTofEffDifferenceHftMC3D_ptVsEtaVsZVtx[j] = new TH3F( *dynamic_cast<TH3F*>( file->Get( TString("effDifferece_HftMinusMC_pion") ) ) ); //ALERT!!!
      mhTofEffDifferenceHftMC3D_ptVsEtaVsZVtx[j]->SetDirectory(0);
    }
    file->Close();
  } else
    cout << "ERROR in Efficiencies when loading the 3D tof eff. corrections file!!!" << endl;
  
  
  
  
  file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) +"tpcEff_tofEff/TofCorrSyst.root");
  if( file->IsOpen() ){
    for(int j=0; j<Util::nDefinedParticles; ++j){
      mhTofEffExtraCorr2D_ptVsEta[j] = new TH2D( *dynamic_cast<TH2D*>( file->Get( TString(j==0?"pion":(j==1?"kaon":"proton"))+TString("_tofEffExtraCorr") ) ) );
      mhTofEffExtraCorr2D_ptVsEta[j]->SetDirectory(0);
      mhTofEffSystErr2D_ptVsEta[j] = new TH2D( *dynamic_cast<TH2D*>( file->Get( TString(j==0?"pion":(j==1?"kaon":"proton"))+TString("_tofEffSystErr") ) ) );
      mhTofEffSystErr2D_ptVsEta[j]->SetDirectory(0);
    }
    file->Close();
  } else
    cout << "ERROR in Efficiencies when loading the 3D tof eff. corrections file!!!" << endl;
  
  
  
  
  
  file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) +"tpcEff_tofEff/tpcSystPileUp.root");
  if( file->IsOpen() ){
    for(int i=0; i<Util::nRunRanges; ++i){
      for(int j=0; j<Util::nSigns; ++j){
        for(int k=0; k<2; ++k){
          mhSystematicsPileUpTPC[i][j][k] = dynamic_cast<TH1D*>( file->Get(Form("histSystPileUpPt%d%d%d",i,j,k)) );
          mhSystematicsPileUpTPC[i][j][k]->SetDirectory(0);
        }
      }
    }
    file->Close();
  } else
    cout << "ERROR in Efficiencies when tpc systematics file!!!" << endl;
  
  
  
  
  
  
  file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) +"tpcEff_tofEff/massAccForFit.root");
  if( file->IsOpen() ){
    for(int i=0; i<Util::nDeltaPhiRanges; ++i){
      mAcceptanceForDeltaPhiCuts_Vs_t1t2[i] = dynamic_cast<TEfficiency*>( file->Get( "mAcceptanceForDeltaPhiCuts_Vs_t1t2"+TString(Form("_deltaPhiBin_%d", i+1)) ) );
      mAcceptanceForDeltaPhiCuts_Vs_t1t2[i]->SetDirectory(0);
      mAcceptanceForDeltaPhiCuts_Vs_t1t2_TH2F[i] = dynamic_cast<TH2F*>( file->Get( "mAcceptanceForDeltaPhiCuts_Vs_t1t2"+TString(Form("_deltaPhiBin_%d_TH2F", i+1)) ) );
      mAcceptanceForDeltaPhiCuts_Vs_t1t2_TH2F[i]->SetDirectory(0);
    }
    file->Close();
  } else
    cout << "ERROR in Efficiencies when tpc systematics file!!!" << endl;
  
  
  
  file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) +"tpcEff_tofEff/massAccForFit_GenEx.root");
  if( file->IsOpen() ){
    for(int i=0; i<Util::nDeltaPhiRanges; ++i){
      mAcceptanceForDeltaPhiCuts_Vs_t1t2_GenEx[i] = dynamic_cast<TEfficiency*>( file->Get( "mAcceptanceForDeltaPhiCuts_Vs_t1t2"+TString(Form("_deltaPhiBin_%d", i+1)) ) );
      mAcceptanceForDeltaPhiCuts_Vs_t1t2_GenEx[i]->SetDirectory(0);
      mAcceptanceForDeltaPhiCuts_Vs_t1t2_GenEx[i]->SetName(Form("%s_GenEx", mAcceptanceForDeltaPhiCuts_Vs_t1t2_GenEx[i]->GetName()));
      mAcceptanceForDeltaPhiCuts_Vs_t1t2_GenEx_TH2F[i] = dynamic_cast<TH2F*>( file->Get( "mAcceptanceForDeltaPhiCuts_Vs_t1t2"+TString(Form("_deltaPhiBin_%d_TH2F", i+1)) ) );
      mAcceptanceForDeltaPhiCuts_Vs_t1t2_GenEx_TH2F[i]->SetDirectory(0);
      mAcceptanceForDeltaPhiCuts_Vs_t1t2_GenEx_TH2F[i]->SetName(Form("%s_GenEx", mAcceptanceForDeltaPhiCuts_Vs_t1t2_GenEx_TH2F[i]->GetName()));
    }
    file->Close();
  } else
    cout << "ERROR in Efficiencies when tpc systematics file!!!" << endl;
  
  
  
  file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) +"tpcEff_tofEff/massAccForFit_DiMe.root");
  if( file->IsOpen() ){
    for(int i=0; i<Util::nDeltaPhiRanges; ++i){
      mAcceptanceForDeltaPhiCuts_Vs_t1t2_DiMe[i] = dynamic_cast<TEfficiency*>( file->Get( "mAcceptanceForDeltaPhiCuts_Vs_t1t2"+TString(Form("_deltaPhiBin_%d", i+1)) ) );
      mAcceptanceForDeltaPhiCuts_Vs_t1t2_DiMe[i]->SetDirectory(0);
      mAcceptanceForDeltaPhiCuts_Vs_t1t2_DiMe[i]->SetName(Form("%s_DiMe", mAcceptanceForDeltaPhiCuts_Vs_t1t2_DiMe[i]->GetName()));
      mAcceptanceForDeltaPhiCuts_Vs_t1t2_DiMe_TH2F[i] = dynamic_cast<TH2F*>( file->Get( "mAcceptanceForDeltaPhiCuts_Vs_t1t2"+TString(Form("_deltaPhiBin_%d_TH2F", i+1)) ) );
      mAcceptanceForDeltaPhiCuts_Vs_t1t2_DiMe_TH2F[i]->SetDirectory(0);
      mAcceptanceForDeltaPhiCuts_Vs_t1t2_DiMe_TH2F[i]->SetName(Form("%s_DiMe", mAcceptanceForDeltaPhiCuts_Vs_t1t2_DiMe_TH2F[i]->GetName()));
    }
    file->Close();
  } else
    cout << "ERROR in Efficiencies when tpc systematics file!!!" << endl;
  
  
  
  
  
  
  file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) +"energyLossCorrections/energyLoss3D_OnePrtAlso.root");
  if( file->IsOpen() ){
    for(int k=0; k<2; ++k){
      TString id2(k==0 ? "SD":"CD");
      for(int i=0; i<Util::nSigns; ++i){
	for(int j=0; j<Util::nDefinedParticles; ++j){
	  for(int l=0; l<20; ++l){ // l numbers zVertex bin
	    TString id1 = Form("%d", l);
	    if( k==0 ){
	      mEnergyLossCorrection3D[i][j][l] = new TH3F( *dynamic_cast<TH3F*>( file->Get("hEnergyLoss3D"+id2+histNames[i][j]+id1) ) );
	      mEnergyLossCorrection3D[i][j][l]->SetDirectory(0);
	    } else{
	      mEnergyLossCorrection3D[i][j][l]->Add( dynamic_cast<TH3F*>( file->Get("hEnergyLoss3D"+id2+histNames[i][j]+id1) ) );
	    }
	  }
	}
      }
    }
    file->Close();
  } else
    cout << "ERROR in Efficiencies when loading the 3D-energy loss corrections file!!!" << endl;

  
  
  
  
  file = TFile::Open( ( (codePath=="") ? codePath : (codePath + "/")) +"deadMaterial/secondariesDeadMaterial.root");
  if( file->IsOpen() ){
    for(int i=0; i<Util::nSigns; ++i){
      for(int j=0; j<Util::nDefinedParticles; ++j){
        mhDeadMaterialEffectTpcRecoEff3D[i][j] = dynamic_cast<TH3F* >( file->Get("hDeadMaterialSDCD"+histNames[i][j]) );
        mhDeadMaterialEffectTpcRecoEff3D[i][j]->SetDirectory(0);
      }
    }
    file->Close();
  } else
    cout << "ERROR in Efficiencies when loading the 3D-dead material secondaries effect file!!!" << endl;

}



Double_t Efficiencies::pidEff(UInt_t pid, Double_t pMax, Double_t pMin) const{
  if(pid >= Util::nDefinedParticles){
    cout << "ERROR in Efficiencies::pidEff(UInt_t pid, Double_t pMax, Double_t pMin) const" << endl;
    return 0.0;
  }
  const TEfficiency* eff = mPidEfficiency[pid][pid];
  double pMaxThis = (pMax > eff->GetPassedHistogram()->GetXaxis()->GetXmax()) ? (eff->GetPassedHistogram()->GetXaxis()->GetXmax()-0.0001) : pMax;
  double pMinThis = (pMin > eff->GetPassedHistogram()->GetYaxis()->GetXmax()) ? (eff->GetPassedHistogram()->GetYaxis()->GetXmax()-0.0001) : pMin;
  int binNumber = eff->FindFixBin(pMaxThis, pMinThis);
  const double pidEffVal = eff->GetEfficiency( binNumber );
  
  return pidEffVal;
}


Double_t Efficiencies::probability_TofL0MultAtLeast1VsClusterSize( unsigned int clusterSize ) const{
  if( clusterSize>6 ) clusterSize = 6;
  return mProbability_TofL0MultAtLeast1VsClusterSize->GetEfficiency( mProbability_TofL0MultAtLeast1VsClusterSize->FindFixBin(clusterSize) );
}

Double_t Efficiencies::probability_TofL0MultAtLeast2VsClusterSize( unsigned int clusterSize ) const{
  if( clusterSize>6 ) clusterSize = 6;
  return mProbability_TofL0MultAtLeast2VsClusterSize->GetEfficiency( mProbability_TofL0MultAtLeast2VsClusterSize->FindFixBin(clusterSize) );
}

Double_t Efficiencies::probability_TofL0MultAtLeast1VsNRecoTofHits_1TofVtx_1TofTrack( unsigned int nRecoHits ) const{
  if( nRecoHits>15 ) nRecoHits = 15;
  return mProbability_TofL0MultAtLeast1VsNRecoTofHits_1TofVtx_1TofTrack->GetEfficiency( mProbability_TofL0MultAtLeast1VsNRecoTofHits_1TofVtx_1TofTrack->FindFixBin(nRecoHits) );
}

Double_t Efficiencies::probability_TofL0MultAtLeast1VsNRecoTofHits_1TofVtx_2TofTracks( unsigned int nRecoHits ) const{
  if( nRecoHits>15 ) nRecoHits = 15;
  return mProbability_TofL0MultAtLeast1VsNRecoTofHits_1TofVtx_2TofTracks->GetEfficiency( mProbability_TofL0MultAtLeast1VsNRecoTofHits_1TofVtx_2TofTracks->FindFixBin(nRecoHits) );
}

Double_t Efficiencies::probability_TofL0MultAtLeast2VsNRecoTofHits_1TofVtx_1TofTrack( unsigned int nRecoHits ) const{
  if( nRecoHits>15 ) nRecoHits = 15;
  return mProbability_TofL0MultAtLeast2VsNRecoTofHits_1TofVtx_1TofTrack->GetEfficiency( mProbability_TofL0MultAtLeast2VsNRecoTofHits_1TofVtx_1TofTrack->FindFixBin(nRecoHits) );
}

Double_t Efficiencies::probability_TofL0MultAtLeast2VsNRecoTofHits_1TofVtx_2TofTracks( unsigned int nRecoHits ) const{
  if( nRecoHits>15 ) nRecoHits = 15;
  return mProbability_TofL0MultAtLeast2VsNRecoTofHits_1TofVtx_2TofTracks->GetEfficiency( mProbability_TofL0MultAtLeast2VsNRecoTofHits_1TofVtx_2TofTracks->FindFixBin(nRecoHits) );
}

Double_t Efficiencies::probability_TofL0MultAtLeast1VsTrackEta_1TofVtx_1TofTrack( double eta ) const{
  if( eta < -1.0 ) eta = -0.99;
  if( eta >  1.0 ) eta =  0.99;
  return mProbability_TofL0MultAtLeast1VsTrackEta_1TofVtx_1TofTrack->GetEfficiency( mProbability_TofL0MultAtLeast1VsTrackEta_1TofVtx_1TofTrack->FindFixBin( eta ) );
}

Double_t Efficiencies::probability_TofL0MultAtLeast2VsTrackEta_1TofVtx_1TofTrack( double eta ) const{
  if( eta < -1.0 ) eta = -0.99;
  if( eta >  1.0 ) eta =  0.99;
  return mProbability_TofL0MultAtLeast2VsTrackEta_1TofVtx_1TofTrack->GetEfficiency( mProbability_TofL0MultAtLeast2VsTrackEta_1TofVtx_1TofTrack->FindFixBin( eta ) );
}


Double_t Efficiencies::tpcEff3DForTF1(Double_t* x, Double_t* p) const{
  UInt_t sign = ((p[0] < 0.5) ? Util::PLUS : Util::MINUS);
  UInt_t pid = ((p[1] < 0.5) ? Util::PION : ((p[1] < 1.5) ? Util::KAON : Util::PROTON) );
  Double_t zVx = p[2];
  Double_t eta = p[3];
  return tpcRecoEff3D(sign, pid, Util::RUN_RANGE_2, zVx, eta, x[0]);
}

Double_t Efficiencies::tofEff3DForTF1(Double_t* x, Double_t* p) const{
  UInt_t sign = ((p[0] < 0.5) ? Util::PLUS : Util::MINUS);
  UInt_t pid = ((p[1] < 0.5) ? Util::PION : ((p[1] < 1.5) ? Util::KAON : Util::PROTON) );
  Double_t zVx = p[2];
  Double_t eta = p[3];
  return tofEff3D(sign, pid, zVx, eta, x[0]);
}





Double_t Efficiencies::onlineAndOfflineVetoEff(UInt_t branchesConf, UInt_t runNumber) const{
  if(branchesConf<Util::nBranchesConfigurations){
    const unsigned int reducedRunNumber = runNumber - 16000000;
    int binNumber = mFullVetoEff[branchesConf]->FindFixBin( reducedRunNumber );
    if( binNumber > 0 && binNumber <= mFullVetoEff[branchesConf]->GetTotalHistogram()->GetNbinsX() ){
      double factor = 1.0;
      #ifdef TEST_GENEX_PION_EFFICIENCIES
      double val = mhEmbeddingVetoCorrectionFactor->GetBinContent( mhEmbeddingVetoCorrectionFactor->FindBin( reducedRunNumber ) );
      double valerr = mhEmbeddingVetoCorrectionFactor->GetBinError( mhEmbeddingVetoCorrectionFactor->FindBin( reducedRunNumber ) );
      if( val>0 && valerr>0 )
        factor = val;
//       cout << (mFullVetoEff[branchesConf]->GetEfficiency( binNumber ) / val) << " "  << mFullVetoEff[branchesConf]->GetEfficiency( binNumber ) << endl;
      #endif
      return mFullVetoEff[branchesConf]->GetEfficiency( binNumber ) / factor;
    } else{ std::cerr << "ERROR in Efficiencies::onlineAndOfflineVetoEff(UInt_t, UInt_t)" << std::endl; return 9e9; }
  } else{ std::cerr << "ERROR in Efficiencies::onlineAndOfflineVetoEff(UInt_t, UInt_t)" << std::endl; return 9e9; }
}



Double_t Efficiencies::purityPid(UInt_t pid, Double_t m2) const{
  if( pid<Util::nDefinedParticles && m2>=0.0 ){
    int binNumber = mhSquaredMassPurity[pid]->FindFixBin( m2 );
    if( m2 > 1.5 )
      binNumber = mhSquaredMassPurity[pid]->GetNbinsX();
    return mhSquaredMassPurity[pid]->GetBinContent( binNumber );
  }else{ std::cerr << "ERROR in Efficiencies::purityPid(UInt_t, Double_t)  m2=" << m2 << std::endl; return 9e9; }
}





Double_t Efficiencies::migrationsAndFakesCorrectionCentralTracks(Double_t pt, Double_t eta) const{
  const double maxPt = mMigrationsAndFakesCorrection_PtEta->GetYaxis()->GetXmax();
  int binNumber = mMigrationsAndFakesCorrection_PtEta->FindBin( eta, pt>maxPt ? 0.99*maxPt : pt );
//   cout << "migrationsAndFakesCorrectionCentralTracks    eta = " << eta << "  pt = " << pt << " factor = " <<  mMigrationsAndFakesCorrection_PtEta->GetBinContent( binNumber ) << endl;
  if( binNumber > 0 && binNumber <= (mMigrationsAndFakesCorrection_PtEta->GetNbinsX()+2)*(mMigrationsAndFakesCorrection_PtEta->GetNbinsY()+2) ){
    return mMigrationsAndFakesCorrection_PtEta->GetBinContent( binNumber );
  }
  else
    return 0.;
}

Double_t Efficiencies::migrationsAndFakesCorrectionForwardProtons(Int_t side, Double_t px, Double_t py) const{
  if( side<Util::nSides && side>=0 ){
    int binNumber = mMigrationsAndFakesCorrection_PxPy_PerSide[side]->FindBin( px, py );
    if( binNumber > 0 && binNumber <= (mMigrationsAndFakesCorrection_PxPy_PerSide[side]->GetNbinsX()+2)*(mMigrationsAndFakesCorrection_PxPy_PerSide[side]->GetNbinsY()+2) ){
      double corr = mMigrationsAndFakesCorrection_PxPy_PerSide[side]->GetBinContent( binNumber );
      if( corr>100 || corr<1e-5 || corr!=corr ){
        cout << "Efficiencies::migrationsAndFakesCorrectionForwardProtons(): corr = " << corr << endl;
        corr = 1.0;
      }
      return corr;
    }
    else{
//       cout << "Efficiencies::migrationsAndFakesCorrectionForwardProtons(): px = " << px << "   py = " << py << "   Returning 1.0 "<< endl;
      return 1.;
    }
  } else{
    cout << "Efficiencies::migrationsAndFakesCorrectionForwardProtons(): error!!!! side = " << side << endl;
    return 0.;
  }
}





Double_t Efficiencies::migrationsAndFakesCorrectionForwardProtons_LimitedMandelstamT(Double_t px, Double_t py) const{
  int binNumber = mMigrationsAndFakesCorrection_PxPy_LimitedMandelstamT->FindBin( px, py );
  if( binNumber > 0 && binNumber <= (mMigrationsAndFakesCorrection_PxPy_LimitedMandelstamT->GetNbinsX()+2)*(mMigrationsAndFakesCorrection_PxPy_LimitedMandelstamT->GetNbinsY()+2) ){
    double corr = mMigrationsAndFakesCorrection_PxPy_LimitedMandelstamT->GetBinContent( binNumber );
    if( corr>100 || corr<1e-5 || corr!=corr ){
      cout << "Efficiencies::migrationsAndFakesCorrectionForwardProtons_LimitedMandelstamT(): corr = " << corr << endl;
      corr = 1.0;
    }
    return corr;
  }
  else return 0.;
}

Double_t Efficiencies::migrationsAndFakesCorrectionForwardProtons_LimitedMandelstamT(Double_t t) const{
  int binNumber = mMigrationsAndFakesCorrection_MandelstamT_LimitedMandelstamT->FindBin( t );
  double corr = mMigrationsAndFakesCorrection_MandelstamT_LimitedMandelstamT->GetBinContent( binNumber );
  return corr;
}

Double_t Efficiencies::migrationsAndFakesCorrectionForwardProtons_LimitedMandelstamT(Int_t br, Double_t t) const{
  if( br<Util::nBranches && br>=0 ){
    int binNumber = mMigrationsAndFakesCorrection_MandelstamT_LimitedMandelstamT_PerBranch[br]->FindBin( t );
    double corr = mMigrationsAndFakesCorrection_MandelstamT_LimitedMandelstamT_PerBranch[br]->GetBinContent( binNumber );
    return corr;
  } else{
    cout << "Efficiencies::migrationsAndFakesCorrectionForwardProtons_LimitedMandelstamT(): error!!!! " << endl;
    return 0.;
  }
}

Double_t Efficiencies::migrationsCorrectionWithinFiducialRegionForwardProtons(UInt_t side, Double_t px, Double_t py) const{
  if(side<Util::nSides){
    int binNumber = mMigrationsWithinFiducialRegionCorrection_PxPy[side]->FindBin( px, py );
    if( binNumber > 0 && binNumber <= (mMigrationsWithinFiducialRegionCorrection_PxPy[side]->GetNbinsX()+2)*(mMigrationsWithinFiducialRegionCorrection_PxPy[side]->GetNbinsY()+2) )
      return mMigrationsWithinFiducialRegionCorrection_PxPy[side]->GetBinContent( binNumber );
    else return 0.;
  } else{ std::cerr << "ERROR in Efficiencies::migrationsCorrectionWithinFiducialRegionForwardProtons(UInt_t, Double_t, Double_t)" << std::endl; return 9e9; }
}


Double_t Efficiencies::trackEffZVtxPxPy(UInt_t branch, Double_t zVtx, Double_t px, Double_t py, Bool_t useCorrection) const{
  if( branch >= Util::nBranches ) return 0.;
  int binNumber = mRpTrackRecoEffiPerBranch[branch]->FindBin( zVtx, px, py );
  if( binNumber > 0 && binNumber <= (mRpTrackRecoEffiPerBranch[branch]->GetNbinsX()+2)*(mRpTrackRecoEffiPerBranch[branch]->GetNbinsY()+2)*(mRpTrackRecoEffiPerBranch[branch]->GetNbinsZ()+2) ){
    double correction = 0.0;
    if( useCorrection ){
      int binX = mRpTrackRecoEffiCorrectionPerBranch[branch]->GetXaxis()->FindBin( px );
      int binY = mRpTrackRecoEffiCorrectionPerBranch[branch]->GetYaxis()->FindBin( py );
      correction = mRpTrackRecoEffiCorrectionPerBranch[branch]->GetBinContent( binX, binY );
    }
    double eff = mRpTrackRecoEffiPerBranch[branch]->GetBinContent( binNumber ) + correction;
    if( eff > 1.0 ) eff = 1.0;
    return eff;
  }
  else
    return 0.;
}



Double_t Efficiencies::rpTrackEffPxPySystematicError(UInt_t branch, Double_t px, Double_t py) const{
  if( branch >= Util::nBranches ) return 0.;
  int binNumber = mRpTrackRecoEffiSystematicsPerBranch[branch]->FindBin( px, py );
  if( binNumber > 0 && binNumber <= (mRpTrackRecoEffiSystematicsPerBranch[branch]->GetNbinsX()+2)*(mRpTrackRecoEffiSystematicsPerBranch[branch]->GetNbinsY()+2) ){
    return mRpTrackRecoEffiSystematicsPerBranch[branch]->GetBinContent( binNumber );
  }
  else
    return 0.;
}


Double_t Efficiencies::rpDeadMatVetoEffPxPySystematicError(UInt_t branch, Double_t px, Double_t py) const{
  if( branch >= Util::nBranches ) return 0.;
  int binNumber = mhDeadMatSystematics[branch]->FindBin( px, py );
  if( binNumber > 0 && binNumber <= (mhDeadMatSystematics[branch]->GetNbinsX()+2)*(mhDeadMatSystematics[branch]->GetNbinsY()+2) ){
    return mhDeadMatSystematics[branch]->GetBinContent( binNumber );
  }
  else
    return 0.;
}


Double_t Efficiencies::trackEffZVtxPxPy_FullRecoEfficiency(UInt_t branch, Double_t zVtx, Double_t px, Double_t py, Bool_t useCorrection) const{
  if( branch >= Util::nBranches ) return 0.;
  int binNumber = mRpTrackRecoEffiPerBranch_FullRecoEfficiency[branch]->FindBin( zVtx, px, py );
  if( binNumber > 0 && binNumber <= (mRpTrackRecoEffiPerBranch_FullRecoEfficiency[branch]->GetNbinsX()+2)*(mRpTrackRecoEffiPerBranch_FullRecoEfficiency[branch]->GetNbinsY()+2)*(mRpTrackRecoEffiPerBranch_FullRecoEfficiency[branch]->GetNbinsZ()+2) ){
//     if( mRpTrackRecoEffiPerBranch_FullRecoEfficiency[branch]->GetBinContent( binNumber ) < 0.0000001 ){
//       int a, b, c;
//       mRpTrackRecoEffiPerBranch_FullRecoEfficiency[branch]->GetBinXYZ( binNumber, a, b, c );
//       cout << "RpTrackRecoEffiPerBranch_FullRecoEfficiency" << endl;
//       cout << "Branch = " << branch << endl;
//       cout << mRpTrackRecoEffiPerBranch_FullRecoEfficiency[branch]->GetXaxis()->GetBinLowEdge( a ) << " < "  << zVtx << " < " << mRpTrackRecoEffiPerBranch_FullRecoEfficiency[branch]->GetXaxis()->GetBinUpEdge( a ) << endl;
//       cout << mRpTrackRecoEffiPerBranch_FullRecoEfficiency[branch]->GetYaxis()->GetBinLowEdge( b ) << " < "  << px << " < " << mRpTrackRecoEffiPerBranch_FullRecoEfficiency[branch]->GetYaxis()->GetBinUpEdge( b ) << endl;
//       cout << mRpTrackRecoEffiPerBranch_FullRecoEfficiency[branch]->GetZaxis()->GetBinLowEdge( c ) << " < "  << py << " < " << mRpTrackRecoEffiPerBranch_FullRecoEfficiency[branch]->GetZaxis()->GetBinUpEdge( c ) << endl;
//       cout << "---------" << endl;
//     }
    double eff = mRpTrackRecoEffiPerBranch_FullRecoEfficiency[branch]->GetBinContent( binNumber );
    if( !(eff>0) ){
      const double binWidthPx = mRpTrackRecoEffiPerBranch_FullRecoEfficiency[branch]->GetYaxis()->GetBinWidth(1);
      const double binWidthPy = mRpTrackRecoEffiPerBranch_FullRecoEfficiency[branch]->GetZaxis()->GetBinWidth(1);
      eff = 0;
      int n = 0;
      for(int i=-1; i<2; ++i)
        for(int j=-1; j<2; ++j){
          if(!i&&!j) continue;
          binNumber = mRpTrackRecoEffiPerBranch_FullRecoEfficiency[branch]->FindBin( zVtx, px+i*binWidthPx, py+j*binWidthPy );
          if( binNumber > 0 && binNumber <= (mRpTrackRecoEffiPerBranch_FullRecoEfficiency[branch]->GetNbinsX()+2)*(mRpTrackRecoEffiPerBranch_FullRecoEfficiency[branch]->GetNbinsY()+2)*(mRpTrackRecoEffiPerBranch_FullRecoEfficiency[branch]->GetNbinsZ()+2) ){
            double effTmp = mRpTrackRecoEffiPerBranch_FullRecoEfficiency[branch]->GetBinContent( binNumber );
            if( effTmp > 0 ){
              eff += effTmp;
              ++n;
            }
          }
        }
      eff /= n;
    } else{
      if( fabs(py)>0.2 && fabs(py)<0.4 && px>0 ){
        int a, b, c;
        mRpTrackRecoEffiPerBranch_FullRecoEfficiency[branch]->GetBinXYZ( binNumber, a, b, c );
        double edges_px[2];
        double edges_py[2];
        edges_px[0] = mRpTrackRecoEffiPerBranch_FullRecoEfficiency[branch]->GetYaxis()->GetBinLowEdge( b );
        edges_py[0] = mRpTrackRecoEffiPerBranch_FullRecoEfficiency[branch]->GetZaxis()->GetBinLowEdge( c );
        edges_px[1] = mRpTrackRecoEffiPerBranch_FullRecoEfficiency[branch]->GetYaxis()->GetBinUpEdge( b );
        edges_py[1] = mRpTrackRecoEffiPerBranch_FullRecoEfficiency[branch]->GetZaxis()->GetBinUpEdge( c );
        double minR = 9999;
        double maxR = -1;
        for(int e1=0; e1<2; ++e1){
          for(int e2=0; e2<2; ++e2){
            const double tmpR = sqrt( pow( edges_px[e1]+0.3, 2) + pow(edges_py[e2], 2) );
            if( minR>tmpR ) minR = tmpR;
            if( maxR<tmpR ) maxR = tmpR;
          }
        }
        if( minR < 0.5 && maxR > 0.5 ){
          const double tmp_px = mRpTrackRecoEffiPerBranch_FullRecoEfficiency[branch]->GetYaxis()->GetBinCenter( b );
          const double tmp_py = mRpTrackRecoEffiPerBranch_FullRecoEfficiency[branch]->GetZaxis()->GetBinCenter( c );
          const double r = sqrt( pow( tmp_px+0.3, 2) + pow(tmp_py, 2) );
          const double pxBinWidth = mRpTrackRecoEffiPerBranch_FullRecoEfficiency[branch]->GetYaxis()->GetBinWidth( b )/2;
          const double pxHalfBinWidth = (sqrt(2.)/2.)*mRpTrackRecoEffiPerBranch_FullRecoEfficiency[branch]->GetYaxis()->GetBinWidth( b );
          const double eff2 = mRpTrackRecoEffiPerBranch_FullRecoEfficiency[branch]->GetBinContent( mRpTrackRecoEffiPerBranch_FullRecoEfficiency[branch]->FindBin( zVtx, px - pxBinWidth, py ) );
          const double f = ((0.5-r)/pxHalfBinWidth + 1.0)/2.;
          eff = f*eff + (1.-f)*eff2;
        }
      }
    }
    
    double correction = 0.0;
    if( useCorrection ){
      int binX = mRpTrackRecoEffiCorrectionPerBranch[branch]->GetXaxis()->FindBin( px );
      int binY = mRpTrackRecoEffiCorrectionPerBranch[branch]->GetYaxis()->FindBin( py );
      correction = mRpTrackRecoEffiCorrectionPerBranch[branch]->GetBinContent( binX, binY );
    }
    
    return eff + correction;
  }
  else
    return 0.;
}

Double_t Efficiencies::trackEffZVtxPxPy_DeadMaterialVetoEff(UInt_t branch, Double_t zVtx, Double_t px, Double_t py, Bool_t useCorrection) const{
  if( branch >= Util::nBranches ) return 0.;
  int binNumber = mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff[branch]->FindBin( zVtx, px, py );
  if( binNumber > 0 && binNumber <= (mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff[branch]->GetNbinsX()+2)*(mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff[branch]->GetNbinsY()+2)*(mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff[branch]->GetNbinsZ()+2) ){
//     if( mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff[branch]->GetBinContent( binNumber ) < 0.0000001 ){
//       int a, b, c;
//       mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff[branch]->GetBinXYZ( binNumber, a, b, c );
//       cout << "RpTrackRecoEffiPerBranchPerTrackType" << endl;
//       cout << "Branch = " << branch << endl;
//       cout << mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff[branch]->GetXaxis()->GetBinLowEdge( a ) << " < "  << zVtx << " < " << mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff[branch]->GetXaxis()->GetBinUpEdge( a ) << endl;
//       cout << mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff[branch]->GetYaxis()->GetBinLowEdge( b ) << " < "  << px << " < " << mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff[branch]->GetYaxis()->GetBinUpEdge( b ) << endl;
//       cout << mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff[branch]->GetZaxis()->GetBinLowEdge( c ) << " < "  << py << " < " << mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff[branch]->GetZaxis()->GetBinUpEdge( c ) << endl;
//       cout << "---------" << endl;
//     }
    double eff = mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff[branch]->GetBinContent( binNumber );
    if( !(eff>0) ){
      const double binWidthPx = mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff[branch]->GetYaxis()->GetBinWidth(1);
      const double binWidthPy = mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff[branch]->GetZaxis()->GetBinWidth(1);
      eff = 0;
      int n = 0;
      for(int i=-1; i<2; ++i)
        for(int j=-1; j<2; ++j){
          if(!i&&!j) continue;
          binNumber = mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff[branch]->FindBin( zVtx, px+i*binWidthPx, py+j*binWidthPy );
          if( binNumber > 0 && binNumber <= (mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff[branch]->GetNbinsX()+2)*(mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff[branch]->GetNbinsY()+2)*(mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff[branch]->GetNbinsZ()+2) ){
            double effTmp = mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff[branch]->GetBinContent( binNumber );
            if( effTmp > 0 ){
              eff += effTmp;
              ++n;
            }
          }
        }
      eff /= n;
    } else{
      if( fabs(py)>0.2 && fabs(py)<0.4 && px>0 ){
        int a, b, c;
        mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff[branch]->GetBinXYZ( binNumber, a, b, c );
        double edges_px[2];
        double edges_py[2];
        edges_px[0] = mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff[branch]->GetYaxis()->GetBinLowEdge( b );
        edges_py[0] = mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff[branch]->GetZaxis()->GetBinLowEdge( c );
        edges_px[1] = mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff[branch]->GetYaxis()->GetBinUpEdge( b );
        edges_py[1] = mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff[branch]->GetZaxis()->GetBinUpEdge( c );
        double minR = 9999;
        double maxR = -1;
        for(int e1=0; e1<2; ++e1){
          for(int e2=0; e2<2; ++e2){
            const double tmpR = sqrt( pow( edges_px[e1]+0.3, 2) + pow(edges_py[e2], 2) );
            if( minR>tmpR ) minR = tmpR;
            if( maxR<tmpR ) maxR = tmpR;
          }
        }
        if( minR < 0.5 && maxR > 0.5 ){
          const double tmp_px = mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff[branch]->GetYaxis()->GetBinCenter( b );
          const double tmp_py = mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff[branch]->GetZaxis()->GetBinCenter( c );
          const double r = sqrt( pow( tmp_px+0.3, 2) + pow(tmp_py, 2) );
          const double pxBinWidth = mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff[branch]->GetYaxis()->GetBinWidth( b )/2;
          const double pxHalfBinWidth = (sqrt(2.)/2.)*mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff[branch]->GetYaxis()->GetBinWidth( b );
          const double eff2 = mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff[branch]->GetBinContent( mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff[branch]->FindBin( zVtx, px - pxBinWidth, py ) );
          const double f = ((0.5-r)/pxHalfBinWidth + 1.0)/2.;
          eff = f*eff + (1.-f)*eff2;
        }
      }
    }
    
    double correction = 1.0;
    if( useCorrection ){
      int binX = mhDeadMatCorrectionMultiplicativeFactor[branch]->GetXaxis()->FindBin( px );
      int binY = mhDeadMatCorrectionMultiplicativeFactor[branch]->GetYaxis()->FindBin( py );
      correction = mhDeadMatCorrectionMultiplicativeFactor[branch]->GetBinContent( binX, binY );
    }
    
    return eff * correction;
    
  }
  else
    return 0.;
}



Double_t Efficiencies::trackEffZVtxPxPy_DeadMaterialVetoEff_ProtonTrackRequired(UInt_t branch, Double_t zVtx, Double_t px, Double_t py, Bool_t useCorrection) const{
  if( branch >= Util::nBranches ) return 0.;
  int binNumber = mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff_ProtonTrackRequired[branch]->FindBin( zVtx, px, py );
  if( binNumber > 0 && binNumber <= (mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff_ProtonTrackRequired[branch]->GetNbinsX()+2)*(mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff_ProtonTrackRequired[branch]->GetNbinsY()+2)*(mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff_ProtonTrackRequired[branch]->GetNbinsZ()+2) ){
    //     if( mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff_ProtonTrackRequired[branch]->GetBinContent( binNumber ) < 0.0000001 ){
    //       int a, b, c;
    //       mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff_ProtonTrackRequired[branch]->GetBinXYZ( binNumber, a, b, c );
    //       cout << "RpTrackRecoEffiPerBranchPerTrackType" << endl;
    //       cout << "Branch = " << branch << endl;
    //       cout << mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff_ProtonTrackRequired[branch]->GetXaxis()->GetBinLowEdge( a ) << " < "  << zVtx << " < " << mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff_ProtonTrackRequired[branch]->GetXaxis()->GetBinUpEdge( a ) << endl;
    //       cout << mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff_ProtonTrackRequired[branch]->GetYaxis()->GetBinLowEdge( b ) << " < "  << px << " < " << mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff_ProtonTrackRequired[branch]->GetYaxis()->GetBinUpEdge( b ) << endl;
    //       cout << mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff_ProtonTrackRequired[branch]->GetZaxis()->GetBinLowEdge( c ) << " < "  << py << " < " << mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff_ProtonTrackRequired[branch]->GetZaxis()->GetBinUpEdge( c ) << endl;
    //       cout << "---------" << endl;
    //     }
    double eff = mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff_ProtonTrackRequired[branch]->GetBinContent( binNumber );
    if( !(eff>0) ){
      const double binWidthPx = mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff_ProtonTrackRequired[branch]->GetYaxis()->GetBinWidth(1);
      const double binWidthPy = mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff_ProtonTrackRequired[branch]->GetZaxis()->GetBinWidth(1);
      eff = 0;
      int n = 0;
      for(int i=-1; i<2; ++i)
        for(int j=-1; j<2; ++j){
          if(!i&&!j) continue;
          binNumber = mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff_ProtonTrackRequired[branch]->FindBin( zVtx, px+i*binWidthPx, py+j*binWidthPy );
          if( binNumber > 0 && binNumber <= (mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff_ProtonTrackRequired[branch]->GetNbinsX()+2)*(mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff_ProtonTrackRequired[branch]->GetNbinsY()+2)*(mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff_ProtonTrackRequired[branch]->GetNbinsZ()+2) ){
            double effTmp = mRpTrackRecoEffiPerBranch_DeadMaterialVetoEff_ProtonTrackRequired[branch]->GetBinContent( binNumber );
            if( effTmp > 0 ){
              eff += effTmp;
              ++n;
            }
          }
        }
        eff /= n;
    }
    
    double correction = 1.0;
    if( useCorrection ){
      int binX = mhDeadMatCorrectionMultiplicativeFactor[branch]->GetXaxis()->FindBin( px );
      int binY = mhDeadMatCorrectionMultiplicativeFactor[branch]->GetYaxis()->FindBin( py );
      correction = mhDeadMatCorrectionMultiplicativeFactor[branch]->GetBinContent( binX, binY );
    }
    
    return eff * correction;
    
  }
  else
    return 0.;
}


Double_t Efficiencies::trackEffTypeZVtxPxPy(UInt_t branch, UInt_t type, Double_t zVtx, Double_t px, Double_t py) const{
  if( branch >= Util::nBranches || type > 1 ) return 0.;
  int binNumber = mRpTrackRecoEffiPerBranchPerTrackType[branch][type]->FindBin( zVtx, px, py );
  if( binNumber > 0 && binNumber <= (mRpTrackRecoEffiPerBranchPerTrackType[branch][type]->GetNbinsX()+2)*(mRpTrackRecoEffiPerBranchPerTrackType[branch][type]->GetNbinsY()+2)*(mRpTrackRecoEffiPerBranchPerTrackType[branch][type]->GetNbinsZ()+2) ){
//     if( mRpTrackRecoEffiPerBranchPerTrackType[branch][type]->GetBinContent( binNumber ) < 0.0000001 ){
//     int a, b, c;
//     mRpTrackRecoEffiPerBranchPerTrackType[branch][type]->GetBinXYZ( binNumber, a, b, c );
//     cout << "Type = " << type << "  branch = " << branch << endl;
//     cout << mRpTrackRecoEffiPerBranchPerTrackType[branch][type]->GetXaxis()->GetBinLowEdge( a ) << " < "  << zVtx << " < " << mRpTrackRecoEffiPerBranchPerTrackType[branch][type]->GetXaxis()->GetBinUpEdge( a ) << endl;
//     cout << mRpTrackRecoEffiPerBranchPerTrackType[branch][type]->GetYaxis()->GetBinLowEdge( b ) << " < "  << px << " < " << mRpTrackRecoEffiPerBranchPerTrackType[branch][type]->GetYaxis()->GetBinUpEdge( b ) << endl;
//     cout << mRpTrackRecoEffiPerBranchPerTrackType[branch][type]->GetZaxis()->GetBinLowEdge( c ) << " < "  << py << " < " << mRpTrackRecoEffiPerBranchPerTrackType[branch][type]->GetZaxis()->GetBinUpEdge( c ) << endl;
//     cout << "---------" << endl;
//     }
    return mRpTrackRecoEffiPerBranchPerTrackType[branch][type]->GetBinContent( binNumber );
  }
  else
    return 0.;
}

TF1* Efficiencies::tpcEff3DTF1(UInt_t sign, UInt_t pid, Double_t zVx, Double_t eta) const{
  if(sign<Util::nSigns && pid<Util::nDefinedParticles){
    TF1 *tempPtr = new TF1( "f", this, &Efficiencies::tpcEff3DForTF1, 0, 20, 4, "Efficiencies","tpcEff3DForTF1");
    tempPtr->FixParameter(0, sign);
    tempPtr->FixParameter(1, pid);
    tempPtr->FixParameter(2, zVx);
    tempPtr->FixParameter(3, eta);
    return tempPtr;
  }
  else{ std::cerr << "ERROR in Efficiencies::tpcEff3DTF1(UInt_t, UInt_t, Double_t, Double_t)" << std::endl; return nullptr; }
}

TF1* Efficiencies::tofEff3DTF1(UInt_t sign, UInt_t pid, Double_t zVx, Double_t eta) const{
  if(sign<Util::nSigns && pid<Util::nDefinedParticles){
    TF1 *tempPtr = new TF1( "f", this, &Efficiencies::tofEff3DForTF1, 0, 20, 4, "Efficiencies","tofEff3DForTF1");
    tempPtr->FixParameter(0, sign);
    tempPtr->FixParameter(1, pid);
    tempPtr->FixParameter(2, zVx);
    tempPtr->FixParameter(3, eta);
    return tempPtr;
  }
  else{ std::cerr << "ERROR in Efficiencies::tofEff3DTF1(UInt_t, UInt_t, Double_t, Double_t)" << std::endl; return nullptr; }
}



Double_t Efficiencies::tpcRecoEff3D(UInt_t sign, UInt_t pid, UInt_t rr, Double_t zVx, Double_t eta, Double_t pt, Bool_t recoCuts) const{
  TH3F *hist = recoCuts ? mhTpcRecoEff3D_withRecoCuts[sign][pid] : mhTpcRecoEff3D[sign][pid][rr];
  double newZVx = (hist->GetXaxis()->GetBinCenter(1) > zVx ? hist->GetXaxis()->GetBinCenter(1)+1e-3 : (hist->GetXaxis()->GetBinCenter(hist->GetNbinsX()) < zVx ? hist->GetXaxis()->GetBinCenter(hist->GetNbinsX())-1e-3 : zVx));
  double newEta = (hist->GetZaxis()->GetBinCenter(1) > eta ? hist->GetZaxis()->GetBinCenter(1)+1e-3 : (hist->GetZaxis()->GetBinCenter(hist->GetNbinsZ()) < eta ? hist->GetZaxis()->GetBinCenter(hist->GetNbinsZ())-1e-3 : eta));
  double newPt = (hist->GetYaxis()->GetBinCenter(1) > (pt + mTpcRecoEff3D_InflationShift[sign][pid]) ? hist->GetYaxis()->GetBinCenter(1)+1e-3 : (hist->GetYaxis()->GetBinCenter(hist->GetNbinsY()) < (pt + mTpcRecoEff3D_InflationShift[sign][pid]) ? hist->GetYaxis()->GetBinCenter(hist->GetNbinsY())-1e-3 : (pt + mTpcRecoEff3D_InflationShift[sign][pid])));
//   if( mTpcRecoEff3D_MultiplicationFactor[sign][pid] * hist->Interpolate( newZVx, newPt, newEta ) < 1e-7 ) cout << "TPC EFF = 0     " << zVx << " " << eta << " " << pt << "     " << hist->Interpolate( newZVx, newPt, newEta ) << endl;
  if( mTpcRecoEff3D_isUnchanged[sign][pid] )
    return mInterpolate ? hist->Interpolate( newZVx, newPt, newEta ) : hist->GetBinContent( hist->FindBin( newZVx, newPt, newEta) );
  else{
    double nominalEff = hist->Interpolate( newZVx, newPt, newEta );
    double highPtEff = ( hist->Interpolate( newZVx, hist->GetYaxis()->GetBinCenter(hist->GetNbinsY())-1e-3, newEta) + 
			 hist->Interpolate( newZVx, hist->GetYaxis()->GetBinCenter(hist->GetNbinsY()-1), newEta ) +
			 hist->Interpolate( newZVx, hist->GetYaxis()->GetBinCenter(hist->GetNbinsY()-2), newEta ) )/3.; // averaging of three last high-pT bins (efficiency plateu) to reduce fluctuations
    double efficiencyScaling = (highPtEff + mTpcRecoEff3D_PlateuShift[sign][pid]) / highPtEff;
    return efficiencyScaling * ( mInterpolate ? hist->Interpolate( newZVx, newPt, newEta ) : hist->GetBinContent( hist->FindBin( newZVx, newPt, newEta) ) );
  }
}


Double_t Efficiencies::tpcRecoEff3D_phi(UInt_t sign, UInt_t pid, UInt_t rr, Double_t phi, Double_t eta, Double_t pt) const{
  TH3F *hist = mhTpcRecoEff3D_Phi[sign][pid][rr];
  double newPhi = (phi + Util::PI/12.) < 3.14 ? (phi + Util::PI/12.) : -3.10;
  double newEta = (hist->GetZaxis()->GetBinCenter(1) > eta ? hist->GetZaxis()->GetBinCenter(1)+1e-3 : (hist->GetZaxis()->GetBinCenter(hist->GetNbinsZ()) < eta ? hist->GetZaxis()->GetBinCenter(hist->GetNbinsZ())-1e-3 : eta));
  double newPt = (hist->GetYaxis()->GetBinCenter(1) > (pt + mTpcRecoEff3D_InflationShift[sign][pid]) ? hist->GetYaxis()->GetBinCenter(1)+1e-3 : (hist->GetYaxis()->GetBinCenter(hist->GetNbinsY()) < (pt + mTpcRecoEff3D_InflationShift[sign][pid]) ? hist->GetYaxis()->GetBinCenter(hist->GetNbinsY())-1e-3 : (pt + mTpcRecoEff3D_InflationShift[sign][pid])));
  //   if( mTpcRecoEff3D_MultiplicationFactor[sign][pid] * hist->Interpolate( newPhi, newPt, newEta ) < 1e-7 ) cout << "TPC EFF = 0     " << zVx << " " << eta << " " << pt << "     " << hist->Interpolate( newPhi, newPt, newEta ) << endl;
  if( mTpcRecoEff3D_isUnchanged[sign][pid] )
    return mInterpolate ? hist->Interpolate( newPhi, newPt, newEta ) : hist->GetBinContent( hist->FindBin( newPhi, newPt, newEta) );
  else{
    double nominalEff = hist->Interpolate( newPhi, newPt, newEta );
    double highPtEff = ( hist->Interpolate( newPhi, hist->GetYaxis()->GetBinCenter(hist->GetNbinsY())-1e-3, newEta) + 
    hist->Interpolate( newPhi, hist->GetYaxis()->GetBinCenter(hist->GetNbinsY()-1), newEta ) +
    hist->Interpolate( newPhi, hist->GetYaxis()->GetBinCenter(hist->GetNbinsY()-2), newEta ) )/3.; // averaging of three last high-pT bins (efficiency plateu) to reduce fluctuations
    double efficiencyScaling = (highPtEff + mTpcRecoEff3D_PlateuShift[sign][pid]) / highPtEff;
    return efficiencyScaling * ( mInterpolate ? hist->Interpolate( newPhi, newPt, newEta ) : hist->GetBinContent( hist->FindBin( newPhi, newPt, newEta) ) );
  }
}


Double_t Efficiencies::tofEff3D(UInt_t sign, UInt_t pid, Double_t zVx, Double_t eta, Double_t pt, Bool_t recoCuts) const{
  TH3F *hist = recoCuts ? mhTofEff3D_withRecoCuts[sign][pid] : mhTofEff3D[sign][pid];
  double newZVx = (hist->GetXaxis()->GetBinCenter(1) > zVx ? hist->GetXaxis()->GetBinCenter(1)+1e-3 : (hist->GetXaxis()->GetBinCenter(hist->GetNbinsX()) < zVx ? hist->GetXaxis()->GetBinCenter(hist->GetNbinsX())-1e-3 : zVx));
  double newEta = (hist->GetZaxis()->GetBinCenter(1) > eta ? hist->GetZaxis()->GetBinCenter(1)+1e-3 : (hist->GetZaxis()->GetBinCenter(hist->GetNbinsZ()) < eta ? hist->GetZaxis()->GetBinCenter(hist->GetNbinsZ())-1e-3 : eta));
  double newPt = (hist->GetYaxis()->GetBinCenter(1) > (pt + mTofEff3D_InflationShift[sign][pid]) ? hist->GetYaxis()->GetBinCenter(1)+1e-3 : (hist->GetYaxis()->GetBinCenter(hist->GetNbinsY()) < (pt + mTofEff3D_InflationShift[sign][pid]) ? hist->GetYaxis()->GetBinCenter(hist->GetNbinsY())-1e-3 : (pt + mTofEff3D_InflationShift[sign][pid])));
//   if( mTofEff3D_MultiplicationFactor[sign][pid] * hist->Interpolate( newZVx, newPt, newEta ) < 1e-7 ) cout << "TOF EFF = 0     " << zVx << " " << eta << " " << pt << "     " << hist->Interpolate( newZVx, newPt, newEta ) << endl;
  if( mTofEff3D_isUnchanged[sign][pid] ){
//     cout << "Eff comparison:\t\tBin Content: " << hist->GetBinContent( hist->FindBin( newZVx, newPt, newEta) ) << "\tStandard ROOT: " << hist->Interpolate( newZVx, newPt, newEta ) << "\tTricubic: " << u_out << " " << v_out << " " << w_out << endl;
    return mInterpolate ? hist->Interpolate( newZVx, newPt, newEta ) : hist->GetBinContent( hist->FindBin( newZVx, newPt, newEta) );
  }
  else{
    double nominalEff = hist->Interpolate( newZVx, newPt, newEta );
    double highPtEff = ( hist->Interpolate( newZVx, hist->GetYaxis()->GetBinCenter(hist->GetNbinsY())-1e-3, newEta) + 
			 hist->Interpolate( newZVx, hist->GetYaxis()->GetBinCenter(hist->GetNbinsY()-1), newEta ) +
			 hist->Interpolate( newZVx, hist->GetYaxis()->GetBinCenter(hist->GetNbinsY()-2), newEta ) )/3.; // averaging of three last high-pT bins (efficiency plateu) to reduce fluctuations
    double efficiencyScaling = (highPtEff + mTofEff3D_PlateuShift[sign][pid]) / highPtEff;
    return efficiencyScaling * ( mInterpolate ? hist->Interpolate( newZVx, newPt, newEta ) : hist->GetBinContent( hist->FindBin( newZVx, newPt, newEta) ) );
  }
}



Double_t Efficiencies::tofEff3D_phi(UInt_t sign, UInt_t pid, Double_t phi, Double_t eta, Double_t pt) const{
  TH3F *hist = mhTofEff3D_Phi[sign][pid];
  double newPhi = (phi + Util::PI/12.) < 3.14 ? (phi + Util::PI/12.) : -3.10;
  double newEta = (hist->GetZaxis()->GetBinCenter(1) > eta ? hist->GetZaxis()->GetBinCenter(1)+1e-3 : (hist->GetZaxis()->GetBinCenter(hist->GetNbinsZ()) < eta ? hist->GetZaxis()->GetBinCenter(hist->GetNbinsZ())-1e-3 : eta));
  double newPt = (hist->GetYaxis()->GetBinCenter(1) > (pt + mTofEff3D_InflationShift[sign][pid]) ? hist->GetYaxis()->GetBinCenter(1)+1e-3 : (hist->GetYaxis()->GetBinCenter(hist->GetNbinsY()) < (pt + mTofEff3D_InflationShift[sign][pid]) ? hist->GetYaxis()->GetBinCenter(hist->GetNbinsY())-1e-3 : (pt + mTofEff3D_InflationShift[sign][pid])));
  //   if( mTofEff3D_MultiplicationFactor[sign][pid] * hist->Interpolate( newPhi, newPt, newEta ) < 1e-7 ) cout << "TOF EFF = 0     " << zVx << " " << eta << " " << pt << "     " << hist->Interpolate( newPhi, newPt, newEta ) << endl;
  if( mTofEff3D_isUnchanged[sign][pid] ){
    //     cout << "Eff comparison:\t\tBin Content: " << hist->GetBinContent( hist->FindBin( newPhi, newPt, newEta) ) << "\tStandard ROOT: " << hist->Interpolate( newPhi, newPt, newEta ) << "\tTricubic: " << u_out << " " << v_out << " " << w_out << endl;
    return mInterpolate ? hist->Interpolate( newPhi, newPt, newEta ) : hist->GetBinContent( hist->FindBin( newPhi, newPt, newEta) );
  }
  else{
    double nominalEff = hist->Interpolate( newPhi, newPt, newEta );
    double highPtEff = ( hist->Interpolate( newPhi, hist->GetYaxis()->GetBinCenter(hist->GetNbinsY())-1e-3, newEta) + 
    hist->Interpolate( newPhi, hist->GetYaxis()->GetBinCenter(hist->GetNbinsY()-1), newEta ) +
    hist->Interpolate( newPhi, hist->GetYaxis()->GetBinCenter(hist->GetNbinsY()-2), newEta ) )/3.; // averaging of three last high-pT bins (efficiency plateu) to reduce fluctuations
    double efficiencyScaling = (highPtEff + mTofEff3D_PlateuShift[sign][pid]) / highPtEff;
    return efficiencyScaling * ( mInterpolate ? hist->Interpolate( newPhi, newPt, newEta ) : hist->GetBinContent( hist->FindBin( newPhi, newPt, newEta) ) );
  }
}


Double_t Efficiencies::tofEff2DCorrection( UInt_t pid, Double_t eta, Double_t pt ) const{
  TH2F *hist = mhTofEffCorrection2D_ptVsEta[pid];
  double newEta = (hist->GetXaxis()->GetBinCenter(1) > eta ? hist->GetXaxis()->GetBinCenter(1)+1e-3 : (hist->GetXaxis()->GetBinCenter(hist->GetNbinsX()) < eta ? hist->GetXaxis()->GetBinCenter(hist->GetNbinsX())-1e-3 : eta));
  double newPt = (hist->GetYaxis()->GetBinCenter(1) > pt ? hist->GetYaxis()->GetBinCenter(1)+1e-3 : (hist->GetYaxis()->GetBinCenter(hist->GetNbinsY()) < pt ? hist->GetYaxis()->GetBinCenter(hist->GetNbinsY())-1e-3 : pt));
  return hist->GetBinContent( hist->FindBin( newEta, newPt ) );
}

Double_t Efficiencies::tofEff2DExtraCorrection( UInt_t pid, Double_t eta, Double_t pt ) const{
  TH2D *hist = mhTofEffExtraCorr2D_ptVsEta[pid];
  double newEta = (hist->GetXaxis()->GetBinCenter(1) > eta ? hist->GetXaxis()->GetBinCenter(1)+1e-3 : (hist->GetXaxis()->GetBinCenter(hist->GetNbinsX()) < eta ? hist->GetXaxis()->GetBinCenter(hist->GetNbinsX())-1e-3 : eta));
  double newPt = (hist->GetYaxis()->GetBinCenter(1) > pt ? hist->GetYaxis()->GetBinCenter(1)+1e-3 : (hist->GetYaxis()->GetBinCenter(hist->GetNbinsY()) < pt ? hist->GetYaxis()->GetBinCenter(hist->GetNbinsY())-1e-3 : pt));
  return hist->GetBinContent( hist->FindBin( newEta, newPt ) );
}

Double_t Efficiencies::tofEff2DSystematicError( UInt_t pid, Double_t eta, Double_t pt ) const{
  TH2D *hist = mhTofEffSystErr2D_ptVsEta[pid];
  double newEta = (hist->GetXaxis()->GetBinCenter(1) > eta ? hist->GetXaxis()->GetBinCenter(1)+1e-3 : (hist->GetXaxis()->GetBinCenter(hist->GetNbinsX()) < eta ? hist->GetXaxis()->GetBinCenter(hist->GetNbinsX())-1e-3 : eta));
  double newPt = (hist->GetYaxis()->GetBinCenter(1) > pt ? hist->GetYaxis()->GetBinCenter(1)+1e-3 : (hist->GetYaxis()->GetBinCenter(hist->GetNbinsY()) < pt ? hist->GetYaxis()->GetBinCenter(hist->GetNbinsY())-1e-3 : pt));
  return hist->GetBinContent( hist->FindBin( newEta, newPt ) );
}

Double_t Efficiencies::tofEffDifferenceHftMC3D( UInt_t pid, Double_t zVx, Double_t eta, Double_t pt ) const{
    TH3F *hist = mhTofEffDifferenceHftMC3D_ptVsEtaVsZVtx[pid];
    double newZVx = (hist->GetXaxis()->GetBinCenter(1) > zVx ? hist->GetXaxis()->GetBinCenter(1)+1e-3 : (hist->GetXaxis()->GetBinCenter(hist->GetNbinsX()) < zVx ? hist->GetXaxis()->GetBinCenter(hist->GetNbinsX())-1e-3 : zVx));
    double newEta = (hist->GetYaxis()->GetBinCenter(1) > eta ? hist->GetYaxis()->GetBinCenter(1)+1e-3 : (hist->GetYaxis()->GetBinCenter(hist->GetNbinsZ()) < eta ? hist->GetYaxis()->GetBinCenter(hist->GetNbinsZ())-1e-3 : eta));
    double newPt = (hist->GetZaxis()->GetBinCenter(1) > pt ? hist->GetZaxis()->GetBinCenter(1)+1e-3 : (hist->GetZaxis()->GetBinCenter(hist->GetNbinsY()) < pt ? hist->GetZaxis()->GetBinCenter(hist->GetNbinsY())-1e-3 : pt));
    return /*mInterpolate ? hist->Interpolate( newZVx, newEta, newPt ) :*/ hist->GetBinContent( hist->FindBin( newZVx, newEta, newPt ) );
}

Double_t Efficiencies::deadMatEffectTpcRecoEff3D(UInt_t sign, UInt_t pid, Double_t zVx, Double_t eta, Double_t pt) const{
  TH3F *hist = mhDeadMaterialEffectTpcRecoEff3D[sign][pid];
  double newZVx = (hist->GetXaxis()->GetBinCenter(1) > zVx ? hist->GetXaxis()->GetBinCenter(1)+1e-3 : (hist->GetXaxis()->GetBinCenter(hist->GetNbinsX()) < zVx ? hist->GetXaxis()->GetBinCenter(hist->GetNbinsX())-1e-3 : zVx));
  double newEta = (hist->GetZaxis()->GetBinCenter(1) > eta ? hist->GetZaxis()->GetBinCenter(1)+1e-3 : (hist->GetZaxis()->GetBinCenter(hist->GetNbinsZ()) < eta ? hist->GetZaxis()->GetBinCenter(hist->GetNbinsZ())-1e-3 : eta));
  double newPt = (hist->GetYaxis()->GetBinCenter(1) > pt ? hist->GetYaxis()->GetBinCenter(1)+1e-3 : (hist->GetYaxis()->GetBinCenter(hist->GetNbinsY()) < pt ? hist->GetYaxis()->GetBinCenter(hist->GetNbinsY())-1e-3 : pt));
  return hist->GetBinContent( hist->FindBin( newZVx, newPt, newEta) );
}


Double_t Efficiencies::getLowHighPileUpSystematics(Double_t pt, /*Double_t eta, */UInt_t runRange, UInt_t charge) const{
  const TH1D* hist = mhSystematicsPileUpTPC[runRange][charge][0];
  const TH1D* hist1 = mhSystematicsPileUpTPC[runRange][charge][1];
  double newPt = (hist->GetXaxis()->GetBinCenter(1) > pt ? hist->GetXaxis()->GetBinCenter(1)+1e-3 : (hist->GetXaxis()->GetBinCenter(hist->GetNbinsX()) < pt ? hist->GetXaxis()->GetBinCenter(hist->GetNbinsX())-1e-3 : pt));
  double value0 = fabs(hist->GetBinContent(hist->GetXaxis()->FindBin(newPt)));
  double value1 = fabs(hist1->GetBinContent(hist1->GetXaxis()->FindBin(newPt)));
  
  return value0 > value1 ? value0 : value1;
}




Double_t Efficiencies::zVertexCutEff(UInt_t fillNumber, Double_t minZVtx, Double_t maxZVtx, Double_t shiftZVtxMean, Double_t shiftZVtxSigma) const{
  int binNumber = mhZVtxMeanVsFill->GetXaxis()->FindBin(fillNumber);
  if( binNumber>0 && binNumber<=mhZVtxMeanVsFill->GetNbinsX() && minZVtx < maxZVtx ){
    const double zVtxMean = mhZVtxMeanVsFill->GetBinContent( binNumber ) + shiftZVtxMean;
    const double zVtxSigma = mhZVtxSigmaVsFill->GetBinContent( binNumber ) + shiftZVtxSigma;
    double eff = 0.5 * (TMath::Erf( (maxZVtx-zVtxMean) / (sqrt(2.0)*zVtxSigma) ) - TMath::Erf( (minZVtx-zVtxMean) / (sqrt(2.0)*zVtxSigma) ));
    return eff;
  } else{
    std::cerr << "ERROR in Efficiencies::zVertexCutEff(UInt_t, Double_t, Double_t, Double_t=0, Double_t=0)" << std::endl; return -1.; 
  }
}

Double_t Efficiencies::zVertexCutEff(Double_t minZVtx, Double_t maxZVtx, Double_t zVtxMean, Double_t zVtxSigma) const{
  double eff = 0.5 * (TMath::Erf( (maxZVtx-zVtxMean) / (sqrt(2.0)*zVtxSigma) ) - TMath::Erf( (minZVtx-zVtxMean) / (sqrt(2.0)*zVtxSigma) ));
  return eff;
}



Double_t Efficiencies::correctedDEdxMC(UInt_t pid, Double_t dEdx, Double_t p) const{
  if(pid<Util::nDefinedParticlesExtended){
    Util *util = Util::instance();
    const double sigma_Data = mDEdxWidthParametersData[pid][0] + mDEdxWidthParametersData[pid][1]*std::exp(-mDEdxWidthParametersData[pid][2]*p) + mDEdxWidthParametersData[pid][3]*std::atan(mDEdxWidthParametersData[pid][4]*(p - mDEdxWidthParametersData[pid][5]));
    const double sigma_MC = mDEdxWidthParametersMC[pid][0] + mDEdxWidthParametersMC[pid][1]*std::exp(-mDEdxWidthParametersMC[pid][2]*p) + mDEdxWidthParametersMC[pid][3]*std::atan(mDEdxWidthParametersMC[pid][4]*(p - mDEdxWidthParametersMC[pid][5]));
    const double meanOffset_Data = mDEdxMeanOffsetParametersData[pid][0] + mDEdxMeanOffsetParametersData[pid][1]*std::exp(-mDEdxMeanOffsetParametersData[pid][2]*p) + mDEdxMeanOffsetParametersData[pid][3]*std::atan(mDEdxMeanOffsetParametersData[pid][4]*(p - mDEdxMeanOffsetParametersData[pid][5]));
    const double meanOffset_MC = mDEdxMeanOffsetParametersMC[pid][0] + mDEdxMeanOffsetParametersMC[pid][1]*std::exp(-mDEdxMeanOffsetParametersMC[pid][2]*p) + mDEdxMeanOffsetParametersMC[pid][3]*std::atan(mDEdxMeanOffsetParametersMC[pid][4]*(p - mDEdxMeanOffsetParametersMC[pid][5]));
    const double BetheBloch_MPV = 1e-6*std::exp( util->bichsel()->dedx()->GetMostProbableZ(std::log10(p/util->massExtended(static_cast<Util::PARTICLE_NAME_EXTENDED>(pid))),1.) );
    const double mean_Data = BetheBloch_MPV - meanOffset_Data;
    const double mean_MC = BetheBloch_MPV - meanOffset_MC;
    const double a = sigma_Data / sigma_MC;
    const double b = mean_Data / std::pow( mean_MC, a );
    return b * std::pow(dEdx, a);
  }
  else{ std::cerr << "ERROR in Efficiencies::correctedDEdxMC(UInt_t, Double_t, Double_t)" << std::endl; return dEdx; }
}


Double_t Efficiencies::correctedNSigmaMC(UInt_t nSigmaPid, UInt_t pid, Double_t dEdx, Double_t p, Double_t dEdxError, Double_t log2dx) const{
  if(pid<Util::nDefinedParticlesExtended && nSigmaPid<Util::nDefinedParticles){
    Util *util = Util::instance();
    const double dEdx_corrected = correctedDEdxMC(pid, dEdx, p);
//     const double sigma_Data = mDEdxWidthParametersData[pid][0] + mDEdxWidthParametersData[pid][1]*std::exp(-mDEdxWidthParametersData[pid][2]*p) + mDEdxWidthParametersData[pid][3]*std::atan(mDEdxWidthParametersData[pid][4]*(p - mDEdxWidthParametersData[pid][5]));
//     const double sigma_MC = mDEdxWidthParametersMC[pid][0] + mDEdxWidthParametersMC[pid][1]*std::exp(-mDEdxWidthParametersMC[pid][2]*p) + mDEdxWidthParametersMC[pid][3]*std::atan(mDEdxWidthParametersMC[pid][4]*(p - mDEdxWidthParametersMC[pid][5]));
    const double BetheBloch_MPV = 1e-6*std::exp( util->bichsel()->dedx()->GetMostProbableZ(std::log10(p/util->mass(static_cast<Util::PARTICLE_NAME>(nSigmaPid))), log2dx) );
    return std::log(dEdx_corrected/BetheBloch_MPV) / (dEdxError /** (sigma_Data/sigma_MC)*/);
  }
  else{ std::cerr << "ERROR in Efficiencies::correctedNSigmaMC(UInt_t, UInt_t, Double_t, Double_t)" << std::endl; return 9999; }
}


Double_t Efficiencies::meanDEdx(UInt_t inputType, UInt_t pid, Double_t p) const{
  if(pid<Util::nDefinedParticlesExtended && inputType<Util::nInputTypes){
    Util *util = Util::instance();
    const double* parametersArray = (inputType==Util::DATA) ? mDEdxMeanOffsetParametersData[pid] : mDEdxMeanOffsetParametersMC[pid];
    double meanOffset = parametersArray[0] + parametersArray[1]*std::exp(-parametersArray[2]*p) + parametersArray[3]*std::atan(parametersArray[4]*(p - parametersArray[5]));
    double BetheBloch_MPV = 1e-6*std::exp( util->bichsel()->dedx()->GetMostProbableZ(std::log10(p/util->massExtended(static_cast<Util::PARTICLE_NAME_EXTENDED>(pid))),1.) );
    return BetheBloch_MPV - meanOffset;
  }
  else{ std::cerr << "ERROR in Efficiencies::meanDEdx(UInt_t, UInt_t, Double_t)" << std::endl; return 9999; }
}

Double_t Efficiencies::sigmaDEdx(UInt_t inputType, UInt_t pid, Double_t p) const{
  if(pid<Util::nDefinedParticlesExtended && inputType<Util::nInputTypes){
    Util *util = Util::instance();
    const double* parametersArray = (inputType==Util::DATA) ? mDEdxWidthParametersData[pid] : mDEdxWidthParametersMC[pid];
    return parametersArray[0] + parametersArray[1]*std::exp(-parametersArray[2]*p) + parametersArray[3]*std::atan(parametersArray[4]*(p - parametersArray[5]));
  }
  else{ std::cerr << "ERROR in Efficiencies::sigmaDEdx(UInt_t, UInt_t, Double_t)" << std::endl; return 9999; }
}


TVector3 Efficiencies::correctedMomentumVec(UInt_t sign, UInt_t pid, Double_t zVtx, const TVector3 &p) const{
  if(sign<Util::nSigns && pid<Util::nDefinedParticles){
    int Vzbin = 19;
    for(int i=-10; i<10; ++i)
      if(zVtx > i*10. && zVtx < (i+1)*10.){
	Vzbin = i+10;
	break;
      }
    const double rawPt = p.Pt();
    const double eta = p.Eta();
    mEnergyLossCorrection3D[sign][pid][Vzbin]->GetXaxis()->SetRange( mEnergyLossCorrection3D[sign][pid][Vzbin]->GetXaxis()->FindBin(rawPt), mEnergyLossCorrection3D[sign][pid][Vzbin]->GetXaxis()->FindBin(rawPt) );
    mEnergyLossCorrection3D[sign][pid][Vzbin]->GetYaxis()->SetRange( mEnergyLossCorrection3D[sign][pid][Vzbin]->GetYaxis()->FindBin(eta), mEnergyLossCorrection3D[sign][pid][Vzbin]->GetYaxis()->FindBin(eta) );
    double corrPt = rawPt - mEnergyLossCorrection3D[sign][pid][Vzbin]->GetMean(3);
    double correctionFactor = corrPt / rawPt;
    return TVector3( correctionFactor*p[0], correctionFactor*p[1], correctionFactor*p[2] );
  }
  else{ std::cerr << "ERROR in Efficiencies::correctedMomentumVec(UInt_t, UInt_t, TVector3 &)" << std::endl; return TVector3(); }
}


Double_t Efficiencies::vertexingEff(Double_t deltaZ0, Bool_t isData) const{
  TEfficiency *efObj = isData ? mVertexingEfficiencyVsDeltaZ0 : mVertexingEfficiencyVsDeltaZ0_MC;
  return efObj->GetEfficiency( efObj->FindFixBin( std::fabs(deltaZ0) ) );
}


Double_t Efficiencies::deltaZ0CutEff(Double_t lowerPt, Bool_t isData) const{
  if( lowerPt < 0.2 ) lowerPt = 0.2001; else if( lowerPt > 2.0 ) lowerPt = 1.999;
  TEfficiency *efObj = isData ? mDeltaZ0CutEfficiencyVsLowerPt : mDeltaZ0CutEfficiencyVsLowerPt_MC;
  return efObj->GetEfficiency( efObj->FindFixBin( lowerPt ) );
}


Double_t Efficiencies::missingPtCutEff(Double_t pmax, Double_t pmin, Bool_t isData) const{
  TEfficiency *effObj = isData ? mPtMissCutEff : mPtMissCutEff_MC;
  const double eff = effObj->GetEfficiency( effObj->FindFixBin( min(pmax, 2.99), min(pmin, 2.99) ) );
  if( eff > 0.01 ) return eff;
  else{ std::cerr << "ERROR in Efficiencies::missingPtCutEff(Double_t, Double_t, Bool_t)" << std::endl; return (isData?0.01:1.0); }
}


Double_t Efficiencies::deltaPhiCutEff(UInt_t deltaPhiBin, Double_t t1, Double_t t2) const{
  if(deltaPhiBin<Util::nDeltaPhiRanges){
    const double eff = mAcceptanceForDeltaPhiCuts_Vs_t1t2[deltaPhiBin]->GetEfficiency( mAcceptanceForDeltaPhiCuts_Vs_t1t2[deltaPhiBin]->FindFixBin( t1, t2 ) );
    return eff;
  }
  else{ std::cerr << "ERROR in Efficiencies::deltaPhiCutEff(UInt_t, Double_t, Double_t)" << std::endl; return 0.00000000001; }
}

Double_t Efficiencies::deltaPhiCutEff_GenEx(UInt_t deltaPhiBin, Double_t t1, Double_t t2) const{
  if(deltaPhiBin<Util::nDeltaPhiRanges){
    const double eff = mAcceptanceForDeltaPhiCuts_Vs_t1t2_GenEx[deltaPhiBin]->GetEfficiency( mAcceptanceForDeltaPhiCuts_Vs_t1t2_GenEx[deltaPhiBin]->FindFixBin( t1, t2 ) );
    return eff;
  }
  else{ std::cerr << "ERROR in Efficiencies::deltaPhiCutEff_GenEx(UInt_t, Double_t, Double_t)" << std::endl; return 0.00000000001; }
}

Double_t Efficiencies::deltaPhiCutEff_DiMe(UInt_t deltaPhiBin, Double_t t1, Double_t t2) const{
  if(deltaPhiBin<Util::nDeltaPhiRanges){
    const double eff = mAcceptanceForDeltaPhiCuts_Vs_t1t2_DiMe[deltaPhiBin]->GetEfficiency( mAcceptanceForDeltaPhiCuts_Vs_t1t2_DiMe[deltaPhiBin]->FindFixBin( t1, t2 ) );
    return eff;
  }
  else{ std::cerr << "ERROR in Efficiencies::deltaPhiCutEff_DiMe(UInt_t, Double_t, Double_t)" << std::endl; return 0.00000000001; }
}


Double_t Efficiencies::goodRecoProtonProbability(UInt_t branchConf) const{
  if(branchConf<Util::nBranchesConfigurations){
      return mhRpCoincidence->GetBinContent( mhRpCoincidence->GetXaxis()->FindBin( static_cast<double>(branchConf) ) );
  }
  else{ std::cerr << "ERROR in Efficiencies::goodRecoProtonProbability(UInt_t)" << std::endl; return 0.0; }
}


Double_t Efficiencies::simultaneousBothSidesProtonLossProbability(UInt_t branchConf) const{
  if(branchConf<Util::nBranchesConfigurations){
      return mhSimultaneousEastWestRpTrackLossProbability->GetBinContent( mhSimultaneousEastWestRpTrackLossProbability->GetXaxis()->FindBin( static_cast<double>(branchConf) ) );
  }
  else{ std::cerr << "ERROR in Efficiencies::simultaneousBothSidesProtonLossProbability(UInt_t)" << std::endl; return 0.0; }
}

Double_t Efficiencies::rpTrackRecoEffEastWestCorrelation(UInt_t branchConf) const{
  if(branchConf<Util::nBranchesConfigurations){
    return mhEastWestRpTrackEffCorrelation->GetBinContent( mhEastWestRpTrackEffCorrelation->GetXaxis()->FindBin( static_cast<double>(branchConf) ) );
  }
  else{ std::cerr << "ERROR in Efficiencies::rpTrackRecoEffEastWestCorrelation(UInt_t)" << std::endl; return 0.0; }
}

Double_t Efficiencies::tpcTotalEffPlusMinusCorrelation(Double_t avrgPt) const{
  return mhPlusMinusTpcTrackTotalEfficiencyCorrelationVsAveragePt->Interpolate( avrgPt );
//   return mhPlusMinusTpcTrackTotalEfficiencyCorrelationVsAveragePt->GetBinContent( mhPlusMinusTpcTrackTotalEfficiencyCorrelationVsAveragePt->FindBin( avrgPt ) );
}

Double_t Efficiencies::tpcRecoEffPlusMinusCorrelation(Double_t avrgPt) const{
  return mhPlusMinusTpcTrackRecoEfficiencyCorrelationVsAveragePt->Interpolate( avrgPt );
//   return mhPlusMinusTpcTrackRecoEfficiencyCorrelationVsAveragePt->GetBinContent( mhPlusMinusTpcTrackRecoEfficiencyCorrelationVsAveragePt->FindBin( avrgPt ) );
}

Double_t Efficiencies::rpDeadMatVetoEffEastWestCorrelation(UInt_t branchConf) const{
  if(branchConf<Util::nBranchesConfigurations){
    return mhEastWestVetoEffCorrelation->GetBinContent( mhEastWestVetoEffCorrelation->GetXaxis()->FindBin( static_cast<double>(branchConf) ) );
  }
  else{ std::cerr << "ERROR in Efficiencies::rpDeadMatVetoEffEastWestCorrelation(UInt_t)" << std::endl; return 0.0; }
}

Double_t Efficiencies::rpTotalEffEastWestCorrelation(UInt_t branchConf) const{
  if(branchConf<Util::nBranchesConfigurations){
    return mhEastWestTotalEffCorrelation->GetBinContent( mhEastWestTotalEffCorrelation->GetXaxis()->FindBin( static_cast<double>(branchConf) ) );
  }
  else{ std::cerr << "ERROR in Efficiencies::rpTotalEffEastWestCorrelation(UInt_t)" << std::endl; return 0.0; }
}

Double_t Efficiencies::rpTotalEffEastWestCorrelation(UInt_t branchConf, Double_t tE, Double_t tW) const{
  if(branchConf<Util::nBranchesConfigurations){
    if( tE > mhEastWestTotalEffCorrelationVsMandelstamT[branchConf]->GetXaxis()->GetXmin() &&
        tE < mhEastWestTotalEffCorrelationVsMandelstamT[branchConf]->GetXaxis()->GetXmax() &&
        tW > mhEastWestTotalEffCorrelationVsMandelstamT[branchConf]->GetYaxis()->GetXmin() &&
        tW < mhEastWestTotalEffCorrelationVsMandelstamT[branchConf]->GetYaxis()->GetXmax()
    )
      return mhEastWestTotalEffCorrelationVsMandelstamT[branchConf]->Interpolate( tE, tW );
    else{
      int bin1 = (tE < mhEastWestTotalEffCorrelationVsMandelstamT[branchConf]->GetXaxis()->GetXmin()) ? 1 : ( (tE > mhEastWestTotalEffCorrelationVsMandelstamT[branchConf]->GetXaxis()->GetXmax()) ? mhEastWestTotalEffCorrelationVsMandelstamT[branchConf]->GetNbinsX() : mhEastWestTotalEffCorrelationVsMandelstamT[branchConf]->GetXaxis()->FindBin(tE) );
      int bin2 = (tW < mhEastWestTotalEffCorrelationVsMandelstamT[branchConf]->GetYaxis()->GetXmin()) ? 1 : ( (tW > mhEastWestTotalEffCorrelationVsMandelstamT[branchConf]->GetYaxis()->GetXmax()) ? mhEastWestTotalEffCorrelationVsMandelstamT[branchConf]->GetNbinsY() : mhEastWestTotalEffCorrelationVsMandelstamT[branchConf]->GetXaxis()->FindBin(tW) );
      return mhEastWestTotalEffCorrelationVsMandelstamT[branchConf]->GetBinContent( bin1, bin2 );
    }
  }
  else{ std::cerr << "ERROR in Efficiencies::rpTotalEffEastWestCorrelation(UInt_t, Double_t)" << std::endl; return 0.0; }
}

Double_t Efficiencies::rpTpcTotalEfficiencyCorrelation(Double_t mass) const{
  return mhRpTpcTotalEfficiencyCorrelationVsInvMass->Interpolate( mass );
}

Double_t Efficiencies::simultaneousBothSidesProtonLossProbability(UInt_t branchConf, Double_t zVtx) const{
  if(branchConf<Util::nBranchesConfigurations){
      return mhSimultaneousEastWestRpTrackLossProbability_zVertex->GetBinContent( mhSimultaneousEastWestRpTrackLossProbability_zVertex->GetXaxis()->FindBin( static_cast<double>(branchConf) ), mhSimultaneousEastWestRpTrackLossProbability_zVertex->GetYaxis()->FindBin( zVtx ) );
  }
  else{ std::cerr << "ERROR in Efficiencies::simultaneousBothSidesProtonLossProbability(UInt_t, Double_t)" << std::endl; return 0.0; }
}


void Efficiencies::multiplyTofEffPlateu(UInt_t sign, UInt_t pid, Double_t multFactor){
  if(sign<Util::nSigns && pid<Util::nDefinedParticles && multFactor>0){
      mTofMatchingEfficiency[sign][pid]->FixParameter(0, multFactor*mTofMatchingEffParameters[sign][pid][0]);
  }
  else{ std::cerr << "ERROR in Efficiencies::multiplyTofEffPlateu(UInt_t, UInt_t, Double_t)" << std::endl; }
}


void Efficiencies::shiftTofEffPlateu(UInt_t sign, UInt_t pid, Double_t shift){
  if(sign<Util::nSigns && pid<Util::nDefinedParticles ){
      mTofMatchingEfficiency[sign][pid]->FixParameter(0, shift + mTofMatchingEffParameters[sign][pid][0]);
  }
  else{ std::cerr << "ERROR in Efficiencies::shiftTofEffPlateu(UInt_t, UInt_t, Double_t)" << std::endl; }
}

void Efficiencies::shiftTpcEffPlateu(UInt_t sign, UInt_t pid, Double_t shift){
  if(sign<Util::nSigns && pid<Util::nDefinedParticles ){
      mTpcReconstructionEfficiency[sign][pid]->FixParameter(0, shift + mTpcRecoEffParameters[sign][pid][0]);
  }
  else{ std::cerr << "ERROR in Efficiencies::shiftTpcEffPlateu(UInt_t, UInt_t, Double_t)" << std::endl; }
}

void Efficiencies::multiplyTpcRecoEffPlateu(UInt_t sign, UInt_t pid, Double_t multFactor){
  if(sign<Util::nSigns && pid<Util::nDefinedParticles && multFactor>0){
      mTpcReconstructionEfficiency[sign][pid]->FixParameter(0, multFactor*mTpcRecoEffParameters[sign][pid][0]);
  }
  else{ std::cerr << "ERROR in Efficiencies::multiplyTpcRecoEffPlateu(UInt_t, UInt_t, Double_t)" << std::endl; }
}


void Efficiencies::shiftTofEffInflation(UInt_t sign, UInt_t pid, Double_t shift){
  if(sign<Util::nSigns && pid<Util::nDefinedParticles){
      mTofMatchingEfficiency[sign][pid]->FixParameter(1, shift + mTofMatchingEffParameters[sign][pid][1]);
  }
  else{ std::cerr << "ERROR in Efficiencies::shiftTofEffInflation(UInt_t, UInt_t, Double_t)" << std::endl; }
}

void Efficiencies::shiftTpcRecoEffInflation(UInt_t sign, UInt_t pid, Double_t shift){
  if(sign<Util::nSigns && pid<Util::nDefinedParticles){
      mTpcReconstructionEfficiency[sign][pid]->FixParameter(1, shift + mTpcRecoEffParameters[sign][pid][1]);
  }
  else{ std::cerr << "ERROR in Efficiencies::shiftTpcRecoEffInflation(UInt_t, UInt_t, Double_t)" << std::endl; }
}

void Efficiencies::resetTofEffToDefault(UInt_t sign, UInt_t pid){
  if(sign<Util::nSigns && pid<Util::nDefinedParticles){
    for(int par=0; par<3; ++par)
      mTofMatchingEfficiency[sign][pid]->FixParameter(par, mTofMatchingEffParameters[sign][pid][par]);
  }
  else{ std::cerr << "ERROR in Efficiencies::resetTofEffToDefault(UInt_t, UInt_t)" << std::endl; }
}

void Efficiencies::resetTpcRecoEffToDefault(UInt_t sign, UInt_t pid){
  if(sign<Util::nSigns && pid<Util::nDefinedParticles){
    for(int par=0; par<3; ++par)
      mTpcReconstructionEfficiency[sign][pid]->FixParameter(par, mTpcRecoEffParameters[sign][pid][par]);
  }
  else{ std::cerr << "ERROR in Efficiencies::resetTpcRecoEffToDefault(UInt_t, UInt_t)" << std::endl; }
}




void Efficiencies::shiftTofEff3DPlateu(UInt_t sign, UInt_t pid, Double_t shift){
  if(sign<Util::nSigns && pid<Util::nDefinedParticles){
    mTofEff3D_PlateuShift[sign][pid] = shift;
    mTofEff3D_isUnchanged[sign][pid] = false;
  }
  else{ std::cerr << "ERROR in Efficiencies::shiftTofEff3DPlateu(UInt_t, UInt_t, Double_t)" << std::endl; }
}

void Efficiencies::multiplyTofEff3DPlateu(UInt_t sign, UInt_t pid, Double_t multFactor){
  if(sign<Util::nSigns && pid<Util::nDefinedParticles && multFactor>0){
    mTofEff3D_MultiplicationFactor[sign][pid] = multFactor;
//     mTofEff3D_isUnchanged[sign][pid] = false;
  }
  else{ std::cerr << "ERROR in Efficiencies::multiplyTofEff3DPlateu(UInt_t, UInt_t, Double_t)" << std::endl; }
}

void Efficiencies::shiftTofEff3DInflation(UInt_t sign, UInt_t pid, Double_t shift){
  if(sign<Util::nSigns && pid<Util::nDefinedParticles){
      mTofEff3D_InflationShift[sign][pid] = shift;
//       mTofEff3D_isUnchanged[sign][pid] = false;
  }
  else{ std::cerr << "ERROR in Efficiencies::shiftTofEff3DInflation(UInt_t, UInt_t, Double_t)" << std::endl; }
}

void Efficiencies::shiftTpcRecoEff3DPlateu(UInt_t sign, UInt_t pid, Double_t shift){
  if(sign<Util::nSigns && pid<Util::nDefinedParticles){
    mTpcRecoEff3D_PlateuShift[sign][pid] = shift;
    mTpcRecoEff3D_isUnchanged[sign][pid] = false;
  }
  else{ std::cerr << "ERROR in Efficiencies::shiftTpcRecoEff3DPlateu(UInt_t, UInt_t, Double_t)" << std::endl; }
}

void Efficiencies::shiftTpcRecoEff3DInflation(UInt_t sign, UInt_t pid, Double_t shift){
  if(sign<Util::nSigns && pid<Util::nDefinedParticles){
    mTpcRecoEff3D_InflationShift[sign][pid] = shift;
//     mTpcRecoEff3D_isUnchanged[sign][pid] = false;
  }
  else{ std::cerr << "ERROR in Efficiencies::shiftTpcRecoEff3DInflation(UInt_t, UInt_t, Double_t)" << std::endl; }
}

void Efficiencies::multiplyTpcRecoEff3DPlateu(UInt_t sign, UInt_t pid, Double_t multFactor){
  if(sign<Util::nSigns && pid<Util::nDefinedParticles && multFactor>0){
    mTpcRecoEff3D_MultiplicationFactor[sign][pid] = multFactor;
//     mTpcRecoEff3D_isUnchanged[sign][pid] = false;
  }
  else{ std::cerr << "ERROR in Efficiencies::multiplyTpcRecoEff3DPlateu(UInt_t, UInt_t, Double_t)" << std::endl; }
}

void Efficiencies::resetTofEff3DToDefault(UInt_t sign, UInt_t pid){
  if(sign<Util::nSigns && pid<Util::nDefinedParticles){
    mTofEff3D_PlateuShift[sign][pid] = 0.0;
    mTofEff3D_MultiplicationFactor[sign][pid] = 1.0;
    mTofEff3D_InflationShift[sign][pid] = 0.0;
    mTofEff3D_isUnchanged[sign][pid] = true;
  }
  else{ std::cerr << "ERROR in Efficiencies::resetTofEff3DToDefault(UInt_t, UInt_t)" << std::endl; }
}

void Efficiencies::resetTpcRecoEff3DToDefault(UInt_t sign, UInt_t pid){
  if(sign<Util::nSigns && pid<Util::nDefinedParticles){
    mTpcRecoEff3D_PlateuShift[sign][pid] = 0.0;
    mTpcRecoEff3D_MultiplicationFactor[sign][pid] = 1.0;
    mTpcRecoEff3D_InflationShift[sign][pid] = 0.0;
    mTpcRecoEff3D_isUnchanged[sign][pid] = true;
  }
  else{ std::cerr << "ERROR in Efficiencies::resetTpcRecoEff3DToDefault(UInt_t, UInt_t)" << std::endl; }
}
