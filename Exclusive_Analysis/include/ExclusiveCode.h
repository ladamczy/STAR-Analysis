#ifndef EXCLUSIVE_CODE_H
#define EXCLUSIVE_CODE_H

#include <iostream>
#include <string>    
#include <utility>
#include <sstream> 
#include <algorithm> 
#include <stdio.h> 
#include <stdlib.h> 
#include <vector> 
#include <fstream> 
#include <cmath> 
#include <cstdlib>
#include <sys/stat.h>
#include <iterator>
#include <ostream>
#include <TRandom.h>
#include <iomanip>
#include <stdexcept>
#include <limits>
#include "TROOT.h"
#include "TSystem.h"
#include "TThread.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include "TProfile.h"
#include <TH2.h> 
#include <TF1.h> 
#include <TF2.h> 
#include <THStack.h> 
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TStyle.h> 
#include <TGraph.h> 
#include <TGraph2D.h> 
#include <TGraphErrors.h> 
#include <TCanvas.h> 
#include <TLegend.h> 
#include <TGaxis.h> 
#include <TString.h> 
#include <TColor.h> 
#include <TLine.h> 
#include <TExec.h> 
#include <TFitResultPtr.h> 
#include <TFitResult.h> 
#include <TLatex.h> 
#include <TMath.h>
#include <TLorentzVector.h>
#include <ROOT/TThreadedObject.hxx>
#include <TTreeReader.h>
#include <ROOT/TTreeProcessorMT.hxx>
#include "StRPEvent.h"
#include "StUPCRpsTrack.h"
#include "StUPCRpsTrackPoint.h"
#include "StUPCEvent.h"
#include "StUPCTrack.h"
#include "StUPCBemcCluster.h"
#include "StUPCVertex.h"
#include "StUPCTofHit.h"
#include "StUPCV0.h"
#include "BeamPosition.h"
#include <StRPEvent.h>




namespace ExclusiveK0K0{

    const double PION_PT = 0.2;
    const double PION_ETA = 0.9;
    const double N_FIT = 20;
    const float MASS_PION = 0.13957061;
    const double MASS_KAON =  497.611/1000.0;
    const double BEAM_ENERGY = 254.867;
    const double MASS_PROTON = 938.272/1000.0;
    const double MASS_LAMBDA = 1.115;

    const double PT_MISS= 0.15; 
    const double DCA_BEAMLINE = 2.5;
    const double DCA_DAUGHTERS = 2.5;
    const double DCA_BEAMLINE_23 = 1.5;
    const double DCA_DAUGHTERS_23 = 1.5;

    const double COS_PNT_ANG = 0.925;

    const double COS_PNT_ANG_2 = 0.95;
    const double DECAY_LENGTH = 3.0;

}

void CheckN1(vector <bool *> vCuts,bool & condition, bool & conditionN1 );

void FindProtons(bool isMC, StRPEvent *rpEvt, StUPCEvent *upcEvt, TLorentzVector & proton1, TLorentzVector & proton2);
void SeparateTracks(StUPCEvent *upcEvt, vector <StUPCTrack const*> &tracksWithTofHit, vector <StUPCTrack const*> &tracksWithoutTofHit, bool isMC,  TH1D* HistNumWithTofTrakcs,TH1D* HistNumWithoutTofTrakcs);

bool ValidNumberOfTofTracks(vector <int> vNums, vector <StUPCTrack const*>  tracksWithTofHit , vector <StUPCTrack const*>  tracksWithoutTofHit ,vector<TH1D *>&HistPtPionWithTof,vector<TH1D *>&HistPtPionWithoutTof, vector<TH1D *>& HistEtaPionWithTof, vector<TH1D *> &HistEtaPionWithoutTof, vector<TH1D *>& HistNfitPionWithTof,vector<TH1D *>  &HistNfitPionWithoutTof  );
bool AreTofTracksGood(vector <StUPCTrack const*>  &tracksWithTofHit);
bool FindTracks( vector <StUPCTrack const*> &tracksWithTofHit, vector <StUPCTrack const*> &tracksWithoutTofHit, vector <StUPCTrack const*> &tracksWithoutTofHitGoodChargeThreeTof);


bool CheckNumberOfClusters(  StUPCEvent *upcEvt,      vector <StUPCTrack const*> &tracksWithTofHit, int &  totalCluster);
void GetBeamPar(StUPCEvent *upcEvt, double *beamPar, bool isMC);
bool FindProtonsAndPions(  StUPCEvent *upcEvt, vector <StUPCTrack const*> &tracksWithTofHit, vector <StUPCTrack const*> &goodTracksWithoutTofHit,vector <StUPCTrack const*> &  v1, vector <StUPCTrack const*> & v2, double * beamPar);
bool FindPions(  StUPCEvent *upcEvt, vector <StUPCTrack const*> &tracksWithTofHit, vector <StUPCTrack const*> &tracksWithoutTofHit, vector <StUPCTrack const*> &goodTracksWithoutTofHit,vector <StUPCTrack const*> &PosNegPionLeadingKaon, vector <StUPCTrack const*> & vPosNegPionSubLeadingKaon, double * beamPar);
bool CheckPtMiss(StUPCV0 &leadingKaon, StUPCV0 &subLeadingKaon, TLorentzVector proton1, TLorentzVector proton2, double & pTmiss);


#endif