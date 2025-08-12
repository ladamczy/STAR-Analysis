#ifndef __CINT__
#include "TROOT.h"
#include "TSystem.h"
#include <iostream>
#include <bitset>
#include <string>    
#include <utility>
#include <sstream> 
#include <algorithm> 
#include <stdio.h> 
#include <stdlib.h> 
#include <vector> 
#include <fstream> 
#include <cmath>
#include <math.h> 
#include <cstdlib>
#include <bitset>
#include <map>
#include <set>
#include <vector>
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TRandom3.h"
#include "TDatime.h"
#include "StarRoot/TUnixTime.h"
#include "StChain.h"
#include "StMessMgr.h"
#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h" 
#include "StMuDSTMaker/COMMON/StMuBTofHit.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuMcTrack.h"
#include "StMuDSTMaker/COMMON/StMuMcVertex.h"
#include "StMuDSTMaker/COMMON/StMuPrimaryVertex.h"
#include "StBTofHeader.h"
#include "StThreeVectorF.hh"
#include "StPhysicalHelixD.hh"
#endif

using namespace std;

// Function implementations
double tofPathLength(const StThreeVector<double>* beginPoint,
                     const StThreeVector<double>* endPoint,
                     const double curvature) {
    if (!beginPoint || !endPoint) {
        return -999.0;
    }
    
    StThreeVector<double> diff = *endPoint - *beginPoint;
    double straightDistance = diff.mag();
    
    if (fabs(curvature) < 1e-6) {
        return straightDistance;
    }
    
    double radius = 1.0 / fabs(curvature);
    double chord = straightDistance;
    
    if (chord > 2.0 * radius) {
        return straightDistance;
    }
    
    double halfAngle = asin(chord / (2.0 * radius));
    double arcLength = 2.0 * radius * halfAngle;
    
    return arcLength;
}

// Additional overloads for compatibility
double tofPathLength(const StThreeVectorF* beginPoint,
                     const StThreeVectorF* endPoint,
                     const double curvature) {
    if (!beginPoint || !endPoint) return -999.0;
    StThreeVector<double> bp(beginPoint->x(), beginPoint->y(), beginPoint->z());
    StThreeVector<double> ep(endPoint->x(), endPoint->y(), endPoint->z());
    return tofPathLength(&bp, &ep, curvature);
}

double tofPathLength(const StThreeVectorD* beginPoint,
                     const StThreeVectorF* endPoint,
                     const double curvature) {
    if (!beginPoint || !endPoint) return -999.0;
    StThreeVector<double> bp(beginPoint->x(), beginPoint->y(), beginPoint->z());
    StThreeVector<double> ep(endPoint->x(), endPoint->y(), endPoint->z());
    return tofPathLength(&bp, &ep, curvature);
}

double tofPathLength(const StThreeVectorF* beginPoint,
                     const StThreeVectorD* endPoint,
                     const double curvature) {
    if (!beginPoint || !endPoint) return -999.0;
    StThreeVector<double> bp(beginPoint->x(), beginPoint->y(), beginPoint->z());
    StThreeVector<double> ep(endPoint->x(), endPoint->y(), endPoint->z());
    return tofPathLength(&bp, &ep, curvature);
}

double tofPathLength(const TVector3 beginPoint,
                     const TVector3 endPoint,
                     const double curvature) {
    StThreeVector<double> bp(beginPoint.X(), beginPoint.Y(), beginPoint.Z());
    StThreeVector<double> ep(endPoint.X(), endPoint.Y(), endPoint.Z());
    return tofPathLength(&bp, &ep, curvature);
}
// Calculate m2TOF directly from TOF information
double calculateM2TOF(const StMuBTofPidTraits &tofPid, const StThreeVectorF &momentum) {
    if(tofPid.matchFlag() == 0) return -999.0;  // No TOF match
    
    double tof = tofPid.timeOfFlight();      // Time of flight (ns)
    double pathLength = tofPid.pathLength(); // Path length (cm)
    double p = momentum.mag();               // Total momentum (GeV/c)
    
    if(tof <= 0 || pathLength <= 0 || p <= 0) return -999.0;
    
    // Constants
    const double c = 29.9792458;  // Speed of light (cm/ns)
    
    // Calculate beta = path/(c*time)
    double beta = pathLength / (c * tof);
    
    if(beta <= 0 || beta >= 1.0) return -999.0;
    
    // Calculate m2 = p2(1/β2 - 1)
    double betaSq = beta * beta;
    double m2 = p * p * (1.0/betaSq - 1.0);
    
    return m2;
}
bool passFiducialCut(double eta, double zvtx) {
    return (eta < -zvtx/250.0 + 0.9) && 
           (eta > -zvtx/250.0 - 0.9) && 
           (fabs(eta) < 0.9);
}

void UnifiedMuDstQA(TString InputFileList, Int_t nFiles = 1, Int_t nEvents = 0, 
                    TString OutputFile = "unified_analysis_output.root", Bool_t isJobMode = false) 
{
    cout << "==============================================================" << endl;
    cout << "     Unified MuDST QA Analysis (Fixed Version)             " << endl;
    cout << "==============================================================" << endl;
    cout << "Input: " << InputFileList << endl;
    cout << "Files: " << nFiles << ", Events: " << nEvents << endl;
    cout << "Output: " << OutputFile << endl;
    cout << "Mode: " << (isJobMode ? "Job Processing" : "Interactive") << endl;
    cout << "==============================================================" << endl;

    // Grid detection and library loading
    Bool_t isGridJob = (gSystem->Getenv("JOBID") != 0) || 
                       (gSystem->Getenv("_CONDOR_JOB_AD") != 0) ||
                       (gSystem->Getenv("SUMS_EXECUTION") != 0) ||
                       isJobMode;
    
    cout << "Grid job detected: " << (isGridJob ? "YES" : "NO") << endl;

    // Library loading
    cout << "Checking STAR library status..." << endl;
    TClass* chainClass = TClass::GetClass("StChain");
    Bool_t librariesLoaded = (chainClass != 0);
    cout << "STAR libraries already loaded: " << (librariesLoaded ? "YES" : "NO") << endl;
    
    if (!librariesLoaded) {
        cout << "Loading STAR libraries..." << endl;
        TString libMacro = "/afs/rhic.bnl.gov/star/packages/SL20c/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C";
        cout << "Loading macro: " << libMacro << endl;
        gROOT->LoadMacro(libMacro.Data());
        loadSharedLibraries();
        cout << "Libraries loaded successfully." << endl;
    }
    
    if (isGridJob) {
        cout << "Applying grid optimizations..." << endl;
        gROOT->SetBatch(kTRUE);
        gSystem->SetFPEMask(0);
    }
    
    cout << "Environment initialized." << endl;

    // Create output file
    TFile *tags_output = new TFile(OutputFile, "recreate");
    if (!tags_output || tags_output->IsZombie()) {
        cerr << "Error: Cannot create output file " << OutputFile << endl;
        return;
    }
    tags_output->cd();

    // Physics constants
    // Constants for analysis
    const double VzMax = 100.0; // Maximum Vz for fiducial cut
    const double PtMax = 3.0; // Maximum transverse momentum // we use it for efficiency ToF and TPC
    
    //pion parameters
    /*const double partMass = 0.13957018; // Pion mass
    const double partMassSq = partMass * partMass;
    const int gePidPlus = 8;  // pi+
    const int gePidMinus = 9; // pi-
    const double PtMin = 0.25; // Minimum transverse momentum
    cout << "Going to analyze Pions (#pi+ and #pi-) " << endl;*/
    
    // kaon parameters
    /*const double partMass = 0.493677; // Kaon mass
    const double partMassSq = partMass * partMass;
    const int gePidPlus = 11;  // K+
    const int gePidMinus = 12; // K-
    const double PtMin = 0.3; // Minimum transverse momentum
    const double m2tofMin = 0.15; // not used
    cout << "Going to analyze Kaons (K+ and K-) " << endl;*/

     // proton parameters 
     const double partMass = 0.938272; // Proton mass
     const double partMassSq = partMass * partMass;
     const int gePidPlus = 14;  // p+
     const int gePidMinus = 15; // p-
     const double PtMin = 0.4; // Minimum transverse momentum
     const double m2tofMin = 0.6; //not used
     cout << "Going to analyze Protons (p+ and p-) " << endl;

    // Random number generator
    TRandom3 gen(0);

    // Set default event count
    if (nEvents == 0) nEvents = 1000;

    // Analysis flow
    TString anaCutName[] = {"All", "MuDst", "McTrcks/vrtx", "VrtxSel", "2recoTrks"};
    enum ANA_CUTS { kAll = 1, kMuDst, kMcTracksAndVerticies, kVertexSelected, kTwoRecoTracks, nAnaCuts };
    TH1D *hAnaFlow = new TH1D("AnalysisFlow", "CutsFlow", nAnaCuts-1, 1, nAnaCuts);
    for(int tb = 1; tb < nAnaCuts; ++tb) 
        hAnaFlow->GetXaxis()->SetBinLabel(tb, anaCutName[tb-1]);
    
    // Track flow - FIXED labels
    TString anaTrackName[] = {"ALL", "MCMATCHED", "MCTRACK", "IDMATCH", "GEPID", "MCVERTEX", "PID", "QATRUTH", "NHITSFIT", "NHITSDEDX", "FLAG", "CHARGE", "TOF"};
    enum TRACK_CUTS { ALL = 1, MCMATCHED, MCTRACK, IDMATCH, GEPID, MCVERTEX, PID, QATRUTH, NHITSFIT, NHITSDEDX, FLAG, CHARGE, TOF, nTracksCuts};
    TH1D *hTrackFlow = new TH1D("TracksFlow", "Track Cuts Flow", nTracksCuts-1, 1, nTracksCuts);
    for(int tb = 1; tb < nTracksCuts; ++tb) 
        hTrackFlow->GetXaxis()->SetBinLabel(tb, anaTrackName[tb-1]);

    cout << "Booking histograms and trees..." << endl;
   
    // Monte Carlo histograms
    TH2F *hNMCVerteciesAndTracks = new TH2F("hNMCVerteciesAndTracks", "", 50, -0.5, 49.5, 50, -0.5, 49.5);
  
    // TOF histograms
    TH1F *hTray = new TH1F("hTray", "", 100, 0, 1000);
    TH1F *hModule = new TH1F("hModule", "", 100, 0, 1000);
    TH1F *hCell = new TH1F("hCell", "", 100, 0, 1000);
    TH1F *hLeadingEdgeTime = new TH1F("hLeadingEdgeTime", "", 100, 0, 60000);
    TH1F *hTrailingEdgeTime = new TH1F("hTrailingEdgeTime", "", 100, 0, 60000);
    TH1F *hTofTime = new TH1F("hTofTime", "", 100, 0, 100);
    TH1F *hTofPathLength = new TH1F("hTofPathLength", "", 100, -2000, 1000);
    TH1F *hNTOFHits = new TH1F("hNTOFHits", "", 100, 0, 100);
    TH1F *hPhiMc = new TH1F("hPhiMc", "Phi of Mc tracks", 20, -TMath::Pi(), TMath::Pi());
    TH1F *hPtMc = new TH1F("hPtMc", "Pt of Mc tracks", 60, 0, 3.0);
    TH1F *hEtaMc = new TH1F("hEtaMc", "Eta of Mc tracks", 20, -1.0, 1.0);
    TH1F *hVzMc = new TH1F("hVzMc", "Vz of MC tracks", 20, -100, 100);
    
    // Reconstructed track histograms
    TH2F *hNVerteciesAndTracks = new TH2F("hNVerteciesAndTracks", "", 100, -0.5, 99.5, 250, -0.5, 249.5);
    TH1F *hNGlobalTracks = new TH1F("hNGlobalTracks", "", 4000, -0.5, 3999.5);
    
    TH1F *hPhi = new TH1F("hPhi", "Phi of matched RC tracks", 20, -TMath::Pi(), TMath::Pi());
    TH1F *hPt = new TH1F("hPt", "Pt of matched RC tracks", 60, 0, 3.0);
    TH1F *hVz = new TH1F("hVz", "Vz of matched RC tracks", 20, -100, 100);
    TH1F *hEta = new TH1F("hEta", "Eta of matched RC tracks", 20, -1.0, 1.0);
    
    TH1F *hPhi_reco = new TH1F("hPhi_reco", "True Phi of matched RC tracks", 20, -TMath::Pi(), TMath::Pi());
    TH1F *hPt_reco = new TH1F("hPt_reco", "True Pt of matched RC tracks", 60, 0, 3.0);
    TH1F *hEta_reco = new TH1F("hEta_reco", "True Eta of matched RC tracks", 20, -1.0, 1.0);
    TH1F *hVz_reco = new TH1F("hVz_reco", "True Vz of matched RC tracks", 20, -100, 100);

    //TPC efficiency numerators (MC truth of TPC reconstructed tracks) - positive
    TH1F *hSelPtMc_P = new TH1F("hSelPtMc_P", "True Pt of selected MC tracks", 60, 0, 3.0);//100, 0.1, 3.1
    TH1F *hSelEtaMc_P = new TH1F("hSelEtaMc_P", "True Eta of selected MC tracks", 20, -1.0, 1.0); //20, -1.0, 1.0
    TH1F *hSelPhiMc_P = new TH1F("hSelPhiMc_P", "True Phi of selected MC tracks", 20, -TMath::Pi(), TMath::Pi());
    TH1F *hSelVzMc_P = new TH1F("hSelVzMc_P", "True Vz of selected MC tracks", 20, -100, 100);
    // TPC efficiency numerators (MC truth of TPC reconstructed tracks) - negative
    TH1F *hSelPtMc_N = new TH1F("hSelPtMc_N", "True Pt of selected MC tracks", 60, 0, 3.0);
    TH1F *hSelEtaMc_N = new TH1F("hSelEtaMc_N", "True Eta of selected MC tracks", 20, -1.0, 1.0);
    TH1F *hSelPhiMc_N = new TH1F("hSelPhiMc_N", "True Phi of selected MC tracks", 20, -TMath::Pi(), TMath::Pi());
    TH1F *hSelVzMc_N = new TH1F("hSelVzMc_N", "True Vz of selected MC tracks", 20, -100, 100); 

    // TPC efficiency numerators and (MC truth of TPC reconstructed tracks) - positive
    TH1F *hSelPtTrue_P = new TH1F("hSelPtTrue_P", "True Pt of selected matched RC tracks", 60, 0, 3.0);
    TH1F *hSelEtaTrue_P = new TH1F("hSelEtaTrue_P", "True Eta of selected matched RC tracks", 20, -1.0, 1.0);
    TH1F *hSelPhiTrue_P = new TH1F("hSelPhiTrue_P", "True Phi of selected matched RC tracks", 20, -TMath::Pi(), TMath::Pi());
    TH1F *hSelVzTrue_P = new TH1F("hSelVzTrue_P", "True Vz of selected matched RC tracks", 20, -100, 100);
    // TPC efficiency numerators and (MC truth of TPC reconstructed tracks) - negative
    TH1F *hSelPtTrue_N = new TH1F("hSelPtTrue_N", "True Pt of selected matched RC tracks", 60, 0, 3.0);
    TH1F *hSelEtaTrue_N = new TH1F("hSelEtaTrue_N", "True Eta of selected matched RC tracks", 20, -1.0, 1.0);
    TH1F *hSelPhiTrue_N = new TH1F("hSelPhiTrue_N", "True Phi of selected matched RC tracks", 20, -TMath::Pi(), TMath::Pi());
    TH1F *hSelVzTrue_N = new TH1F("hSelVzTrue_N", "True Vz of selected matched RC tracks", 20, -100, 100);

    // TOF efficiency numerators (MC truth of TPC+TOF matched tracks) - positive
    TH1F *hPtTpcTofMatched_P = new TH1F("hPtTpcTofMatched_P", "True Pt reconstructed TPC+TOF matched", 60, 0, 3.0);
    TH1F *hEtaTpcTofMatched_P = new TH1F("hEtaTpcTofMatched_P", "True Eta reconstructed TPC+TOF matched", 20, -1.0, 1.0);
    TH1F *hPhiTpcTofMatched_P = new TH1F("hPhiTpcTofMatched_P", "True Phi reconstructed TPC+TOF matched", 20, -TMath::Pi(), TMath::Pi());
    TH1F *hVzTpcTofMatched_P = new TH1F("hVzTpcTofMatched_P", "True Vz reconstructed TPC+TOF matched", 20, -100, 100);

    // TOF efficiency numerators (MC truth of TPC+TOF matched tracks) - negative
    TH1F *hPtTpcTofMatched_N = new TH1F("hPtTpcTofMatched_N", "True Pt reconstructed TPC+TOF matched", 60, 0, 3.0);
    TH1F *hEtaTpcTofMatched_N = new TH1F("hEtaTpcTofMatched_N", "True Eta reconstructed TPC+TOF matched", 20, -1.0, 1.0);
    TH1F *hPhiTpcTofMatched_N = new TH1F("hPhiTpcTofMatched_N", "True Phi reconstructed TPC+TOF matched", 20, -TMath::Pi(), TMath::Pi());
    TH1F *hVzTpcTofMatched_N = new TH1F("hVzTpcTofMatched_N", "True Vz reconstructed TPC+TOF matched", 20, -100, 100);
    
    // ToF efficiency denominator and (MC truth of TPC reconstructed tracks) - positive
    TH1F *hSelPtReco_P = new TH1F("hSelPtReco_P", "Reco Pt of selected matched RC tracks", 60, 0, 3.0);
    TH1F *hSelEtaReco_P = new TH1F("hSelEtaReco_P", "Reco Eta of selected matched RC tracks", 20, -1.0, 1.0);
    TH1F *hSelPhiReco_P = new TH1F("hSelPhiReco_P", "Reco Phi of selected matched RC tracks", 20, -TMath::Pi(), TMath::Pi());
    TH1F *hSelVzReco_P = new TH1F("hSelVzReco_P", "Reco Vz of selected matched RC tracks", 20, -100, 100);
    // ToF efficiency denominator and (MC truth of TPC reconstructed tracks) - negative
    TH1F *hSelPtReco_N = new TH1F("hSelPtReco_N", "Reco Pt of selected matched RC tracks", 60, 0, 3.0);
    TH1F *hSelEtaReco_N = new TH1F("hSelEtaReco_N", "Reco Eta of selected matched RC tracks", 20, -1.0, 1.0);
    TH1F *hSelPhiReco_N = new TH1F("hSelPhiReco_N", "Reco Phi of selected matched RC tracks", 20, -TMath::Pi(), TMath::Pi());
    TH1F *hSelVzReco_N = new TH1F("hSelVzReco_N", "Reco Vz of selected matched RC tracks", 20, -100, 100);
    
    // =============================================================================
    // TRUE 3D HISTOGRAMS FOR PROPER TEfficiency CONSTRUCTION
    // =============================================================================
    
    // Define consistent binning scheme for all 3D histograms
    // Using same binning as your 1D histograms but in 3D arrays
    Double_t ptBins[] = {0.20, 0.25, 0.30, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                         1.0, 1.1, 1.2, 1.3, 1.4, 1.6, 2.0, 3.0};
    Int_t nPtBins = sizeof(ptBins)/sizeof(ptBins[0]) - 1;
    
    Double_t etaBins[] = {-1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,
                          0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    Int_t nEtaBins = sizeof(etaBins)/sizeof(etaBins[0]) - 1;
    
    Double_t vzBins[] = {-100, -90, -80, -70, -60, -50, -40, -30, -20, -10,
                         0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    Int_t nVzBins = sizeof(vzBins)/sizeof(vzBins[0]) - 1;
    
    // Numerator: MC truth of TPC reconstructed tracks
    TH3F *h3D_TPC_num_P = new TH3F("h3D_TPC_num_P", "TPC 3D Numerator (+)", nPtBins, ptBins, nEtaBins, etaBins, nVzBins, vzBins);
    h3D_TPC_num_P->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h3D_TPC_num_P->GetYaxis()->SetTitle("#eta");
    h3D_TPC_num_P->GetZaxis()->SetTitle("V_{z} (cm)");
    
    // Denominator: MC truth of all generated tracks in acceptance
    TH3F *h3D_TPC_den_P = new TH3F("h3D_TPC_den_P", "TPC 3D Denominator (+)", nPtBins, ptBins, nEtaBins, etaBins, nVzBins, vzBins);
    h3D_TPC_den_P->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h3D_TPC_den_P->GetYaxis()->SetTitle("#eta");
    h3D_TPC_den_P->GetZaxis()->SetTitle("V_{z} (cm)");
    
    // TPC 3D Efficiency Histograms - NEGATIVE PARTICLES
    TH3F *h3D_TPC_num_N = new TH3F("h3D_TPC_num_N", "TPC 3D Numerator (-)", nPtBins, ptBins, nEtaBins, etaBins, nVzBins, vzBins);
    h3D_TPC_num_N->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h3D_TPC_num_N->GetYaxis()->SetTitle("#eta");
    h3D_TPC_num_N->GetZaxis()->SetTitle("V_{z} (cm)");
    
    TH3F *h3D_TPC_den_N = new TH3F("h3D_TPC_den_N", "TPC 3D Denominator (-)", nPtBins, ptBins, nEtaBins, etaBins, nVzBins, vzBins);
    h3D_TPC_den_N->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h3D_TPC_den_N->GetYaxis()->SetTitle("#eta");
    h3D_TPC_den_N->GetZaxis()->SetTitle("V_{z} (cm)");
    
    // TOF 3D Efficiency Histograms - POSITIVE PARTICLES
    // Numerator: TPC+TOF matched tracks (using reconstructed variables)
    TH3F *h3D_TOF_num_P = new TH3F("h3D_TOF_num_P", "TOF 3D Numerator (+)",nPtBins, ptBins, nEtaBins, etaBins, nVzBins, vzBins);
    h3D_TOF_num_P->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h3D_TOF_num_P->GetYaxis()->SetTitle("#eta");
    h3D_TOF_num_P->GetZaxis()->SetTitle("V_{z} (cm)");
    
    // Denominator: TPC reconstructed tracks (using reconstructed variables)
    TH3F *h3D_TOF_den_P = new TH3F("h3D_TOF_den_P", "TOF 3D Denominator (+)", nPtBins, ptBins, nEtaBins, etaBins, nVzBins, vzBins);
    h3D_TOF_den_P->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h3D_TOF_den_P->GetYaxis()->SetTitle("#eta");
    h3D_TOF_den_P->GetZaxis()->SetTitle("V_{z} (cm)");
    
    // TOF 3D Efficiency Histograms - NEGATIVE PARTICLES
    TH3F *h3D_TOF_num_N = new TH3F("h3D_TOF_num_N", "TOF 3D Numerator (-)", nPtBins, ptBins, nEtaBins, etaBins, nVzBins, vzBins);
    h3D_TOF_num_N->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h3D_TOF_num_N->GetYaxis()->SetTitle("#eta");
    h3D_TOF_num_N->GetZaxis()->SetTitle("V_{z} (cm)");
    
    TH3F *h3D_TOF_den_N = new TH3F("h3D_TOF_den_N", "TOF 3D Denominator (-)", nPtBins, ptBins, nEtaBins, etaBins, nVzBins, vzBins);
    h3D_TOF_den_N->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h3D_TOF_den_N->GetYaxis()->SetTitle("#eta");
    h3D_TOF_den_N->GetZaxis()->SetTitle("V_{z} (cm)");
    
    cout << "✅ 3D efficiency histograms booked successfully:" << endl;

    // CORRECTED error counters - only count actual errors
    UInt_t mcFromWrongVertex = 0, mcArrayInconsistency = 0, idTruthMissmatched = 0;
    UInt_t wrongGID = 0, noMCinfo = 0, noMuDst = 0;
    
    // Track counters for embedded data
    UInt_t backgroundTracks = 0;  // idTruth = 0 (normal!)
    UInt_t mcMatchedTracks = 0;   // idTruth > 0 (MC matched)
    UInt_t badIdTruth = 0;        // idTruth < 0 or > NoMcTracks (actual errors)

    // MC counters
    UInt_t totalMcParticles = 0;
    UInt_t selectedMcParticles = 0;

    // Tree variables
    const unsigned int nSigns = 2;
    double mTofTimeOrig[nSigns], mTofLengthOrig[nSigns], mTofTime[nSigns], mTofLength[nSigns];

    TTree* mTree = new TTree("recTree", "recTree");
    for (int i = 0; i < nSigns; ++i) {
        mTree->Branch(Form("tofTimeInNs%i", i), &mTofTime[i]);
        mTree->Branch(Form("tofLengthInCm%i", i), &mTofLength[i]);
        mTree->Branch(Form("tofLengthOrig%i", i), &mTofLengthOrig[i]);
        mTree->Branch(Form("tofTimeOrig%i", i), &mTofTimeOrig[i]);
    }
    
    cout << "Histograms and trees booked successfully." << endl;

    // Setup STAR analysis chain
    StChain* chain = new StChain;
    StMuDstMaker* muDstMaker = new StMuDstMaker(0, 0, "", InputFileList, "MuDst", nFiles);
    
    // Configure data branches
    muDstMaker->SetStatus("*", 0);
    muDstMaker->SetStatus("MuEvent", 1);
    muDstMaker->SetStatus("PrimaryVertices", 1);
    muDstMaker->SetStatus("PrimaryTracks", 1);
    muDstMaker->SetStatus("GlobalTracks", 1);
    muDstMaker->SetStatus("MCAll", 1);
    muDstMaker->SetStatus("BTofHeader", 1);
    muDstMaker->SetStatus("BTofHit", 1);
    muDstMaker->SetStatus("Event", 1);
    muDstMaker->SetStatus("McEvent", 1);
    muDstMaker->SetStatus("TofHit", 1);
    muDstMaker->SetStatus("TofData", 1);
    muDstMaker->SetStatus("TofRawData", 1);
    muDstMaker->SetStatus("BTofRawHit", 1);
    muDstMaker->SetDebug(1);
    
    // Initialize the chain
    Int_t iInit = chain->Init();
    if (iInit) {
        chain->FatalErr(iInit, "on init");
        return;
    }

    // Process events
    cout << "Processing " << nEvents << " events..." << endl;
    
    Int_t istat = 0, i = 1;
    Int_t progressStep = (nEvents / 20 > 1) ? nEvents / 20 : 1;
    
    while (i <= nEvents && istat != 2) {
        if (i % progressStep == 0 && !isJobMode) {
            cout << "Processing event " << i << " (" << (100*i/nEvents) << "%)" << endl;
        }
        
        chain->Clear();
        istat = chain->Make(i);
        
        if (istat == 2) {
            cout << "Last event processed. Status = " << istat << endl;
            break;
        }
        if (istat == 3) {
            cout << "Error event processed. Status = " << istat << endl;
            i++;
            continue;
        }
        
        i++;
        
        if (istat != kStOK) continue;

        hAnaFlow->Fill(kAll);
        
        // Get MuDST data
        StMuDst* mMuDst = muDstMaker->muDst();
        if (!mMuDst) {
            noMuDst++;
            continue;
        }
        
        hAnaFlow->Fill(kMuDst);
        
        // Get MC arrays
        TClonesArray *MuMcVertices = mMuDst->mcArray(0);
        TClonesArray *MuMcTracks = mMuDst->mcArray(1);
        if (!MuMcVertices || !MuMcTracks) {
            noMCinfo++;
            continue;
        }
        
        Int_t NoMuMcVertices = MuMcVertices->GetEntriesFast();
        Int_t NoMuMcTracks = MuMcTracks->GetEntriesFast();
        
        if (!NoMuMcVertices || !NoMuMcTracks) {
            noMCinfo++;
            continue;
        }
        
        // FIXED: Fill this bin when MC info is available
        hAnaFlow->Fill(kMcTracksAndVerticies);
        
        Int_t nMc = 0;
        
        // Process MC tracks
        for (Int_t itrk = 0; itrk < NoMuMcTracks; itrk++) {
            StMuMcTrack *mcTrack = (StMuMcTrack*)MuMcTracks->UncheckedAt(itrk);
            if (!mcTrack) continue;
            
            const int Gid = mcTrack->GePid();
            
            // Check if it's a particle first
            if (Gid == gePidPlus || Gid == gePidMinus) {
                totalMcParticles++;
                
                // Select only from primary vertex
                Int_t IdVx = mcTrack->IdVx();
                if (IdVx != 1) {
                    mcFromWrongVertex++;
                    continue;
                }
                
                // Get MC vertex for this track
                StMuMcVertex* mcVertex = (StMuMcVertex*)MuMcVertices->UncheckedAt(IdVx-1);
                if (!mcVertex) continue;
                
                // Fill all MC histograms with TRUE variables
                double truePt = mcTrack->Pxyz().perp();
                double truePhi = mcTrack->Pxyz().phi();
                double trueEta = mcTrack->Pxyz().pseudoRapidity();
                double trueVz = mcVertex->XyzV().z();
                
                hPtMc->Fill(truePt);
                hPhiMc->Fill(truePhi);
                hEtaMc->Fill(trueEta);
                hVzMc->Fill(trueVz);  // 
                //these Mc were for TPC
                bool isPositive = (Gid == gePidPlus);   // +
                bool isNegative = (Gid == gePidMinus);  // -

                //TPC efficiency denominator
                if (fabs(trueVz) < VzMax && passFiducialCut(trueEta, trueVz)) {  // For pt distribution: cut on nominal eta and vz
                    if (isPositive) hSelPtMc_P->Fill(truePt);
                    if (isNegative) hSelPtMc_N->Fill(truePt);
                }
                if (truePt > PtMin && truePt < PtMax && fabs(trueVz) < VzMax) { //&& passFiducialCut(trueEta, trueVz)) {  // For eta distribution: cut on nominal pt and vz
                    if (isPositive) hSelEtaMc_P->Fill(trueEta);
                    if (isNegative) hSelEtaMc_N->Fill(trueEta);
                }
                if (truePt > PtMin && truePt < PtMax && passFiducialCut(trueEta, trueVz)) {  // For vz distribution: cut on nominal pt and eta
                    if (isPositive) hSelVzMc_P->Fill(trueVz);
                    if (isNegative) hSelVzMc_N->Fill(trueVz);
                }
                if (truePt > PtMin && truePt < PtMax && fabs(trueVz) < VzMax && passFiducialCut(trueEta, trueVz)) {  // For phi distribution: cut on nominal pt, eta and vz
                    if (isPositive) hSelPhiMc_P->Fill(truePhi);
                    if (isNegative) hSelPhiMc_N->Fill(truePhi);
                }
                // 3D TPC Deno
                if (isPositive) h3D_TPC_den_P->Fill(truePt, trueEta, trueVz);
                if (isNegative) h3D_TPC_den_N->Fill(truePt, trueEta, trueVz);
              
                selectedMcParticles++;
                
            } else {
                wrongGID++;
            }
        }
        
        hNMCVerteciesAndTracks->Fill(NoMuMcVertices, nMc);
        
        // Vertex selection
        //removing vetex dependancy
        /*int const originalVertexId = mMuDst->currentVertexIndex();
        StMuPrimaryVertex* selectedVertex = 0;
        
        mMuDst->setVertexIndex(0);
        selectedVertex = mMuDst->primaryVertex();
        
        if (!selectedVertex) { 
            mMuDst->setVertexIndex(originalVertexId);
            continue;
        }*/
        
        hAnaFlow->Fill(kVertexSelected);
        //Global Track information
        hNGlobalTracks->Fill(mMuDst->numberOfGlobalTracks());

        // Get track information
        TObjArray* tracks = mMuDst->globalTracks();

        int nTracks = mMuDst->numberOfPrimaryTracks();
        int nVertecies = mMuDst->numberOfPrimaryVertices();
        hNVerteciesAndTracks->Fill(nVertecies, nTracks);
        
        // TOF information
        Int_t nTofHits = mMuDst->numberOfBTofHit();
        hNTOFHits->Fill(nTofHits);
        for (Int_t iTof = 0; iTof < nTofHits; iTof++) {
            StMuBTofHit *tofHit = mMuDst->btofHit(iTof);
            if (tofHit) {
                hTray->Fill(tofHit->tray());
                hModule->Fill(tofHit->module());
                hCell->Fill(tofHit->cell());
                hLeadingEdgeTime->Fill(tofHit->leadingEdgeTime());
                hTrailingEdgeTime->Fill(tofHit->trailingEdgeTime());
            }
        }
        TObjArrayIter GetTracks(tracks);
        StMuTrack* ptrack;
        int nRecoTracks = 0;
        
        // Initialize TOF arrays
        for (int j = 0; j < nSigns; j++) {
            mTofTime[j] = -999;
            mTofLength[j] = -999;
            mTofTimeOrig[j] = -999;
            mTofLengthOrig[j] = -999;
        }
        
        // FIXED: Process reconstructed tracks for embedded data
        while ((ptrack = (StMuTrack*)GetTracks.Next())) {
            hTrackFlow->Fill(ALL);

            // Get ID Truth value
            Int_t idTruthValue = ptrack->idTruth();
            
            // FIXED: For embedded data, any idTruth outside valid range is background
            if (idTruthValue <= 0 || idTruthValue > NoMuMcTracks) {
                backgroundTracks++;
                continue;  // Skip background tracks - NOT an error!
            }
            
            // We have a valid MC-matched track
            mcMatchedTracks++;
            hTrackFlow->Fill(MCMATCHED);
            
            // Access MC track safely
            StMuMcTrack *mcTrack = (StMuMcTrack*)MuMcTracks->UncheckedAt(idTruthValue - 1);
            if (!mcTrack) {
                mcArrayInconsistency++;
                continue;
            }
            hTrackFlow->Fill(MCTRACK);
            
            // Verify ID matching
            if (mcTrack->Id() != idTruthValue) {
                idTruthMissmatched++;
                continue;
            }
            hTrackFlow->Fill(IDMATCH);
            
            // Check particle ID
            if(mcTrack->GePid() != gePidPlus && mcTrack->GePid() != gePidMinus) 
                continue;
            hTrackFlow->Fill(GEPID);
            
            // Check vertex
            if(mcTrack->IdVx() != 1) {
                continue;
            }
            hTrackFlow->Fill(MCVERTEX);
            
            // PID selection
            //if(ptrack->nSigmaPion() > 3.0)   //did not use for pions and got agreement
            //    continue;

            //if(fabs(ptrack->nSigmaKaon()) > 3.0)
            //    continue;

            //if(fabs(ptrack->nSigmaProton()) > 3.0)
            //    continue;

            hTrackFlow->Fill(PID);
            
            // Quality cuts
            if(ptrack->qaTruth() < 95.) 
                continue;
            hTrackFlow->Fill(QATRUTH);
            
            if(ptrack->nHitsFit() < 20)
                continue;
            hTrackFlow->Fill(NHITSFIT);
            
            if(ptrack->nHitsDedx() < 15)
                continue;
            hTrackFlow->Fill(NHITSDEDX);
            
            if(ptrack->flag() <= 0)
                continue;
            hTrackFlow->Fill(FLAG);
            
            if(abs(ptrack->charge()) != 1) 
                continue;
            hTrackFlow->Fill(CHARGE);
            
            // TOF matching
            const StMuBTofPidTraits &tofPid = ptrack->btofPidTraits();
            //if(tofPid.matchFlag() == 0)  //for TPC efficiency we commented so that no bias
            //    continue;  
            hTrackFlow->Fill(TOF);

            StThreeVectorF p = ptrack->p();

            // Calculate m2TOF ------------------------------------// did not work
            //did not use for pions
            /*double m2tof = calculateM2TOF(tofPid, p);
            
            // Apply kaon m2TOF cut
            if(m2tof < m2tofMin) {  //
                continue;
            }*/
            
            // TOF calculations
            int charge = ptrack->charge() > 0 ? 0 : 1;

            hTofTime->Fill(tofPid.timeOfFlight());
            hTofPathLength->Fill(tofPid.pathLength());
          
            double tofTime = (tofPid.pathLength()/299792458.0) * sqrt(1 + ((partMassSq)/p.mag2())) * pow(10.0, 7.0);
            tofTime += gen.Gaus(0.0, 0.1);
            
            mTofTimeOrig[charge] = tofPid.timeOfFlight();
            mTofLengthOrig[charge] = tofPid.pathLength();
            mTofTime[charge] = tofTime;          
            
            // Path length calculation
            // Old (vertex-dependent):
            // StPhysicalHelixD p1Helix = ptrack->helix();
            // StThreeVector<double> vtx = selectedVertex->position(); 
            
            // New (vertex-independent):
            StPhysicalHelixD p1Helix = ptrack->helix();  // Use track origin directly

            double x,y,z;
            x = (double)ptrack->helix().origin().x();
            y = (double)ptrack->helix().origin().y();
            z = (double)ptrack->helix().origin().z();
            StThreeVector<double> bp(x,y,z);

            x = (double)tofPid.position().x();
            y = (double)tofPid.position().y();
            z = (double)tofPid.position().z();
            StThreeVector<double> ep(x,y,z);
          
            // Old (vertex-dependent):
            // p1Helix.moveOrigin(p1Helix.pathLength(vtx));
            // mTofLength[charge] = tofPathLength(&bp, &ep, ptrack->helix().curvature());
            
            // New (vertex-independent):
            mTofLength[charge] = tofPathLength(
                &bp, 
                &ep, 
                ptrack->helix().curvature()
            );

            // Get reconstructed variables
            StThreeVectorF pReco = ptrack->p();
            double recoEta = pReco.pseudoRapidity();
            double recoPhi = pReco.phi();
            double recoPt = pReco.perp();
              
            // Get TRUE variables from MC track
            double truePt = mcTrack->Pxyz().perp();
            double truePhi = mcTrack->Pxyz().phi();
            double trueEta = mcTrack->Pxyz().pseudoRapidity();
            StMuMcVertex* mcVertex = (StMuMcVertex*)MuMcVertices->UncheckedAt(mcTrack->IdVx()-1);
            double trueVz = mcVertex ? mcVertex->XyzV().z() : -999.0;

            // Add track-MC association
            double deltaEtaPhi = sqrt(pow(trueEta - recoEta, 2) + pow(truePhi - recoPhi, 2));
            if(deltaEtaPhi > 0.15)
                continue;
            
            // Fill QA histograms with RECONSTRUCTED variables
            hPhi->Fill(pReco.phi());
            hPt->Fill(pReco.perp());
            hEta->Fill(pReco.pseudoRapidity());
            //hVz->Fill(selectedVertex->position().z());

            // Fill TRUE variable histograms for all matched tracks
            hPhi_reco->Fill(truePhi);
            hPt_reco->Fill(truePt);
            hEta_reco->Fill(trueEta);
            hVz_reco->Fill(trueVz);

            const int Gid = mcTrack->GePid();
            bool isPositive = (Gid == gePidPlus);   // +
            bool isNegative = (Gid == gePidMinus);  // -
          
            // TPC efficiency numerators (MC truth of TPC reconstructed tracks) 
            if (fabs(trueVz) < VzMax && passFiducialCut(trueEta, trueVz)) {  //For pT efficiency: cut on eta and vz
                if (isPositive) hSelPtTrue_P->Fill(truePt);
                if (isNegative) hSelPtTrue_N->Fill(truePt);
            }
            if (truePt > PtMin && truePt < PtMax && fabs(trueVz) < VzMax) {   // For eta efficiency: cut on pt and vz
                if (isPositive) hSelEtaTrue_P->Fill(trueEta);
                if (isNegative) hSelEtaTrue_N->Fill(trueEta);
            }
            if (truePt > PtMin && truePt < PtMax && passFiducialCut(trueEta, trueVz)) {   // For vz efficiency: cut on pt and eta
                if (isPositive) hSelVzTrue_P->Fill(trueVz);
                if (isNegative) hSelVzTrue_N->Fill(trueVz);
            }
            if (truePt > PtMin && truePt < PtMax && fabs(trueVz) < VzMax && passFiducialCut(trueEta, trueVz)) {   // For phi efficiency: cut on pt, eta and vz
                if (isPositive) hSelPhiTrue_P->Fill(truePhi);
                if (isNegative) hSelPhiTrue_N->Fill(truePhi);
            }
            //3D TPC num
            if (isPositive) h3D_TPC_num_P->Fill(truePt, trueEta, trueVz);
            if (isNegative) h3D_TPC_num_N->Fill(truePt, trueEta, trueVz);

            // ToF efficiency denominator 
            if (fabs(trueVz) < VzMax && passFiducialCut(recoEta, trueVz)) {  // For pT efficiency: cut on eta and vz
                if (isPositive) hSelPtReco_P->Fill(recoPt);
                if (isNegative) hSelPtReco_N->Fill(recoPt);
            } 
            if (recoPt > PtMin && recoPt < PtMax && fabs(trueVz) < VzMax) {   // For eta efficiency: cut on pt and vz
                if (isPositive) hSelEtaReco_P->Fill(recoEta);
                if (isNegative) hSelEtaReco_N->Fill(recoEta);
            } 
            if (recoPt > PtMin && recoPt < PtMax && passFiducialCut(recoEta, trueVz)) {   // For vz efficiency: cut on pt and eta
                if (isPositive) hSelVzReco_P->Fill(trueVz);
                if (isNegative) hSelVzReco_N->Fill(trueVz);
            } 
            if (recoPt > PtMin && recoPt < PtMax && fabs(trueVz) < VzMax && passFiducialCut(recoEta, trueVz)) {  // For phi efficiency: cut on pt, eta and vz
                if (isPositive) hSelPhiReco_P->Fill(recoPhi);
                if (isNegative) hSelPhiReco_N->Fill(recoPhi);
            }
            // 3D TOF deno
            if (isPositive) h3D_TOF_den_P->Fill(recoPt, recoEta, trueVz);
            if (isNegative) h3D_TOF_den_N->Fill(recoPt, recoEta, trueVz);

            // ToF efficiency numerators
            if(tofPid.matchFlag() != 0) {
                if (fabs(trueVz) < VzMax && passFiducialCut(recoEta, trueVz)) {  //// For pt efficiency: cut on RECO eta and TRUE vz, fill RECO pt
                    if (isPositive) hPtTpcTofMatched_P->Fill(recoPt);  //
                    if (isNegative) hPtTpcTofMatched_N->Fill(recoPt);
                }  
                if (recoPt > PtMin && recoPt < PtMax && fabs(trueVz) < VzMax) {   ////For eta efficiency: cut on RECO pt and TRUE vz, fill RECO eta
                    if (isPositive) hEtaTpcTofMatched_P->Fill(recoEta);  //  
                    if (isNegative) hEtaTpcTofMatched_N->Fill(recoEta);
                } 
                if (recoPt > PtMin && recoPt < PtMax && passFiducialCut(recoEta, trueVz)) {  // For vz efficiency: cut on RECO pt and RECO eta, fill TRUE vz
                    if (isPositive) hVzTpcTofMatched_P->Fill(trueVz);  // 
                    if (isNegative) hVzTpcTofMatched_N->Fill(trueVz);
                }  
                if (recoPt > PtMin && recoPt < PtMax && fabs(trueVz) < VzMax && passFiducialCut(recoEta, trueVz)) {  //For phi efficiency: cut on RECO pt, RECO eta and TRUE vz, fill RECO phi
                    if (isPositive) hPhiTpcTofMatched_P->Fill(recoPhi);  //
                    if (isNegative) hPhiTpcTofMatched_N->Fill(recoPhi);
                }
                //3D TOF num
                if (isPositive) h3D_TOF_num_P->Fill(recoPt, recoEta, trueVz);
                if (isNegative) h3D_TOF_num_N->Fill(recoPt, recoEta, trueVz);

            }
            
            nRecoTracks++;
        }
        
        // Require exactly 2 reconstructed tracks (one positive, one negative)
        if(nRecoTracks != nSigns)
            continue;

        hAnaFlow->Fill(kTwoRecoTracks);
        mTree->Fill();
    }
    
    cout << "Event processing completed." << endl;
    
    /*// Calculate efficiency 
    cout << "Calculating efficiency histograms using TRUE variables..." << endl;
    
    TH1F *hEffPt_P = new TH1F("hEffPt_P", "Efficiency in pt bins", 60, 0, 3.0);
    TH1F *hEffVz_P = new TH1F("hEffVz_P", "Efficiency in Vz bins", 20, -100, 100);
    TH1F *hEffEta_P = new TH1F("hEffEta_P", "Efficiency in eta bins", 20, -1.0, 1.0);
    TH1F *hEffPhi_P = new TH1F("hEffPhi_P", "Efficiency in phi bins", 20, -TMath::Pi(), TMath::Pi());
   
    TH1F *hEffPt_N = new TH1F("hEffPt_N", "Efficiency in pt bins", 60, 0, 3.0);
    TH1F *hEffVz_N = new TH1F("hEffVz_N", "Efficiency in Vz bins", 20, -100, 100);
    TH1F *hEffEta_N = new TH1F("hEffEta_N", "Efficiency in eta bins", 20, -1.0, 1.0);
    TH1F *hEffPhi_N = new TH1F("hEffPhi_N", "Efficiency in phi bins", 20, -TMath::Pi(), TMath::Pi());

    cout << "\n=== HISTOGRAM ENTRIES CHECK ===" << endl;
    cout << "MC Truth histograms:" << endl;
    cout << "  hSelPtMc_P: " << hSelPtMc_P->GetEntries() << endl;
    cout << "  hSelEtaMc_P: " << hSelEtaMc_P->GetEntries() << endl;
    cout << "  hSelPhiMc_P: " << hSelPhiMc_P->GetEntries() << endl;
    cout << "  hSelVzMc_P: " << hSelVzMc_P->GetEntries() << endl;
    
    cout << "Reconstructed True histograms:" << endl;
    cout << "  hSelPtTrue_P: " << hSelPtTrue_P->GetEntries() << endl;
    cout << "  hSelEtaTrue_P: " << hSelEtaTrue_P->GetEntries() << endl;
    cout << "  hSelPhiTrue_P: " << hSelPhiTrue_P->GetEntries() << endl;
    cout << "  hSelVzTrue_P: " << hSelVzTrue_P->GetEntries() << endl;
    cout << "================================" << endl;
    
    // Calculate efficiency: Reco(true vars) / MC(true vars)
    if (hSelPtMc_P->GetEntries() > 0) {
        hEffPt_P->Divide(hSelPtTrue_P, hSelPtMc_P, 1, 1, "B");  // ✅ Correct
        cout << "Pt efficiency calculated using true variables" << endl;
    } else {
        cout << "WARNING: No MC truth entries for Pt efficiency!" << endl;
    }
    
    if (hSelVzMc_P->GetEntries() > 0) {  // ✅ FIXED: Use hSelVzMc_P, not hVzMc
        hEffVz_P->Divide(hSelVzTrue_P, hSelVzMc_P, 1, 1, "B");  // ✅ Correct
        // Cap efficiency at 1.0
        for (int bin = 1; bin <= hEffVz_P->GetNbinsX(); bin++) {
            if (hEffVz_P->GetBinContent(bin) > 1.0) {
                hEffVz_P->SetBinContent(bin, 1.0);
                hEffVz_P->SetBinError(bin, 0.0);
            }
        }
        cout << "Vz efficiency calculated using true variables" << endl;
    } else {
        cout << "WARNING: No MC truth entries for Vz efficiency!" << endl;
    }
    
    if (hSelEtaMc_P->GetEntries() > 0) {
        hEffEta_P->Divide(hSelEtaTrue_P, hSelEtaMc_P, 1, 1, "B");  // ✅ Correct
        cout << "Eta efficiency calculated using true variables" << endl;
    } else {
        cout << "WARNING: No MC truth entries for Eta efficiency!" << endl;
    }
    
    if (hSelPhiMc_P->GetEntries() > 0) {
        hEffPhi_P->Divide(hSelPhiTrue_P, hSelPhiMc_P, 1, 1, "B");  // ✅ Correct
        // Cap efficiency at 1.0
        for (int bin = 1; bin <= hEffPhi_P->GetNbinsX(); bin++) {
            if (hEffPhi_P->GetBinContent(bin) > 1.0) {
                hEffPhi_P->SetBinContent(bin, 1.0);
                hEffPhi_P->SetBinError(bin, 0.0);
            }
        }
        cout << "Phi efficiency calculated using true variables" << endl;
    } else {
        cout << "WARNING: No MC truth entries for Phi efficiency!" << endl;
    }
    hEffPt_N->Divide(hSelPtTrue_N, hSelPtMc_N, 1, 1, "B");
    hEffVz_N->Divide(hSelVzTrue_N, hSelVzMc_N, 1, 1, "B");
    hEffEta_N->Divide(hSelEtaTrue_N, hSelEtaMc_N, 1, 1, "B");
    hEffPhi_N->Divide(hSelPhiTrue_N, hSelPhiMc_N, 1, 1, "B");

    cout << "\n=== EFFICIENCY VERIFICATION ===" << endl;
    for (int bin = 1; bin <= 10; bin++) {  // Check first 10 bins
        double mcEntries = hSelPtMc_P->GetBinContent(bin);
        double recoEntries = hSelPtTrue_P->GetBinContent(bin);
        double efficiency = hEffPt_P->GetBinContent(bin);
        double expectedEff = (mcEntries > 0) ? recoEntries / mcEntries : 0;
        
        cout << "Bin " << bin << " (pT=" << hEffPt_P->GetBinCenter(bin) << "): ";
        cout << "MC=" << mcEntries << ", Reco=" << recoEntries;
        cout << ", Eff=" << efficiency << ", Expected=" << expectedEff << endl;
        
        if (fabs(efficiency - expectedEff) > 0.001 && mcEntries > 0) {
            cout << "  *** ERROR: Efficiency mismatch! ***" << endl;
        }
    }
    cout << "===============================" << endl;*/
    // Finalize chain
    if (i > 1) {
        cout << "Finishing chain..." << endl;
        chain->Finish();
        delete chain;
        chain = 0;
    }
    
    // Write output
    cout << "Writing output file..." << endl;
    tags_output->cd();
    
    // Write all histograms
    hAnaFlow->Write();
    hTrackFlow->Write();

    hNMCVerteciesAndTracks->Write();
    hNVerteciesAndTracks->Write();
    hNGlobalTracks->Write();

    hPhiMc->Write();
    hPtMc->Write();
    hEtaMc->Write();
    hVzMc->Write();

    hPhi->Write();
    hPt->Write();
    hEta->Write();
    //hVz->Write();

    hSelPtMc_P->Write();
    hSelEtaMc_P->Write();
    hSelPhiMc_P->Write();
    hSelVzMc_P->Write();    
    hSelPtTrue_P->Write();
    hSelEtaTrue_P->Write();
    hSelPhiTrue_P->Write();
    hSelVzTrue_P->Write();

    hSelPtReco_P->Write();
    hSelEtaReco_P->Write();
    hSelPhiReco_P->Write();
    hSelVzReco_P->Write();
    hSelPtReco_N->Write();  
    hSelEtaReco_N->Write();
    hSelPhiReco_N->Write();
    hSelVzReco_N->Write();
    
    hSelPtMc_N->Write(); 
    hSelEtaMc_N->Write();
    hSelPhiMc_N->Write();
    hSelVzMc_N->Write();
    hSelPtTrue_N->Write();  
    hSelEtaTrue_N->Write();
    hSelPhiTrue_N->Write();
    hSelVzTrue_N->Write();

    hPtTpcTofMatched_P->Write();
    hEtaTpcTofMatched_P->Write();
    hPhiTpcTofMatched_P->Write();
    hVzTpcTofMatched_P->Write();
    hPtTpcTofMatched_N->Write();
    hEtaTpcTofMatched_N->Write();
    hPhiTpcTofMatched_N->Write();
    hVzTpcTofMatched_N->Write();

    hPhi_reco->Write();
    hPt_reco->Write();
    hEta_reco->Write();
    hVz_reco->Write();

    // Write 3D efficiency histograms
    h3D_TPC_num_P->Write();
    h3D_TPC_den_P->Write();
    h3D_TPC_num_N->Write();
    h3D_TPC_den_N->Write();
    h3D_TOF_num_P->Write();
    h3D_TOF_den_P->Write();
    h3D_TOF_num_N->Write();
    h3D_TOF_den_N->Write();
    
    cout << "✅ 3D efficiency histograms written to file" << endl;

    // Write TOF histograms
    hTray->Write();
    hModule->Write();
    hCell->Write();
    hLeadingEdgeTime->Write();
    hTrailingEdgeTime->Write();
    hTofTime->Write();
    hTofPathLength->Write();
    hNTOFHits->Write();
    
    // Write tree
    mTree->Write();
    
    // Close output file
    tags_output->Close();
    cout << "Output file written successfully." << endl;
    
    // Calculate actual errors (not including background tracks)
    UInt_t actualErrors = mcArrayInconsistency + idTruthMissmatched + badIdTruth;
    
    cout << "\n==============================================================" << endl;
    cout << "                  ANALYSIS SUMMARY" << endl;
    cout << "==============================================================" << endl;
    
    cout << "\nEmbedded Data Track Summary:" << endl;
    cout << "  Background tracks (idTruth≤0): " << backgroundTracks << " (NORMAL)" << endl;
    cout << "  MC-matched tracks (idTruth>0): " << mcMatchedTracks << endl;
    cout << "  Bad idTruth values: " << badIdTruth << endl;
    
    cout << "\nMC Particle Summary:" << endl;
    cout << "  Total MC particles: " << totalMcParticles << endl;
    cout << "  Selected MC particles (|eta|<0.9): " << selectedMcParticles << endl;
    
    if (actualErrors > 0) {
        cout << "\nActual Error Summary:" << endl;
        cout << "  Missing muDst: " << noMuDst << endl;
        cout << "  Missing MC info: " << noMCinfo << endl;
        cout << "  Wrong GID of MC particle: " << wrongGID << endl;
        cout << "  MC originated from wrong vertex: " << mcFromWrongVertex << endl;
        cout << "  MC array inconsistency: " << mcArrayInconsistency << endl;
        cout << "  Bad ID truth values: " << badIdTruth << endl;
        cout << "  ID truth mismatched: " << idTruthMissmatched << endl;
        cout << "  Total ACTUAL errors: " << actualErrors << endl;
    } else {
        cout << "\n✅ No actual errors found! Analysis is working correctly." << endl;
    }
    
    cout << "\nOutput Information:" << endl;
    cout << "  Output file: " << OutputFile << endl;
    cout << "  Tree entries: " << mTree->GetEntries() << endl;
    cout << "  Analysis flow: All bins should now be filled!" << endl;
    cout << "  Track flow: MC-matched tracks should progress through cuts!" << endl;

    cout << "\n==============================================================" << endl;
}
