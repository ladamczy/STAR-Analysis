// -------------------------------------------------------------------========
// EventAnalyzer.cxx
// 
// Minimal code to demonstrate RP data quality and SC 7 correlation issues
// for K_S K_S exclusive production analysis
//
// Purpose: Show that dataset lacks exclusive K_S K_S events due to
//          kinematic correlation failures (SC 7.3)
// -------------------------------------------------------------------========

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "StRPEvent.h"
#include "StUPCRpsTrack.h"
#include "StUPCEvent.h"
#include "StUPCTrack.h"
#include "StUPCV0.h"
#include "StUPCVertex.h" 

using namespace std;

// -------------------------------------------------------------------========
// CONSTANTS
// -------------------------------------------------------------------========
const float MASS_PION = 0.13957061;  // GeV/c^2
const float BEAM_MOMENTUM = 255.0;   // GeV/c for 510 GeV collision
const float SQRT_S = 510.0;          // GeV

// Selection cuts (matching thesis)
const double CUT_VTX_Z = 80.0;       // cm
const double CUT_ETA = 0.9;
const double CUT_PT = 0.2;           // GeV/c
const int CUT_NHITS_FIT = 20;
const int CUT_NHITS_DEDX = 15;

// V0 topology cuts (SC 5)
const double CUT_DCA_BEAMLINE = 2.5; // cm
const double CUT_DCA_DAUGHTERS = 2.5;// cm
const double CUT_DECAY_LENGTH = 3.0; // cm
const double CUT_COS_POINTING = 0.925;

// K_S mass window
const double MASS_KAON_LOW = 0.44;   // GeV/c^2
const double MASS_KAON_HIGH = 0.54;  // GeV/c^2

// SC 7 cuts (exclusive kinematics) - LOOSENED for demonstration
const double CUT_XI_MIN = 0.0001;
const double CUT_XI_MAX = 0.5;
const double CUT_CORRELATION1 = 0.5;  // |m_KK/sqrt(s*xi_E*xi_W) - 1| < 0.5
const double CUT_CORRELATION2 = 2.0;  // |y_KK - 0.5*ln(xi_E/xi_W)| < 2.0
const double CUT_CORRELATION3 = 0.5;  // |m_KK*exp(±y_KK)/sqrt(s*xi) - 1| < 0.5

// -------------------------------------------------------------------========
// MAIN ANALYSIS
// -------------------------------------------------------------------========
int main(int argc, char* argv[]) {
    
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " <input_file_list> <output_file>" << endl;
        return 1;
    }
    
    cout << "\n============================================" << endl;
    cout << "   EVENT ANALYZER - RP DIAGNOSTIC TOOL" << endl;
    cout << "============================================\n" << endl;
    
    // -----------------------------------------------------------------------
    // Open input files
    // -----------------------------------------------------------------------
    ifstream inputFileList(argv[1]);
    if (!inputFileList) {
        cerr << "ERROR: Cannot open input file list: " << argv[1] << endl;
        return 1;
    }
    
    TChain* chain = new TChain("mUPCTree");
    string inputFileName;
    while (getline(inputFileList, inputFileName)) {
        chain->Add(inputFileName.c_str());
    }
    inputFileList.close();
    
    cout << "Total entries in chain: " << chain->GetEntries() << endl;
    
    // -----------------------------------------------------------------------
    // Set up branches
    // -----------------------------------------------------------------------
    static StUPCEvent* upcEvt = nullptr;
    static StRPEvent* rpEvt = nullptr;
    
    chain->SetBranchAddress("mUPCEvent", &upcEvt);
    chain->GetEntry(0);
    
    bool isMC = (upcEvt->getRunNumber() == 1);
    if (!isMC) {
        chain->SetBranchAddress("correctedRpEvent", &rpEvt);
    }
    
    cout << "Running on: " << (isMC ? "MC" : "Real Data") << endl;
    cout << "\n============================================\n" << endl;
    
    // -----------------------------------------------------------------------
    // Create output file and histograms
    // -----------------------------------------------------------------------
    TFile* outFile = TFile::Open(argv[2], "RECREATE");
    
    // RP quality histograms
    TH1D* hRPProtonMomentum = new TH1D("hRPProtonMomentum", "RP Proton Momentum;p (GeV/c);Counts", 100, 200, 260);
    TH1D* hRPXi = new TH1D("hRPXi", "RP #xi distribution;#xi;Counts", 100, -0.1, 0.2);
    TH2D* hRPXiEvsXiW = new TH2D("hRPXiEvsXiW", "RP #xi_{E} vs #xi_{W};#xi_{E};#xi_{W}", 100, -0.1, 0.2, 100, -0.1, 0.2);
    
    // K_S K_S mass and kinematics
    TH1D* hMassKK = new TH1D("hMassKK", "K_{S}K_{S} Invariant Mass (Before SC 7);M_{KK} (GeV/c^{2});Counts", 100, 0.8, 2.0);
    TH1D* hMassKK_AfterSC7 = new TH1D("hMassKK_AfterSC7", "K_{S}K_{S} Invariant Mass (After SC 7);M_{KK} (GeV/c^{2});Counts", 100, 0.8, 2.0);
    
    // SC 7 correlation distributions
    TH1D* hCorrelation1 = new TH1D("hCorrelation1", "SC 7.1: m_{KK}/#sqrt{s#xi_{E}#xi_{W}};Correlation Value;Counts", 100, 0, 50);
    TH1D* hCorrelation2 = new TH1D("hCorrelation2", "SC 7.2: y_{KK} - 0.5ln(#xi_{E}/#xi_{W});Correlation Value;Counts", 100, -5, 5);
    TH1D* hCorrelation3a = new TH1D("hCorrelation3a", "SC 7.3a: m_{KK}exp(y)/#sqrt{s#xi_{W}};Correlation Value;Counts", 100, 0, 2);
    TH1D* hCorrelation3b = new TH1D("hCorrelation3b", "SC 7.3b: m_{KK}exp(-y)/#sqrt{s#xi_{E}};Correlation Value;Counts", 100, 0, 2);
    
    // 2D correlation plots
    TH2D* hXiSumVsMKK = new TH2D("hXiSumVsMKK", "#xi_{E}+#xi_{W} vs M_{KK};M_{KK} (GeV/c^{2});#xi_{E}+#xi_{W}", 50, 0.8, 2.0, 100, 0, 0.3);
    TH2D* hCorr1VsCorr3a = new TH2D("hCorr1VsCorr3a", "Correlation 1 vs 3a;Corr1;Corr3a", 50, 0, 5, 50, 0, 2);
    
    // SC 8 & 9 diagnostic histograms
    TH1D* hAvgZ = new TH1D("hAvgZ", "Average Decay Z Position;#LT z #GT (cm);Counts", 100, -200, 200);
    TH1D* hDeltaZ = new TH1D("hDeltaZ", "Z Distance Between Vertices;|#Delta z| (cm);Counts", 100, 0, 100);
    TH1D* hCosAngle1 = new TH1D("hCosAngle1", "cos(#theta) Leading K_{S};cos(#theta);Counts", 100, -1.1, 1.1);
    TH1D* hCosAngle2 = new TH1D("hCosAngle2", "cos(#theta) Subleading K_{S};cos(#theta);Counts", 100, -1.1, 1.1);

    // Cut flow histogram
    TH1D* hCutFlow = new TH1D("hCutFlow", "Cut Flow;Cut Stage;Events", 18, -0.5, 17.5);
    hCutFlow->GetXaxis()->SetBinLabel(1, "All Events");
    hCutFlow->GetXaxis()->SetBinLabel(2, "N_{vtx}=1");
    hCutFlow->GetXaxis()->SetBinLabel(3, "|z_{vtx}|<80");
    hCutFlow->GetXaxis()->SetBinLabel(4, "N_{TOF}=4");
    hCutFlow->GetXaxis()->SetBinLabel(5, "#Sigma q=0");
    hCutFlow->GetXaxis()->SetBinLabel(6, "|#eta|<0.9");
    hCutFlow->GetXaxis()->SetBinLabel(7, "p_{T}>0.2");
    hCutFlow->GetXaxis()->SetBinLabel(8, "N_{fit}#geq20");
    hCutFlow->GetXaxis()->SetBinLabel(9, "N_{dE/dx}#geq15");
    hCutFlow->GetXaxis()->SetBinLabel(10, "V0 Topology");
    hCutFlow->GetXaxis()->SetBinLabel(11, "K_{S} Mass");
    hCutFlow->GetXaxis()->SetBinLabel(12, "#xi Range");
    hCutFlow->GetXaxis()->SetBinLabel(13, "SC 7.1");
    hCutFlow->GetXaxis()->SetBinLabel(14, "SC 7.2");
    hCutFlow->GetXaxis()->SetBinLabel(15, "SC 7.3");
    hCutFlow->GetXaxis()->SetBinLabel(16, "SC 8.1");  
    hCutFlow->GetXaxis()->SetBinLabel(17, "SC 8.2");  
    hCutFlow->GetXaxis()->SetBinLabel(18, "SC 9");    
    
    // -----------------------------------------------------------------------
    // Event counters
    // -----------------------------------------------------------------------
    int nTotal = 0;
    int nWithRP = 0;
    int nRP2Tracks = 0;
    int nKSKS = 0;
    int nXiRange = 0;
    int nPassSC71 = 0;
    int nPassSC72 = 0;
    int nPassSC73 = 0;

    int nPassSC81 = 0;  
    int nPassSC82 = 0;  
    int nPassSC9 = 0;    
    
    int nNegativeXi = 0;
    int nLargeXi = 0;
    int nReasonableXi = 0;
    
    // -----------------------------------------------------------------------
    // EVENT LOOP
    // -----------------------------------------------------------------------
    Long64_t nEntries = chain->GetEntries();
    cout << "Starting event loop over " << nEntries << " events..." << endl;
    
    for (Long64_t iEvent = 0; iEvent < nEntries; ++iEvent) {
        
        if (iEvent % 100000 == 0) {
            cout << "Processing event " << iEvent << "/" << nEntries << endl;
        }
        
        chain->GetEntry(iEvent);
        
        hCutFlow->Fill(0); // All events
        nTotal++;
        
        // -------------------------------------------------------------------
        // Basic event selection
        // -------------------------------------------------------------------
        
        // Single vertex
        if (upcEvt->getNumberOfVertices() != 1) continue;
        hCutFlow->Fill(1);
        
        // Vertex position
        if (fabs(upcEvt->getVertex(0)->getPosZ()) > CUT_VTX_Z) continue;
        hCutFlow->Fill(2);
        
        // Four TOF-matched tracks
        vector<const StUPCTrack*> tofTracks;
        for (int iTrk = 0; iTrk < upcEvt->getNumberOfTracks(); ++iTrk) {
            if (upcEvt->getTrack(iTrk)->getFlag(StUPCTrack::kTof)) {
                tofTracks.push_back(upcEvt->getTrack(iTrk));
            }
        }
        if (tofTracks.size() != 4) continue;
        hCutFlow->Fill(3);
        
        // Total charge = 0
        int totalCharge = 0;
        for (auto trk : tofTracks) totalCharge += trk->getCharge();
        if (totalCharge != 0) continue;
        hCutFlow->Fill(4);
        
        // Eta cut
        bool etaOK = true;
        for (auto trk : tofTracks) {
            if (fabs(trk->getEta()) > CUT_ETA) {
                etaOK = false;
                break;
            }
        }
        if (!etaOK) continue;
        hCutFlow->Fill(5);
        
        // pT cut
        bool ptOK = true;
        for (auto trk : tofTracks) {
            if (trk->getPt() < CUT_PT) {
                ptOK = false;
                break;
            }
        }
        if (!ptOK) continue;
        hCutFlow->Fill(6);
        
        // NhitsFit cut
        bool nFitOK = true;
        for (auto trk : tofTracks) {
            if (trk->getNhitsFit() < CUT_NHITS_FIT) {
                nFitOK = false;
                break;
            }
        }
        if (!nFitOK) continue;
        hCutFlow->Fill(7);
        
        // NhitsDEdx cut
        bool nDedxOK = true;
        for (auto trk : tofTracks) {
            if (trk->getNhitsDEdx() < CUT_NHITS_DEDX) {
                nDedxOK = false;
                break;
            }
        }
        if (!nDedxOK) continue;
        hCutFlow->Fill(8);
        
        // -------------------------------------------------------------------
        // Separate positive and negative tracks
        // -------------------------------------------------------------------
        vector<const StUPCTrack*> posTracks, negTracks;
        for (auto trk : tofTracks) {
            if (trk->getCharge() > 0) posTracks.push_back(trk);
            else negTracks.push_back(trk);
        }
        
        if (posTracks.size() != 2 || negTracks.size() != 2) continue;
        
        // -------------------------------------------------------------------
        // Reconstruct V0 candidates
        // -------------------------------------------------------------------
        TVector3 beamVec(0, 0, 0);
        double beamLine[4] = {0, 0, 0, 0};
        
        StUPCV0 v0_00(posTracks[0], negTracks[0], MASS_PION, MASS_PION, 1, 1, 
                      beamVec, beamLine, upcEvt->getMagneticField(), isMC);
        StUPCV0 v0_01(posTracks[0], negTracks[1], MASS_PION, MASS_PION, 1, 1, 
                      beamVec, beamLine, upcEvt->getMagneticField(), isMC);
        StUPCV0 v0_10(posTracks[1], negTracks[0], MASS_PION, MASS_PION, 1, 1, 
                      beamVec, beamLine, upcEvt->getMagneticField(), isMC);
        StUPCV0 v0_11(posTracks[1], negTracks[1], MASS_PION, MASS_PION, 1, 1, 
                      beamVec, beamLine, upcEvt->getMagneticField(), isMC);
        
        // Choose best pairing (minimize distance from K_S mass)
        double dist1 = sqrt(pow(v0_00.m() - 0.4976, 2) + pow(v0_11.m() - 0.4976, 2));
        double dist2 = sqrt(pow(v0_01.m() - 0.4976, 2) + pow(v0_10.m() - 0.4976, 2));
        
        StUPCV0* ks1 = nullptr;
        StUPCV0* ks2 = nullptr;
        
        if (dist1 < dist2) {
            ks1 = &v0_00;
            ks2 = &v0_11;
        } else {
            ks1 = &v0_01;
            ks2 = &v0_10;
        }
        
        // Ensure ks1 has higher pT (leading)
        if (ks2->pt() > ks1->pt()) {
            StUPCV0* temp = ks1;
            ks1 = ks2;
            ks2 = temp;
        }
        
        // -------------------------------------------------------------------
        // V0 topology cuts (SC 5)
        // -------------------------------------------------------------------
        bool topologyOK = true;
        
        // DCA to beamline
        if (ks1->DCABeamLine() > CUT_DCA_BEAMLINE) topologyOK = false;
        if (ks2->DCABeamLine() > CUT_DCA_BEAMLINE) topologyOK = false;
        
        // DCA between daughters
        if (ks1->dcaDaughters() > CUT_DCA_DAUGHTERS) topologyOK = false;
        if (ks2->dcaDaughters() > CUT_DCA_DAUGHTERS) topologyOK = false;
        
        // Decay length OR pointing angle
        bool ks1OK = (ks1->decayLength() >= CUT_DECAY_LENGTH) || 
                     (ks1->pointingAngle() >= CUT_COS_POINTING);
        bool ks2OK = (ks2->decayLength() >= CUT_DECAY_LENGTH) || 
                     (ks2->pointingAngle() >= CUT_COS_POINTING);
        
        if (!ks1OK || !ks2OK) topologyOK = false;
        
        if (!topologyOK) continue;
        hCutFlow->Fill(9);
        
        // -------------------------------------------------------------------
        // K_S mass window
        // -------------------------------------------------------------------
        if (ks1->m() < MASS_KAON_LOW || ks1->m() > MASS_KAON_HIGH) continue;
        if (ks2->m() < MASS_KAON_LOW || ks2->m() > MASS_KAON_HIGH) continue;
        hCutFlow->Fill(10);
        
        nKSKS++;
        
        // -------------------------------------------------------------------
        // K_S K_S system kinematics
        // -------------------------------------------------------------------
        TLorentzVector lv_ks1 = ks1->lorentzVector();
        TLorentzVector lv_ks2 = ks2->lorentzVector();
        TLorentzVector lv_kk = lv_ks1 + lv_ks2;
        
        double mKK = lv_kk.M();
        double yKK = lv_kk.Rapidity();
        
        hMassKK->Fill(mKK);
        
        // -------------------------------------------------------------------
        // Read Roman Pot data
        // -------------------------------------------------------------------
        if (isMC || !rpEvt) continue;
        
        nWithRP++;
        
        if (rpEvt->getNumberOfTracks() != 2) continue;
        nRP2Tracks++;
        
        // Get RP proton kinematics
        StUPCRpsTrack* trk_E = rpEvt->getTrack(0);
        StUPCRpsTrack* trk_W = rpEvt->getTrack(1);
        
        trk_E->setEvent(rpEvt);
        trk_W->setEvent(rpEvt);
        
        TVector3 p_E = trk_E->pVec();
        TVector3 p_W = trk_W->pVec();
        
        double xiE = trk_E->xi(BEAM_MOMENTUM);
        double xiW = trk_W->xi(BEAM_MOMENTUM);
        
        // Fill RP quality histograms
        hRPProtonMomentum->Fill(p_E.Mag());
        hRPProtonMomentum->Fill(p_W.Mag());
        hRPXi->Fill(xiE);
        hRPXi->Fill(xiW);
        hRPXiEvsXiW->Fill(xiE, xiW);
        
        // Track xi quality
        if (xiE < 0) nNegativeXi++;
        else if (xiE > 0.1) nLargeXi++;
        else if (xiE >= 0.001 && xiE <= 0.1) nReasonableXi++;
        
        if (xiW < 0) nNegativeXi++;
        else if (xiW > 0.1) nLargeXi++;
        else if (xiW >= 0.001 && xiW <= 0.1) nReasonableXi++;
        
        // -------------------------------------------------------------------
        // SC 7: Kinematic correlations for exclusive production
        // -------------------------------------------------------------------
        
        // SC 7.0: Xi range check
        if (xiE < CUT_XI_MIN || xiE > CUT_XI_MAX) continue;
        if (xiW < CUT_XI_MIN || xiW > CUT_XI_MAX) continue;
        hCutFlow->Fill(11);
        nXiRange++;
        
        // Calculate correlations
        double correlation1 = mKK / sqrt(SQRT_S * xiE * xiW);
        double correlation2 = yKK - 0.5 * log(xiE / xiW);
        double correlation3a = mKK * exp(yKK) / sqrt(SQRT_S * xiW);
        double correlation3b = mKK * exp(-yKK) / sqrt(SQRT_S * xiE);
        
        // Fill correlation histograms (BEFORE cuts)
        hCorrelation1->Fill(correlation1);
        hCorrelation2->Fill(correlation2);
        hCorrelation3a->Fill(correlation3a);
        hCorrelation3b->Fill(correlation3b);
        
        hXiSumVsMKK->Fill(mKK, xiE + xiW);
        hCorr1VsCorr3a->Fill(correlation1, correlation3a);
        
        // SC 7.1: Energy-momentum conservation
        if (fabs(correlation1 - 1.0) > CUT_CORRELATION1) continue;
        hCutFlow->Fill(12);
        nPassSC71++;
        
        // SC 7.2: Rapidity correlation
        if (fabs(correlation2) > CUT_CORRELATION2) continue;
        hCutFlow->Fill(13);
        nPassSC72++;
        
        // SC 7.3: Individual proton-system correlations
        if (fabs(correlation3a - 1.0) > CUT_CORRELATION3 || 
            fabs(correlation3b - 1.0) > CUT_CORRELATION3) continue;
        hCutFlow->Fill(14);
        nPassSC73++;
        
        // -------------------------------------------------------------------
        // SC 8: Decay Vertex Positions
        // -------------------------------------------------------------------
        
        TVector3 vtx1 = ks1->decayVertex();
        TVector3 vtx2 = ks2->decayVertex();
        
        // SC 8.1: Average z position |<z>| < 80 cm
        double avgZ = (vtx1.Z() + vtx2.Z()) / 2.0;
        hAvgZ->Fill(avgZ);
        if (fabs(avgZ) > 80.0) continue;
        hCutFlow->Fill(15);
        nPassSC81++;
        
        // SC 8.2: Distance between vertices |Δz| < 15 cm
        double deltaZ = fabs(vtx1.Z() - vtx2.Z());
        hDeltaZ->Fill(deltaZ);
        if (deltaZ > 15.0) continue;
        hCutFlow->Fill(16);
        nPassSC82++;
        
        // -------------------------------------------------------------------
        // SC 9: Angle Between Pions in K_S Rest Frame
        // -------------------------------------------------------------------
        
        // Get daughter momenta for ks1 (leading)
        TVector3 mom_pos1, mom_neg1;
        const StUPCTrack* pos1 = (ks1 == &v0_00 || ks1 == &v0_01) ? posTracks[0] : posTracks[1];
        const StUPCTrack* neg1 = (ks1 == &v0_00 || ks1 == &v0_10) ? negTracks[0] : negTracks[1];
        pos1->getMomentum(mom_pos1);
        neg1->getMomentum(mom_neg1);
        
        TLorentzVector pi_plus1, pi_minus1;
        pi_plus1.SetXYZM(mom_pos1.X(), mom_pos1.Y(), mom_pos1.Z(), MASS_PION);
        pi_minus1.SetXYZM(mom_neg1.X(), mom_neg1.Y(), mom_neg1.Z(), MASS_PION);
        
        // Boost to K_S rest frame
        TVector3 boost1 = lv_ks1.BoostVector();
        pi_plus1.Boost(-boost1);
        pi_minus1.Boost(-boost1);
        
        // Calculate angle
        double cosAngle1 = pi_plus1.Vect().Dot(pi_minus1.Vect()) / 
                          (pi_plus1.Vect().Mag() * pi_minus1.Vect().Mag());
        
        // Get daughter momenta for ks2 (subleading)
        TVector3 mom_pos2, mom_neg2;
        const StUPCTrack* pos2 = (ks2 == &v0_00 || ks2 == &v0_01) ? posTracks[0] : posTracks[1];
        const StUPCTrack* neg2 = (ks2 == &v0_00 || ks2 == &v0_10) ? negTracks[0] : negTracks[1];
        pos2->getMomentum(mom_pos2);
        neg2->getMomentum(mom_neg2);
        
        TLorentzVector pi_plus2, pi_minus2;
        pi_plus2.SetXYZM(mom_pos2.X(), mom_pos2.Y(), mom_pos2.Z(), MASS_PION);
        pi_minus2.SetXYZM(mom_neg2.X(), mom_neg2.Y(), mom_neg2.Z(), MASS_PION);
        
        // Boost to K_S rest frame
        TVector3 boost2 = lv_ks2.BoostVector();
        pi_plus2.Boost(-boost2);
        pi_minus2.Boost(-boost2);
        
        // Calculate angle
        double cosAngle2 = pi_plus2.Vect().Dot(pi_minus2.Vect()) / 
                          (pi_plus2.Vect().Mag() * pi_minus2.Vect().Mag());

        hCosAngle1->Fill(cosAngle1);
        hCosAngle2->Fill(cosAngle2);                  
        
        // SC 9: Require pions to be back-to-back (cos(θ) > -0.8)
        // Real K_S decays have cos(θ) ≈ -1 (180 degrees)
        if (cosAngle1 > -0.8 || cosAngle2 > -0.8) continue;
        hCutFlow->Fill(17);
        nPassSC9++;
        
        // -------------------------------------------------------------------
        // EVENT PASSED ALL CUTS!
        // -------------------------------------------------------------------
        hMassKK_AfterSC7->Fill(mKK);
        
    } // End event loop
    
    // -----------------------------------------------------------------------
    // Print summary
    // -----------------------------------------------------------------------
    cout << "\n============================================" << endl;
    cout << "         ANALYSIS SUMMARY" << endl;
    cout << "============================================" << endl;
    cout << "Total events processed:     " << nTotal << endl;
    cout << "K_S K_S candidates (pre-SC7): " << nKSKS << endl;
    cout << "\n--- Roman Pot Statistics ---" << endl;
    cout << "Events with RP data:        " << nWithRP << " (" << 100.0*nWithRP/nKSKS << "%)" << endl;
    cout << "Events with 2 RP tracks:    " << nRP2Tracks << " (" << 100.0*nRP2Tracks/nWithRP << "%)" << endl;
    
    int nTotalXi = nNegativeXi + nLargeXi + nReasonableXi;
    cout << "\n--- RP Xi Quality ---" << endl;
    cout << "Total xi measurements:      " << nTotalXi << endl;
    cout << "Negative xi:                " << nNegativeXi << " (" << 100.0*nNegativeXi/nTotalXi << "%)" << endl;
    cout << "Too large xi (>0.1):        " << nLargeXi << " (" << 100.0*nLargeXi/nTotalXi << "%)" << endl;
    cout << "Reasonable xi (0.001-0.1):  " << nReasonableXi << " (" << 100.0*nReasonableXi/nTotalXi << "%)" << endl;
    
    cout << "\n--- SC 7 Cut Flow ---" << endl;
    cout << "K_S K_S with RP (2 tracks): " << nRP2Tracks << endl;
    cout << "After xi range check:       " << nXiRange << " (" << 100.0*nXiRange/nRP2Tracks << "%)" << endl;
    cout << "After SC 7.1:               " << nPassSC71 << " (" << 100.0*nPassSC71/nXiRange << "%)" << endl;
    cout << "After SC 7.2:               " << nPassSC72 << " (" << 100.0*nPassSC72/nXiRange << "%)" << endl;
    cout << "After SC 7.3:               " << nPassSC73 << " (" << 100.0*nPassSC73/nXiRange << "%)" << endl;
    
    cout << "\n--- SC 8 & 9 Cut Flow ---" << endl;
    cout << "After SC 8.1 (avg z):       " << nPassSC81 << " (" << 100.0*nPassSC81/nPassSC73 << "%)" << endl;
    cout << "After SC 8.2 (delta z):     " << nPassSC82 << " (" << 100.0*nPassSC82/nPassSC73 << "%)" << endl;
    cout << "After SC 9 (pi angle):      " << nPassSC9 << " (" << 100.0*nPassSC9/nPassSC73 << "%)" << endl;
    
    cout << "\n============================================" << endl;
    cout << "CONCLUSION:" << endl;
    if (nPassSC9 == 0) {
        cout << "❌ NO exclusive K_S K_S events found!" << endl;
        if (nPassSC73 == 0) {
            cout << "   SC 7.3 rejects all candidates." << endl;
        } else {
            cout << "   SC 7 passed " << nPassSC73 << " events," << endl;
            cout << "   but SC 8/9 rejected them all." << endl;
        }
        cout << "   This indicates the dataset does NOT" << endl;
        cout << "   contain exclusive K_S K_S production." << endl;
    } else {
        cout << "✓ Found " << nPassSC9 << " exclusive candidates!" << endl;
    }
    cout << "============================================\n" << endl;

    
    // -----------------------------------------------------------------------
    // Write output and cleanup
    // -----------------------------------------------------------------------
    outFile->cd();
    hCutFlow->Write();
    hRPProtonMomentum->Write();
    hRPXi->Write();
    hRPXiEvsXiW->Write();
    hMassKK->Write();
    hMassKK_AfterSC7->Write();
    hCorrelation1->Write();
    hCorrelation2->Write();
    hCorrelation3a->Write();
    hCorrelation3b->Write();
    hXiSumVsMKK->Write();
    hCorr1VsCorr3a->Write();
    hAvgZ->Write();
    hDeltaZ->Write();
    hCosAngle1->Write();
    hCosAngle2->Write();

    
    outFile->Close();
    
    cout << "Output written to: " << argv[2] << endl;
    cout << "\nDone!" << endl;
    
    return 0;
}
