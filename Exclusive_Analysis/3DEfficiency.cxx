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
#include <TH3F.h>  // Added for 3D histograms
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
#include <Afterburner.h>
#include <StRPEvent.h>
#include <map>

using namespace std;

// Function prototype for finding protons
void FindProtons(bool isMC, StRPEvent *rpEvt, StUPCEvent *upcEvt, TLorentzVector & proton1, TLorentzVector & proton2);

bool passFiducialCut(double eta, double zvtx) {
    return (eta < -zvtx/250.0 + 0.9) && 
           (eta > -zvtx/250.0 - 0.9) && 
           (fabs(eta) < 0.9);
}

int main(int argc, char** argv)  
{
    // Particle masses and beam energy definitions
    float massPion = 0.13957061; // Mass of pion in GeV
    double massKaon = 497.611 / 1000.0; // Mass of kaon in GeV
    double massProton = 938.272 / 1000.0; // Mass of proton in GeV

    const double VzMax = 100.0; // Maximum Vz for fiducial cut
    const double PtMax = 3.0; // Maximum transverse momentum // we use it for efficiency ToF and TPC
    const double PtMin = 0.25; // Minimum transverse momentum
    
    //cut constants
    const double kDeltaRCut = 0.15;  //0.15

    // Open input file containing paths to data files
    ifstream inputFilePathList(argv[1]);
    if (!inputFilePathList) 
    {
        cerr << "Failed to open input file." << std::endl;
        return 1;
    }

    // Create TChain to read UPC tree from multiple files
    TChain *chain = new TChain("mUPCTree"); 
    string inputFileName;
    while (std::getline(inputFilePathList, inputFileName))
    {
        chain->Add(inputFileName.c_str());
    }
    inputFilePathList.close();

    // Initialize variables for event processing
    bool isMC = 0; // Flag to check if data is Monte Carlo (MC)
    static StUPCEvent *upcEvt = nullptr; // Pointer to UPC event
    static StRPEvent *correctedRpEvent = nullptr; // Pointer to corrected RP event

    // Load offset file for corrections (not applicable) //
    //LoadOffsetFile("../data/OffSetsCorrectionsRun17.list", mCorrection);

    // Set branch addresses for UPC event and check if it's MC
    chain->SetBranchAddress("mUPCEvent", &upcEvt);
    chain->GetEntry(0);

    if (upcEvt->getRunNumber() == 1)
    {
        isMC = 1;
    }
    // Set branch address for corrected RP event if not MC
    if (isMC == 0)
    {
        chain->SetBranchAddress("correctedRpEvent", &correctedRpEvent);
    }
    // =============================================================================
    // 1D HISTOGRAM BINNING (consistent with UnifiedMuDstQA.C)
    // =============================================================================
    //TPC efficiency denominator (MC truth of TPC reconstructed tracks) - positive
    TH1F *hSelPtMc_P = new TH1F("hSelPtMc_P", "True Pt of selected MC tracks", 60, 0, 3.0);//100, 0.1, 3.1
    TH1F *hSelEtaMc_P = new TH1F("hSelEtaMc_P", "True Eta of selected MC tracks", 20, -1.0, 1.0); //20, -1.0, 1.0
    TH1F *hSelPhiMc_P = new TH1F("hSelPhiMc_P", "True Phi of selected MC tracks", 20, -TMath::Pi(), TMath::Pi());
    TH1F *hSelVzMc_P = new TH1F("hSelVzMc_P", "True Vz of selected MC tracks", 20, -100, 100);
    // TPC efficiency denominator (MC truth of TPC reconstructed tracks) - negative
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
    // 3D HISTOGRAM BINNING (consistent with UnifiedMuDstQA.C)
    // =============================================================================
    
    // Define consistent binning scheme for all 3D histograms
    Double_t ptBins[] = {0.20, 0.25, 0.30, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                         1.0, 1.1, 1.2, 1.3, 1.4, 1.6, 2.0, 3.0};
    Int_t nPtBins = sizeof(ptBins)/sizeof(ptBins[0]) - 1;
    
    Double_t etaBins[] = {-1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,
                          0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    Int_t nEtaBins = sizeof(etaBins)/sizeof(etaBins[0]) - 1;
    
    Double_t vzBins[] = {-100, -90, -80, -70, -60, -50, -40, -30, -20, -10,
                         0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    Int_t nVzBins = sizeof(vzBins)/sizeof(vzBins[0]) - 1;

    // =============================================================================
    // 3D EFFICIENCY HISTOGRAMS (SEPARATED BY CHARGE)
    // =============================================================================
    
    // TPC 3D Efficiency Histograms - POSITIVE PARTICLES
    // 1. True pion+ pT, eta and Vz (denominator)
    TH3F *h3D_TPC_TruePions_P = new TH3F("h3D_TPC_TruePions_P", "True Pions 3D (+)", 
                                         nPtBins, ptBins, nEtaBins, etaBins, nVzBins, vzBins);
    h3D_TPC_TruePions_P->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h3D_TPC_TruePions_P->GetYaxis()->SetTitle("#eta");
    h3D_TPC_TruePions_P->GetZaxis()->SetTitle("V_{z} (cm)");
    
    // 2. Reconstructed pion+ pT, eta and Vz matched with true pion (numerator)
    TH3F *h3D_TPC_RecoMatchedPions_P = new TH3F("h3D_TPC_RecoMatchedPions_P", "Reconstructed Matched Pions 3D (+)", 
                                                nPtBins, ptBins, nEtaBins, etaBins, nVzBins, vzBins);
    h3D_TPC_RecoMatchedPions_P->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h3D_TPC_RecoMatchedPions_P->GetYaxis()->SetTitle("#eta");
    h3D_TPC_RecoMatchedPions_P->GetZaxis()->SetTitle("V_{z} (cm)");
    
    // TPC 3D Efficiency Histograms - NEGATIVE PARTICLES  
    TH3F *h3D_TPC_TruePions_N = new TH3F("h3D_TPC_TruePions_N", "True Pions 3D (-)", 
                                         nPtBins, ptBins, nEtaBins, etaBins, nVzBins, vzBins);
    h3D_TPC_TruePions_N->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h3D_TPC_TruePions_N->GetYaxis()->SetTitle("#eta");
    h3D_TPC_TruePions_N->GetZaxis()->SetTitle("V_{z} (cm)");
    
    TH3F *h3D_TPC_RecoMatchedPions_N = new TH3F("h3D_TPC_RecoMatchedPions_N", "Reconstructed Matched Pions 3D (-)", 
                                                nPtBins, ptBins, nEtaBins, etaBins, nVzBins, vzBins);
    h3D_TPC_RecoMatchedPions_N->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h3D_TPC_RecoMatchedPions_N->GetYaxis()->SetTitle("#eta");
    h3D_TPC_RecoMatchedPions_N->GetZaxis()->SetTitle("V_{z} (cm)");
    
    // TOF 3D Efficiency Histograms - POSITIVE PARTICLES
    // 3. Reconstructed pion+ pT, eta and Vz matched with true pion (denominator for TOF)
    TH3F *h3D_TOF_RecoMatchedPions_P = new TH3F("h3D_TOF_RecoMatchedPions_P", "TOF Denominator: Reconstructed Matched Pions 3D (+)", 
                                                nPtBins, ptBins, nEtaBins, etaBins, nVzBins, vzBins);
    h3D_TOF_RecoMatchedPions_P->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h3D_TOF_RecoMatchedPions_P->GetYaxis()->SetTitle("#eta");
    h3D_TOF_RecoMatchedPions_P->GetZaxis()->SetTitle("V_{z} (cm)");
    
    // 4. Reconstructed pion+ pT,eta and Vz matched with ToF hit (numerator for TOF)
    TH3F *h3D_TOF_RecoMatchedPionsWithTOF_P = new TH3F("h3D_TOF_RecoMatchedPionsWithTOF_P", "TOF Numerator: Reconstructed Matched Pions with TOF 3D (+)", 
                                                       nPtBins, ptBins, nEtaBins, etaBins, nVzBins, vzBins);
    h3D_TOF_RecoMatchedPionsWithTOF_P->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h3D_TOF_RecoMatchedPionsWithTOF_P->GetYaxis()->SetTitle("#eta");
    h3D_TOF_RecoMatchedPionsWithTOF_P->GetZaxis()->SetTitle("V_{z} (cm)");
    
    // TOF 3D Efficiency Histograms - NEGATIVE PARTICLES
    TH3F *h3D_TOF_RecoMatchedPions_N = new TH3F("h3D_TOF_RecoMatchedPions_N", "TOF Denominator: Reconstructed Matched Pions 3D (-)", 
                                                nPtBins, ptBins, nEtaBins, etaBins, nVzBins, vzBins);
    h3D_TOF_RecoMatchedPions_N->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h3D_TOF_RecoMatchedPions_N->GetYaxis()->SetTitle("#eta");
    h3D_TOF_RecoMatchedPions_N->GetZaxis()->SetTitle("V_{z} (cm)");
    
    TH3F *h3D_TOF_RecoMatchedPionsWithTOF_N = new TH3F("h3D_TOF_RecoMatchedPionsWithTOF_N", "TOF Numerator: Reconstructed Matched Pions with TOF 3D (-)", 
                                                       nPtBins, ptBins, nEtaBins, etaBins, nVzBins, vzBins);
    h3D_TOF_RecoMatchedPionsWithTOF_N->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h3D_TOF_RecoMatchedPionsWithTOF_N->GetYaxis()->SetTitle("#eta");
    h3D_TOF_RecoMatchedPionsWithTOF_N->GetZaxis()->SetTitle("V_{z} (cm)");

    //cut flow histograms
    TH1D* HistEventCutFlow = new TH1D("HistEventCutFlow", ";Cut;Events", 5, 0, 5);
    HistEventCutFlow->GetXaxis()->SetBinLabel(1, "All Events");
    HistEventCutFlow->GetXaxis()->SetBinLabel(2, "Trigger Cuts");
    HistEventCutFlow->GetXaxis()->SetBinLabel(3, "Vertex Cuts");
    HistEventCutFlow->GetXaxis()->SetBinLabel(4, ">= 3 TOF Tracks");
    HistEventCutFlow->GetXaxis()->SetBinLabel(5, "Final Selection");

    TH1D* HistTrackCutFlow = new TH1D("HistTrackCutFlow", ";Cut;Tracks", 6, 0, 6);
    HistTrackCutFlow->GetXaxis()->SetBinLabel(1, "All Tracks");
    HistTrackCutFlow->GetXaxis()->SetBinLabel(2, "|#eta| < 0.9");
    HistTrackCutFlow->GetXaxis()->SetBinLabel(3, "p_{T} > 0.2 GeV");
    HistTrackCutFlow->GetXaxis()->SetBinLabel(4, "NhitsFit >= 20");
    HistTrackCutFlow->GetXaxis()->SetBinLabel(5, "kV0 Flag");
    HistTrackCutFlow->GetXaxis()->SetBinLabel(6, "TOF Matched");
    
    cout << "3D efficiency histograms booked successfully:" << endl;
    cout << "  - TPC True Pions 3D (separated by charge)" << endl;
    cout << "  - TPC Reconstructed Matched Pions 3D (separated by charge)" << endl;
    cout << "  - TOF Reconstructed Matched Pions 3D (separated by charge)" << endl;

    // Matching quality monitoring histograms
    TH1F* hMatchingDeltaR = new TH1F("hMatchingDeltaR", "DeltaR of matched tracks;#DeltaR;Entries", 100, 0, 0.5);
    TH1F* hMatchingDeltaPt = new TH1F("hMatchingDeltaPt", "Relative pT difference;#Delta p_{T}/p_{T};Entries", 100, -0.5, 0.5);
    TH1F* hMatchingDeltaEta = new TH1F("hMatchingDeltaEta", "Eta difference;#Delta#eta;Entries", 100, -0.2, 0.2);
    TH1F* hMatchingDeltaPhi = new TH1F("hMatchingDeltaPhi", "Phi difference;#Delta#phi;Entries", 100, -0.2, 0.2);
    TH1F* hMatchingCharge = new TH1F("hMatchingCharge", "Charge matching;0=wrong, 1=correct;Entries", 2, -0.5, 1.5);
    TH1F* hMatchingScore = new TH1F("hMatchingScore", "Combined matching score;Score;Entries", 100, 0, 3.0);
    TH1F* hMatchingEfficiency = new TH1F("hMatchingEfficiency", "Matching efficiency per event;Efficiency;Events", 100, 0, 1.0);
    
    // Separate matching quality by charge
    TH2F* hMatchingDeltaR_vs_Pt_Pos = new TH2F("hMatchingDeltaR_vs_Pt_Pos", 
        "DeltaR vs pT (positive);p_{T} (GeV/c);#DeltaR", 60, 0, 3.0, 100, 0, 0.5);
    TH2F* hMatchingDeltaR_vs_Pt_Neg = new TH2F("hMatchingDeltaR_vs_Pt_Neg", 
        "DeltaR vs pT (negative);p_{T} (GeV/c);#DeltaR", 60, 0, 3.0, 100, 0, 0.5);
    
    // Track purity checks
    TH1F* hUnmatchedMCPions = new TH1F("hUnmatchedMCPions", "pT of unmatched MC pions;p_{T} (GeV/c);Entries", 60, 0, 3.0);
    TH1F* hUnmatchedRecoTracks = new TH1F("hUnmatchedRecoTracks", "pT of unmatched reco tracks;p_{T} (GeV/c);Entries", 60, 0, 3.0);
    
    
    // Counters for statistics
    int totalTruePionsPos = 0, totalTruePionsNeg = 0;
    int matchedPionsPos = 0, matchedPionsNeg = 0;
    int chargeMatchedPionsPos = 0, chargeMatchedPionsNeg = 0;


    // Loop over all entries in the TChain
    for (Long64_t i = 0; i < chain->GetEntries(); ++i) 
    {
        if (i%1000000 == 0)  
        {
            cout << i << "/" <<  chain->GetEntries() << endl;
        }

        chain->GetEntry(i);

        // Fill event cut flow: All events
        HistEventCutFlow->Fill(0);

        // Select tracks based on quality criteria
        vector <StUPCTrack const*> tracksWithTofHit;
        vector <StUPCTrack const*> tracksWithSmallDca;
        vector <StUPCTrack const*> tracksForTrueLevel;
        
        // SELECT tracks with kV0 flag, absEta < 0.9 and pT > 0.15 GeV and Nhit >= 20
        for (Int_t j = 0; j<upcEvt->getNumberOfTracks(); j++)
		{
            if (isMC == 0)
            {
                if ( upcEvt->getTrack(j)->getFlag(StUPCTrack::kV0) ) //(abs(upcEvt->getTrack(j)->getEta()) <= 0.9) and (upcEvt->getTrack(j)->getPt() >= 0.2) and (upcEvt->getTrack(j)->getNhitsFit() >= 20) and
                {
                    tracksWithSmallDca.push_back(upcEvt->getTrack(j));
                    if(upcEvt->getTrack(j)->getFlag(StUPCTrack::kTof))
                    {
                        tracksWithTofHit.push_back(upcEvt->getTrack(j));
                    } 
                }
            }
            else if (isMC == 1)
            {
                if ( upcEvt->getTrack(j)->getNhitsFit() >= 20 ) // removed (abs(upcEvt->getTrack(j)->getEta()) <= 0.9) and (upcEvt->getTrack(j)->getPt() >= 0.25) and
                {
                    tracksWithSmallDca.push_back(upcEvt->getTrack(j));
                    tracksForTrueLevel.push_back(upcEvt->getTrack(j));
                    if(upcEvt->getTrack(j)->getFlag(StUPCTrack::kTof))
                    {
                        tracksWithTofHit.push_back(upcEvt->getTrack(j));
                    } 
                }
            }
		} 
        
        
        // ============================
        // True Level Analysis (MC Only) with 3D Histograms
        // ============================
        /*if (isMC==1) 
        {
            // Extract true pions from MC truth information
            TParticle* particle;
            vector <TParticle *> vPions; // Vector to store true pions
            for (int k = 0; k <  upcEvt->getNumberOfMCParticles(); k++)
            {
                particle = upcEvt->getMCParticle(k);
                if (particle->GetPDG()->PdgCode() == 211 or particle->GetPDG()->PdgCode() == -211)
                {
                    vPions.push_back(particle);                    
                }
            }

            TLorentzVector productionVertex; 
            TLorentzVector lorentzTrack;        
            // Loop over true pions and match them to reconstructed tracks
            for (int j = 0; j < vPions.size(); j++) 
            {
                 // Apply vertex and kinematic cuts 
                TLorentzVector lorentzVector;
                lorentzVector.SetPxPyPzE(vPions[j]->Px(), vPions[j]->Py(), vPions[j]->Pz(), vPions[j]->Energy()); 
                double pT = sqrt(pow(vPions[j]->Px(),2) + pow(vPions[j]->Py(),2)) ;
                double eta = lorentzVector.Eta();
                vPions[j]->ProductionVertex(productionVertex);

                double truthVertexR = sqrt(pow(productionVertex.X(),2) + pow(productionVertex.Y(),2)); // Transverse vertex position
                double truthVertexZ = productionVertex.Z(); // Longitudinal vertex position

                // Determine charge from PDG code
                int pdgCode = vPions[j]->GetPDG()->PdgCode();
                bool isPositive = (pdgCode == 211);   // π+
                bool isNegative = (pdgCode == -211);  // π-

                // Apply cuts on the production vertex. // //
                if (abs(truthVertexZ) > VzMax || pT < PtMin || abs(eta) > 0.9 || truthVertexR > 3.0 || !passFiducialCut(eta, truthVertexZ))
                {
                    continue; // // Skip pions outside acceptance or failing fiducial cut
                }
               //if (truthVertexR > 3.0) continue;
               //if (!passFiducialCut(eta, truthVertexZ)) continue;

                //TPC efficiency denominator
                //if (fabs(truthVertexZ) < VzMax && passFiducialCut(eta, truthVertexZ)) {  // For pt distribution: cut on nominal eta and vz
                //    if (isPositive) hSelPtMc_P->Fill(pT);
                //    if (isNegative) hSelPtMc_N->Fill(pT);
                //}
                //if (pT > PtMin && pT < PtMax && fabs(truthVertexZ) < VzMax) { //&& passFiducialCut(eta, truthVertexZ)) {  // For eta distribution: cut on nominal pt and vz
                //    if (isPositive) hSelEtaMc_P->Fill(eta);
                //    if (isNegative) hSelEtaMc_N->Fill(eta);
                //}
                //if (pT > PtMin && pT < PtMax && passFiducialCut(eta, truthVertexZ)) {  // For vz distribution: cut on nominal pt and eta
                //    if (isPositive) hSelVzMc_P->Fill(truthVertexZ);
                //    if (isNegative) hSelVzMc_N->Fill(truthVertexZ);
                //}
                
                // Fill 3D histogram for true pions (TPC denominator) separated by charge
                if (isPositive) h3D_TPC_TruePions_P->Fill(pT, eta, truthVertexZ);
                if (isNegative) h3D_TPC_TruePions_N->Fill(pT, eta, truthVertexZ);
                
                // Match reconstructed tracks to true pions
                vector <double> distr;
                vector <StUPCTrack const *> matchedTPCpions;
                if (tracksForTrueLevel.size() == 0) {break;}
                for (int m = 0; m < tracksForTrueLevel.size(); m++) 
                {
                    TLorentzVector pion, track;
                    pion.SetPxPyPzE(vPions[j]->Px(), vPions[j]->Py(), vPions[j]->Pz(), vPions[j]->Energy());
                    tracksForTrueLevel[m]->getLorentzVector(track, massPion);
                    double r = pion.DeltaR(track);  // Spatial distance in η-φ space
                    matchedTPCpions.push_back(tracksForTrueLevel[m]);
                    distr.push_back(r);
                }
   
                // Find the minimum distance and corresponding index
                auto it = std::min_element(distr.begin(), distr.end());
                int index = std::distance(distr.begin(), it);              

                if (distr[index] < kDeltaRCut) 
                {  
                    // Get reconstructed track variables
                    double recoPt = matchedTPCpions[index]->getPt();
                    double recoEta = matchedTPCpions[index]->getEta();
                    // Use true Vz for 3D histograms (consistent with UnifiedMuDstQA.C)
                    //double recoVz = truthVertexZ;
                    
                    // TPC efficiency numerators (MC truth of TPC reconstructed tracks) 
                    //if (fabs(truthVertexZ) < VzMax && passFiducialCut(eta, truthVertexZ)) {  //For pT efficiency: cut on eta and vz
                    //    if (isPositive) hSelPtTrue_P->Fill(pT);
                    //    if (isNegative) hSelPtTrue_N->Fill(pT);
                    //}
                    //if (pT > PtMin && pT < PtMax && fabs(truthVertexZ) < VzMax) {   // For eta efficiency: cut on pt and vz
                    //    if (isPositive) hSelEtaTrue_P->Fill(eta);
                    //    if (isNegative) hSelEtaTrue_N->Fill(eta);
                    //}
                    //if (pT > PtMin && pT < PtMax && passFiducialCut(eta, truthVertexZ)) {   // For vz efficiency: cut on pt and eta
                    //    if (isPositive) hSelVzTrue_P->Fill(truthVertexZ);
                    //    if (isNegative) hSelVzTrue_N->Fill(truthVertexZ);
                    //}

                    // ToF efficiency denominator 
                    //if (fabs(truthVertexZ) < VzMax && passFiducialCut(recoEta, truthVertexZ)) {  // For pT efficiency: cut on eta and vz
                    //    if (isPositive) hSelPtReco_P->Fill(recoPt);
                    //    if (isNegative) hSelPtReco_N->Fill(recoPt);
                    //} 
                    //if (recoPt > PtMin && recoPt < PtMax && fabs(truthVertexZ) < VzMax) {   // For eta efficiency: cut on pt and vz
                    //    if (isPositive) hSelEtaReco_P->Fill(recoEta);
                    //    if (isNegative) hSelEtaReco_N->Fill(recoEta);
                    //} 
                    //if (recoPt > PtMin && recoPt < PtMax && passFiducialCut(recoEta, truthVertexZ)) {   // For vz efficiency: cut on pt and eta
                    //    if (isPositive) hSelVzReco_P->Fill(truthVertexZ);
                    //    if (isNegative) hSelVzReco_N->Fill(truthVertexZ);
                    //}
                    
                    // Fill 3D histogram for reconstructed matched pions (TPC numerator) separated by charge
                    if (isPositive) h3D_TPC_RecoMatchedPions_P->Fill(pT, eta, truthVertexZ);//(recoPt, recoEta, recoVz);
                    if (isNegative) h3D_TPC_RecoMatchedPions_N->Fill(pT, eta, truthVertexZ);//(recoPt, recoEta, recoVz);
                    
                    // Fill 3D histogram for TOF denominator (all reconstructed matched pions) separated by charge
                    if (isPositive) h3D_TOF_RecoMatchedPions_P->Fill(recoPt, recoEta, truthVertexZ);
                    if (isNegative) h3D_TOF_RecoMatchedPions_N->Fill(recoPt, recoEta, truthVertexZ);
                

                    if (matchedTPCpions[index]->getFlag(StUPCTrack::kTof)) 
                    {
                        // ToF efficiency numerators
                        //if (fabs(truthVertexZ) < VzMax && passFiducialCut(recoEta, truthVertexZ)) {  //// For pt efficiency: cut on RECO eta and TRUE vz, fill RECO pt
                        //    if (isPositive) hPtTpcTofMatched_P->Fill(recoPt);  //
                        //    if (isNegative) hPtTpcTofMatched_N->Fill(recoPt);
                        //}  
                        //if (recoPt > PtMin && recoPt < PtMax && fabs(truthVertexZ) < VzMax) {   ////For eta efficiency: cut on RECO pt and TRUE vz, fill RECO eta
                        //    if (isPositive) hEtaTpcTofMatched_P->Fill(recoEta);  //  
                        //    if (isNegative) hEtaTpcTofMatched_N->Fill(recoEta);
                        //} 
                        //if (recoPt > PtMin && recoPt < PtMax && passFiducialCut(recoEta, truthVertexZ)) {  // For vz efficiency: cut on RECO pt and RECO eta, fill TRUE vz
                        //    if (isPositive) hVzTpcTofMatched_P->Fill(truthVertexZ);  // 
                        //    if (isNegative) hVzTpcTofMatched_N->Fill(truthVertexZ);
                        //}
                        
                        // Fill 3D histogram for TOF numerator (reconstructed matched pions with TOF) separated by charge
                        if (isPositive) h3D_TOF_RecoMatchedPionsWithTOF_P->Fill(recoPt, recoEta, truthVertexZ);
                        if (isNegative) h3D_TOF_RecoMatchedPionsWithTOF_N->Fill(recoPt, recoEta, truthVertexZ);
                        
                    }
                    
                    // Remove the matched track to avoid double-counting
                    tracksForTrueLevel.erase(tracksForTrueLevel.begin()+index);
                }    
            }            
        }
        */
        
        // MATCHING ALGORITHM CONSISTENT WITH UnifiedMuDstQA.C
        
        /*  if (isMC==1) 
        {
            // Debug counters
            static int eventsProcessed = 0;
            static int totalValidMCPions = 0;
            static int totalMatches = 0;
            static int eventsWithTracks = 0;
            
            eventsProcessed++;
            
            // STEP 1: Collect reco tracks with MINIMAL cuts (same philosophy as before)
            vector<StUPCTrack const*> allRecoTracks;
            
            for (Int_t j = 0; j < upcEvt->getNumberOfTracks(); j++)
            {
                StUPCTrack const* track = upcEvt->getTrack(j);
                
                // Very basic cuts for track collection
                if (track->getNhitsFit() >= 20 &&           //10
                    fabs(track->getEta()) <= 0.9 &&         //1.2
                    track->getPt() >= 0.25)                  //0.1
                {
                    allRecoTracks.push_back(track);
                }
            }
            
            if (allRecoTracks.size() == 0) {
                continue; // Skip events with no tracks
            }
            
            eventsWithTracks++;
            
            // STEP 2: Collect valid MC pions (same as before)
            vector<TParticle*> validMCPions;
            
            for (int k = 0; k < upcEvt->getNumberOfMCParticles(); k++)
            {
                TParticle* particle = upcEvt->getMCParticle(k);
                int pdgCode = particle->GetPDG()->PdgCode();
                
                if (pdgCode != 211 && pdgCode != -211) continue;
                
                TLorentzVector pionTrueVec;
                pionTrueVec.SetPxPyPzE(particle->Px(), particle->Py(), 
                                    particle->Pz(), particle->Energy());
                
                TLorentzVector productionVertex;
                particle->ProductionVertex(productionVertex);
                double truthVertexZ = productionVertex.Z();
                double truthVertexR = sqrt(pow(productionVertex.X(),2) + pow(productionVertex.Y(),2));
                
                // Apply physics cuts
                if (fabs(truthVertexZ) > VzMax || 
                    pionTrueVec.Pt() < PtMin || 
                    fabs(pionTrueVec.Eta()) > 0.9 || 
                    truthVertexR > 3.0 ||
                    !passFiducialCut(pionTrueVec.Eta(), truthVertexZ)) {
                    continue;
                }
                
                validMCPions.push_back(particle);
                totalValidMCPions++;
                
                // Fill denominator histograms
                if (pdgCode == 211) {
                    h3D_TPC_TruePions_P->Fill(pionTrueVec.Pt(), pionTrueVec.Eta(), truthVertexZ);
                } else {
                    h3D_TPC_TruePions_N->Fill(pionTrueVec.Pt(), pionTrueVec.Eta(), truthVertexZ);
                }
            }
            
            // STEP 3: Build matching candidates (same as before)
            struct MatchCandidate {
                TParticle* mcPion;
                StUPCTrack const* recoTrack;
                double deltaR;
                double deltaPt;
                double score;
                bool chargeMatch;
                
                void calculateScore() {
                    score = deltaR * 3.0 + deltaPt * 1.0 + (chargeMatch ? 0.0 : 2.0);
                }
            };
            
            vector<MatchCandidate> allCandidates;
            
            for (TParticle* mcPion : validMCPions) 
            {
                TLorentzVector mcVec;
                mcVec.SetPxPyPzE(mcPion->Px(), mcPion->Py(), mcPion->Pz(), mcPion->Energy());
                
                int mcCharge = (mcPion->GetPDG()->PdgCode() > 0) ? +1 : -1;
                
                for (StUPCTrack const* recoTrack : allRecoTracks) 
                {
                    TLorentzVector recoVec;
                    recoTrack->getLorentzVector(recoVec, massPion);
                    
                    double deltaEta = recoVec.Eta() - mcVec.Eta();
                    double deltaPhi = TVector2::Phi_mpi_pi(recoVec.Phi() - mcVec.Phi());
                    double deltaR = sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
                    double deltaPt = fabs((recoVec.Pt() - mcVec.Pt()) / mcVec.Pt());
                    
                    if (deltaR > 1.2 || deltaPt > 1.5) continue;
                    
                    MatchCandidate candidate;
                    candidate.mcPion = mcPion;
                    candidate.recoTrack = recoTrack;
                    candidate.deltaR = deltaR;
                    candidate.deltaPt = deltaPt;
                    candidate.chargeMatch = (recoTrack->getCharge() == mcCharge);
                    candidate.calculateScore();
                    
                    allCandidates.push_back(candidate);
                }
            }
            
            // Sort by matching score
            std::sort(allCandidates.begin(), allCandidates.end(), 
                    [](const MatchCandidate& a, const MatchCandidate& b) {
                        return a.score < b.score;
                    });
            
            // STEP 4: Apply UnifiedMuDstQA.C style cuts AFTER matching
            std::set<TParticle*> usedMCPions;
            std::set<StUPCTrack const*> usedRecoTracks;
            
            const double finalDeltaRCut = 0.3; //0.3;   // Tighter than 0.6
            const double finalDeltaPtCut = 0.15; //0.15; // Tighter than 0.3
            const bool requireChargeMatch = true;
            
            for (const auto& candidate : allCandidates) {
                // Check if already used
                if (usedMCPions.count(candidate.mcPion)) continue;
                if (usedRecoTracks.count(candidate.recoTrack)) continue;
                
                // Apply matching quality cuts
                if (candidate.deltaR > finalDeltaRCut) continue;
                if (candidate.deltaPt > finalDeltaPtCut) continue;
                if (requireChargeMatch && !candidate.chargeMatch) continue;
                
                StUPCTrack const* recoTrack = candidate.recoTrack;

                
                // APPLY EXACT UnifiedMuDstQA.C CUTS FOR PIONS
                // Basic kinematic cuts
                if (fabs(recoTrack->getEta()) > 0.9) continue;
                if (recoTrack->getPt() < 0.25) continue;                // PtMin for pions
                if (recoTrack->getNhitsFit() < 20) continue;
                
                // Additional cuts from UnifiedMuDstQA.C:
                // NOTE: These may not be available in StUPCTrack - check your framework
                // if (recoTrack->getNhitsDedx() < 15) continue;        // dE/dx quality - CHECK IF AVAILABLE
                // if (recoTrack->getFlag() <= 0) continue;             // Track flag - CHECK IF AVAILABLE
                if (abs(recoTrack->getCharge()) != 1) continue;         // Valid charge
                
                // CRITICAL MISSING CUT: qaTruth() is NOT available in StUPCTrack
                // This explains efficiency differences with UnifiedMuDstQA.C
                // UnifiedMuDstQA.C: if(ptrack->qaTruth() < 95.) continue;
                
                // UPC-specific cut (not in UnifiedMuDstQA.C):
                // For non-MC data, apply kV0 flag
                if (!isMC && !recoTrack->getFlag(StUPCTrack::kV0)) continue;
                
                // For non-MC data, apply kV0 flag
                if (!isMC && !recoTrack->getFlag(StUPCTrack::kV0)) continue;
                
                // Accept this match
                usedMCPions.insert(candidate.mcPion);
                usedRecoTracks.insert(candidate.recoTrack);
                totalMatches++;
                
                // Fill efficiency numerators (same as before)
                TParticle* mcPion = candidate.mcPion;
                
                TLorentzVector productionVertex;
                mcPion->ProductionVertex(productionVertex);
                double truthVertexZ = productionVertex.Z();
                
                TLorentzVector mcVec;
                mcVec.SetPxPyPzE(mcPion->Px(), mcPion->Py(), mcPion->Pz(), mcPion->Energy());
                
                double recoPt = recoTrack->getPt();
                double recoEta = recoTrack->getEta();
                double recoVz = truthVertexZ; // Use MC truth for vertex
                
                int pdgCode = mcPion->GetPDG()->PdgCode();
                
                if (pdgCode == 211) {
                    h3D_TPC_RecoMatchedPions_P->Fill(mcVec.Pt(), mcVec.Eta(), truthVertexZ);
                    h3D_TOF_RecoMatchedPions_P->Fill(recoPt, recoEta, recoVz);
                    
                    if (recoTrack->getFlag(StUPCTrack::kTof)) {
                        h3D_TOF_RecoMatchedPionsWithTOF_P->Fill(recoPt, recoEta, recoVz);
                    }
                }
                
                if (pdgCode == -211) {
                    h3D_TPC_RecoMatchedPions_N->Fill(mcVec.Pt(), mcVec.Eta(), truthVertexZ);
                    h3D_TOF_RecoMatchedPions_N->Fill(recoPt, recoEta, recoVz);
                    
                    if (recoTrack->getFlag(StUPCTrack::kTof)) {
                        h3D_TOF_RecoMatchedPionsWithTOF_N->Fill(recoPt, recoEta, recoVz);
                    }
                }
                
                // Fill monitoring histograms
                hMatchingDeltaR->Fill(candidate.deltaR);
                hMatchingDeltaPt->Fill(candidate.deltaPt);
                hMatchingCharge->Fill(candidate.chargeMatch ? 1.0 : 0.0);
            }
            
            // Report progress
            if (eventsProcessed % 5000 == 0) {
                double matchingEff = (double)totalMatches / max(totalValidMCPions, 1);
                double eventUtilization = (double)eventsWithTracks / eventsProcessed;
                
                cout << "\n=== CONSISTENT MATCHING REPORT at Event " << i << " ===" << endl;
                cout << "Events processed: " << eventsProcessed << endl;
                cout << "Events with tracks: " << eventsWithTracks 
                    << " (" << 100.0*eventUtilization << "%)" << endl;
                cout << "Valid MC pions: " << totalValidMCPions << endl;
                cout << "Total matches: " << totalMatches << endl;
                cout << "OVERALL MATCHING EFFICIENCY: " << 100.0*matchingEff << "%" << endl;
                cout << "Matches this event: " << usedMCPions.size() << endl;
                cout << "============================================\n" << endl;
            }
        }

        */

        //trying
        int tofDenominatorCount_P = 0;
        int tofNumeratorCount_P = 0;
        int tofDenominatorCount_N = 0;
        int tofNumeratorCount_N = 0;

        if (isMC == 1) 
        {
            TLorentzVector productionVertex;  // 
             TLorentzVector lorentzTrack;     // 
            
            // Separate true pions by charge
            vector <TParticle *> vPositivePions;
            vector <TParticle *> vNegativePions;
            
            for (int k = 0; k < upcEvt->getNumberOfMCParticles(); k++)
            {
                TParticle* particle = upcEvt->getMCParticle(k);
                if (particle->GetPDG()->PdgCode() == 211)       // pi+
                    vPositivePions.push_back(particle);
                else if (particle->GetPDG()->PdgCode() == -211) // pi-
                    vNegativePions.push_back(particle);
            }
            
            // Separate reconstructed tracks by charge
            vector <StUPCTrack const *> positiveRecoTracks;
            vector <StUPCTrack const *> negativeRecoTracks;
            
            for (int m = 0; m < tracksForTrueLevel.size(); m++) 
            {
                if (tracksForTrueLevel[m]->getCharge() > 0)
                    positiveRecoTracks.push_back(tracksForTrueLevel[m]);
                else if (tracksForTrueLevel[m]->getCharge() < 0)
                    negativeRecoTracks.push_back(tracksForTrueLevel[m]);
            }
            
            // Match positive pions to positive tracks only
            for (int j = 0; j < vPositivePions.size(); j++) 
            {
                // Apply your existing cuts on true pion
                TLorentzVector lorentzVector;
                lorentzVector.SetPxPyPzE(vPositivePions[j]->Px(), vPositivePions[j]->Py(), 
                                        vPositivePions[j]->Pz(), vPositivePions[j]->Energy()); 
                double pT = sqrt(pow(vPositivePions[j]->Px(),2) + pow(vPositivePions[j]->Py(),2));
                double eta = lorentzVector.Eta();
                
                vPositivePions[j]->ProductionVertex(productionVertex);
                double truthVertexZ = productionVertex.Z();
                double truthVertexR = sqrt(pow(productionVertex.X(),2) + pow(productionVertex.Y(),2));
                
                // Apply your existing cuts
                if (abs(truthVertexZ) > VzMax || pT < PtMin || abs(eta) > 0.9 || 
                    truthVertexR > 3.0 || !passFiducialCut(eta, truthVertexZ)) {
                    continue;
                }
                
                // Fill true positive pion histogram
                h3D_TPC_TruePions_P->Fill(pT, eta, truthVertexZ);
                totalTruePionsPos++; // counter
                
                // Match only to positive reconstructed tracks
                double minDeltaR = 999.0;
                int bestMatchIndex = -1;
                
                for (int m = 0; m < positiveRecoTracks.size(); m++) 
                {
                    TLorentzVector pion, track;
                    pion.SetPxPyPzE(vPositivePions[j]->Px(), vPositivePions[j]->Py(), 
                                vPositivePions[j]->Pz(), vPositivePions[j]->Energy());
                    positiveRecoTracks[m]->getLorentzVector(track, massPion);
                    
                    double deltaR = pion.DeltaR(track);
                    if (deltaR < minDeltaR) {
                        minDeltaR = deltaR;
                        bestMatchIndex = m;
                    }
                }
                
                // Check if match is good enough
                if (bestMatchIndex >= 0 && minDeltaR < kDeltaRCut) 
                {
                    matchedPionsPos++;
                    chargeMatchedPionsPos++; // Always matches now by construction

                    // Fill histograms
                    h3D_TPC_RecoMatchedPions_P->Fill(pT, eta, truthVertexZ);
                    
                    double recoPt = positiveRecoTracks[bestMatchIndex]->getPt();
                    double recoEta = positiveRecoTracks[bestMatchIndex]->getEta();
                    h3D_TOF_RecoMatchedPions_P->Fill(recoPt, recoEta, truthVertexZ);
                    tofDenominatorCount_P++;
                    
                    if (positiveRecoTracks[bestMatchIndex]->getFlag(StUPCTrack::kTof)) {
                        h3D_TOF_RecoMatchedPionsWithTOF_P->Fill(recoPt, recoEta, truthVertexZ);
                        tofNumeratorCount_P++;
                    }
                    
                    // Remove the matched track to avoid double-counting
                    //positiveRecoTracks.erase(positiveRecoTracks.begin() + bestMatchIndex);
                    // Instead, mark as used to prevent double-matching
                    positiveRecoTracks[bestMatchIndex] = nullptr;
                }
            }
            
            // Match negative pions to negative tracks only
            for (int j = 0; j < vNegativePions.size(); j++) 
            {
                // Similar logic for negative pions...
                TLorentzVector lorentzVector;
                lorentzVector.SetPxPyPzE(vNegativePions[j]->Px(), vNegativePions[j]->Py(), 
                                        vNegativePions[j]->Pz(), vNegativePions[j]->Energy()); 
                double pT = sqrt(pow(vNegativePions[j]->Px(),2) + pow(vNegativePions[j]->Py(),2));
                double eta = lorentzVector.Eta();
                
                vNegativePions[j]->ProductionVertex(productionVertex);
                double truthVertexZ = productionVertex.Z();
                double truthVertexR = sqrt(pow(productionVertex.X(),2) + pow(productionVertex.Y(),2));
                
                if (abs(truthVertexZ) > VzMax || pT < PtMin || abs(eta) > 0.9 || 
                    truthVertexR > 3.0 || !passFiducialCut(eta, truthVertexZ)) {
                    continue;
                }
                
                h3D_TPC_TruePions_N->Fill(pT, eta, truthVertexZ);
                totalTruePionsNeg++;
                
                // Match only to negative reconstructed tracks
                double minDeltaR = 999.0;
                int bestMatchIndex = -1;
                
                for (int m = 0; m < negativeRecoTracks.size(); m++) 
                {
                    TLorentzVector pion, track;
                    pion.SetPxPyPzE(vNegativePions[j]->Px(), vNegativePions[j]->Py(), 
                                vNegativePions[j]->Pz(), vNegativePions[j]->Energy());
                    negativeRecoTracks[m]->getLorentzVector(track, massPion);
                    
                    double deltaR = pion.DeltaR(track);
                    if (deltaR < minDeltaR) {
                        minDeltaR = deltaR;
                        bestMatchIndex = m;
                    }
                }
                
                if (bestMatchIndex >= 0 && minDeltaR < kDeltaRCut) 
                {
                    matchedPionsNeg++;
                    chargeMatchedPionsNeg++;
                    
                    h3D_TPC_RecoMatchedPions_N->Fill(pT, eta, truthVertexZ);
                    
                    double recoPt = negativeRecoTracks[bestMatchIndex]->getPt();
                    double recoEta = negativeRecoTracks[bestMatchIndex]->getEta();
                    h3D_TOF_RecoMatchedPions_N->Fill(recoPt, recoEta, truthVertexZ);
                    tofDenominatorCount_N++;
                    
                    if (negativeRecoTracks[bestMatchIndex]->getFlag(StUPCTrack::kTof)) {
                        h3D_TOF_RecoMatchedPionsWithTOF_N->Fill(recoPt, recoEta, truthVertexZ);
                        tofNumeratorCount_N++;
                    }
                    
                    //negativeRecoTracks.erase(negativeRecoTracks.begin() + bestMatchIndex);
                    // Instead, mark as used to prevent double-matching
                    negativeRecoTracks[bestMatchIndex] = nullptr;
                }
            }
        }



        // Fill event cut flow: Events with a reconstructed primary vertex
        HistEventCutFlow->Fill(1);

        // Fill event cut flow: >= 2 tracks with small DCA
        if (tracksWithSmallDca.size() < 2) 
        {
            continue;
        }
        HistEventCutFlow->Fill(2);         
        
        // Skip events with fewer than 3 TOF-matched tracks
		if (tracksWithTofHit.size() < 3) 
		{
			continue;
		}
        
        // Fill event cut flow: >= 3 TOF-matched tracks
        HistEventCutFlow->Fill(3);

         // Fill event cut flow: Final selection
         HistEventCutFlow->Fill(4);

         //Track flow
        for (Int_t j = 0; j < upcEvt->getNumberOfTracks(); j++)
        {
            StUPCTrack const* track = upcEvt->getTrack(j);

            // Track cut flow: All tracks
            HistTrackCutFlow->Fill(0);

            // Track cut flow: |eta| < 0.9
            if (abs(track->getEta()) < 0.9)
            {
                HistTrackCutFlow->Fill(1);
                // Track cut flow: pT > 0.2 GeV
                if (track->getPt() >= 0.25)
                {
                    HistTrackCutFlow->Fill(2);
                    // Track cut flow: NhitsFit >= 20
                    if (track->getNhitsFit() >= 20)
                    {
                        HistTrackCutFlow->Fill(3);
                        // Track cut flow: kV0 flag
                        if (track->getFlag(StUPCTrack::kV0))
                        {
                            HistTrackCutFlow->Fill(4);
                            // Track cut flow: TOF matched
                            if (track->getFlag(StUPCTrack::kTof))
                            {
                                HistTrackCutFlow->Fill(5);
                            }
                        }
                    }  
                }
            } 
        }
    }
    // At the end of your event loop
    cout << "Matching Statistics:" << endl;
    cout << "True π+: " << totalTruePionsPos << ", Matched: " << matchedPionsPos 
        << ", Charge matched: " << chargeMatchedPionsPos << endl;
    cout << "True π-: " << totalTruePionsNeg << ", Matched: " << matchedPionsNeg 
        << ", Charge matched: " << chargeMatchedPionsNeg << endl;
    
    // TOF matching stats for positive pions
    cout << "Positive TOF Denominator entries: " << tofDenominatorCount_P << endl;
    cout << "Positive TOF Numerator entries: " << tofNumeratorCount_P << endl;
    cout << "h3D_TOF_RecoMatchedPions_P integral: "  << h3D_TOF_RecoMatchedPions_P->Integral() << endl;
    cout << "h3D_TOF_RecoMatchedPionsWithTOF_P integral: "  << h3D_TOF_RecoMatchedPionsWithTOF_P->Integral() << endl;   
    
    // TOF matching stats for negative pions
    cout << "Negative TOF Denominator entries: " << tofDenominatorCount_N << endl;
    cout << "Negative TOF Numerator entries: " << tofNumeratorCount_N << endl;
    cout << "h3D_TOF_RecoMatchedPions_N integral: "  << h3D_TOF_RecoMatchedPions_N->Integral() << endl; 
    cout << "h3D_TOF_RecoMatchedPionsWithTOF_N integral: "  << h3D_TOF_RecoMatchedPionsWithTOF_N->Integral() << endl;

    // Write histograms to output file
    TFile *outfile = TFile::Open(argv[2], "recreate");

    HistEventCutFlow->Write();
    HistTrackCutFlow->Write();
    
    // =============================================================================
    // WRITE 3D EFFICIENCY HISTOGRAMS (SEPARATED BY CHARGE)
    // =============================================================================
    
    // Write 3D efficiency histograms for positive particles
    h3D_TPC_TruePions_P->Write();
    h3D_TPC_RecoMatchedPions_P->Write();
    h3D_TOF_RecoMatchedPions_P->Write();
    h3D_TOF_RecoMatchedPionsWithTOF_P->Write();
    
    // Write 3D efficiency histograms for negative particles
    h3D_TPC_TruePions_N->Write();
    h3D_TPC_RecoMatchedPions_N->Write();
    h3D_TOF_RecoMatchedPions_N->Write();
    h3D_TOF_RecoMatchedPionsWithTOF_N->Write();

    /*hSelPtMc_P->Write();
    hSelEtaMc_P->Write();
    hSelVzMc_P->Write();    
    hSelPtTrue_P->Write();
    hSelEtaTrue_P->Write();
    hSelVzTrue_P->Write();

    hSelPtReco_P->Write();
    hSelEtaReco_P->Write();
    hSelVzReco_P->Write();
    hSelPtReco_N->Write();  
    hSelEtaReco_N->Write();
    hSelVzReco_N->Write();
    
    hSelPtMc_N->Write(); 
    hSelEtaMc_N->Write();
    hSelVzMc_N->Write();
    hSelPtTrue_N->Write();  
    hSelEtaTrue_N->Write();
    hSelVzTrue_N->Write();

    hPtTpcTofMatched_P->Write();
    hEtaTpcTofMatched_P->Write();
    hVzTpcTofMatched_P->Write();
    hPtTpcTofMatched_N->Write();
    hEtaTpcTofMatched_N->Write();
    hVzTpcTofMatched_N->Write();*/
    
    cout << "Done" << endl;
  
    
    outfile->Close();

    return 0;
}
