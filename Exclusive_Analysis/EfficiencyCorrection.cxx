#include "ExclusiveCode.h"
#include "SystematicCuts.h"
#include <TH3F.h>
#include <TH1F.h>
#include <TFile.h>
#include <algorithm>
#include <iostream>

using namespace std;

// ===================================================================
// EXACT WEIGHT LOGIC FROM K0K0_Analysis.C
// ===================================================================
double GetK0TrackWeight(StUPCTrack const* trk, double vz, TH3F* hTpcEff, TH3F* hTofEff, bool hasTof, int& rejCounter) {
    double pt = trk->getPt();
    double eta = trk->getEta();
    
    int binX = hTpcEff->GetXaxis()->FindBin(pt);
    int binY = hTpcEff->GetYaxis()->FindBin(eta);
    int binZ = hTpcEff->GetZaxis()->FindBin(vz);
    
    double eff_TPC = hTpcEff->GetBinContent(binX, binY, binZ);
    double eff_TOF_given_TPC = hTofEff->GetBinContent(binX, binY, binZ);

    // Apply safety caps
    if (eff_TPC > 0.98) eff_TPC = 0.98;
    if (eff_TOF_given_TPC > 0.98) eff_TOF_given_TPC = 0.98;

    // Return 0 weight if in a detector dead zone (instead of 'continue')
    if (eff_TPC < 0.01) return 0.0;
    
    const double kMinEfficiency = 0.01;
    
    if (eff_TPC <= kMinEfficiency) return 0.0;
    
    if (eff_TOF_given_TPC > 0.95) {
        rejCounter++;
        return 0.0; // Reject
    }

    double final_weight = 1.0;
    if (hasTof) {
        if (eff_TOF_given_TPC > kMinEfficiency) {
            final_weight = 1.0 / (eff_TPC * eff_TOF_given_TPC);
        } else {
            hasTof = false;
        }
    }
    
    if (!hasTof) {
        final_weight = 1.0 / (eff_TPC * (1.0 - eff_TOF_given_TPC));
    }
    
    return final_weight;
}

bool passFiducialCut(double eta, double zvtx) {
    return (eta < -zvtx/250.0 + 0.9) && 
           (eta > -zvtx/250.0 - 0.9) && 
           (fabs(eta) < 0.9);
}

int main(int argc, char** argv)  
{
    if (argc < 3) {
        cerr << "Usage: ./EfficiencyCorrection <input_list> <output.root>" << endl;
        return 1;
    }

    cout << "==============================================================" << endl;
    cout << "  UPC Data & MC Efficiency Correction - K0K0 Style" << endl;
    cout << "==============================================================" << endl;

    //some constants & config
    const double VzMax = 80.0; 

    CutConfig config("nominal");

    // ===================================================================
    // LOAD EFFICIENCY HISTOGRAMS
    // ===================================================================
    string effFilePath = "~/Downloads/SPK0K0StylePions_April16_0.root";//SPK0K0StylePions_March13_0
    TFile* fEff = TFile::Open(effFilePath.c_str(), "READ");
    if (!fEff || fEff->IsZombie()) {
        cerr << "Error: Cannot open efficiency file " << effFilePath << endl;
        return 1;
    }

    /*TH3F* h3D_TPC_Eff_P = (TH3F*)((TH3F*)fEff->Get("h3D_TPC_RecoMatchedPions_P"))->Clone("h3D_TPC_Eff_P");
    h3D_TPC_Eff_P->Divide((TH3F*)fEff->Get("h3D_TPC_TruePions_P"));

    TH3F* h3D_TPC_Eff_N = (TH3F*)((TH3F*)fEff->Get("h3D_TPC_RecoMatchedPions_N"))->Clone("h3D_TPC_Eff_N");
    h3D_TPC_Eff_N->Divide((TH3F*)fEff->Get("h3D_TPC_TruePions_N"));

    TH3F* h3D_TOF_Eff_P = (TH3F*)((TH3F*)fEff->Get("h3D_TOF_RecoMatchedPionsWithTOF_P"))->Clone("h3D_TOF_Eff_P");
    h3D_TOF_Eff_P->Divide((TH3F*)fEff->Get("h3D_TOF_RecoMatchedPions_P"));

    TH3F* h3D_TOF_Eff_N = (TH3F*)((TH3F*)fEff->Get("h3D_TOF_RecoMatchedPionsWithTOF_N"))->Clone("h3D_TOF_Eff_N");
    h3D_TOF_Eff_N->Divide((TH3F*)fEff->Get("h3D_TOF_RecoMatchedPions_N"));*/
   
    
    TH3F* h3D_TPC_Eff_P = (TH3F*)((TH3F*)fEff->Get("h3D_TPC_RecoMatchedParticles_P"))->Clone("h3D_TPC_Eff_P");
    h3D_TPC_Eff_P->Divide((TH3F*)fEff->Get("h3D_TPC_TrueParticles_P"));
    TH3F* h3D_TPC_Eff_N = (TH3F*)((TH3F*)fEff->Get("h3D_TPC_RecoMatchedParticles_N"))->Clone("h3D_TPC_Eff_N");
    h3D_TPC_Eff_N->Divide((TH3F*)fEff->Get("h3D_TPC_TrueParticles_N"));
    TH3F* h3D_TOF_Eff_P = (TH3F*)((TH3F*)fEff->Get("h3D_TOF_RecoMatchedParticlesWithTOF_P"))->Clone("h3D_TOF_Eff_P");
    h3D_TOF_Eff_P->Divide((TH3F*)fEff->Get("h3D_TOF_RecoMatchedParticlesWithTOF_P"));
    TH3F* h3D_TOF_Eff_N = (TH3F*)((TH3F*)fEff->Get("h3D_TOF_RecoMatchedParticlesWithTOF_N"))->Clone("h3D_TOF_Eff_N");
    h3D_TOF_Eff_N->Divide((TH3F*)fEff->Get("h3D_TOF_RecoMatchedParticlesWithTOF_N"));
    
    
    cout << "✅ Efficiency maps successfully loaded & divided." << endl;

    // Print efficiency comparison
    cout << "pi+ avg TPC efficiency: " << h3D_TPC_Eff_P->GetMean() << endl;
    cout << "pi- avg TPC efficiency: " << h3D_TPC_Eff_N->GetMean() << endl;
    cout << "pi+ avg TOF efficiency: " << h3D_TOF_Eff_P->GetMean() << endl;
    cout << "pi- avg TOF efficiency: " << h3D_TOF_Eff_N->GetMean() << endl;

    // ===================================================================
    // SETUP CHAIN
    // ===================================================================
    ifstream inputFilePathList(argv[1]);
    TChain *chain = new TChain("mUPCTree");
    string inputFileName;
    while (getline(inputFilePathList, inputFileName)) chain->Add(inputFileName.c_str());
    inputFilePathList.close();

    bool isMC = 0;
    static StUPCEvent *upcEvt = 0x0;
    static StRPEvent  *correctedRpEvent = 0x0;
    chain->SetBranchAddress("mUPCEvent", &upcEvt);
    chain->GetEntry(0);
    if (upcEvt->getRunNumber() == 1) isMC = 1;
    if (isMC == 0) chain->SetBranchAddress("correctedRpEvent", &correctedRpEvent);
    
    //long long nEntries = 100000;//chain->GetEntries();
    long long nEntries = chain->GetEntries();
    cout << "Processing " << nEntries << " events. Mode: " << (isMC ? "MC" : "Real Data") << endl;

    TFile *outfile = TFile::Open(argv[2], "recreate");

    // ===================================================================
    // 1D TRACK LEVEL VERIFICATION HISTOGRAMS (pT, eta, Vz)
    // ===================================================================
    Double_t ptBins[] = {0.20, 0.25, 0.30, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6, 2.0, 3.0};
    Int_t nPtBins = sizeof(ptBins)/sizeof(ptBins[0]) - 1;
    Double_t etaBins[] = {-1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    Int_t nEtaBins = sizeof(etaBins)/sizeof(etaBins[0]) - 1;
    Double_t vzBins[] = {-100, -90, -80, -70, -60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    Int_t nVzBins = sizeof(vzBins)/sizeof(vzBins[0]) - 1;

    // RAW 1D Distributions
    TH1F *h1D_RawPt_P  = new TH1F("h1D_RawPt_P",  "Raw p_{T} (#pi^{+});p_{T} (GeV/c);Counts", nPtBins, ptBins);
    TH1F *h1D_RawEta_P = new TH1F("h1D_RawEta_P", "Raw #eta (#pi^{+});#eta;Counts", nEtaBins, etaBins);
    TH1F *h1D_RawVz_P  = new TH1F("h1D_RawVz_P",  "Raw V_{z} (#pi^{+});V_{z} (cm);Counts", nVzBins, vzBins);
    TH1F *h1D_RawPt_N  = new TH1F("h1D_RawPt_N",  "Raw p_{T} (#pi^{-});p_{T} (GeV/c);Counts", nPtBins, ptBins);
    TH1F *h1D_RawEta_N = new TH1F("h1D_RawEta_N", "Raw #eta (#pi^{-});#eta;Counts", nEtaBins, etaBins);
    TH1F *h1D_RawVz_N  = new TH1F("h1D_RawVz_N",  "Raw V_{z} (#pi^{-});V_{z} (cm);Counts", nVzBins, vzBins);

    // CORRECTED 1D Distributions
    TH1F *h1D_CorrectedPt_P  = new TH1F("h1D_CorrectedPt_P",  "Corrected p_{T} (#pi^{+});p_{T} (GeV/c);Counts", nPtBins, ptBins); h1D_CorrectedPt_P->Sumw2();
    TH1F *h1D_CorrectedEta_P = new TH1F("h1D_CorrectedEta_P", "Corrected #eta (#pi^{+});#eta;Counts", nEtaBins, etaBins); h1D_CorrectedEta_P->Sumw2();
    TH1F *h1D_CorrectedVz_P  = new TH1F("h1D_CorrectedVz_P",  "Corrected V_{z} (#pi^{+});V_{z} (cm);Counts", nVzBins, vzBins); h1D_CorrectedVz_P->Sumw2();
    TH1F *h1D_CorrectedPt_N  = new TH1F("h1D_CorrectedPt_N",  "Corrected p_{T} (#pi^{-});p_{T} (GeV/c);Counts", nPtBins, ptBins); h1D_CorrectedPt_N->Sumw2();
    TH1F *h1D_CorrectedEta_N = new TH1F("h1D_CorrectedEta_N", "Corrected #eta (#pi^{-});#eta;Counts", nEtaBins, etaBins); h1D_CorrectedEta_N->Sumw2();
    TH1F *h1D_CorrectedVz_N  = new TH1F("h1D_CorrectedVz_N",  "Corrected V_{z} (#pi^{-});V_{z} (cm);Counts", nVzBins, vzBins); h1D_CorrectedVz_N->Sumw2();

    //Weighted Missing pT Histogram
    TH1F *hPtMiss_Corrected = new TH1F("hPtMiss_Corrected", "Corrected Missing p_{T};p_{T}^{miss} [GeV/c];Counts", 50, 0.0, 1.0);
    hPtMiss_Corrected->Sumw2();

    // ===================================================================
    // 4-PION SYSTEM HISTOGRAMS
    // ===================================================================
    Double_t massBins4pi[] = {0.3, 0.6, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.5, 3.8};
    Int_t nMassBins4pi = sizeof(massBins4pi)/sizeof(massBins4pi[0]) - 1;
    const int nPtBins4pi = 10;
    double ptBinEdges4pi[nPtBins4pi + 1] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0};
    // 4-TOF SAMPLE HISTOGRAMS
    TH1F *h1D_Reco_InvMass_Raw_4pi_4TOF = new TH1F("h1D_Reco_InvMass_Raw_4pi_4TOF", "Raw Data 4#pi Mass;M_{4#pi} (GeV/c^{2});Counts", 21, 0.9, 3.0);
    TH1F *h1D_Reco_InvMass_Corrected_4pi_4TOF = new TH1F("h1D_Reco_InvMass_Corrected_4pi_4TOF", "Corrected Data 4#pi Mass;M_{4#pi} (GeV/c^{2});Counts", 21, 0.9, 3.0);
    h1D_Reco_InvMass_Corrected_4pi_4TOF->Sumw2();
    TH1F *h1D_Reco_Pt_Raw_4pi_4TOF = new TH1F("h1D_Reco_Pt_Raw_4pi_4TOF", "Raw Data 4#pi p_{T};p_{T,K_{0}K_{0}} (GeV/c);Counts", nPtBins4pi, ptBinEdges4pi);
    TH1F *h1D_Reco_Pt_Corrected_4pi_4TOF = new TH1F("h1D_Reco_Pt_Corrected_4pi_4TOF", "Corrected Data 4#pi p_{T};p_{T,K_{0}K_{0}} (GeV/c);Counts", nPtBins4pi, ptBinEdges4pi);
    h1D_Reco_Pt_Corrected_4pi_4TOF->Sumw2(); 
    TH1F *h1D_Reco_Y_Raw_4pi_4TOF = new TH1F("h1D_Reco_Y_Raw_4pi_4TOF", "Raw Data 4#pi Rapidity;y_{K_{0}K_{0}};Counts", 40, -1.0, 1.0);
    TH1F *h1D_Reco_Y_Corrected_4pi_4TOF = new TH1F("h1D_Reco_Y_Corrected_4pi_4TOF", "Corrected Data 4#pi Rapidity;y_{K_{0}K_{0}};Counts", 40, -1.0, 1.0);
    h1D_Reco_Y_Corrected_4pi_4TOF->Sumw2();
    // 3-TOF SAMPLE HISTOGRAMS
    TH1F *h1D_Reco_InvMass_Raw_4pi_3TOF = new TH1F("h1D_Reco_InvMass_Raw_4pi_3TOF", "Raw Data 4#pi Mass (3 TOF);M_{4#pi} (GeV/c^{2});Counts", 21, 0.9, 3.0);
    TH1F *h1D_Reco_InvMass_Corrected_4pi_3TOF = new TH1F("h1D_Reco_InvMass_Corrected_4pi_3TOF", "Corrected Data 4#pi Mass (3 TOF);M_{4#pi} (GeV/c^{2});Counts", 21, 0.9, 3.0);
    h1D_Reco_InvMass_Corrected_4pi_3TOF->Sumw2();
    TH1F *h1D_Reco_Pt_Raw_4pi_3TOF = new TH1F("h1D_Reco_Pt_Raw_4pi_3TOF", "Raw Data 4#pi p_{T} (3 TOF);p_{T,4#pi} (GeV/c);Counts", nPtBins4pi, ptBinEdges4pi);
    TH1F *h1D_Reco_Pt_Corrected_4pi_3TOF = new TH1F("h1D_Reco_Pt_Corrected_4pi_3TOF", "Corrected Data 4#pi p_{T} (3 TOF);p_{T,4#pi} (GeV/c);Counts", nPtBins4pi, ptBinEdges4pi);
    h1D_Reco_Pt_Corrected_4pi_3TOF->Sumw2(); 
    TH1F *h1D_Reco_Y_Raw_4pi_3TOF = new TH1F("h1D_Reco_Y_Raw_4pi_3TOF", "Raw Data 4#pi Rapidity (3 TOF);y_{4#pi};Counts", 40, -1.0, 1.0);
    TH1F *h1D_Reco_Y_Corrected_4pi_3TOF = new TH1F("h1D_Reco_Y_Corrected_4pi_3TOF", "Corrected Data 4#pi Rapidity (3 TOF);y_{4#pi};Counts", 40, -1.0, 1.0);
    h1D_Reco_Y_Corrected_4pi_3TOF->Sumw2();
    // 2-TOF SAMPLE HISTOGRAMS
    TH1F *h1D_Reco_InvMass_Raw_4pi_2TOF = new TH1F("h1D_Reco_InvMass_Raw_4pi_2TOF", "Raw Data 4#pi Mass (2 TOF);M_{4#pi} (GeV/c^{2});Counts", 21, 0.9, 3.0);
    TH1F *h1D_Reco_InvMass_Corrected_4pi_2TOF = new TH1F("h1D_Reco_InvMass_Corrected_4pi_2TOF", "Corrected Data 4#pi Mass (2 TOF);M_{4#pi} (GeV/c^{2});Counts", 21, 0.9, 3.0);
    h1D_Reco_InvMass_Corrected_4pi_2TOF->Sumw2();
    TH1F *h1D_Reco_Pt_Raw_4pi_2TOF = new TH1F("h1D_Reco_Pt_Raw_4pi_2TOF", "Raw Data 4#pi p_{T} (2 TOF);p_{T,4#pi} (GeV/c);Counts", nPtBins4pi, ptBinEdges4pi);
    TH1F *h1D_Reco_Pt_Corrected_4pi_2TOF = new TH1F("h1D_Reco_Pt_Corrected_4pi_2TOF", "Corrected Data 4#pi p_{T} (2 TOF);p_{T,4#pi} (GeV/c);Counts", nPtBins4pi, ptBinEdges4pi);
    h1D_Reco_Pt_Corrected_4pi_2TOF->Sumw2(); 
    TH1F *h1D_Reco_Y_Raw_4pi_2TOF = new TH1F("h1D_Reco_Y_Raw_4pi_2TOF", "Raw Data 4#pi Rapidity (2 TOF);y_{4#pi};Counts", 40, -1.0, 1.0);
    TH1F *h1D_Reco_Y_Corrected_4pi_2TOF = new TH1F("h1D_Reco_Y_Corrected_4pi_2TOF", "Corrected Data 4#pi Rapidity (2 TOF);y_{4#pi};Counts", 40, -1.0, 1.0);
    h1D_Reco_Y_Corrected_4pi_2TOF->Sumw2();

    // -------------------------------------------------------------------
    TH1D* d1 = new TH1D("d1","",1,0,1); TH1D* d2 = new TH1D("d2","",1,0,1);
    vector<TH1D*> dVec;
    for(int k = 0; k < 4; k++) dVec.push_back(new TH1D(Form("dVec_%d",k),"",1,0,1));

    int count_4pi_Total = 0, count_4pi_4TOF = 0, count_4pi_3TOF = 0, count_4pi_2TOF_symmetric = 0, count_4pi_Other = 0;
    int nRejected_HighTofEff_P = 0, nRejected_HighTofEff_N = 0;

    // Counters for the final summary statement
    int nEvents_4TOF = 0;
    int nEvents_3TOF = 0;
    int nEvents_2TOF = 0;
    double expectedYield_4TOF = 0.0;
    double expectedYield_3TOF = 0.0;
    double expectedYield_2TOF = 0.0;
    int count_K0si_Other = 0;

    // ===================================================================
    // EVENT LOOP
    // ===================================================================
    for (Long64_t i = 0; i < nEntries; ++i)
    {
        chain->GetEntry(i);
        if (i % 10000 == 0 && i > 0) cout << "Processing event " << i << " (" << (100*i/nEntries) << "%)" << endl;

        // Safely extract vertex for Data vs MC
        if (isMC == 0) {
            if (upcEvt->getNumberOfVertices() != 1) continue;
            if (abs(upcEvt->getVertex(0)->getPosZ()) >= VzMax) continue;
        } else {
            if (upcEvt->getNumberOfVertices() < 1) continue; // MC fallback safety
        }
        double eventVz = upcEvt->getVertex(0)->getPosZ();

        // -------------------------------------------------------------------
        // 1D INCLUSIVE TRACK EFFICIENCY VERIFICATION (STRICT K0K0 CUTS)
        // -------------------------------------------------------------------
        if (abs(eventVz) <= VzMax) { 
            for (int trkIdx = 0; trkIdx < upcEvt->getNumberOfTracks(); trkIdx++) {
                StUPCTrack const* track = upcEvt->getTrack(trkIdx);
                if (!track) continue;
                
                double pt = track->getPt();
                double eta = track->getEta();
                
                // STRICT K0K0 MAP CUTS: pT > 0.25 and diagonal fiducial cut
                if (track->getNhitsFit() >= 20 && pt >= 0.2 && pt <= 3.0 && passFiducialCut(eta, eventVz) ) {// 
                    
                    bool hasTof = track->getFlag(StUPCTrack::kTof);
                    short charge = track->getCharge();
                    
                    if (charge > 0) {
                        double wP = GetK0TrackWeight(track, eventVz, h3D_TPC_Eff_P, h3D_TOF_Eff_P, hasTof, nRejected_HighTofEff_P);
                        if (wP > 0) {
                            h1D_RawPt_P->Fill(pt);        h1D_CorrectedPt_P->Fill(pt, wP);
                            h1D_RawEta_P->Fill(eta);      h1D_CorrectedEta_P->Fill(eta, wP);
                            h1D_RawVz_P->Fill(eventVz);   h1D_CorrectedVz_P->Fill(eventVz, wP);
                        }
                    } else if (charge < 0) {
                        double wN = GetK0TrackWeight(track, eventVz, h3D_TPC_Eff_N, h3D_TOF_Eff_N, hasTof, nRejected_HighTofEff_N);
                        if (wN > 0) {
                            h1D_RawPt_N->Fill(pt);        h1D_CorrectedPt_N->Fill(pt, wN);
                            h1D_RawEta_N->Fill(eta);      h1D_CorrectedEta_N->Fill(eta, wN);
                            h1D_RawVz_N->Fill(eventVz);   h1D_CorrectedVz_N->Fill(eventVz, wN);
                        }
                    }
                }
            }
        }
        // -------------------------------------------------------------------

        TLorentzVector protonE, protonW;
        FindProtons(isMC, correctedRpEvent, upcEvt, protonE, protonW);

        vector<StUPCTrack const*> tracksWithTofHit, tracksWithoutTofHit;
        SeparateTracks(upcEvt, tracksWithTofHit, tracksWithoutTofHit, isMC, d1, d2);

        vector<int> N_TOF_SELECTED = {2,3,4,5};
        if (!ValidNumberOfTofTracks(N_TOF_SELECTED, tracksWithTofHit, tracksWithoutTofHit, dVec, dVec, dVec, dVec, dVec, dVec)) continue;
        if (!AreTofTracksGood(tracksWithTofHit)) continue;

        vector<StUPCTrack const*> goodTracksWithoutTofHit;
        if (!FindTracks(tracksWithTofHit, tracksWithoutTofHit, goodTracksWithoutTofHit)) continue;

        double beamPar[4] = {};
        GetBeamPar(upcEvt, beamPar, isMC);

        vector<StUPCTrack const*> vPosNegPionLeadingKaon;
        vector<StUPCTrack const*> vPosNegPionSubLeadingKaon;
        if (!FindPions(upcEvt, tracksWithTofHit, tracksWithoutTofHit, goodTracksWithoutTofHit,
                       vPosNegPionLeadingKaon, vPosNegPionSubLeadingKaon, beamPar)) continue;

        TVector3 const tryVec(0,0,0);
        StUPCV0 leadingKaon(vPosNegPionLeadingKaon[0], vPosNegPionLeadingKaon[1], ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PION, 1, 1, tryVec, beamPar, upcEvt->getMagneticField(), true);
        StUPCV0 subLeadingKaon(vPosNegPionSubLeadingKaon[0], vPosNegPionSubLeadingKaon[1], ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PION, 1, 1, tryVec, beamPar, upcEvt->getMagneticField(), true);

        // --- FULL EXCLUSIVITY CUTS (SC4 - SC9) ---
        double lkMass = leadingKaon.m();
        double skMass = subLeadingKaon.m();
        if (!(lkMass > config.massWinLow && lkMass < config.massWinHigh && skMass > config.massWinLow && skMass < config.massWinHigh)) continue;
        
        double pTmiss;
        if (!CheckPtMiss(leadingKaon, subLeadingKaon, protonE, protonW, pTmiss) || pTmiss > config.ptMissMax) continue;

        int totalCluster = 0;
        CheckNumberOfClusters(upcEvt, tracksWithTofHit, totalCluster);
        if (totalCluster > (int)(2 * tracksWithTofHit.size() + 1)) continue;

        double dcaBeamLeadingCut = (tracksWithTofHit.size() >= 3) ? config.dcaBeamline4 : config.dcaBeamline23;
        double dcaDauLeadingCut  = (tracksWithTofHit.size() >= 3) ? config.dcaDaughters4 : config.dcaDaughters23;
        double dcaBeamSubLeadingCut = (tracksWithTofHit.size() >= 4) ? config.dcaBeamline4 : config.dcaBeamline23;
        double dcaDauSubLeadingCut  = (tracksWithTofHit.size() >= 4) ? config.dcaDaughters4 : config.dcaDaughters23;

        if (leadingKaon.dcaDaughters() > dcaDauLeadingCut || subLeadingKaon.dcaDaughters() > dcaDauSubLeadingCut) continue;
        if (leadingKaon.DCABeamLine() > dcaBeamLeadingCut || subLeadingKaon.DCABeamLine() > dcaBeamSubLeadingCut) continue;
        
        bool paPass = false;
        if (tracksWithTofHit.size() == 5 || tracksWithTofHit.size() == 4) {
            paPass = (leadingKaon.decayLengthHypo() <= config.decayLength || leadingKaon.pointingAngleHypo() >= config.cosPA) &&
                     (subLeadingKaon.decayLengthHypo() <= config.decayLength || subLeadingKaon.pointingAngleHypo() >= config.cosPA);
        } else if (tracksWithTofHit.size() == 3) {
            paPass = (leadingKaon.decayLengthHypo() <= config.decayLength || leadingKaon.pointingAngleHypo() >= config.cosPA) &&
                     subLeadingKaon.pointingAngleHypo() >= config.cosPA23;
        } else if (tracksWithTofHit.size() == 2) {
            paPass = (leadingKaon.pointingAngleHypo() >= config.cosPA23 && subLeadingKaon.pointingAngleHypo() >= config.cosPA23);
        }
        if (!paPass) continue;

        // -------------------------------------------------------------------
        // CORRECTED isMC LOGIC FOR ROMAN POTS (Matches EventAnalyzer.cxx)
        // -------------------------------------------------------------------
        double ksiE = 0, ksiW = 0;
        if (isMC == 0) {
            ksiE = correctedRpEvent->getTrack(0)->xi(ExclusiveK0K0::BEAM_ENERGY);
            ksiW = correctedRpEvent->getTrack(1)->xi(ExclusiveK0K0::BEAM_ENERGY);
        } else {
            ksiE = (255. - protonE.E()) / 255.;
            ksiW = (255. - protonW.E()) / 255.;
        }
        double sumProtonMomentumX = protonE.X() + protonW.X();
        double sumProtonMomentumY = protonE.Y() + protonW.Y();
        if (!(abs(ksiE) >= 0.007 || abs(ksiW) >= 0.007 || abs(sumProtonMomentumX) >= 0.1 || abs(sumProtonMomentumY) >= 0.1)) continue;

        TLorentzVector K0K0 = leadingKaon.lorentzVector() + subLeadingKaon.lorentzVector();
        double massK0K0 = K0K0.M();

        double ksiENew = (ksiE < 0.006) ? 0.003 : ksiE;
        double ksiWNew = (ksiW < 0.006) ? 0.003 : ksiW;
        if (abs(massK0K0/510.0 - sqrt(ksiENew*ksiWNew)) >= config.ksiCorrMax) continue;
        if (abs(K0K0.Rapidity() - 0.5*log(ksiENew/ksiWNew)) >= config.etaCorrMax) continue;
        if (abs(massK0K0/510.0 * exp(-K0K0.Rapidity()) - ksiE) >= config.separateCorrMax || 
            abs(massK0K0/510.0 * exp(+K0K0.Rapidity()) - ksiW) >= config.separateCorrMax) continue;

        double zdiff = leadingKaon.decayVertex().Z() - subLeadingKaon.decayVertex().Z();
        double zmean = (leadingKaon.decayVertex().Z() + subLeadingKaon.decayVertex().Z()) / 2.0;
        if (abs(zmean) > config.zmeanMax || abs(zdiff) > config.zdiffMax) continue;

        if (abs(leadingKaon.cosThetaStar()) > config.cosThetaStarMax || abs(subLeadingKaon.cosThetaStar()) > config.cosThetaStarMax) continue;

        // ===================================================================
        // PASSED ALL CUTS -> APPLY EFFICIENCY WEIGHTS FOR 4-PION EVENT
        // ===================================================================
        auto hasTof = [&](StUPCTrack const* trk) {
            return std::find(tracksWithTofHit.begin(), tracksWithTofHit.end(), trk) != tracksWithTofHit.end();
        };

        int dummyRej = 0;
        double w_LK_P = GetK0TrackWeight(vPosNegPionLeadingKaon[0], eventVz, h3D_TPC_Eff_P, h3D_TOF_Eff_P, hasTof(vPosNegPionLeadingKaon[0]), dummyRej);
        double w_LK_N = GetK0TrackWeight(vPosNegPionLeadingKaon[1], eventVz, h3D_TPC_Eff_N, h3D_TOF_Eff_N, hasTof(vPosNegPionLeadingKaon[1]), dummyRej);
        double w_SK_P = GetK0TrackWeight(vPosNegPionSubLeadingKaon[0], eventVz, h3D_TPC_Eff_P, h3D_TOF_Eff_P, hasTof(vPosNegPionSubLeadingKaon[0]), dummyRej);
        double w_SK_N = GetK0TrackWeight(vPosNegPionSubLeadingKaon[1], eventVz, h3D_TPC_Eff_N, h3D_TOF_Eff_N, hasTof(vPosNegPionSubLeadingKaon[1]), dummyRej);

        double eventWeight = w_LK_P * w_LK_N * w_SK_P * w_SK_N;
        if (eventWeight <= 0.0) continue; 

        int totalTOF = hasTof(vPosNegPionLeadingKaon[0]) + hasTof(vPosNegPionLeadingKaon[1]) + hasTof(vPosNegPionSubLeadingKaon[0]) + hasTof(vPosNegPionSubLeadingKaon[1]);

        count_4pi_Total++;
        // ===================================================================
        // FILL HISTOGRAMS BASED ON TOPOLOGY
        // ===================================================================
        if (totalTOF == 4) {
            nEvents_4TOF++;                   
            expectedYield_4TOF += eventWeight; 
            
            // Fill 4-TOF plots
            h1D_Reco_InvMass_Raw_4pi_4TOF->Fill(massK0K0);
            h1D_Reco_InvMass_Corrected_4pi_4TOF->Fill(massK0K0, eventWeight);
            h1D_Reco_Pt_Raw_4pi_4TOF->Fill(K0K0.Pt());
            h1D_Reco_Pt_Corrected_4pi_4TOF->Fill(K0K0.Pt(), eventWeight);
            h1D_Reco_Y_Raw_4pi_4TOF->Fill(K0K0.Rapidity());
            h1D_Reco_Y_Corrected_4pi_4TOF->Fill(K0K0.Rapidity(), eventWeight);
            
            // Missing pT (Usually kept for 4-TOF, but you can add it to others if needed)
            hPtMiss_Corrected->Fill(pTmiss, eventWeight);    
        }
        else if (totalTOF == 3) {
            nEvents_3TOF++;
            double weight_3TOF = eventWeight / 4.0;
            expectedYield_3TOF += eventWeight;
            // Fill 3-TOF plots
            h1D_Reco_InvMass_Raw_4pi_3TOF->Fill(massK0K0);
            h1D_Reco_InvMass_Corrected_4pi_3TOF->Fill(massK0K0, weight_3TOF);
            h1D_Reco_Pt_Raw_4pi_3TOF->Fill(K0K0.Pt());
            h1D_Reco_Pt_Corrected_4pi_3TOF->Fill(K0K0.Pt(), weight_3TOF);
            h1D_Reco_Y_Raw_4pi_3TOF->Fill(K0K0.Rapidity());
            h1D_Reco_Y_Corrected_4pi_3TOF->Fill(K0K0.Rapidity(), weight_3TOF);
        }
        else if (totalTOF == 2) {
            nEvents_2TOF++; 
            double weight_2TOF = eventWeight / 4.0;
            expectedYield_2TOF += eventWeight;
            // Fill 2-TOF plots
            h1D_Reco_InvMass_Raw_4pi_2TOF->Fill(massK0K0);
            h1D_Reco_InvMass_Corrected_4pi_2TOF->Fill(massK0K0, weight_2TOF);
            h1D_Reco_Pt_Raw_4pi_2TOF->Fill(K0K0.Pt());
            h1D_Reco_Pt_Corrected_4pi_2TOF->Fill(K0K0.Pt(), weight_2TOF);
            h1D_Reco_Y_Raw_4pi_2TOF->Fill(K0K0.Rapidity());
            h1D_Reco_Y_Corrected_4pi_2TOF->Fill(K0K0.Rapidity(), weight_2TOF);
        }
        else {
            count_K0si_Other++;
        }

    } // end event loop

    // ===================================================================
    // SUMMARY PRINTOUT
    // ===================================================================
    cout << "\n==============================================================" << endl;
    cout << " DATA/MC CORRECTION SUMMARY" << endl;
    cout << "==============================================================" << endl;

    cout << "\n=== 1D SINGLE TRACK INCLUSIVE CHECK ===" << endl;
    double rawTracksP = h1D_RawPt_P->Integral();
    double corrTracksP = h1D_CorrectedPt_P->Integral();
    double rawTracksN = h1D_RawPt_N->Integral();
    double corrTracksN = h1D_CorrectedPt_N->Integral();
    
    cout << "Total Raw Pi+ Tracks:       " << rawTracksP << endl;
    cout << "Total Corrected Pi+ Tracks: " << corrTracksP << endl;
    if (rawTracksP > 0) cout << ">>> AVERAGE Pi+ TRACK WEIGHT: " << corrTracksP / rawTracksP << " <<<" << endl;

    cout << "\nTotal Raw Pi- Tracks:       " << rawTracksN << endl;
    cout << "Total Corrected Pi- Tracks: " << corrTracksN << endl;
    if (rawTracksN > 0) cout << ">>> AVERAGE Pi- TRACK WEIGHT: " << corrTracksN / rawTracksN << " <<<" << endl;

    cout << "\n=== EVENT TOPOLOGY & YIELDS ===" << endl;
    cout << "4-TOF Sample:" << endl;
    cout << "  Raw Data Events found:         " << nEvents_4TOF << endl;
    cout << "  Corrected Expected Yield:      " << expectedYield_4TOF << endl;
    if (nEvents_4TOF > 0) {
        cout << "  -> Avg 4-TOF Event Weight:     " << expectedYield_4TOF / nEvents_4TOF << endl;
    }

    cout << "\n3-TOF Sample:" << endl;
    cout << "  Raw Data Events found:         " << nEvents_3TOF << endl;
    cout << "  Corrected Expected Yield:      " << expectedYield_3TOF << endl;
    if (nEvents_3TOF > 0) {
        cout << "  -> Avg 3-TOF Event Weight:     " << expectedYield_3TOF / nEvents_3TOF << endl;
    }
    
    cout << "\n2-TOF Sample:" << endl;
    cout << "  Raw Data Events found:         " << nEvents_2TOF << endl;
    cout << "  Corrected Expected Yield:      " << expectedYield_2TOF << endl;
    if (nEvents_2TOF > 0) {
        cout << "  -> Avg 2-TOF Event Weight:     " << expectedYield_2TOF / nEvents_2TOF << endl;
    }

    cout << "\n=== PHYSICS OBSERVABLES CHECK (4-TOF Base) ===" << endl;
    if (h1D_Reco_InvMass_Corrected_4pi_4TOF->GetEntries() > 0) {
        cout << "4-Pion Mass Mean: " << h1D_Reco_InvMass_Corrected_4pi_4TOF->GetMean() << " GeV/c^2" << endl;
        cout << "Missing pT Mean:  " << hPtMiss_Corrected->GetMean() << " GeV/c" << endl;
    }

    outfile->cd();
    // Write 1D verifications
    h1D_RawPt_P->Write();  h1D_CorrectedPt_P->Write();
    h1D_RawEta_P->Write(); h1D_CorrectedEta_P->Write();
    h1D_RawVz_P->Write();  h1D_CorrectedVz_P->Write();
    h1D_RawPt_N->Write();  h1D_CorrectedPt_N->Write();
    h1D_RawEta_N->Write(); h1D_CorrectedEta_N->Write();
    h1D_RawVz_N->Write();  h1D_CorrectedVz_N->Write();
    // Write Final 4p output
    h1D_Reco_InvMass_Raw_4pi_4TOF->Write();
    h1D_Reco_InvMass_Corrected_4pi_4TOF->Write();
    h1D_Reco_Pt_Raw_4pi_4TOF->Write();
    h1D_Reco_Pt_Corrected_4pi_4TOF->Write();
    h1D_Reco_Y_Raw_4pi_4TOF->Write();
    h1D_Reco_Y_Corrected_4pi_4TOF->Write();
    hPtMiss_Corrected->Write();
    // Write 3-TOF output
    h1D_Reco_InvMass_Raw_4pi_3TOF->Write();
    h1D_Reco_InvMass_Corrected_4pi_3TOF->Write();
    h1D_Reco_Pt_Raw_4pi_3TOF->Write();
    h1D_Reco_Pt_Corrected_4pi_3TOF->Write();
    h1D_Reco_Y_Raw_4pi_3TOF->Write();
    h1D_Reco_Y_Corrected_4pi_3TOF->Write();
    // Write 2-TOF output
    h1D_Reco_InvMass_Raw_4pi_2TOF->Write();
    h1D_Reco_InvMass_Corrected_4pi_2TOF->Write();
    h1D_Reco_Pt_Raw_4pi_2TOF->Write();
    h1D_Reco_Pt_Corrected_4pi_2TOF->Write();
    h1D_Reco_Y_Raw_4pi_2TOF->Write();
    h1D_Reco_Y_Corrected_4pi_2TOF->Write();

    // Clean up
    delete d1; delete d2;
    for(auto d : dVec) delete d;
    
    fEff->Close();
    outfile->Close();
    cout << "\n✅ Done. Check the 1D & 4pi histograms in: " << argv[2] << endl;
    return 0;
}