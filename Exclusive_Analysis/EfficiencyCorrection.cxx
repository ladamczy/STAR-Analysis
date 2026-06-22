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

    // -------------------------------------------------------------
    // THE REALISTIC PHYSICS CAPS
    // -------------------------------------------------------------
    // 1. Cap TPC efficiency to prevent 1/0.001 = 1000x explosions
    if (eff_TPC > 0.90) eff_TPC = 0.90; 
    
    // 2. The most critical cap: Stop the TOF miss penalty from exploding.
    // If we cap at 0.85, the maximum miss penalty is 1 / (1 - 0.85) = 6.6
    if (eff_TOF_given_TPC > 0.85) eff_TOF_given_TPC = 0.85; 

    // Return 0 weight if in a true detector dead zone
    const double kMinEfficiency = 0.05; // Ignore tracks with < 5% efficiency
    if (eff_TPC <= kMinEfficiency) return 0.0;
    // -------------------------------------------------------------

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
    const int NFitHitsCut = 20; //17;//20 // 22;  check in ExclusiveCode.h

    CutConfig config("nominal");

    // ===================================================================
    // LOAD EFFICIENCY HISTOGRAMS
    // ===================================================================
    string effFilePath = "~/Downloads/RootFiles/SPK0K0StylePions_April20_0.root";//SPK0K0StylePions_April23_nHist22_0.root SPK0K0StylePions_April20_0.root//SPK0K0StylePions_April24_nHist17_0.root//SPK0K0StylePions_April20_0.root//SPK0K0StylePions_April16_0.root //SPK0K0StylePions_March13_0  SPK0K0StylePions_April16_0
    TFile* fEff = TFile::Open(effFilePath.c_str(), "READ");
    if (!fEff || fEff->IsZombie()) {
        cerr << "Error: Cannot open efficiency file " << effFilePath << endl;
        return 1;
    }

    TH3F* h3D_TPC_Eff_P = (TH3F*)((TH3F*)fEff->Get("h3D_TPC_RecoMatchedPions_P"))->Clone("h3D_TPC_Eff_P");
    h3D_TPC_Eff_P->Divide((TH3F*)fEff->Get("h3D_TPC_TruePions_P"));

    TH3F* h3D_TPC_Eff_N = (TH3F*)((TH3F*)fEff->Get("h3D_TPC_RecoMatchedPions_N"))->Clone("h3D_TPC_Eff_N");
    h3D_TPC_Eff_N->Divide((TH3F*)fEff->Get("h3D_TPC_TruePions_N"));

    TH3F* h3D_TOF_Eff_P = (TH3F*)((TH3F*)fEff->Get("h3D_TOF_RecoMatchedPionsWithTOF_P"))->Clone("h3D_TOF_Eff_P");
    h3D_TOF_Eff_P->Divide((TH3F*)fEff->Get("h3D_TOF_RecoMatchedPions_P"));

    TH3F* h3D_TOF_Eff_N = (TH3F*)((TH3F*)fEff->Get("h3D_TOF_RecoMatchedPionsWithTOF_N"))->Clone("h3D_TOF_Eff_N");
    h3D_TOF_Eff_N->Divide((TH3F*)fEff->Get("h3D_TOF_RecoMatchedPions_N"));
    /*
    TH3F* h3D_TPC_Eff_P = (TH3F*)((TH3F*)fEff->Get("h3D_TPC_RecoMatchedParticles_P"))->Clone("h3D_TPC_Eff_P");
    h3D_TPC_Eff_P->Divide((TH3F*)fEff->Get("h3D_TPC_TrueParticles_P"));
    TH3F* h3D_TPC_Eff_N = (TH3F*)((TH3F*)fEff->Get("h3D_TPC_RecoMatchedParticles_N"))->Clone("h3D_TPC_Eff_N");
    h3D_TPC_Eff_N->Divide((TH3F*)fEff->Get("h3D_TPC_TrueParticles_N"));
    TH3F* h3D_TOF_Eff_P = (TH3F*)((TH3F*)fEff->Get("h3D_TOF_RecoMatchedParticlesWithTOF_P"))->Clone("h3D_TOF_Eff_P");
    h3D_TOF_Eff_P->Divide((TH3F*)fEff->Get("h3D_TOF_RecoMatchedParticlesWithTOF_P"));
    TH3F* h3D_TOF_Eff_N = (TH3F*)((TH3F*)fEff->Get("h3D_TOF_RecoMatchedParticlesWithTOF_N"))->Clone("h3D_TOF_Eff_N");
    h3D_TOF_Eff_N->Divide((TH3F*)fEff->Get("h3D_TOF_RecoMatchedParticlesWithTOF_N"));
    */
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
    Double_t ptBins[] = {0.1, 0.15, 0.20, 0.25, 0.30, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6, 2.0, 3.0};
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
    TH1F *hPtMiss_Raw = new TH1F("hPtMiss_Raw", "Raw Missing p_{T};p_{T}^{miss} [GeV/c];Counts", 50, 0.0, 1.0);
    TH1F *hPtMiss_Corrected = new TH1F("hPtMiss_Corrected", "Corrected Missing p_{T};p_{T}^{miss} [GeV/c];Counts", 50, 0.0, 1.0);
    hPtMiss_Corrected->Sumw2();

    // ===================================================================
    // 4-PION SYSTEM HISTOGRAMS
    // ===================================================================
    //Double_t massBins4pi[] = {0.2, 0.4, 0.6, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.5, 3.8};
    //Int_t nMassBins4pi = sizeof(massBins4pi)/sizeof(massBins4pi[0]) - 1;
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
    TH1F *h1D_Reco_PtMiss_Raw_4pi_4TOF = new TH1F("h1D_Reco_PtMiss_Raw_4pi_4TOF", "Raw Data 4#pi Missing p_{T};p_{T}^{miss} (GeV/c);Counts", 50, 0.0, 1.0);
    TH1F *h1D_Reco_PtMiss_Corrected_4pi_4TOF = new TH1F("h1D_Reco_PtMiss_Corrected_4pi_4TOF", "Corrected Data 4#pi Missing p_{T};p_{T}^{miss} (GeV/c);Counts", 50, 0.0, 1.0);
    h1D_Reco_PtMiss_Corrected_4pi_4TOF->Sumw2();
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
    TH1F *h1D_Reco_PtMiss_Raw_4pi_3TOF = new TH1F("h1D_Reco_PtMiss_Raw_4pi_3TOF", "Raw Data 4#pi Missing p_{T};p_{T}^{miss} (GeV/c);Counts", 50, 0.0, 1.0);
    TH1F *h1D_Reco_PtMiss_Corrected_4pi_3TOF = new TH1F("h1D_Reco_PtMiss_Corrected_4pi_3TOF", "Corrected Data 4#pi Missing p_{T};p_{T}^{miss} (GeV/c);Counts", 50, 0.0, 1.0);
    h1D_Reco_PtMiss_Corrected_4pi_3TOF->Sumw2();
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
    TH1F *h1D_Reco_PtMiss_Raw_4pi_2TOF = new TH1F("h1D_Reco_PtMiss_Raw_4pi_2TOF", "Raw Data 4#pi Missing p_{T};p_{T}^{miss} (GeV/c);Counts", 50, 0.0, 1.0);
    TH1F *h1D_Reco_PtMiss_Corrected_4pi_2TOF = new TH1F("h1D_Reco_PtMiss_Corrected_4pi_2TOF", "Corrected Data 4#pi Missing p_{T};p_{T}^{miss} (GeV/c);Counts", 50, 0.0, 1.0);
    h1D_Reco_PtMiss_Corrected_4pi_2TOF->Sumw2();

    TH2D* h2_PtMiss_Vs_Mass = new TH2D("h2_PtMiss_Vs_Mass", "Missing p_{T} vs Mass;M(K_{S}^{0}K_{S}^{0}) [GeV/c^{2}];p_{T}^{miss} [GeV/c]", 21, 0.9, 3.0, 50, 0.0, 1.0);
    TH2D* h2_PtMiss_Vs_Mass_Corrected = new TH2D("h2_PtMiss_Vs_Mass_Corrected", "Corrected Missing p_{T} vs Mass;M(K_{S}^{0}K_{S}^{0}) [GeV/c^{2}];p_{T}^{miss} [GeV/c]", 21, 0.9, 3.0, 50, 0.0, 1.0);   
    h2_PtMiss_Vs_Mass_Corrected->Sumw2();

    TH2D* h2_PtMiss_Vs_Mass_2TOF = new TH2D("h2_PtMiss_Vs_Mass_2TOF", "Missing p_{T} vs Mass (2 TOF);M(K_{S}^{0}K_{S}^{0}) [GeV/c^{2}];p_{T}^{miss} [GeV/c]", 21, 0.9, 3.0, 50, 0.0, 1.0);
    TH2D* h2_PtMiss_Vs_Mass_2TOF_Corrected = new TH2D("h2_PtMiss_Vs_Mass_2TOF_Corrected", "Corrected Missing p_{T} vs Mass (2 TOF);M(K_{S}^{0}K_{S}^{0}) [GeV/c^{2}];p_{T}^{miss} [GeV/c]", 21, 0.9, 3.0, 50, 0.0, 1.0);
    h2_PtMiss_Vs_Mass_2TOF_Corrected->Sumw2();

    TH2D* h2_PtMiss_Vs_Mass_3TOF = new TH2D("h2_PtMiss_Vs_Mass_3TOF", "Missing p_{T} vs Mass (3 TOF);M(K_{S}^{0}K_{S}^{0}) [GeV/c^{2}];p_{T}^{miss} [GeV/c]", 21, 0.9, 3.0, 50, 0.0, 1.0);
    TH2D* h2_PtMiss_Vs_Mass_3TOF_Corrected = new TH2D("h2_PtMiss_Vs_Mass_3TOF_Corrected", "Corrected Missing p_{T} vs Mass (3 TOF);M(K_{S}^{0}K_{S}^{0}) [GeV/c^{2}];p_{T}^{miss} [GeV/c]", 21, 0.9, 3.0, 50, 0.0, 1.0);
    h2_PtMiss_Vs_Mass_3TOF_Corrected->Sumw2();

    TH2D* h2_PtMiss_Vs_Mass_4TOF = new TH2D("h2_PtMiss_Vs_Mass_4TOF", "Missing p_{T} vs Mass (4 TOF);M(K_{S}^{0}K_{S}^{0}) [GeV/c^{2}];p_{T}^{miss} [GeV/c]", 21, 0.9, 3.0, 50, 0.0, 1.0);
    TH2D* h2_PtMiss_Vs_Mass_4TOF_Corrected = new TH2D("h2_PtMiss_Vs_Mass_4TOF_Corrected", "Corrected Missing p_{T} vs Mass (4 TOF);M(K_{S}^{0}K_{S}^{0}) [GeV/c^{2}];p_{T}^{miss} [GeV/c]", 21, 0.9, 3.0, 50, 0.0, 1.0);
    h2_PtMiss_Vs_Mass_4TOF_Corrected->Sumw2();

    // 2-TOF LIKE-SIGN (LS) BACKGROUND HISTOGRAM
    TH1D* h1D_Reco_InvMass_Raw_LS_4pi_2TOF = new TH1D("h1D_Reco_InvMass_Raw_LS_4pi_2TOF", "Like-Sign Raw Background (2 TOF);M_{4#pi} (GeV/c^{2});Counts", 21, 0.9, 3.0); 
    TH1D* h1D_Reco_InvMass_Corrected_LS_4pi_2TOF = new TH1D("h1D_Reco_InvMass_Corrected_LS_4pi_2TOF", "Like-Sign Corrected Background (2 TOF);M_{4#pi} (GeV/c^{2});Counts", 21, 0.9, 3.0);
    h1D_Reco_InvMass_Corrected_LS_4pi_2TOF->Sumw2();

    TH2D* h2_PtMiss_Vs_Mass_Raw_4TOF_LS = new TH2D("h2_PtMiss_Vs_Mass_Raw_4TOF_LS", "LS Missing p_{T} vs Mass (4 TOF);M(K_{S}^{0}K_{S}^{0}) [GeV/c^{2}];p_{T}^{miss} [GeV/c]", 21, 0.9, 3.0, 50, 0.0, 1.0);
    TH2D* h2_PtMiss_Vs_Mass_Corrected_4TOF_LS = new TH2D("h2_PtMiss_Vs_Mass_Corrected_4TOF_LS", "LS Missing p_{T} vs Mass (4 TOF);M(K_{S}^{0}K_{S}^{0}) [GeV/c^{2}];p_{T}^{miss} [GeV/c]", 21, 0.9, 3.0, 50, 0.0, 1.0);
    h2_PtMiss_Vs_Mass_Corrected_4TOF_LS->Sumw2();
    TH2D* h2_PtMiss_Vs_Mass_Raw_3TOF_LS = new TH2D("h2_PtMiss_Vs_Mass_Raw_3TOF_LS", "LS Missing p_{T} vs Mass (3 TOF);M(K_{S}^{0}K_{S}^{0}) [GeV/c^{2}];p_{T}^{miss} [GeV/c]", 21, 0.9, 3.0, 50, 0.0, 1.0);
    TH2D* h2_PtMiss_Vs_Mass_Corrected_3TOF_LS = new TH2D("h2_PtMiss_Vs_Mass_Corrected_3TOF_LS", "LS Missing p_{T} vs Mass (3 TOF);M(K_{S}^{0}K_{S}^{0}) [GeV/c^{2}];p_{T}^{miss} [GeV/c]", 21, 0.9, 3.0, 50, 0.0, 1.0);
    h2_PtMiss_Vs_Mass_Corrected_3TOF_LS->Sumw2();
    TH2D* h2_PtMiss_Vs_Mass_Raw_2TOF_LS = new TH2D("h2_PtMiss_Vs_Mass_Raw_2TOF_LS", "LS Missing p_{T} vs Mass (2 TOF);M(K_{S}^{0}K_{S}^{0}) [GeV/c^{2}];p_{T}^{miss} [GeV/c]", 21, 0.9, 3.0, 50, 0.0, 1.0);
    TH2D* h2_PtMiss_Vs_Mass_Corrected_2TOF_LS = new TH2D("h2_PtMiss_Vs_Mass_Corrected_2TOF_LS", "LS Missing p_{T} vs Mass (2 TOF);M(K_{S}^{0}K_{S}^{0}) [GeV/c^{2}];p_{T}^{miss} [GeV/c]", 21, 0.9, 3.0, 50, 0.0, 1.0);
    h2_PtMiss_Vs_Mass_Corrected_2TOF_LS->Sumw2();

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

    double expectedYield_LS_2TOF = 0.0;
    int nEvents_LS_2TOF = 0;

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
                
                // STRICT K0K0 MAP CUTS: pT > 0.2 and diagonal fiducial cut
                if (track->getNhitsFit() >= NFitHitsCut && pt >= 0.2 && pt <= 3.0 && passFiducialCut(eta, eventVz) ) {// 
                    
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
        
    
        // START OF LIKE-SIGN (LS) BACKGROUND EVALUATION
        // 1. Map the tracks: [0] is positive, [1] is negative
        StUPCTrack const* pi_plus_1  = vPosNegPionLeadingKaon[0];
        StUPCTrack const* pi_minus_1 = vPosNegPionLeadingKaon[1];
        StUPCTrack const* pi_plus_2  = vPosNegPionSubLeadingKaon[0];
        StUPCTrack const* pi_minus_2 = vPosNegPionSubLeadingKaon[1];

        // 2. Build the Fake Like-Sign Kaons (++, --)
        StUPCV0 fakeKaonPlus(pi_plus_1, pi_plus_2, ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PION, 1, 1, tryVec, beamPar, upcEvt->getMagneticField(), true);
        StUPCV0 fakeKaonMinus(pi_minus_1, pi_minus_2, ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PION, 1, 1, tryVec, beamPar, upcEvt->getMagneticField(), true);

        // 3. Independent Boolean to track LS success (EXCLUDING pTmiss)
        bool passLSCuts = true;

        if (!(fakeKaonPlus.m() > config.massWinLow && fakeKaonPlus.m() < config.massWinHigh && 
              fakeKaonMinus.m() > config.massWinLow && fakeKaonMinus.m() < config.massWinHigh)) passLSCuts = false;

        double dcaBeamLeadingCutLS = (tracksWithTofHit.size() >= 3) ? config.dcaBeamline4 : config.dcaBeamline23;
        double dcaDauLeadingCutLS  = (tracksWithTofHit.size() >= 3) ? config.dcaDaughters4 : config.dcaDaughters23;
        double dcaBeamSubLeadingCutLS = (tracksWithTofHit.size() >= 4) ? config.dcaBeamline4 : config.dcaBeamline23;
        double dcaDauSubLeadingCutLS  = (tracksWithTofHit.size() >= 4) ? config.dcaDaughters4 : config.dcaDaughters23;

        if (fakeKaonPlus.dcaDaughters() > dcaDauLeadingCutLS || fakeKaonMinus.dcaDaughters() > dcaDauSubLeadingCutLS) passLSCuts = false;
        if (fakeKaonPlus.DCABeamLine() > dcaBeamLeadingCutLS || fakeKaonMinus.DCABeamLine() > dcaBeamSubLeadingCutLS) passLSCuts = false;

        bool lsPaPass = false;
        if (tracksWithTofHit.size() == 5 || tracksWithTofHit.size() == 4) {
            lsPaPass = (fakeKaonPlus.decayLengthHypo() <= config.decayLength || fakeKaonPlus.pointingAngleHypo() >= config.cosPA) &&
                       (fakeKaonMinus.decayLengthHypo() <= config.decayLength || fakeKaonMinus.pointingAngleHypo() >= config.cosPA);
        } else if (tracksWithTofHit.size() == 3) {
            lsPaPass = (fakeKaonPlus.decayLengthHypo() <= config.decayLength || fakeKaonPlus.pointingAngleHypo() >= config.cosPA) &&
                       fakeKaonMinus.pointingAngleHypo() >= config.cosPA23;
        } else if (tracksWithTofHit.size() == 2) {
            lsPaPass = (fakeKaonPlus.pointingAngleHypo() >= config.cosPA23 && fakeKaonMinus.pointingAngleHypo() >= config.cosPA23);
        }
        if (!lsPaPass) passLSCuts = false;

        // CALCULATE pTmiss BUT DO NOT CUT YET
        double pTmiss_LS;
        CheckPtMiss(fakeKaonPlus, fakeKaonMinus, protonE, protonW, pTmiss_LS);

        double zdiff_LS = fakeKaonPlus.decayVertex().Z() - fakeKaonMinus.decayVertex().Z();
        double zmean_LS = (fakeKaonPlus.decayVertex().Z() + fakeKaonMinus.decayVertex().Z()) / 2.0;
        if (abs(zmean_LS) > config.zmeanMax || abs(zdiff_LS) > config.zdiffMax) passLSCuts = false;

        if (abs(fakeKaonPlus.cosThetaStar()) > config.cosThetaStarMax || abs(fakeKaonMinus.cosThetaStar()) > config.cosThetaStarMax) passLSCuts = false;

        // 4. Fill 2D and 1D Histograms
        if (passLSCuts) {
            auto hasTofLS = [&](StUPCTrack const* trk) {
                return std::find(tracksWithTofHit.begin(), tracksWithTofHit.end(), trk) != tracksWithTofHit.end();
            };

            int dummyRejLS = 0;
            double w_LK_P_LS = GetK0TrackWeight(pi_plus_1, eventVz, h3D_TPC_Eff_P, h3D_TOF_Eff_P, hasTofLS(pi_plus_1), dummyRejLS);
            double w_LK_N_LS = GetK0TrackWeight(pi_minus_1, eventVz, h3D_TPC_Eff_N, h3D_TOF_Eff_N, hasTofLS(pi_minus_1), dummyRejLS);
            double w_SK_P_LS = GetK0TrackWeight(pi_plus_2, eventVz, h3D_TPC_Eff_P, h3D_TOF_Eff_P, hasTofLS(pi_plus_2), dummyRejLS);
            double w_SK_N_LS = GetK0TrackWeight(pi_minus_2, eventVz, h3D_TPC_Eff_N, h3D_TOF_Eff_N, hasTofLS(pi_minus_2), dummyRejLS);

            double lsEventWeight = w_LK_P_LS * w_LK_N_LS * w_SK_P_LS * w_SK_N_LS;

            if (lsEventWeight > 0.0) {
                int nTofLK_LS = hasTofLS(pi_plus_1) + hasTofLS(pi_minus_1);
                int nTofSK_LS = hasTofLS(pi_plus_2) + hasTofLS(pi_minus_2);
                int totalTOF_LS = nTofLK_LS + nTofSK_LS;
                TLorentzVector K0K0_LS = fakeKaonPlus.lorentzVector() + fakeKaonMinus.lorentzVector();
                double massK0K0_LS = K0K0_LS.M();

                // FILL 2D HISTOGRAMS (No pTmiss cut applied yet)
                if (totalTOF_LS == 4) {
                    h2_PtMiss_Vs_Mass_Corrected_4TOF_LS->Fill(massK0K0_LS, pTmiss_LS, lsEventWeight);
                    h2_PtMiss_Vs_Mass_Raw_4TOF_LS->Fill(massK0K0_LS, pTmiss_LS);
                } else if (totalTOF_LS == 3) {
                    double weight_3TOF_LS = lsEventWeight / 4.0;
                    h2_PtMiss_Vs_Mass_Corrected_3TOF_LS->Fill(massK0K0_LS, pTmiss_LS, weight_3TOF_LS);
                    h2_PtMiss_Vs_Mass_Raw_3TOF_LS->Fill(massK0K0_LS, pTmiss_LS);
                } else if (totalTOF_LS == 2 && nTofLK_LS == 1 && nTofSK_LS == 1) {
                    double weight_2TOF_LS = lsEventWeight / 4.0;
                    h2_PtMiss_Vs_Mass_Corrected_2TOF_LS->Fill(massK0K0_LS, pTmiss_LS, weight_2TOF_LS);
                    h2_PtMiss_Vs_Mass_Raw_2TOF_LS->Fill(massK0K0_LS, pTmiss_LS);

                    // FILL 1D HISTOGRAM AND TRACKERS (Apply pTmiss cut here!)
                    if (pTmiss_LS <= config.ptMissMax) {
                        h1D_Reco_InvMass_Raw_LS_4pi_2TOF->Fill(massK0K0_LS);
                        h1D_Reco_InvMass_Corrected_LS_4pi_2TOF->Fill(massK0K0_LS, weight_2TOF_LS);
                        expectedYield_LS_2TOF += weight_2TOF_LS;
                        nEvents_LS_2TOF++;
                    }
                }
            }
        }

        // END OF LIKE-SIGN (LS) BACKGROUND EVALUATION 

        // --- FULL EXCLUSIVITY CUTS (SC4 - SC9) ---
        double lkMass = leadingKaon.m();
        double skMass = subLeadingKaon.m();
        if (!(lkMass > config.massWinLow && lkMass < config.massWinHigh && skMass > config.massWinLow && skMass < config.massWinHigh)) continue;
        
        double pTmiss;
        //if (!CheckPtMiss(leadingKaon, subLeadingKaon, protonE, protonW, pTmiss) || pTmiss > config.ptMissMax) continue; 
        // Just calculate it. Ignore the boolean return and DO NOT 'continue' yet!
        CheckPtMiss(leadingKaon, subLeadingKaon, protonE, protonW, pTmiss);

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

        //int totalTOF = hasTof(vPosNegPionLeadingKaon[0]) + hasTof(vPosNegPionLeadingKaon[1]) + hasTof(vPosNegPionSubLeadingKaon[0]) + hasTof(vPosNegPionSubLeadingKaon[1]);
        
        int nTofLK = hasTof(vPosNegPionLeadingKaon[0]) + hasTof(vPosNegPionLeadingKaon[1]);
        int nTofSK = hasTof(vPosNegPionSubLeadingKaon[0]) + hasTof(vPosNegPionSubLeadingKaon[1]);
        int totalTOF = nTofLK + nTofSK;

        // Require exactly 1 TOF hit from the leading kaon AND 1 from the subleading kaon
        bool isSymmetric2TOF = (totalTOF == 2 && nTofLK == 1 && nTofSK == 1);

        //if (massK0K0 > 1.4 && massK0K0 < 1.9) { 
        hPtMiss_Raw->Fill(pTmiss);
        hPtMiss_Corrected->Fill(pTmiss, eventWeight); 
        if (totalTOF == 4) {
            h1D_Reco_PtMiss_Raw_4pi_4TOF->Fill(pTmiss);
            h1D_Reco_PtMiss_Corrected_4pi_4TOF->Fill(pTmiss, eventWeight);
            h2_PtMiss_Vs_Mass_4TOF->Fill(massK0K0, pTmiss);
            h2_PtMiss_Vs_Mass_4TOF_Corrected->Fill(massK0K0, pTmiss, eventWeight);
        } else if (totalTOF == 3) {
            double weight_3TOF = eventWeight / 4.0;
            h1D_Reco_PtMiss_Raw_4pi_3TOF->Fill(pTmiss);
            h1D_Reco_PtMiss_Corrected_4pi_3TOF->Fill(pTmiss, weight_3TOF);
            h2_PtMiss_Vs_Mass_3TOF->Fill(massK0K0, pTmiss);
            h2_PtMiss_Vs_Mass_3TOF_Corrected->Fill(massK0K0, pTmiss, weight_3TOF);
        } else if (totalTOF == 2) {
            double weight_2TOF = eventWeight / 4.0;
            h1D_Reco_PtMiss_Raw_4pi_2TOF->Fill(pTmiss);
            h1D_Reco_PtMiss_Corrected_4pi_2TOF->Fill(pTmiss, weight_2TOF);
            h2_PtMiss_Vs_Mass_2TOF->Fill(massK0K0, pTmiss);
            h2_PtMiss_Vs_Mass_2TOF_Corrected->Fill(massK0K0, pTmiss, weight_2TOF);  
        }
        //}
        h2_PtMiss_Vs_Mass->Fill(massK0K0, pTmiss);
        h2_PtMiss_Vs_Mass_Corrected->Fill(massK0K0, pTmiss, eventWeight);

        if (pTmiss > config.ptMissMax) continue;

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
        //else if (totalTOF == 2) {
        else if (isSymmetric2TOF) {
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
        cout << "Missing pT Mean:  " << h1D_Reco_PtMiss_Corrected_4pi_4TOF->GetMean() << " GeV/c" << endl;
    }

    cout << "\n=== 2-TOF SAMPLE DIAGNOSTICS (LIKE-SIGN METHOD) ===" << endl;
    cout << "Opposite-Sign (OS) Candidates (Signal + Background):" << endl;
    cout << "  Raw OS Events:                 " << nEvents_2TOF << endl;
    cout << "  Corrected OS Yield:            " << expectedYield_2TOF << endl;
    
    cout << "\nLike-Sign (LS) Candidates (Pure Combinatorial Background):" << endl;
    cout << "  Raw LS Events:                 " << nEvents_LS_2TOF << endl;
    cout << "  Corrected LS Yield:            " << expectedYield_LS_2TOF << endl;

    cout << "\nNet Extraction:" << endl;
    double netYield_2TOF = expectedYield_2TOF - expectedYield_LS_2TOF;
    double bgFraction_2TOF = (expectedYield_2TOF > 0) ? (expectedYield_LS_2TOF / expectedYield_2TOF) * 100.0 : 0.0;
    
    cout << "  Calculated Net Signal:         " << netYield_2TOF << endl;
    cout << "  Combinatorial Background %:    " << bgFraction_2TOF << " %" << endl;

    if (bgFraction_2TOF > 100.0) {
        cout << "  [!] WARNING: Background exceeds signal. Check LS scaling or acceptance." << endl;
    } else if (nEvents_LS_2TOF == 0 && nEvents_2TOF > 0) {
        cout << "  [!] WARNING: Zero LS background found. Check StUPCV0 charge logic." << endl;
    } else {
        cout << "  [+] SUCCESS: Like-Sign background successfully modeled." << endl;
    }

    outfile->cd();
    // Write 1D verifications
    h1D_RawPt_P->Write();  h1D_CorrectedPt_P->Write();
    h1D_RawEta_P->Write(); h1D_CorrectedEta_P->Write();
    h1D_RawVz_P->Write();  h1D_CorrectedVz_P->Write();
    h1D_RawPt_N->Write();  h1D_CorrectedPt_N->Write();
    h1D_RawEta_N->Write(); h1D_CorrectedEta_N->Write();
    h1D_RawVz_N->Write();  h1D_CorrectedVz_N->Write();
    hPtMiss_Raw->Write();
    hPtMiss_Corrected->Write();
    h2_PtMiss_Vs_Mass->Write();
    h2_PtMiss_Vs_Mass_Corrected->Write();
    h2_PtMiss_Vs_Mass_4TOF->Write();
    h2_PtMiss_Vs_Mass_4TOF_Corrected->Write();
    h2_PtMiss_Vs_Mass_3TOF->Write();
    h2_PtMiss_Vs_Mass_3TOF_Corrected->Write();
    h2_PtMiss_Vs_Mass_2TOF->Write();
    h2_PtMiss_Vs_Mass_2TOF_Corrected->Write();
    // Write Final 4p output
    h1D_Reco_InvMass_Raw_4pi_4TOF->Write();
    h1D_Reco_InvMass_Corrected_4pi_4TOF->Write();
    h1D_Reco_Pt_Raw_4pi_4TOF->Write();
    h1D_Reco_Pt_Corrected_4pi_4TOF->Write();
    h1D_Reco_Y_Raw_4pi_4TOF->Write();
    h1D_Reco_Y_Corrected_4pi_4TOF->Write();
    h1D_Reco_PtMiss_Raw_4pi_4TOF->Write();
    h1D_Reco_PtMiss_Corrected_4pi_4TOF->Write();
    // Write 3-TOF output
    h1D_Reco_InvMass_Raw_4pi_3TOF->Write();
    h1D_Reco_InvMass_Corrected_4pi_3TOF->Write();
    h1D_Reco_Pt_Raw_4pi_3TOF->Write();
    h1D_Reco_Pt_Corrected_4pi_3TOF->Write();
    h1D_Reco_Y_Raw_4pi_3TOF->Write();
    h1D_Reco_Y_Corrected_4pi_3TOF->Write();
    h1D_Reco_PtMiss_Raw_4pi_3TOF->Write();
    h1D_Reco_PtMiss_Corrected_4pi_3TOF->Write();
    // Write 2-TOF output
    h1D_Reco_InvMass_Raw_4pi_2TOF->Write();
    h1D_Reco_InvMass_Corrected_4pi_2TOF->Write();
    h1D_Reco_Pt_Raw_4pi_2TOF->Write();
    h1D_Reco_Pt_Corrected_4pi_2TOF->Write();
    h1D_Reco_Y_Raw_4pi_2TOF->Write();
    h1D_Reco_Y_Corrected_4pi_2TOF->Write();
    h1D_Reco_PtMiss_Raw_4pi_2TOF->Write();
    h1D_Reco_PtMiss_Corrected_4pi_2TOF->Write();
    //LS
    h1D_Reco_InvMass_Raw_LS_4pi_2TOF->Write();
    h1D_Reco_InvMass_Corrected_LS_4pi_2TOF->Write();
    h2_PtMiss_Vs_Mass_Raw_4TOF_LS->Write();
    h2_PtMiss_Vs_Mass_Corrected_4TOF_LS->Write();
    h2_PtMiss_Vs_Mass_Raw_3TOF_LS->Write();
    h2_PtMiss_Vs_Mass_Corrected_3TOF_LS->Write();
    h2_PtMiss_Vs_Mass_Raw_2TOF_LS->Write();
    h2_PtMiss_Vs_Mass_Corrected_2TOF_LS->Write();
    // Clean up
    delete d1; delete d2;
    for(auto d : dVec) delete d;
    
    fEff->Close();
    outfile->Close();
    cout << "\n✅ Done. Check the 1D & 4pi histograms in: " << argv[2] << endl;
    return 0;
}
