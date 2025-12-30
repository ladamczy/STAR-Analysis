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
#include <TVector3.h>
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
#include "StEfficiencyCorrector3D.h"

int getXiBin(double xi);
bool LambdaCut(const StUPCV0& L, char cut_type='0');
double getRapidity(const double pt, const double eta, const double phi, const double m);
using namespace std;

int main(int argc, char** argv)  
{
    //loading efficiency

    // TFile *effFile = TFile::Open("Efficiency/SingleParticleMC_ProtonsAug12_0.root", "READ");
    // if (!file || effFile->IsZombie()) {
    //     cerr << "ERROR: Cannot open ROOT file" << endl;
    //     return 1;
    // }
    // TH3F* h3D_TPC_num_plus = (TH3F*)effFile->Get("h3D_TPC_num_P");
    // TH3F* h3D_TPC_den_plus = (TH3F*)effFile->Get("h3D_TPC_den_P");
    // TH3F* h3D_TOF_num_plus = (TH3F*)effFile->Get("h3D_TOF_num_P");
    // TH3F* h3D_TOF_den_plus = (TH3F*)effFile->Get("h3D_TOF_den_P");

    // TH3F* h3D_TPC_num_minus = (TH3F*)effFile->Get("h3D_TPC_num_N");
    // TH3F* h3D_TPC_den_minus = (TH3F*)effFile->Get("h3D_TPC_den_N");
    // TH3F* h3D_TOF_num_minus = (TH3F*)effFile->Get("h3D_TOF_num_N");
    // TH3F* h3D_TOF_den_minus = (TH3F*)effFile->Get("h3D_TOF_den_N");

    // TH3F* tpcEfficiency_plus = (TH3F*)h3D_TPC_num_plus->Clone("tpcEfficiency_plus");
    // TH3F* tofEfficiency_plus = (TH3F*)h3D_TOF_num_plus->Clone("tofEfficiency_plus");

    // TH3F* tpcEfficiency_minus = (TH3F*)h3D_TPC_num_plus->Clone("tpcEfficiency_minus");
    // TH3F* tofEfficiency_minus = (TH3F*)h3D_TOF_num_plus->Clone("tofEfficiency_minus");

    // tpcEfficiency_plus->Divide(h3D_TPC_num_plus, h3D_TPC_den_plus, 1, 1, "B");
    // tofEfficiency_plus->Divide(h3D_TOF_num_plus, h3D_TOF_den_plus, 1, 1, "B"); //efficiency

    // tpcEfficiency_minus->Divide(h3D_TPC_num_minus, h3D_TPC_den_minus, 1, 1, "B");
    // tofEfficiency_minus->Divide(h3D_TOF_num_minus, h3D_TOF_den_minus, 1, 1, "B");

    // StEfficiencyCorrector3D corrector_plus;
    // StEfficiencyCorrector3D corrector_minus;

    // corrector_plus.setTpcEfficiency(tpcEfficiency_plus, +1, 2, true);
    // corrector_plus.setTofEfficiency(tofEfficiency_plus, +1, 2, true);

    // corrector_minus.setTpcEfficiency(tpcEfficiency_minus, -1, 2, true);
    // corrector_minus.setTofEfficiency(tofEfficiency_minus, -1, 2, true);

//    ifstream inputFilePathList(argv[1]);
//    if (!inputFilePathList) 
//
//    {
//        cerr << "Failed to open input file." << std::endl;
//        return 1;
//    }

    TChain *chain = new TChain("mUPCTree"); 

    string inputFileName;
    string outputFileName;

    std::string run = argv[1];
    inputFileName="/data1/sd_star_2017/presel_"+run+".root";      
    outputFileName="Run/"+run+".root";       

    chain->AddFile(inputFileName.c_str());

    static StUPCEvent * upcEvt = 0x0;
    static StRPEvent  * rpEvt = 0x0;

    chain->SetBranchAddress("mUPCEvent", &upcEvt);
    chain->SetBranchAddress("correctedRpEvent",  &rpEvt);


    TFile *outfile = TFile::Open(outputFileName.c_str(), "recreate"); 

    TDirectory* dirNSigma = outfile->mkdir("dEdx_Nsigma");
    TDirectory* dirProtons = outfile->mkdir("Protons");
    TDirectory* dirLambdas = outfile->mkdir("Lambdas");
    TDirectory* proton_lambda_distances = outfile->mkdir("proton_lambda_distances");

    outfile->cd();
    TH1D* HistNumOfPrimaryTracksToF = new TH1D("HistNumOfPrimaryTracksToF", "; Num of tracks; # events", 50 ,0, 50);
    TH1D* HistNumOfPToFWest = new TH1D("HistNumOfPToFWest", "; Num of tracks when proton on west; # events", 50 ,0, 50);
    TH1D* HistNumOfPToFEast = new TH1D("HistNumOfPToFEast", "; Num of tracks when proton on east; # events", 50 ,0, 50);
    TH1D* HistXiProtonWest = new TH1D("HistXiProtonWest", "xi of the proton on west;xi; # events", 100, -0.1, 1);
    TH1D* HistXiProtonEast = new TH1D("HistXiProtonEast", "xi of the proton on east;xi; # events", 100, -0.1, 1);
    TH1D* HistLogXiProtonWest = new TH1D("HistLogXiProtonWest", "logxi of the proton on west; logxi; # events", 100, -6, 0);
    TH1D* HistLogXiProtonEast = new TH1D("HistLogXiProtonEast", "logxi of the proton on east; logxi; # events", 100, -6, 0);
    TH1D* HistPtProtonWest = new TH1D("HistPtProtonWest", "pT of the proton on west; pT [GeV]; # events", 100, 0, 4);
    TH1D* HistPtProtonEast = new TH1D("HistPtProtonEast", "pT of the proton on west; pT [GeV]; pT [GeV]; # events", 100, 0, 4);
    TH1D* HistEtaProtonWest = new TH1D("HistEtaProtonWest", "Eta of the proton on west; #eta; # events", 100, -10, 10);
    TH1D* HistEtaProtonEast = new TH1D("HistEtaProtonEast", "Eta of the proton on east; #eta; # events", 100, -10, 10);
    TH1D* HistPtTracksWest = new TH1D("HistPtTracksWest", "pT of the reconstructed tracks when proton on west; pT [GeV]; # events", 100, 0, 4);
    TH1D* HistPtTracksEast = new TH1D("HistPtTracksEast", "pT of the reconstructed tracks when proton on east; pT [GeV]; # events", 100, 0, 4);
    TH1D* HistEtaTracksWest = new TH1D("HistEtaTracksWest", "Eta of the reconstructed tracks when proton on west; #eta; # events", 20, -1, 1);
    TH1D* HistEtaTracksEast = new TH1D("HistEtaTracksEast", "Eta of the reconstructed tracks when proton on east; #eta; # events", 20, -1, 1);
    TH1D* HistPtTracksWestCut = new TH1D("HistPtTracksWestCut", "pT of the reconstructed tracks when proton on west; pT [GeV]; # events", 100, 0, 4);
    TH1D* HistPtTracksEastCut = new TH1D("HistPtTracksEastCut", "pT of the reconstructed tracks when proton on east; pT [GeV]; # events", 100, 0, 4);
    TH1D* HistEtaTracksWestCut = new TH1D("HistEtaTracksWestCut", "Eta of the reconstructed tracks when proton on west; #eta; # events", 20, -1, 1);
    TH1D* HistEtaTracksEastCut = new TH1D("HistEtaTracksEastCut", "Eta of the reconstructed tracks when proton on east; #eta; # events", 20, -1, 1);
    
    TH2F* hNSigmaPiPlus = new TH2F("hNSigmaPiPlus",  ";p_{T} [GeV/c];n#sigma^{#pi^{+}}", 300, 0.0, 3.0, 200, -50.0, 50.0); //300 bins 0.0-3.0 for pT and 200 bins -50.0-50.0 for n_sigma
    TH2F* hNSigmaPiMinus = new TH2F("hNSigmaPiMinus", ";p_{T} [GeV/c];n#sigma^{#pi^{-}}", 300, 0.0, 3.0, 200, -50.0, 50.0);
    TH2F* hNSigmaKPlus = new TH2F("hNSigmaKPlus",   ";p_{T} [GeV/c];n#sigma^{K^{+}}", 300, 0.0, 3.0, 200, -50.0, 50.0);
    TH2F* hNSigmaKMinus = new TH2F("hNSigmaKMinus",  ";p_{T} [GeV/c];n#sigma^{K^{-}}", 300, 0.0, 3.0, 200, -50.0, 50.0);
    TH2F* hNSigmaPPlus = new TH2F("hNSigmaPPlus",   ";p_{T} [GeV/c];n#sigma^{p^{+}}", 300, 0.0, 3.0, 200, -50.0, 50.0);
    TH2F* hNSigmaPMinus = new TH2F("hNSigmaPMinus",  ";p_{T} [GeV/c];n#sigma^{p^{-}}", 300, 0.0, 3.0, 200, -50.0, 50.0);

    TH2F* hdEdx = new TH2F("hdEdx",  ";q #times p [GeV/c]; dE/dx [keV/cm]", 600, -5.0, 5.0, 600, 0, 100);


    // analysis of protons and antiprotons created

    TH1D* hPtProtonEastX = new TH1D("hPtProtonEastX","pT of protons (East); pT [GeV]; # events",60,0,3);
    TH1D* hPtProtonWestX = new TH1D("hPtProtonWestX","pT of protons (West); pT [GeV]; # events",60,0,3);
    TH1D* hPtProtonX = new TH1D("hPtProtonX","pT of protons; pT [GeV]; # events",60,0,3);

    TH1D* hPtNProtonEastX = new TH1D("hPtNProtonEastX","pT of antiprotons (East); pT [GeV]; # events",60,0,3);
    TH1D* hPtNProtonWestX = new TH1D("hPtNProtonWestX","pT of antiprotons (West); pT [GeV]; # events",60,0,3);
    TH1D* hPtNProtonX = new TH1D("hPtNProtonX","pT of antiprotons; pT [GeV]; # events",60,0,3);

    // eta (40 bins, −1–1)
    TH1D* hEtaProtonEastX = new TH1D("hEtaProtonEastX","#eta of protons (East); #eta; # events",20,-1,1);
    TH1D* hEtaProtonWestX = new TH1D("hEtaProtonWestX","#eta of protons (West); #eta; # events",20,-1,1);
    TH1D* hEtaProtonX = new TH1D("hEtaProtonX","#eta of proton; #eta; # events",20,-1,1);

    TH1D* hEtaNProtonEastX = new TH1D("hEtaNProtonEastX","#eta of antiprotons (East); #eta; # events",20,-1,1);
    TH1D* hEtaNProtonWestX = new TH1D("hEtaNProtonWestX","#eta of antiprotons (West); #eta; # events",20,-1,1);
    TH1D* hEtaNProtonX = new TH1D("hEtaNProtonX","#eta of antiprotons; #eta; # events",20,-1,1);

    TH1D* hRapidityProtonEastX  = new TH1D("hRapidityProtonEastX","#it{y} of protons (East); #it{y}; # events",20,-1,1);
    TH1D* hRapidityProtonWestX  = new TH1D("hRapidityProtonWestX","#it{y} of protons (West); #it{y}; # events",20,-1,1);
    TH1D* hRapidityProtonX      = new TH1D("hRapidityProtonX","#it{y} of protons; #it{y}; # events",20,-1,1);

    TH1D* hRapidityNProtonEastX = new TH1D("hRapidityNProtonEastX","#it{y} of antiprotons (East); #it{y}; # events",20,-1,1);
    TH1D* hRapidityNProtonWestX = new TH1D("hRapidityNProtonWestX","#it{y} of antiprotons (West); #it{y}; # events",20,-1,1);
    TH1D* hRapidityNProtonX     = new TH1D("hRapidityNProtonX","#it{y} of antiprotons; #it{y}; # events",20,-1,1);

    TH1D* hMultiplicityProtonEastX = new TH1D("hMultiplicityProtonEastX", "Multiplicity of protons (East); N_{p}; # events", 30, 0, 30);
    TH1D* hMultiplicityProtonWestX = new TH1D("hMultiplicityProtonWestX", "Multiplicity of protons (West); N_{p}; # events", 30, 0, 30);
    TH1D* hMultiplicityProtonX = new TH1D("hMultiplicityProtonX", "Multiplicity of protons; N_{p}; # events", 30, 0, 30);

    TH1D* hMultiplicityNProtonEastX = new TH1D("hMultiplicityNProtonEastX", "Multiplicity of antiprotons (East); N_{#bar{p}}; # events", 30, 0, 30);
    TH1D* hMultiplicityNProtonWestX = new TH1D("hMultiplicityNProtonWestX", "Multiplicity of antiprotons (West); N_{#bar{p}}; # events", 30, 0, 30);
    TH1D* hMultiplicityNProtonX = new TH1D("hMultiplicityNProtonX", "Multiplicity of antiprotons; N_{#bar{p}}; # events", 30, 0, 30);

    TH1D* hProton_DCA_E = new TH1D("hProton_DCA_E", "DCA of protons (East); DCA [cm]; #events", 100, 0, 40);
    TH1D* hProton_DCA_W = new TH1D("hProton_DCA_W", "DCA of protons (West); DCA [cm]; #events", 100, 0, 40);
    TH1D* hAntiproton_DCA_E = new TH1D("hAntiProton_DCA_E", "DCA of antiprotons (East); DCA [cm]; #events", 100, 0, 40);
    TH1D* hAntiproton_DCA_W = new TH1D("hAntiProton_DCA_W", "DCA of antiprotons (West); DCA [cm]; #events", 100, 0, 40);

    TH1D* hProton_DCA_E_pm = new TH1D("hProton_DCA_E_pm", "DCA of protons (East); DCA [cm]; #events", 100, 0, 40);
    TH1D* hProton_DCA_W_pm = new TH1D("hProton_DCA_W_pm", "DCA of protons (West); DCA [cm]; #events", 100, 0, 40);
    TH1D* hAntiproton_DCA_E_pm = new TH1D("hAntiProton_DCA_E_pm", "DCA of antiprotons (East); DCA [cm]; #events", 100, 0, 40);
    TH1D* hAntiproton_DCA_W_pm = new TH1D("hAntiProton_DCA_W_pm", "DCA of antiprotons (West); DCA [cm]; #events", 100, 0, 40);

    // TH1D* hProton_Primaries_E = new TH1D("hProton_Primaries_E", "Multiplicities of other primaries (East); N_{primaries}; #events", 20, 0, 20);
    // TH1D* hProton_Primaries_W = new TH1D("hProton_Primaries_W", "Multiplicities of other primaries (West); N_{primaries}; #events", 20, 0, 20);
    // TH1D* hAntiproton_Primaries_E = new TH1D("hAntiProton_Primaries_E", "Multiplicities of other primaries (East); N_{primaries}; #events", 20, 0, 20);
    // TH1D* hAntiproton_Primaries_W = new TH1D("hAntiProton_Primaries_W", "Multiplicities of other primaries (West); N_{primaries}; #events", 20, 0, 20);

    TH1D* hPAntiP_Primaries_E = new TH1D("hPAntiP_Primaries_E", "Multiplicities of other primaries (East); N_{primaries}; #events", 20, 0, 20);
    TH1D* hPAntiP_Primaries_W = new TH1D("hPAntiP_Primaries_W", "Multiplicities of other primaries (West); N_{primaries}; #events", 20, 0, 20);

    TH1D* hVz_check = new TH1D("Vz check", "Vz check", 200, -400, 400);

    TH3F* h3D_protons = new TH3F("h3D_protons",
                             "Protons: pT vs #eta vs Vz; p_{T} [GeV/c]; #eta; V_{z} [mm]",
                             60, 0.0, 3.0,      // pT bins
                             40, -1.0, 1.0,     // eta bins
                             80, -100.0, 100.0);  // Vz bins

    TH3F* h3D_antiprotons = new TH3F("h3D_antiprotons",
                            "Antiprotons: pT vs #eta vs Vz; p_{T} [GeV/c]; #eta; V_{z} [mm]",
                            60, 0.0, 3.0,      // pT bins
                            40, -1.0, 1.0,     // eta bins
                            80, -100.0, 100.0);  // Vz bins

    // making 3D histograms in terms of East and West side of forward proton -> eta changes
    TH3F* h3D_protonsEast = new TH3F("h3D_protonsEast",
                             "Protons east: pT vs #eta vs Vz; p_{T} [GeV/c]; #eta; V_{z} [mm]",
                             60, 0.0, 3.0,      // pT bins
                             40, -1.0, 1.0,     // eta bins
                             80, -100.0, 100.0);  // Vz bins

    TH3F* h3D_antiprotonsEast = new TH3F("h3D_antiprotonsEast",
                            "Antiprotons east: pT vs #eta vs Vz; p_{T} [GeV/c]; #eta; V_{z} [mm]",
                            60, 0.0, 3.0,      // pT bins
                            40, -1.0, 1.0,     // eta bins
                            80, -100.0, 100.0);  // Vz bins

    TH3F* h3D_protonsWest = new TH3F("h3D_protonsWest",
                             "Protons west: pT vs #eta vs Vz; p_{T} [GeV/c]; #eta; V_{z} [mm]",
                             60, 0.0, 3.0,      // pT bins
                             40, -1.0, 1.0,    // eta bins
                             80, -100.0, 100.0);  // Vz bins

    TH3F* h3D_antiprotonsWest = new TH3F("h3D_antiprotonsWest",
                            "Antiprotons west: pT vs #eta vs Vz; p_{T} [GeV/c]; #eta; V_{z} [mm]",
                            60, 0.0, 3.0,      // pT bins
                            40, -1.0, 1.0,    // eta bins
                            80, -100.0, 100.0);  // Vz bins


    TH1D* hEta_xi0_proton_E = new TH1D("hEta_xi0_proton_E", ";#eta;counts", 20, -1, 1);
    TH1D* hEta_xi1_proton_E = new TH1D("hEta_xi1_proton_E", ";#eta;counts", 20, -1, 1);
    TH1D* hEta_xi2_proton_E = new TH1D("hEta_xi2_proton_E", ";#eta;counts", 20, -1, 1);
    TH1D* hEta_xi3_proton_E = new TH1D("hEta_xi3_proton_E", ";#eta;counts", 20, -1, 1);
    TH1D* hEta_xi4_proton_E = new TH1D("hEta_xi4_proton_E", ";#eta;counts", 20, -1, 1);

    TH1D* hEta_xi0_proton_W = new TH1D("hEta_xi0_proton_W", ";#eta;counts", 20, -1, 1);
    TH1D* hEta_xi1_proton_W = new TH1D("hEta_xi1_proton_W", ";#eta;counts", 20, -1, 1);
    TH1D* hEta_xi2_proton_W = new TH1D("hEta_xi2_proton_W", ";#eta;counts", 20, -1, 1);
    TH1D* hEta_xi3_proton_W = new TH1D("hEta_xi3_proton_W", ";#eta;counts", 20, -1, 1);
    TH1D* hEta_xi4_proton_W = new TH1D("hEta_xi4_proton_W", ";#eta;counts", 20, -1, 1);

    TH1D* hEta_xi0_proton_C = new TH1D("hEta_xi0_proton_C", ";#eta;counts", 20, -1, 1);
    TH1D* hEta_xi1_proton_C = new TH1D("hEta_xi1_proton_C", ";#eta;counts", 20, -1, 1);
    TH1D* hEta_xi2_proton_C = new TH1D("hEta_xi2_proton_C", ";#eta;counts", 20, -1, 1);
    TH1D* hEta_xi3_proton_C = new TH1D("hEta_xi3_proton_C", ";#eta;counts", 20, -1, 1);
    TH1D* hEta_xi4_proton_C = new TH1D("hEta_xi4_proton_C", ";#eta;counts", 20, -1, 1);

    TH1D* hEta_xi0_antiproton_E = new TH1D("hEta_xi0_antiproton_E", ";#eta;counts", 20, -1, 1);
    TH1D* hEta_xi1_antiproton_E = new TH1D("hEta_xi1_antiproton_E", ";#eta;counts", 20, -1, 1);
    TH1D* hEta_xi2_antiproton_E = new TH1D("hEta_xi2_antiproton_E", ";#eta;counts", 20, -1, 1);
    TH1D* hEta_xi3_antiproton_E = new TH1D("hEta_xi3_antiproton_E", ";#eta;counts", 20, -1, 1);
    TH1D* hEta_xi4_antiproton_E = new TH1D("hEta_xi4_antiproton_E", ";#eta;counts", 20, -1, 1);

    TH1D* hEta_xi0_antiproton_W = new TH1D("hEta_xi0_antiproton_W", ";#eta;counts", 20, -1, 1);
    TH1D* hEta_xi1_antiproton_W = new TH1D("hEta_xi1_antiproton_W", ";#eta;counts", 20, -1, 1);
    TH1D* hEta_xi2_antiproton_W = new TH1D("hEta_xi2_antiproton_W", ";#eta;counts", 20, -1, 1);
    TH1D* hEta_xi3_antiproton_W = new TH1D("hEta_xi3_antiproton_W", ";#eta;counts", 20, -1, 1);
    TH1D* hEta_xi4_antiproton_W = new TH1D("hEta_xi4_antiproton_W", ";#eta;counts", 20, -1, 1);

    TH1D* hEta_xi0_antiproton_C = new TH1D("hEta_xi0_antiproton_C", ";#eta;counts", 20, -1, 1);
    TH1D* hEta_xi1_antiproton_C = new TH1D("hEta_xi1_antiproton_C", ";#eta;counts", 20, -1, 1);
    TH1D* hEta_xi2_antiproton_C = new TH1D("hEta_xi2_antiproton_C", ";#eta;counts", 20, -1, 1);
    TH1D* hEta_xi3_antiproton_C = new TH1D("hEta_xi3_antiproton_C", ";#eta;counts", 20, -1, 1);
    TH1D* hEta_xi4_antiproton_C = new TH1D("hEta_xi4_antiproton_C", ";#eta;counts", 20, -1, 1);



    TH1D* hRapidity_xi0_proton_E = new TH1D("hRapidity_xi0_proton_E", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_proton_E = new TH1D("hRapidity_xi1_proton_E", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_proton_E = new TH1D("hRapidity_xi2_proton_E", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_proton_E = new TH1D("hRapidity_xi3_proton_E", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_proton_E = new TH1D("hRapidity_xi4_proton_E", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_proton_W = new TH1D("hRapidity_xi0_proton_W", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_proton_W = new TH1D("hRapidity_xi1_proton_W", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_proton_W = new TH1D("hRapidity_xi2_proton_W", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_proton_W = new TH1D("hRapidity_xi3_proton_W", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_proton_W = new TH1D("hRapidity_xi4_proton_W", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_proton_C = new TH1D("hRapidity_xi0_proton_C", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_proton_C = new TH1D("hRapidity_xi1_proton_C", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_proton_C = new TH1D("hRapidity_xi2_proton_C", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_proton_C = new TH1D("hRapidity_xi3_proton_C", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_proton_C = new TH1D("hRapidity_xi4_proton_C", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_antiproton_E = new TH1D("hRapidity_xi0_antiproton_E", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_antiproton_E = new TH1D("hRapidity_xi1_antiproton_E", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_antiproton_E = new TH1D("hRapidity_xi2_antiproton_E", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_antiproton_E = new TH1D("hRapidity_xi3_antiproton_E", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_antiproton_E = new TH1D("hRapidity_xi4_antiproton_E", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_antiproton_W = new TH1D("hRapidity_xi0_antiproton_W", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_antiproton_W = new TH1D("hRapidity_xi1_antiproton_W", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_antiproton_W = new TH1D("hRapidity_xi2_antiproton_W", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_antiproton_W = new TH1D("hRapidity_xi3_antiproton_W", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_antiproton_W = new TH1D("hRapidity_xi4_antiproton_W", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_antiproton_C = new TH1D("hRapidity_xi0_antiproton_C", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_antiproton_C = new TH1D("hRapidity_xi1_antiproton_C", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_antiproton_C = new TH1D("hRapidity_xi2_antiproton_C", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_antiproton_C = new TH1D("hRapidity_xi3_antiproton_C", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_antiproton_C = new TH1D("hRapidity_xi4_antiproton_C", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_lambda_E = new TH1D("hRapidity_xi0_lambda_E", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_lambda_E = new TH1D("hRapidity_xi1_lambda_E", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_lambda_E = new TH1D("hRapidity_xi2_lambda_E", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_lambda_E = new TH1D("hRapidity_xi3_lambda_E", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_lambda_E = new TH1D("hRapidity_xi4_lambda_E", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_lambda_W = new TH1D("hRapidity_xi0_lambda_W", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_lambda_W = new TH1D("hRapidity_xi1_lambda_W", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_lambda_W = new TH1D("hRapidity_xi2_lambda_W", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_lambda_W = new TH1D("hRapidity_xi3_lambda_W", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_lambda_W = new TH1D("hRapidity_xi4_lambda_W", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_lambda_C = new TH1D("hRapidity_xi0_lambda_C", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_lambda_C = new TH1D("hRapidity_xi1_lambda_C", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_lambda_C = new TH1D("hRapidity_xi2_lambda_C", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_lambda_C = new TH1D("hRapidity_xi3_lambda_C", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_lambda_C = new TH1D("hRapidity_xi4_lambda_C", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_antilambda_E = new TH1D("hRapidity_xi0_antilambda_E", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_antilambda_E = new TH1D("hRapidity_xi1_antilambda_E", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_antilambda_E = new TH1D("hRapidity_xi2_antilambda_E", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_antilambda_E = new TH1D("hRapidity_xi3_antilambda_E", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_antilambda_E = new TH1D("hRapidity_xi4_antilambda_E", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_antilambda_W = new TH1D("hRapidity_xi0_antilambda_W", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_antilambda_W = new TH1D("hRapidity_xi1_antilambda_W", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_antilambda_W = new TH1D("hRapidity_xi2_antilambda_W", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_antilambda_W = new TH1D("hRapidity_xi3_antilambda_W", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_antilambda_W = new TH1D("hRapidity_xi4_antilambda_W", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_antilambda_C = new TH1D("hRapidity_xi0_antilambda_C", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_antilambda_C = new TH1D("hRapidity_xi1_antilambda_C", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_antilambda_C = new TH1D("hRapidity_xi2_antilambda_C", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_antilambda_C = new TH1D("hRapidity_xi3_antilambda_C", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_antilambda_C = new TH1D("hRapidity_xi4_antilambda_C", ";y;counts", 20, -1, 1);


    TH1D* hLambda_DCA_W=new TH1D("hLambda_DCA_W","#Lambda DCA daughters (West);DCA [cm];Events",100,0,40);
    TH1D* hLambda_DCA_E=new TH1D("hLambda_DCA_E","#Lambda DCA daughters (East);DCA [cm];Events",100,0,40);
    TH1D* hLambda_DCABeamLine_W=new TH1D("hLambda_DCABeamLine_W","#Lambda DCA to Beamline (West);DCA [cm];Events",100,0,40);
    TH1D* hLambda_DCABeamLine_E=new TH1D("hLambda_DCABeamLine_E","#Lambda DCA to Beamline (East);DCA [cm];Events",100,0,40);
    TH1D* hLambda_PointingAngle_W=new TH1D("hLambda_PointingAngle_W","#Lambda Pointing angle (West);cos(#theta);Events",100,-1.1,1.1);
    TH1D* hLambda_PointingAngle_E=new TH1D("hLambda_PointingAngle_E","#Lambda Pointing angle (East);cos(#theta);Events",100,-1.1,1.1);
    TH1D* hLambda_DecayLength_W=new TH1D("hLambda_DecayLength_W","#Lambda Decay length (West);L [cm];Events",100,0,100);
    TH1D* hLambda_DecayLength_E=new TH1D("hLambda_DecayLength_E","#Lambda Decay length (East);L [cm];Events",100,0,100);
    TH1D* hLambda_Mass_W=new TH1D("hLambda_Mass_W","#Lambda invariant mass (West);m [GeV/c^{2}];Events",100,1.05,1.25);
    TH1D* hLambda_Mass_E=new TH1D("hLambda_Mass_E","#Lambda invariant mass (East);m [GeV/c^{2}];Events",100,1.05,1.25);
    TH1D* hLambda_Mass_Background_Decay0Length_W=new TH1D("hLambda_Mass_Background_Decay0Length_W","#Lambda invariant mass background #lambda<3 (West);m [GeV/c^{2}];Events",100,1.05,1.25);
    TH1D* hLambda_Mass_Background_Decay0Length_E=new TH1D("hLambda_Mass_Background_Decay0Length_E","#Lambda invariant mass background #lambda<3 (East);m [GeV/c^{2}];Events",100,1.05,1.25);
    TH1D* hLambda_Mass_Background_Pointing0Angle_E=new TH1D("hLambda_Mass_Background_Pointing0Angle_E","#Lambda invariant mass background cos(#theta)<0.925 (East);m [GeV/c^{2}];Events",100,1.05,1.25);
    TH1D* hLambda_Mass_Background_Pointing0Angle_W=new TH1D("hLambda_Mass_Background_Pointing0Angle_W","#Lambda invariant mass background cos(#theta)<0.925 (West);m [GeV/c^{2}];Events",100,1.05,1.25);
    TH1D* hLambda_Mass_Background_DecayLength_W=new TH1D("hLambda_Mass_Background_DecayLength_W","#Lambda invariant mass background 3<DecayLength<5 (West);m [GeV/c^{2}];Events",100,1.05,1.25);
    TH1D* hLambda_Mass_Background_PointingAngle_W=new TH1D("hLambda_Mass_Background_PointingAngle_W","#Lambda invariant mass background 0.925<cos(#theta)<0.99 (West);m [GeV/c^{2}];Events",100,1.05,1.25);
    TH1D* hLambda_Mass_Background_DecayLength_E=new TH1D("hLambda_Mass_Background_DecayLength_E","#Lambda invariant mass background 3<DecayLength<5 (East);m [GeV/c^{2}];Events",100,1.05,1.25);
    TH1D* hLambda_Mass_Background_PointingAngle_E=new TH1D("hLambda_Mass_Background_PointingAngle_E","#Lambda invariant mass background 0.925<cos(#theta)<0.99 (East);m [GeV/c^{2}];Events",100,1.05,1.25);
    TH1D* hLambda_Mass_Background_Neg_PointingAngle_E=new TH1D("hLambda_Mass_Background_Neg_PointingAngle_E","#Lambda invariant mass background cos(#theta)<-0.925 (East);m [GeV/c^{2}];Events",100,1.05,1.25);
    TH1D* hLambda_Mass_Background_Neg_PointingAngle_W=new TH1D("hLambda_Mass_Background_Neg_PointingAngle_W","#Lambda invariant mass background cos(#theta)<-0.925 (West);m [GeV/c^{2}];Events",100,1.05,1.25);
    TH1D* hPT_Lambda_E = new TH1D("pT_lambda_E","p_{T} of #Lambda (East);p_{T} [GeV/c];Counts",60,0,3);
    TH1D* hPT_Lambda_W = new TH1D("pT_lambda_W","p_{T} of #Lambda (West);p_{T} [GeV/c];Counts",60,0,3);
    TH1D* hEta_Lambda_E = new TH1D("eta_lambda_E","#eta of #Lambda (East);#eta;Counts",20,-1,1);
    TH1D* hEta_Lambda_W = new TH1D("eta_lambda_W","#eta of #Lambda (West);#eta;Counts",20,-1,1);
    TH1D* hRapidity_Lambda_E = new TH1D("rapidity_lambda_E","y of #Lambda (East);y;Counts",20,-1,1);
    TH1D* hRapidity_Lambda_W = new TH1D("rapidity_lambda_W","y of #Lambda (West);y;Counts",20,-1,1);

    TH1D* hAntiLambda_DCA_W=new TH1D("hAntiLambda_DCA_W","#bar{#Lambda} DCA daughters (West);DCA [cm];Events",100,0,40);
    TH1D* hAntiLambda_DCA_E=new TH1D("hAntiLambda_DCA_E","#bar{#Lambda} DCA daughters (East);DCA [cm];Events",100,0,40);
    TH1D* hAntiLambda_DCABeamLine_W=new TH1D("hAntiLambda_DCABeamLine_W","#bar{#Lambda} DCA to Beamline (West);DCA [cm];Events",100,0,40);
    TH1D* hAntiLambda_DCABeamLine_E=new TH1D("hAntiLambda_DCABeamLine_E","#bar{#Lambda} DCA to Beamline (East);DCA [cm];Events",100,0,40);
    TH1D* hAntiLambda_PointingAngle_W=new TH1D("hAntiLambda_PointingAngle_W","#bar{#Lambda} Pointing angle (West);cos(#theta);Events",100,-1.1,1.1);
    TH1D* hAntiLambda_PointingAngle_E=new TH1D("hAntiLambda_PointingAngle_E","#bar{#Lambda} Pointing angle (East);cos(#theta);Events",100,-1.1,1.1);
    TH1D* hAntiLambda_DecayLength_W=new TH1D("hAntiLambda_DecayLength_W","#bar{#Lambda} Decay length (West);L [cm];Events",100,0,100);
    TH1D* hAntiLambda_DecayLength_E=new TH1D("hAntiLambda_DecayLength_E","#bar{#Lambda} Decay length (East);L [cm];Events",100,0,100);
    TH1D* hAntiLambda_Mass_W=new TH1D("hAntiLambda_Mass_W","#bar{#Lambda} invariant mass (West);m [GeV/c^{2}];Events",100,1.05,1.25);
    TH1D* hAntiLambda_Mass_E=new TH1D("hAntiLambda_Mass_E","#bar{#Lambda} invariant mass (East);m [GeV/c^{2}];Events",100,1.05,1.25);
    TH1D* hAntiLambda_Mass_Background_Decay0Length_W=new TH1D("hAntiLambda_Mass_Background_Decay0Length_W","#bar{#Lambda} invariant mass background #lambda<3 (West);m [GeV/c^{2}];Events",100,1.05,1.25);
    TH1D* hAntiLambda_Mass_Background_Decay0Length_E=new TH1D("hAntiLambda_Mass_Background_Decay0Length_E","#bar{#Lambda} invariant mass background #lambda<3 (East);m [GeV/c^{2}];Events",100,1.05,1.25);
    TH1D* hAntiLambda_Mass_Background_Pointing0Angle_E=new TH1D("hAntiLambda_Mass_Background_Pointing0Angle_E","#bar{#Lambda} invariant mass background cos(#theta)<0.925 (East);m [GeV/c^{2}];Events",100,1.05,1.25);
    TH1D* hAntiLambda_Mass_Background_Pointing0Angle_W=new TH1D("hAntiLambda_Mass_Background_Pointing0Angle_W","#bar{#Lambda} invariant mass background cos(#theta)<0.925 (West);m [GeV/c^{2}];Events",100,1.05,1.25);
    TH1D* hAntiLambda_Mass_Background_DecayLength_W=new TH1D("hAntiLambda_Mass_Background_DecayLength_W","#bar{#Lambda} invariant mass background 3<DecayLength<5 (West);m [GeV/c^{2}];Events",100,1.05,1.25);
    TH1D* hAntiLambda_Mass_Background_PointingAngle_W=new TH1D("hAntiLambda_Mass_Background_PointingAngle_W","#bar{#Lambda} invariant mass background 0.925<cos(#theta)<0.99 (West);m [GeV/c^{2}];Events",100,1.05,1.25);
    TH1D* hAntiLambda_Mass_Background_DecayLength_E=new TH1D("hAntiLambda_Mass_Background_DecayLength_E","#bar{#Lambda} invariant mass background 3<DecayLength<5 (East);m [GeV/c^{2}];Events",100,1.05,1.25);
    TH1D* hAntiLambda_Mass_Background_PointingAngle_E=new TH1D("hAntiLambda_Mass_Background_PointingAngle_E","#bar{#Lambda} invariant mass background 0.925<cos(#theta)<0.99 (East);m [GeV/c^{2}];Events",100,1.05,1.25);
    TH1D* hAntiLambda_Mass_Background_Neg_PointingAngle_E=new TH1D("hAntiLambda_Mass_Background_Neg_PointingAngle_E","#bar{#Lambda} invariant mass background cos(#theta)<-0.925 (East);m [GeV/c^{2}];Events",100,1.05,1.25);
    TH1D* hAntiLambda_Mass_Background_Neg_PointingAngle_W=new TH1D("hAntiLambda_Mass_Background_Neg_PointingAngle_W","#bar{#Lambda} invariant mass background cos(#theta)<-0.925 (West);m [GeV/c^{2}];Events",100,1.05,1.25);
    TH1D* hPT_antiLambda_E = new TH1D("pT_antiLambda_E","p_{T} of #bar{#Lambda} (East);p_{T} [GeV/c];Counts",60,0,3);
    TH1D* hPT_antiLambda_W = new TH1D("pT_antiLambda_W","p_{T} of #bar{#Lambda} (West);p_{T} [GeV/c];Counts",60,0,3);
    TH1D* hEta_antiLambda_E = new TH1D("eta_antiLambda_E","#eta of #bar{#Lambda} (East);#eta;Counts",20,-1,1);
    TH1D* hEta_antiLambda_W = new TH1D("eta_antiLambda_W","#eta of #bar{#Lambda} (West);#eta;Counts",20,-1,1);
    TH1D* hRapidity_antiLambda_E = new TH1D("rapidity_antiLambda_E","y of #bar{#Lambda} (East);y;Counts",20,-1,1);
    TH1D* hRapidity_antiLambda_W = new TH1D("rapidity_antiLambda_W","y of #bar{#Lambda} (West);y;Counts",20,-1,1);

    TH1D* hDistanceBeam_proton=new TH1D("hDistanceBeam_proton","Distance from beam for protons;#it{y};# events",160,-8,8);
    TH1D* hDistanceBeam_antiproton=new TH1D("hDistanceBeam_antiproton","Distance from beam for antiprotons;#it{y};# events",160,-8,8);
    TH1D* hDistanceBeam_Lambda=new TH1D("hDistanceBeam_Lambda","Distance from beam for #Lambda;#it{y};# events",160,-8,8);
    TH1D* hDistanceBeam_antiLambda=new TH1D("hDistanceBeam_antiLambda","Distance from beam for #bar{#Lambda};#it{y};# events",160,-8,8);

    TH1D* hDistanceEdge_proton=new TH1D("hDistanceEdge_proton","Distance to edge for protons;#it{y};# events",80,-8,8);
    TH1D* hDistanceEdge_antiproton=new TH1D("hDistanceEdge_antiproton","Distance to edge for antiprotons;#it{y};# events",80,-8,8);
    TH1D* hDistanceEdge_Lambda=new TH1D("hDistanceEdge_Lambda","Distance to edge for #Lambda;#it{y};# events",80,-8,8);
    TH1D* hDistanceEdge_antiLambda=new TH1D("hDistanceEdge_antiLambda","Distance to edge for #bar{#Lambda};#it{y};# events",80,-8,8);


    TH1D* hDistanceBeam_proton_E=new TH1D("hDistanceBeam_proton_E","Distance from beam for protons East;#it{y};# events",160,-8,8);
    TH1D* hDistanceBeam_proton_W=new TH1D("hDistanceBeam_proton_W","Distance from beam for protons West;#it{y};# events",160,-8,8);
    TH1D* hDistanceBeam_antiproton_E=new TH1D("hDistanceBeam_antiproton_E","Distance from beam for antiprotons East;#it{y};# events",160,-8,8);
    TH1D* hDistanceBeam_antiproton_W=new TH1D("hDistanceBeam_antiproton_W","Distance from beam for antiprotons West;#it{y};# events",160,-8,8);
    TH1D* hDistanceBeam_Lambda_E=new TH1D("hDistanceBeam_Lambda_E","Distance from beam for #Lambda East;#it{y};# events",160,-8,8);
    TH1D* hDistanceBeam_Lambda_W=new TH1D("hDistanceBeam_Lambda_W","Distance from beam for #Lambda West;#it{y};# events",160,-8,8);
    TH1D* hDistanceBeam_antiLambda_E=new TH1D("hDistanceBeam_antiLambda_E","Distance from beam for #bar{#Lambda} East;#it{y};# events",160,-8,8);
    TH1D* hDistanceBeam_antiLambda_W=new TH1D("hDistanceBeam_antiLambda_W","Distance from beam for #bar{#Lambda} West;#it{y};# events",160,-8,8);

    TH1D* hDistanceEdge_proton_E=new TH1D("hDistanceEdge_proton_E","Distance to edge for protons East;#it{y};# events",80,-8,8);
    TH1D* hDistanceEdge_proton_W=new TH1D("hDistanceEdge_proton_W","Distance to edge for protons West;#it{y};# events",80,-8,8);
    TH1D* hDistanceEdge_antiproton_E=new TH1D("hDistanceEdge_antiproton_E","Distance to edge for antiprotons East;#it{y};# events",80,-8,8);
    TH1D* hDistanceEdge_antiproton_W=new TH1D("hDistanceEdge_antiproton_W","Distance to edge for antiprotons West;#it{y};# events",80,-8,8);
    TH1D* hDistanceEdge_Lambda_E=new TH1D("hDistanceEdge_Lambda_E","Distance to edge for #Lambda East;#it{y};# events",80,-8,8);
    TH1D* hDistanceEdge_Lambda_W=new TH1D("hDistanceEdge_Lambda_W","Distance to edge for #Lambda West;#it{y};# events",80,-8,8);
    TH1D* hDistanceEdge_antiLambda_E=new TH1D("hDistanceEdge_antiLambda_E","Distance to edge for #bar{#Lambda} East;#it{y};# events",80,-8,8);
    TH1D* hDistanceEdge_antiLambda_W=new TH1D("hDistanceEdge_antiLambda_W","Distance to edge for #bar{#Lambda} West;#it{y};# events",80,-8,8);

    double massPion = 0.13957061;
    double massProton = 0.93827;
    double massLambda = 1.115638;
    

    for (Long64_t i = 0; i < chain->GetEntries(); ++i) {
        chain->GetEntry(i);
        if( upcEvt->getNumberOfVertices() != 1) continue;
        if (rpEvt->getNumberOfTracks() != 1 ) continue;

        //   cout << upcEvt->getEventNumber() << " " 
        //   cout << upcEvt->getNPrimVertices() <<  endl; //1, 2, 3

        StUPCRpsTrack *proton = rpEvt->getTrack(0);

        bool isEast = (proton->branch() < 2);
        bool isWest = (proton->branch() > 1);

        double xi_proton = proton->xi(254.867);
        int b=getXiBin(xi_proton);

        if (xi_proton<1e-7) continue;

        double y_beam = 0.0;
        double y_edge = 0.0;

        if (isWest) { //positive rapidity
            y_beam = -6.30;
            y_edge = 6.30 + log(xi_proton);
        } 
        else { //negative rapidity
            y_beam = 6.30;
            y_edge = -6.30 - log(xi_proton);
        }

        int NumOfPrimaryTracksToF = 0;

        int MultiplicityProtonX = 0;
        int MultiplicityNProtonX=0;
        int MultiplicityProtonEastX=0;
        int MultiplicityNProtonEastX=0;
        int MultiplicityProtonWestX=0;
        int MultiplicityNProtonWestX=0;

        //vz in cm
        double vz = upcEvt->getVertex(0)->getPosZ();
        hVz_check->Fill(vz);
        if(fabs(vz)>100.0) continue; //vz cut

        int NumOfGoodPrimaryTracks=0;

        for(int j=0; j<upcEvt->getNumberOfTracks(); j++) {
            if (upcEvt->getTrack(j)->getFlag(StUPCTrack::kPrimary)
                && upcEvt->getTrack(j)->getFlag(StUPCTrack::kTof)
                && upcEvt->getTrack(j)->getNhits()>20
                && upcEvt->getTrack(j)->getPt()>0.2
                && fabs(upcEvt->getTrack(j)->getEta()) < 0.9
                && upcEvt->getTrack(j)->getEta() < -(vz/250.0) + 0.9
                && upcEvt->getTrack(j)->getEta() > -(vz/250.0) - 0.9) {
                    NumOfGoodPrimaryTracks++;
                }
        }

        if (NumOfGoodPrimaryTracks<2) continue;
        int proton_high_p_candidate=0;

        for (int j = 0; j < upcEvt->getNumberOfTracks(); j++) {
            if( upcEvt->getTrack(j)->getNhits()>20  
                && upcEvt->getTrack(j)->getFlag(StUPCTrack::kTof) 
                && upcEvt->getTrack(j)->getFlag(StUPCTrack::kPrimary)) {

                if  ( upcEvt->getTrack(j)->getPt() > 0.2 &&
                    ( fabs(upcEvt->getTrack(j)->getEta()) < 0.9 &&
                    upcEvt->getTrack(j)->getEta() < -(vz/250.0) + 0.9 &&
                    upcEvt->getTrack(j)->getEta() > -(vz/250.0) - 0.9 ) ){ // fiducial cuts

                    NumOfPrimaryTracksToF++;
                    TVector3 momentum;
                    upcEvt->getTrack(j)->getMomentum(momentum);
                    double p = momentum.Mag();

                    if (isEast) {
                        HistPtTracksEastCut->Fill(upcEvt->getTrack(j)->getPt());
                        HistEtaTracksEastCut->Fill(upcEvt->getTrack(j)->getEta());
                    }
                    if (isWest) {
                        HistPtTracksWestCut->Fill(upcEvt->getTrack(j)->getPt());
                        HistEtaTracksWestCut->Fill(upcEvt->getTrack(j)->getEta());
                    }
                    if((fabs(upcEvt->getTrack(j)->getNSigmasTPCProton())<3) && (p<0.9)) { //proton further selection based on nsigma and p
                        double pt  = upcEvt->getTrack(j)->getPt();
                        double eta = upcEvt->getTrack(j)->getEta();
                        double dca = upcEvt->getTrack(j)->getDcaXY();
                        double rapidity = getRapidity(pt, eta, upcEvt->getTrack(j)->getPhi(), massProton);
                        double distance_from_beam = rapidity - y_beam;
                        double distance_from_edge = y_edge-rapidity;
                        if(upcEvt->getTrack(j)->getCharge()>0){ //proton
                            hDistanceBeam_proton->Fill(distance_from_beam);
                            hDistanceEdge_proton->Fill(distance_from_edge);
                            h3D_protons->Fill(pt, eta, vz);
                            hPtProtonX->Fill(pt);
                            hEtaProtonX->Fill(eta);
                            hRapidityProtonX->Fill(rapidity);
                            MultiplicityProtonX++;
                            if(isEast) {
                                h3D_protonsEast->Fill(pt, eta, vz);
                                hPtProtonEastX->Fill(pt);
                                hEtaProtonEastX->Fill(eta);
                                hRapidityProtonEastX->Fill(rapidity);
                                hDistanceBeam_proton_E->Fill(distance_from_beam);
                                hDistanceEdge_proton_E->Fill(distance_from_edge);
                                hProton_DCA_E->Fill(dca);
                                MultiplicityProtonEastX++;

                                if (b == 0) { hEta_xi0_proton_E->Fill(eta); hRapidity_xi0_proton_E->Fill(rapidity); }
                                if (b == 1) { hEta_xi1_proton_E->Fill(eta); hRapidity_xi1_proton_E->Fill(rapidity); }
                                if (b == 2) { hEta_xi2_proton_E->Fill(eta); hRapidity_xi2_proton_E->Fill(rapidity); }
                                if (b == 3) { hEta_xi3_proton_E->Fill(eta); hRapidity_xi3_proton_E->Fill(rapidity); }
                                if (b == 4) { hEta_xi4_proton_E->Fill(eta); hRapidity_xi4_proton_E->Fill(rapidity); }

                                if (b == 0) { hEta_xi0_proton_C->Fill(-eta); hRapidity_xi0_proton_C->Fill(-rapidity); }
                                if (b == 1) { hEta_xi1_proton_C->Fill(-eta); hRapidity_xi1_proton_C->Fill(-rapidity); }
                                if (b == 2) { hEta_xi2_proton_C->Fill(-eta); hRapidity_xi2_proton_C->Fill(-rapidity); }
                                if (b == 3) { hEta_xi3_proton_C->Fill(-eta); hRapidity_xi3_proton_C->Fill(-rapidity); }
                                if (b == 4) { hEta_xi4_proton_C->Fill(-eta); hRapidity_xi4_proton_C->Fill(-rapidity); }

                            }
                            if(isWest) {
                                h3D_protonsWest->Fill(pt, eta, vz);
                                hPtProtonWestX->Fill(pt);
                                hEtaProtonWestX->Fill(eta);
                                hRapidityProtonWestX->Fill(rapidity);
                                hDistanceBeam_proton_W->Fill(distance_from_beam);
                                hDistanceEdge_proton_W->Fill(distance_from_edge);
                                hProton_DCA_W->Fill(dca);
                                MultiplicityProtonWestX++;

                                if (b == 0) { hEta_xi0_proton_W->Fill(eta); hRapidity_xi0_proton_W->Fill(rapidity); }
                                if (b == 1) { hEta_xi1_proton_W->Fill(eta); hRapidity_xi1_proton_W->Fill(rapidity); }
                                if (b == 2) { hEta_xi2_proton_W->Fill(eta); hRapidity_xi2_proton_W->Fill(rapidity); }
                                if (b == 3) { hEta_xi3_proton_W->Fill(eta); hRapidity_xi3_proton_W->Fill(rapidity); }
                                if (b == 4) { hEta_xi4_proton_W->Fill(eta); hRapidity_xi4_proton_W->Fill(rapidity); }

                                if (b == 0) { hEta_xi0_proton_C->Fill(eta); hRapidity_xi0_proton_C->Fill(rapidity); }
                                if (b == 1) { hEta_xi1_proton_C->Fill(eta); hRapidity_xi1_proton_C->Fill(rapidity); }
                                if (b == 2) { hEta_xi2_proton_C->Fill(eta); hRapidity_xi2_proton_C->Fill(rapidity); }
                                if (b == 3) { hEta_xi3_proton_C->Fill(eta); hRapidity_xi3_proton_C->Fill(rapidity); }
                                if (b == 4) { hEta_xi4_proton_C->Fill(eta); hRapidity_xi4_proton_C->Fill(rapidity); }

                            }

                        } 
                        if(upcEvt->getTrack(j)->getCharge()<0){
                            h3D_antiprotons->Fill(pt, eta, vz);
                            hPtNProtonX->Fill(pt);
                            hEtaNProtonX->Fill(eta);
                            hRapidityNProtonX->Fill(rapidity);
                            MultiplicityNProtonX++;
                            hDistanceBeam_antiproton->Fill(distance_from_beam);
                            hDistanceEdge_antiproton->Fill(distance_from_edge);
                            if(isEast) {
                                h3D_antiprotonsEast->Fill(pt, eta, vz);
                                hPtNProtonEastX->Fill(pt);
                                hEtaNProtonEastX->Fill(eta);
                                hRapidityNProtonEastX->Fill(rapidity);
                                hDistanceBeam_antiproton_E->Fill(distance_from_beam);
                                hDistanceEdge_antiproton_E->Fill(distance_from_edge);
                                hAntiproton_DCA_E->Fill(dca);
                                MultiplicityNProtonEastX++;

                                if (b == 0) { hEta_xi0_antiproton_E->Fill(eta); hRapidity_xi0_antiproton_E->Fill(rapidity); }
                                if (b == 1) { hEta_xi1_antiproton_E->Fill(eta); hRapidity_xi1_antiproton_E->Fill(rapidity); }
                                if (b == 2) { hEta_xi2_antiproton_E->Fill(eta); hRapidity_xi2_antiproton_E->Fill(rapidity); }
                                if (b == 3) { hEta_xi3_antiproton_E->Fill(eta); hRapidity_xi3_antiproton_E->Fill(rapidity); }
                                if (b == 4) { hEta_xi4_antiproton_E->Fill(eta); hRapidity_xi4_antiproton_E->Fill(rapidity); }

                                if (b == 0) { hEta_xi0_antiproton_C->Fill(-eta); hRapidity_xi0_antiproton_C->Fill(-rapidity); }
                                if (b == 1) { hEta_xi1_antiproton_C->Fill(-eta); hRapidity_xi1_antiproton_C->Fill(-rapidity); }
                                if (b == 2) { hEta_xi2_antiproton_C->Fill(-eta); hRapidity_xi2_antiproton_C->Fill(-rapidity); }
                                if (b == 3) { hEta_xi3_antiproton_C->Fill(-eta); hRapidity_xi3_antiproton_C->Fill(-rapidity); }
                                if (b == 4) { hEta_xi4_antiproton_C->Fill(-eta); hRapidity_xi4_antiproton_C->Fill(-rapidity); }
                            }
                            if(isWest) {
                                h3D_antiprotonsWest->Fill(pt, eta, vz);
                                hPtNProtonWestX->Fill(pt);
                                hEtaNProtonWestX->Fill(eta);
                                hRapidityNProtonWestX->Fill(rapidity);
                                hDistanceBeam_antiproton_W->Fill(distance_from_beam);
                                hDistanceEdge_antiproton_W->Fill(distance_from_edge);
                                hAntiproton_DCA_W->Fill(dca);
                                MultiplicityNProtonWestX++;

                                if (b == 0) { hEta_xi0_antiproton_W->Fill(eta); hRapidity_xi0_antiproton_W->Fill(rapidity); }
                                if (b == 1) { hEta_xi1_antiproton_W->Fill(eta); hRapidity_xi1_antiproton_W->Fill(rapidity); }
                                if (b == 2) { hEta_xi2_antiproton_W->Fill(eta); hRapidity_xi2_antiproton_W->Fill(rapidity); }
                                if (b == 3) { hEta_xi3_antiproton_W->Fill(eta); hRapidity_xi3_antiproton_W->Fill(rapidity); }
                                if (b == 4) { hEta_xi4_antiproton_W->Fill(eta); hRapidity_xi4_antiproton_W->Fill(rapidity); }

                                if (b == 0) { hEta_xi0_antiproton_C->Fill(eta); hRapidity_xi0_antiproton_C->Fill(rapidity); }
                                if (b == 1) { hEta_xi1_antiproton_C->Fill(eta); hRapidity_xi1_antiproton_C->Fill(rapidity); }
                                if (b == 2) { hEta_xi2_antiproton_C->Fill(eta); hRapidity_xi2_antiproton_C->Fill(rapidity); }
                                if (b == 3) { hEta_xi3_antiproton_C->Fill(eta); hRapidity_xi3_antiproton_C->Fill(rapidity); }
                                if (b == 4) { hEta_xi4_antiproton_C->Fill(eta); hRapidity_xi4_antiproton_C->Fill(rapidity); }

                            }
                        }
                    }

                    if((fabs(upcEvt->getTrack(j)->getNSigmasTPCProton())<3) && (p>0.9)) { //proton antiproton candidates for p>0.9

                        double DCA = upcEvt->getTrack(j)->getDcaXY();
                        proton_high_p_candidate++;

                        if (upcEvt->getTrack(j)->getCharge()>0) { // proton
                            if(isEast) {
                                hProton_DCA_E_pm->Fill(DCA);
                            } else if (isWest) {
                                hProton_DCA_W_pm->Fill(DCA);
                            }
                        } else if (upcEvt->getTrack(j)->getCharge()<0) { //antiproton
                            if(isEast) {
                                hAntiproton_DCA_E_pm->Fill(DCA);
                            } else if (isWest) {
                                hAntiproton_DCA_W_pm->Fill(DCA);
                            }
                        }
                    } // end of proton antiproton candidates as p>0.9
                }
                if (fabs(upcEvt->getTrack(j)->getEta()) < 0.9 &&
                    upcEvt->getTrack(j)->getEta() < -(vz/250.0) + 0.9 &&
                    upcEvt->getTrack(j)->getEta() > -(vz/250.0) - 0.9) {
                    if (isEast) HistPtTracksEast->Fill(upcEvt->getTrack(j)->getPt());
                    if (isWest) HistPtTracksWest->Fill(upcEvt->getTrack(j)->getPt());
                    if(upcEvt->getTrack(j)->getCharge()>0) { //positive charge
                        hNSigmaPiPlus->Fill(upcEvt->getTrack(j)->getPt(), upcEvt->getTrack(j)->getNSigmasTPCPion());
                        hNSigmaKPlus->Fill(upcEvt->getTrack(j)->getPt(), upcEvt->getTrack(j)->getNSigmasTPCKaon());
                        hNSigmaPPlus->Fill(upcEvt->getTrack(j)->getPt(), upcEvt->getTrack(j)->getNSigmasTPCProton());
                    } else { //negative charge
                            hNSigmaPiMinus->Fill(upcEvt->getTrack(j)->getPt(), upcEvt->getTrack(j)->getNSigmasTPCPion());
                            hNSigmaKMinus->Fill(upcEvt->getTrack(j)->getPt(), upcEvt->getTrack(j)->getNSigmasTPCKaon());
                            hNSigmaPMinus->Fill(upcEvt->getTrack(j)->getPt(), upcEvt->getTrack(j)->getNSigmasTPCProton());
                    }
                    TVector3 momentum_01;
                    upcEvt->getTrack(j)->getMomentum(momentum_01);
                    double pq = (momentum_01.Mag())*(upcEvt->getTrack(j)->getCharge());
                    hdEdx->Fill(pq, upcEvt->getTrack(j)->getDEdxSignal()*1e6); //GeV originally, keV now
                } 
                if (upcEvt->getTrack(j)->getPt()>0.2) {
                    if (isEast) HistEtaTracksEast->Fill(upcEvt->getTrack(j)->getEta());
                    if (isWest) HistEtaTracksWest->Fill(upcEvt->getTrack(j)->getEta()); 
                }
            }

            if(upcEvt->getTrack(j)->getFlag(StUPCTrack::kTof)
                && upcEvt->getTrack(j)->getFlag(StUPCTrack::kV0)
                && upcEvt->getTrack(j)->getNhits()>20 
                && upcEvt->getTrack(j)->getPt()>0.2
                && fabs(upcEvt->getTrack(j)->getEta())<0.9) {
                        
                StUPCTrack* track1=upcEvt->getTrack(j);

                for(size_t k=j+1; k<upcEvt->getNumberOfTracks(); k++) {

                    if (upcEvt->getTrack(k)->getFlag(StUPCTrack::kTof)
                    && upcEvt->getTrack(k)->getFlag(StUPCTrack::kV0)
                    && upcEvt->getTrack(k)->getNhits()>20 
                    && upcEvt->getTrack(k)->getPt()>0.2
                    && fabs(upcEvt->getTrack(k)->getEta())<0.9) {

                        StUPCTrack* track2=upcEvt->getTrack(k);

                        if (track1->getCharge() == track2->getCharge()) continue;

                        double bField = upcEvt->getMagneticField();
                        double beamline[4]; 
                        beamline[0] = upcEvt->getBeamXPosition();
                        beamline[2] = upcEvt->getBeamXSlope();
                        beamline[1] = upcEvt->getBeamYPosition();
                        beamline[3] = upcEvt->getBeamYSlope();

                        const TVector3 PrimVrtx(upcEvt->getVertex(0)->getPosX(), upcEvt->getVertex(0)->getPosY(), upcEvt->getVertex(0)->getPosZ());

                        StUPCV0 L01(track1,track2, massPion, massProton, j, k, PrimVrtx, beamline, bField, false);
                        if ( abs(L01.m()-1.115) < 0.02 ) { //checking hypothesis of lambda/lambdabar, when track1 is a pion and track2 is a proton
                            double rapidity = getRapidity(L01.pt(), L01.eta(), L01.phi(), massLambda);
                            double distance_from_beam = rapidity - y_beam;
                            double distance_from_edge = y_edge-rapidity;
                            if(track2->getCharge()>0) { //if hypothesis is true, then if proton.charge()>0 then its lambda
                                if (LambdaCut(L01)) {
                                    hDistanceBeam_Lambda->Fill(distance_from_beam);
                                    hDistanceEdge_Lambda->Fill(distance_from_edge);
                                    if(isWest) {hDistanceBeam_Lambda_W->Fill(distance_from_beam);hDistanceEdge_Lambda_W->Fill(distance_from_edge);}
                                    if(isEast) {hDistanceBeam_Lambda_E->Fill(distance_from_beam);hDistanceEdge_Lambda_E->Fill(distance_from_edge);}
                                }
                                if(isWest) {
                                    if (LambdaCut(L01, 'd')) hLambda_DCA_W->Fill(L01.dcaDaughters());
                                    if (LambdaCut(L01, 'b')) hLambda_DCABeamLine_W->Fill(L01.DCABeamLine());
                                    if (LambdaCut(L01, 'a')) hLambda_PointingAngle_W->Fill(std::cos(L01.pointingAngle()));
                                    if (LambdaCut(L01, 'l')) hLambda_DecayLength_W->Fill(L01.decayLength());
                                    if (LambdaCut(L01, 'm')) hLambda_Mass_W->Fill(L01.m());
                                    if (LambdaCut(L01, 'n') && L01.decayLength()<3 && std::cos(L01.pointingAngle())>0.925) hLambda_Mass_Background_Decay0Length_W->Fill(L01.m());
                                    if (LambdaCut(L01, 'n') && L01.decayLength()>3 &&std::cos(L01.pointingAngle())<0.925) hLambda_Mass_Background_Pointing0Angle_W->Fill(L01.m());
                                    if (LambdaCut(L01, 'n') && L01.decayLength()>3 && L01.decayLength()<5 && std::cos(L01.pointingAngle())>0.99) hLambda_Mass_Background_DecayLength_W->Fill(L01.m());
                                    if (LambdaCut(L01, 'n') && std::cos(L01.pointingAngle())>0.925 && std::cos(L01.pointingAngle())<0.99 && std::cos(L01.pointingAngle())>5) hLambda_Mass_Background_PointingAngle_W->Fill(L01.m());
                                    if (LambdaCut(L01, 'n') && std::cos(L01.pointingAngle())<-0.925 && L01.decayLength()>5) hLambda_Mass_Background_Neg_PointingAngle_W->Fill(L01.m());
                                    if (LambdaCut(L01, 'p')) hPT_Lambda_W->Fill(L01.pt());
                                    if (LambdaCut(L01)) {hEta_Lambda_W->Fill(L01.eta()); hRapidity_Lambda_W->Fill(rapidity); }
                                    if (LambdaCut(L01)) {
                                        if (b == 0) { hRapidity_xi0_lambda_W->Fill(rapidity); }
                                        if (b == 1) { hRapidity_xi1_lambda_W->Fill(rapidity); }
                                        if (b == 2) { hRapidity_xi2_lambda_W->Fill(rapidity); }
                                        if (b == 3) { hRapidity_xi3_lambda_W->Fill(rapidity); }
                                        if (b == 4) { hRapidity_xi4_lambda_W->Fill(rapidity); }

                                        if (b == 0) { hRapidity_xi0_lambda_C->Fill(rapidity); }
                                        if (b == 1) { hRapidity_xi1_lambda_C->Fill(rapidity); }
                                        if (b == 2) { hRapidity_xi2_lambda_C->Fill(rapidity); }
                                        if (b == 3) { hRapidity_xi3_lambda_C->Fill(rapidity); }
                                        if (b == 4) { hRapidity_xi4_lambda_C->Fill(rapidity); }
                                    }

                                } else if(isEast) {
                                    if (LambdaCut(L01, 'd')) hLambda_DCA_E->Fill(L01.dcaDaughters());
                                    if (LambdaCut(L01, 'b')) hLambda_DCABeamLine_E->Fill(L01.DCABeamLine());
                                    if (LambdaCut(L01, 'a')) hLambda_PointingAngle_E->Fill(std::cos(L01.pointingAngle()));
                                    if (LambdaCut(L01, 'l')) hLambda_DecayLength_E->Fill(L01.decayLength());
                                    if (LambdaCut(L01, 'm')) hLambda_Mass_E->Fill(L01.m());
                                    if (LambdaCut(L01, 'n') && L01.decayLength()<3 && std::cos(L01.pointingAngle())>0.925) hLambda_Mass_Background_Decay0Length_E->Fill(L01.m());
                                    if (LambdaCut(L01, 'n') && L01.decayLength()>3 && std::cos(L01.pointingAngle())<0.925) hLambda_Mass_Background_Pointing0Angle_E->Fill(L01.m());
                                    if (LambdaCut(L01, 'n') && L01.decayLength()>3 && L01.decayLength()<5 && std::cos(L01.pointingAngle())>0.99) hLambda_Mass_Background_DecayLength_E->Fill(L01.m());
                                    if (LambdaCut(L01, 'n') && std::cos(L01.pointingAngle())>0.925 && std::cos(L01.pointingAngle())<0.99 && L01.decayLength()>5) hLambda_Mass_Background_PointingAngle_E->Fill(L01.m());
                                    if (LambdaCut(L01, 'n') && std::cos(L01.pointingAngle())<-0.925 && L01.decayLength()>5) hLambda_Mass_Background_Neg_PointingAngle_E->Fill(L01.m());
                                    if (LambdaCut(L01, 'p')) hPT_Lambda_E->Fill(L01.pt());
                                    if (LambdaCut(L01)) {hEta_Lambda_E->Fill(L01.eta()); hRapidity_Lambda_E->Fill(rapidity); }
                                    if (LambdaCut(L01)) {
                                        if (b == 0) { hRapidity_xi0_lambda_E->Fill(rapidity); }
                                        if (b == 1) { hRapidity_xi1_lambda_E->Fill(rapidity); }
                                        if (b == 2) { hRapidity_xi2_lambda_E->Fill(rapidity); }
                                        if (b == 3) { hRapidity_xi3_lambda_E->Fill(rapidity); }
                                        if (b == 4) { hRapidity_xi4_lambda_E->Fill(rapidity); }

                                        if (b == 0) { hRapidity_xi0_lambda_C->Fill(-rapidity); }
                                        if (b == 1) { hRapidity_xi1_lambda_C->Fill(-rapidity); }
                                        if (b == 2) { hRapidity_xi2_lambda_C->Fill(-rapidity); }
                                        if (b == 3) { hRapidity_xi3_lambda_C->Fill(-rapidity); }
                                        if (b == 4) { hRapidity_xi4_lambda_C->Fill(-rapidity); }
                                    }

                                }
                            }
                            else { //lambdabar
                                if (LambdaCut(L01)) {
                                    hDistanceBeam_antiLambda->Fill(distance_from_beam);
                                    hDistanceEdge_antiLambda->Fill(distance_from_edge);
                                    if(isWest) {hDistanceBeam_antiLambda_W->Fill(distance_from_beam);hDistanceEdge_antiLambda_W->Fill(distance_from_edge);}
                                    if(isEast) {hDistanceBeam_antiLambda_E->Fill(distance_from_beam);hDistanceEdge_antiLambda_E->Fill(distance_from_edge);}
                                }
                                if(isWest) {
                                    if (LambdaCut(L01, 'd')) hAntiLambda_DCA_W->Fill(L01.dcaDaughters());
                                    if (LambdaCut(L01, 'b')) hAntiLambda_DCABeamLine_W->Fill(L01.DCABeamLine());
                                    if (LambdaCut(L01, 'a')) hAntiLambda_PointingAngle_W->Fill(std::cos(L01.pointingAngle()));
                                    if (LambdaCut(L01, 'l')) hAntiLambda_DecayLength_W->Fill(L01.decayLength());
                                    if (LambdaCut(L01, 'm')) hAntiLambda_Mass_W->Fill(L01.m());
                                    if (LambdaCut(L01, 'n') && L01.decayLength()<3 && std::cos(L01.pointingAngle())>0.925) hAntiLambda_Mass_Background_Decay0Length_W->Fill(L01.m());
                                    if (LambdaCut(L01, 'n') && L01.decayLength()>3 && std::cos(L01.pointingAngle())<0.925) hAntiLambda_Mass_Background_Pointing0Angle_W->Fill(L01.m());
                                    if (LambdaCut(L01, 'n') && L01.decayLength()>3 && L01.decayLength()<5 && std::cos(L01.pointingAngle())>0.99) hAntiLambda_Mass_Background_DecayLength_W->Fill(L01.m());
                                    if (LambdaCut(L01, 'n') && std::cos(L01.pointingAngle())>0.925 && std::cos(L01.pointingAngle())<0.99 && L01.decayLength()>5) hAntiLambda_Mass_Background_PointingAngle_W->Fill(L01.m());
                                    if (LambdaCut(L01, 'n') && std::cos(L01.pointingAngle())<-0.925 && L01.decayLength()>5) hAntiLambda_Mass_Background_Neg_PointingAngle_W->Fill(L01.m());
                                    if (LambdaCut(L01, 'p')) hPT_antiLambda_W->Fill(L01.pt());
                                    if (LambdaCut(L01)) {hEta_antiLambda_W->Fill(L01.eta()); hRapidity_antiLambda_W->Fill(rapidity); }
                                    if(LambdaCut(L01)) {
                                        if (b == 0) { hRapidity_xi0_antilambda_W->Fill(rapidity); }
                                        if (b == 1) { hRapidity_xi1_antilambda_W->Fill(rapidity); }
                                        if (b == 2) { hRapidity_xi2_antilambda_W->Fill(rapidity); }
                                        if (b == 3) { hRapidity_xi3_antilambda_W->Fill(rapidity); }
                                        if (b == 4) { hRapidity_xi4_antilambda_W->Fill(rapidity); }

                                        if (b == 0) { hRapidity_xi0_antilambda_C->Fill(rapidity); }
                                        if (b == 1) { hRapidity_xi1_antilambda_C->Fill(rapidity); }
                                        if (b == 2) { hRapidity_xi2_antilambda_C->Fill(rapidity); }
                                        if (b == 3) { hRapidity_xi3_antilambda_C->Fill(rapidity); }
                                        if (b == 4) { hRapidity_xi4_antilambda_C->Fill(rapidity); }
                                    }
                                } else if (isEast) {
                                    if (LambdaCut(L01, 'd')) hAntiLambda_DCA_E->Fill(L01.dcaDaughters());
                                    if (LambdaCut(L01, 'b')) hAntiLambda_DCABeamLine_E->Fill(L01.DCABeamLine());
                                    if (LambdaCut(L01, 'a')) hAntiLambda_PointingAngle_E->Fill(std::cos(L01.pointingAngle()));
                                    if (LambdaCut(L01, 'l')) hAntiLambda_DecayLength_E->Fill(L01.decayLength());
                                    if (LambdaCut(L01, 'm')) hAntiLambda_Mass_E->Fill(L01.m());
                                    if (LambdaCut(L01, 'n') && L01.decayLength()<3 && std::cos(L01.pointingAngle())>0.925) hAntiLambda_Mass_Background_Decay0Length_E->Fill(L01.m());
                                    if (LambdaCut(L01, 'n') && L01.decayLength()>3 && std::cos(L01.pointingAngle())<0.925) hAntiLambda_Mass_Background_Pointing0Angle_E->Fill(L01.m());
                                    if (LambdaCut(L01, 'n') && L01.decayLength()>3 && L01.decayLength()<5 && std::cos(L01.pointingAngle())>0.99) hAntiLambda_Mass_Background_DecayLength_E->Fill(L01.m());
                                    if (LambdaCut(L01, 'n') && std::cos(L01.pointingAngle())>0.925 && std::cos(L01.pointingAngle())<0.99 && L01.decayLength()>5) hAntiLambda_Mass_Background_PointingAngle_E->Fill(L01.m());
                                    if (LambdaCut(L01, 'n') && std::cos(L01.pointingAngle())<-0.925 && L01.decayLength()>5) hAntiLambda_Mass_Background_Neg_PointingAngle_E->Fill(L01.m());
                                    if (LambdaCut(L01, 'p')) hPT_antiLambda_E->Fill(L01.pt());
                                    if (LambdaCut(L01)) {hEta_antiLambda_E->Fill(L01.eta()); hRapidity_antiLambda_E->Fill(rapidity); }
                                    if(LambdaCut(L01)) {
                                        if (b == 0) { hRapidity_xi0_antilambda_E->Fill(rapidity); }
                                        if (b == 1) { hRapidity_xi1_antilambda_E->Fill(rapidity); }
                                        if (b == 2) { hRapidity_xi2_antilambda_E->Fill(rapidity); }
                                        if (b == 3) { hRapidity_xi3_antilambda_E->Fill(rapidity); }
                                        if (b == 4) { hRapidity_xi4_antilambda_E->Fill(rapidity); }

                                        if (b == 0) { hRapidity_xi0_antilambda_C->Fill(-rapidity); }
                                        if (b == 1) { hRapidity_xi1_antilambda_C->Fill(-rapidity); }
                                        if (b == 2) { hRapidity_xi2_antilambda_C->Fill(-rapidity); }
                                        if (b == 3) { hRapidity_xi3_antilambda_C->Fill(-rapidity); }
                                        if (b == 4) { hRapidity_xi4_antilambda_C->Fill(-rapidity); }
                                    }
                                }
                            }
                        }
                        StUPCV0 L02(track1,track2, massProton, massPion, j, k, PrimVrtx, beamline, bField, false);
                        if ( abs(L02.m()-1.115) < 0.02 ) {
                            double rapidity = getRapidity(L02.pt(), L02.eta(), L02.phi(), massLambda);
                            double distance_from_beam = rapidity - y_beam;
                            double distance_from_edge = y_edge-rapidity;
                            if(track1->getCharge()>0) { //if true, then its lambda
                                if (LambdaCut(L02)) {
                                    hDistanceBeam_Lambda->Fill(distance_from_beam);
                                    hDistanceEdge_Lambda->Fill(distance_from_edge);
                                    if(isWest) {hDistanceBeam_Lambda_W->Fill(distance_from_beam);hDistanceEdge_Lambda_W->Fill(distance_from_edge);}
                                    if(isEast) {hDistanceBeam_Lambda_E->Fill(distance_from_beam);hDistanceEdge_Lambda_E->Fill(distance_from_edge);}
                                }
                                if(isWest) {
                                    if (LambdaCut(L02, 'd')) hLambda_DCA_W->Fill(L02.dcaDaughters());
                                    if (LambdaCut(L02, 'b')) hLambda_DCABeamLine_W->Fill(L02.DCABeamLine());
                                    if (LambdaCut(L02, 'a')) hLambda_PointingAngle_W->Fill(std::cos(L02.pointingAngle()));
                                    if (LambdaCut(L02, 'l')) hLambda_DecayLength_W->Fill(L02.decayLength());
                                    if (LambdaCut(L02, 'm')) hLambda_Mass_W->Fill(L02.m());
                                    if (LambdaCut(L02, 'n') && L02.decayLength()<3 && std::cos(L02.pointingAngle())>0.925) hLambda_Mass_Background_Decay0Length_W->Fill(L02.m());
                                    if (LambdaCut(L02, 'n') && L02.decayLength()>3 && std::cos(L02.pointingAngle())<0.925) hLambda_Mass_Background_Pointing0Angle_W->Fill(L02.m());
                                    if (LambdaCut(L02, 'n') && L02.decayLength()>3 && L02.decayLength()<5 && std::cos(L02.pointingAngle())>0.99) hLambda_Mass_Background_DecayLength_W->Fill(L02.m());
                                    if (LambdaCut(L02, 'n') && std::cos(L02.pointingAngle())>0.925 && std::cos(L02.pointingAngle())<0.99 && L02.decayLength()>5) hLambda_Mass_Background_PointingAngle_W->Fill(L02.m());
                                    if (LambdaCut(L02, 'n') && std::cos(L02.pointingAngle())<-0.925 && L02.decayLength()>5) hLambda_Mass_Background_Neg_PointingAngle_W->Fill(L02.m());
                                    if (LambdaCut(L02, 'p')) hPT_Lambda_W->Fill(L02.pt());
                                    if (LambdaCut(L02)) {hEta_Lambda_W->Fill(L02.eta()); hRapidity_Lambda_W->Fill(rapidity); }
                                    if (LambdaCut(L02)) {
                                        if (b == 0) { hRapidity_xi0_lambda_W->Fill(rapidity); }
                                        if (b == 1) { hRapidity_xi1_lambda_W->Fill(rapidity); }
                                        if (b == 2) { hRapidity_xi2_lambda_W->Fill(rapidity); }
                                        if (b == 3) { hRapidity_xi3_lambda_W->Fill(rapidity); }
                                        if (b == 4) { hRapidity_xi4_lambda_W->Fill(rapidity); }

                                        if (b == 0) { hRapidity_xi0_lambda_C->Fill(rapidity); }
                                        if (b == 1) { hRapidity_xi1_lambda_C->Fill(rapidity); }
                                        if (b == 2) { hRapidity_xi2_lambda_C->Fill(rapidity); }
                                        if (b == 3) { hRapidity_xi3_lambda_C->Fill(rapidity); }
                                        if (b == 4) { hRapidity_xi4_lambda_C->Fill(rapidity); }
                                    }
                                } else if(isEast) {
                                    if (LambdaCut(L02, 'd')) hLambda_DCA_E->Fill(L02.dcaDaughters());
                                    if (LambdaCut(L02, 'b')) hLambda_DCABeamLine_E->Fill(L02.DCABeamLine());
                                    if (LambdaCut(L02, 'a')) hLambda_PointingAngle_E->Fill(std::cos(L02.pointingAngle()));
                                    if (LambdaCut(L02, 'l')) hLambda_DecayLength_E->Fill(L02.decayLength());
                                    if (LambdaCut(L02, 'm')) hLambda_Mass_E->Fill(L02.m());
                                    if (LambdaCut(L02, 'n') && L02.decayLength()<3 && std::cos(L02.pointingAngle())>0.925) hLambda_Mass_Background_Decay0Length_E->Fill(L02.m());
                                    if (LambdaCut(L02, 'n') && L02.decayLength()>3 && std::cos(L02.pointingAngle())<0.925) hLambda_Mass_Background_Pointing0Angle_E->Fill(L02.m());
                                    if (LambdaCut(L02, 'n') && L02.decayLength()>3 && L02.decayLength()<5 && std::cos(L02.pointingAngle())>0.99) hLambda_Mass_Background_DecayLength_E->Fill(L02.m());
                                    if (LambdaCut(L02, 'n') && std::cos(L02.pointingAngle())>0.925 && std::cos(L02.pointingAngle())<0.99 && L02.decayLength()>5) hLambda_Mass_Background_PointingAngle_E->Fill(L02.m());
                                    if (LambdaCut(L02, 'n') && std::cos(L02.pointingAngle())<-0.925 && L02.decayLength()>5) hLambda_Mass_Background_Neg_PointingAngle_E->Fill(L02.m());
                                    if (LambdaCut(L02, 'p')) hPT_Lambda_E->Fill(L02.pt());
                                    if (LambdaCut(L02)) {hEta_Lambda_E->Fill(L02.eta()); hRapidity_Lambda_E->Fill(rapidity); }
                                    if (LambdaCut(L02)) {
                                        if (b == 0) { hRapidity_xi0_lambda_E->Fill(rapidity); }
                                        if (b == 1) { hRapidity_xi1_lambda_E->Fill(rapidity); }
                                        if (b == 2) { hRapidity_xi2_lambda_E->Fill(rapidity); }
                                        if (b == 3) { hRapidity_xi3_lambda_E->Fill(rapidity); }
                                        if (b == 4) { hRapidity_xi4_lambda_E->Fill(rapidity); }

                                        if (b == 0) { hRapidity_xi0_lambda_C->Fill(-rapidity); }
                                        if (b == 1) { hRapidity_xi1_lambda_C->Fill(-rapidity); }
                                        if (b == 2) { hRapidity_xi2_lambda_C->Fill(-rapidity); }
                                        if (b == 3) { hRapidity_xi3_lambda_C->Fill(-rapidity); }
                                        if (b == 4) { hRapidity_xi4_lambda_C->Fill(-rapidity); }
                                    }
                                }
                            }
                            else { //lambdabar
                                if (LambdaCut(L02)) {
                                    hDistanceBeam_antiLambda->Fill(distance_from_beam);
                                    hDistanceEdge_antiLambda->Fill(distance_from_edge);
                                    if(isWest) {hDistanceBeam_antiLambda_W->Fill(distance_from_beam);hDistanceEdge_antiLambda_W->Fill(distance_from_edge);}
                                    if(isEast) {hDistanceBeam_antiLambda_E->Fill(distance_from_beam);hDistanceEdge_antiLambda_E->Fill(distance_from_edge);}
                                }
                                if(isWest) {
                                    if (LambdaCut(L02, 'd')) hAntiLambda_DCA_W->Fill(L02.dcaDaughters());
                                    if (LambdaCut(L02, 'b')) hAntiLambda_DCABeamLine_W->Fill(L02.DCABeamLine());
                                    if (LambdaCut(L02, 'a')) hAntiLambda_PointingAngle_W->Fill(std::cos(L02.pointingAngle()));
                                    if (LambdaCut(L02, 'l')) hAntiLambda_DecayLength_W->Fill(L02.decayLength());
                                    if (LambdaCut(L02, 'm')) hAntiLambda_Mass_W->Fill(L02.m());
                                    if (LambdaCut(L02, 'n') && L02.decayLength()<3 && std::cos(L02.pointingAngle())>0.925) hAntiLambda_Mass_Background_Decay0Length_W->Fill(L02.m());
                                    if (LambdaCut(L02, 'n') && L02.decayLength()>3 && std::cos(L02.pointingAngle())<0.925) hAntiLambda_Mass_Background_Pointing0Angle_W->Fill(L02.m());
                                    if (LambdaCut(L02, 'n') && L02.decayLength()>3 && L02.decayLength()<5 && std::cos(L02.pointingAngle())>0.99) hAntiLambda_Mass_Background_DecayLength_W->Fill(L02.m());
                                    if (LambdaCut(L02, 'n') && std::cos(L02.pointingAngle())>0.925 && std::cos(L02.pointingAngle())<0.99 && L02.decayLength()>5) hAntiLambda_Mass_Background_PointingAngle_W->Fill(L02.m());
                                    if (LambdaCut(L02, 'n') && std::cos(L02.pointingAngle())<-0.925 && L02.decayLength()>5) hAntiLambda_Mass_Background_Neg_PointingAngle_W->Fill(L02.m());
                                    if (LambdaCut(L02, 'p')) hPT_antiLambda_W->Fill(L02.pt());
                                    if (LambdaCut(L02)) {hEta_antiLambda_W->Fill(L02.eta()); hRapidity_antiLambda_W->Fill(rapidity); }
                                    if(LambdaCut(L02)) {
                                        if (b == 0) { hRapidity_xi0_antilambda_W->Fill(rapidity); }
                                        if (b == 1) { hRapidity_xi1_antilambda_W->Fill(rapidity); }
                                        if (b == 2) { hRapidity_xi2_antilambda_W->Fill(rapidity); }
                                        if (b == 3) { hRapidity_xi3_antilambda_W->Fill(rapidity); }
                                        if (b == 4) { hRapidity_xi4_antilambda_W->Fill(rapidity); }

                                        if (b == 0) { hRapidity_xi0_antilambda_C->Fill(rapidity); }
                                        if (b == 1) { hRapidity_xi1_antilambda_C->Fill(rapidity); }
                                        if (b == 2) { hRapidity_xi2_antilambda_C->Fill(rapidity); }
                                        if (b == 3) { hRapidity_xi3_antilambda_C->Fill(rapidity); }
                                        if (b == 4) { hRapidity_xi4_antilambda_C->Fill(rapidity); }
                                    }
                                } else if (isEast) {
                                    if (LambdaCut(L02, 'd')) hAntiLambda_DCA_E->Fill(L02.dcaDaughters());
                                    if (LambdaCut(L02, 'b')) hAntiLambda_DCABeamLine_E->Fill(L02.DCABeamLine());
                                    if (LambdaCut(L02, 'a')) hAntiLambda_PointingAngle_E->Fill(std::cos(L02.pointingAngle()));
                                    if (LambdaCut(L02, 'l')) hAntiLambda_DecayLength_E->Fill(L02.decayLength());
                                    if (LambdaCut(L02, 'm')) hAntiLambda_Mass_E->Fill(L02.m());
                                    if (LambdaCut(L02, 'n') && L02.decayLength()<3 && std::cos(L02.pointingAngle())>0.925) hAntiLambda_Mass_Background_Decay0Length_E->Fill(L02.m());
                                    if (LambdaCut(L02, 'n') && L02.decayLength()>3 && std::cos(L02.pointingAngle())<0.925) hAntiLambda_Mass_Background_Pointing0Angle_E->Fill(L02.m());
                                    if (LambdaCut(L02, 'n') && L02.decayLength()>3 && L02.decayLength()<5 && std::cos(L02.pointingAngle())>0.99) hAntiLambda_Mass_Background_DecayLength_E->Fill(L02.m());
                                    if (LambdaCut(L02, 'n') && std::cos(L02.pointingAngle())>0.925 && std::cos(L02.pointingAngle())<0.99 && L02.decayLength()>5) hAntiLambda_Mass_Background_PointingAngle_E->Fill(L02.m());
                                    if (LambdaCut(L02, 'n') && std::cos(L02.pointingAngle())<-0.925 && L02.decayLength()>5) hAntiLambda_Mass_Background_Neg_PointingAngle_E->Fill(L02.m());
                                    if (LambdaCut(L02, 'p')) hPT_antiLambda_E->Fill(L02.pt());
                                    if (LambdaCut(L02)) {hEta_antiLambda_E->Fill(L02.eta()); hRapidity_antiLambda_E->Fill(rapidity); }
                                    if(LambdaCut(L02)) {
                                        if (b == 0) { hRapidity_xi0_antilambda_E->Fill(rapidity); }
                                        if (b == 1) { hRapidity_xi1_antilambda_E->Fill(rapidity); }
                                        if (b == 2) { hRapidity_xi2_antilambda_E->Fill(rapidity); }
                                        if (b == 3) { hRapidity_xi3_antilambda_E->Fill(rapidity); }
                                        if (b == 4) { hRapidity_xi4_antilambda_E->Fill(rapidity); }

                                        if (b == 0) { hRapidity_xi0_antilambda_C->Fill(-rapidity); }
                                        if (b == 1) { hRapidity_xi1_antilambda_C->Fill(-rapidity); }
                                        if (b == 2) { hRapidity_xi2_antilambda_C->Fill(-rapidity); }
                                        if (b == 3) { hRapidity_xi3_antilambda_C->Fill(-rapidity); }
                                        if (b == 4) { hRapidity_xi4_antilambda_C->Fill(-rapidity); }
                                    }
                                }
                            } // charge if-else
                        }  // end of second hypothesis loop
                    } // second loop for lambda selection
                } // second tracks loop
            } // end of loop for lambda selection
        } // end of loop for tracks in an event
       
        HistNumOfPrimaryTracksToF->Fill(NumOfPrimaryTracksToF);
        hMultiplicityProtonX->Fill(MultiplicityProtonX);
        hMultiplicityNProtonX->Fill(MultiplicityNProtonX);

        if (isEast) {
            HistNumOfPToFEast->Fill(NumOfPrimaryTracksToF);
            hMultiplicityProtonEastX->Fill(MultiplicityProtonEastX);
            hMultiplicityNProtonEastX->Fill(MultiplicityNProtonEastX);
            if(proton_high_p_candidate!=0) {
                hPAntiP_Primaries_W->Fill(NumOfGoodPrimaryTracks-proton_high_p_candidate);
            }
        }
        if (isWest) {
            HistNumOfPToFWest->Fill(NumOfPrimaryTracksToF);
            hMultiplicityProtonWestX->Fill(MultiplicityProtonWestX);
            hMultiplicityNProtonWestX->Fill(MultiplicityNProtonWestX);
            if(proton_high_p_candidate!=0) { 
                hPAntiP_Primaries_E->Fill(NumOfGoodPrimaryTracks-proton_high_p_candidate);
            }
        }

       if (proton->branch() < 2)  { //east
       HistXiProtonEast->Fill(proton->xi(254.867));
       HistLogXiProtonEast->Fill(log10(proton->xi(254.867)));
       HistPtProtonEast->Fill(proton->pt());
       HistEtaProtonEast->Fill(proton->eta());
       }
       if (proton->branch() > 1)  { //west
       HistXiProtonWest->Fill(proton->xi(254.867));
       HistLogXiProtonWest->Fill(log10(proton->xi(254.867)));
       HistPtProtonWest->Fill(proton->pt());
       HistEtaProtonWest->Fill(proton->eta());
       }

    } // end of loop for events
    
    HistNumOfPrimaryTracksToF->Write();
    HistNumOfPToFWest->Write();
    HistNumOfPToFEast->Write();

    HistXiProtonWest->Write();
    HistXiProtonEast->Write();

    HistLogXiProtonWest->Write();
    HistLogXiProtonEast->Write();

    HistPtProtonWest->Write();
    HistPtProtonEast->Write();

    HistEtaProtonWest->Write();
    HistEtaProtonEast->Write();

    HistPtTracksWest->Write();
    HistPtTracksEast->Write();

    HistEtaTracksWest->Write();
    HistEtaTracksEast->Write();

    HistPtTracksWestCut->Write();
    HistPtTracksEastCut->Write();

    HistEtaTracksWestCut->Write();
    HistEtaTracksEastCut->Write();

    dirNSigma->cd();

    hNSigmaPiPlus->Write();
    hNSigmaPiMinus->Write();
    hNSigmaKPlus->Write();
    hNSigmaKMinus->Write();
    hNSigmaPPlus->Write();
    hNSigmaPMinus->Write();

    hdEdx->Write();

    dirProtons->cd();

    hPtProtonEastX->Write();
    hPtProtonWestX->Write();
    hPtProtonX->Write();

    hPtNProtonEastX->Write();
    hPtNProtonWestX->Write();
    hPtNProtonX->Write();

    hEtaProtonEastX->Write();
    hEtaProtonWestX->Write();
    hEtaProtonX->Write();
    hEtaNProtonEastX->Write();
    hEtaNProtonWestX->Write();
    hEtaNProtonX->Write();

    hRapidityProtonEastX->Write();
    hRapidityProtonWestX->Write();
    hRapidityProtonX->Write();

    hRapidityNProtonEastX->Write();
    hRapidityNProtonWestX->Write();
    hRapidityNProtonX->Write();


    hMultiplicityProtonEastX->Write();
    hMultiplicityProtonWestX->Write();
    hMultiplicityProtonX->Write();

    hMultiplicityNProtonEastX->Write();
    hMultiplicityNProtonWestX->Write();
    hMultiplicityNProtonX->Write();
    hVz_check->Write();

    h3D_protons->Write();
    h3D_antiprotons->Write();
    h3D_protonsEast->Write();
    h3D_antiprotonsEast->Write();
    h3D_protonsWest->Write();
    h3D_antiprotonsWest->Write();

    hEta_xi0_proton_E->Write();
    hEta_xi1_proton_E->Write();
    hEta_xi2_proton_E->Write();
    hEta_xi3_proton_E->Write();
    hEta_xi4_proton_E->Write();

    hEta_xi0_proton_W->Write();
    hEta_xi1_proton_W->Write();
    hEta_xi2_proton_W->Write();
    hEta_xi3_proton_W->Write();
    hEta_xi4_proton_W->Write();

    hEta_xi0_proton_C->Write();
    hEta_xi1_proton_C->Write();
    hEta_xi2_proton_C->Write();
    hEta_xi3_proton_C->Write();
    hEta_xi4_proton_C->Write();

    hEta_xi0_antiproton_E->Write();
    hEta_xi1_antiproton_E->Write();
    hEta_xi2_antiproton_E->Write();
    hEta_xi3_antiproton_E->Write();
    hEta_xi4_antiproton_E->Write();

    hEta_xi0_antiproton_W->Write();
    hEta_xi1_antiproton_W->Write();
    hEta_xi2_antiproton_W->Write();
    hEta_xi3_antiproton_W->Write();
    hEta_xi4_antiproton_W->Write();

    hEta_xi0_antiproton_C->Write();
    hEta_xi1_antiproton_C->Write();
    hEta_xi2_antiproton_C->Write();
    hEta_xi3_antiproton_C->Write();
    hEta_xi4_antiproton_C->Write();

    hRapidity_xi0_proton_E->Write();
    hRapidity_xi1_proton_E->Write();
    hRapidity_xi2_proton_E->Write();
    hRapidity_xi3_proton_E->Write();
    hRapidity_xi4_proton_E->Write();

    hRapidity_xi0_proton_W->Write();
    hRapidity_xi1_proton_W->Write();
    hRapidity_xi2_proton_W->Write();
    hRapidity_xi3_proton_W->Write();
    hRapidity_xi4_proton_W->Write();

    hRapidity_xi0_proton_C->Write();
    hRapidity_xi1_proton_C->Write();
    hRapidity_xi2_proton_C->Write();
    hRapidity_xi3_proton_C->Write();
    hRapidity_xi4_proton_C->Write();

    hRapidity_xi0_antiproton_E->Write();
    hRapidity_xi1_antiproton_E->Write();
    hRapidity_xi2_antiproton_E->Write();
    hRapidity_xi3_antiproton_E->Write();
    hRapidity_xi4_antiproton_E->Write();

    hRapidity_xi0_antiproton_W->Write();
    hRapidity_xi1_antiproton_W->Write();
    hRapidity_xi2_antiproton_W->Write();
    hRapidity_xi3_antiproton_W->Write();
    hRapidity_xi4_antiproton_W->Write();

    hRapidity_xi0_antiproton_C->Write();
    hRapidity_xi1_antiproton_C->Write();
    hRapidity_xi2_antiproton_C->Write();
    hRapidity_xi3_antiproton_C->Write();
    hRapidity_xi4_antiproton_C->Write();

    hProton_DCA_E->Write();
    hProton_DCA_W->Write();
    hAntiproton_DCA_E->Write();
    hAntiproton_DCA_W->Write();

    hProton_DCA_E_pm->Write();
    hProton_DCA_W_pm->Write();
    hAntiproton_DCA_E_pm->Write();
    hAntiproton_DCA_W_pm->Write();

    // hProton_Primaries_E->Write();
    // hProton_Primaries_W->Write();
    // hAntiproton_Primaries_E->Write();
    // hAntiproton_Primaries_W->Write();

    hPAntiP_Primaries_E->Write();
    hPAntiP_Primaries_W->Write();

    dirLambdas->cd();

    hLambda_DCA_W->Write();
    hLambda_DCA_E->Write();
    hLambda_DCABeamLine_W->Write();
    hLambda_DCABeamLine_E->Write();
    hLambda_PointingAngle_W->Write();
    hLambda_PointingAngle_E->Write();
    hLambda_DecayLength_W->Write();
    hLambda_DecayLength_E->Write();
    hLambda_Mass_W->Write();
    hLambda_Mass_E->Write();
    hLambda_Mass_Background_Decay0Length_W->Write();
    hLambda_Mass_Background_Decay0Length_E->Write();
    hLambda_Mass_Background_DecayLength_W->Write();
    hLambda_Mass_Background_PointingAngle_W->Write();
    hLambda_Mass_Background_DecayLength_E->Write();
    hLambda_Mass_Background_PointingAngle_E->Write();
    hPT_Lambda_E->Write();
    hPT_Lambda_W->Write();
    hEta_Lambda_E->Write();
    hEta_Lambda_W->Write();

    hAntiLambda_DCA_W->Write();
    hAntiLambda_DCA_E->Write();
    hAntiLambda_DCABeamLine_W->Write();
    hAntiLambda_DCABeamLine_E->Write();
    hAntiLambda_PointingAngle_W->Write();
    hAntiLambda_PointingAngle_E->Write();
    hAntiLambda_DecayLength_W->Write();
    hAntiLambda_DecayLength_E->Write();
    hAntiLambda_Mass_W->Write();
    hAntiLambda_Mass_E->Write();
    hAntiLambda_Mass_Background_Decay0Length_W->Write();
    hAntiLambda_Mass_Background_Decay0Length_E->Write();
    hAntiLambda_Mass_Background_DecayLength_W->Write();
    hAntiLambda_Mass_Background_PointingAngle_W->Write();
    hAntiLambda_Mass_Background_DecayLength_E->Write();
    hAntiLambda_Mass_Background_PointingAngle_E->Write();
    hPT_antiLambda_E->Write();
    hPT_antiLambda_W->Write();
    hEta_antiLambda_E->Write();
    hEta_antiLambda_W->Write();

    hLambda_Mass_Background_Neg_PointingAngle_E->Write();
    hLambda_Mass_Background_Neg_PointingAngle_W->Write();
    hAntiLambda_Mass_Background_Neg_PointingAngle_E->Write();
    hAntiLambda_Mass_Background_Neg_PointingAngle_W->Write();
    hLambda_Mass_Background_Pointing0Angle_E->Write();
    hLambda_Mass_Background_Pointing0Angle_W->Write();
    hAntiLambda_Mass_Background_Pointing0Angle_E->Write();
    hAntiLambda_Mass_Background_Pointing0Angle_W->Write();

    hRapidity_Lambda_E->Write();
    hRapidity_Lambda_W->Write();
    hRapidity_antiLambda_E->Write();
    hRapidity_antiLambda_W->Write();

    hRapidity_xi0_lambda_E->Write();
    hRapidity_xi1_lambda_E->Write();
    hRapidity_xi2_lambda_E->Write();
    hRapidity_xi3_lambda_E->Write();
    hRapidity_xi4_lambda_E->Write();

    hRapidity_xi0_lambda_W->Write();
    hRapidity_xi1_lambda_W->Write();
    hRapidity_xi2_lambda_W->Write();
    hRapidity_xi3_lambda_W->Write();
    hRapidity_xi4_lambda_W->Write();

    hRapidity_xi0_lambda_C->Write();
    hRapidity_xi1_lambda_C->Write();
    hRapidity_xi2_lambda_C->Write();
    hRapidity_xi3_lambda_C->Write();
    hRapidity_xi4_lambda_C->Write();

    hRapidity_xi0_antilambda_E->Write();
    hRapidity_xi1_antilambda_E->Write();
    hRapidity_xi2_antilambda_E->Write();
    hRapidity_xi3_antilambda_E->Write();
    hRapidity_xi4_antilambda_E->Write();

    hRapidity_xi0_antilambda_W->Write();
    hRapidity_xi1_antilambda_W->Write();
    hRapidity_xi2_antilambda_W->Write();
    hRapidity_xi3_antilambda_W->Write();
    hRapidity_xi4_antilambda_W->Write();

    hRapidity_xi0_antilambda_C->Write();
    hRapidity_xi1_antilambda_C->Write();
    hRapidity_xi2_antilambda_C->Write();
    hRapidity_xi3_antilambda_C->Write();
    hRapidity_xi4_antilambda_C->Write();

    proton_lambda_distances->cd();

    hDistanceBeam_proton->Write();
    hDistanceBeam_antiproton->Write();
    hDistanceBeam_Lambda->Write();
    hDistanceBeam_antiLambda->Write();

    hDistanceEdge_proton->Write();
    hDistanceEdge_antiproton->Write();
    hDistanceEdge_Lambda->Write();
    hDistanceEdge_antiLambda->Write();

    hDistanceBeam_proton_E->Write();
    hDistanceBeam_proton_W->Write();
    hDistanceBeam_antiproton_E->Write();
    hDistanceBeam_antiproton_W->Write();
    hDistanceBeam_Lambda_E->Write();
    hDistanceBeam_Lambda_W->Write();
    hDistanceBeam_antiLambda_E->Write();
    hDistanceBeam_antiLambda_W->Write();

    hDistanceEdge_proton_E->Write();
    hDistanceEdge_proton_W->Write();
    hDistanceEdge_antiproton_E->Write();
    hDistanceEdge_antiproton_W->Write();
    hDistanceEdge_Lambda_E->Write();
    hDistanceEdge_Lambda_W->Write();
    hDistanceEdge_antiLambda_E->Write();
    hDistanceEdge_antiLambda_W->Write();

    outfile->Close();
    return 0;
}
 
int getXiBin(double xi) {
    const int nXi = 5;
    double xi_low[nXi]  = {0.0,   0.02, 0.05, 0.1, 0.2};
    double xi_high[nXi] = {0.02, 0.05, 0.1, 0.2, 0.4};
    for (int i = 0; i < nXi; i++) {
        if (xi >= xi_low[i] && xi < xi_high[i])
            return i;
    }
    return -1; // poza zakresem
}

bool LambdaCut(const StUPCV0& L, char cut_type) {
    bool mass_cut = (L.m() > (1.1157-0.009675) && L.m() < 1.1157+0.009675); //1.115683
    bool dcaDaughters_cut = L.dcaDaughters()<1.5;
    bool PointingAngle_cut = std::cos(L.pointingAngle())>0.99;
    bool dcaBeamline_cut=L.DCABeamLine()<1.5;
    bool decayLength_cut=L.decayLength()>5;
    bool pt_cut=L.pt()>0.8;
    if(cut_type == 'd') { //without daughtersDCA cut
        return (PointingAngle_cut && dcaBeamline_cut && mass_cut && decayLength_cut && pt_cut);
    } else if(cut_type == 'a') { //without pointing angle cut
        return (dcaDaughters_cut&& dcaBeamline_cut && mass_cut && decayLength_cut && pt_cut);
    } else if(cut_type == 'b') { // without dcaBeamline cut
        return (dcaDaughters_cut && PointingAngle_cut && mass_cut && decayLength_cut && pt_cut);
    } else if(cut_type == 'm') { //without mass cut
        return (dcaDaughters_cut && PointingAngle_cut && dcaBeamline_cut && decayLength_cut && pt_cut);
    } else if(cut_type == 'l') { //without decay length cut
        return (dcaDaughters_cut && PointingAngle_cut && dcaBeamline_cut && mass_cut && pt_cut);
    } else if(cut_type == 'n') { //noise
        return (dcaDaughters_cut && dcaBeamline_cut && pt_cut);
    } else if(cut_type == 'p') { //without pt cut
        return (dcaDaughters_cut && PointingAngle_cut && dcaBeamline_cut && mass_cut && decayLength_cut);
    }
    else { //all cuts
        return (dcaDaughters_cut && PointingAngle_cut && dcaBeamline_cut && mass_cut && decayLength_cut && pt_cut);
    }
}

double getRapidity(const double pt, const double eta, const double phi, const double m) {
    TLorentzVector v;
    v.SetPtEtaPhiM(pt, eta, phi,m);
    return v.Rapidity();
}