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
#include "StPicoHelix.h"
#include "dEdxParameterization.h"

// double TOF_pions=0.2658*3;
// double TOF_kaons=0.2934*3;
// double TOF_protons=0.3869*3;

// double TOFKAONS_pions=0.202*3;
// double TOFKAONS_kaons=0.272*3;
// double TOFKAONS_protons=0.312*3;

// double TOFPIONS_pions=0.160*3;
// double TOFPIONS_kaons=0.262*3;
// double TOFPIONS_protons=0.262*3;

double TOF_pions=0.224*3;
double TOF_kaons=0.250*3;
double TOF_protons=0.300*3;

double TOFKAONS_pions=0.181*3;
double TOFKAONS_kaons=0.220*3;
double TOFKAONS_protons=0.283*3;

double TOFPIONS_pions=0.148*3;
double TOFPIONS_kaons=0.187*3;
double TOFPIONS_protons=0.250*3;

int getXiBin(double xi);
bool LambdaCut(const StUPCV0& L, char cut_type='0');
double getRapidity(const double pt, const double eta, const double phi, const double m);
int getMomentumBin(double p);
double GetTheoreticaldEdx(double p, double massHypothesis, double dx = 4.0);
using namespace std;

int main(int argc, char** argv)  
{
    // dEdxParameterization* dedxParam = new dEdxParameterization("P10", 0, 0, 0, 1, 1);

    TChain *chain = new TChain("mUPCTree"); 

    string inputFileName;
    string outputFileName;

    std::string run = argv[1];
    inputFileName="/data2/sd_star_2017/presel_"+run+".root";      
    outputFileName="Run/"+run+".root";       

    chain->AddFile(inputFileName.c_str());

    static StUPCEvent * upcEvt = 0x0;
    static StRPEvent  * rpEvt = 0x0;

    chain->SetBranchAddress("mUPCEvent", &upcEvt);
    chain->SetBranchAddress("correctedRpEvent",  &rpEvt);


    TFile *outfile = TFile::Open(outputFileName.c_str(), "recreate"); 

    TDirectory* dirNSigma = outfile->mkdir("dEdx_Nsigma");
    TDirectory* dirProtons = outfile->mkdir("Protons");
    TDirectory* dirKaons = outfile->mkdir("Kaons");
    TDirectory* dirPions = outfile->mkdir("Pions");
    TDirectory* protons_lower_momenta = dirProtons->mkdir("protons_lower_momenta");
    TDirectory* protons_mid_momenta = dirProtons->mkdir("protons_mid_momenta");
    TDirectory* protons_higher_momenta = dirProtons->mkdir("protons_higher_momenta");
    TDirectory* kaons_lower_momenta = dirKaons->mkdir("kaons_lower_momenta");
    TDirectory* kaons_mid_momenta = dirKaons->mkdir("kaons_mid_momenta");
    TDirectory* kaons_higher_momenta = dirKaons->mkdir("kaons_higher_momenta");
    TDirectory* pions_lower_momenta = dirPions->mkdir("pions_lower_momenta");
    TDirectory* pions_mid_momenta = dirPions->mkdir("pions_mid_momenta");
    TDirectory* pions_higher_momenta = dirPions->mkdir("pions_higher_momenta");
    TDirectory* dirLambdas = outfile->mkdir("Lambdas");
    TDirectory* proton_lambda_distances = outfile->mkdir("proton_lambda_distances");
    TDirectory* TOF_Analysis = outfile->mkdir("TOF_analysis");
    TDirectory* deltat_protons_lower_momenta = TOF_Analysis->mkdir("deltat_protons_lower_momenta");
    TDirectory* deltat_protons_mid_momenta = TOF_Analysis->mkdir("deltat_protons_mid_momenta");
    TDirectory* deltat_protons_higher_momenta = TOF_Analysis->mkdir("deltat_protons_higher_momenta");
    TDirectory* deltat_kaons_lower_momenta = TOF_Analysis->mkdir("deltat_kaons_lower_momenta");
    TDirectory* deltat_kaons_mid_momenta = TOF_Analysis->mkdir("deltat_kaons_mid_momenta");
    TDirectory* deltat_kaons_higher_momenta = TOF_Analysis->mkdir("deltat_kaons_higher_momenta");
    TDirectory* deltat_pions_lower_momenta = TOF_Analysis->mkdir("deltat_pions_lower_momenta");
    TDirectory* deltat_pions_mid_momenta = TOF_Analysis->mkdir("deltat_pions_mid_momenta");
    TDirectory* deltat_pions_higher_momenta = TOF_Analysis->mkdir("deltat_pions_higher_momenta");
    TDirectory* deltat_anchors = TOF_Analysis->mkdir("deltat_anchors");

    outfile->cd();
    TH1D* HistNumOfPrimaryTracksToF = new TH1D("HistNumOfPrimaryTracksToF", "; Num of tracks; # events", 50 ,0, 50);
    TH1D* HistNumOfPToFWest = new TH1D("HistNumOfPToFWest", "; Num of tracks when proton on west; # events", 50 ,0, 50);
    TH1D* HistNumOfPToFEast = new TH1D("HistNumOfPToFEast", "; Num of tracks when proton on east; # events", 50 ,0, 50);
    TH1D* HistXiProtonWest = new TH1D("HistXiProtonWest", "xi of the proton on west;xi; # events", 100, -0.1, 1);
    TH1D* HistXiProtonEast = new TH1D("HistXiProtonEast", "xi of the proton on east;xi; # events", 100, -0.1, 1);
    TH1D* HistLogXiProtonWest = new TH1D("HistLogXiProtonWest", "logxi of the proton on west; logxi; # events", 100, -6, 0);
    TH1D* HistLogXiProtonEast = new TH1D("HistLogXiProtonEast", "logxi of the proton on east; logxi; # events", 100, -6, 0);
    TH1D* HistPtProtonWest = new TH1D("HistPtProtonWest", "pT of the proton on west; pT [GeV]; # events", 100, 0, 4);
    TH1D* HistPtProtonEast = new TH1D("HistPtProtonEast", "pT of the proton on east; pT [GeV]; # events", 100, 0, 4);
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

    TH1D* hPtProtonEastX_high_pt = new TH1D("hPtProtonEastX_high_pt","pT of protons (East); pT [GeV]; # events",60,0,3);
    TH1D* hPtProtonWestX_high_pt = new TH1D("hPtProtonWestX_high_pt","pT of protons (West); pT [GeV]; # events",60,0,3);
    TH1D* hPtProtonX_high_pt = new TH1D("hPtProtonX_high_pt","pT of protons; pT [GeV]; # events",60,0,3);

    TH1D* hPtNProtonEastX_high_pt = new TH1D("hPtNProtonEastX_high_pt","pT of antiprotons (East); pT [GeV]; # events",60,0,3);
    TH1D* hPtNProtonWestX_high_pt = new TH1D("hPtNProtonWestX_high_pt","pT of antiprotons (West); pT [GeV]; # events",60,0,3);
    TH1D* hPtNProtonX_high_pt = new TH1D("hPtNProtonX_high_pt","pT of antiprotons; pT [GeV]; # events",60,0,3);

    TH1D* hPtProtonEastX_mid_pt = new TH1D("hPtProtonEastX_mid_pt","pT of protons (East); pT [GeV]; # events",60,0,3);
    TH1D* hPtProtonWestX_mid_pt = new TH1D("hPtProtonWestX_mid_pt","pT of protons (West); pT [GeV]; # events",60,0,3);
    TH1D* hPtProtonX_mid_pt = new TH1D("hPtProtonX_mid_pt","pT of protons; pT [GeV]; # events",60,0,3);

    TH1D* hPtNProtonEastX_mid_pt = new TH1D("hPtNProtonEastX_mid_pt","pT of antiprotons (East); pT [GeV]; # events",60,0,3);
    TH1D* hPtNProtonWestX_mid_pt = new TH1D("hPtNProtonWestX_mid_pt","pT of antiprotons (West); pT [GeV]; # events",60,0,3);
    TH1D* hPtNProtonX_mid_pt = new TH1D("hPtNProtonX_mid_pt","pT of antiprotons; pT [GeV]; # events",60,0,3);

    TH1D* hPtProtonEastX = new TH1D("hPtProtonEastX","pT of protons (East); pT [GeV]; # events",60,0,3);
    TH1D* hPtProtonWestX = new TH1D("hPtProtonWestX","pT of protons (West); pT [GeV]; # events",60,0,3);
    TH1D* hPtProtonX = new TH1D("hPtProtonX","pT of protons; pT [GeV]; # events",60,0,3);

    TH1D* hPtNProtonEastX = new TH1D("hPtNProtonEastX","pT of antiprotons (East); pT [GeV]; # events",60,0,3);
    TH1D* hPtNProtonWestX = new TH1D("hPtNProtonWestX","pT of antiprotons (West); pT [GeV]; # events",60,0,3);
    TH1D* hPtNProtonX = new TH1D("hPtNProtonX","pT of antiprotons; pT [GeV]; # events",60,0,3);

    // eta (40 bins, −1–1)

    TH1D* hEtaProtonEastX_high_pt = new TH1D("hEtaProtonEastX_high_pt","#eta of protons (East); #eta; # events",20,-1,1);
    TH1D* hEtaProtonWestX_high_pt = new TH1D("hEtaProtonWestX_high_pt","#eta of protons (West); #eta; # events",20,-1,1);
    TH1D* hEtaProtonX_high_pt = new TH1D("hEtaProtonX_high_pt","#eta of proton; #eta; # events",20,-1,1);

    TH1D* hEtaNProtonEastX_high_pt = new TH1D("hEtaNProtonEastX_high_pt","#eta of antiprotons (East); #eta; # events",20,-1,1);
    TH1D* hEtaNProtonWestX_high_pt = new TH1D("hEtaNProtonWestX_high_pt","#eta of antiprotons (West); #eta; # events",20,-1,1);
    TH1D* hEtaNProtonX_high_pt = new TH1D("hEtaNProtonX_high_pt","#eta of antiprotons; #eta; # events",20,-1,1);

    TH1D* hEtaProtonEastX_mid_pt = new TH1D("hEtaProtonEastX_mid_pt","#eta of protons (East); #eta; # events",20,-1,1);
    TH1D* hEtaProtonWestX_mid_pt = new TH1D("hEtaProtonWestX_mid_pt","#eta of protons (West); #eta; # events",20,-1,1);
    TH1D* hEtaProtonX_mid_pt = new TH1D("hEtaProtonX_mid_pt","#eta of proton; #eta; # events",20,-1,1);

    TH1D* hEtaNProtonEastX_mid_pt = new TH1D("hEtaNProtonEastX_mid_pt","#eta of antiprotons (East); #eta; # events",20,-1,1);
    TH1D* hEtaNProtonWestX_mid_pt = new TH1D("hEtaNProtonWestX_mid_pt","#eta of antiprotons (West); #eta; # events",20,-1,1);
    TH1D* hEtaNProtonX_mid_pt = new TH1D("hEtaNProtonX_mid_pt","#eta of antiprotons; #eta; # events",20,-1,1);

    TH1D* hEtaProtonEastX = new TH1D("hEtaProtonEastX","#eta of protons (East); #eta; # events",20,-1,1);
    TH1D* hEtaProtonWestX = new TH1D("hEtaProtonWestX","#eta of protons (West); #eta; # events",20,-1,1);
    TH1D* hEtaProtonX = new TH1D("hEtaProtonX","#eta of proton; #eta; # events",20,-1,1);

    TH1D* hEtaNProtonEastX = new TH1D("hEtaNProtonEastX","#eta of antiprotons (East); #eta; # events",20,-1,1);
    TH1D* hEtaNProtonWestX = new TH1D("hEtaNProtonWestX","#eta of antiprotons (West); #eta; # events",20,-1,1);
    TH1D* hEtaNProtonX = new TH1D("hEtaNProtonX","#eta of antiprotons; #eta; # events",20,-1,1);

    TH1D* hRapidityProtonEastX_high_pt  = new TH1D("hRapidityProtonEastX_high_pt","#it{y} of protons (East); #it{y}; # events",20,-1,1);
    TH1D* hRapidityProtonWestX_high_pt  = new TH1D("hRapidityProtonWestX_high_pt","#it{y} of protons (West); #it{y}; # events",20,-1,1);
    TH1D* hRapidityProtonX_high_pt      = new TH1D("hRapidityProtonX_high_pt","#it{y} of protons; #it{y}; # events",20,-1,1);

    TH1D* hRapidityNProtonEastX_high_pt = new TH1D("hRapidityNProtonEastX_high_pt","#it{y} of antiprotons (East); #it{y}; # events",20,-1,1);
    TH1D* hRapidityNProtonWestX_high_pt = new TH1D("hRapidityNProtonWestX_high_pt","#it{y} of antiprotons (West); #it{y}; # events",20,-1,1);
    TH1D* hRapidityNProtonX_high_pt     = new TH1D("hRapidityNProtonX_high_pt","#it{y} of antiprotons; #it{y}; # events",20,-1,1);

    TH1D* hRapidityProtonEastX_mid_pt  = new TH1D("hRapidityProtonEastX_mid_pt","#it{y} of protons (East); #it{y}; # events",20,-1,1);
    TH1D* hRapidityProtonWestX_mid_pt  = new TH1D("hRapidityProtonWestX_mid_pt","#it{y} of protons (West); #it{y}; # events",20,-1,1);
    TH1D* hRapidityProtonX_mid_pt      = new TH1D("hRapidityProtonX_mid_pt","#it{y} of protons; #it{y}; # events",20,-1,1);

    TH1D* hRapidityNProtonEastX_mid_pt = new TH1D("hRapidityNProtonEastX_mid_pt","#it{y} of antiprotons (East); #it{y}; # events",20,-1,1);
    TH1D* hRapidityNProtonWestX_mid_pt = new TH1D("hRapidityNProtonWestX_mid_pt","#it{y} of antiprotons (West); #it{y}; # events",20,-1,1);
    TH1D* hRapidityNProtonX_mid_pt     = new TH1D("hRapidityNProtonX_mid_pt","#it{y} of antiprotons; #it{y}; # events",20,-1,1);

    TH1D* hRapidityProtonEastX  = new TH1D("hRapidityProtonEastX","#it{y} of protons (East); #it{y}; # events",20,-1,1);
    TH1D* hRapidityProtonWestX  = new TH1D("hRapidityProtonWestX","#it{y} of protons (West); #it{y}; # events",20,-1,1);
    TH1D* hRapidityProtonX      = new TH1D("hRapidityProtonX","#it{y} of protons; #it{y}; # events",20,-1,1);

    TH1D* hRapidityNProtonEastX = new TH1D("hRapidityNProtonEastX","#it{y} of antiprotons (East); #it{y}; # events",20,-1,1);
    TH1D* hRapidityNProtonWestX = new TH1D("hRapidityNProtonWestX","#it{y} of antiprotons (West); #it{y}; # events",20,-1,1);
    TH1D* hRapidityNProtonX     = new TH1D("hRapidityNProtonX","#it{y} of antiprotons; #it{y}; # events",20,-1,1);

    TH1D* hMultiplicityProtonEastX_high_pt = new TH1D("hMultiplicityProtonEastX_high_pt", "Multiplicity of protons (East); N_{p}; # events", 30, 0, 30);
    TH1D* hMultiplicityProtonWestX_high_pt = new TH1D("hMultiplicityProtonWestX_high_pt", "Multiplicity of protons (West); N_{p}; # events", 30, 0, 30);
    TH1D* hMultiplicityProtonX_high_pt = new TH1D("hMultiplicityProtonX_high_pt", "Multiplicity of protons; N_{p}; # events", 30, 0, 30);

    TH1D* hMultiplicityNProtonEastX_high_pt = new TH1D("hMultiplicityNProtonEastX_high_pt", "Multiplicity of antiprotons (East); N_{#bar{p}}; # events", 30, 0, 30);
    TH1D* hMultiplicityNProtonWestX_high_pt = new TH1D("hMultiplicityNProtonWestX_high_pt", "Multiplicity of antiprotons (West); N_{#bar{p}}; # events", 30, 0, 30);
    TH1D* hMultiplicityNProtonX_high_pt = new TH1D("hMultiplicityNProtonX_high_pt", "Multiplicity of antiprotons; N_{#bar{p}}; # events", 30, 0, 30);

    TH1D* hMultiplicityProtonEastX_mid_pt = new TH1D("hMultiplicityProtonEastX_mid_pt", "Multiplicity of protons (East); N_{p}; # events", 30, 0, 30);
    TH1D* hMultiplicityProtonWestX_mid_pt = new TH1D("hMultiplicityProtonWestX_mid_pt", "Multiplicity of protons (West); N_{p}; # events", 30, 0, 30);
    TH1D* hMultiplicityProtonX_mid_pt = new TH1D("hMultiplicityProtonX_mid_pt", "Multiplicity of protons; N_{p}; # events", 30, 0, 30);

    TH1D* hMultiplicityNProtonEastX_mid_pt = new TH1D("hMultiplicityNProtonEastX_mid_pt", "Multiplicity of antiprotons (East); N_{#bar{p}}; # events", 30, 0, 30);
    TH1D* hMultiplicityNProtonWestX_mid_pt = new TH1D("hMultiplicityNProtonWestX_mid_pt", "Multiplicity of antiprotons (West); N_{#bar{p}}; # events", 30, 0, 30);
    TH1D* hMultiplicityNProtonX_mid_pt = new TH1D("hMultiplicityNProtonX_mid_pt", "Multiplicity of antiprotons; N_{#bar{p}}; # events", 30, 0, 30);

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

    TH1D* hProton_DCA_E_high_pt = new TH1D("hProton_DCA_E_high_pt", "DCA of protons (East); DCA [cm]; #events", 100, 0, 40);
    TH1D* hProton_DCA_W_high_pt = new TH1D("hProton_DCA_W_high_pt", "DCA of protons (West); DCA [cm]; #events", 100, 0, 40);
    TH1D* hAntiproton_DCA_E_high_pt = new TH1D("hAntiProton_DCA_E_high_pt", "DCA of antiprotons (East); DCA [cm]; #events", 100, 0, 40);
    TH1D* hAntiproton_DCA_W_high_pt = new TH1D("hAntiProton_DCA_W_high_pt", "DCA of antiprotons (West); DCA [cm]; #events", 100, 0, 40);

    TH1D* hProton_DCA_E_mid_pt = new TH1D("hProton_DCA_E_mid_pt", "DCA of protons (East); DCA [cm]; #events", 100, 0, 40);
    TH1D* hProton_DCA_W_mid_pt = new TH1D("hProton_DCA_W_mid_pt", "DCA of protons (West); DCA [cm]; #events", 100, 0, 40);
    TH1D* hAntiproton_DCA_E_mid_pt = new TH1D("hAntiProton_DCA_E_mid_pt", "DCA of antiprotons (East); DCA [cm]; #events", 100, 0, 40);
    TH1D* hAntiproton_DCA_W_mid_pt = new TH1D("hAntiProton_DCA_W_mid_pt", "DCA of antiprotons (West); DCA [cm]; #events", 100, 0, 40);

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

    TH1D* hRapidity_xi0_proton_E_mid_pt = new TH1D("hRapidity_xi0_proton_E_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_proton_E_mid_pt = new TH1D("hRapidity_xi1_proton_E_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_proton_E_mid_pt = new TH1D("hRapidity_xi2_proton_E_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_proton_E_mid_pt = new TH1D("hRapidity_xi3_proton_E_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_proton_E_mid_pt = new TH1D("hRapidity_xi4_proton_E_mid_pt", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_proton_W_mid_pt = new TH1D("hRapidity_xi0_proton_W_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_proton_W_mid_pt = new TH1D("hRapidity_xi1_proton_W_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_proton_W_mid_pt = new TH1D("hRapidity_xi2_proton_W_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_proton_W_mid_pt = new TH1D("hRapidity_xi3_proton_W_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_proton_W_mid_pt = new TH1D("hRapidity_xi4_proton_W_mid_pt", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_proton_C_mid_pt = new TH1D("hRapidity_xi0_proton_C_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_proton_C_mid_pt = new TH1D("hRapidity_xi1_proton_C_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_proton_C_mid_pt = new TH1D("hRapidity_xi2_proton_C_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_proton_C_mid_pt = new TH1D("hRapidity_xi3_proton_C_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_proton_C_mid_pt = new TH1D("hRapidity_xi4_proton_C_mid_pt", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_antiproton_E_mid_pt = new TH1D("hRapidity_xi0_antiproton_E_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_antiproton_E_mid_pt = new TH1D("hRapidity_xi1_antiproton_E_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_antiproton_E_mid_pt = new TH1D("hRapidity_xi2_antiproton_E_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_antiproton_E_mid_pt = new TH1D("hRapidity_xi3_antiproton_E_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_antiproton_E_mid_pt = new TH1D("hRapidity_xi4_antiproton_E_mid_pt", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_antiproton_W_mid_pt = new TH1D("hRapidity_xi0_antiproton_W_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_antiproton_W_mid_pt = new TH1D("hRapidity_xi1_antiproton_W_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_antiproton_W_mid_pt = new TH1D("hRapidity_xi2_antiproton_W_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_antiproton_W_mid_pt = new TH1D("hRapidity_xi3_antiproton_W_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_antiproton_W_mid_pt = new TH1D("hRapidity_xi4_antiproton_W_mid_pt", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_antiproton_C_mid_pt = new TH1D("hRapidity_xi0_antiproton_C_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_antiproton_C_mid_pt = new TH1D("hRapidity_xi1_antiproton_C_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_antiproton_C_mid_pt = new TH1D("hRapidity_xi2_antiproton_C_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_antiproton_C_mid_pt = new TH1D("hRapidity_xi3_antiproton_C_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_antiproton_C_mid_pt = new TH1D("hRapidity_xi4_antiproton_C_mid_pt", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_proton_E_high_pt = new TH1D("hRapidity_xi0_proton_E_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_proton_E_high_pt = new TH1D("hRapidity_xi1_proton_E_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_proton_E_high_pt = new TH1D("hRapidity_xi2_proton_E_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_proton_E_high_pt = new TH1D("hRapidity_xi3_proton_E_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_proton_E_high_pt = new TH1D("hRapidity_xi4_proton_E_high_pt", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_proton_W_high_pt = new TH1D("hRapidity_xi0_proton_W_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_proton_W_high_pt = new TH1D("hRapidity_xi1_proton_W_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_proton_W_high_pt = new TH1D("hRapidity_xi2_proton_W_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_proton_W_high_pt = new TH1D("hRapidity_xi3_proton_W_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_proton_W_high_pt = new TH1D("hRapidity_xi4_proton_W_high_pt", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_proton_C_high_pt = new TH1D("hRapidity_xi0_proton_C_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_proton_C_high_pt = new TH1D("hRapidity_xi1_proton_C_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_proton_C_high_pt = new TH1D("hRapidity_xi2_proton_C_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_proton_C_high_pt = new TH1D("hRapidity_xi3_proton_C_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_proton_C_high_pt = new TH1D("hRapidity_xi4_proton_C_high_pt", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_antiproton_E_high_pt = new TH1D("hRapidity_xi0_antiproton_E_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_antiproton_E_high_pt = new TH1D("hRapidity_xi1_antiproton_E_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_antiproton_E_high_pt = new TH1D("hRapidity_xi2_antiproton_E_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_antiproton_E_high_pt = new TH1D("hRapidity_xi3_antiproton_E_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_antiproton_E_high_pt = new TH1D("hRapidity_xi4_antiproton_E_high_pt", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_antiproton_W_high_pt = new TH1D("hRapidity_xi0_antiproton_W_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_antiproton_W_high_pt = new TH1D("hRapidity_xi1_antiproton_W_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_antiproton_W_high_pt = new TH1D("hRapidity_xi2_antiproton_W_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_antiproton_W_high_pt = new TH1D("hRapidity_xi3_antiproton_W_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_antiproton_W_high_pt = new TH1D("hRapidity_xi4_antiproton_W_high_pt", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_antiproton_C_high_pt = new TH1D("hRapidity_xi0_antiproton_C_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_antiproton_C_high_pt = new TH1D("hRapidity_xi1_antiproton_C_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_antiproton_C_high_pt = new TH1D("hRapidity_xi2_antiproton_C_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_antiproton_C_high_pt = new TH1D("hRapidity_xi3_antiproton_C_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_antiproton_C_high_pt = new TH1D("hRapidity_xi4_antiproton_C_high_pt", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_kaon_plus_E = new TH1D("hRapidity_xi0_kaon_plus_E", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_kaon_plus_E = new TH1D("hRapidity_xi1_kaon_plus_E", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_kaon_plus_E = new TH1D("hRapidity_xi2_kaon_plus_E", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_kaon_plus_E = new TH1D("hRapidity_xi3_kaon_plus_E", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_kaon_plus_E = new TH1D("hRapidity_xi4_kaon_plus_E", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_kaon_plus_W = new TH1D("hRapidity_xi0_kaon_plus_W", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_kaon_plus_W = new TH1D("hRapidity_xi1_kaon_plus_W", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_kaon_plus_W = new TH1D("hRapidity_xi2_kaon_plus_W", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_kaon_plus_W = new TH1D("hRapidity_xi3_kaon_plus_W", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_kaon_plus_W = new TH1D("hRapidity_xi4_kaon_plus_W", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_kaon_plus_C = new TH1D("hRapidity_xi0_kaon_plus_C", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_kaon_plus_C = new TH1D("hRapidity_xi1_kaon_plus_C", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_kaon_plus_C = new TH1D("hRapidity_xi2_kaon_plus_C", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_kaon_plus_C = new TH1D("hRapidity_xi3_kaon_plus_C", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_kaon_plus_C = new TH1D("hRapidity_xi4_kaon_plus_C", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_kaon_minus_E = new TH1D("hRapidity_xi0_kaon_minus_E", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_kaon_minus_E = new TH1D("hRapidity_xi1_kaon_minus_E", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_kaon_minus_E = new TH1D("hRapidity_xi2_kaon_minus_E", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_kaon_minus_E = new TH1D("hRapidity_xi3_kaon_minus_E", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_kaon_minus_E = new TH1D("hRapidity_xi4_kaon_minus_E", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_kaon_minus_W = new TH1D("hRapidity_xi0_kaon_minus_W", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_kaon_minus_W = new TH1D("hRapidity_xi1_kaon_minus_W", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_kaon_minus_W = new TH1D("hRapidity_xi2_kaon_minus_W", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_kaon_minus_W = new TH1D("hRapidity_xi3_kaon_minus_W", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_kaon_minus_W = new TH1D("hRapidity_xi4_kaon_minus_W", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_kaon_minus_C = new TH1D("hRapidity_xi0_kaon_minus_C", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_kaon_minus_C = new TH1D("hRapidity_xi1_kaon_minus_C", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_kaon_minus_C = new TH1D("hRapidity_xi2_kaon_minus_C", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_kaon_minus_C = new TH1D("hRapidity_xi3_kaon_minus_C", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_kaon_minus_C = new TH1D("hRapidity_xi4_kaon_minus_C", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_kaon_plus_E_mid_pt = new TH1D("hRapidity_xi0_kaon_plus_E_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_kaon_plus_E_mid_pt = new TH1D("hRapidity_xi1_kaon_plus_E_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_kaon_plus_E_mid_pt = new TH1D("hRapidity_xi2_kaon_plus_E_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_kaon_plus_E_mid_pt = new TH1D("hRapidity_xi3_kaon_plus_E_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_kaon_plus_E_mid_pt = new TH1D("hRapidity_xi4_kaon_plus_E_mid_pt", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_kaon_plus_W_mid_pt = new TH1D("hRapidity_xi0_kaon_plus_W_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_kaon_plus_W_mid_pt = new TH1D("hRapidity_xi1_kaon_plus_W_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_kaon_plus_W_mid_pt = new TH1D("hRapidity_xi2_kaon_plus_W_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_kaon_plus_W_mid_pt = new TH1D("hRapidity_xi3_kaon_plus_W_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_kaon_plus_W_mid_pt = new TH1D("hRapidity_xi4_kaon_plus_W_mid_pt", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_kaon_plus_C_mid_pt = new TH1D("hRapidity_xi0_kaon_plus_C_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_kaon_plus_C_mid_pt = new TH1D("hRapidity_xi1_kaon_plus_C_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_kaon_plus_C_mid_pt = new TH1D("hRapidity_xi2_kaon_plus_C_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_kaon_plus_C_mid_pt = new TH1D("hRapidity_xi3_kaon_plus_C_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_kaon_plus_C_mid_pt = new TH1D("hRapidity_xi4_kaon_plus_C_mid_pt", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_kaon_minus_E_mid_pt = new TH1D("hRapidity_xi0_kaon_minus_E_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_kaon_minus_E_mid_pt = new TH1D("hRapidity_xi1_kaon_minus_E_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_kaon_minus_E_mid_pt = new TH1D("hRapidity_xi2_kaon_minus_E_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_kaon_minus_E_mid_pt = new TH1D("hRapidity_xi3_kaon_minus_E_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_kaon_minus_E_mid_pt = new TH1D("hRapidity_xi4_kaon_minus_E_mid_pt", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_kaon_minus_W_mid_pt = new TH1D("hRapidity_xi0_kaon_minus_W_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_kaon_minus_W_mid_pt = new TH1D("hRapidity_xi1_kaon_minus_W_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_kaon_minus_W_mid_pt = new TH1D("hRapidity_xi2_kaon_minus_W_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_kaon_minus_W_mid_pt = new TH1D("hRapidity_xi3_kaon_minus_W_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_kaon_minus_W_mid_pt = new TH1D("hRapidity_xi4_kaon_minus_W_mid_pt", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_kaon_minus_C_mid_pt = new TH1D("hRapidity_xi0_kaon_minus_C_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_kaon_minus_C_mid_pt = new TH1D("hRapidity_xi1_kaon_minus_C_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_kaon_minus_C_mid_pt = new TH1D("hRapidity_xi2_kaon_minus_C_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_kaon_minus_C_mid_pt = new TH1D("hRapidity_xi3_kaon_minus_C_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_kaon_minus_C_mid_pt = new TH1D("hRapidity_xi4_kaon_minus_C_mid_pt", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_kaon_plus_E_high_pt = new TH1D("hRapidity_xi0_kaon_plus_E_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_kaon_plus_E_high_pt = new TH1D("hRapidity_xi1_kaon_plus_E_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_kaon_plus_E_high_pt = new TH1D("hRapidity_xi2_kaon_plus_E_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_kaon_plus_E_high_pt = new TH1D("hRapidity_xi3_kaon_plus_E_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_kaon_plus_E_high_pt = new TH1D("hRapidity_xi4_kaon_plus_E_high_pt", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_kaon_plus_W_high_pt = new TH1D("hRapidity_xi0_kaon_plus_W_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_kaon_plus_W_high_pt = new TH1D("hRapidity_xi1_kaon_plus_W_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_kaon_plus_W_high_pt = new TH1D("hRapidity_xi2_kaon_plus_W_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_kaon_plus_W_high_pt = new TH1D("hRapidity_xi3_kaon_plus_W_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_kaon_plus_W_high_pt = new TH1D("hRapidity_xi4_kaon_plus_W_high_pt", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_kaon_plus_C_high_pt = new TH1D("hRapidity_xi0_kaon_plus_C_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_kaon_plus_C_high_pt = new TH1D("hRapidity_xi1_kaon_plus_C_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_kaon_plus_C_high_pt = new TH1D("hRapidity_xi2_kaon_plus_C_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_kaon_plus_C_high_pt = new TH1D("hRapidity_xi3_kaon_plus_C_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_kaon_plus_C_high_pt = new TH1D("hRapidity_xi4_kaon_plus_C_high_pt", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_kaon_minus_E_high_pt = new TH1D("hRapidity_xi0_kaon_minus_E_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_kaon_minus_E_high_pt = new TH1D("hRapidity_xi1_kaon_minus_E_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_kaon_minus_E_high_pt = new TH1D("hRapidity_xi2_kaon_minus_E_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_kaon_minus_E_high_pt = new TH1D("hRapidity_xi3_kaon_minus_E_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_kaon_minus_E_high_pt = new TH1D("hRapidity_xi4_kaon_minus_E_high_pt", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_kaon_minus_W_high_pt = new TH1D("hRapidity_xi0_kaon_minus_W_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_kaon_minus_W_high_pt = new TH1D("hRapidity_xi1_kaon_minus_W_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_kaon_minus_W_high_pt = new TH1D("hRapidity_xi2_kaon_minus_W_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_kaon_minus_W_high_pt = new TH1D("hRapidity_xi3_kaon_minus_W_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_kaon_minus_W_high_pt = new TH1D("hRapidity_xi4_kaon_minus_W_high_pt", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_kaon_minus_C_high_pt = new TH1D("hRapidity_xi0_kaon_minus_C_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_kaon_minus_C_high_pt = new TH1D("hRapidity_xi1_kaon_minus_C_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_kaon_minus_C_high_pt = new TH1D("hRapidity_xi2_kaon_minus_C_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_kaon_minus_C_high_pt = new TH1D("hRapidity_xi3_kaon_minus_C_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_kaon_minus_C_high_pt = new TH1D("hRapidity_xi4_kaon_minus_C_high_pt", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_pion_plus_E = new TH1D("hRapidity_xi0_pion_plus_E", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_pion_plus_E = new TH1D("hRapidity_xi1_pion_plus_E", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_pion_plus_E = new TH1D("hRapidity_xi2_pion_plus_E", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_pion_plus_E = new TH1D("hRapidity_xi3_pion_plus_E", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_pion_plus_E = new TH1D("hRapidity_xi4_pion_plus_E", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_pion_plus_W = new TH1D("hRapidity_xi0_pion_plus_W", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_pion_plus_W = new TH1D("hRapidity_xi1_pion_plus_W", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_pion_plus_W = new TH1D("hRapidity_xi2_pion_plus_W", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_pion_plus_W = new TH1D("hRapidity_xi3_pion_plus_W", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_pion_plus_W = new TH1D("hRapidity_xi4_pion_plus_W", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_pion_plus_C = new TH1D("hRapidity_xi0_pion_plus_C", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_pion_plus_C = new TH1D("hRapidity_xi1_pion_plus_C", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_pion_plus_C = new TH1D("hRapidity_xi2_pion_plus_C", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_pion_plus_C = new TH1D("hRapidity_xi3_pion_plus_C", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_pion_plus_C = new TH1D("hRapidity_xi4_pion_plus_C", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_pion_minus_E = new TH1D("hRapidity_xi0_pion_minus_E", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_pion_minus_E = new TH1D("hRapidity_xi1_pion_minus_E", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_pion_minus_E = new TH1D("hRapidity_xi2_pion_minus_E", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_pion_minus_E = new TH1D("hRapidity_xi3_pion_minus_E", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_pion_minus_E = new TH1D("hRapidity_xi4_pion_minus_E", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_pion_minus_W = new TH1D("hRapidity_xi0_pion_minus_W", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_pion_minus_W = new TH1D("hRapidity_xi1_pion_minus_W", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_pion_minus_W = new TH1D("hRapidity_xi2_pion_minus_W", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_pion_minus_W = new TH1D("hRapidity_xi3_pion_minus_W", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_pion_minus_W = new TH1D("hRapidity_xi4_pion_minus_W", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_pion_minus_C = new TH1D("hRapidity_xi0_pion_minus_C", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_pion_minus_C = new TH1D("hRapidity_xi1_pion_minus_C", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_pion_minus_C = new TH1D("hRapidity_xi2_pion_minus_C", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_pion_minus_C = new TH1D("hRapidity_xi3_pion_minus_C", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_pion_minus_C = new TH1D("hRapidity_xi4_pion_minus_C", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_pion_plus_E_mid_pt = new TH1D("hRapidity_xi0_pion_plus_E_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_pion_plus_E_mid_pt = new TH1D("hRapidity_xi1_pion_plus_E_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_pion_plus_E_mid_pt = new TH1D("hRapidity_xi2_pion_plus_E_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_pion_plus_E_mid_pt = new TH1D("hRapidity_xi3_pion_plus_E_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_pion_plus_E_mid_pt = new TH1D("hRapidity_xi4_pion_plus_E_mid_pt", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_pion_plus_W_mid_pt = new TH1D("hRapidity_xi0_pion_plus_W_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_pion_plus_W_mid_pt = new TH1D("hRapidity_xi1_pion_plus_W_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_pion_plus_W_mid_pt = new TH1D("hRapidity_xi2_pion_plus_W_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_pion_plus_W_mid_pt = new TH1D("hRapidity_xi3_pion_plus_W_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_pion_plus_W_mid_pt = new TH1D("hRapidity_xi4_pion_plus_W_mid_pt", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_pion_plus_C_mid_pt = new TH1D("hRapidity_xi0_pion_plus_C_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_pion_plus_C_mid_pt = new TH1D("hRapidity_xi1_pion_plus_C_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_pion_plus_C_mid_pt = new TH1D("hRapidity_xi2_pion_plus_C_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_pion_plus_C_mid_pt = new TH1D("hRapidity_xi3_pion_plus_C_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_pion_plus_C_mid_pt = new TH1D("hRapidity_xi4_pion_plus_C_mid_pt", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_pion_minus_E_mid_pt = new TH1D("hRapidity_xi0_pion_minus_E_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_pion_minus_E_mid_pt = new TH1D("hRapidity_xi1_pion_minus_E_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_pion_minus_E_mid_pt = new TH1D("hRapidity_xi2_pion_minus_E_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_pion_minus_E_mid_pt = new TH1D("hRapidity_xi3_pion_minus_E_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_pion_minus_E_mid_pt = new TH1D("hRapidity_xi4_pion_minus_E_mid_pt", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_pion_minus_W_mid_pt = new TH1D("hRapidity_xi0_pion_minus_W_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_pion_minus_W_mid_pt = new TH1D("hRapidity_xi1_pion_minus_W_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_pion_minus_W_mid_pt = new TH1D("hRapidity_xi2_pion_minus_W_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_pion_minus_W_mid_pt = new TH1D("hRapidity_xi3_pion_minus_W_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_pion_minus_W_mid_pt = new TH1D("hRapidity_xi4_pion_minus_W_mid_pt", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_pion_minus_C_mid_pt = new TH1D("hRapidity_xi0_pion_minus_C_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_pion_minus_C_mid_pt = new TH1D("hRapidity_xi1_pion_minus_C_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_pion_minus_C_mid_pt = new TH1D("hRapidity_xi2_pion_minus_C_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_pion_minus_C_mid_pt = new TH1D("hRapidity_xi3_pion_minus_C_mid_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_pion_minus_C_mid_pt = new TH1D("hRapidity_xi4_pion_minus_C_mid_pt", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_pion_plus_E_high_pt = new TH1D("hRapidity_xi0_pion_plus_E_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_pion_plus_E_high_pt = new TH1D("hRapidity_xi1_pion_plus_E_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_pion_plus_E_high_pt = new TH1D("hRapidity_xi2_pion_plus_E_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_pion_plus_E_high_pt = new TH1D("hRapidity_xi3_pion_plus_E_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_pion_plus_E_high_pt = new TH1D("hRapidity_xi4_pion_plus_E_high_pt", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_pion_plus_W_high_pt = new TH1D("hRapidity_xi0_pion_plus_W_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_pion_plus_W_high_pt = new TH1D("hRapidity_xi1_pion_plus_W_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_pion_plus_W_high_pt = new TH1D("hRapidity_xi2_pion_plus_W_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_pion_plus_W_high_pt = new TH1D("hRapidity_xi3_pion_plus_W_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_pion_plus_W_high_pt = new TH1D("hRapidity_xi4_pion_plus_W_high_pt", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_pion_plus_C_high_pt = new TH1D("hRapidity_xi0_pion_plus_C_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_pion_plus_C_high_pt = new TH1D("hRapidity_xi1_pion_plus_C_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_pion_plus_C_high_pt = new TH1D("hRapidity_xi2_pion_plus_C_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_pion_plus_C_high_pt = new TH1D("hRapidity_xi3_pion_plus_C_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_pion_plus_C_high_pt = new TH1D("hRapidity_xi4_pion_plus_C_high_pt", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_pion_minus_E_high_pt = new TH1D("hRapidity_xi0_pion_minus_E_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_pion_minus_E_high_pt = new TH1D("hRapidity_xi1_pion_minus_E_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_pion_minus_E_high_pt = new TH1D("hRapidity_xi2_pion_minus_E_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_pion_minus_E_high_pt = new TH1D("hRapidity_xi3_pion_minus_E_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_pion_minus_E_high_pt = new TH1D("hRapidity_xi4_pion_minus_E_high_pt", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_pion_minus_W_high_pt = new TH1D("hRapidity_xi0_pion_minus_W_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_pion_minus_W_high_pt = new TH1D("hRapidity_xi1_pion_minus_W_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_pion_minus_W_high_pt = new TH1D("hRapidity_xi2_pion_minus_W_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_pion_minus_W_high_pt = new TH1D("hRapidity_xi3_pion_minus_W_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_pion_minus_W_high_pt = new TH1D("hRapidity_xi4_pion_minus_W_high_pt", ";y;counts", 20, -1, 1);

    TH1D* hRapidity_xi0_pion_minus_C_high_pt = new TH1D("hRapidity_xi0_pion_minus_C_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi1_pion_minus_C_high_pt = new TH1D("hRapidity_xi1_pion_minus_C_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi2_pion_minus_C_high_pt = new TH1D("hRapidity_xi2_pion_minus_C_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi3_pion_minus_C_high_pt = new TH1D("hRapidity_xi3_pion_minus_C_high_pt", ";y;counts", 20, -1, 1);
    TH1D* hRapidity_xi4_pion_minus_C_high_pt = new TH1D("hRapidity_xi4_pion_minus_C_high_pt", ";y;counts", 20, -1, 1);

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
    TH1D* hDistanceBeam_proton_high_pt=new TH1D("hDistanceBeam_proton_high_pt","Distance from beam for protons;#it{y};# events",160,-8,8);
    TH1D* hDistanceBeam_antiproton_high_pt=new TH1D("hDistanceBeam_antiproton_high_pt","Distance from beam for antiprotons;#it{y};# events",160,-8,8);
    TH1D* hDistanceBeam_proton_mid_pt = new TH1D("hDistanceBeam_proton_mid_pt","Distance from beam for protons;#it{y};# events",160,-8,8);
    TH1D* hDistanceBeam_antiproton_mid_pt = new TH1D("hDistanceBeam_antiproton_mid_pt","Distance from beam for antiprotons;#it{y};# events",160,-8,8);
    TH1D* hDistanceBeam_Lambda=new TH1D("hDistanceBeam_Lambda","Distance from beam for #Lambda;#it{y};# events",160,-8,8);
    TH1D* hDistanceBeam_antiLambda=new TH1D("hDistanceBeam_antiLambda","Distance from beam for #bar{#Lambda};#it{y};# events",160,-8,8);

    TH1D* hDistanceEdge_proton=new TH1D("hDistanceEdge_proton","Distance to edge for protons;#it{y};# events",80,-8,8);
    TH1D* hDistanceEdge_antiproton=new TH1D("hDistanceEdge_antiproton","Distance to edge for antiprotons;#it{y};# events",80,-8,8);
    TH1D* hDistanceEdge_proton_high_pt=new TH1D("hDistanceEdge_proton_high_pt","Distance to edge for protons;#it{y};# events",80,-8,8);
    TH1D* hDistanceEdge_antiproton_high_pt=new TH1D("hDistanceEdge_antiproton_high_pt","Distance to edge for antiprotons;#it{y};# events",80,-8,8);
    TH1D* hDistanceEdge_proton_mid_pt = new TH1D("hDistanceEdge_proton_mid_pt","Distance to edge for protons;#it{y};# events",80,-8,8);
    TH1D* hDistanceEdge_antiproton_mid_pt = new TH1D("hDistanceEdge_antiproton_mid_pt","Distance to edge for antiprotons;#it{y};# events",80,-8,8);
    TH1D* hDistanceEdge_Lambda=new TH1D("hDistanceEdge_Lambda","Distance to edge for #Lambda;#it{y};# events",80,-8,8);
    TH1D* hDistanceEdge_antiLambda=new TH1D("hDistanceEdge_antiLambda","Distance to edge for #bar{#Lambda};#it{y};# events",80,-8,8);


    TH1D* hDistanceBeam_proton_E=new TH1D("hDistanceBeam_proton_E","Distance from beam for protons East;#it{y};# events",160,-8,8);
    TH1D* hDistanceBeam_proton_W=new TH1D("hDistanceBeam_proton_W","Distance from beam for protons West;#it{y};# events",160,-8,8);
    TH1D* hDistanceBeam_antiproton_E=new TH1D("hDistanceBeam_antiproton_E","Distance from beam for antiprotons East;#it{y};# events",160,-8,8);
    TH1D* hDistanceBeam_antiproton_W=new TH1D("hDistanceBeam_antiproton_W","Distance from beam for antiprotons West;#it{y};# events",160,-8,8);
    TH1D* hDistanceBeam_proton_E_high_pt=new TH1D("hDistanceBeam_proton_E_high_pt","Distance from beam for protons East;#it{y};# events",160,-8,8);
    TH1D* hDistanceBeam_proton_W_high_pt=new TH1D("hDistanceBeam_proton_W_high_pt","Distance from beam for protons West;#it{y};# events",160,-8,8);
    TH1D* hDistanceBeam_antiproton_E_high_pt=new TH1D("hDistanceBeam_antiproton_E_high_pt","Distance from beam for antiprotons East;#it{y};# events",160,-8,8);
    TH1D* hDistanceBeam_antiproton_W_high_pt=new TH1D("hDistanceBeam_antiproton_W_high_pt","Distance from beam for antiprotons West;#it{y};# events",160,-8,8);
    TH1D* hDistanceBeam_proton_E_mid_pt = new TH1D("hDistanceBeam_proton_E_mid_pt","Distance from beam for protons East;#it{y};# events",160,-8,8);
    TH1D* hDistanceBeam_proton_W_mid_pt = new TH1D("hDistanceBeam_proton_W_mid_pt","Distance from beam for protons West;#it{y};# events",160,-8,8);
    TH1D* hDistanceBeam_antiproton_E_mid_pt = new TH1D("hDistanceBeam_antiproton_E_mid_pt","Distance from beam for antiprotons East;#it{y};# events",160,-8,8);
    TH1D* hDistanceBeam_antiproton_W_mid_pt = new TH1D("hDistanceBeam_antiproton_W_mid_pt","Distance from beam for antiprotons West;#it{y};# events",160,-8,8);
    TH1D* hDistanceBeam_Lambda_E=new TH1D("hDistanceBeam_Lambda_E","Distance from beam for #Lambda East;#it{y};# events",160,-8,8);
    TH1D* hDistanceBeam_Lambda_W=new TH1D("hDistanceBeam_Lambda_W","Distance from beam for #Lambda West;#it{y};# events",160,-8,8);
    TH1D* hDistanceBeam_antiLambda_E=new TH1D("hDistanceBeam_antiLambda_E","Distance from beam for #bar{#Lambda} East;#it{y};# events",160,-8,8);
    TH1D* hDistanceBeam_antiLambda_W=new TH1D("hDistanceBeam_antiLambda_W","Distance from beam for #bar{#Lambda} West;#it{y};# events",160,-8,8);

    TH1D* hDistanceEdge_proton_E=new TH1D("hDistanceEdge_proton_E","Distance to edge for protons East;#it{y};# events",80,-8,8);
    TH1D* hDistanceEdge_proton_W=new TH1D("hDistanceEdge_proton_W","Distance to edge for protons West;#it{y};# events",80,-8,8);
    TH1D* hDistanceEdge_antiproton_E=new TH1D("hDistanceEdge_antiproton_E","Distance to edge for antiprotons East;#it{y};# events",80,-8,8);
    TH1D* hDistanceEdge_antiproton_W=new TH1D("hDistanceEdge_antiproton_W","Distance to edge for antiprotons West;#it{y};# events",80,-8,8);
    TH1D* hDistanceEdge_proton_E_high_pt=new TH1D("hDistanceEdge_proton_E_high_pt","Distance to edge for protons East;#it{y};# events",80,-8,8);
    TH1D* hDistanceEdge_proton_W_high_pt=new TH1D("hDistanceEdge_proton_W_high_pt","Distance to edge for protons West;#it{y};# events",80,-8,8);
    TH1D* hDistanceEdge_antiproton_E_high_pt=new TH1D("hDistanceEdge_antiproton_E_high_pt","Distance to edge for antiprotons East;#it{y};# events",80,-8,8);
    TH1D* hDistanceEdge_antiproton_W_high_pt=new TH1D("hDistanceEdge_antiproton_W_high_pt","Distance to edge for antiprotons West;#it{y};# events",80,-8,8);
    TH1D* hDistanceEdge_proton_E_mid_pt = new TH1D("hDistanceEdge_proton_E_mid_pt","Distance to edge for protons East;#it{y};# events",80,-8,8);
    TH1D* hDistanceEdge_proton_W_mid_pt = new TH1D("hDistanceEdge_proton_W_mid_pt","Distance to edge for protons West;#it{y};# events",80,-8,8);
    TH1D* hDistanceEdge_antiproton_E_mid_pt = new TH1D("hDistanceEdge_antiproton_E_mid_pt","Distance to edge for antiprotons East;#it{y};# events",80,-8,8);
    TH1D* hDistanceEdge_antiproton_W_mid_pt = new TH1D("hDistanceEdge_antiproton_W_mid_pt","Distance to edge for antiprotons West;#it{y};# events",80,-8,8);
    TH1D* hDistanceEdge_Lambda_E=new TH1D("hDistanceEdge_Lambda_E","Distance to edge for #Lambda East;#it{y};# events",80,-8,8);
    TH1D* hDistanceEdge_Lambda_W=new TH1D("hDistanceEdge_Lambda_W","Distance to edge for #Lambda West;#it{y};# events",80,-8,8);
    TH1D* hDistanceEdge_antiLambda_E=new TH1D("hDistanceEdge_antiLambda_E","Distance to edge for #bar{#Lambda} East;#it{y};# events",80,-8,8);
    TH1D* hDistanceEdge_antiLambda_W=new TH1D("hDistanceEdge_antiLambda_W","Distance to edge for #bar{#Lambda} West;#it{y};# events",80,-8,8);

    TH1D* hNumOfParticlesPassSelectionWest=new TH1D("hNumOfParticlesPassSelectionWest","hNumOfParticlesPassSelectionWest;N;# events",10,0,10);
    TH1D* hNumOfParticlesPassSelectionEast=new TH1D("hNumOfParticlesPassSelectionEast","hNumOfParticlesPassSelectionEast;N;# events",10,0,10);

    TH1D* hEventsPassedSelection=new TH1D("hEventsPassedSelection","hEventsPassedSelection;N;# events",5,0,5);

    TH2D* hDeltaT1vsT2_protons = new TH2D("hDeltaT1vsT2_protons", "#Delta t_{1} vs #Delta t_{2} (protons);#Delta t_{1} [ns] ;#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_antiprotons = new TH2D("hDeltaT1vsT2_antiprotons", "#Delta t_{1} vs #Delta t_{2} (antiprotons);#Delta t_{1} [ns];#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_ppbar = new TH2D("hDeltaT1vsT2_ppbar", "#Delta t_{1} vs #Delta t_{2};#Delta t_{1} [ns] ;#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT5_protons = new TH2D("hDeltaT1vsT5_protons", "#Delta t_{1} vs #Delta t_{5} (protons);#Delta t_{1} [ns] ;#Delta t_{5} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT5_antiprotons  = new TH2D("hDeltaT1vsT5_antiprotons", "#Delta t_{1} vs #Delta t_{5} (antiprotons) _mid_p;#Delta t_{1} [ns] ;#Delta t_{5} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_protons_east = new TH2D("hDeltaT1vsT2_protons_east", "#Delta t_{1} vs #Delta t_{2} (protons east);#Delta t_{1} [ns] ;#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_antiprotons_east = new TH2D("hDeltaT1vsT2_antiprotons_east", "#Delta t_{1} vs #Delta t_{2} (antiprotons east);#Delta t_{1} [ns];#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_protons_west = new TH2D("hDeltaT1vsT2_protons_west", "#Delta t_{1} vs #Delta t_{2} (protons west);#Delta t_{1} [ns] ;#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_antiprotons_west = new TH2D("hDeltaT1vsT2_antiprotons_west", "#Delta t_{1} vs #Delta t_{2} (antiprotons west);#Delta t_{1} [ns];#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);

    TH2D* hDeltaT1vsT2_kaons_p = new TH2D("hDeltaT1vsT2_kaons_p", "#Delta t_{1} vs #Delta t_{2} (kaons_p);#Delta t_{1} [ns] ;#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_kaons = new TH2D("hDeltaT1vsT2_kaons", "#Delta t_{1} vs #Delta t_{2};#Delta t_{1} [ns] ;#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT5_kaons_p = new TH2D("hDeltaT1vsT5_kaons_p", "#Delta t_{1} vs #Delta t_{5} (kaons_p);#Delta t_{1} [ns] ;#Delta t_{5} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_kaons_p_east = new TH2D("hDeltaT1vsT2_kaons_p_east", "#Delta t_{1} vs #Delta t_{2} (kaons_p east);#Delta t_{1} [ns] ;#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_kaons_p_west = new TH2D("hDeltaT1vsT2_kaons_p_west", "#Delta t_{1} vs #Delta t_{2} (kaons_p west);#Delta t_{1} [ns] ;#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_kaons_n = new TH2D("hDeltaT1vsT2_kaons_n", "#Delta t_{1} vs #Delta t_{2} (kaons_n);#Delta t_{1} [ns];#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT5_kaons_n  = new TH2D("hDeltaT1vsT5_kaons_n", "#Delta t_{1} vs #Delta t_{5} (kaons_n) _mid_p;#Delta t_{1} [ns] ;#Delta t_{5} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_kaons_n_east = new TH2D("hDeltaT1vsT2_kaons_n_east", "#Delta t_{1} vs #Delta t_{2} (kaons_n east);#Delta t_{1} [ns];#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_kaons_n_west = new TH2D("hDeltaT1vsT2_kaons_n_west", "#Delta t_{1} vs #Delta t_{2} (kaons_n west);#Delta t_{1} [ns];#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);

    TH2D* hDeltaT1vsT2_pions_p = new TH2D("hDeltaT1vsT2_pions_p", "#Delta t_{1} vs #Delta t_{2} (pions_p);#Delta t_{1} [ns] ;#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_pions = new TH2D("hDeltaT1vsT2_pions", "#Delta t_{1} vs #Delta t_{2};#Delta t_{1} [ns] ;#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT5_pions_p = new TH2D("hDeltaT1vsT5_pions_p", "#Delta t_{1} vs #Delta t_{5} (pions_p);#Delta t_{1} [ns] ;#Delta t_{5} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_pions_p_east = new TH2D("hDeltaT1vsT2_pions_p_east", "#Delta t_{1} vs #Delta t_{2} (pions_p east);#Delta t_{1} [ns] ;#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_pions_p_west = new TH2D("hDeltaT1vsT2_pions_p_west", "#Delta t_{1} vs #Delta t_{2} (pions_p west);#Delta t_{1} [ns] ;#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_pions_n = new TH2D("hDeltaT1vsT2_pions_n", "#Delta t_{1} vs #Delta t_{2} (pions_n);#Delta t_{1} [ns];#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT5_pions_n  = new TH2D("hDeltaT1vsT5_pions_n", "#Delta t_{1} vs #Delta t_{5} (pions_n) _mid_p;#Delta t_{1} [ns] ;#Delta t_{5} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_pions_n_east = new TH2D("hDeltaT1vsT2_pions_n_east", "#Delta t_{1} vs #Delta t_{2} (pions_n east);#Delta t_{1} [ns];#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_pions_n_west = new TH2D("hDeltaT1vsT2_pions_n_west", "#Delta t_{1} vs #Delta t_{2} (pions_n west);#Delta t_{1} [ns];#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);

    TH1D* hDeltaT_pionpion=new TH1D("hDeltaT_pion_pion","#Delta t pion-pion hypothesis;#Delta t [ns];# events",400,-20.0,20.0);
    TH1D* hDeltaT_kaonpion=new TH1D("hDeltaT_kaon_pion","#Delta t kaon-pion hypothesis;#Delta t [ns];# events",400,-20.0,20.0);
    TH1D* hDeltaT_protonproton=new TH1D("hDeltaT_proton_proton","#Delta t proton-proton hypothesis;#Delta t [ns];# events",400,-20.0,20.0);

    TH2D* hDeltaT1vsT2_protons_high_p = new TH2D("hDeltaT1vsT2_protons_high_p", "#Delta t_{1} vs #Delta t_{2} (protons) _high_p;#Delta t_{1} [ns] ;#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_ppbar_high_p = new TH2D("hDeltaT1vsT2_ppbar_high_p", "#Delta t_{1} vs #Delta t_{2};#Delta t_{1} [ns] ;#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT5_protons_high_p  = new TH2D("hDeltaT1vsT5_protons_high_p", "#Delta t_{1} vs #Delta t_{5} (protons) _high_p;#Delta t_{1} [ns] ;#Delta t_{5} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT5_antiprotons_high_p  = new TH2D("hDeltaT1vsT5_antiprotons_high_p", "#Delta t_{1} vs #Delta t_{5} (antiprotons) _mid_p;#Delta t_{1} [ns] ;#Delta t_{5} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_antiprotons_high_p = new TH2D("hDeltaT1vsT2_antiprotons_high_p", "#Delta t_{1} vs #Delta t_{2} (antiprotons) _high_p;#Delta t_{1} [ns];#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_protons_east_high_p = new TH2D("hDeltaT1vsT2_protons_east_high_p", "#Delta t_{1} vs #Delta t_{2} (protons east) _high_p;#Delta t_{1} [ns] ;#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_antiprotons_east_high_p = new TH2D("hDeltaT1vsT2_antiprotons_east_high_p", "#Delta t_{1} vs #Delta t_{2} (antiprotons east) _high_p;#Delta t_{1} [ns];#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_protons_west_high_p = new TH2D("hDeltaT1vsT2_protons_west_high_p", "#Delta t_{1} vs #Delta t_{2} (protons west) _high_p;#Delta t_{1} [ns] ;#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_antiprotons_west_high_p = new TH2D("hDeltaT1vsT2_antiprotons_west_high_p", "#Delta t_{1} vs #Delta t_{2} (antiprotons west) _high_p;#Delta t_{1} [ns];#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);

    TH2D* hDeltaT1vsT2_kaons_p_high_p = new TH2D("hDeltaT1vsT2_kaons_p_high_p", "#Delta t_{1} vs #Delta t_{2} (kaons_p) _high_p;#Delta t_{1} [ns] ;#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_kaons_high_p = new TH2D("hDeltaT1vsT2_kaons_high_p", "#Delta t_{1} vs #Delta t_{2} _high_p;#Delta t_{1} [ns] ;#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT5_kaons_p_high_p = new TH2D("hDeltaT1vsT5_kaons_p_high_p", "#Delta t_{1} vs #Delta t_{5} (kaons_p) _high_p;#Delta t_{1} [ns] ;#Delta t_{5} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT5_kaons_n_high_p = new TH2D("hDeltaT1vsT5_kaons_n_high_p", "#Delta t_{1} vs #Delta t_{5} (kaons_n) _high_p;#Delta t_{1} [ns] ;#Delta t_{5} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_kaons_n_high_p = new TH2D("hDeltaT1vsT2_kaons_n_high_p", "#Delta t_{1} vs #Delta t_{2} (kaons_n) _high_p;#Delta t_{1} [ns];#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_kaons_p_east_high_p = new TH2D("hDeltaT1vsT2_kaons_p_east_high_p", "#Delta t_{1} vs #Delta t_{2} (kaons_p east) _high_p;#Delta t_{1} [ns] ;#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_kaons_n_east_high_p = new TH2D("hDeltaT1vsT2_kaons_n_east_high_p", "#Delta t_{1} vs #Delta t_{2} (kaons_n east) _high_p;#Delta t_{1} [ns];#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_kaons_p_west_high_p = new TH2D("hDeltaT1vsT2_kaons_p_west_high_p", "#Delta t_{1} vs #Delta t_{2} (kaons_p west) _high_p;#Delta t_{1} [ns] ;#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_kaons_n_west_high_p = new TH2D("hDeltaT1vsT2_kaons_n_west_high_p", "#Delta t_{1} vs #Delta t_{2} (kaons_n west) _high_p;#Delta t_{1} [ns];#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);

    TH2D* hDeltaT1vsT2_pions_p_high_p = new TH2D("hDeltaT1vsT2_pions_p_high_p", "#Delta t_{1} vs #Delta t_{2} (pions_p) _high_p;#Delta t_{1} [ns] ;#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_pions_high_p = new TH2D("hDeltaT1vsT2_pions_high_p", "#Delta t_{1} vs #Delta t_{2} _high_p;#Delta t_{1} [ns] ;#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT5_pions_p_high_p = new TH2D("hDeltaT1vsT5_pions_p_high_p", "#Delta t_{1} vs #Delta t_{5} (pions_p) _high_p;#Delta t_{1} [ns] ;#Delta t_{5} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT5_pions_n_high_p = new TH2D("hDeltaT1vsT5_pions_n_high_p", "#Delta t_{1} vs #Delta t_{5} (pions_n) _high_p;#Delta t_{1} [ns] ;#Delta t_{5} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_pions_n_high_p = new TH2D("hDeltaT1vsT2_pions_n_high_p", "#Delta t_{1} vs #Delta t_{2} (pions_n) _high_p;#Delta t_{1} [ns];#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_pions_p_east_high_p = new TH2D("hDeltaT1vsT2_pions_p_east_high_p", "#Delta t_{1} vs #Delta t_{2} (pions_p east) _high_p;#Delta t_{1} [ns] ;#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_pions_n_east_high_p = new TH2D("hDeltaT1vsT2_pions_n_east_high_p", "#Delta t_{1} vs #Delta t_{2} (pions_n east) _high_p;#Delta t_{1} [ns];#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_pions_p_west_high_p = new TH2D("hDeltaT1vsT2_pions_p_west_high_p", "#Delta t_{1} vs #Delta t_{2} (pions_p west) _high_p;#Delta t_{1} [ns] ;#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_pions_n_west_high_p = new TH2D("hDeltaT1vsT2_pions_n_west_high_p", "#Delta t_{1} vs #Delta t_{2} (pions_n west) _high_p;#Delta t_{1} [ns];#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    
    TH1D* hDeltaT_pionpion_high_p = new TH1D("hDeltaT_pion_pion_high_p", "#Delta t pion-pion hypothesis _high_p;#Delta t [ns];# events", 400, -20.0, 20.0);
    TH1D* hDeltaT_kaonpion_high_p = new TH1D("hDeltaT_kaon_pion_high_p", "#Delta t kaon-pion hypothesis _high_p;#Delta t [ns];# events", 400, -20.0, 20.0);
    TH1D* hDeltaT_protonproton_high_p = new TH1D("hDeltaT_proton_proton_high_p", "#Delta t proton-proton hypothesis _high_p;#Delta t [ns];# events", 400, -20.0, 20.0);

    TH2D* hDeltaT1vsT2_protons_mid_p = new TH2D("hDeltaT1vsT2_protons_mid_p", "#Delta t_{1} vs #Delta t_{2} (protons) _mid_p;#Delta t_{1} [ns] ;#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_ppbar_mid_p = new TH2D("hDeltaT1vsT2_ppbar_mid_p", "#Delta t_{1} vs #Delta t_{2};#Delta t_{1} [ns] ;#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT5_protons_mid_p  = new TH2D("hDeltaT1vsT5_protons_mid_p", "#Delta t_{1} vs #Delta t_{5} (protons) _mid_p;#Delta t_{1} [ns] ;#Delta t_{5} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT5_antiprotons_mid_p  = new TH2D("hDeltaT1vsT5_antiprotons_mid_p", "#Delta t_{1} vs #Delta t_{5} (antiprotons) _mid_p;#Delta t_{1} [ns] ;#Delta t_{5} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_antiprotons_mid_p = new TH2D("hDeltaT1vsT2_antiprotons_mid_p", "#Delta t_{1} vs #Delta t_{2} (antiprotons) _mid_p;#Delta t_{1} [ns];#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_protons_east_mid_p = new TH2D("hDeltaT1vsT2_protons_east_mid_p", "#Delta t_{1} vs #Delta t_{2} (protons east) _mid_p;#Delta t_{1} [ns] ;#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_antiprotons_east_mid_p = new TH2D("hDeltaT1vsT2_antiprotons_east_mid_p", "#Delta t_{1} vs #Delta t_{2} (antiprotons east) _mid_p;#Delta t_{1} [ns];#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_protons_west_mid_p = new TH2D("hDeltaT1vsT2_protons_west_mid_p", "#Delta t_{1} vs #Delta t_{2} (protons west) _mid_p;#Delta t_{1} [ns] ;#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_antiprotons_west_mid_p = new TH2D("hDeltaT1vsT2_antiprotons_west_mid_p", "#Delta t_{1} vs #Delta t_{2} (antiprotons west) _mid_p;#Delta t_{1} [ns];#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);

    TH2D* hDeltaT1vsT2_kaons_p_mid_p = new TH2D("hDeltaT1vsT2_kaons_p_mid_p", "#Delta t_{1} vs #Delta t_{2} (kaons_p) _mid_p;#Delta t_{1} [ns] ;#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_kaons_mid_p = new TH2D("hDeltaT1vsT2_kaons_mid_p", "#Delta t_{1} vs #Delta t_{2} _mid_p;#Delta t_{1} [ns] ;#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT5_kaons_p_mid_p = new TH2D("hDeltaT1vsT5_kaons_p_mid_p", "#Delta t_{1} vs #Delta t_{5} (kaons_p) _mid_p;#Delta t_{1} [ns] ;#Delta t_{5} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT5_kaons_n_mid_p = new TH2D("hDeltaT1vsT5_kaons_n_mid_p", "#Delta t_{1} vs #Delta t_{5} (kaons_n) _mid_p;#Delta t_{1} [ns] ;#Delta t_{5} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_kaons_n_mid_p = new TH2D("hDeltaT1vsT2_kaons_n_mid_p", "#Delta t_{1} vs #Delta t_{2} (kaons_n) _mid_p;#Delta t_{1} [ns];#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_kaons_p_east_mid_p = new TH2D("hDeltaT1vsT2_kaons_p_east_mid_p", "#Delta t_{1} vs #Delta t_{2} (kaons_p east) _mid_p;#Delta t_{1} [ns] ;#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_kaons_n_east_mid_p = new TH2D("hDeltaT1vsT2_kaons_n_east_mid_p", "#Delta t_{1} vs #Delta t_{2} (kaons_n east) _mid_p;#Delta t_{1} [ns];#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_kaons_p_west_mid_p = new TH2D("hDeltaT1vsT2_kaons_p_west_mid_p", "#Delta t_{1} vs #Delta t_{2} (kaons_p west) _mid_p;#Delta t_{1} [ns] ;#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_kaons_n_west_mid_p = new TH2D("hDeltaT1vsT2_kaons_n_west_mid_p", "#Delta t_{1} vs #Delta t_{2} (kaons_n west) _mid_p;#Delta t_{1} [ns];#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);

    TH2D* hDeltaT1vsT2_pions_p_mid_p = new TH2D("hDeltaT1vsT2_pions_p_mid_p", "#Delta t_{1} vs #Delta t_{2} (pions_p) _mid_p;#Delta t_{1} [ns] ;#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_pions_mid_p = new TH2D("hDeltaT1vsT2_pions_mid_p", "#Delta t_{1} vs #Delta t_{2} _mid_p;#Delta t_{1} [ns] ;#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT5_pions_p_mid_p = new TH2D("hDeltaT1vsT5_pions_p_mid_p", "#Delta t_{1} vs #Delta t_{5} (pions_p) _mid_p;#Delta t_{1} [ns] ;#Delta t_{5} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT5_pions_n_mid_p = new TH2D("hDeltaT1vsT5_pions_n_mid_p", "#Delta t_{1} vs #Delta t_{5} (pions_n) _mid_p;#Delta t_{1} [ns] ;#Delta t_{5} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_pions_n_mid_p = new TH2D("hDeltaT1vsT2_pions_n_mid_p", "#Delta t_{1} vs #Delta t_{2} (pions_n) _mid_p;#Delta t_{1} [ns];#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_pions_p_east_mid_p = new TH2D("hDeltaT1vsT2_pions_p_east_mid_p", "#Delta t_{1} vs #Delta t_{2} (pions_p east) _mid_p;#Delta t_{1} [ns] ;#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_pions_n_east_mid_p = new TH2D("hDeltaT1vsT2_pions_n_east_mid_p", "#Delta t_{1} vs #Delta t_{2} (pions_n east) _mid_p;#Delta t_{1} [ns];#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_pions_p_west_mid_p = new TH2D("hDeltaT1vsT2_pions_p_west_mid_p", "#Delta t_{1} vs #Delta t_{2} (pions_p west) _mid_p;#Delta t_{1} [ns] ;#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);
    TH2D* hDeltaT1vsT2_pions_n_west_mid_p = new TH2D("hDeltaT1vsT2_pions_n_west_mid_p", "#Delta t_{1} vs #Delta t_{2} (pions_n west) _mid_p;#Delta t_{1} [ns];#Delta t_{2} [ns]", 200, -20.0, 20.0, 200, -20.0, 20.0);

    TH1D* hDeltaT_pionpion_mid_p = new TH1D("hDeltaT_pion_pion_mid_p", "#Delta t pion-pion hypothesis _mid_p;#Delta t [ns];# events", 400, -20.0, 20.0);
    TH1D* hDeltaT_kaonpion_mid_p = new TH1D("hDeltaT_kaon_pion_mid_p", "#Delta t kaon-pion hypothesis _mid_p;#Delta t [ns];# events", 400, -20.0, 20.0);
    TH1D* hDeltaT_protonproton_mid_p = new TH1D("hDeltaT_proton_proton_mid_p", "#Delta t proton-proton hypothesis _mid_p;#Delta t [ns];# events", 400, -20.0, 20.0);

    TH1D* hDeltaT_global_anchor_Kppi_high_pt = new TH1D("hDeltaT_global_anchor_Kppi_high_pt", "#Delta t K^{+} - pion (anchor) (high p_{T});#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_KpK_high_pt = new TH1D("hDeltaT_global_anchor_KpK_high_pt", "#Delta t K^{+} - kaon (anchor) (high p_{T});#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_Kpp_high_pt = new TH1D("hDeltaT_global_anchor_Kpp_high_pt", "#Delta t K^{+} - proton (anchor) (high p_{T});#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_Knpi_high_pt = new TH1D("hDeltaT_global_anchor_Knpi_high_pt", "#Delta t K^{-} - pion (anchor) (high p_{T});#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_KnK_high_pt = new TH1D("hDeltaT_global_anchor_KnK_high_pt", "#Delta t K^{-} - kaon (anchor) (high p_{T});#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_Knp_high_pt = new TH1D("hDeltaT_global_anchor_Knp_high_pt", "#Delta t K^{-} - proton (anchor) (high p_{T});#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_Kppi_mid_pt = new TH1D("hDeltaT_global_anchor_Kppi_mid_pt", "#Delta t K^{+} - pion (anchor) (mid p_{T});#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_KpK_mid_pt = new TH1D("hDeltaT_global_anchor_KpK_mid_pt", "#Delta t K^{+} - kaon (anchor) (mid p_{T});#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_Kpp_mid_pt = new TH1D("hDeltaT_global_anchor_Kpp_mid_pt", "#Delta t K^{+} - proton (anchor) (mid p_{T});#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_Knpi_mid_pt = new TH1D("hDeltaT_global_anchor_Knpi_mid_pt", "#Delta t K^{-} - pion (anchor) (mid p_{T});#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_KnK_mid_pt = new TH1D("hDeltaT_global_anchor_KnK_mid_pt", "#Delta t K^{-} - kaon (anchor) (mid p_{T});#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_Knp_mid_pt = new TH1D("hDeltaT_global_anchor_Knp_mid_pt", "#Delta t K^{-} - proton (anchor) (mid p_{T});#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_Kppi = new TH1D("hDeltaT_global_anchor_Kppi", "#Delta t K^{+} - pion (anchor);#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_KpK = new TH1D("hDeltaT_global_anchor_KpK", "#Delta t K^{+} - kaon (anchor);#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_Kpp = new TH1D("hDeltaT_global_anchor_Kpp", "#Delta t K^{+} - proton (anchor);# Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_Knpi = new TH1D("hDeltaT_global_anchor_Knpi", "#Delta t K^{-} - pion (anchor);#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_KnK = new TH1D("hDeltaT_global_anchor_KnK", "#Delta t K^{-} - kaon (anchor);#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_Knp = new TH1D("hDeltaT_global_anchor_Knp", "#Delta t K^{-} - proton (anchor);#Delta t [ns];# events", 80, -4.0, 4.0);

    TH1D* hDeltaT_global_anchor_Pippi_high_pt = new TH1D("hDeltaT_global_anchor_Pippi_high_pt", "#Delta t #pi^{+} - pion (anchor) (high p_{T});#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_PipK_high_pt = new TH1D("hDeltaT_global_anchor_PipK_high_pt", "#Delta t #pi^{+} - kaon (anchor) (high p_{T});#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_Pipp_high_pt = new TH1D("hDeltaT_global_anchor_Pipp_high_pt", "#Delta t #pi^{+} - proton (anchor) (high p_{T});#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_Pinpi_high_pt = new TH1D("hDeltaT_global_anchor_Pinpi_high_pt", "#Delta t #pi^{-} - pion (anchor) (high p_{T});#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_PinK_high_pt = new TH1D("hDeltaT_global_anchor_PinK_high_pt", "#Delta t #pi^{-} - kaon (anchor) (high p_{T});#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_Pinp_high_pt = new TH1D("hDeltaT_global_anchor_Pinp_high_pt", "#Delta t #pi^{-} - proton (anchor) (high p_{T});#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_Pippi_mid_pt = new TH1D("hDeltaT_global_anchor_Pippi_mid_pt", "# Delta t #pi^{+} - pion (anchor) (mid p_{T});# Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_PipK_mid_pt = new TH1D("hDeltaT_global_anchor_PipK_mid_pt", "# Delta t #pi^{+} - kaon (anchor) (mid p_{T});# Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_Pipp_mid_pt = new TH1D("hDeltaT_global_anchor_Pipp_mid_pt", "# Delta t #pi^{+} - proton (anchor) (mid p_{T});# Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_Pinpi_mid_pt = new TH1D("hDeltaT_global_anchor_Pinpi_mid_pt", "#Delta t #pi^{-} - pion (anchor) (mid p_{T});#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_PinK_mid_pt = new TH1D("hDeltaT_global_anchor_PinK_mid_pt", "#Delta t #pi^{-} - kaon (anchor) (mid p_{T});#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_Pinp_mid_pt = new TH1D("hDeltaT_global_anchor_Pinp_mid_pt", "#Delta t #pi^{-} - proton (anchor) (mid p_{T});#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_Pippi = new TH1D("hDeltaT_global_anchor_Pippi", "#Delta t #pi^{+} - pion (anchor);#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_PipK = new TH1D("hDeltaT_global_anchor_PipK", "#Delta t #pi^{+} - kaon (anchor);#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_Pipp = new TH1D("hDeltaT_global_anchor_Pipp", "#Delta t #pi^{+} - proton (anchor);#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_Pinpi = new TH1D("hDeltaT_global_anchor_Pinpi", "#Delta t #pi^{-} - pion (anchor);#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_PinK = new TH1D("hDeltaT_global_anchor_PinK", "#Delta t #pi^{-} - kaon (anchor);#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_Pinp = new TH1D("hDeltaT_global_anchor_Pinp", "#Delta t #pi^{-} - proton (anchor);#Delta t [ns];# events", 80, -4.0, 4.0);

    TH1D* hDeltaT_global_anchor_Pppi_high_pt = new TH1D("hDeltaT_global_anchor_Pppi_high_pt", "#Delta t p - pion (anchor) (high p_{T});#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_PpK_high_pt = new TH1D("hDeltaT_global_anchor_PpK_high_pt", "#Delta t p - kaon (anchor) (high p_{T});#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_Ppp_high_pt = new TH1D("hDeltaT_global_anchor_Ppp_high_pt", "#Delta t p - proton (anchor) (high p_{T});#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_Pnpi_high_pt = new TH1D("hDeltaT_global_anchor_Pnpi_high_pt", "#Delta t #bar{p} - pion (anchor) (high p_{T});#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_PnK_high_pt = new TH1D("hDeltaT_global_anchor_PnK_high_pt", "#Delta t #bar{p} - kaon (anchor) (high p_{T});#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_Pnp_high_pt = new TH1D("hDeltaT_global_anchor_Pnp_high_pt", "#Delta t #bar{p} - proton (anchor) (high p_{T});#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_Pppi_mid_pt = new TH1D("hDeltaT_global_anchor_Pppi_mid_pt", "#Delta t p - pion (anchor) (mid p_{T});#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_PpK_mid_pt = new TH1D("hDeltaT_global_anchor_PpK_mid_pt", "#Delta t p - kaon (anchor) (mid p_{T});#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_Ppp_mid_pt = new TH1D("hDeltaT_global_anchor_Ppp_mid_pt", "#Delta t p - proton (anchor) (mid p_{T});#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_Pnpi_mid_pt = new TH1D("hDeltaT_global_anchor_Pnpi_mid_pt", "#Delta t #bar{p} - pion (anchor) (mid p_{T});#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_PnK_mid_pt = new TH1D("hDeltaT_global_anchor_PnK_mid_pt", "#Delta t #bar{p} - kaon (anchor) (mid p_{T});#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_Pnp_mid_pt = new TH1D("hDeltaT_global_anchor_Pnp_mid_pt", "#Delta t #bar{p} - proton (anchor) (mid p_{T});#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_Pppi = new TH1D("hDeltaT_global_anchor_Pppi", "#Delta t p - pion (anchor);#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_PpK = new TH1D("hDeltaT_global_anchor_PpK", "#Delta t p - kaon (anchor);#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_Ppp = new TH1D("hDeltaT_global_anchor_Ppp", "#Delta t p - proton (anchor);#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_Pnpi = new TH1D("hDeltaT_global_anchor_Pnpi", "#Delta t #bar{p} - pion (anchor);#Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_PnK = new TH1D("hDeltaT_global_anchor_PnK", "#Delta t #bar{p} - kaon (anchor);# Delta t [ns];# events", 80, -4.0, 4.0);
    TH1D* hDeltaT_global_anchor_Pnp = new TH1D("hDeltaT_global_anchor_Pnp", "#Delta t #bar{p} - proton (anchor);# Delta t [ns];# events", 80, -4.0, 4.0);

    TH1D* hTPClength = new TH1D("hTPClength", "track length;length [cm];# events", 500, -500.0, 500.0);
    TH1D* hTPCdx = new TH1D("hTPCdx", "TPC track dx;dx [cm];# events", 120, -20.0, 20.0);
    TH2D* hdxVsP = new TH2D("hdxVsP", "dx vs p; p [GeV/c]; dx [cm]", 100, 0, 5, 100, 0, 20.0);

    const std::vector<double> dx_targets = {3.3, 3.4, 3.5, 3.6, 3.7};
    const int nDx = dx_targets.size();
    TH2D* hRatioVsP_pion[nDx];

    for (int i = 0; i < nDx; ++i) {
        TString name  = Form("hRatioVsP_pion_dx_%.2f", dx_targets[i]);
        TString title = Form("dE/dx Ratio (Pion) for dx = %.2f cm;p [GeV/c]; (dE/dx)_{meas} / (dE/dx)_{Bichsel}", dx_targets[i]);
        
        hRatioVsP_pion[i] = new TH2D(name, title, 200, 0.15, 5.0, 200, 0.5, 1.5);
    }

    TH2D* hRatioVsP_kaon[nDx];
    TH2D* hRatioVsP_proton[nDx];

    for (int i = 0; i < nDx; ++i) {
        hRatioVsP_kaon[i] = new TH2D(
            Form("hRatioVsP_kaon_dx_%.2f", dx_targets[i]),
            Form("dE/dx Ratio (Kaon) for dx = %.2f cm;p [GeV/c]; (dE/dx)_{meas} / (dE/dx)_{Bichsel}", dx_targets[i]),
            200, 0.15, 5.0, 200, 0.5, 1.5
        );

        hRatioVsP_proton[i] = new TH2D(
            Form("hRatioVsP_proton_dx_%.2f", dx_targets[i]),
            Form("dE/dx Ratio (Proton) for dx = %.2f cm;p [GeV/c]; (dE/dx)_{meas} / (dE/dx)_{Bichsel}", dx_targets[i]),
            200, 0.15, 5.0, 200, 0.5, 1.5
        );
    }

    TH2D* hRatioTOFVsP_pion[nDx];
    TH2D* hRatioTOFVsP_kaon[nDx];
    TH2D* hRatioTOFVsP_proton[nDx];

    for (int i = 0; i < nDx; ++i) {
        hRatioTOFVsP_pion[i] = new TH2D(
            Form("hRatioTOFVsP_pion_dx_%.2f", dx_targets[i]),
            Form("dE/dx Ratio TOF (Pion) for dx = %.2f cm;p [GeV/c]; (dE/dx)_{meas} / (dE/dx)_{Bichsel}", dx_targets[i]),
            200, 0.15, 5.0, 200, 0.5, 1.5
        );

        hRatioTOFVsP_kaon[i] = new TH2D(
            Form("hRatioTOFVsP_kaon_dx_%.2f", dx_targets[i]),
            Form("dE/dx Ratio TOF (Kaon) for dx = %.2f cm;p [GeV/c]; (dE/dx)_{meas} / (dE/dx)_{Bichsel}", dx_targets[i]),
            200, 0.15, 5.0, 200, 0.5, 1.5
        );

        hRatioTOFVsP_proton[i] = new TH2D(
            Form("hRatioTOFVsP_proton_dx_%.2f", dx_targets[i]),
            Form("dE/dx Ratio TOF (Proton) for dx = %.2f cm;p [GeV/c]; (dE/dx)_{meas} / (dE/dx)_{Bichsel}", dx_targets[i]),
            200, 0.15, 5.0, 200, 0.5, 1.5
        );
    }



    std::vector<std::string> particles = {"Kp", "Kn", "Pip", "Pin", "Pp", "Pn"};
    std::vector<std::string> particle_titles = {"K^{+}", "K^{-}", "#pi^{+}", "#pi^{-}", "p", "#bar{p}"};

    std::vector<std::string> partners = {"pi", "K", "p"};
    std::vector<std::string> partner_titles = {"pion", "kaon", "proton"};

    std::vector<std::string> mom_bins = {"_p1", "_p2", "_p3", "_p4", "_p5", "_p6", "_p7", "_p8", "_p9", "_p10"};
    std::vector<std::string> mom_titles = {" (p1)", " (p2)", " (p3)", " (p4)", " (p5)", " (p6)", " (p7)", " (p8)", " (p9)", " (p10)"};

    std::vector<std::vector<std::vector<TH1D*>>> hDeltaT(
        particles.size(), 
        std::vector<std::vector<TH1D*>>(partners.size(), std::vector<TH1D*>(mom_bins.size()))
    );

    for (size_t i = 0; i < particles.size(); ++i) {
        for (size_t j = 0; j < partners.size(); ++j) {
            for (size_t k = 0; k < mom_bins.size(); ++k) {
                
                // Generuje np. "hDeltaT_global_anchor_Kppi_p1"
                std::string name = "hDeltaT_global_anchor_" + particles[i] + partners[j] + mom_bins[k];
                
                // Generuje np. "#Delta t K^{+} - pion (anchor) (p1);#Delta t [ns];# events"
                std::string title = "#Delta t " + particle_titles[i] + " - " + partner_titles[j] + " (anchor)" + mom_titles[k] + ";#Delta t [ns];# events";
                
                hDeltaT[i][j][k] = new TH1D(name.c_str(), title.c_str(), 80, -4.0, 4.0);
            }
        }
    }

    TH1D* hDeltaT_lambda_proton_pion=new TH1D("hDeltaT_lambda_proton_pion","#Delta t lambda proton-pion hypothesis;#Delta t [ns];# events",4000,-20.0,20.0);

    double massPion = 0.13957061;
    double massKaon = 0.49367;
    double massProton = 0.93827;
    double massLambda = 1.115638;

    int EventsPassedSelection=0;

    std::vector<int> goodTrackIdx;
    std::vector<int> proton_Idx;
    std::vector<int> pion_Idx;
    std::vector<int> kaon_Idx;

    Bool_t TreeIsEast = false;
    Double_t xi_proton = -999.0;
    Int_t nCharged = 0;

    std::vector<double> track_p;
    std::vector<double> track_pt;
    std::vector<double> track_eta;
    Double_t event_vz=-999.0;
    
    std::vector<double> track_nSigmaPion;
    std::vector<double> track_nSigmaKaon;
    std::vector<double> track_nSigmaProton;
    
    std::vector<double> track_dTofPion;
    std::vector<double> track_dTofKaon;
    std::vector<double> track_dTofProton;

    std::vector<int> track_charge;
    

    TTree *outTree = new TTree("outTree", "Tree");

    outTree->Branch("TreeIsEast", &TreeIsEast, "TreeIsEast/O");
    outTree->Branch("xi_proton", &xi_proton, "xi_proton/D");
    outTree->Branch("nCharged", &nCharged, "nCharged/I");

    outTree->Branch("track_p", &track_p);
    outTree->Branch("track_pt", &track_pt);
    outTree->Branch("track_eta", &track_eta);
    outTree->Branch("track_charge", &track_charge);
    outTree->Branch("event_vz", &event_vz, "event_vz/D");

    outTree->Branch("track_nSigmaPion", &track_nSigmaPion);
    outTree->Branch("track_nSigmaKaon", &track_nSigmaKaon);
    outTree->Branch("track_nSigmaProton", &track_nSigmaProton);

    outTree->Branch("track_dTofPion", &track_dTofPion);
    outTree->Branch("track_dTofKaon", &track_dTofKaon);
    outTree->Branch("track_dTofProton", &track_dTofProton);

    for (Long64_t i = 0; i < chain->GetEntries(); ++i) {
        chain->GetEntry(i);
        if( upcEvt->getNumberOfVertices() != 1) continue;
        if (rpEvt->getNumberOfTracks() != 1 ) continue;

        goodTrackIdx.clear();
        proton_Idx.clear();
        pion_Idx.clear();
        kaon_Idx.clear();


        track_p.clear();
        track_pt.clear();
        track_eta.clear();
        track_nSigmaPion.clear();
        track_nSigmaKaon.clear();
        track_nSigmaProton.clear();
        track_dTofPion.clear();
        track_dTofKaon.clear();
        track_dTofProton.clear();
        track_charge.clear();
        //   cout << upcEvt->getEventNumber() << " " 
        //   cout << upcEvt->getNPrimVertices() <<  endl; //1, 2, 3

        StUPCRpsTrack *proton = rpEvt->getTrack(0);

        bool isEast = (proton->branch() < 2);
        bool isWest = (proton->branch() > 1);
        xi_proton = proton->xi(254.867);

        if (proton->branch() < 2)  { //east
            HistXiProtonEast->Fill(proton->xi(254.867));
            HistLogXiProtonEast->Fill(log(proton->xi(254.867)));
            HistPtProtonEast->Fill(proton->pt());
            HistEtaProtonEast->Fill(proton->eta());
        }
        if (proton->branch() > 1)  { //west
            HistXiProtonWest->Fill(proton->xi(254.867));
            HistLogXiProtonWest->Fill(log(proton->xi(254.867)));
            HistPtProtonWest->Fill(proton->pt());
            HistEtaProtonWest->Fill(proton->eta());
        }
        int b=getXiBin(xi_proton);

        //if (xi_proton<1e-7) continue;

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

        int MultiplicityProtonX_mid_pt=0;
        int MultiplicityProtonWestX_mid_pt=0;
        int MultiplicityProtonEastX_mid_pt=0;
        int MultiplicityNProtonX_mid_pt=0;
        int MultiplicityNProtonWestX_mid_pt=0;
        int MultiplicityNProtonEastX_mid_pt=0;

        int MultiplicityProtonX_high_pt=0;
        int MultiplicityProtonWestX_high_pt=0;
        int MultiplicityProtonEastX_high_pt=0;
        int MultiplicityNProtonX_high_pt=0;
        int MultiplicityNProtonWestX_high_pt=0;
        int MultiplicityNProtonEastX_high_pt=0;

        //vz in cm
        double vz = upcEvt->getVertex(0)->getPosZ();
        hVz_check->Fill(vz);
        event_vz = vz;
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

        nCharged = NumOfGoodPrimaryTracks;

        if (NumOfGoodPrimaryTracks<2) continue;

        for (int j = 0; j < upcEvt->getNumberOfTracks(); j++) {
            if( upcEvt->getTrack(j)->getNhits()>20  
                && upcEvt->getTrack(j)->getFlag(StUPCTrack::kTof) ) {

                if  ( upcEvt->getTrack(j)->getPt() > 0.2 &&
                    ( fabs(upcEvt->getTrack(j)->getEta()) < 0.9 &&
                    upcEvt->getTrack(j)->getEta() < -(vz/250.0) + 0.9 &&
                    upcEvt->getTrack(j)->getEta() > -(vz/250.0) - 0.9 ) ){ // fiducial cuts
                    if(isWest) {
                        hNumOfParticlesPassSelectionWest->Fill(1);
                    } else {
                        hNumOfParticlesPassSelectionEast->Fill(1);
                    }
                }
            }
        }

        int proton_high_p_candidate=0;
        int proton_mid_p_candidate=0;


        for (int j = 0; j < upcEvt->getNumberOfTracks(); j++) {
            if( upcEvt->getTrack(j)->getNhits()>20  
                && upcEvt->getTrack(j)->getFlag(StUPCTrack::kTof) 
                && upcEvt->getTrack(j)->getFlag(StUPCTrack::kPrimary)) {

                if  ( upcEvt->getTrack(j)->getPt() > 0.2 &&  // pT, eta and vz/eta cut
                    ( fabs(upcEvt->getTrack(j)->getEta()) < 0.9 &&
                    upcEvt->getTrack(j)->getEta() < -(vz/250.0) + 0.9 &&
                    upcEvt->getTrack(j)->getEta() > -(vz/250.0) - 0.9 ) ){
                    NumOfPrimaryTracksToF++;
                    TVector3 momentum;
                    upcEvt->getTrack(j)->getMomentum(momentum);
                    double p = momentum.Mag();

                    goodTrackIdx.push_back(j);

                    if (isEast) {
                        HistPtTracksEastCut->Fill(upcEvt->getTrack(j)->getPt());
                        HistEtaTracksEastCut->Fill(upcEvt->getTrack(j)->getEta());
                    }
                    if (isWest) {
                        HistPtTracksWestCut->Fill(upcEvt->getTrack(j)->getPt());
                        HistEtaTracksWestCut->Fill(upcEvt->getTrack(j)->getEta());
                    }

                    double pt  = upcEvt->getTrack(j)->getPt();
                    double eta = upcEvt->getTrack(j)->getEta();
                    double dca = upcEvt->getTrack(j)->getDcaXY();

                    if((fabs(upcEvt->getTrack(j)->getNSigmasTPCPion())<3) && (upcEvt->getTrack(j)->getNhitsDEdx()>=15)) {
                        double rapidity = getRapidity(pt, eta, upcEvt->getTrack(j)->getPhi(), massPion);
                        if(upcEvt->getTrack(j)->getCharge()>0){ //pion+
                            pion_Idx.push_back(j);
                        } 
                        if(upcEvt->getTrack(j)->getCharge()<0){ //pion-
                            pion_Idx.push_back(j);
                        }
                        if(fabs(upcEvt->getTrack(j)->getNSigmasTPCPion())<0.25){
                            double dedx_meas=upcEvt->getTrack(j)->getDEdxSignal()*1e6;
                            for (int i = 0; i < nDx; ++i) {
                                double dx_val = dx_targets[i];
                                double dedx_theo = GetTheoreticaldEdx(p, massPion, dx_val); 
                                double ratio = dedx_meas/dedx_theo;
                                hRatioVsP_pion[i]->Fill(p, ratio);
                            }
                        }
                    }


                    if((fabs(upcEvt->getTrack(j)->getNSigmasTPCKaon()) < 3) && (upcEvt->getTrack(j)->getNhitsDEdx() >= 15)) {
                        double rapidity = getRapidity(pt, eta, upcEvt->getTrack(j)->getPhi(), massKaon);

                        if(upcEvt->getTrack(j)->getCharge() > 0){ // kaon+
                            kaon_Idx.push_back(j);
                        } 
                        if(upcEvt->getTrack(j)->getCharge() < 0){ // kaon-
                            kaon_Idx.push_back(j);
                        }
                        if(fabs(upcEvt->getTrack(j)->getNSigmasTPCKaon())<0.25){
                            double dedx_meas=upcEvt->getTrack(j)->getDEdxSignal()*1e6;
                            for (int i = 0; i < nDx; ++i) {
                                double dx_val = dx_targets[i];
                                double dedx_theo = GetTheoreticaldEdx(p, massKaon, dx_val); 
                                double ratio = dedx_meas/dedx_theo;
                                hRatioVsP_kaon[i]->Fill(p, ratio);
                            }
                        }
                    }


                    if((fabs(upcEvt->getTrack(j)->getNSigmasTPCProton())<3) && (upcEvt->getTrack(j)->getNhitsDEdx()>=15)) { //proton further selection based on nsigma and p
                        double rapidity = getRapidity(pt, eta, upcEvt->getTrack(j)->getPhi(), massProton);
                        double distance_from_beam = rapidity - y_beam;
                        double distance_from_edge = y_edge-rapidity;
                        if(upcEvt->getTrack(j)->getCharge()>0){ //proton
                            proton_Idx.push_back(j);
                        } 
                        if(upcEvt->getTrack(j)->getCharge()<0){ //antiproton
                            proton_Idx.push_back(j);
                        }
                        if(fabs(upcEvt->getTrack(j)->getNSigmasTPCProton())<0.25){
                            double dedx_meas=upcEvt->getTrack(j)->getDEdxSignal()*1e6;
                            for (int i = 0; i < nDx; ++i) {
                                double dx_val = dx_targets[i];
                                double dedx_theo = GetTheoreticaldEdx(p, massProton, dx_val); 
                                double ratio = dedx_meas/dedx_theo;
                                hRatioVsP_proton[i]->Fill(p, ratio);
                            }
                        }
                    }
                }
                if (fabs(upcEvt->getTrack(j)->getEta()) < 0.9 &&  //tracks && nsigmas 2D plots
                    upcEvt->getTrack(j)->getEta() < -(vz/250.0) + 0.9 &&
                    upcEvt->getTrack(j)->getEta() > -(vz/250.0) - 0.9) {
                    TVector3 momentum;
                    upcEvt->getTrack(j)->getMomentum(momentum);
                    double p = momentum.Mag();
                    if (isEast) HistPtTracksEast->Fill(upcEvt->getTrack(j)->getPt());
                    if (isWest) HistPtTracksWest->Fill(upcEvt->getTrack(j)->getPt());
                    if (upcEvt->getTrack(j)->getNhitsDEdx()>=15) {
                        if(upcEvt->getTrack(j)->getCharge()>0) { //positive charge
                            hNSigmaPiPlus->Fill(p, upcEvt->getTrack(j)->getNSigmasTPCPion());
                            hNSigmaKPlus->Fill(p, upcEvt->getTrack(j)->getNSigmasTPCKaon());
                            hNSigmaPPlus->Fill(p, upcEvt->getTrack(j)->getNSigmasTPCProton());
                        } else { //negative charge
                            hNSigmaPiMinus->Fill(p, upcEvt->getTrack(j)->getNSigmasTPCPion());
                            hNSigmaKMinus->Fill(p, upcEvt->getTrack(j)->getNSigmasTPCKaon());
                            hNSigmaPMinus->Fill(p, upcEvt->getTrack(j)->getNSigmasTPCProton());
                        }
                        TVector3 momentum_01;
                        upcEvt->getTrack(j)->getMomentum(momentum_01);
                        double pq = (momentum_01.Mag())*(upcEvt->getTrack(j)->getCharge());
                        hdEdx->Fill(pq, upcEvt->getTrack(j)->getDEdxSignal()*1e6); //GeV originally, keV now

                        StPicoHelix helix(upcEvt->getTrack(j)->getCurvature(), upcEvt->getTrack(j)->getDipAngle(),upcEvt->getTrack(j)->getPhase(), upcEvt->getTrack(j)->getOrigin(), -upcEvt->getTrack(j)->getCharge()/abs(upcEvt->getTrack(j)->getCharge()));
                        auto pathTof = helix.pathLength(222.0 * centimeter);
                        auto pathOut = helix.pathLength(200.0 * centimeter);
                        auto pathIn  = helix.pathLength(50.0  * centimeter);

                        double tofData = upcEvt->getTrack(j)->getTofPathLength();

                        // 3. Sprawdzamy, które rozwiązanie z pary (.first czy .second) jest bliższe realnemu TOF
                        double diff1 = std::abs(pathTof.first  - tofData);
                        double diff2 = std::abs(pathTof.second - tofData);

                        double sOut = (diff1 < diff2) ? pathOut.first : pathOut.second;
                        double sIn  = (diff1 < diff2) ? pathIn.first  : pathIn.second;

                        // 4. Wynikowa długość w TPC oraz dx
                        double TPClength = std::abs(sOut - sIn);
                        double dx = TPClength / upcEvt->getTrack(j)->getNhitsDEdx();

                        hdxVsP->Fill(momentum_01.Mag(),dx);

                        hTPClength->Fill(TPClength);
                        hTPCdx->Fill(dx);
                    }
                } 
                if (upcEvt->getTrack(j)->getPt()>0.2) {
                    if (isEast) HistEtaTracksEast->Fill(upcEvt->getTrack(j)->getEta());
                    if (isWest) HistEtaTracksWest->Fill(upcEvt->getTrack(j)->getEta()); 
                }
            }

            if(upcEvt->getTrack(j)->getFlag(StUPCTrack::kTof) //lambda selection - first daughter track
                && upcEvt->getTrack(j)->getFlag(StUPCTrack::kV0)
                && upcEvt->getTrack(j)->getNhits()>20 
                && upcEvt->getTrack(j)->getPt()>0.2
                && fabs(upcEvt->getTrack(j)->getEta())<0.9) {
                        
                StUPCTrack* track1=upcEvt->getTrack(j);

                for(size_t k=j+1; k<upcEvt->getNumberOfTracks(); k++) { //lambda selection - second daughter track

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
                            
                            double t0_proton = track2->getT0(massProton);
                            double t0_pion = track1->getT0(massPion);
                            double delta_t_lambda = t0_proton-t0_pion;
                            hDeltaT_lambda_proton_pion->Fill(delta_t_lambda);


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
                                    if (LambdaCut(L01, 'n') && std::cos(L01.pointingAngle())>0.925 && std::cos(L01.pointingAngle())<0.99 && L01.decayLength()>5) hLambda_Mass_Background_PointingAngle_W->Fill(L01.m());
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

                            double t0_proton = track1->getT0(massProton);
                            double t0_pion = track2->getT0(massPion);
                            double delta_t_lambda = t0_proton-t0_pion;
                            hDeltaT_lambda_proton_pion->Fill(delta_t_lambda);

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
        } // end of loop for tracks in an eventd

        int global_ref_idx_dedx = -1;
        double maksymalna_separacja = -1.0;
        double wyznaczony_t0 = 0.0;
        bool znaleziono_anchor = false;
        int best_pid = -1;

        for (int i_trk : goodTrackIdx) { //for loop for finding the best track for PID reference (anchor track)
            StUPCTrack *trk = upcEvt->getTrack(i_trk);
            if (!trk) continue;
        
            double chi2_pi = pow(trk->getNSigmasTPCPion(), 2);
            double chi2_K  = pow(trk->getNSigmasTPCKaon(), 2);
            double chi2_p  = pow(trk->getNSigmasTPCProton(), 2);
            std::pair<double, int> chi2_list[3] = {
                {chi2_pi, 0},
                {chi2_K,  1},
                {chi2_p,  2}
            };
            std::sort(chi2_list, chi2_list + 3);
        
            double current_best_chi2 = chi2_list[0].first;
            double delta_chi2 = chi2_list[1].first - chi2_list[0].first;
            if (delta_chi2 > maksymalna_separacja) {
                maksymalna_separacja = delta_chi2;
                global_ref_idx_dedx = i_trk;
                znaleziono_anchor = true;
        
                best_pid = chi2_list[0].second;
                if (best_pid == 0)      wyznaczony_t0 = trk->getT0(massPion);
                else if (best_pid == 1) wyznaczony_t0 = trk->getT0(massKaon);
                else if (best_pid == 2) wyznaczony_t0 = trk->getT0(massProton);
            }
        }

        if(isWest) hEventsPassedSelection->Fill(1);
        if(isEast) hEventsPassedSelection->Fill(2);
        TreeIsEast = isEast;

        for (int i_trk : goodTrackIdx) { //tree filling with nsigmas from TOF and TPC for all good tracks in the event
            StUPCTrack *trk = upcEvt->getTrack(i_trk);
            if (!trk) continue;
            if ((upcEvt->getTrack(i_trk)->getNhitsDEdx()<15)) continue;

            TVector3 momentum;
            trk->getMomentum(momentum);
            double p = momentum.Mag();

            track_p.push_back(p);
            track_pt.push_back(trk->getPt());
            track_eta.push_back(trk->getEta());
            track_charge.push_back(trk->getCharge());
            track_nSigmaPion.push_back(trk->getNSigmasTPCPion());
            track_nSigmaKaon.push_back(trk->getNSigmasTPCKaon());
            track_nSigmaProton.push_back(trk->getNSigmasTPCProton());
            
            double tof_limit_pi = 1.0;
            double tof_limit_K  = 1.0;
            double tof_limit_p  = 1.0;

            if (best_pid == 0) { // Ref = PION (_pions)
                tof_limit_pi = TOFPIONS_pions;
                tof_limit_K  = TOFKAONS_pions;
                tof_limit_p  = TOF_pions;
            } else if (best_pid == 1) { // Ref = KAON (_kaons)
                tof_limit_pi = TOFPIONS_kaons;
                tof_limit_K  = TOFKAONS_kaons;
                tof_limit_p  = TOF_kaons;
            } else if (best_pid == 2) { 
                tof_limit_pi = TOFPIONS_protons;
                tof_limit_K  = TOFKAONS_protons;
                tof_limit_p  = TOF_protons;
            }
            double dt_pi = trk->getT0(massPion) - wyznaczony_t0;
            track_dTofPion.push_back(dt_pi / (tof_limit_pi / 3.0));

            double dt_K = trk->getT0(massKaon) - wyznaczony_t0;
            track_dTofKaon.push_back(dt_K / (tof_limit_K / 3.0));

            double dt_p = trk->getT0(massProton) - wyznaczony_t0;
            track_dTofProton.push_back(dt_p / (tof_limit_p / 3.0));
        }
        nCharged = track_p.size();
        outTree->Fill();
        if (isWest && isEast) {
            std::cout << " [PROBLEM] Event ma ZARÓWNO isWest=1 JAK I isEast=1! " << std::endl;
        }

        for (int i_k : kaon_Idx) {
            TVector3 mom_k;
            upcEvt->getTrack(i_k)->getMomentum(mom_k);
            double p_k = mom_k.Mag();

            enum MomentumCat { LOW_P, MID_P, HIGH_P };
            MomentumCat pCat = LOW_P;
            if (p_k >= 0.5 && p_k <= 0.9) pCat = MID_P;
            else if (p_k > 0.9)           pCat = HIGH_P;

            bool isPositive = (upcEvt->getTrack(i_k)->getCharge() > 0);
            if (znaleziono_anchor && global_ref_idx_dedx != i_k) {
                double t0_kaon_hip = upcEvt->getTrack(i_k)->getT0(massKaon);
                double delta_t_global = t0_kaon_hip - wyznaczony_t0;

                int k_bin = getMomentumBin(p_k);
                if (k_bin != -1) {
                    int i_part = -1;
                    if      (particles[0] == "Kp"  && isPositive) i_part = 0;
                    else if (particles[1] == "Kn"  && !isPositive) i_part = 1;

                    int j_part = -1;
                    double chi2_pi = pow(upcEvt->getTrack(global_ref_idx_dedx)->getNSigmasTPCPion(), 2);
                    double chi2_K  = pow(upcEvt->getTrack(global_ref_idx_dedx)->getNSigmasTPCKaon(), 2);
                    double chi2_p  = pow(upcEvt->getTrack(global_ref_idx_dedx)->getNSigmasTPCProton(), 2);

                    if (chi2_pi <= chi2_K && chi2_pi <= chi2_p) j_part = 0;
                    else if (chi2_K <= chi2_pi && chi2_K <= chi2_p) j_part = 1;
                    else j_part = 2;

                    if (i_part != -1 && j_part != -1) {
                        hDeltaT[i_part][j_part][k_bin]->Fill(delta_t_global);

                        double min_delta=delta_t_global;
                        double tof_limit = TOFKAONS_pions;
                        if (j_part == 1) tof_limit = TOFKAONS_kaons;
                        else if (j_part == 2) tof_limit = TOFKAONS_protons;

                        if (fabs(min_delta) / (tof_limit / 3.0) < 0.25 && (fabs(upcEvt->getTrack(i_k)->getNSigmasTPCKaon())<0.25)) {
                            double dedx_meas=upcEvt->getTrack(i_k)->getDEdxSignal()*1e6;
                            TVector3 mom;
                            upcEvt->getTrack(i_k)->getMomentum(mom);
                            double p=mom.Mag();
                            for (int i = 0; i < nDx; ++i) {
                                double dx_val = dx_targets[i];
                                double dedx_theo = GetTheoreticaldEdx(p, massKaon, dx_val); 
                                double ratio = dedx_meas/dedx_theo;
                                hRatioTOFVsP_kaon[i]->Fill(p, ratio);
                            }
                        }
                    }
                }

                double chi2_pi = pow(upcEvt->getTrack(global_ref_idx_dedx)->getNSigmasTPCPion(), 2);
                double chi2_K  = pow(upcEvt->getTrack(global_ref_idx_dedx)->getNSigmasTPCKaon(), 2);
                double chi2_p  = pow(upcEvt->getTrack(global_ref_idx_dedx)->getNSigmasTPCProton(), 2);

                TH1* hDeltaT_global = nullptr;

                if (chi2_pi <= chi2_K && chi2_pi <= chi2_p) {
                    if (isPositive) {
                        if (pCat == LOW_P)       hDeltaT_global = hDeltaT_global_anchor_Kppi;
                        else if (pCat == MID_P)  hDeltaT_global = hDeltaT_global_anchor_Kppi_mid_pt;
                        else                     hDeltaT_global = hDeltaT_global_anchor_Kppi_high_pt;
                    } else {
                        if (pCat == LOW_P)       hDeltaT_global = hDeltaT_global_anchor_Knpi;
                        else if (pCat == MID_P)  hDeltaT_global = hDeltaT_global_anchor_Knpi_mid_pt;
                        else                     hDeltaT_global = hDeltaT_global_anchor_Knpi_high_pt;
                    }
                } 
                else if (chi2_K <= chi2_pi && chi2_K <= chi2_p) {
                    if (isPositive) {
                        if (pCat == LOW_P)       hDeltaT_global = hDeltaT_global_anchor_KpK;
                        else if (pCat == MID_P)  hDeltaT_global = hDeltaT_global_anchor_KpK_mid_pt;
                        else                     hDeltaT_global = hDeltaT_global_anchor_KpK_high_pt;
                    } else {
                        if (pCat == LOW_P)       hDeltaT_global = hDeltaT_global_anchor_KnK;
                        else if (pCat == MID_P)  hDeltaT_global = hDeltaT_global_anchor_KnK_mid_pt;
                        else                     hDeltaT_global = hDeltaT_global_anchor_KnK_high_pt;
                    }
                } 
                else {
                    if (isPositive) {
                        if (pCat == LOW_P)       hDeltaT_global = hDeltaT_global_anchor_Kpp;
                        else if (pCat == MID_P)  hDeltaT_global = hDeltaT_global_anchor_Kpp_mid_pt;
                        else                     hDeltaT_global = hDeltaT_global_anchor_Kpp_high_pt;
                    } else {
                        if (pCat == LOW_P)       hDeltaT_global = hDeltaT_global_anchor_Knp;
                        else if (pCat == MID_P)  hDeltaT_global = hDeltaT_global_anchor_Knp_mid_pt;
                        else                     hDeltaT_global = hDeltaT_global_anchor_Knp_high_pt;
                    }
                }

                if (hDeltaT_global) hDeltaT_global->Fill(delta_t_global);
            }

            int ref_idx = -1;
            double min_p = 9999.0;

            for (int i_trk : goodTrackIdx) {
                if (i_trk == i_k) continue;
                TVector3 mom;
                upcEvt->getTrack(i_trk)->getMomentum(mom);

                if (mom.Mag() < min_p) {
                    min_p = mom.Mag();
                    ref_idx = i_trk;
                }
            }

            if (ref_idx >= 0) {
                double t0_kaon_hip = upcEvt->getTrack(i_k)->getT0(massKaon);
                double t0_pion     = upcEvt->getTrack(ref_idx)->getT0(massPion);
                double t0_kaon     = upcEvt->getTrack(ref_idx)->getT0(massKaon);
                double t0_proton   = upcEvt->getTrack(ref_idx)->getT0(massProton);

                double delta_t1 = t0_kaon_hip - t0_pion;
                double delta_t2 = t0_kaon_hip - t0_kaon;
                double delta_t5 = t0_kaon_hip - t0_proton;

                TH2* hT1vsT2_charge = nullptr;
                TH2* hT1vsT2_all    = nullptr;
                TH2* hT1vsT5_charge = nullptr;

                if (isPositive) {
                    if (pCat == LOW_P) {
                        hT1vsT2_charge = hDeltaT1vsT2_kaons_p;
                        hT1vsT2_all    = hDeltaT1vsT2_kaons;
                        hT1vsT5_charge = hDeltaT1vsT5_kaons_p;
                    } else if (pCat == MID_P) {
                        hT1vsT2_charge = hDeltaT1vsT2_kaons_p_mid_p;
                        hT1vsT2_all    = hDeltaT1vsT2_kaons_mid_p;
                        hT1vsT5_charge = hDeltaT1vsT5_kaons_p_mid_p;
                    } else {
                        hT1vsT2_charge = hDeltaT1vsT2_kaons_p_high_p;
                        hT1vsT2_all    = hDeltaT1vsT2_kaons_high_p;
                        hT1vsT5_charge = hDeltaT1vsT5_kaons_p_high_p;
                    }
                } else {
                    if (pCat == LOW_P) {
                        hT1vsT2_charge = hDeltaT1vsT2_kaons_n;
                        hT1vsT2_all    = hDeltaT1vsT2_kaons;
                        hT1vsT5_charge = hDeltaT1vsT5_kaons_n;
                    } else if (pCat == MID_P) {
                        hT1vsT2_charge = hDeltaT1vsT2_kaons_n_mid_p;
                        hT1vsT2_all    = hDeltaT1vsT2_kaons_mid_p;
                        hT1vsT5_charge = hDeltaT1vsT5_kaons_n_mid_p;
                    } else {
                        hT1vsT2_charge = hDeltaT1vsT2_kaons_n_high_p;
                        hT1vsT2_all    = hDeltaT1vsT2_kaons_high_p;
                        hT1vsT5_charge = hDeltaT1vsT5_kaons_n_high_p;
                    }
                }

                hT1vsT2_charge->Fill(delta_t1, delta_t2);
                hT1vsT2_all->Fill(delta_t1, delta_t2);
                hT1vsT5_charge->Fill(delta_t1, delta_t5);

                if (!(fabs(delta_t1) < TOFKAONS_pions || fabs(delta_t2) < TOFKAONS_kaons || fabs(delta_t5) < TOFKAONS_protons)) continue;

                double pt       = upcEvt->getTrack(i_k)->getPt();
                double eta      = upcEvt->getTrack(i_k)->getEta();
                double rapidity = getRapidity(pt, eta, upcEvt->getTrack(i_k)->getPhi(), massKaon);

                TH1* hRap_EW = nullptr;
                TH1* hRap_C  = nullptr;

                if (isEast) {
                    if (isPositive) {
                        if (pCat == LOW_P) {
                            if (b == 0) { hRap_EW = hRapidity_xi0_kaon_plus_E; hRap_C = hRapidity_xi0_kaon_plus_C; }
                            if (b == 1) { hRap_EW = hRapidity_xi1_kaon_plus_E; hRap_C = hRapidity_xi1_kaon_plus_C; }
                            if (b == 2) { hRap_EW = hRapidity_xi2_kaon_plus_E; hRap_C = hRapidity_xi2_kaon_plus_C; }
                            if (b == 3) { hRap_EW = hRapidity_xi3_kaon_plus_E; hRap_C = hRapidity_xi3_kaon_plus_C; }
                            if (b == 4) { hRap_EW = hRapidity_xi4_kaon_plus_E; hRap_C = hRapidity_xi4_kaon_plus_C; }
                        } else if (pCat == MID_P) {
                            if (b == 0) { hRap_EW = hRapidity_xi0_kaon_plus_E_mid_pt; hRap_C = hRapidity_xi0_kaon_plus_C_mid_pt; }
                            if (b == 1) { hRap_EW = hRapidity_xi1_kaon_plus_E_mid_pt; hRap_C = hRapidity_xi1_kaon_plus_C_mid_pt; }
                            if (b == 2) { hRap_EW = hRapidity_xi2_kaon_plus_E_mid_pt; hRap_C = hRapidity_xi2_kaon_plus_C_mid_pt; }
                            if (b == 3) { hRap_EW = hRapidity_xi3_kaon_plus_E_mid_pt; hRap_C = hRapidity_xi3_kaon_plus_C_mid_pt; }
                            if (b == 4) { hRap_EW = hRapidity_xi4_kaon_plus_E_mid_pt; hRap_C = hRapidity_xi4_kaon_plus_C_mid_pt; }
                        } else {
                            if (b == 0) { hRap_EW = hRapidity_xi0_kaon_plus_E_high_pt; hRap_C = hRapidity_xi0_kaon_plus_C_high_pt; }
                            if (b == 1) { hRap_EW = hRapidity_xi1_kaon_plus_E_high_pt; hRap_C = hRapidity_xi1_kaon_plus_C_high_pt; }
                            if (b == 2) { hRap_EW = hRapidity_xi2_kaon_plus_E_high_pt; hRap_C = hRapidity_xi2_kaon_plus_C_high_pt; }
                            if (b == 3) { hRap_EW = hRapidity_xi3_kaon_plus_E_high_pt; hRap_C = hRapidity_xi3_kaon_plus_C_high_pt; }
                            if (b == 4) { hRap_EW = hRapidity_xi4_kaon_plus_E_high_pt; hRap_C = hRapidity_xi4_kaon_plus_C_high_pt; }
                        }
                    } else { // Negative
                        if (pCat == LOW_P) {
                            if (b == 0) { hRap_EW = hRapidity_xi0_kaon_minus_E; hRap_C = hRapidity_xi0_kaon_minus_C; }
                            if (b == 1) { hRap_EW = hRapidity_xi1_kaon_minus_E; hRap_C = hRapidity_xi1_kaon_minus_C; }
                            if (b == 2) { hRap_EW = hRapidity_xi2_kaon_minus_E; hRap_C = hRapidity_xi2_kaon_minus_C; }
                            if (b == 3) { hRap_EW = hRapidity_xi3_kaon_minus_E; hRap_C = hRapidity_xi3_kaon_minus_C; }
                            if (b == 4) { hRap_EW = hRapidity_xi4_kaon_minus_E; hRap_C = hRapidity_xi4_kaon_minus_C; }
                        } else if (pCat == MID_P) {
                            if (b == 0) { hRap_EW = hRapidity_xi0_kaon_minus_E_mid_pt; hRap_C = hRapidity_xi0_kaon_minus_C_mid_pt; }
                            if (b == 1) { hRap_EW = hRapidity_xi1_kaon_minus_E_mid_pt; hRap_C = hRapidity_xi1_kaon_minus_C_mid_pt; }
                            if (b == 2) { hRap_EW = hRapidity_xi2_kaon_minus_E_mid_pt; hRap_C = hRapidity_xi2_kaon_minus_C_mid_pt; }
                            if (b == 3) { hRap_EW = hRapidity_xi3_kaon_minus_E_mid_pt; hRap_C = hRapidity_xi3_kaon_minus_C_mid_pt; }
                            if (b == 4) { hRap_EW = hRapidity_xi4_kaon_minus_E_mid_pt; hRap_C = hRapidity_xi4_kaon_minus_C_mid_pt; }
                        } else {
                            if (b == 0) { hRap_EW = hRapidity_xi0_kaon_minus_E_high_pt; hRap_C = hRapidity_xi0_kaon_minus_C_high_pt; }
                            if (b == 1) { hRap_EW = hRapidity_xi1_kaon_minus_E_high_pt; hRap_C = hRapidity_xi1_kaon_minus_C_high_pt; }
                            if (b == 2) { hRap_EW = hRapidity_xi2_kaon_minus_E_high_pt; hRap_C = hRapidity_xi2_kaon_minus_C_high_pt; }
                            if (b == 3) { hRap_EW = hRapidity_xi3_kaon_minus_E_high_pt; hRap_C = hRapidity_xi3_kaon_minus_C_high_pt; }
                            if (b == 4) { hRap_EW = hRapidity_xi4_kaon_minus_E_high_pt; hRap_C = hRapidity_xi4_kaon_minus_C_high_pt; }
                        }
                    }
                }
                
                if (isWest) {
                    if (isPositive) {
                        if (pCat == LOW_P) {
                            if (b == 0) { hRap_EW = hRapidity_xi0_kaon_plus_W; hRap_C = hRapidity_xi0_kaon_plus_C; }
                            if (b == 1) { hRap_EW = hRapidity_xi1_kaon_plus_W; hRap_C = hRapidity_xi1_kaon_plus_C; }
                            if (b == 2) { hRap_EW = hRapidity_xi2_kaon_plus_W; hRap_C = hRapidity_xi2_kaon_plus_C; }
                            if (b == 3) { hRap_EW = hRapidity_xi3_kaon_plus_W; hRap_C = hRapidity_xi3_kaon_plus_C; }
                            if (b == 4) { hRap_EW = hRapidity_xi4_kaon_plus_W; hRap_C = hRapidity_xi4_kaon_plus_C; }
                        } else if (pCat == MID_P) {
                            if (b == 0) { hRap_EW = hRapidity_xi0_kaon_plus_W_mid_pt; hRap_C = hRapidity_xi0_kaon_plus_C_mid_pt; }
                            if (b == 1) { hRap_EW = hRapidity_xi1_kaon_plus_W_mid_pt; hRap_C = hRapidity_xi1_kaon_plus_C_mid_pt; }
                            if (b == 2) { hRap_EW = hRapidity_xi2_kaon_plus_W_mid_pt; hRap_C = hRapidity_xi2_kaon_plus_C_mid_pt; }
                            if (b == 3) { hRap_EW = hRapidity_xi3_kaon_plus_W_mid_pt; hRap_C = hRapidity_xi3_kaon_plus_C_mid_pt; }
                            if (b == 4) { hRap_EW = hRapidity_xi4_kaon_plus_W_mid_pt; hRap_C = hRapidity_xi4_kaon_plus_C_mid_pt; }
                        } else {
                            if (b == 0) { hRap_EW = hRapidity_xi0_kaon_plus_W_high_pt; hRap_C = hRapidity_xi0_kaon_plus_C_high_pt; }
                            if (b == 1) { hRap_EW = hRapidity_xi1_kaon_plus_W_high_pt; hRap_C = hRapidity_xi1_kaon_plus_C_high_pt; }
                            if (b == 2) { hRap_EW = hRapidity_xi2_kaon_plus_W_high_pt; hRap_C = hRapidity_xi2_kaon_plus_C_high_pt; }
                            if (b == 3) { hRap_EW = hRapidity_xi3_kaon_plus_W_high_pt; hRap_C = hRapidity_xi3_kaon_plus_C_high_pt; }
                            if (b == 4) { hRap_EW = hRapidity_xi4_kaon_plus_W_high_pt; hRap_C = hRapidity_xi4_kaon_plus_C_high_pt; }
                        }
                    } else { // Negative
                        if (pCat == LOW_P) {
                            if (b == 0) { hRap_EW = hRapidity_xi0_kaon_minus_W; hRap_C = hRapidity_xi0_kaon_minus_C; }
                            if (b == 1) { hRap_EW = hRapidity_xi1_kaon_minus_W; hRap_C = hRapidity_xi1_kaon_minus_C; }
                            if (b == 2) { hRap_EW = hRapidity_xi2_kaon_minus_W; hRap_C = hRapidity_xi2_kaon_minus_C; }
                            if (b == 3) { hRap_EW = hRapidity_xi3_kaon_minus_W; hRap_C = hRapidity_xi3_kaon_minus_C; }
                            if (b == 4) { hRap_EW = hRapidity_xi4_kaon_minus_W; hRap_C = hRapidity_xi4_kaon_minus_C; }
                        } else if (pCat == MID_P) {
                            if (b == 0) { hRap_EW = hRapidity_xi0_kaon_minus_W_mid_pt; hRap_C = hRapidity_xi0_kaon_minus_C_mid_pt; }
                            if (b == 1) { hRap_EW = hRapidity_xi1_kaon_minus_W_mid_pt; hRap_C = hRapidity_xi1_kaon_minus_C_mid_pt; }
                            if (b == 2) { hRap_EW = hRapidity_xi2_kaon_minus_W_mid_pt; hRap_C = hRapidity_xi2_kaon_minus_C_mid_pt; }
                            if (b == 3) { hRap_EW = hRapidity_xi3_kaon_minus_W_mid_pt; hRap_C = hRapidity_xi3_kaon_minus_C_mid_pt; }
                            if (b == 4) { hRap_EW = hRapidity_xi4_kaon_minus_W_mid_pt; hRap_C = hRapidity_xi4_kaon_minus_C_mid_pt; }
                        } else {
                            if (b == 0) { hRap_EW = hRapidity_xi0_kaon_minus_W_high_pt; hRap_C = hRapidity_xi0_kaon_minus_C_high_pt; }
                            if (b == 1) { hRap_EW = hRapidity_xi1_kaon_minus_W_high_pt; hRap_C = hRapidity_xi1_kaon_minus_C_high_pt; }
                            if (b == 2) { hRap_EW = hRapidity_xi2_kaon_minus_W_high_pt; hRap_C = hRapidity_xi2_kaon_minus_C_high_pt; }
                            if (b == 3) { hRap_EW = hRapidity_xi3_kaon_minus_W_high_pt; hRap_C = hRapidity_xi3_kaon_minus_C_high_pt; }
                            if (b == 4) { hRap_EW = hRapidity_xi4_kaon_minus_W_high_pt; hRap_C = hRapidity_xi4_kaon_minus_C_high_pt; }
                        }
                    }
                }
                if (hRap_EW) hRap_EW->Fill(rapidity);
                if (hRap_C)  hRap_C->Fill(isEast ? -rapidity : rapidity);
            }
        }

        for (int i_p : pion_Idx) {
            TVector3 mom_p;
            upcEvt->getTrack(i_p)->getMomentum(mom_p);
            double p_p = mom_p.Mag();

            enum MomentumCat { LOW_P, MID_P, HIGH_P };
            MomentumCat pCat = LOW_P;
            if (p_p >= 0.5 && p_p <= 0.9)      pCat = MID_P;
            else if (p_p > 0.9)                pCat = HIGH_P;

            bool isPositive = (upcEvt->getTrack(i_p)->getCharge() > 0);

            if (znaleziono_anchor && global_ref_idx_dedx != i_p) {
                double t0_pion_hip = upcEvt->getTrack(i_p)->getT0(massPion);
                double delta_t_global = t0_pion_hip - wyznaczony_t0;

                int k_bin = getMomentumBin(p_p);
                if (k_bin != -1) {
                    int i_part = -1; 
                    // Indeksy w particles: 0:Kp, 1:Kn, 2:Pip, 3:Pin, 4:Pp, 5:Pn
                    if      (particles[2] == "Pip" && isPositive) i_part = 2;
                    else if (particles[3] == "Pin" && !isPositive) i_part = 3;

                    int j_part = -1; // partnerzy: 0:pion, 1:kaon, 2:proton
                    double chi2_pi = pow(upcEvt->getTrack(global_ref_idx_dedx)->getNSigmasTPCPion(), 2);
                    double chi2_K  = pow(upcEvt->getTrack(global_ref_idx_dedx)->getNSigmasTPCKaon(), 2);
                    double chi2_p  = pow(upcEvt->getTrack(global_ref_idx_dedx)->getNSigmasTPCProton(), 2);

                    if      (chi2_pi <= chi2_K && chi2_pi <= chi2_p) j_part = 0;
                    else if (chi2_K <= chi2_pi && chi2_K <= chi2_p) j_part = 1;
                    else    j_part = 2;

                    if (i_part != -1 && j_part != -1) {
                        hDeltaT[i_part][j_part][k_bin]->Fill(delta_t_global);

                        double min_delta=delta_t_global;
                        double tof_limit = TOFPIONS_pions;
                        if (j_part == 1) tof_limit = TOFPIONS_kaons;
                        else if (j_part == 2) tof_limit = TOFPIONS_protons;

                        if (fabs(min_delta) / (tof_limit / 3.0) < 0.25 && (fabs(upcEvt->getTrack(i_p)->getNSigmasTPCPion())<0.25)) {
                                double dedx_meas=upcEvt->getTrack(i_p)->getDEdxSignal()*1e6;
                                TVector3 mom;
                                upcEvt->getTrack(i_p)->getMomentum(mom);
                                double p=mom.Mag();
                                for (int i = 0; i < nDx; ++i) {
                                    double dx_val = dx_targets[i];
                                    double dedx_theo = GetTheoreticaldEdx(p, massPion, dx_val); 
                                    double ratio = dedx_meas/dedx_theo;
                                    hRatioTOFVsP_pion[i]->Fill(p, ratio);
                                }
                            }
                    }
                }

                double chi2_pi = pow(upcEvt->getTrack(global_ref_idx_dedx)->getNSigmasTPCPion(), 2);
                double chi2_K  = pow(upcEvt->getTrack(global_ref_idx_dedx)->getNSigmasTPCKaon(), 2);
                double chi2_p  = pow(upcEvt->getTrack(global_ref_idx_dedx)->getNSigmasTPCProton(), 2);

                TH1* hDeltaT_global = nullptr;

                if (chi2_pi <= chi2_K && chi2_pi <= chi2_p) {
                    if (isPositive) {
                        if (pCat == LOW_P)       hDeltaT_global = hDeltaT_global_anchor_Pippi;
                        else if (pCat == MID_P)  hDeltaT_global = hDeltaT_global_anchor_Pippi_mid_pt;
                        else                     hDeltaT_global = hDeltaT_global_anchor_Pippi_high_pt;
                    } else {
                        if (pCat == LOW_P)       hDeltaT_global = hDeltaT_global_anchor_Pinpi;
                        else if (pCat == MID_P)  hDeltaT_global = hDeltaT_global_anchor_Pinpi_mid_pt;
                        else                     hDeltaT_global = hDeltaT_global_anchor_Pinpi_high_pt;
                    }
                } 
                else if (chi2_K <= chi2_pi && chi2_K <= chi2_p) {
                    if (isPositive) {
                        if (pCat == LOW_P)       hDeltaT_global = hDeltaT_global_anchor_PipK;
                        else if (pCat == MID_P)  hDeltaT_global = hDeltaT_global_anchor_PipK_mid_pt;
                        else                     hDeltaT_global = hDeltaT_global_anchor_PipK_high_pt;
                    } else {
                        if (pCat == LOW_P)       hDeltaT_global = hDeltaT_global_anchor_PinK;
                        else if (pCat == MID_P)  hDeltaT_global = hDeltaT_global_anchor_PinK_mid_pt;
                        else                     hDeltaT_global = hDeltaT_global_anchor_PinK_high_pt;
                    }
                } 
                else {
                    if (isPositive) {
                        if (pCat == LOW_P)       hDeltaT_global = hDeltaT_global_anchor_Pipp;
                        else if (pCat == MID_P)  hDeltaT_global = hDeltaT_global_anchor_Pipp_mid_pt;
                        else                     hDeltaT_global = hDeltaT_global_anchor_Pipp_high_pt;
                    } else {
                        if (pCat == LOW_P)       hDeltaT_global = hDeltaT_global_anchor_Pinp;
                        else if (pCat == MID_P)  hDeltaT_global = hDeltaT_global_anchor_Pinp_mid_pt;
                        else                     hDeltaT_global = hDeltaT_global_anchor_Pinp_high_pt;
                    }
                }

                if (hDeltaT_global) hDeltaT_global->Fill(delta_t_global);
            }
            int ref_idx = -1;
            double min_p = 9999.0;

            for (int i_trk : goodTrackIdx) {
                if (i_trk == i_p) continue;
                TVector3 mom;
                upcEvt->getTrack(i_trk)->getMomentum(mom);

                if (mom.Mag() < min_p) {
                    min_p = mom.Mag();
                    ref_idx = i_trk;
                }
            }

            if (ref_idx >= 0) {
                double t0_pion_hip = upcEvt->getTrack(i_p)->getT0(massPion);
                double t0_pion     = upcEvt->getTrack(ref_idx)->getT0(massPion);
                double t0_kaon     = upcEvt->getTrack(ref_idx)->getT0(massKaon);
                double t0_proton   = upcEvt->getTrack(ref_idx)->getT0(massProton);

                double delta_t1 = t0_pion_hip - t0_pion;
                double delta_t2 = t0_pion_hip - t0_kaon;
                double delta_t5 = t0_pion_hip - t0_proton;

                TH2* hT1vsT2_charge = nullptr;
                TH2* hT1vsT2_all    = nullptr;
                TH2* hT1vsT5_charge = nullptr;

                if (isPositive) {
                    if (pCat == LOW_P) {
                        hT1vsT2_charge = hDeltaT1vsT2_pions_p;
                        hT1vsT2_all    = hDeltaT1vsT2_pions;
                        hT1vsT5_charge = hDeltaT1vsT5_pions_p;
                    } else if (pCat == MID_P) {
                        hT1vsT2_charge = hDeltaT1vsT2_pions_p_mid_p;
                        hT1vsT2_all    = hDeltaT1vsT2_pions_mid_p;
                        hT1vsT5_charge = hDeltaT1vsT5_pions_p_mid_p;
                    } else {
                        hT1vsT2_charge = hDeltaT1vsT2_pions_p_high_p; 
                        hT1vsT2_all    = hDeltaT1vsT2_pions_high_p;
                        hT1vsT5_charge = hDeltaT1vsT5_pions_p_high_p;
                    }
                } else { // Negative
                    if (pCat == LOW_P) {
                        hT1vsT2_charge = hDeltaT1vsT2_pions_n;
                        hT1vsT2_all    = hDeltaT1vsT2_pions;
                        hT1vsT5_charge = hDeltaT1vsT5_pions_n;
                    } else if (pCat == MID_P) {
                        hT1vsT2_charge = hDeltaT1vsT2_pions_n_mid_p;
                        hT1vsT2_all    = hDeltaT1vsT2_pions_mid_p;
                        hT1vsT5_charge = hDeltaT1vsT5_pions_n_mid_p;
                    } else {
                        hT1vsT2_charge = hDeltaT1vsT2_pions_n_high_p;
                        hT1vsT2_all    = hDeltaT1vsT2_pions_high_p;
                        hT1vsT5_charge = hDeltaT1vsT5_pions_n_high_p;
                    }
                }

                if (hT1vsT2_charge) hT1vsT2_charge->Fill(delta_t1, delta_t2);
                if (hT1vsT2_all)    hT1vsT2_all->Fill(delta_t1, delta_t2);
                if (hT1vsT5_charge) hT1vsT5_charge->Fill(delta_t1, delta_t5);

                if (!(fabs(delta_t1) < TOFPIONS_pions || fabs(delta_t2) < TOFPIONS_kaons || fabs(delta_t5) < TOFPIONS_protons)) continue;

                // double abs1 = std::fabs(delta_t1);
                // double abs2 = std::fabs(delta_t2);
                // double abs5 = std::fabs(delta_t5);
                
                // double min_delta = abs1;
                // double tof_limit = TOFPIONS_pions;
                
                // if (abs2 < min_delta) { min_delta = abs2; tof_limit = TOFPIONS_kaons; }
                // if (abs5 < min_delta) { min_delta = abs5; tof_limit = TOFPIONS_protons; }
                
                // if (min_delta / (tof_limit / 3.0) < 0.25 && (fabs(upcEvt->getTrack(i_p)->getNSigmasTPCPion())<0.25)) {
                //     double dedx_meas=upcEvt->getTrack(i_p)->getDEdxSignal()*1e6;
                //     TVector3 mom;
                //     upcEvt->getTrack(i_p)->getMomentum(mom);
                //     double p=mom.Mag();
                //     for (int i = 0; i < nDx; ++i) {
                //         double dx_val = dx_targets[i];
                //         double dedx_theo = GetTheoreticaldEdx(p, massPion, dx_val); 
                //         double ratio = dedx_meas/dedx_theo;
                //         hRatioTOFVsP_pion[i]->Fill(p, ratio);
                //     }
                // }

                double pt       = upcEvt->getTrack(i_p)->getPt();
                double eta      = upcEvt->getTrack(i_p)->getEta();
                double rapidity = getRapidity(pt, eta, upcEvt->getTrack(i_p)->getPhi(), massPion); 

                TH1* hRap_EW = nullptr;
                TH1* hRap_C  = nullptr;

                if (isEast) {
                    if (isPositive) {
                        if (pCat == LOW_P) {
                            if (b == 0) { hRap_EW = hRapidity_xi0_pion_plus_E; hRap_C = hRapidity_xi0_pion_plus_C; }
                            if (b == 1) { hRap_EW = hRapidity_xi1_pion_plus_E; hRap_C = hRapidity_xi1_pion_plus_C; }
                            if (b == 2) { hRap_EW = hRapidity_xi2_pion_plus_E; hRap_C = hRapidity_xi2_pion_plus_C; }
                            if (b == 3) { hRap_EW = hRapidity_xi3_pion_plus_E; hRap_C = hRapidity_xi3_pion_plus_C; }
                            if (b == 4) { hRap_EW = hRapidity_xi4_pion_plus_E; hRap_C = hRapidity_xi4_pion_plus_C; }
                        } else if (pCat == MID_P) {
                            if (b == 0) { hRap_EW = hRapidity_xi0_pion_plus_E_mid_pt; hRap_C = hRapidity_xi0_pion_plus_C_mid_pt; }
                            if (b == 1) { hRap_EW = hRapidity_xi1_pion_plus_E_mid_pt; hRap_C = hRapidity_xi1_pion_plus_C_mid_pt; }
                            if (b == 2) { hRap_EW = hRapidity_xi2_pion_plus_E_mid_pt; hRap_C = hRapidity_xi2_pion_plus_C_mid_pt; }
                            if (b == 3) { hRap_EW = hRapidity_xi3_pion_plus_E_mid_pt; hRap_C = hRapidity_xi3_pion_plus_C_mid_pt; }
                            if (b == 4) { hRap_EW = hRapidity_xi4_pion_plus_E_mid_pt; hRap_C = hRapidity_xi4_pion_plus_C_mid_pt; }
                        } else { // HIGH_P
                            if (b == 0) { hRap_EW = hRapidity_xi0_pion_plus_E_high_pt; hRap_C = hRapidity_xi0_pion_plus_C_high_pt; }
                            if (b == 1) { hRap_EW = hRapidity_xi1_pion_plus_E_high_pt; hRap_C = hRapidity_xi1_pion_plus_C_high_pt; }
                            if (b == 2) { hRap_EW = hRapidity_xi2_pion_plus_E_high_pt; hRap_C = hRapidity_xi2_pion_plus_C_high_pt; }
                            if (b == 3) { hRap_EW = hRapidity_xi3_pion_plus_E_high_pt; hRap_C = hRapidity_xi3_pion_plus_C_high_pt; }
                            if (b == 4) { hRap_EW = hRapidity_xi4_pion_plus_E_high_pt; hRap_C = hRapidity_xi4_pion_plus_C_high_pt; }
                        }
                    } else { // Negative
                        if (pCat == LOW_P) {
                            if (b == 0) { hRap_EW = hRapidity_xi0_pion_minus_E; hRap_C = hRapidity_xi0_pion_minus_C; }
                            if (b == 1) { hRap_EW = hRapidity_xi1_pion_minus_E; hRap_C = hRapidity_xi1_pion_minus_C; }
                            if (b == 2) { hRap_EW = hRapidity_xi2_pion_minus_E; hRap_C = hRapidity_xi2_pion_minus_C; }
                            if (b == 3) { hRap_EW = hRapidity_xi3_pion_minus_E; hRap_C = hRapidity_xi3_pion_minus_C; }
                            if (b == 4) { hRap_EW = hRapidity_xi4_pion_minus_E; hRap_C = hRapidity_xi4_pion_minus_C; }
                        } else if (pCat == MID_P) {
                            if (b == 0) { hRap_EW = hRapidity_xi0_pion_minus_E_mid_pt; hRap_C = hRapidity_xi0_pion_minus_C_mid_pt; }
                            if (b == 1) { hRap_EW = hRapidity_xi1_pion_minus_E_mid_pt; hRap_C = hRapidity_xi1_pion_minus_C_mid_pt; }
                            if (b == 2) { hRap_EW = hRapidity_xi2_pion_minus_E_mid_pt; hRap_C = hRapidity_xi2_pion_minus_C_mid_pt; }
                            if (b == 3) { hRap_EW = hRapidity_xi3_pion_minus_E_mid_pt; hRap_C = hRapidity_xi3_pion_minus_C_mid_pt; }
                            if (b == 4) { hRap_EW = hRapidity_xi4_pion_minus_E_mid_pt; hRap_C = hRapidity_xi4_pion_minus_C_mid_pt; }
                        } else { // HIGH_P
                            if (b == 0) { hRap_EW = hRapidity_xi0_pion_minus_E_high_pt; hRap_C = hRapidity_xi0_pion_minus_C_high_pt; }
                            if (b == 1) { hRap_EW = hRapidity_xi1_pion_minus_E_high_pt; hRap_C = hRapidity_xi1_pion_minus_C_high_pt; }
                            if (b == 2) { hRap_EW = hRapidity_xi2_pion_minus_E_high_pt; hRap_C = hRapidity_xi2_pion_minus_C_high_pt; }
                            if (b == 3) { hRap_EW = hRapidity_xi3_pion_minus_E_high_pt; hRap_C = hRapidity_xi3_pion_minus_C_high_pt; }
                            if (b == 4) { hRap_EW = hRapidity_xi4_pion_minus_E_high_pt; hRap_C = hRapidity_xi4_pion_minus_C_high_pt; }
                        }
                    }
                }
                
                if (isWest) {
                    if (isPositive) {
                        if (pCat == LOW_P) {
                            if (b == 0) { hRap_EW = hRapidity_xi0_pion_plus_W; hRap_C = hRapidity_xi0_pion_plus_C; }
                            if (b == 1) { hRap_EW = hRapidity_xi1_pion_plus_W; hRap_C = hRapidity_xi1_pion_plus_C; }
                            if (b == 2) { hRap_EW = hRapidity_xi2_pion_plus_W; hRap_C = hRapidity_xi2_pion_plus_C; }
                            if (b == 3) { hRap_EW = hRapidity_xi3_pion_plus_W; hRap_C = hRapidity_xi3_pion_plus_C; }
                            if (b == 4) { hRap_EW = hRapidity_xi4_pion_plus_W; hRap_C = hRapidity_xi4_pion_plus_C; }
                        } else if (pCat == MID_P) {
                            if (b == 0) { hRap_EW = hRapidity_xi0_pion_plus_W_mid_pt; hRap_C = hRapidity_xi0_pion_plus_C_mid_pt; }
                            if (b == 1) { hRap_EW = hRapidity_xi1_pion_plus_W_mid_pt; hRap_C = hRapidity_xi1_pion_plus_C_mid_pt; }
                            if (b == 2) { hRap_EW = hRapidity_xi2_pion_plus_W_mid_pt; hRap_C = hRapidity_xi2_pion_plus_C_mid_pt; }
                            if (b == 3) { hRap_EW = hRapidity_xi3_pion_plus_W_mid_pt; hRap_C = hRapidity_xi3_pion_plus_C_mid_pt; }
                            if (b == 4) { hRap_EW = hRapidity_xi4_pion_plus_W_mid_pt; hRap_C = hRapidity_xi4_pion_plus_C_mid_pt; }
                        } else { // HIGH_P
                            if (b == 0) { hRap_EW = hRapidity_xi0_pion_plus_W_high_pt; hRap_C = hRapidity_xi0_pion_plus_C_high_pt; }
                            if (b == 1) { hRap_EW = hRapidity_xi1_pion_plus_W_high_pt; hRap_C = hRapidity_xi1_pion_plus_C_high_pt; }
                            if (b == 2) { hRap_EW = hRapidity_xi2_pion_plus_W_high_pt; hRap_C = hRapidity_xi2_pion_plus_C_high_pt; }
                            if (b == 3) { hRap_EW = hRapidity_xi3_pion_plus_W_high_pt; hRap_C = hRapidity_xi3_pion_plus_C_high_pt; }
                            if (b == 4) { hRap_EW = hRapidity_xi4_pion_plus_W_high_pt; hRap_C = hRapidity_xi4_pion_plus_C_high_pt; }
                        }
                    } else { // Negative
                        if (pCat == LOW_P) {
                            if (b == 0) { hRap_EW = hRapidity_xi0_pion_minus_W; hRap_C = hRapidity_xi0_pion_minus_C; }
                            if (b == 1) { hRap_EW = hRapidity_xi1_pion_minus_W; hRap_C = hRapidity_xi1_pion_minus_C; }
                            if (b == 2) { hRap_EW = hRapidity_xi2_pion_minus_W; hRap_C = hRapidity_xi2_pion_minus_C; }
                            if (b == 3) { hRap_EW = hRapidity_xi3_pion_minus_W; hRap_C = hRapidity_xi3_pion_minus_C; }
                            if (b == 4) { hRap_EW = hRapidity_xi4_pion_minus_W; hRap_C = hRapidity_xi4_pion_minus_C; }
                        } else if (pCat == MID_P) {
                            if (b == 0) { hRap_EW = hRapidity_xi0_pion_minus_W_mid_pt; hRap_C = hRapidity_xi0_pion_minus_C_mid_pt; }
                            if (b == 1) { hRap_EW = hRapidity_xi1_pion_minus_W_mid_pt; hRap_C = hRapidity_xi1_pion_minus_C_mid_pt; }
                            if (b == 2) { hRap_EW = hRapidity_xi2_pion_minus_W_mid_pt; hRap_C = hRapidity_xi2_pion_minus_C_mid_pt; }
                            if (b == 3) { hRap_EW = hRapidity_xi3_pion_minus_W_mid_pt; hRap_C = hRapidity_xi3_pion_minus_C_mid_pt; }
                            if (b == 4) { hRap_EW = hRapidity_xi4_pion_minus_W_mid_pt; hRap_C = hRapidity_xi4_pion_minus_C_mid_pt; }
                        } else { // HIGH_P
                            if (b == 0) { hRap_EW = hRapidity_xi0_pion_minus_W_high_pt; hRap_C = hRapidity_xi0_pion_minus_C_high_pt; }
                            if (b == 1) { hRap_EW = hRapidity_xi1_pion_minus_W_high_pt; hRap_C = hRapidity_xi1_pion_minus_C_high_pt; }
                            if (b == 2) { hRap_EW = hRapidity_xi2_pion_minus_W_high_pt; hRap_C = hRapidity_xi2_pion_minus_C_high_pt; }
                            if (b == 3) { hRap_EW = hRapidity_xi3_pion_minus_W_high_pt; hRap_C = hRapidity_xi3_pion_minus_C_high_pt; }
                            if (b == 4) { hRap_EW = hRapidity_xi4_pion_minus_W_high_pt; hRap_C = hRapidity_xi4_pion_minus_C_high_pt; }
                        }
                    }
                }

                if (hRap_EW) hRap_EW->Fill(rapidity);
                if (hRap_C)  hRap_C->Fill(isEast ? -rapidity : rapidity);
            }
        }

        for (int i_p : proton_Idx) {
            TVector3 mom_p;
            upcEvt->getTrack(i_p)->getMomentum(mom_p);
            double p_p = mom_p.Mag();

            enum MomentumCat { LOW_P, MID_P, HIGH_P };
            MomentumCat pCat = LOW_P;
            if (p_p >= 0.5 && p_p <= 0.9)      pCat = MID_P;
            else if (p_p > 0.9)                pCat = HIGH_P;

            bool isPositive = (upcEvt->getTrack(i_p)->getCharge() > 0);

            if (znaleziono_anchor && global_ref_idx_dedx != i_p) {
                double t0_proton_hip = upcEvt->getTrack(i_p)->getT0(massProton);
                double delta_t_global = t0_proton_hip - wyznaczony_t0;

                int k_bin = getMomentumBin(p_p);
                if (k_bin != -1) {
                    int i_part = -1;
                    if      (particles[4] == "Pp"  && isPositive) i_part = 4;
                    else if (particles[5] == "Pn"  && !isPositive) i_part = 5;

                    int j_part = -1;
                    double chi2_pi = pow(upcEvt->getTrack(global_ref_idx_dedx)->getNSigmasTPCPion(), 2);
                    double chi2_K  = pow(upcEvt->getTrack(global_ref_idx_dedx)->getNSigmasTPCKaon(), 2);
                    double chi2_p  = pow(upcEvt->getTrack(global_ref_idx_dedx)->getNSigmasTPCProton(), 2);

                    if (chi2_pi <= chi2_K && chi2_pi <= chi2_p) j_part = 0;
                    else if (chi2_K <= chi2_pi && chi2_K <= chi2_p) j_part = 1;
                    else j_part = 2;

                    if (i_part != -1 && j_part != -1) {
                        hDeltaT[i_part][j_part][k_bin]->Fill(delta_t_global);

                        double min_delta=delta_t_global;
                        double tof_limit = TOF_pions;
                        if (j_part == 1) tof_limit = TOF_kaons;
                        else if (j_part == 2) tof_limit = TOF_protons;

                        if (fabs(min_delta) / (tof_limit / 3.0) < 0.25 && (fabs(upcEvt->getTrack(i_p)->getNSigmasTPCProton())<0.25)) {
                            double dedx_meas=upcEvt->getTrack(i_p)->getDEdxSignal()*1e6;
                            TVector3 mom;
                            upcEvt->getTrack(i_p)->getMomentum(mom);
                            double p=mom.Mag();
                            for (int i = 0; i < nDx; ++i) {
                                double dx_val = dx_targets[i];
                                double dedx_theo = GetTheoreticaldEdx(p, massProton, dx_val); 
                                double ratio = dedx_meas/dedx_theo;
                                hRatioTOFVsP_proton[i]->Fill(p, ratio);
                            }
                        }
                    }
                }

                double chi2_pi = pow(upcEvt->getTrack(global_ref_idx_dedx)->getNSigmasTPCPion(), 2);
                double chi2_K  = pow(upcEvt->getTrack(global_ref_idx_dedx)->getNSigmasTPCKaon(), 2);
                double chi2_p  = pow(upcEvt->getTrack(global_ref_idx_dedx)->getNSigmasTPCProton(), 2);

                TH1* hDeltaT_global = nullptr;

                if (chi2_pi <= chi2_K && chi2_pi <= chi2_p) {
                    if (isPositive) {
                        if (pCat == LOW_P)       hDeltaT_global = hDeltaT_global_anchor_Pppi;
                        else if (pCat == MID_P)  hDeltaT_global = hDeltaT_global_anchor_Pppi_mid_pt;
                        else                     hDeltaT_global = hDeltaT_global_anchor_Pppi_high_pt;
                    } else {
                        if (pCat == LOW_P)       hDeltaT_global = hDeltaT_global_anchor_Pnpi;
                        else if (pCat == MID_P)  hDeltaT_global = hDeltaT_global_anchor_Pnpi_mid_pt;
                        else                     hDeltaT_global = hDeltaT_global_anchor_Pnpi_high_pt;
                    }
                } 
                else if (chi2_K <= chi2_pi && chi2_K <= chi2_p) {
                    if (isPositive) {
                        if (pCat == LOW_P)       hDeltaT_global = hDeltaT_global_anchor_PpK;
                        else if (pCat == MID_P)  hDeltaT_global = hDeltaT_global_anchor_PpK_mid_pt;
                        else                     hDeltaT_global = hDeltaT_global_anchor_PpK_high_pt;
                    } else {
                        if (pCat == LOW_P)       hDeltaT_global = hDeltaT_global_anchor_PnK;
                        else if (pCat == MID_P)  hDeltaT_global = hDeltaT_global_anchor_PnK_mid_pt;
                        else                     hDeltaT_global = hDeltaT_global_anchor_PnK_high_pt;
                    }
                } 
                else {
                    if (isPositive) {
                        if (pCat == LOW_P)       hDeltaT_global = hDeltaT_global_anchor_Ppp;
                        else if (pCat == MID_P)  hDeltaT_global = hDeltaT_global_anchor_Ppp_mid_pt;
                        else                     hDeltaT_global = hDeltaT_global_anchor_Ppp_high_pt;
                    } else {
                        if (pCat == LOW_P)       hDeltaT_global = hDeltaT_global_anchor_Pnp;
                        else if (pCat == MID_P)  hDeltaT_global = hDeltaT_global_anchor_Pnp_mid_pt;
                        else                     hDeltaT_global = hDeltaT_global_anchor_Pnp_high_pt;
                    }
                }

                if (hDeltaT_global) hDeltaT_global->Fill(delta_t_global);
            }

            int ref_idx = -1;
            double min_p = 9999.0;

            for (int i_trk : goodTrackIdx) {
                if (i_trk == i_p) continue;
                TVector3 mom;
                upcEvt->getTrack(i_trk)->getMomentum(mom);

                if (mom.Mag() < min_p) {
                    min_p = mom.Mag();
                    ref_idx = i_trk;
                }
            }

            if (ref_idx >= 0) {
                double t0_nucleon      = upcEvt->getTrack(i_p)->getT0(massProton);
                double t0_pion         = upcEvt->getTrack(ref_idx)->getT0(massPion);
                double t0_kaon         = upcEvt->getTrack(ref_idx)->getT0(massKaon);
                double t0_test_nucleon = upcEvt->getTrack(ref_idx)->getT0(massProton);

                double delta_t1 = t0_nucleon - t0_pion;
                double delta_t2 = t0_nucleon - t0_kaon;
                double delta_t5 = t0_nucleon - t0_test_nucleon;

                TH2* hT1vsT2_charge      = nullptr;
                TH2* hT1vsT2_charge_west = nullptr;
                TH2* hT1vsT2_charge_east = nullptr;
                TH2* hT1vsT5_charge      = nullptr;
                TH2* hT1vsT2_ppbar_all   = nullptr;

                if (isPositive) {
                    if (pCat == LOW_P) {
                        hT1vsT2_charge      = hDeltaT1vsT2_protons;
                        hT1vsT2_charge_west = hDeltaT1vsT2_protons_west;
                        hT1vsT2_charge_east = hDeltaT1vsT2_protons_east;
                        hT1vsT5_charge      = hDeltaT1vsT5_protons;
                        hT1vsT2_ppbar_all   = hDeltaT1vsT2_ppbar;
                    } else if (pCat == MID_P) {
                        hT1vsT2_charge      = hDeltaT1vsT2_protons_mid_p;
                        hT1vsT2_charge_west = hDeltaT1vsT2_protons_west_mid_p;
                        hT1vsT2_charge_east = hDeltaT1vsT2_protons_east_mid_p;
                        hT1vsT5_charge      = hDeltaT1vsT5_protons_mid_p;
                        hT1vsT2_ppbar_all   = hDeltaT1vsT2_ppbar_mid_p;
                    } else { // HIGH_P
                        hT1vsT2_charge      = hDeltaT1vsT2_protons_high_p;
                        hT1vsT2_charge_west = hDeltaT1vsT2_protons_west_high_p;
                        hT1vsT2_charge_east = hDeltaT1vsT2_protons_east_high_p;
                        hT1vsT5_charge      = hDeltaT1vsT5_protons_high_p;
                        hT1vsT2_ppbar_all   = hDeltaT1vsT2_ppbar_high_p;
                    }
                } else { // Negative
                    if (pCat == LOW_P) {
                        hT1vsT2_charge      = hDeltaT1vsT2_antiprotons;
                        hT1vsT2_charge_west = hDeltaT1vsT2_antiprotons_west;
                        hT1vsT2_charge_east = hDeltaT1vsT2_antiprotons_east;
                        hT1vsT5_charge      = hDeltaT1vsT5_antiprotons;
                        hT1vsT2_ppbar_all   = hDeltaT1vsT2_ppbar;
                    } else if (pCat == MID_P) {
                        hT1vsT2_charge      = hDeltaT1vsT2_antiprotons_mid_p;
                        hT1vsT2_charge_west = hDeltaT1vsT2_antiprotons_west_mid_p;
                        hT1vsT2_charge_east = hDeltaT1vsT2_antiprotons_east_mid_p;
                        hT1vsT5_charge      = hDeltaT1vsT5_antiprotons_mid_p;
                        hT1vsT2_ppbar_all   = hDeltaT1vsT2_ppbar_mid_p;
                    } else { // HIGH_P
                        hT1vsT2_charge      = hDeltaT1vsT2_antiprotons_high_p;
                        hT1vsT2_charge_west = hDeltaT1vsT2_antiprotons_west_high_p;
                        hT1vsT2_charge_east = hDeltaT1vsT2_antiprotons_east_high_p;
                        hT1vsT5_charge      = hDeltaT1vsT5_antiprotons_high_p;
                        hT1vsT2_ppbar_all   = hDeltaT1vsT2_ppbar_high_p;
                    }
                }

                if (hT1vsT2_charge)    hT1vsT2_charge->Fill(delta_t1, delta_t2);
                if (hT1vsT5_charge)    hT1vsT5_charge->Fill(delta_t1, delta_t5);
                if (hT1vsT2_ppbar_all) hT1vsT2_ppbar_all->Fill(delta_t1, delta_t2);

                if (isWest) {
                    if (hT1vsT2_charge_west) hT1vsT2_charge_west->Fill(delta_t1, delta_t2);
                } else {
                    if (hT1vsT2_charge_east) hT1vsT2_charge_east->Fill(delta_t1, delta_t2);
                }

                TH1* hDT_pp   = nullptr;
                TH1* hDT_pipi = nullptr;
                TH1* hDT_Kpi  = nullptr;

                if (pCat == LOW_P) {
                    hDT_pp   = hDeltaT_protonproton;
                    hDT_pipi = hDeltaT_pionpion;
                    hDT_Kpi  = hDeltaT_kaonpion;
                } else if (pCat == MID_P) {
                    hDT_pp   = hDeltaT_protonproton_mid_p;
                    hDT_pipi = hDeltaT_pionpion_mid_p;
                    hDT_Kpi  = hDeltaT_kaonpion_mid_p;
                } else { // HIGH_P
                    hDT_pp   = hDeltaT_protonproton_high_p;
                    hDT_pipi = hDeltaT_pionpion_high_p;
                    hDT_Kpi  = hDeltaT_kaonpion_high_p;
                }

                if ((fabs(delta_t1) > TOF_pions && fabs(delta_t2) > TOF_kaons)) { 
                    if (hDT_pp) hDT_pp->Fill(delta_t5);
                
                    if ((fabs(delta_t5) > TOF_protons)) {
                        double t0_pion_rejected = upcEvt->getTrack(i_p)->getT0(massPion);
                        double t0_kaon_rejected = upcEvt->getTrack(i_p)->getT0(massKaon);
                        
                        double delta_t3 = t0_pion_rejected - t0_pion;
                        double delta_t4 = t0_kaon_rejected - t0_pion;
                        
                        if (hDT_pipi) hDT_pipi->Fill(delta_t3);
                        if (hDT_Kpi)  hDT_Kpi->Fill(delta_t4);
                    }
                }

                if (!(fabs(delta_t1) < TOF_pions || fabs(delta_t2) < TOF_kaons || fabs(delta_t5) < TOF_protons)) continue;

                // double abs1 = std::fabs(delta_t1);
                // double abs2 = std::fabs(delta_t2);
                // double abs5 = std::fabs(delta_t5);
                
                // double min_delta = abs1;
                // double tof_limit = TOF_pions;
                
                // if (abs2 < min_delta) { min_delta = abs2; tof_limit = TOF_kaons; }
                // if (abs5 < min_delta) { min_delta = abs5; tof_limit = TOF_protons; }
                

                // if (min_delta / (tof_limit / 3.0) < 0.25 && (fabs(upcEvt->getTrack(i_p)->getNSigmasTPCProton())<0.25)) {
                //     double dedx_meas=upcEvt->getTrack(i_p)->getDEdxSignal()*1e6;
                //     TVector3 mom;
                //     upcEvt->getTrack(i_p)->getMomentum(mom);
                //     double p=mom.Mag();
                //     for (int i = 0; i < nDx; ++i) {
                //         double dx_val = dx_targets[i];
                //         double dedx_theo = GetTheoreticaldEdx(p, massProton, dx_val); 
                //         double ratio = dedx_meas/dedx_theo;
                //         hRatioTOFVsP_proton[i]->Fill(p, ratio);
                //     }
                // }
                
                double pt       = upcEvt->getTrack(i_p)->getPt();
                double eta      = upcEvt->getTrack(i_p)->getEta();
                double dca      = upcEvt->getTrack(i_p)->getDcaXY();
                double rapidity = getRapidity(pt, eta, upcEvt->getTrack(i_p)->getPhi(), massProton);
                
                double distance_from_beam = rapidity - y_beam;
                double distance_from_edge = y_edge - rapidity;

                TH3* h3D = nullptr;
                TH1* hPt = nullptr; TH1* hEta = nullptr; TH1* hRap = nullptr;
                TH1* hDistBeam = nullptr; TH1* hDistEdge = nullptr;

                TH3* h3D_E = nullptr;
                TH1* hPt_E = nullptr; TH1* hEta_E = nullptr; TH1* hRap_E = nullptr;
                TH1* hDistBeam_E = nullptr; TH1* hDistEdge_E = nullptr; TH1* hDCA_E = nullptr;

                TH3* h3D_W = nullptr;
                TH1* hPt_W = nullptr; TH1* hEta_W = nullptr; TH1* hRap_W = nullptr;
                TH1* hDistBeam_W = nullptr; TH1* hDistEdge_W = nullptr; TH1* hDCA_W = nullptr;

                TH1* hEta_xi_EW[5] = {nullptr};
                TH1* hRap_xi_EW[5] = {nullptr};
                TH1* hEta_xi_C[5]  = {nullptr};
                TH1* hRap_xi_C[5]  = {nullptr};

                if (isPositive) {
                    if (pCat == LOW_P) {
                        h3D = h3D_protons; hPt = hPtProtonX; hEta = hEtaProtonX; hRap = hRapidityProtonX;
                        hDistBeam = hDistanceBeam_proton; hDistEdge = hDistanceEdge_proton;
                        MultiplicityProtonX++;

                        if (isEast) {
                            h3D_E = h3D_protonsEast; hPt_E = hPtProtonEastX; hEta_E = hEtaProtonEastX; hRap_E = hRapidityProtonEastX;
                            hDistBeam_E = hDistanceBeam_proton_E; hDistEdge_E = hDistanceEdge_proton_E; hDCA_E = hProton_DCA_E;
                            MultiplicityProtonEastX++;

                            hEta_xi_EW[0] = hEta_xi0_proton_E; hRap_xi_EW[0] = hRapidity_xi0_proton_E;
                            hEta_xi_EW[1] = hEta_xi1_proton_E; hRap_xi_EW[1] = hRapidity_xi1_proton_E;
                            hEta_xi_EW[2] = hEta_xi2_proton_E; hRap_xi_EW[2] = hRapidity_xi2_proton_E;
                            hEta_xi_EW[3] = hEta_xi3_proton_E; hRap_xi_EW[3] = hRapidity_xi3_proton_E;
                            hEta_xi_EW[4] = hEta_xi4_proton_E; hRap_xi_EW[4] = hRapidity_xi4_proton_E;

                            hEta_xi_C[0] = hEta_xi0_proton_C; hRap_xi_C[0] = hRapidity_xi0_proton_C;
                            hEta_xi_C[1] = hEta_xi1_proton_C; hRap_xi_C[1] = hRapidity_xi1_proton_C;
                            hEta_xi_C[2] = hEta_xi2_proton_C; hRap_xi_C[2] = hRapidity_xi2_proton_C;
                            hEta_xi_C[3] = hEta_xi3_proton_C; hRap_xi_C[3] = hRapidity_xi3_proton_C;
                            hEta_xi_C[4] = hEta_xi4_proton_C; hRap_xi_C[4] = hRapidity_xi4_proton_C;
                        }
                        if (isWest) {
                            h3D_W = h3D_protonsWest; hPt_W = hPtProtonWestX; hEta_W = hEtaProtonWestX; hRap_W = hRapidityProtonWestX;
                            hDistBeam_W = hDistanceBeam_proton_W; hDistEdge_W = hDistanceEdge_proton_W; hDCA_W = hProton_DCA_W;
                            MultiplicityProtonWestX++;

                            hEta_xi_EW[0] = hEta_xi0_proton_W; hRap_xi_EW[0] = hRapidity_xi0_proton_W;
                            hEta_xi_EW[1] = hEta_xi1_proton_W; hRap_xi_EW[1] = hRapidity_xi1_proton_W;
                            hEta_xi_EW[2] = hEta_xi2_proton_W; hRap_xi_EW[2] = hRapidity_xi2_proton_W;
                            hEta_xi_EW[3] = hEta_xi3_proton_W; hRap_xi_EW[3] = hRapidity_xi3_proton_W;
                            hEta_xi_EW[4] = hEta_xi4_proton_W; hRap_xi_EW[4] = hRapidity_xi4_proton_W;

                            hEta_xi_C[0] = hEta_xi0_proton_C; hRap_xi_C[0] = hRapidity_xi0_proton_C;
                            hEta_xi_C[1] = hEta_xi1_proton_C; hRap_xi_C[1] = hRapidity_xi1_proton_C;
                            hEta_xi_C[2] = hEta_xi2_proton_C; hRap_xi_C[2] = hRapidity_xi2_proton_C;
                            hEta_xi_C[3] = hEta_xi3_proton_C; hRap_xi_C[3] = hRapidity_xi3_proton_C;
                            hEta_xi_C[4] = hEta_xi4_proton_C; hRap_xi_C[4] = hRapidity_xi4_proton_C;
                        }

                    } else if (pCat == MID_P) {
                        hPt = hPtProtonX_mid_pt; hEta = hEtaProtonX_mid_pt; hRap = hRapidityProtonX_mid_pt;
                        hDistBeam = hDistanceBeam_proton_mid_pt; hDistEdge = hDistanceEdge_proton_mid_pt;
                        MultiplicityProtonX_mid_pt++;

                        if (isEast) {
                            hPt_E = hPtProtonEastX_mid_pt; hEta_E = hEtaProtonEastX_mid_pt; hRap_E = hRapidityProtonEastX_mid_pt;
                            hDistBeam_E = hDistanceBeam_proton_E_mid_pt; hDistEdge_E = hDistanceEdge_proton_E_mid_pt; hDCA_E = hProton_DCA_E_mid_pt;
                            MultiplicityProtonEastX_mid_pt++;
                            
                            hRap_xi_EW[0] = hRapidity_xi0_proton_E_mid_pt; hRap_xi_C[0] = hRapidity_xi0_proton_C_mid_pt;
                            hRap_xi_EW[1] = hRapidity_xi1_proton_E_mid_pt; hRap_xi_C[1] = hRapidity_xi1_proton_C_mid_pt;
                            hRap_xi_EW[2] = hRapidity_xi2_proton_E_mid_pt; hRap_xi_C[2] = hRapidity_xi2_proton_C_mid_pt;
                            hRap_xi_EW[3] = hRapidity_xi3_proton_E_mid_pt; hRap_xi_C[3] = hRapidity_xi3_proton_C_mid_pt;
                            hRap_xi_EW[4] = hRapidity_xi4_proton_E_mid_pt; hRap_xi_C[4] = hRapidity_xi4_proton_C_mid_pt;
                        }
                        if (isWest) {
                            hPt_W = hPtProtonWestX_mid_pt; hEta_W = hEtaProtonWestX_mid_pt; hRap_W = hRapidityProtonWestX_mid_pt;
                            hDistBeam_W = hDistanceBeam_proton_W_mid_pt; hDistEdge_W = hDistanceEdge_proton_W_mid_pt; hDCA_W = hProton_DCA_W_mid_pt;
                            MultiplicityProtonWestX_mid_pt++;
                            
                            hRap_xi_EW[0] = hRapidity_xi0_proton_W_mid_pt; hRap_xi_C[0] = hRapidity_xi0_proton_C_mid_pt;
                            hRap_xi_EW[1] = hRapidity_xi1_proton_W_mid_pt; hRap_xi_C[1] = hRapidity_xi1_proton_C_mid_pt;
                            hRap_xi_EW[2] = hRapidity_xi2_proton_W_mid_pt; hRap_xi_C[2] = hRapidity_xi2_proton_C_mid_pt;
                            hRap_xi_EW[3] = hRapidity_xi3_proton_W_mid_pt; hRap_xi_C[3] = hRapidity_xi3_proton_C_mid_pt;
                            hRap_xi_EW[4] = hRapidity_xi4_proton_W_mid_pt; hRap_xi_C[4] = hRapidity_xi4_proton_C_mid_pt;
                        }

                    } else { // HIGH_P
                        hPt = hPtProtonX_high_pt; hEta = hEtaProtonX_high_pt; hRap = hRapidityProtonX_high_pt;
                        hDistBeam = hDistanceBeam_proton_high_pt; hDistEdge = hDistanceEdge_proton_high_pt;
                        MultiplicityProtonX_high_pt++;

                        if (isEast) {
                            hPt_E = hPtProtonEastX_high_pt; hEta_E = hEtaProtonEastX_high_pt; hRap_E = hRapidityProtonEastX_high_pt;
                            hDistBeam_E = hDistanceBeam_proton_E_high_pt; hDistEdge_E = hDistanceEdge_proton_E_high_pt; hDCA_E = hProton_DCA_E_high_pt;
                            MultiplicityProtonEastX_high_pt++;
                            
                            hRap_xi_EW[0] = hRapidity_xi0_proton_E_high_pt; hRap_xi_C[0] = hRapidity_xi0_proton_C_high_pt;
                            hRap_xi_EW[1] = hRapidity_xi1_proton_E_high_pt; hRap_xi_C[1] = hRapidity_xi1_proton_C_high_pt;
                            hRap_xi_EW[2] = hRapidity_xi2_proton_E_high_pt; hRap_xi_C[2] = hRapidity_xi2_proton_C_high_pt;
                            hRap_xi_EW[3] = hRapidity_xi3_proton_E_high_pt; hRap_xi_C[3] = hRapidity_xi3_proton_C_high_pt;
                            hRap_xi_EW[4] = hRapidity_xi4_proton_E_high_pt; hRap_xi_C[4] = hRapidity_xi4_proton_C_high_pt;
                        }
                        if (isWest) {
                            hPt_W = hPtProtonWestX_high_pt; hEta_W = hEtaProtonWestX_high_pt; hRap_W = hRapidityProtonWestX_high_pt;
                            hDistBeam_W = hDistanceBeam_proton_W_high_pt; hDistEdge_W = hDistanceEdge_proton_W_high_pt; hDCA_W = hProton_DCA_W_high_pt;
                            MultiplicityProtonWestX_high_pt++;

                            hRap_xi_EW[0] = hRapidity_xi0_proton_W_high_pt; hRap_xi_C[0] = hRapidity_xi0_proton_C_high_pt;
                            hRap_xi_EW[1] = hRapidity_xi1_proton_W_high_pt; hRap_xi_C[1] = hRapidity_xi1_proton_C_high_pt;
                            hRap_xi_EW[2] = hRapidity_xi2_proton_W_high_pt; hRap_xi_C[2] = hRapidity_xi2_proton_C_high_pt;
                            hRap_xi_EW[3] = hRapidity_xi3_proton_W_high_pt; hRap_xi_C[3] = hRapidity_xi3_proton_C_high_pt;
                            hRap_xi_EW[4] = hRapidity_xi4_proton_W_high_pt; hRap_xi_C[4] = hRapidity_xi4_proton_C_high_pt;
                        }
                    }
                } else { // ANTIPROTONS
                    if (pCat == LOW_P) {
                        h3D = h3D_antiprotons; hPt = hPtNProtonX; hEta = hEtaNProtonX; hRap = hRapidityNProtonX;
                        hDistBeam = hDistanceBeam_antiproton; hDistEdge = hDistanceEdge_antiproton;
                        MultiplicityNProtonX++;

                        if (isEast) {
                            h3D_E = h3D_antiprotonsEast; hPt_E = hPtNProtonEastX; hEta_E = hEtaNProtonEastX; hRap_E = hRapidityNProtonEastX;
                            hDistBeam_E = hDistanceBeam_antiproton_E; hDistEdge_E = hDistanceEdge_antiproton_E; hDCA_E = hAntiproton_DCA_E;
                            MultiplicityNProtonEastX++;
                            
                            hEta_xi_EW[0] = hEta_xi0_antiproton_E; hRap_xi_EW[0] = hRapidity_xi0_antiproton_E;
                            hEta_xi_EW[1] = hEta_xi1_antiproton_E; hRap_xi_EW[1] = hRapidity_xi1_antiproton_E;
                            hEta_xi_EW[2] = hEta_xi2_antiproton_E; hRap_xi_EW[2] = hRapidity_xi2_antiproton_E;
                            hEta_xi_EW[3] = hEta_xi3_antiproton_E; hRap_xi_EW[3] = hRapidity_xi3_antiproton_E;
                            hEta_xi_EW[4] = hEta_xi4_antiproton_E; hRap_xi_EW[4] = hRapidity_xi4_antiproton_E;

                            hEta_xi_C[0] = hEta_xi0_antiproton_C; hRap_xi_C[0] = hRapidity_xi0_antiproton_C;
                            hEta_xi_C[1] = hEta_xi1_antiproton_C; hRap_xi_C[1] = hRapidity_xi1_antiproton_C;
                            hEta_xi_C[2] = hEta_xi2_antiproton_C; hRap_xi_C[2] = hRapidity_xi2_antiproton_C;
                            hEta_xi_C[3] = hEta_xi3_antiproton_C; hRap_xi_C[3] = hRapidity_xi3_antiproton_C;
                            hEta_xi_C[4] = hEta_xi4_antiproton_C; hRap_xi_C[4] = hRapidity_xi4_antiproton_C;
                        }
                        if (isWest) {
                            h3D_W = h3D_antiprotonsWest; hPt_W = hPtNProtonWestX; hEta_W = hEtaNProtonWestX; hRap_W = hRapidityNProtonWestX;
                            hDistBeam_W = hDistanceBeam_antiproton_W; hDistEdge_W = hDistanceEdge_antiproton_W; hDCA_W = hAntiproton_DCA_W;
                            MultiplicityNProtonWestX++;
                            
                            hEta_xi_EW[0] = hEta_xi0_antiproton_W; hRap_xi_EW[0] = hRapidity_xi0_antiproton_W;
                            hEta_xi_EW[1] = hEta_xi1_antiproton_W; hRap_xi_EW[1] = hRapidity_xi1_antiproton_W;
                            hEta_xi_EW[2] = hEta_xi2_antiproton_W; hRap_xi_EW[2] = hRapidity_xi2_antiproton_W;
                            hEta_xi_EW[3] = hEta_xi3_antiproton_W; hRap_xi_EW[3] = hRapidity_xi3_antiproton_W;
                            hEta_xi_EW[4] = hEta_xi4_antiproton_W; hRap_xi_EW[4] = hRapidity_xi4_antiproton_W;
                            
                            hEta_xi_C[0] = hEta_xi0_antiproton_C; hRap_xi_C[0] = hRapidity_xi0_antiproton_C;
                            hEta_xi_C[1] = hEta_xi1_antiproton_C; hRap_xi_C[1] = hRapidity_xi1_antiproton_C;
                            hEta_xi_C[2] = hEta_xi2_antiproton_C; hRap_xi_C[2] = hRapidity_xi2_antiproton_C;
                            hEta_xi_C[3] = hEta_xi3_antiproton_C; hRap_xi_C[3] = hRapidity_xi3_antiproton_C;
                            hEta_xi_C[4] = hEta_xi4_antiproton_C; hRap_xi_C[4] = hRapidity_xi4_antiproton_C;
                        }

                    } else if (pCat == MID_P) {
                        hPt = hPtNProtonX_mid_pt; hEta = hEtaNProtonX_mid_pt; hRap = hRapidityNProtonX_mid_pt;
                        hDistBeam = hDistanceBeam_antiproton_mid_pt; hDistEdge = hDistanceEdge_antiproton_mid_pt;
                        MultiplicityNProtonX_mid_pt++;

                        if (isEast) {
                            hPt_E = hPtNProtonEastX_mid_pt; hEta_E = hEtaNProtonEastX_mid_pt; hRap_E = hRapidityNProtonEastX_mid_pt;
                            hDistBeam_E = hDistanceBeam_antiproton_E_mid_pt; hDistEdge_E = hDistanceEdge_antiproton_E_mid_pt; hDCA_E = hAntiproton_DCA_E_mid_pt;
                            MultiplicityNProtonEastX_mid_pt++;
                            
                            hRap_xi_EW[0] = hRapidity_xi0_antiproton_E_mid_pt; hRap_xi_C[0] = hRapidity_xi0_antiproton_C_mid_pt;
                            hRap_xi_EW[1] = hRapidity_xi1_antiproton_E_mid_pt; hRap_xi_C[1] = hRapidity_xi1_antiproton_C_mid_pt;
                            hRap_xi_EW[2] = hRapidity_xi2_antiproton_E_mid_pt; hRap_xi_C[2] = hRapidity_xi2_antiproton_C_mid_pt;
                            hRap_xi_EW[3] = hRapidity_xi3_antiproton_E_mid_pt; hRap_xi_C[3] = hRapidity_xi3_antiproton_C_mid_pt;
                            hRap_xi_EW[4] = hRapidity_xi4_antiproton_E_mid_pt; hRap_xi_C[4] = hRapidity_xi4_antiproton_C_mid_pt;
                        }
                        if (isWest) {
                            hPt_W = hPtNProtonWestX_mid_pt; hEta_W = hEtaNProtonWestX_mid_pt; hRap_W = hRapidityNProtonWestX_mid_pt;
                            hDistBeam_W = hDistanceBeam_antiproton_W_mid_pt; hDistEdge_W = hDistanceEdge_antiproton_W_mid_pt; hDCA_W = hAntiproton_DCA_W_mid_pt;
                            MultiplicityNProtonWestX_mid_pt++;
                            
                            hRap_xi_EW[0] = hRapidity_xi0_antiproton_W_mid_pt; hRap_xi_C[0] = hRapidity_xi0_antiproton_C_mid_pt;
                            hRap_xi_EW[1] = hRapidity_xi1_antiproton_W_mid_pt; hRap_xi_C[1] = hRapidity_xi1_antiproton_C_mid_pt;
                            hRap_xi_EW[2] = hRapidity_xi2_antiproton_W_mid_pt; hRap_xi_C[2] = hRapidity_xi2_antiproton_C_mid_pt;
                            hRap_xi_EW[3] = hRapidity_xi3_antiproton_W_mid_pt; hRap_xi_C[3] = hRapidity_xi3_antiproton_C_mid_pt;
                            hRap_xi_EW[4] = hRapidity_xi4_antiproton_W_mid_pt; hRap_xi_C[4] = hRapidity_xi4_antiproton_C_mid_pt;
                        }

                    } else { // HIGH_P
                        hPt = hPtNProtonX_high_pt; hEta = hEtaNProtonX_high_pt; hRap = hRapidityNProtonX_high_pt;
                        hDistBeam = hDistanceBeam_antiproton_high_pt; hDistEdge = hDistanceEdge_antiproton_high_pt;
                        MultiplicityNProtonX_high_pt++;

                        if (isEast) {
                            hPt_E = hPtNProtonEastX_high_pt; hEta_E = hEtaNProtonEastX_high_pt; hRap_E = hRapidityNProtonEastX_high_pt;
                            hDistBeam_E = hDistanceBeam_antiproton_E_high_pt; hDistEdge_E = hDistanceEdge_antiproton_E_high_pt; hDCA_E = hAntiproton_DCA_E_high_pt;
                            MultiplicityNProtonEastX_high_pt++;

                            hRap_xi_EW[0] = hRapidity_xi0_antiproton_E_high_pt; hRap_xi_C[0] = hRapidity_xi0_antiproton_C_high_pt;
                            hRap_xi_EW[1] = hRapidity_xi1_antiproton_E_high_pt; hRap_xi_C[1] = hRapidity_xi1_antiproton_C_high_pt;
                            hRap_xi_EW[2] = hRapidity_xi2_antiproton_E_high_pt; hRap_xi_C[2] = hRapidity_xi2_antiproton_C_high_pt;
                            hRap_xi_EW[3] = hRapidity_xi3_antiproton_E_high_pt; hRap_xi_C[3] = hRapidity_xi3_antiproton_C_high_pt;
                            hRap_xi_EW[4] = hRapidity_xi4_antiproton_E_high_pt; hRap_xi_C[4] = hRapidity_xi4_antiproton_C_high_pt;
                        }
                        if (isWest) {
                            hPt_W = hPtNProtonWestX_high_pt; hEta_W = hEtaNProtonWestX_high_pt; hRap_W = hRapidityNProtonWestX_high_pt;
                            hDistBeam_W = hDistanceBeam_antiproton_W_high_pt; hDistEdge_W = hDistanceEdge_antiproton_W_high_pt; hDCA_W = hAntiproton_DCA_W_high_pt;
                            MultiplicityNProtonWestX_high_pt++;

                            hRap_xi_EW[0] = hRapidity_xi0_antiproton_W_high_pt; hRap_xi_C[0] = hRapidity_xi0_antiproton_C_high_pt;
                            hRap_xi_EW[1] = hRapidity_xi1_antiproton_W_high_pt; hRap_xi_C[1] = hRapidity_xi1_antiproton_C_high_pt;
                            hRap_xi_EW[2] = hRapidity_xi2_antiproton_W_high_pt; hRap_xi_C[2] = hRapidity_xi2_antiproton_C_high_pt;
                            hRap_xi_EW[3] = hRapidity_xi3_antiproton_W_high_pt; hRap_xi_C[3] = hRapidity_xi3_antiproton_C_high_pt;
                            hRap_xi_EW[4] = hRapidity_xi4_antiproton_W_high_pt; hRap_xi_C[4] = hRapidity_xi4_antiproton_C_high_pt;
                        }
                    }
                }

                if (h3D)       h3D->Fill(pt, eta, vz);
                if (hPt)       hPt->Fill(pt);
                if (hEta)      hEta->Fill(eta);
                if (hRap)      hRap->Fill(rapidity);
                if (hDistBeam) hDistBeam->Fill(distance_from_beam);
                if (hDistEdge) hDistEdge->Fill(distance_from_edge);

                if (isEast) {
                    if (h3D_E)       h3D_E->Fill(pt, eta, vz);
                    if (hPt_E)       hPt_E->Fill(pt);
                    if (hEta_E)      hEta_E->Fill(eta);
                    if (hRap_E)      hRap_E->Fill(rapidity);
                    if (hDistBeam_E) hDistBeam_E->Fill(distance_from_beam);
                    if (hDistEdge_E) hDistEdge_E->Fill(distance_from_edge);
                    if (hDCA_E)      hDCA_E->Fill(dca);

                    if (b >= 0 && b <= 4) {
                        if (hEta_xi_EW[b]) hEta_xi_EW[b]->Fill(eta);
                        if (hRap_xi_EW[b]) hRap_xi_EW[b]->Fill(rapidity);
                        if (hEta_xi_C[b])  hEta_xi_C[b]->Fill(-eta);
                        if (hRap_xi_C[b])  hRap_xi_C[b]->Fill(-rapidity);
                    }
                }

                if (isWest) {
                    if (h3D_W)       h3D_W->Fill(pt, eta, vz);
                    if (hPt_W)       hPt_W->Fill(pt);
                    if (hEta_W)      hEta_W->Fill(eta);
                    if (hRap_W)      hRap_W->Fill(rapidity);
                    if (hDistBeam_W) hDistBeam_W->Fill(distance_from_beam);
                    if (hDistEdge_W) hDistEdge_W->Fill(distance_from_edge);
                    if (hDCA_W)      hDCA_W->Fill(dca);

                    if (b >= 0 && b <= 4) {
                        if (hEta_xi_EW[b]) hEta_xi_EW[b]->Fill(eta);
                        if (hRap_xi_EW[b]) hRap_xi_EW[b]->Fill(rapidity);
                        if (hEta_xi_C[b])  hEta_xi_C[b]->Fill(eta);
                        if (hRap_xi_C[b])  hRap_xi_C[b]->Fill(rapidity);
                    }
                }
            }
        }
            
        HistNumOfPrimaryTracksToF->Fill(NumOfPrimaryTracksToF);
        hMultiplicityProtonX->Fill(MultiplicityProtonX);
        hMultiplicityNProtonX->Fill(MultiplicityNProtonX);
        hMultiplicityProtonX_high_pt->Fill(MultiplicityProtonX_high_pt);
        hMultiplicityNProtonX_high_pt->Fill(MultiplicityNProtonX_high_pt);

        if (isEast) {
            HistNumOfPToFEast->Fill(NumOfPrimaryTracksToF);
            hMultiplicityProtonEastX->Fill(MultiplicityProtonEastX);
            hMultiplicityNProtonEastX->Fill(MultiplicityNProtonEastX);
            hMultiplicityProtonEastX_high_pt->Fill(MultiplicityProtonEastX_high_pt);
            hMultiplicityNProtonEastX_high_pt->Fill(MultiplicityNProtonEastX_high_pt);
            if(proton_high_p_candidate!=0) {
                hPAntiP_Primaries_W->Fill(NumOfGoodPrimaryTracks-proton_high_p_candidate);
            }
        }
        if (isWest) {
            HistNumOfPToFWest->Fill(NumOfPrimaryTracksToF);
            hMultiplicityProtonWestX->Fill(MultiplicityProtonWestX);
            hMultiplicityNProtonWestX->Fill(MultiplicityNProtonWestX);
            hMultiplicityProtonWestX_high_pt->Fill(MultiplicityProtonWestX_high_pt);
            hMultiplicityNProtonWestX_high_pt->Fill(MultiplicityNProtonWestX_high_pt);
            if(proton_high_p_candidate!=0) { 
                hPAntiP_Primaries_E->Fill(NumOfGoodPrimaryTracks-proton_high_p_candidate);
            }
        }

    //    if (proton->branch() < 2)  { //east
    //    HistXiProtonEast->Fill(proton->xi(254.867));
    //    HistLogXiProtonEast->Fill(log(proton->xi(254.867)));
    //    HistPtProtonEast->Fill(proton->pt());
    //    HistEtaProtonEast->Fill(proton->eta());
    //    }
    //    if (proton->branch() > 1)  { //west
    //    HistXiProtonWest->Fill(proton->xi(254.867));
    //    HistLogXiProtonWest->Fill(log(proton->xi(254.867)));
    //    HistPtProtonWest->Fill(proton->pt());
    //    HistEtaProtonWest->Fill(proton->eta());
    //    }

    } // end of loop for events


    hEventsPassedSelection->Write();
    hTPClength->Write();
    hTPCdx->Write();
    hdxVsP->Write();

    for (int i = 0; i < nDx; ++i) {
        if (hRatioVsP_pion[i])   hRatioVsP_pion[i]->Write();
        if (hRatioVsP_kaon[i])   hRatioVsP_kaon[i]->Write();
        if (hRatioVsP_proton[i]) hRatioVsP_proton[i]->Write();
    }

    for (int i = 0; i < nDx; ++i) {
        if (hRatioTOFVsP_pion[i])   hRatioTOFVsP_pion[i]->Write();
        if (hRatioTOFVsP_kaon[i])   hRatioTOFVsP_kaon[i]->Write();
        if (hRatioTOFVsP_proton[i]) hRatioTOFVsP_proton[i]->Write();
    }
    
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
    protons_lower_momenta->cd();

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

    // hProton_Primaries_E->Write();
    // hProton_Primaries_W->Write();
    // hAntiproton_Primaries_E->Write();
    // hAntiproton_Primaries_W->Write();

    hPAntiP_Primaries_E->Write();
    hPAntiP_Primaries_W->Write();

    protons_mid_momenta->cd();

    hPtProtonEastX_mid_pt->Write();
    hPtProtonWestX_mid_pt->Write();
    hPtProtonX_mid_pt->Write();

    hPtNProtonEastX_mid_pt->Write();
    hPtNProtonWestX_mid_pt->Write();
    hPtNProtonX_mid_pt->Write();

    hEtaProtonEastX_mid_pt->Write();
    hEtaProtonWestX_mid_pt->Write();
    hEtaProtonX_mid_pt->Write();

    hEtaNProtonEastX_mid_pt->Write();
    hEtaNProtonWestX_mid_pt->Write();
    hEtaNProtonX_mid_pt->Write();

    hRapidityProtonEastX_mid_pt->Write();
    hRapidityProtonWestX_mid_pt->Write();
    hRapidityProtonX_mid_pt->Write();

    hRapidityNProtonEastX_mid_pt->Write();
    hRapidityNProtonWestX_mid_pt->Write();
    hRapidityNProtonX_mid_pt->Write();

    hMultiplicityProtonEastX_mid_pt->Write();
    hMultiplicityProtonWestX_mid_pt->Write();
    hMultiplicityProtonX_mid_pt->Write();

    hMultiplicityNProtonEastX_mid_pt->Write();
    hMultiplicityNProtonWestX_mid_pt->Write();
    hMultiplicityNProtonX_mid_pt->Write();

    hProton_DCA_E_mid_pt->Write();
    hProton_DCA_W_mid_pt->Write();
    hAntiproton_DCA_E_mid_pt->Write();
    hAntiproton_DCA_W_mid_pt->Write();

    hRapidity_xi0_proton_E_mid_pt->Write();
    hRapidity_xi1_proton_E_mid_pt->Write();
    hRapidity_xi2_proton_E_mid_pt->Write();
    hRapidity_xi3_proton_E_mid_pt->Write();
    hRapidity_xi4_proton_E_mid_pt->Write();

    hRapidity_xi0_proton_W_mid_pt->Write();
    hRapidity_xi1_proton_W_mid_pt->Write();
    hRapidity_xi2_proton_W_mid_pt->Write();
    hRapidity_xi3_proton_W_mid_pt->Write();
    hRapidity_xi4_proton_W_mid_pt->Write();

    hRapidity_xi0_proton_C_mid_pt->Write();
    hRapidity_xi1_proton_C_mid_pt->Write();
    hRapidity_xi2_proton_C_mid_pt->Write();
    hRapidity_xi3_proton_C_mid_pt->Write();
    hRapidity_xi4_proton_C_mid_pt->Write();

    hRapidity_xi0_antiproton_E_mid_pt->Write();
    hRapidity_xi1_antiproton_E_mid_pt->Write();
    hRapidity_xi2_antiproton_E_mid_pt->Write();
    hRapidity_xi3_antiproton_E_mid_pt->Write();
    hRapidity_xi4_antiproton_E_mid_pt->Write();

    hRapidity_xi0_antiproton_W_mid_pt->Write();
    hRapidity_xi1_antiproton_W_mid_pt->Write();
    hRapidity_xi2_antiproton_W_mid_pt->Write();
    hRapidity_xi3_antiproton_W_mid_pt->Write();
    hRapidity_xi4_antiproton_W_mid_pt->Write();

    hRapidity_xi0_antiproton_C_mid_pt->Write();
    hRapidity_xi1_antiproton_C_mid_pt->Write();
    hRapidity_xi2_antiproton_C_mid_pt->Write();
    hRapidity_xi3_antiproton_C_mid_pt->Write();
    hRapidity_xi4_antiproton_C_mid_pt->Write();

    protons_higher_momenta->cd();

    hPtProtonEastX_high_pt->Write();
    hPtProtonWestX_high_pt->Write();
    hPtProtonX_high_pt->Write();

    hPtNProtonEastX_high_pt->Write();
    hPtNProtonWestX_high_pt->Write();
    hPtNProtonX_high_pt->Write();

    hEtaProtonEastX_high_pt->Write();
    hEtaProtonWestX_high_pt->Write();
    hEtaProtonX_high_pt->Write();
    hEtaNProtonEastX_high_pt->Write();
    hEtaNProtonWestX_high_pt->Write();
    hEtaNProtonX_high_pt->Write();

    hRapidityProtonEastX_high_pt->Write();
    hRapidityProtonWestX_high_pt->Write();
    hRapidityProtonX_high_pt->Write();

    hRapidityNProtonEastX_high_pt->Write();
    hRapidityNProtonWestX_high_pt->Write();
    hRapidityNProtonX_high_pt->Write();


    hMultiplicityProtonEastX_high_pt->Write();
    hMultiplicityProtonWestX_high_pt->Write();
    hMultiplicityProtonX_high_pt->Write();

    hMultiplicityNProtonEastX_high_pt->Write();
    hMultiplicityNProtonWestX_high_pt->Write();
    hMultiplicityNProtonX_high_pt->Write();

    hProton_DCA_E_high_pt->Write();
    hProton_DCA_W_high_pt->Write();
    hAntiproton_DCA_E_high_pt->Write();
    hAntiproton_DCA_W_high_pt->Write();

    hProton_DCA_E_pm->Write();
    hProton_DCA_W_pm->Write();
    hAntiproton_DCA_E_pm->Write();
    hAntiproton_DCA_W_pm->Write();

    
    hRapidity_xi0_proton_E_high_pt->Write();
    hRapidity_xi1_proton_E_high_pt->Write();
    hRapidity_xi2_proton_E_high_pt->Write();
    hRapidity_xi3_proton_E_high_pt->Write();
    hRapidity_xi4_proton_E_high_pt->Write();

    hRapidity_xi0_proton_W_high_pt->Write();
    hRapidity_xi1_proton_W_high_pt->Write();
    hRapidity_xi2_proton_W_high_pt->Write();
    hRapidity_xi3_proton_W_high_pt->Write();
    hRapidity_xi4_proton_W_high_pt->Write();

    hRapidity_xi0_proton_C_high_pt->Write();
    hRapidity_xi1_proton_C_high_pt->Write();
    hRapidity_xi2_proton_C_high_pt->Write();
    hRapidity_xi3_proton_C_high_pt->Write();
    hRapidity_xi4_proton_C_high_pt->Write();

    hRapidity_xi0_antiproton_E_high_pt->Write();
    hRapidity_xi1_antiproton_E_high_pt->Write();
    hRapidity_xi2_antiproton_E_high_pt->Write();
    hRapidity_xi3_antiproton_E_high_pt->Write();
    hRapidity_xi4_antiproton_E_high_pt->Write();

    hRapidity_xi0_antiproton_W_high_pt->Write();
    hRapidity_xi1_antiproton_W_high_pt->Write();
    hRapidity_xi2_antiproton_W_high_pt->Write();
    hRapidity_xi3_antiproton_W_high_pt->Write();
    hRapidity_xi4_antiproton_W_high_pt->Write();

    hRapidity_xi0_antiproton_C_high_pt->Write();
    hRapidity_xi1_antiproton_C_high_pt->Write();
    hRapidity_xi2_antiproton_C_high_pt->Write();
    hRapidity_xi3_antiproton_C_high_pt->Write();
    hRapidity_xi4_antiproton_C_high_pt->Write();

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

    hDistanceBeam_proton_mid_pt->Write();
    hDistanceBeam_antiproton_mid_pt->Write();

    hDistanceEdge_proton_mid_pt->Write();
    hDistanceEdge_antiproton_mid_pt->Write();

    hDistanceBeam_proton_E_mid_pt->Write();
    hDistanceBeam_proton_W_mid_pt->Write();
    hDistanceBeam_antiproton_E_mid_pt->Write();
    hDistanceBeam_antiproton_W_mid_pt->Write();

    hDistanceEdge_proton_E_mid_pt->Write();
    hDistanceEdge_proton_W_mid_pt->Write();
    hDistanceEdge_antiproton_E_mid_pt->Write();
    hDistanceEdge_antiproton_W_mid_pt->Write();

    hDistanceBeam_proton_high_pt->Write();
    hDistanceBeam_antiproton_high_pt->Write();

    hDistanceEdge_proton_high_pt->Write();
    hDistanceEdge_antiproton_high_pt->Write();

    hDistanceBeam_proton_E_high_pt->Write();
    hDistanceBeam_proton_W_high_pt->Write();
    hDistanceBeam_antiproton_E_high_pt->Write();
    hDistanceBeam_antiproton_W_high_pt->Write();

    hDistanceEdge_proton_E_high_pt->Write();
    hDistanceEdge_proton_W_high_pt->Write();
    hDistanceEdge_antiproton_E_high_pt->Write();
    hDistanceEdge_antiproton_W_high_pt->Write();

    hNumOfParticlesPassSelectionWest->Write();
    hNumOfParticlesPassSelectionEast->Write();

    TOF_Analysis->cd();

    hDeltaT_lambda_proton_pion->Write();

    deltat_anchors->cd();

    hDeltaT_global_anchor_Kppi_high_pt->Write();
    hDeltaT_global_anchor_KpK_high_pt->Write();
    hDeltaT_global_anchor_Kpp_high_pt->Write();
    hDeltaT_global_anchor_Knpi_high_pt->Write();
    hDeltaT_global_anchor_KnK_high_pt->Write();
    hDeltaT_global_anchor_Knp_high_pt->Write();
    hDeltaT_global_anchor_Kppi_mid_pt->Write();
    hDeltaT_global_anchor_KpK_mid_pt->Write();
    hDeltaT_global_anchor_Kpp_mid_pt->Write();
    hDeltaT_global_anchor_Knpi_mid_pt->Write();
    hDeltaT_global_anchor_KnK_mid_pt->Write();
    hDeltaT_global_anchor_Knp_mid_pt->Write();
    hDeltaT_global_anchor_Kppi->Write();
    hDeltaT_global_anchor_KpK->Write();
    hDeltaT_global_anchor_Kpp->Write();
    hDeltaT_global_anchor_Knpi->Write();
    hDeltaT_global_anchor_KnK->Write();
    hDeltaT_global_anchor_Knp->Write();

    hDeltaT_global_anchor_Pippi_high_pt->Write();
    hDeltaT_global_anchor_PipK_high_pt->Write();
    hDeltaT_global_anchor_Pipp_high_pt->Write();
    hDeltaT_global_anchor_Pinpi_high_pt->Write();
    hDeltaT_global_anchor_PinK_high_pt->Write();
    hDeltaT_global_anchor_Pinp_high_pt->Write();
    hDeltaT_global_anchor_Pippi_mid_pt->Write();
    hDeltaT_global_anchor_PipK_mid_pt->Write();
    hDeltaT_global_anchor_Pipp_mid_pt->Write();
    hDeltaT_global_anchor_Pinpi_mid_pt->Write();
    hDeltaT_global_anchor_PinK_mid_pt->Write();
    hDeltaT_global_anchor_Pinp_mid_pt->Write();
    hDeltaT_global_anchor_Pippi->Write();
    hDeltaT_global_anchor_PipK->Write();
    hDeltaT_global_anchor_Pipp->Write();
    hDeltaT_global_anchor_Pinpi->Write();
    hDeltaT_global_anchor_PinK->Write();
    hDeltaT_global_anchor_Pinp->Write();

    hDeltaT_global_anchor_Pppi_high_pt->Write();
    hDeltaT_global_anchor_PpK_high_pt->Write();
    hDeltaT_global_anchor_Ppp_high_pt->Write();
    hDeltaT_global_anchor_Pnpi_high_pt->Write();
    hDeltaT_global_anchor_PnK_high_pt->Write();
    hDeltaT_global_anchor_Pnp_high_pt->Write();
    hDeltaT_global_anchor_Pppi_mid_pt->Write();
    hDeltaT_global_anchor_PpK_mid_pt->Write();
    hDeltaT_global_anchor_Ppp_mid_pt->Write();
    hDeltaT_global_anchor_Pnpi_mid_pt->Write();
    hDeltaT_global_anchor_PnK_mid_pt->Write();
    hDeltaT_global_anchor_Pnp_mid_pt->Write();
    hDeltaT_global_anchor_Pppi->Write();
    hDeltaT_global_anchor_PpK->Write();
    hDeltaT_global_anchor_Ppp->Write();
    hDeltaT_global_anchor_Pnpi->Write();
    hDeltaT_global_anchor_PnK->Write();
    hDeltaT_global_anchor_Pnp->Write();

    for (size_t i = 0; i < particles.size(); ++i) {
        for (size_t j = 0; j < partners.size(); ++j) {
            for (size_t k = 0; k < mom_bins.size(); ++k) {
                if (hDeltaT[i][j][k]) {
                    hDeltaT[i][j][k]->Write();
                }
            }
        }
    }

    deltat_protons_lower_momenta->cd();

    hDeltaT1vsT2_protons->Write();
    hDeltaT1vsT2_ppbar->Write();
    hDeltaT1vsT5_protons->Write();
    hDeltaT1vsT5_antiprotons->Write();
    hDeltaT1vsT2_antiprotons->Write();
    hDeltaT1vsT2_protons_east->Write();
    hDeltaT1vsT2_antiprotons_east->Write();
    hDeltaT1vsT2_protons_west->Write();
    hDeltaT1vsT2_antiprotons_west->Write();

    hDeltaT_pionpion->Write();
    hDeltaT_kaonpion->Write();
    hDeltaT_protonproton->Write();

    deltat_protons_mid_momenta->cd();

    hDeltaT1vsT5_protons_mid_p->Write();
    hDeltaT1vsT5_antiprotons_mid_p->Write();
    hDeltaT1vsT2_ppbar_mid_p->Write();
    hDeltaT1vsT2_protons_mid_p->Write();
    hDeltaT1vsT2_antiprotons_mid_p->Write();

    hDeltaT1vsT2_protons_east_mid_p->Write();
    hDeltaT1vsT2_antiprotons_east_mid_p->Write();

    hDeltaT1vsT2_protons_west_mid_p->Write();
    hDeltaT1vsT2_antiprotons_west_mid_p->Write();

    hDeltaT_pionpion_mid_p->Write();
    hDeltaT_kaonpion_mid_p->Write();
    hDeltaT_protonproton_mid_p->Write();

    deltat_protons_higher_momenta->cd();
    hDeltaT1vsT5_protons_high_p->Write();
    hDeltaT1vsT5_antiprotons_high_p->Write();
    hDeltaT1vsT2_ppbar_high_p->Write();
    hDeltaT1vsT2_protons_high_p->Write();
    hDeltaT1vsT2_antiprotons_high_p->Write();
    hDeltaT1vsT2_protons_east_high_p->Write();
    hDeltaT1vsT2_antiprotons_east_high_p->Write();
    hDeltaT1vsT2_protons_west_high_p->Write();
    hDeltaT1vsT2_antiprotons_west_high_p->Write();

    hDeltaT_pionpion_high_p->Write();
    hDeltaT_kaonpion_high_p->Write();
    hDeltaT_protonproton_high_p->Write();

    kaons_lower_momenta->cd();

    hRapidity_xi0_kaon_plus_E->Write();
    hRapidity_xi1_kaon_plus_E->Write();
    hRapidity_xi2_kaon_plus_E->Write();
    hRapidity_xi3_kaon_plus_E->Write();
    hRapidity_xi4_kaon_plus_E->Write();

    hRapidity_xi0_kaon_plus_W->Write();
    hRapidity_xi1_kaon_plus_W->Write();
    hRapidity_xi2_kaon_plus_W->Write();
    hRapidity_xi3_kaon_plus_W->Write();
    hRapidity_xi4_kaon_plus_W->Write();

    hRapidity_xi0_kaon_plus_C->Write();
    hRapidity_xi1_kaon_plus_C->Write();
    hRapidity_xi2_kaon_plus_C->Write();
    hRapidity_xi3_kaon_plus_C->Write();
    hRapidity_xi4_kaon_plus_C->Write();

    hRapidity_xi0_kaon_minus_E->Write();
    hRapidity_xi1_kaon_minus_E->Write();
    hRapidity_xi2_kaon_minus_E->Write();
    hRapidity_xi3_kaon_minus_E->Write();
    hRapidity_xi4_kaon_minus_E->Write();

    hRapidity_xi0_kaon_minus_W->Write();
    hRapidity_xi1_kaon_minus_W->Write();
    hRapidity_xi2_kaon_minus_W->Write();
    hRapidity_xi3_kaon_minus_W->Write();
    hRapidity_xi4_kaon_minus_W->Write();

    hRapidity_xi0_kaon_minus_C->Write();
    hRapidity_xi1_kaon_minus_C->Write();
    hRapidity_xi2_kaon_minus_C->Write();
    hRapidity_xi3_kaon_minus_C->Write();
    hRapidity_xi4_kaon_minus_C->Write();

    deltat_kaons_lower_momenta->cd();

    hDeltaT1vsT2_kaons_p->Write();
    hDeltaT1vsT2_kaons->Write();
    hDeltaT1vsT5_kaons_p->Write();
    hDeltaT1vsT2_kaons_p_east->Write();
    hDeltaT1vsT2_kaons_p_west->Write();

    hDeltaT1vsT2_kaons_n->Write();
    hDeltaT1vsT5_kaons_n->Write();
    hDeltaT1vsT2_kaons_n_east->Write();
    hDeltaT1vsT2_kaons_n_west->Write();

    kaons_mid_momenta->cd();

    hRapidity_xi0_kaon_plus_E_mid_pt->Write();
    hRapidity_xi1_kaon_plus_E_mid_pt->Write();
    hRapidity_xi2_kaon_plus_E_mid_pt->Write();
    hRapidity_xi3_kaon_plus_E_mid_pt->Write();
    hRapidity_xi4_kaon_plus_E_mid_pt->Write();

    hRapidity_xi0_kaon_plus_W_mid_pt->Write();
    hRapidity_xi1_kaon_plus_W_mid_pt->Write();
    hRapidity_xi2_kaon_plus_W_mid_pt->Write();
    hRapidity_xi3_kaon_plus_W_mid_pt->Write();
    hRapidity_xi4_kaon_plus_W_mid_pt->Write();

    hRapidity_xi0_kaon_plus_C_mid_pt->Write();
    hRapidity_xi1_kaon_plus_C_mid_pt->Write();
    hRapidity_xi2_kaon_plus_C_mid_pt->Write();
    hRapidity_xi3_kaon_plus_C_mid_pt->Write();
    hRapidity_xi4_kaon_plus_C_mid_pt->Write();

    hRapidity_xi0_kaon_minus_E_mid_pt->Write();
    hRapidity_xi1_kaon_minus_E_mid_pt->Write();
    hRapidity_xi2_kaon_minus_E_mid_pt->Write();
    hRapidity_xi3_kaon_minus_E_mid_pt->Write();
    hRapidity_xi4_kaon_minus_E_mid_pt->Write();

    hRapidity_xi0_kaon_minus_W_mid_pt->Write();
    hRapidity_xi1_kaon_minus_W_mid_pt->Write();
    hRapidity_xi2_kaon_minus_W_mid_pt->Write();
    hRapidity_xi3_kaon_minus_W_mid_pt->Write();
    hRapidity_xi4_kaon_minus_W_mid_pt->Write();

    hRapidity_xi0_kaon_minus_C_mid_pt->Write();
    hRapidity_xi1_kaon_minus_C_mid_pt->Write();
    hRapidity_xi2_kaon_minus_C_mid_pt->Write();
    hRapidity_xi3_kaon_minus_C_mid_pt->Write();
    hRapidity_xi4_kaon_minus_C_mid_pt->Write();

    deltat_kaons_mid_momenta->cd();

    hDeltaT1vsT2_kaons_p_mid_p->Write();
    hDeltaT1vsT2_kaons_mid_p->Write();
    hDeltaT1vsT5_kaons_p_mid_p->Write();
    hDeltaT1vsT5_kaons_n_mid_p->Write();
    hDeltaT1vsT2_kaons_n_mid_p->Write();
    hDeltaT1vsT2_kaons_p_east_mid_p->Write();
    hDeltaT1vsT2_kaons_n_east_mid_p->Write();
    hDeltaT1vsT2_kaons_p_west_mid_p->Write();
    hDeltaT1vsT2_kaons_n_west_mid_p->Write();

    kaons_higher_momenta->cd();

    hRapidity_xi0_kaon_plus_E_high_pt->Write();
    hRapidity_xi1_kaon_plus_E_high_pt->Write();
    hRapidity_xi2_kaon_plus_E_high_pt->Write();
    hRapidity_xi3_kaon_plus_E_high_pt->Write();
    hRapidity_xi4_kaon_plus_E_high_pt->Write();

    hRapidity_xi0_kaon_plus_W_high_pt->Write();
    hRapidity_xi1_kaon_plus_W_high_pt->Write();
    hRapidity_xi2_kaon_plus_W_high_pt->Write();
    hRapidity_xi3_kaon_plus_W_high_pt->Write();
    hRapidity_xi4_kaon_plus_W_high_pt->Write();

    hRapidity_xi0_kaon_plus_C_high_pt->Write();
    hRapidity_xi1_kaon_plus_C_high_pt->Write();
    hRapidity_xi2_kaon_plus_C_high_pt->Write();
    hRapidity_xi3_kaon_plus_C_high_pt->Write();
    hRapidity_xi4_kaon_plus_C_high_pt->Write();

    // Kaon Minus (High Pt)
    hRapidity_xi0_kaon_minus_E_high_pt->Write();
    hRapidity_xi1_kaon_minus_E_high_pt->Write();
    hRapidity_xi2_kaon_minus_E_high_pt->Write();
    hRapidity_xi3_kaon_minus_E_high_pt->Write();
    hRapidity_xi4_kaon_minus_E_high_pt->Write();

    hRapidity_xi0_kaon_minus_W_high_pt->Write();
    hRapidity_xi1_kaon_minus_W_high_pt->Write();
    hRapidity_xi2_kaon_minus_W_high_pt->Write();
    hRapidity_xi3_kaon_minus_W_high_pt->Write();
    hRapidity_xi4_kaon_minus_W_high_pt->Write();

    hRapidity_xi0_kaon_minus_C_high_pt->Write();
    hRapidity_xi1_kaon_minus_C_high_pt->Write();
    hRapidity_xi2_kaon_minus_C_high_pt->Write();
    hRapidity_xi3_kaon_minus_C_high_pt->Write();
    hRapidity_xi4_kaon_minus_C_high_pt->Write();

    deltat_kaons_higher_momenta->cd();

    hDeltaT1vsT2_kaons_p_high_p->Write();
    hDeltaT1vsT2_kaons_high_p->Write();
    hDeltaT1vsT5_kaons_p_high_p->Write();
    hDeltaT1vsT5_kaons_n_high_p->Write();
    hDeltaT1vsT2_kaons_n_high_p->Write();
    hDeltaT1vsT2_kaons_p_east_high_p->Write();
    hDeltaT1vsT2_kaons_n_east_high_p->Write();
    hDeltaT1vsT2_kaons_p_west_high_p->Write();
    hDeltaT1vsT2_kaons_n_west_high_p->Write();

    pions_lower_momenta->cd();

    hRapidity_xi0_pion_plus_E->Write();
    hRapidity_xi1_pion_plus_E->Write();
    hRapidity_xi2_pion_plus_E->Write();
    hRapidity_xi3_pion_plus_E->Write();
    hRapidity_xi4_pion_plus_E->Write();

    hRapidity_xi0_pion_plus_W->Write();
    hRapidity_xi1_pion_plus_W->Write();
    hRapidity_xi2_pion_plus_W->Write();
    hRapidity_xi3_pion_plus_W->Write();
    hRapidity_xi4_pion_plus_W->Write();

    hRapidity_xi0_pion_plus_C->Write();
    hRapidity_xi1_pion_plus_C->Write();
    hRapidity_xi2_pion_plus_C->Write();
    hRapidity_xi3_pion_plus_C->Write();
    hRapidity_xi4_pion_plus_C->Write();

    hRapidity_xi0_pion_minus_E->Write();
    hRapidity_xi1_pion_minus_E->Write();
    hRapidity_xi2_pion_minus_E->Write();
    hRapidity_xi3_pion_minus_E->Write();
    hRapidity_xi4_pion_minus_E->Write();

    hRapidity_xi0_pion_minus_W->Write();
    hRapidity_xi1_pion_minus_W->Write();
    hRapidity_xi2_pion_minus_W->Write();
    hRapidity_xi3_pion_minus_W->Write();
    hRapidity_xi4_pion_minus_W->Write();

    hRapidity_xi0_pion_minus_C->Write();
    hRapidity_xi1_pion_minus_C->Write();
    hRapidity_xi2_pion_minus_C->Write();
    hRapidity_xi3_pion_minus_C->Write();
    hRapidity_xi4_pion_minus_C->Write();

    deltat_pions_lower_momenta->cd();

    hDeltaT1vsT2_pions_p->Write();
    hDeltaT1vsT2_pions->Write();
    hDeltaT1vsT5_pions_p->Write();
    hDeltaT1vsT2_pions_p_east->Write();
    hDeltaT1vsT2_pions_p_west->Write();

    hDeltaT1vsT2_pions_n->Write();
    hDeltaT1vsT5_pions_n->Write();
    hDeltaT1vsT2_pions_n_east->Write();
    hDeltaT1vsT2_pions_n_west->Write();

    pions_mid_momenta->cd();

    hRapidity_xi0_pion_plus_E_mid_pt->Write();
    hRapidity_xi1_pion_plus_E_mid_pt->Write();
    hRapidity_xi2_pion_plus_E_mid_pt->Write();
    hRapidity_xi3_pion_plus_E_mid_pt->Write();
    hRapidity_xi4_pion_plus_E_mid_pt->Write();

    hRapidity_xi0_pion_plus_W_mid_pt->Write();
    hRapidity_xi1_pion_plus_W_mid_pt->Write();
    hRapidity_xi2_pion_plus_W_mid_pt->Write();
    hRapidity_xi3_pion_plus_W_mid_pt->Write();
    hRapidity_xi4_pion_plus_W_mid_pt->Write();

    hRapidity_xi0_pion_plus_C_mid_pt->Write();
    hRapidity_xi1_pion_plus_C_mid_pt->Write();
    hRapidity_xi2_pion_plus_C_mid_pt->Write();
    hRapidity_xi3_pion_plus_C_mid_pt->Write();
    hRapidity_xi4_pion_plus_C_mid_pt->Write();

    hRapidity_xi0_pion_minus_E_mid_pt->Write();
    hRapidity_xi1_pion_minus_E_mid_pt->Write();
    hRapidity_xi2_pion_minus_E_mid_pt->Write();
    hRapidity_xi3_pion_minus_E_mid_pt->Write();
    hRapidity_xi4_pion_minus_E_mid_pt->Write();

    hRapidity_xi0_pion_minus_W_mid_pt->Write();
    hRapidity_xi1_pion_minus_W_mid_pt->Write();
    hRapidity_xi2_pion_minus_W_mid_pt->Write();
    hRapidity_xi3_pion_minus_W_mid_pt->Write();
    hRapidity_xi4_pion_minus_W_mid_pt->Write();

    hRapidity_xi0_pion_minus_C_mid_pt->Write();
    hRapidity_xi1_pion_minus_C_mid_pt->Write();
    hRapidity_xi2_pion_minus_C_mid_pt->Write();
    hRapidity_xi3_pion_minus_C_mid_pt->Write();
    hRapidity_xi4_pion_minus_C_mid_pt->Write();

    deltat_pions_mid_momenta->cd();

    hDeltaT1vsT2_pions_p_mid_p->Write();
    hDeltaT1vsT2_pions_mid_p->Write();
    hDeltaT1vsT5_pions_p_mid_p->Write();
    hDeltaT1vsT5_pions_n_mid_p->Write();
    hDeltaT1vsT2_pions_n_mid_p->Write();
    hDeltaT1vsT2_pions_p_east_mid_p->Write();
    hDeltaT1vsT2_pions_n_east_mid_p->Write();
    hDeltaT1vsT2_pions_p_west_mid_p->Write();
    hDeltaT1vsT2_pions_n_west_mid_p->Write();

    pions_higher_momenta->cd();

    hRapidity_xi0_pion_plus_E_high_pt->Write();
    hRapidity_xi1_pion_plus_E_high_pt->Write();
    hRapidity_xi2_pion_plus_E_high_pt->Write();
    hRapidity_xi3_pion_plus_E_high_pt->Write();
    hRapidity_xi4_pion_plus_E_high_pt->Write();

    hRapidity_xi0_pion_plus_W_high_pt->Write();
    hRapidity_xi1_pion_plus_W_high_pt->Write();
    hRapidity_xi2_pion_plus_W_high_pt->Write();
    hRapidity_xi3_pion_plus_W_high_pt->Write();
    hRapidity_xi4_pion_plus_W_high_pt->Write();

    hRapidity_xi0_pion_plus_C_high_pt->Write();
    hRapidity_xi1_pion_plus_C_high_pt->Write();
    hRapidity_xi2_pion_plus_C_high_pt->Write();
    hRapidity_xi3_pion_plus_C_high_pt->Write();
    hRapidity_xi4_pion_plus_C_high_pt->Write();

    hRapidity_xi0_pion_minus_E_high_pt->Write();
    hRapidity_xi1_pion_minus_E_high_pt->Write();
    hRapidity_xi2_pion_minus_E_high_pt->Write();
    hRapidity_xi3_pion_minus_E_high_pt->Write();
    hRapidity_xi4_pion_minus_E_high_pt->Write();

    hRapidity_xi0_pion_minus_W_high_pt->Write();
    hRapidity_xi1_pion_minus_W_high_pt->Write();
    hRapidity_xi2_pion_minus_W_high_pt->Write();
    hRapidity_xi3_pion_minus_W_high_pt->Write();
    hRapidity_xi4_pion_minus_W_high_pt->Write();

    hRapidity_xi0_pion_minus_C_high_pt->Write();
    hRapidity_xi1_pion_minus_C_high_pt->Write();
    hRapidity_xi2_pion_minus_C_high_pt->Write();
    hRapidity_xi3_pion_minus_C_high_pt->Write();
    hRapidity_xi4_pion_minus_C_high_pt->Write();

    deltat_pions_higher_momenta->cd();

    hDeltaT1vsT2_pions_p_high_p->Write();
    hDeltaT1vsT2_pions_high_p->Write();
    hDeltaT1vsT5_pions_p_high_p->Write();
    hDeltaT1vsT5_pions_n_high_p->Write();
    hDeltaT1vsT2_pions_n_high_p->Write();
    hDeltaT1vsT2_pions_p_east_high_p->Write();
    hDeltaT1vsT2_pions_n_east_high_p->Write();
    hDeltaT1vsT2_pions_p_west_high_p->Write();
    hDeltaT1vsT2_pions_n_west_high_p->Write();

    outTree->Write("", TObject::kOverwrite);
    Long64_t treeWestCount = outTree->GetEntries("TreeIsEast == 0"); // lub "TreeIsWest == 1"
    double histWestCount   = hEventsPassedSelection->GetBinContent(2); // bin 2 dla isWest

    std::cout << "Liczba przypadków West w Tree: " << treeWestCount << std::endl;
    std::cout << "Liczba przypadków West w Hist (bin 2): " << histWestCount << std::endl;
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
    bool mass_cut = (L.m() > (1.1157-0.00487) && L.m() < 1.1157+0.00487); //1.115683
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

int getMomentumBin(double p) {
    if (p >= 0.2 && p < 0.3) return 0; // p1
    if (p >= 0.3 && p < 0.4) return 1; // p2
    if (p >= 0.4 && p < 0.5) return 2; // p3
    // Od 0.5 skok co 0.2
    if (p >= 0.5 && p < 0.7) return 3; // p4
    if (p >= 0.7 && p < 0.9) return 4; // p5
    if (p >= 0.9 && p < 1.1) return 5; // p6
    if (p >= 1.1 && p < 1.3) return 6; // p7
    if (p >= 1.3 && p < 1.5) return 7; // p8
    if (p >= 1.5 && p < 1.7) return 8; // p9
    if (p >= 1.7)            return 9;
    return -1;
}

double GetTheoreticaldEdx(double p, double massHypothesis, double dx) {
    static dEdxParameterization dEdxModel("p10");

    double bg = p / massHypothesis;
    double log10_bg = TMath::Log10(bg);
    double log2_dx = TMath::Log2(dx);

    return dEdxModel.GetI70(log10_bg, log2_dx);
}