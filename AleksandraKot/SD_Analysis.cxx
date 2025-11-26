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
    TH1D* HistEtaTracksWest = new TH1D("HistEtaTracksWest", "Eta of the reconstructed tracks when proton on west; #eta; # events", 100, -3, 3);
    TH1D* HistEtaTracksEast = new TH1D("HistEtaTracksEast", "Eta of the reconstructed tracks when proton on east; #eta; # events", 100, -3, 3);
    TH1D* HistPtTracksWestCut = new TH1D("HistPtTracksWestCut", "pT of the reconstructed tracks when proton on west; pT [GeV]; # events", 100, 0, 4);
    TH1D* HistPtTracksEastCut = new TH1D("HistPtTracksEastCut", "pT of the reconstructed tracks when proton on east; pT [GeV]; # events", 100, 0, 4);
    TH1D* HistEtaTracksWestCut = new TH1D("HistEtaTracksWestCut", "Eta of the reconstructed tracks when proton on west; #eta; # events", 100, -3, 3);
    TH1D* HistEtaTracksEastCut = new TH1D("HistEtaTracksEastCut", "Eta of the reconstructed tracks when proton on east; #eta; # events", 100, -3, 3);
    
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
    TH1D* hEtaProtonEastX = new TH1D("hEtaProtonEastX","#eta of protons (East); #eta; # events",40,-1,1);
    TH1D* hEtaProtonWestX = new TH1D("hEtaProtonWestX","#eta of protons (West); #eta; # events",40,-1,1);
    TH1D* hEtaProtonX = new TH1D("hEtaProtonX","#eta of proton; #eta; # events",40,-1,1);

    TH1D* hEtaNProtonEastX = new TH1D("hEtaNProtonEastX","#eta of antiprotons (East); #eta; # events",40,-1,1);
    TH1D* hEtaNProtonWestX = new TH1D("hEtaNProtonWestX","#eta of antiprotons (West); #eta; # events",40,-1,1);
    TH1D* hEtaNProtonX = new TH1D("hEtaNProtonX","#eta of antiprotons; #eta; # events",40,-1,1);

    TH1D* hMultiplicityProtonEastX = new TH1D("hMultiplicityProtonEastX", "Multiplicity of protons (East); N_{p}; # events", 30, 0, 30);
    TH1D* hMultiplicityProtonWestX = new TH1D("hMultiplicityProtonWestX", "Multiplicity of protons (West); N_{p}; # events", 30, 0, 30);
    TH1D* hMultiplicityProtonX = new TH1D("hMultiplicityProtonX", "Multiplicity of protons; N_{p}; # events", 30, 0, 30);

    TH1D* hMultiplicityNProtonEastX = new TH1D("hMultiplicityNProtonEastX", "Multiplicity of antiprotons (East); N_{#bar{p}}; # events", 30, 0, 30);
    TH1D* hMultiplicityNProtonWestX = new TH1D("hMultiplicityNProtonWestX", "Multiplicity of antiprotons (West); N_{#bar{p}}; # events", 30, 0, 30);
    TH1D* hMultiplicityNProtonX = new TH1D("hMultiplicityNProtonX", "Multiplicity of antiprotons; N_{#bar{p}}; # events", 30, 0, 30);

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


    

    for (Long64_t i = 0; i < chain->GetEntries(); ++i) {
        chain->GetEntry(i);
        if( upcEvt->getNumberOfVertices() != 1) continue;
        if (rpEvt->getNumberOfTracks() != 1 ) continue;

        //   cout << upcEvt->getEventNumber() << " " 
        //        << upcEvt->getNPrimVertices() <<  endl; //1, 2, 3

        StUPCRpsTrack *proton = rpEvt->getTrack(0);

        bool isEast = (proton->branch() < 2);
        bool isWest = (proton->branch() > 1);

        double xi_proton = proton->xi(254.867);
        int b=getXiBin(xi_proton);

        int NumOfPrimaryTracksToF = 0;

        int MultiplicityProtonX = 0;
        int MultiplicityNProtonX=0;
        int MultiplicityProtonEastX=0;
        int MultiplicityNProtonEastX=0;
        int MultiplicityProtonWestX=0;
        int MultiplicityNProtonWestX=0;

        //vz in mm
        double vz = upcEvt->getVertex(0)->getPosZ();
        hVz_check->Fill(vz);
        if(fabs(vz)>100.0) continue; //vz cut

        for (int j = 0; j < upcEvt->getNumberOfTracks(); j++) {
            if(   upcEvt->getTrack(j)->getNhits()>15  
                && upcEvt->getTrack(j)->getFlag(StUPCTrack::kTof) 
                && upcEvt->getTrack(j)->getFlag(StUPCTrack::kPrimary)) {

                if (upcEvt->getTrack(j)->getPt()>0.2  && abs(upcEvt->getTrack(j)->getEta())<0.9){ // fiducial cuts
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
                        if(upcEvt->getTrack(j)->getCharge()>0){
                            h3D_protons->Fill(pt, eta, vz);
                            hPtProtonX->Fill(pt);
                            hEtaProtonX->Fill(eta);
                            MultiplicityProtonX++;
                            if(isEast) {
                                h3D_protonsEast->Fill(pt, eta, vz);
                                hPtProtonEastX->Fill(pt);
                                hEtaProtonEastX->Fill(eta);
                                MultiplicityProtonEastX++;

                                if (b == 0) hEta_xi0_proton_E->Fill(eta);
                                if (b == 1) hEta_xi1_proton_E->Fill(eta);
                                if (b == 2) hEta_xi2_proton_E->Fill(eta);
                                if (b == 3) hEta_xi3_proton_E->Fill(eta);
                                if (b == 4) hEta_xi4_proton_E->Fill(eta);

                                if (b == 0) hEta_xi0_proton_C->Fill(-eta);
                                if (b == 1) hEta_xi1_proton_C->Fill(-eta);
                                if (b == 2) hEta_xi2_proton_C->Fill(-eta);
                                if (b == 3) hEta_xi3_proton_C->Fill(-eta);
                                if (b == 4) hEta_xi4_proton_C->Fill(-eta);

                            }
                            if(isWest) {
                                h3D_protonsWest->Fill(pt, eta, vz);
                                hPtProtonWestX->Fill(pt);
                                hEtaProtonWestX->Fill(eta);
                                MultiplicityProtonWestX++;

                                if (b == 0) hEta_xi0_proton_W->Fill(eta);
                                if (b == 1) hEta_xi1_proton_W->Fill(eta);
                                if (b == 2) hEta_xi2_proton_W->Fill(eta);
                                if (b == 3) hEta_xi3_proton_W->Fill(eta);
                                if (b == 4) hEta_xi4_proton_W->Fill(eta);

                                if (b == 0) hEta_xi0_proton_C->Fill(eta);
                                if (b == 1) hEta_xi1_proton_C->Fill(eta);
                                if (b == 2) hEta_xi2_proton_C->Fill(eta);
                                if (b == 3) hEta_xi3_proton_C->Fill(eta);
                                if (b == 4) hEta_xi4_proton_C->Fill(eta);
                            }

                        } 
                        if(upcEvt->getTrack(j)->getCharge()<0){
                            h3D_antiprotons->Fill(pt, eta, vz);
                            hPtNProtonX->Fill(pt);
                            hEtaNProtonX->Fill(eta);
                            MultiplicityNProtonX++;
                            if(isEast) {
                                h3D_antiprotonsEast->Fill(pt, eta, vz);
                                hPtNProtonEastX->Fill(pt);
                                hEtaNProtonEastX->Fill(eta);
                                MultiplicityNProtonEastX++;

                                if (b == 0) hEta_xi0_antiproton_E->Fill(eta);
                                if (b == 1) hEta_xi1_antiproton_E->Fill(eta);
                                if (b == 2) hEta_xi2_antiproton_E->Fill(eta);
                                if (b == 3) hEta_xi3_antiproton_E->Fill(eta);
                                if (b == 4) hEta_xi4_antiproton_E->Fill(eta);

                                if (b == 0) hEta_xi0_antiproton_C->Fill(-eta);
                                if (b == 1) hEta_xi1_antiproton_C->Fill(-eta);
                                if (b == 2) hEta_xi2_antiproton_C->Fill(-eta);
                                if (b == 3) hEta_xi3_antiproton_C->Fill(-eta);
                                if (b == 4) hEta_xi4_antiproton_C->Fill(-eta);
                            }
                            if(isWest) {
                                h3D_antiprotonsWest->Fill(pt, eta, vz);
                                hPtNProtonWestX->Fill(pt);
                                hEtaNProtonWestX->Fill(eta);
                                MultiplicityNProtonWestX++;

                                if (b == 0) hEta_xi0_antiproton_W->Fill(eta);
                                if (b == 1) hEta_xi1_antiproton_W->Fill(eta);
                                if (b == 2) hEta_xi2_antiproton_W->Fill(eta);
                                if (b == 3) hEta_xi3_antiproton_W->Fill(eta);
                                if (b == 4) hEta_xi4_antiproton_W->Fill(eta);

                                if (b == 0) hEta_xi0_antiproton_C->Fill(eta);
                                if (b == 1) hEta_xi1_antiproton_C->Fill(eta);
                                if (b == 2) hEta_xi2_antiproton_C->Fill(eta);
                                if (b == 3) hEta_xi3_antiproton_C->Fill(eta);
                                if (b == 4) hEta_xi4_antiproton_C->Fill(eta);
                            }
                        }
                    }
                }

                if (abs(upcEvt->getTrack(j)->getEta())<0.9) {
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
                    double pq = (upcEvt->getTrack(j)->getPt())*(upcEvt->getTrack(j)->getCharge());
                    hdEdx->Fill(pq, upcEvt->getTrack(j)->getDEdxSignal()*1e6); //GeV originally, keV now
                } 
                if (upcEvt->getTrack(j)->getPt()>0.2) {
                    if (isEast) HistEtaTracksEast->Fill(upcEvt->getTrack(j)->getEta());
                    if (isWest) HistEtaTracksWest->Fill(upcEvt->getTrack(j)->getEta()); 
                }
            }
        }
       
        HistNumOfPrimaryTracksToF->Fill(NumOfPrimaryTracksToF);
        hMultiplicityProtonX->Fill(MultiplicityProtonX);
        hMultiplicityNProtonX->Fill(MultiplicityNProtonX);

        if (isEast) {
            HistNumOfPToFEast->Fill(NumOfPrimaryTracksToF);
            hMultiplicityProtonEastX->Fill(MultiplicityProtonEastX);
            hMultiplicityNProtonEastX->Fill(MultiplicityNProtonEastX);
        }
        if (isWest) {
            HistNumOfPToFWest->Fill(NumOfPrimaryTracksToF);
            hMultiplicityProtonWestX->Fill(MultiplicityProtonWestX);
            hMultiplicityNProtonWestX->Fill(MultiplicityNProtonWestX);
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