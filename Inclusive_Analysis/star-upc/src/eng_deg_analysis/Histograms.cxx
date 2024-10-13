// Run by: ./Ana file.list output/
// e.g. ./Ana /star/u/truhlar/star-upcDst/build/run17.list ./output/


// Table of RP indecies and names
// RP_ID   0,    1,    2,   3,   4,   5,   6, 7
// RP_name E1U, E1D, E2U, E2D, W1U, W1D, W2U, W2D

// c++ headers
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

// ROOT headers
#include "TROOT.h"
#include "TSystem.h"
#include "TThread.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include <TH2.h> 
#include <TF1.h> 
#include <TF2.h> 
#include <THStack.h> 
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

// picoDst headers
#include "StRPEvent.h"
#include "StUPCRpsTrack.h"
#include "StUPCRpsTrackPoint.h"
#include "StUPCEvent.h"
#include "StUPCTrack.h"
#include "StUPCBemcCluster.h"
#include "StUPCVertex.h"
#include "StUPCTofHit.h"

using namespace std;

// enums are very usefull 
enum { kAll = 1, kCPT,  kRP, kOneVertex, kTPCTOF, 
    kTotQ, kMax};
enum SIDE {E = 0, East = 0, W = 1, West = 1, nSides};
enum PARTICLES {Pion = 0, Kaon = 1, Proton = 2, nParticles};
enum BRANCH_ID { EU, ED, WU, WD, nBranches };
enum RP_ID {E1U, E1D, E2U, E2D, W1U, W1D, W2U, W2D, nRomanPots};

const double beamMomentum = 255;
const double c = 299792458.0;
const double particleMass[nParticles] = { 0.13957, 0.493677, 0.93827 }; // pion, kaon, proton in GeV /c^2 
const int nTriggers = 17;
const int triggerID[] = { 570209, 570219, 570229, 570701, 570702, 570703, 570704, 570705, 
                  570709, 570711, 570712, 570719, 590701, 590703, 590705, 590708, 590709};
// 570702 RP_UPC // 570712 RP_UPC // 570703 RP_SDT // 570709 RP_ET // 570719 RP_ET // 570701 RP_CPT2 // 570711 RP_CPT2 // 570705 RP_CPT2noBBCL // 570704 RP_Zerobias // 590703 RP_SDT // 590709 RP_ET // 590701 RP_CPT2 // 590705 RP_CPT2noBBCL // 590708 RP_CPTnoBBCL // 570209 JPsi*HTTP // 570219 JPsi*HTTP // 570229 JPsi*HTTP
const int CEPtriggers[] = { 570701, 570705,       570711, 590701, 590705,       590708};
//                          RP_CPT2 RP_CPT2noBBCL RP_CPT2 RP_CPT2 RP_CPT2noBBCL RP_CPTnoBBCL

bool ConnectInput(int argc, char** argv, TChain* fileChain);
bool CheckTriggers(StUPCEvent* localupcEvt);
long long GetFileSize(string filename);

//_____________________________________________________________________________
int main(int argc, char** argv) 
{
    int nthreads = 2;
    if(argc==4){
        nthreads = atoi(argv[3]);
    }
    cout<<"Program is running on "<<nthreads<<" threads"<<endl;
    ROOT::EnableThreadSafety();
    ROOT::EnableImplicitMT(nthreads); //turn on multicore processing
    ROOT::EnableThreadSafety();
    
    TChain* upcChain = new TChain("mUPCTree");    //chain with files to iterate through

    const string& outputFolder = argv[2];

    if(!ConnectInput(argc, argv, upcChain))
    {
        cout << "Wrong input parameters..." << endl; 
        return 1;
    }

    // Define the function that will process a subrange of the tree.
    // The function must receive only one parameter, a TTreeReader,
    // and it must be thread safe. To enforce the latter requirement,
    // TThreadedObject histograms will be used.

    //2 photomultipliers in each roman pot
    ROOT::TThreadedObject<TH2I> pmBefore("pmBeforeHist", "Photomultiplier coactivation before ADC cuts", nRomanPots*2, 0, nRomanPots*2, nRomanPots*2, 0, nRomanPots*2);
    ROOT::TThreadedObject<TH2I> pmAfter("pmAfterHist", "Photomultiplier coactivation after ADC cuts", nRomanPots*2, 0, nRomanPots*2, nRomanPots*2, 0, nRomanPots*2);
    TString brName[nBranches] = {TString("EU"), TString("ED"), TString("WU"), TString("WD")};
    TString rpName[nRomanPots] = {TString("E1U"), TString("E1D"), TString("E2U"), TString("E2D"), TString("W1U"), TString("W1D"), TString("W2U"), TString("W2D")};
    TString pmName[nRomanPots*2] = {TString("E1U1"), TString("E1U2"), TString("E1D1"), TString("E1D2"), TString("E2U1"), TString("E2U2"), TString("E2D1"), TString("E2D2"), TString("W1U1"), TString("W1U2"), TString("W1D1"), TString("W1D2"), TString("W2U1"), TString("W2U2"), TString("W2D1"), TString("W2D2")};
    for(int i=0; i<nRomanPots*2; ++i){
        //because 0 is illegal bin number
        pmBefore->GetXaxis()->SetBinLabel(i+1, pmName[i]);
        pmBefore->GetYaxis()->SetBinLabel(i+1, pmName[i]);
        pmAfter->GetXaxis()->SetBinLabel(i+1, pmName[i]);
        pmAfter->GetYaxis()->SetBinLabel(i+1, pmName[i]);
    }
    //histograms for measuring plane efficiency, there are 4 planes per roman pot
    ROOT::TThreadedObject<TH1D> planes3("planes3Hist", "planes3Hist", nRomanPots*4, 0, nRomanPots*4);
    ROOT::TThreadedObject<TH1D> planes4("planes4Hist", "planes4Hist", nRomanPots*4, 0, nRomanPots*4);
    //histograms for adc & tac
    ROOT::TThreadedObject<TH2I> adch("adch", "ADC values;PM Id", nRomanPots*2, 0, nRomanPots*2, 3200, 0, 3200);       //max value was around that
    ROOT::TThreadedObject<TH2I> tach("tach", "TAC values;PM Id", nRomanPots*2, 0, nRomanPots*2, 2400, 0, 2400);       //max value was around that
    for(int i=0; i<nRomanPots*2; ++i){
        //because 0 is illegal bin number
        adch->GetXaxis()->SetBinLabel(i+1, pmName[i]);
        tach->GetXaxis()->SetBinLabel(i+1, pmName[i]);
    }

    //histograms, with wide range, RP
    ROOT::TThreadedObject<TH1D> phiHist("phiHist", "#phi angle of reconstructed protons", 1000, -3.1415926535, 3.1415926535);
    ROOT::TThreadedObject<TH1D> xiHist("xiHist", "#xi variable of reconstructed protons", 160000, -50, 1);
    ROOT::TThreadedObject<TH1D> tHist("tHist", "Four-momentum transfer", 240000, -23900, 100);
    ROOT::TThreadedObject<TH1D> EHist("EHist", "Energy of reconstructed protons", 240000, 0, 480000);
    ROOT::TThreadedObject<TH1D> ptHist("ptHist", "Transverse momentum of reconstructed protons", 192000, 0, 9600);
    ROOT::TThreadedObject<TH1D> etaHist("etaHist", "#eta variable of reconstructed protons", 2000, -20, 20);
    ROOT::TThreadedObject<TH2D> tphiHist("tphiHist", "t(#phi) dependency of reconstructed protons", 1000, -3.1415926535, 3.1415926535, 25000, -4900, 100);
    ROOT::TThreadedObject<TH2D> thetaphiHist("thetaphiHist", "#theta(#phi) dependency of reconstructed protons", 1000, -3.1415926535, 3.1415926535, 1000, 0, 0.05);
    ROOT::TThreadedObject<TH2D> etaphiHist("etaphiHist", "#eta(#phi) dependency of reconstructed protons", 1000, -3.1415926535, 3.1415926535, 2000, -20, 20);
    ROOT::TThreadedObject<TH2D> EdphiHist("EdphiHist", "Energy of reconstructed protons vs polar angle between them", 1000, -3.1415926535*2, 3.1415926535*2, 4000, 0, 2000);
    ROOT::TThreadedObject<TH2D> NdphiHist("NdphiHist", "Number of created particles vs polar angle between protons", 1000, -3.1415926535*2, 3.1415926535*2, 15, 0, 15);
    ROOT::TThreadedObject<TH1D> EprecHist("EprecHist", "Energy of reconstructed protons (around beam energy)", 300000, 240, 270);
    ROOT::TThreadedObject<TH1D> EeltrigHist("EeltrigHist", "Elastic trigger activation", 240000, 0, 480000);
    ROOT::TThreadedObject<TH2D> pxpyHist("pxpyHist", "Transverse momentum of reconstructed protons, by axis", 2000, -100, 100, 2000, -100, 100);
    ROOT::TThreadedObject<TH2D> pxpyprecHist("pxpyprecHist", "Transverse momentum of reconstructed protons, by axis around (0,0)", 400, -2, 2, 400, -2, 2);
    ROOT::TThreadedObject<TH2D> pxpyfidHist("pxpyfidHist", "Sum of transverse momenta of reconstructed protons, by axis", 2000, -100, 100, 2000, -100, 100);
    ROOT::TThreadedObject<TH2D> pxpyfidprecHist("pxpyfidprecHist", "Sum of transverse momenta of reconstructed protons, by axis around (0,0)", 400, -2, 2, 400, -2, 2);
    ROOT::TThreadedObject<TH2D> pxpyfidmoreprecHist("pxpyfidmoreprecHist", "Sum of transverse momenta of reconstructed protons, by axis around (0,0), very precise", 600, -0.3, 0.3, 600, -0.3, 0.3);
    ROOT::TThreadedObject<TH2D> pxpylackHist("pxpylackHist", "Sum of transverse momenta of all produced particles, by axis", 2000, -100, 100, 2000, -100, 100);
    ROOT::TThreadedObject<TH2D> pxpylackprecHist("pxpylackprecHist", "Sum of transverse momenta of all produced particles, by axis around (0,0)", 400, -2, 2, 400, -2, 2);
    ROOT::TThreadedObject<TH2D> pxpylackmoreprecHist("pxpylackmoreprecHist", "Sum of transverse momenta of all produced particles, by axis around (0,0), very precise", 600, -0.3, 0.3, 600, -0.3, 0.3);

    //histograms, with wide range, UPC
    ROOT::TThreadedObject<TH1D> MHist("MHist", "m_{inv} of reconstructed kaon candidates", 250000, 0, 50);
    ROOT::TThreadedObject<TH1D> MbcgHist("MbcgHist", "Background m_{inv} from pairs of particles", 250000, 0, 50);
    ROOT::TThreadedObject<TH1D> MforThesisHist("MforThesisHist", "m_{inv} of reconstructed kaon candidates;m [GeV/c^{2}]", 3000, 0, 3);
    ROOT::TThreadedObject<TH1D> MbcgforThesisHist("MbcgforThesisHist", "Background m_{inv};m [GeV/c^{2}", 3000, 0, 3);
    ROOT::TThreadedObject<TH1D> MforMatchHist("MforMatchHist", "m_{inv} of reconstructed kaon candidates;m [GeV/c^{2}]", 300, 0.46, 0.52);
    ROOT::TThreadedObject<TH1D> MbcgforMatchHist("MbcgforMatchHist", "Background m_{inv};m [GeV/c^{2}", 300, 0.46, 0.52);
    ROOT::TThreadedObject<TH1D> MrawHist("MrawHist", "m_{inv} of detected particles", 250000, 0, 2500);
    ROOT::TThreadedObject<TH1D> NparHist("NparHist", "Number of detected particles", 3, 0, 3);
    NparHist->GetXaxis()->SetBinLabel(1, "positive");
    NparHist->GetXaxis()->SetBinLabel(2, "all");
    NparHist->GetXaxis()->SetBinLabel(3, "negative");
    ROOT::TThreadedObject<TH2D> ProtonidmassHist("ProtonidmassHist", "Proton #frac{dE}{dx} identification vs TOF mass", 500, 0, 50, 700, -50, 20);
    ROOT::TThreadedObject<TH2D> KaonidmassHist("KaonidmassHist", "Kaon #frac{dE}{dx} identification vs TOF mass", 500, 0, 50, 600, -30, 30);
    ROOT::TThreadedObject<TH2D> PionidmassHist("PionidmassHist", "Pion #frac{dE}{dx} identification vs TOF mass", 500, 0, 50, 400, -5, 35);
    ROOT::TThreadedObject<TH2D> ElectronidmassHist("ElectronidmassHist", "Electron #frac{dE}{dx} identification vs TOF mass", 500, 0, 50, 400, -10, 30);
    ROOT::TThreadedObject<TH2D> dEdxidmassHist("dEdxidmassHist", "#frac{dE}{dx} identification vs TOF mass", 5000, 0, 5, 4, 0, 4);
    dEdxidmassHist->GetYaxis()->SetBinLabel(1, "Proton");
    dEdxidmassHist->GetYaxis()->SetBinLabel(2, "Kaon");
    dEdxidmassHist->GetYaxis()->SetBinLabel(3, "Pion");
    dEdxidmassHist->GetYaxis()->SetBinLabel(4, "Electron");
    ROOT::TThreadedObject<TH2D> pxpxHist("pxpxHist", "X-axis momentum of reconstructed protons vs particles from central production;", 2000, -1000, 1000, 2000, -10, 10);
    ROOT::TThreadedObject<TH2D> pypyHist("pypyHist", "Y-axis momentum of reconstructed protons vs particles from central production", 2000, -1000, 1000, 2000, -10, 10);
    ROOT::TThreadedObject<TH1D> TOFtimeHist("TOFtimeHist", "Time given by TOF detector", 55000, 0, 55000);
    ROOT::TThreadedObject<TH1D> TOFpathHist("TOFpathHist", "Path length given by TOF detector", 4000, 200, 600);
    ROOT::TThreadedObject<TH1D> TOFvelHist("TOFvelHist", "Velocity given by TOF detector", 30000, 0, 3e8);
    ROOT::TThreadedObject<TH2D> m2TOFpqHist("m2TOFpqHist", "m_{TOF}^{2} vs pq;pq;m^{2}", 1000, -2.5, 2.5, 1000, 0, 2500);
    ROOT::TThreadedObject<TH2D> m2TOFpqprecHist("m2TOFpqprecHist", "m_{TOF}^{2} vs pq for small pq values;pq;m^{2}", 1000, -2.5, 2.5, 1000, 0, 4);

    ROOT::TThreadedObject<TH2D> m2TOFpathHist("m2TOFpathHist", "m_{TOF}^{2} vs TOF path lenght;L;m^{2}", 1000, 200, 600, 1000, 0, 50);
    ROOT::TThreadedObject<TH2D> m2TOFtimeHist("m2TOFtimeHist", "m_{TOF}^{2} vs TOF flight time;t;m^{2}", 550, 0, 55000, 1000, 0, 50);
    ROOT::TThreadedObject<TH2D> TOFpathtimeHist("TOFpathtimeHist", "TOF path lenght vs TOF flight time;t;L", 550, 0, 55000, 1000, 200, 600);
    ROOT::TThreadedObject<TH2D> pTOFtimeHist("pTOFtimeHist", "momentum vs TOF flight time;t;p", 550, 0, 55000, 1000, 0, 1.5);
    ROOT::TThreadedObject<TH2D> pTOFpathHist("pTOFpathHist", "momentum vs TOF path lenght;L;p", 1000, 200, 600, 1000, 0, 1.5);
    
    ROOT::TThreadedObject<TH2D> dEdxpqHist("dEdxpqHist", "dE/dx vs pq;pq;dE/dx", 1000, -10, 10, 1000, 0, 4e-5);
    ROOT::TThreadedObject<TH1I> NverHist("NverHist", "Number of vertices", 8, 0, 8);
    ROOT::TThreadedObject<TH1I> NprimverHist("NprimverHist", "Number of primary vertices", 8, 0, 8);
    ROOT::TThreadedObject<TH1D> NTOFHist("NTOFHist", "Number of detected particles with and without TOF tag", 2, 0, 2);
    NTOFHist->GetXaxis()->SetBinLabel(1, "TOF unknown");
    NTOFHist->GetXaxis()->SetBinLabel(2, "TOF tagged");


    //standalone histograms

    //histograms for each branch & filling loops
    ROOT::TThreadedObject<TH2D> tphiHisttab[nBranches];
    ROOT::TThreadedObject<TH2D> thetaphiHisttab[nBranches];
    ROOT::TThreadedObject<TH2D> etaphiHisttab[nBranches];
    for (int i = 0; i < nBranches; i++){
        string name = "tphiHist" + string(brName[i].Data());    string title = "Four-momentum transfer (#phi angle dependant) for branch " + string(brName[i].Data());
        //https://stackoverflow.com/questions/2520843/operator-new-for-array-of-class-without-default-constructor
        //no idea why but it works
        new(&tphiHisttab[i])ROOT::TThreadedObject<TH2D>(name.c_str(), title.c_str(), 1000, -3.1415926535, 3.1415926535, 25000, -4900, 100);
        name = "thetaphiHist" + string(brName[i].Data());    title = "#theta(#phi) dependency of reconstructed protons for branch " + string(brName[i].Data());
        new(&thetaphiHisttab[i])ROOT::TThreadedObject<TH2D>(name.c_str(), title.c_str(), 1000, -3.1415926535, 3.1415926535, 1000, 0, 0.05);
        name = "etaphiHist" + string(brName[i].Data());    title = "#eta(#phi) dependency of reconstructed protons for branch " + string(brName[i].Data());
        new(&etaphiHisttab[i])ROOT::TThreadedObject<TH2D>(name.c_str(), title.c_str(), 1000, -3.1415926535, 3.1415926535, 1000, -15, 15);
    }

    //histograms for each RP & filling loops
    ROOT::TThreadedObject<TH2D> yxHisttab[nRomanPots];
    for (int i = 0; i < nRomanPots; i++){
        string name = "yxHist" + string(rpName[i].Data());  string title = "Intensity of detected hits for Roman Pot " + string(rpName[i].Data());
        new(&yxHisttab[i])ROOT::TThreadedObject<TH2D>(name.c_str(), title.c_str(), 1200, -0.06, 0.06, 2000, -0.1, 0.1);
    }

    auto myFunction = [&](TFile* myFile) {
        //test if tree is not empty
        TFile* tempFile = new TFile(myFile->GetTitle());
        TTree* tempTree = (TTree*)tempFile->Get("mUPCTree");
        if(tempTree->GetEntries()==0){
            delete tempFile;
            return 0;
        }
        //creating a reader and all stuff
        TTreeReader myReader(tempTree);
        TTreeReaderValue<StRPEvent> StRPEventInstance(myReader, "mRPEvent");
        TTreeReaderValue<StUPCEvent> StUPCEventInstance(myReader, "mUPCEvent");
        //variable initialization
        myReader.Next();
        //for performance, RP
        auto pmBeforeLocal = pmBefore.Get();
        auto pmAfterLocal = pmAfter.Get();
        auto planes3Local = planes3.Get();
        auto planes4Local = planes4.Get();
        auto adchLocal = adch.Get();
        auto tachLocal = tach.Get();
        auto phiHistLocal = phiHist.Get();
        auto xiHistLocal = xiHist.Get();
        auto tHistLocal = tHist.Get();
        auto EHistLocal = EHist.Get();
        auto ptHistLocal = ptHist.Get();
        auto etaHistLocal = etaHist.Get();
        auto tphiHistLocal = tphiHist.Get();
        auto thetaphiHistLocal = thetaphiHist.Get();
        auto etaphiHistLocal = etaphiHist.Get();
        auto EdphiHistLocal = EdphiHist.Get();
        auto NdphiHistLocal = NdphiHist.Get();
        auto EprecHistLocal = EprecHist.Get();
        auto EeltrigHistLocal = EeltrigHist.Get();
        auto pxpyHistLocal = pxpyHist.Get();
        auto pxpyprecHistLocal = pxpyprecHist.Get();
        auto pxpyfidHistLocal = pxpyfidHist.Get();
        auto pxpyfidprecHistLocal = pxpyfidprecHist.Get();
        auto pxpyfidmoreprecHistLocal = pxpyfidmoreprecHist.Get();
        auto pxpylackHistLocal = pxpylackHist.Get();
        auto pxpylackprecHistLocal = pxpylackprecHist.Get();
        auto pxpylackmoreprecHistLocal = pxpylackmoreprecHist.Get();
        //UPC
        auto MHistLocal = MHist.Get();
        auto MbcgHistLocal = MbcgHist.Get();
        auto MforThesisHistLocal = MforThesisHist.Get();
        auto MbcgforThesisHistLocal = MbcgforThesisHist.Get();
        auto MforMatchHistLocal = MforMatchHist.Get();
        auto MbcgforMatchHistLocal = MbcgforMatchHist.Get();
        auto MrawHistLocal = MrawHist.Get();
        auto NparHistLocal = NparHist.Get();
        auto ProtonidmassHistLocal = ProtonidmassHist.Get();
        auto KaonidmassHistLocal = KaonidmassHist.Get();
        auto PionidmassHistLocal = PionidmassHist.Get();
        auto ElectronidmassHistLocal = ElectronidmassHist.Get();
        auto dEdxidmassHistLocal = dEdxidmassHist.Get();
        auto pxpxHistLocal = pxpxHist.Get();
        auto pypyHistLocal = pypyHist.Get();
        auto TOFtimeHistLocal = TOFtimeHist.Get();
        auto TOFpathHistLocal = TOFpathHist.Get();
        auto TOFvelHistLocal = TOFvelHist.Get();
        auto m2TOFpqHistLocal = m2TOFpqHist.Get();
        auto m2TOFpqprecHistLocal = m2TOFpqprecHist.Get();

        auto m2TOFpathHistLocal = m2TOFpathHist.Get();
        auto m2TOFtimeHistLocal = m2TOFtimeHist.Get();
        auto TOFpathtimeHistLocal = TOFpathtimeHist.Get();
        auto pTOFtimeHistLocal = pTOFtimeHist.Get();
        auto pTOFpathHistLocal = pTOFpathHist.Get();

        auto dEdxpqHistLocal = dEdxpqHist.Get();
        auto NverHistLocal = NverHist.Get();
        auto NprimverHistLocal = NprimverHist.Get();
        auto NTOFHistLocal = NTOFHist.Get();
        //for performance, tables
        //branches
        std::shared_ptr<TH2D> tphiHisttabLocal[nBranches];
        std::shared_ptr<TH2D> thetaphiHisttabLocal[nBranches];
        std::shared_ptr<TH2D> etaphiHisttabLocal[nBranches];
        for (int i = 0; i < nBranches; i++)
        {
            tphiHisttabLocal[i] = tphiHisttab[i].Get();
            thetaphiHisttabLocal[i] = thetaphiHisttab[i].Get();
            etaphiHisttabLocal[i] = etaphiHisttab[i].Get();
        }
        //RPs
        std::shared_ptr<TH2D> yxHisttabLocal[nRomanPots];
        for (int i = 0; i < nRomanPots; i++)
        {
            yxHisttabLocal[i] = yxHisttab[i].Get();
        }

        //local variables
        vector<int> pmNumbersBeforeLocal;
        vector<int> pmNumbersAfterLocal;
        vector<TLorentzVector> posPart;
        vector<TLorentzVector> negPart;
        TLorentzVector trackVector;
        TClonesArray* trackPointsLocal;
        TClonesArray* tracksLocal;
        StUPCRpsTrackPoint* trackPointInstanceLocal;
        StUPCRpsTrack* trackInstanceLocal;
        StUPCTrack* trackInstanceLocalUPC;
        bool isWorthy;
        do{
            //checking pm activation
            isWorthy = false;
            for (int i = 0; i < nRomanPots*2; i++)
            {
                adchLocal->Fill(i, StRPEventInstance->adc(i/2, i%2));
                tachLocal->Fill(i, StRPEventInstance->tac(i/2, i%2));
                //i/2 & i%2 - clever trick to number pms appropriately, deleted tac:  || StRPEventInstance->tac(i/2,i%2)!=0
                if(StRPEventInstance->adc(i/2,i%2)!=0){
                    pmNumbersBeforeLocal.push_back(i);
                }
                //checking if event activated anything and is worthy of further analysis
                if(StRPEventInstance->adc(i/2, i%2)>=0 && StRPEventInstance->adc(i/2, i%2)<65536){
                    isWorthy = true;
                }
            }

            //filtering out inapropriate adc values
            if(!isWorthy){
                continue;
            }

            //checking pm activation again
            for (int i = 0; i < nRomanPots*2; i++)
            {
                //i/2 & i%2 - clever trick to number pms appropriately
                if(StRPEventInstance->adc(i/2,i%2)!=0 || StRPEventInstance->tac(i/2,i%2)!=0){
                    pmNumbersAfterLocal.push_back(i);
                }
            }

            //checking planes efficiency & other things with single points
            trackPointsLocal = StRPEventInstance->getTrackPoints();
            for (int i = 0; i < trackPointsLocal->GetEntries(); i++){
                trackPointInstanceLocal = (StUPCRpsTrackPoint*)trackPointsLocal->ConstructedAt(i);

                //filling things RP-specific
                yxHisttabLocal[trackPointInstanceLocal->rpId()]->Fill(trackPointInstanceLocal->x(), trackPointInstanceLocal->y());

                //planes efficiency
                if(trackPointInstanceLocal->planesUsed()<3){
                    continue;
                }
                //goes through all planes and adds each detection to respective plane
                for (int j = 0; j < 4; j++){
                    //if there's 3 planes it fills the one where there's lack of signal, and not the other ones
                    //other ones don't count in efficiency of other planes
                    if(trackPointInstanceLocal->clusterId(j)<0 && trackPointInstanceLocal->planesUsed()==3){
                        planes4Local->Fill(trackPointInstanceLocal->rpId()*4 + j, 1.0);
                    }
                    //fills all 4 planes if particle went through 4 of them
                    if(trackPointInstanceLocal->planesUsed()==4){
                        planes4Local->Fill(trackPointInstanceLocal->rpId()*4 + j, 1.0);
                        planes3Local->Fill(trackPointInstanceLocal->rpId()*4 + j, 1.0);
                    }
                }
            }

            //filling all sorts of histograms, and also the one with clusters
            //UPC
            Double_t vel;
            TLorentzVector UPCsum = {0,0,0,0};
            for (int i = 0; i < StUPCEventInstance->getNumberOfTracks(); i++){
                trackInstanceLocalUPC = StUPCEventInstance->getTrack(i);
                //checking if the track is a good quality one (a pion, from primary vertex, through TOF etc)
                //TODO id part
                if(!trackInstanceLocalUPC->getFlag(StUPCTrack::kPrimary) || !trackInstanceLocalUPC->getFlag(StUPCTrack::kTof))
                    continue;
                if(trackInstanceLocalUPC->getTofPathLength()<0 || trackInstanceLocalUPC->getTofTime()<0)
                    continue;
                //transformation from cm/10ps to m/s, it's a factor of 1e9
                vel = trackInstanceLocalUPC->getTofPathLength()/trackInstanceLocalUPC->getTofTime()*1e9;
                //the only difference will be a very simple cut, so it's filled before to check how many are superluminal
                TOFvelHistLocal->Fill(vel);
                double nsigmaTab[4];
                nsigmaTab[0] = abs(trackInstanceLocalUPC->getNSigmasTPCProton());
                nsigmaTab[1] = abs(trackInstanceLocalUPC->getNSigmasTPCKaon());
                nsigmaTab[2] = abs(trackInstanceLocalUPC->getNSigmasTPCPion());
                nsigmaTab[3] = abs(trackInstanceLocalUPC->getNSigmasTPCElectron());
                //if checks positive gets written into arrays
                trackInstanceLocalUPC->getLorentzVector(trackVector, particleMass[Pion]);
                //particle recognition (only pions allowed)
                if(abs(trackInstanceLocalUPC->getNSigmasTPCPion())>3){
                    continue;
                }
                //changed because without t0 TOF mass doesn't work
                if(vel<c){
                    NTOFHistLocal->Fill(1);
                    //divided by v and gamma in place
                    //we use GeV/c and GeV, so we have to multiply by c
                    Double_t tempmass = trackVector.P()/vel*c*sqrt(1-(vel/c)*(vel/c));
                    MrawHistLocal->Fill(tempmass);
                    ProtonidmassHistLocal->Fill(tempmass, trackInstanceLocalUPC->getNSigmasTPCProton());
                    KaonidmassHistLocal->Fill(tempmass, trackInstanceLocalUPC->getNSigmasTPCKaon());
                    PionidmassHistLocal->Fill(tempmass, trackInstanceLocalUPC->getNSigmasTPCPion());
                    ElectronidmassHistLocal->Fill(tempmass, trackInstanceLocalUPC->getNSigmasTPCElectron());
                    if(trackVector.P()<1)
                        dEdxidmassHistLocal->Fill(tempmass, distance(nsigmaTab, min_element(nsigmaTab, nsigmaTab + 4)));
                    m2TOFpqHistLocal->Fill(trackVector.P()*trackInstanceLocalUPC->getCharge(), tempmass*tempmass);
                    m2TOFpqprecHistLocal->Fill(trackVector.P()*trackInstanceLocalUPC->getCharge(), tempmass*tempmass);
                    m2TOFpathHistLocal->Fill(trackInstanceLocalUPC->getTofPathLength(), tempmass*tempmass);
                    m2TOFtimeHistLocal->Fill(trackInstanceLocalUPC->getTofTime(), tempmass*tempmass);
                }
                UPCsum += trackVector;
                if(trackInstanceLocalUPC->getCharge()==1){
                    posPart.push_back(trackVector);
                }else if (trackInstanceLocalUPC->getCharge()==-1){
                    negPart.push_back(trackVector);
                }else{
                    cout<<"Neutral Particle or weirdly charged"<<endl;
                }
                TOFtimeHistLocal->Fill(trackInstanceLocalUPC->getTofTime());
                TOFpathHistLocal->Fill(trackInstanceLocalUPC->getTofPathLength());
                TOFpathtimeHistLocal->Fill(trackInstanceLocalUPC->getTofTime(), trackInstanceLocalUPC->getTofPathLength());
                pTOFtimeHistLocal->Fill(trackInstanceLocalUPC->getTofTime(), trackVector.P());
                pTOFpathHistLocal->Fill(trackInstanceLocalUPC->getTofPathLength(), trackVector.P());
                dEdxpqHistLocal->Fill(trackVector.P()*trackInstanceLocalUPC->getCharge(), trackInstanceLocalUPC->getDEdxSignal());
            }
            NTOFHistLocal->Fill(0.0, StUPCEventInstance->getNumberOfTracks());
            NverHistLocal->Fill(StUPCEventInstance->getNumberOfVertices());
            NprimverHistLocal->Fill(StUPCEventInstance->getNPrimVertices());

            //filling all sorts of histograms, and also the one with clusters
            //RP
            TLorentzVector RPsum = {0,0,0,0};
            if (StUPCEventInstance->isTrigger(570709) || StUPCEventInstance->isTrigger(570719) || StUPCEventInstance->isTrigger(590709))
                EeltrigHistLocal->Fill(trackVector.E());
                
            tracksLocal = StRPEventInstance->getTracks();
            for (int i = 0; i < tracksLocal->GetEntries(); i++){
                trackInstanceLocal = (StUPCRpsTrack*)tracksLocal->ConstructedAt(i);
                phiHistLocal->Fill(trackInstanceLocal->phi());
                xiHistLocal->Fill(trackInstanceLocal->xi(beamMomentum));
                tHistLocal->Fill(trackInstanceLocal->t(beamMomentum));
                trackVector.SetVectM(trackInstanceLocal->pVec(), particleMass[Proton]);
                EHistLocal->Fill(trackVector.E());
                ptHistLocal->Fill(trackInstanceLocal->pt());
                etaHistLocal->Fill(trackInstanceLocal->eta());
                tphiHistLocal->Fill(trackInstanceLocal->phi(), trackInstanceLocal->t(beamMomentum));
                thetaphiHistLocal->Fill(trackInstanceLocal->phi(), trackInstanceLocal->theta());
                etaphiHistLocal->Fill(trackInstanceLocal->phi(), trackInstanceLocal->eta());
                EprecHistLocal->Fill(trackVector.E());
                pxpyHistLocal->Fill(trackVector.X(), trackVector.Y());
                pxpyprecHistLocal->Fill(trackVector.X(), trackVector.Y());
                RPsum += trackVector;
                //filling things branch-specific
                tphiHisttabLocal[trackInstanceLocal->branch()]->Fill(trackInstanceLocal->phi(), trackInstanceLocal->t(beamMomentum));
                thetaphiHisttabLocal[trackInstanceLocal->branch()]->Fill(trackInstanceLocal->phi(), trackInstanceLocal->theta());
                etaphiHisttabLocal[trackInstanceLocal->branch()]->Fill(trackInstanceLocal->phi(), trackInstanceLocal->eta());
                //filling things that need two or more tracks
                for (int j = 0; j < i; j++)
                {
                    EdphiHistLocal->Fill(trackInstanceLocal->phi()-((StUPCRpsTrack*)tracksLocal->ConstructedAt(j))->phi(), trackVector.E());
                    NdphiHistLocal->Fill(trackInstanceLocal->phi()-((StUPCRpsTrack*)tracksLocal->ConstructedAt(j))->phi(), negPart.size()+posPart.size());
                }
                //checking if there are particles with any coordinate zero
                if(trackVector.X()==0 || trackVector.Y()==0 || trackVector.Z()==0){
                    cout<<"-----------------------------------------------"<<endl;
                    cout<<"X: "<<trackVector.X()<<endl;
                    cout<<"Y: "<<trackVector.Y()<<endl;
                    cout<<"Z: "<<trackVector.Z()<<endl;
                }
            }

            //momentum comparison
            //be careful - possible empty RPs
            if(tracksLocal->GetEntries()>0){
                pxpxHistLocal->Fill(RPsum.X(), UPCsum.X());
                pypyHistLocal->Fill(RPsum.Y(), UPCsum.Y());
                pxpyfidHistLocal->Fill(RPsum.X(), RPsum.Y());
                pxpyfidprecHistLocal->Fill(RPsum.X(), RPsum.Y());
                pxpyfidmoreprecHistLocal->Fill(RPsum.X(), RPsum.Y());
                pxpylackHistLocal->Fill(RPsum.X()+UPCsum.X(), RPsum.Y()+UPCsum.Y());
                pxpylackprecHistLocal->Fill(RPsum.X()+UPCsum.X(), RPsum.Y()+UPCsum.Y());
                pxpylackmoreprecHistLocal->Fill(RPsum.X()+UPCsum.X(), RPsum.Y()+UPCsum.Y());
            }

            TLorentzVector tempSum = {0, 0, 0, 0};

            //invariant mass
            for (long unsigned int posi = 0; posi < posPart.size(); posi++)
            {
                for (long unsigned int negi = 0; negi < negPart.size(); negi++)
                {
                    tempSum = posPart[posi] + negPart[negi];
                    MHistLocal->Fill(tempSum.M());
                    MforThesisHistLocal->Fill(tempSum.M());
                    MforMatchHistLocal->Fill(tempSum.M());
                }
            }
            //background mass (pairs of the same charge)
            for (long unsigned int posi = 0; posi < posPart.size(); posi++)
            {
                for (long unsigned int posi2 = 0; posi2 < posi; posi2++)
                {
                    tempSum = posPart[posi] + posPart[posi2];
                    MbcgHistLocal->Fill(tempSum.M());
                    MbcgforThesisHistLocal->Fill(tempSum.M());
                    MbcgforMatchHistLocal->Fill(tempSum.M());
                }
            }
            for (long unsigned int negi = 0; negi < negPart.size(); negi++)
            {
                for (long unsigned int negi2 = 0; negi2 < negi; negi2++)
                {
                    tempSum = negPart[negi] + negPart[negi2];
                    MbcgHistLocal->Fill(tempSum.M());
                    MbcgforThesisHistLocal->Fill(tempSum.M());
                    MbcgforMatchHistLocal->Fill(tempSum.M());
                }
            }
            //how many charged particles are there
            NparHistLocal->Fill((Double_t)0, (Double_t)posPart.size());
            NparHistLocal->Fill((Double_t)1, (Double_t)posPart.size() + (Double_t)negPart.size());
            NparHistLocal->Fill((Double_t)2, (Double_t)negPart.size());

            //filling histograms of coincidence
            for (long unsigned int i = 0; i < pmNumbersBeforeLocal.size(); i++){
                for (long unsigned int j = i; j < pmNumbersBeforeLocal.size(); j++){
                    pmBeforeLocal->Fill(pmNumbersBeforeLocal[i], pmNumbersBeforeLocal[j]);
                    //for a nice, symmetric graph
                    pmBeforeLocal->Fill(pmNumbersBeforeLocal[j], pmNumbersBeforeLocal[i]);
                }
            }
            for (long unsigned int i = 0; i < pmNumbersAfterLocal.size(); i++){
                for (long unsigned int j = i; j < pmNumbersAfterLocal.size(); j++){
                    pmAfterLocal->Fill(pmNumbersAfterLocal[i], pmNumbersAfterLocal[j]);
                    //for a nice, symmetric graph
                    pmAfterLocal->Fill(pmNumbersAfterLocal[j], pmNumbersAfterLocal[i]);
                }
            }

            pmNumbersBeforeLocal.clear();
            pmNumbersAfterLocal.clear();
            posPart.clear();
            negPart.clear();
            trackPointsLocal->Clear();
            tracksLocal->Clear();
        }while (myReader.Next());

        cout<<"Finished operation on input file "<<myFile->GetTitle()<<endl;

        delete tempFile;

        return 0;
    };

    auto redFunction = [](const std::vector<int> &mapV)
    {
        return std::accumulate(mapV.begin(), mapV.end(), 0);
    };

    vector<TFile*> listOfFiles;
    //creating the list of TFile*
    TObjArray* tempList = upcChain->GetListOfFiles();
    for (int i = 0; i < tempList->GetEntries(); i++)
    {
        listOfFiles.push_back((TFile*)(tempList->At(i)));
    }

    // Create a TreeProcessor: specify the file and the tree in it
    ROOT::TThreadExecutor TreeProcessor(nthreads);
    // Launch the parallel processing of the tree
    TreeProcessor.MapReduce(myFunction, listOfFiles, redFunction);
    // Use the TThreadedObject::Merge method to merge the thread private tree
    // into the final result
    //changes from shared pointer into normal one
    //RP
    auto pmBeforeFinal = pmBefore.Merge();
    auto pmAfterFinal = pmAfter.Merge();
    auto planes3Final = planes3.Merge();
    auto planes4Final = planes4.Merge();
    auto adchFinal = adch.Merge();
    auto tachFinal = tach.Merge();
    auto phiHistFinal = phiHist.Merge();
    auto xiHistFinal = xiHist.Merge();
    auto tHistFinal = tHist.Merge();
    auto EHistFinal = EHist.Merge();
    auto ptHistFinal = ptHist.Merge();
    auto etaHistFinal = etaHist.Merge();
    auto tphiHistFinal = tphiHist.Merge();
    auto thetaphiHistFinal = thetaphiHist.Merge();
    auto etaphiHistFinal = etaphiHist.Merge();
    auto EdphiHistFinal = EdphiHist.Merge();
    auto NdphiHistFinal = NdphiHist.Merge();
    auto EprecHistFinal = EprecHist.Merge();
    auto EeltrigHistFinal = EeltrigHist.Merge();
    auto pxpyHistFinal = pxpyHist.Merge();
    auto pxpyprecHistFinal = pxpyprecHist.Merge();
    auto pxpyfidHistFinal = pxpyfidHist.Merge();
    auto pxpyfidprecHistFinal = pxpyfidprecHist.Merge();
    auto pxpyfidmoreprecHistFinal = pxpyfidmoreprecHist.Merge();
    auto pxpylackHistFinal = pxpylackHist.Merge();
    auto pxpylackprecHistFinal = pxpylackprecHist.Merge();
    auto pxpylackmoreprecHistFinal = pxpylackmoreprecHist.Merge();
    //UPC
    auto MHistFinal = MHist.Merge();
    auto MforThesisHistFinal = MforThesisHist.Merge();
    auto MbcgforThesisHistFinal = MbcgforThesisHist.Merge();
    auto MforMatchHistFinal = MforMatchHist.Merge();
    auto MbcgforMatchHistFinal = MbcgforMatchHist.Merge();
    auto MbcgHistFinal = MbcgHist.Merge();
    auto MrawHistFinal = MrawHist.Merge();
    auto NparHistFinal = NparHist.Merge();
    auto ProtonidmassHistFinal = ProtonidmassHist.Merge();
    auto KaonidmassHistFinal = KaonidmassHist.Merge();
    auto PionidmassHistFinal = PionidmassHist.Merge();
    auto ElectronidmassHistFinal = ElectronidmassHist.Merge();
    auto dEdxidmassHistFinal = dEdxidmassHist.Merge();
    auto pxpxHistFinal = pxpxHist.Merge();
    auto pypyHistFinal = pypyHist.Merge();
    auto TOFtimeHistFinal = TOFtimeHist.Merge();
    auto TOFpathHistFinal = TOFpathHist.Merge();
    auto TOFvelHistFinal = TOFvelHist.Merge();
    auto m2TOFpqHistFinal = m2TOFpqHist.Merge();
    auto m2TOFpqprecHistFinal = m2TOFpqprecHist.Merge();

    auto m2TOFpathHistFinal = m2TOFpathHist.Merge();
    auto m2TOFtimeHistFinal = m2TOFtimeHist.Merge();
    auto TOFpathtimeHistFinal = TOFpathtimeHist.Merge();
    auto pTOFtimeHistFinal = pTOFtimeHist.Merge();
    auto pTOFpathHistFinal = pTOFpathHist.Merge();

    auto dEdxpqHistFinal = dEdxpqHist.Merge();
    auto NverHistFinal = NverHist.Merge();
    auto NprimverHistFinal = NprimverHist.Merge();
    auto NTOFHistFinal = NTOFHist.Merge();

    //merging tables specific for branches now
    std::shared_ptr<TH2D> tphiHisttabFinal[nBranches];
    std::shared_ptr<TH2D> thetaphiHisttabFinal[nBranches];
    std::shared_ptr<TH2D> etaphiHisttabFinal[nBranches];
    for (int i = 0; i < nBranches; i++)
    {
        tphiHisttabFinal[i] = tphiHisttab[i].Merge();
    }
    for (int i = 0; i < nBranches; i++)
    {
        thetaphiHisttabFinal[i] = thetaphiHisttab[i].Merge();
    }
    for (int i = 0; i < nBranches; i++)
    {
        etaphiHisttabFinal[i] = etaphiHisttab[i].Merge();
    }
    //merging tables specific for RP now
    std::shared_ptr<TH2D> yxHisttabFinal[nRomanPots];
    for (int i = 0; i < nRomanPots; i++)
    {
        yxHisttabFinal[i] = yxHisttab[i].Merge();
    } 

    //making the planes efficiency histogram final
    planes3Final->Divide(planes4Final.get());

    //setting up a tree & output file
    string outfileName = outputFolder + "AnaOutput_Hist.root";
    cout<<"Created output file "<<outfileName<<endl;
    TFile* outputFileHist = TFile::Open(outfileName.c_str(), "recreate");

    //writing to file
    //RP
    outputFileHist->cd();
    pmBeforeFinal->Write();
    pmAfterFinal->Write();
    planes3Final->Write();
    adchFinal->Write();
    tachFinal->Write();
    phiHistFinal->Write();
    xiHistFinal->Write();
    tHistFinal->Write();
    EHistFinal->Write();
    ptHistFinal->Write();
    etaHistFinal->Write();
    tphiHistFinal->Write();
    thetaphiHistFinal->Write();
    etaphiHistFinal->Write();
    EdphiHistFinal->Write();
    NdphiHistFinal->Write();
    EprecHistFinal->Write();
    EeltrigHistFinal->Write();
    pxpyHistFinal->Write();
    pxpyprecHistFinal->Write();
    pxpyfidHistFinal->Write();
    pxpyfidprecHistFinal->Write();
    pxpyfidmoreprecHistFinal->Write();
    pxpylackHistFinal->Write();
    pxpylackprecHistFinal->Write();
    pxpylackmoreprecHistFinal->Write();
    //UPC
    MHistFinal->Write();
    MforThesisHistFinal->Write();
    MbcgforThesisHistFinal->Write();
    MforMatchHistFinal->Write();
    MbcgforMatchHistFinal->Write();
    MbcgHistFinal->Write();
    MrawHistFinal->Write();
    NparHistFinal->Write();
    ProtonidmassHistFinal->Write();
    KaonidmassHistFinal->Write();
    PionidmassHistFinal->Write();
    ElectronidmassHistFinal->Write();
    dEdxidmassHistFinal->Write();
    pxpxHistFinal->Write();
    pypyHistFinal->Write();
    TOFtimeHistFinal->Write();
    TOFpathHistFinal->Write();
    TOFvelHistFinal->Write();
    m2TOFpqHistFinal->Write();
    m2TOFpqprecHistFinal->Write();

    m2TOFpathHistFinal->Write();
    m2TOFtimeHistFinal->Write();
    TOFpathtimeHistFinal->Write();
    pTOFtimeHistFinal->Write();
    pTOFpathHistFinal->Write();

    dEdxpqHistFinal->Write();
    NverHistFinal->Write();
    NprimverHistFinal->Write();
    NTOFHistFinal->Write();
    //writing branch tables to file
    for (int i = 0; i < nBranches; i++)
    {
        tphiHisttabFinal[i]->Write();
    }
    for (int i = 0; i < nBranches; i++)
    {
        thetaphiHisttabFinal[i]->Write();
    }
    for (int i = 0; i < nBranches; i++)
    {
        etaphiHisttabFinal[i]->Write();
    }
    //writing RP tables to file
    for (int i = 0; i < nRomanPots; i++)
    {
        yxHisttabFinal[i]->Write();
    }

    outputFileHist->Close();

    cout<<"Finished processing "<<endl;
    cout<<"Ending Analysis... GOOD BYE!"<<endl;
    return 0;
}//main

bool ConnectInput(int argc, char** argv, TChain* fileChain) 
{
    int fileId = -1;
    string line;
    int lineId=0;

    const string& input = argv[1];
    cout<<"Using list "<<input<<endl;
    if(input.find(".list") != string::npos )
    {
        cout << "Input from chain" << endl;
        ifstream instr(input.c_str());
        if (!instr.is_open())
        {
            cout<< "Couldn't open: "<<input.c_str()<<endl;
            return false;
        }
        //for testing if file exists
        TFile* infile;

        while(getline(instr, line)) 
        {
            if(fileId==lineId || fileId== -1)
            {
                fileChain->AddFile(line.c_str());
                infile = TFile::Open(line.c_str(), "read");
                if(!infile)
                {
                    cout<< "Couldn't open: "<<line.c_str()<<endl;
                    return false;
                }
                infile->Close();
            }
            lineId++;
        }
        instr.close();
    }

    return true;
}//ConnectInput

bool CheckTriggers(StUPCEvent* localupcEvt)
{

    bool CPTtrigger = false;
    for(int var = 0; var < nTriggers; ++var)
    {
        if(localupcEvt->isTrigger(triggerID[var]))
        {
            //Checked if it is CPT trigger
            for (int i = 0; i < *(&CEPtriggers + 1) - CEPtriggers; ++i)
                if(triggerID[var] == CEPtriggers[i])
                    CPTtrigger=true;
        }
    }

    return CPTtrigger;
}

long long GetFileSize(string filename){
    struct stat64 stat_buf;
    int rc = stat64(filename.c_str(), &stat_buf);
    return rc==0 ? stat_buf.st_size : -1;
}