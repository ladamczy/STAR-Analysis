// Run by: ./Ana file.list output.root
// e.g. ./Ana /star/u/truhlar/star-upcDst/build/run17.list ./output.root


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

// my headers
#include "ProcessingInsideLoop.h"
#include "ProcessingOutsideLoop.h"

using namespace std;

// enums are very usefull 
enum{
    kAll = 1, kCPT, kRP, kOneVertex, kTPCTOF,
    kTotQ, kMax
};
enum SIDE{ E = 0, East = 0, W = 1, West = 1, nSides };
enum PARTICLES{ Pion = 0, Kaon = 1, Proton = 2, nParticles };
enum BRANCH_ID{ EU, ED, WU, WD, nBranches };
enum RP_ID{ E1U, E1D, E2U, E2D, W1U, W1D, W2U, W2D, nRomanPots };

const double particleMass[nParticles] = { 0.13957, 0.493677, 0.93827 }; // pion, kaon, proton in GeV /c^2 
const int nTriggers = 17;
const int triggerID[] = { 570209, 570219, 570229, 570701, 570702, 570703, 570704, 570705,
                  570709, 570711, 570712, 570719, 590701, 590703, 590705, 590708, 590709 };
// 570702 RP_UPC // 570712 RP_UPC // 570703 RP_SDT // 570709 RP_ET // 570719 RP_ET // 570701 RP_CPT2 // 570711 RP_CPT2 // 570705 RP_CPT2noBBCL // 570704 RP_Zerobias // 590703 RP_SDT // 590709 RP_ET // 590701 RP_CPT2 // 590705 RP_CPT2noBBCL // 590708 RP_CPTnoBBCL // 570209 JPsi*HTTP // 570219 JPsi*HTTP // 570229 JPsi*HTTP
const int CEPtriggers[] = { 570701, 570705, 570711, 590701, 590705, 590708 };

bool ConnectInput(int argc, char **argv, TChain *fileChain);
bool CheckTriggers(StUPCEvent *localupcEvt);
long long GetFileSize(string filename);

//_____________________________________________________________________________
int main(int argc, char **argv){
    int nthreads = 2;
    if(argc==4){
        nthreads = atoi(argv[3]);
    }
    cout<<"Program is running on "<<nthreads<<" threads"<<endl;
    ROOT::EnableThreadSafety();
    ROOT::EnableImplicitMT(nthreads); //turn on multicore processing
    ROOT::EnableThreadSafety();

    TChain *upcChain = new TChain("mUPCTree");    //chain with files to iterate through

    if(!ConnectInput(argc, argv, upcChain)){
        cout<<"Wrong input parameters..."<<endl;
        return 1;
    }

    const string &outputFolder = argv[2];

    //HISTOGRAMS
    //trigger histograms
    ROOT::TThreadedObject<TH1D> trigHist("trigHist", ";;events", 4, 0, 4);
    trigHist->GetXaxis()->SetBinLabel(1, "RP_CPT2 start");
    trigHist->GetXaxis()->SetBinLabel(2, "RP_CPT2noBBCL start");
    trigHist->GetXaxis()->SetBinLabel(3, "RP_CPT2 end");
    trigHist->GetXaxis()->SetBinLabel(4, "RP_CPT2noBBCL endS");
    //one proton on each side
    ROOT::TThreadedObject<TH1D> oneProtonEUHist("oneProtonEUHist", ";;events", 15, 0, 15);
    ROOT::TThreadedObject<TH1D> oneProtonEDHist("oneProtonEDHist", ";;events", 15, 0, 15);
    ROOT::TThreadedObject<TH1D> oneProtonWUHist("oneProtonWUHist", ";;events", 15, 0, 15);
    ROOT::TThreadedObject<TH1D> oneProtonWDHist("oneProtonWDHist", ";;events", 15, 0, 15);
    //histograms for measuring plane efficiency, there are 4 planes per roman pot and difference between 3 and 4 planes
    ROOT::TThreadedObject<TH1D> planes3("planes3Hist", "planes3Hist", nRomanPots*4, 0, nRomanPots*4);
    ROOT::TThreadedObject<TH1D> planes4("planes4Hist", "planes4Hist", nRomanPots*4, 0, nRomanPots*4);
    ROOT::TThreadedObject<TH2D> pxpy3EHist("pxpy3EHist", ";p_{x} [GeV];p_{y} [GeV]", 400, -2, 2, 400, -2, 2);
    ROOT::TThreadedObject<TH2D> pxpy4EHist("pxpy4EHist", ";p_{x} [GeV];p_{y} [GeV]", 400, -2, 2, 400, -2, 2);
    ROOT::TThreadedObject<TH2D> pxpy3WHist("pxpy3WHist", ";p_{x} [GeV];p_{y} [GeV]", 400, -2, 2, 400, -2, 2);
    ROOT::TThreadedObject<TH2D> pxpy4WHist("pxpy4WHist", ";p_{x} [GeV];p_{y} [GeV]", 400, -2, 2, 400, -2, 2);
    //exactly one primary vertex            number of vertices and prime vertices
    ROOT::TThreadedObject<TH1D> NverHist("NverHist", ";Vertices;events", 8, 0, 8);
    ROOT::TThreadedObject<TH1D> NprimverHist("NprimverHist", ";Primary vertices;events", 8, 0, 8);
    //Vertex middle position, TPC RP
    ROOT::TThreadedObject<TH1D> verZHist("verZHist", ";VP_{zTPC} [cm];events", 120, -240, 240);
    ROOT::TThreadedObject<TH1D> verZRPHist("verZRPHist", ";VP_{zRP} [cm];events", 120, -240, 240);
    //vertex position difference, 1D and 2D
    ROOT::TThreadedObject<TH1D> verDeltaZHist("verDeltaZHist", ";VP_{zTPC}-VP_{zRP} [cm];events", 240, -240, 240);
    ROOT::TThreadedObject<TH2D> verZRPZTPCHist("verZRPZTPCHist", ";VP_{zRP} [cm];VP_{zTPC} [cm]", 240, -240, 240, 240, -240, 240);
    //fiducial cuts and checking where was such hard cut
    ROOT::TThreadedObject<TH2D> pxpyEHist("pxpyEHist", ";p_{x} [GeV];p_{y} [GeV]", 300, -1.5, 1.5, 300, -1.5, 1.5);
    ROOT::TThreadedObject<TH2D> pxpyWHist("pxpyWHist", ";p_{x} [GeV];p_{y} [GeV]", 300, -1.5, 1.5, 300, -1.5, 1.5);
    ROOT::TThreadedObject<TH2D> pxpyEOnePrimaryHist("pxpyEOnePrimaryHist", ";p_{x} [GeV];p_{y} [GeV]", 300, -1.5, 1.5, 300, -1.5, 1.5);
    ROOT::TThreadedObject<TH2D> pxpyWOnePrimaryHist("pxpyWOnePrimaryHist", ";p_{x} [GeV];p_{y} [GeV]", 300, -1.5, 1.5, 300, -1.5, 1.5);
    ROOT::TThreadedObject<TH2D> pxpyEPrimaryMiddleHist("pxpyEPrimaryMiddleHist", ";p_{x} [GeV];p_{y} [GeV]", 300, -1.5, 1.5, 300, -1.5, 1.5);
    ROOT::TThreadedObject<TH2D> pxpyWPrimaryMiddleHist("pxpyWPrimaryMiddleHist", ";p_{x} [GeV];p_{y} [GeV]", 300, -1.5, 1.5, 300, -1.5, 1.5);
    //elastic collisions
    ROOT::TThreadedObject<TH2D> pxpySumHist("pxpySumHist", ";#Sigma p_{x} [GeV];#Sigma p_{y} [GeV]", 400, -2, 2, 400, -2, 2);
    //TOF match         particle number             particle number without TOF checking
    int nEvents = 10;
    ROOT::TThreadedObject<TH2D> NpartHist("NpartHist", ";;Number of produced TOF-matched tracks", nEvents, 0, nEvents, 40, 0, 40);
    NpartHist->GetXaxis()->SetBinLabel(1, "All");
    NpartHist->GetXaxis()->SetBinLabel(2, "CPT cuts");
    NpartHist->GetXaxis()->SetBinLabel(3, "#splitline{one proton}{each side}");
    NpartHist->GetXaxis()->SetBinLabel(4, "#splitline{8 planes}{used in RPs}");
    NpartHist->GetXaxis()->SetBinLabel(5, "one primary vertex");
    NpartHist->GetXaxis()->SetBinLabel(6, "#splitline{primary vertex}{in central position}");
    NpartHist->GetXaxis()->SetBinLabel(7, "#splitline{primary vertex}{unambiguous position}");
    NpartHist->GetXaxis()->SetBinLabel(8, "Fiducial cuts");
    NpartHist->GetXaxis()->SetBinLabel(9, "#splitline{Elastic}{collision cuts}");
    NpartHist->GetXaxis()->SetBinLabel(10, "#splitline{central tracks}{quality check}");
    ROOT::TThreadedObject<TH2D> NallpartHist("NallpartHist", ";;Number of tracks", nEvents, 0, nEvents, 40, 0, 40);
    NallpartHist->GetXaxis()->SetBinLabel(1, "All");
    NallpartHist->GetXaxis()->SetBinLabel(2, "CPT cuts");
    NallpartHist->GetXaxis()->SetBinLabel(3, "#splitline{one proton}{each side}");
    NallpartHist->GetXaxis()->SetBinLabel(4, "#splitline{8 planes}{used in RPs}");
    NallpartHist->GetXaxis()->SetBinLabel(5, "one primary vertex");
    NallpartHist->GetXaxis()->SetBinLabel(6, "#splitline{primary vertex}{in central position}");
    NallpartHist->GetXaxis()->SetBinLabel(7, "#splitline{primary vertex}{unambiguous position}");
    NallpartHist->GetXaxis()->SetBinLabel(8, "Fiducial cuts");
    NallpartHist->GetXaxis()->SetBinLabel(9, "#splitline{Elastic}{collision cuts}");
    NallpartHist->GetXaxis()->SetBinLabel(10, "#splitline{central tracks}{quality check}");
    //proper track reconstruction and energy evaluation
    ROOT::TThreadedObject<TH1D> NfitHitsHist("NfitHitsHist", ";N_{TPCfit}^{hits};Number of tracks", 50, 0, 50);
    ROOT::TThreadedObject<TH1D> NdEdxHitsHist("NdEdxHitsHist", ";N_{dE/dx}^{hits};Number of tracks", 50, 0, 50);
    //DCA from vertex
    ROOT::TThreadedObject<TH1D> DCAHist("DCAHist", ";R [cm];Number of tracks", 200, 0, 4);
    //pseudorapidity & momentum
    ROOT::TThreadedObject<TH1D> ptHist("ptHist", ";p_{T} [GeV];Number of tracks", 200, 0, 2);
    ROOT::TThreadedObject<TH1D> etaHist("etaHist", ";#eta;Number of tracks", 100, -2.5, 2.5);
    //particle identification
    ROOT::TThreadedObject<TH2D> dEdxpqHist("dEdxpqHist", ";pq [GeV];dE/dx [GeV/cm]", 1000, -4, 4, 1000, 0, 4e-5);
    //invariant mass calculation
    //mass histograms
    ROOT::TThreadedObject<TH1D> MforThesisHist("MforThesisHist", ";m_{inv} [GeV];Number of pairs", 300, 0, 3);
    ROOT::TThreadedObject<TH1D> MforMatchHist("MforMatchHist", ";m_{inv} [GeV];Number of pairs", 100, 0.42, 0.56);
    ROOT::TThreadedObject<TH1D> MLambdaforMatchHist("MLambdaforMatchHist", ";m_{inv} [GeV];Number of pairs", 100, 1.06, 1.16);
    ROOT::TThreadedObject<TH1D> MpipiforThesisHist("MpipiforThesisHist", ";m_{#pi^{+}#pi^{-}} [GeV];Number of pairs", 300, 0, 3);
    ROOT::TThreadedObject<TH1D> MpipiforMatchHist("MpipiforMatchHist", ";m_{#pi^{+}#pi^{-}} [GeV];Number of pairs", 100, 0.42, 0.56);
    ROOT::TThreadedObject<TH1D> MpipiLambdaforMatchHist("MpipiLambdaforMatchHist", ";m_{#pi^{+}#pi^{-}} [GeV];Number of pairs", 100, 1.06, 1.16);
    ROOT::TThreadedObject<TH1D> MppiforThesisHist("MppiforThesisHist", ";m_{p^{#pm}#pi^{#mp}} [GeV];Number of pairs", 300, 0, 3);
    ROOT::TThreadedObject<TH1D> MppiforMatchHist("MppiforMatchHist", ";m_{p^{#pm}#pi^{#mp}} [GeV];Number of pairs", 100, 0.42, 0.56);
    ROOT::TThreadedObject<TH1D> MppiLambdaforMatchHist("MppiLambdaforMatchHist", ";m_{p^{#pm}#pi^{#mp}} [GeV];Number of pairs", 100, 1.06, 1.16);
    ROOT::TThreadedObject<TH1D> MKKforThesisHist("MKKforThesisHist", ";m_{K^{#pm}K^{#mp}} [GeV];Number of pairs", 300, 0, 3);
    ROOT::TThreadedObject<TH1D> MKKforMatchHist("MKKforMatchHist", ";m_{K^{#pm}K^{#mp}} [GeV];Number of pairs", 100, 0.42, 0.56);
    ROOT::TThreadedObject<TH1D> MKKLambdaforMatchHist("MKKLambdaforMatchHist", ";m_{K^{#pm}K^{#mp}} [GeV];Number of pairs", 100, 1.06, 1.16);

    //analysis results
    //histogram for measuring event number
    ROOT::TThreadedObject<TH1D> events("events", ";;events", nEvents, 0, nEvents);
    events->GetXaxis()->SetBinLabel(1, "All");
    events->GetXaxis()->SetBinLabel(2, "CPT cuts");
    events->GetXaxis()->SetBinLabel(3, "#splitline{one proton}{each side}");
    events->GetXaxis()->SetBinLabel(4, "#splitline{8 planes}{used in RPs}");
    events->GetXaxis()->SetBinLabel(5, "one primary vertex");
    events->GetXaxis()->SetBinLabel(6, "#splitline{primary vertex}{in central position}");
    events->GetXaxis()->SetBinLabel(7, "#splitline{primary vertex}{unambiguous position}");
    events->GetXaxis()->SetBinLabel(8, "Fiducial cuts");
    events->GetXaxis()->SetBinLabel(9, "#splitline{Elastic}{collision cuts}");
    events->GetXaxis()->SetBinLabel(10, "#splitline{central tracks}{quality check}");

    //additional histograms
    ROOT::TThreadedObject<TH1D> MppluspiminusforMatchHist("MppluspiminusforMatchHist", ";m_{p^{+}#pi^{-}} [GeV];Number of pairs", 100, 0.42, 0.56);
    ROOT::TThreadedObject<TH1D> MpminuspiplusforMatchHist("MpminuspiplusforMatchHist", ";m_{p^{-}#pi^{+}} [GeV];Number of pairs", 100, 0.42, 0.56);
    ROOT::TThreadedObject<TH1D> phiForEtaHist("phiForEtaHist", ";#phi [rad];Number of tracks", 100, -3.14159, 3.14159);
    ROOT::TThreadedObject<TH1D> phiForKHist("phiForKHist", ";#phi [rad];Number of tracks", 100, -3.14159, 3.14159);
    ROOT::TThreadedObject<TH1D> etaForKHist("etaForKHist", ";#phi [rad];Number of tracks", 100, -2.5, 2.5);
    ROOT::TThreadedObject<TH1D> phiForLambdaHist("phiForLambdaHist", ";#phi [rad];Number of tracks", 100, -3.14159, 3.14159);
    ROOT::TThreadedObject<TH1D> etaForLambdaHist("etaForLambdaHist", ";#phi [rad];Number of tracks", 100, -2.5, 2.5);

    ProcessingOutsideLoop outsideprocessing;
    outsideprocessing.AddHistogram(TH1D("MpipiforMatchHist_test", ";m_{#pi^{+}#pi^{-}} [GeV];Number of pairs", 100, 0.42, 0.56));

    // Define the function that will process a subrange of the tree.
    // The function must receive only one parameter, a TTreeReader,
    // and it must be thread safe. To enforce the latter requirement,
    // TThreadedObject histograms will be used.
    //but maybe later
    auto myFunction = [&](TFile *myFile){
        //test if tree is not empty
        TFile *tempFile = new TFile(myFile->GetTitle());
        TTree *tempTree = (TTree *)tempFile->Get("mUPCTree");
        if(tempTree->GetEntries()==0){
            delete tempFile;
            return 0;
        }
        //creating a reader and all stuff
        TTreeReader myReader(tempTree);
        TTreeReaderValue<StUPCEvent> StUPCEventInstance(myReader, "mUPCEvent");
        TTreeReaderValue<StRPEvent> StRPEventInstance(myReader, "mRPEvent");
        //variable initialization
        myReader.Next();
        //histograms
        auto trigHistLocal = trigHist.Get();
        auto oneProtonEUHistLocal = oneProtonEUHist.Get();
        auto oneProtonEDHistLocal = oneProtonEDHist.Get();
        auto oneProtonWUHistLocal = oneProtonWUHist.Get();
        auto oneProtonWDHistLocal = oneProtonWDHist.Get();
        auto planes3Local = planes3.Get();
        auto planes4Local = planes4.Get();
        auto pxpy3EHistLocal = pxpy3EHist.Get();
        auto pxpy4EHistLocal = pxpy4EHist.Get();
        auto pxpy3WHistLocal = pxpy3WHist.Get();
        auto pxpy4WHistLocal = pxpy4WHist.Get();
        auto NverHistLocal = NverHist.Get();
        auto NprimverHistLocal = NprimverHist.Get();
        auto verZHistLocal = verZHist.Get();
        auto verZRPHistLocal = verZRPHist.Get();
        auto verDeltaZHistLocal = verDeltaZHist.Get();
        auto verZRPZTPCHistLocal = verZRPZTPCHist.Get();
        auto pxpyEHistLocal = pxpyEHist.Get();
        auto pxpyWHistLocal = pxpyWHist.Get();
        auto pxpyEOnePrimaryHistLocal = pxpyEOnePrimaryHist.Get();
        auto pxpyWOnePrimaryHistLocal = pxpyWOnePrimaryHist.Get();
        auto pxpyEPrimaryMiddleHistLocal = pxpyEPrimaryMiddleHist.Get();
        auto pxpyWPrimaryMiddleHistLocal = pxpyWPrimaryMiddleHist.Get();
        auto pxpySumHistLocal = pxpySumHist.Get();
        auto NpartHistLocal = NpartHist.Get();
        auto NallpartHistLocal = NallpartHist.Get();
        auto NfitHitsHistLocal = NfitHitsHist.Get();
        auto NdEdxHitsHistLocal = NdEdxHitsHist.Get();
        auto DCAHistLocal = DCAHist.Get();
        auto ptHistLocal = ptHist.Get();
        auto etaHistLocal = etaHist.Get();
        auto dEdxpqHistLocal = dEdxpqHist.Get();
        auto MforThesisHistLocal = MforThesisHist.Get();
        auto MforMatchHistLocal = MforMatchHist.Get();
        auto MLambdaforMatchHistLocal = MLambdaforMatchHist.Get();
        auto MpipiforThesisHistLocal = MpipiforThesisHist.Get();
        auto MpipiforMatchHistLocal = MpipiforMatchHist.Get();
        auto MpipiLambdaforMatchHistLocal = MpipiLambdaforMatchHist.Get();
        auto MppiforThesisHistLocal = MppiforThesisHist.Get();
        auto MppiforMatchHistLocal = MppiforMatchHist.Get();
        auto MppiLambdaforMatchHistLocal = MppiLambdaforMatchHist.Get();
        auto MKKforThesisHistLocal = MKKforThesisHist.Get();
        auto MKKforMatchHistLocal = MKKforMatchHist.Get();
        auto MKKLambdaforMatchHistLocal = MKKLambdaforMatchHist.Get();
        auto eventsLocal = events.Get();
        auto MppluspiminusforMatchHistLocal = MppluspiminusforMatchHist.Get();
        auto MpminuspiplusforMatchHistLocal = MpminuspiplusforMatchHist.Get();
        auto phiForEtaHistLocal = phiForEtaHist.Get();
        auto phiForKHistLocal = phiForKHist.Get();
        auto etaForKHistLocal = etaForKHist.Get();
        auto phiForLambdaHistLocal = phiForLambdaHist.Get();
        auto etaForLambdaHistLocal = etaForLambdaHist.Get();

        ProcessingInsideLoop insideprocessing;
        insideprocessing.GetLocalHistograms(&outsideprocessing);

        int filtered_entries = 0;
        //for changing branch address
        StUPCEvent *tempUPCpointer = StUPCEventInstance.Get();
        StRPEvent *tempRPpointer = StRPEventInstance.Get();

        TClonesArray *trackPointsLocal;
        StUPCRpsTrackPoint *trackPointInstanceLocal;
        TLorentzVector trackVector;
        Double_t verDeltaZ;
        TVector3 pVector;
        vector<TLorentzVector> posPion;
        vector<TLorentzVector> negPion;
        vector<TLorentzVector> posKaon;
        vector<TLorentzVector> negKaon;
        vector<TLorentzVector> posProton;
        vector<TLorentzVector> negProton;
        int whatParticleIsThis;
        do{
            tempUPCpointer = StUPCEventInstance.Get();
            tempRPpointer = StRPEventInstance.Get();

            eventsLocal->Fill(0.0);
            int numberOfOkayTracks = 0;
            for(Int_t i = 0; i<tempUPCpointer->getNumberOfTracks(); i++){
                if(tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kTof)){
                    numberOfOkayTracks++;
                }
            }
            NpartHistLocal->Fill(0.0, numberOfOkayTracks);
            NallpartHistLocal->Fill(0.0, tempUPCpointer->getNumberOfTracks());
            //TESTS & HISTOGRAMS
            //triggers
            if(!CheckTriggers(StUPCEventInstance.Get())){
                continue;
            }
            numberOfOkayTracks = 0;
            for(Int_t i = 0; i<tempUPCpointer->getNumberOfTracks(); i++){
                if(tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kTof)){
                    numberOfOkayTracks++;
                }
            }
            eventsLocal->Fill(1);
            NpartHistLocal->Fill(1, numberOfOkayTracks);
            NallpartHistLocal->Fill(1, tempUPCpointer->getNumberOfTracks());
            for(int i = 0; i<6; i++){
                if(tempUPCpointer->isTrigger(CEPtriggers[i])){
                    switch(CEPtriggers[i]){
                    case 570701:
                    case 570711:
                    case 590701:
                        trigHist->Fill(0);
                        break;

                    case 570705:
                    case 590705:
                    case 590708:
                        trigHist->Fill(1);
                        break;

                    default:
                        break;
                    }

                }
            }

            //one proton each side
            int numberOfTracksPerSide[nSides] = { 0, 0 };
            int numberOfTracksPerBranch[nBranches] = { 0, 0, 0, 0 };
            for(unsigned int k = 0; k<tempRPpointer->getNumberOfTracks(); ++k){
                // Get pointer to k-th track in Roman Pot data collection
                StUPCRpsTrack *trk = tempRPpointer->getTrack(k);
                trk->setEvent(tempRPpointer);
                // Get ID of a branch in which this k-th track was reconstructed
                int j = trk->branch();
                int side = j<2 ? E : W;
                numberOfTracksPerSide[side]++;
                numberOfTracksPerBranch[j]++;
            }
            oneProtonEUHistLocal->Fill(numberOfTracksPerBranch[EU]);
            oneProtonEDHistLocal->Fill(numberOfTracksPerBranch[ED]);
            oneProtonWUHistLocal->Fill(numberOfTracksPerBranch[WU]);
            oneProtonWDHistLocal->Fill(numberOfTracksPerBranch[WD]);
            if(numberOfTracksPerSide[0]!=1||numberOfTracksPerSide[1]!=1){
                continue;
            }
            eventsLocal->Fill(2);
            numberOfOkayTracks = 0;
            for(Int_t i = 0; i<tempUPCpointer->getNumberOfTracks(); i++){
                if(tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kTof)){
                    numberOfOkayTracks++;
                }
            }
            NpartHistLocal->Fill(2, numberOfOkayTracks);
            NallpartHistLocal->Fill(2, tempUPCpointer->getNumberOfTracks());

            //checking planes efficiency & other things with single points
            trackPointsLocal = StRPEventInstance->getTrackPoints();
            for(int i = 0; i<trackPointsLocal->GetEntries(); i++){
                trackPointInstanceLocal = (StUPCRpsTrackPoint *)trackPointsLocal->ConstructedAt(i);

                //planes efficiency
                if(trackPointInstanceLocal->planesUsed()<3){
                    continue;
                }
                //goes through all planes and adds each detection to respective plane
                for(int j = 0; j<4; j++){
                    //if there's 3 planes it fills the one where there's lack of signal, and not the other ones
                    //other ones don't count in efficiency of other planes
                    if(trackPointInstanceLocal->clusterId(j)<0&&trackPointInstanceLocal->planesUsed()==3){
                        planes4Local->Fill(trackPointInstanceLocal->rpId()*4+j, 1.0);
                    }
                    //fills all 4 planes if particle went through 4 of them
                    if(trackPointInstanceLocal->planesUsed()==4){
                        planes4Local->Fill(trackPointInstanceLocal->rpId()*4+j, 1.0);
                        planes3Local->Fill(trackPointInstanceLocal->rpId()*4+j, 1.0);
                    }
                }
            }
            // trackPointsLocal->Clear();

            //checking momentum histogram for at least 3 planes
            for(unsigned int k = 0; k<tempRPpointer->getNumberOfTracks(); ++k){
                // Get pointer to k-th track in Roman Pot data collection
                StUPCRpsTrack *trk = tempRPpointer->getTrack(k);
                trk->setEvent(tempRPpointer);
                //there were problems with apparently not having track point like, entirely???
                //so the first is check point if they do have them
                //and then if points are of good quality
                if(trk->getTrackPoint(0)==nullptr||trk->getTrackPoint(1)==nullptr){
                    continue;
                }
                //check if track has at least 3 of 4 RP planes used
                if(trk->getTrackPoint(0)->planesUsed()<3||trk->getTrackPoint(1)->planesUsed()<3){
                    continue;
                }
                //East
                if(trk->branch()<2){
                    pxpy3EHistLocal->Fill(trk->pVec().X(), trk->pVec().Y());
                }
                //West
                if(trk->branch()>=2){
                    pxpy3WHistLocal->Fill(trk->pVec().X(), trk->pVec().Y());
                }
            }


            //four planes on each TrackPoint
            bool areAllPlanesPresent = true;
            for(unsigned int k = 0; k<tempRPpointer->getNumberOfTracks(); ++k){
                // Get pointer to k-th track in Roman Pot data collection
                StUPCRpsTrack *trk = tempRPpointer->getTrack(k);
                trk->setEvent(tempRPpointer);
                // test if all 8 planes are present
                if(trk->planesUsed()!=8){
                    areAllPlanesPresent = false;
                    break;
                }
            }
            if(!areAllPlanesPresent){
                continue;
            }
            eventsLocal->Fill(3);
            numberOfOkayTracks = 0;
            for(Int_t i = 0; i<tempUPCpointer->getNumberOfTracks(); i++){
                if(tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kTof)){
                    numberOfOkayTracks++;
                }
            }
            NpartHistLocal->Fill(3, numberOfOkayTracks);
            NallpartHistLocal->Fill(3, tempUPCpointer->getNumberOfTracks());
            //checking momentum histogram for 4 planes
            for(unsigned int k = 0; k<tempRPpointer->getNumberOfTracks(); ++k){
                // Get pointer to k-th track in Roman Pot data collection
                StUPCRpsTrack *trk = tempRPpointer->getTrack(k);
                trk->setEvent(tempRPpointer);
                //East
                if(trk->branch()<2){
                    pxpy4EHistLocal->Fill(trk->pVec().X(), trk->pVec().Y());
                }
                //West
                if(trk->branch()>=2){
                    pxpy4WHistLocal->Fill(trk->pVec().X(), trk->pVec().Y());
                }
                trackVector.SetVectM(trk->pVec(), particleMass[Proton]);
            }

            //exactly one primary vertex to which TOF tracks match to
            NverHistLocal->Fill(tempUPCpointer->getNumberOfVertices());
            NprimverHistLocal->Fill(tempUPCpointer->getNPrimVertices());
            vector<Int_t> primaryVertices;
            Int_t VertexId;
            for(Int_t i = 0; i<tempUPCpointer->getNumberOfTracks(); i++){
                VertexId = tempUPCpointer->getTrack(i)->getVertexId();
                if(tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kTof)&&tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kPrimary)&&!(find(primaryVertices.begin(), primaryVertices.end(), VertexId)!=primaryVertices.end())){
                    primaryVertices.push_back(VertexId);
                }
            }
            if(primaryVertices.size()!=1){
                continue;
            }
            VertexId = primaryVertices[0];
            eventsLocal->Fill(4);
            numberOfOkayTracks = 0;
            for(Int_t i = 0; i<tempUPCpointer->getNumberOfTracks(); i++){
                if(tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kTof)){
                    numberOfOkayTracks++;
                }
            }
            NpartHistLocal->Fill(4, numberOfOkayTracks);
            NallpartHistLocal->Fill(4, tempUPCpointer->getNumberOfTracks());
            for(unsigned int k = 0; k<tempRPpointer->getNumberOfTracks(); ++k){
                // Get pointer to k-th track in Roman Pot data collection
                StUPCRpsTrack *trk = tempRPpointer->getTrack(k);
                trk->setEvent(tempRPpointer);
                //East
                if(trk->branch()<2){
                    pxpyEOnePrimaryHistLocal->Fill(trk->pVec().X(), trk->pVec().Y());
                }
                //West
                if(trk->branch()>=2){
                    pxpyWOnePrimaryHistLocal->Fill(trk->pVec().X(), trk->pVec().Y());
                }
            }

            //primary vertex is placed within 80cm of the centre
            verZHistLocal->Fill(tempUPCpointer->getVertexId(VertexId)->getPosZ());
            if(abs(tempUPCpointer->getVertexId(VertexId)->getPosZ())>=80){
                continue;
            }
            eventsLocal->Fill(5);
            numberOfOkayTracks = 0;
            for(Int_t i = 0; i<tempUPCpointer->getNumberOfTracks(); i++){
                if(tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kTof)){
                    numberOfOkayTracks++;
                }
            }
            NpartHistLocal->Fill(5, numberOfOkayTracks);
            NallpartHistLocal->Fill(5, tempUPCpointer->getNumberOfTracks());
            for(unsigned int k = 0; k<tempRPpointer->getNumberOfTracks(); ++k){
                // Get pointer to k-th track in Roman Pot data collection
                StUPCRpsTrack *trk = tempRPpointer->getTrack(k);
                trk->setEvent(tempRPpointer);
                //East
                if(trk->branch()<2){
                    pxpyEPrimaryMiddleHistLocal->Fill(trk->pVec().X(), trk->pVec().Y());
                }
                //West
                if(trk->branch()>=2){
                    pxpyWPrimaryMiddleHistLocal->Fill(trk->pVec().X(), trk->pVec().Y());
                }
            }

            //vertexes from TPC and RP are not too far away (36cm)
            bool areBothRPTracksWithTime = true;
            Double_t time = 0;
            for(unsigned int k = 0; k<tempRPpointer->getNumberOfTracks(); ++k){
                // Get pointer to k-th track in Roman Pot data collection
                StUPCRpsTrack *trk = tempRPpointer->getTrack(k);
                trk->setEvent(tempRPpointer);
                // Get ID of a branch in which this k-th track was reconstructed
                int j = trk->branch();
                //deletes those tracks which do not have confirmed time
                if(trk->time()<0){
                    areBothRPTracksWithTime = false;
                    break;
                }
                int side = j<2 ? E : W;
                //to not make the loop go twice, also west side is the "+" one
                if(side==E){
                    time -= trk->time();
                }
                if(side==W){
                    time += trk->time();
                }
            }
            if(!areBothRPTracksWithTime){
                continue;
            }
            //*100 at the end it to convert from m to cm
            verDeltaZ = tempUPCpointer->getVertexId(VertexId)->getPosZ()-299792458.0*time/2*100;
            verDeltaZHistLocal->Fill(verDeltaZ);
            verZRPHistLocal->Fill(299792458.0*time/2*100);
            verZRPZTPCHistLocal->Fill(299792458.0*time/2*100, tempUPCpointer->getVertexId(VertexId)->getPosZ());
            if(abs(verDeltaZ)>36){
                continue;
            }
            eventsLocal->Fill(6);
            numberOfOkayTracks = 0;
            for(Int_t i = 0; i<tempUPCpointer->getNumberOfTracks(); i++){
                if(tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kTof)){
                    numberOfOkayTracks++;
                }
            }
            NpartHistLocal->Fill(6, numberOfOkayTracks);
            NallpartHistLocal->Fill(6, tempUPCpointer->getNumberOfTracks());

            //fiducial cuts
            for(unsigned int k = 0; k<tempRPpointer->getNumberOfTracks(); ++k){
                // Get pointer to k-th track in Roman Pot data collection
                StUPCRpsTrack *trk = tempRPpointer->getTrack(k);
                trk->setEvent(tempRPpointer);
                //East
                if(trk->branch()<2){
                    pxpyEHistLocal->Fill(trk->pVec().X(), trk->pVec().Y());
                }
                //West
                if(trk->branch()>=2){
                    pxpyWHistLocal->Fill(trk->pVec().X(), trk->pVec().Y());
                }
            }
            bool passedFiducialCuts = true;
            //px barrier
            Double_t pxbarrier = -0.27;
            //py low nand high barrier
            Double_t pylowbarrier = 0.4;
            Double_t pyhighbarrier = 0.8;
            //p circle parameters
            Double_t pxcenter = -0.6;
            Double_t pycenter = 0;
            Double_t pradius = 1.1;
            for(unsigned int k = 0; k<tempRPpointer->getNumberOfTracks(); ++k){
                // Get pointer to k-th track in Roman Pot data collection
                StUPCRpsTrack *trk = tempRPpointer->getTrack(k);
                trk->setEvent(tempRPpointer);
                pVector = trk->pVec();
                bool f1 = pow(pVector.X()-pxcenter, 2)+pow(pVector.Y()-pycenter, 2)<pradius*pradius;
                bool f2 = pylowbarrier<abs(pVector.Y())&&abs(pVector.Y())<pyhighbarrier;
                bool f3 = pxbarrier<pVector.X();
                if(!(f1&&f2&&f3)){
                    passedFiducialCuts = false;
                    break;
                }
            }
            if(!passedFiducialCuts){
                continue;
            }
            eventsLocal->Fill(7);
            numberOfOkayTracks = 0;
            for(Int_t i = 0; i<tempUPCpointer->getNumberOfTracks(); i++){
                if(tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kTof)){
                    numberOfOkayTracks++;
                }
            }
            NpartHistLocal->Fill(7, numberOfOkayTracks);
            NallpartHistLocal->Fill(7, tempUPCpointer->getNumberOfTracks());

            //elastic cuts
            TLorentzVector RPsum = { 0,0,0,0 };
            //checking total momentum
            for(unsigned int k = 0; k<tempRPpointer->getNumberOfTracks(); ++k){
                // Get pointer to k-th track in Roman Pot data collection
                StUPCRpsTrack *trk = tempRPpointer->getTrack(k);
                trk->setEvent(tempRPpointer);
                trackVector.SetVectM(trk->pVec(), particleMass[Proton]);
                RPsum += trackVector;
            }
            pxpySumHistLocal->Fill(RPsum.X(), RPsum.Y());
            //xposition of the gauss
            Double_t xGauss = -0.035;
            //yposition of the gauss
            Double_t yGauss = 0.01;
            //radius of the gauss
            //was 0.07, 0.1 FOR TEST
            Double_t rGauss = 0.07;
            if(pow(RPsum.X()-xGauss, 2)+pow(RPsum.Y()-yGauss, 2)<rGauss*rGauss){
                continue;
            }
            eventsLocal->Fill(8);
            numberOfOkayTracks = 0;
            for(Int_t i = 0; i<tempUPCpointer->getNumberOfTracks(); i++){
                if(tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kTof)){
                    numberOfOkayTracks++;
                }
            }
            NpartHistLocal->Fill(8, numberOfOkayTracks);
            NallpartHistLocal->Fill(8, tempUPCpointer->getNumberOfTracks());

            //quality test of central tracks
            //1. Number of hits
            //2. DCA in R and z (the second one also doubles as spread limiter)
            //3. kinematic range for optimal measurement
            int TracksValid = 0;
            TLorentzVector UPCsum = { 0,0,0,0 };
            for(Int_t i = 0; i<tempUPCpointer->getNumberOfTracks(); i++){
                //basic tests
                if(!tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kTof)){
                    continue;
                }
                if(tempUPCpointer->getTrack(i)->getTofPathLength()<=0||tempUPCpointer->getTrack(i)->getTofTime()<=0){
                    continue;
                }
                //quality
                NfitHitsHistLocal->Fill(tempUPCpointer->getTrack(i)->getNhitsFit());
                NdEdxHitsHistLocal->Fill(tempUPCpointer->getTrack(i)->getNhitsDEdx());
                if(tempUPCpointer->getTrack(i)->getNhitsDEdx()<15||tempUPCpointer->getTrack(i)->getNhitsFit()<25){
                    continue;
                }
                DCAHistLocal->Fill(sqrt(pow(tempUPCpointer->getTrack(i)->getDcaXY(), 2)+pow(tempUPCpointer->getTrack(i)->getDcaZ(), 2)));
                if((pow(tempUPCpointer->getTrack(i)->getDcaXY(), 2)+pow(tempUPCpointer->getTrack(i)->getDcaZ(), 2))<1&&tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kPrimary)){
                    continue;
                }
                ptHistLocal->Fill(tempUPCpointer->getTrack(i)->getPt());
                etaHistLocal->Fill(tempUPCpointer->getTrack(i)->getEta());
                phiForEtaHistLocal->Fill(tempUPCpointer->getTrack(i)->getPhi());
                if(tempUPCpointer->getTrack(i)->getPt()<=0.2||abs(tempUPCpointer->getTrack(i)->getEta())>=0.7){
                    continue;
                }
                tempUPCpointer->getTrack(i)->getLorentzVector(trackVector, particleMass[Pion]);
                dEdxpqHistLocal->Fill(trackVector.P()*tempUPCpointer->getTrack(i)->getCharge(), tempUPCpointer->getTrack(i)->getDEdxSignal());
                if(abs(tempUPCpointer->getTrack(i)->getNSigmasTPCProton())<3){
                    tempUPCpointer->getTrack(i)->getLorentzVector(trackVector, particleMass[Proton]);
                    whatParticleIsThis = Proton;
                } else if(abs(tempUPCpointer->getTrack(i)->getNSigmasTPCKaon())<3){
                    tempUPCpointer->getTrack(i)->getLorentzVector(trackVector, particleMass[Kaon]);
                    whatParticleIsThis = Kaon;
                } else if(abs(tempUPCpointer->getTrack(i)->getNSigmasTPCPion())<3){
                    tempUPCpointer->getTrack(i)->getLorentzVector(trackVector, particleMass[Pion]);
                    whatParticleIsThis = Pion;
                } else{
                    continue;
                }
                //at this point track is Valid(TM)
                TracksValid++;
                UPCsum += trackVector;
                if(tempUPCpointer->getTrack(i)->getCharge()==1){
                    switch(whatParticleIsThis){
                    case Proton:
                        posProton.push_back(trackVector);
                        break;
                    case Kaon:
                        posKaon.push_back(trackVector);
                        break;
                    case Pion:
                        posPion.push_back(trackVector);
                        break;
                    default:
                        break;
                    }
                } else if(tempUPCpointer->getTrack(i)->getCharge()==-1){
                    switch(whatParticleIsThis){
                    case Proton:
                        negProton.push_back(trackVector);
                        break;
                    case Kaon:
                        negKaon.push_back(trackVector);
                        break;
                    case Pion:
                        negPion.push_back(trackVector);
                        break;
                    default:
                        break;
                    }
                }
            }
            if(TracksValid==0){
                continue;
            }

            //triggers at the end
            for(int i = 0; i<6; i++){
                if(tempUPCpointer->isTrigger(CEPtriggers[i])){
                    switch(CEPtriggers[i]){
                    case 570701:
                    case 570711:
                    case 590701:
                        trigHist->Fill(2);
                        break;

                    case 570705:
                    case 590705:
                    case 590708:
                        trigHist->Fill(3);
                        break;

                    default:
                        break;
                    }

                }
            }
            //more histograms
            eventsLocal->Fill(9);
            NpartHistLocal->Fill(9, TracksValid);
            NallpartHistLocal->Fill(9, tempUPCpointer->getNumberOfTracks());


            //for mass
            TLorentzVector tempSum = { 0, 0, 0, 0 };
            //invariant mass, p+ pi-
            for(long unsigned int posi = 0; posi<posProton.size(); posi++){
                for(long unsigned int negi = 0; negi<negPion.size(); negi++){
                    tempSum = posProton[posi]+negPion[negi];
                    MforThesisHistLocal->Fill(tempSum.M());
                    MforMatchHistLocal->Fill(tempSum.M());
                    MLambdaforMatchHistLocal->Fill(tempSum.M());
                    //specific
                    MppiforThesisHistLocal->Fill(tempSum.M());
                    MppiforMatchHistLocal->Fill(tempSum.M());
                    MppiLambdaforMatchHistLocal->Fill(tempSum.M());
                    //additional
                    if(tempSum.M()>1.11&&tempSum.M()<1.12){
                        phiForLambdaHistLocal->Fill(tempSum.Phi());
                        etaForLambdaHistLocal->Fill(tempSum.Eta());
                    }
                    MppluspiminusforMatchHistLocal->Fill(tempSum.M());
                }
            }
            //invariant mass, pi+ p-
            for(long unsigned int posi = 0; posi<posPion.size(); posi++){
                for(long unsigned int negi = 0; negi<negProton.size(); negi++){
                    tempSum = posPion[posi]+negProton[negi];
                    MforThesisHistLocal->Fill(tempSum.M());
                    MforMatchHistLocal->Fill(tempSum.M());
                    MLambdaforMatchHistLocal->Fill(tempSum.M());
                    //specific
                    MppiforThesisHistLocal->Fill(tempSum.M());
                    MppiforMatchHistLocal->Fill(tempSum.M());
                    MppiLambdaforMatchHistLocal->Fill(tempSum.M());
                    //additional
                    if(tempSum.M()>1.11&&tempSum.M()<1.12){
                        phiForLambdaHistLocal->Fill(tempSum.Phi());
                        etaForLambdaHistLocal->Fill(tempSum.Eta());
                    }
                    MpminuspiplusforMatchHistLocal->Fill(tempSum.M());
                }
            }
            //invariant mass, pi+ pi-
            for(long unsigned int posi = 0; posi<posPion.size(); posi++){
                for(long unsigned int negi = 0; negi<negPion.size(); negi++){
                    tempSum = posPion[posi]+negPion[negi];
                    MforThesisHistLocal->Fill(tempSum.M());
                    MforMatchHistLocal->Fill(tempSum.M());
                    MLambdaforMatchHistLocal->Fill(tempSum.M());
                    //specific
                    MpipiforThesisHistLocal->Fill(tempSum.M());
                    MpipiforMatchHistLocal->Fill(tempSum.M());
                    MpipiLambdaforMatchHistLocal->Fill(tempSum.M());
                    //additional
                    if(tempSum.M()>0.492&&tempSum.M()<0.502){
                        phiForKHistLocal->Fill(tempSum.Phi());
                        etaForKHistLocal->Fill(tempSum.Eta());
                    }
                    insideprocessing.Fill(0, tempSum.M());
                }
            }
            //invariant mass, K+ K-
            for(long unsigned int posi = 0; posi<posKaon.size(); posi++){
                for(long unsigned int negi = 0; negi<negKaon.size(); negi++){
                    tempSum = posKaon[posi]+negKaon[negi];
                    MforThesisHistLocal->Fill(tempSum.M());
                    MforMatchHistLocal->Fill(tempSum.M());
                    MLambdaforMatchHistLocal->Fill(tempSum.M());
                    //specific
                    MKKforThesisHistLocal->Fill(tempSum.M());
                    MKKforMatchHistLocal->Fill(tempSum.M());
                    MKKLambdaforMatchHistLocal->Fill(tempSum.M());
                }
            }

            filtered_entries++;

            posPion.clear();
            negPion.clear();
            posKaon.clear();
            negKaon.clear();
            posProton.clear();
            negProton.clear();
        } while(myReader.Next());

        //waiting for file opening to check if there were any filtered entries
        if(filtered_entries==0){
            cout<<"Finished operation on file "<<myFile->GetTitle()<<endl;
            cout<<"Analyzed "<<tempTree->GetEntries()<<" entries"<<endl;
            cout<<"There were 0 filtered entries"<<endl;

            delete tempFile;
            return 0;
        }

        cout<<"Finished operation on file "<<myFile->GetTitle()<<endl;
        cout<<"Analyzed "<<tempTree->GetEntries()<<" entries"<<endl;
        cout<<"Filtered "<<filtered_entries<<" entries"<<endl;

        delete tempFile;

        return filtered_entries;
        };

    auto redFunction = [](const std::vector<int> &mapV){
        return std::accumulate(mapV.begin(), mapV.end(), 0);
        };

    int filtered_entries = 0;
    vector<TFile *> listOfFiles;
    //creating the list of TFile*
    TObjArray *tempList = upcChain->GetListOfFiles();
    for(int i = 0; i<tempList->GetEntries(); i++){
        listOfFiles.push_back((TFile *)(tempList->At(i)));
    }

    // Create a TreeProcessor: specify the file and the tree in it
    ROOT::TThreadExecutor TreeProcessor(nthreads);
    // Launch the parallel processing of the tree
    filtered_entries = TreeProcessor.MapReduce(myFunction, listOfFiles, redFunction);
    // Use the TThreadedObject::Merge method to merge the thread private tree
    // into the final result

    auto trigHistFinal = trigHist.Merge();
    auto oneProtonEUHistFinal = oneProtonEUHist.Merge();
    auto oneProtonEDHistFinal = oneProtonEDHist.Merge();
    auto oneProtonWUHistFinal = oneProtonWUHist.Merge();
    auto oneProtonWDHistFinal = oneProtonWDHist.Merge();
    auto planes3Final = planes3.Merge();
    auto planes4Final = planes4.Merge();
    auto pxpy3EHistFinal = pxpy3EHist.Merge();
    auto pxpy4EHistFinal = pxpy4EHist.Merge();
    auto pxpy3WHistFinal = pxpy3WHist.Merge();
    auto pxpy4WHistFinal = pxpy4WHist.Merge();
    auto NverHistFinal = NverHist.Merge();
    auto NprimverHistFinal = NprimverHist.Merge();
    auto verZHistFinal = verZHist.Merge();
    auto verZRPHistFinal = verZRPHist.Merge();
    auto verDeltaZHistFinal = verDeltaZHist.Merge();
    auto verZRPZTPCHistFinal = verZRPZTPCHist.Merge();
    auto pxpyEHistFinal = pxpyEHist.Merge();
    auto pxpyWHistFinal = pxpyWHist.Merge();
    auto pxpyEOnePrimaryHistFinal = pxpyEOnePrimaryHist.Merge();
    auto pxpyWOnePrimaryHistFinal = pxpyWOnePrimaryHist.Merge();
    auto pxpyEPrimaryMiddleHistFinal = pxpyEPrimaryMiddleHist.Merge();
    auto pxpyWPrimaryMiddleHistFinal = pxpyWPrimaryMiddleHist.Merge();
    auto pxpySumHistFinal = pxpySumHist.Merge();
    auto NpartHistFinal = NpartHist.Merge();
    auto NallpartHistFinal = NallpartHist.Merge();
    auto NfitHitsHistFinal = NfitHitsHist.Merge();
    auto NdEdxHitsHistFinal = NdEdxHitsHist.Merge();
    auto DCAHistFinal = DCAHist.Merge();
    auto ptHistFinal = ptHist.Merge();
    auto etaHistFinal = etaHist.Merge();
    auto dEdxpqHistFinal = dEdxpqHist.Merge();
    auto MforThesisHistFinal = MforThesisHist.Merge();
    auto MforMatchHistFinal = MforMatchHist.Merge();
    auto MLambdaforMatchHistFinal = MLambdaforMatchHist.Merge();
    auto MpipiforThesisHistFinal = MpipiforThesisHist.Merge();
    auto MpipiforMatchHistFinal = MpipiforMatchHist.Merge();
    auto MpipiLambdaforMatchHistFinal = MpipiLambdaforMatchHist.Merge();
    auto MppiforThesisHistFinal = MppiforThesisHist.Merge();
    auto MppiforMatchHistFinal = MppiforMatchHist.Merge();
    auto MppiLambdaforMatchHistFinal = MppiLambdaforMatchHist.Merge();
    auto MKKforThesisHistFinal = MKKforThesisHist.Merge();
    auto MKKforMatchHistFinal = MKKforMatchHist.Merge();
    auto MKKLambdaforMatchHistFinal = MKKLambdaforMatchHist.Merge();
    auto eventsFinal = events.Merge();
    auto MppluspiminusforMatchHistFinal = MppluspiminusforMatchHist.Merge();
    auto MpminuspiplusforMatchHistFinal = MpminuspiplusforMatchHist.Merge();
    auto phiForEtaHistFinal = phiForEtaHist.Merge();
    auto phiForKHistFinal = phiForKHist.Merge();
    auto etaForKHistFinal = etaForKHist.Merge();
    auto phiForLambdaHistFinal = phiForLambdaHist.Merge();
    auto etaForLambdaHistFinal = etaForLambdaHist.Merge();

    outsideprocessing.Merge();

    //setting up a tree & output file
    string outfileName = outputFolder+"AnaOutput_OnlyHistForOverleaf_with_new_histogram_handling.root";
    cout<<"Created output file "<<outfileName<<endl;
    TFile *outputFileHist = TFile::Open(outfileName.c_str(), "recreate");

    outputFileHist->cd();

    trigHistFinal->Write();
    oneProtonEUHistFinal->Write();
    oneProtonEDHistFinal->Write();
    oneProtonWUHistFinal->Write();
    oneProtonWDHistFinal->Write();
    planes3Final->Write();
    planes4Final->Write();
    pxpy3EHistFinal->Write();
    pxpy4EHistFinal->Write();
    pxpy3WHistFinal->Write();
    pxpy4WHistFinal->Write();
    NverHistFinal->Write();
    NprimverHistFinal->Write();
    verZHistFinal->Write();
    verZRPHistFinal->Write();
    verDeltaZHistFinal->Write();
    verZRPZTPCHistFinal->Write();
    pxpyEHistFinal->Write();
    pxpyWHistFinal->Write();
    pxpyEOnePrimaryHistFinal->Write();
    pxpyWOnePrimaryHistFinal->Write();
    pxpyEPrimaryMiddleHistFinal->Write();
    pxpyWPrimaryMiddleHistFinal->Write();
    pxpySumHistFinal->Write();
    NpartHistFinal->Write();
    NallpartHistFinal->Write();
    NfitHitsHistFinal->Write();
    NdEdxHitsHistFinal->Write();
    DCAHistFinal->Write();
    ptHistFinal->Write();
    etaHistFinal->Write();
    dEdxpqHistFinal->Write();
    MforThesisHistFinal->Write();
    MforMatchHistFinal->Write();
    MLambdaforMatchHistFinal->Write();
    MpipiforThesisHistFinal->Write();
    MpipiforMatchHistFinal->Write();
    MpipiLambdaforMatchHistFinal->Write();
    MppiforThesisHistFinal->Write();
    MppiforMatchHistFinal->Write();
    MppiLambdaforMatchHistFinal->Write();
    MKKforThesisHistFinal->Write();
    MKKforMatchHistFinal->Write();
    MKKLambdaforMatchHistFinal->Write();
    eventsFinal->Write();
    MppluspiminusforMatchHistFinal->Write();
    MpminuspiplusforMatchHistFinal->Write();
    phiForEtaHistFinal->Write();
    phiForKHistFinal->Write();
    etaForKHistFinal->Write();
    phiForLambdaHistFinal->Write();
    etaForLambdaHistFinal->Write();

    outsideprocessing.SaveToFile(outputFileHist);

    outputFileHist->Close();

    cout<<"Finished processing "<<endl;
    cout<<"Analyzed total "<<upcChain->GetEntries()<<" entries"<<endl;
    cout<<"Filtered total "<<filtered_entries<<" entries"<<endl;
    cout<<"Ending Analysis... GOOD BYE!"<<endl;
    return 0;
}//main

bool ConnectInput(int argc, char **argv, TChain *fileChain){
    int fileId = -1;
    string line;
    int lineId = 0;

    const string &input = argv[1];
    cout<<"Using list "<<input<<endl;
    if(input.find(".list")!=string::npos){
        cout<<"Input from chain"<<endl;
        ifstream instr(input.c_str());
        if(!instr.is_open()){
            cout<<"Couldn't open: "<<input.c_str()<<endl;
            return false;
        }
        //for testing if file exists
        TFile *infile;

        while(getline(instr, line)){
            if(fileId==lineId||fileId==-1){
                fileChain->AddFile(line.c_str());
                infile = TFile::Open(line.c_str(), "read");
                if(!infile){
                    cout<<"Couldn't open: "<<line.c_str()<<endl;
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

bool CheckTriggers(StUPCEvent *localupcEvt){

    bool CPTtrigger = false;
    for(int var = 0; var<nTriggers; ++var){
        if(localupcEvt->isTrigger(triggerID[var])){
            //Checked if it is CPT trigger
            for(int i = 0; i<*(&CEPtriggers+1)-CEPtriggers; ++i)
                if(triggerID[var]==CEPtriggers[i])
                    CPTtrigger = true;
        }
    }

    return CPTtrigger;
}

long long GetFileSize(string filename){
    struct stat64 stat_buf;
    int rc = stat64(filename.c_str(), &stat_buf);
    return rc==0 ? stat_buf.st_size : -1;
}