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

using namespace std;

// enums are very usefull 
enum { kAll = 1, kCPT,  kRP, kOneVertex, kTPCTOF, 
    kTotQ, kMax};
enum SIDE {E = 0, East = 0, W = 1, West = 1, nSides};
enum PARTICLES {Pion = 0, Kaon = 1, Proton = 2, nParticles};
enum BRANCH_ID { EU, ED, WU, WD, nBranches };
enum RP_ID {E1U, E1D, E2U, E2D, W1U, W1D, W2U, W2D, nRomanPots};

const double particleMass[nParticles] = { 0.13957, 0.493677, 0.93827 }; // pion, kaon, proton in GeV /c^2 
const int nTriggers = 17;
const int triggerID[] = { 570209, 570219, 570229, 570701, 570702, 570703, 570704, 570705, 
                  570709, 570711, 570712, 570719, 590701, 590703, 590705, 590708, 590709};
// 570702 RP_UPC // 570712 RP_UPC // 570703 RP_SDT // 570709 RP_ET // 570719 RP_ET // 570701 RP_CPT2 // 570711 RP_CPT2 // 570705 RP_CPT2noBBCL // 570704 RP_Zerobias // 590703 RP_SDT // 590709 RP_ET // 590701 RP_CPT2 // 590705 RP_CPT2noBBCL // 590708 RP_CPTnoBBCL // 570209 JPsi*HTTP // 570219 JPsi*HTTP // 570229 JPsi*HTTP
const int CEPtriggers[] = { 570701, 570705, 570711, 590701, 590705, 590708};

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

    if(!ConnectInput(argc, argv, upcChain))
    {
        cout << "Wrong input parameters..." << endl; 
        return 1;
    }

    const string& outputFolder = argv[2];

    //HISTOGRAMS
    //histograms for measuring plane efficiency, there are 4 planes per roman pot
    ROOT::TThreadedObject<TH1D> planes3("planes3Hist", "planes3Hist", nRomanPots*4, 0, nRomanPots*4);
    ROOT::TThreadedObject<TH1D> planes4("planes4Hist", "planes4Hist", nRomanPots*4, 0, nRomanPots*4);
    //histogram for measuring event number
    int nEvents = 8;
    ROOT::TThreadedObject<TH1D> events("events", ";;events", nEvents, 0, nEvents);
    events->GetXaxis()->SetBinLabel(1, "All");
    events->GetXaxis()->SetBinLabel(2, "CPT cuts");
    events->GetXaxis()->SetBinLabel(3, "#splitline{one proton}{each side}");
    events->GetXaxis()->SetBinLabel(4, "#splitline{8 planes}{used in RPs}");
    events->GetXaxis()->SetBinLabel(5, "one primary vertex");
    events->GetXaxis()->SetBinLabel(6, "#splitline{primary vertex}{in central position}");
    events->GetXaxis()->SetBinLabel(7, "#splitline{primary vertex}{unambiguous position}");
    events->GetXaxis()->SetBinLabel(8, "#splitline{central tracks}{quality check}");
    //particle number
    ROOT::TThreadedObject<TH2D> NpartHist("NpartHist", ";;Number of created TOF-matched particles", nEvents, 0, nEvents, 40, 0, 40);
    NpartHist->GetXaxis()->SetBinLabel(1, "All");
    NpartHist->GetXaxis()->SetBinLabel(2, "CPT cuts");
    NpartHist->GetXaxis()->SetBinLabel(3, "#splitline{one proton}{each side}");
    NpartHist->GetXaxis()->SetBinLabel(4, "#splitline{8 planes}{used in RPs}");
    NpartHist->GetXaxis()->SetBinLabel(5, "one primary vertex");
    NpartHist->GetXaxis()->SetBinLabel(6, "#splitline{primary vertex}{in central position}");
    NpartHist->GetXaxis()->SetBinLabel(7, "#splitline{primary vertex}{unambiguous position}");
    NpartHist->GetXaxis()->SetBinLabel(8, "#splitline{central tracks}{quality check}");
    //particle number without TOF checking
    ROOT::TThreadedObject<TH2D> NallpartHist("NallpartHist", ";;Number of created TOF-matched particles", nEvents, 0, nEvents, 40, 0, 40);
    NallpartHist->GetXaxis()->SetBinLabel(1, "All");
    NallpartHist->GetXaxis()->SetBinLabel(2, "CPT cuts");
    NallpartHist->GetXaxis()->SetBinLabel(3, "#splitline{one proton}{each side}");
    NallpartHist->GetXaxis()->SetBinLabel(4, "#splitline{8 planes}{used in RPs}");
    NallpartHist->GetXaxis()->SetBinLabel(5, "one primary vertex");
    NallpartHist->GetXaxis()->SetBinLabel(6, "#splitline{primary vertex}{in central position}");
    NallpartHist->GetXaxis()->SetBinLabel(7, "#splitline{primary vertex}{unambiguous position}");
    NallpartHist->GetXaxis()->SetBinLabel(8, "#splitline{central tracks}{quality check}");
    //momentum, after 3 and 4 planes, with repect to sides
    ROOT::TThreadedObject<TH2D> pxpy3EHist("pxpy3EHist", ";p_{x} [GeV];p_{y} [GeV]", 400, -2, 2, 400, -2, 2);
    ROOT::TThreadedObject<TH2D> pxpy4EHist("pxpy4EHist", ";p_{x} [GeV];p_{y} [GeV]", 400, -2, 2, 400, -2, 2);
    ROOT::TThreadedObject<TH2D> pxpy3WHist("pxpy3WHist", ";p_{x} [GeV];p_{y} [GeV]", 400, -2, 2, 400, -2, 2);
    ROOT::TThreadedObject<TH2D> pxpy4WHist("pxpy4WHist", ";p_{x} [GeV];p_{y} [GeV]", 400, -2, 2, 400, -2, 2);
    //px vs px and py vs py
    ROOT::TThreadedObject<TH2D> pxpxHist("pxpxHist", ";p_{xRP} [GeV];p_{xTPC} [GeV]", 800, -4, 4, 800, -4, 4);
    ROOT::TThreadedObject<TH2D> pypyHist("pypyHist", ";p_{yRP} [GeV];p_{yTPC} [GeV]", 800, -4, 4, 800, -4, 4);
    //vertex z position, TPC
    ROOT::TThreadedObject<TH1D> verZHist("verZHist", ";VP_{zTPC} [cm];events", 120, -240, 240);
    //vertex z position, RP
    ROOT::TThreadedObject<TH1D> verZRPHist("verZRPHist", ";VP_{zRP} [cm];events", 120, -240, 240);
    //vertex position difference
    ROOT::TThreadedObject<TH1D> verDeltaZHist("verDeltaZHist", ";VP_{zTPC}-VP_{zRP} [cm];events", 240, -240, 240);
    //correlation vertex position TPC RP
    ROOT::TThreadedObject<TH2D> verZRPZTPCHist("verZRPZTPCHist", ";VP_{zRP} [cm];VP_{zTPC} [cm]", 240, -240, 240, 240, -240, 240);
    //final one for fiducial cuts
    ROOT::TThreadedObject<TH2D> pxpyEHist("pxpyEHist", ";p_{x} [GeV];p_{y} [GeV]", 300, -1.5, 1.5, 300, -1.5, 1.5);
    ROOT::TThreadedObject<TH2D> pxpyWHist("pxpyWHist", ";p_{x} [GeV];p_{y} [GeV]", 300, -1.5, 1.5, 300, -1.5, 1.5);
    //sum of proton momenta
    ROOT::TThreadedObject<TH2D> pxpySumEHist("pxpySumEHist", ";#Sigma p_{x} [GeV];#Sigma p_{y} [GeV]", 400, -2, 2, 400, -2, 2);
    //lacking momentum
    ROOT::TThreadedObject<TH2D> pxpylackHist("pxpylackHist", ";#Sigma p_{x} [GeV];#Sigma p_{y} [GeV]", 400, -2, 2, 400, -2, 2);
    //dE/dx dependant of p
    ROOT::TThreadedObject<TH2D> dEdxpqHist("dEdxpqHist", ";pq [GeV];dE/dx [GeV/cm]", 1000, -4, 4, 1000, 0, 4e-5);
    //number of vertices and prime vertices
    ROOT::TThreadedObject<TH1D> NverHist("NverHist", ";Vertices;events", 8, 0, 8);
    ROOT::TThreadedObject<TH1D> NprimverHist("NprimverHist", ";Primary vertices;events", 8, 0, 8);
    //NfitHits and NdEdxHits
    ROOT::TThreadedObject<TH1D> NfitHitsHist("NfitHitsHist", ";N_{TPCfit};tracks", 50, 0, 50);
    ROOT::TThreadedObject<TH1D> NdEdxHitsHist("NdEdxHitsHist", ";N_{dE/dx};tracks", 50, 0, 50);
    //DCAR and DCAZ
    ROOT::TThreadedObject<TH1D> DCARHist("DCARHist", ";R [cm];tracks", 160, 0, 3.2);
    ROOT::TThreadedObject<TH1D> DCAZHist("DCAZHist", ";z [cm];tracks", 160, -4, 4);
    //pt and eta
    ROOT::TThreadedObject<TH1D> ptHist("ptHist", ";p_{T} [GeV];tracks", 100, 0, 1);
    ROOT::TThreadedObject<TH1D> etaHist("etaHist", ";#eta;tracks", 100, -2.5, 2.5);
    //trigger histograms
    ROOT::TThreadedObject<TH1D> trigBeginningHist("trigBeginningHist", ";;events", 6, 0, 6);
    trigBeginningHist->GetXaxis()->SetBinLabel(1, "570701");
    trigBeginningHist->GetXaxis()->SetBinLabel(2, "570705");
    trigBeginningHist->GetXaxis()->SetBinLabel(3, "570711");
    trigBeginningHist->GetXaxis()->SetBinLabel(4, "590701");
    trigBeginningHist->GetXaxis()->SetBinLabel(5, "590705");
    trigBeginningHist->GetXaxis()->SetBinLabel(6, "590708");
    ROOT::TThreadedObject<TH1D> trigEndingHist("trigEndingHist", ";;events", 6, 0, 6);
    trigEndingHist->GetXaxis()->SetBinLabel(1, "570701");
    trigEndingHist->GetXaxis()->SetBinLabel(2, "570705");
    trigEndingHist->GetXaxis()->SetBinLabel(3, "570711");
    trigEndingHist->GetXaxis()->SetBinLabel(4, "590701");
    trigEndingHist->GetXaxis()->SetBinLabel(5, "590705");
    trigEndingHist->GetXaxis()->SetBinLabel(6, "590708");
    //mass histograms
    ROOT::TThreadedObject<TH1D> MforThesisHist("MforThesisHist", ";m_{inv} [GeV];", 300, 0, 3);
    ROOT::TThreadedObject<TH1D> MbcgforThesisHist("MbcgforThesisHist", ";m_{inv} [GeV];", 300, 0, 3);
    ROOT::TThreadedObject<TH1D> MforMatchHist("MforMatchHist", ";m_{inv} [GeV];", 100, 0.42, 0.56);
    ROOT::TThreadedObject<TH1D> MbcgforMatchHist("MbcgforMatchHist", ";m_{inv} [GeV];", 100, 0.42, 0.56);
    //DCAR, DCAZ and deltaDCAZ for pions that create K0
    ROOT::TThreadedObject<TH1D> pionDCARHist("pionDCARHist", ";R [cm];tracks", 160, 0, 1.5);
    ROOT::TThreadedObject<TH1D> pionDCAZHist("pionDCAZHist", ";z [cm];tracks", 160, -1, 1);
    ROOT::TThreadedObject<TH1D> piondeltaDCAZHist("piondeltaDCAZHist", ";#Delta z [cm];pion pairs", 160, -2, 2);
    //VPz, DCAR, DCAZ, dTOF for the worst cases
    ROOT::TThreadedObject<TH1D> verZBadHist("verZBadHist", ";VP_{zTPC} [cm];events", 120, -240, 240);
    ROOT::TThreadedObject<TH1D> DCARBadHist("DCARBadHist", ";R [cm];tracks", 160, 0, 3.2);
    ROOT::TThreadedObject<TH1D> DCAZBadHist("DCAZBadHist", ";z [cm];tracks", 160, -4, 4);
    ROOT::TThreadedObject<TH1D> dTOFHist("dTOFHist", ";;events", 3, 0, 3);
    dTOFHist->GetXaxis()->SetBinLabel(1, "Went through");
    dTOFHist->GetXaxis()->SetBinLabel(2, "Got filtered, but not by TOF");
    dTOFHist->GetXaxis()->SetBinLabel(3, "Got filtered by TOF");

    // Define the function that will process a subrange of the tree.
    // The function must receive only one parameter, a TTreeReader,
    // and it must be thread safe. To enforce the latter requirement,
    // TThreadedObject histograms will be used.
    //but maybe later
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
        TTreeReaderValue<StUPCEvent> StUPCEventInstance(myReader, "mUPCEvent");
        TTreeReaderValue<StRPEvent> StRPEventInstance(myReader, "mRPEvent");
        //variable initialization
        myReader.Next();
        //histograms
        auto planes3Local = planes3.Get();
        auto planes4Local = planes4.Get();
        auto eventsLocal = events.Get();
        auto NpartHistLocal = NpartHist.Get();
        auto NallpartHistLocal = NallpartHist.Get();
        auto pxpy3EHistLocal = pxpy3EHist.Get();
        auto pxpy4EHistLocal = pxpy4EHist.Get();
        auto pxpy3WHistLocal = pxpy3WHist.Get();
        auto pxpy4WHistLocal = pxpy4WHist.Get();
        auto pxpxHistLocal = pxpxHist.Get();
        auto pypyHistLocal = pypyHist.Get();
        auto verZHistLocal = verZHist.Get();
        auto verZRPHistLocal = verZRPHist.Get();
        auto verZRPZTPCHistLocal = verZRPZTPCHist.Get();
        auto verDeltaZHistLocal = verDeltaZHist.Get();
        auto pxpyEHistLocal = pxpyEHist.Get();
        auto pxpyWHistLocal = pxpyWHist.Get();
        auto pxpySumEHistLocal = pxpySumEHist.Get();
        auto pxpylackHistLocal = pxpylackHist.Get();
        auto dEdxpqHistLocal = dEdxpqHist.Get();
        auto NverHistLocal = NverHist.Get();
        auto NprimverHistLocal = NprimverHist.Get();
        auto NfitHitsHistLocal = NfitHitsHist.Get();
        auto NdEdxHitsHistLocal = NdEdxHitsHist.Get();
        auto DCARHistLocal = DCARHist.Get();
        auto DCAZHistLocal = DCAZHist.Get();
        auto ptHistLocal = ptHist.Get();
        auto etaHistLocal = etaHist.Get();
        auto trigBeginningHistLocal = trigBeginningHist.Get();
        auto trigEndingHistLocal = trigEndingHist.Get();
        auto MforThesisHistLocal = MforThesisHist.Get();
        auto MbcgforThesisHistLocal = MbcgforThesisHist.Get();
        auto MforMatchHistLocal = MforMatchHist.Get();
        auto MbcgforMatchHistLocal = MbcgforMatchHist.Get();
        auto pionDCARHistLocal = pionDCARHist.Get();
        auto pionDCAZHistLocal = pionDCAZHist.Get();
        auto piondeltaDCAZHistLocal = piondeltaDCAZHist.Get();
        auto verZBadHistLocal = verZBadHist.Get();
        auto DCARBadHistLocal = DCARBadHist.Get();
        auto DCAZBadHistLocal = DCAZBadHist.Get();
        auto dTOFHistLocal = dTOFHist.Get();
        //setting up a tree & output file
        string fileName = string(myFile->GetTitle());
        fileName = fileName.substr(fileName.find_last_of("/\\")+1);
        string outfileName = outputFolder + fileName;
        cout<<"Created output file "<<outfileName<<endl;
        TFile* outputFile = TFile::Open(outfileName.c_str(), "recreate");
        TTree* mUPCTree = new TTree("mUPCTree", "mUPCTree");
        int filtered_entries = 0;
        //setting up branches
        mUPCTree->Branch("mUPCEvent", StUPCEventInstance.Get());
        mUPCTree->Branch("mRPEvent", StRPEventInstance.Get());
        //for changing branch address
        StUPCEvent* tempUPCpointer = StUPCEventInstance.Get();
        StRPEvent* tempRPpointer = StRPEventInstance.Get();

        TClonesArray* trackPointsLocal;
        StUPCRpsTrackPoint* trackPointInstanceLocal;
        TLorentzVector trackVector;
        Double_t verDeltaZ;
        vector<TLorentzVector> posPart;
        vector<TLorentzVector> negPart;
        vector<StUPCTrack*> posPion;
        vector<StUPCTrack*> negPion;
        do{
            tempUPCpointer = StUPCEventInstance.Get();
            tempRPpointer = StRPEventInstance.Get();

            eventsLocal->Fill(0.0);
            int numberOfOkayTracks = 0;
            for (Int_t i = 0; i < tempUPCpointer->getNumberOfTracks(); i++){
                if(tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kTof) && tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kPrimary)){
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
            for (Int_t i = 0; i < tempUPCpointer->getNumberOfTracks(); i++){
                if(tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kTof) && tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kPrimary)){
                    numberOfOkayTracks++;
                }
            }
            eventsLocal->Fill(1);
            NpartHistLocal->Fill(1, numberOfOkayTracks);
            NallpartHistLocal->Fill(1, tempUPCpointer->getNumberOfTracks());
            for (int i = 0; i < 6; i++){
                if(tempUPCpointer->isTrigger(CEPtriggers[i])){
                    trigBeginningHistLocal->Fill(i);
                }
            }
            
            //one proton each side
            int numberOfTracksPerSide[nSides] = {0, 0};
            for(unsigned int k = 0; k < tempRPpointer->getNumberOfTracks(); ++k)
            {
                // Get pointer to k-th track in Roman Pot data collection
                StUPCRpsTrack *trk = tempRPpointer->getTrack(k);
                trk->setEvent(tempRPpointer);
                // Get ID of a branch in which this k-th track was reconstructed
                int j = trk->branch();
                int side = j<2 ? E : W;
                numberOfTracksPerSide[side]++;
                //test to shorten the loop
                if(numberOfTracksPerSide[0]>1 || numberOfTracksPerSide[1]>1)
                    break;
            }
            if(numberOfTracksPerSide[0]!=1 || numberOfTracksPerSide[1]!=1){
                continue;
            }
            eventsLocal->Fill(2);
            numberOfOkayTracks = 0;
            for (Int_t i = 0; i < tempUPCpointer->getNumberOfTracks(); i++){
                if(tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kTof) && tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kPrimary)){
                    numberOfOkayTracks++;
                }
            }
            NpartHistLocal->Fill(2, numberOfOkayTracks);
            NallpartHistLocal->Fill(2, tempUPCpointer->getNumberOfTracks());

            //checking planes efficiency & other things with single points
            trackPointsLocal = StRPEventInstance->getTrackPoints();
            for (int i = 0; i < trackPointsLocal->GetEntries(); i++){
                trackPointInstanceLocal = (StUPCRpsTrackPoint*)trackPointsLocal->ConstructedAt(i);

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
            // trackPointsLocal->Clear();

            //checking momentum histogram for at least 3 planes
            for(unsigned int k = 0; k < tempRPpointer->getNumberOfTracks(); ++k)
            {
                // Get pointer to k-th track in Roman Pot data collection
                StUPCRpsTrack *trk = tempRPpointer->getTrack(k);
                trk->setEvent(tempRPpointer);
                //there were problems with apparently not having track point like, entirely???
                //so the first is check point if they do have them
                //and then if points are of good quality
                if(trk->getTrackPoint(0)==nullptr || trk->getTrackPoint(1)==nullptr){
                    continue;
                }
                //check if track has at least 3 of 4 RP planes used
                if(trk->getTrackPoint(0)->planesUsed()<3 || trk->getTrackPoint(1)->planesUsed()<3){
                    continue;
                }
                //East
                if(trk->branch()<2){
                    pxpy3EHistLocal->Fill(trk->pVec().X(),trk->pVec().Y());
                }
                //West
                if(trk->branch()>=2){
                    pxpy3WHistLocal->Fill(trk->pVec().X(),trk->pVec().Y());
                }
            }

            TLorentzVector RPsum = {0,0,0,0};
            //four planes on each TrackPoint
            bool areAllPlanesPresent = true;
            for(unsigned int k = 0; k < tempRPpointer->getNumberOfTracks(); ++k)
            {
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
            for (Int_t i = 0; i < tempUPCpointer->getNumberOfTracks(); i++){
                if(tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kTof) && tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kPrimary)){
                    numberOfOkayTracks++;
                }
            }
            NpartHistLocal->Fill(3, numberOfOkayTracks);
            NallpartHistLocal->Fill(3, tempUPCpointer->getNumberOfTracks());
            //checking momentum histogram for 4 planes
            for(unsigned int k = 0; k < tempRPpointer->getNumberOfTracks(); ++k)
            {
                // Get pointer to k-th track in Roman Pot data collection
                StUPCRpsTrack *trk = tempRPpointer->getTrack(k);
                trk->setEvent(tempRPpointer);
                //East
                if(trk->branch()<2){
                    pxpy4EHistLocal->Fill(trk->pVec().X(),trk->pVec().Y());
                }
                //West
                if(trk->branch()>=2){
                    pxpy4WHistLocal->Fill(trk->pVec().X(),trk->pVec().Y());
                }
                trackVector.SetVectM(trk->pVec(), particleMass[Proton]);
                RPsum += trackVector;
            }

            //exactly one primary vertex to which TOF tracks match to
            NverHistLocal->Fill(tempUPCpointer->getNumberOfVertices());
            NprimverHistLocal->Fill(tempUPCpointer->getNPrimVertices());
            vector<Int_t> primaryVertices;
            Int_t VertexId;
            for (Int_t i = 0; i < tempUPCpointer->getNumberOfTracks(); i++){
                VertexId = tempUPCpointer->getTrack(i)->getVertexId();
                if(tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kTof) && tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kPrimary) && !(find(primaryVertices.begin(), primaryVertices.end(), VertexId)!=primaryVertices.end())){
                    primaryVertices.push_back(VertexId);
                }
            }
            if(primaryVertices.size()!=1){
                continue;
            }
            VertexId = primaryVertices[0];
            eventsLocal->Fill(4);
            numberOfOkayTracks = 0;
            for (Int_t i = 0; i < tempUPCpointer->getNumberOfTracks(); i++){
                if(tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kTof) && tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kPrimary)){
                    numberOfOkayTracks++;
                }
            }
            NpartHistLocal->Fill(4, numberOfOkayTracks);
            NallpartHistLocal->Fill(4, tempUPCpointer->getNumberOfTracks());

            //primary vertex is placed within 80cm of the centre
            verZHistLocal->Fill(tempUPCpointer->getVertexId(VertexId)->getPosZ());
            if(abs(tempUPCpointer->getVertexId(VertexId)->getPosZ())>=80){
                continue;
            }
            eventsLocal->Fill(5);
            numberOfOkayTracks = 0;
            for (Int_t i = 0; i < tempUPCpointer->getNumberOfTracks(); i++){
                if(tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kTof) && tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kPrimary)){
                    numberOfOkayTracks++;
                }
            }
            NpartHistLocal->Fill(5, numberOfOkayTracks);
            NallpartHistLocal->Fill(5, tempUPCpointer->getNumberOfTracks());

            //vertexes from TPC and RP are not too far away (36cm)
            bool areBothRPTracksWithTime = true;
            Double_t time = 0;
            for(unsigned int k = 0; k < tempRPpointer->getNumberOfTracks(); ++k)
            {
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
            for (Int_t i = 0; i < tempUPCpointer->getNumberOfTracks(); i++){
                if(tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kTof) && tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kPrimary)){
                    numberOfOkayTracks++;
                }
            }
            NpartHistLocal->Fill(6, numberOfOkayTracks);
            NallpartHistLocal->Fill(6, tempUPCpointer->getNumberOfTracks());

            //if passed all tests we may copy
            mUPCTree->SetBranchAddress("mUPCEvent", &tempUPCpointer);
            mUPCTree->SetBranchAddress("mRPEvent", &tempRPpointer);
            mUPCTree->Fill();
            filtered_entries++;

            //triggers at the end
            for (int i = 0; i < 6; i++){
                if(tempUPCpointer->isTrigger(CEPtriggers[i])){
                    trigEndingHistLocal->Fill(i);
                }
            }

            //only for histograms of triggers, it won't be saved
            //quality test of central tracks
            //1. Number of hits
            //2. DCA in R and z (the second one also doubles as spread limiter)
            //3. kinematic range for optimal measurement
            int TracksValid = 0;
            TLorentzVector UPCsum = {0,0,0,0};
            for (Int_t i = 0; i < tempUPCpointer->getNumberOfTracks(); i++){
                //basic tests
                if(!tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kPrimary) || !tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kTof)){
                    continue;
                }
                if(tempUPCpointer->getTrack(i)->getTofPathLength()<=0 || tempUPCpointer->getTrack(i)->getTofTime()<=0){
                    continue;
                }
                //quality
                NfitHitsHistLocal->Fill(tempUPCpointer->getTrack(i)->getNhitsFit());
                NdEdxHitsHistLocal->Fill(tempUPCpointer->getTrack(i)->getNhitsDEdx());
                if(tempUPCpointer->getTrack(i)->getNhitsDEdx()<15 || tempUPCpointer->getTrack(i)->getNhitsFit()<25){
                    continue;
                }
                DCARHistLocal->Fill(tempUPCpointer->getTrack(i)->getDcaXY());
                DCAZHistLocal->Fill(tempUPCpointer->getTrack(i)->getDcaZ());
                if(abs(tempUPCpointer->getTrack(i)->getDcaXY())>=1.5 || abs(tempUPCpointer->getTrack(i)->getDcaZ())>=1){
                    continue;
                }
                ptHistLocal->Fill(tempUPCpointer->getTrack(i)->getPt());
                etaHistLocal->Fill(tempUPCpointer->getTrack(i)->getEta());
                if(tempUPCpointer->getTrack(i)->getPt()<=0.2 || abs(tempUPCpointer->getTrack(i)->getEta())>=0.7){
                    continue;
                }
                tempUPCpointer->getTrack(i)->getLorentzVector(trackVector, particleMass[Pion]);
                dEdxpqHistLocal->Fill(trackVector.P()*tempUPCpointer->getTrack(i)->getCharge(), tempUPCpointer->getTrack(i)->getDEdxSignal());
                if(abs(tempUPCpointer->getTrack(i)->getNSigmasTPCPion())>3){
                    continue;
                }
                //at this point track is Valid(TM)
                TracksValid++;
                UPCsum += trackVector;
                if(tempUPCpointer->getTrack(i)->getCharge()==1){
                    posPart.push_back(trackVector);
                    posPion.push_back(tempUPCpointer->getTrack(i));
                }else if (tempUPCpointer->getTrack(i)->getCharge()==-1){
                    negPart.push_back(trackVector);
                    negPion.push_back(tempUPCpointer->getTrack(i));
                }
                pionDCARHistLocal->Fill(tempUPCpointer->getTrack(i)->getDcaXY());
                pionDCAZHistLocal->Fill(tempUPCpointer->getTrack(i)->getDcaZ());
            }
            if(TracksValid==0){
                int howManyAreTOFMatched = 0;
                verZBadHistLocal->Fill(tempUPCpointer->getVertexId(VertexId)->getPosZ());
                for (Int_t i = 0; i < tempUPCpointer->getNumberOfTracks(); i++){
                    DCARBadHistLocal->Fill(tempUPCpointer->getTrack(i)->getDcaXY());
                    DCAZBadHistLocal->Fill(tempUPCpointer->getTrack(i)->getDcaZ());
                    if(tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kPrimary) && tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kTof)){
                        howManyAreTOFMatched++;
                    }
                }
                if(howManyAreTOFMatched==0){
                    dTOFHistLocal->Fill(2);
                }else{
                    dTOFHistLocal->Fill(1);
                }
                continue;
            }
            //more histograms
            dTOFHistLocal->Fill(0);
            eventsLocal->Fill(7);
            NpartHistLocal->Fill(7, TracksValid);
            NallpartHistLocal->Fill(7, tempUPCpointer->getNumberOfTracks());
            pxpxHistLocal->Fill(RPsum.X(), UPCsum.X());
            pypyHistLocal->Fill(RPsum.Y(), UPCsum.Y());
            for(unsigned int k = 0; k < tempRPpointer->getNumberOfTracks(); ++k)
            {
                // Get pointer to k-th track in Roman Pot data collection
                StUPCRpsTrack *trk = tempRPpointer->getTrack(k);
                trk->setEvent(tempRPpointer);
                //East
                if(trk->branch()<2){
                    pxpyEHistLocal->Fill(trk->pVec().X(),trk->pVec().Y());
                }
                //West
                if(trk->branch()>=2){
                    pxpyWHistLocal->Fill(trk->pVec().X(),trk->pVec().Y());
                }
            }
            pxpySumEHistLocal->Fill(RPsum.X(), RPsum.Y());
            pxpylackHistLocal->Fill(RPsum.X()+UPCsum.X(), RPsum.Y()+UPCsum.Y());
            //for mass
            TLorentzVector tempSum = {0, 0, 0, 0};
            //invariant mass
            for (long unsigned int posi = 0; posi < posPart.size(); posi++)
            {
                for (long unsigned int negi = 0; negi < negPart.size(); negi++)
                {
                    tempSum = posPart[posi] + negPart[negi];
                    MforThesisHistLocal->Fill(tempSum.M());
                    MforMatchHistLocal->Fill(tempSum.M());
                    piondeltaDCAZHistLocal->Fill(posPion[posi]->getDcaZ()-negPion[negi]->getDcaZ());
                }
            }
            //background mass (pairs of the same charge)
            for (long unsigned int posi = 0; posi < posPart.size(); posi++)
            {
                for (long unsigned int posi2 = 0; posi2 < posi; posi2++)
                {
                    tempSum = posPart[posi] + posPart[posi2];
                    MbcgforThesisHistLocal->Fill(tempSum.M());
                    MbcgforMatchHistLocal->Fill(tempSum.M());
                }
            }
            for (long unsigned int negi = 0; negi < negPart.size(); negi++)
            {
                for (long unsigned int negi2 = 0; negi2 < negi; negi2++)
                {
                    tempSum = negPart[negi] + negPart[negi2];
                    MbcgforThesisHistLocal->Fill(tempSum.M());
                    MbcgforMatchHistLocal->Fill(tempSum.M());
                }
            }

            posPart.clear();
            negPart.clear();
            posPion.clear();
            negPion.clear();
        }while (myReader.Next());

        //waiting for file opening to check if there were any filtered entries
        if(filtered_entries==0)
        {
            cout<<"Finished operation on output file "<<outfileName<<endl;
            cout<<"Analyzed "<<tempTree->GetEntries()<<" entries"<<endl;
            cout<<"There were 0 filtered entries, the file will be deleted"<<endl;

            delete tempFile;
            mUPCTree->Delete();

            outputFile->cd();
            outputFile->Close();
            gSystem->Unlink(outfileName.c_str());
            return 0;
        }

        outputFile->cd();
        mUPCTree->Write();
        outputFile->Close();

        cout<<"Finished operation on output file "<<outfileName<<endl;
        cout<<"Analyzed "<<tempTree->GetEntries()<<" entries"<<endl;
        cout<<"Filtered "<<filtered_entries<<" entries"<<endl;

        delete tempFile;

        return filtered_entries;
    };

    auto redFunction = [](const std::vector<int> &mapV)
    {
        return std::accumulate(mapV.begin(), mapV.end(), 0);
    };

    int filtered_entries = 0;
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
    filtered_entries = TreeProcessor.MapReduce(myFunction, listOfFiles, redFunction);
    // Use the TThreadedObject::Merge method to merge the thread private tree
    // into the final result

    auto planes3Final = planes3.Merge();
    auto planes4Final = planes4.Merge();
    auto eventsFinal = events.Merge();
    auto NpartHistFinal = NpartHist.Merge();
    auto NallpartHistFinal = NallpartHist.Merge();
    auto pxpy3EHistFinal = pxpy3EHist.Merge();
    auto pxpy4EHistFinal = pxpy4EHist.Merge();
    auto pxpy3WHistFinal = pxpy3WHist.Merge();
    auto pxpy4WHistFinal = pxpy4WHist.Merge();
    auto pxpxHistFinal = pxpxHist.Merge();
    auto pypyHistFinal = pypyHist.Merge();
    auto verZHistFinal = verZHist.Merge();
    auto verZRPHistFinal = verZRPHist.Merge();
    auto verZRPZTPCHistFinal = verZRPZTPCHist.Merge();
    auto verDeltaZHistFinal = verDeltaZHist.Merge();
    auto pxpyEHistFinal = pxpyEHist.Merge();
    auto pxpyWHistFinal = pxpyWHist.Merge();
    auto pxpySumEHistFinal = pxpySumEHist.Merge();
    auto pxpylackHistFinal = pxpylackHist.Merge();
    auto dEdxpqHistFinal = dEdxpqHist.Merge();
    auto NverHistFinal = NverHist.Merge();
    auto NprimverHistFinal = NprimverHist.Merge();
    auto NfitHitsHistFinal = NfitHitsHist.Merge();
    auto NdEdxHitsHistFinal = NdEdxHitsHist.Merge();
    auto DCARHistFinal = DCARHist.Merge();
    auto DCAZHistFinal = DCAZHist.Merge();
    auto ptHistFinal = ptHist.Merge();
    auto etaHistFinal = etaHist.Merge();
    auto trigBeginningHistFinal = trigBeginningHist.Merge();
    auto trigEndingHistFinal = trigEndingHist.Merge();
    auto MforThesisHistFinal = MforThesisHist.Merge();
    auto MbcgforThesisHistFinal = MbcgforThesisHist.Merge();
    auto MforMatchHistFinal = MforMatchHist.Merge();
    auto MbcgforMatchHistFinal = MbcgforMatchHist.Merge();
    auto pionDCARHistFinal = pionDCARHist.Merge();
    auto pionDCAZHistFinal = pionDCAZHist.Merge();
    auto piondeltaDCAZHistFinal = piondeltaDCAZHist.Merge();
    auto verZBadHistFinal = verZBadHist.Merge();
    auto DCARBadHistFinal = DCARBadHist.Merge();
    auto DCAZBadHistFinal = DCAZBadHist.Merge();
    auto dTOFHistFinal = dTOFHist.Merge();

    //setting up a tree & output file
    string outfileName = outputFolder + "AnaOutput_HistForOverleaf.root";
    cout<<"Created output file "<<outfileName<<endl;
    TFile* outputFileHist = TFile::Open(outfileName.c_str(), "recreate");

    outputFileHist->cd();
    planes3Final->Write();
    planes4Final->Write();
    eventsFinal->Write();
    NpartHistFinal->Write();
    NallpartHistFinal->Write();
    pxpy3EHistFinal->Write();
    pxpy4EHistFinal->Write();
    pxpy3WHistFinal->Write();
    pxpy4WHistFinal->Write();
    pxpxHistFinal->Write();
    pypyHistFinal->Write();
    verZHistFinal->Write();
    verZRPHistFinal->Write();
    verZRPZTPCHistFinal->Write();
    verDeltaZHistFinal->Write();
    pxpyEHistFinal->Write();
    pxpyWHistFinal->Write();
    pxpySumEHistFinal->Write();
    pxpylackHistFinal->Write();
    dEdxpqHistFinal->Write();
    NverHistFinal->Write();
    NprimverHistFinal->Write();
    NfitHitsHistFinal->Write();
    NdEdxHitsHistFinal->Write();
    DCARHistFinal->Write();
    DCAZHistFinal->Write();
    ptHistFinal->Write();
    etaHistFinal->Write();
    trigBeginningHistFinal->Write();
    trigEndingHistFinal->Write();
    MforThesisHistFinal->Write();
    MbcgforThesisHistFinal->Write();
    MforMatchHistFinal->Write();
    MbcgforMatchHistFinal->Write();
    pionDCARHistFinal->Write();
    pionDCAZHistFinal->Write();
    piondeltaDCAZHistFinal->Write();
    verZBadHistFinal->Write();
    DCARBadHistFinal->Write();
    DCAZBadHistFinal->Write();
    dTOFHistFinal->Write();

    outputFileHist->Close();

    cout<<"Finished processing "<<endl;
    cout<<"Analyzed total "<<upcChain->GetEntries()<<" entries"<<endl;
    cout<<"Filtered total "<<filtered_entries<<" entries"<<endl;
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