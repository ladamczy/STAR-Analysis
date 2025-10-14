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
void DCASweep(TChain* inputChain, TH2D* outputHist, double testedDCA, int nthreads);
void VertexCheck(TChain* inputChain, TH1D* outputHistTab1D[], int nhists1D, TH2D* outputHistTab2D[], int nhists2D, double testedDCA, int nthreads);

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
    if(!ConnectInput(argc, argv, upcChain)){
        cout << "Wrong input parameters..." << endl; 
        return 1;
    }

    const string& outputFolder = argv[2];
    //KONIEC WSTĘPNYCH*****************************************************************************************

    //VertexCheck DCA>0****************************************************************************************

    //setting up a tree & output file
    string outfileName = outputFolder + "AnaOutput_VertexCheck_DCA0.root";
    cout<<"Created output file "<<outfileName<<endl;
    TFile* outputFileHist = TFile::Open(outfileName.c_str(), "recreate");
    TH1D* histtab1D[2];
    TH2D* histtab2D[4];
    histtab1D[0] = new TH1D("kPrimaryTestAfter", "ratio of primary tracks after DCA test;m_{#pi^{+}#pi^{-}} [GeV];% of tracks with kPrimary", 100, 0.42, 0.56);
    histtab1D[1] = new TH1D("kPrimaryTestAfterBackground", "ratio of primary tracks before DCA test;m_{#pi^{+}#pi^{-}} [GeV];% of tracks with kPrimary", 100, 0.42, 0.56);
    histtab2D[0] = new TH2D("ZTrackDifference", "m_{inv} vs #Delta DCA_{Z};m_{#pi^{+}#pi^{-}} [GeV];#Delta DCA_{Z} [cm]", 100, 0.42, 0.56, 50, 0, 5);
    histtab2D[1] = new TH2D("XVertexError", "m_{inv} vs Vertex position error (X axis);m_{#pi^{+}#pi^{-}} [GeV];#pos_{X} [cm]", 100, 0.42, 0.56, 50, 0, 0.05);
    histtab2D[2] = new TH2D("YVertexError", "m_{inv} vs Vertex position error (Y axis);m_{#pi^{+}#pi^{-}} [GeV];#pos_{Y} [cm]", 100, 0.42, 0.56, 50, 0, 0.05);
    histtab2D[3] = new TH2D("ZVertexError", "m_{inv} vs Vertex position error (Z axis);m_{#pi^{+}#pi^{-}} [GeV];#pos_{Z} [cm]", 100, 0.42, 0.56, 100, 0, 1);

    VertexCheck(upcChain, histtab1D, 2, histtab2D, 4, 0., nthreads);
    
    outputFileHist->cd();
    for (int i = 0; i < 4; i++){
        histtab2D[i]->Write();
    }
    histtab1D[0]->Divide(histtab1D[1]);
    histtab1D[0]->Write();
    
    outputFileHist->Close();

    cout<<"Finished processing VertexCheck"<<endl;
    cout<<"Analyzed total "<<upcChain->GetEntries()<<" entries"<<endl;    

    //VertexCheck DCA<1****************************************************************************************

    //setting up a tree & output file
    // string outfileName = outputFolder + "AnaOutput_VertexCheck_DCA1.root";
    outfileName = outputFolder + "AnaOutput_VertexCheck_DCA1.root";
    cout<<"Created output file "<<outfileName<<endl;
    // TFile* outputFileHist = TFile::Open(outfileName.c_str(), "recreate");
    outputFileHist = TFile::Open(outfileName.c_str(), "recreate");
    // TH1D* histtab1D[2];
    // TH2D* histtab2D[4];
    histtab1D[0] = new TH1D("kPrimaryTestAfter", "ratio of primary tracks after DCA test;m_{#pi^{+}#pi^{-}} [GeV];% of tracks with kPrimary", 100, 0.42, 0.56);
    histtab1D[1] = new TH1D("kPrimaryTestAfterBackground", "ratio of primary tracks before DCA test;m_{#pi^{+}#pi^{-}} [GeV];% of tracks with kPrimary", 100, 0.42, 0.56);
    histtab2D[0] = new TH2D("ZTrackDifference", "m_{inv} vs #Delta DCA_{Z};m_{#pi^{+}#pi^{-}} [GeV];#Delta DCA_{Z} [cm]", 100, 0.42, 0.56, 50, 0, 5);
    histtab2D[1] = new TH2D("XVertexError", "m_{inv} vs Vertex position error (X axis);m_{#pi^{+}#pi^{-}} [GeV];#pos_{X} [cm]", 100, 0.42, 0.56, 50, 0, 0.05);
    histtab2D[2] = new TH2D("YVertexError", "m_{inv} vs Vertex position error (Y axis);m_{#pi^{+}#pi^{-}} [GeV];#pos_{Y} [cm]", 100, 0.42, 0.56, 50, 0, 0.05);
    histtab2D[3] = new TH2D("ZVertexError", "m_{inv} vs Vertex position error (Z axis);m_{#pi^{+}#pi^{-}} [GeV];#pos_{Z} [cm]", 100, 0.42, 0.56, 100, 0, 1);

    VertexCheck(upcChain, histtab1D, 2, histtab2D, 4, -1., nthreads);
    
    outputFileHist->cd();
    for (int i = 0; i < 4; i++){
        histtab2D[i]->Write();
    }
    histtab1D[0]->Divide(histtab1D[1]);
    histtab1D[0]->Write();
    
    outputFileHist->Close();

    cout<<"Finished processing VertexCheck"<<endl;
    cout<<"Analyzed total "<<upcChain->GetEntries()<<" entries"<<endl;    

    return 0;


    //DCASweep***************************************************************************
    int bins = 100;
    double lowerEdge = 0;
    double upperEdge = 10;

    //setting up a tree & output file
    outfileName = outputFolder + "AnaOutput_DCASweep.root";
    cout<<"Created output file "<<outfileName<<endl;
    outputFileHist = TFile::Open(outfileName.c_str(), "recreate");
    TH2D* MinvDCAsweepHist = new TH2D("MinvDCAsweepHist", "m_{inv} vs DCA sweep;m_{#pi^{+}#pi^{-}} [GeV];DCA [cm]", 100, 0.42, 0.56, bins, lowerEdge, upperEdge);

    for (int i = 0; i < bins; i++){
        DCASweep(upcChain, MinvDCAsweepHist, (double)i*(upperEdge-lowerEdge)/bins+lowerEdge, nthreads);
    }
    
    outputFileHist->cd();
    MinvDCAsweepHist->Write();
    outputFileHist->Close();

    cout<<"Finished processing DCASweep"<<endl;
    cout<<"Analyzed total "<<upcChain->GetEntries()<<" entries"<<endl;
    cout<<"Ending Analysis... GOOD BYE!"<<endl;
    return 0;
}//main

bool ConnectInput(int argc, char** argv, TChain* fileChain) {
    int fileId = -1;
    string line;
    int lineId=0;

    const string& input = argv[1];
    cout<<"Using list "<<input<<endl;
    if(input.find(".list") != string::npos ){
        cout << "Input from chain" << endl;
        ifstream instr(input.c_str());
        if (!instr.is_open()){
            cout<< "Couldn't open: "<<input.c_str()<<endl;
            return false;
        }
        //for testing if file exists
        TFile* infile;

        while(getline(instr, line)){
            if(fileId==lineId || fileId== -1){
                fileChain->AddFile(line.c_str());
                infile = TFile::Open(line.c_str(), "read");
                if(!infile){
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

bool CheckTriggers(StUPCEvent* localupcEvt){

    bool CPTtrigger = false;
    for(int var = 0; var < nTriggers; ++var){
        if(localupcEvt->isTrigger(triggerID[var])){
            //Checked if it is CPT trigger
            for (int i = 0; i < *(&CEPtriggers + 1) - CEPtriggers; ++i)
                if(triggerID[var] == CEPtriggers[i])
                    CPTtrigger=true;
        }
    }
    return CPTtrigger;
}

void DCASweep(TChain* inputChain, TH2D* outputHist, double testedDCA, int nthreads){

    ROOT::EnableThreadSafety();
    ROOT::EnableImplicitMT(nthreads); //turn on multicore processing
    ROOT::EnableThreadSafety();

    ROOT::TThreadedObject<TH2D> TempHist(outputHist->GetName(), outputHist->GetTitle(), outputHist->GetNbinsX(), outputHist->GetXaxis()->GetXmin(), outputHist->GetXaxis()->GetXmax(), outputHist->GetNbinsY(), outputHist->GetYaxis()->GetXmin(), outputHist->GetYaxis()->GetXmax());

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
        // for future reference, if needed
        // TTreeReaderValue<StRPEvent> StRPEventInstance(myReader, "mRPEvent");
        //variable initialization
        myReader.Next();
        //histogram
        auto TempHistLocal = TempHist.Get();

        int filtered_entries = 0;
        //for changing branch address
        StUPCEvent* tempUPCpointer = StUPCEventInstance.Get();
        // for future reference, if needed
        // StRPEvent* tempRPpointer = StRPEventInstance.Get();

        TLorentzVector trackVector;
        vector<TLorentzVector> posPion;
        vector<TLorentzVector> negPion;
        int whatParticleIsThis;
        do{
            tempUPCpointer = StUPCEventInstance.Get();
            // for future reference, if needed
            // tempRPpointer = StRPEventInstance.Get();

            //quality test of central tracks
            //1. Number of hits
            //2. DCA in R and z (the second one also doubles as spread limiter)
            //3. kinematic range for optimal measurement
            int TracksValid = 0;
            TLorentzVector UPCsum = {0,0,0,0};
            for (Int_t i = 0; i < tempUPCpointer->getNumberOfTracks(); i++){
                //basic tests
                if(!tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kTof)){
                    continue;
                }
                if(tempUPCpointer->getTrack(i)->getTofPathLength()<=0 || tempUPCpointer->getTrack(i)->getTofTime()<=0){
                    continue;
                }
                //quality
                if(tempUPCpointer->getTrack(i)->getNhitsDEdx()<15 || tempUPCpointer->getTrack(i)->getNhitsFit()<25){
                    continue;
                }
                if((pow(tempUPCpointer->getTrack(i)->getDcaXY(), 2) + pow(tempUPCpointer->getTrack(i)->getDcaZ(), 2))<=testedDCA*testedDCA && tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kPrimary)){
                    continue;
                }
                if(tempUPCpointer->getTrack(i)->getPt()<=0.2 || abs(tempUPCpointer->getTrack(i)->getEta())>=0.7){
                    continue;
                }
                tempUPCpointer->getTrack(i)->getLorentzVector(trackVector, particleMass[Pion]);
                if(abs(tempUPCpointer->getTrack(i)->getNSigmasTPCProton())<3){
                    tempUPCpointer->getTrack(i)->getLorentzVector(trackVector, particleMass[Proton]);
                    whatParticleIsThis = Proton;
                }else if (abs(tempUPCpointer->getTrack(i)->getNSigmasTPCKaon())<3){
                    tempUPCpointer->getTrack(i)->getLorentzVector(trackVector, particleMass[Kaon]);
                    whatParticleIsThis = Kaon;
                }else if (abs(tempUPCpointer->getTrack(i)->getNSigmasTPCPion())<3){
                    tempUPCpointer->getTrack(i)->getLorentzVector(trackVector, particleMass[Pion]);
                    whatParticleIsThis = Pion;
                }else{
                    continue;
                }
                //at this point track is Valid(TM)
                TracksValid++;
                UPCsum += trackVector;
                if(tempUPCpointer->getTrack(i)->getCharge()==1){
                    switch (whatParticleIsThis){
                    case Pion:
                        posPion.push_back(trackVector);
                        break;                    
                    default:
                        break;
                    }
                }else if (tempUPCpointer->getTrack(i)->getCharge()==-1){
                    switch (whatParticleIsThis){
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

            //for mass
            TLorentzVector tempSum = {0, 0, 0, 0};
            //invariant mass, pi+ pi-
            for (long unsigned int posi = 0; posi < posPion.size(); posi++){
                for (long unsigned int negi = 0; negi < negPion.size(); negi++)
                {
                    tempSum = posPion[posi] + negPion[negi];
                    TempHistLocal->Fill(tempSum.M(), testedDCA);
                }
            }

            filtered_entries++;

            posPion.clear();
            negPion.clear();
        }while (myReader.Next());

        //waiting for file opening to check if there were any filtered entries
        if(filtered_entries==0)
        {
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

    auto redFunction = [](const std::vector<int> &mapV)
    {
        return std::accumulate(mapV.begin(), mapV.end(), 0);
    };

    vector<TFile*> listOfFiles;
    //creating the list of TFile*
    TObjArray* tempList = inputChain->GetListOfFiles();
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
    auto TempHistFinal = TempHist.Merge();

    outputHist->Add(TempHistFinal.get());
}//DCASweep

// histtab1D[0] = new TH1D("kPrimaryTestAfter", "ratio of primary tracks after DCA test;m_{#pi^{+}#pi^{-}} [GeV];% of tracks with kPrimary", 100, 0.42, 0.56);
// histtab1D[1] = new TH1D("kPrimaryTestAfterBackground", "ratio of primary tracks before DCA test;m_{#pi^{+}#pi^{-}} [GeV];% of tracks with kPrimary", 100, 0.42, 0.56);
// histtab2D[0] = new TH2D("ZTrackDifference", "m_{inv} vs #Delta DCA_{Z};m_{#pi^{+}#pi^{-}} [GeV];#Delta DCA_{Z} [cm]", 100, 0.42, 0.56, 100, 0, 10);
// histtab2D[1] = new TH2D("XVertexError", "m_{inv} vs Vertex position error (X axis);m_{#pi^{+}#pi^{-}} [GeV];#pos_{X} [cm]", 100, 0.42, 0.56, 200, -10, 10);
// histtab2D[2] = new TH2D("YVertexError", "m_{inv} vs Vertex position error (Y axis);m_{#pi^{+}#pi^{-}} [GeV];#pos_{Y} [cm]", 100, 0.42, 0.56, 200, -10, 10);
// histtab2D[3] = new TH2D("ZVertexError", "m_{inv} vs Vertex position error (Z axis);m_{#pi^{+}#pi^{-}} [GeV];#pos_{Z} [cm]", 100, 0.42, 0.56, 200, -10, 10);

void VertexCheck(TChain* inputChain, TH1D* outputHistTab1D[], int nhists1D, TH2D* outputHistTab2D[], int nhists2D, double testedDCA, int nthreads){

    ROOT::EnableThreadSafety();
    ROOT::EnableImplicitMT(nthreads); //turn on multicore processing
    ROOT::EnableThreadSafety();

    ROOT::TThreadedObject<TH1D> TempHistTab1D[nhists1D];
    ROOT::TThreadedObject<TH2D> TempHistTab2D[nhists2D];
    for (int i = 0; i < nhists1D; i++){
        new(&TempHistTab1D[i]) ROOT::TThreadedObject<TH1D>(outputHistTab1D[i]->GetName(), outputHistTab1D[i]->GetTitle(), outputHistTab1D[i]->GetNbinsX(), outputHistTab1D[i]->GetXaxis()->GetXmin(), outputHistTab1D[i]->GetXaxis()->GetXmax());
    }
    for (int i = 0; i < nhists2D; i++){
        new(&TempHistTab2D[i]) ROOT::TThreadedObject<TH2D>(outputHistTab2D[i]->GetName(), outputHistTab2D[i]->GetTitle(), outputHistTab2D[i]->GetNbinsX(), outputHistTab2D[i]->GetXaxis()->GetXmin(), outputHistTab2D[i]->GetXaxis()->GetXmax(), outputHistTab2D[i]->GetNbinsY(), outputHistTab2D[i]->GetYaxis()->GetXmin(), outputHistTab2D[i]->GetYaxis()->GetXmax());
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
        TTreeReaderValue<StUPCEvent> StUPCEventInstance(myReader, "mUPCEvent");
        // for future reference, if needed
        // TTreeReaderValue<StRPEvent> StRPEventInstance(myReader, "mRPEvent");
        //variable initialization
        myReader.Next();
        //histograms
        shared_ptr<TH1D> TempHistLocalTab1D[nhists1D];
        shared_ptr<TH2D> TempHistLocalTab2D[nhists2D];
        for (int i = 0; i < nhists1D; i++){
            TempHistLocalTab1D[i] = TempHistTab1D[i].Get();
        }
        for (int i = 0; i < nhists2D; i++){
            TempHistLocalTab2D[i] = TempHistTab2D[i].Get();
        }

        int filtered_entries = 0;
        //for changing branch address
        StUPCEvent* tempUPCpointer = StUPCEventInstance.Get();
        // for future reference, if needed
        // StRPEvent* tempRPpointer = StRPEventInstance.Get();

        TLorentzVector trackVector;
        vector<TLorentzVector> posPion;
        vector<TLorentzVector> negPion;
        vector<bool> posPionWasPrimary;
        vector<bool> negPionWasPrimary;
        vector<float> posPionDCAZ;
        vector<float> negPionDCAZ;
        vector<StUPCVertex*> posPionVertex;
        vector<StUPCVertex*> negPionVertex;
        int whatParticleIsThis;
        do{
            tempUPCpointer = StUPCEventInstance.Get();
            // for future reference, if needed
            // tempRPpointer = StRPEventInstance.Get();

            //quality test of central tracks
            //1. Number of hits
            //2. DCA in R and z (the second one also doubles as spread limiter)
            //3. kinematic range for optimal measurement
            int TracksValid = 0;
            TLorentzVector UPCsum = {0,0,0,0};
            for (Int_t i = 0; i < tempUPCpointer->getNumberOfTracks(); i++){
                //basic tests
                if(!tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kTof)){
                    continue;
                }
                if(tempUPCpointer->getTrack(i)->getTofPathLength()<=0 || tempUPCpointer->getTrack(i)->getTofTime()<=0){
                    continue;
                }
                //quality
                if(tempUPCpointer->getTrack(i)->getNhitsDEdx()<15 || tempUPCpointer->getTrack(i)->getNhitsFit()<25){
                    continue;
                }
                if(tempUPCpointer->getTrack(i)->getPt()<=0.2 || abs(tempUPCpointer->getTrack(i)->getEta())>=0.7){
                    continue;
                }
                if(testedDCA<0){
                    if((pow(tempUPCpointer->getTrack(i)->getDcaXY(), 2) + pow(tempUPCpointer->getTrack(i)->getDcaZ(), 2))>testedDCA*testedDCA && tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kPrimary)){
                        continue;
                    }                    
                }else{
                    if((pow(tempUPCpointer->getTrack(i)->getDcaXY(), 2) + pow(tempUPCpointer->getTrack(i)->getDcaZ(), 2))<=testedDCA*testedDCA && tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kPrimary)){
                        continue;
                    }
                }
                tempUPCpointer->getTrack(i)->getLorentzVector(trackVector, particleMass[Pion]);
                if(abs(tempUPCpointer->getTrack(i)->getNSigmasTPCProton())<3){
                    tempUPCpointer->getTrack(i)->getLorentzVector(trackVector, particleMass[Proton]);
                    whatParticleIsThis = Proton;
                }else if (abs(tempUPCpointer->getTrack(i)->getNSigmasTPCKaon())<3){
                    tempUPCpointer->getTrack(i)->getLorentzVector(trackVector, particleMass[Kaon]);
                    whatParticleIsThis = Kaon;
                }else if (abs(tempUPCpointer->getTrack(i)->getNSigmasTPCPion())<3){
                    tempUPCpointer->getTrack(i)->getLorentzVector(trackVector, particleMass[Pion]);
                    whatParticleIsThis = Pion;
                }else{
                    continue;
                }
                //at this point track is Valid(TM)
                TracksValid++;
                UPCsum += trackVector;
                if(tempUPCpointer->getTrack(i)->getCharge()==1){
                    switch (whatParticleIsThis){
                    case Pion:
                        posPion.push_back(trackVector);
                        posPionWasPrimary.push_back(tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kPrimary));
                        posPionDCAZ.push_back(tempUPCpointer->getTrack(i)->getDcaZ());
                        posPionVertex.push_back(tempUPCpointer->getTrack(i)->getVertex());
                        break;                    
                    default:
                        break;
                    }
                }else if (tempUPCpointer->getTrack(i)->getCharge()==-1){
                    switch (whatParticleIsThis){
                    case Pion:
                        negPion.push_back(trackVector);
                        negPionWasPrimary.push_back(tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kPrimary));
                        negPionDCAZ.push_back(tempUPCpointer->getTrack(i)->getDcaZ());
                        negPionVertex.push_back(tempUPCpointer->getTrack(i)->getVertex());
                        break;                    
                    default:
                        break;
                    }
                }
            }
            if(TracksValid==0){
                continue;
            }

            //for mass
            TLorentzVector tempSum = {0, 0, 0, 0};
            //invariant mass, pi+ pi-
            for (long unsigned int posi = 0; posi < posPion.size(); posi++){
                for (long unsigned int negi = 0; negi < negPion.size(); negi++)
                {
                    tempSum = posPion[posi] + negPion[negi];
                    if(posPionWasPrimary[posi] && negPionWasPrimary[negi]){
                        TempHistLocalTab1D[0]->Fill(tempSum.M());
                    }
                    TempHistLocalTab1D[1]->Fill(tempSum.M());
                    TempHistLocalTab2D[0]->Fill(tempSum.M(), abs(posPionDCAZ[posi]-negPionDCAZ[negi]));
                    TempHistLocalTab2D[1]->Fill(tempSum.M(), posPionVertex[posi]->getErrX());
                    TempHistLocalTab2D[1]->Fill(tempSum.M(), negPionVertex[negi]->getErrX());
                    TempHistLocalTab2D[2]->Fill(tempSum.M(), posPionVertex[posi]->getErrY());
                    TempHistLocalTab2D[2]->Fill(tempSum.M(), negPionVertex[negi]->getErrY());
                    TempHistLocalTab2D[3]->Fill(tempSum.M(), posPionVertex[posi]->getErrZ());
                    TempHistLocalTab2D[3]->Fill(tempSum.M(), negPionVertex[negi]->getErrZ());
                }
            }

            filtered_entries++;

            posPion.clear();
            negPion.clear();
            posPionWasPrimary.clear();
            negPionWasPrimary.clear();
            posPionDCAZ.clear();
            negPionDCAZ.clear();
            posPionVertex.clear();
            negPionVertex.clear();
        }while (myReader.Next());

        //waiting for file opening to check if there were any filtered entries
        if(filtered_entries==0)
        {
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

    auto redFunction = [](const std::vector<int> &mapV)
    {
        return std::accumulate(mapV.begin(), mapV.end(), 0);
    };

    vector<TFile*> listOfFiles;
    //creating the list of TFile*
    TObjArray* tempList = inputChain->GetListOfFiles();
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
    shared_ptr<TH1D> TempHistFinalTab1D[nhists1D];
    shared_ptr<TH2D> TempHistFinalTab2D[nhists2D];
    for (int i = 0; i < nhists1D; i++){
        TempHistFinalTab1D[i] = TempHistTab1D[i].Merge();
    }
    for (int i = 0; i < nhists2D; i++){
        TempHistFinalTab2D[i] = TempHistTab2D[i].Merge();
    }
    //uzupełnianie końca
    for (int i = 0; i < nhists1D; i++){
        outputHistTab1D[i]->Add(TempHistFinalTab1D[i].get());
    }
    for (int i = 0; i < nhists2D; i++){
        outputHistTab2D[i]->Add(TempHistFinalTab2D[i].get());
    }
}