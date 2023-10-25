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

const double particleMass[nParticles] = { 0.13957, 0.497611, 0.93827}; // pion, kaon, proton in GeV /c^2 
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
    //histogram for measuring event number
    int nEvents = 3;
    ROOT::TThreadedObject<TH1D> eventsFid("eventsFid", ";;eventsFid", nEvents, 0, nEvents);
    eventsFid->GetXaxis()->SetBinLabel(1, "Before cuts");
    eventsFid->GetXaxis()->SetBinLabel(2, "Fiducial cuts");
    eventsFid->GetXaxis()->SetBinLabel(3, "#splitline{Elimination of}{ellastic collisions}");
    //particle number
    ROOT::TThreadedObject<TH2D> NpartFidHist("NpartFidHist", ";;Number of created TOF-matched particles", nEvents, 0, nEvents, 40, 0, 40);
    NpartFidHist->GetXaxis()->SetBinLabel(1, "Before cuts");
    NpartFidHist->GetXaxis()->SetBinLabel(2, "Fiducial cuts");
    NpartFidHist->GetXaxis()->SetBinLabel(3, "#splitline{Elimination of}{ellastic collisions}");
    //particle number without TOF checking
    ROOT::TThreadedObject<TH2D> NallpartFidHist("NallpartFidHist", ";;Number of created TOF-matched particles", nEvents, 0, nEvents, 40, 0, 40);
    NallpartFidHist->GetXaxis()->SetBinLabel(1, "Before cuts");
    NallpartFidHist->GetXaxis()->SetBinLabel(2, "Fiducial cuts");
    NallpartFidHist->GetXaxis()->SetBinLabel(3, "#splitline{Elimination of}{ellastic collisions}");
    //px vs px and py vs py
    ROOT::TThreadedObject<TH2D> pxpxFidHist("pxpxFidHist", ";p_{xRP} [GeV];p_{xTPC} [GeV]", 800, -4, 4, 800, -4, 4);
    ROOT::TThreadedObject<TH2D> pypyFidHist("pypyFidHist", ";p_{yRP} [GeV];p_{yTPC} [GeV]", 800, -4, 4, 800, -4, 4);
    //final one for fiducial cuts
    ROOT::TThreadedObject<TH2D> pxpyEFidHist("pxpyEFidHist", ";p_{x} [GeV];p_{y} [GeV]", 300, -1.5, 1.5, 300, -1.5, 1.5);
    ROOT::TThreadedObject<TH2D> pxpyWFidHist("pxpyWFidHist", ";p_{x} [GeV];p_{y} [GeV]", 300, -1.5, 1.5, 300, -1.5, 1.5);
    //sum of proton momenta
    ROOT::TThreadedObject<TH2D> pxpySumEFidHist("pxpySumEFidHist", ";#Sigma p_{x} [GeV];#Sigma p_{y} [GeV]", 400, -2, 2, 400, -2, 2);
    //lacking momentum
    ROOT::TThreadedObject<TH2D> pxpylackFidHist("pxpylackFidHist", ";#Sigma p_{x} [GeV];#Sigma p_{y} [GeV]", 400, -2, 2, 400, -2, 2);
    //trigger histograms
    ROOT::TThreadedObject<TH1D> trigBeginningFidHist("trigBeginningFidHist", ";;events", 6, 0, 6);
    trigBeginningFidHist->GetXaxis()->SetBinLabel(1, "570701");
    trigBeginningFidHist->GetXaxis()->SetBinLabel(2, "570705");
    trigBeginningFidHist->GetXaxis()->SetBinLabel(3, "570711");
    trigBeginningFidHist->GetXaxis()->SetBinLabel(4, "590701");
    trigBeginningFidHist->GetXaxis()->SetBinLabel(5, "590705");
    trigBeginningFidHist->GetXaxis()->SetBinLabel(6, "590708");
    ROOT::TThreadedObject<TH1D> trigEndingFidHist("trigEndingFidHist", ";;events", 6, 0, 6);
    trigEndingFidHist->GetXaxis()->SetBinLabel(1, "570701");
    trigEndingFidHist->GetXaxis()->SetBinLabel(2, "570705");
    trigEndingFidHist->GetXaxis()->SetBinLabel(3, "570711");
    trigEndingFidHist->GetXaxis()->SetBinLabel(4, "590701");
    trigEndingFidHist->GetXaxis()->SetBinLabel(5, "590705");
    trigEndingFidHist->GetXaxis()->SetBinLabel(6, "590708");
    //mass histograms
    ROOT::TThreadedObject<TH1D> MforThesisFidHist("MforThesisFidHist", ";m_{inv} [GeV];", 300, 0, 3);
    ROOT::TThreadedObject<TH1D> MbcgforThesisFidHist("MbcgforThesisFidHist", ";m_{inv} [GeV];", 300, 0, 3);
    ROOT::TThreadedObject<TH1D> MforMatchFidHist("MforMatchFidHist", ";m_{inv} [GeV];", 100, 0.42, 0.56);
    ROOT::TThreadedObject<TH1D> MbcgforMatchFidHist("MbcgforMatchFidHist", ";m_{inv} [GeV];", 100, 0.42, 0.56);
    //Energy before cuts
    ROOT::TThreadedObject<TH1D> EbeforeFidHist("EbeforeFidHist", ";E_{p} [GeV];", 200, 100, 300);
    //Energy after fiducial cuts
    ROOT::TThreadedObject<TH1D> EafterfiducialFidHist("EafterfiducialFidHist", ";E_{p} [GeV];", 200, 100, 300);
    //Energy after both cuts
    ROOT::TThreadedObject<TH1D> EafterallFidHist("EafterallFidHist", ";E_{p} [GeV];", 200, 100, 300);

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
        auto eventsFidLocal = eventsFid.Get();
        auto NpartFidHistLocal = NpartFidHist.Get();
        auto NallpartFidHistLocal = NallpartFidHist.Get();
        auto pxpxFidHistLocal = pxpxFidHist.Get();
        auto pypyFidHistLocal = pypyFidHist.Get();
        auto pxpyEFidHistLocal = pxpyEFidHist.Get();
        auto pxpyWFidHistLocal = pxpyWFidHist.Get();
        auto pxpySumEFidHistLocal = pxpySumEFidHist.Get();
        auto pxpylackFidHistLocal = pxpylackFidHist.Get();
        auto trigBeginningFidHistLocal = trigBeginningFidHist.Get();
        auto trigEndingFidHistLocal = trigEndingFidHist.Get();
        auto MforThesisFidHistLocal = MforThesisFidHist.Get();
        auto MbcgforThesisFidHistLocal = MbcgforThesisFidHist.Get();
        auto MforMatchFidHistLocal = MforMatchFidHist.Get();
        auto MbcgforMatchFidHistLocal = MbcgforMatchFidHist.Get();
        auto EbeforeFidHistLocal = EbeforeFidHist.Get();
        auto EafterfiducialFidHistLocal = EafterfiducialFidHist.Get();
        auto EafterallFidHistLocal = EafterallFidHist.Get();
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

        TLorentzVector trackVector;
        TVector3 pVector;
        vector<TLorentzVector> posPart;
        vector<TLorentzVector> negPart;
        vector<StUPCTrack*> posPion;
        vector<StUPCTrack*> negPion;
        bool passedFiducialCuts;
        bool passedElasticCuts;
        do{
            tempUPCpointer = StUPCEventInstance.Get();
            tempRPpointer = StRPEventInstance.Get();
            passedFiducialCuts = true;
            passedElasticCuts = true;

            //triggers at the beginning
            for (int i = 0; i < 6; i++){
                if(tempUPCpointer->isTrigger(CEPtriggers[i])){
                    trigBeginningFidHistLocal->Fill(i);
                }
            }
            //TESTS & HISTOGRAMS

            //fiducial cuts
            //px barrier
            Double_t pxbarrier = -0.27;
            //py low nand high barrier
            Double_t pylowbarrier = 0.4;
            Double_t pyhighbarrier = 0.8;
            //p circle parameters
            Double_t pxcenter = -0.6;
            Double_t pycenter = 0;
            Double_t pradius = 1.1;
            for(unsigned int k = 0; k < tempRPpointer->getNumberOfTracks(); ++k){
                // Get pointer to k-th track in Roman Pot data collection
                StUPCRpsTrack *trk = tempRPpointer->getTrack(k);
                trk->setEvent(tempRPpointer);
                pVector = trk->pVec();
                bool f1 = pow(pVector.X()-pxcenter, 2) + pow(pVector.Y()-pycenter, 2) < pradius*pradius;
                bool f2 = pylowbarrier < abs(pVector.Y()) && abs(pVector.Y()) < pyhighbarrier;
                bool f3 = pxbarrier < pVector.X();
                if(!(f1 && f2 && f3)){
                    passedFiducialCuts = false;
                    break;
                }
            }

            //elastic cuts
            TLorentzVector RPsum = {0,0,0,0};
            //checking total momentum
            for(unsigned int k = 0; k < tempRPpointer->getNumberOfTracks(); ++k){
                // Get pointer to k-th track in Roman Pot data collection
                StUPCRpsTrack *trk = tempRPpointer->getTrack(k);
                trk->setEvent(tempRPpointer);
                trackVector.SetVectM(trk->pVec(), particleMass[Proton]);
                RPsum += trackVector;
            }
            //xposition of the gauss
            Double_t xGauss = -0.035;
            //yposition of the gauss
            Double_t yGauss = 0.01;
            //radius of the gauss
            //was 0.07, 0.1 FOR TEST
            Double_t rGauss = 0.07;
            if (pow(RPsum.X()-xGauss, 2) + pow(RPsum.Y()-yGauss, 2) < rGauss*rGauss){
                passedElasticCuts = false;
            }

            //if passed all tests we may copy
            if(passedFiducialCuts && passedElasticCuts){
                mUPCTree->SetBranchAddress("mUPCEvent", &tempUPCpointer);
                mUPCTree->SetBranchAddress("mRPEvent", &tempRPpointer);
                mUPCTree->Fill();
                filtered_entries++;
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
                //TEST NA NIE PRIMARY
                // if(!tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kPrimary) || !tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kTof)){
                //     continue;
                // }
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
                // if(abs(tempUPCpointer->getTrack(i)->getDcaXY())>=1.5 || abs(tempUPCpointer->getTrack(i)->getDcaZ())>=1){
                //     continue;
                // }
                if((pow(tempUPCpointer->getTrack(i)->getDcaXY(), 2) + pow(tempUPCpointer->getTrack(i)->getDcaZ(), 2))<1 && tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kPrimary)){
                    continue;
                }
                //lambda cut
                if((pow(tempUPCpointer->getTrack(i)->getDcaXY(), 2) + pow(tempUPCpointer->getTrack(i)->getDcaZ(), 2))>4 && tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kPrimary)){
                    continue;
                }
                if(tempUPCpointer->getTrack(i)->getPt()<=0.2 || abs(tempUPCpointer->getTrack(i)->getEta())>=0.7){
                    continue;
                }
                tempUPCpointer->getTrack(i)->getLorentzVector(trackVector, particleMass[Pion]);
                if(abs(tempUPCpointer->getTrack(i)->getNSigmasTPCPion())>3){
                    continue;
                }
                //at this point track is Valid(TM)
                TracksValid++;
                UPCsum += trackVector;
                //because structure of the program changed
                if(!(passedFiducialCuts && passedElasticCuts)){
                    continue;
                }
                if(tempUPCpointer->getTrack(i)->getCharge()==1){
                    posPart.push_back(trackVector);
                    posPion.push_back(tempUPCpointer->getTrack(i));
                }else if (tempUPCpointer->getTrack(i)->getCharge()==-1){
                    negPart.push_back(trackVector);
                    negPion.push_back(tempUPCpointer->getTrack(i));
                }
            }
            if(TracksValid==0){
                continue;
            }
            //more histograms
            //after only quality cuts
            eventsFidLocal->Fill(0.0);
            int numberOfOkayTracks = 0;
            for (Int_t i = 0; i < tempUPCpointer->getNumberOfTracks(); i++){
                if(tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kTof) && tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kPrimary)){
                    numberOfOkayTracks++;
                }
            }
            NpartFidHistLocal->Fill(0.0, numberOfOkayTracks);
            NallpartFidHistLocal->Fill(0.0, tempUPCpointer->getNumberOfTracks());
            for(unsigned int k = 0; k < tempRPpointer->getNumberOfTracks(); ++k)
            {
                // Get pointer to k-th track in Roman Pot data collection
                StUPCRpsTrack *trk = tempRPpointer->getTrack(k);
                trk->setEvent(tempRPpointer);
                trackVector.SetVectM(tempRPpointer->getTrack(k)->pVec(), particleMass[Proton]);
                EbeforeFidHistLocal->Fill(trackVector.E());
            }
            

            //after fiducial cuts
            if(!passedFiducialCuts){
                continue;
            }
            eventsFidLocal->Fill(1);
            numberOfOkayTracks = 0;
            for (Int_t i = 0; i < tempUPCpointer->getNumberOfTracks(); i++){
                if(tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kTof) && tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kPrimary)){
                    numberOfOkayTracks++;
                }
            }
            NpartFidHistLocal->Fill(1, numberOfOkayTracks);
            NallpartFidHistLocal->Fill(1, tempUPCpointer->getNumberOfTracks());
            for(unsigned int k = 0; k < tempRPpointer->getNumberOfTracks(); ++k)
            {
                // Get pointer to k-th track in Roman Pot data collection
                StUPCRpsTrack *trk = tempRPpointer->getTrack(k);
                trk->setEvent(tempRPpointer);
                trackVector.SetVectM(tempRPpointer->getTrack(k)->pVec(), particleMass[Proton]);
                EafterfiducialFidHistLocal->Fill(trackVector.E());
            }

            //after also elastic cuts
            if(!passedElasticCuts){
                continue;
            }
            //triggers at the end
            for (int i = 0; i < 6; i++){
                if(tempUPCpointer->isTrigger(CEPtriggers[i])){
                    trigEndingFidHistLocal->Fill(i);
                }
            }
            eventsFidLocal->Fill(2);
            NpartFidHistLocal->Fill(2, TracksValid);
            NallpartFidHistLocal->Fill(2, tempUPCpointer->getNumberOfTracks());
            for(unsigned int k = 0; k < tempRPpointer->getNumberOfTracks(); ++k)
            {
                // Get pointer to k-th track in Roman Pot data collection
                StUPCRpsTrack *trk = tempRPpointer->getTrack(k);
                trk->setEvent(tempRPpointer);
                trackVector.SetVectM(tempRPpointer->getTrack(k)->pVec(), particleMass[Proton]);
                EafterallFidHistLocal->Fill(trackVector.E());
            }
            
            pxpxFidHistLocal->Fill(RPsum.X(), UPCsum.X());
            pypyFidHistLocal->Fill(RPsum.Y(), UPCsum.Y());
            for(unsigned int k = 0; k < tempRPpointer->getNumberOfTracks(); ++k)
            {
                // Get pointer to k-th track in Roman Pot data collection
                StUPCRpsTrack *trk = tempRPpointer->getTrack(k);
                trk->setEvent(tempRPpointer);
                //East
                if(trk->branch()<2){
                    pxpyEFidHistLocal->Fill(trk->pVec().X(),trk->pVec().Y());
                }
                //West
                if(trk->branch()>=2){
                    pxpyWFidHistLocal->Fill(trk->pVec().X(),trk->pVec().Y());
                }
            }
            pxpySumEFidHistLocal->Fill(RPsum.X(), RPsum.Y());
            pxpylackFidHistLocal->Fill(RPsum.X()+UPCsum.X(), RPsum.Y()+UPCsum.Y());
            //for mass
            TLorentzVector tempSum = {0, 0, 0, 0};
            //invariant mass
            for (long unsigned int posi = 0; posi < posPart.size(); posi++)
            {
                for (long unsigned int negi = 0; negi < negPart.size(); negi++)
                {
                    tempSum = posPart[posi] + negPart[negi];
                    MforThesisFidHistLocal->Fill(tempSum.M());
                    MforMatchFidHistLocal->Fill(tempSum.M());
                }
            }
            //background mass (pairs of the same charge)
            for (long unsigned int posi = 0; posi < posPart.size(); posi++)
            {
                for (long unsigned int posi2 = 0; posi2 < posi; posi2++)
                {
                    tempSum = posPart[posi] + posPart[posi2];
                    MbcgforThesisFidHistLocal->Fill(tempSum.M());
                    MbcgforMatchFidHistLocal->Fill(tempSum.M());
                }
            }
            for (long unsigned int negi = 0; negi < negPart.size(); negi++)
            {
                for (long unsigned int negi2 = 0; negi2 < negi; negi2++)
                {
                    tempSum = negPart[negi] + negPart[negi2];
                    MbcgforThesisFidHistLocal->Fill(tempSum.M());
                    MbcgforMatchFidHistLocal->Fill(tempSum.M());
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

    auto eventsFidFinal = eventsFid.Merge();
    auto NpartFidHistFinal = NpartFidHist.Merge();
    auto NallpartFidHistFinal = NallpartFidHist.Merge();
    auto pxpxFidHistFinal = pxpxFidHist.Merge();
    auto pypyFidHistFinal = pypyFidHist.Merge();
    auto pxpyEFidHistFinal = pxpyEFidHist.Merge();
    auto pxpyWFidHistFinal = pxpyWFidHist.Merge();
    auto pxpySumEFidHistFinal = pxpySumEFidHist.Merge();
    auto pxpylackFidHistFinal = pxpylackFidHist.Merge();
    auto trigBeginningFidHistFinal = trigBeginningFidHist.Merge();
    auto trigEndingFidHistFinal = trigEndingFidHist.Merge();
    auto MforThesisFidHistFinal = MforThesisFidHist.Merge();
    auto MbcgforThesisFidHistFinal = MbcgforThesisFidHist.Merge();
    auto MforMatchFidHistFinal = MforMatchFidHist.Merge();
    auto MbcgforMatchFidHistFinal = MbcgforMatchFidHist.Merge();
    auto EbeforeFidHistFinal = EbeforeFidHist.Merge();
    auto EafterfiducialFidHistFinal = EafterfiducialFidHist.Merge();
    auto EafterallFidHistFinal = EafterallFidHist.Merge();

    //setting up a tree & output file
    string outfileName = outputFolder + "AnaOutput_FidHistForOverleafv2.root";
    cout<<"Created output file "<<outfileName<<endl;
    TFile* outputFileFidHist = TFile::Open(outfileName.c_str(), "recreate");

    outputFileFidHist->cd();

    eventsFidFinal->Write();
    NpartFidHistFinal->Write();
    NallpartFidHistFinal->Write();
    pxpxFidHistFinal->Write();
    pypyFidHistFinal->Write();
    pxpyEFidHistFinal->Write();
    pxpyWFidHistFinal->Write();
    pxpySumEFidHistFinal->Write();
    pxpylackFidHistFinal->Write();
    trigBeginningFidHistFinal->Write();
    trigEndingFidHistFinal->Write();
    MforThesisFidHistFinal->Write();
    MbcgforThesisFidHistFinal->Write();
    MforMatchFidHistFinal->Write();
    MbcgforMatchFidHistFinal->Write();
    EbeforeFidHistFinal->Write();
    EafterfiducialFidHistFinal->Write();
    EafterallFidHistFinal->Write();

    outputFileFidHist->Close();

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