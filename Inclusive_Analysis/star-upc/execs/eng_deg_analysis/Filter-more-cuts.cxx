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
void FillNewStUPCEvent(StUPCEvent* oldEvent, StUPCEvent* newEvent);
void FillNewStUPCTrack(StUPCTrack* oldTrack, StUPCTrack* newTrack);
long long GetFileSize(string filename);

//_____________________________________________________________________________
int main(int argc, char** argv) 
{
    cout<<"Program is running on 1 thread, because of incompatibility of sorts with multithreading"<<endl;
    ROOT::EnableThreadSafety();
    
    TChain* upcChain = new TChain("mUPCTree");    //chain with files to iterate through

    if(!ConnectInput(argc, argv, upcChain))
    {
        cout << "Wrong input parameters..." << endl; 
        return 1;
    }

    const string& outputFolder = argv[2];

    int filtered_entries = 0;
    vector<TFile*> listOfFiles;
    //creating the list of TFile*
    TObjArray* tempList = upcChain->GetListOfFiles();
    for (int i = 0; i < tempList->GetEntries(); i++)
    {
        listOfFiles.push_back((TFile*)(tempList->At(i)));
    }

    //for some reason declaring it in a loop or perhaps a TreeReader is a no-go
    StUPCEvent* tempUPCEvent = new StUPCEvent();

    //loop for filtering
    for(unsigned int nFile = 0; nFile < listOfFiles.size(); nFile++){
        //test if tree is not empty
        TFile* tempFile = new TFile(listOfFiles[nFile]->GetTitle());
        TTree* tempTree = (TTree*)tempFile->Get("mUPCTree");
        if(tempTree->GetEntries()==0){
            delete tempFile;
            continue;;
        }
        //creating a reader and all stuff
        TTreeReader myReader(tempTree);
        TTreeReaderValue<StUPCEvent> StUPCEventInstance(myReader, "mUPCEvent");
        TTreeReaderValue<StRPEvent> StRPEventInstance(myReader, "mRPEvent");
        //variable initialization
        myReader.Next();
        //setting up a tree & output file
        string fileName = string(listOfFiles[nFile]->GetTitle());
        fileName = fileName.substr(fileName.find_last_of("/\\")+1);
        string outfileName = outputFolder + fileName;
        cout<<"Created output file "<<outfileName<<endl;
        TFile* outputFile = TFile::Open(outfileName.c_str(), "recreate");
        TTree* mUPCTree = new TTree("mUPCTree", "mUPCTree");
        int filtered_entries_inside = 0;
        //setting up branches
        mUPCTree->Branch("mUPCEvent", StUPCEventInstance.Get());
        mUPCTree->Branch("mRPEvent", StRPEventInstance.Get());
        //for changing branch address
        StUPCEvent* tempUPCpointer = StUPCEventInstance.Get();
        StRPEvent* tempRPpointer = StRPEventInstance.Get();
        //for filtering UPC tracks

        do{
            tempUPCpointer = StUPCEventInstance.Get();
            tempRPpointer = StRPEventInstance.Get();
            //exactly one primary vertex to which TOF tracks match to
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
            //primary vertex is placed within 80cm of the centre
            if(abs(tempUPCpointer->getVertexId(VertexId)->getPosZ())>=80){
                continue;
            }
            //additional check if vertexes from TPC and RP are not too far away (36cm)
            int numberOfTracksPerSide[nSides] = {0, 0};
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
                    continue;
                }
                int side = j<2 ? E : W;
                numberOfTracksPerSide[side]++;
                //to not make the loop go twice, also west side is the "+" one
                if(side==E){
                    time -= trk->time();
                }
                if(side==W){
                    time += trk->time();
                }
            }            
            //test if there's exactly one track and if the vertex positions are not too far (36cm)
            //*100 at the end it to convert from m to cm
            if(numberOfTracksPerSide[0]!=1 || numberOfTracksPerSide[1]!=1 || abs(tempUPCpointer->getVertexId(VertexId)->getPosZ()-299792458.0*time/2*100)>36){
                continue;
            }

            //we create dummy StUPCEvent to copy only needed tracks and preinitialize it
            tempUPCEvent->clearEvent();
            FillNewStUPCEvent(tempUPCpointer,tempUPCEvent);
            //some variables that don't need to be initialised every loop
            //quality criteria of all UPC tracks
            //1. Number of hits
            //2. DCA in R and z (the second one also doubles as spread limiter)
            //3. kinematic range for optimal measurement
            //we copy only those worthy
            Int_t leftUPCTracks = 0;
            for (Int_t i = 0; i < tempUPCpointer->getNumberOfTracks(); i++){
                if(tempUPCpointer->getTrack(i)->getNhitsDEdx()<15 || tempUPCpointer->getTrack(i)->getNhitsFit()<25){
                    continue;
                }
                if(abs(tempUPCpointer->getTrack(i)->getDcaXY())>=1.5 || abs(tempUPCpointer->getTrack(i)->getDcaZ())>=1){
                    continue;
                }
                if(tempUPCpointer->getTrack(i)->getPt()<=0.2 || abs(tempUPCpointer->getTrack(i)->getEta())>=0.7){
                    continue;
                }
                FillNewStUPCTrack(tempUPCpointer->getTrack(i), tempUPCEvent->addTrack());
                leftUPCTracks++;
            }
            //if there is nothing left, we do not pass the event
            if(leftUPCTracks==0){
                continue;
            }
            
            //if all tests passed then you may copy
            mUPCTree->SetBranchAddress("mUPCEvent", &tempUPCEvent);
            mUPCTree->SetBranchAddress("mRPEvent", &tempRPpointer);
            mUPCTree->Fill();
            filtered_entries_inside++;
        }while (myReader.Next());

        //waiting for file opening to check if there were any filtered entries
        if(filtered_entries_inside==0)
        {
            cout<<"Finished operation on output file "<<outfileName<<endl;
            cout<<"Analyzed "<<tempTree->GetEntries()<<" entries"<<endl;
            cout<<"There were 0 filtered entries, the file will be deleted"<<endl;

            delete tempFile;
            mUPCTree->Delete();

            outputFile->cd();
            outputFile->Close();
            gSystem->Unlink(outfileName.c_str());
            continue;
        }

        outputFile->cd();
        mUPCTree->Write();
        outputFile->Close();

        cout<<"Finished operation on output file "<<outfileName<<endl;
        cout<<"Analyzed "<<tempTree->GetEntries()<<" entries"<<endl;
        cout<<"Filtered "<<filtered_entries_inside<<" entries"<<endl;

        delete tempFile;

        filtered_entries += filtered_entries_inside;
    };

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

void FillNewStUPCEvent(StUPCEvent* oldEvent, StUPCEvent* newEvent){
    //6 because that is no. of triggers
    for (int i = 0; i < 6; i++){
        if(newEvent->isTrigger(CEPtriggers[i]))
            newEvent->addTriggerId(CEPtriggers[i]);
    }   
    newEvent->setRunNumber(oldEvent->getRunNumber());
    newEvent->setEventNumber(oldEvent->getEventNumber());
    newEvent->setFillNumber(oldEvent->getFillNumber());
    newEvent->setNGlobTracks(oldEvent->getNGlobTracks());
    newEvent->setNPrimTracks(oldEvent->getNPrimTracks());
    newEvent->setNPrimVertices(oldEvent->getNPrimVertices());
    // all the arrays are probably unnecessary
    //maybe the TOF one isn't, but it won't matter anymore anyway
}

void FillNewStUPCTrack(StUPCTrack* oldTrack, StUPCTrack* newTrack){
    for (size_t i = 0; i < 4; i++)
    {
        if(oldTrack->getFlag(StUPCTrack::Flag(i)))
            newTrack->setFlag(StUPCTrack::Flag(i));
    }
    newTrack->setPtEtaPhi(oldTrack->getPt(), oldTrack->getEta(), oldTrack->getPhi());
    newTrack->setCurvatureDipAnglePhase(oldTrack->getCurvature(), oldTrack->getDipAngle(), oldTrack->getPhase());
    newTrack->setOrigin(oldTrack->getOrigin());
    newTrack->setDcaXY(oldTrack->getDcaXY());
    newTrack->setDcaZ(oldTrack->getDcaZ());
    newTrack->setCharge(oldTrack->getCharge());
    newTrack->setNhits(oldTrack->getNhits());
    newTrack->setNhitsFit(oldTrack->getNhitsFit());
    newTrack->setChi2(oldTrack->getChi2());
    newTrack->setNhitsDEdx(oldTrack->getNhitsDEdx());
    newTrack->setDEdxSignal(oldTrack->getDEdxSignal());
    for (size_t i = 0; i < 4; i++)
    {
        newTrack->setNSigmasTPC(StUPCTrack::Part(i), oldTrack->getNSigmasTPC(StUPCTrack::Part(i)));
    }
    newTrack->setTofTime(oldTrack->getTofTime());
    newTrack->setTofPathLength(oldTrack->getTofPathLength());
    newTrack->setVertexId(oldTrack->getVertexId());
    newTrack->setEvent(oldTrack->getEvent());

}

long long GetFileSize(string filename){
    struct stat64 stat_buf;
    int rc = stat64(filename.c_str(), &stat_buf);
    return rc==0 ? stat_buf.st_size : -1;
}