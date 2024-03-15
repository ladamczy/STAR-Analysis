// Run by: ./Ana file.list output/folder/ no_of_cores
// e.g. ~/STAR-Analysis/Inclusive_Analysis/star-upc/build/bin/Preselection_CPT_2p_WE_3of4planes_fiducial /run/media/adam/OneTouch/starlist.list /run/media/adam/OneTouch/star_data_1st_cleaning/ 8


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

//my headers
#include "UsefulThings.h"

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

const double particleMass[nParticles] = { 0.13957, 0.497611, 0.93827 }; // pion, kaon, proton in GeV /c^2 

//_____________________________________________________________________________
int main(int argc, char **argv){
    int nthreads = 2;
    if(argc==4){
        nthreads = atoi(argv[3]);
    }
    cout<<"Program is running on "<<nthreads<<" threads"<<endl;
    ROOT::EnableThreadSafety();
    ROOT::EnableImplicitMT(nthreads); //turn on multicore processing

    TChain *upcChain = new TChain("mUPCTree");    //chain with files to iterate through

    if(!ConnectInput(argc, argv, upcChain)){
        cout<<"Wrong input parameters..."<<endl;
        return 1;
    }

    const string &outputFolder = argv[2];

    auto myFunction = [&](TFile *myFile){
        //test if tree is not empty
        TFile *tempFile = TFile::Open(myFile->GetTitle());
        TTree *tempTree = (TTree *)tempFile->Get("mUPCTree");
        string fileName = string(myFile->GetTitle());
        fileName = fileName.substr(fileName.find_last_of("/\\")+1);
        if(tempTree->GetEntries()==0){
            cout<<"Input file "<<fileName<<" has 0 entries"<<endl;
            return 0;
        }

        //creating a reader and all stuff
        TTreeReader myReader(tempTree);
        TTreeReaderValue<StUPCEvent> StUPCEventInstance(myReader, "mUPCEvent");
        TTreeReaderValue<StRPEvent> StRPEventInstance(myReader, "mRPEvent");
        //setting up a tree & output files
        string outfileNames[4];
        TFile *outputFiles[4];
        TTree *mUPCTrees[4];
        outfileNames[0] = outputFolder+"570701upto18083025/"+fileName;
        outfileNames[1] = outputFolder+"570704upto18083025/"+fileName;
        outfileNames[2] = outputFolder+"570705over18083025/"+fileName;
        outfileNames[3] = outputFolder+"570704over18083025/"+fileName;
        for(int i = 0; i<4; i++){
            outputFiles[i] = TFile::Open(outfileNames[i].c_str(), "recreate");
            mUPCTrees[i] = new TTree("mUPCTree", "mUPCTree");
        }
        int filtered_entries[4] = { 0, 0, 0, 0 };
        //value initialization
        myReader.Next();
        //for changing branch address
        StUPCEvent *tempUPCpointer = StUPCEventInstance.Get();
        StRPEvent *tempRPpointer = StRPEventInstance.Get();
        //setting up branches
        for(int i = 0; i<4; i++){
            mUPCTrees[i]->Branch("mUPCEvent", tempUPCpointer);
            mUPCTrees[i]->Branch("mRPEvent", tempRPpointer);
            mUPCTrees[i]->SetBranchAddress("mUPCEvent", &tempUPCpointer);
            mUPCTrees[i]->SetBranchAddress("mRPEvent", &tempRPpointer);
        }
        // actual copying
        do{
            //for some reason it *needs* to be here, God knows why
            tempUPCpointer = StUPCEventInstance.Get();
            tempRPpointer = StRPEventInstance.Get();

            //testing for all trigger and run number combinations
            if(tempUPCpointer->isTrigger(570701)&&tempUPCpointer->getRunNumber()<=18083025){
                mUPCTrees[0]->Fill();
                filtered_entries[0]++;
            }
            if(tempUPCpointer->isTrigger(570704)&&tempUPCpointer->getRunNumber()<=18083025){
                mUPCTrees[1]->Fill();
                filtered_entries[1]++;
            }
            if(tempUPCpointer->isTrigger(570705)&&tempUPCpointer->getRunNumber()>18083025){
                mUPCTrees[2]->Fill();
                filtered_entries[2]++;
            }
            if(tempUPCpointer->isTrigger(570704)&&tempUPCpointer->getRunNumber()>18083025){
                mUPCTrees[3]->Fill();
                filtered_entries[3]++;
            }

        } while(myReader.Next());

        //waiting for file opening to check if there were any filtered entries
        int total_filtered_entries = 0;
        for(int i = 0; i<4; i++){
            if(filtered_entries[i]==0){
                mUPCTrees[i]->Delete();
                outputFiles[i]->cd();
                outputFiles[i]->Close();
                gSystem->Unlink(outfileNames[i].c_str());
            }
        }

        for(int i = 0; i<4; i++){
            total_filtered_entries += filtered_entries[i];
            if(filtered_entries[i]!=0){
                outputFiles[i]->cd();
                mUPCTrees[i]->Write();
                outputFiles[i]->Close();
            }
        }

        if(total_filtered_entries==0){
            cout<<"Finished operation on input file "<<fileName<<" with "<<tempTree->GetEntries()<<" entries and 0 filtered entries, all resulting files will be deleted"<<endl;
            tempFile->Close();
            return 0;
        }

        cout<<"Finished operation on input file "<<fileName<<" with "<<tempTree->GetEntries()<<" entries and "<<total_filtered_entries<<" filtered entries"<<endl;
        tempFile->Close();

        return total_filtered_entries;
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

    cout<<"Finished processing "<<endl;
    cout<<"Analyzed total "<<upcChain->GetEntries()<<" entries"<<endl;
    cout<<"Filtered total "<<filtered_entries<<" entries"<<endl;
    cout<<"Ending Analysis... GOOD BYE!"<<endl;
    return 0;
}