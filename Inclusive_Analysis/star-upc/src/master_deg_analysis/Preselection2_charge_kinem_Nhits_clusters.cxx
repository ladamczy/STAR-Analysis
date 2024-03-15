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
        //setting up a tree & output file
        string outfileName = outputFolder+fileName;
        TFile *outputFile = TFile::Open(outfileName.c_str(), "recreate");
        TTree *mUPCTree = new TTree("mUPCTree", "mUPCTree");
        int filtered_entries = 0;
        //value initialization
        myReader.Next();
        //for changing branch address
        StUPCEvent *tempUPCpointer = StUPCEventInstance.Get();
        StRPEvent *tempRPpointer = StRPEventInstance.Get();
        //setting up branches
        mUPCTree->Branch("mUPCEvent", tempUPCpointer);
        mUPCTree->Branch("mRPEvent", tempRPpointer);
        mUPCTree->SetBranchAddress("mUPCEvent", &tempUPCpointer);
        mUPCTree->SetBranchAddress("mRPEvent", &tempRPpointer);

        // additional helpful variables
        std::vector<StUPCTrack *> vector_Track;

        // actual copying
        do{
            //for some reason it *needs* to be here, God knows why
            tempUPCpointer = StUPCEventInstance.Get();
            tempRPpointer = StRPEventInstance.Get();

            //clearing vectors and such
            vector_Track.clear();

            //tests
            //at least one pair of opposite signs
            int totalCharge = 0;
            for(int i = 0; i<vector_Track.size(); i++){
                totalCharge += vector_Track[i]->getCharge();
            }
            if(abs(totalCharge)==vector_Track.size()){
                continue;
            }
            //kinematic range
            bool areAllTracksInRange = true;
            for(long unsigned int i = 0; i<vector_Track.size(); i++){
                if(!(abs(vector_Track[i]->getEta())<0.9&&vector_Track[i]->getPt()>0.2)){
                    areAllTracksInRange = false;
                    break;
                }
            }
            if(!areAllTracksInRange){
                continue;
            }
            //number of detection hits 
            bool areAllTracksWithEnoughHits = true;
            for(long unsigned int i = 0; i<vector_Track.size(); i++){
                if(!(vector_Track[i]->getNhitsFit()>25&&vector_Track[i]->getNhitsDEdx()>15)){
                    areAllTracksWithEnoughHits = false;
                    break;
                }
            }
            if(!areAllTracksWithEnoughHits){
                continue;
            }

            //TODO selekcja Patrycji
            // SELECTION: number of cluster 
            int isNumberOfTofClusterSmall = 0;
            Int_t nTofHits = tempUPCpointer->getNumberOfHits();
            vector <Int_t> vTray;
            vector <Int_t> vTrayUniqueVector;
            vector <Int_t> vMmodule;

            //filling with trays and modules (a pair for every hit)
            for(Int_t i = 0; i<nTofHits; i++){
                vTray.push_back(Int_t(tempUPCpointer->getHit(i)->getTray()));
                vTrayUniqueVector.push_back(Int_t(tempUPCpointer->getHit(i)->getTray()));
                vMmodule.push_back(Int_t(tempUPCpointer->getHit(i)->getModule()));
            }

            //making vTrayUniqueVector unique 
            sort(vTrayUniqueVector.begin(), vTrayUniqueVector.end());
            auto last = unique(vTrayUniqueVector.begin(), vTrayUniqueVector.end());
            vTrayUniqueVector.erase(last, vTrayUniqueVector.end());


            vector <vector <Int_t>> vModuleUniqueTray;
            //for every unique tray we add modules from hits with this unique tray
            for(long unsigned int i = 0; i<vTrayUniqueVector.size(); i++){
                vector <Int_t> vModuleUnique;
                for(int j = 0; j<nTofHits; j++){
                    if(vTrayUniqueVector[i]==vTray[j]){
                        vModuleUnique.push_back(vMmodule[j]);
                    }
                }
                vModuleUniqueTray.push_back(vModuleUnique);
                vModuleUnique.clear();
            }

            int totalCluster = 0;
            for(long unsigned int i = 0; i<vModuleUniqueTray.size(); i++){
                //from modules of unique trays we make unique modules of unique trays
                vector <Int_t> vec = vModuleUniqueTray[i];
                sort(vec.begin(), vec.end());
                auto last = unique(vec.begin(), vec.end());
                vec.erase(last, vec.end());

                //if there is one unique module for that unique tray for that hit we increase counter by 1
                if(vec.size()==1){
                    totalCluster += 1;
                }

                //if there is more than one unique module per hit per tray we
                //
                for(long unsigned int j = 0; j<vec.size()-1; j++){
                    Int_t modNum = vec[j];
                    int diff = 1;
                    int num = 0;
                    for(long unsigned int z = j+1; z<vec.size(); z++){
                        //if the module is the next module by diff
                        //you increase the difference and check the next pair
                        if(modNum+diff==vec[z]){
                            num += 0;
                            diff += 1;

                            if(z==j+1){
                                num += 1;
                            }
                        }

                        else if(j==(vec.size()-2)&&vec[vec.size()-2]+1!=vec[vec.size()-1]){
                            num += 2;
                            continue;
                        }

                        else{
                            num += 1;
                        }

                        j += (diff);
                    }
                    totalCluster += num;
                }
            }


            if(totalCluster<=9){
                isNumberOfTofClusterSmall = 1;
            }
            if(!isNumberOfTofClusterSmall){
                continue;
            }


            //end of tests

            //filling
            mUPCTree->Fill();
            filtered_entries++;
        } while(myReader.Next());

        //waiting for file opening to check if there were any filtered entries
        if(filtered_entries==0){
            cout<<"Finished operation on output file "<<outfileName<<" with "<<tempTree->GetEntries()<<" entries and 0 filtered entries, the file will be deleted"<<endl;

            tempFile->Close();
            mUPCTree->Delete();

            outputFile->cd();
            outputFile->Close();
            gSystem->Unlink(outfileName.c_str());
            return 0;
        }

        outputFile->cd();
        mUPCTree->Write();
        outputFile->Close();

        cout<<"Finished operation on output file "<<outfileName<<" with "<<tempTree->GetEntries()<<" entries and "<<filtered_entries<<" filtered entries"<<endl;
        tempFile->Close();

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

    cout<<"Finished processing "<<endl;
    cout<<"Analyzed total "<<upcChain->GetEntries()<<" entries"<<endl;
    cout<<"Filtered total "<<filtered_entries<<" entries"<<endl;
    cout<<"Ending Analysis... GOOD BYE!"<<endl;
    return 0;
}