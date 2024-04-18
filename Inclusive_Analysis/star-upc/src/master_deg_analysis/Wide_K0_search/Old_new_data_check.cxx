//cpp headers

//ROOT headers
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
#include "StUPCV0.h"
#include "BeamPosition.h"

//my headers
#include "UsefulThings.h"
#include "ProcessingInsideLoop.h"
#include "ProcessingOutsideLoop.h"

enum{
    kAll = 1, kCPT, kRP, kOneVertex, kTPCTOF,
    kTotQ, kMax
};
enum SIDE{ E = 0, East = 0, W = 1, West = 1, nSides };
enum PARTICLES{ Pion = 0, Kaon = 1, Proton = 2, nParticles };
const double particleMass[nParticles] = { 0.13957, 0.497611, 0.93827 }; // pion, kaon, proton in GeV /c^2 
enum BRANCH_ID{ EU, ED, WU, WD, nBranches };
enum RP_ID{ E1U, E1D, E2U, E2D, W1U, W1D, W2U, W2D, nRomanPots };

bool CustomConnectInput(int arg, char **argv, TChain *fileChain);
// bool tuple_sort(tuple<int, int, int> t1, tuple<int, int, int> t2);

int main(int argc, char **argv){

    // int nthreads = 2;
    // if(argc==4){
    //     nthreads = atoi(argv[3]);
    // }
    // cout<<"Program is running on "<<nthreads<<" threads"<<endl;
    ROOT::EnableThreadSafety();
    //actually i'm not sure if it's needed here
    // ROOT::EnableImplicitMT(nthreads); //turn on multicore processing

    //preparing inputs
    TChain *upcChainOld = new TChain("mUPCTree");
    TChain *upcChainNew = new TChain("mUPCTree");
    if(CustomConnectInput(1, argv, upcChainOld)){
        cout<<"All old files connected"<<endl;
    }
    if(CustomConnectInput(2, argv, upcChainNew)){
        cout<<"All new files connected"<<endl;
    }

    //tables
    std::vector<tuple<int, int, int>> id_old;
    std::vector<tuple<int, int, int>> id_new;
    std::vector<tuple<int, int, int>> id_difference;

    //setting up variables
    StUPCEvent *tempUPCpointer = new StUPCEvent();
    // StRPEvent *tempRPpointer = nullptr;

    //old data
    //setting up branches
    upcChainOld->Branch("mUPCEvent", tempUPCpointer);
    // upcChainOld->Branch("mRPEvent", tempRPpointer);
    upcChainOld->SetBranchAddress("mUPCEvent", &tempUPCpointer);
    // upcChainOld->SetBranchAddress("mRPEvent", &tempRPpointer);
    //fill loop
    std::cout<<"Old data: "<<upcChainOld->GetEntries()<<std::endl;
    for(size_t i = 0; i<upcChainOld->GetEntries(); i++){
        upcChainOld->GetEntry(i);
        // if(i%100000==0){
        //     std::cout<<i<<std::endl;
        // }
        id_old.push_back(tuple(tempUPCpointer->getFillNumber(), tempUPCpointer->getRunNumber(), tempUPCpointer->getEventNumber()));
        std::cout<<tempUPCpointer->getFillNumber()<<" "<<tempUPCpointer->getRunNumber()<<" "<<tempUPCpointer->getEventNumber()<<std::endl;
    }
    //new data
    //setting up branches
    upcChainNew->Branch("mUPCEvent", tempUPCpointer);
    // upcChainNew->Branch("mRPEvent", tempRPpointer);
    upcChainNew->SetBranchAddress("mUPCEvent", &tempUPCpointer);
    // upcChainNew->SetBranchAddress("mRPEvent", &tempRPpointer);
    //fill loop
    std::cout<<"New data: "<<upcChainNew->GetEntries()<<std::endl;
    for(size_t i = 0; i<upcChainNew->GetEntries(); i++){
        upcChainNew->GetEntry(i);
        // if(i%100000==0){
        //     std::cout<<i<<std::endl;
        // }
        id_new.push_back(tuple(tempUPCpointer->getFillNumber(), tempUPCpointer->getRunNumber(), tempUPCpointer->getEventNumber()));
        std::cout<<tempUPCpointer->getFillNumber()<<" "<<tempUPCpointer->getRunNumber()<<" "<<tempUPCpointer->getEventNumber()<<std::endl;
    }

    //old - new
    // std::set_difference(id_old.begin(), id_old.end(), id_new.begin(), id_new.end(), std::back_inserter(id_difference));
    // std::cout<<"Difference old - new"<<std::endl;
    // std::cout<<id_difference.size()<<std::endl;
    // for(size_t i = 0; i<id_difference.size(); i++){
    //     std::cout<<std::get<0>(id_difference[i])<<" "<<std::get<1>(id_difference[i])<<" "<<std::get<2>(id_difference[i])<<std::endl;
    // }
    // std::sort(id_difference.begin(), id_difference.end());
    // id_difference.erase(unique(id_difference.begin(), id_difference.end()), id_difference.end());
    // std::cout<<"Difference old - new, deduplicated"<<std::endl;
    // std::cout<<id_difference.size()<<std::endl;
    // for(size_t i = 0; i<id_difference.size(); i++){
    //     std::cout<<std::get<0>(id_difference[i])<<" "<<std::get<1>(id_difference[i])<<" "<<std::get<2>(id_difference[i])<<std::endl;
    // }
    // id_difference.clear();
    // //new - old
    // std::set_difference(id_new.begin(), id_new.end(), id_old.begin(), id_old.end(), std::back_inserter(id_difference));
    // std::cout<<"Difference new - old"<<std::endl;
    // std::cout<<id_difference.size()<<std::endl;
    // for(size_t i = 0; i<id_difference.size(); i++){
    //     std::cout<<std::get<0>(id_difference[i])<<" "<<std::get<1>(id_difference[i])<<" "<<std::get<2>(id_difference[i])<<std::endl;
    // }
    // std::sort(id_difference.begin(), id_difference.end());
    // id_difference.erase(unique(id_difference.begin(), id_difference.end()), id_difference.end());
    // std::cout<<"Difference new - old, deduplicated"<<std::endl;
    // std::cout<<id_difference.size()<<std::endl;
    // for(size_t i = 0; i<id_difference.size(); i++){
    //     std::cout<<std::get<0>(id_difference[i])<<" "<<std::get<1>(id_difference[i])<<" "<<std::get<2>(id_difference[i])<<std::endl;
    // }
}

bool CustomConnectInput(int arg, char **argv, TChain *fileChain){
    int fileId = -1;
    string line;
    int lineId = 0;
    //for testing if file exists
    TFile *infile;

    const string &input = argv[arg];
    if(input.find(".list")!=string::npos){
        cout<<"Using list "<<input<<endl;
        ifstream instr(input.c_str());
        if(!instr.is_open()){
            cout<<"Couldn't open: "<<input.c_str()<<endl;
            return false;
        }


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
    } else if(input.find(".root")!=string::npos){
        cout<<"Using root file "<<input<<endl;
        fileChain->AddFile(input.c_str());
        infile = TFile::Open(input.c_str(), "read");
        if(!infile){
            cout<<"Couldn't open: "<<input.c_str()<<endl;
            return false;
        }
        infile->Close();
    }

    cout<<"There are "<<fileChain->GetEntries()<<" entries"<<endl;
    return true;
}

// bool tuple_sort(tuple<int, int, int> t1, tuple<int, int, int> t2){
//     return std::get<0>(t1)<std::get<0>(t2);
// }