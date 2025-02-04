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
const double particleMass[nParticles] = { 0.13957, 0.493677, 0.93827 }; // pion, kaon, proton in GeV /c^2 
enum BRANCH_ID{ EU, ED, WU, WD, nBranches };
enum RP_ID{ E1U, E1D, E2U, E2D, W1U, W1D, W2U, W2D, nRomanPots };

//creating custom structure
//for storing the data
struct DataHolder{
    int fill, run, event;
    std::vector<double> RP_px, RP_py, TPC_pt, TPC_eta, TPC_phi;
    bool nonPrimaryWithNoAdditionalFlags;

    DataHolder(int fill, int run, int event, std::vector<double> RP_px, std::vector<double> RP_py,
        std::vector<double>TPC_pt, std::vector<double>TPC_eta, std::vector<double>TPC_phi):
        fill(fill), run(run), event(event), RP_px(RP_px), RP_py(RP_py), TPC_pt(TPC_pt),
        TPC_eta(TPC_eta), TPC_phi(TPC_phi){}

    bool operator==(const DataHolder& rhs){
        return std::tie(fill, run, event, RP_px, RP_py, TPC_pt, TPC_eta, TPC_phi)==std::tie(rhs.fill, rhs.run, rhs.event, rhs.RP_px, rhs.RP_py, rhs.TPC_pt, rhs.TPC_eta, rhs.TPC_phi);
    }
    bool operator<(const DataHolder& rhs){
        return std::tie(fill, run, event, RP_px, RP_py, TPC_pt, TPC_eta, TPC_phi)<std::tie(rhs.fill, rhs.run, rhs.event, rhs.RP_px, rhs.RP_py, rhs.TPC_pt, rhs.TPC_eta, rhs.TPC_phi);
    }
    bool operator>(const DataHolder& rhs){
        return std::tie(fill, run, event, RP_px, RP_py, TPC_pt, TPC_eta, TPC_phi)>std::tie(rhs.fill, rhs.run, rhs.event, rhs.RP_px, rhs.RP_py, rhs.TPC_pt, rhs.TPC_eta, rhs.TPC_phi);
    }
};

bool CustomConnectInput(int arg, char **argv, TChain *fileChain);
void FillingFunction(TChain* fileChain, StUPCEvent* tempUPCpointer, StRPEvent* tempRPpointer, std::vector<DataHolder>* dataStorage, std::string dataName);
void DifferentiatingFunction(std::vector<DataHolder>* dataFirst, std::vector<DataHolder>* dataSecond, std::string name);

int main(int argc, char **argv){

    int nthreads = 2;
    if(argc==4){
        nthreads = atoi(argv[3]);
    }
    cout<<"Program is running on "<<nthreads<<" threads"<<endl;
    ROOT::EnableThreadSafety();
    ROOT::EnableImplicitMT(nthreads); //turn on multicore processing

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
    std::vector<DataHolder> data_old;
    std::vector<DataHolder> data_new;

    //setting up variables
    StUPCEvent* tempUPCpointer = new StUPCEvent();
    StRPEvent* tempRPpointer = new StRPEvent();

    FillingFunction(upcChainOld, tempUPCpointer, tempRPpointer, &data_old, "Old data");
    FillingFunction(upcChainNew, tempUPCpointer, tempRPpointer, &data_new, "New data");

    std::sort(data_old.begin(), data_old.end());
    std::sort(data_new.begin(), data_new.end());

    DifferentiatingFunction(&data_old, &data_new, "old - new");
    DifferentiatingFunction(&data_new, &data_old, "new - old");

    //other check:
    //if all the tracks with no kPrimary have either kV0 or kCEP
    printf("Old data, if there exist some tracks with no kPrimary flag:\n");
    for(size_t i = 0; i<data_old.size(); i++){
        if(data_old[i].nonPrimaryWithNoAdditionalFlags){
            printf("%d %d %d\n", data_old[i].fill, data_old[i].run, data_old[i].event);
        }
    }
    printf("New data, if there exist some tracks with no kPrimary and no kV0 or kCEP flags:\n");
    for(size_t i = 0; i<data_new.size(); i++){
        if(data_new[i].nonPrimaryWithNoAdditionalFlags){
            printf("%d %d %d\n", data_new[i].fill, data_new[i].run, data_new[i].event);
        }
    }
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

void FillingFunction(TChain* fileChain, StUPCEvent* tempUPCpointer, StRPEvent* tempRPpointer, std::vector<DataHolder>* dataStorage, std::string dataName){
    //setting up branches
    fileChain->SetBranchAddress("mUPCEvent", &tempUPCpointer);
    fileChain->SetBranchAddress("mRPEvent", &tempRPpointer);
    //setting up other variables
    std::vector<double> RP_px, RP_py, TPC_pt, TPC_eta, TPC_phi;
    bool NoAdditionalFlags = false;
    //fill loop
    std::cout<<dataName+": "<<fileChain->GetEntries()<<std::endl;
    for(Long64_t i = 0; i<fileChain->GetEntries(); i++){
        fileChain->GetEntry(i);
        if(i%100000==0){
            std::cout<<i<<std::endl;
        }
        for(size_t j = 0; j<tempRPpointer->getNumberOfTracks(); j++){
            RP_px.push_back(tempRPpointer->getTrack(j)->pVec().X());
            RP_py.push_back(tempRPpointer->getTrack(j)->pVec().Y());
        }
        for(size_t j = 0; j<tempUPCpointer->getNumberOfTracks(); j++){
            if(tempUPCpointer->getTrack(j)->getFlag(StUPCTrack::kPrimary)){
                TPC_pt.push_back(tempUPCpointer->getTrack(j)->getPt());
                TPC_eta.push_back(tempUPCpointer->getTrack(j)->getEta());
                TPC_phi.push_back(tempUPCpointer->getTrack(j)->getPhi());
            } else if(!tempUPCpointer->getTrack(j)->getFlag(StUPCTrack::kV0)&&!tempUPCpointer->getTrack(j)->getFlag(StUPCTrack::kCEP)){
                NoAdditionalFlags = true;
            }
        }
        std::sort(TPC_pt.begin(), TPC_pt.end());
        std::sort(TPC_eta.begin(), TPC_eta.end());
        std::sort(TPC_phi.begin(), TPC_phi.end());

        //filling the data
        dataStorage->push_back(DataHolder(tempUPCpointer->getFillNumber(), tempUPCpointer->getRunNumber(), tempUPCpointer->getEventNumber(),
            RP_px, RP_py, TPC_pt, TPC_eta, TPC_phi));
        dataStorage->back().nonPrimaryWithNoAdditionalFlags = NoAdditionalFlags;
        //cleaning
        RP_px.clear();
        RP_py.clear();
        TPC_pt.clear();
        TPC_eta.clear();
        TPC_phi.clear();
        NoAdditionalFlags = false;
    }
}

void DifferentiatingFunction(std::vector<DataHolder>* dataFirst, std::vector<DataHolder>* dataSecond, std::string name){
    std::vector<DataHolder> data_difference;
    std::set_difference(dataFirst->begin(), dataFirst->end(), dataSecond->begin(), dataSecond->end(), std::back_inserter(data_difference));
    printf("Difference %s: %ld\n", name.c_str(), data_difference.size());
    data_difference.erase(unique(data_difference.begin(), data_difference.end()), data_difference.end());
    printf("Difference %s, deduplicated: %ld\nList of differing entries:\n", name.c_str(), data_difference.size());
    for(size_t i = 0; i<data_difference.size(); i++){
        printf("%d %d %d\n", data_difference[i].fill, data_difference[i].run, data_difference[i].event);
    }
    data_difference.clear();
}