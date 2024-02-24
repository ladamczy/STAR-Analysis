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
#include <Afterburner.h>

//my headers
#include "UsefulThings.h"
#include "ProcessingInsideLoop.h"
#include "ProcessingOutsideLoop.h"

enum{
    kAll = 1, kCPT, kRP, kOneVertex, kTPCTOF,
    kTotQ, kMax
};
// enum SIDE{ E = 0, East = 0, W = 1, West = 1, nSides };
// enum PARTICLES{ Pion = 0, Kaon = 1, Proton = 2, nParticles };
const double particleMass[nParticles] = { 0.13957, 0.497611, 0.93827 }; // pion, kaon, proton in GeV /c^2 
// enum BRANCH_ID{ EU, ED, WU, WD, nBranches };
// enum RP_ID{ E1U, E1D, E2U, E2D, W1U, W1D, W2U, W2D, nRomanPots };

int main(int argc, char **argv){

    int nthreads = 2;
    if(argc==4){
        nthreads = atoi(argv[3]);
    }
    cout<<"Program is running on "<<nthreads<<" threads"<<endl;
    ROOT::EnableThreadSafety();
    //actually i'm not sure if it's needed here
    // ROOT::EnableImplicitMT(nthreads); //turn on multicore processing

    //preparing input & output
    TChain *upcChain = new TChain("mUPCTree");
    if(ConnectInput(argc, argv, upcChain)){
        cout<<"All files connected"<<endl;
    }
    const string &outputFolder = argv[2];


    //histograms
    ProcessingOutsideLoop outsideprocessing;
    outsideprocessing.AddHistogram(TH1D("BBCL_E_lower", "BBCL East up to run 18083025 inclusive", 65536, 0, 65536));
    outsideprocessing.AddHistogram(TH1D("BBCL_W_lower", "BBCL West up to run 18083025 inclusive", 65536, 0, 65536));
    outsideprocessing.AddHistogram(TH1D("BBCS_E_lower", "BBCS East up to run 18083025 inclusive", 65536, 0, 65536));
    outsideprocessing.AddHistogram(TH1D("BBCS_W_lower", "BBCS West up to run 18083025 inclusive", 65536, 0, 65536));
    outsideprocessing.AddHistogram(TH1D("BBCL_E_upper", "BBCL East after run 18083026 inclusive", 65536, 0, 65536));
    outsideprocessing.AddHistogram(TH1D("BBCL_W_upper", "BBCL West after run 18083026 inclusive", 65536, 0, 65536));
    outsideprocessing.AddHistogram(TH1D("BBCS_E_upper", "BBCS East after run 18083026 inclusive", 65536, 0, 65536));
    outsideprocessing.AddHistogram(TH1D("BBCS_W_upper", "BBCS West after run 18083026 inclusive", 65536, 0, 65536));

    //processing
    //defining TreeProcessor
    ROOT::TTreeProcessorMT TreeProc(*upcChain, nthreads);

    //defining processing function
    auto myFunction = [&](TTreeReader &myReader){
        //getting values from TChain, in-loop histogram initialization
        TTreeReaderValue<StUPCEvent> StUPCEventInstance(myReader, "mUPCEvent");
        TTreeReaderValue<StRPEvent> StRPEventInstance(myReader, "mRPEvent");
        ProcessingInsideLoop insideprocessing;
        StUPCEvent *tempUPCpointer;
        StRPEvent *tempRPpointer;
        insideprocessing.GetLocalHistograms(&outsideprocessing);

        //helpful variables

        //actual loop
        while(myReader.Next()){
            //in a TTree, it *would* be constant, in TChain however not necessarily
            tempUPCpointer = StUPCEventInstance.Get();
            // tempRPpointer = StRPEventInstance.Get();
            if(tempUPCpointer->getRunNumber()<=18083025){
                insideprocessing.Fill("BBCL_E_lower", tempUPCpointer->getBBCLargeEast());
                insideprocessing.Fill("BBCL_W_lower", tempUPCpointer->getBBCLargeWest());
                insideprocessing.Fill("BBCS_E_lower", tempUPCpointer->getBBCSmallEast());
                insideprocessing.Fill("BBCS_W_lower", tempUPCpointer->getBBCSmallWest());
            } else if(tempUPCpointer->getRunNumber()>18083025){
                insideprocessing.Fill("BBCL_E_upper", tempUPCpointer->getBBCLargeEast());
                insideprocessing.Fill("BBCL_W_upper", tempUPCpointer->getBBCLargeWest());
                insideprocessing.Fill("BBCS_E_upper", tempUPCpointer->getBBCSmallEast());
                insideprocessing.Fill("BBCS_W_upper", tempUPCpointer->getBBCSmallWest());
            }

        }
        return 0;
        };

    TreeProc.Process(myFunction);

    outsideprocessing.Merge();

    //setting up a tree & output file
    string path = string(argv[0]);
    string outfileName;
    if(outputFolder.find(".root")!=std::string::npos){
        outfileName = outputFolder;
    } else{
        outfileName = outputFolder+"AnaOutput_"+path.substr(path.find_last_of("/\\")+1)+".root";
    }
    cout<<"Created output file "<<outfileName<<endl;
    TFile *outputFileHist = TFile::Open(outfileName.c_str(), "recreate");
    outsideprocessing.SaveToFile(outputFileHist);
    outputFileHist->Close();

    return 0;
}