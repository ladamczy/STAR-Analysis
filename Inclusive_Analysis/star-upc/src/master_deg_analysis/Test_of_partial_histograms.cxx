//cpp headers
#include <algorithm>

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
#include "Afterburner.h"

//my headers
#include "UsefulThings.h"
#include "PartialHistograms.h"

enum{
    kAll = 1, kCPT, kRP, kOneVertex, kTPCTOF,
    kTotQ, kMax
};
// enum SIDE{ E = 0, East = 0, W = 1, West = 1, nSides };
// enum PARTICLES{ Pion = 0, Kaon = 1, Proton = 2, nParticles };
const double particleMass[nParticles] = { 0.13957, 0.497611, 0.93827 }; // pion, kaon, proton in GeV /c^2 
// enum BRANCH_ID{ EU, ED, WU, WD, nBranches };
// enum RP_ID{ E1U, E1D, E2U, E2D, W1U, W1D, W2U, W2D, nRomanPots };

//necesssary cuts
bool eventAndRunCut(StRPEvent *, StUPCEvent *, std::vector<int> *, std::vector<int> *);
//histogram cuts
bool oneVertex(StRPEvent *, StUPCEvent *, TH1 *, std::vector<int> *, std::vector<int> *);
bool zeroCut(StRPEvent *, StUPCEvent *, TH1 *, std::vector<int> *, std::vector<int> *);

//global variables
double kaonMassWindowPresentationLow = 0.46;
double kaonMassWindowPresentationHigh = 0.53;
//Afterburner, yet again
vector<vector<double>> beamData = ReadFillPositionData("STAR-Analysis/share/Run7PolarizationWithPosition.csv");

int main(int argc, char **argv){
    //preparing input & output
    TChain *upcChain = new TChain("mUPCTree");
    if(ConnectInput(argc, argv, upcChain)){
        cout<<"All files connected"<<endl;
    }
    const string &outputFolder = argv[2];

    TTreeReader myReader(upcChain);
    TTreeReaderValue<StUPCEvent> StUPCEventInstance(myReader, "mUPCEvent");
    TTreeReaderValue<StRPEvent> StRPEventInstance(myReader, "mRPEvent");
    StUPCEvent *tempUPCpointer;
    StRPEvent *tempRPpointer;
    PartialHistograms histObject;

    //Afterburner things
    LoadOffsetFile("STAR-Analysis/share/OffSetsCorrectionsRun17.list", mCorrection);

    //histograms
    TH1D hist1("hist1", "primary vertex check", 10, 0, 10);
    TH1D hist2("hist2", "primary vertex check", 10, 0, 10);

    //necessary cuts
    histObject.AddNecessaryCut(eventAndRunCut, true);
    // histObject.AddNecessaryCut(protonCuts, true);
    // histObject.AddNecessaryCut(OppositeCharges, true);
    //histogram cuts
    histObject.AddHistogramCut(oneVertex, &hist1, true);
    histObject.AddHistogramCut(zeroCut, &hist2, true);

    int counter = 0;
    while(myReader.Next()){
        counter++;
        if(counter%100000==0){
            cout<<"Analysed entry no "<<counter<<endl;
        }
        //in a TTree, it *would* be constant, in TChain however not necessarily
        tempUPCpointer = StUPCEventInstance.Get();
        // tempRPpointer = StRPEventInstance.Get();
        //modified for afterburner
        tempRPpointer = new StRPEvent(*StRPEventInstance.Get());
        tempRPpointer->clearEvent();
        runAfterburner(StRPEventInstance.Get(), tempRPpointer, tempUPCpointer->getRunNumber());

        //actual  processing
        histObject.SetEventPointers(tempRPpointer, tempUPCpointer);
        histObject.ProcessEvent();

        //final in-loop cleaning after Afterburner
        delete tempRPpointer;
    }

    //setting up a tree & output file
    string path = string(argv[0]);
    string outfileName = outputFolder+"AnaOutput_"+path.substr(path.find_last_of("/\\")+1)+".root";
    cout<<"Created output file "<<outfileName<<endl;
    TFile *outputFileHist = TFile::Open(outfileName.c_str(), "recreate");

    //histograms
    outputFileHist->cd();
    hist1.Write();
    hist2.Write();
    outputFileHist->Close();
    return 0;
}

//necessary cuts
bool eventAndRunCut(StRPEvent *RPEvent, StUPCEvent *UPCEvent, std::vector<int> *RPTrackIDs, std::vector<int> *UPCTrackIDs){
    if(UPCEvent->isTrigger(570701)){
        return true;
    }
    return false;
}


////////////////
//histogram cuts
////////////////



bool oneVertex(StRPEvent *RPEvent, StUPCEvent *UPCEvent, TH1 *hist, std::vector<int> *RPTrackIDs, std::vector<int> *UPCTrackIDs){
    //histogram
    if(hist!=nullptr){
        hist->Fill(UPCEvent->getNPrimVertices());
        return false;
    }
    //cut
    if(UPCEvent->getNPrimVertices()!=1){
        return false;
    }
    return true;
}

bool zeroCut(StRPEvent *RPEvent, StUPCEvent *UPCEvent, TH1 *hist, std::vector<int> *RPTrackIDs, std::vector<int> *UPCTrackIDs){
    //histogram
    if(hist!=nullptr){
        hist->Fill(UPCEvent->getNPrimVertices());
        return false;
    }
    //cut    
    return true;
}

// bool a(StRPEvent *RPEvent, StUPCEvent *UPCEvent, std::vector<int> *RPTrackIDs, std::vector<int> *UPCTrackIDs){

// }