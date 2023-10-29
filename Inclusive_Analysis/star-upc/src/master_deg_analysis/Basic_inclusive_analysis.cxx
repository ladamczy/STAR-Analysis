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

//my headers
#include "UsefulThings.h"
#include "ProcessingInsideLoop.h"
#include "ProcessingOutsideLoop.h"

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
    outsideprocessing.AddHistogram(TH1D("MpipiforMatchHist", "K^{0}_{S} mass;m_{#pi^{+}#pi^{-}} [GeV];Number of pairs", 100, 0.42, 0.56));
    int triggers[] = { 570701, 570705, 570711, 590701, 590705, 590708 };
    outsideprocessing.AddHistogram(TH1D("triggerHist", "Data triggers;Trigger ID;Number of events", 6, 0, 6));
    for(int i = 0;i<6;i++){
        outsideprocessing.GetPointer1D(1)->GetXaxis()->SetBinLabel(i+1, to_string(triggers[i]).c_str());
    }


    //processing
    //defining TreeProcessor
    ROOT::TTreeProcessorMT TreeProc(*upcChain, nthreads);

    //defining processing function
    auto myFunction = [&](TTreeReader &myReader){
        //getting values from TChain, in-loop histogram initialization
        TTreeReaderValue<StUPCEvent> StUPCEventInstance(myReader, "mUPCEvent");
        TTreeReaderValue<StRPEvent> StRPEventInstance(myReader, "mRPEvent");
        ProcessingInsideLoop insideprocessing;
        StUPCEvent *localUPCEvent;
        StRPEvent *localRPEvent;
        insideprocessing.GetLocalHistograms(&outsideprocessing);

        //actual loop
        while(myReader.Next()){
            //in a TTree, it *would* be constant, in TChain however not necessarily
            localUPCEvent = StUPCEventInstance.Get();
            localRPEvent = StRPEventInstance.Get();

            //cuts & histogram filling
            //used triggers histograms
            for(int i = 0;i<6;i++){
                if(localUPCEvent->isTrigger(triggers[i])){
                    insideprocessing.Fill(1, i+1);
                }
            }
        }
        return 0;
        };

    TreeProc.Process(myFunction);

    outsideprocessing.Merge();

    //setting up a tree & output file
    string path = string(argv[0]);
    string outfileName = outputFolder+"AnaOutput_"+path.substr(path.find_last_of("/\\")+1)+".root";
    cout<<"Created output file "<<outfileName<<endl;
    TFile *outputFileHist = TFile::Open(outfileName.c_str(), "recreate");
    outsideprocessing.SaveToFile(outputFileHist);
    outputFileHist->Close();

    return 0;
}
