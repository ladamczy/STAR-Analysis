#include <iostream>

//ROOT headers
#include <ROOT/TThreadedObject.hxx>
#include <TTreeReader.h>
#include <ROOT/TTreeProcessorMT.hxx>

//STAR headers
#include <StarGenEvent.h>
#include <StarGenParticle.h>

//my headers
#include <ProcessingInsideLoop.h>
#include <ProcessingOutsideLoop.h>
#include <UsefulThings.h>

int main(int argc, char* argv[]){

    //argv[1] - input file
    //argv[2] - output folder
    //argv[3] - #of cores

    int nthreads = 1;
    if(argc==4){
        nthreads = atoi(argv[3]);
    }

    cout<<"Program is running on "<<nthreads<<" threads"<<endl;
    ROOT::EnableThreadSafety();
    //actually i'm not sure if it's needed here
    // ROOT::EnableImplicitMT(nthreads); //turn on multicore processing

    //preparing input & output
    TChain* eventFiles = new TChain("genevents");
    if(ConnectInput(argc, argv, eventFiles)){
        cout<<"All files connected"<<endl;
    }
    const string& outputFolder = argv[2];

    //histograms
    ProcessingOutsideLoop outsideprocessing;
    //mixing TOF between events proved to be a failure
    //mass histograms with TOF first and mixing event pairs after
    outsideprocessing.AddHistogram(TH1D("TestEnergy", "Enegy of simulated particles;E [GeV];Number of particles", 200, 0.0, 2.0));

    //processing
    //defining TreeProcessor
    ROOT::TTreeProcessorMT TreeProc(*eventFiles, nthreads);

    //defining processing function
    auto myFunction = [&](TTreeReader& myReader){
        //getting values from TChain, in-loop histogram initialization
        TTreeReaderValue<StarGenEvent> tempStarGenEventInstance(myReader, "primaryEvent");
        ProcessingInsideLoop insideprocessing;
        StarGenEvent* tempStarGenEventpointer;
        insideprocessing.GetLocalHistograms(&outsideprocessing);

        while(myReader.Next()){
            tempStarGenEventpointer = tempStarGenEventInstance.Get();
            //below is the loop

            //test case - get energy of every particle that was left
            for(size_t particle_index = 0; particle_index<tempStarGenEventpointer->GetNumberOfParticles(); particle_index++){
                insideprocessing.Fill("TestEnergy", (*tempStarGenEventpointer)[particle_index]->GetEnergy());
            }

            //event loop finish
        }

        //lambda finish
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
        outfileName = outputFolder+"SimOutput_"+path.substr(path.find_last_of("/\\")+1)+".root";
    }
    cout<<"Created output file "<<outfileName<<endl;
    TFile* outputFileHist = TFile::Open(outfileName.c_str(), "recreate");
    outsideprocessing.SaveToFile(outputFileHist);
    outputFileHist->Close();

    return 0;
}
