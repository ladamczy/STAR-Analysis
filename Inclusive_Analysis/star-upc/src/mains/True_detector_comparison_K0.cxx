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
#include <filesystem>

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

//repeating things
#include "UsefulThings.h"

using namespace std;

// enums are very usefull 
enum { kAll = 1, kCPT,  kRP, kOneVertex, kTPCTOF, 
    kTotQ, kMax};
enum SIDE {E = 0, East = 0, W = 1, West = 1, nSides};
enum PARTICLES {Pion = 0, Kaon = 1, Proton = 2, nParticles};
enum BRANCH_ID { EU, ED, WU, WD, nBranches };
enum RP_ID {E1U, E1D, E2U, E2D, W1U, W1D, W2U, W2D, nRomanPots};

const double particleMass[nParticles] = { 0.13957, 0.497611, 0.93827}; // pion, kaon, proton in GeV /c^2 

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

    if(!ConnectInput(argc, argv, upcChain)){
        cout << "Wrong input parameters..." << endl; 
        return 1;
    }

    const string& outputFolder = argv[2];

    //HISTOGRAMS
    ROOT::TThreadedObject<TH1D> metricbetween = ROOT::TThreadedObject<TH1D>("metricbetween", "Distance between simulation and event in #phi-#eta space;;events", 18, 0, TMath::Pi());
    
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
        //unavailable
        // TTreeReaderValue<StRPEvent> StRPEventInstance(myReader, "mRPEvent");
        //variable initialization
        myReader.Next();
        //histograms
        std::shared_ptr<TH1D> metricbetweenLocal;

        metricbetweenLocal = metricbetween.Get();

        int filtered_entries = 0;
        //for changing branch address
        StUPCEvent* tempUPCpointer = StUPCEventInstance.Get();
        //unavailable
        // StRPEvent* tempRPpointer = StRPEventInstance.Get();

        int total_charged_MC_particles = 0;
        do{
            tempUPCpointer = StUPCEventInstance.Get();
            vector<TParticle*> true_particles;
            vector<StUPCTrack*> detected_particles;
            vector<double> metric_difference;
            vector<bool> was_used(tempUPCpointer->getNumberOfTracks(), false);

            //assigning and identifying
            for (int i = 0; i < tempUPCpointer->getNumberOfTracks(); i++){
                //processing
                int best_choice = -1;
                double best_metric = 100.0; //absurdly high
                for (int j = 0; j < tempUPCpointer->getNumberOfMCParticles(); j++){
                    //counting particles present and eliminating neutral ones
                    if (int(round(tempUPCpointer->getMCParticle(j)->GetPDG()->Charge()))!=0){
                        total_charged_MC_particles++;
                    }else{
                        continue;
                    }
                    //check metric and match
                    double temp_metric = MetricCheck(tempUPCpointer->getMCParticle(j), tempUPCpointer->getTrack(i));
                    if(temp_metric<0 || was_used[j]){
                        continue;
                    }else if(temp_metric<best_metric){
                        best_choice = j;
                        best_metric = temp_metric;
                        was_used[j] = true;
                    }
                }
                if(best_metric<100.0 && best_choice>=0){
                    true_particles.push_back(tempUPCpointer->getMCParticle(best_choice));
                    detected_particles.push_back(tempUPCpointer->getTrack(i));
                    metric_difference.push_back(best_metric);
                }
            }

            //things to do with identified particles
            for (int i = 0; i < metric_difference.size(); i++){
                metricbetweenLocal->Fill(metric_difference[i]);
                filtered_entries++;
            }
            
            
            true_particles.clear();
            detected_particles.clear();
            metric_difference.clear();
        }while (myReader.Next());

        //waiting for file opening to check if there were any filtered entries
        if(filtered_entries==0){
            cout<<"Finished operation on file "<<myFile->GetTitle()<<endl;
            cout<<"Analyzed "<<tempTree->GetEntries()<<" entries"<<endl;
            cout<<"There were 0 filtered entries"<<endl;

            delete tempFile;
            return 0;
        }

        cout<<"Finished operation on file "<<myFile->GetTitle()<<endl;
        cout<<"Analyzed "<<tempTree->GetEntries()<<" entries"<<endl;
        cout<<"Matched "<<filtered_entries<<" out of "<<total_charged_MC_particles<<endl;

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
    std::shared_ptr<TH1D> metricbetweenFinal;

    metricbetweenFinal = metricbetween.Merge();

    metricbetweenFinal->SetMinimum(0);

    //setting up a tree & output file
    string path = string(argv[0]);
    string outfileName = outputFolder + "AnaOutput_" + path.substr(path.find_last_of("/\\") + 1) + ".root";
    cout<<"Created output file "<<outfileName<<endl;
    TFile* outputFileHist = TFile::Open(outfileName.c_str(), "recreate");

    outputFileHist->cd();
    metricbetweenFinal->Write();
    outputFileHist->Close();

    cout<<"Finished processing "<<endl;
    cout<<"Analyzed total "<<upcChain->GetEntries()<<" entries"<<endl;
    cout<<"Filtered total "<<filtered_entries<<" entries"<<endl;
    cout<<"Ending Analysis... GOOD BYE!"<<endl;
    return 0;
}//main