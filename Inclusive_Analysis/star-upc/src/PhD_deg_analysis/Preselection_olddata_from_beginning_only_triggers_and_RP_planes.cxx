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
#include <Afterburner.h>

//my headers
#include "UsefulThings.h"

using namespace std;

// enums are very usefull 
enum{
    kAll = 1, kCPT, kRP, kOneVertex, kTPCTOF,
    kTotQ, kMax
};
// enum SIDE{ E = 0, East = 0, W = 1, West = 1, nSides };
// enum PARTICLES{ Pion = 0, Kaon = 1, Proton = 2, nParticles };
// enum BRANCH_ID{ EU, ED, WU, WD, nBranches };
// enum RP_ID{ E1U, E1D, E2U, E2D, W1U, W1D, W2U, W2D, nRomanPots };
// 570702 RP_UPC
// 570712 RP_UPC
// 570703 RP_SDT
// 570709 RP_ET
// 570719 RP_ET
// 570701 RP_CPT2
// 570711 RP_CPT2
// 570705 RP_CPT2noBBCL
// 570704 RP_Zerobias
// 590703 RP_SDT
// 590709 RP_ET
// 590701 RP_CPT2
// 590705 RP_CPT2noBBCL
// 590708 RP_CPTnoBBCL
// 570209 JPsi*HTTP
// 570219 JPsi*HTTP
// 570229 JPsi*HTTP


const double particleMass[nParticles] = { 0.13957, 0.493677, 0.93827 }; // pion, kaon, proton in GeV /c^2 
bool IsInXiElasticSpot(StUPCRpsTrack*, StUPCRpsTrack*);
bool IsInMomElasticSpot(StUPCRpsTrack*, StUPCRpsTrack*);

//_____________________________________________________________________________
int main(int argc, char** argv){
    int nthreads = 2;
    if(argc==4){
        nthreads = atoi(argv[3]);
    }
    std::cout<<"Program is running on "<<nthreads<<" threads"<<endl;
    ROOT::EnableThreadSafety();
    ROOT::EnableImplicitMT(nthreads); //turn on multicore processing

    TChain* upcChain = new TChain("mUPCTree");    //chain with files to iterate through

    if(!ConnectInput(argc, argv, upcChain)){
        std::cout<<"Wrong input parameters..."<<endl;
        return 1;
    }

    const string& outputFolder = argv[2];

    //Afterburner things
    LoadOffsetFile("STAR-Analysis/share/OffSetsCorrectionsRun17.list", mCorrection);

    auto myFunction = [&](TFile* myFile){
        //test if tree is not empty
        TFile* tempFile = TFile::Open(myFile->GetTitle());
        TTree* tempTree = (TTree*)tempFile->Get("mUPCTree");
        string fileName = string(myFile->GetTitle());
        fileName = fileName.substr(fileName.find_last_of("/\\")+1);
        if(tempTree->GetEntries()==0){
            std::cout<<"Input file "<<fileName<<" has 0 entries"<<endl;
            return 0;
        }

        //creating a reader and all stuff
        TTreeReader myReader(tempTree);
        TTreeReaderValue<StUPCEvent> StUPCEventInstance(myReader, "mUPCEvent");
        TTreeReaderValue<StRPEvent> StRPEventInstance(myReader, "mRPEvent");
        //setting up a tree & output file
        string outfileName = outputFolder+fileName;
        TFile* outputFile = TFile::Open(outfileName.c_str(), "recreate");
        TTree* mUPCTree = new TTree("mUPCTree", "mUPCTree");
        int filtered_entries = 0;
        //value initialization
        myReader.Next();
        //for changing branch address
        StUPCEvent* tempUPCpointer = StUPCEventInstance.Get();
        StRPEvent* tempRPpointer = StRPEventInstance.Get();
        //setting up branches
        mUPCTree->Branch("mUPCEvent", tempUPCpointer);
        mUPCTree->Branch("mRPEvent", tempRPpointer);
        mUPCTree->SetBranchAddress("mUPCEvent", &tempUPCpointer);
        mUPCTree->SetBranchAddress("mRPEvent", &tempRPpointer);

        // additional helpful variables
        bool IsFirstLoop = true;
        double firstBranch, secondBranch;
        bool goodQuality;

        // actual copying
        do{
            goodQuality = true;
            //for some reason it *needs* to be here, God knows why
            tempUPCpointer = StUPCEventInstance.Get();
            if(!IsFirstLoop){
                //final in-loop cleaning after Afterburner has to happen here cause we might not get to it later on
                delete tempRPpointer;
            }
            IsFirstLoop = false;
            tempRPpointer = new StRPEvent(*StRPEventInstance.Get());
            tempRPpointer->clearEvent();
            runAfterburner(StRPEventInstance.Get(), tempRPpointer, tempUPCpointer->getRunNumber());
            //cause it changed address
            mUPCTree->SetBranchAddress("mRPEvent", &tempRPpointer);

            //tests
            //trigger 570701 (no BBCL veto, minimal scaling)
            if(!tempUPCpointer->isTrigger(570701) or !(tempUPCpointer->getRunNumber()<=18083025)){
                continue;
            }
            //2 tracks
            if(tempRPpointer->getNumberOfTracks()!=2){
                continue;
            }
            //1 track east, 1 track west (neat trick - assigning negative to east by
            //substracting 1.5, and if after multiplying they are <0, they are from opposite sides
            firstBranch = tempRPpointer->getTrack(0)->branch();
            secondBranch = tempRPpointer->getTrack(1)->branch();
            if((firstBranch-1.5)*(secondBranch-1.5)>0){
                continue;
            }
            //at least 3 out of 4 planes on both and both should have both RPs hit
            //also fiducial
            for(unsigned int k = 0; k<tempRPpointer->getNumberOfTracks(); ++k){
                // Get pointer to k-th track in Roman Pot data collection
                StUPCRpsTrack* trk = tempRPpointer->getTrack(k);
                trk->setEvent(tempRPpointer);
                //there were problems with apparently not having track point like, entirely???
                //so the first is check point if they do have them
                //and then if points are of good quality
                if(trk->getTrackPoint(0)==nullptr||trk->getTrackPoint(1)==nullptr){
                    goodQuality = false;
                    break;
                }
                //check if track has at least 3 of 4 RP planes used
                if(trk->getTrackPoint(0)->planesUsed()<3||trk->getTrackPoint(1)->planesUsed()<3){
                    goodQuality = false;
                    break;
                }
            }
            if(!goodQuality){ continue; }
            //non-elastic
            if(IsInXiElasticSpot(tempRPpointer->getTrack(0), tempRPpointer->getTrack(1)) or IsInMomElasticSpot(tempRPpointer->getTrack(0), tempRPpointer->getTrack(1))){
                continue;
            }
            //end of tests

            //filling
            if(goodQuality){
                mUPCTree->Fill();
                filtered_entries++;
            }

        } while(myReader.Next());

        //waiting for file opening to check if there were any filtered entries
        if(filtered_entries==0){
            std::cout<<"Finished operation on output file "<<outfileName<<" with "<<tempTree->GetEntries()<<" entries and 0 filtered entries, the file will be deleted"<<endl;

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

        std::cout<<"Finished operation on output file "<<outfileName<<" with "<<tempTree->GetEntries()<<" entries and "<<filtered_entries<<" filtered entries"<<endl;
        tempFile->Close();

        return filtered_entries;
    };

    auto redFunction = [](const std::vector<int>& mapV){
        return std::accumulate(mapV.begin(), mapV.end(), 0);
    };

    int filtered_entries = 0;
    vector<TFile*> listOfFiles;
    //creating the list of TFile*
    TObjArray* tempList = upcChain->GetListOfFiles();
    for(int i = 0; i<tempList->GetEntries(); i++){
        listOfFiles.push_back((TFile*)(tempList->At(i)));
    }

    // Create a TreeProcessor: specify the file and the tree in it
    ROOT::TThreadExecutor TreeProcessor(nthreads);
    // Launch the parallel processing of the tree
    filtered_entries = TreeProcessor.MapReduce(myFunction, listOfFiles, redFunction);

    std::cout<<"Finished processing "<<endl;
    std::cout<<"Analyzed total "<<upcChain->GetEntries()<<" entries"<<endl;
    std::cout<<"Filtered total "<<filtered_entries<<" entries"<<endl;
    std::cout<<"Ending Analysis... GOOD BYE!"<<endl;
    return 0;
}

//NON MAIN
bool IsInXiElasticSpot(StUPCRpsTrack* east, StUPCRpsTrack* west){
    //after Afterburner
    double x_0 = -4.48170e-04;
    double sigma_x = 1.79095e-03;
    double y_0 = -8.04898e-04;
    double sigma_y = 2.12035e-03;
    return pow((east->xi(beamMomentum)-x_0)/sigma_x, 2)+pow((west->xi(beamMomentum)-y_0)/sigma_y, 2)<3*3;
}

bool IsInMomElasticSpot(StUPCRpsTrack* east, StUPCRpsTrack* west){
    //after Afterburner
    double x_0 = 5.06472e-03;
    double sigma_x = 3.42004e-02;
    double y_0 = 5.98219e-04;
    double sigma_y = 3.15726e-02;
    double x = east->pVec().X()+west->pVec().X();
    double y = east->pVec().Y()+west->pVec().Y();
    return pow((x-x_0)/sigma_x, 2)+pow((y-y_0)/sigma_y, 2)<3*3;
}