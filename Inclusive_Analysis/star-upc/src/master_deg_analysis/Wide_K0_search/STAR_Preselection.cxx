//  Run by:
//  ./Ana full/path/to/file.root full/path/to/output/folder/
//  or
//  ./Ana full/path/to/file.root full/path/to/output/folder/file2.root

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

const double particleMass[nParticles] = { 0.13957, 0.497611, 0.93827 }; // pion, kaon, proton in GeV /c^2 
const double beamMomentum = 254.867;

bool hasEnding(std::string const &, std::string const &);
bool IsInXiElasticSpot(StUPCRpsTrack *, StUPCRpsTrack *);
bool IsInMomElasticSpot(StUPCRpsTrack *, StUPCRpsTrack *);

//_____________________________________________________________________________
int main(int argc, char **argv){

    const string &input = argv[1];
    const string &outputPath = argv[2];

    TFile *infile = TFile::Open(input.c_str(), "read");
    //check if file exists
    if(!infile){
        cout<<"Couldn't open: "<<input.c_str()<<endl;
        return 1;
    }
    //extract file name
    string fileName = input.substr(input.find_last_of("/\\")+1);

    //Afterburner things
    LoadOffsetFile("STAR-Analysis/share/OffSetsCorrectionsRun17.list", mCorrection);

    //creating a TTree and all stuff
    TTree *inputTree = (TTree *)infile->Get("mUPCTree");
    //setting up a tree & output file
    string outfileName;
    if(hasEnding(outputPath, "/")){
        outfileName = outputPath+fileName;
    } else if(hasEnding(outputPath, ".root")){
        outfileName = outputPath;
    }
    TFile *outputFile = TFile::Open(outfileName.c_str(), "recreate");
    TTree *mUPCTree = new TTree("mUPCTree", "mUPCTree");
    int filtered_entries = 0;

    // additional helpful variables
    StUPCEvent *tempUPCpointer = nullptr;
    StRPEvent *tempRPpointerBeforeAfterburner = nullptr;
    StRPEvent *tempRPpointer;
    // bool IsCPTAndRun;
    // bool IsCPTnoBBCLAndRun;
    bool goodQuality;
    double firstBranch, secondBranch;
    // bool f1, f2, f3;
    // double px, py;
    bool IsFirstLoop = true;

    //setting up branches
    inputTree->Branch("mUPCEvent", tempUPCpointer);
    inputTree->Branch("mRPEvent", tempRPpointerBeforeAfterburner);
    inputTree->SetBranchAddress("mUPCEvent", &tempUPCpointer);
    inputTree->SetBranchAddress("mRPEvent", &tempRPpointerBeforeAfterburner);
    mUPCTree->Branch("mUPCEvent", tempUPCpointer);
    mUPCTree->Branch("mRPEvent", tempRPpointerBeforeAfterburner);

    // actual copying
    for(Long64_t i = 0; i<inputTree->GetEntries(); i++){
        inputTree->GetEntry(i);
        if(!IsFirstLoop){
            //final in-loop cleaning after Afterburner has to happen here cause we might not get to it
            //because of multiple "continue;"
            delete tempRPpointer;
        }
        if(IsFirstLoop){
            IsFirstLoop = false;
        }
        // IsCPTAndRun = false;
        // IsCPTnoBBCLAndRun = false;
        goodQuality = true;
        //Afterburner fixes
        tempRPpointer = new StRPEvent(*tempRPpointerBeforeAfterburner);
        tempRPpointer->clearEvent();
        runAfterburner(tempRPpointerBeforeAfterburner, tempRPpointer, tempUPCpointer->getRunNumber());
        //cause it changed address
        mUPCTree->SetBranchAddress("mRPEvent", &tempRPpointer);


        //tests
        //run & trigger combination test
        if(tempUPCpointer->isTrigger(570701)&&tempUPCpointer->getRunNumber()<=18083025){
            // IsCPTAndRun = true;
        } else if(tempUPCpointer->isTrigger(570705)&&tempUPCpointer->getRunNumber()>18083025){
            // IsCPTnoBBCLAndRun = true;
        } else{
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
            StUPCRpsTrack *trk = tempRPpointer->getTrack(k);
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

            //fiducial
            // px = trk->pVec().X();
            // py = trk->pVec().Y();
            // f1 = (0.4<abs(py)&&abs(py)<0.8);
            // f2 = (-0.27<px);
            // f3 = (pow(px+0.6, 2)+pow(py, 2)<1.25);
            // if(!(f1&&f2&&f3)){
            //     goodQuality = false;
            //     break;
            // }

            //Xi cut
            // if(IsCPTAndRun){
            //     if(trk->xi(beamMomentum)<=0.005 or trk->xi(beamMomentum)>=0.2){
            //         goodQuality = false;
            //         break;
            //     }
            // } else if(IsCPTnoBBCLAndRun){
            //     if(trk->xi(beamMomentum)<=0.005 or trk->xi(beamMomentum)>=0.08){
            //         goodQuality = false;
            //         break;
            //     }
            // } else{
            //     goodQuality = false;
            //     break;
            // }
        }
        if(!goodQuality){ continue; }
        //non-elastic
        // if(IsInXiElasticSpot(tempRPpointer->getTrack(0), tempRPpointer->getTrack(1)) or IsInMomElasticSpot(tempRPpointer->getTrack(0), tempRPpointer->getTrack(1))){
        //     continue;
        // }
        //UPC tests
        //at least 2 good tracks of opposite signs
        int nOfGoodTracks = 0;
        int chargeOfGoodTracks = 0;
        for(int i = 0; i<tempUPCpointer->getNumberOfTracks(); i++){
            StUPCTrack *tmptrk = tempUPCpointer->getTrack(i);
            if(tmptrk->getFlag(StUPCTrack::kTof)&&abs(tmptrk->getEta())<1.1&&tmptrk->getPt()>0.15&&tmptrk->getNhitsFit()>15&&(!tmptrk->getFlag(StUPCTrack::kV0))&&(!tmptrk->getFlag(StUPCTrack::kCEP))){
                nOfGoodTracks++;
                chargeOfGoodTracks += tmptrk->getCharge();
            }
        }
        if(nOfGoodTracks<2 or abs(chargeOfGoodTracks)==nOfGoodTracks){
            continue;
        }

        //end of tests

        //filling
        if(goodQuality){
            mUPCTree->Fill();
            filtered_entries++;
        }
    }

    //waiting for file opening to check if there were any filtered entries
    if(filtered_entries==0){
        cout<<"Finished operation on output file "<<outfileName<<" with "<<inputTree->GetEntries()<<" entries and 0 filtered entries, the file will be deleted"<<endl;

        infile->Close();
        mUPCTree->Delete();

        outputFile->cd();
        outputFile->Close();
        gSystem->Unlink(outfileName.c_str());
        return 0;
    }

    outputFile->cd();
    mUPCTree->Write();
    outputFile->Close();

    cout<<"Finished operation on output file "<<outfileName<<" with "<<inputTree->GetEntries()<<" entries and "<<filtered_entries<<" filtered entries"<<endl;
    infile->Close();

    return 0;
}

//NON MAIN
bool hasEnding(std::string const &fullString, std::string const &ending){
    if(fullString.length()>=ending.length()){
        return (0==fullString.compare(fullString.length()-ending.length(), ending.length(), ending));
    } else{
        return false;
    }
}

// bool CheckValue(ROOT::Internal::TTreeReaderValueBase &value){
//     if(value.GetSetupStatus()<0){
//         std::cerr<<"Error "<<value.GetSetupStatus()
//             <<"setting up reader for "<<value.GetBranchName()<<'\n';
//         return false;
//     }
//     return true;
// }

bool IsInXiElasticSpot(StUPCRpsTrack *east, StUPCRpsTrack *west){
    //before Afterburner
    // double x_0 = 4.30588e-03;
    // double sigma_x = 2.02340e-03;
    // double y_0 = 1.72097e-03;
    // double sigma_y = 2.26638e-03;
    //after Afterburner
    double x_0 = -4.48170e-04;
    double sigma_x = 1.79095e-03;
    double y_0 = -8.04898e-04;
    double sigma_y = 2.12035e-03;
    return pow((east->xi(beamMomentum)-x_0)/sigma_x, 2)+pow((west->xi(beamMomentum)-y_0)/sigma_y, 2)<3*3;
}

bool IsInMomElasticSpot(StUPCRpsTrack *east, StUPCRpsTrack *west){
    //before Afterburner
    // double x_0 = -3.82151e-02;
    // double sigma_x = 3.67545e-02;
    // double y_0 = 1.98348e-03;
    // double sigma_y = 3.40440e-02;
    //after Afterburner
    double x_0 = 5.06472e-03;
    double sigma_x = 3.42004e-02;
    double y_0 = 5.98219e-04;
    double sigma_y = 3.15726e-02;
    double x = east->pVec().X()+west->pVec().X();
    double y = east->pVec().Y()+west->pVec().Y();
    return pow((x-x_0)/sigma_x, 2)+pow((y-y_0)/sigma_y, 2)<3*3;
}