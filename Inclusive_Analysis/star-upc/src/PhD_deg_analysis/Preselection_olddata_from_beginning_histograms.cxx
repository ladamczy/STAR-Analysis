//cpp headers
#include <map>

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
// const double particleMass[nParticles] = { 0.13957, 0.493677, 0.93827 }; // pion, kaon, proton in GeV /c^2 
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

bool IsInXiElasticSpot(StUPCRpsTrack*, StUPCRpsTrack*);
bool IsInMomElasticSpot(StUPCRpsTrack*, StUPCRpsTrack*);

int main(int argc, char** argv){

    int nthreads = 1;
    if(argc==4){
        nthreads = atoi(argv[3]);
    }

    cout<<"Program is running on "<<nthreads<<" threads"<<endl;
    ROOT::EnableThreadSafety();
    //actually i'm not sure if it's needed here
    // ROOT::EnableImplicitMT(nthreads); //turn on multicore processing

    //preparing input & output
    TChain* upcChain = new TChain("mUPCTree");
    if(ConnectInput(argc, argv, upcChain)){
        cout<<"All files connected"<<endl;
    }
    const string& outputFolder = argv[2];

    //histograms
    ProcessingOutsideLoop outsideprocessing;

    //RP_FIDUCIAL
    outsideprocessing.AddHistogram(TH2D("RP_FIDUCIAL_east", "RP_FIDUCIAL_east", 200, -1., 1., 300, -1.5, 1.5));
    outsideprocessing.AddHistogram(TH2D("RP_FIDUCIAL_west", "RP_FIDUCIAL_west", 200, -1., 1., 300, -1.5, 1.5));

    //processing
    //defining TreeProcessor
    ROOT::TTreeProcessorMT TreeProc(*upcChain, nthreads);

    //other things
    int eventsProcessed = 0;

    //Afterburner things
    LoadOffsetFile("STAR-Analysis/share/OffSetsCorrectionsRun17.list", mCorrection);

    //defining processing function
    auto myFunction = [&](TTreeReader& myReader){
        //getting values from TChain, in-loop histogram initialization
        TTreeReaderValue<StUPCEvent> StUPCEventInstance(myReader, "mUPCEvent");
        TTreeReaderValue<StRPEvent> StRPEventInstance(myReader, "mRPEvent");
        ProcessingInsideLoop insideprocessing;
        StUPCEvent* tempUPCpointer;
        StRPEvent* tempRPpointer;
        insideprocessing.GetLocalHistograms(&outsideprocessing);

        //helpful variables
        bool goodQuality;
        double firstBranch, secondBranch;
        bool f1, f2, f3;
        double px, py;
        bool IsFirstLoop = true;

        //actual loop
        while(myReader.Next()){
            //in a TTree, it *would* be constant, in TChain however not necessarily
            tempUPCpointer = StUPCEventInstance.Get();
            if(!IsFirstLoop){
                //final in-loop cleaning after Afterburner has to happen here cause we might not get to it later on
                delete tempRPpointer;
            }
            IsFirstLoop = false;
            tempRPpointer = new StRPEvent(*StRPEventInstance.Get());
            tempRPpointer->clearEvent();
            runAfterburner(StRPEventInstance.Get(), tempRPpointer, tempUPCpointer->getRunNumber());

            //cause I want to see what's going on
            if(eventsProcessed%10000==0){
                cout<<"Processed "<<eventsProcessed<<" events"<<endl;
            }
            eventsProcessed++;

            //loop  cleaning
            goodQuality = true;

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

                //fiducial
                //we measure them so they're off
                // px = trk->pVec().X();
                // py = trk->pVec().Y();
                // f1 = (0.4<abs(py)&&abs(py)<0.8);
                // f2 = (-0.27<px);
                // f3 = (pow(px+0.6, 2)+pow(py, 2)<1.25);
                // if(!(f1&&f2&&f3)){
                //     goodQuality = false;
                //     break;
                // }

                //Xi cut (without BBCL veto)
                if(trk->xi(beamMomentum)<=0.005 or trk->xi(beamMomentum)>=0.2){
                    goodQuality = false;
                    break;
                }
            }
            if(!goodQuality){ continue; }
            //non-elastic
            if(IsInXiElasticSpot(tempRPpointer->getTrack(0), tempRPpointer->getTrack(1)) or IsInMomElasticSpot(tempRPpointer->getTrack(0), tempRPpointer->getTrack(1))){
                continue;
            }

            //histogram filling
            if(firstBranch<1.5){
                //0th is east
                insideprocessing.Fill("RP_FIDUCIAL_east", tempRPpointer->getTrack(0)->pVec().X(), tempRPpointer->getTrack(0)->pVec().Y());
                insideprocessing.Fill("RP_FIDUCIAL_west", tempRPpointer->getTrack(1)->pVec().X(), tempRPpointer->getTrack(1)->pVec().Y());
            } else{
                insideprocessing.Fill("RP_FIDUCIAL_west", tempRPpointer->getTrack(0)->pVec().X(), tempRPpointer->getTrack(0)->pVec().Y());
                insideprocessing.Fill("RP_FIDUCIAL_east", tempRPpointer->getTrack(1)->pVec().X(), tempRPpointer->getTrack(1)->pVec().Y());
            }

            //lambda finish
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
    TFile* outputFileHist = TFile::Open(outfileName.c_str(), "recreate");
    outsideprocessing.SaveToFile(outputFileHist);
    outputFileHist->Close();

    return 0;
}

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