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
#include <TEfficiency.h>
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
#include "StUPCV0.h"
#include "BeamPosition.h"
#include "StPicoPhysicalHelix.h"

//my headers
#include "UsefulThings.h"
#include "ProcessingInsideLoop.h"
#include "ProcessingOutsideLoop.h"

using namespace std;

// enums are very usefull 
enum{
    kAll = 1, kCPT, kRP, kOneVertex, kTPCTOF,
    kTotQ, kMax
};
// enum SIDE{ E = 0, East = 0, W = 1, West = 1, nSides };
enum PARTICLES{ Pion = 0, Kaon = 1, Proton = 2, nParticles };
// enum BRANCH_ID{ EU, ED, WU, WD, nBranches };
// enum RP_ID{ E1U, E1D, E2U, E2D, W1U, W1D, W2U, W2D, nRomanPots };
// enum ERRORS{ OK, RP_PLANES, RP_FIDUCIAL, RP_XI, RP_ELASTIC, TRACKS_TOF, TRACKS_PT, TRACKS_ETA, TRACKS_NHITS, PAIRS_DCADAUGHTERS, PAIRS_DCABEAMLINE, PAIRS_DECAYLENGTH };
enum ERRORS{ OK, RP_PLANES, RP_FIDUCIAL, RP_XI, RP_ELASTIC, TRACKS_TOF, TRACKS_PT, TRACKS_ETA, TRACKS_NHITS, PAIRS_DCADAUGHTERS, PAIRS_DCABEAMLINE, PAIRS_DECAYLENGTH };

const double particleMass[nParticles] = { 0.13957, 0.493677, 0.93827 }; // pion, kaon, proton in GeV /c^2 
bool IsInXiElasticSpot(StUPCRpsTrack *, StUPCRpsTrack *);
bool IsInMomElasticSpot(StUPCRpsTrack *, StUPCRpsTrack *);
bool tuple_sort(tuple<double, int, int>, tuple<double, int, int>);
bool tuple_sort_big(tuple<double, int, int, bool>, tuple<double, int, int, bool>);
bool isPi(StUPCTrack *);
bool isProton(StUPCTrack *);
double deltaT(TVector3 decayVertex, StUPCTrack *track1, StUPCTrack *track2, double m1, double m2);
double RP_PV_z(StRPEvent *rpevent);
//central detector track test functions
bool kToFFlag(StUPCTrack *tempTrack);
bool pTRange(StUPCTrack *tempTrack);
bool etaRange(StUPCTrack *tempTrack);
bool NHitsNumber(StUPCTrack *tempTrack);
bool emptyTest(StUPCTrack *tempTrack);
//central detector pair test function
bool dcaDaughters(StUPCV0 *tempParticle);
bool dcaBeamline(StUPCV0 *tempParticle);
bool decayLengthPointingAngle(StUPCV0 *tempParticle);
bool emptyParticleTest(StUPCV0 *tempParticle);

int main(int argc, char *argv[]){
    int nthreads = 1;
    enum dataType{ old_data, new_data_old_tracks, new_data_new_tracks };
    int inputDataType = -1;
    if(argc==5){
        nthreads = atoi(argv[3]);
        inputDataType = atoi(argv[4]);
    } else if(argc==4){
        inputDataType = atoi(argv[3]);
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

    //useful constants
    //k0
    double kaonMassWindowNarrowLow = 0.48;
    double kaonMassWindowNarrowHigh = 0.52;
    int kaonMassWindowNarrowBins = 40;
    double kaonMassWindowWideLow = 0.44;
    double kaonMassWindowWideHigh = 0.54;
    int kaonMassWindowWideBins = 100;
    //lambda
    double lambdaMassWindowNarrowLow = 1.09;
    double lambdaMassWindowNarrowHigh = 1.13;
    int lambdaMassWindowNarrowBins = 40;
    double lambdaMassWindowWideLow = 1.06;
    double lambdaMassWindowWideHigh = 1.16;
    int lambdaMassWindowWideBins = 100;
    vector<vector<double>> beamData = ReadFillPositionData("STAR-Analysis/share/Run7PolarizationWithPosition.csv");

    //histograms
    ProcessingOutsideLoop outsideprocessing;
    //RP_PLANES
    outsideprocessing.AddHistogram(TH2D("RP_PLANES_east", "RP_PLANES_east", 100, -0.5, 0.5, 200, -1, 1));
    outsideprocessing.AddHistogram(TH2D("RP_PLANES_west", "RP_PLANES_west", 100, -0.5, 0.5, 200, -1, 1));
    //RP_FIDUCIAL
    outsideprocessing.AddHistogram(TH2D("RP_FIDUCIAL_east", "RP_FIDUCIAL_east", 200, -1., 1., 300, -1.5, 1.5));
    outsideprocessing.AddHistogram(TH2D("RP_FIDUCIAL_west", "RP_FIDUCIAL_west", 200, -1., 1., 300, -1.5, 1.5));
    //RP_XI
    outsideprocessing.AddHistogram(TH1D("RP_XI_east_CPT2", "RP_XI_east_CPT2", 600, -0.05, 0.25));
    outsideprocessing.AddHistogram(TH1D("RP_XI_west_CPT2", "RP_XI_west_CPT2", 600, -0.05, 0.25));
    outsideprocessing.AddHistogram(TH1D("RP_XI_east_CPT2noBBCL", "RP_XI_east_CPT2noBBCL", 600, -0.05, 0.25));
    outsideprocessing.AddHistogram(TH1D("RP_XI_west_CPT2noBBCL", "RP_XI_west_CPT2noBBCL", 600, -0.05, 0.25));
    //RP_ELASTIC
    outsideprocessing.AddHistogram(TH2D("RP_ELASTIC_xi", "RP_ELASTIC_xi", 60, -0.01, 0.05, 60, -0.01, 0.05));
    outsideprocessing.AddHistogram(TH2D("RP_ELASTIC_theta", "RP_ELASTIC_theta", 60, -3e-3, 3e-3, 60, -2e-3, 2e-3));
    outsideprocessing.AddHistogram(TH2D("RP_ELASTIC_p", "RP_ELASTIC_p", 80, -0.6, 1., 50, -0.5, 0.5));
    //TRACKS_TOF
    outsideprocessing.AddHistogram(TH1D("TRACKS_TOF", "TRACKS_TOF", 40, 0, 40));
    outsideprocessing.AddHistogram(TH1D("TRACKS_TOF_OK", "TRACKS_TOF_OK", 40, 0, 40));
    //TRACKS_PT
    outsideprocessing.AddHistogram(TH1D("TRACKS_PT", "TRACKS_PT", 300, 0, 3));
    //TRACKS_ETA
    outsideprocessing.AddHistogram(TH1D("TRACKS_ETA", "TRACKS_ETA", 240, -1.2, 1.2));
    //TRACKS_NHITS
    outsideprocessing.AddHistogram(TH1D("TRACKS_NHITS", "TRACKS_NHITS", 60, 0, 60));
    //PAIRS_DCADAUGHTERS
    outsideprocessing.AddHistogram(TH1D("PAIRS_DCADAUGHTERS_K0", "PAIRS_DCADAUGHTERS_K0", 120, 0, 6));
    outsideprocessing.AddHistogram(TH1D("PAIRS_DCADAUGHTERS_Lambda0", "PAIRS_DCADAUGHTERS_Lambda0", 120, 0, 6));
    //PAIRS_DCABEAMLINE
    outsideprocessing.AddHistogram(TH1D("PAIRS_DCABEAMLINE_K0", "PAIRS_DCABEAMLINE_K0", 120, 0, 6));
    outsideprocessing.AddHistogram(TH1D("PAIRS_DCABEAMLINE_Lambda0", "PAIRS_DCABEAMLINE_Lambda0", 120, 0, 6));
    //PAIRS_DECAYLENGTH
    outsideprocessing.AddHistogram(TH1D("PAIRS_DECAYLENGTH_K0_LEN", "PAIRS_DECAYLENGTH_K0_LEN", 240, 0, 12));
    outsideprocessing.AddHistogram(TH1D("PAIRS_DECAYLENGTH_Lambda0_LEN", "PAIRS_DECAYLENGTH_Lambda0_LEN", 240, 0, 12));
    outsideprocessing.AddHistogram(TH1D("PAIRS_DECAYLENGTH_K0_ANGLE", "PAIRS_DECAYLENGTH_K0_ANGLE", 100, 0.9, 1.));
    outsideprocessing.AddHistogram(TH1D("PAIRS_DECAYLENGTH_Lambda0_ANGLE", "PAIRS_DECAYLENGTH_Lambda0_ANGLE", 100, 0.9, 1.));
    outsideprocessing.AddHistogram(TH2D("PAIRS_DECAYLENGTH_K0", "PAIRS_DECAYLENGTH_K0", 240, 0, 12, 100, 0.9, 1));
    outsideprocessing.AddHistogram(TH2D("PAIRS_DECAYLENGTH_Lambda0", "PAIRS_DECAYLENGTH_Lambda0", 240, 0, 12, 100, 0.9, 1));
    outsideprocessing.AddHistogram(TH2D("PAIRS_DECAYLENGTH_K0_BIGGER_RANGE", "PAIRS_DECAYLENGTH_K0_BIGGER_RANGE", 240, 0, 12, 1000, 0., 1));
    outsideprocessing.AddHistogram(TH2D("PAIRS_DECAYLENGTH_Lambda0_BIGGER_RANGE", "PAIRS_DECAYLENGTH_Lambda0_BIGGER_RANGE", 240, 0, 12, 1000, 0., 1));
    //REMOVING_DUPLICATES
    outsideprocessing.AddHistogram(TH1D("REMOVING_DUPLICATES_K0_BEFORE", "REMOVING_DUPLICATES_K0_BEFORE", kaonMassWindowWideBins, kaonMassWindowWideLow, kaonMassWindowWideHigh));
    outsideprocessing.AddHistogram(TH1D("REMOVING_DUPLICATES_K0_AFTER", "REMOVING_DUPLICATES_K0_AFTER", kaonMassWindowWideBins, kaonMassWindowWideLow, kaonMassWindowWideHigh));
    outsideprocessing.AddHistogram(TH1D("REMOVING_DUPLICATES_Lambda0_BEFORE", "REMOVING_DUPLICATES_Lambda0_BEFORE", lambdaMassWindowWideBins, lambdaMassWindowWideLow, lambdaMassWindowWideHigh));
    outsideprocessing.AddHistogram(TH1D("REMOVING_DUPLICATES_Lambda0_AFTER", "REMOVING_DUPLICATES_Lambda0_AFTER", lambdaMassWindowWideBins, lambdaMassWindowWideLow, lambdaMassWindowWideHigh));

    //processing
    //defining TreeProcessor
    ROOT::TTreeProcessorMT TreeProc(*upcChain, nthreads);

    int eventsProcessed = 0;

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
        int badCut;
        double firstBranch, secondBranch;
        bool f1, f2, f3;
        double px, py;
        double xiLowerLimit, xiUpperLimit;
        std::vector<StUPCTrack *> vector_Track_positive;
        std::vector<StUPCTrack *> vector_Track_negative;
        StUPCTrack *tempTrack;
        TVector3 vertexPrimary;
        std::vector<double> tempBeamVector;
        double beamValues[4];
        StUPCV0 *tempParticle;
        std::vector<tuple<double, int, int>> vector_K0_pairs;
        vector_K0_pairs.reserve(100);
        std::vector<tuple<double, int, int, bool>> vector_Lambda_pairs;
        vector_Lambda_pairs.reserve(100);
        std::vector<int> vector_all;
        vector_all.reserve(100);
        bool didItStartWithOK = true;
        bool isTrackOkay = true;
        bool isParticleOkay = true;

        bool (*centralTracksCutsTab[])(StUPCTrack *) = { kToFFlag, pTRange, etaRange, NHitsNumber };
        int centralTracksCutsTabLen = 4;
        bool (*centralV0sCutsTab[])(StUPCV0 *) = { dcaDaughters, dcaBeamline, decayLengthPointingAngle };
        int centralV0sCutsTabLen = 3;

        StUPCRpsTrack *eastTrack;
        StUPCRpsTrack *westTrack;

        //actual loop
        while(myReader.Next()){
            //in a TTree, it *would* be constant, in TChain however not necessarily
            tempUPCpointer = StUPCEventInstance.Get();
            tempRPpointer = StRPEventInstance.Get();

            //cleaning the loop
            badCut = OK;

            //cause I want to see what's going on
            if(eventsProcessed%10000==0){
                cout<<"Processed "<<eventsProcessed<<" events"<<endl;
            }
            eventsProcessed++;

            //cuts already taken care of during preselection:
            //triggers
            //2 protons, 1 on each side
            //two stations, with >=3 planes on each side

            //cuts remaining
            //RP planes
            for(unsigned int k = 0; k<tempRPpointer->getNumberOfTracks(); ++k){
                // Get pointer to k-th track in Roman Pot data collection
                StUPCRpsTrack *trk = tempRPpointer->getTrack(k);
                trk->setEvent(tempRPpointer);
                if(trk->getTrackPoint(0)->planesUsed()<3||trk->getTrackPoint(1)->planesUsed()<3){
                    badCut = RP_PLANES;
                    break;
                }
            }
            //fiducial region
            for(unsigned int k = 0; k<tempRPpointer->getNumberOfTracks(); ++k){
                // Get pointer to k-th track in Roman Pot data collection
                StUPCRpsTrack *trk = tempRPpointer->getTrack(k);
                trk->setEvent(tempRPpointer);
                px = trk->pVec().X();
                py = trk->pVec().Y();
                f1 = (0.4<abs(py)&&abs(py)<0.8);
                f2 = (-0.27<px);
                f3 = (pow(px+0.6, 2)+pow(py, 2)<1.25);
                if(!(f1&&f2&&f3)){
                    if(badCut!=OK){
                        goto EndLoop;
                    }
                    badCut = RP_FIDUCIAL;
                    break;
                }
            }
            //Xi test
            //can do that because of the preselection
            if(tempUPCpointer->isTrigger(570701)&&tempUPCpointer->getRunNumber()<=18083025){
                xiLowerLimit = 0.005;
                xiUpperLimit = 0.2;
            } else if(tempUPCpointer->isTrigger(570705)&&tempUPCpointer->getRunNumber()>18083025){
                xiLowerLimit = 0.005;
                xiUpperLimit = 0.08;
            }
            for(unsigned int k = 0; k<tempRPpointer->getNumberOfTracks(); ++k){
                // Get pointer to k-th track in Roman Pot data collection
                StUPCRpsTrack *trk = tempRPpointer->getTrack(k);
                trk->setEvent(tempRPpointer);
                if(trk->xi(beamMomentum)<=xiLowerLimit or trk->xi(beamMomentum)>=xiUpperLimit){
                    if(badCut!=OK){
                        goto EndLoop;
                    }
                    badCut = RP_XI;
                    break;
                }
            }
            //elastic collisions
            if(IsInXiElasticSpot(tempRPpointer->getTrack(0), tempRPpointer->getTrack(1)) or IsInMomElasticSpot(tempRPpointer->getTrack(0), tempRPpointer->getTrack(1))){
                if(badCut!=OK){
                    goto EndLoop;
                }
                badCut = RP_ELASTIC;
            }
            //central detector stuff
            didItStartWithOK = (badCut==OK);
            for(int cut = 0; cut<centralTracksCutsTabLen+centralV0sCutsTabLen+1; cut++){
                //cleaning the loop
                vector_Track_positive.clear();
                vector_Track_negative.clear();
                vector_K0_pairs.clear();
                vector_Lambda_pairs.clear();
                vector_all.clear();
                //piece to omit all cuts if one was already "used" earlier
                if(!didItStartWithOK){
                    cut = centralTracksCutsTabLen+centralV0sCutsTabLen;
                }
                //piece to go through all cuts, including the "okay" one
                if(didItStartWithOK&&cut!=centralTracksCutsTabLen+centralV0sCutsTabLen){
                    //setting badCut to the cut used
                    //RP_ELASTIC is the last before central cuts
                    //change if necessary
                    badCut = RP_ELASTIC+1+cut;
                }
                if(didItStartWithOK&&cut==centralTracksCutsTabLen+centralV0sCutsTabLen){
                    badCut = OK;
                }
                //cuts
                for(int track = 0; track<tempUPCpointer->getNumberOfTracks(); track++){
                    tempTrack = tempUPCpointer->getTrack(track);
                    //mandatory ones
                    if(inputDataType==new_data_old_tracks){
                        //kV0 false and kCEP false
                        if(tempTrack->getFlag(StUPCTrack::kV0)||tempTrack->getFlag(StUPCTrack::kCEP)){
                            continue;
                        }
                    } else if(inputDataType==new_data_new_tracks){
                        //kV0 true
                        if(!tempTrack->getFlag(StUPCTrack::kV0)){
                            continue;
                        }
                    }
                    //optional ones
                    isTrackOkay = true;
                    for(int trackCut = 0; trackCut<centralTracksCutsTabLen; trackCut++){
                        if(trackCut==cut)
                            continue;
                        if(!centralTracksCutsTab[trackCut](tempTrack)){
                            isTrackOkay = false;
                            break;
                        }
                    }
                    //filling the vectors
                    if(isTrackOkay&&tempTrack->getCharge()>0){
                        vector_Track_positive.push_back(tempTrack);
                    } else if(isTrackOkay){
                        vector_Track_negative.push_back(tempTrack);
                    }
                }
                //working on pairs
                switch(inputDataType){
                case old_data:
                    vertexPrimary = { tempUPCpointer->getVertex(0)->getPosX(), tempUPCpointer->getVertex(0)->getPosY(), tempUPCpointer->getVertex(0)->getPosZ() };
                    tempBeamVector = FindPosition(tempUPCpointer->getFillNumber(), vertexPrimary.Z(), beamData[0], beamData[1], beamData[2], beamData[3], beamData[4], beamData[5], beamData[6], beamData[7], beamData[8]);
                    beamValues[0] = tempBeamVector[0];
                    beamValues[1] = tempBeamVector[1];
                    beamValues[2] = tempBeamVector[2];
                    beamValues[3] = tempBeamVector[3];
                    break;
                case new_data_old_tracks:
                case new_data_new_tracks:
                    beamValues[0] = tempUPCpointer->getBeamXPosition();
                    beamValues[1] = tempUPCpointer->getBeamYPosition();
                    beamValues[2] = tempUPCpointer->getBeamXSlope();
                    beamValues[3] = tempUPCpointer->getBeamYSlope();
                    vertexPrimary = { beamValues[0]+RP_PV_z(tempRPpointer)*beamValues[2], beamValues[1]+RP_PV_z(tempRPpointer)*beamValues[3], RP_PV_z(tempRPpointer) };
                    break;
                default:
                    break;
                }
                //actual loop for K0
                for(long unsigned int i = 0; i<vector_Track_positive.size(); i++){
                    for(long unsigned int j = 0; j<vector_Track_negative.size(); j++){
                        // tempParticle = new StUPCV0(vector_Track[i], vector_Track[j], particleMass[0], particleMass[0], 1, 1, vertexPrimary, beamValues, tempUPCpointer->getMagneticField(), false);
                        tempParticle = new StUPCV0(vector_Track_positive[i], vector_Track_negative[j], particleMass[0], particleMass[0], 1, 1, { 0,0,0 }, beamValues, tempUPCpointer->getMagneticField(), false);
                        //tests before filling
                        isParticleOkay = true;
                        for(int particleCut = 0; particleCut<centralV0sCutsTabLen; particleCut++){
                            if(particleCut+centralTracksCutsTabLen==cut)
                                continue;
                            if(!centralV0sCutsTab[particleCut](tempParticle)){
                                isParticleOkay = false;
                                break;
                            }
                        }
                        //filling
                        double mK0candidate = tempParticle->m();
                        if(mK0candidate>kaonMassWindowWideLow&&mK0candidate<kaonMassWindowWideHigh&&isParticleOkay){
                            vector_K0_pairs.push_back(tuple(tempParticle->dcaDaughters(), i, j));
                        }
                        //finishing
                        delete tempParticle;
                    }
                }
                //actual loop for Lambda (1/2)
                //positive is a proton
                for(long unsigned int i = 0; i<vector_Track_positive.size(); i++){
                    for(long unsigned int j = 0; j<vector_Track_negative.size(); j++){
                        tempParticle = new StUPCV0(vector_Track_positive[i], vector_Track_negative[j], particleMass[2], particleMass[0], 1, 1, { 0,0,0 }, beamValues, tempUPCpointer->getMagneticField(), false);
                        //tests if accept the particle
                        //tests before filling
                        isParticleOkay = true;
                        for(int particleCut = 0; particleCut<centralV0sCutsTabLen; particleCut++){
                            if(particleCut+centralTracksCutsTabLen==cut)
                                continue;
                            if(!centralV0sCutsTab[particleCut](tempParticle)){
                                isParticleOkay = false;
                                break;
                            }
                        }
                        //filling
                        double mLambdacandidate = tempParticle->m();
                        if(mLambdacandidate>lambdaMassWindowWideLow&&mLambdacandidate<lambdaMassWindowWideHigh&&isParticleOkay){
                            vector_Lambda_pairs.push_back(tuple(tempParticle->dcaDaughters(), i, j, true));
                        }
                        //finishing
                        delete tempParticle;
                    }
                }
                //actual loop for Lambda (2/2)
                //negative is an antiproton
                for(long unsigned int i = 0; i<vector_Track_positive.size(); i++){
                    for(long unsigned int j = 0; j<vector_Track_negative.size(); j++){
                        tempParticle = new StUPCV0(vector_Track_positive[i], vector_Track_negative[j], particleMass[0], particleMass[2], 1, 1, { 0,0,0 }, beamValues, tempUPCpointer->getMagneticField(), false);
                        //tests if accept the particle
                        //tests before filling
                        isParticleOkay = true;
                        for(int particleCut = 0; particleCut<centralV0sCutsTabLen; particleCut++){
                            if(particleCut+centralTracksCutsTabLen==cut)
                                continue;
                            if(!centralV0sCutsTab[particleCut](tempParticle)){
                                isParticleOkay = false;
                                break;
                            }
                        }
                        //filling
                        double mLambdacandidate = tempParticle->m();
                        if(mLambdacandidate>lambdaMassWindowWideLow&&mLambdacandidate<lambdaMassWindowWideHigh){
                            vector_Lambda_pairs.push_back(tuple(tempParticle->dcaDaughters(), i, j, false));
                        }
                        //finishing
                        delete tempParticle;
                    }
                }
                //check particle deduplication results
                if(badCut==OK){
                    for(size_t i = 0; i<vector_K0_pairs.size(); i++){
                        tempParticle = new StUPCV0(vector_Track_positive[std::get<1>(vector_K0_pairs[i])], vector_Track_negative[std::get<2>(vector_K0_pairs[i])], particleMass[0], particleMass[0], 1, 1, { 0,0,0 }, beamValues, tempUPCpointer->getMagneticField(), false);
                        insideprocessing.Fill("REMOVING_DUPLICATES_K0_BEFORE", tempParticle->m());
                        delete tempParticle;
                    }
                    for(size_t i = 0; i<vector_Lambda_pairs.size(); i++){
                        if(std::get<3>(vector_Lambda_pairs[i])){
                            tempParticle = new StUPCV0(vector_Track_positive[std::get<1>(vector_Lambda_pairs[i])], vector_Track_negative[std::get<2>(vector_Lambda_pairs[i])], particleMass[2], particleMass[0], 1, 1, { 0,0,0 }, beamValues, tempUPCpointer->getMagneticField(), false);
                        } else{
                            tempParticle = new StUPCV0(vector_Track_positive[std::get<1>(vector_Lambda_pairs[i])], vector_Track_negative[std::get<2>(vector_Lambda_pairs[i])], particleMass[0], particleMass[2], 1, 1, { 0,0,0 }, beamValues, tempUPCpointer->getMagneticField(), false);
                        }
                        insideprocessing.Fill("REMOVING_DUPLICATES_Lambda0_BEFORE", tempParticle->m());
                        delete tempParticle;
                    }
                }
                //get rid of pairs with shared particles
                //we sort in ascending order of DCAdaugters
                std::sort(vector_K0_pairs.begin(), vector_K0_pairs.end(), tuple_sort);
                std::sort(vector_Lambda_pairs.begin(), vector_Lambda_pairs.end(), tuple_sort_big);
                //removing duplicates from both vectors, even cross-appearing
                //for vector_all, i for vector_K0_pairs will be i+1, and j for vector_Lambda_pairs will be -(j+1)
                int counter_K0 = 0;
                int counter_Lambda = 0;
                double pos_part_index_i, neg_part_index_i, pos_part_index_j, neg_part_index_j;
                for(int counter = 0; counter<vector_K0_pairs.size()+vector_Lambda_pairs.size(); counter++){
                    //if counter comes to an end of either vector, just copy the rest of the other one
                    if(counter_K0>=vector_K0_pairs.size()){
                        vector_all.push_back(-counter_Lambda-1);
                        counter_Lambda++;
                        continue;
                    }
                    if(counter_Lambda>=vector_Lambda_pairs.size()){
                        vector_all.push_back(counter_K0+1);
                        counter_K0++;
                        continue;
                    }
                    //if that didn't happen, compare normally 
                    if(std::get<0>(vector_K0_pairs[counter_K0])<std::get<0>(vector_Lambda_pairs[counter_Lambda])){
                        vector_all.push_back(counter_K0+1);
                        counter_K0++;
                    } else{
                        vector_all.push_back(-counter_Lambda-1);
                        counter_Lambda++;
                    }
                }
                //now erase repeating particles
                for(int i = 0; i<int(vector_all.size())-1; i++){
                    for(size_t j = i+1; j<vector_all.size(); j++){
                        if(vector_all[i]>0){
                            pos_part_index_i = std::get<1>(vector_K0_pairs[vector_all[i]-1]);
                            neg_part_index_i = std::get<2>(vector_K0_pairs[vector_all[i]-1]);
                        } else{
                            pos_part_index_i = std::get<1>(vector_Lambda_pairs[-vector_all[i]-1]);
                            neg_part_index_i = std::get<2>(vector_Lambda_pairs[-vector_all[i]-1]);
                        }
                        if(vector_all[j]>0){
                            pos_part_index_j = std::get<1>(vector_K0_pairs[vector_all[j]-1]);
                            neg_part_index_j = std::get<2>(vector_K0_pairs[vector_all[j]-1]);
                        } else{
                            pos_part_index_j = std::get<1>(vector_Lambda_pairs[-vector_all[j]-1]);
                            neg_part_index_j = std::get<2>(vector_Lambda_pairs[-vector_all[j]-1]);
                        }
                        //removing the index and value from all vectors, including vector_all
                        if(pos_part_index_i==pos_part_index_j or neg_part_index_i==neg_part_index_j){
                            //first erasing from particle-dependent vectors
                            if(vector_all[j]>0){
                                vector_K0_pairs.erase(vector_K0_pairs.begin()+vector_all[j]-1);
                            } else{
                                vector_Lambda_pairs.erase(vector_Lambda_pairs.begin()-vector_all[j]-1);
                            }
                            //then from allvector
                            vector_all.erase(vector_all.begin()+j);
                            j--;
                        }
                    }
                }
                //now fill the histograms (only if there are actually any particles)
                if(vector_K0_pairs.size()==0&&vector_Lambda_pairs.size()==0)
                    continue;
                //check which one to fill
                switch(badCut){
                case TRACKS_TOF:
                    insideprocessing.Fill("TRACKS_TOF", vector_Track_negative.size()+vector_Track_positive.size());
                    break;
                case TRACKS_PT:
                    for(size_t i = 0; i<vector_Track_negative.size(); i++){
                        insideprocessing.Fill("TRACKS_PT", vector_Track_negative[i]->getPt());
                    }
                    for(size_t i = 0; i<vector_Track_positive.size(); i++){
                        insideprocessing.Fill("TRACKS_PT", vector_Track_positive[i]->getPt());
                    }
                    break;
                case TRACKS_ETA:
                    for(size_t i = 0; i<vector_Track_negative.size(); i++){
                        insideprocessing.Fill("TRACKS_ETA", vector_Track_negative[i]->getEta());
                    }
                    for(size_t i = 0; i<vector_Track_positive.size(); i++){
                        insideprocessing.Fill("TRACKS_ETA", vector_Track_positive[i]->getEta());
                    }
                    break;
                case TRACKS_NHITS:
                    for(size_t i = 0; i<vector_Track_negative.size(); i++){
                        insideprocessing.Fill("TRACKS_NHITS", vector_Track_negative[i]->getNhits());
                    }
                    for(size_t i = 0; i<vector_Track_positive.size(); i++){
                        insideprocessing.Fill("TRACKS_NHITS", vector_Track_positive[i]->getNhits());
                    }
                    break;
                case PAIRS_DCADAUGHTERS:
                    for(size_t i = 0; i<vector_K0_pairs.size(); i++){
                        tempParticle = new StUPCV0(vector_Track_positive[std::get<1>(vector_K0_pairs[i])], vector_Track_negative[std::get<2>(vector_K0_pairs[i])], particleMass[0], particleMass[0], 1, 1, { 0,0,0 }, beamValues, tempUPCpointer->getMagneticField(), false);
                        insideprocessing.Fill("PAIRS_DCADAUGHTERS_K0", tempParticle->dcaDaughters());
                        delete tempParticle;
                    }
                    for(size_t i = 0; i<vector_Lambda_pairs.size(); i++){
                        if(std::get<3>(vector_Lambda_pairs[i])){
                            tempParticle = new StUPCV0(vector_Track_positive[std::get<1>(vector_Lambda_pairs[i])], vector_Track_negative[std::get<2>(vector_Lambda_pairs[i])], particleMass[2], particleMass[0], 1, 1, { 0,0,0 }, beamValues, tempUPCpointer->getMagneticField(), false);
                        } else{
                            tempParticle = new StUPCV0(vector_Track_positive[std::get<1>(vector_Lambda_pairs[i])], vector_Track_negative[std::get<2>(vector_Lambda_pairs[i])], particleMass[0], particleMass[2], 1, 1, { 0,0,0 }, beamValues, tempUPCpointer->getMagneticField(), false);
                        }
                        insideprocessing.Fill("PAIRS_DCADAUGHTERS_Lambda0", tempParticle->dcaDaughters());
                        delete tempParticle;
                    }
                    break;
                case PAIRS_DCABEAMLINE:
                    for(size_t i = 0; i<vector_K0_pairs.size(); i++){
                        tempParticle = new StUPCV0(vector_Track_positive[std::get<1>(vector_K0_pairs[i])], vector_Track_negative[std::get<2>(vector_K0_pairs[i])], particleMass[0], particleMass[0], 1, 1, { 0,0,0 }, beamValues, tempUPCpointer->getMagneticField(), false);
                        insideprocessing.Fill("PAIRS_DCABEAMLINE_K0", tempParticle->DCABeamLine());
                        delete tempParticle;
                    }
                    for(size_t i = 0; i<vector_Lambda_pairs.size(); i++){
                        if(std::get<3>(vector_Lambda_pairs[i])){
                            tempParticle = new StUPCV0(vector_Track_positive[std::get<1>(vector_Lambda_pairs[i])], vector_Track_negative[std::get<2>(vector_Lambda_pairs[i])], particleMass[2], particleMass[0], 1, 1, { 0,0,0 }, beamValues, tempUPCpointer->getMagneticField(), false);
                        } else{
                            tempParticle = new StUPCV0(vector_Track_positive[std::get<1>(vector_Lambda_pairs[i])], vector_Track_negative[std::get<2>(vector_Lambda_pairs[i])], particleMass[0], particleMass[2], 1, 1, { 0,0,0 }, beamValues, tempUPCpointer->getMagneticField(), false);
                        }
                        insideprocessing.Fill("PAIRS_DCABEAMLINE_Lambda0", tempParticle->DCABeamLine());
                        delete tempParticle;
                    }
                    break;
                case PAIRS_DECAYLENGTH:
                    for(size_t i = 0; i<vector_K0_pairs.size(); i++){
                        tempParticle = new StUPCV0(vector_Track_positive[std::get<1>(vector_K0_pairs[i])], vector_Track_negative[std::get<2>(vector_K0_pairs[i])], particleMass[0], particleMass[0], 1, 1, { 0,0,0 }, beamValues, tempUPCpointer->getMagneticField(), false);
                        insideprocessing.Fill("PAIRS_DECAYLENGTH_K0_LEN", tempParticle->decayLengthHypo());
                        insideprocessing.Fill("PAIRS_DECAYLENGTH_K0_ANGLE", tempParticle->pointingAngleHypo());
                        insideprocessing.Fill("PAIRS_DECAYLENGTH_K0", tempParticle->decayLengthHypo(), tempParticle->pointingAngleHypo());
                        insideprocessing.Fill("PAIRS_DECAYLENGTH_K0_BIGGER_RANGE", tempParticle->decayLengthHypo(), tempParticle->pointingAngleHypo());
                        delete tempParticle;
                    }
                    for(size_t i = 0; i<vector_Lambda_pairs.size(); i++){
                        if(std::get<3>(vector_Lambda_pairs[i])){
                            tempParticle = new StUPCV0(vector_Track_positive[std::get<1>(vector_Lambda_pairs[i])], vector_Track_negative[std::get<2>(vector_Lambda_pairs[i])], particleMass[2], particleMass[0], 1, 1, { 0,0,0 }, beamValues, tempUPCpointer->getMagneticField(), false);
                        } else{
                            tempParticle = new StUPCV0(vector_Track_positive[std::get<1>(vector_Lambda_pairs[i])], vector_Track_negative[std::get<2>(vector_Lambda_pairs[i])], particleMass[0], particleMass[2], 1, 1, { 0,0,0 }, beamValues, tempUPCpointer->getMagneticField(), false);
                        }
                        insideprocessing.Fill("PAIRS_DECAYLENGTH_Lambda0_LEN", tempParticle->decayLengthHypo());
                        insideprocessing.Fill("PAIRS_DECAYLENGTH_Lambda0_ANGLE", tempParticle->pointingAngleHypo());
                        insideprocessing.Fill("PAIRS_DECAYLENGTH_Lambda0", tempParticle->decayLengthHypo(), tempParticle->pointingAngleHypo());
                        insideprocessing.Fill("PAIRS_DECAYLENGTH_Lambda0_BIGGER_RANGE", tempParticle->decayLengthHypo(), tempParticle->pointingAngleHypo());
                        delete tempParticle;
                    }
                    break;
                case OK:
                    insideprocessing.Fill("TRACKS_TOF_OK", vector_Track_negative.size()+vector_Track_positive.size());
                    for(size_t i = 0; i<vector_K0_pairs.size(); i++){
                        tempParticle = new StUPCV0(vector_Track_positive[std::get<1>(vector_K0_pairs[i])], vector_Track_negative[std::get<2>(vector_K0_pairs[i])], particleMass[0], particleMass[0], 1, 1, { 0,0,0 }, beamValues, tempUPCpointer->getMagneticField(), false);
                        insideprocessing.Fill("REMOVING_DUPLICATES_K0_AFTER", tempParticle->m());
                        delete tempParticle;
                    }
                    for(size_t i = 0; i<vector_Lambda_pairs.size(); i++){
                        if(std::get<3>(vector_Lambda_pairs[i])){
                            tempParticle = new StUPCV0(vector_Track_positive[std::get<1>(vector_Lambda_pairs[i])], vector_Track_negative[std::get<2>(vector_Lambda_pairs[i])], particleMass[2], particleMass[0], 1, 1, { 0,0,0 }, beamValues, tempUPCpointer->getMagneticField(), false);
                        } else{
                            tempParticle = new StUPCV0(vector_Track_positive[std::get<1>(vector_Lambda_pairs[i])], vector_Track_negative[std::get<2>(vector_Lambda_pairs[i])], particleMass[0], particleMass[2], 1, 1, { 0,0,0 }, beamValues, tempUPCpointer->getMagneticField(), false);
                        }
                        insideprocessing.Fill("REMOVING_DUPLICATES_Lambda0_AFTER", tempParticle->m());
                        delete tempParticle;
                    }
                    break;
                default:
                    //means there's RP cut, do not interfere then
                    break;
                }
            }
            if(vector_K0_pairs.size()==0&&vector_Lambda_pairs.size()==0)
                goto EndLoop;
            //filling the histograms based on which cut was violated 
            if(tempRPpointer->getTrack(0)->pVec().Z()<0){
                eastTrack = tempRPpointer->getTrack(0);
                westTrack = tempRPpointer->getTrack(1);
            } else{
                eastTrack = tempRPpointer->getTrack(1);
                westTrack = tempRPpointer->getTrack(0);
            }
            switch(badCut){
            case RP_PLANES:
                insideprocessing.Fill("RP_PLANES_east", eastTrack->pVec().X(), eastTrack->pVec().Y());
                insideprocessing.Fill("RP_PLANES_west", westTrack->pVec().X(), westTrack->pVec().Y());
                break;
            case RP_FIDUCIAL:
                insideprocessing.Fill("RP_FIDUCIAL_east", eastTrack->pVec().X(), eastTrack->pVec().Y());
                insideprocessing.Fill("RP_FIDUCIAL_west", westTrack->pVec().X(), westTrack->pVec().Y());
                break;
            case RP_XI:
                if(tempUPCpointer->isTrigger(570701)&&tempUPCpointer->getRunNumber()<=18083025){
                    insideprocessing.Fill("RP_XI_east_CPT2", eastTrack->xi(beamMomentum));
                    insideprocessing.Fill("RP_XI_west_CPT2", westTrack->xi(beamMomentum));
                } else if(tempUPCpointer->isTrigger(570705)&&tempUPCpointer->getRunNumber()>18083025){
                    insideprocessing.Fill("RP_XI_east_CPT2noBBCL", eastTrack->xi(beamMomentum));
                    insideprocessing.Fill("RP_XI_west_CPT2noBBCL", westTrack->xi(beamMomentum));
                }
                break;
            case RP_ELASTIC:
                insideprocessing.Fill("RP_ELASTIC_xi", eastTrack->xi(beamMomentum), westTrack->xi(beamMomentum));
                insideprocessing.Fill("RP_ELASTIC_theta", eastTrack->theta(StUPCRpsTrack::rpsAngleThetaX)+westTrack->theta(StUPCRpsTrack::rpsAngleThetaX), eastTrack->theta(StUPCRpsTrack::rpsAngleThetaY)+westTrack->theta(StUPCRpsTrack::rpsAngleThetaY));
                insideprocessing.Fill("RP_ELASTIC_p", eastTrack->pVec().X()+westTrack->pVec().X(), eastTrack->pVec().Y()+westTrack->pVec().Y());
                break;
            case OK:
                insideprocessing.Fill("RP_PLANES_east", eastTrack->pVec().X(), eastTrack->pVec().Y());
                insideprocessing.Fill("RP_PLANES_west", westTrack->pVec().X(), westTrack->pVec().Y());
                insideprocessing.Fill("RP_FIDUCIAL_east", eastTrack->pVec().X(), eastTrack->pVec().Y());
                insideprocessing.Fill("RP_FIDUCIAL_west", westTrack->pVec().X(), westTrack->pVec().Y());
                if(tempUPCpointer->isTrigger(570701)&&tempUPCpointer->getRunNumber()<=18083025){
                    insideprocessing.Fill("RP_XI_east_CPT2", eastTrack->xi(beamMomentum));
                    insideprocessing.Fill("RP_XI_west_CPT2", westTrack->xi(beamMomentum));
                } else if(tempUPCpointer->isTrigger(570705)&&tempUPCpointer->getRunNumber()>18083025){
                    insideprocessing.Fill("RP_XI_east_CPT2noBBCL", eastTrack->xi(beamMomentum));
                    insideprocessing.Fill("RP_XI_west_CPT2noBBCL", westTrack->xi(beamMomentum));
                }
                insideprocessing.Fill("RP_ELASTIC_xi", eastTrack->xi(beamMomentum), westTrack->xi(beamMomentum));
                insideprocessing.Fill("RP_ELASTIC_theta", eastTrack->theta(StUPCRpsTrack::rpsAngleThetaX)+westTrack->theta(StUPCRpsTrack::rpsAngleThetaX), eastTrack->theta(StUPCRpsTrack::rpsAngleThetaY)+westTrack->theta(StUPCRpsTrack::rpsAngleThetaY));
                insideprocessing.Fill("RP_ELASTIC_p", eastTrack->pVec().X()+westTrack->pVec().X(), eastTrack->pVec().Y()+westTrack->pVec().Y());
                break;
            default:
                //used when cuts are those on central state
                //so no histograms get filled
                break;
            }

            //endpoint for twice failed events
        EndLoop:
            continue;
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
        if(inputDataType==old_data){
            outfileName = outputFolder+"AnaOutput_"+path.substr(path.find_last_of("/\\")+1)+".root";
        } else if(inputDataType==new_data_old_tracks){
            outfileName = outputFolder+"AnaOutput_"+path.substr(path.find_last_of("/\\")+1)+"_new_data_old_tracks.root";
        } else if(inputDataType==new_data_new_tracks){
            outfileName = outputFolder+"AnaOutput_"+path.substr(path.find_last_of("/\\")+1)+"_new_data_new_tracks.root";
        } else{
            outfileName = outputFolder+"AnaOutput_"+path.substr(path.find_last_of("/\\")+1)+"_unspecified.root";
        }
    }
    cout<<"Created output file "<<outfileName<<endl;
    TFile *outputFileHist = TFile::Open(outfileName.c_str(), "recreate");
    outsideprocessing.SaveToFile(outputFileHist);
    outputFileHist->Close();

    return 0;

    return 0;
}

//NON MAIN
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


bool tuple_sort(tuple<double, int, int> t1, tuple<double, int, int> t2){
    return std::get<0>(t1)<std::get<0>(t2);
}

bool tuple_sort_big(tuple<double, int, int, bool> t1, tuple<double, int, int, bool> t2){
    return std::get<0>(t1)<std::get<0>(t2);
}

bool isPi(StUPCTrack *track){
    if(track->getNhitsDEdx()<20){
        return false;
    }
    return abs(track->getNSigmasTPCPion())<3;
}

bool isProton(StUPCTrack *track){
    if(track->getNhitsDEdx()<20){
        return false;
    }
    return abs(track->getNSigmasTPCProton())<3;
}

double deltaT(TVector3 decayVertex, StUPCTrack *track1, StUPCTrack *track2, double m1, double m2){
    double middleOfToF = 214.4;
    StPicoPhysicalHelix p1Helix = StPicoPhysicalHelix(track1->getCurvature(),
        track1->getDipAngle(),
        track1->getPhase(),
        track1->getOrigin(),
        track1->getCharge());
    StPicoPhysicalHelix p2Helix = StPicoPhysicalHelix(track2->getCurvature(),
        track2->getDipAngle(),
        track2->getPhase(),
        track2->getOrigin(),
        track2->getCharge());

    //my way
    p1Helix.moveOrigin(p1Helix.pathLength(decayVertex));
    p2Helix.moveOrigin(p2Helix.pathLength(decayVertex));
    std::pair<double, double> stemp = p1Helix.pathLength(middleOfToF);
    double s1 = (stemp.first<=0) ? stemp.second : min(stemp.first, stemp.second);
    stemp = p2Helix.pathLength(middleOfToF);
    double s2 = (stemp.first<=0) ? stemp.second : min(stemp.first, stemp.second);
    TLorentzVector vectTemp;
    double t1, t2;
    //to change BoostVector from v/c to v in cm/ToF time
    //ToFlst in ps/ToFtime:
    double ToFlst = 24;
    //m/s = 100cm/1e12ps = 100cm/(1e12ps/ToFlst) = 1e-10cm*ToFlst
    double coeff = 299792458.0*1e-10*ToFlst;
    track1->getLorentzVector(vectTemp, m1);
    t1 = track1->getTofTime()-s1/(vectTemp.BoostVector().Mag()*coeff);
    track2->getLorentzVector(vectTemp, m2);
    t2 = track2->getTofTime()-s2/(vectTemp.BoostVector().Mag()*coeff);

    //ToF way
    // TLorentzVector vectTemp;
    // double t1, t2;
    // track1->getLorentzVector(vectTemp, m1);
    // t1 = track1->getTofTime()-track1->getTofPathLength()/vectTemp.BoostVector().Mag();
    // track2->getLorentzVector(vectTemp, m2);
    // t2 = track2->getTofTime()-track2->getTofPathLength()/vectTemp.BoostVector().Mag();

    // std::cout<<"AAAAAAAAAAAAAAAAAAAAA"<<std::endl;
    // std::cout<<track2->getTofTime()<<" "<<s2/(vectTemp.BoostVector().Mag()*coeff)<<" "<<t2<<" "<<t1-t2<<std::endl;
    return t1-t2;

    // track1->getLorentzVector(vectTemp, m1);
    // double temp = s1/(vectTemp.BoostVector().Mag()*299792458.0*1e-10);
    // track2->getLorentzVector(vectTemp, m2);
    // temp -= s2/(vectTemp.BoostVector().Mag()*299792458.0*1e-10);
    // if(track1->getTofTime()-track2->getTofTime()==0)
    //     printf("Same: %f\n", track1->getTofTime());
    // return (track1->getTofTime()-track2->getTofTime())/temp;
}

double RP_PV_z(StRPEvent *rpevent){
    StUPCRpsTrack *p_plus = rpevent->getTrack(0);
    p_plus->setEvent(rpevent);
    StUPCRpsTrack *p_minus;
    if(p_plus->getTrackPoint(0)->z()>0){
        p_plus = rpevent->getTrack(0);
        p_minus = rpevent->getTrack(1);
        p_plus->setEvent(rpevent);
        p_minus->setEvent(rpevent);
    } else{
        p_plus = rpevent->getTrack(1);
        p_minus = rpevent->getTrack(0);
        p_plus->setEvent(rpevent);
        p_minus->setEvent(rpevent);
    }
    double v1, v2, d1, d2, t1, t2;
    TLorentzVector temp;
    t1 = p_plus->time();
    t2 = p_minus->time();
    // printf("%lf %lf\n", t1*1e9, t2*1e9);
    d1 = 0.5*(p_plus->getTrackPoint(0)->z()+p_plus->getTrackPoint(1)->z());
    d2 = -0.5*(p_minus->getTrackPoint(0)->z()+p_minus->getTrackPoint(1)->z());
    temp.SetVectM(p_plus->pVec(), particleMass[Proton]);
    v1 = temp.BoostVector().Z()*299792458.0*100; //na cm/s
    temp.SetVectM(p_minus->pVec(), particleMass[Proton]);
    v2 = -temp.BoostVector().Z()*299792458.0*100; //na cm/s

    return (t1-t2-(d1/v1-d2/v2))/(1/v1+1/v2);
}

//CENTRAL DETECTOR TEST FUNCTIONS
bool kToFFlag(StUPCTrack *tempTrack){
    return tempTrack->getFlag(StUPCTrack::kTof);
}

bool pTRange(StUPCTrack *tempTrack){
    return tempTrack->getPt()>0.2;
}

bool etaRange(StUPCTrack *tempTrack){
    return abs(tempTrack->getEta())<0.9;
}

bool NHitsNumber(StUPCTrack *tempTrack){
    return tempTrack->getNhits()>20;
}

bool emptyTest(StUPCTrack *tempTrack){
    return true;
}

//CENTRAL DETECTOR PAIR TESTS
bool dcaDaughters(StUPCV0 *tempParticle){
    return tempParticle->dcaDaughters()<=2.5;

}

bool dcaBeamline(StUPCV0 *tempParticle){
    return tempParticle->DCABeamLine()<=2.5;
}

bool decayLengthPointingAngle(StUPCV0 *tempParticle){
    return tempParticle->pointingAngleHypo()>0.925 or tempParticle->decayLengthHypo()<3.0;
}

bool emptyParticleTest(StUPCV0 *tempParticle){
    return true;
}
