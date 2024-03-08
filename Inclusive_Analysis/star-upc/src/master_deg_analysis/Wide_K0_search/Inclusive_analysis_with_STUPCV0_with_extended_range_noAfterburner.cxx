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
const double particleMass[nParticles] = { 0.13957, 0.497611, 0.93827 }; // pion, kaon, proton in GeV /c^2 
enum BRANCH_ID{ EU, ED, WU, WD, nBranches };
enum RP_ID{ E1U, E1D, E2U, E2D, W1U, W1D, W2U, W2D, nRomanPots };

bool IsInXiElasticSpot(StUPCRpsTrack *, StUPCRpsTrack *);
bool IsInMomElasticSpot(StUPCRpsTrack *, StUPCRpsTrack *);

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

    //useful constants
    double kaonMassWindowNarrowLow = 0.48;
    double kaonMassWindowNarrowHigh = 0.52;
    int kaonMassWindowNarrowBins = 40;
    double kaonMassWindowWideLow = 0.44;
    double kaonMassWindowWideHigh = 0.54;
    int kaonMassWindowWideBins = 100;
    vector<vector<double>> beamData = ReadFillPositionData("STAR-Analysis/share/Run7PolarizationWithPosition.csv");

    //histograms
    ProcessingOutsideLoop outsideprocessing;
    outsideprocessing.AddHistogram(TH1D("MpipiNarrow", "K^{0}_{S} mass in narrow range;m_{#pi^{+}#pi^{-}} [GeV];Number of pairs", kaonMassWindowNarrowBins, kaonMassWindowNarrowLow, kaonMassWindowNarrowHigh));
    outsideprocessing.AddHistogram(TH1D("MpipiWide", "K^{0}_{S} mass in wide range;m_{#pi^{+}#pi^{-}} [GeV];Number of pairs", kaonMassWindowWideBins, kaonMassWindowWideLow, kaonMassWindowWideHigh));
    outsideprocessing.AddHistogram(TH1D("NotherTracks", "Number of tracks not used in K0 reconstruction", 20, 0, 20));

    // int triggers[] = { 570701, 570705, 570711, 590701, 590705, 590708 };
    // outsideprocessing.AddHistogram(TH1D("triggerHist", "Data triggers;Trigger ID;Number of events", 6, 0, 6));
    // for(int i = 0;i<6;i++){
    //     outsideprocessing.GetPointer1D(1)->GetXaxis()->SetBinLabel(i+1, to_string(triggers[i]).c_str());
    // }

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
        std::vector<StUPCTrack *> vector_Track_positive;
        std::vector<StUPCTrack *> vector_Track_negative;
        StUPCTrack *tempTrack;
        TVector3 vertexPrimary;
        std::vector<double> tempBeamVector;
        double beamValues[4];
        StUPCV0 *tempParticle;
        std::vector<int> vector_usedForK0_pair;
        vector_usedForK0_pair.reserve(100);
        std::vector<double> vector_usedForK0_DCAdaugters;
        vector_usedForK0_DCAdaugters.reserve(100);
        std::vector<int> vector_usedForK0_mass;
        vector_usedForK0_mass.reserve(100);

        //actual loop
        while(myReader.Next()){
            //in a TTree, it *would* be constant, in TChain however not necessarily
            tempUPCpointer = StUPCEventInstance.Get();
            tempRPpointer = StRPEventInstance.Get();

            //cleaning the loop
            vector_Track_positive.clear();
            vector_Track_negative.clear();
            vector_usedForK0_pair.clear();
            vector_usedForK0_DCAdaugters.clear();

            //cause I want to see what's going on
            eventsProcessed++;
            if(eventsProcessed%10000<nthreads){
                cout<<"Processed "<<eventsProcessed<<" events"<<endl;
            }

            //cuts & histogram filling
            //selecting tracks matching criteria:
            //TOF
            //pt & eta 
            //Nhits
            for(int i = 0; i<tempUPCpointer->getNumberOfTracks(); i++){
                tempTrack = tempUPCpointer->getTrack(i);
                if(!tempTrack->getFlag(StUPCTrack::kTof)){
                    continue;
                }
                if(tempTrack->getPt()<=0.2 or abs(tempTrack->getEta())>=0.9){
                    continue;
                }
                if(tempTrack->getNhits()<=20){
                    continue;
                }
                if(tempTrack->getCharge()>0){
                    vector_Track_positive.push_back(tempTrack);
                } else{
                    vector_Track_negative.push_back(tempTrack);
                }
            }
            for(size_t i = 0; i<vector_Track_positive.size(); i++){
                vector_usedForK0_pair[i] = -1;
                vector_usedForK0_DCAdaugters[i] = 0;
                vector_usedForK0_mass[i] = -1;
            }
            //selecting all K0s matching criteria:
            // DCAdaugter<=2.5cm
            // DCAbeamline<=2.5cm
            // decayLenghth<3cm or cos(pointing angle)>0.925
            vertexPrimary = { tempUPCpointer->getVertex(0)->getPosX(), tempUPCpointer->getVertex(0)->getPosY(), tempUPCpointer->getVertex(0)->getPosZ() };
            tempBeamVector = FindPosition(tempUPCpointer->getFillNumber(), vertexPrimary.Z(), beamData[0], beamData[1], beamData[2], beamData[3], beamData[4], beamData[5], beamData[6], beamData[7], beamData[8]);
            beamValues[0] = tempBeamVector[0];
            beamValues[1] = tempBeamVector[1];
            beamValues[2] = tempBeamVector[2];
            beamValues[3] = tempBeamVector[3];
            //actual loop
            for(long unsigned int i = 0; i<vector_Track_positive.size(); i++){
                for(long unsigned int j = 0; j<vector_Track_negative.size(); j++){
                    // tempParticle = new StUPCV0(vector_Track[i], vector_Track[j], particleMass[0], particleMass[0], 1, 1, vertexPrimary, beamValues, tempUPCpointer->getMagneticField(), true);
                    tempParticle = new StUPCV0(vector_Track_positive[i], vector_Track_negative[j], particleMass[0], particleMass[0], 1, 1, { 0,0,0 }, beamValues, tempUPCpointer->getMagneticField(), true);
                    //tests if accept the particle
                    //test unnecessary due to previous costraints
                    // bool K0test1 = vector_Track[i]->getCharge()*vector_Track[j]->getCharge()<0;
                    bool K0test2 = tempParticle->dcaDaughters()<=2.5;
                    bool K0test3 = tempParticle->DCABeamLine()<=2.5;
                    bool K0test4 = (tempParticle->pointingAngleHypo()>0.925 or tempParticle->decayLength()<3.0);
                    if(!(K0test2&&K0test3&&K0test4)){
                        continue;
                    }
                    //filling
                    double mK0candidate = tempParticle->m();
                    if(mK0candidate>kaonMassWindowNarrowLow&&mK0candidate<kaonMassWindowNarrowHigh){
                        //check if other K0 shares a pi+
                        if(vector_usedForK0_pair[i]>=0){
                            if(vector_usedForK0_DCAdaugters[i]>tempParticle->dcaDaughters()){
                                vector_usedForK0_pair[i] = j;
                                vector_usedForK0_DCAdaugters[i] = tempParticle->dcaDaughters();
                                vector_usedForK0_mass[i] = mK0candidate;
                            }
                            //then maybe it shares pi-?
                        } else if(std::find(vector_usedForK0_pair.begin(), vector_usedForK0_pair.end(), j)!=vector_usedForK0_pair.end()){
                            int index = std::find(vector_usedForK0_pair.begin(), vector_usedForK0_pair.end(), j)-vector_usedForK0_pair.begin();
                            if(vector_usedForK0_DCAdaugters[index]>tempParticle->dcaDaughters()){
                                vector_usedForK0_pair[index] = -1;
                                vector_usedForK0_DCAdaugters[index] = 0;
                                vector_usedForK0_mass[index] = -1;
                                vector_usedForK0_pair[i] = j;
                                vector_usedForK0_DCAdaugters[i] = tempParticle->dcaDaughters();
                                vector_usedForK0_mass[i] = mK0candidate;
                            }
                            //if a new K0 shares nothing, we simply note it down
                        } else{
                            vector_usedForK0_pair[i] = j;
                            vector_usedForK0_DCAdaugters[i] = tempParticle->dcaDaughters();
                            vector_usedForK0_mass[i] = mK0candidate;
                        }
                        // insideprocessing.Fill("MpipiNarrow", mK0candidate);
                    }
                    if(mK0candidate>kaonMassWindowWideLow&&mK0candidate<kaonMassWindowWideHigh){
                        insideprocessing.Fill("MpipiWide", mK0candidate);
                    }
                    //finishing
                    delete tempParticle;
                }
            }
            //loop after all remaining K0
            for(size_t i = 0; i<vector_usedForK0_pair.size(); i++){
                if(vector_usedForK0_mass[i]>0){
                    insideprocessing.Fill("MpipiNarrow", vector_usedForK0_mass[i]);
                }
            }
            //sum of not used particles
            int numberOfTracksNotUsed = vector_Track_positive.size()+vector_Track_negative.size();
            for(size_t i = 0; i<vector_usedForK0_pair.size(); i++){
                if(vector_usedForK0_pair[i]>=0){
                    numberOfTracksNotUsed -= 2;
                }
            }
            insideprocessing.Fill("NotherTracks", numberOfTracksNotUsed);
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