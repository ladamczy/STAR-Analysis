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

bool tuple_sort(tuple<double, int, int>, tuple<double, int, int>);
bool tuple_sort_big(tuple<double, int, int, bool>, tuple<double, int, int, bool>);
bool isPi(StUPCTrack *);
bool isProton(StUPCTrack *);

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
    outsideprocessing.AddHistogram(TH1D("MpipiNarrow", "K^{0}_{S} mass in narrow range;m_{#pi^{+}#pi^{-}} [GeV];Number of pairs", kaonMassWindowNarrowBins, kaonMassWindowNarrowLow, kaonMassWindowNarrowHigh));
    outsideprocessing.AddHistogram(TH1D("MpipiWide", "K^{0}_{S} mass in wide range;m_{#pi^{+}#pi^{-}} [GeV];Number of pairs", kaonMassWindowWideBins, kaonMassWindowWideLow, kaonMassWindowWideHigh));
    outsideprocessing.AddHistogram(TH1D("MpipiVeryWide", "Pion pair mass in a very wide range;m_{#pi^{+}#pi^{-}} [GeV];Number of pairs", 500, 0, 5));
    outsideprocessing.AddHistogram(TH1D("MpipiNarrowWithPidEcut", "K^{0}_{S} mass in narrow range with dE/dx cuts;m_{#pi^{+}#pi^{-}} [GeV];Number of pairs", kaonMassWindowNarrowBins, kaonMassWindowNarrowLow, kaonMassWindowNarrowHigh));
    outsideprocessing.AddHistogram(TH1D("MpipiVeryWideWithPidEcut", "Pion pair mass in a very wide range with dE/dx cuts;m_{#pi^{+}#pi^{-}} [GeV];Number of pairs", 500, 0, 5));
    outsideprocessing.AddHistogram(TH1D("NotherTracks", "Number of tracks excluding K^{0}_{S} daughters", 20, 0, 20));
    outsideprocessing.AddHistogram(TH2D("VertexMultiplicity", "Number of good quality tracks attributed to a vertex;vertex ID;tracks number", 10, 0, 10, 15, 0, 15));
    outsideprocessing.AddHistogram(TH2D("VertexMultiplicityWithoutK0", "Number of good quality tracks attributed to a vertex except those making a K0;vertex ID;tracks number", 10, 0, 10, 15, 0, 15));
    outsideprocessing.AddHistogram(TH1D("K0multiplicity", "Number of K^{0}_{S} detected", 5, 0, 5));
    outsideprocessing.AddHistogram(TH1D("PVV0K0dist", "PV - V0K0 z axis distance;d[cm]", 200, -10, 10));
    outsideprocessing.AddHistogram(TH1D("invdecaylenghthHypocut", "Mass of pion pair with decaylenghthHypo>3cm cut", kaonMassWindowWideBins, kaonMassWindowWideLow, kaonMassWindowWideHigh));
    outsideprocessing.AddHistogram(TH1D("K0PVdistance", "Distance between K0 line-of-flight and PV;d[cm]", 100, 0, 100));
    outsideprocessing.AddHistogram(TH2D("VertexIdvsPos", "Vertex ID vs vertex list position;position;ID", 10, 0, 10, 10, 0, 10));
    outsideprocessing.AddHistogram(TH2D("VertexPrimvsAll", "Number of primary vertices vs number of all vertices;all;primary", 10, 0, 10, 10, 0, 10));
    //labda
    outsideprocessing.AddHistogram(TH1D("MppiNarrow", "#Lambda^{0} mass in narrow range;m_{p^{#pm}#pi^{#mp}} [GeV];Number of pairs", lambdaMassWindowNarrowBins, lambdaMassWindowNarrowLow, lambdaMassWindowNarrowHigh));
    outsideprocessing.AddHistogram(TH1D("MppiWide", "#Lambda^{0} mass in wide range;m_{p^{#pm}#pi^{#mp}} [GeV];Number of pairs", lambdaMassWindowWideBins, lambdaMassWindowWideLow, lambdaMassWindowWideHigh));
    outsideprocessing.AddHistogram(TH1D("MppiVeryWide", "Pion pair mass in a very wide range;m_{p^{#pm}#pi^{#mp}} [GeV];Number of pairs", 500, 0, 5));
    outsideprocessing.AddHistogram(TH1D("MppiNarrowWithPidEcut", "#Lambda^{0} mass in narrow range with dE/dx cuts;m_{p^{#pm}#pi^{#mp}} [GeV];Number of pairs", lambdaMassWindowNarrowBins, lambdaMassWindowNarrowLow, lambdaMassWindowNarrowHigh));
    outsideprocessing.AddHistogram(TH1D("MppiVeryWideWithPidEcut", "Pion pair mass in a very wide range with dE/dx cuts;m_{p^{#pm}#pi^{#mp}} [GeV];Number of pairs", 500, 0, 5));

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
        // StRPEvent *tempRPpointer;
        insideprocessing.GetLocalHistograms(&outsideprocessing);

        //helpful variables
        std::vector<StUPCTrack *> vector_Track_positive;
        std::vector<StUPCTrack *> vector_Track_negative;
        StUPCTrack *tempTrack;
        TVector3 vertexPrimary;
        std::vector<double> tempBeamVector;
        double beamValues[4];
        StUPCV0 *tempParticle;
        std::vector<tuple<double, int, int>> vector_K0_pairs;
        vector_K0_pairs.reserve(100);
        std::vector<int> numberOfTracksTiedToVertex(100, 0);
        std::vector<tuple<double, int, int, bool>> vector_Lambda_pairs;
        vector_Lambda_pairs.reserve(100);

        //actual loop
        while(myReader.Next()){
            //in a TTree, it *would* be constant, in TChain however not necessarily
            tempUPCpointer = StUPCEventInstance.Get();
            // tempRPpointer = StRPEventInstance.Get();

            //test


            //cleaning the loop
            vector_Track_positive.clear();
            vector_Track_negative.clear();
            vector_K0_pairs.clear();
            vector_Lambda_pairs.clear();
            std::fill(numberOfTracksTiedToVertex.begin(), numberOfTracksTiedToVertex.end(), 0);

            //cause I want to see what's going on
            if(eventsProcessed%10000<nthreads){
                cout<<"Processed "<<eventsProcessed<<" events"<<endl;
            }
            eventsProcessed++;

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
            //testing amount of tracks tied to a vertex
            for(size_t i = 0; i<vector_Track_positive.size(); i++){
                numberOfTracksTiedToVertex[vector_Track_positive[i]->getVertexId()]++;
            }
            for(size_t i = 0; i<vector_Track_negative.size(); i++){
                numberOfTracksTiedToVertex[vector_Track_negative[i]->getVertexId()]++;
            }
            for(size_t i = 0; i<numberOfTracksTiedToVertex.size(); i++){
                if(numberOfTracksTiedToVertex[i]!=0)
                    insideprocessing.Fill("VertexMultiplicity", i, numberOfTracksTiedToVertex[i]);
            }
            //selecting all K0s matching criteria:
            // DCAdaugter<=2.5cm
            // DCAbeamline<=2.5cm
            // decayLenghth<3cm or cos(pointing angle)>0.925
            for(int i = 0; i<tempUPCpointer->getNumberOfVertices(); i++){
                insideprocessing.Fill("VertexIdvsPos", i, tempUPCpointer->getVertex(i)->getId());
            }
            insideprocessing.Fill("VertexPrimvsAll", tempUPCpointer->getNumberOfVertices(), tempUPCpointer->getNPrimVertices());
            vertexPrimary = { tempUPCpointer->getVertex(0)->getPosX(), tempUPCpointer->getVertex(0)->getPosY(), tempUPCpointer->getVertex(0)->getPosZ() };
            tempBeamVector = FindPosition(tempUPCpointer->getFillNumber(), vertexPrimary.Z(), beamData[0], beamData[1], beamData[2], beamData[3], beamData[4], beamData[5], beamData[6], beamData[7], beamData[8]);
            beamValues[0] = tempBeamVector[0];
            beamValues[1] = tempBeamVector[1];
            beamValues[2] = tempBeamVector[2];
            beamValues[3] = tempBeamVector[3];
            //actual loop for K0
            for(long unsigned int i = 0; i<vector_Track_positive.size(); i++){
                for(long unsigned int j = 0; j<vector_Track_negative.size(); j++){
                    // tempParticle = new StUPCV0(vector_Track[i], vector_Track[j], particleMass[0], particleMass[0], 1, 1, vertexPrimary, beamValues, tempUPCpointer->getMagneticField(), false);
                    tempParticle = new StUPCV0(vector_Track_positive[i], vector_Track_negative[j], particleMass[0], particleMass[0], 1, 1, { 0,0,0 }, beamValues, tempUPCpointer->getMagneticField(), false);
                    //tests if accept the particle
                    //test unnecessary due to previous costraints
                    // bool K0test1 = vector_Track[i]->getCharge()*vector_Track[j]->getCharge()<0;
                    bool K0test2 = tempParticle->dcaDaughters()<=2.5;
                    bool K0test3 = tempParticle->DCABeamLine()<=2.5;
                    bool K0test4 = tempParticle->pointingAngleHypo()>0.925;
                    bool K0test5 = tempParticle->decayLengthHypo()<3.0;
                    //tests before filling
                    double mK0candidate = tempParticle->m();
                    if(K0test2&&K0test3&&(K0test4||(!K0test5))){
                        insideprocessing.Fill("invdecaylenghthHypocut", mK0candidate);
                    }
                    if(!(K0test2&&K0test3&&(K0test4||K0test5))){
                        continue;
                    }
                    //filling
                    if(mK0candidate>kaonMassWindowNarrowLow&&mK0candidate<kaonMassWindowNarrowHigh){
                        vector_K0_pairs.push_back(tuple(tempParticle->dcaDaughters(), i, j));
                    }
                    if(mK0candidate>kaonMassWindowWideLow&&mK0candidate<kaonMassWindowWideHigh){
                        insideprocessing.Fill("MpipiWide", mK0candidate);
                    }
                    insideprocessing.Fill("MpipiVeryWide", mK0candidate);
                    if(isPi(vector_Track_positive[i])&&isPi(vector_Track_negative[j])){
                        insideprocessing.Fill("MpipiVeryWideWithPidEcut", mK0candidate);
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
                    //test unnecessary due to previous costraints
                    // bool K0test1 = vector_Track[i]->getCharge()*vector_Track[j]->getCharge()<0;
                    bool Lambdatest2 = tempParticle->dcaDaughters()<=2.5;
                    bool Lambdatest3 = tempParticle->DCABeamLine()<=2.5;
                    bool Lambdatest4 = tempParticle->pointingAngleHypo()>0.925;
                    bool Lambdatest5 = tempParticle->decayLengthHypo()<3.0;
                    //tests before filling
                    double mLambdacandidate = tempParticle->m();
                    if(Lambdatest2&&Lambdatest3&&(Lambdatest4||(!Lambdatest5))){
                        // insideprocessing.Fill("invdecaylenghthHypocut", mLambdacandidate);
                    }
                    if(!(Lambdatest2&&Lambdatest3&&(Lambdatest4||Lambdatest5))){
                        continue;
                    }
                    //filling
                    if(mLambdacandidate>lambdaMassWindowNarrowLow&&mLambdacandidate<lambdaMassWindowNarrowHigh){
                        vector_Lambda_pairs.push_back(tuple(tempParticle->dcaDaughters(), i, j, true));
                    }
                    if(mLambdacandidate>lambdaMassWindowWideLow&&mLambdacandidate<lambdaMassWindowWideHigh){
                        insideprocessing.Fill("MppiWide", mLambdacandidate);
                    }
                    insideprocessing.Fill("MppiVeryWide", mLambdacandidate);
                    if(isProton(vector_Track_positive[i])&&isPi(vector_Track_negative[j])){
                        insideprocessing.Fill("MppiVeryWideWithPidEcut", mLambdacandidate);
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
                    //test unnecessary due to previous costraints
                    // bool K0test1 = vector_Track[i]->getCharge()*vector_Track[j]->getCharge()<0;
                    bool Lambdatest2 = tempParticle->dcaDaughters()<=2.5;
                    bool Lambdatest3 = tempParticle->DCABeamLine()<=2.5;
                    bool Lambdatest4 = tempParticle->pointingAngleHypo()>0.925;
                    bool Lambdatest5 = tempParticle->decayLengthHypo()<3.0;
                    //tests before filling
                    double mLambdacandidate = tempParticle->m();
                    if(Lambdatest2&&Lambdatest3&&(Lambdatest4||(!Lambdatest5))){
                        // insideprocessing.Fill("invdecaylenghthHypocut", mLambdacandidate);
                    }
                    if(!(Lambdatest2&&Lambdatest3&&(Lambdatest4||Lambdatest5))){
                        continue;
                    }
                    //filling
                    if(mLambdacandidate>lambdaMassWindowNarrowLow&&mLambdacandidate<lambdaMassWindowNarrowHigh){
                        vector_Lambda_pairs.push_back(tuple(tempParticle->dcaDaughters(), i, j, false));
                    }
                    if(mLambdacandidate>lambdaMassWindowWideLow&&mLambdacandidate<lambdaMassWindowWideHigh){
                        insideprocessing.Fill("MppiWide", mLambdacandidate);
                    }
                    insideprocessing.Fill("MppiVeryWide", mLambdacandidate);
                    if(isPi(vector_Track_positive[i])&&isProton(vector_Track_negative[j])){
                        insideprocessing.Fill("MppiVeryWideWithPidEcut", mLambdacandidate);
                    }
                    //finishing
                    delete tempParticle;
                }
            }
            //get rid of pairs with shared particles
            //we sort in ascending order of DCAdaugters
            std::sort(vector_K0_pairs.begin(), vector_K0_pairs.end(), tuple_sort);
            std::sort(vector_Lambda_pairs.begin(), vector_Lambda_pairs.end(), tuple_sort_big);
            //first, we check for repeating pi+ in K0 candidates
            for(int i = 0; i<int(vector_K0_pairs.size())-1; i++){
                for(size_t j = i+1; j<vector_K0_pairs.size(); j++){
                    if(std::get<1>(vector_K0_pairs[i])==std::get<1>(vector_K0_pairs[j])){
                        vector_K0_pairs.erase(vector_K0_pairs.begin()+j);
                        j--;
                    }
                }
            }
            //then, for repeating pi- in K0 candidates
            for(int i = 0; i<int(vector_K0_pairs.size())-1; i++){
                for(size_t j = i+1; j<vector_K0_pairs.size(); j++){
                    if(std::get<2>(vector_K0_pairs[i])==std::get<2>(vector_K0_pairs[j])){
                        vector_K0_pairs.erase(vector_K0_pairs.begin()+j);
                        j--;
                    }
                }
            }
            //and then we do this for positive tracks of Lambda daughters
            for(int i = 0; i<int(vector_Lambda_pairs.size())-1; i++){
                for(size_t j = i+1; j<vector_Lambda_pairs.size(); j++){
                    if(std::get<1>(vector_Lambda_pairs[i])==std::get<1>(vector_Lambda_pairs[j])){
                        vector_Lambda_pairs.erase(vector_Lambda_pairs.begin()+j);
                        j--;
                    }
                }
            }
            //and then for negative
            for(int i = 0; i<int(vector_Lambda_pairs.size())-1; i++){
                for(size_t j = i+1; j<vector_Lambda_pairs.size(); j++){
                    if(std::get<2>(vector_Lambda_pairs[i])==std::get<2>(vector_Lambda_pairs[j])){
                        vector_Lambda_pairs.erase(vector_Lambda_pairs.begin()+j);
                        j--;
                    }
                }
            }
            //loop after all remaining K0
            for(size_t i = 0; i<vector_K0_pairs.size(); i++){
                tempParticle = new StUPCV0(vector_Track_positive[std::get<1>(vector_K0_pairs[i])], vector_Track_negative[std::get<2>(vector_K0_pairs[i])], particleMass[0], particleMass[0], 1, 1, { 0,0,0 }, beamValues, tempUPCpointer->getMagneticField(), false);
                insideprocessing.Fill("MpipiNarrow", tempParticle->m());
                insideprocessing.Fill("PVV0K0dist", tempUPCpointer->getVertex(0)->getPosZ()-tempParticle->prodVertexHypo().Z());
                insideprocessing.Fill("K0PVdistance", tempParticle->DcaToPrimaryVertex());
                if(isPi(vector_Track_positive[std::get<1>(vector_K0_pairs[i])])&&isPi(vector_Track_negative[std::get<2>(vector_K0_pairs[i])])){
                    insideprocessing.Fill("MpipiNarrowWithPidEcut", tempParticle->m());
                }
                delete tempParticle;
            }
            insideprocessing.Fill("K0multiplicity", vector_K0_pairs.size());
            //sum of not used particles
            if(vector_K0_pairs.size()>0){
                int numberOfTracksNotUsed = vector_Track_positive.size()+vector_Track_negative.size()-2*vector_K0_pairs.size();
                insideprocessing.Fill("NotherTracks", numberOfTracksNotUsed);
            }
            //testing amount of tracks tied to a vertex minus those used for K0
            for(size_t i = 0; i<vector_K0_pairs.size(); i++){
                numberOfTracksTiedToVertex[vector_Track_positive[std::get<1>(vector_K0_pairs[i])]->getVertexId()]--;
                numberOfTracksTiedToVertex[vector_Track_negative[std::get<2>(vector_K0_pairs[i])]->getVertexId()]--;
            }
            for(size_t i = 0; i<numberOfTracksTiedToVertex.size(); i++){
                if(numberOfTracksTiedToVertex[i]!=0)
                    insideprocessing.Fill("VertexMultiplicityWithoutK0", i, numberOfTracksTiedToVertex[i]);
            }

            //loop after all remaining Lambda
            for(size_t i = 0; i<vector_Lambda_pairs.size(); i++){
                if(std::get<3>(vector_Lambda_pairs[i])){
                    tempParticle = new StUPCV0(vector_Track_positive[std::get<1>(vector_Lambda_pairs[i])], vector_Track_negative[std::get<2>(vector_Lambda_pairs[i])], particleMass[2], particleMass[0], 1, 1, { 0,0,0 }, beamValues, tempUPCpointer->getMagneticField(), false);
                } else{
                    tempParticle = new StUPCV0(vector_Track_positive[std::get<1>(vector_Lambda_pairs[i])], vector_Track_negative[std::get<2>(vector_Lambda_pairs[i])], particleMass[0], particleMass[2], 1, 1, { 0,0,0 }, beamValues, tempUPCpointer->getMagneticField(), false);
                }
                insideprocessing.Fill("MppiNarrow", tempParticle->m());
                // insideprocessing.Fill("PVV0Lambdadist", tempUPCpointer->getVertex(0)->getPosZ()-tempParticle->prodVertexHypo().Z());
                // insideprocessing.Fill("LambdaPVdistance", tempParticle->DcaToPrimaryVertex());
                if(std::get<3>(vector_Lambda_pairs[i])&&isProton(vector_Track_positive[std::get<1>(vector_Lambda_pairs[i])])&&isPi(vector_Track_negative[std::get<2>(vector_Lambda_pairs[i])])){
                    insideprocessing.Fill("MppiNarrowWithPidEcut", tempParticle->m());
                } else if((!std::get<3>(vector_Lambda_pairs[i]))&&isPi(vector_Track_positive[std::get<1>(vector_Lambda_pairs[i])])&&isProton(vector_Track_negative[std::get<2>(vector_Lambda_pairs[i])])){
                    insideprocessing.Fill("MppiNarrowWithPidEcut", tempParticle->m());
                }

                delete tempParticle;
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