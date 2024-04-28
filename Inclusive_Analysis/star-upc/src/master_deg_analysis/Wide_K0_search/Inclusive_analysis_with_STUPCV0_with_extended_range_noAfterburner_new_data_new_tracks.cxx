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
#include "StPicoPhysicalHelix.h"

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
double deltaT(TVector3 decayVertex, StUPCTrack *track1, StUPCTrack *track2, double m1, double m2);

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

    //histograms
    ProcessingOutsideLoop outsideprocessing;
    outsideprocessing.AddHistogram(TH1D("MpipiNarrow", "K^{0}_{S} mass in narrow range;m_{#pi^{+}#pi^{-}} [GeV];Number of pairs", kaonMassWindowNarrowBins, kaonMassWindowNarrowLow, kaonMassWindowNarrowHigh));
    outsideprocessing.AddHistogram(TH1D("MpipiWide", "K^{0}_{S} mass in wide range;m_{#pi^{+}#pi^{-}} [GeV];Number of pairs", kaonMassWindowWideBins, kaonMassWindowWideLow, kaonMassWindowWideHigh));
    outsideprocessing.AddHistogram(TH1D("MpipiVeryWide", "Pion pair mass in a very wide range;m_{#pi^{+}#pi^{-}} [GeV];Number of pairs", 500, 0, 5));
    outsideprocessing.AddHistogram(TH1D("MpipiNarrowWithPidEcut", "K^{0}_{S} mass in narrow range with dE/dx cuts;m_{#pi^{+}#pi^{-}} [GeV];Number of pairs", kaonMassWindowNarrowBins, kaonMassWindowNarrowLow, kaonMassWindowNarrowHigh));
    outsideprocessing.AddHistogram(TH1D("MpipiWideWithPidEcut", "K^{0}_{S} mass in a very wide range with dE/dx cuts;m_{#pi^{+}#pi^{-}} [GeV];Number of pairs", kaonMassWindowWideBins, kaonMassWindowWideLow, kaonMassWindowWideHigh));
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
    outsideprocessing.AddHistogram(TH1D("MppiWideWithPidEcut", "#Lambda^{0} mass in wide range with dE/dx cuts;m_{p^{#pm}#pi^{#mp}} [GeV];Number of pairs", lambdaMassWindowWideBins, lambdaMassWindowWideLow, lambdaMassWindowWideHigh));
    outsideprocessing.AddHistogram(TH1D("MppiVeryWideWithPidEcut", "Pion pair mass in a very wide range with dE/dx cuts;m_{p^{#pm}#pi^{#mp}} [GeV];Number of pairs", 500, 0, 5));

    outsideprocessing.AddHistogram(TH1D("etaK0", "K^{0}_{S} events density by #eta;#eta;", 60, -1.5, 1.5));
    outsideprocessing.AddHistogram(TH1D("etapiK0", "#pi from K^{0}_{S} decay density by #eta;#eta;", 60, -1.5, 1.5));
    outsideprocessing.AddHistogram(TH1D("etaLambda", "#Lambda^{0} events density by #eta;#eta;", 60, -1.5, 1.5));
    outsideprocessing.AddHistogram(TH1D("etapiLambda", "#pi from #Lambda^{0} decay density by #eta;#eta;", 60, -1.5, 1.5));
    outsideprocessing.AddHistogram(TH1D("phiK0", "K^{0}_{S} events density by #phi;#phi;", 63, -TMath::Pi(), TMath::Pi()));
    outsideprocessing.AddHistogram(TH1D("phipiK0", "#pi from K^{0}_{S} decay density by #phi;#phi;", 63, -TMath::Pi(), TMath::Pi()));
    outsideprocessing.AddHistogram(TH1D("phiLambda", "#Lambda^{0} events density by #phi;#phi;", 63, -TMath::Pi(), TMath::Pi()));
    outsideprocessing.AddHistogram(TH1D("phipiLambda", "#pi from #Lambda^{0} decay density by #phi;#phi;", 63, -TMath::Pi(), TMath::Pi()));

    outsideprocessing.AddHistogram(TH1D("MpipiWideMissing", "K^{0}_{S} mass in wide range which doesn't pass the criteria;m_{#pi^{+}#pi^{-}} [GeV];Number of pairs", kaonMassWindowWideBins, kaonMassWindowWideLow, kaonMassWindowWideHigh));
    outsideprocessing.AddHistogram(TH1D("MppiWideMissing", "#Lambda^{0} mass in wide range which doesn't pass the criteria;m_{p^{#pm}#pi^{#mp}} [GeV];Number of pairs", lambdaMassWindowWideBins, lambdaMassWindowWideLow, lambdaMassWindowWideHigh));
    outsideprocessing.AddHistogram(TH1D("MpipiWideCheck", "K^{0}_{S} mass in wide range with additional cuts;m_{#pi^{+}#pi^{-}} [GeV];Number of pairs", kaonMassWindowWideBins, kaonMassWindowWideLow, kaonMassWindowWideHigh));
    outsideprocessing.AddHistogram(TH1D("MppiWideCheck", "#Lambda^{0} mass in wide range with additional cuts;m_{p^{#pm}#pi^{#mp}} [GeV];Number of pairs", lambdaMassWindowWideBins, lambdaMassWindowWideLow, lambdaMassWindowWideHigh));

    outsideprocessing.AddHistogram(TH1D("K0dcaDaughtersSignal", "K0dcaDaughtersSignal", 50, 0, 2.5));
    outsideprocessing.AddHistogram(TH1D("K0dcaDaughtersBackground", "K0dcaDaughtersBackground", 50, 0, 2.5));
    outsideprocessing.AddHistogram(TH1D("LambdadcaDaughtersSignal", "LambdadcaDaughtersSignal", 50, 0, 2.5));
    outsideprocessing.AddHistogram(TH1D("LambdadcaDaughtersBackground", "LambdadcaDaughtersBackground", 50, 0, 2.5));
    outsideprocessing.AddHistogram(TH1D("K0dcaBeamlineSignal", "K0dcaBeamlineSignal", 50, 0, 2.5));
    outsideprocessing.AddHistogram(TH1D("K0dcaBeamlineBackground", "K0dcaBeamlineBackground", 50, 0, 2.5));
    outsideprocessing.AddHistogram(TH1D("LambdadcaBeamlineSignal", "LambdadcaBeamlineSignal", 50, 0, 2.5));
    outsideprocessing.AddHistogram(TH1D("LambdadcaBeamlineBackground", "LambdadcaBeamlineBackground", 50, 0, 2.5));
    outsideprocessing.AddHistogram(TH1D("K0pointingAngleHypoSignal", "K0pointingAngleHypoSignal", 50, 0.8, 1));
    outsideprocessing.AddHistogram(TH1D("K0pointingAngleHypoBackground", "K0pointingAngleHypoBackground", 50, 0.8, 1));
    outsideprocessing.AddHistogram(TH1D("LambdapointingAngleHypoSignal", "LambdapointingAngleHypoSignal", 50, 0.8, 1));
    outsideprocessing.AddHistogram(TH1D("LambdapointingAngleHypoBackground", "LambdapointingAngleHypoBackground", 50, 0.8, 1));
    outsideprocessing.AddHistogram(TH1D("K0pointingAngleHypoSignalDetailed", "K0pointingAngleHypoSignalDetailed", 50, 0.98, 1));
    outsideprocessing.AddHistogram(TH1D("K0pointingAngleHypoBackgroundDetailed", "K0pointingAngleHypoBackgroundDetailed", 50, 0.98, 1));
    outsideprocessing.AddHistogram(TH1D("LambdapointingAngleHypoSignalDetailed", "LambdapointingAngleHypoSignalDetailed", 50, 0.98, 1));
    outsideprocessing.AddHistogram(TH1D("LambdapointingAngleHypoBackgroundDetailed", "LambdapointingAngleHypoBackgroundDetailed", 50, 0.98, 1));
    outsideprocessing.AddHistogram(TH1D("K0decayLengthHypoSignal", "K0decayLengthHypoSignal", 50, 0, 10));
    outsideprocessing.AddHistogram(TH1D("K0decayLengthHypoBackground", "K0decayLengthHypoBackground", 50, 0, 10));
    outsideprocessing.AddHistogram(TH1D("LambdadecayLengthHypoSignal", "LambdadecayLengthHypoSignal", 50, 0, 10));
    outsideprocessing.AddHistogram(TH1D("LambdadecayLengthHypoBackground", "LambdadecayLengthHypoBackground", 50, 0, 10));
    outsideprocessing.AddHistogram(TH1D("K0cosThetaStarSignal", "K0cosThetaStarSignal", 50, -1, 1));
    outsideprocessing.AddHistogram(TH1D("K0cosThetaStarBackground", "K0cosThetaStarBackground", 50, -1, 1));
    outsideprocessing.AddHistogram(TH1D("LambdacosThetaStarSignal", "LambdacosThetaStarSignal", 50, -1, 1));
    outsideprocessing.AddHistogram(TH1D("LambdacosThetaStarBackground", "LambdacosThetaStarBackground", 50, -1, 1));

    outsideprocessing.AddHistogram(TH1D("K0ToFDecayTime", "K0ToFDecayTime", 40, -2000000, 2000000));
    outsideprocessing.AddHistogram(TH1D("LambdaToFDecayTimeSignal", "LambdaToFDecayTimeSignal", 40, -2000000, 2000000));
    outsideprocessing.AddHistogram(TH1D("LambdaToFDecayTimeBackground", "LambdaToFDecayTimeBackground", 40, -2000000, 2000000));

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
        std::vector<int> vector_all;
        vector_all.reserve(100);

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
            vector_all.clear();
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
                //tests for flags of new data
                //new data - kV0 true
                if(!tempTrack->getFlag(StUPCTrack::kV0)){
                    continue;
                }
                //rest of the tests
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
            // vertexPrimary = { tempUPCpointer->getVertex(0)->getPosX(), tempUPCpointer->getVertex(0)->getPosY(), tempUPCpointer->getVertex(0)->getPosZ() };
            // beamValues[0] = tempUPCpointer->getBeamXPosition();
            // beamValues[1] = tempUPCpointer->getBeamYPosition();
            // beamValues[2] = tempUPCpointer->getBeamXSlope();
            // beamValues[3] = tempUPCpointer->getBeamYSlope();
            beamValues[0] = 0;
            beamValues[1] = 0;
            beamValues[2] = 0;
            beamValues[3] = 0;
            //actual loop
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
                    if(!(K0test2&&K0test3&&(K0test4||(!K0test5)))){
                        insideprocessing.Fill("invdecaylenghthHypocut", mK0candidate);
                    }
                    if(!(K0test2&&K0test3&&(K0test4||K0test5))){
                        insideprocessing.Fill("MpipiWideMissing", mK0candidate);
                        continue;
                    }
                    //filling
                    if(mK0candidate>kaonMassWindowWideLow&&mK0candidate<kaonMassWindowWideHigh){
                        vector_K0_pairs.push_back(tuple(tempParticle->dcaDaughters(), i, j));
                    }
                    insideprocessing.Fill("MpipiVeryWide", tempParticle->m());
                    if(isPi(vector_Track_positive[i])&&isPi(vector_Track_negative[j])){
                        insideprocessing.Fill("MpipiVeryWideWithPidEcut", tempParticle->m());
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
                        insideprocessing.Fill("MppiWideMissing", mLambdacandidate);
                        continue;
                    }
                    //filling
                    if(mLambdacandidate>lambdaMassWindowWideLow&&mLambdacandidate<lambdaMassWindowWideHigh){
                        vector_Lambda_pairs.push_back(tuple(tempParticle->dcaDaughters(), i, j, true));
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
                        insideprocessing.Fill("MppiWideMissing", mLambdacandidate);
                        continue;
                    }
                    //filling
                    if(mLambdacandidate>lambdaMassWindowWideLow&&mLambdacandidate<lambdaMassWindowWideHigh){
                        vector_Lambda_pairs.push_back(tuple(tempParticle->dcaDaughters(), i, j, false));
                    }
                    insideprocessing.Fill("MppiVeryWide", mLambdacandidate);
                    if(isProton(vector_Track_positive[i])&&isPi(vector_Track_negative[j])){
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
                        }
                        else{
                            vector_Lambda_pairs.erase(vector_Lambda_pairs.begin()-vector_all[j]-1);
                        }
                        //then from allvector
                        vector_all.erase(vector_all.begin()+j);
                        j--;
                    }
                }
            }
            //loop after all remaining K0
            for(size_t i = 0; i<vector_K0_pairs.size(); i++){
                tempParticle = new StUPCV0(vector_Track_positive[std::get<1>(vector_K0_pairs[i])], vector_Track_negative[std::get<2>(vector_K0_pairs[i])], particleMass[0], particleMass[0], 1, 1, { 0,0,0 }, beamValues, tempUPCpointer->getMagneticField(), false);
                if(tempParticle->m()>kaonMassWindowNarrowLow&&tempParticle->m()<kaonMassWindowNarrowHigh){
                    insideprocessing.Fill("MpipiNarrow", tempParticle->m());
                }
                insideprocessing.Fill("MpipiWide", tempParticle->m());
                insideprocessing.Fill("K0ToFDecayTime", deltaT(tempParticle->decayVertex(), vector_Track_positive[std::get<1>(vector_K0_pairs[i])], vector_Track_negative[std::get<2>(vector_K0_pairs[i])], particleMass[0], particleMass[0]));
                if(isPi(vector_Track_positive[std::get<1>(vector_K0_pairs[i])])&&isPi(vector_Track_negative[std::get<2>(vector_K0_pairs[i])])){
                    if(tempParticle->m()>kaonMassWindowNarrowLow&&tempParticle->m()<kaonMassWindowNarrowHigh){
                        insideprocessing.Fill("MpipiNarrowWithPidEcut", tempParticle->m());
                    }
                    insideprocessing.Fill("MpipiWideWithPidEcut", tempParticle->m());
                }
                // insideprocessing.Fill("PVV0K0dist", tempUPCpointer->getVertex(0)->getPosZ()-tempParticle->prodVertexHypo().Z());
                insideprocessing.Fill("K0PVdistance", tempParticle->DcaToPrimaryVertex());
                insideprocessing.Fill("etaK0", tempParticle->eta());
                insideprocessing.Fill("phiK0", tempParticle->phi());
                insideprocessing.Fill("etapiK0", vector_Track_positive[std::get<1>(vector_K0_pairs[i])]->getEta());
                insideprocessing.Fill("etapiK0", vector_Track_negative[std::get<2>(vector_K0_pairs[i])]->getEta());
                insideprocessing.Fill("phipiK0", vector_Track_positive[std::get<1>(vector_K0_pairs[i])]->getPhi());
                insideprocessing.Fill("phipiK0", vector_Track_negative[std::get<2>(vector_K0_pairs[i])]->getPhi());

                if(tempParticle->pointingAngleHypo()>0.99){
                    insideprocessing.Fill("MpipiWideCheck", tempParticle->m());
                }

                //signal&background
                if(tempParticle->m()>0.49&&tempParticle->m()<0.51){
                    insideprocessing.Fill("K0dcaDaughtersSignal", tempParticle->dcaDaughters());
                    insideprocessing.Fill("K0dcaBeamlineSignal", tempParticle->DCABeamLine());
                    insideprocessing.Fill("K0pointingAngleHypoSignal", tempParticle->pointingAngleHypo());
                    insideprocessing.Fill("K0pointingAngleHypoSignalDetailed", tempParticle->pointingAngleHypo());
                    insideprocessing.Fill("K0decayLengthHypoSignal", tempParticle->decayLengthHypo());
                    insideprocessing.Fill("K0cosThetaStarSignal", tempParticle->cosThetaStar());
                } else{
                    insideprocessing.Fill("K0dcaDaughtersBackground", tempParticle->dcaDaughters());
                    insideprocessing.Fill("K0dcaBeamlineBackground", tempParticle->DCABeamLine());
                    insideprocessing.Fill("K0pointingAngleHypoBackground", tempParticle->pointingAngleHypo());
                    insideprocessing.Fill("K0pointingAngleHypoBackgroundDetailed", tempParticle->pointingAngleHypo());
                    insideprocessing.Fill("K0decayLengthHypoBackground", tempParticle->decayLengthHypo());
                    insideprocessing.Fill("K0cosThetaStarBackground", tempParticle->cosThetaStar());
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
                if(tempParticle->m()>lambdaMassWindowNarrowLow&&tempParticle->m()<lambdaMassWindowNarrowHigh){
                    insideprocessing.Fill("MppiNarrow", tempParticle->m());
                }
                insideprocessing.Fill("MppiWide", tempParticle->m());
                if(std::get<3>(vector_Lambda_pairs[i])&&isProton(vector_Track_positive[std::get<1>(vector_Lambda_pairs[i])])&&isPi(vector_Track_negative[std::get<2>(vector_Lambda_pairs[i])])){
                    if(tempParticle->m()>lambdaMassWindowNarrowLow&&tempParticle->m()<lambdaMassWindowNarrowHigh){
                        insideprocessing.Fill("MppiNarrowWithPidEcut", tempParticle->m());
                    }
                    insideprocessing.Fill("MppiWideWithPidEcut", tempParticle->m());
                } else if((!std::get<3>(vector_Lambda_pairs[i]))&&isPi(vector_Track_positive[std::get<1>(vector_Lambda_pairs[i])])&&isProton(vector_Track_negative[std::get<2>(vector_Lambda_pairs[i])])){
                    if(tempParticle->m()>lambdaMassWindowNarrowLow&&tempParticle->m()<lambdaMassWindowNarrowHigh){
                        insideprocessing.Fill("MppiNarrowWithPidEcut", tempParticle->m());
                    }
                    insideprocessing.Fill("MppiWideWithPidEcut", tempParticle->m());
                }
                // insideprocessing.Fill("PVV0Lambdadist", tempUPCpointer->getVertex(0)->getPosZ()-tempParticle->prodVertexHypo().Z());
                // insideprocessing.Fill("LambdaPVdistance", tempParticle->DcaToPrimaryVertex());
                insideprocessing.Fill("etaLambda", tempParticle->eta());
                insideprocessing.Fill("phiLambda", tempParticle->phi());
                if(std::get<3>(vector_Lambda_pairs[i])){
                    insideprocessing.Fill("etapiLambda", vector_Track_negative[std::get<2>(vector_Lambda_pairs[i])]->getEta());
                    insideprocessing.Fill("phipiLambda", vector_Track_negative[std::get<2>(vector_Lambda_pairs[i])]->getPhi());

                } else{
                    insideprocessing.Fill("etapiLambda", vector_Track_positive[std::get<1>(vector_Lambda_pairs[i])]->getEta());
                    insideprocessing.Fill("phipiLambda", vector_Track_positive[std::get<1>(vector_Lambda_pairs[i])]->getPhi());
                }
                //signal & background stuff
                bool temptest1 = std::get<3>(vector_Lambda_pairs[i])&&isProton(vector_Track_positive[std::get<1>(vector_Lambda_pairs[i])])&&isPi(vector_Track_negative[std::get<2>(vector_Lambda_pairs[i])]);
                bool temptest2 = (!std::get<3>(vector_Lambda_pairs[i]))&&isPi(vector_Track_positive[std::get<1>(vector_Lambda_pairs[i])])&&isProton(vector_Track_negative[std::get<2>(vector_Lambda_pairs[i])]);
                if((temptest1 or temptest2)&&(tempParticle->m()>1.11&&tempParticle->m()<1.12)){
                    insideprocessing.Fill("LambdadcaDaughtersSignal", tempParticle->dcaDaughters());
                    insideprocessing.Fill("LambdadcaBeamlineSignal", tempParticle->DCABeamLine());
                    insideprocessing.Fill("LambdapointingAngleHypoSignal", tempParticle->pointingAngleHypo());
                    insideprocessing.Fill("LambdapointingAngleHypoSignalDetailed", tempParticle->pointingAngleHypo());
                    insideprocessing.Fill("LambdadecayLengthHypoSignal", tempParticle->decayLengthHypo());
                    insideprocessing.Fill("LambdacosThetaStarSignal", tempParticle->cosThetaStar());
                    if(temptest1)
                        insideprocessing.Fill("LambdaToFDecayTimeSignal", deltaT(tempParticle->decayVertex(), vector_Track_positive[std::get<1>(vector_Lambda_pairs[i])], vector_Track_negative[std::get<2>(vector_Lambda_pairs[i])], particleMass[2], particleMass[0]));
                    else if(temptest2){
                        insideprocessing.Fill("LambdaToFDecayTimeSignal", deltaT(tempParticle->decayVertex(), vector_Track_positive[std::get<1>(vector_Lambda_pairs[i])], vector_Track_negative[std::get<2>(vector_Lambda_pairs[i])], particleMass[0], particleMass[2]));
                    }
                } else{
                    insideprocessing.Fill("LambdadcaDaughtersBackground", tempParticle->dcaDaughters());
                    insideprocessing.Fill("LambdadcaBeamlineBackground", tempParticle->DCABeamLine());
                    insideprocessing.Fill("LambdapointingAngleHypoBackground", tempParticle->pointingAngleHypo());
                    insideprocessing.Fill("LambdapointingAngleHypoBackgroundDetailed", tempParticle->pointingAngleHypo());
                    insideprocessing.Fill("LambdadecayLengthHypoBackground", tempParticle->decayLengthHypo());
                    insideprocessing.Fill("LambdacosThetaStarBackground", tempParticle->cosThetaStar());
                    if(temptest1)
                        insideprocessing.Fill("LambdaToFDecayTimeBackground", deltaT(tempParticle->decayVertex(), vector_Track_positive[std::get<1>(vector_Lambda_pairs[i])], vector_Track_negative[std::get<2>(vector_Lambda_pairs[i])], particleMass[2], particleMass[0]));
                    else if(temptest2){
                        insideprocessing.Fill("LambdaToFDecayTimeBackground", deltaT(tempParticle->decayVertex(), vector_Track_positive[std::get<1>(vector_Lambda_pairs[i])], vector_Track_negative[std::get<2>(vector_Lambda_pairs[i])], particleMass[0], particleMass[2]));
                    }
                }

                if(tempParticle->pointingAngleHypo()>0.99&&(temptest1 or temptest2)){
                    insideprocessing.Fill("MppiWideCheck", tempParticle->m());
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
    double ToFlst = 1e-3;//was 24
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