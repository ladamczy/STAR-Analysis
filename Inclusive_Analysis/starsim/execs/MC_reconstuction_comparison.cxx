#include <iostream>

//ROOT headers
#include <ROOT/TThreadedObject.hxx>
#include <TTreeReader.h>
#include <ROOT/TTreeProcessorMT.hxx>
#include <TGraph.h>
#include <TCanvas.h>
#include <TParticlePDG.h>
#include <TLine.h>
#include <TEfficiency.h>

//STAR headers
#include <StarGenEvent.h>
#include <StarGenParticle.h>

// picoDst headers
#include "StRPEvent.h"
#include "StUPCRpsTrack.h"
#include "StUPCRpsTrackPoint.h"
#include "StUPCEvent.h"
#include "StUPCTrack.h"
#include "StUPCBemcCluster.h"
#include "StUPCVertex.h"
#include "StUPCTofHit.h"
#include "StPicoPhysicalHelix.h"

//my headers
#include <ProcessingInsideLoop.h>
#include <ProcessingOutsideLoop.h>
#include <UsefulThings.h>

void PrintBigger(TParticle* input, std::string additional_stuff = "");
double distanceToBeamline(TVector3 point, double x0 = -0.126857, double y0 = -0.127072, double dxdz = -0.000945, double dydz = 0.000440);
double DeltaTFix(double deltaT, TParticle* particle1, TParticle* particle2, StUPCEvent* upcEvent);
double DeltaTFixSimilarModules(double deltaT, TParticle* particle1, TParticle* particle2, StUPCEvent* upcEvent);

int main(int argc, char* argv[]){

    //argv[1] - input file
    //argv[2] - output folder/file
    //argv[3] - PDG id/codename of particle tested
    //argv[4] - number of cores
    //argv[5] - event to print (optional)
    if(argc>6||argc<4){
        printf("Invalid number of arguments. Proper argument usage:\n");
        printf("argv[1] - input file\n");
        printf("argv[2] - output folder/file\n");
        printf("argv[3] - PDG id/codename of particle tested\n");
        printf("argv[4] - number of cores\n");
        printf("argv[5] - event to print (optional)\n");
        printf("\n");
        printf("Particle\tPDG code\tcodename\n");
        printf("K0S\t\t310\t\tK0S\n");
        printf("Lambda0\t\t3122\t\tLambda0\n");
        printf("Lambda0bar\t-3122\t\tLambda0bar\n");
        printf("K*(892)\t\t313\t\tKstar\n");
        printf("K*(892)bar\t-313\t\tKstarbar\n");
        printf("phi(1020)\t333\t\tphi\n");
        return 0;
    }

    //Useful IDs
    const int K0SPDGid = 310;
    const int LambdaPDGid = 3122;
    const int LambdabarPDGid = -3122;
    const int KstarPDGid = 313;
    const int KstarbarPDGid = -313;
    const int phiPDGid = 333;
    const int piplusPDGid = 211;
    const int piminusPDGid = -211;
    const int KplusPDGid = 321;
    const int KminusPDGid = -321;
    const int pplusPDGid = 2212;
    const int pminusPDGid = -2212;
    //masses
    std::map<int, double> massmap = { {K0SPDGid, 0.497611},
                                    {LambdaPDGid, 1.115683},
                                    {LambdabarPDGid, 1.115683},
                                    {KstarPDGid, 0.89167},
                                    {KstarbarPDGid, 0.89167},
                                    {phiPDGid, 1.019461},
                                    {piplusPDGid, 0.139570},
                                    {piminusPDGid, 0.139570},
                                    {KplusPDGid, 0.493677},
                                    {KminusPDGid, 0.493677},
                                    {pplusPDGid, 0.938272},
                                    {pminusPDGid, 0.938272} };

    int PDGmain, eventToPrint = -1;
    int PDGpositive, PDGnegative;
    if(argc==6){
        eventToPrint = atoi(argv[5]);
    }
    //atoi(not number) returns 0
    //so there is a nice check if argv[3] is a PDG id or a codename
    PDGmain = atoi(argv[3]);
    //consult the map for the key of particle codename
    //and fill PDGmain with proper id
    if(PDGmain==0){
        //why is there no switch for strings, I will never fathom
        //instead i have std::map with custom comparison function for c strings
        auto PDGmap = std::map<const char*, int, std::function<bool(const char*, const char*)>>{
            [](const char* a, const char* b){
                return strcmp(a,b)<0;
            }
        };
        PDGmap = { {"K0S", K0SPDGid},
            {"Lambda0", LambdaPDGid},
            {"Lambda0bar", LambdabarPDGid},
            {"Kstar", KstarPDGid},
            {"Kstatbar", KstarbarPDGid},
            {"phi", phiPDGid} };
        if(PDGmap.count(argv[3])==0){
            printf("This codename (%s) is not implemented/invalid\n", argv[3]);
            return 1;
        } else{
            PDGmain = PDGmap[argv[3]];
        }
    }

    printf("Chosen main particle PDG id: %d\n", PDGmain);

    //setting particle PDG number and decay products
    switch(PDGmain){
    case K0SPDGid:
        PDGpositive = piplusPDGid;
        PDGnegative = piminusPDGid;
        break;
    case LambdaPDGid:
        PDGpositive = pplusPDGid;
        PDGnegative = piminusPDGid;
        break;
    case LambdabarPDGid:
        PDGpositive = piplusPDGid;
        PDGnegative = pminusPDGid;
        break;
    case KstarPDGid:
        PDGpositive = KplusPDGid;
        PDGnegative = piminusPDGid;
        break;
    case KstarbarPDGid:
        PDGpositive = piplusPDGid;
        PDGnegative = KminusPDGid;
        break;
    case phiPDGid:
        PDGpositive = KplusPDGid;
        PDGnegative = KminusPDGid;
        break;
    default:
        printf("This PDG number is not implemented\n");
        return 1;
    }

    //multithread processing preparation
    int nthreads = 1;
    if(argc>4){
        nthreads = atoi(argv[4]);
    }
    cout<<"Program is running on "<<nthreads<<" threads"<<endl;
    ROOT::EnableThreadSafety();
    //actually i'm not sure if it's needed here
    ROOT::EnableImplicitMT(nthreads); //turn on multicore processing

    //preparing input & output
    TChain* eventFiles = new TChain("mUPCTree");
    if(ConnectInput(argc, argv, eventFiles)){
        cout<<"All files connected"<<endl;
    }
    const string& outputFolder = argv[2];

    //setting IDs and other phrases to customise the histograms

    //histograms
    ProcessingOutsideLoop outsideprocessing;

    //histograms for efficiency
    outsideprocessing.AddHistogram(TH2D("K0STotal", "#eta vs p_{T} of K^{0}_{S};#eta;p_{T} [GeV/c]", 60, -3.0, 3.0, 40, 0, 2));
    outsideprocessing.AddHistogram(TH2D("K0SAfterMCpionCuts", "#eta vs p_{T} of K^{0}_{S} after MC fiducial cuts;#eta;p_{T} [GeV/c]", 60, -3.0, 3.0, 40, 0, 2));
    outsideprocessing.AddHistogram(TH2D("K0SAfterTPCpionCuts", "#eta vs p_{T} of K^{0}_{S} after MC cuts and TPC reconstruction;#eta;p_{T} [GeV/c]", 60, -3.0, 3.0, 40, 0, 2));

    //rest of the histograms
    outsideprocessing.AddHistogram(TH1D("ToFhitTimeDifferenceExtremelyWide", "Time difference between pair of ToF hits;#Delta t [ns];pairs", 100, -10, 10));
    outsideprocessing.AddHistogram(TH1D("ToFhitTimeDifferenceSimilarPositionExtremelyWide", "Time difference between pair of ToF hits to similar modules;#Delta t [ns];pairs", 100, -10, 10));
    outsideprocessing.AddHistogram(TH1D("ToFhitTimeDifferenceSamePositionExtremelyWide", "Time difference between pair of ToF hits to the same modules;#Delta t [ns];pairs", 100, -10, 10));
    outsideprocessing.AddHistogram(TH1D("ToFhitTimeDifference", "Time difference between pair of ToF hits;#Delta t [ns];pairs", 300, -1.5, 1.5));
    outsideprocessing.AddHistogram(TH1D("ToFhitTimeDifferenceSimilarPosition", "Time difference between pair of ToF hits to similar modules;#Delta t [ns];pairs", 300, -1.5, 1.5));
    outsideprocessing.AddHistogram(TH1D("ToFhitTimeDifferenceSamePosition", "Time difference between pair of ToF hits to the same modules;#Delta t [ns];pairs", 300, -1.5, 1.5));
    outsideprocessing.AddHistogram(TH2D("TracksVsTOFHits", "Number of tracks with kToF flag vs number of ToF hits;tracks;ToF hits", 50, 0, 50, 50, 0, 50));
    outsideprocessing.AddHistogram(TH1D("FlowOfEvents", "Independent event checks (of MC particles));;events", 1, 0, 1));
    outsideprocessing.AddHistogram(TH1D("K0SIndex", "Index of K^{0}_{S};index;particles", 30, 0, 30));
    outsideprocessing.AddHistogram(TH1D("K0Sdetectability", "How many K^{0}_{S} are detectable", 1, 0, 1));
    outsideprocessing.AddHistogram(TH2D("K0SdecayProductsKinematics", "#eta vs p_{T} of #pi^{#pm} from K^{0}_{S} decay;#eta;p_{T} [GeV/c]", 60, -3.0, 3.0, 40, 0, 2));
    outsideprocessing.AddHistogram(TH3D("NumberOfPions", "Number of pions in event;positive;negative;neutral", 5, 0, 5, 5, 0, 5, 5, 0, 5));
    outsideprocessing.AddHistogram(TH1D("MCParticles", "MC number of particles;particles;events", 25, 0, 25));
    outsideprocessing.AddHistogram(TH1D("TPCParticles", "TPC number of particles;particles;events", 25, 0, 25));
    outsideprocessing.AddHistogram(TH1D("DifferenceParticles", "(MC - TPC) difference in number of particles;difference;events", 20, -10, 10));
    outsideprocessing.AddHistogram(TH1D("Distance", "Distance MC-TPC in #eta-#phi space;distance;track pairs", 120, 0, 6));
    outsideprocessing.AddHistogram(TH1D("DistanceCloser", "Distance MC-TPC in #eta-#phi space;distance;track pairs", 100, 0, 0.5));
    outsideprocessing.AddHistogram(TH1D("MpipiTPC", "#pi^{+}#pi^{-} pair mass;m_{#pi^{+}#pi^{-}} [GeV/c^{2}];pairs", 40, 0.4, 0.6));
    outsideprocessing.AddHistogram(TH1D("MpipiTPCExtremelyWide", "#pi^{+}#pi^{-} pair mass;m_{#pi^{+}#pi^{-}} [GeV/c^{2}];pairs", 300, 0.0, 3.0));
    outsideprocessing.AddHistogram(TH1D("MpipiMC", "#pi^{+}#pi^{-} pair mass (only K^{0}_{S} decay products));m_{#pi^{+}#pi^{-}} [GeV/c^{2}];pairs", 40, 0.4, 0.6));
    outsideprocessing.AddHistogram(TH1D("MpipiMCExtremelyWide", "#pi^{+}#pi^{-} pair mass (only K^{0}_{S} decay products));m_{#pi^{+}#pi^{-}} [GeV/c^{2}];pairs", 300, 0.0, 3.0));
    outsideprocessing.AddHistogram(TH2D("MpipiMCExtremelyWideParticlesCreated", "Particles created in vertex vs #pi^{+}#pi^{-} pair mass (only K^{0}_{S} decay products));m_{#pi^{+}#pi^{-}} [GeV/c^{2}];pairs", 300, 0.0, 3.0, 30, 2, 32));
    outsideprocessing.AddHistogram(TH1D("MpipiMCExtremelyWideParticlesNames", "Particles created in vertex vs #pi^{+}#pi^{-} pair mass (only K^{0}_{S} decay products));m_{#pi^{+}#pi^{-}} [GeV/c^{2}];pairs", 1, 0, 1));
    outsideprocessing.AddHistogram(TH1D("MpipiFlow", "#pi^{+}#pi^{-} pairs after MC cuts;;pairs", 1, 0, 1));
    outsideprocessing.AddHistogram(TH2D("MpipiPairs", "particle pairs;positive;negative", 1, 0, 1, 1, 0, 1));
    outsideprocessing.AddHistogram(TH2D("MpipiMothers", "particle number;positive;negative", 20, 0, 20, 20, 0, 20));
    outsideprocessing.AddHistogram(TH1D("MpipiMotherName", "#pi^{+}#pi^{-} pair mothers;;mothers", 1, 0, 1));
    outsideprocessing.AddHistogram(TH1D("MpipiDeltaTExtremelyWide", "Time between pion decay;#Delta t_{0} [ns];events", 100, -5, 5));
    outsideprocessing.AddHistogram(TH1D("MpipiDeltaTFixCheckExtremelyWide", "Time correction between pion decay;#Delta t_{0} [ns];events", 100, -5, 5));
    outsideprocessing.AddHistogram(TH1D("MpipiDeltaTAfterFixExtremelyWide", "Corrected time between pion decay;#Delta t_{0} [ns];events", 100, -5, 5));
    outsideprocessing.AddHistogram(TH1D("MpipiDeltaTSimilarTOFModulesFixCheckExtremelyWide", "Time correction between pion decay;#Delta t_{0} [ns];events", 100, -5, 5));
    outsideprocessing.AddHistogram(TH1D("MpipiDeltaTSimilarTOFModulesAfterFixExtremelyWide", "Corrected time between pion decay;#Delta t_{0} [ns];events", 100, -5, 5));
    outsideprocessing.AddHistogram(TH1D("MpipiDeltaT", "Time between pion decay;#Delta t_{0} [ns];events", 300, -1.5, 1.5));
    outsideprocessing.AddHistogram(TH1D("MpipiDeltaTFixCheck", "Time correction between pion decay;#Delta t_{0} [ns];events", 300, -1.5, 1.5));
    outsideprocessing.AddHistogram(TH1D("MpipiDeltaTAfterFix", "Corrected time between pion decay;#Delta t_{0} [ns];events", 500, -2.5, 2.5));
    outsideprocessing.AddHistogram(TH1D("MpipiDeltaTSimilarTOFModulesFixCheck", "Time correction between pion decay;#Delta t_{0} [ns];events", 300, -1.5, 1.5));
    outsideprocessing.AddHistogram(TH1D("MpipiDeltaTSimilarTOFModulesAfterFix", "Corrected time between pion decay;#Delta t_{0} [ns];events", 500, -2.5, 2.5));
    outsideprocessing.AddHistogram(TH2D("MpipiDeltaTvsEvent", "Time between pion decay;#Delta t_{0} [ns];event number", 100, -5, 5, 1, 0, 1));
    outsideprocessing.AddHistogram(TH2D("DeltaTvsNumberOfMCParticles", "Time between pion decay vs number of MC particles;#Delta t_{0} [ns];particles", 100, -5, 5, 10, 0, 100));
    outsideprocessing.AddHistogram(TH1D("MpipiDeltaTOnlyTwoTracksDetected", "Time between pion decay;#Delta t_{0} [ns];events", 100, -5, 5));
    outsideprocessing.AddHistogram(TH1D("MpipiDeltaTMoreThanTwoTracksDetected", "Time between pion decay;#Delta t_{0} [ns];events", 100, -5, 5));

    outsideprocessing.AddHistogram(TH1D("MpipiAfterDeltaT", "#pi^{+}#pi^{-} pair mass after #Delta t_{0} cut;m_{#pi^{+}#pi^{-}} [GeV/c^{2}];pairs", 40, 0.4, 0.6));
    outsideprocessing.AddHistogram(TH1D("MpipiAfterDeltaTExtremelyWide", "#pi^{+}#pi^{-} pair mass after #Delta t_{0} cut;m_{#pi^{+}#pi^{-}} [GeV/c^{2}];pairs", 300, 0.0, 3.0));
    outsideprocessing.AddHistogram(TH2D("MpipiAfterDeltaTMass", "#pi^{+}#pi^{-} pair mass after #Delta t_{0} cut;m_{#pi^{+}#pi^{-}} [GeV/c^{2}];pairs", 40, 0.4, 0.6, 100, -5, 5));

    outsideprocessing.AddHistogram(TH1D("MpipiAfterDeltaTNotPassed", "#pi^{+}#pi^{-} pair mass after failed #Delta t_{0} cut;m_{#pi^{+}#pi^{-}} [GeV/c^{2}];pairs", 40, 0.4, 0.6));
    outsideprocessing.AddHistogram(TH1D("MpipiAfterDeltaTNotPassedExtremelyWide", "#pi^{+}#pi^{-} pair mass after failed #Delta t_{0} cut;m_{#pi^{+}#pi^{-}} [GeV/c^{2}];pairs", 300, 0.0, 3.0));
    outsideprocessing.AddHistogram(TH2D("MpipiAfterDeltaTNotPassedMass", "#pi^{+}#pi^{-} pair mass after failed #Delta t_{0} cut;m_{#pi^{+}#pi^{-}} [GeV/c^{2}];pairs", 40, 0.4, 0.6, 50, -5, 5));
    outsideprocessing.AddHistogram(TH1D("MpipiAfterDeltaTNotPassedLengthOfFlightPosX", "Flight x dimension", 120, -6, 6));
    outsideprocessing.AddHistogram(TH1D("MpipiAfterDeltaTNotPassedLengthOfFlightPosY", "Flight y dimension", 120, -6, 6));
    outsideprocessing.AddHistogram(TH1D("MpipiAfterDeltaTNotPassedLengthOfFlightPosZ", "Flight z dimension", 120, -6, 6));
    outsideprocessing.AddHistogram(TH1D("MpipiAfterDeltaTNotPassedLengthOfFlightTime", "Flight time", 100, 0, 1e-6));

    outsideprocessing.AddHistogram(TH1D("MpipiAfterDeltaTNotPassedLengthOfFlightPosDistance", "Flight distance", 100, 0, 10));
    outsideprocessing.AddHistogram(TH1D("MpipiAfterDeltaTPassedPionTOFLength", "TOF path length of #pi^{#pm} from K^{0}_{S} decay", 40, 200, 400));
    outsideprocessing.AddHistogram(TH1D("MpipiAfterDeltaTNotPassedPionTOFLength", "TOF path length of #pi^{#pm} from K^{0}_{S} decay", 40, 200, 400));
    outsideprocessing.AddHistogram(TH1D("MpipiAfterDeltaTPassedPionTOFTime", "TOF time of #pi^{#pm} from K^{0}_{S} decay", 100, 0, 50));
    outsideprocessing.AddHistogram(TH1D("MpipiAfterDeltaTNotPassedPionTOFTime", "TOF time of #pi^{#pm} from K^{0}_{S} decay", 100, 0, 50));
    outsideprocessing.AddHistogram(TH1D("MpipiAfterDeltaTPassedPionMass", "#pi^{#pm} pair mass from K^{0}_{S} decay;m_{#pi^{#pm}} [GeV/c^{2}];pairs", 100, 0., 1.0));
    outsideprocessing.AddHistogram(TH1D("MpipiAfterDeltaTNotPassedPionMass", "#pi^{#pm} pair mass from K^{0}_{S} decay;m_{#pi^{#pm}} [GeV/c^{2}];pairs", 100, 0., 1.0));
    outsideprocessing.AddHistogram(TH1D("MpipiAfterDeltaTPassedK0SDecayVertexBeamlineDistance", "Distance from K^{0}_{S} decay vertex to the beamline;d [cm];events", 50, 0, 20));
    outsideprocessing.AddHistogram(TH1D("MpipiAfterDeltaTNotPassedK0SDecayVertexBeamlineDistance", "Distance from K^{0}_{S} decay vertex to the beamline;d [cm];events", 50, 0, 20));
    outsideprocessing.AddHistogram(TH1D("MpipiAfterDeltaTPassedPionAngle", "#pi^{#pm} pair angle from K^{0}_{S} decay;angle [rad];pairs", 63, 0., 2*TMath::Pi()));
    outsideprocessing.AddHistogram(TH1D("MpipiAfterDeltaTNotPassedPionAngle", "#pi^{#pm} pair angle from K^{0}_{S} decay;angle [rad];pairs", 63, 0., 2*TMath::Pi()));
    outsideprocessing.AddHistogram(TH1D("TOFLeadingEdgeAfterDeltaTPassed", "TOF modules leading edge time;time [ns];number of modules", 100, 0, 20));
    outsideprocessing.AddHistogram(TH1D("TOFLeadingEdgeAfterDeltaTNotPassed", "TOF modules leading edge time;time [ns];number of modules", 100, 0, 20));


    outsideprocessing.AddHistogram(TH1D("MpipiMCnotK0SMother", "#pi^{+}#pi^{-} pair mass (everything except K^{0}_{S} mother verification));m_{#pi^{+}#pi^{-}} [GeV/c^{2}];pairs", 40, 0.4, 0.6));
    outsideprocessing.AddHistogram(TH1D("MpipiMCnotK0SMotherExtremelyWide", "#pi^{+}#pi^{-} pair mass (everything except K^{0}_{S} mother verification));m_{#pi^{+}#pi^{-}} [GeV/c^{2}];pairs", 300, 0.0, 3.0));

    outsideprocessing.AddHistogram(TH1D("DistanceAbnormal", "Distance MC-TPC in #eta-#phi space;distance;track pairs", 100, 0, 0.5));
    outsideprocessing.AddHistogram(TH1D("DistanceNormal", "Distance MC-TPC in #eta-#phi space;distance;track pairs", 100, 0, 0.5));
    outsideprocessing.AddHistogram(TH2D("TOFHitPositionAbnormal", "TOF hit position;z [cm];#phi [rad]", 40, -200, 200, 40, -TMath::Pi(), TMath::Pi()));
    outsideprocessing.AddHistogram(TH2D("TOFHitPositionNormal", "TOF hit position;z [cm];#phi [rad]", 40, -200, 200, 40, -TMath::Pi(), TMath::Pi()));
    outsideprocessing.AddHistogram(TH1D("TOFTimeAbnormal", "TOF time of #pi^{#pm} from K^{0}_{S} decay", 100, 0, 50));
    outsideprocessing.AddHistogram(TH1D("TOFTimeNormal", "TOF time of #pi^{#pm} from K^{0}_{S} decay", 100, 0, 50));
    outsideprocessing.AddHistogram(TH2D("TracksEtaPhiAbnormal", "#eta vs #phi of detected products;#eta;#phi", 40, -1.0, 1.0, 40, -TMath::Pi(), TMath::Pi()));
    outsideprocessing.AddHistogram(TH2D("TracksEtaPhiNormal", "#eta vs #phi of detected products;#eta;#phi", 40, -1.0, 1.0, 40, -TMath::Pi(), TMath::Pi()));
    outsideprocessing.AddHistogram(TH1D("MatchingAbnormal", "Matching of tracks using #eta-#phi space distance;;tracks", 1, 0, 1));
    outsideprocessing.AddHistogram(TH1D("MatchingNormal", "Matching of tracks using #eta-#phi space distance;;tracks", 1, 0, 1));
    outsideprocessing.GetPointer1D("MatchingAbnormal")->Fill("Correct", 0.0);
    outsideprocessing.GetPointer1D("MatchingAbnormal")->Fill("Incorrect", 0.0);
    outsideprocessing.GetPointer1D("MatchingNormal")->Fill("Correct", 0.0);
    outsideprocessing.GetPointer1D("MatchingNormal")->Fill("Incorrect", 0.0);
    outsideprocessing.AddHistogram(TH1D("MpipiDeltaTAllTPCTracks", "Time between pion decay;#Delta t_{0} [ns];events", 100, -5, 5));
    //TODO add tof cluster finding

    //processing
    //defining TreeProcessor
    ROOT::TTreeProcessorMT TreeProc(*eventFiles, nthreads);
    //defining processing function
    auto myFunction = [&](TTreeReader& myReader){
        //getting values from TChain, in-loop histogram initialization
        TTreeReaderValue<StUPCEvent> StUPCEventInstance(myReader, "mUPCEvent");
        // TTreeReaderValue<StRPEvent> StRPEventInstance(myReader, "mRPEvent");
        ProcessingInsideLoop insideprocessing;
        StUPCEvent* tempUPCpointer;
        // StRPEvent* tempRPpointer;
        insideprocessing.GetLocalHistograms(&outsideprocessing);

        std::vector<TParticle*> positiveMC;
        std::vector<TParticle*> negativeMC;
        std::vector<StUPCTrack*> positiveTrack;
        std::vector<StUPCTrack*> negativeTrack;
        std::vector<int> chosen_MC_list_positive;
        std::vector<int> chosen_MC_list_negative;
        std::vector<int> chosen_MC_list_copy;

        while(myReader.Next()){
            tempUPCpointer = StUPCEventInstance.Get();
            // tempRPpointer = StRPEventInstance.Get();

            //cleaning
            int tracksWithToFflag = 0;
            positiveMC.clear();
            negativeMC.clear();
            positiveTrack.clear();
            negativeTrack.clear();
            chosen_MC_list_positive.clear();
            chosen_MC_list_negative.clear();
            //below is the loop

            //filtering MC particles
            for(size_t i = 0; i<tempUPCpointer->getNumberOfMCParticles(); i++){
                TParticle* temp = tempUPCpointer->getMCParticle(i);
                //filtering
                if(fabs(temp->Eta())>0.9||temp->Pt()<0.2){
                    continue;
                }
                if(temp->GetPDG()->Charge()==0||abs(temp->GetPdgCode())<10){
                    //if not charged or gluon/quark
                    continue;
                }
                //filling
                if(temp->GetPDG()->Charge()>0){
                    positiveMC.push_back(temp);
                } else{
                    negativeMC.push_back(temp);
                }
            }

            //filtering tracks
            for(int i = 0; i<tempUPCpointer->getNumberOfTracks(); i++){
                StUPCTrack* tempTrack = tempUPCpointer->getTrack(i);
                //TODO check with & without TOF
                if(!tempTrack->getFlag(StUPCTrack::kTof)){
                    continue;
                }
                tracksWithToFflag++;
                if(tempTrack->getTofPathLength()<=0){
                    continue;
                }
                if(tempTrack->getTofTime()<=0){
                    continue;
                }
                // ALL TRACKS ARE PRIMARY
                // if(tempTrack->getFlag(StUPCTrack::kPrimary)){
                //     continue;
                // }
                //filtering
                if(fabs(tempTrack->getEta())>0.9||tempTrack->getPt()<0.2){
                    continue;
                }
                if(tempTrack->getNhitsFit()<=20||tempTrack->getNhitsDEdx()<15){
                    continue;
                }
                //filling
                if(tempTrack->getCharge()>0){
                    positiveTrack.push_back(tempTrack);
                } else{
                    negativeTrack.push_back(tempTrack);
                }
            }

            for(int i = 0; i<tempUPCpointer->getNumberOfHits()-1; i++){
                for(int j = i+1; j<tempUPCpointer->getNumberOfHits(); j++){
                    insideprocessing.Fill("ToFhitTimeDifference", tempUPCpointer->getHit(i)->getLeadingEdgeTime()-tempUPCpointer->getHit(j)->getLeadingEdgeTime());
                    insideprocessing.Fill("ToFhitTimeDifferenceExtremelyWide", tempUPCpointer->getHit(i)->getLeadingEdgeTime()-tempUPCpointer->getHit(j)->getLeadingEdgeTime());
                    int currentTray = tempUPCpointer->getHit(i)->getTray();
                    int currentModule = tempUPCpointer->getHit(i)->getModule();
                    int TOFmoduleDifference = abs(currentTray-tempUPCpointer->getHit(j)->getTray())+abs(currentModule-tempUPCpointer->getHit(j)->getModule());
                    if(TOFmoduleDifference==0){
                        insideprocessing.Fill("ToFhitTimeDifferenceSamePosition", tempUPCpointer->getHit(i)->getLeadingEdgeTime()-tempUPCpointer->getHit(j)->getLeadingEdgeTime());
                        insideprocessing.Fill("ToFhitTimeDifferenceSamePositionExtremelyWide", tempUPCpointer->getHit(i)->getLeadingEdgeTime()-tempUPCpointer->getHit(j)->getLeadingEdgeTime());
                    }
                    if(TOFmoduleDifference==1){
                        insideprocessing.Fill("ToFhitTimeDifferenceSimilarPosition", tempUPCpointer->getHit(i)->getLeadingEdgeTime()-tempUPCpointer->getHit(j)->getLeadingEdgeTime());
                        insideprocessing.Fill("ToFhitTimeDifferenceSimilarPositionExtremelyWide", tempUPCpointer->getHit(i)->getLeadingEdgeTime()-tempUPCpointer->getHit(j)->getLeadingEdgeTime());
                    }
                }
            }

            bool duplicated_matches_exist = false;

            //histograms
            insideprocessing.Fill("TracksVsTOFHits", tracksWithToFflag, tempUPCpointer->getNumberOfHits());
            //early histogram filling, for proper order
            insideprocessing.Fill("FlowOfEvents", "All events", 1.0);
            insideprocessing.Fill("FlowOfEvents", "2x K^{0}_{S}", 0.0);
            insideprocessing.Fill("FlowOfEvents", "4x#pi^{#pm/0}", 0.0);
            insideprocessing.Fill("FlowOfEvents", "2x#pi^{+}, 2x#pi^{-}", 0.0);
            insideprocessing.Fill("K0Sdetectability", "All K^{0}_{S}", 0.0);
            insideprocessing.Fill("K0Sdetectability", "Decayed into #pi^{+}#pi^{-} pair", 0.0);
            insideprocessing.Fill("K0Sdetectability", "Both MC pions in p_{T}-#eta range", 0.0);
            int n_of_K0S = 0;
            int n_of_piplus = 0;
            int n_of_piminus = 0;
            int n_of_pizero = 0;
            for(size_t i = 0; i<tempUPCpointer->getNumberOfMCParticles(); i++){
                if(abs(tempUPCpointer->getMCParticle(i)->GetPdgCode())==PDGmain){
                    n_of_K0S++;
                    insideprocessing.Fill("K0SIndex", i);
                    insideprocessing.Fill("K0STotal", tempUPCpointer->getMCParticle(i)->Eta(), tempUPCpointer->getMCParticle(i)->Pt());
                    //K0Sdetectability things
                    insideprocessing.Fill("K0Sdetectability", 0., 1.0);
                    int decayVertex = tempUPCpointer->getMCParticle(i)->GetFirstDaughter();
                    bool hasPiPlus = false, hasPiMinus = false;
                    bool isPiPlusDetectable = false, isPiMinusDetectable = false;
                    for(size_t j = 0; j<tempUPCpointer->getNumberOfMCParticles(); j++){
                        TParticle* temp = tempUPCpointer->getMCParticle(j);
                        if(temp->GetFirstMother()!=decayVertex){
                            continue;
                        }
                        //when we check we actually have a particle correlated
                        //with the proper vertex, we can continue analysis
                        if(temp->GetPdgCode()==PDGpositive){
                            hasPiPlus = true;
                            if(fabs(temp->Eta())<=0.9&&temp->Pt()>=0.2){
                                isPiPlusDetectable = true;
                            }
                            insideprocessing.Fill("K0SdecayProductsKinematics", temp->Eta(), temp->Pt());
                        }
                        if(temp->GetPdgCode()==PDGnegative){
                            hasPiMinus = true;
                            if(fabs(temp->Eta())<=0.9&&temp->Pt()>=0.2){
                                isPiMinusDetectable = true;
                            }
                            insideprocessing.Fill("K0SdecayProductsKinematics", temp->Eta(), temp->Pt());
                        }
                    }
                    if(hasPiPlus&&hasPiMinus){
                        insideprocessing.Fill("K0Sdetectability", 1, 1.0);
                    }
                    if(isPiPlusDetectable&&isPiMinusDetectable){
                        insideprocessing.Fill("K0Sdetectability", 2, 1.0);
                        insideprocessing.Fill("K0SAfterMCpionCuts", tempUPCpointer->getMCParticle(i)->Eta(), tempUPCpointer->getMCParticle(i)->Pt());
                    }
                }
                if(tempUPCpointer->getMCParticle(i)->GetPdgCode()==PDGpositive)
                    n_of_piplus++;
                if(tempUPCpointer->getMCParticle(i)->GetPdgCode()==PDGnegative)
                    n_of_piminus++;
                if(abs(tempUPCpointer->getMCParticle(i)->GetPdgCode())==111)//pi0 special case
                    n_of_pizero++;
            }
            insideprocessing.Fill("NumberOfPions", n_of_piplus, n_of_piminus, n_of_pizero);
            if(n_of_K0S==2)
                insideprocessing.Fill("FlowOfEvents", "2x K^{0}_{S}", 1.0);
            if(n_of_piplus==2&&n_of_piminus==2)
                insideprocessing.Fill("FlowOfEvents", "2x#pi^{+}, 2x#pi^{-}", 1.0);
            if(n_of_piplus+n_of_piminus+n_of_pizero==4)
                insideprocessing.Fill("FlowOfEvents", "4x#pi^{#pm/0}", 1.0);
            //the rest
            insideprocessing.Fill("MCParticles", positiveMC.size()+negativeMC.size());
            insideprocessing.Fill("TPCParticles", positiveTrack.size()+negativeTrack.size());
            insideprocessing.Fill("DifferenceParticles", int(positiveMC.size()+negativeMC.size())-int(positiveTrack.size()+negativeTrack.size()));
            //loops for fitting the tracks
            if(myReader.GetCurrentEntry()==eventToPrint){
                printf("Event %d:\n", myReader.GetCurrentEntry());
                printf("Positive:\n");
            }
            //positive
            for(size_t i = 0; i<positiveTrack.size(); i++){
                double min_distance = 100;
                int min_distance_id = -1;
                for(size_t j = 0; j<positiveMC.size(); j++){
                    //getting distance in eta-phi space
                    TVector3 vecMC, vecTPC;
                    vecMC.SetXYZ(positiveMC[j]->Px(), positiveMC[j]->Py(), positiveMC[j]->Pz());
                    positiveTrack[i]->getMomentum(vecTPC);
                    double etaphi_distance = vecMC.DrEtaPhi(vecTPC);
                    //distance cutoff
                    if(etaphi_distance>0.1){
                        continue;
                    }
                    insideprocessing.Fill("Distance", etaphi_distance);
                    insideprocessing.Fill("DistanceCloser", etaphi_distance);
                    //checking if the distance to j-th MC track is the smallest
                    min_distance = min(min_distance, etaphi_distance);
                    if(min_distance==etaphi_distance){
                        min_distance_id = j;
                    }
                }
                //checking for which MC track we have the nearest match
                chosen_MC_list_positive.push_back(min_distance_id);
                if(myReader.GetCurrentEntry()==eventToPrint){
                    //finding global index of MC particle
                    for(size_t MC_particle_index = 0; MC_particle_index<tempUPCpointer->getNumberOfMCParticles(); MC_particle_index++){
                        if(tempUPCpointer->getMCParticle(MC_particle_index)==positiveMC[min_distance_id]){
                            printf("TPC TruthId:\t%d,\t\tMC array index:\t%d\n", positiveTrack[i]->getIdTruth(), MC_particle_index);
                        }
                    }
                }
            }
            chosen_MC_list_copy = chosen_MC_list_positive;
            //removing all the -1
            chosen_MC_list_copy.erase(std::remove(chosen_MC_list_copy.begin(), chosen_MC_list_copy.end(), -1), chosen_MC_list_copy.end());
            std::sort(chosen_MC_list_copy.begin(), chosen_MC_list_copy.end());
            const auto duplicate_positive = std::adjacent_find(chosen_MC_list_copy.begin(), chosen_MC_list_copy.end());
            if(duplicate_positive!=chosen_MC_list_copy.end()){
                printf("Event with conflict (positive) = %d\n", myReader.GetCurrentEntry());
                for(size_t index = 0; index<chosen_MC_list_positive.size(); index++){
                    printf("%d, ", chosen_MC_list_positive[index]);
                }
                printf("\n");
                duplicated_matches_exist = true;
            }
            if(myReader.GetCurrentEntry()==eventToPrint){
                printf("Negative:\n");
            }
            //negative
            for(size_t i = 0; i<negativeTrack.size(); i++){
                double min_distance = 100;
                int min_distance_id = -1;
                for(size_t j = 0; j<negativeMC.size(); j++){
                    //getting distance in eta-phi space
                    TVector3 vecMC, vecTPC;
                    vecMC.SetXYZ(negativeMC[j]->Px(), negativeMC[j]->Py(), negativeMC[j]->Pz());
                    negativeTrack[i]->getMomentum(vecTPC);
                    double etaphi_distance = vecMC.DrEtaPhi(vecTPC);
                    //distance cutoff
                    if(etaphi_distance>0.1){
                        continue;
                    }
                    insideprocessing.Fill("Distance", etaphi_distance);
                    insideprocessing.Fill("DistanceCloser", etaphi_distance);
                    //checking if the distance to j-th MC track is the smallest
                    min_distance = min(min_distance, etaphi_distance);
                    if(min_distance==etaphi_distance){
                        min_distance_id = j;
                    }
                }
                //checking for which MC track we have the nearest match
                chosen_MC_list_negative.push_back(min_distance_id);
                if(myReader.GetCurrentEntry()==eventToPrint){
                    //finding global index of MC particle
                    for(size_t MC_particle_index = 0; MC_particle_index<tempUPCpointer->getNumberOfMCParticles(); MC_particle_index++){
                        if(tempUPCpointer->getMCParticle(MC_particle_index)==negativeMC[min_distance_id]){
                            printf("TPC TruthId:\t%d,\t\tMC array index:\t%d\n", negativeTrack[i]->getIdTruth(), MC_particle_index);
                        }
                    }
                }
            }
            chosen_MC_list_copy = chosen_MC_list_negative;
            //removing all the -1
            chosen_MC_list_copy.erase(std::remove(chosen_MC_list_copy.begin(), chosen_MC_list_copy.end(), -1), chosen_MC_list_copy.end());
            std::sort(chosen_MC_list_copy.begin(), chosen_MC_list_copy.end());
            const auto duplicate_negative = std::adjacent_find(chosen_MC_list_copy.begin(), chosen_MC_list_copy.end());
            if(duplicate_negative!=chosen_MC_list_copy.end()){
                printf("Event with conflict (negative) = %d\n", myReader.GetCurrentEntry());
                for(size_t index = 0; index<chosen_MC_list_negative.size(); index++){
                    printf("%d, ", chosen_MC_list_negative[index]);
                }
                printf("\n");
                duplicated_matches_exist = true;
            }

            //checking the mass-fit of pairs when no duplicates
            if(!duplicated_matches_exist){
                for(size_t i = 0; i<positiveTrack.size(); i++){
                    for(size_t j = 0; j<negativeTrack.size(); j++){
                        //check if the track is associated with anything
                        if(chosen_MC_list_positive[i]<0||chosen_MC_list_negative[j]<0){
                            continue;
                        }
                        //tracks
                        TLorentzVector posTrack, negTrack;
                        positiveTrack[i]->getLorentzVector(posTrack, massmap[PDGpositive]);
                        negativeTrack[j]->getLorentzVector(negTrack, massmap[PDGnegative]);
                        insideprocessing.Fill("MpipiTPC", (posTrack+negTrack).M());
                        insideprocessing.Fill("MpipiTPCExtremelyWide", (posTrack+negTrack).M());
                        //MC particles
                        //choosing MC particle associated with these particular tracks
                        int posPDG = positiveMC[chosen_MC_list_positive[i]]->GetPdgCode();
                        int posProductionVertex = positiveMC[chosen_MC_list_positive[i]]->GetFirstMother();
                        int negPDG = negativeMC[chosen_MC_list_negative[j]]->GetPdgCode();
                        int negProductionVertex = negativeMC[chosen_MC_list_negative[j]]->GetFirstMother();
                        //loop for finding mother particle
                        int posMotherPDG = 0;
                        int MotherParticleIndex = -1;
                        for(size_t MCindex = 0; MCindex<tempUPCpointer->getNumberOfMCParticles(); MCindex++){
                            if(tempUPCpointer->getMCParticle(MCindex)->GetFirstDaughter()==posProductionVertex){
                                posMotherPDG = tempUPCpointer->getMCParticle(MCindex)->GetPdgCode();
                                MotherParticleIndex = MCindex;
                                break;
                            }
                        }
                        //both pions, same mother, mother is K0S
                        insideprocessing.Fill("MpipiFlow", "TPC", 1.0);
                        insideprocessing.Fill("MpipiPairs", positiveMC[chosen_MC_list_positive[i]]->GetPDG()->GetName(), negativeMC[chosen_MC_list_negative[j]]->GetPDG()->GetName(), 1.0);
                        if(posPDG==PDGpositive&&negPDG==PDGnegative){
                            insideprocessing.Fill("MpipiFlow", "Pion pair", 1.0);
                            insideprocessing.Fill("MpipiMothers", posProductionVertex, negProductionVertex);
                            if(posProductionVertex==negProductionVertex){
                                insideprocessing.Fill("MpipiFlow", "Same mother", 1.0);
                                for(size_t MCindex = 0; MCindex<tempUPCpointer->getNumberOfMCParticles(); MCindex++){
                                    if(tempUPCpointer->getMCParticle(MCindex)->GetFirstDaughter()==posProductionVertex){
                                        insideprocessing.Fill("MpipiMotherName", tempUPCpointer->getMCParticle(MCindex)->GetPDG()->GetName(), 1.0);
                                        break;
                                    }
                                }
                                if(abs(posMotherPDG)==PDGmain){
                                    insideprocessing.Fill("MpipiFlow", "K^{0}_{S} mother", 1.0);
                                    insideprocessing.Fill("MpipiMC", (posTrack+negTrack).M());
                                    insideprocessing.Fill("MpipiMCExtremelyWide", (posTrack+negTrack).M());
                                    int particlesFromVertex = 0;
                                    for(size_t MCparticle = 0; MCparticle<tempUPCpointer->getNumberOfMCParticles(); MCparticle++){
                                        if(tempUPCpointer->getMCParticle(MCparticle)->GetFirstMother()==posProductionVertex){
                                            particlesFromVertex++;
                                        }
                                    }
                                    insideprocessing.Fill("MpipiMCExtremelyWideParticlesCreated", (posTrack+negTrack).M(), particlesFromVertex);
                                    if(particlesFromVertex>2){
                                        for(size_t MCparticle = 0; MCparticle<tempUPCpointer->getNumberOfMCParticles(); MCparticle++){
                                            if(tempUPCpointer->getMCParticle(MCparticle)->GetFirstMother()==posProductionVertex){
                                                insideprocessing.Fill("MpipiMCExtremelyWideParticlesNames", tempUPCpointer->getMCParticle(MCparticle)->GetPDG()->GetName(), 1.);
                                            }
                                        }
                                    }
                                    insideprocessing.Fill("K0SAfterTPCpionCuts", tempUPCpointer->getMCParticle(MotherParticleIndex)->Eta(), tempUPCpointer->getMCParticle(MotherParticleIndex)->Pt());
                                    //stuff to do with deltaT
                                    double deltaT = DeltaT0(positiveTrack[i], negativeTrack[j], massmap[PDGpositive], massmap[PDGnegative]);
                                    insideprocessing.Fill("MpipiDeltaTExtremelyWide", deltaT);
                                    insideprocessing.Fill("MpipiDeltaTAfterFixExtremelyWide", DeltaTFix(deltaT, positiveMC[chosen_MC_list_positive[i]], negativeMC[chosen_MC_list_negative[j]], tempUPCpointer));
                                    insideprocessing.Fill("MpipiDeltaTFixCheckExtremelyWide", DeltaTFix(deltaT, positiveMC[chosen_MC_list_positive[i]], negativeMC[chosen_MC_list_negative[j]], tempUPCpointer)-deltaT);
                                    insideprocessing.Fill("MpipiDeltaTSimilarTOFModulesAfterFixExtremelyWide", DeltaTFixSimilarModules(deltaT, positiveMC[chosen_MC_list_positive[i]], negativeMC[chosen_MC_list_negative[j]], tempUPCpointer));
                                    insideprocessing.Fill("MpipiDeltaTSimilarTOFModulesFixCheckExtremelyWide", DeltaTFixSimilarModules(deltaT, positiveMC[chosen_MC_list_positive[i]], negativeMC[chosen_MC_list_negative[j]], tempUPCpointer)-deltaT);
                                    insideprocessing.Fill("MpipiDeltaT", deltaT);
                                    insideprocessing.Fill("MpipiDeltaTAfterFix", DeltaTFix(deltaT, positiveMC[chosen_MC_list_positive[i]], negativeMC[chosen_MC_list_negative[j]], tempUPCpointer));
                                    insideprocessing.Fill("MpipiDeltaTFixCheck", DeltaTFix(deltaT, positiveMC[chosen_MC_list_positive[i]], negativeMC[chosen_MC_list_negative[j]], tempUPCpointer)-deltaT);
                                    insideprocessing.Fill("MpipiDeltaTSimilarTOFModulesAfterFix", DeltaTFixSimilarModules(deltaT, positiveMC[chosen_MC_list_positive[i]], negativeMC[chosen_MC_list_negative[j]], tempUPCpointer));
                                    insideprocessing.Fill("MpipiDeltaTSimilarTOFModulesFixCheck", DeltaTFixSimilarModules(deltaT, positiveMC[chosen_MC_list_positive[i]], negativeMC[chosen_MC_list_negative[j]], tempUPCpointer)-deltaT);
                                    insideprocessing.Fill("MpipiDeltaTvsEvent", deltaT, to_string(myReader.GetCurrentEntry()).c_str(), 1.);
                                    if(positiveTrack.size()==1&&negativeTrack.size()==1){
                                        insideprocessing.Fill("MpipiDeltaTOnlyTwoTracksDetected", deltaT);
                                    } else{
                                        //if no tracks either positive or negative
                                        //then no loops, so the only alternative to 1&1
                                        //is more tracks, so no checking necessary
                                        insideprocessing.Fill("MpipiDeltaTMoreThanTwoTracksDetected", deltaT);
                                    }
                                    insideprocessing.Fill("DeltaTvsNumberOfMCParticles", deltaT, tempUPCpointer->getNumberOfMCParticles());
                                    //value nicked from .txt file with actual data
                                    //replaced with a better cut of 0.6
                                    // if(fabs(deltaT)<3*0.13124272289253383){
                                    if(fabs(deltaT)<0.6){
                                        insideprocessing.Fill("MpipiAfterDeltaT", (posTrack+negTrack).M());
                                        insideprocessing.Fill("MpipiAfterDeltaTExtremelyWide", (posTrack+negTrack).M());
                                        insideprocessing.Fill("MpipiAfterDeltaTMass", (posTrack+negTrack).M(), deltaT);
                                        insideprocessing.Fill("MpipiAfterDeltaTPassedPionTOFLength", positiveTrack[i]->getTofPathLength());
                                        insideprocessing.Fill("MpipiAfterDeltaTPassedPionTOFLength", negativeTrack[j]->getTofPathLength());
                                        insideprocessing.Fill("MpipiAfterDeltaTPassedPionTOFTime", positiveTrack[i]->getTofTime());
                                        insideprocessing.Fill("MpipiAfterDeltaTPassedPionTOFTime", negativeTrack[j]->getTofTime());
                                        insideprocessing.Fill("MpipiAfterDeltaTPassedPionMass", sqrt(M2TOF(positiveTrack[i], negativeTrack[j])));
                                        TLorentzVector K0SdecayVertex;
                                        positiveMC[chosen_MC_list_positive[i]]->ProductionVertex(K0SdecayVertex);
                                        insideprocessing.Fill("MpipiAfterDeltaTPassedK0SDecayVertexBeamlineDistance", distanceToBeamline(K0SdecayVertex.Vect()));
                                        insideprocessing.Fill("MpipiAfterDeltaTPassedPionAngle", posTrack.Angle(negTrack.Vect()));
                                        for(size_t TOFhit = 0; TOFhit<tempUPCpointer->getNumberOfHits(); TOFhit++){
                                            insideprocessing.Fill("TOFLeadingEdgeAfterDeltaTPassed", tempUPCpointer->getHit(TOFhit)->getLeadingEdgeTime());
                                        }
                                    } else{
                                        insideprocessing.Fill("MpipiAfterDeltaTNotPassed", (posTrack+negTrack).M());
                                        insideprocessing.Fill("MpipiAfterDeltaTNotPassedExtremelyWide", (posTrack+negTrack).M());
                                        insideprocessing.Fill("MpipiAfterDeltaTNotPassedMass", (posTrack+negTrack).M(), deltaT);
                                        TParticle* mother = tempUPCpointer->getMCParticle(MotherParticleIndex);
                                        TLorentzVector MotherProductionVertex, MotherDecayVertex;
                                        mother->ProductionVertex(MotherProductionVertex);
                                        //mother decay vertex is daughter production vertex
                                        positiveMC[chosen_MC_list_positive[i]]->ProductionVertex(MotherDecayVertex);
                                        TLorentzVector Length = MotherDecayVertex-MotherProductionVertex;
                                        insideprocessing.Fill("MpipiAfterDeltaTNotPassedLengthOfFlightPosX", Length.X());
                                        insideprocessing.Fill("MpipiAfterDeltaTNotPassedLengthOfFlightPosY", Length.Y());
                                        insideprocessing.Fill("MpipiAfterDeltaTNotPassedLengthOfFlightPosZ", Length.Z());
                                        insideprocessing.Fill("MpipiAfterDeltaTNotPassedLengthOfFlightTime", Length.T());
                                        insideprocessing.Fill("MpipiAfterDeltaTNotPassedLengthOfFlightPosDistance", Length.P());
                                        insideprocessing.Fill("MpipiAfterDeltaTNotPassedPionTOFLength", positiveTrack[i]->getTofPathLength());
                                        insideprocessing.Fill("MpipiAfterDeltaTNotPassedPionTOFLength", negativeTrack[j]->getTofPathLength());
                                        insideprocessing.Fill("MpipiAfterDeltaTNotPassedPionTOFTime", positiveTrack[i]->getTofTime());
                                        insideprocessing.Fill("MpipiAfterDeltaTNotPassedPionTOFTime", negativeTrack[j]->getTofTime());
                                        insideprocessing.Fill("MpipiAfterDeltaTNotPassedPionMass", sqrt(M2TOF(positiveTrack[i], negativeTrack[j])));
                                        TLorentzVector K0SdecayVertex;
                                        positiveMC[chosen_MC_list_positive[i]]->ProductionVertex(K0SdecayVertex);
                                        insideprocessing.Fill("MpipiAfterDeltaTNotPassedK0SDecayVertexBeamlineDistance", distanceToBeamline(K0SdecayVertex.Vect()));
                                        insideprocessing.Fill("MpipiAfterDeltaTNotPassedPionAngle", posTrack.Angle(negTrack.Vect()));
                                        for(size_t TOFhit = 0; TOFhit<tempUPCpointer->getNumberOfHits(); TOFhit++){
                                            insideprocessing.Fill("TOFLeadingEdgeAfterDeltaTNotPassed", tempUPCpointer->getHit(TOFhit)->getLeadingEdgeTime());
                                        }
                                        //part for checking separate tracks
                                        //setting temp values to not switch them around
                                        //assumption is that to abnormal tracks was added 1ns, but it might be wrong
                                        //and as we always have t0+ - t0-, then for deltaT>0 the abnormal track is the positive one
                                        StUPCTrack* normalTrack, * abnormalTrack;
                                        TLorentzVector normalTrackVector, abnormalTrackVector;
                                        TParticle* normalMCparticle, * abnormalMCParticle;
                                        if(deltaT>0){
                                            abnormalTrack = positiveTrack[i];
                                            normalTrack = negativeTrack[j];
                                            abnormalTrackVector = posTrack;
                                            normalTrackVector = negTrack;
                                            abnormalMCParticle = positiveMC[chosen_MC_list_positive[i]];
                                            normalMCparticle = negativeMC[chosen_MC_list_negative[j]];
                                        } else{
                                            abnormalTrack = negativeTrack[j];
                                            normalTrack = positiveTrack[i];
                                            abnormalTrackVector = negTrack;
                                            normalTrackVector = posTrack;
                                            abnormalMCParticle = negativeMC[chosen_MC_list_negative[j]];
                                            normalMCparticle = positiveMC[chosen_MC_list_positive[i]];
                                        }
                                        //histograms of distance in eta-phi space
                                        TVector3 vecMCTempAbnormal, vecMCTempNormal;
                                        vecMCTempAbnormal.SetXYZ(abnormalMCParticle->Px(), abnormalMCParticle->Py(), abnormalMCParticle->Pz());
                                        vecMCTempNormal.SetXYZ(normalMCparticle->Px(), normalMCparticle->Py(), normalMCparticle->Pz());
                                        insideprocessing.Fill("DistanceAbnormal", abnormalTrackVector.Vect().DrEtaPhi(vecMCTempAbnormal));
                                        insideprocessing.Fill("DistanceNormal", normalTrackVector.Vect().DrEtaPhi(vecMCTempNormal));
                                        //histograms of ToF cluster position
                                        double middleOfToF = 214.4;
                                        StPicoPhysicalHelix p1Helix = StPicoPhysicalHelix(abnormalTrack->getCurvature(),
                                            abnormalTrack->getDipAngle(),
                                            abnormalTrack->getPhase(),
                                            abnormalTrack->getOrigin(),
                                            abnormalTrack->getCharge());
                                        StPicoPhysicalHelix p2Helix = StPicoPhysicalHelix(normalTrack->getCurvature(),
                                            normalTrack->getDipAngle(),
                                            normalTrack->getPhase(),
                                            normalTrack->getOrigin(),
                                            normalTrack->getCharge());
                                        std::pair<double, double> stemp = p1Helix.pathLength(middleOfToF);
                                        double s1 = (stemp.first<=0) ? stemp.second : min(stemp.first, stemp.second);
                                        TVector3 TOFhitPosition = p1Helix.at(s1);
                                        insideprocessing.Fill("TOFHitPositionAbnormal", TOFhitPosition.Z(), TOFhitPosition.Phi());
                                        stemp = p2Helix.pathLength(middleOfToF);
                                        double s2 = (stemp.first<=0) ? stemp.second : min(stemp.first, stemp.second);
                                        TOFhitPosition = p2Helix.at(s2);
                                        insideprocessing.Fill("TOFHitPositionNormal", TOFhitPosition.Z(), TOFhitPosition.Phi());
                                        //the rest of the histograms
                                        insideprocessing.Fill("TOFTimeAbnormal", abnormalTrack->getTofTime());
                                        insideprocessing.Fill("TOFTimeNormal", normalTrack->getTofTime());
                                        insideprocessing.Fill("TracksEtaPhiAbnormal", vecMCTempAbnormal.Eta(), vecMCTempAbnormal.Phi());
                                        insideprocessing.Fill("TracksEtaPhiNormal", vecMCTempNormal.Eta(), vecMCTempNormal.Phi());
                                        for(size_t MC_particle_index = 0; MC_particle_index<tempUPCpointer->getNumberOfMCParticles(); MC_particle_index++){
                                            //looking for match for abnormal track
                                            if(tempUPCpointer->getMCParticle(MC_particle_index)==abnormalMCParticle){
                                                if(abnormalTrack->getIdTruth()==MC_particle_index+1){
                                                    insideprocessing.Fill("MatchingAbnormal", "Correct", 1.0);
                                                } else{
                                                    insideprocessing.Fill("MatchingAbnormal", "Incorrect", 1.0);
                                                }
                                            }
                                            //looking for match for normal track
                                            if(tempUPCpointer->getMCParticle(MC_particle_index)==normalMCparticle){
                                                if(normalTrack->getIdTruth()==MC_particle_index+1){
                                                    insideprocessing.Fill("MatchingNormal", "Correct", 1.0);
                                                } else{
                                                    insideprocessing.Fill("MatchingNormal", "Incorrect", 1.0);
                                                }
                                            }
                                        }
                                    }
                                } else{
                                    insideprocessing.Fill("MpipiMCnotK0SMother", (posTrack+negTrack).M());
                                    insideprocessing.Fill("MpipiMCnotK0SMotherExtremelyWide", (posTrack+negTrack).M());
                                }
                            }
                        }
                    }
                }
            }

            //checking things for all TPC tracks
            for(size_t i = 0; i<positiveTrack.size(); i++){
                for(size_t j = 0; j<negativeTrack.size(); j++){
                    TLorentzVector posTrack, negTrack;
                    positiveTrack[i]->getLorentzVector(posTrack, massmap[PDGpositive]);
                    negativeTrack[j]->getLorentzVector(negTrack, massmap[PDGnegative]);
                    double deltaT = DeltaT0(positiveTrack[i], negativeTrack[j], massmap[PDGpositive], massmap[PDGnegative]);
                    //filling histograms
                    insideprocessing.Fill("MpipiDeltaTAllTPCTracks", deltaT);
                }
            }


            //special part where one event is drawn
            if(myReader.GetCurrentEntry()==eventToPrint){
                //statistics
                printf("Number of MC particles (excluding protons) before filter: %d\n", tempUPCpointer->getNumberOfMCParticles()-2);
                printf("Number of MC particles (excluding protons) after filter: %d\n", positiveMC.size()+negativeMC.size());
                printf("Number of tracks (excluding protons) before filter: %d\n", tempUPCpointer->getNumberOfTracks());
                printf("Number of tracks (excluding protons) after filter: %d\n", positiveTrack.size()+negativeTrack.size());
                for(size_t MCindex = 0; MCindex<tempUPCpointer->getNumberOfMCParticles(); MCindex++){
                    //we only know IdTruth because there is great care taken not to change the order of particles
                    //it is NOT written into upcDst file!!!
                    PrintBigger(tempUPCpointer->getMCParticle(MCindex), "\tIdTruth:\t"+to_string(MCindex+1));
                }


                //drawing
                TCanvas c1("c1", "c1", 1200, 800);

                //MC particles graph
                TGraph ParticlesMC(0);
                ParticlesMC.SetNameTitle("MC", "MC;#eta;#phi");
                ParticlesMC.SetMarkerStyle(20);
                ParticlesMC.SetMarkerSize(2);
                ParticlesMC.SetMarkerColor(4);
                printf("MC positive:\n");
                for(size_t particle_index = 0; particle_index<positiveMC.size(); particle_index++){
                    TParticle* temp = positiveMC[particle_index];
                    TVector3 tempVec(temp->Px(), temp->Py(), temp->Pz());
                    ParticlesMC.AddPoint(tempVec.Eta(), tempVec.Phi());
                    printf("Eta:\t%f,\tPhi:\t%f\n", tempVec.Eta(), tempVec.Phi());
                }
                printf("MC negative:\n");
                for(size_t particle_index = 0; particle_index<negativeMC.size(); particle_index++){
                    TParticle* temp = negativeMC[particle_index];
                    TVector3 tempVec(temp->Px(), temp->Py(), temp->Pz());
                    ParticlesMC.AddPoint(tempVec.Eta(), tempVec.Phi());
                    printf("Eta:\t%f,\tPhi:\t%f\n", tempVec.Eta(), tempVec.Phi());
                }

                //TPC particles graph
                TGraph ParticlesTPC(0);
                ParticlesTPC.SetNameTitle("TPC", "TPC;#eta;#phi");
                ParticlesTPC.SetMarkerStyle(21);
                ParticlesTPC.SetMarkerSize(1.5);
                ParticlesTPC.SetMarkerColor(2);
                printf("Track positive:\n");
                for(int particle_index = 0; particle_index<positiveTrack.size(); particle_index++){
                    StUPCTrack* tempTrack = positiveTrack[particle_index];
                    TVector3 tempVec;
                    tempTrack->getMomentum(tempVec);
                    ParticlesTPC.AddPoint(tempVec.Eta(), tempVec.Phi());
                    printf("Eta:\t%f,\tPhi:\t%f\n", tempVec.Eta(), tempVec.Phi());
                }
                printf("Track negative:\n");
                for(int particle_index = 0; particle_index<negativeTrack.size(); particle_index++){
                    StUPCTrack* tempTrack = negativeTrack[particle_index];
                    TVector3 tempVec;
                    tempTrack->getMomentum(tempVec);
                    ParticlesTPC.AddPoint(tempVec.Eta(), tempVec.Phi());
                    printf("Eta:\t%f,\tPhi:\t%f\n", tempVec.Eta(), tempVec.Phi());
                }

                //drawing and saving canvas (with protection against empty graphs)
                bool drawnMC = false;
                if(ParticlesMC.GetN()){
                    ParticlesMC.Draw("ap");
                    ParticlesMC.GetXaxis()->SetLimits(-1.0, 1.0);
                    ParticlesMC.GetHistogram()->SetMinimum(-TMath::Pi());
                    ParticlesMC.GetHistogram()->SetMaximum(TMath::Pi());
                    drawnMC = true;
                }
                if(ParticlesTPC.GetN()){
                    if(drawnMC){
                        ParticlesTPC.Draw("same p");
                    } else{
                        ParticlesTPC.Draw("ap");
                        ParticlesTPC.GetXaxis()->SetLimits(-1.0, 1.0);
                        ParticlesTPC.GetHistogram()->SetMinimum(-TMath::Pi());
                        ParticlesTPC.GetHistogram()->SetMaximum(TMath::Pi());
                    }
                }
                c1.BuildLegend();
                c1.SetTitle("Matching test");
                c1.SaveAs("Matching.png");
            }

        }
        //event loop finish

        //lambda finish
        return 0;
    };

    TreeProc.Process(myFunction);

    //merging and tidying up
    outsideprocessing.Merge();
    outsideprocessing.GetPointerAfterMerge1D("FlowOfEvents")->LabelsDeflate();
    outsideprocessing.GetPointerAfterMerge1D("K0Sdetectability")->LabelsDeflate();
    outsideprocessing.GetPointerAfterMerge1D("MpipiMCExtremelyWideParticlesNames")->LabelsDeflate();
    outsideprocessing.GetPointerAfterMerge1D("MpipiFlow")->LabelsDeflate();
    outsideprocessing.GetPointerAfterMerge2D("MpipiPairs")->LabelsDeflate("X");
    outsideprocessing.GetPointerAfterMerge2D("MpipiPairs")->LabelsDeflate("Y");
    outsideprocessing.GetPointerAfterMerge1D("MpipiMotherName")->LabelsDeflate();
    outsideprocessing.GetPointerAfterMerge2D("MpipiDeltaTvsEvent")->LabelsDeflate("Y");
    outsideprocessing.GetPointerAfterMerge1D("MatchingAbnormal")->LabelsDeflate();
    outsideprocessing.GetPointerAfterMerge1D("MatchingNormal")->LabelsDeflate();

    //setting up a tree & output file
    string path = string(argv[0]);
    string outfileName;
    if(outputFolder.find(".root")!=std::string::npos){
        outfileName = outputFolder;
    } else{
        outfileName = outputFolder+"SimOutput_"+path.substr(path.find_last_of("/\\")+1)+".root";
    }
    cout<<"Created output file "<<outfileName<<endl;

    //saving histograms
    TFile* outputFileHist = TFile::Open(outfileName.c_str(), "recreate");
    //first, custom ones (efficiency)
    TEfficiency K0SMCEfficiency(*outsideprocessing.GetPointerAfterMerge2D("K0SAfterMCpionCuts"), *outsideprocessing.GetPointerAfterMerge2D("K0STotal"));
    K0SMCEfficiency.SetNameTitle("K0SMCEfficiency", "K^{0}_{S} MC detection efficiency (judged by pions inside fiducial region)");
    TEfficiency K0STPCEfficiency(*outsideprocessing.GetPointerAfterMerge2D("K0SAfterTPCpionCuts"), *outsideprocessing.GetPointerAfterMerge2D("K0SAfterMCpionCuts"));
    K0STPCEfficiency.SetNameTitle("K0STPCEfficiency", "K^{0}_{S} TPC detection efficiency (judged by TPC reconstruction of pions already in fiducial cuts)");
    TEfficiency K0STotalEfficiency(*outsideprocessing.GetPointerAfterMerge2D("K0SAfterTPCpionCuts"), *outsideprocessing.GetPointerAfterMerge2D("K0STotal"));
    K0STotalEfficiency.SetNameTitle("K0STotalEfficiency", "K^{0}_{S} total detection efficiency (judged by MC pions inside fiducial region, properly reconstructed)");
    K0SMCEfficiency.Write();
    K0STPCEfficiency.Write();
    K0STotalEfficiency.Write();
    //then the rest
    outsideprocessing.SaveToFile(outputFileHist);
    outputFileHist->Close();

    return 0;
}

void PrintBigger(TParticle* input, std::string additional_stuff){
    Printf("TParticle: %-13s  p: %8f %8f %8f \tVertex: %8e %8e %8e \tProd. Vertex: %5d %5d \tDecay Vertex: %5d \tTOF tray:%5d \tTOF module:%5d%s",
        input->GetName(), input->Px(), input->Py(), input->Pz(), input->Vx(), input->Vy(), input->Vz(),
        input->GetFirstMother(), input->GetSecondMother(), input->GetFirstDaughter(),
        input->GetLastDaughter()/100, input->GetLastDaughter()%100, additional_stuff.c_str());
}

double distanceToBeamline(TVector3 point, double x0, double y0, double dxdz, double dydz){
    TVector3 beamlinePoint(x0, y0, 0.);
    TVector3 beamlineDirection(dxdz, dydz, 1.);
    //https://math.stackexchange.com/questions/2353288/point-to-line-distance-in-3d-using-cross-product
    return beamlineDirection.Cross(beamlinePoint-point).Mag()/beamlineDirection.Mag();
}

double DeltaTFix(double deltaT, TParticle* particle1, TParticle* particle2, StUPCEvent* upcEvent){
    std::map<std::pair<int, int>, double> TOFTimingMap;
    std::map<std::pair<int, int>, double> TOFTimingDifferenceMap;
    //tof maps filling
    for(int i = 0; i<upcEvent->getNumberOfHits(); i++){
        int currentTray = upcEvent->getHit(i)->getTray();
        int currentModule = upcEvent->getHit(i)->getModule();
        double currentLeadingEdge = upcEvent->getHit(i)->getLeadingEdgeTime();
        //insertion; if trying to insert a second hit in the same position (so emplace gives "false")
        //a later hit is inserted and a difference of those is recorded
        if(!TOFTimingMap.emplace(std::pair<int, int>(currentTray, currentModule), currentLeadingEdge).second){
            //checking if the tof hit is bigger than the one currently stored
            if(TOFTimingMap[std::pair<int, int>(currentTray, currentModule)]<currentLeadingEdge){
                TOFTimingDifferenceMap[std::pair<int, int>(currentTray, currentModule)] += currentLeadingEdge-TOFTimingMap[std::pair<int, int>(currentTray, currentModule)];
                TOFTimingMap[std::pair<int, int>(currentTray, currentModule)] = currentLeadingEdge;
            } else{
                TOFTimingDifferenceMap[std::pair<int, int>(currentTray, currentModule)] -= currentLeadingEdge-TOFTimingMap[std::pair<int, int>(currentTray, currentModule)];
            }
        } else{
            TOFTimingDifferenceMap.emplace(std::pair<int, int>(currentTray, currentModule), 0.);
        }
    }
    int tray1 = particle1->GetLastDaughter()/100;
    int module1 = particle1->GetLastDaughter()%100;
    int tray2 = particle2->GetLastDaughter()/100;
    int module2 = particle2->GetLastDaughter()%100;
    return deltaT+TOFTimingDifferenceMap[std::pair<int, int>(tray1, module1)]-TOFTimingDifferenceMap[std::pair<int, int>(tray2, module2)];
}

double DeltaTFixSimilarModules(double deltaT, TParticle* particle1, TParticle* particle2, StUPCEvent* upcEvent){
    std::map<std::pair<int, int>, double> TOFTimingMap;
    std::map<std::pair<int, int>, double> TOFTimingMapLowerTray;
    std::map<std::pair<int, int>, double> TOFTimingMapUpperTray;
    std::map<std::pair<int, int>, double> TOFTimingMapLowerModule;
    std::map<std::pair<int, int>, double> TOFTimingMapUpperModule;
    //tof main map filling
    for(int i = 0; i<upcEvent->getNumberOfHits(); i++){
        int currentTray = upcEvent->getHit(i)->getTray();
        int currentModule = upcEvent->getHit(i)->getModule();
        double currentLeadingEdge = upcEvent->getHit(i)->getLeadingEdgeTime();
        //insertion; if trying to insert a second hit in the same position (so emplace gives "false")
        //an earlier hit is kept
        if(!TOFTimingMap.emplace(std::pair<int, int>(currentTray, currentModule), currentLeadingEdge).second){
            //checking if the tof hit is bigger than the one currently stored
            if(TOFTimingMap[std::pair<int, int>(currentTray, currentModule)]>currentLeadingEdge){
                TOFTimingMap[std::pair<int, int>(currentTray, currentModule)] = currentLeadingEdge;
            }
        }
    }
    //tof neighbour maps filling
    for(int i = 0; i<upcEvent->getNumberOfHits(); i++){
        int currentTray = upcEvent->getHit(i)->getTray();
        int currentModule = upcEvent->getHit(i)->getModule();
        double currentLeadingEdge = upcEvent->getHit(i)->getLeadingEdgeTime();
        //going througn alleged neighbours and checking if the hit can be a neighbour to someone
        if(TOFTimingMap.count(std::pair<int, int>(currentTray+1, currentModule))==1){
            TOFTimingMapLowerTray.emplace(std::pair<int, int>(currentTray, currentModule), currentLeadingEdge);
        }
        if(TOFTimingMap.count(std::pair<int, int>(currentTray-1, currentModule))==1){
            TOFTimingMapUpperTray.emplace(std::pair<int, int>(currentTray, currentModule), currentLeadingEdge);
        }
        if(TOFTimingMap.count(std::pair<int, int>(currentTray, currentModule+1))==1){
            TOFTimingMapLowerModule.emplace(std::pair<int, int>(currentTray, currentModule), currentLeadingEdge);
        }
        if(TOFTimingMap.count(std::pair<int, int>(currentTray, currentModule-1))==1){
            TOFTimingMapUpperModule.emplace(std::pair<int, int>(currentTray, currentModule), currentLeadingEdge);
        }
    }
    double deltaTfix = 0;
    //checking particle1 times
    int TOFtray = particle1->GetLastDaughter()/100;
    int TOFmodule = particle1->GetLastDaughter()%100;
    std::pair<int, int> TOFhit(TOFtray, TOFmodule);
    if(TOFTimingMapLowerModule.count(TOFhit)==1){
        deltaTfix += TOFTimingMapLowerModule[TOFhit]-TOFTimingMap[TOFhit];
    } else if(TOFTimingMapUpperModule.count(TOFhit)==1){
        deltaTfix += TOFTimingMapUpperModule[TOFhit]-TOFTimingMap[TOFhit];
    } else if(TOFTimingMapLowerTray.count(TOFhit)==1){
        deltaTfix += TOFTimingMapLowerTray[TOFhit]-TOFTimingMap[TOFhit];
    } else if(TOFTimingMapUpperTray.count(TOFhit)==1){
        deltaTfix += TOFTimingMapUpperTray[TOFhit]-TOFTimingMap[TOFhit];
    }
    //checking particle2 times
    TOFtray = particle2->GetLastDaughter()/100;
    TOFmodule = particle2->GetLastDaughter()%100;
    TOFhit = std::pair<int, int>(TOFtray, TOFmodule);
    if(TOFTimingMapLowerModule.count(TOFhit)==1){
        deltaTfix -= TOFTimingMapLowerModule[TOFhit]-TOFTimingMap[TOFhit];
    } else if(TOFTimingMapUpperModule.count(TOFhit)==1){
        deltaTfix -= TOFTimingMapUpperModule[TOFhit]-TOFTimingMap[TOFhit];
    } else if(TOFTimingMapLowerTray.count(TOFhit)==1){
        deltaTfix -= TOFTimingMapLowerTray[TOFhit]-TOFTimingMap[TOFhit];
    } else if(TOFTimingMapUpperTray.count(TOFhit)==1){
        deltaTfix -= TOFTimingMapUpperTray[TOFhit]-TOFTimingMap[TOFhit];
    }

    return deltaT+deltaTfix;
}

