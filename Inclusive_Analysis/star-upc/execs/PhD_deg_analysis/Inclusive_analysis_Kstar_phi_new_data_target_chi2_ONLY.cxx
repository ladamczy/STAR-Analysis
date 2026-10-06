//cpp headers
#include <map>
#include <fstream>
#include <deque>

//ROOT headers
#include <ROOT/TThreadedObject.hxx>
#include <TTreeReader.h>
#include <ROOT/TTreeProcessorMT.hxx>
#include <TRandom3.h>
#include <TFileCollection.h>
#include <THashList.h>
#include <TFileInfo.h>
#include <TSystem.h>

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
#include "StEfficiencyCorrector3D.h"

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
const double particleMass[nParticles] = { 0.13957, 0.493677, 0.93827 }; // pion, kaon, proton in GeV /c^2 
enum BRANCH_ID{ EU, ED, WU, WD, nBranches };
enum RP_ID{ E1U, E1D, E2U, E2D, W1U, W1D, W2U, W2D, nRomanPots };
enum SUSPECTED_PARTICLES{ K0S, Lambda, Kstar, Phi };
string particleNicks[EXTENDED_PARTICLES::nParticlesExtended] = { "e", "pi", "K", "p" };

enum CHARGE{ positive = 0, negative = 1, nCharge = 2 };
TFile* getNewestFile(std::string folder, std::string filename);
StEfficiencyCorrector3D* getInitialisedEfficiencyCorrector(std::string folder, CHARGE charge, PARTICLES particle);

int main(int argc, char** argv){

    //argv:
    //1 - .root file or .list list
    //2 - output folder
    //3 - (optional) number of cores
    //4 - (optional) number of events to look back for (default 1 to speed things up, 1000 gives good results for mixed background)
    //5 - (optional) folder with efficiency files for Sneha's efficiency correction (named SPPion*.root/SPProton*.root/SPKaon*.root)

    int nthreads = 1;
    int total_events_in_queue = 1000;
    std::string efficiency_folder = "";
    switch(argc){
    case 6:
        nthreads = atoi(argv[3]);
        total_events_in_queue = atoi(argv[4]);
        efficiency_folder = argv[5];
        break;
    case 5:
        //if the last one is a path, then we have # of cores & path to efficiency corrections
        //if both are numbers, we have # of cores & number of previous events
        if(atoi(argv[4])==0){
            nthreads = atoi(argv[3]);
            efficiency_folder = argv[4];
        } else{
            nthreads = atoi(argv[3]);
            total_events_in_queue = atoi(argv[4]);
        }
        break;
    case 4:
        //if it is the path, then we have path to efficiency corrections
        //if it is the number, we have number of cores
        if(atoi(argv[3])==0){
            efficiency_folder = argv[3];
        } else{
            nthreads = atoi(argv[3]);
        }
        break;
    default:
        printf("Invalid number of arguments. Proper argument usage:\n");
        printf("argv[1] - .root file or .list list\n");
        printf("argv[2] - output folder\n");
        printf("argv[3] - (optional) number of cores (default 1)\n");
        printf("argv[4] - (optional) number of events to look back for (default 1000 - gives good results for mixed background)\n");
        printf("argv[5] - (optional) folder for Sneha's efficiency correction files(named SPPion*.root/SPProton*.root/SPKaon*.root)\n");
        return 1;
        break;
    }
    //summary
    printf("Program is running on %d threads\n", nthreads);
    printf("Previous events to combine for mixed background: %d\n", total_events_in_queue);
    printf("Efficiency corrections: %s\n", efficiency_folder.length()!=0 ? efficiency_folder.c_str() : "Not used");
    //because number of events includes the current one, we need to add 1
    total_events_in_queue++;

    //setting up efficiency correction
    StEfficiencyCorrector3D* totalEfficiencyCorrector[nCharge][nParticles];
    bool loadedCorrectorCorrectly = true;
    for(size_t i = 0; i<nCharge; i++){
        for(size_t j = 0; j<nParticles; j++){
            totalEfficiencyCorrector[i][j] = getInitialisedEfficiencyCorrector(efficiency_folder, static_cast<CHARGE>(i), static_cast<PARTICLES>(j));
            if(totalEfficiencyCorrector[i][j]==nullptr){
                loadedCorrectorCorrectly = false;
                break;
            }
        }
        if(!loadedCorrectorCorrectly)
            break;
    }
    printf("Corrections initialised %s\n", loadedCorrectorCorrectly ? "successfully!" : "unsuccessfully! Falling back to default weight of 1.");
    //a nice wrapper
    auto correction_coefficient = [&](int total_pair_charge, PARTICLES positive_id, PARTICLES negative_id, StUPCTrack* positive_track, StUPCTrack* negative_track, double Vz_positive, double Vz_negative){
        //fallback if something broke
        if(!loadedCorrectorCorrectly)
            return 1.;
        //if success, we gotta load all the variables
        double pt1 = positive_track->getPt();
        double pt2 = negative_track->getPt();
        double eta1 = positive_track->getEta();
        double eta2 = negative_track->getEta();

        //choosing correct signs for detected tracks
        CHARGE first_track_charge_enum, second_track_charge_enum;
        int first_track_charge, second_track_charge;
        if(total_pair_charge==0){
            first_track_charge_enum = positive;
            first_track_charge = 1;
            second_track_charge_enum = negative;
            second_track_charge = -1;
        } else if(total_pair_charge>0){
            first_track_charge_enum = positive;
            first_track_charge = 1;
            second_track_charge_enum = positive;
            second_track_charge = 1;
        } else if(total_pair_charge<0){
            first_track_charge_enum = negative;
            first_track_charge = -1;
            second_track_charge_enum = negative;
            second_track_charge = -1;
        } else{
            return 1.;
        }

        //calculating correction
        double efficiency = 1.;
        //positive
        switch(positive_id){
        case Pion:
            efficiency *= totalEfficiencyCorrector[first_track_charge_enum][positive_id]->getCombinedEfficiency(eta1, pt1, Vz_positive, first_track_charge, StEfficiencyCorrector3D::PION);
            break;
        case Kaon:
            efficiency *= totalEfficiencyCorrector[first_track_charge_enum][positive_id]->getCombinedEfficiency(eta1, pt1, Vz_positive, first_track_charge, StEfficiencyCorrector3D::KAON);
            break;
        case Proton:
            efficiency *= totalEfficiencyCorrector[first_track_charge_enum][positive_id]->getCombinedEfficiency(eta1, pt1, Vz_positive, first_track_charge, StEfficiencyCorrector3D::PROTON);
            break;
        default:
            return 1.;
            break;
        }
        //negative
        switch(negative_id){
        case Pion:
            efficiency *= totalEfficiencyCorrector[second_track_charge_enum][negative_id]->getCombinedEfficiency(eta2, pt2, Vz_negative, second_track_charge, StEfficiencyCorrector3D::PION);
            break;
        case Kaon:
            efficiency *= totalEfficiencyCorrector[second_track_charge_enum][negative_id]->getCombinedEfficiency(eta2, pt2, Vz_negative, second_track_charge, StEfficiencyCorrector3D::KAON);
            break;
        case Proton:
            efficiency *= totalEfficiencyCorrector[second_track_charge_enum][negative_id]->getCombinedEfficiency(eta2, pt2, Vz_negative, second_track_charge, StEfficiencyCorrector3D::PROTON);
            break;
        default:
            return 1.;
            break;
        }

        return 1./efficiency;
    };

    ROOT::EnableThreadSafety();
    //actually i'm not sure if it's needed here
    ROOT::EnableImplicitMT(nthreads); //turn on multicore processing

    //preparing input & output
    TChain* upcChain = new TChain("mUPCTree");
    if(ConnectInput(argc, argv, upcChain)){
        cout<<"All files connected"<<endl;
    }
    const string& outputFolder = argv[2];

    //adding chi2 histograms
    //file with sigma values:
    ifstream sigmaFile;
    string line;
    char sigmaName[10];
    double sigmaValue;
    map<string, double> sigmaMap;
    sigmaFile.open("STAR-Analysis/Inclusive_Analysis/star-upc/AnaOutput_Inclusive_analysis_Kstar_phi_new_data_TOF_tests_sigmaValues.txt");
    while(getline(sigmaFile, line)){
        sscanf(line.c_str(), "%s\t\t%lf", sigmaName, &sigmaValue);
        sigmaMap.insert({ string(sigmaName), sigmaValue });
    }
    sigmaFile.close();

    //pairs that are actively looked for
    std::vector<std::string> pairTab = { "Kpi", "piK", "ppi", "pip", "KK", "pipi", "pp" };

    //histograms
    ProcessingOutsideLoop outsideprocessing;
    outsideprocessing.AddHistogram(TH1D("pairInfoSignal", "", 1, 0, 1));
    outsideprocessing.AddHistogram(TH1D("pairInfoBackgroundSameSign", "", 1, 0, 1));
    outsideprocessing.AddHistogram(TH1D("pairInfoBackgroundTrackRotation", "", 1, 0, 1));
    outsideprocessing.AddHistogram(TH1D("pairInfoBackgroundRandomTrackRotation", "", 1, 0, 1));

    outsideprocessing.AddHistogram(TH1D("Flowchart", "", 1, 0, 1));

    //mass histograms (signal)
    outsideprocessing.AddHistogram(TH1D("MKpiChi2", ";m_{K^{+}#pi^{-}} [GeV/c^{2}];Number of pairs", 200, 0.5, 2.0));
    outsideprocessing.AddHistogram(TH1D("MpiKChi2", ";m_{K^{-}#pi^{+}} [GeV/c^{2}];Number of pairs", 200, 0.5, 2.0));
    outsideprocessing.AddHistogram(TH1D("MppiChi2", ";m_{p^{+}#pi^{-}} [GeV/c^{2}];Number of pairs", 500, 1.0, 2.5));
    outsideprocessing.AddHistogram(TH1D("MpipChi2", ";m_{p^{-}#pi^{+}} [GeV/c^{2}];Number of pairs", 500, 1.0, 2.5));
    outsideprocessing.AddHistogram(TH1D("MKKChi2", ";m_{K^{+}K^{-}} [GeV/c^{2}];Number of pairs", 50, 0.99, 1.05));
    outsideprocessing.AddHistogram(TH1D("MpipiChi2", ";m_{#pi^{+}#pi^{-}} [GeV/c^{2}];Number of pairs", 600, 0.2, 1.4));
    outsideprocessing.AddHistogram(TH1D("MppChi2", ";m_{p^{+}p^{-}} [GeV/c^{2}];Number of pairs", 500, 1.5, 3.5));
    //closer histograms (signal)
    outsideprocessing.AddHistogram(TH1D("MKKChi2Close", ";m_{K^{+}K^{-}} [GeV/c^{2}];Number of pairs", 50, 0.99, 1.05));
    outsideprocessing.AddHistogram(TH1D("MKpiChi2Close", ";m_{K^{+}#pi^{-}} [GeV/c^{2}];Number of pairs", 50, 0.7, 1.1));
    outsideprocessing.AddHistogram(TH1D("MpiKChi2Close", ";m_{K^{-}#pi^{+}} [GeV/c^{2}];Number of pairs", 50, 0.7, 1.1));
    //adding mass histograms grouped by category (signal)
    getCategoryHistograms(outsideprocessing, pairTab);

    //mass histograms (background, same sign)
    outsideprocessing.AddHistogram(TH1D("MKpiChi2BcgSameSign", ";m_{K^{+}#pi^{-}} [GeV/c^{2}];Number of pairs", 200, 0.5, 2.0));
    outsideprocessing.AddHistogram(TH1D("MpiKChi2BcgSameSign", ";m_{K^{-}#pi^{+}} [GeV/c^{2}];Number of pairs", 200, 0.5, 2.0));
    outsideprocessing.AddHistogram(TH1D("MppiChi2BcgSameSign", ";m_{p^{+}#pi^{-}} [GeV/c^{2}];Number of pairs", 500, 1.0, 2.5));
    outsideprocessing.AddHistogram(TH1D("MpipChi2BcgSameSign", ";m_{p^{-}#pi^{+}} [GeV/c^{2}];Number of pairs", 500, 1.0, 2.5));
    outsideprocessing.AddHistogram(TH1D("MKKChi2BcgSameSign", ";m_{K^{+}K^{-}} [GeV/c^{2}];Number of pairs", 50, 0.99, 1.05));
    outsideprocessing.AddHistogram(TH1D("MpipiChi2BcgSameSign", ";m_{#pi^{+}#pi^{-}} [GeV/c^{2}];Number of pairs", 600, 0.2, 1.4));
    outsideprocessing.AddHistogram(TH1D("MppChi2BcgSameSign", ";m_{p^{+}p^{-}} [GeV/c^{2}];Number of pairs", 500, 1.5, 3.5));
    //closer histograms (background, same sign)
    outsideprocessing.AddHistogram(TH1D("MKKChi2BcgSameSignClose", ";m_{K^{+}K^{-}} [GeV/c^{2}];Number of pairs", 50, 0.99, 1.05));
    outsideprocessing.AddHistogram(TH1D("MKpiChi2BcgSameSignClose", ";m_{K^{+}#pi^{-}} [GeV/c^{2}];Number of pairs", 50, 0.7, 1.1));
    outsideprocessing.AddHistogram(TH1D("MpiKChi2BcgSameSignClose", ";m_{K^{-}#pi^{+}} [GeV/c^{2}];Number of pairs", 50, 0.7, 1.1));
    //adding mass histograms grouped by category (background, same sign)
    getCategoryHistograms(outsideprocessing, pairTab, "BcgSameSign");

    //mass histograms (background, track rotation)
    outsideprocessing.AddHistogram(TH1D("MKpiChi2BcgTrackRotation", ";m_{K^{+}#pi^{-}} [GeV/c^{2}];Number of pairs", 200, 0.5, 2.0));
    outsideprocessing.AddHistogram(TH1D("MpiKChi2BcgTrackRotation", ";m_{K^{-}#pi^{+}} [GeV/c^{2}];Number of pairs", 200, 0.5, 2.0));
    outsideprocessing.AddHistogram(TH1D("MppiChi2BcgTrackRotation", ";m_{p^{+}#pi^{-}} [GeV/c^{2}];Number of pairs", 500, 1.0, 2.5));
    outsideprocessing.AddHistogram(TH1D("MpipChi2BcgTrackRotation", ";m_{p^{-}#pi^{+}} [GeV/c^{2}];Number of pairs", 500, 1.0, 2.5));
    outsideprocessing.AddHistogram(TH1D("MKKChi2BcgTrackRotation", ";m_{K^{+}K^{-}} [GeV/c^{2}];Number of pairs", 50, 0.99, 1.05));
    outsideprocessing.AddHistogram(TH1D("MpipiChi2BcgTrackRotation", ";m_{#pi^{+}#pi^{-}} [GeV/c^{2}];Number of pairs", 600, 0.2, 1.4));
    outsideprocessing.AddHistogram(TH1D("MppChi2BcgTrackRotation", ";m_{p^{+}p^{-}} [GeV/c^{2}];Number of pairs", 500, 1.5, 3.5));
    //closer histograms (background, track rotation)
    outsideprocessing.AddHistogram(TH1D("MKKChi2BcgTrackRotationClose", ";m_{K^{+}K^{-}} [GeV/c^{2}];Number of pairs", 50, 0.99, 1.05));
    outsideprocessing.AddHistogram(TH1D("MKpiChi2BcgTrackRotationClose", ";m_{K^{+}#pi^{-}} [GeV/c^{2}];Number of pairs", 50, 0.7, 1.1));
    outsideprocessing.AddHistogram(TH1D("MpiKChi2BcgTrackRotationClose", ";m_{K^{-}#pi^{+}} [GeV/c^{2}];Number of pairs", 50, 0.7, 1.1));
    //adding mass histograms grouped by category (background, track rotation)
    getCategoryHistograms(outsideprocessing, pairTab, "BcgTrackRotation");

    //mass histograms (background, random track rotation)
    outsideprocessing.AddHistogram(TH1D("MKpiChi2BcgRandomTrackRotation", ";m_{K^{+}#pi^{-}} [GeV/c^{2}];Number of pairs", 200, 0.5, 2.0));
    outsideprocessing.AddHistogram(TH1D("MpiKChi2BcgRandomTrackRotation", ";m_{K^{-}#pi^{+}} [GeV/c^{2}];Number of pairs", 200, 0.5, 2.0));
    outsideprocessing.AddHistogram(TH1D("MppiChi2BcgRandomTrackRotation", ";m_{p^{+}#pi^{-}} [GeV/c^{2}];Number of pairs", 500, 1.0, 2.5));
    outsideprocessing.AddHistogram(TH1D("MpipChi2BcgRandomTrackRotation", ";m_{p^{-}#pi^{+}} [GeV/c^{2}];Number of pairs", 500, 1.0, 2.5));
    outsideprocessing.AddHistogram(TH1D("MKKChi2BcgRandomTrackRotation", ";m_{K^{+}K^{-}} [GeV/c^{2}];Number of pairs", 50, 0.99, 1.05));
    outsideprocessing.AddHistogram(TH1D("MpipiChi2BcgRandomTrackRotation", ";m_{#pi^{+}#pi^{-}} [GeV/c^{2}];Number of pairs", 600, 0.2, 1.4));
    outsideprocessing.AddHistogram(TH1D("MppChi2BcgRandomTrackRotation", ";m_{p^{+}p^{-}} [GeV/c^{2}];Number of pairs", 500, 1.5, 3.5));
    //closer histograms (background, random track rotation)
    outsideprocessing.AddHistogram(TH1D("MKKChi2BcgRandomTrackRotationClose", ";m_{K^{+}K^{-}} [GeV/c^{2}];Number of pairs", 50, 0.99, 1.05));
    outsideprocessing.AddHistogram(TH1D("MKpiChi2BcgRandomTrackRotationClose", ";m_{K^{+}#pi^{-}} [GeV/c^{2}];Number of pairs", 50, 0.7, 1.1));
    outsideprocessing.AddHistogram(TH1D("MpiKChi2BcgRandomTrackRotationClose", ";m_{K^{-}#pi^{+}} [GeV/c^{2}];Number of pairs", 50, 0.7, 1.1));
    //adding mass histograms grouped by category (background, random track rotation)
    getCategoryHistograms(outsideprocessing, pairTab, "BcgRandomTrackRotation");

    //mass histograms (background, mixed events)
    outsideprocessing.AddHistogram(TH1D("MKpiChi2BcgMixedEvent", ";m_{K^{+}#pi^{-}} [GeV/c^{2}];Number of pairs", 200, 0.5, 2.0));
    outsideprocessing.AddHistogram(TH1D("MpiKChi2BcgMixedEvent", ";m_{K^{-}#pi^{+}} [GeV/c^{2}];Number of pairs", 200, 0.5, 2.0));
    outsideprocessing.AddHistogram(TH1D("MppiChi2BcgMixedEvent", ";m_{p^{+}#pi^{-}} [GeV/c^{2}];Number of pairs", 500, 1.0, 2.5));
    outsideprocessing.AddHistogram(TH1D("MpipChi2BcgMixedEvent", ";m_{p^{-}#pi^{+}} [GeV/c^{2}];Number of pairs", 500, 1.0, 2.5));
    outsideprocessing.AddHistogram(TH1D("MKKChi2BcgMixedEvent", ";m_{K^{+}K^{-}} [GeV/c^{2}];Number of pairs", 50, 0.99, 1.05));
    outsideprocessing.AddHistogram(TH1D("MpipiChi2BcgMixedEvent", ";m_{#pi^{+}#pi^{-}} [GeV/c^{2}];Number of pairs", 600, 0.2, 1.4));
    outsideprocessing.AddHistogram(TH1D("MppChi2BcgMixedEvent", ";m_{p^{+}p^{-}} [GeV/c^{2}];Number of pairs", 500, 1.5, 3.5));
    //closer histograms (background, mixed events)
    outsideprocessing.AddHistogram(TH1D("MKKChi2BcgMixedEventClose", ";m_{K^{+}K^{-}} [GeV/c^{2}];Number of pairs", 50, 0.99, 1.05));
    outsideprocessing.AddHistogram(TH1D("MKpiChi2BcgMixedEventClose", ";m_{K^{+}#pi^{-}} [GeV/c^{2}];Number of pairs", 50, 0.7, 1.1));
    outsideprocessing.AddHistogram(TH1D("MpiKChi2BcgMixedEventClose", ";m_{K^{-}#pi^{+}} [GeV/c^{2}];Number of pairs", 50, 0.7, 1.1));
    //adding mass histograms grouped by category (background, mixed events)
    getCategoryHistograms(outsideprocessing, pairTab, "BcgMixedEvent");

    //mass histograms (background, mixed events same sign)
    outsideprocessing.AddHistogram(TH1D("MKpiChi2BcgMixedEventSameSign", ";m_{K^{+}#pi^{-}} [GeV/c^{2}];Number of pairs", 200, 0.5, 2.0));
    outsideprocessing.AddHistogram(TH1D("MpiKChi2BcgMixedEventSameSign", ";m_{K^{-}#pi^{+}} [GeV/c^{2}];Number of pairs", 200, 0.5, 2.0));
    outsideprocessing.AddHistogram(TH1D("MppiChi2BcgMixedEventSameSign", ";m_{p^{+}#pi^{-}} [GeV/c^{2}];Number of pairs", 500, 1.0, 2.5));
    outsideprocessing.AddHistogram(TH1D("MpipChi2BcgMixedEventSameSign", ";m_{p^{-}#pi^{+}} [GeV/c^{2}];Number of pairs", 500, 1.0, 2.5));
    outsideprocessing.AddHistogram(TH1D("MKKChi2BcgMixedEventSameSign", ";m_{K^{+}K^{-}} [GeV/c^{2}];Number of pairs", 50, 0.99, 1.05));
    outsideprocessing.AddHistogram(TH1D("MpipiChi2BcgMixedEventSameSign", ";m_{#pi^{+}#pi^{-}} [GeV/c^{2}];Number of pairs", 600, 0.2, 1.4));
    outsideprocessing.AddHistogram(TH1D("MppChi2BcgMixedEventSameSign", ";m_{p^{+}p^{-}} [GeV/c^{2}];Number of pairs", 500, 1.5, 3.5));
    //closer histograms (background, mixed events)
    outsideprocessing.AddHistogram(TH1D("MKKChi2BcgMixedEventSameSignClose", ";m_{K^{+}K^{-}} [GeV/c^{2}];Number of pairs", 50, 0.99, 1.05));
    outsideprocessing.AddHistogram(TH1D("MKpiChi2BcgMixedEventSameSignClose", ";m_{K^{+}#pi^{-}} [GeV/c^{2}];Number of pairs", 50, 0.7, 1.1));
    outsideprocessing.AddHistogram(TH1D("MpiKChi2BcgMixedEventSameSignClose", ";m_{K^{-}#pi^{+}} [GeV/c^{2}];Number of pairs", 50, 0.7, 1.1));
    //adding mass histograms grouped by category (background, mixed events same sign)
    getCategoryHistograms(outsideprocessing, pairTab, "BcgMixedEventSameSign");


    //other mass histograms
    outsideprocessing.AddHistogram(TH1D("MKKSuspiciousPeakTestedAsPionPair", ";m_{#pi^{+}#pi^{-}} [GeV/c^{2}];Number of pairs", 400, 0.25, 0.65));
    outsideprocessing.AddHistogram(TH1D("MKKSuspiciousPeakTestedAsPionPairNeighbourhood", ";m_{#pi^{+}#pi^{-}} [GeV/c^{2}];Number of pairs", 400, 0.25, 0.65));
    outsideprocessing.AddHistogram(TH1D("MKKSuspiciousPeakTestedWithStrictChi2LessThan3", ";m_{K^{+}K^{-}} [GeV/c^{2}];Number of pairs", 500, 0.9, 2.4));
    outsideprocessing.AddHistogram(TH1D("MKKSuspiciousPeakTestedWithStrictChi2LessThan1", ";m_{K^{+}K^{-}} [GeV/c^{2}];Number of pairs", 500, 0.9, 2.4));

    //other *other* histograms
    outsideprocessing.AddHistogram(TH2D("dEdxTPCEnergyLoss", "dE/dx TPC Energy loss distribution;pq [GeV/c];dE/dx [keV/cm]", 300, -3, 3, 100, 0, 40));
    outsideprocessing.AddHistogram(TH1D("PrimaryVertexPositionX", ";x [cm];Number of events", 160, -0.8, 0.8));
    outsideprocessing.AddHistogram(TH1D("PrimaryVertexPositionY", ";y [cm];Number of events", 160, -0.8, 0.8));
    outsideprocessing.AddHistogram(TH1D("PrimaryVertexPositionZ", ";z [cm];Number of events", 400, -200., 200));
    outsideprocessing.AddHistogram(TH1D("PrimaryVertexPositionTotal", ";R [cm];Number of events", 400, 0., 400));
    outsideprocessing.AddHistogram(TH1D("MxTotal", ";M_{X} [GeV/c^{2}];Number of events", 400, 0., 20.));
    outsideprocessing.AddHistogram(TH1D("XiWTotal", ";#xi_{W};Number of events", 400, -0.05, 0.35));
    outsideprocessing.AddHistogram(TH1D("XiETotal", ";#xi_{E};Number of events", 400, -0.05, 0.35));
    outsideprocessing.AddHistogram(TH2D("XiBothTotal", ";#xi_{W};#xi_{E}", 400, -0.05, 0.35, 400, -0.05, 0.35));
    //histograms for eventual cutting of mixed events
    outsideprocessing.AddHistogram(TH1D("PrimaryVertexPositionDifferenceX", ";#Delta x [cm];Number of event pairs", 160, -0.8, 0.8));
    outsideprocessing.AddHistogram(TH1D("PrimaryVertexPositionDifferenceY", ";#Delta y [cm];Number of event pairs", 160, -0.8, 0.8));
    outsideprocessing.AddHistogram(TH1D("PrimaryVertexPositionDifferenceZ", ";#Delta z [cm];Number of event pairs", 400, -200., 200));
    outsideprocessing.AddHistogram(TH1D("PrimaryVertexPositionDifferenceTotal", ";#Delta R [cm];Number of event pairs", 400, 0., 400));
    outsideprocessing.AddHistogram(TH1D("MxDifference", ";#Delta M_{X} [GeV/c^{2}];Number of event pairs", 400, -10., 10));
    outsideprocessing.AddHistogram(TH1D("MxDifferenceClose", ";#Delta M_{X} [GeV/c^{2}];Number of event pairs", 400, -0.2, 0.2));
    outsideprocessing.AddHistogram(TH1D("XiWDifference", ";#Delta #xi_{W};Number of event pairs", 400, -0.1, 0.1));
    outsideprocessing.AddHistogram(TH1D("XiEDifference", ";#Delta #xi_{E};Number of event pairs", 400, -0.1, 0.1));
    outsideprocessing.AddHistogram(TH2D("XiBothDifference", ";#Delta #xi_{W};#Delta #xi_{E}", 400, -0.1, 0.1, 400, -0.1, 0.1));

    for(size_t i = 0; i<outsideprocessing.GetNumberOfHistograms(); i++){
        if(&outsideprocessing.GetPointer1D(i)!=nullptr){
            printf("%s\n", outsideprocessing.GetPointer1D(i)->GetName());
        } else if(&outsideprocessing.GetPointer2D(i)!=nullptr){
            printf("%s\n", outsideprocessing.GetPointer2D(i)->GetName());
        } else if(&outsideprocessing.GetPointer3D(i)!=nullptr){
            printf("%s\n", outsideprocessing.GetPointer3D(i)->GetName());
        } else{
            printf("What.\n");
        }
    }


    //processing
    //defining TreeProcessor
    ROOT::TTreeProcessorMT TreeProc(*upcChain, nthreads);

    //other things
    int eventsProcessed = 0;

    //defining processing function
    auto myFunction = [&](TTreeReader& myReader){
        //getting values from TChain, in-loop histogram initialization
        TTreeReaderValue<StUPCEvent> StUPCEventInstance(myReader, "mUPCEvent");
        TTreeReaderValue<StRPEvent> StRPEventInstance(myReader, "mRPEvent");
        ProcessingInsideLoop insideprocessing;
        StUPCEvent* tempUPCpointer;
        StRPEvent* tempRPpointer;
        insideprocessing.GetLocalHistograms(&outsideprocessing);

        //filling pairInfo histograms in correct order of bins
        for(auto&& histlastname:{ "Signal", "BackgroundSameSign", "BackgroundTrackRotation", "BackgroundRandomTrackRotation" }){
            std::string histfirstname = "pairInfo";
            insideprocessing.Fill((histfirstname+histlastname).c_str(), "OK", 0.0);
            insideprocessing.Fill((histfirstname+histlastname).c_str(), "TOF wrong", 0.0);
            insideprocessing.Fill((histfirstname+histlastname).c_str(), "dEdx wrong", 0.0);
            insideprocessing.Fill((histfirstname+histlastname).c_str(), "Both wrong", 0.0);
        }

        //helpful variables
        std::vector<StUPCTrack*> vector_Track_positive;
        std::vector<StUPCTrack*> vector_Track_negative;
        StUPCTrack* tempTrack;
        TLorentzVector positive_track;
        TLorentzVector negative_track;
        TLorentzVector positive_track2;
        TLorentzVector negative_track2;
        double mass, chi2pipi, chi2Kpi, eta, pT, phi, correction;
        map<string, double> chi2Map;
        bool isdEdxOk, isTOFOk;
        string tempPairName;
        TRandom3 random_generator;

        //queues for vectors of tracks from previous events assigned to particular pairs
        std::deque<std::vector<StUPCTrack*>> queue_of_previous_vector_Tracks_Kpi_positive;
        std::deque<std::vector<StUPCTrack*>> queue_of_previous_vector_Tracks_piK_positive;
        std::deque<std::vector<StUPCTrack*>> queue_of_previous_vector_Tracks_ppi_positive;
        std::deque<std::vector<StUPCTrack*>> queue_of_previous_vector_Tracks_pip_positive;
        std::deque<std::vector<StUPCTrack*>> queue_of_previous_vector_Tracks_KK_positive;
        std::deque<std::vector<StUPCTrack*>> queue_of_previous_vector_Tracks_pipi_positive;
        std::deque<std::vector<StUPCTrack*>> queue_of_previous_vector_Tracks_pp_positive;
        std::deque<std::vector<StUPCTrack*>> queue_of_previous_vector_Tracks_Kpi_negative;
        std::deque<std::vector<StUPCTrack*>> queue_of_previous_vector_Tracks_piK_negative;
        std::deque<std::vector<StUPCTrack*>> queue_of_previous_vector_Tracks_ppi_negative;
        std::deque<std::vector<StUPCTrack*>> queue_of_previous_vector_Tracks_pip_negative;
        std::deque<std::vector<StUPCTrack*>> queue_of_previous_vector_Tracks_KK_negative;
        std::deque<std::vector<StUPCTrack*>> queue_of_previous_vector_Tracks_pipi_negative;
        std::deque<std::vector<StUPCTrack*>> queue_of_previous_vector_Tracks_pp_negative;
        //positions of primary vertex from previous events
        std::deque<TVector3> queue_of_previous_PV_positions;
        //total central mass
        std::deque<double> queue_of_previous_Mx;
        //xi of protons
        std::deque<double> queue_of_previous_Xi_W;
        std::deque<double> queue_of_previous_Xi_E;

        //actual loop
        while(myReader.Next()){
            //in a TTree, it *would* be constant, in TChain however not necessarily
            tempUPCpointer = StUPCEventInstance.Get();
            tempRPpointer = StRPEventInstance.Get();

            //cleaning the loop
            vector_Track_positive.clear();
            vector_Track_negative.clear();

            //cause I want to see what's going on
            if(eventsProcessed%10000==0){
                cout<<"Processed "<<eventsProcessed<<" events"<<endl;
            }
            eventsProcessed++;

            //additional cuts that normally are used in only-RP-cuts examples
            //cause the data i used isn't properly filtered
            //making sure i don't get garbage from zerobias trigger
            insideprocessing.Fill("Flowchart", "All", 1.0);
            if(tempUPCpointer->isTrigger(570704)){
                continue;
            }
            insideprocessing.Fill("Flowchart", "Not zerobias", 1.0);
            //at least 2 good tracks
            int nOfGoodTracks = 0;
            for(int i = 0; i<tempUPCpointer->getNumberOfTracks(); i++){
                StUPCTrack* tmptrk = tempUPCpointer->getTrack(i);
                if(tmptrk->getFlag(StUPCTrack::kTof)&&tmptrk->getFlag(StUPCTrack::kPrimary)&&fabs(tmptrk->getEta())<0.9&&tmptrk->getPt()>0.2){
                    nOfGoodTracks++;
                }
            }
            if(nOfGoodTracks<2){
                continue;
            }
            insideprocessing.Fill("Flowchart", "#geq 2 ToF tracks", 1.0);
            //exactly one vertex
            if(tempUPCpointer->getNumberOfVertices()!=1){
                continue;
            }
            insideprocessing.Fill("Flowchart", "One vertex", 1.0);

            //cuts & histogram filling
            //selecting tracks matching criteria:
            //TOF
            //pt & eta 
            //Nhits

            //reusing good track counter
            nOfGoodTracks = 0;
            for(int i = 0; i<tempUPCpointer->getNumberOfTracks(); i++){
                tempTrack = tempUPCpointer->getTrack(i);
                if(!tempTrack->getFlag(StUPCTrack::kTof)){
                    continue;
                }
                if(!tempTrack->getFlag(StUPCTrack::kPrimary)){
                    continue;
                }
                if(tempTrack->getPt()<=0.2 or fabs(tempTrack->getEta())>=0.9){
                    continue;
                }
                if(tempTrack->getNhits()<=20){
                    continue;
                }
                //grouping properly reconstructed particles by their charge
                if(tempTrack->getCharge()>0){
                    vector_Track_positive.push_back(tempTrack);
                } else{
                    vector_Track_negative.push_back(tempTrack);
                }
                nOfGoodTracks++;
                //filling dE/dx (with the same cuts as in MC_signal_id etc.)
                if(tempTrack->getNhitsDEdx()>=15&&tempTrack->getTofPathLength()>0&&tempTrack->getTofTime()>0){
                    TVector3 temp_momentum;
                    tempTrack->getMomentum(temp_momentum);
                    insideprocessing.Fill("dEdxTPCEnergyLoss", temp_momentum.Mag()*tempTrack->getCharge(), tempTrack->getDEdxSignal()*1e6);
                }
            }

            //filling a chi2 map with keys for all the possibilities
            for(auto const& imap:sigmaMap){
                chi2Map.insert({ imap.first, 0. });
            }

            //filling other one-in-an-event histograms
            TVector3 PV_position(tempUPCpointer->getVertex(0)->getPosX(), tempUPCpointer->getVertex(0)->getPosY(), tempUPCpointer->getVertex(0)->getPosZ());
            insideprocessing.Fill("PrimaryVertexPositionX", PV_position.X());
            insideprocessing.Fill("PrimaryVertexPositionY", PV_position.Y());
            insideprocessing.Fill("PrimaryVertexPositionZ", PV_position.Z());
            insideprocessing.Fill("PrimaryVertexPositionTotal", PV_position.Mag());
            double Mx = tempRPpointer->getTrack(0)->xi(beamMomentum)*tempRPpointer->getTrack(1)->xi(beamMomentum)*510.;
            insideprocessing.Fill("MxTotal", Mx);
            double Xi_W, Xi_E;
            if(tempRPpointer->getTrack(0)->branch()==WU||tempRPpointer->getTrack(0)->branch()==WD){
                Xi_W = tempRPpointer->getTrack(0)->xi(beamMomentum);
                Xi_E = tempRPpointer->getTrack(1)->xi(beamMomentum);
            } else{
                Xi_W = tempRPpointer->getTrack(1)->xi(beamMomentum);
                Xi_E = tempRPpointer->getTrack(0)->xi(beamMomentum);
            }
            insideprocessing.Fill("XiWTotal", Xi_W);
            insideprocessing.Fill("XiETotal", Xi_E);
            insideprocessing.Fill("XiBothTotal", Xi_W, Xi_E);

            //########## SIGNAL EXTRACTION ###############

            //loop through identified particles (signal)
            for(long unsigned int i = 0; i<vector_Track_positive.size(); i++){
                for(long unsigned int j = 0; j<vector_Track_negative.size(); j++){
                    isdEdxOk = (vector_Track_positive[i]->getNhitsDEdx()>=15)&&(vector_Track_negative[j]->getNhitsDEdx()>=15);
                    isTOFOk = (vector_Track_positive[i]->getTofPathLength()>0)&&(vector_Track_positive[i]->getTofTime()>0)&&(vector_Track_negative[j]->getTofPathLength()>0)&&(vector_Track_negative[j]->getTofTime()>0);
                    if(isdEdxOk&&isTOFOk){
                        insideprocessing.Fill("pairInfoSignal", "OK", 1.0);
                    } else if(isdEdxOk&&!isTOFOk){
                        insideprocessing.Fill("pairInfoSignal", "TOF wrong", 1.0);
                        continue;
                    } else if(!isdEdxOk&&isTOFOk){
                        insideprocessing.Fill("pairInfoSignal", "dEdx wrong", 1.0);
                        continue;
                    } else if(!isdEdxOk&&!isTOFOk){
                        insideprocessing.Fill("pairInfoSignal", "Both wrong", 1.0);
                        continue;
                    }

                    //chi2
                    for(size_t pos = 0; pos<nParticlesExtended; pos++){
                        for(size_t neg = 0; neg<nParticlesExtended; neg++){
                            tempPairName = particleNicks[pos]+"_"+particleNicks[neg];
                            chi2Map[tempPairName] = getChi2(vector_Track_positive[i], vector_Track_negative[j], pos, neg, sigmaMap[tempPairName]);
                        }
                    }

                    //mass tests on different pairs
                    if(chi2Map["K_pi"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Kaon]);
                        vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Pion]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(0, Kaon, Pion, vector_Track_positive[i], vector_Track_negative[j], PV_position.Z(), PV_position.Z());
                        insideprocessing.Fill("MKpiChi2", mass, correction);
                        insideprocessing.Fill("MKpiChi2Close", mass, correction);
                        insideprocessing.Fill("MKpiChi2eta", mass, eta, correction);
                        insideprocessing.Fill("MKpiChi2pT", mass, pT, correction);
                    }
                    if(chi2Map["pi_K"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                        vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Kaon]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(0, Pion, Kaon, vector_Track_positive[i], vector_Track_negative[j], PV_position.Z(), PV_position.Z());
                        insideprocessing.Fill("MpiKChi2", mass, correction);
                        insideprocessing.Fill("MpiKChi2Close", mass, correction);
                        insideprocessing.Fill("MpiKChi2eta", mass, eta, correction);
                        insideprocessing.Fill("MpiKChi2pT", mass, pT, correction);
                    }
                    if(chi2Map["p_pi"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Proton]);
                        vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Pion]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(0, Proton, Pion, vector_Track_positive[i], vector_Track_negative[j], PV_position.Z(), PV_position.Z());
                        insideprocessing.Fill("MppiChi2", mass, correction);
                        insideprocessing.Fill("MppiChi2eta", mass, eta, correction);
                        insideprocessing.Fill("MppiChi2pT", mass, pT, correction);
                    }
                    if(chi2Map["pi_p"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                        vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Proton]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(0, Pion, Proton, vector_Track_positive[i], vector_Track_negative[j], PV_position.Z(), PV_position.Z());
                        insideprocessing.Fill("MpipChi2", mass, correction);
                        insideprocessing.Fill("MpipChi2eta", mass, eta, correction);
                        insideprocessing.Fill("MpipChi2pT", mass, pT, correction);
                    }
                    if(chi2Map["K_K"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Kaon]);
                        vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Kaon]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(0, Kaon, Kaon, vector_Track_positive[i], vector_Track_negative[j], PV_position.Z(), PV_position.Z());
                        insideprocessing.Fill("MKKChi2", mass, correction);
                        insideprocessing.Fill("MKKChi2Close", mass, correction);
                        insideprocessing.Fill("MKKChi2eta", mass, eta, correction);
                        insideprocessing.Fill("MKKChi2pT", mass, pT, correction);
                        //test of suspicious peak and its neighbourhood
                        if(chi2Map["K_K"]<3){
                            insideprocessing.Fill("MKKSuspiciousPeakTestedWithStrictChi2LessThan3", mass, correction);
                        }
                        if(chi2Map["K_K"]<1){
                            insideprocessing.Fill("MKKSuspiciousPeakTestedWithStrictChi2LessThan1", mass, correction);
                        }
                        correction = correction_coefficient(0, Pion, Pion, vector_Track_positive[i], vector_Track_negative[j], PV_position.Z(), PV_position.Z());
                        if(1.06<mass&&mass<1.08){
                            vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                            vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Pion]);
                            insideprocessing.Fill("MKKSuspiciousPeakTestedAsPionPair", (positive_track+negative_track).M(), correction);
                        } else if(1.05<mass&&mass<1.09){
                            vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                            vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Pion]);
                            insideprocessing.Fill("MKKSuspiciousPeakTestedAsPionPairNeighbourhood", (positive_track+negative_track).M(), correction);
                        }
                    }
                    if(chi2Map["pi_pi"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                        vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Pion]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(0, Pion, Pion, vector_Track_positive[i], vector_Track_negative[j], PV_position.Z(), PV_position.Z());
                        insideprocessing.Fill("MpipiChi2", mass);
                        insideprocessing.Fill("MpipiChi2eta", mass, eta);
                        insideprocessing.Fill("MpipiChi2pT", mass, pT);
                    }
                    if(chi2Map["p_p"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Proton]);
                        vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Proton]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(0, Proton, Proton, vector_Track_positive[i], vector_Track_negative[j], PV_position.Z(), PV_position.Z());
                        insideprocessing.Fill("MppChi2", mass);
                        insideprocessing.Fill("MppChi2eta", mass, eta);
                        insideprocessing.Fill("MppChi2pT", mass, pT);
                    }
                }
            }

            //########## BACKGROUND EXTRACTION (SAME-SIGN) ###############

            //loop through identified particles (background, same-sign positive)
            for(long int i = 0; i+1<vector_Track_positive.size(); i++){
                for(long int j = i+1; j<vector_Track_positive.size(); j++){
                    isdEdxOk = (vector_Track_positive[i]->getNhitsDEdx()>=15)&&(vector_Track_positive[j]->getNhitsDEdx()>=15);
                    isTOFOk = (vector_Track_positive[i]->getTofPathLength()>0)&&(vector_Track_positive[i]->getTofTime()>0)&&(vector_Track_positive[j]->getTofPathLength()>0)&&(vector_Track_positive[j]->getTofTime()>0);
                    if(isdEdxOk&&isTOFOk){
                        insideprocessing.Fill("pairInfoBackgroundSameSign", "OK", 1.0);
                    } else if(isdEdxOk&&!isTOFOk){
                        insideprocessing.Fill("pairInfoBackgroundSameSign", "TOF wrong", 1.0);
                        continue;
                    } else if(!isdEdxOk&&isTOFOk){
                        insideprocessing.Fill("pairInfoBackgroundSameSign", "dEdx wrong", 1.0);
                        continue;
                    } else if(!isdEdxOk&&!isTOFOk){
                        insideprocessing.Fill("pairInfoBackgroundSameSign", "Both wrong", 1.0);
                        continue;
                    }

                    //chi2
                    for(size_t pos1 = 0; pos1<nParticlesExtended; pos1++){
                        for(size_t pos2 = 0; pos2<nParticlesExtended; pos2++){
                            tempPairName = particleNicks[pos1]+"_"+particleNicks[pos2];
                            chi2Map[tempPairName] = getChi2(vector_Track_positive[i], vector_Track_positive[j], pos1, pos2, sigmaMap[tempPairName]);
                        }
                    }

                    //mass tests on different pairs
                    if(chi2Map["K_pi"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Kaon]);
                        vector_Track_positive[j]->getLorentzVector(positive_track2, particleMass[Pion]);
                        mass = (positive_track+positive_track2).M();
                        eta = (positive_track+positive_track2).Eta();
                        pT = (positive_track+positive_track2).Pt();
                        correction = correction_coefficient(2, Kaon, Pion, vector_Track_positive[i], vector_Track_positive[j], PV_position.Z(), PV_position.Z());
                        insideprocessing.Fill("MKpiChi2BcgSameSign", mass, correction);
                        insideprocessing.Fill("MKpiChi2BcgSameSignClose", mass, correction);
                        insideprocessing.Fill("MKpiChi2BcgSameSigneta", mass, eta, correction);
                        insideprocessing.Fill("MKpiChi2BcgSameSignpT", mass, pT, correction);
                    }
                    if(chi2Map["pi_K"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                        vector_Track_positive[j]->getLorentzVector(positive_track2, particleMass[Kaon]);
                        mass = (positive_track+positive_track2).M();
                        eta = (positive_track+positive_track2).Eta();
                        pT = (positive_track+positive_track2).Pt();
                        correction = correction_coefficient(2, Pion, Kaon, vector_Track_positive[i], vector_Track_positive[j], PV_position.Z(), PV_position.Z());
                        insideprocessing.Fill("MpiKChi2BcgSameSign", mass, correction);
                        insideprocessing.Fill("MpiKChi2BcgSameSignClose", mass, correction);
                        insideprocessing.Fill("MpiKChi2BcgSameSigneta", mass, eta, correction);
                        insideprocessing.Fill("MpiKChi2BcgSameSignpT", mass, pT, correction);
                    }
                    if(chi2Map["p_pi"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Proton]);
                        vector_Track_positive[j]->getLorentzVector(positive_track2, particleMass[Pion]);
                        mass = (positive_track+positive_track2).M();
                        eta = (positive_track+positive_track2).Eta();
                        pT = (positive_track+positive_track2).Pt();
                        correction = correction_coefficient(2, Proton, Pion, vector_Track_positive[i], vector_Track_positive[j], PV_position.Z(), PV_position.Z());
                        insideprocessing.Fill("MppiChi2BcgSameSign", mass, correction);
                        insideprocessing.Fill("MppiChi2BcgSameSigneta", mass, eta, correction);
                        insideprocessing.Fill("MppiChi2BcgSameSignpT", mass, pT, correction);
                    }
                    if(chi2Map["pi_p"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                        vector_Track_positive[j]->getLorentzVector(positive_track2, particleMass[Proton]);
                        mass = (positive_track+positive_track2).M();
                        eta = (positive_track+positive_track2).Eta();
                        pT = (positive_track+positive_track2).Pt();
                        correction = correction_coefficient(2, Pion, Proton, vector_Track_positive[i], vector_Track_positive[j], PV_position.Z(), PV_position.Z());
                        insideprocessing.Fill("MpipChi2BcgSameSign", mass, correction);
                        insideprocessing.Fill("MpipChi2BcgSameSigneta", mass, eta, correction);
                        insideprocessing.Fill("MpipChi2BcgSameSignpT", mass, pT, correction);
                    }
                    if(chi2Map["K_K"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Kaon]);
                        vector_Track_positive[j]->getLorentzVector(positive_track2, particleMass[Kaon]);
                        mass = (positive_track+positive_track2).M();
                        eta = (positive_track+positive_track2).Eta();
                        pT = (positive_track+positive_track2).Pt();
                        correction = correction_coefficient(2, Kaon, Kaon, vector_Track_positive[i], vector_Track_positive[j], PV_position.Z(), PV_position.Z());
                        insideprocessing.Fill("MKKChi2BcgSameSign", mass, correction);
                        insideprocessing.Fill("MKKChi2BcgSameSignClose", mass, correction);
                        insideprocessing.Fill("MKKChi2BcgSameSigneta", mass, eta, correction);
                        insideprocessing.Fill("MKKChi2BcgSameSignpT", mass, pT, correction);
                    }
                    if(chi2Map["pi_pi"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                        vector_Track_positive[j]->getLorentzVector(positive_track2, particleMass[Pion]);
                        mass = (positive_track+positive_track2).M();
                        eta = (positive_track+positive_track2).Eta();
                        pT = (positive_track+positive_track2).Pt();
                        correction = correction_coefficient(2, Pion, Pion, vector_Track_positive[i], vector_Track_positive[j], PV_position.Z(), PV_position.Z());
                        insideprocessing.Fill("MpipiChi2BcgSameSign", mass, correction);
                        insideprocessing.Fill("MpipiChi2BcgSameSigneta", mass, eta, correction);
                        insideprocessing.Fill("MpipiChi2BcgSameSignpT", mass, pT, correction);
                    }
                    if(chi2Map["p_p"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Proton]);
                        vector_Track_positive[j]->getLorentzVector(positive_track2, particleMass[Proton]);
                        mass = (positive_track+positive_track2).M();
                        eta = (positive_track+positive_track2).Eta();
                        pT = (positive_track+positive_track2).Pt();
                        correction = correction_coefficient(2, Proton, Proton, vector_Track_positive[i], vector_Track_positive[j], PV_position.Z(), PV_position.Z());
                        insideprocessing.Fill("MppChi2BcgSameSign", mass, correction);
                        insideprocessing.Fill("MppChi2BcgSameSigneta", mass, eta, correction);
                        insideprocessing.Fill("MppChi2BcgSameSignpT", mass, pT, correction);
                    }
                }
            }
            //loop through identified particles (background, same-sign negative)
            for(long int i = 0; i+1<vector_Track_negative.size(); i++){
                for(long int j = i+1; j<vector_Track_negative.size(); j++){
                    isdEdxOk = (vector_Track_negative[i]->getNhitsDEdx()>=15)&&(vector_Track_negative[j]->getNhitsDEdx()>=15);
                    isTOFOk = (vector_Track_negative[i]->getTofPathLength()>0)&&(vector_Track_negative[i]->getTofTime()>0)&&(vector_Track_negative[j]->getTofPathLength()>0)&&(vector_Track_negative[j]->getTofTime()>0);
                    if(isdEdxOk&&isTOFOk){
                        insideprocessing.Fill("pairInfoBackgroundSameSign", "OK", 1.0);
                    } else if(isdEdxOk&&!isTOFOk){
                        insideprocessing.Fill("pairInfoBackgroundSameSign", "TOF wrong", 1.0);
                        continue;
                    } else if(!isdEdxOk&&isTOFOk){
                        insideprocessing.Fill("pairInfoBackgroundSameSign", "dEdx wrong", 1.0);
                        continue;
                    } else if(!isdEdxOk&&!isTOFOk){
                        insideprocessing.Fill("pairInfoBackgroundSameSign", "Both wrong", 1.0);
                        continue;
                    }

                    //chi2
                    for(size_t neg1 = 0; neg1<nParticlesExtended; neg1++){
                        for(size_t neg2 = 0; neg2<nParticlesExtended; neg2++){
                            tempPairName = particleNicks[neg1]+"_"+particleNicks[neg2];
                            chi2Map[tempPairName] = getChi2(vector_Track_negative[i], vector_Track_negative[j], neg1, neg2, sigmaMap[tempPairName]);
                        }
                    }

                    //mass tests on different pairs
                    if(chi2Map["K_pi"]<9){
                        vector_Track_negative[i]->getLorentzVector(negative_track, particleMass[Kaon]);
                        vector_Track_negative[j]->getLorentzVector(negative_track2, particleMass[Pion]);
                        mass = (negative_track+negative_track2).M();
                        eta = (negative_track+negative_track2).Eta();
                        pT = (negative_track+negative_track2).Pt();
                        correction = correction_coefficient(-2, Kaon, Pion, vector_Track_negative[i], vector_Track_negative[j], PV_position.Z(), PV_position.Z());
                        insideprocessing.Fill("MKpiChi2BcgSameSigneta", mass, eta, correction);
                        insideprocessing.Fill("MKpiChi2BcgSameSign", mass, correction);
                        insideprocessing.Fill("MKpiChi2BcgSameSignClose", mass, correction);
                        insideprocessing.Fill("MKpiChi2BcgSameSignpT", mass, pT, correction);
                    }
                    if(chi2Map["pi_K"]<9){
                        vector_Track_negative[i]->getLorentzVector(negative_track, particleMass[Pion]);
                        vector_Track_negative[j]->getLorentzVector(negative_track2, particleMass[Kaon]);
                        mass = (negative_track+negative_track2).M();
                        eta = (negative_track+negative_track2).Eta();
                        pT = (negative_track+negative_track2).Pt();
                        correction = correction_coefficient(-2, Pion, Kaon, vector_Track_negative[i], vector_Track_negative[j], PV_position.Z(), PV_position.Z());
                        insideprocessing.Fill("MpiKChi2BcgSameSign", mass, correction);
                        insideprocessing.Fill("MpiKChi2BcgSameSignClose", mass, correction);
                        insideprocessing.Fill("MpiKChi2BcgSameSigneta", mass, eta, correction);
                        insideprocessing.Fill("MpiKChi2BcgSameSignpT", mass, pT, correction);
                    }
                    if(chi2Map["p_pi"]<9){
                        vector_Track_negative[i]->getLorentzVector(negative_track, particleMass[Proton]);
                        vector_Track_negative[j]->getLorentzVector(negative_track2, particleMass[Pion]);
                        mass = (negative_track+negative_track2).M();
                        eta = (negative_track+negative_track2).Eta();
                        pT = (negative_track+negative_track2).Pt();
                        correction = correction_coefficient(-2, Proton, Pion, vector_Track_negative[i], vector_Track_negative[j], PV_position.Z(), PV_position.Z());
                        insideprocessing.Fill("MppiChi2BcgSameSign", mass, correction);
                        insideprocessing.Fill("MppiChi2BcgSameSigneta", mass, eta, correction);
                        insideprocessing.Fill("MppiChi2BcgSameSignpT", mass, pT, correction);
                    }
                    if(chi2Map["pi_p"]<9){
                        vector_Track_negative[i]->getLorentzVector(negative_track, particleMass[Pion]);
                        vector_Track_negative[j]->getLorentzVector(negative_track2, particleMass[Proton]);
                        mass = (negative_track+negative_track2).M();
                        eta = (negative_track+negative_track2).Eta();
                        pT = (negative_track+negative_track2).Pt();
                        correction = correction_coefficient(-2, Pion, Proton, vector_Track_negative[i], vector_Track_negative[j], PV_position.Z(), PV_position.Z());
                        insideprocessing.Fill("MpipChi2BcgSameSign", mass, correction);
                        insideprocessing.Fill("MpipChi2BcgSameSigneta", mass, eta, correction);
                        insideprocessing.Fill("MpipChi2BcgSameSignpT", mass, pT, correction);
                    }
                    if(chi2Map["K_K"]<9){
                        vector_Track_negative[i]->getLorentzVector(negative_track, particleMass[Kaon]);
                        vector_Track_negative[j]->getLorentzVector(negative_track2, particleMass[Kaon]);
                        mass = (negative_track+negative_track2).M();
                        eta = (negative_track+negative_track2).Eta();
                        pT = (negative_track+negative_track2).Pt();
                        correction = correction_coefficient(-2, Kaon, Kaon, vector_Track_negative[i], vector_Track_negative[j], PV_position.Z(), PV_position.Z());
                        insideprocessing.Fill("MKKChi2BcgSameSign", mass, correction);
                        insideprocessing.Fill("MKKChi2BcgSameSignClose", mass, correction);
                        insideprocessing.Fill("MKKChi2BcgSameSigneta", mass, eta, correction);
                        insideprocessing.Fill("MKKChi2BcgSameSignpT", mass, pT, correction);
                    }
                    if(chi2Map["pi_pi"]<9){
                        vector_Track_negative[i]->getLorentzVector(negative_track, particleMass[Pion]);
                        vector_Track_negative[j]->getLorentzVector(negative_track2, particleMass[Pion]);
                        mass = (negative_track+negative_track2).M();
                        eta = (negative_track+negative_track2).Eta();
                        pT = (negative_track+negative_track2).Pt();
                        correction = correction_coefficient(-2, Pion, Pion, vector_Track_negative[i], vector_Track_negative[j], PV_position.Z(), PV_position.Z());
                        insideprocessing.Fill("MpipiChi2BcgSameSign", mass, correction);
                        insideprocessing.Fill("MpipiChi2BcgSameSigneta", mass, eta, correction);
                        insideprocessing.Fill("MpipiChi2BcgSameSignpT", mass, pT, correction);
                    }
                    if(chi2Map["p_p"]<9){
                        vector_Track_negative[i]->getLorentzVector(negative_track, particleMass[Proton]);
                        vector_Track_negative[j]->getLorentzVector(negative_track2, particleMass[Proton]);
                        mass = (negative_track+negative_track2).M();
                        eta = (negative_track+negative_track2).Eta();
                        pT = (negative_track+negative_track2).Pt();
                        correction = correction_coefficient(-2, Proton, Proton, vector_Track_negative[i], vector_Track_negative[j], PV_position.Z(), PV_position.Z());
                        insideprocessing.Fill("MppChi2BcgSameSign", mass, correction);
                        insideprocessing.Fill("MppChi2BcgSameSigneta", mass, eta, correction);
                        insideprocessing.Fill("MppChi2BcgSameSignpT", mass, pT, correction);
                    }
                }
            }

            //########## BACKGROUND EXTRACTION (TRACK ROTATION) ###############

            //loop through identified particles (track rotation)
            for(long unsigned int i = 0; i<vector_Track_positive.size(); i++){
                for(long unsigned int j = 0; j<vector_Track_negative.size(); j++){
                    isdEdxOk = (vector_Track_positive[i]->getNhitsDEdx()>=15)&&(vector_Track_negative[j]->getNhitsDEdx()>=15);
                    isTOFOk = (vector_Track_positive[i]->getTofPathLength()>0)&&(vector_Track_positive[i]->getTofTime()>0)&&(vector_Track_negative[j]->getTofPathLength()>0)&&(vector_Track_negative[j]->getTofTime()>0);
                    if(isdEdxOk&&isTOFOk){
                        insideprocessing.Fill("pairInfoBackgroundTrackRotation", "OK", 1.0);
                    } else if(isdEdxOk&&!isTOFOk){
                        insideprocessing.Fill("pairInfoBackgroundTrackRotation", "TOF wrong", 1.0);
                        continue;
                    } else if(!isdEdxOk&&isTOFOk){
                        insideprocessing.Fill("pairInfoBackgroundTrackRotation", "dEdx wrong", 1.0);
                        continue;
                    } else if(!isdEdxOk&&!isTOFOk){
                        insideprocessing.Fill("pairInfoBackgroundTrackRotation", "Both wrong", 1.0);
                        continue;
                    }

                    //chi2
                    for(size_t pos = 0; pos<nParticlesExtended; pos++){
                        for(size_t neg = 0; neg<nParticlesExtended; neg++){
                            tempPairName = particleNicks[pos]+"_"+particleNicks[neg];
                            chi2Map[tempPairName] = getChi2(vector_Track_positive[i], vector_Track_negative[j], pos, neg, sigmaMap[tempPairName]);
                        }
                    }

                    //mass tests on different pairs
                    if(chi2Map["K_pi"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Kaon]);
                        vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Pion]);
                        //track rotation (rotating only the negative one)
                        negative_track.RotateZ(TMath::Pi());
                        //the rest as usual
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(0, Kaon, Pion, vector_Track_positive[i], vector_Track_negative[j], PV_position.Z(), PV_position.Z());
                        insideprocessing.Fill("MKpiChi2BcgTrackRotation", mass, correction);
                        insideprocessing.Fill("MKpiChi2BcgTrackRotationClose", mass, correction);
                        insideprocessing.Fill("MKpiChi2BcgTrackRotationeta", mass, eta, correction);
                        insideprocessing.Fill("MKpiChi2BcgTrackRotationpT", mass, pT, correction);
                    }
                    if(chi2Map["pi_K"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                        vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Kaon]);
                        //track rotation (rotating only the negative one)
                        negative_track.RotateZ(TMath::Pi());
                        //the rest as usual
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(0, Pion, Kaon, vector_Track_positive[i], vector_Track_negative[j], PV_position.Z(), PV_position.Z());
                        insideprocessing.Fill("MpiKChi2BcgTrackRotation", mass, correction);
                        insideprocessing.Fill("MpiKChi2BcgTrackRotationClose", mass, correction);
                        insideprocessing.Fill("MpiKChi2BcgTrackRotationeta", mass, eta, correction);
                        insideprocessing.Fill("MpiKChi2BcgTrackRotationpT", mass, pT, correction);
                    }
                    if(chi2Map["p_pi"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Proton]);
                        vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Pion]);
                        //track rotation (rotating only the negative one)
                        negative_track.RotateZ(TMath::Pi());
                        //the rest as usual
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(0, Proton, Pion, vector_Track_positive[i], vector_Track_negative[j], PV_position.Z(), PV_position.Z());
                        insideprocessing.Fill("MppiChi2BcgTrackRotation", mass, correction);
                        insideprocessing.Fill("MppiChi2BcgTrackRotationeta", mass, eta, correction);
                        insideprocessing.Fill("MppiChi2BcgTrackRotationpT", mass, pT, correction);
                    }
                    if(chi2Map["pi_p"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                        vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Proton]);
                        //track rotation (rotating only the negative one)
                        negative_track.RotateZ(TMath::Pi());
                        //the rest as usual
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(0, Pion, Proton, vector_Track_positive[i], vector_Track_negative[j], PV_position.Z(), PV_position.Z());
                        insideprocessing.Fill("MpipChi2BcgTrackRotation", mass, correction);
                        insideprocessing.Fill("MpipChi2BcgTrackRotationeta", mass, eta, correction);
                        insideprocessing.Fill("MpipChi2BcgTrackRotationpT", mass, pT, correction);
                    }
                    if(chi2Map["K_K"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Kaon]);
                        vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Kaon]);
                        //track rotation (rotating only the negative one)
                        negative_track.RotateZ(TMath::Pi());
                        //the rest as usual
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(0, Kaon, Kaon, vector_Track_positive[i], vector_Track_negative[j], PV_position.Z(), PV_position.Z());
                        insideprocessing.Fill("MKKChi2BcgTrackRotation", mass, correction);
                        insideprocessing.Fill("MKKChi2BcgTrackRotationClose", mass, correction);
                        insideprocessing.Fill("MKKChi2BcgTrackRotationeta", mass, eta, correction);
                        insideprocessing.Fill("MKKChi2BcgTrackRotationpT", mass, pT, correction);
                    }
                    if(chi2Map["pi_pi"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                        vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Pion]);
                        //track rotation (rotating only the negative one)
                        negative_track.RotateZ(TMath::Pi());
                        //the rest as usual
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(0, Pion, Pion, vector_Track_positive[i], vector_Track_negative[j], PV_position.Z(), PV_position.Z());
                        insideprocessing.Fill("MpipiChi2BcgTrackRotation", mass, correction);
                        insideprocessing.Fill("MpipiChi2BcgTrackRotationeta", mass, eta, correction);
                        insideprocessing.Fill("MpipiChi2BcgTrackRotationpT", mass, pT, correction);
                    }
                    if(chi2Map["p_p"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Proton]);
                        vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Proton]);
                        //track rotation (rotating only the negative one)
                        negative_track.RotateZ(TMath::Pi());
                        //the rest as usual
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(0, Proton, Proton, vector_Track_positive[i], vector_Track_negative[j], PV_position.Z(), PV_position.Z());
                        insideprocessing.Fill("MppChi2BcgTrackRotation", mass, correction);
                        insideprocessing.Fill("MppChi2BcgTrackRotationeta", mass, eta, correction);
                        insideprocessing.Fill("MppChi2BcgTrackRotationpT", mass, pT, correction);
                    }
                }
            }

            //########## BACKGROUND EXTRACTION (RANDOM TRACK ROTATION) ###############

            //loop through identified particles (track rotation)
            for(long unsigned int i = 0; i<vector_Track_positive.size(); i++){
                for(long unsigned int j = 0; j<vector_Track_negative.size(); j++){
                    isdEdxOk = (vector_Track_positive[i]->getNhitsDEdx()>=15)&&(vector_Track_negative[j]->getNhitsDEdx()>=15);
                    isTOFOk = (vector_Track_positive[i]->getTofPathLength()>0)&&(vector_Track_positive[i]->getTofTime()>0)&&(vector_Track_negative[j]->getTofPathLength()>0)&&(vector_Track_negative[j]->getTofTime()>0);
                    if(isdEdxOk&&isTOFOk){
                        insideprocessing.Fill("pairInfoBackgroundRandomTrackRotation", "OK", 1.0);
                    } else if(isdEdxOk&&!isTOFOk){
                        insideprocessing.Fill("pairInfoBackgroundRandomTrackRotation", "TOF wrong", 1.0);
                        continue;
                    } else if(!isdEdxOk&&isTOFOk){
                        insideprocessing.Fill("pairInfoBackgroundRandomTrackRotation", "dEdx wrong", 1.0);
                        continue;
                    } else if(!isdEdxOk&&!isTOFOk){
                        insideprocessing.Fill("pairInfoBackgroundRandomTrackRotation", "Both wrong", 1.0);
                        continue;
                    }

                    //chi2
                    for(size_t pos = 0; pos<nParticlesExtended; pos++){
                        for(size_t neg = 0; neg<nParticlesExtended; neg++){
                            tempPairName = particleNicks[pos]+"_"+particleNicks[neg];
                            chi2Map[tempPairName] = getChi2(vector_Track_positive[i], vector_Track_negative[j], pos, neg, sigmaMap[tempPairName]);
                        }
                    }

                    //mass tests on different pairs
                    if(chi2Map["K_pi"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Kaon]);
                        vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Pion]);
                        //track rotation (rotating only the negative one)
                        negative_track.RotateZ(random_generator.Uniform(2*TMath::Pi()));
                        //the rest as usual
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(0, Kaon, Pion, vector_Track_positive[i], vector_Track_negative[j], PV_position.Z(), PV_position.Z());
                        insideprocessing.Fill("MKpiChi2BcgRandomTrackRotation", mass, correction);
                        insideprocessing.Fill("MKpiChi2BcgRandomTrackRotationClose", mass, correction);
                        insideprocessing.Fill("MKpiChi2BcgRandomTrackRotationeta", mass, eta, correction);
                        insideprocessing.Fill("MKpiChi2BcgRandomTrackRotationpT", mass, pT, correction);
                    }
                    if(chi2Map["pi_K"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                        vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Kaon]);
                        //track rotation (rotating only the negative one)
                        negative_track.RotateZ(random_generator.Uniform(2*TMath::Pi()));
                        //the rest as usual
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(0, Pion, Kaon, vector_Track_positive[i], vector_Track_negative[j], PV_position.Z(), PV_position.Z());
                        insideprocessing.Fill("MpiKChi2BcgRandomTrackRotation", mass, correction);
                        insideprocessing.Fill("MpiKChi2BcgRandomTrackRotationClose", mass, correction);
                        insideprocessing.Fill("MpiKChi2BcgRandomTrackRotationeta", mass, eta, correction);
                        insideprocessing.Fill("MpiKChi2BcgRandomTrackRotationpT", mass, pT, correction);
                    }
                    if(chi2Map["p_pi"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Proton]);
                        vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Pion]);
                        //track rotation (rotating only the negative one)
                        negative_track.RotateZ(random_generator.Uniform(2*TMath::Pi()));
                        //the rest as usual
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(0, Proton, Pion, vector_Track_positive[i], vector_Track_negative[j], PV_position.Z(), PV_position.Z());
                        insideprocessing.Fill("MppiChi2BcgRandomTrackRotation", mass, correction);
                        insideprocessing.Fill("MppiChi2BcgRandomTrackRotationeta", mass, eta, correction);
                        insideprocessing.Fill("MppiChi2BcgRandomTrackRotationpT", mass, pT, correction);
                    }
                    if(chi2Map["pi_p"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                        vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Proton]);
                        //track rotation (rotating only the negative one)
                        negative_track.RotateZ(random_generator.Uniform(2*TMath::Pi()));
                        //the rest as usual
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(0, Pion, Proton, vector_Track_positive[i], vector_Track_negative[j], PV_position.Z(), PV_position.Z());
                        insideprocessing.Fill("MpipChi2BcgRandomTrackRotation", mass, correction);
                        insideprocessing.Fill("MpipChi2BcgRandomTrackRotationeta", mass, eta, correction);
                        insideprocessing.Fill("MpipChi2BcgRandomTrackRotationpT", mass, pT, correction);
                    }
                    if(chi2Map["K_K"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Kaon]);
                        vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Kaon]);
                        //track rotation (rotating only the negative one)
                        negative_track.RotateZ(random_generator.Uniform(2*TMath::Pi()));
                        //the rest as usual
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(0, Kaon, Kaon, vector_Track_positive[i], vector_Track_negative[j], PV_position.Z(), PV_position.Z());
                        insideprocessing.Fill("MKKChi2BcgRandomTrackRotation", mass, correction);
                        insideprocessing.Fill("MKKChi2BcgRandomTrackRotationClose", mass, correction);
                        insideprocessing.Fill("MKKChi2BcgRandomTrackRotationeta", mass, eta, correction);
                        insideprocessing.Fill("MKKChi2BcgRandomTrackRotationpT", mass, pT, correction);
                    }
                    if(chi2Map["pi_pi"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                        vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Pion]);
                        //track rotation (rotating only the negative one)
                        negative_track.RotateZ(random_generator.Uniform(2*TMath::Pi()));
                        //the rest as usual
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(0, Pion, Pion, vector_Track_positive[i], vector_Track_negative[j], PV_position.Z(), PV_position.Z());
                        insideprocessing.Fill("MpipiChi2BcgRandomTrackRotation", mass, correction);
                        insideprocessing.Fill("MpipiChi2BcgRandomTrackRotationeta", mass, eta, correction);
                        insideprocessing.Fill("MpipiChi2BcgRandomTrackRotationpT", mass, pT, correction);
                    }
                    if(chi2Map["p_p"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Proton]);
                        vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Proton]);
                        //track rotation (rotating only the negative one)
                        negative_track.RotateZ(random_generator.Uniform(2*TMath::Pi()));
                        //the rest as usual
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(0, Proton, Proton, vector_Track_positive[i], vector_Track_negative[j], PV_position.Z(), PV_position.Z());
                        insideprocessing.Fill("MppChi2BcgRandomTrackRotation", mass, correction);
                        insideprocessing.Fill("MppChi2BcgRandomTrackRotationeta", mass, eta, correction);
                        insideprocessing.Fill("MppChi2BcgRandomTrackRotationpT", mass, pT, correction);
                    }
                }
            }

            //########## BACKGROUND EXTRACTION (MIXED EVENT TECHNIQUES) ###############

            //TODO: have a mechanism to remove possible duplicates
            //as in, positive track matching to two negatives gets counted twice

            //picking pairs good for comparing with previous ones (and saving as ones)
            //we load the queue with the copy of the current event
            queue_of_previous_vector_Tracks_Kpi_positive.emplace_back();
            queue_of_previous_vector_Tracks_piK_positive.emplace_back();
            queue_of_previous_vector_Tracks_ppi_positive.emplace_back();
            queue_of_previous_vector_Tracks_pip_positive.emplace_back();
            queue_of_previous_vector_Tracks_KK_positive.emplace_back();
            queue_of_previous_vector_Tracks_pipi_positive.emplace_back();
            queue_of_previous_vector_Tracks_pp_positive.emplace_back();
            queue_of_previous_vector_Tracks_Kpi_negative.emplace_back();
            queue_of_previous_vector_Tracks_piK_negative.emplace_back();
            queue_of_previous_vector_Tracks_ppi_negative.emplace_back();
            queue_of_previous_vector_Tracks_pip_negative.emplace_back();
            queue_of_previous_vector_Tracks_KK_negative.emplace_back();
            queue_of_previous_vector_Tracks_pipi_negative.emplace_back();
            queue_of_previous_vector_Tracks_pp_negative.emplace_back();
            for(long unsigned int i = 0; i<vector_Track_positive.size(); i++){
                for(long unsigned int j = 0; j<vector_Track_negative.size(); j++){
                    isdEdxOk = (vector_Track_positive[i]->getNhitsDEdx()>=15)&&(vector_Track_negative[j]->getNhitsDEdx()>=15);
                    isTOFOk = (vector_Track_positive[i]->getTofPathLength()>0)&&(vector_Track_positive[i]->getTofTime()>0)&&(vector_Track_negative[j]->getTofPathLength()>0)&&(vector_Track_negative[j]->getTofTime()>0);
                    if(isdEdxOk&&isTOFOk){
                        ///pair is ok for harvesting
                    } else{
                        continue;
                    }

                    //chi2
                    for(size_t pos = 0; pos<nParticlesExtended; pos++){
                        for(size_t neg = 0; neg<nParticlesExtended; neg++){
                            tempPairName = particleNicks[pos]+"_"+particleNicks[neg];
                            chi2Map[tempPairName] = getChi2(vector_Track_positive[i], vector_Track_negative[j], pos, neg, sigmaMap[tempPairName]);
                        }
                    }

                    //useful variables
                    double pt, eta, phi;
                    //chi2 tests on different pairs
                    if(chi2Map["K_pi"]<9){
                        //positive
                        queue_of_previous_vector_Tracks_Kpi_positive.back().push_back(new StUPCTrack());
                        vector_Track_positive[i]->getPtEtaPhi(pt, eta, phi);
                        queue_of_previous_vector_Tracks_Kpi_positive.back().back()->setPtEtaPhi(pt, eta, phi);
                        queue_of_previous_vector_Tracks_Kpi_positive.back().back()->setNhitsDEdx(vector_Track_positive[i]->getNhitsDEdx());
                        queue_of_previous_vector_Tracks_Kpi_positive.back().back()->setTofPathLength(vector_Track_positive[i]->getTofPathLength());
                        queue_of_previous_vector_Tracks_Kpi_positive.back().back()->setTofTime(vector_Track_positive[i]->getTofTime());
                        for(size_t part = 0; part<nParticlesExtended; part++){
                            queue_of_previous_vector_Tracks_Kpi_positive.back().back()->setNSigmasTPC(static_cast<StUPCTrack::Part>(part), vector_Track_positive[i]->getNSigmasTPC(static_cast<StUPCTrack::Part>(part)));
                        }
                        //negative
                        queue_of_previous_vector_Tracks_Kpi_negative.back().push_back(new StUPCTrack());
                        vector_Track_negative[j]->getPtEtaPhi(pt, eta, phi);
                        queue_of_previous_vector_Tracks_Kpi_negative.back().back()->setPtEtaPhi(pt, eta, phi);
                        queue_of_previous_vector_Tracks_Kpi_negative.back().back()->setNhitsDEdx(vector_Track_negative[j]->getNhitsDEdx());
                        queue_of_previous_vector_Tracks_Kpi_negative.back().back()->setTofPathLength(vector_Track_negative[j]->getTofPathLength());
                        queue_of_previous_vector_Tracks_Kpi_negative.back().back()->setTofTime(vector_Track_negative[j]->getTofTime());
                        for(size_t part = 0; part<nParticlesExtended; part++){
                            queue_of_previous_vector_Tracks_Kpi_negative.back().back()->setNSigmasTPC(static_cast<StUPCTrack::Part>(part), vector_Track_negative[j]->getNSigmasTPC(static_cast<StUPCTrack::Part>(part)));
                        }
                    }
                    if(chi2Map["pi_K"]<9){
                        //positive
                        queue_of_previous_vector_Tracks_piK_positive.back().push_back(new StUPCTrack());
                        vector_Track_positive[i]->getPtEtaPhi(pt, eta, phi);
                        queue_of_previous_vector_Tracks_piK_positive.back().back()->setPtEtaPhi(pt, eta, phi);
                        queue_of_previous_vector_Tracks_piK_positive.back().back()->setNhitsDEdx(vector_Track_positive[i]->getNhitsDEdx());
                        queue_of_previous_vector_Tracks_piK_positive.back().back()->setTofPathLength(vector_Track_positive[i]->getTofPathLength());
                        queue_of_previous_vector_Tracks_piK_positive.back().back()->setTofTime(vector_Track_positive[i]->getTofTime());
                        for(size_t part = 0; part<nParticlesExtended; part++){
                            queue_of_previous_vector_Tracks_piK_positive.back().back()->setNSigmasTPC(static_cast<StUPCTrack::Part>(part), vector_Track_positive[i]->getNSigmasTPC(static_cast<StUPCTrack::Part>(part)));
                        }
                        //negative
                        queue_of_previous_vector_Tracks_piK_negative.back().push_back(new StUPCTrack());
                        vector_Track_negative[j]->getPtEtaPhi(pt, eta, phi);
                        queue_of_previous_vector_Tracks_piK_negative.back().back()->setPtEtaPhi(pt, eta, phi);
                        queue_of_previous_vector_Tracks_piK_negative.back().back()->setNhitsDEdx(vector_Track_negative[j]->getNhitsDEdx());
                        queue_of_previous_vector_Tracks_piK_negative.back().back()->setTofPathLength(vector_Track_negative[j]->getTofPathLength());
                        queue_of_previous_vector_Tracks_piK_negative.back().back()->setTofTime(vector_Track_negative[j]->getTofTime());
                        for(size_t part = 0; part<nParticlesExtended; part++){
                            queue_of_previous_vector_Tracks_piK_negative.back().back()->setNSigmasTPC(static_cast<StUPCTrack::Part>(part), vector_Track_negative[j]->getNSigmasTPC(static_cast<StUPCTrack::Part>(part)));
                        }
                    }
                    if(chi2Map["p_pi"]<9){
                        //positive
                        queue_of_previous_vector_Tracks_ppi_positive.back().push_back(new StUPCTrack());
                        vector_Track_positive[i]->getPtEtaPhi(pt, eta, phi);
                        queue_of_previous_vector_Tracks_ppi_positive.back().back()->setPtEtaPhi(pt, eta, phi);
                        queue_of_previous_vector_Tracks_ppi_positive.back().back()->setNhitsDEdx(vector_Track_positive[i]->getNhitsDEdx());
                        queue_of_previous_vector_Tracks_ppi_positive.back().back()->setTofPathLength(vector_Track_positive[i]->getTofPathLength());
                        queue_of_previous_vector_Tracks_ppi_positive.back().back()->setTofTime(vector_Track_positive[i]->getTofTime());
                        for(size_t part = 0; part<nParticlesExtended; part++){
                            queue_of_previous_vector_Tracks_ppi_positive.back().back()->setNSigmasTPC(static_cast<StUPCTrack::Part>(part), vector_Track_positive[i]->getNSigmasTPC(static_cast<StUPCTrack::Part>(part)));
                        }
                        //negative
                        queue_of_previous_vector_Tracks_ppi_negative.back().push_back(new StUPCTrack());
                        vector_Track_negative[j]->getPtEtaPhi(pt, eta, phi);
                        queue_of_previous_vector_Tracks_ppi_negative.back().back()->setPtEtaPhi(pt, eta, phi);
                        queue_of_previous_vector_Tracks_ppi_negative.back().back()->setNhitsDEdx(vector_Track_negative[j]->getNhitsDEdx());
                        queue_of_previous_vector_Tracks_ppi_negative.back().back()->setTofPathLength(vector_Track_negative[j]->getTofPathLength());
                        queue_of_previous_vector_Tracks_ppi_negative.back().back()->setTofTime(vector_Track_negative[j]->getTofTime());
                        for(size_t part = 0; part<nParticlesExtended; part++){
                            queue_of_previous_vector_Tracks_ppi_negative.back().back()->setNSigmasTPC(static_cast<StUPCTrack::Part>(part), vector_Track_negative[j]->getNSigmasTPC(static_cast<StUPCTrack::Part>(part)));
                        }
                    }
                    if(chi2Map["pi_p"]<9){
                        //positive
                        queue_of_previous_vector_Tracks_pip_positive.back().push_back(new StUPCTrack());
                        vector_Track_positive[i]->getPtEtaPhi(pt, eta, phi);
                        queue_of_previous_vector_Tracks_pip_positive.back().back()->setPtEtaPhi(pt, eta, phi);
                        queue_of_previous_vector_Tracks_pip_positive.back().back()->setNhitsDEdx(vector_Track_positive[i]->getNhitsDEdx());
                        queue_of_previous_vector_Tracks_pip_positive.back().back()->setTofPathLength(vector_Track_positive[i]->getTofPathLength());
                        queue_of_previous_vector_Tracks_pip_positive.back().back()->setTofTime(vector_Track_positive[i]->getTofTime());
                        for(size_t part = 0; part<nParticlesExtended; part++){
                            queue_of_previous_vector_Tracks_pip_positive.back().back()->setNSigmasTPC(static_cast<StUPCTrack::Part>(part), vector_Track_positive[i]->getNSigmasTPC(static_cast<StUPCTrack::Part>(part)));
                        }
                        //negative
                        queue_of_previous_vector_Tracks_pip_negative.back().push_back(new StUPCTrack());
                        vector_Track_negative[j]->getPtEtaPhi(pt, eta, phi);
                        queue_of_previous_vector_Tracks_pip_negative.back().back()->setPtEtaPhi(pt, eta, phi);
                        queue_of_previous_vector_Tracks_pip_negative.back().back()->setNhitsDEdx(vector_Track_negative[j]->getNhitsDEdx());
                        queue_of_previous_vector_Tracks_pip_negative.back().back()->setTofPathLength(vector_Track_negative[j]->getTofPathLength());
                        queue_of_previous_vector_Tracks_pip_negative.back().back()->setTofTime(vector_Track_negative[j]->getTofTime());
                        for(size_t part = 0; part<nParticlesExtended; part++){
                            queue_of_previous_vector_Tracks_pip_negative.back().back()->setNSigmasTPC(static_cast<StUPCTrack::Part>(part), vector_Track_negative[j]->getNSigmasTPC(static_cast<StUPCTrack::Part>(part)));
                        }
                    }
                    if(chi2Map["K_K"]<9){
                        //positive
                        queue_of_previous_vector_Tracks_KK_positive.back().push_back(new StUPCTrack());
                        vector_Track_positive[i]->getPtEtaPhi(pt, eta, phi);
                        queue_of_previous_vector_Tracks_KK_positive.back().back()->setPtEtaPhi(pt, eta, phi);
                        queue_of_previous_vector_Tracks_KK_positive.back().back()->setNhitsDEdx(vector_Track_positive[i]->getNhitsDEdx());
                        queue_of_previous_vector_Tracks_KK_positive.back().back()->setTofPathLength(vector_Track_positive[i]->getTofPathLength());
                        queue_of_previous_vector_Tracks_KK_positive.back().back()->setTofTime(vector_Track_positive[i]->getTofTime());
                        for(size_t part = 0; part<nParticlesExtended; part++){
                            queue_of_previous_vector_Tracks_KK_positive.back().back()->setNSigmasTPC(static_cast<StUPCTrack::Part>(part), vector_Track_positive[i]->getNSigmasTPC(static_cast<StUPCTrack::Part>(part)));
                        }
                        //negative
                        queue_of_previous_vector_Tracks_KK_negative.back().push_back(new StUPCTrack());
                        vector_Track_negative[j]->getPtEtaPhi(pt, eta, phi);
                        queue_of_previous_vector_Tracks_KK_negative.back().back()->setPtEtaPhi(pt, eta, phi);
                        queue_of_previous_vector_Tracks_KK_negative.back().back()->setNhitsDEdx(vector_Track_negative[j]->getNhitsDEdx());
                        queue_of_previous_vector_Tracks_KK_negative.back().back()->setTofPathLength(vector_Track_negative[j]->getTofPathLength());
                        queue_of_previous_vector_Tracks_KK_negative.back().back()->setTofTime(vector_Track_negative[j]->getTofTime());
                        for(size_t part = 0; part<nParticlesExtended; part++){
                            queue_of_previous_vector_Tracks_KK_negative.back().back()->setNSigmasTPC(static_cast<StUPCTrack::Part>(part), vector_Track_negative[j]->getNSigmasTPC(static_cast<StUPCTrack::Part>(part)));
                        }
                    }
                    if(chi2Map["pi_pi"]<9){
                        //positive
                        queue_of_previous_vector_Tracks_pipi_positive.back().push_back(new StUPCTrack());
                        vector_Track_positive[i]->getPtEtaPhi(pt, eta, phi);
                        queue_of_previous_vector_Tracks_pipi_positive.back().back()->setPtEtaPhi(pt, eta, phi);
                        queue_of_previous_vector_Tracks_pipi_positive.back().back()->setNhitsDEdx(vector_Track_positive[i]->getNhitsDEdx());
                        queue_of_previous_vector_Tracks_pipi_positive.back().back()->setTofPathLength(vector_Track_positive[i]->getTofPathLength());
                        queue_of_previous_vector_Tracks_pipi_positive.back().back()->setTofTime(vector_Track_positive[i]->getTofTime());
                        for(size_t part = 0; part<nParticlesExtended; part++){
                            queue_of_previous_vector_Tracks_pipi_positive.back().back()->setNSigmasTPC(static_cast<StUPCTrack::Part>(part), vector_Track_positive[i]->getNSigmasTPC(static_cast<StUPCTrack::Part>(part)));
                        }
                        //negative
                        queue_of_previous_vector_Tracks_pipi_negative.back().push_back(new StUPCTrack());
                        vector_Track_negative[j]->getPtEtaPhi(pt, eta, phi);
                        queue_of_previous_vector_Tracks_pipi_negative.back().back()->setPtEtaPhi(pt, eta, phi);
                        queue_of_previous_vector_Tracks_pipi_negative.back().back()->setNhitsDEdx(vector_Track_negative[j]->getNhitsDEdx());
                        queue_of_previous_vector_Tracks_pipi_negative.back().back()->setTofPathLength(vector_Track_negative[j]->getTofPathLength());
                        queue_of_previous_vector_Tracks_pipi_negative.back().back()->setTofTime(vector_Track_negative[j]->getTofTime());
                        for(size_t part = 0; part<nParticlesExtended; part++){
                            queue_of_previous_vector_Tracks_pipi_negative.back().back()->setNSigmasTPC(static_cast<StUPCTrack::Part>(part), vector_Track_negative[j]->getNSigmasTPC(static_cast<StUPCTrack::Part>(part)));
                        }
                    }
                    if(chi2Map["p_p"]<9){
                        //positive
                        queue_of_previous_vector_Tracks_pp_positive.back().push_back(new StUPCTrack());
                        vector_Track_positive[i]->getPtEtaPhi(pt, eta, phi);
                        queue_of_previous_vector_Tracks_pp_positive.back().back()->setPtEtaPhi(pt, eta, phi);
                        queue_of_previous_vector_Tracks_pp_positive.back().back()->setNhitsDEdx(vector_Track_positive[i]->getNhitsDEdx());
                        queue_of_previous_vector_Tracks_pp_positive.back().back()->setTofPathLength(vector_Track_positive[i]->getTofPathLength());
                        queue_of_previous_vector_Tracks_pp_positive.back().back()->setTofTime(vector_Track_positive[i]->getTofTime());
                        for(size_t part = 0; part<nParticlesExtended; part++){
                            queue_of_previous_vector_Tracks_pp_positive.back().back()->setNSigmasTPC(static_cast<StUPCTrack::Part>(part), vector_Track_positive[i]->getNSigmasTPC(static_cast<StUPCTrack::Part>(part)));
                        }
                        //negative
                        queue_of_previous_vector_Tracks_pp_negative.back().push_back(new StUPCTrack());
                        vector_Track_negative[j]->getPtEtaPhi(pt, eta, phi);
                        queue_of_previous_vector_Tracks_pp_negative.back().back()->setPtEtaPhi(pt, eta, phi);
                        queue_of_previous_vector_Tracks_pp_negative.back().back()->setNhitsDEdx(vector_Track_negative[j]->getNhitsDEdx());
                        queue_of_previous_vector_Tracks_pp_negative.back().back()->setTofPathLength(vector_Track_negative[j]->getTofPathLength());
                        queue_of_previous_vector_Tracks_pp_negative.back().back()->setTofTime(vector_Track_negative[j]->getTofTime());
                        for(size_t part = 0; part<nParticlesExtended; part++){
                            queue_of_previous_vector_Tracks_pp_negative.back().back()->setNSigmasTPC(static_cast<StUPCTrack::Part>(part), vector_Track_negative[j]->getNSigmasTPC(static_cast<StUPCTrack::Part>(part)));
                        }
                    }
                }
            }
            //noting current position of primary vertex and other values
            queue_of_previous_PV_positions.emplace_back();
            queue_of_previous_PV_positions.back().SetXYZ(tempUPCpointer->getVertex(0)->getPosX(), tempUPCpointer->getVertex(0)->getPosY(), tempUPCpointer->getVertex(0)->getPosZ());
            queue_of_previous_Mx.emplace_back(tempRPpointer->getTrack(0)->xi(beamMomentum)* tempRPpointer->getTrack(1)->xi(beamMomentum)*510.);
            if(tempRPpointer->getTrack(0)->branch()==WU||tempRPpointer->getTrack(0)->branch()==WD){
                queue_of_previous_Xi_W.emplace_back(tempRPpointer->getTrack(0)->xi(beamMomentum));
                queue_of_previous_Xi_E.emplace_back(tempRPpointer->getTrack(1)->xi(beamMomentum));
            } else{
                queue_of_previous_Xi_W.emplace_back(tempRPpointer->getTrack(1)->xi(beamMomentum));
                queue_of_previous_Xi_E.emplace_back(tempRPpointer->getTrack(0)->xi(beamMomentum));
            }

            //useful variables
            double pt, eta, phi;
            //matching this event (on the back, "last" one, n-1) with previous ones (0 to n-2)
            for(size_t evt = 0; evt<min(total_events_in_queue, (int)queue_of_previous_vector_Tracks_Kpi_positive.size())-1; evt++){
                //checking the difference in vertex position
                TVector3 difference = queue_of_previous_PV_positions.back()-queue_of_previous_PV_positions[evt];
                insideprocessing.Fill("PrimaryVertexPositionDifferenceX", difference.X());
                insideprocessing.Fill("PrimaryVertexPositionDifferenceY", difference.Y());
                insideprocessing.Fill("PrimaryVertexPositionDifferenceZ", difference.Z());
                insideprocessing.Fill("PrimaryVertexPositionDifferenceTotal", difference.Mag());
                //if the difference is bigger than (~3sigma of TPC resolution(3cm)) value of 10cm, we skip this event
                if(fabs(difference.Z())>10){
                    continue;
                }
                //check difference in Mx and xi
                double Mx_difference = queue_of_previous_Mx.back()-queue_of_previous_Mx[evt];
                insideprocessing.Fill("MxDifference", Mx_difference);
                insideprocessing.Fill("MxDifferenceClose", Mx_difference);
                double Xi_W_difference = queue_of_previous_Xi_W.back()-queue_of_previous_Xi_W[evt];
                double Xi_E_difference = queue_of_previous_Xi_E.back()-queue_of_previous_Xi_E[evt];
                insideprocessing.Fill("XiWDifference", Xi_W_difference);
                insideprocessing.Fill("XiEDifference", Xi_E_difference);
                insideprocessing.Fill("XiBothDifference", Xi_W_difference, Xi_E_difference);
                //if the difference is bigger than 1 GeV (arbitrarily chosen), we skip this event
                //currently commented because it did not seem to have much impact, and made comparison to MC difficult
                // if(fabs(Mx_difference)>1.){
                //     continue;
                // }


                //MIXED-EVENT DIFFERENT-SIGN
                //############################################################################
                //this positive, past negative
                //Kpi
                for(long unsigned int i = 0; i<queue_of_previous_vector_Tracks_Kpi_positive.back().size(); i++){
                    for(long unsigned int j = 0; j<queue_of_previous_vector_Tracks_Kpi_negative[evt].size(); j++){
                        queue_of_previous_vector_Tracks_Kpi_positive.back()[i]->getLorentzVector(positive_track, particleMass[Kaon]);
                        queue_of_previous_vector_Tracks_Kpi_negative[evt][j]->getLorentzVector(negative_track, particleMass[Pion]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(0, Kaon, Pion, queue_of_previous_vector_Tracks_Kpi_positive.back()[i], queue_of_previous_vector_Tracks_Kpi_negative[evt][j], queue_of_previous_PV_positions.back().Z(), queue_of_previous_PV_positions[evt].Z());
                        insideprocessing.Fill("MKpiChi2BcgMixedEvent", mass, correction);
                        insideprocessing.Fill("MKpiChi2BcgMixedEventClose", mass, correction);
                        insideprocessing.Fill("MKpiChi2BcgMixedEventeta", mass, eta, correction);
                        insideprocessing.Fill("MKpiChi2BcgMixedEventpT", mass, pT, correction);
                    }
                }
                //piK
                for(long unsigned int i = 0; i<queue_of_previous_vector_Tracks_piK_positive.back().size(); i++){
                    for(long unsigned int j = 0; j<queue_of_previous_vector_Tracks_piK_negative[evt].size(); j++){
                        queue_of_previous_vector_Tracks_piK_positive.back()[i]->getLorentzVector(positive_track, particleMass[Pion]);
                        queue_of_previous_vector_Tracks_piK_negative[evt][j]->getLorentzVector(negative_track, particleMass[Kaon]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(0, Pion, Kaon, queue_of_previous_vector_Tracks_piK_positive.back()[i], queue_of_previous_vector_Tracks_piK_negative[evt][j], queue_of_previous_PV_positions.back().Z(), queue_of_previous_PV_positions[evt].Z());
                        insideprocessing.Fill("MpiKChi2BcgMixedEvent", mass, correction);
                        insideprocessing.Fill("MpiKChi2BcgMixedEventClose", mass, correction);
                        insideprocessing.Fill("MpiKChi2BcgMixedEventeta", mass, eta, correction);
                        insideprocessing.Fill("MpiKChi2BcgMixedEventpT", mass, pT, correction);
                    }
                }
                //ppi
                for(long unsigned int i = 0; i<queue_of_previous_vector_Tracks_ppi_positive.back().size(); i++){
                    for(long unsigned int j = 0; j<queue_of_previous_vector_Tracks_ppi_negative[evt].size(); j++){
                        queue_of_previous_vector_Tracks_ppi_positive.back()[i]->getLorentzVector(positive_track, particleMass[Proton]);
                        queue_of_previous_vector_Tracks_ppi_negative[evt][j]->getLorentzVector(negative_track, particleMass[Pion]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(0, Proton, Pion, queue_of_previous_vector_Tracks_ppi_positive.back()[i], queue_of_previous_vector_Tracks_ppi_negative[evt][j], queue_of_previous_PV_positions.back().Z(), queue_of_previous_PV_positions[evt].Z());
                        insideprocessing.Fill("MppiChi2BcgMixedEvent", mass, correction);
                        insideprocessing.Fill("MppiChi2BcgMixedEventeta", mass, eta, correction);
                        insideprocessing.Fill("MppiChi2BcgMixedEventpT", mass, pT, correction);
                    }
                }
                //pip
                for(long unsigned int i = 0; i<queue_of_previous_vector_Tracks_pip_positive.back().size(); i++){
                    for(long unsigned int j = 0; j<queue_of_previous_vector_Tracks_pip_negative[evt].size(); j++){
                        queue_of_previous_vector_Tracks_pip_positive.back()[i]->getLorentzVector(positive_track, particleMass[Pion]);
                        queue_of_previous_vector_Tracks_pip_negative[evt][j]->getLorentzVector(negative_track, particleMass[Proton]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(0, Pion, Proton, queue_of_previous_vector_Tracks_pip_positive.back()[i], queue_of_previous_vector_Tracks_pip_negative[evt][j], queue_of_previous_PV_positions.back().Z(), queue_of_previous_PV_positions[evt].Z());
                        insideprocessing.Fill("MpipChi2BcgMixedEvent", mass, correction);
                        insideprocessing.Fill("MpipChi2BcgMixedEventeta", mass, eta, correction);
                        insideprocessing.Fill("MpipChi2BcgMixedEventpT", mass, pT, correction);
                    }
                }
                //KK
                for(long unsigned int i = 0; i<queue_of_previous_vector_Tracks_KK_positive.back().size(); i++){
                    for(long unsigned int j = 0; j<queue_of_previous_vector_Tracks_KK_negative[evt].size(); j++){
                        queue_of_previous_vector_Tracks_KK_positive.back()[i]->getLorentzVector(positive_track, particleMass[Kaon]);
                        queue_of_previous_vector_Tracks_KK_negative[evt][j]->getLorentzVector(negative_track, particleMass[Kaon]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(0, Kaon, Kaon, queue_of_previous_vector_Tracks_KK_positive.back()[i], queue_of_previous_vector_Tracks_KK_negative[evt][j], queue_of_previous_PV_positions.back().Z(), queue_of_previous_PV_positions[evt].Z());
                        insideprocessing.Fill("MKKChi2BcgMixedEvent", mass, correction);
                        insideprocessing.Fill("MKKChi2BcgMixedEventClose", mass, correction);
                        insideprocessing.Fill("MKKChi2BcgMixedEventeta", mass, eta, correction);
                        insideprocessing.Fill("MKKChi2BcgMixedEventpT", mass, pT, correction);
                    }
                }
                //pipi
                for(long unsigned int i = 0; i<queue_of_previous_vector_Tracks_pipi_positive.back().size(); i++){
                    for(long unsigned int j = 0; j<queue_of_previous_vector_Tracks_pipi_negative[evt].size(); j++){
                        queue_of_previous_vector_Tracks_pipi_positive.back()[i]->getLorentzVector(positive_track, particleMass[Pion]);
                        queue_of_previous_vector_Tracks_pipi_negative[evt][j]->getLorentzVector(negative_track, particleMass[Pion]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(0, Pion, Pion, queue_of_previous_vector_Tracks_pipi_positive.back()[i], queue_of_previous_vector_Tracks_pipi_negative[evt][j], queue_of_previous_PV_positions.back().Z(), queue_of_previous_PV_positions[evt].Z());
                        insideprocessing.Fill("MpipiChi2BcgMixedEvent", mass, correction);
                        insideprocessing.Fill("MpipiChi2BcgMixedEventeta", mass, eta, correction);
                        insideprocessing.Fill("MpipiChi2BcgMixedEventpT", mass, pT, correction);
                    }
                }
                //pp
                for(long unsigned int i = 0; i<queue_of_previous_vector_Tracks_pp_positive.back().size(); i++){
                    for(long unsigned int j = 0; j<queue_of_previous_vector_Tracks_pp_negative[evt].size(); j++){
                        queue_of_previous_vector_Tracks_pp_positive.back()[i]->getLorentzVector(positive_track, particleMass[Proton]);
                        queue_of_previous_vector_Tracks_pp_negative[evt][j]->getLorentzVector(negative_track, particleMass[Proton]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(0, Proton, Proton, queue_of_previous_vector_Tracks_pp_positive.back()[i], queue_of_previous_vector_Tracks_pp_negative[evt][j], queue_of_previous_PV_positions.back().Z(), queue_of_previous_PV_positions[evt].Z());
                        insideprocessing.Fill("MppChi2BcgMixedEvent", mass, correction);
                        insideprocessing.Fill("MppChi2BcgMixedEventeta", mass, eta, correction);
                        insideprocessing.Fill("MppChi2BcgMixedEventpT", mass, pT, correction);
                    }
                }

                //MIXED-EVENT DIFFERENT-SIGN
                //############################################################################
                //this negative, past positive
                //Kpi
                for(long unsigned int i = 0; i<queue_of_previous_vector_Tracks_Kpi_positive[evt].size(); i++){
                    for(long unsigned int j = 0; j<queue_of_previous_vector_Tracks_Kpi_negative.back().size(); j++){
                        queue_of_previous_vector_Tracks_Kpi_positive[evt][i]->getLorentzVector(positive_track, particleMass[Kaon]);
                        queue_of_previous_vector_Tracks_Kpi_negative.back()[j]->getLorentzVector(negative_track, particleMass[Pion]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(0, Kaon, Pion, queue_of_previous_vector_Tracks_Kpi_positive[evt][i], queue_of_previous_vector_Tracks_Kpi_negative.back()[j], queue_of_previous_PV_positions[evt].Z(), queue_of_previous_PV_positions.back().Z());
                        insideprocessing.Fill("MKpiChi2BcgMixedEvent", mass, correction);
                        insideprocessing.Fill("MKpiChi2BcgMixedEventClose", mass, correction);
                        insideprocessing.Fill("MKpiChi2BcgMixedEventeta", mass, eta, correction);
                        insideprocessing.Fill("MKpiChi2BcgMixedEventpT", mass, pT, correction);
                    }
                }
                //piK
                for(long unsigned int i = 0; i<queue_of_previous_vector_Tracks_piK_positive[evt].size(); i++){
                    for(long unsigned int j = 0; j<queue_of_previous_vector_Tracks_piK_negative.back().size(); j++){
                        queue_of_previous_vector_Tracks_piK_positive[evt][i]->getLorentzVector(positive_track, particleMass[Pion]);
                        queue_of_previous_vector_Tracks_piK_negative.back()[j]->getLorentzVector(negative_track, particleMass[Kaon]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(0, Pion, Kaon, queue_of_previous_vector_Tracks_piK_positive[evt][i], queue_of_previous_vector_Tracks_piK_negative.back()[j], queue_of_previous_PV_positions[evt].Z(), queue_of_previous_PV_positions.back().Z());
                        insideprocessing.Fill("MpiKChi2BcgMixedEvent", mass, correction);
                        insideprocessing.Fill("MpiKChi2BcgMixedEventClose", mass, correction);
                        insideprocessing.Fill("MpiKChi2BcgMixedEventeta", mass, eta, correction);
                        insideprocessing.Fill("MpiKChi2BcgMixedEventpT", mass, pT, correction);
                    }
                }
                //ppi
                for(long unsigned int i = 0; i<queue_of_previous_vector_Tracks_ppi_positive[evt].size(); i++){
                    for(long unsigned int j = 0; j<queue_of_previous_vector_Tracks_ppi_negative.back().size(); j++){
                        queue_of_previous_vector_Tracks_ppi_positive[evt][i]->getLorentzVector(positive_track, particleMass[Proton]);
                        queue_of_previous_vector_Tracks_ppi_negative.back()[j]->getLorentzVector(negative_track, particleMass[Pion]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(0, Proton, Pion, queue_of_previous_vector_Tracks_ppi_positive[evt][i], queue_of_previous_vector_Tracks_ppi_negative.back()[j], queue_of_previous_PV_positions[evt].Z(), queue_of_previous_PV_positions.back().Z());
                        insideprocessing.Fill("MppiChi2BcgMixedEvent", mass, correction);
                        insideprocessing.Fill("MppiChi2BcgMixedEventeta", mass, eta, correction);
                        insideprocessing.Fill("MppiChi2BcgMixedEventpT", mass, pT, correction);
                    }
                }
                //pip
                for(long unsigned int i = 0; i<queue_of_previous_vector_Tracks_pip_positive[evt].size(); i++){
                    for(long unsigned int j = 0; j<queue_of_previous_vector_Tracks_pip_negative.back().size(); j++){
                        queue_of_previous_vector_Tracks_pip_positive[evt][i]->getLorentzVector(positive_track, particleMass[Pion]);
                        queue_of_previous_vector_Tracks_pip_negative.back()[j]->getLorentzVector(negative_track, particleMass[Proton]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(0, Pion, Proton, queue_of_previous_vector_Tracks_pip_positive[evt][i], queue_of_previous_vector_Tracks_pip_negative.back()[j], queue_of_previous_PV_positions[evt].Z(), queue_of_previous_PV_positions.back().Z());
                        insideprocessing.Fill("MpipChi2BcgMixedEvent", mass, correction);
                        insideprocessing.Fill("MpipChi2BcgMixedEventeta", mass, eta, correction);
                        insideprocessing.Fill("MpipChi2BcgMixedEventpT", mass, pT, correction);
                    }
                }
                //KK
                for(long unsigned int i = 0; i<queue_of_previous_vector_Tracks_KK_positive[evt].size(); i++){
                    for(long unsigned int j = 0; j<queue_of_previous_vector_Tracks_KK_negative.back().size(); j++){
                        queue_of_previous_vector_Tracks_KK_positive[evt][i]->getLorentzVector(positive_track, particleMass[Kaon]);
                        queue_of_previous_vector_Tracks_KK_negative.back()[j]->getLorentzVector(negative_track, particleMass[Kaon]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(0, Kaon, Kaon, queue_of_previous_vector_Tracks_KK_positive[evt][i], queue_of_previous_vector_Tracks_KK_negative.back()[j], queue_of_previous_PV_positions[evt].Z(), queue_of_previous_PV_positions.back().Z());
                        insideprocessing.Fill("MKKChi2BcgMixedEvent", mass, correction);
                        insideprocessing.Fill("MKKChi2BcgMixedEventClose", mass, correction);
                        insideprocessing.Fill("MKKChi2BcgMixedEventeta", mass, eta, correction);
                        insideprocessing.Fill("MKKChi2BcgMixedEventpT", mass, pT, correction);
                    }
                }
                //pipi
                for(long unsigned int i = 0; i<queue_of_previous_vector_Tracks_pipi_positive[evt].size(); i++){
                    for(long unsigned int j = 0; j<queue_of_previous_vector_Tracks_pipi_negative.back().size(); j++){
                        queue_of_previous_vector_Tracks_pipi_positive[evt][i]->getLorentzVector(positive_track, particleMass[Pion]);
                        queue_of_previous_vector_Tracks_pipi_negative.back()[j]->getLorentzVector(negative_track, particleMass[Pion]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(0, Pion, Pion, queue_of_previous_vector_Tracks_pipi_positive[evt][i], queue_of_previous_vector_Tracks_pipi_negative.back()[j], queue_of_previous_PV_positions[evt].Z(), queue_of_previous_PV_positions.back().Z());
                        insideprocessing.Fill("MpipiChi2BcgMixedEvent", mass, correction);
                        insideprocessing.Fill("MpipiChi2BcgMixedEventeta", mass, eta, correction);
                        insideprocessing.Fill("MpipiChi2BcgMixedEventpT", mass, pT, correction);
                    }
                }
                //pp
                for(long unsigned int i = 0; i<queue_of_previous_vector_Tracks_pp_positive[evt].size(); i++){
                    for(long unsigned int j = 0; j<queue_of_previous_vector_Tracks_pp_negative.back().size(); j++){
                        queue_of_previous_vector_Tracks_pp_positive[evt][i]->getLorentzVector(positive_track, particleMass[Proton]);
                        queue_of_previous_vector_Tracks_pp_negative.back()[j]->getLorentzVector(negative_track, particleMass[Proton]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(0, Proton, Proton, queue_of_previous_vector_Tracks_pp_positive[evt][i], queue_of_previous_vector_Tracks_pp_negative.back()[j], queue_of_previous_PV_positions[evt].Z(), queue_of_previous_PV_positions.back().Z());
                        insideprocessing.Fill("MppChi2BcgMixedEvent", mass, correction);
                        insideprocessing.Fill("MppChi2BcgMixedEventeta", mass, eta, correction);
                        insideprocessing.Fill("MppChi2BcgMixedEventpT", mass, pT, correction);
                    }
                }

                //MIXED-EVENT SAME-SIGN
                //THERE ARE NO POSITIVE AND NEGATIVE TRACKS IN ONE 
                //IT IS JUST A REUSED NAME
                //############################################################################
                //both positive
                //Kpi/piK
                for(long unsigned int i = 0; i<queue_of_previous_vector_Tracks_Kpi_positive[evt].size(); i++){
                    for(long unsigned int j = 0; j<queue_of_previous_vector_Tracks_piK_positive.back().size(); j++){
                        queue_of_previous_vector_Tracks_Kpi_positive[evt][i]->getLorentzVector(positive_track, particleMass[Kaon]);
                        queue_of_previous_vector_Tracks_piK_positive.back()[j]->getLorentzVector(negative_track, particleMass[Pion]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(2, Kaon, Pion, queue_of_previous_vector_Tracks_Kpi_positive[evt][i], queue_of_previous_vector_Tracks_piK_positive.back()[j], queue_of_previous_PV_positions[evt].Z(), queue_of_previous_PV_positions.back().Z());
                        //Kpi
                        insideprocessing.Fill("MKpiChi2BcgMixedEventSameSign", mass, correction);
                        insideprocessing.Fill("MKpiChi2BcgMixedEventSameSignClose", mass, correction);
                        insideprocessing.Fill("MKpiChi2BcgMixedEventSameSigneta", mass, eta, correction);
                        insideprocessing.Fill("MKpiChi2BcgMixedEventSameSignpT", mass, pT, correction);
                        //piK
                        insideprocessing.Fill("MpiKChi2BcgMixedEventSameSign", mass, correction);
                        insideprocessing.Fill("MpiKChi2BcgMixedEventSameSignClose", mass, correction);
                        insideprocessing.Fill("MpiKChi2BcgMixedEventSameSigneta", mass, eta, correction);
                        insideprocessing.Fill("MpiKChi2BcgMixedEventSameSignpT", mass, pT, correction);
                    }
                }
                // for(long unsigned int i = 0; i<queue_of_previous_vector_Tracks_Kpi_positive.back().size(); i++){
                //     for(long unsigned int j = 0; j<queue_of_previous_vector_Tracks_piK_positive[evt].size(); j++){
                //         queue_of_previous_vector_Tracks_Kpi_positive.back()[i]->getLorentzVector(positive_track, particleMass[Kaon]);
                //         queue_of_previous_vector_Tracks_piK_positive[evt][j]->getLorentzVector(negative_track, particleMass[Pion]);
                //         mass = (positive_track+negative_track).M();
                //         eta = (positive_track+negative_track).Eta();
                //         pT = (positive_track+negative_track).Pt();
                //         correction = correction_coefficient(2, Kaon, Pion, queue_of_previous_vector_Tracks_Kpi_positive.back()[i], queue_of_previous_vector_Tracks_piK_positive[evt][j], queue_of_previous_PV_positions.back().Z(), queue_of_previous_PV_positions[evt].Z());
                //         //Kpi
                //         insideprocessing.Fill("MKpiChi2BcgMixedEventSameSign", mass, correction);
                //         insideprocessing.Fill("MKpiChi2BcgMixedEventSameSignClose", mass, correction);
                //         insideprocessing.Fill("MKpiChi2BcgMixedEventSameSigneta", mass, eta, correction);
                //         insideprocessing.Fill("MKpiChi2BcgMixedEventSameSignpT", mass, pT, correction);
                //         //piK
                //         insideprocessing.Fill("MpiKChi2BcgMixedEventSameSign", mass, correction);
                //         insideprocessing.Fill("MpiKChi2BcgMixedEventSameSignClose", mass, correction);
                //         insideprocessing.Fill("MpiKChi2BcgMixedEventSameSigneta", mass, eta, correction);
                //         insideprocessing.Fill("MpiKChi2BcgMixedEventSameSignpT", mass, pT, correction);
                //     }
                // }
                //ppi/pip
                for(long unsigned int i = 0; i<queue_of_previous_vector_Tracks_ppi_positive[evt].size(); i++){
                    for(long unsigned int j = 0; j<queue_of_previous_vector_Tracks_pip_positive.back().size(); j++){
                        queue_of_previous_vector_Tracks_ppi_positive[evt][i]->getLorentzVector(positive_track, particleMass[Proton]);
                        queue_of_previous_vector_Tracks_pip_positive.back()[j]->getLorentzVector(negative_track, particleMass[Pion]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(2, Proton, Pion, queue_of_previous_vector_Tracks_ppi_positive[evt][i], queue_of_previous_vector_Tracks_pip_positive.back()[j], queue_of_previous_PV_positions[evt].Z(), queue_of_previous_PV_positions.back().Z());
                        //ppi
                        insideprocessing.Fill("MppiChi2BcgMixedEventSameSign", mass, correction);
                        insideprocessing.Fill("MppiChi2BcgMixedEventSameSigneta", mass, eta, correction);
                        insideprocessing.Fill("MppiChi2BcgMixedEventSameSignpT", mass, pT, correction);
                        //pip
                        insideprocessing.Fill("MpipChi2BcgMixedEventSameSign", mass, correction);
                        insideprocessing.Fill("MpipChi2BcgMixedEventSameSigneta", mass, eta, correction);
                        insideprocessing.Fill("MpipChi2BcgMixedEventSameSignpT", mass, pT, correction);
                    }
                }
                // for(long unsigned int i = 0; i<queue_of_previous_vector_Tracks_ppi_positive.back().size(); i++){
                //     for(long unsigned int j = 0; j<queue_of_previous_vector_Tracks_pip_positive[evt].size(); j++){
                //         queue_of_previous_vector_Tracks_ppi_positive.back()[i]->getLorentzVector(positive_track, particleMass[Proton]);
                //         queue_of_previous_vector_Tracks_pip_positive[evt][j]->getLorentzVector(negative_track, particleMass[Pion]);
                //         mass = (positive_track+negative_track).M();
                //         eta = (positive_track+negative_track).Eta();
                //         pT = (positive_track+negative_track).Pt();
                //         correction = correction_coefficient(2, Proton, Pion, queue_of_previous_vector_Tracks_ppi_positive.back()[i], queue_of_previous_vector_Tracks_pip_positive[evt][j], queue_of_previous_PV_positions.back().Z(), queue_of_previous_PV_positions[evt].Z());
                //         //ppi
                //         insideprocessing.Fill("MppiChi2BcgMixedEventSameSign", mass, correction);
                //         insideprocessing.Fill("MppiChi2BcgMixedEventSameSigneta", mass, eta, correction);
                //         insideprocessing.Fill("MppiChi2BcgMixedEventSameSignpT", mass, pT, correction);
                //         //pip
                //         insideprocessing.Fill("MpipChi2BcgMixedEventSameSign", mass, correction);
                //         insideprocessing.Fill("MpipChi2BcgMixedEventSameSigneta", mass, eta, correction);
                //         insideprocessing.Fill("MpipChi2BcgMixedEventSameSignpT", mass, pT, correction);
                //     }
                // }
                //KK
                for(long unsigned int i = 0; i<queue_of_previous_vector_Tracks_KK_positive[evt].size(); i++){
                    for(long unsigned int j = 0; j<queue_of_previous_vector_Tracks_KK_positive.back().size(); j++){
                        queue_of_previous_vector_Tracks_KK_positive[evt][i]->getLorentzVector(positive_track, particleMass[Kaon]);
                        queue_of_previous_vector_Tracks_KK_positive.back()[j]->getLorentzVector(negative_track, particleMass[Kaon]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(2, Kaon, Kaon, queue_of_previous_vector_Tracks_KK_positive[evt][i], queue_of_previous_vector_Tracks_KK_positive.back()[j], queue_of_previous_PV_positions[evt].Z(), queue_of_previous_PV_positions.back().Z());
                        insideprocessing.Fill("MKKChi2BcgMixedEventSameSign", mass, correction);
                        insideprocessing.Fill("MKKChi2BcgMixedEventSameSignClose", mass, correction);
                        insideprocessing.Fill("MKKChi2BcgMixedEventSameSigneta", mass, eta, correction);
                        insideprocessing.Fill("MKKChi2BcgMixedEventSameSignpT", mass, pT, correction);
                    }
                }
                //pipi
                for(long unsigned int i = 0; i<queue_of_previous_vector_Tracks_pipi_positive[evt].size(); i++){
                    for(long unsigned int j = 0; j<queue_of_previous_vector_Tracks_pipi_positive.back().size(); j++){
                        queue_of_previous_vector_Tracks_pipi_positive[evt][i]->getLorentzVector(positive_track, particleMass[Pion]);
                        queue_of_previous_vector_Tracks_pipi_positive.back()[j]->getLorentzVector(negative_track, particleMass[Pion]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(2, Pion, Pion, queue_of_previous_vector_Tracks_pipi_positive[evt][i], queue_of_previous_vector_Tracks_pipi_positive.back()[j], queue_of_previous_PV_positions[evt].Z(), queue_of_previous_PV_positions.back().Z());
                        insideprocessing.Fill("MpipiChi2BcgMixedEventSameSign", mass, correction);
                        insideprocessing.Fill("MpipiChi2BcgMixedEventSameSigneta", mass, eta, correction);
                        insideprocessing.Fill("MpipiChi2BcgMixedEventSameSignpT", mass, pT, correction);
                    }
                }
                //pp
                for(long unsigned int i = 0; i<queue_of_previous_vector_Tracks_pp_positive[evt].size(); i++){
                    for(long unsigned int j = 0; j<queue_of_previous_vector_Tracks_pp_positive.back().size(); j++){
                        queue_of_previous_vector_Tracks_pp_positive[evt][i]->getLorentzVector(positive_track, particleMass[Proton]);
                        queue_of_previous_vector_Tracks_pp_positive.back()[j]->getLorentzVector(negative_track, particleMass[Proton]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(2, Proton, Proton, queue_of_previous_vector_Tracks_pp_positive[evt][i], queue_of_previous_vector_Tracks_pp_positive.back()[j], queue_of_previous_PV_positions[evt].Z(), queue_of_previous_PV_positions.back().Z());
                        insideprocessing.Fill("MppChi2BcgMixedEventSameSign", mass, correction);
                        insideprocessing.Fill("MppChi2BcgMixedEventSameSigneta", mass, eta, correction);
                        insideprocessing.Fill("MppChi2BcgMixedEventSameSignpT", mass, pT, correction);
                    }
                }

                //MIXED-EVENT SAME-SIGN
                //############################################################################
                //both negative
                //Kpi/piK
                for(long unsigned int i = 0; i<queue_of_previous_vector_Tracks_Kpi_negative[evt].size(); i++){
                    for(long unsigned int j = 0; j<queue_of_previous_vector_Tracks_piK_negative.back().size(); j++){
                        queue_of_previous_vector_Tracks_Kpi_negative[evt][i]->getLorentzVector(positive_track, particleMass[Kaon]);
                        queue_of_previous_vector_Tracks_piK_negative.back()[j]->getLorentzVector(negative_track, particleMass[Pion]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(-2, Kaon, Pion, queue_of_previous_vector_Tracks_Kpi_negative[evt][i], queue_of_previous_vector_Tracks_piK_negative.back()[j], queue_of_previous_PV_positions[evt].Z(), queue_of_previous_PV_positions.back().Z());
                        //Kpi
                        insideprocessing.Fill("MKpiChi2BcgMixedEventSameSign", mass, correction);
                        insideprocessing.Fill("MKpiChi2BcgMixedEventSameSignClose", mass, correction);
                        insideprocessing.Fill("MKpiChi2BcgMixedEventSameSigneta", mass, eta, correction);
                        insideprocessing.Fill("MKpiChi2BcgMixedEventSameSignpT", mass, pT, correction);
                        //piK
                        insideprocessing.Fill("MpiKChi2BcgMixedEventSameSign", mass, correction);
                        insideprocessing.Fill("MpiKChi2BcgMixedEventSameSignClose", mass, correction);
                        insideprocessing.Fill("MpiKChi2BcgMixedEventSameSigneta", mass, eta, correction);
                        insideprocessing.Fill("MpiKChi2BcgMixedEventSameSignpT", mass, pT, correction);
                    }
                }
                // for(long unsigned int i = 0; i<queue_of_previous_vector_Tracks_Kpi_negative.back().size(); i++){
                //     for(long unsigned int j = 0; j<queue_of_previous_vector_Tracks_piK_negative[evt].size(); j++){
                //         queue_of_previous_vector_Tracks_Kpi_negative.back()[i]->getLorentzVector(positive_track, particleMass[Kaon]);
                //         queue_of_previous_vector_Tracks_piK_negative[evt][j]->getLorentzVector(negative_track, particleMass[Pion]);
                //         mass = (positive_track+negative_track).M();
                //         eta = (positive_track+negative_track).Eta();
                //         pT = (positive_track+negative_track).Pt();
                //         correction = correction_coefficient(-2, Kaon, Pion, queue_of_previous_vector_Tracks_Kpi_negative.back()[i], queue_of_previous_vector_Tracks_piK_negative[evt][j], queue_of_previous_PV_positions.back().Z(), queue_of_previous_PV_positions[evt].Z());
                //         //Kpi
                //         insideprocessing.Fill("MKpiChi2BcgMixedEventSameSign", mass, correction);
                //         insideprocessing.Fill("MKpiChi2BcgMixedEventSameSignClose", mass, correction);
                //         insideprocessing.Fill("MKpiChi2BcgMixedEventSameSigneta", mass, eta, correction);
                //         insideprocessing.Fill("MKpiChi2BcgMixedEventSameSignpT", mass, pT, correction);
                //         //piK
                //         insideprocessing.Fill("MpiKChi2BcgMixedEventSameSign", mass, correction);
                //         insideprocessing.Fill("MpiKChi2BcgMixedEventSameSignClose", mass, correction);
                //         insideprocessing.Fill("MpiKChi2BcgMixedEventSameSigneta", mass, eta, correction);
                //         insideprocessing.Fill("MpiKChi2BcgMixedEventSameSignpT", mass, pT, correction);
                //     }
                // }
                //ppi/pip
                for(long unsigned int i = 0; i<queue_of_previous_vector_Tracks_ppi_negative[evt].size(); i++){
                    for(long unsigned int j = 0; j<queue_of_previous_vector_Tracks_pip_negative.back().size(); j++){
                        queue_of_previous_vector_Tracks_ppi_negative[evt][i]->getLorentzVector(positive_track, particleMass[Proton]);
                        queue_of_previous_vector_Tracks_pip_negative.back()[j]->getLorentzVector(negative_track, particleMass[Pion]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(-2, Proton, Pion, queue_of_previous_vector_Tracks_ppi_negative[evt][i], queue_of_previous_vector_Tracks_pip_negative.back()[j], queue_of_previous_PV_positions[evt].Z(), queue_of_previous_PV_positions.back().Z());
                        //ppi
                        insideprocessing.Fill("MppiChi2BcgMixedEventSameSign", mass, correction);
                        insideprocessing.Fill("MppiChi2BcgMixedEventSameSigneta", mass, eta, correction);
                        insideprocessing.Fill("MppiChi2BcgMixedEventSameSignpT", mass, pT, correction);
                        //pip
                        insideprocessing.Fill("MpipChi2BcgMixedEventSameSign", mass, correction);
                        insideprocessing.Fill("MpipChi2BcgMixedEventSameSigneta", mass, eta, correction);
                        insideprocessing.Fill("MpipChi2BcgMixedEventSameSignpT", mass, pT, correction);
                    }
                }
                // for(long unsigned int i = 0; i<queue_of_previous_vector_Tracks_ppi_negative.back().size(); i++){
                //     for(long unsigned int j = 0; j<queue_of_previous_vector_Tracks_pip_negative[evt].size(); j++){
                //         queue_of_previous_vector_Tracks_ppi_negative.back()[i]->getLorentzVector(positive_track, particleMass[Proton]);
                //         queue_of_previous_vector_Tracks_pip_negative[evt][j]->getLorentzVector(negative_track, particleMass[Pion]);
                //         mass = (positive_track+negative_track).M();
                //         eta = (positive_track+negative_track).Eta();
                //         pT = (positive_track+negative_track).Pt();
                //         correction = correction_coefficient(-2, Proton, Pion, queue_of_previous_vector_Tracks_ppi_negative.back()[i], queue_of_previous_vector_Tracks_pip_negative[evt][j], queue_of_previous_PV_positions.back().Z(), queue_of_previous_PV_positions[evt].Z());
                //         //ppi
                //         insideprocessing.Fill("MppiChi2BcgMixedEventSameSign", mass, correction);
                //         insideprocessing.Fill("MppiChi2BcgMixedEventSameSigneta", mass, eta, correction);
                //         insideprocessing.Fill("MppiChi2BcgMixedEventSameSignpT", mass, pT, correction);
                //         //pip
                //         insideprocessing.Fill("MpipChi2BcgMixedEventSameSign", mass, correction);
                //         insideprocessing.Fill("MpipChi2BcgMixedEventSameSigneta", mass, eta, correction);
                //         insideprocessing.Fill("MpipChi2BcgMixedEventSameSignpT", mass, pT, correction);
                //     }
                // }
                //KK
                for(long unsigned int i = 0; i<queue_of_previous_vector_Tracks_KK_negative[evt].size(); i++){
                    for(long unsigned int j = 0; j<queue_of_previous_vector_Tracks_KK_negative.back().size(); j++){
                        queue_of_previous_vector_Tracks_KK_negative[evt][i]->getLorentzVector(positive_track, particleMass[Kaon]);
                        queue_of_previous_vector_Tracks_KK_negative.back()[j]->getLorentzVector(negative_track, particleMass[Kaon]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(-2, Kaon, Kaon, queue_of_previous_vector_Tracks_KK_negative[evt][i], queue_of_previous_vector_Tracks_KK_negative.back()[j], queue_of_previous_PV_positions[evt].Z(), queue_of_previous_PV_positions.back().Z());
                        insideprocessing.Fill("MKKChi2BcgMixedEventSameSign", mass, correction);
                        insideprocessing.Fill("MKKChi2BcgMixedEventSameSignClose", mass, correction);
                        insideprocessing.Fill("MKKChi2BcgMixedEventSameSigneta", mass, eta, correction);
                        insideprocessing.Fill("MKKChi2BcgMixedEventSameSignpT", mass, pT, correction);
                    }
                }
                //pipi
                for(long unsigned int i = 0; i<queue_of_previous_vector_Tracks_pipi_negative[evt].size(); i++){
                    for(long unsigned int j = 0; j<queue_of_previous_vector_Tracks_pipi_negative.back().size(); j++){
                        queue_of_previous_vector_Tracks_pipi_negative[evt][i]->getLorentzVector(positive_track, particleMass[Pion]);
                        queue_of_previous_vector_Tracks_pipi_negative.back()[j]->getLorentzVector(negative_track, particleMass[Pion]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(-2, Pion, Pion, queue_of_previous_vector_Tracks_pipi_negative[evt][i], queue_of_previous_vector_Tracks_pipi_negative.back()[j], queue_of_previous_PV_positions[evt].Z(), queue_of_previous_PV_positions.back().Z());
                        insideprocessing.Fill("MpipiChi2BcgMixedEventSameSign", mass, correction);
                        insideprocessing.Fill("MpipiChi2BcgMixedEventSameSigneta", mass, eta, correction);
                        insideprocessing.Fill("MpipiChi2BcgMixedEventSameSignpT", mass, pT, correction);
                    }
                }
                //pp
                for(long unsigned int i = 0; i<queue_of_previous_vector_Tracks_pp_negative[evt].size(); i++){
                    for(long unsigned int j = 0; j<queue_of_previous_vector_Tracks_pp_negative.back().size(); j++){
                        queue_of_previous_vector_Tracks_pp_negative[evt][i]->getLorentzVector(positive_track, particleMass[Proton]);
                        queue_of_previous_vector_Tracks_pp_negative.back()[j]->getLorentzVector(negative_track, particleMass[Proton]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        correction = correction_coefficient(-2, Proton, Proton, queue_of_previous_vector_Tracks_pp_negative[evt][i], queue_of_previous_vector_Tracks_pp_negative.back()[j], queue_of_previous_PV_positions[evt].Z(), queue_of_previous_PV_positions.back().Z());
                        insideprocessing.Fill("MppChi2BcgMixedEventSameSign", mass, correction);
                        insideprocessing.Fill("MppChi2BcgMixedEventSameSigneta", mass, eta, correction);
                        insideprocessing.Fill("MppChi2BcgMixedEventSameSignpT", mass, pT, correction);
                    }
                }
            }

            //if the recall limit has been reached, we pop the oldest one
            //(to check if its happening, we check the size of a random pair queue, as all of them have the same length)
            if(queue_of_previous_vector_Tracks_Kpi_positive.size()>total_events_in_queue){
                //Kpi
                for(size_t i = 0; i<queue_of_previous_vector_Tracks_Kpi_positive.front().size(); i++){
                    delete queue_of_previous_vector_Tracks_Kpi_positive.front()[i];
                }
                queue_of_previous_vector_Tracks_Kpi_positive.pop_front();
                for(size_t i = 0; i<queue_of_previous_vector_Tracks_Kpi_negative.front().size(); i++){
                    delete queue_of_previous_vector_Tracks_Kpi_negative.front()[i];
                }
                queue_of_previous_vector_Tracks_Kpi_negative.pop_front();
                //piK
                for(size_t i = 0; i<queue_of_previous_vector_Tracks_piK_positive.front().size(); i++){
                    delete queue_of_previous_vector_Tracks_piK_positive.front()[i];
                }
                queue_of_previous_vector_Tracks_piK_positive.pop_front();
                for(size_t i = 0; i<queue_of_previous_vector_Tracks_piK_negative.front().size(); i++){
                    delete queue_of_previous_vector_Tracks_piK_negative.front()[i];
                }
                queue_of_previous_vector_Tracks_piK_negative.pop_front();
                //ppi
                for(size_t i = 0; i<queue_of_previous_vector_Tracks_ppi_positive.front().size(); i++){
                    delete queue_of_previous_vector_Tracks_ppi_positive.front()[i];
                }
                queue_of_previous_vector_Tracks_ppi_positive.pop_front();
                for(size_t i = 0; i<queue_of_previous_vector_Tracks_ppi_negative.front().size(); i++){
                    delete queue_of_previous_vector_Tracks_ppi_negative.front()[i];
                }
                queue_of_previous_vector_Tracks_ppi_negative.pop_front();
                //pip
                for(size_t i = 0; i<queue_of_previous_vector_Tracks_pip_positive.front().size(); i++){
                    delete queue_of_previous_vector_Tracks_pip_positive.front()[i];
                }
                queue_of_previous_vector_Tracks_pip_positive.pop_front();
                for(size_t i = 0; i<queue_of_previous_vector_Tracks_pip_negative.front().size(); i++){
                    delete queue_of_previous_vector_Tracks_pip_negative.front()[i];
                }
                queue_of_previous_vector_Tracks_pip_negative.pop_front();
                //KK
                for(size_t i = 0; i<queue_of_previous_vector_Tracks_KK_positive.front().size(); i++){
                    delete queue_of_previous_vector_Tracks_KK_positive.front()[i];
                }
                queue_of_previous_vector_Tracks_KK_positive.pop_front();
                for(size_t i = 0; i<queue_of_previous_vector_Tracks_KK_negative.front().size(); i++){
                    delete queue_of_previous_vector_Tracks_KK_negative.front()[i];
                }
                queue_of_previous_vector_Tracks_KK_negative.pop_front();
                //pipi
                for(size_t i = 0; i<queue_of_previous_vector_Tracks_pipi_positive.front().size(); i++){
                    delete queue_of_previous_vector_Tracks_pipi_positive.front()[i];
                }
                queue_of_previous_vector_Tracks_pipi_positive.pop_front();
                for(size_t i = 0; i<queue_of_previous_vector_Tracks_pipi_negative.front().size(); i++){
                    delete queue_of_previous_vector_Tracks_pipi_negative.front()[i];
                }
                queue_of_previous_vector_Tracks_pipi_negative.pop_front();
                //pp
                for(size_t i = 0; i<queue_of_previous_vector_Tracks_pp_positive.front().size(); i++){
                    delete queue_of_previous_vector_Tracks_pp_positive.front()[i];
                }
                queue_of_previous_vector_Tracks_pp_positive.pop_front();
                for(size_t i = 0; i<queue_of_previous_vector_Tracks_pp_negative.front().size(); i++){
                    delete queue_of_previous_vector_Tracks_pp_negative.front()[i];
                }
                queue_of_previous_vector_Tracks_pp_negative.pop_front();
                //primary vertex position
                queue_of_previous_PV_positions.pop_front();
            }

            //lambda finish
        }
        return 0;
    };

    TreeProc.Process(myFunction);

    //merging and tidying up
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
    //saving multicore histgrams and a special one
    outsideprocessing.SaveToFile(outputFileHist);

    outputFileHist->Close();

    return 0;
}

TFile* getNewestFile(std::string folder, std::string filename){
    //mounting path
    std::string full_path;
    if(folder.find_last_of("/")==folder.length()-1){
        full_path = folder+filename;
    } else{
        full_path = folder+"/"+filename;
    }
    //necessary to expand "~" into /home/adam/ or whatever will it be
    full_path = gSystem->ExpandPathName(full_path.c_str());

    //making file collection
    TFileCollection collection;
    int files_found = collection.Add(full_path.c_str());
    std::string result;
    Long_t id, size, flags, modtime = 0, new_modtime;
    switch(files_found){
    case 0:
        //no files
        printf("Search for %s returned 0 results\n", full_path.c_str());
        return nullptr;
        break;
    case 1:
        //exactly one file
        result = static_cast<TFileInfo*>(collection.GetList()->At(0))->GetCurrentUrl()->GetFile();
        printf("Search for %s returned 1 result:\n%s\n", full_path.c_str(), result.c_str());
        return TFile::Open(result.c_str());
        break;
    default:
        //many files; need to iterate through the file collection
        printf("Search for %s returned %d results:\n", full_path.c_str(), collection.GetList()->GetEntries());
        for(auto&& i:*collection.GetList()){
            gSystem->GetPathInfo(static_cast<TFileInfo*>(i)->GetCurrentUrl()->GetFile(), &id, &size, &flags, &new_modtime);
            printf("%s\n", static_cast<TFileInfo*>(i)->GetCurrentUrl()->GetFile());
            if(modtime<new_modtime){
                modtime = new_modtime;
                result = static_cast<TFileInfo*>(i)->GetCurrentUrl()->GetFile();
            }
        }
        printf("Final result chosen:\n%s\n", result.c_str());
        return TFile::Open(result.c_str());
        break;
    }

    //just in case, the program should NOT be even here
    return nullptr;
}

StEfficiencyCorrector3D* getInitialisedEfficiencyCorrector(std::string folder, CHARGE charge, PARTICLES particle){
    //loading file
    TFile* correction_file = nullptr;
    switch(particle){
    case Pion:
        correction_file = getNewestFile(folder, "SPPion*.root");
        break;
    case Kaon:
        correction_file = getNewestFile(folder, "SPKaon*.root");
        break;
    case Proton:
        correction_file = getNewestFile(folder, "SPProton*.root");
        break;
    default:
        return nullptr;
        break;
    }
    //if file not loaded, return nullptr
    if(correction_file==nullptr)
        return nullptr;
    printf("Successfully loaded efficiency correction file\n");

    //getting efficiency histograms out
    TH3F* TPC_num = nullptr;
    TH3F* TPC_den = nullptr;
    TH3F* TOF_num = nullptr;
    TH3F* TOF_den = nullptr;
    switch(charge){
    case positive:
        TPC_num = (TH3F*)correction_file->Get("h3D_TPC_RecoMatched_P");
        TPC_den = (TH3F*)correction_file->Get("h3D_TPC_True_P");
        TOF_num = (TH3F*)correction_file->Get("h3D_TOF_RecoMatchedWithTOF_P");
        TOF_den = (TH3F*)correction_file->Get("h3D_TOF_RecoMatched_P");
        break;
    case negative:
        TPC_num = (TH3F*)correction_file->Get("h3D_TPC_RecoMatched_N");
        TPC_den = (TH3F*)correction_file->Get("h3D_TPC_True_N");
        TOF_num = (TH3F*)correction_file->Get("h3D_TOF_RecoMatchedWithTOF_N");
        TOF_den = (TH3F*)correction_file->Get("h3D_TOF_RecoMatched_N");
        break;
    default:
        return nullptr;
        break;
    }
    if(TPC_num!=nullptr&&TPC_den!=nullptr&&TOF_num!=nullptr&&TOF_den!=nullptr){
        printf("Successfully loaded histograms\n");
    } else{
        printf("Something went wrong with loading histograms!\n");
        return nullptr;
    }

    //calculating efficiency
    TH3F* tpcEfficiency = (TH3F*)TPC_num->Clone("tpcEfficiency");
    TH3F* tofEfficiency = (TH3F*)TOF_num->Clone("tofEfficiency");

    //binomial division for proper error handling - taken directly from Sneha's example
    tpcEfficiency->Divide(TPC_num, TPC_den, 1, 1, "B");
    tofEfficiency->Divide(TOF_num, TOF_den, 1, 1, "B");

    //setting efficiencies
    StEfficiencyCorrector3D* efficiency_corrector = new StEfficiencyCorrector3D();
    bool cloneHist = true;  // Make internal copy, I guess for when the function goes out of scope?
    switch(particle){
    case Pion:
        efficiency_corrector->setTpcEfficiency(tpcEfficiency, charge==positive ? 1 : -1, StEfficiencyCorrector3D::PION, cloneHist);
        efficiency_corrector->setTofEfficiency(tofEfficiency, charge==positive ? 1 : -1, StEfficiencyCorrector3D::PION, cloneHist);
        break;
    case Kaon:
        efficiency_corrector->setTpcEfficiency(tpcEfficiency, charge==positive ? 1 : -1, StEfficiencyCorrector3D::KAON, cloneHist);
        efficiency_corrector->setTofEfficiency(tofEfficiency, charge==positive ? 1 : -1, StEfficiencyCorrector3D::KAON, cloneHist);
        break;
    case Proton:
        efficiency_corrector->setTpcEfficiency(tpcEfficiency, charge==positive ? 1 : -1, StEfficiencyCorrector3D::PROTON, cloneHist);
        efficiency_corrector->setTofEfficiency(tofEfficiency, charge==positive ? 1 : -1, StEfficiencyCorrector3D::PROTON, cloneHist);
        break;
    default:
        delete efficiency_corrector;
        return nullptr;
        break;
    }
    printf("Successfully set up efficiencies\n");

    //cleaning up
    //uncommenting it causes segfault because of nullptr somewhere down the line
    // correction_file->Close();
    // delete correction_file;

    return efficiency_corrector;
}