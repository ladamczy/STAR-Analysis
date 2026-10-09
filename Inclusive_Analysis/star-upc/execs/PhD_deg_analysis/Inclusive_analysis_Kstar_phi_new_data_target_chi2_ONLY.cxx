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
std::pair<PARTICLES, PARTICLES> getPairPID(std::string pair_with_underscore);
std::string reversePairAndRemoveUnderscore(std::string pair_with_underscore);

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
    printf("Efficiency corrections folder: %s\n", efficiency_folder.length()!=0 ? efficiency_folder.c_str() : "Not used");
    //because number of events includes the current one, we need to add 1
    total_events_in_queue++;

    //setting up efficiency correction
    double min_efficiency = 1e-4;
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
    //efficiency histograms
    TH1D efficiency1("efficiency1", "Particle reconstruction efficiency (0.0-1.0);efficiency;particles", 100, 0., 1.);
    TH1D efficiency2("efficiency2", "Particle reconstruction efficiency (10^{-2}-10^{-4});efficiency;particles", 100, 1e-2, 1e-4);
    TH1D efficiency3("efficiency3", "Particle reconstruction efficiency (10^{-4}-10^{-6});efficiency;particles", 100, 1e-4, 1e-6);
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

        efficiency1.Fill(efficiency);
        efficiency2.Fill(efficiency);
        efficiency3.Fill(efficiency);

        if(efficiency<min_efficiency)
            return 1./min_efficiency;
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
    std::vector<std::string> pairTabWithUnderscores = { "K_pi", "pi_K", "p_pi", "pi_p", "K_K", "pi_pi", "p_p" };
    std::vector<std::string> pairTab;
    std::vector<std::pair<PARTICLES, PARTICLES>> pairTabPID;
    for(auto&& temp_pair:pairTabWithUnderscores){
        std::string temp_pair_copy = temp_pair;
        temp_pair_copy.erase(std::remove(temp_pair_copy.begin(), temp_pair_copy.end(), '_'), temp_pair_copy.end());
        pairTab.push_back(temp_pair_copy);
        pairTabPID.push_back(getPairPID(temp_pair));
    }


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
        std::map<std::string, std::deque<std::vector<StUPCTrack*>>> map_of_queue_of_previous_vector_Tracks_positive;
        std::map<std::string, std::deque<std::vector<StUPCTrack*>>> map_of_queue_of_previous_vector_Tracks_negative;
        for(auto&& temp_pair_no_underscores:pairTab){
            map_of_queue_of_previous_vector_Tracks_positive[temp_pair_no_underscores];
            map_of_queue_of_previous_vector_Tracks_negative[temp_pair_no_underscores];
        }

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

                    for(size_t pair_number = 0; pair_number<pairTabWithUnderscores.size(); pair_number++){
                        if(chi2Map[pairTabWithUnderscores[pair_number]]<9){
                            //getting all the indentificators out
                            std::pair<PARTICLES, PARTICLES> pairPID = pairTabPID[pair_number];
                            std::string temp_pair = pairTab[pair_number];

                            //calculating physics (general)
                            vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[pairPID.first]);
                            vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[pairPID.second]);
                            mass = (positive_track+negative_track).M();
                            eta = (positive_track+negative_track).Eta();
                            pT = (positive_track+negative_track).Pt();
                            correction = correction_coefficient(0, pairPID.first, pairPID.second, vector_Track_positive[i], vector_Track_negative[j], PV_position.Z(), PV_position.Z());
                            insideprocessing.Fill(("M"+temp_pair+"Chi2").c_str(), mass, correction);
                            insideprocessing.Fill(("M"+temp_pair+"Chi2eta").c_str(), mass, eta, correction);
                            insideprocessing.Fill(("M"+temp_pair+"Chi2pT").c_str(), mass, pT, correction);

                            //calculating physics (specific)
                            if(temp_pair=="Kpi"||temp_pair=="piK"||temp_pair=="KK"){
                                insideprocessing.Fill(("M"+temp_pair+"Chi2Close").c_str(), mass, correction);
                            }
                            if(temp_pair=="KK"){
                                //test of suspicious peak and its neighbourhood
                                if(chi2Map["K_K"]<3){
                                    insideprocessing.Fill("MKKSuspiciousPeakTestedWithStrictChi2LessThan3", mass, correction);
                                }
                                if(chi2Map["K_K"]<1){
                                    insideprocessing.Fill("MKKSuspiciousPeakTestedWithStrictChi2LessThan1", mass, correction);
                                }
                                //!!!!!!!!!!
                                //CORRECTION CHANGE!!!!!!!!
                                //!!!!!!!!!!
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
                            //end of chi2 test
                        }
                        //end of pair loop
                    }
                    //end of innermost vector pair loop
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

                    for(size_t pair_number = 0; pair_number<pairTabWithUnderscores.size(); pair_number++){
                        if(chi2Map[pairTabWithUnderscores[pair_number]]<9){
                            //getting all the indentificators out
                            std::pair<PARTICLES, PARTICLES> pairPID = pairTabPID[pair_number];
                            std::string temp_pair = pairTab[pair_number];

                            //calculating physics (general)
                            vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[pairPID.first]);
                            vector_Track_positive[j]->getLorentzVector(positive_track2, particleMass[pairPID.second]);
                            mass = (positive_track+positive_track2).M();
                            eta = (positive_track+positive_track2).Eta();
                            pT = (positive_track+positive_track2).Pt();
                            correction = correction_coefficient(2, pairPID.first, pairPID.second, vector_Track_positive[i], vector_Track_positive[j], PV_position.Z(), PV_position.Z());
                            insideprocessing.Fill(("M"+temp_pair+"Chi2BcgSameSign").c_str(), mass, correction);
                            insideprocessing.Fill(("M"+temp_pair+"Chi2BcgSameSigneta").c_str(), mass, eta, correction);
                            insideprocessing.Fill(("M"+temp_pair+"Chi2BcgSameSignpT").c_str(), mass, pT, correction);

                            //calculating physics (specific)
                            if(temp_pair=="Kpi"||temp_pair=="piK"||temp_pair=="KK"){
                                insideprocessing.Fill(("M"+temp_pair+"Chi2BcgSameSignClose").c_str(), mass, correction);
                            }
                            //end of chi2 test
                        }
                        //end of pair loop
                    }
                    //end of innermost vector pair loop
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

                    for(size_t pair_number = 0; pair_number<pairTabWithUnderscores.size(); pair_number++){
                        if(chi2Map[pairTabWithUnderscores[pair_number]]<9){
                            //getting all the indentificators out
                            std::pair<PARTICLES, PARTICLES> pairPID = pairTabPID[pair_number];
                            std::string temp_pair = pairTab[pair_number];

                            //calculating physics (general)
                            vector_Track_negative[i]->getLorentzVector(negative_track, particleMass[pairPID.first]);
                            vector_Track_negative[j]->getLorentzVector(negative_track2, particleMass[pairPID.second]);
                            mass = (negative_track+negative_track2).M();
                            eta = (negative_track+negative_track2).Eta();
                            pT = (negative_track+negative_track2).Pt();
                            correction = correction_coefficient(2, pairPID.first, pairPID.second, vector_Track_negative[i], vector_Track_negative[j], PV_position.Z(), PV_position.Z());
                            insideprocessing.Fill(("M"+temp_pair+"Chi2BcgSameSign").c_str(), mass, correction);
                            insideprocessing.Fill(("M"+temp_pair+"Chi2BcgSameSigneta").c_str(), mass, eta, correction);
                            insideprocessing.Fill(("M"+temp_pair+"Chi2BcgSameSignpT").c_str(), mass, pT, correction);

                            //calculating physics (specific)
                            if(temp_pair=="Kpi"||temp_pair=="piK"||temp_pair=="KK"){
                                insideprocessing.Fill(("M"+temp_pair+"Chi2BcgSameSignClose").c_str(), mass, correction);
                            }
                            //end of chi2 test
                        }
                        //end of pair loop
                    }
                    //end of innermost vector pair loop
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

                    for(size_t pair_number = 0; pair_number<pairTabWithUnderscores.size(); pair_number++){
                        if(chi2Map[pairTabWithUnderscores[pair_number]]<9){
                            //getting all the indentificators out
                            std::pair<PARTICLES, PARTICLES> pairPID = pairTabPID[pair_number];
                            std::string temp_pair = pairTab[pair_number];

                            //calculating physics (general)
                            vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[pairPID.first]);
                            vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[pairPID.second]);
                            //track rotation (rotating only the negative one)
                            negative_track.RotateZ(TMath::Pi());
                            //the rest as usual
                            mass = (positive_track+negative_track).M();
                            eta = (positive_track+negative_track).Eta();
                            pT = (positive_track+negative_track).Pt();
                            correction = correction_coefficient(0, pairPID.first, pairPID.second, vector_Track_positive[i], vector_Track_negative[j], PV_position.Z(), PV_position.Z());
                            insideprocessing.Fill(("M"+temp_pair+"Chi2BcgTrackRotation").c_str(), mass, correction);
                            insideprocessing.Fill(("M"+temp_pair+"Chi2BcgTrackRotationeta").c_str(), mass, eta, correction);
                            insideprocessing.Fill(("M"+temp_pair+"Chi2BcgTrackRotationpT").c_str(), mass, pT, correction);

                            //calculating physics (specific)
                            if(temp_pair=="Kpi"||temp_pair=="piK"||temp_pair=="KK"){
                                insideprocessing.Fill(("M"+temp_pair+"Chi2BcgTrackRotationClose").c_str(), mass, correction);
                            }
                            //end of chi2 test
                        }
                        //end of pair loop
                    }
                    //end of innermost vector pair loop
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


                    for(size_t pair_number = 0; pair_number<pairTabWithUnderscores.size(); pair_number++){
                        if(chi2Map[pairTabWithUnderscores[pair_number]]<9){
                            //getting all the indentificators out
                            std::pair<PARTICLES, PARTICLES> pairPID = pairTabPID[pair_number];
                            std::string temp_pair = pairTab[pair_number];

                            //calculating physics (general)
                            vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[pairPID.first]);
                            vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[pairPID.second]);
                            //track rotation (rotating only the negative one)
                            negative_track.RotateZ(random_generator.Uniform(2*TMath::Pi()));
                            //the rest as usual
                            mass = (positive_track+negative_track).M();
                            eta = (positive_track+negative_track).Eta();
                            pT = (positive_track+negative_track).Pt();
                            correction = correction_coefficient(0, pairPID.first, pairPID.second, vector_Track_positive[i], vector_Track_negative[j], PV_position.Z(), PV_position.Z());
                            insideprocessing.Fill(("M"+temp_pair+"Chi2BcgRandomTrackRotation").c_str(), mass, correction);
                            insideprocessing.Fill(("M"+temp_pair+"Chi2BcgRandomTrackRotationeta").c_str(), mass, eta, correction);
                            insideprocessing.Fill(("M"+temp_pair+"Chi2BcgRandomTrackRotationpT").c_str(), mass, pT, correction);

                            //calculating physics (specific)
                            if(temp_pair=="Kpi"||temp_pair=="piK"||temp_pair=="KK"){
                                insideprocessing.Fill(("M"+temp_pair+"Chi2BcgRandomTrackRotationClose").c_str(), mass, correction);
                            }
                            //end of chi2 test
                        }
                        //end of pair loop
                    }
                    //end of innermost vector pair loop
                }
            }

            //########## BACKGROUND EXTRACTION (MIXED EVENT TECHNIQUES) ###############

            //TODO: have a mechanism to remove possible duplicates
            //as in, positive track matching to two negatives gets counted twice

            //picking pairs good for comparing with previous ones (and saving as ones)
            //we load the queue with the copy of the current event
            for(auto&& temp_pair_no_underscores:pairTab){
                map_of_queue_of_previous_vector_Tracks_positive[temp_pair_no_underscores].emplace_back();
                map_of_queue_of_previous_vector_Tracks_negative[temp_pair_no_underscores].emplace_back();
            }

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


                    for(size_t pair_number = 0; pair_number<pairTabWithUnderscores.size(); pair_number++){
                        if(chi2Map[pairTabWithUnderscores[pair_number]]<9){
                            //getting all the indentificators out
                            std::string temp_pair = pairTab[pair_number];

                            //positive
                            map_of_queue_of_previous_vector_Tracks_positive[temp_pair].back().push_back(new StUPCTrack());
                            vector_Track_positive[i]->getPtEtaPhi(pt, eta, phi);
                            map_of_queue_of_previous_vector_Tracks_positive[temp_pair].back().back()->setPtEtaPhi(pt, eta, phi);
                            map_of_queue_of_previous_vector_Tracks_positive[temp_pair].back().back()->setNhitsDEdx(vector_Track_positive[i]->getNhitsDEdx());
                            map_of_queue_of_previous_vector_Tracks_positive[temp_pair].back().back()->setTofPathLength(vector_Track_positive[i]->getTofPathLength());
                            map_of_queue_of_previous_vector_Tracks_positive[temp_pair].back().back()->setTofTime(vector_Track_positive[i]->getTofTime());
                            for(size_t part = 0; part<nParticlesExtended; part++){
                                map_of_queue_of_previous_vector_Tracks_positive[temp_pair].back().back()->setNSigmasTPC(static_cast<StUPCTrack::Part>(part), vector_Track_positive[i]->getNSigmasTPC(static_cast<StUPCTrack::Part>(part)));
                            }
                            //negative
                            map_of_queue_of_previous_vector_Tracks_negative[temp_pair].back().push_back(new StUPCTrack());
                            vector_Track_negative[j]->getPtEtaPhi(pt, eta, phi);
                            map_of_queue_of_previous_vector_Tracks_negative[temp_pair].back().back()->setPtEtaPhi(pt, eta, phi);
                            map_of_queue_of_previous_vector_Tracks_negative[temp_pair].back().back()->setNhitsDEdx(vector_Track_negative[j]->getNhitsDEdx());
                            map_of_queue_of_previous_vector_Tracks_negative[temp_pair].back().back()->setTofPathLength(vector_Track_negative[j]->getTofPathLength());
                            map_of_queue_of_previous_vector_Tracks_negative[temp_pair].back().back()->setTofTime(vector_Track_negative[j]->getTofTime());
                            for(size_t part = 0; part<nParticlesExtended; part++){
                                map_of_queue_of_previous_vector_Tracks_negative[temp_pair].back().back()->setNSigmasTPC(static_cast<StUPCTrack::Part>(part), vector_Track_negative[j]->getNSigmasTPC(static_cast<StUPCTrack::Part>(part)));
                            }

                            //end of filling current event
                        }
                        //end of pair loop
                    }
                    //end of innermost vector pair loop
                }
            }
            //noting current position of primary vertex and other values
            queue_of_previous_PV_positions.emplace_back(tempUPCpointer->getVertex(0)->getPosX(), tempUPCpointer->getVertex(0)->getPosY(), tempUPCpointer->getVertex(0)->getPosZ());
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
            //################################################################################
            //matching this event (on the back, "last" one, n-1) with previous ones (0 to n-2)
            //################################################################################
            //(to check it, we check the size of a random pair queue, as all of them have the same length)
            for(size_t evt = 0; evt<min(total_events_in_queue, (int)map_of_queue_of_previous_vector_Tracks_positive["Kpi"].size())-1; evt++){
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

                //################### MIXED-EVENT DIFFERENT-SIGN ##################

                for(size_t pair_number = 0; pair_number<pairTabWithUnderscores.size(); pair_number++){
                    //getting all the indentificators out
                    std::pair<PARTICLES, PARTICLES> pairPID = pairTabPID[pair_number];
                    std::string temp_pair = pairTab[pair_number];

                    //LOOPS
                    //this positive, past negative
                    for(long unsigned int i = 0; i<map_of_queue_of_previous_vector_Tracks_positive[temp_pair].back().size(); i++){
                        for(long unsigned int j = 0; j<map_of_queue_of_previous_vector_Tracks_negative[temp_pair][evt].size(); j++){
                            map_of_queue_of_previous_vector_Tracks_positive[temp_pair].back()[i]->getLorentzVector(positive_track, particleMass[pairPID.first]);
                            map_of_queue_of_previous_vector_Tracks_negative[temp_pair][evt][j]->getLorentzVector(negative_track, particleMass[pairPID.second]);
                            mass = (positive_track+negative_track).M();
                            eta = (positive_track+negative_track).Eta();
                            pT = (positive_track+negative_track).Pt();
                            correction = correction_coefficient(0, pairPID.first, pairPID.second, map_of_queue_of_previous_vector_Tracks_positive[temp_pair].back()[i], map_of_queue_of_previous_vector_Tracks_negative[temp_pair][evt][j], queue_of_previous_PV_positions.back().Z(), queue_of_previous_PV_positions[evt].Z());
                            //filling (general)
                            insideprocessing.Fill(("M"+temp_pair+"Chi2BcgMixedEvent").c_str(), mass, correction);
                            insideprocessing.Fill(("M"+temp_pair+"Chi2BcgMixedEventeta").c_str(), mass, eta, correction);
                            insideprocessing.Fill(("M"+temp_pair+"Chi2BcgMixedEventpT").c_str(), mass, pT, correction);
                            //filling (specific)
                            if(temp_pair=="Kpi"||temp_pair=="piK"||temp_pair=="KK"){
                                insideprocessing.Fill(("M"+temp_pair+"Chi2BcgMixedEventClose").c_str(), mass, correction);
                            }
                            //end of pairing
                        }
                    }
                    //this negative, past positive
                    for(long unsigned int i = 0; i<map_of_queue_of_previous_vector_Tracks_positive[temp_pair][evt].size(); i++){
                        for(long unsigned int j = 0; j<map_of_queue_of_previous_vector_Tracks_negative[temp_pair].back().size(); j++){
                            map_of_queue_of_previous_vector_Tracks_positive[temp_pair][evt][i]->getLorentzVector(positive_track, particleMass[pairPID.first]);
                            map_of_queue_of_previous_vector_Tracks_negative[temp_pair].back()[j]->getLorentzVector(negative_track, particleMass[pairPID.second]);
                            mass = (positive_track+negative_track).M();
                            eta = (positive_track+negative_track).Eta();
                            pT = (positive_track+negative_track).Pt();
                            correction = correction_coefficient(0, pairPID.first, pairPID.second, map_of_queue_of_previous_vector_Tracks_positive[temp_pair][evt][i], map_of_queue_of_previous_vector_Tracks_negative[temp_pair].back()[j], queue_of_previous_PV_positions[evt].Z(), queue_of_previous_PV_positions.back().Z());
                            //filling (general)
                            insideprocessing.Fill(("M"+temp_pair+"Chi2BcgMixedEvent").c_str(), mass, correction);
                            insideprocessing.Fill(("M"+temp_pair+"Chi2BcgMixedEventeta").c_str(), mass, eta, correction);
                            insideprocessing.Fill(("M"+temp_pair+"Chi2BcgMixedEventpT").c_str(), mass, pT, correction);
                            //filling (specific)
                            if(temp_pair=="Kpi"||temp_pair=="piK"||temp_pair=="KK"){
                                insideprocessing.Fill(("M"+temp_pair+"Chi2BcgMixedEventClose").c_str(), mass, correction);
                            }
                            //end of pairing
                        }
                    }
                    //end of particle loops
                }

                //####################### MIXED-EVENT SAME-SIGN ############################

                for(size_t pair_number = 0; pair_number<pairTabWithUnderscores.size(); pair_number++){
                    //getting all the indentificators out
                    std::pair<PARTICLES, PARTICLES> pairPID = pairTabPID[pair_number];
                    std::string temp_pair = pairTab[pair_number];
                    std::string reverse_temp_pair = reversePairAndRemoveUnderscore(pairTabWithUnderscores[pair_number]);

                    //LOOPS
                    //both positive
                    for(long unsigned int i = 0; i<map_of_queue_of_previous_vector_Tracks_positive[temp_pair][evt].size(); i++){
                        for(long unsigned int j = 0; j<map_of_queue_of_previous_vector_Tracks_positive[reverse_temp_pair].back().size(); j++){
                            map_of_queue_of_previous_vector_Tracks_positive[temp_pair][evt][i]->getLorentzVector(positive_track, particleMass[pairPID.first]);
                            map_of_queue_of_previous_vector_Tracks_positive[reverse_temp_pair].back()[j]->getLorentzVector(positive_track2, particleMass[pairPID.second]);
                            mass = (positive_track+positive_track2).M();
                            eta = (positive_track+positive_track2).Eta();
                            pT = (positive_track+positive_track2).Pt();
                            correction = correction_coefficient(2, pairPID.first, pairPID.second, map_of_queue_of_previous_vector_Tracks_positive[temp_pair][evt][i], map_of_queue_of_previous_vector_Tracks_positive[reverse_temp_pair].back()[j], queue_of_previous_PV_positions[evt].Z(), queue_of_previous_PV_positions.back().Z());
                            //filling (general)
                            insideprocessing.Fill(("M"+temp_pair+"Chi2BcgMixedEventSameSign").c_str(), mass, correction);
                            insideprocessing.Fill(("M"+temp_pair+"Chi2BcgMixedEventSameSigneta").c_str(), mass, eta, correction);
                            insideprocessing.Fill(("M"+temp_pair+"Chi2BcgMixedEventSameSignpT").c_str(), mass, pT, correction);
                            //filling (specific)
                            if(temp_pair=="Kpi"||temp_pair=="piK"||temp_pair=="KK"){
                                insideprocessing.Fill(("M"+temp_pair+"Chi2BcgMixedEventSameSignClose").c_str(), mass, correction);
                            }
                            //end of pairing
                        }
                    }
                    //both negative
                    for(long unsigned int i = 0; i<map_of_queue_of_previous_vector_Tracks_negative[temp_pair][evt].size(); i++){
                        for(long unsigned int j = 0; j<map_of_queue_of_previous_vector_Tracks_negative[reverse_temp_pair].back().size(); j++){
                            map_of_queue_of_previous_vector_Tracks_negative[temp_pair][evt][i]->getLorentzVector(negative_track, particleMass[pairPID.first]);
                            map_of_queue_of_previous_vector_Tracks_negative[reverse_temp_pair].back()[j]->getLorentzVector(negative_track2, particleMass[pairPID.second]);
                            mass = (negative_track+negative_track2).M();
                            eta = (negative_track+negative_track2).Eta();
                            pT = (negative_track+negative_track2).Pt();
                            correction = correction_coefficient(-2, pairPID.first, pairPID.second, map_of_queue_of_previous_vector_Tracks_negative[temp_pair][evt][i], map_of_queue_of_previous_vector_Tracks_negative[reverse_temp_pair].back()[j], queue_of_previous_PV_positions[evt].Z(), queue_of_previous_PV_positions.back().Z());
                            //filling (general)
                            insideprocessing.Fill(("M"+temp_pair+"Chi2BcgMixedEventSameSign").c_str(), mass, correction);
                            insideprocessing.Fill(("M"+temp_pair+"Chi2BcgMixedEventSameSigneta").c_str(), mass, eta, correction);
                            insideprocessing.Fill(("M"+temp_pair+"Chi2BcgMixedEventSameSignpT").c_str(), mass, pT, correction);
                            //filling (specific)
                            if(temp_pair=="Kpi"||temp_pair=="piK"||temp_pair=="KK"){
                                insideprocessing.Fill(("M"+temp_pair+"Chi2BcgMixedEventSameSignClose").c_str(), mass, correction);
                            }
                            //end of pairing
                        }
                    }
                    //if both particles are different, then we need to have crossover in both directions
                    //i.e. for K+pi+, we need:
                    //current K+   previous pi+
                    //current pi+  previous K+
                    //while for same pairs it would be:
                    //current K+   previous K+
                    //current K+  previous K+
                    //which obviously does not make sense
                    if(temp_pair==reverse_temp_pair){
                        continue;
                    }
                    //both positive
                    for(long unsigned int i = 0; i<map_of_queue_of_previous_vector_Tracks_positive[temp_pair].back().size(); i++){
                        for(long unsigned int j = 0; j<map_of_queue_of_previous_vector_Tracks_positive[reverse_temp_pair][evt].size(); j++){
                            map_of_queue_of_previous_vector_Tracks_positive[temp_pair].back()[i]->getLorentzVector(positive_track, particleMass[pairPID.first]);
                            map_of_queue_of_previous_vector_Tracks_positive[reverse_temp_pair][evt][j]->getLorentzVector(positive_track2, particleMass[pairPID.second]);
                            mass = (positive_track+positive_track2).M();
                            eta = (positive_track+positive_track2).Eta();
                            pT = (positive_track+positive_track2).Pt();
                            correction = correction_coefficient(2, pairPID.first, pairPID.second, map_of_queue_of_previous_vector_Tracks_positive[temp_pair].back()[i], map_of_queue_of_previous_vector_Tracks_positive[reverse_temp_pair][evt][j], queue_of_previous_PV_positions.back().Z(), queue_of_previous_PV_positions[evt].Z());
                            //filling (general)
                            insideprocessing.Fill(("M"+temp_pair+"Chi2BcgMixedEventSameSign").c_str(), mass, correction);
                            insideprocessing.Fill(("M"+temp_pair+"Chi2BcgMixedEventSameSigneta").c_str(), mass, eta, correction);
                            insideprocessing.Fill(("M"+temp_pair+"Chi2BcgMixedEventSameSignpT").c_str(), mass, pT, correction);
                            //filling (specific)
                            if(temp_pair=="Kpi"||temp_pair=="piK"||temp_pair=="KK"){
                                insideprocessing.Fill(("M"+temp_pair+"Chi2BcgMixedEventSameSignClose").c_str(), mass, correction);
                            }
                            //end of pairing
                        }
                    }
                    //both negative
                    for(long unsigned int i = 0; i<map_of_queue_of_previous_vector_Tracks_negative[temp_pair].back().size(); i++){
                        for(long unsigned int j = 0; j<map_of_queue_of_previous_vector_Tracks_negative[reverse_temp_pair][evt].size(); j++){
                            map_of_queue_of_previous_vector_Tracks_negative[temp_pair].back()[i]->getLorentzVector(negative_track, particleMass[pairPID.first]);
                            map_of_queue_of_previous_vector_Tracks_negative[reverse_temp_pair][evt][j]->getLorentzVector(negative_track2, particleMass[pairPID.second]);
                            mass = (negative_track+negative_track2).M();
                            eta = (negative_track+negative_track2).Eta();
                            pT = (negative_track+negative_track2).Pt();
                            correction = correction_coefficient(-2, pairPID.first, pairPID.second, map_of_queue_of_previous_vector_Tracks_negative[temp_pair].back()[i], map_of_queue_of_previous_vector_Tracks_negative[reverse_temp_pair][evt][j], queue_of_previous_PV_positions.back().Z(), queue_of_previous_PV_positions[evt].Z());
                            //filling (general)
                            insideprocessing.Fill(("M"+temp_pair+"Chi2BcgMixedEventSameSign").c_str(), mass, correction);
                            insideprocessing.Fill(("M"+temp_pair+"Chi2BcgMixedEventSameSigneta").c_str(), mass, eta, correction);
                            insideprocessing.Fill(("M"+temp_pair+"Chi2BcgMixedEventSameSignpT").c_str(), mass, pT, correction);
                            //filling (specific)
                            if(temp_pair=="Kpi"||temp_pair=="piK"||temp_pair=="KK"){
                                insideprocessing.Fill(("M"+temp_pair+"Chi2BcgMixedEventSameSignClose").c_str(), mass, correction);
                            }
                            //end of pairing
                        }
                    }
                    //end of THE LOOPS
                }
                //end of previous event loop
            }
            //if the recall limit has been reached, we pop the oldest one
            //(to check if its happening, we check the size of a random pair queue, as all of them have the same length)
            if(map_of_queue_of_previous_vector_Tracks_positive["Kpi"].size()>total_events_in_queue){
                for(auto&& temp_pair_no_underscores:pairTab){
                    for(size_t i = 0; i<map_of_queue_of_previous_vector_Tracks_positive[temp_pair_no_underscores].front().size(); i++){
                        delete map_of_queue_of_previous_vector_Tracks_positive[temp_pair_no_underscores].front()[i];
                    }
                    map_of_queue_of_previous_vector_Tracks_positive[temp_pair_no_underscores].pop_front();
                    for(size_t i = 0; i<map_of_queue_of_previous_vector_Tracks_negative[temp_pair_no_underscores].front().size(); i++){
                        delete map_of_queue_of_previous_vector_Tracks_negative[temp_pair_no_underscores].front()[i];
                    }
                    map_of_queue_of_previous_vector_Tracks_negative[temp_pair_no_underscores].pop_front();
                }
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
    outputFileHist->cd();
    efficiency1.Write();
    efficiency2.Write();
    efficiency3.Write();
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

std::pair<PARTICLES, PARTICLES> getPairPID(std::string pair_with_underscore){
    int underscore_position = pair_with_underscore.find("_");
    std::string firstPID = pair_with_underscore.substr(0, underscore_position);
    std::string secondPID = pair_with_underscore.substr(underscore_position+1, pair_with_underscore.size()-1-underscore_position);
    //example:
    //0123
    //K_pi
    //underscore_position: 1
    //first string: from 0, length 1
    //second string: from 2 (1+1), length 2 (4-1-1)
    std::pair<PARTICLES, PARTICLES> output;
    if(firstPID=="pi"){
        output.first = Pion;
    } else if(firstPID=="K"){
        output.first = Kaon;
    } else if(firstPID=="p"){
        output.first = Proton;
    } else{
        output.first = nParticles;
    }
    if(secondPID=="pi"){
        output.second = Pion;
    } else if(secondPID=="K"){
        output.second = Kaon;
    } else if(secondPID=="p"){
        output.second = Proton;
    } else{
        output.second = nParticles;
    }
    return output;
}

std::string reversePairAndRemoveUnderscore(std::string pair_with_underscore){
    int underscore_position = pair_with_underscore.find("_");
    std::string firstPID = pair_with_underscore.substr(0, underscore_position);
    std::string secondPID = pair_with_underscore.substr(underscore_position+1, pair_with_underscore.size()-1-underscore_position);
    return secondPID+firstPID;
}