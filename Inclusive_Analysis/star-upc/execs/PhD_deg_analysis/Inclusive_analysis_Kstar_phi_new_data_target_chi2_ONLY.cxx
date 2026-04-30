//cpp headers
#include <map>
#include <fstream>
#include <deque>

//ROOT headers
#include <ROOT/TThreadedObject.hxx>
#include <TTreeReader.h>
#include <ROOT/TTreeProcessorMT.hxx>
#include <TRandom3.h>

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
#include "UsefulThings.h"
#include "ProcessingInsideLoop.h"
#include "ProcessingOutsideLoop.h"

enum{
    kAll = 1, kCPT, kRP, kOneVertex, kTPCTOF,
    kTotQ, kMax
};
enum SIDE{ E = 0, East = 0, W = 1, West = 1, nSides };
enum PARTICLES{ Pion = 0, Kaon = 1, Proton = 2, nParticles };
enum EXTENDED_PARTICLES{ ExtElectron = 0, ExtPion = 1, ExtKaon = 2, ExtProton = 3, nParticlesExtended };
const double particleMass[nParticles] = { 0.13957, 0.493677, 0.93827 }; // pion, kaon, proton in GeV /c^2 
const double particleMassExtended[nParticlesExtended] = { 0.000510999, 0.13957, 0.493677, 0.93827 }; // electron, pion, kaon, proton in GeV /c^2 
enum BRANCH_ID{ EU, ED, WU, WD, nBranches };
enum RP_ID{ E1U, E1D, E2U, E2D, W1U, W1D, W2U, W2D, nRomanPots };
enum SUSPECTED_PARTICLES{ K0S, Lambda, Kstar, Phi };
string particleNicks[nParticlesExtended] = { "e", "pi", "K", "p" };

double getChi2(StUPCTrack* positive, StUPCTrack* negative, int positiveId, int negativeId, double sigmaT);

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
    outsideprocessing.AddHistogram(TH1D("pairInfoSignal", "", 1, 0, 1));
    outsideprocessing.AddHistogram(TH1D("pairInfoBackgroundSameSign", "", 1, 0, 1));
    outsideprocessing.AddHistogram(TH1D("pairInfoBackgroundTrackRotation", "", 1, 0, 1));
    outsideprocessing.AddHistogram(TH1D("pairInfoBackgroundRandomTrackRotation", "", 1, 0, 1));

    outsideprocessing.AddHistogram(TH1D("Flowchart", "", 1, 0, 1));

    //adding chi2 histograms
    //file with sigma values:
    ifstream sigmaFile;
    string line;
    char sigmaName[10];
    double sigmaValue;
    map<string, double> sigmaMap;
    sigmaFile.open("STAR-Analysis/Inclusive_Analysis/star-upc/AnaOutput_Inclusive_analysis_Kstar_phi_old_data_TOF_tests_sigmaValues.txt");
    while(getline(sigmaFile, line)){
        sscanf(line.c_str(), "%s\t\t%lf", sigmaName, &sigmaValue);
        sigmaMap.insert({ string(sigmaName), sigmaValue });
    }
    sigmaFile.close();

    std::vector<std::string> pairTab = { "Kpi", "piK", "ppi", "pip", "KK", "pipi", "pp" };

    //mass histograms (signal)
    outsideprocessing.AddHistogram(TH1D("MKpiChi2", ";m_{K^{+}#pi^{-}} [GeV];Number of pairs", 200, 0.5, 2.0));
    outsideprocessing.AddHistogram(TH1D("MpiKChi2", ";m_{#pi^{+}K^{-}} [GeV];Number of pairs", 200, 0.5, 2.0));
    outsideprocessing.AddHistogram(TH1D("MppiChi2", ";m_{p^{+}#pi^{-}} [GeV];Number of pairs", 500, 1.0, 2.5));
    outsideprocessing.AddHistogram(TH1D("MpipChi2", ";m_{#pi^{+}p^{-}} [GeV];Number of pairs", 500, 1.0, 2.5));
    outsideprocessing.AddHistogram(TH1D("MKKChi2", ";m_{K^{+}K^{-}} [GeV];Number of pairs", 500, 0.9, 2.4));
    outsideprocessing.AddHistogram(TH1D("MpipiChi2", ";m_{#pi^{+}#pi^{-}} [GeV];Number of pairs", 600, 0.2, 1.4));
    outsideprocessing.AddHistogram(TH1D("MppChi2", ";m_{p^{+}p^{-}} [GeV];Number of pairs", 500, 1.5, 3.5));
    //closer histograms (signal)
    outsideprocessing.AddHistogram(TH1D("MKKChi2Close", ";m_{K^{+}K^{-}} [GeV];Number of pairs", 50, 0.99, 1.05));
    outsideprocessing.AddHistogram(TH1D("MKpiChi2Close", ";m_{K^{+}#pi^{-}} [GeV];Number of pairs", 50, 0.7, 1.1));
    outsideprocessing.AddHistogram(TH1D("MpiKChi2Close", ";m_{#pi^{+}K^{-}} [GeV];Number of pairs", 50, 0.7, 1.1));
    //adding mass histograms grouped by category (signal)
    getCategoryHistograms(outsideprocessing, pairTab);

    //mass histograms (background, same sign)
    outsideprocessing.AddHistogram(TH1D("MKpiChi2BcgSameSign", ";m_{K^{+}#pi^{-}} [GeV];Number of pairs", 200, 0.5, 2.0));
    outsideprocessing.AddHistogram(TH1D("MpiKChi2BcgSameSign", ";m_{#pi^{+}K^{-}} [GeV];Number of pairs", 200, 0.5, 2.0));
    outsideprocessing.AddHistogram(TH1D("MppiChi2BcgSameSign", ";m_{p^{+}#pi^{-}} [GeV];Number of pairs", 500, 1.0, 2.5));
    outsideprocessing.AddHistogram(TH1D("MpipChi2BcgSameSign", ";m_{#pi^{+}p^{-}} [GeV];Number of pairs", 500, 1.0, 2.5));
    outsideprocessing.AddHistogram(TH1D("MKKChi2BcgSameSign", ";m_{K^{+}K^{-}} [GeV];Number of pairs", 500, 0.9, 2.4));
    outsideprocessing.AddHistogram(TH1D("MpipiChi2BcgSameSign", ";m_{#pi^{+}#pi^{-}} [GeV];Number of pairs", 600, 0.2, 1.4));
    outsideprocessing.AddHistogram(TH1D("MppChi2BcgSameSign", ";m_{p^{+}p^{-}} [GeV];Number of pairs", 500, 1.5, 3.5));
    //closer histograms (background, same sign)
    outsideprocessing.AddHistogram(TH1D("MKKChi2BcgSameSignClose", ";m_{K^{+}K^{-}} [GeV];Number of pairs", 50, 0.99, 1.05));
    outsideprocessing.AddHistogram(TH1D("MKpiChi2BcgSameSignClose", ";m_{K^{+}#pi^{-}} [GeV];Number of pairs", 50, 0.7, 1.1));
    outsideprocessing.AddHistogram(TH1D("MpiKChi2BcgSameSignClose", ";m_{#pi^{+}K^{-}} [GeV];Number of pairs", 50, 0.7, 1.1));
    //adding mass histograms grouped by category (background, same sign)
    getCategoryHistograms(outsideprocessing, pairTab, "BcgSameSign");

    //mass histograms (background, track rotation)
    outsideprocessing.AddHistogram(TH1D("MKpiChi2BcgTrackRotation", ";m_{K^{+}#pi^{-}} [GeV];Number of pairs", 200, 0.5, 2.0));
    outsideprocessing.AddHistogram(TH1D("MpiKChi2BcgTrackRotation", ";m_{#pi^{+}K^{-}} [GeV];Number of pairs", 200, 0.5, 2.0));
    outsideprocessing.AddHistogram(TH1D("MppiChi2BcgTrackRotation", ";m_{p^{+}#pi^{-}} [GeV];Number of pairs", 500, 1.0, 2.5));
    outsideprocessing.AddHistogram(TH1D("MpipChi2BcgTrackRotation", ";m_{#pi^{+}p^{-}} [GeV];Number of pairs", 500, 1.0, 2.5));
    outsideprocessing.AddHistogram(TH1D("MKKChi2BcgTrackRotation", ";m_{K^{+}K^{-}} [GeV];Number of pairs", 500, 0.9, 2.4));
    outsideprocessing.AddHistogram(TH1D("MpipiChi2BcgTrackRotation", ";m_{#pi^{+}#pi^{-}} [GeV];Number of pairs", 600, 0.2, 1.4));
    outsideprocessing.AddHistogram(TH1D("MppChi2BcgTrackRotation", ";m_{p^{+}p^{-}} [GeV];Number of pairs", 500, 1.5, 3.5));
    //closer histograms (background, track rotation)
    outsideprocessing.AddHistogram(TH1D("MKKChi2BcgTrackRotationClose", ";m_{K^{+}K^{-}} [GeV];Number of pairs", 50, 0.99, 1.05));
    outsideprocessing.AddHistogram(TH1D("MKpiChi2BcgTrackRotationClose", ";m_{K^{+}#pi^{-}} [GeV];Number of pairs", 50, 0.7, 1.1));
    outsideprocessing.AddHistogram(TH1D("MpiKChi2BcgTrackRotationClose", ";m_{#pi^{+}K^{-}} [GeV];Number of pairs", 50, 0.7, 1.1));
    //adding mass histograms grouped by category (background, track rotation)
    getCategoryHistograms(outsideprocessing, pairTab, "BcgTrackRotation");

    //mass histograms (background, random track rotation)
    outsideprocessing.AddHistogram(TH1D("MKpiChi2BcgRandomTrackRotation", ";m_{K^{+}#pi^{-}} [GeV];Number of pairs", 200, 0.5, 2.0));
    outsideprocessing.AddHistogram(TH1D("MpiKChi2BcgRandomTrackRotation", ";m_{#pi^{+}K^{-}} [GeV];Number of pairs", 200, 0.5, 2.0));
    outsideprocessing.AddHistogram(TH1D("MppiChi2BcgRandomTrackRotation", ";m_{p^{+}#pi^{-}} [GeV];Number of pairs", 500, 1.0, 2.5));
    outsideprocessing.AddHistogram(TH1D("MpipChi2BcgRandomTrackRotation", ";m_{#pi^{+}p^{-}} [GeV];Number of pairs", 500, 1.0, 2.5));
    outsideprocessing.AddHistogram(TH1D("MKKChi2BcgRandomTrackRotation", ";m_{K^{+}K^{-}} [GeV];Number of pairs", 500, 0.9, 2.4));
    outsideprocessing.AddHistogram(TH1D("MpipiChi2BcgRandomTrackRotation", ";m_{#pi^{+}#pi^{-}} [GeV];Number of pairs", 600, 0.2, 1.4));
    outsideprocessing.AddHistogram(TH1D("MppChi2BcgRandomTrackRotation", ";m_{p^{+}p^{-}} [GeV];Number of pairs", 500, 1.5, 3.5));
    //closer histograms (background, random track rotation)
    outsideprocessing.AddHistogram(TH1D("MKKChi2BcgRandomTrackRotationClose", ";m_{K^{+}K^{-}} [GeV];Number of pairs", 50, 0.99, 1.05));
    outsideprocessing.AddHistogram(TH1D("MKpiChi2BcgRandomTrackRotationClose", ";m_{K^{+}#pi^{-}} [GeV];Number of pairs", 50, 0.7, 1.1));
    outsideprocessing.AddHistogram(TH1D("MpiKChi2BcgRandomTrackRotationClose", ";m_{#pi^{+}K^{-}} [GeV];Number of pairs", 50, 0.7, 1.1));
    //adding mass histograms grouped by category (background, random track rotation)
    getCategoryHistograms(outsideprocessing, pairTab, "BcgRandomTrackRotation");

    //mass histograms (background, mixed events)
    outsideprocessing.AddHistogram(TH1D("MKpiChi2BcgMixedEvent", ";m_{K^{+}#pi^{-}} [GeV];Number of pairs", 200, 0.5, 2.0));
    outsideprocessing.AddHistogram(TH1D("MpiKChi2BcgMixedEvent", ";m_{#pi^{+}K^{-}} [GeV];Number of pairs", 200, 0.5, 2.0));
    outsideprocessing.AddHistogram(TH1D("MppiChi2BcgMixedEvent", ";m_{p^{+}#pi^{-}} [GeV];Number of pairs", 500, 1.0, 2.5));
    outsideprocessing.AddHistogram(TH1D("MpipChi2BcgMixedEvent", ";m_{#pi^{+}p^{-}} [GeV];Number of pairs", 500, 1.0, 2.5));
    outsideprocessing.AddHistogram(TH1D("MKKChi2BcgMixedEvent", ";m_{K^{+}K^{-}} [GeV];Number of pairs", 500, 0.9, 2.4));
    outsideprocessing.AddHistogram(TH1D("MpipiChi2BcgMixedEvent", ";m_{#pi^{+}#pi^{-}} [GeV];Number of pairs", 600, 0.2, 1.4));
    outsideprocessing.AddHistogram(TH1D("MppChi2BcgMixedEvent", ";m_{p^{+}p^{-}} [GeV];Number of pairs", 500, 1.5, 3.5));
    //closer histograms (background, mixed events)
    outsideprocessing.AddHistogram(TH1D("MKKChi2BcgMixedEventClose", ";m_{K^{+}K^{-}} [GeV];Number of pairs", 50, 0.99, 1.05));
    outsideprocessing.AddHistogram(TH1D("MKpiChi2BcgMixedEventClose", ";m_{K^{+}#pi^{-}} [GeV];Number of pairs", 50, 0.7, 1.1));
    outsideprocessing.AddHistogram(TH1D("MpiKChi2BcgMixedEventClose", ";m_{#pi^{+}K^{-}} [GeV];Number of pairs", 50, 0.7, 1.1));
    //adding mass histograms grouped by category (background, mixed events)
    getCategoryHistograms(outsideprocessing, pairTab, "BcgMixedEvent");


    //other histograms
    outsideprocessing.AddHistogram(TH1D("MKKSuspiciousPeakTestedAsPionPair", ";m_{#pi^{+}#pi^{-}} [GeV];Number of pairs", 400, 0.25, 0.65));
    outsideprocessing.AddHistogram(TH1D("MKKSuspiciousPeakTestedAsPionPairNeighbourhood", ";m_{#pi^{+}#pi^{-}} [GeV];Number of pairs", 400, 0.25, 0.65));
    outsideprocessing.AddHistogram(TH1D("MKKSuspiciousPeakTestedWithStrictChi2LessThan3", ";m_{K^{+}K^{-}} [GeV];Number of pairs", 500, 0.9, 2.4));
    outsideprocessing.AddHistogram(TH1D("MKKSuspiciousPeakTestedWithStrictChi2LessThan1", ";m_{K^{+}K^{-}} [GeV];Number of pairs", 500, 0.9, 2.4));

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
        // TTreeReaderValue<StRPEvent> StRPEventInstance(myReader, "mRPEvent");
        ProcessingInsideLoop insideprocessing;
        StUPCEvent* tempUPCpointer;
        // StRPEvent* tempRPpointer;
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
        double mass, chi2pipi, chi2Kpi, eta, pT, phi;
        map<string, double> chi2Map;
        bool isdEdxOk, isTOFOk;
        string tempPairName;
        TRandom3 random_generator;

        //queues for vectors of tracks from previous events assigned to particular pairs
        // std::deque<std::vector<StUPCTrack*>> queue_of_previous_vector_Tracks_Kpi_positive;
        // std::deque<std::vector<StUPCTrack*>> queue_of_previous_vector_Tracks_piK_positive;
        // std::deque<std::vector<StUPCTrack*>> queue_of_previous_vector_Tracks_ppi_positive;
        // std::deque<std::vector<StUPCTrack*>> queue_of_previous_vector_Tracks_pip_positive;
        // std::deque<std::vector<StUPCTrack*>> queue_of_previous_vector_Tracks_KK_positive;
        // std::deque<std::vector<StUPCTrack*>> queue_of_previous_vector_Tracks_pipi_positive;
        // std::deque<std::vector<StUPCTrack*>> queue_of_previous_vector_Tracks_pp_positive;
        // std::deque<std::vector<StUPCTrack*>> queue_of_previous_vector_Tracks_Kpi_negative;
        // std::deque<std::vector<StUPCTrack*>> queue_of_previous_vector_Tracks_piK_negative;
        // std::deque<std::vector<StUPCTrack*>> queue_of_previous_vector_Tracks_ppi_negative;
        // std::deque<std::vector<StUPCTrack*>> queue_of_previous_vector_Tracks_pip_negative;
        // std::deque<std::vector<StUPCTrack*>> queue_of_previous_vector_Tracks_KK_negative;
        // std::deque<std::vector<StUPCTrack*>> queue_of_previous_vector_Tracks_pipi_negative;
        // std::deque<std::vector<StUPCTrack*>> queue_of_previous_vector_Tracks_pp_negative;

        //TODO
        //Teraz mamy wszystkie pary cząstek i nie mamy dobrego tła mixed events, jak ograniczymy pary to dobrego tła mieć nie będziemy na 100%
        //funkcja na dopasowywanie rejonu dopasowywania tła

        std::deque<std::vector<StUPCTrack*>> queue_of_previous_vector_Tracks_positive;
        std::deque<std::vector<StUPCTrack*>> queue_of_previous_vector_Tracks_negative;

        //parameters
        const int previous_events_in_queue = 1000;

        //actual loop
        while(myReader.Next()){
            //in a TTree, it *would* be constant, in TChain however not necessarily
            tempUPCpointer = StUPCEventInstance.Get();
            // tempRPpointer = StRPEventInstance.Get();

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
                if(tempTrack->getCharge()>0){
                    vector_Track_positive.push_back(tempTrack);
                } else{
                    vector_Track_negative.push_back(tempTrack);
                }
                nOfGoodTracks++;
            }

            //filling a chi2 map with keys for all the possibilities
            for(auto const& imap:sigmaMap){
                chi2Map.insert({ imap.first, 0. });
            }

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
                        insideprocessing.Fill("MKpiChi2", mass);
                        insideprocessing.Fill("MKpiChi2Close", mass);
                        insideprocessing.Fill("MKpiChi2eta", mass, eta);
                        insideprocessing.Fill("MKpiChi2pT", mass, pT);
                    }
                    if(chi2Map["pi_K"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                        vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Kaon]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        insideprocessing.Fill("MpiKChi2", mass);
                        insideprocessing.Fill("MpiKChi2Close", mass);
                        insideprocessing.Fill("MpiKChi2eta", mass, eta);
                        insideprocessing.Fill("MpiKChi2pT", mass, pT);
                    }
                    if(chi2Map["p_pi"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Proton]);
                        vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Pion]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        insideprocessing.Fill("MppiChi2", mass);
                        insideprocessing.Fill("MppiChi2eta", mass, eta);
                        insideprocessing.Fill("MppiChi2pT", mass, pT);
                    }
                    if(chi2Map["pi_p"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                        vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Proton]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        insideprocessing.Fill("MpipChi2", mass);
                        insideprocessing.Fill("MpipChi2eta", mass, eta);
                        insideprocessing.Fill("MpipChi2pT", mass, pT);
                    }
                    if(chi2Map["K_K"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Kaon]);
                        vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Kaon]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        insideprocessing.Fill("MKKChi2", mass);
                        insideprocessing.Fill("MKKChi2Close", mass);
                        insideprocessing.Fill("MKKChi2eta", mass, eta);
                        insideprocessing.Fill("MKKChi2pT", mass, pT);
                        //test of suspicious peak and its neighbourhood
                        if(chi2Map["K_K"]<3){
                            insideprocessing.Fill("MKKSuspiciousPeakTestedWithStrictChi2LessThan3", mass);
                        }
                        if(chi2Map["K_K"]<1){
                            insideprocessing.Fill("MKKSuspiciousPeakTestedWithStrictChi2LessThan1", mass);
                        }
                        if(1.06<mass&&mass<1.08){
                            vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                            vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Pion]);
                            mass = (positive_track+negative_track).M();
                            insideprocessing.Fill("MKKSuspiciousPeakTestedAsPionPair", mass);
                        }
                        if((1.05<mass&&mass<1.06)||(1.08<mass&&mass<1.09)){
                            vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                            vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Pion]);
                            mass = (positive_track+negative_track).M();
                            insideprocessing.Fill("MKKSuspiciousPeakTestedAsPionPairNeighbourhood", mass);
                        }
                    }
                    if(chi2Map["pi_pi"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                        vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Pion]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
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
                        insideprocessing.Fill("MKpiChi2BcgSameSign", mass);
                        insideprocessing.Fill("MKpiChi2BcgSameSignClose", mass);
                        insideprocessing.Fill("MKpiChi2BcgSameSigneta", mass, eta);
                        insideprocessing.Fill("MKpiChi2BcgSameSignpT", mass, pT);
                    }
                    if(chi2Map["pi_K"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                        vector_Track_positive[j]->getLorentzVector(positive_track2, particleMass[Kaon]);
                        mass = (positive_track+positive_track2).M();
                        eta = (positive_track+positive_track2).Eta();
                        pT = (positive_track+positive_track2).Pt();
                        insideprocessing.Fill("MpiKChi2BcgSameSign", mass);
                        insideprocessing.Fill("MpiKChi2BcgSameSignClose", mass);
                        insideprocessing.Fill("MpiKChi2BcgSameSigneta", mass, eta);
                        insideprocessing.Fill("MpiKChi2BcgSameSignpT", mass, pT);
                    }
                    if(chi2Map["p_pi"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Proton]);
                        vector_Track_positive[j]->getLorentzVector(positive_track2, particleMass[Pion]);
                        mass = (positive_track+positive_track2).M();
                        eta = (positive_track+positive_track2).Eta();
                        pT = (positive_track+positive_track2).Pt();
                        insideprocessing.Fill("MppiChi2BcgSameSign", mass);
                        insideprocessing.Fill("MppiChi2BcgSameSigneta", mass, eta);
                        insideprocessing.Fill("MppiChi2BcgSameSignpT", mass, pT);
                    }
                    if(chi2Map["pi_p"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                        vector_Track_positive[j]->getLorentzVector(positive_track2, particleMass[Proton]);
                        mass = (positive_track+positive_track2).M();
                        eta = (positive_track+positive_track2).Eta();
                        pT = (positive_track+positive_track2).Pt();
                        insideprocessing.Fill("MpipChi2BcgSameSign", mass);
                        insideprocessing.Fill("MpipChi2BcgSameSigneta", mass, eta);
                        insideprocessing.Fill("MpipChi2BcgSameSignpT", mass, pT);
                    }
                    if(chi2Map["K_K"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Kaon]);
                        vector_Track_positive[j]->getLorentzVector(positive_track2, particleMass[Kaon]);
                        mass = (positive_track+positive_track2).M();
                        eta = (positive_track+positive_track2).Eta();
                        pT = (positive_track+positive_track2).Pt();
                        insideprocessing.Fill("MKKChi2BcgSameSign", mass);
                        insideprocessing.Fill("MKKChi2BcgSameSignClose", mass);
                        insideprocessing.Fill("MKKChi2BcgSameSigneta", mass, eta);
                        insideprocessing.Fill("MKKChi2BcgSameSignpT", mass, pT);
                    }
                    if(chi2Map["pi_pi"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                        vector_Track_positive[j]->getLorentzVector(positive_track2, particleMass[Pion]);
                        mass = (positive_track+positive_track2).M();
                        eta = (positive_track+positive_track2).Eta();
                        pT = (positive_track+positive_track2).Pt();
                        insideprocessing.Fill("MpipiChi2BcgSameSign", mass);
                        insideprocessing.Fill("MpipiChi2BcgSameSigneta", mass, eta);
                        insideprocessing.Fill("MpipiChi2BcgSameSignpT", mass, pT);
                    }
                    if(chi2Map["p_p"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Proton]);
                        vector_Track_positive[j]->getLorentzVector(positive_track2, particleMass[Proton]);
                        mass = (positive_track+positive_track2).M();
                        eta = (positive_track+positive_track2).Eta();
                        pT = (positive_track+positive_track2).Pt();
                        insideprocessing.Fill("MppChi2BcgSameSign", mass);
                        insideprocessing.Fill("MppChi2BcgSameSigneta", mass, eta);
                        insideprocessing.Fill("MppChi2BcgSameSignpT", mass, pT);
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
                        insideprocessing.Fill("MKpiChi2BcgSameSigneta", mass, eta);
                        insideprocessing.Fill("MKpiChi2BcgSameSign", mass);
                        insideprocessing.Fill("MKpiChi2BcgSameSignClose", mass);
                        insideprocessing.Fill("MKpiChi2BcgSameSignpT", mass, pT);
                    }
                    if(chi2Map["pi_K"]<9){
                        vector_Track_negative[i]->getLorentzVector(negative_track, particleMass[Pion]);
                        vector_Track_negative[j]->getLorentzVector(negative_track2, particleMass[Kaon]);
                        mass = (negative_track+negative_track2).M();
                        eta = (negative_track+negative_track2).Eta();
                        pT = (negative_track+negative_track2).Pt();
                        insideprocessing.Fill("MpiKChi2BcgSameSign", mass);
                        insideprocessing.Fill("MpiKChi2BcgSameSignClose", mass);
                        insideprocessing.Fill("MpiKChi2BcgSameSigneta", mass, eta);
                        insideprocessing.Fill("MpiKChi2BcgSameSignpT", mass, pT);
                    }
                    if(chi2Map["p_pi"]<9){
                        vector_Track_negative[i]->getLorentzVector(negative_track, particleMass[Proton]);
                        vector_Track_negative[j]->getLorentzVector(negative_track2, particleMass[Pion]);
                        mass = (negative_track+negative_track2).M();
                        eta = (negative_track+negative_track2).Eta();
                        pT = (negative_track+negative_track2).Pt();
                        insideprocessing.Fill("MppiChi2BcgSameSign", mass);
                        insideprocessing.Fill("MppiChi2BcgSameSigneta", mass, eta);
                        insideprocessing.Fill("MppiChi2BcgSameSignpT", mass, pT);
                    }
                    if(chi2Map["pi_p"]<9){
                        vector_Track_negative[i]->getLorentzVector(negative_track, particleMass[Pion]);
                        vector_Track_negative[j]->getLorentzVector(negative_track2, particleMass[Proton]);
                        mass = (negative_track+negative_track2).M();
                        eta = (negative_track+negative_track2).Eta();
                        pT = (negative_track+negative_track2).Pt();
                        insideprocessing.Fill("MpipChi2BcgSameSign", mass);
                        insideprocessing.Fill("MpipChi2BcgSameSigneta", mass, eta);
                        insideprocessing.Fill("MpipChi2BcgSameSignpT", mass, pT);
                    }
                    if(chi2Map["K_K"]<9){
                        vector_Track_negative[i]->getLorentzVector(negative_track, particleMass[Kaon]);
                        vector_Track_negative[j]->getLorentzVector(negative_track2, particleMass[Kaon]);
                        mass = (negative_track+negative_track2).M();
                        eta = (negative_track+negative_track2).Eta();
                        pT = (negative_track+negative_track2).Pt();
                        insideprocessing.Fill("MKKChi2BcgSameSign", mass);
                        insideprocessing.Fill("MKKChi2BcgSameSignClose", mass);
                        insideprocessing.Fill("MKKChi2BcgSameSigneta", mass, eta);
                        insideprocessing.Fill("MKKChi2BcgSameSignpT", mass, pT);
                    }
                    if(chi2Map["pi_pi"]<9){
                        vector_Track_negative[i]->getLorentzVector(negative_track, particleMass[Pion]);
                        vector_Track_negative[j]->getLorentzVector(negative_track2, particleMass[Pion]);
                        mass = (negative_track+negative_track2).M();
                        eta = (negative_track+negative_track2).Eta();
                        pT = (negative_track+negative_track2).Pt();
                        insideprocessing.Fill("MpipiChi2BcgSameSign", mass);
                        insideprocessing.Fill("MpipiChi2BcgSameSigneta", mass, eta);
                        insideprocessing.Fill("MpipiChi2BcgSameSignpT", mass, pT);
                    }
                    if(chi2Map["p_p"]<9){
                        vector_Track_negative[i]->getLorentzVector(negative_track, particleMass[Proton]);
                        vector_Track_negative[j]->getLorentzVector(negative_track2, particleMass[Proton]);
                        mass = (negative_track+negative_track2).M();
                        eta = (negative_track+negative_track2).Eta();
                        pT = (negative_track+negative_track2).Pt();
                        insideprocessing.Fill("MppChi2BcgSameSign", mass);
                        insideprocessing.Fill("MppChi2BcgSameSigneta", mass, eta);
                        insideprocessing.Fill("MppChi2BcgSameSignpT", mass, pT);
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
                        insideprocessing.Fill("MKpiChi2BcgTrackRotation", mass);
                        insideprocessing.Fill("MKpiChi2BcgTrackRotationClose", mass);
                        insideprocessing.Fill("MKpiChi2BcgTrackRotationeta", mass, eta);
                        insideprocessing.Fill("MKpiChi2BcgTrackRotationpT", mass, pT);
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
                        insideprocessing.Fill("MpiKChi2BcgTrackRotation", mass);
                        insideprocessing.Fill("MpiKChi2BcgTrackRotationClose", mass);
                        insideprocessing.Fill("MpiKChi2BcgTrackRotationeta", mass, eta);
                        insideprocessing.Fill("MpiKChi2BcgTrackRotationpT", mass, pT);
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
                        insideprocessing.Fill("MppiChi2BcgTrackRotation", mass);
                        insideprocessing.Fill("MppiChi2BcgTrackRotationeta", mass, eta);
                        insideprocessing.Fill("MppiChi2BcgTrackRotationpT", mass, pT);
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
                        insideprocessing.Fill("MpipChi2BcgTrackRotation", mass);
                        insideprocessing.Fill("MpipChi2BcgTrackRotationeta", mass, eta);
                        insideprocessing.Fill("MpipChi2BcgTrackRotationpT", mass, pT);
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
                        insideprocessing.Fill("MKKChi2BcgTrackRotation", mass);
                        insideprocessing.Fill("MKKChi2BcgTrackRotationClose", mass);
                        insideprocessing.Fill("MKKChi2BcgTrackRotationeta", mass, eta);
                        insideprocessing.Fill("MKKChi2BcgTrackRotationpT", mass, pT);
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
                        insideprocessing.Fill("MpipiChi2BcgTrackRotation", mass);
                        insideprocessing.Fill("MpipiChi2BcgTrackRotationeta", mass, eta);
                        insideprocessing.Fill("MpipiChi2BcgTrackRotationpT", mass, pT);
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
                        insideprocessing.Fill("MppChi2BcgTrackRotation", mass);
                        insideprocessing.Fill("MppChi2BcgTrackRotationeta", mass, eta);
                        insideprocessing.Fill("MppChi2BcgTrackRotationpT", mass, pT);
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
                        insideprocessing.Fill("MKpiChi2BcgRandomTrackRotation", mass);
                        insideprocessing.Fill("MKpiChi2BcgRandomTrackRotationClose", mass);
                        insideprocessing.Fill("MKpiChi2BcgRandomTrackRotationeta", mass, eta);
                        insideprocessing.Fill("MKpiChi2BcgRandomTrackRotationpT", mass, pT);
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
                        insideprocessing.Fill("MpiKChi2BcgRandomTrackRotation", mass);
                        insideprocessing.Fill("MpiKChi2BcgRandomTrackRotationClose", mass);
                        insideprocessing.Fill("MpiKChi2BcgRandomTrackRotationeta", mass, eta);
                        insideprocessing.Fill("MpiKChi2BcgRandomTrackRotationpT", mass, pT);
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
                        insideprocessing.Fill("MppiChi2BcgRandomTrackRotation", mass);
                        insideprocessing.Fill("MppiChi2BcgRandomTrackRotationeta", mass, eta);
                        insideprocessing.Fill("MppiChi2BcgRandomTrackRotationpT", mass, pT);
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
                        insideprocessing.Fill("MpipChi2BcgRandomTrackRotation", mass);
                        insideprocessing.Fill("MpipChi2BcgRandomTrackRotationeta", mass, eta);
                        insideprocessing.Fill("MpipChi2BcgRandomTrackRotationpT", mass, pT);
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
                        insideprocessing.Fill("MKKChi2BcgRandomTrackRotation", mass);
                        insideprocessing.Fill("MKKChi2BcgRandomTrackRotationClose", mass);
                        insideprocessing.Fill("MKKChi2BcgRandomTrackRotationeta", mass, eta);
                        insideprocessing.Fill("MKKChi2BcgRandomTrackRotationpT", mass, pT);
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
                        insideprocessing.Fill("MpipiChi2BcgRandomTrackRotation", mass);
                        insideprocessing.Fill("MpipiChi2BcgRandomTrackRotationeta", mass, eta);
                        insideprocessing.Fill("MpipiChi2BcgRandomTrackRotationpT", mass, pT);
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
                        insideprocessing.Fill("MppChi2BcgRandomTrackRotation", mass);
                        insideprocessing.Fill("MppChi2BcgRandomTrackRotationeta", mass, eta);
                        insideprocessing.Fill("MppChi2BcgRandomTrackRotationpT", mass, pT);
                    }
                }
            }

            //########## BACKGROUND EXTRACTION (MIXED EVENT) ###############

            //connecting current tracks with previous tracks
            for(size_t past_event = 0; past_event<queue_of_previous_vector_Tracks_positive.size(); past_event++){
                //current +, previous -
                for(long unsigned int i = 0; i<vector_Track_positive.size(); i++){
                    for(long unsigned int j = 0; j<queue_of_previous_vector_Tracks_negative[past_event].size(); j++){
                        isdEdxOk = (vector_Track_positive[i]->getNhitsDEdx()>=15)&&(queue_of_previous_vector_Tracks_negative[past_event][j]->getNhitsDEdx()>=15);
                        isTOFOk = (vector_Track_positive[i]->getTofPathLength()>0)&&(vector_Track_positive[i]->getTofTime()>0)&&(queue_of_previous_vector_Tracks_negative[past_event][j]->getTofPathLength()>0)&&(queue_of_previous_vector_Tracks_negative[past_event][j]->getTofTime()>0);
                        if(isdEdxOk&&isTOFOk){
                            //tracks are okay
                        } else{
                            continue;
                        }

                        //chi2
                        for(size_t pos = 0; pos<nParticlesExtended; pos++){
                            for(size_t neg = 0; neg<nParticlesExtended; neg++){
                                tempPairName = particleNicks[pos]+"_"+particleNicks[neg];
                                chi2Map[tempPairName] = getChi2(vector_Track_positive[i], queue_of_previous_vector_Tracks_negative[past_event][j], pos, neg, sigmaMap[tempPairName]);
                            }
                        }

                        //track rotation (rotating only the negative one)

                        //mass tests on different pairs
                        if(chi2Map["K_pi"]<9){
                            vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Kaon]);
                            queue_of_previous_vector_Tracks_negative[past_event][j]->getLorentzVector(negative_track, particleMass[Pion]);
                            //track rotation (rotating only the negative one)
                            negative_track.RotateZ(TMath::Pi());
                            //the rest as usual
                            mass = (positive_track+negative_track).M();
                            eta = (positive_track+negative_track).Eta();
                            pT = (positive_track+negative_track).Pt();
                            insideprocessing.Fill("MKpiChi2BcgMixedEvent", mass);
                            insideprocessing.Fill("MKpiChi2BcgMixedEventClose", mass);
                            insideprocessing.Fill("MKpiChi2BcgMixedEventeta", mass, eta);
                            insideprocessing.Fill("MKpiChi2BcgMixedEventpT", mass, pT);
                        }
                        if(chi2Map["pi_K"]<9){
                            vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                            queue_of_previous_vector_Tracks_negative[past_event][j]->getLorentzVector(negative_track, particleMass[Kaon]);
                            //track rotation (rotating only the negative one)
                            negative_track.RotateZ(TMath::Pi());
                            //the rest as usual
                            mass = (positive_track+negative_track).M();
                            eta = (positive_track+negative_track).Eta();
                            pT = (positive_track+negative_track).Pt();
                            insideprocessing.Fill("MpiKChi2BcgMixedEvent", mass);
                            insideprocessing.Fill("MpiKChi2BcgMixedEventClose", mass);
                            insideprocessing.Fill("MpiKChi2BcgMixedEventeta", mass, eta);
                            insideprocessing.Fill("MpiKChi2BcgMixedEventpT", mass, pT);
                        }
                        if(chi2Map["p_pi"]<9){
                            vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Proton]);
                            queue_of_previous_vector_Tracks_negative[past_event][j]->getLorentzVector(negative_track, particleMass[Pion]);
                            //track rotation (rotating only the negative one)
                            negative_track.RotateZ(TMath::Pi());
                            //the rest as usual
                            mass = (positive_track+negative_track).M();
                            eta = (positive_track+negative_track).Eta();
                            pT = (positive_track+negative_track).Pt();
                            insideprocessing.Fill("MppiChi2BcgMixedEvent", mass);
                            insideprocessing.Fill("MppiChi2BcgMixedEventeta", mass, eta);
                            insideprocessing.Fill("MppiChi2BcgMixedEventpT", mass, pT);
                        }
                        if(chi2Map["pi_p"]<9){
                            vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                            queue_of_previous_vector_Tracks_negative[past_event][j]->getLorentzVector(negative_track, particleMass[Proton]);
                            //track rotation (rotating only the negative one)
                            negative_track.RotateZ(TMath::Pi());
                            //the rest as usual
                            mass = (positive_track+negative_track).M();
                            eta = (positive_track+negative_track).Eta();
                            pT = (positive_track+negative_track).Pt();
                            insideprocessing.Fill("MpipChi2BcgMixedEvent", mass);
                            insideprocessing.Fill("MpipChi2BcgMixedEventeta", mass, eta);
                            insideprocessing.Fill("MpipChi2BcgMixedEventpT", mass, pT);
                        }
                        if(chi2Map["K_K"]<9){
                            vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Kaon]);
                            queue_of_previous_vector_Tracks_negative[past_event][j]->getLorentzVector(negative_track, particleMass[Kaon]);
                            //track rotation (rotating only the negative one)
                            negative_track.RotateZ(TMath::Pi());
                            //the rest as usual
                            mass = (positive_track+negative_track).M();
                            eta = (positive_track+negative_track).Eta();
                            pT = (positive_track+negative_track).Pt();
                            insideprocessing.Fill("MKKChi2BcgMixedEvent", mass);
                            insideprocessing.Fill("MKKChi2BcgMixedEventClose", mass);
                            insideprocessing.Fill("MKKChi2BcgMixedEventeta", mass, eta);
                            insideprocessing.Fill("MKKChi2BcgMixedEventpT", mass, pT);
                        }
                        if(chi2Map["pi_pi"]<9){
                            vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                            queue_of_previous_vector_Tracks_negative[past_event][j]->getLorentzVector(negative_track, particleMass[Pion]);
                            //track rotation (rotating only the negative one)
                            negative_track.RotateZ(TMath::Pi());
                            //the rest as usual
                            mass = (positive_track+negative_track).M();
                            eta = (positive_track+negative_track).Eta();
                            pT = (positive_track+negative_track).Pt();
                            insideprocessing.Fill("MpipiChi2BcgMixedEvent", mass);
                            insideprocessing.Fill("MpipiChi2BcgMixedEventeta", mass, eta);
                            insideprocessing.Fill("MpipiChi2BcgMixedEventpT", mass, pT);
                        }
                        if(chi2Map["p_p"]<9){
                            vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Proton]);
                            queue_of_previous_vector_Tracks_negative[past_event][j]->getLorentzVector(negative_track, particleMass[Proton]);
                            //track rotation (rotating only the negative one)
                            negative_track.RotateZ(TMath::Pi());
                            //the rest as usual
                            mass = (positive_track+negative_track).M();
                            eta = (positive_track+negative_track).Eta();
                            pT = (positive_track+negative_track).Pt();
                            insideprocessing.Fill("MppChi2BcgMixedEvent", mass);
                            insideprocessing.Fill("MppChi2BcgMixedEventeta", mass, eta);
                            insideprocessing.Fill("MppChi2BcgMixedEventpT", mass, pT);
                        }
                    }
                }
                //previous +, current -
                for(long unsigned int i = 0; i<queue_of_previous_vector_Tracks_positive[past_event].size(); i++){
                    for(long unsigned int j = 0; j<vector_Track_negative.size(); j++){
                        isdEdxOk = (queue_of_previous_vector_Tracks_positive[past_event][i]->getNhitsDEdx()>=15)&&(vector_Track_negative[j]->getNhitsDEdx()>=15);
                        isTOFOk = (queue_of_previous_vector_Tracks_positive[past_event][i]->getTofPathLength()>0)&&(queue_of_previous_vector_Tracks_positive[past_event][i]->getTofTime()>0)&&(vector_Track_negative[j]->getTofPathLength()>0)&&(vector_Track_negative[j]->getTofTime()>0);
                        if(isdEdxOk&&isTOFOk){
                            //tracks are okay
                        } else{
                            continue;
                        }

                        //chi2
                        for(size_t pos = 0; pos<nParticlesExtended; pos++){
                            for(size_t neg = 0; neg<nParticlesExtended; neg++){
                                tempPairName = particleNicks[pos]+"_"+particleNicks[neg];
                                chi2Map[tempPairName] = getChi2(queue_of_previous_vector_Tracks_positive[past_event][i], vector_Track_negative[j], pos, neg, sigmaMap[tempPairName]);
                            }
                        }

                        //track rotation (rotating only the negative one)

                        //mass tests on different pairs
                        if(chi2Map["K_pi"]<9){
                            queue_of_previous_vector_Tracks_positive[past_event][i]->getLorentzVector(positive_track, particleMass[Kaon]);
                            vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Pion]);
                            //track rotation (rotating only the negative one)
                            negative_track.RotateZ(TMath::Pi());
                            //the rest as usual
                            mass = (positive_track+negative_track).M();
                            eta = (positive_track+negative_track).Eta();
                            pT = (positive_track+negative_track).Pt();
                            insideprocessing.Fill("MKpiChi2BcgMixedEvent", mass);
                            insideprocessing.Fill("MKpiChi2BcgMixedEventClose", mass);
                            insideprocessing.Fill("MKpiChi2BcgMixedEventeta", mass, eta);
                            insideprocessing.Fill("MKpiChi2BcgMixedEventpT", mass, pT);
                        }
                        if(chi2Map["pi_K"]<9){
                            queue_of_previous_vector_Tracks_positive[past_event][i]->getLorentzVector(positive_track, particleMass[Pion]);
                            vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Kaon]);
                            //track rotation (rotating only the negative one)
                            negative_track.RotateZ(TMath::Pi());
                            //the rest as usual
                            mass = (positive_track+negative_track).M();
                            eta = (positive_track+negative_track).Eta();
                            pT = (positive_track+negative_track).Pt();
                            insideprocessing.Fill("MpiKChi2BcgMixedEvent", mass);
                            insideprocessing.Fill("MpiKChi2BcgMixedEventClose", mass);
                            insideprocessing.Fill("MpiKChi2BcgMixedEventeta", mass, eta);
                            insideprocessing.Fill("MpiKChi2BcgMixedEventpT", mass, pT);
                        }
                        if(chi2Map["p_pi"]<9){
                            queue_of_previous_vector_Tracks_positive[past_event][i]->getLorentzVector(positive_track, particleMass[Proton]);
                            vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Pion]);
                            //track rotation (rotating only the negative one)
                            negative_track.RotateZ(TMath::Pi());
                            //the rest as usual
                            mass = (positive_track+negative_track).M();
                            eta = (positive_track+negative_track).Eta();
                            pT = (positive_track+negative_track).Pt();
                            insideprocessing.Fill("MppiChi2BcgMixedEvent", mass);
                            insideprocessing.Fill("MppiChi2BcgMixedEventeta", mass, eta);
                            insideprocessing.Fill("MppiChi2BcgMixedEventpT", mass, pT);
                        }
                        if(chi2Map["pi_p"]<9){
                            queue_of_previous_vector_Tracks_positive[past_event][i]->getLorentzVector(positive_track, particleMass[Pion]);
                            vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Proton]);
                            //track rotation (rotating only the negative one)
                            negative_track.RotateZ(TMath::Pi());
                            //the rest as usual
                            mass = (positive_track+negative_track).M();
                            eta = (positive_track+negative_track).Eta();
                            pT = (positive_track+negative_track).Pt();
                            insideprocessing.Fill("MpipChi2BcgMixedEvent", mass);
                            insideprocessing.Fill("MpipChi2BcgMixedEventeta", mass, eta);
                            insideprocessing.Fill("MpipChi2BcgMixedEventpT", mass, pT);
                        }
                        if(chi2Map["K_K"]<9){
                            queue_of_previous_vector_Tracks_positive[past_event][i]->getLorentzVector(positive_track, particleMass[Kaon]);
                            vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Kaon]);
                            //track rotation (rotating only the negative one)
                            negative_track.RotateZ(TMath::Pi());
                            //the rest as usual
                            mass = (positive_track+negative_track).M();
                            eta = (positive_track+negative_track).Eta();
                            pT = (positive_track+negative_track).Pt();
                            insideprocessing.Fill("MKKChi2BcgMixedEvent", mass);
                            insideprocessing.Fill("MKKChi2BcgMixedEventClose", mass);
                            insideprocessing.Fill("MKKChi2BcgMixedEventeta", mass, eta);
                            insideprocessing.Fill("MKKChi2BcgMixedEventpT", mass, pT);
                        }
                        if(chi2Map["pi_pi"]<9){
                            queue_of_previous_vector_Tracks_positive[past_event][i]->getLorentzVector(positive_track, particleMass[Pion]);
                            vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Pion]);
                            //track rotation (rotating only the negative one)
                            negative_track.RotateZ(TMath::Pi());
                            //the rest as usual
                            mass = (positive_track+negative_track).M();
                            eta = (positive_track+negative_track).Eta();
                            pT = (positive_track+negative_track).Pt();
                            insideprocessing.Fill("MpipiChi2BcgMixedEvent", mass);
                            insideprocessing.Fill("MpipiChi2BcgMixedEventeta", mass, eta);
                            insideprocessing.Fill("MpipiChi2BcgMixedEventpT", mass, pT);
                        }
                        if(chi2Map["p_p"]<9){
                            queue_of_previous_vector_Tracks_positive[past_event][i]->getLorentzVector(positive_track, particleMass[Proton]);
                            vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Proton]);
                            //track rotation (rotating only the negative one)
                            negative_track.RotateZ(TMath::Pi());
                            //the rest as usual
                            mass = (positive_track+negative_track).M();
                            eta = (positive_track+negative_track).Eta();
                            pT = (positive_track+negative_track).Pt();
                            insideprocessing.Fill("MppChi2BcgMixedEvent", mass);
                            insideprocessing.Fill("MppChi2BcgMixedEventeta", mass, eta);
                            insideprocessing.Fill("MppChi2BcgMixedEventpT", mass, pT);
                        }
                    }
                }
            }
            //moving current tracks to previous tracks container
            //if the recall limit has been reached, we pop the oldest one
            if(queue_of_previous_vector_Tracks_positive.size()==previous_events_in_queue){
                for(size_t i = 0; i<queue_of_previous_vector_Tracks_positive.front().size(); i++){
                    delete queue_of_previous_vector_Tracks_positive.front()[i];
                }
                queue_of_previous_vector_Tracks_positive.pop_front();
                for(size_t i = 0; i<queue_of_previous_vector_Tracks_negative.front().size(); i++){
                    delete queue_of_previous_vector_Tracks_negative.front()[i];
                }
                queue_of_previous_vector_Tracks_negative.pop_front();
            }
            //after (eventual) popping of the oldest one, we load the queue with the copy of the current event one
            queue_of_previous_vector_Tracks_positive.emplace_back();
            for(size_t i = 0; i<vector_Track_positive.size(); i++){
                queue_of_previous_vector_Tracks_positive.back().push_back(new StUPCTrack());
                double pt, eta, phi;
                vector_Track_positive[i]->getPtEtaPhi(pt, eta, phi);
                queue_of_previous_vector_Tracks_positive.back()[i]->setPtEtaPhi(pt, eta, phi);
                queue_of_previous_vector_Tracks_positive.back()[i]->setNhitsDEdx(vector_Track_positive[i]->getNhitsDEdx());
                queue_of_previous_vector_Tracks_positive.back()[i]->setTofPathLength(vector_Track_positive[i]->getTofPathLength());
                queue_of_previous_vector_Tracks_positive.back()[i]->setTofTime(vector_Track_positive[i]->getTofTime());
                for(size_t part = 0; part<nParticlesExtended; part++){
                    queue_of_previous_vector_Tracks_positive.back()[i]->setNSigmasTPC(static_cast<StUPCTrack::Part>(part), vector_Track_positive[i]->getNSigmasTPC(static_cast<StUPCTrack::Part>(part)));
                }
            }
            queue_of_previous_vector_Tracks_negative.emplace_back();
            for(size_t i = 0; i<vector_Track_negative.size(); i++){
                queue_of_previous_vector_Tracks_negative.back().push_back(new StUPCTrack());
                double pt, eta, phi;
                vector_Track_negative[i]->getPtEtaPhi(pt, eta, phi);
                queue_of_previous_vector_Tracks_negative.back()[i]->setPtEtaPhi(pt, eta, phi);
                queue_of_previous_vector_Tracks_negative.back()[i]->setNhitsDEdx(vector_Track_negative[i]->getNhitsDEdx());
                queue_of_previous_vector_Tracks_negative.back()[i]->setTofPathLength(vector_Track_negative[i]->getTofPathLength());
                queue_of_previous_vector_Tracks_negative.back()[i]->setTofTime(vector_Track_negative[i]->getTofTime());
                for(size_t part = 0; part<nParticlesExtended; part++){
                    queue_of_previous_vector_Tracks_negative.back()[i]->setNSigmasTPC(static_cast<StUPCTrack::Part>(part), vector_Track_negative[i]->getNSigmasTPC(static_cast<StUPCTrack::Part>(part)));
                }
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

double getChi2(StUPCTrack* positive, StUPCTrack* negative, int positiveId, int negativeId, double sigmaT){
    double sigma1 = pow(positive->getNSigmasTPC(static_cast<StUPCTrack::Part>(positiveId)), 2);
    double sigma2 = pow(negative->getNSigmasTPC(static_cast<StUPCTrack::Part>(negativeId)), 2);
    double sigma3 = pow(DeltaT0(positive, negative, particleMassExtended[positiveId], particleMassExtended[negativeId])/sigmaT, 2);
    return sigma1+sigma2+sigma3;
}