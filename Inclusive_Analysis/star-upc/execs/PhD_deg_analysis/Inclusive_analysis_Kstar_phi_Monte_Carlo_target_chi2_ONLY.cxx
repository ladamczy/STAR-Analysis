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

TParticle* getMother(StUPCEvent* inputUPCpointer, TParticle* particle1, TParticle* particle2);
double getChi2(StUPCTrack* positive, StUPCTrack* negative, int positiveId, int negativeId, double sigmaT);

int main(int argc, char** argv){

    //argv:
    //1 - .root file or .list list
    //2 - output folder
    //3 - number of cores

    int nthreads = 1;
    if(argc>=4){
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
    //mass histograms (signal) with identification of the mothers
    outsideprocessing.AddHistogram(TH2D("MKpiMothersChi2", ";m_{K^{+}#pi^{-}} [GeV];Mother symbol", 200, 0.5, 2.0, 1, 0, 1));
    outsideprocessing.AddHistogram(TH2D("MpiKMothersChi2", ";m_{#pi^{+}K^{-}} [GeV];Mother symbol", 200, 0.5, 2.0, 1, 0, 1));
    outsideprocessing.AddHistogram(TH2D("MppiMothersChi2", ";m_{p^{+}#pi^{-}} [GeV];Mother symbol", 500, 1.0, 2.5, 1, 0, 1));
    outsideprocessing.AddHistogram(TH2D("MpipMothersChi2", ";m_{#pi^{+}p^{-}} [GeV];Mother symbol", 500, 1.0, 2.5, 1, 0, 1));
    outsideprocessing.AddHistogram(TH2D("MKKMothersChi2", ";m_{K^{+}K^{-}} [GeV];Mother symbol", 500, 0.9, 2.4, 1, 0, 1));
    outsideprocessing.AddHistogram(TH2D("MpipiMothersChi2", ";m_{#pi^{+}#pi^{-}} [GeV];Mother symbol", 600, 0.2, 1.4, 1, 0, 1));
    outsideprocessing.AddHistogram(TH2D("MppMothersChi2", ";m_{p^{+}p^{-}} [GeV];Mother symbol", 500, 1.5, 3.5, 1, 0, 1));

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

        //filling pairInfo histogram in correct order of bins
        for(auto&& histlastname:{ "Signal" }){
            std::string histfirstname = "pairInfo";
            insideprocessing.Fill((histfirstname+histlastname).c_str(), "OK", 0.0);
            insideprocessing.Fill((histfirstname+histlastname).c_str(), "TOF wrong", 0.0);
            insideprocessing.Fill((histfirstname+histlastname).c_str(), "dEdx wrong", 0.0);
            insideprocessing.Fill((histfirstname+histlastname).c_str(), "Both wrong", 0.0);
        }

        //helpful variables
        std::vector<StUPCTrack*> vector_Track_positive;
        std::vector<StUPCTrack*> vector_Track_negative;
        std::vector<TParticle*> vector_MC_Track_positive;
        std::vector<TParticle*> vector_MC_Track_negative;
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

        //actual loop
        while(myReader.Next()){
            //in a TTree, it *would* be constant, in TChain however not necessarily
            tempUPCpointer = StUPCEventInstance.Get();
            // tempRPpointer = StRPEventInstance.Get();

            //cleaning the loop
            vector_Track_positive.clear();
            vector_Track_negative.clear();
            vector_MC_Track_positive.clear();
            vector_MC_Track_negative.clear();
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
                    vector_MC_Track_positive.push_back(tempUPCpointer->getMCParticle(tempTrack->getIdTruth()-1));
                } else{
                    vector_Track_negative.push_back(tempTrack);
                    vector_MC_Track_negative.push_back(tempUPCpointer->getMCParticle(tempTrack->getIdTruth()-1));
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
                        std::string motherName = "nonresonant";
                        if(getMother(tempUPCpointer, vector_MC_Track_positive[i], vector_MC_Track_negative[j])!=nullptr){
                            motherName = std::string(getMother(tempUPCpointer, vector_MC_Track_positive[i], vector_MC_Track_negative[j])->GetName());
                        }
                        insideprocessing.Fill("MKpiMothersChi2", mass, motherName.c_str(), 1.0);
                    }
                    if(chi2Map["pi_K"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                        vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Kaon]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        insideprocessing.Fill("MpiKChi2", mass);
                        std::string motherName = "nonresonant";
                        if(getMother(tempUPCpointer, vector_MC_Track_positive[i], vector_MC_Track_negative[j])!=nullptr){
                            motherName = std::string(getMother(tempUPCpointer, vector_MC_Track_positive[i], vector_MC_Track_negative[j])->GetName());
                        }
                        insideprocessing.Fill("MpiKMothersChi2", mass, motherName.c_str(), 1.0);
                    }
                    if(chi2Map["p_pi"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Proton]);
                        vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Pion]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        insideprocessing.Fill("MppiChi2", mass);
                        std::string motherName = "nonresonant";
                        if(getMother(tempUPCpointer, vector_MC_Track_positive[i], vector_MC_Track_negative[j])!=nullptr){
                            motherName = std::string(getMother(tempUPCpointer, vector_MC_Track_positive[i], vector_MC_Track_negative[j])->GetName());
                        }
                        insideprocessing.Fill("MppiMothersChi2", mass, motherName.c_str(), 1.0);
                    }
                    if(chi2Map["pi_p"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                        vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Proton]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        insideprocessing.Fill("MpipChi2", mass);
                        std::string motherName = "nonresonant";
                        if(getMother(tempUPCpointer, vector_MC_Track_positive[i], vector_MC_Track_negative[j])!=nullptr){
                            motherName = std::string(getMother(tempUPCpointer, vector_MC_Track_positive[i], vector_MC_Track_negative[j])->GetName());
                        }
                        insideprocessing.Fill("MpipMothersChi2", mass, motherName.c_str(), 1.0);
                    }
                    if(chi2Map["K_K"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Kaon]);
                        vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Kaon]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        insideprocessing.Fill("MKKChi2", mass);
                        std::string motherName = "nonresonant";
                        if(getMother(tempUPCpointer, vector_MC_Track_positive[i], vector_MC_Track_negative[j])!=nullptr){
                            motherName = std::string(getMother(tempUPCpointer, vector_MC_Track_positive[i], vector_MC_Track_negative[j])->GetName());
                        }
                        insideprocessing.Fill("MKKMothersChi2", mass, motherName.c_str(), 1.0);
                        //test of suspicious peak and its neighbourhood
                        // if(chi2Map["K_K"]<3){
                        //     insideprocessing.Fill("MKKSuspiciousPeakTestedWithStrictChi2LessThan3", mass);
                        // }
                        // if(chi2Map["K_K"]<1){
                        //     insideprocessing.Fill("MKKSuspiciousPeakTestedWithStrictChi2LessThan1", mass);
                        // }
                        // if(1.06<mass&&mass<1.08){
                        //     vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                        //     vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Pion]);
                        //     mass = (positive_track+negative_track).M();
                        //     insideprocessing.Fill("MKKSuspiciousPeakTestedAsPionPair", mass);
                        // }
                        // if((1.05<mass&&mass<1.06)||(1.08<mass&&mass<1.09)){
                        //     vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                        //     vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Pion]);
                        //     mass = (positive_track+negative_track).M();
                        //     insideprocessing.Fill("MKKSuspiciousPeakTestedAsPionPairNeighbourhood", mass);
                        // }
                    }
                    if(chi2Map["pi_pi"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                        vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Pion]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        insideprocessing.Fill("MpipiChi2", mass);
                        std::string motherName = "nonresonant";
                        if(getMother(tempUPCpointer, vector_MC_Track_positive[i], vector_MC_Track_negative[j])!=nullptr){
                            motherName = std::string(getMother(tempUPCpointer, vector_MC_Track_positive[i], vector_MC_Track_negative[j])->GetName());
                        }
                        insideprocessing.Fill("MpipiMothersChi2", mass, motherName.c_str(), 1.0);
                    }
                    if(chi2Map["p_p"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Proton]);
                        vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Proton]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        insideprocessing.Fill("MppChi2", mass);
                        std::string motherName = "nonresonant";
                        if(getMother(tempUPCpointer, vector_MC_Track_positive[i], vector_MC_Track_negative[j])!=nullptr){
                            motherName = std::string(getMother(tempUPCpointer, vector_MC_Track_positive[i], vector_MC_Track_negative[j])->GetName());
                        }
                        insideprocessing.Fill("MppMothersChi2", mass, motherName.c_str(), 1.0);
                    }
                }
            }

            //lambda finish
        }
        return 0;
    };

    TreeProc.Process(myFunction);

    //merging and tidying up
    outsideprocessing.Merge();
    for(auto&& histpair:pairTab){
        outsideprocessing.GetPointerAfterMerge2D(("M"+histpair+"MothersChi2").c_str())->LabelsDeflate("Y");
    }

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

TParticle* getMother(StUPCEvent* inputUPCpointer, TParticle* particle1, TParticle* particle2){
    int particle1ProductionVertex = particle1->GetFirstMother();
    int particle2ProductionVertex = particle2->GetFirstMother();
    if(particle1ProductionVertex!=particle2ProductionVertex){
        return nullptr;
    }
    //loop through all MC particles in search of main particle daughters
    for(size_t MC_mother = 0; MC_mother<inputUPCpointer->getNumberOfMCParticles(); MC_mother++){
        //mother and daughter here refer to vertex numbers!!!
        int motherDecayVertex = inputUPCpointer->getMCParticle(MC_mother)->GetFirstDaughter();
        //if we found a daughter (a particle whose production vertex is its mothers' decay vertex) we save its MC container number
        if(motherDecayVertex==particle1ProductionVertex){
            return inputUPCpointer->getMCParticle(MC_mother);
        }
    }
    return nullptr;
}

double getChi2(StUPCTrack* positive, StUPCTrack* negative, int positiveId, int negativeId, double sigmaT){
    double sigma1 = pow(positive->getNSigmasTPC(static_cast<StUPCTrack::Part>(positiveId)), 2);
    double sigma2 = pow(negative->getNSigmasTPC(static_cast<StUPCTrack::Part>(negativeId)), 2);
    double sigma3 = pow(DeltaT0(positive, negative, particleMassExtended[positiveId], particleMassExtended[negativeId])/sigmaT, 2);
    return sigma1+sigma2+sigma3;
}