//cpp headers
#include <map>
#include <fstream>

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
    outsideprocessing.AddHistogram(TH1D("pairInfoBackground", "", 1, 0, 1));
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
    //mass histograms (background)
    outsideprocessing.AddHistogram(TH1D("MKpiChi2Bcg", ";m_{K#pi} [GeV];Number of pairs", 200, 0.5, 2.0));
    outsideprocessing.AddHistogram(TH1D("MpiKChi2Bcg", ";m_{#piK} [GeV];Number of pairs", 200, 0.5, 2.0));
    outsideprocessing.AddHistogram(TH1D("MppiChi2Bcg", ";m_{p#pi} [GeV];Number of pairs", 500, 1.0, 2.5));
    outsideprocessing.AddHistogram(TH1D("MpipChi2Bcg", ";m_{#pip} [GeV];Number of pairs", 500, 1.0, 2.5));
    outsideprocessing.AddHistogram(TH1D("MKKChi2Bcg", ";m_{KK} [GeV];Number of pairs", 500, 0.9, 2.4));
    outsideprocessing.AddHistogram(TH1D("MpipiChi2Bcg", ";m_{#pi#pi} [GeV];Number of pairs", 600, 0.2, 1.4));
    outsideprocessing.AddHistogram(TH1D("MppChi2Bcg", ";m_{pp} [GeV];Number of pairs", 500, 1.5, 3.5));
    //closer histograms (background)
    outsideprocessing.AddHistogram(TH1D("MKKChi2BcgClose", ";m_{KK} [GeV];Number of pairs", 50, 0.99, 1.05));
    outsideprocessing.AddHistogram(TH1D("MKpiChi2BcgClose", ";m_{K#pi} [GeV];Number of pairs", 50, 0.7, 1.1));
    outsideprocessing.AddHistogram(TH1D("MpiKChi2BcgClose", ";m_{#piK} [GeV];Number of pairs", 50, 0.7, 1.1));
    //adding mass histograms grouped by category (background)
    getCategoryHistograms(outsideprocessing, pairTab, "Bcg");

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
            if(tempUPCpointer->isTrigger(570704)){
                continue;
            }
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
            //exactly one vertex
            if(tempUPCpointer->getNumberOfVertices()!=1){
                continue;
            }

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

            //loop through identified particles (background, same-sign positive)
            for(long int i = 0; i+1<vector_Track_positive.size(); i++){
                for(long int j = i+1; j<vector_Track_positive.size(); j++){
                    isdEdxOk = (vector_Track_positive[i]->getNhitsDEdx()>=15)&&(vector_Track_positive[j]->getNhitsDEdx()>=15);
                    isTOFOk = (vector_Track_positive[i]->getTofPathLength()>0)&&(vector_Track_positive[i]->getTofTime()>0)&&(vector_Track_positive[j]->getTofPathLength()>0)&&(vector_Track_positive[j]->getTofTime()>0);
                    if(isdEdxOk&&isTOFOk){
                        insideprocessing.Fill("pairInfoBackground", "OK", 1.0);
                    } else if(isdEdxOk&&!isTOFOk){
                        insideprocessing.Fill("pairInfoBackground", "TOF wrong", 1.0);
                        continue;
                    } else if(!isdEdxOk&&isTOFOk){
                        insideprocessing.Fill("pairInfoBackground", "dEdx wrong", 1.0);
                        continue;
                    } else if(!isdEdxOk&&!isTOFOk){
                        insideprocessing.Fill("pairInfoBackground", "Both wrong", 1.0);
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
                        insideprocessing.Fill("MKpiChi2Bcg", mass);
                        insideprocessing.Fill("MKpiChi2BcgClose", mass);
                        insideprocessing.Fill("MKpiChi2Bcgeta", mass, eta);
                        insideprocessing.Fill("MKpiChi2BcgpT", mass, pT);
                    }
                    if(chi2Map["pi_K"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                        vector_Track_positive[j]->getLorentzVector(positive_track2, particleMass[Kaon]);
                        mass = (positive_track+positive_track2).M();
                        eta = (positive_track+positive_track2).Eta();
                        pT = (positive_track+positive_track2).Pt();
                        insideprocessing.Fill("MpiKChi2Bcg", mass);
                        insideprocessing.Fill("MpiKChi2BcgClose", mass);
                        insideprocessing.Fill("MpiKChi2Bcgeta", mass, eta);
                        insideprocessing.Fill("MpiKChi2BcgpT", mass, pT);
                    }
                    if(chi2Map["p_pi"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Proton]);
                        vector_Track_positive[j]->getLorentzVector(positive_track2, particleMass[Pion]);
                        mass = (positive_track+positive_track2).M();
                        eta = (positive_track+positive_track2).Eta();
                        pT = (positive_track+positive_track2).Pt();
                        insideprocessing.Fill("MppiChi2Bcg", mass);
                        insideprocessing.Fill("MppiChi2Bcgeta", mass, eta);
                        insideprocessing.Fill("MppiChi2BcgpT", mass, pT);
                    }
                    if(chi2Map["pi_p"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                        vector_Track_positive[j]->getLorentzVector(positive_track2, particleMass[Proton]);
                        mass = (positive_track+positive_track2).M();
                        eta = (positive_track+positive_track2).Eta();
                        pT = (positive_track+positive_track2).Pt();
                        insideprocessing.Fill("MpipChi2Bcg", mass);
                        insideprocessing.Fill("MpipChi2Bcgeta", mass, eta);
                        insideprocessing.Fill("MpipChi2BcgpT", mass, pT);
                    }
                    if(chi2Map["K_K"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Kaon]);
                        vector_Track_positive[j]->getLorentzVector(positive_track2, particleMass[Kaon]);
                        mass = (positive_track+positive_track2).M();
                        eta = (positive_track+positive_track2).Eta();
                        pT = (positive_track+positive_track2).Pt();
                        insideprocessing.Fill("MKKChi2Bcg", mass);
                        insideprocessing.Fill("MKKChi2BcgClose", mass);
                        insideprocessing.Fill("MKKChi2Bcgeta", mass, eta);
                        insideprocessing.Fill("MKKChi2BcgpT", mass, pT);
                    }
                    if(chi2Map["pi_pi"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                        vector_Track_positive[j]->getLorentzVector(positive_track2, particleMass[Pion]);
                        mass = (positive_track+positive_track2).M();
                        eta = (positive_track+positive_track2).Eta();
                        pT = (positive_track+positive_track2).Pt();
                        insideprocessing.Fill("MpipiChi2Bcg", mass);
                        insideprocessing.Fill("MpipiChi2Bcgeta", mass, eta);
                        insideprocessing.Fill("MpipiChi2BcgpT", mass, pT);
                    }
                    if(chi2Map["p_p"]<9){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Proton]);
                        vector_Track_positive[j]->getLorentzVector(positive_track2, particleMass[Proton]);
                        mass = (positive_track+positive_track2).M();
                        eta = (positive_track+positive_track2).Eta();
                        pT = (positive_track+positive_track2).Pt();
                        insideprocessing.Fill("MppChi2Bcg", mass);
                        insideprocessing.Fill("MppChi2Bcgeta", mass, eta);
                        insideprocessing.Fill("MppChi2BcgpT", mass, pT);
                    }
                }
            }
            //loop through identified particles (background, same-sign negative)
            for(long int i = 0; i+1<vector_Track_negative.size(); i++){
                for(long int j = i+1; j<vector_Track_negative.size(); j++){
                    isdEdxOk = (vector_Track_negative[i]->getNhitsDEdx()>=15)&&(vector_Track_negative[j]->getNhitsDEdx()>=15);
                    isTOFOk = (vector_Track_negative[i]->getTofPathLength()>0)&&(vector_Track_negative[i]->getTofTime()>0)&&(vector_Track_negative[j]->getTofPathLength()>0)&&(vector_Track_negative[j]->getTofTime()>0);
                    if(isdEdxOk&&isTOFOk){
                        insideprocessing.Fill("pairInfoBackground", "OK", 1.0);
                    } else if(isdEdxOk&&!isTOFOk){
                        insideprocessing.Fill("pairInfoBackground", "TOF wrong", 1.0);
                        continue;
                    } else if(!isdEdxOk&&isTOFOk){
                        insideprocessing.Fill("pairInfoBackground", "dEdx wrong", 1.0);
                        continue;
                    } else if(!isdEdxOk&&!isTOFOk){
                        insideprocessing.Fill("pairInfoBackground", "Both wrong", 1.0);
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
                        insideprocessing.Fill("MKpiChi2Bcgeta", mass, eta);
                        insideprocessing.Fill("MKpiChi2Bcg", mass);
                        insideprocessing.Fill("MKpiChi2BcgClose", mass);
                        insideprocessing.Fill("MKpiChi2BcgpT", mass, pT);
                    }
                    if(chi2Map["pi_K"]<9){
                        vector_Track_negative[i]->getLorentzVector(negative_track, particleMass[Pion]);
                        vector_Track_negative[j]->getLorentzVector(negative_track2, particleMass[Kaon]);
                        mass = (negative_track+negative_track2).M();
                        eta = (negative_track+negative_track2).Eta();
                        pT = (negative_track+negative_track2).Pt();
                        insideprocessing.Fill("MpiKChi2Bcg", mass);
                        insideprocessing.Fill("MpiKChi2BcgClose", mass);
                        insideprocessing.Fill("MpiKChi2Bcgeta", mass, eta);
                        insideprocessing.Fill("MpiKChi2BcgpT", mass, pT);
                    }
                    if(chi2Map["p_pi"]<9){
                        vector_Track_negative[i]->getLorentzVector(negative_track, particleMass[Proton]);
                        vector_Track_negative[j]->getLorentzVector(negative_track2, particleMass[Pion]);
                        mass = (negative_track+negative_track2).M();
                        eta = (negative_track+negative_track2).Eta();
                        pT = (negative_track+negative_track2).Pt();
                        insideprocessing.Fill("MppiChi2Bcg", mass);
                        insideprocessing.Fill("MppiChi2Bcgeta", mass, eta);
                        insideprocessing.Fill("MppiChi2BcgpT", mass, pT);
                    }
                    if(chi2Map["pi_p"]<9){
                        vector_Track_negative[i]->getLorentzVector(negative_track, particleMass[Pion]);
                        vector_Track_negative[j]->getLorentzVector(negative_track2, particleMass[Proton]);
                        mass = (negative_track+negative_track2).M();
                        eta = (negative_track+negative_track2).Eta();
                        pT = (negative_track+negative_track2).Pt();
                        insideprocessing.Fill("MpipChi2Bcg", mass);
                        insideprocessing.Fill("MpipChi2Bcgeta", mass, eta);
                        insideprocessing.Fill("MpipChi2BcgpT", mass, pT);
                    }
                    if(chi2Map["K_K"]<9){
                        vector_Track_negative[i]->getLorentzVector(negative_track, particleMass[Kaon]);
                        vector_Track_negative[j]->getLorentzVector(negative_track2, particleMass[Kaon]);
                        mass = (negative_track+negative_track2).M();
                        eta = (negative_track+negative_track2).Eta();
                        pT = (negative_track+negative_track2).Pt();
                        insideprocessing.Fill("MKKChi2Bcg", mass);
                        insideprocessing.Fill("MKKChi2BcgClose", mass);
                        insideprocessing.Fill("MKKChi2Bcgeta", mass, eta);
                        insideprocessing.Fill("MKKChi2BcgpT", mass, pT);
                    }
                    if(chi2Map["pi_pi"]<9){
                        vector_Track_negative[i]->getLorentzVector(negative_track, particleMass[Pion]);
                        vector_Track_negative[j]->getLorentzVector(negative_track2, particleMass[Pion]);
                        mass = (negative_track+negative_track2).M();
                        eta = (negative_track+negative_track2).Eta();
                        pT = (negative_track+negative_track2).Pt();
                        insideprocessing.Fill("MpipiChi2Bcg", mass);
                        insideprocessing.Fill("MpipiChi2Bcgeta", mass, eta);
                        insideprocessing.Fill("MpipiChi2BcgpT", mass, pT);
                    }
                    if(chi2Map["p_p"]<9){
                        vector_Track_negative[i]->getLorentzVector(negative_track, particleMass[Proton]);
                        vector_Track_negative[j]->getLorentzVector(negative_track2, particleMass[Proton]);
                        mass = (negative_track+negative_track2).M();
                        eta = (negative_track+negative_track2).Eta();
                        pT = (negative_track+negative_track2).Pt();
                        insideprocessing.Fill("MppChi2Bcg", mass);
                        insideprocessing.Fill("MppChi2Bcgeta", mass, eta);
                        insideprocessing.Fill("MppChi2BcgpT", mass, pT);
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