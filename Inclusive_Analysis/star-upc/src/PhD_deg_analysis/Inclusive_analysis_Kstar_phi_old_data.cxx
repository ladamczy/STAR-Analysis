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

int main(int argc, char **argv){

    int nthreads = 1;
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

    //histograms
    ProcessingOutsideLoop outsideprocessing;
    int n_ptBins = 9;
    double ptBins[] = { 0,0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8 };
    int n_etaBins = 10;
    double etaBins[] = { -1., -0.8, -0.6, -0.4, -0.2, 0., 0.2, 0.4, 0.6, 0.8, 1.0 };
    //phi, Kstar, Lambda, K0S
    outsideprocessing.AddHistogram(TH1D("MKKWide", ";m_{K^{+}K^{-}} [GeV];Number of pairs", 500, 0.0, 5.0));
    outsideprocessing.AddHistogram(TH1D("MKKWidedEdx", ";m_{K^{+}K^{-}} [GeV];Number of pairs", 500, 0.0, 5.0));
    outsideprocessing.AddHistogram(TH1D("MKKNarrow", ";m_{K^{+}K^{-}} [GeV];Number of pairs", 100, 0.9, 1.2));
    outsideprocessing.AddHistogram(TH1D("MKKNarrowdEdx", ";m_{K^{+}K^{-}} [GeV];Number of pairs", 100, 0.9, 1.2));
    outsideprocessing.AddHistogram(TH1D("MKpiWide", ";m_{K^{#pm}#pi^{#mp}} [GeV];Number of pairs", 500, 0.0, 5.0));
    outsideprocessing.AddHistogram(TH1D("MKpiWidedEdx", ";m_{K^{#pm}#pi^{#mp}} [GeV];Number of pairs", 500, 0.0, 5.0));
    outsideprocessing.AddHistogram(TH1D("MKpiNarrow", ";m_{K^{#pm}#pi^{#mp}} [GeV];Number of pairs", 100, 0.8, 1.0));
    outsideprocessing.AddHistogram(TH1D("MKpiNarrowdEdx", ";m_{K^{#pm}#pi^{#mp}} [GeV];Number of pairs", 100, 0.8, 1.0));
    outsideprocessing.AddHistogram(TH1D("MppiWide", ";m_{p^{#pm}#pi^{#mp}} [GeV];Number of pairs", 500, 0.0, 5.0));
    outsideprocessing.AddHistogram(TH1D("MppiWidedEdx", ";m_{p^{#pm}#pi^{#mp}} [GeV];Number of pairs", 500, 0.0, 5.0));
    outsideprocessing.AddHistogram(TH1D("MppiNarrow", ";m_{p^{#pm}#pi^{#mp}} [GeV];Number of pairs", 100, 0.9, 1.3));
    outsideprocessing.AddHistogram(TH1D("MppiNarrowdEdx", ";m_{p^{#pm}#pi^{#mp}} [GeV];Number of pairs", 100, 0.9, 1.3));
    outsideprocessing.AddHistogram(TH1D("MpipiWide", ";m_{#pi^{#pm}#pi^{#mp}} [GeV];Number of pairs", 500, 0.0, 5.0));
    outsideprocessing.AddHistogram(TH1D("MpipiWidedEdx", ";m_{#pi^{#pm}#pi^{#mp}} [GeV];Number of pairs", 500, 0.0, 5.0));
    outsideprocessing.AddHistogram(TH1D("MpipiNarrow", ";m_{#pi^{#pm}#pi^{#mp}} [GeV];Number of pairs", 100, 0.4, 0.6));
    outsideprocessing.AddHistogram(TH1D("MpipiNarrowdEdx", ";m_{#pi^{#pm}#pi^{#mp}} [GeV];Number of pairs", 100, 0.4, 0.6));
    outsideprocessing.AddHistogram(TH1D("MKKWideNoVeto", ";m_{K^{+}K^{-}} [GeV];Number of pairs", 500, 0.0, 5.0));
    outsideprocessing.AddHistogram(TH1D("MKKWidedEdxNoVeto", ";m_{K^{+}K^{-}} [GeV];Number of pairs", 500, 0.0, 5.0));
    outsideprocessing.AddHistogram(TH1D("MKKNarrowNoVeto", ";m_{K^{+}K^{-}} [GeV];Number of pairs", 100, 0.9, 1.2));
    outsideprocessing.AddHistogram(TH1D("MKKNarrowdEdxNoVeto", ";m_{K^{+}K^{-}} [GeV];Number of pairs", 100, 0.9, 1.2));
    outsideprocessing.AddHistogram(TH1D("MKpiWideNoVeto", ";m_{K^{#pm}#pi^{#mp}} [GeV];Number of pairs", 500, 0.0, 5.0));
    outsideprocessing.AddHistogram(TH1D("MKpiWidedEdxNoVeto", ";m_{K^{#pm}#pi^{#mp}} [GeV];Number of pairs", 500, 0.0, 5.0));
    outsideprocessing.AddHistogram(TH1D("MKpiNarrowNoVeto", ";m_{K^{#pm}#pi^{#mp}} [GeV];Number of pairs", 100, 0.8, 1.0));
    outsideprocessing.AddHistogram(TH1D("MKpiNarrowdEdxNoVeto", ";m_{K^{#pm}#pi^{#mp}} [GeV];Number of pairs", 100, 0.8, 1.0));
    //background
    outsideprocessing.AddHistogram(TH1D("MKKWideBackground", ";m_{K^{#pm}K^{#pm}} [GeV];Number of pairs", 500, 0.0, 5.0));
    outsideprocessing.AddHistogram(TH1D("MKKWidedEdxBackground", ";m_{K^{#pm}K^{#pm}} [GeV];Number of pairs", 500, 0.0, 5.0));
    outsideprocessing.AddHistogram(TH1D("MKKNarrowBackground", ";m_{K^{#pm}K^{#pm}} [GeV];Number of pairs", 100, 0.9, 1.2));
    outsideprocessing.AddHistogram(TH1D("MKKNarrowdEdxBackground", ";m_{K^{#pm}K^{#pm}} [GeV];Number of pairs", 100, 0.9, 1.2));
    outsideprocessing.AddHistogram(TH1D("MKpiWideBackground", ";m_{K^{#pm}#pi^{#pm}} [GeV];Number of pairs", 500, 0.0, 5.0));
    outsideprocessing.AddHistogram(TH1D("MKpiWidedEdxBackground", ";m_{K^{#pm}#pi^{#pm}} [GeV];Number of pairs", 500, 0.0, 5.0));
    outsideprocessing.AddHistogram(TH1D("MKpiNarrowBackground", ";m_{K^{#pm}#pi^{#pm}} [GeV];Number of pairs", 100, 0.8, 1.0));
    outsideprocessing.AddHistogram(TH1D("MKpiNarrowdEdxBackground", ";m_{K^{#pm}#pi^{#pm}} [GeV];Number of pairs", 100, 0.8, 1.0));
    //other
    outsideprocessing.AddHistogram(TH2D("Mphipt2DHist", "Mphipt2DHist", 50, 0.9, 1.1, n_ptBins, ptBins));
    outsideprocessing.AddHistogram(TH2D("Mphieta2DHist", "Mphieta2DHist", 50, 0.9, 1.1, n_etaBins, etaBins));

    //processing
    //defining TreeProcessor
    ROOT::TTreeProcessorMT TreeProc(*upcChain, nthreads);

    //other things
    int eventsProcessed = 0;
    double K0WindowLow = 0.48;
    double K0WindowHigh = 0.51;
    double Lambda0WindowLow = 1.05;
    double Lambda0WindowHigh = 1.15;

    //defining processing function
    auto myFunction = [&](TTreeReader &myReader){
        //getting values from TChain, in-loop histogram initialization
        TTreeReaderValue<StUPCEvent> StUPCEventInstance(myReader, "mUPCEvent");
        // TTreeReaderValue<StRPEvent> StRPEventInstance(myReader, "mRPEvent");
        ProcessingInsideLoop insideprocessing;
        StUPCEvent *tempUPCpointer;
        // StRPEvent *tempRPpointer;
        insideprocessing.GetLocalHistograms(&outsideprocessing);

        //helpful variables
        std::vector<StUPCTrack *> vector_Track_positive;
        std::vector<StUPCTrack *> vector_Track_negative;
        int vector_Track_positive_dEdx;
        int vector_Track_negative_dEdx;
        int vector_bcg_track_1_dEdx;
        int vector_bcg_track_2_dEdx;
        StUPCTrack *tempTrack;
        TLorentzVector positive_track;
        TLorentzVector negative_track;
        TLorentzVector bcg_track_1;
        TLorentzVector bcg_track_2;
        TVector3 tempMomentum;
        double mass;

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
                if(!tempTrack->getFlag(StUPCTrack::kPrimary)){
                    continue;
                }
                if(tempTrack->getPt()<=0.2 or abs(tempTrack->getEta())>=0.7){
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

            //loop through particles
            for(long unsigned int i = 0; i<vector_Track_positive.size(); i++){
                for(long unsigned int j = 0; j<vector_Track_negative.size(); j++){
                    //no veto
                    //Kstar
                    vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Kaon]);
                    vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Pion]);
                    mass = (positive_track+negative_track).M();
                    insideprocessing.Fill("MKpiWideNoVeto", mass);
                    insideprocessing.Fill("MKpiNarrowNoVeto", mass);
                    vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                    vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Kaon]);
                    mass = (positive_track+negative_track).M();
                    insideprocessing.Fill("MKpiWideNoVeto", mass);
                    insideprocessing.Fill("MKpiNarrowNoVeto", mass);
                    //phi
                    vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Kaon]);
                    vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Kaon]);
                    mass = (positive_track+negative_track).M();
                    insideprocessing.Fill("MKKWideNoVeto", mass);
                    insideprocessing.Fill("MKKNarrowNoVeto", mass);

                    //normal, with veto
                    //K0S (for vetoing)
                    vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                    vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Pion]);
                    mass = (positive_track+negative_track).M();
                    insideprocessing.Fill("MpipiWide", mass);
                    insideprocessing.Fill("MpipiNarrow", mass);
                    if(K0WindowLow<mass&&mass<K0WindowHigh){
                        continue;
                    }
                    //Lambda0 (for vetoing)
                    vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Proton]);
                    vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Pion]);
                    mass = (positive_track+negative_track).M();
                    insideprocessing.Fill("MppiWide", mass);
                    insideprocessing.Fill("MppiNarrow", mass);
                    if(Lambda0WindowLow<mass&&mass<Lambda0WindowHigh){
                        continue;
                    }
                    vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                    vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Proton]);
                    mass = (positive_track+negative_track).M();
                    insideprocessing.Fill("MppiWide", mass);
                    insideprocessing.Fill("MppiNarrow", mass);
                    if(Lambda0WindowLow<mass&&mass<Lambda0WindowHigh){
                        continue;
                    }
                    //Kstar
                    vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Kaon]);
                    vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Pion]);
                    mass = (positive_track+negative_track).M();
                    insideprocessing.Fill("MKpiWide", mass);
                    insideprocessing.Fill("MKpiNarrow", mass);
                    vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                    vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Kaon]);
                    mass = (positive_track+negative_track).M();
                    insideprocessing.Fill("MKpiWide", mass);
                    insideprocessing.Fill("MKpiNarrow", mass);
                    //phi
                    vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Kaon]);
                    vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Kaon]);
                    mass = (positive_track+negative_track).M();
                    insideprocessing.Fill("MKKWide", mass);
                    insideprocessing.Fill("MKKNarrow", mass);
                    if(mass<1.080){
                        insideprocessing.Fill("Mphipt2DHist", mass, (positive_track+negative_track).Pt());
                        insideprocessing.Fill("Mphieta2DHist", mass, (positive_track+negative_track).Eta());
                    }
                }
            }

            //loop through identified particles
            for(long unsigned int i = 0; i<vector_Track_positive.size(); i++){
                for(long unsigned int j = 0; j<vector_Track_negative.size(); j++){
                    if(vector_Track_positive[i]->getNhitsDEdx()<15 or vector_Track_negative[j]->getNhitsDEdx()<15){
                        continue;
                    }
                    //positive
                    vector_Track_positive_dEdx = Pion;
                    vector_Track_positive[i]->getMomentum(tempMomentum);
                    if(abs(vector_Track_positive[i]->getNSigmasTPCKaon())<3&&tempMomentum.Mag()<0.4){
                        vector_Track_positive_dEdx = Kaon;
                    } else if(abs(vector_Track_positive[i]->getNSigmasTPCProton())<3&&tempMomentum.Mag()<0.8){
                        vector_Track_positive_dEdx = Proton;
                    }
                    vector_Track_negative_dEdx = Pion;
                    vector_Track_negative[j]->getMomentum(tempMomentum);
                    if(abs(vector_Track_negative[j]->getNSigmasTPCKaon())<3&&tempMomentum.Mag()<0.4){
                        vector_Track_negative_dEdx = Kaon;
                    } else if(abs(vector_Track_negative[j]->getNSigmasTPCProton())<3&&tempMomentum.Mag()<0.8){
                        vector_Track_negative_dEdx = Proton;
                    }
                    //identified
                    vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[vector_Track_positive_dEdx]);
                    vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[vector_Track_negative_dEdx]);
                    mass = (positive_track+negative_track).M();

                    //no veto
                    //KK pair
                    if(vector_Track_positive_dEdx==Kaon&&vector_Track_negative_dEdx==Kaon){
                        insideprocessing.Fill("MKKWidedEdxNoVeto", mass);
                        insideprocessing.Fill("MKKNarrowdEdxNoVeto", mass);
                    }
                    //Kpi pair
                    if((vector_Track_positive_dEdx==Kaon&&vector_Track_negative_dEdx==Pion)or(vector_Track_positive_dEdx==Pion&&vector_Track_negative_dEdx==Kaon)){
                        insideprocessing.Fill("MKpiWidedEdxNoVeto", mass);
                        insideprocessing.Fill("MKpiNarrowdEdxNoVeto", mass);
                    }

                    //with veto
                    if(vector_Track_positive_dEdx==Pion&&vector_Track_negative_dEdx==Pion){
                        insideprocessing.Fill("MpipiWidedEdx", mass);
                        insideprocessing.Fill("MpipiNarrowdEdx", mass);
                        //cut for K0S
                        if(K0WindowLow<mass&&mass<K0WindowHigh){
                            continue;
                        }
                    }
                    if((vector_Track_positive_dEdx==Proton&&vector_Track_negative_dEdx==Pion)or(vector_Track_positive_dEdx==Pion&&vector_Track_negative_dEdx==Proton)){
                        insideprocessing.Fill("MppiWidedEdx", mass);
                        insideprocessing.Fill("MppiNarrowdEdx", mass);
                        //cut for Lambda0
                        if(Lambda0WindowLow<mass&&mass<Lambda0WindowHigh){
                            continue;
                        }
                    }
                    //KK pair
                    if(vector_Track_positive_dEdx==Kaon&&vector_Track_negative_dEdx==Kaon){
                        insideprocessing.Fill("MKKWidedEdx", mass);
                        insideprocessing.Fill("MKKNarrowdEdx", mass);
                    }
                    //Kpi pair
                    if((vector_Track_positive_dEdx==Kaon&&vector_Track_negative_dEdx==Pion)or(vector_Track_positive_dEdx==Pion&&vector_Track_negative_dEdx==Kaon)){
                        insideprocessing.Fill("MKpiWidedEdx", mass);
                        insideprocessing.Fill("MKpiNarrowdEdx", mass);
                    }
                }
            }

            //loop through background - positive
            for(long unsigned int i = 0; i<vector_Track_positive.size(); i++){
                for(long unsigned int j = i+1; j<vector_Track_positive.size(); j++){
                    //normal
                    //Kpi
                    vector_Track_positive[i]->getLorentzVector(bcg_track_1, particleMass[Kaon]);
                    vector_Track_positive[j]->getLorentzVector(bcg_track_2, particleMass[Pion]);
                    mass = (bcg_track_1+bcg_track_2).M();
                    insideprocessing.Fill("MKpiWideBackground", mass);
                    insideprocessing.Fill("MKpiNarrowBackground", mass);
                    //KK
                    vector_Track_positive[i]->getLorentzVector(bcg_track_1, particleMass[Kaon]);
                    vector_Track_positive[j]->getLorentzVector(bcg_track_2, particleMass[Kaon]);
                    mass = (bcg_track_1+bcg_track_2).M();
                    insideprocessing.Fill("MKKWideBackground", mass);
                    insideprocessing.Fill("MKKNarrowBackground", mass);

                    //dEdx
                    if(vector_Track_positive[i]->getNhitsDEdx()<15 or vector_Track_positive[j]->getNhitsDEdx()<15){
                        continue;
                    }
                    vector_bcg_track_1_dEdx = Pion;
                    vector_Track_positive[i]->getMomentum(tempMomentum);
                    if(abs(vector_Track_positive[i]->getNSigmasTPCKaon())<3&&tempMomentum.Mag()<0.4){
                        vector_bcg_track_1_dEdx = Kaon;
                    } else if(abs(vector_Track_positive[i]->getNSigmasTPCProton())<3&&tempMomentum.Mag()<0.8){
                        vector_bcg_track_1_dEdx = Proton;
                    }
                    vector_bcg_track_2_dEdx = Pion;
                    vector_Track_positive[j]->getMomentum(tempMomentum);
                    if(abs(vector_Track_positive[j]->getNSigmasTPCKaon())<3&&tempMomentum.Mag()<0.4){
                        vector_bcg_track_2_dEdx = Kaon;
                    } else if(abs(vector_Track_positive[j]->getNSigmasTPCProton())<3&&tempMomentum.Mag()<0.8){
                        vector_bcg_track_2_dEdx = Proton;
                    }
                    vector_Track_positive[i]->getLorentzVector(bcg_track_1, particleMass[vector_bcg_track_1_dEdx]);
                    vector_Track_positive[j]->getLorentzVector(bcg_track_2, particleMass[vector_bcg_track_2_dEdx]);
                    mass = (bcg_track_1+bcg_track_2).M();
                    //Kpi
                    if((vector_bcg_track_1_dEdx==Kaon&&vector_bcg_track_2_dEdx==Pion)or(vector_bcg_track_1_dEdx==Pion&&vector_bcg_track_2_dEdx==Kaon)){
                        insideprocessing.Fill("MKpiWidedEdxBackground", mass);
                        insideprocessing.Fill("MKpiNarrowdEdxBackground", mass);
                    }
                    //KK
                    if(vector_bcg_track_1_dEdx==Kaon&&vector_bcg_track_2_dEdx==Kaon){
                        insideprocessing.Fill("MKKWidedEdxBackground", mass);
                        insideprocessing.Fill("MKKNarrowdEdxBackground", mass);
                    }
                }
            }

            //loop through background - negative
            for(long unsigned int i = 0; i<vector_Track_negative.size(); i++){
                for(long unsigned int j = i+1; j<vector_Track_negative.size(); j++){
                    //Kpi
                    vector_Track_negative[i]->getLorentzVector(bcg_track_1, particleMass[Kaon]);
                    vector_Track_negative[j]->getLorentzVector(bcg_track_2, particleMass[Pion]);
                    mass = (bcg_track_1+bcg_track_2).M();
                    insideprocessing.Fill("MKpiWideBackground", mass);
                    insideprocessing.Fill("MKpiNarrowBackground", mass);
                    //KK
                    vector_Track_negative[i]->getLorentzVector(bcg_track_1, particleMass[Kaon]);
                    vector_Track_negative[j]->getLorentzVector(bcg_track_2, particleMass[Kaon]);
                    mass = (bcg_track_1+bcg_track_2).M();
                    insideprocessing.Fill("MKKWideBackground", mass);
                    insideprocessing.Fill("MKKNarrowBackground", mass);

                    //dEdx
                    if(vector_Track_negative[i]->getNhitsDEdx()<15 or vector_Track_negative[j]->getNhitsDEdx()<15){
                        continue;
                    }
                    vector_bcg_track_1_dEdx = Pion;
                    vector_Track_negative[i]->getMomentum(tempMomentum);
                    if(abs(vector_Track_negative[i]->getNSigmasTPCKaon())<3&&tempMomentum.Mag()<0.4){
                        vector_bcg_track_1_dEdx = Kaon;
                    } else if(abs(vector_Track_negative[i]->getNSigmasTPCProton())<3&&tempMomentum.Mag()<0.8){
                        vector_bcg_track_1_dEdx = Proton;
                    }
                    vector_bcg_track_2_dEdx = Pion;
                    vector_Track_negative[j]->getMomentum(tempMomentum);
                    if(abs(vector_Track_negative[j]->getNSigmasTPCKaon())<3&&tempMomentum.Mag()<0.4){
                        vector_bcg_track_2_dEdx = Kaon;
                    } else if(abs(vector_Track_negative[j]->getNSigmasTPCProton())<3&&tempMomentum.Mag()<0.8){
                        vector_bcg_track_2_dEdx = Proton;
                    }
                    vector_Track_negative[i]->getLorentzVector(bcg_track_1, particleMass[vector_bcg_track_1_dEdx]);
                    vector_Track_negative[j]->getLorentzVector(bcg_track_2, particleMass[vector_bcg_track_2_dEdx]);
                    mass = (bcg_track_1+bcg_track_2).M();
                    //Kpi
                    if((vector_bcg_track_1_dEdx==Kaon&&vector_bcg_track_2_dEdx==Pion)or(vector_bcg_track_1_dEdx==Pion&&vector_bcg_track_2_dEdx==Kaon)){
                        insideprocessing.Fill("MKpiWidedEdxBackground", mass);
                        insideprocessing.Fill("MKpiNarrowdEdxBackground", mass);
                    }
                    //KK
                    if(vector_bcg_track_1_dEdx==Kaon&&vector_bcg_track_2_dEdx==Kaon){
                        insideprocessing.Fill("MKKWidedEdxBackground", mass);
                        insideprocessing.Fill("MKKNarrowdEdxBackground", mass);
                    }
                }
            }

            //lambda finish
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