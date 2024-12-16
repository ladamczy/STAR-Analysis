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
    outsideprocessing.AddHistogram(TH1D("deltaT0Kpi", ";#Deltat_{0} [ns];Number of pairs", 100, -10, 10));
    outsideprocessing.AddHistogram(TH1D("deltaT0KpiNarrow", ";#Deltat_{0} [ns];Number of pairs", 100, -1, 1));
    outsideprocessing.AddHistogram(TH1D("MKpiTest", ";m_{K^{#pm}#pi^{#mp}} [GeV];Number of pairs", 100, 0.7, 1.1));
    outsideprocessing.AddHistogram(TH1D("deltaT0KKNarrow", ";#Deltat_{0} [ns];Number of pairs", 100, -1, 1));
    outsideprocessing.AddHistogram(TH1D("MKKTest", ";m_{K^{+}K^{-}} [GeV];Number of pairs", 100, 0.8, 1.2));
    outsideprocessing.AddHistogram(TH1D("MKKTest2", ";m_{TOF} [GeV];Number of pairs", 100, 0.3, 0.7));
    outsideprocessing.AddHistogram(TH1D("MKKTest3", ";m_{K^{+}K^{-}} [GeV];Number of pairs", 100, 0.8, 1.2));

    //background
    outsideprocessing.AddHistogram(TH1D("MKKWideBackground", ";m_{K^{#pm}K^{#pm}} [GeV];Number of pairs", 500, 0.0, 5.0));
    outsideprocessing.AddHistogram(TH1D("MKKWidedEdxBackground", ";m_{K^{#pm}K^{#pm}} [GeV];Number of pairs", 500, 0.0, 5.0));
    outsideprocessing.AddHistogram(TH1D("MKKNarrowBackground", ";m_{K^{#pm}K^{#pm}} [GeV];Number of pairs", 100, 0.9, 1.2));
    outsideprocessing.AddHistogram(TH1D("MKKNarrowdEdxBackground", ";m_{K^{#pm}K^{#pm}} [GeV];Number of pairs", 100, 0.9, 1.2));
    outsideprocessing.AddHistogram(TH1D("MKpiWideBackground", ";m_{K^{#pm}#pi^{#pm}} [GeV];Number of pairs", 500, 0.0, 5.0));
    outsideprocessing.AddHistogram(TH1D("MKpiWidedEdxBackground", ";m_{K^{#pm}#pi^{#pm}} [GeV];Number of pairs", 500, 0.0, 5.0));
    outsideprocessing.AddHistogram(TH1D("MKpiNarrowBackground", ";m_{K^{#pm}#pi^{#pm}} [GeV];Number of pairs", 100, 0.8, 1.0));
    outsideprocessing.AddHistogram(TH1D("MKpiNarrowdEdxBackground", ";m_{K^{#pm}#pi^{#pm}} [GeV];Number of pairs", 100, 0.8, 1.0));
    //phi(1020) histograms
    outsideprocessing.AddHistogram(TH2D("Mphipt2DHist", "Mphipt2DHist", 50, 0.9, 1.1, n_ptBins, ptBins));
    outsideprocessing.AddHistogram(TH2D("Mphieta2DHist", "Mphieta2DHist", 50, 0.9, 1.1, n_etaBins, etaBins));
    //other
    outsideprocessing.AddHistogram(TH2D("MissingpXY", "MissingpXY", 100, -5., 5., 100, -5., 5.));
    outsideprocessing.AddHistogram(TH1D("MissingpT", "MissingpT", 100, 0., 1.));
    outsideprocessing.AddHistogram(TH1D("MissingpTWide", "MissingpT", 250, 0., 5.));

    outsideprocessing.AddHistogram(TH1D("MissingpTtests", "MissingpT", 100, 0., 1.));
    outsideprocessing.AddHistogram(TH1D("MissingpZ", "MissingpZ", 100, -5, 5.));
    outsideprocessing.AddHistogram(TH1D("MKKExtraNarrowNoVeto", ";m_{K^{+}K^{-}} [GeV];Number of pairs", 140, 0.98, 1.05));
    outsideprocessing.AddHistogram(TH2D("xi", ";#xi_{E};#xi_{W}", 150, -0.05, 0.25, 150, -0.05, 0.25));

    outsideprocessing.AddHistogram(TH1D("t0piTest", ";t-t_{0} [ns]", 100, 0, 50));
    outsideprocessing.AddHistogram(TH1D("t0piDifference", ";t_{0, 2nd}-t_{0, 1st} [ns]", 100, 0, 5));
    outsideprocessing.AddHistogram(TH1D("deltaT0pipi", ";t_{0}^{+}-t_{0}^{-} [ns]", 200, -20, 20));
    outsideprocessing.AddHistogram(TH1D("deltaT0Kpi2", ";t_{0}^{+}-t_{0}^{-} [ns]", 200, -20, 20));
    outsideprocessing.AddHistogram(TH1D("deltaT0KK", ";t_{0}^{+}-t_{0}^{-} [ns]", 200, -20, 20));
    outsideprocessing.AddHistogram(TH1D("deltaT0ppi", ";t_{0}^{+}-t_{0}^{-} [ns]", 200, -20, 20));
    TH1D temp("choice", ";choice", 16, 0, 16);
    // int pairchoice = ppiPair*8+pipiPair*4+KpiPair*2+KKPair;
    temp.GetXaxis()->SetBinLabel(1, "Nothing");
    temp.GetXaxis()->SetBinLabel(2, "KK");
    temp.GetXaxis()->SetBinLabel(3, "K#pi");
    temp.GetXaxis()->SetBinLabel(4, "K#pi+KK");
    temp.GetXaxis()->SetBinLabel(5, "#pi#pi");
    temp.GetXaxis()->SetBinLabel(6, "#pi#pi+KK");
    temp.GetXaxis()->SetBinLabel(7, "#pi#pi+Kpi");
    temp.GetXaxis()->SetBinLabel(8, "#pi#pi+Kpi+KK");
    temp.GetXaxis()->SetBinLabel(9, "p#pi");
    temp.GetXaxis()->SetBinLabel(11, "p#pi+K#pi");
    temp.GetXaxis()->SetBinLabel(15, "p#pi+KK+K#pi");
    temp.GetXaxis()->SetBinLabel(16, "All");
    outsideprocessing.AddHistogram(temp);
    outsideprocessing.AddHistogram(TH1D("MKKNarrowChoice", ";m_{K^{#pm}K^{#pm}} [GeV];Number of pairs", 100, 0.9, 1.7));
    outsideprocessing.AddHistogram(TH1D("MKpiNarrowChoice", ";m_{K^{#pm}#pi^{#mp}} [GeV];Number of pairs", 25, 0.8, 1.0));
    outsideprocessing.AddHistogram(TH1D("MppiNarrowChoice", ";m_{K^{#pm}#pi^{#mp}} [GeV];Number of pairs", 75, 0.9, 1.5));
    outsideprocessing.AddHistogram(TH1D("MpipiChoice", ";m_{#pi^{#pm}#pi^{#mp}} [GeV];Number of pairs", 100, 0., 2.));
    outsideprocessing.AddHistogram(TH1D("MpipiNoChoice", ";m_{#pi^{#pm}#pi^{#mp}} [GeV];Number of pairs", 100, 0., 2.));
    outsideprocessing.AddHistogram(TH1D("M2TOFpipiKK", ";m_{TOF}^{2} [GeV^{2}];Number of pairs", 100, -0.25, 0.75));
    outsideprocessing.AddHistogram(TH2D("chi2pipiKK", ";n#sigma_{#pi};n#sigma_{K}", 50, 0, 10, 50, 0, 10));

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
        TTreeReaderValue<StRPEvent> StRPEventInstance(myReader, "mRPEvent");
        ProcessingInsideLoop insideprocessing;
        StUPCEvent *tempUPCpointer;
        StRPEvent* tempRPpointer;
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
        TVector3 totalMomentum;
        double mass;
        double deltaT0;
        double m2TOF;
        std::vector<double> vector_t0_all_pions;
        double min_t0;

        //actual loop
        while(myReader.Next()){
            //in a TTree, it *would* be constant, in TChain however not necessarily
            tempUPCpointer = StUPCEventInstance.Get();
            tempRPpointer = StRPEventInstance.Get();

            //cleaning the loop
            vector_Track_positive.clear();
            vector_Track_negative.clear();
            totalMomentum = { 0, 0, 0 };

            //cause I want to see what's going on
            if(eventsProcessed%10000==0){
                cout<<"Processed "<<eventsProcessed<<" events"<<endl;
            }
            eventsProcessed++;

            //additional cuts that normally are used in only-RP-cuts examples
            //cause the data i used isn't properly filtered
            //at least 2 good tracks
            int nOfGoodTracks = 0;
            for(int i = 0; i<tempUPCpointer->getNumberOfTracks(); i++){
                StUPCTrack* tmptrk = tempUPCpointer->getTrack(i);
                if(tmptrk->getFlag(StUPCTrack::kTof)&&abs(tmptrk->getEta())<0.7&&tmptrk->getPt()>0.2){
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
            //also missing momenta fill
            tempMomentum.SetXYZ(0, 0, 0);
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
                tempTrack->getMomentum(tempMomentum);
                totalMomentum += tempMomentum;
            }

            //missing pt
            totalMomentum += tempRPpointer->getTrack(0)->pVec();
            totalMomentum += tempRPpointer->getTrack(1)->pVec();
            insideprocessing.Fill("MissingpT", totalMomentum.Pt());
            insideprocessing.Fill("MissingpTWide", totalMomentum.Pt());
            insideprocessing.Fill("MissingpXY", totalMomentum.X(), totalMomentum.Y());
            insideprocessing.Fill("MissingpZ", totalMomentum.Z());
            //tests to bring back the peak at lowpT
            totalMomentum = tempRPpointer->getTrack(0)->pVec();
            totalMomentum += tempRPpointer->getTrack(1)->pVec();
            for(long unsigned int i = 0; i<vector_Track_positive.size(); i++){
                if(vector_Track_positive[i]->getNhitsDEdx()>=15){
                    vector_Track_positive[i]->getMomentum(tempMomentum);
                    totalMomentum += tempMomentum;
                }
            }
            for(long unsigned int j = 0; j<vector_Track_negative.size(); j++){
                if(vector_Track_negative[j]->getNhitsDEdx()>=15){
                    vector_Track_negative[j]->getMomentum(tempMomentum);
                    totalMomentum += tempMomentum;
                }
            }
            insideprocessing.Fill("MissingpTtests", totalMomentum.Pt());
            if(tempRPpointer->getTrack(0)->pVec().Z()>0){
                insideprocessing.Fill("xi", tempRPpointer->getTrack(1)->xi(254.867), tempRPpointer->getTrack(0)->xi(254.867));
            } else{
                insideprocessing.Fill("xi", tempRPpointer->getTrack(0)->xi(254.867), tempRPpointer->getTrack(1)->xi(254.867));
            }

            //loop through particles
            for(long unsigned int i = 0; i<vector_Track_positive.size(); i++){
                for(long unsigned int j = 0; j<vector_Track_negative.size(); j++){
                    //no veto
                    //Kstar K+pi-
                    vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Kaon]);
                    vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Pion]);
                    mass = (positive_track+negative_track).M();
                    insideprocessing.Fill("MKpiWideNoVeto", mass);
                    insideprocessing.Fill("MKpiNarrowNoVeto", mass);
                    deltaT0 = DeltaT0(vector_Track_positive[i], vector_Track_negative[j], particleMass[Kaon], particleMass[Pion]);
                    insideprocessing.Fill("deltaT0Kpi", deltaT0);
                    insideprocessing.Fill("deltaT0KpiNarrow", deltaT0);
                    if(abs(deltaT0)<0.35){
                        insideprocessing.Fill("MKpiTest", mass);
                    }
                    //Kstar K-pi+
                    vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                    vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Kaon]);
                    mass = (positive_track+negative_track).M();
                    insideprocessing.Fill("MKpiWideNoVeto", mass);
                    insideprocessing.Fill("MKpiNarrowNoVeto", mass);
                    deltaT0 = DeltaT0(vector_Track_positive[i], vector_Track_negative[j], particleMass[Pion], particleMass[Kaon]);
                    insideprocessing.Fill("deltaT0Kpi", deltaT0);
                    insideprocessing.Fill("deltaT0KpiNarrow", deltaT0);
                    if(abs(deltaT0)<0.35){
                        insideprocessing.Fill("MKpiTest", mass);
                    }
                    //phi
                    vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Kaon]);
                    vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Kaon]);
                    mass = (positive_track+negative_track).M();
                    insideprocessing.Fill("MKKWideNoVeto", mass);
                    insideprocessing.Fill("MKKNarrowNoVeto", mass);
                    insideprocessing.Fill("MKKExtraNarrowNoVeto", mass);
                    deltaT0 = DeltaT0(vector_Track_positive[i], vector_Track_negative[j], particleMass[Kaon], particleMass[Kaon]);
                    insideprocessing.Fill("deltaT0KKNarrow", deltaT0);
                    if(abs(deltaT0)<0.35){
                        insideprocessing.Fill("MKKTest", mass);
                    }
                    m2TOF = M2TOF(vector_Track_positive[i], vector_Track_negative[j]);
                    insideprocessing.Fill("MKKTest2", sqrt(m2TOF));
                    if(m2TOF<0.53&&m2TOF>0.47){
                        insideprocessing.Fill("MKKTest3", mass);
                    }

                    //other things
                    if(vector_Track_positive[i]->getTofPathLength()>0&&vector_Track_positive[i]->getTofTime()>0&&vector_Track_negative[j]->getTofPathLength()>0&&vector_Track_negative[j]->getTofTime()){
                        insideprocessing.Fill("deltaT0pipi", DeltaT0(vector_Track_positive[i], vector_Track_negative[j], particleMass[Pion], particleMass[Pion]));
                        insideprocessing.Fill("deltaT0Kpi2", DeltaT0(vector_Track_positive[i], vector_Track_negative[j], particleMass[Pion], particleMass[Kaon]));
                        insideprocessing.Fill("deltaT0Kpi2", DeltaT0(vector_Track_positive[i], vector_Track_negative[j], particleMass[Kaon], particleMass[Pion]));
                        insideprocessing.Fill("deltaT0KK", DeltaT0(vector_Track_positive[i], vector_Track_negative[j], particleMass[Kaon], particleMass[Kaon]));
                        insideprocessing.Fill("deltaT0ppi", DeltaT0(vector_Track_positive[i], vector_Track_negative[j], particleMass[Pion], particleMass[Proton]));
                        insideprocessing.Fill("deltaT0ppi", DeltaT0(vector_Track_positive[i], vector_Track_negative[j], particleMass[Proton], particleMass[Pion]));
                        double t0cutoff = 0.6;
                        bool pipiPair = abs(DeltaT0(vector_Track_positive[i], vector_Track_negative[j], particleMass[Pion], particleMass[Pion]))<t0cutoff;
                        bool KpiPair = abs(DeltaT0(vector_Track_positive[i], vector_Track_negative[j], particleMass[Pion], particleMass[Kaon]))<t0cutoff or abs(DeltaT0(vector_Track_positive[i], vector_Track_negative[j], particleMass[Kaon], particleMass[Pion]))<t0cutoff;
                        bool KKPair = abs(DeltaT0(vector_Track_positive[i], vector_Track_negative[j], particleMass[Kaon], particleMass[Kaon]))<t0cutoff;
                        bool ppiPair = abs(DeltaT0(vector_Track_positive[i], vector_Track_negative[j], particleMass[Pion], particleMass[Proton]))<t0cutoff or abs(DeltaT0(vector_Track_positive[i], vector_Track_negative[j], particleMass[Proton], particleMass[Pion]))<t0cutoff;
                        int pairchoice = ppiPair*8+pipiPair*4+KpiPair*2+KKPair;
                        if(pairchoice==5){
                            insideprocessing.Fill("M2TOFpipiKK", M2TOF(vector_Track_positive[i], vector_Track_negative[j]));
                            double chi2pion = pow(vector_Track_positive[i]->getNSigmasTPCPion(), 2)+pow(vector_Track_negative[j]->getNSigmasTPCPion(), 2);
                            double chi2kaon = pow(vector_Track_positive[i]->getNSigmasTPCKaon(), 2)+pow(vector_Track_negative[j]->getNSigmasTPCKaon(), 2);
                            insideprocessing.Fill("chi2pipiKK", sqrt(chi2pion), sqrt(chi2kaon));
                        }
                        insideprocessing.Fill("choice", pairchoice);
                        if(pairchoice==8){
                            if(abs(DeltaT0(vector_Track_positive[i], vector_Track_negative[j], particleMass[Pion], particleMass[Proton]))<abs(DeltaT0(vector_Track_positive[i], vector_Track_negative[j], particleMass[Proton], particleMass[Pion]))){
                                vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                                vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Proton]);
                            } else{
                                vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Proton]);
                                vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Pion]);
                            }
                            mass = (positive_track+negative_track).M();
                            insideprocessing.Fill("MppiNarrowChoice", mass);
                        } else if(pairchoice==4){
                            vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                            vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Pion]);
                            mass = (positive_track+negative_track).M();
                            insideprocessing.Fill("MpipiChoice", mass);
                        } else if(pairchoice==2){
                            if(abs(DeltaT0(vector_Track_positive[i], vector_Track_negative[j], particleMass[Pion], particleMass[Kaon]))<abs(DeltaT0(vector_Track_positive[i], vector_Track_negative[j], particleMass[Kaon], particleMass[Pion]))){
                                vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                                vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Kaon]);
                            } else{
                                vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Kaon]);
                                vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Pion]);
                            }
                            mass = (positive_track+negative_track).M();
                            insideprocessing.Fill("MKpiNarrowChoice", mass);
                        } else if(pairchoice==1){
                            vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Kaon]);
                            vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Kaon]);
                            mass = (positive_track+negative_track).M();
                            insideprocessing.Fill("MKKNarrowChoice", mass);
                        } else if(pairchoice==0){
                            vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                            vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Pion]);
                            mass = (positive_track+negative_track).M();
                            insideprocessing.Fill("MpipiNoChoice", mass);
                        }
                    }

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

            vector_t0_all_pions.clear();
            for(long unsigned int i = 0; i<vector_Track_positive.size(); i++){
                if(vector_Track_positive[i]->getTofPathLength()<=0 or vector_Track_positive[i]->getTofTime()<=0){
                    vector_t0_all_pions.push_back(1e6);
                } else{
                    vector_t0_all_pions.push_back(t0(vector_Track_positive[i], particleMass[Pion]));
                }
            }
            for(long unsigned int j = 0; j<vector_Track_negative.size(); j++){
                if(vector_Track_negative[j]->getTofPathLength()<=0 or vector_Track_negative[j]->getTofTime()<=0){
                    vector_t0_all_pions.push_back(1e6);
                } else{
                    vector_t0_all_pions.push_back(t0(vector_Track_negative[j], particleMass[Pion]));
                }
            }
            min_t0 = *min_element(vector_t0_all_pions.begin(), vector_t0_all_pions.end());
            for(size_t i = 0; i<vector_t0_all_pions.size(); i++){
                insideprocessing.Fill("t0piTest", vector_t0_all_pions[i]-min_t0);
            }
            if(vector_t0_all_pions.size()>=2){
                sort(vector_t0_all_pions.begin(), vector_t0_all_pions.end());
                insideprocessing.Fill("t0piDifference", vector_t0_all_pions[1]-vector_t0_all_pions[0]);
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