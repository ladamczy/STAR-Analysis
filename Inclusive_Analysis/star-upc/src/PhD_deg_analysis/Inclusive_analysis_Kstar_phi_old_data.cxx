//cpp headers
#include <map>

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
    outsideprocessing.AddHistogram(TH2D("chi2pipivsKpi1", ";#pi^{+}#pi^{-};K^{+}#pi^{-}", 100, 0, 100, 100, 0, 100));
    outsideprocessing.AddHistogram(TH1D("MKpiChi2Narrow", ";m_{K^{#pm}#pi^{#mp}} [GeV];Number of pairs", 50, 0.8, 1.0));
    outsideprocessing.AddHistogram(TH1D("MKpiChi2Wide", ";m_{K^{#pm}#pi^{#mp}} [GeV];Number of pairs", 150, 0.0, 1.5));
    outsideprocessing.AddHistogram(TH1D("MKpiChi2WideTest3", ";m_{K^{#pm}#pi^{#mp}} [GeV];Number of pairs", 150, 0.0, 1.5));
    outsideprocessing.AddHistogram(TH1D("MKpiChi2WideTest6", ";m_{K^{#pm}#pi^{#mp}} [GeV];Number of pairs", 150, 0.0, 1.5));
    outsideprocessing.AddHistogram(TH1D("MKpiChi2WideTest9", ";m_{K^{#pm}#pi^{#mp}} [GeV];Number of pairs", 150, 0.0, 1.5));
    outsideprocessing.AddHistogram(TH1D("MKpiChi2WideTest12", ";m_{K^{#pm}#pi^{#mp}} [GeV];Number of pairs", 150, 0.0, 1.5));

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
    auto myFunction = [&](TTreeReader& myReader){
        //getting values from TChain, in-loop histogram initialization
        TTreeReaderValue<StUPCEvent> StUPCEventInstance(myReader, "mUPCEvent");
        // TTreeReaderValue<StRPEvent> StRPEventInstance(myReader, "mRPEvent");
        // TTreeReaderValue<StRPEvent> StRPEventInstance(myReader, "mRPEvent");
        ProcessingInsideLoop insideprocessing;
        StUPCEvent* tempUPCpointer;
        // StRPEvent* tempRPpointer;
        // StRPEvent* tempRPpointer;
        insideprocessing.GetLocalHistograms(&outsideprocessing);

        //helpful variables
        std::vector<StUPCTrack*> vector_Track_positive;
        std::vector<StUPCTrack*> vector_Track_negative;
        int vector_Track_positive_dEdx;
        int vector_Track_negative_dEdx;
        int vector_bcg_track_1_dEdx;
        int vector_bcg_track_2_dEdx;
        StUPCTrack* tempTrack;
        TLorentzVector positive_track;
        TLorentzVector negative_track;
        TLorentzVector bcg_track_1;
        TLorentzVector bcg_track_2;
        double mass, chi1, chi2;

        //actual loop
        while(myReader.Next()){
            //in a TTree, it *would* be constant, in TChain however not necessarily
            tempUPCpointer = StUPCEventInstance.Get();
            // tempRPpointer = StRPEventInstance.Get();
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

            //loop through identified particles
            for(long unsigned int i = 0; i<vector_Track_positive.size(); i++){
                for(long unsigned int j = 0; j<vector_Track_negative.size(); j++){
                    if(vector_Track_positive[i]->getNhitsDEdx()<15 or vector_Track_negative[j]->getNhitsDEdx()<15){
                        continue;
                    }
                    //chi2
                    vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Kaon]);
                    vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Pion]);
                    mass = (positive_track+negative_track).M();
                    chi1 = pow(vector_Track_positive[i]->getNSigmasTPCPion(), 2)+pow(vector_Track_negative[j]->getNSigmasTPCPion(), 2)+pow(DeltaT0(vector_Track_positive[i], vector_Track_negative[j], particleMass[Pion], particleMass[Pion])/0.131243, 2);
                    chi2 = pow(vector_Track_positive[i]->getNSigmasTPCKaon(), 2)+pow(vector_Track_negative[j]->getNSigmasTPCPion(), 2)+pow(DeltaT0(vector_Track_positive[i], vector_Track_negative[j], particleMass[Kaon], particleMass[Pion])/0.125116, 2);
                    insideprocessing.Fill("chi2pipivsKpi1", chi1, chi2);
                    if(chi2<9&&chi1>9){
                        insideprocessing.Fill("MKpiChi2Narrow", mass);
                        insideprocessing.Fill("MKpiChi2Wide", mass);
                    }
                    if(chi2<3&&chi1>3){
                        insideprocessing.Fill("MKpiChi2WideTest3", mass);
                    }
                    if(chi2<6&&chi1>6){
                        insideprocessing.Fill("MKpiChi2WideTest6", mass);
                    }
                    if(chi2<9&&chi1>9){
                        insideprocessing.Fill("MKpiChi2WideTest9", mass);
                    }
                    if(chi2<12&&chi1>12){
                        insideprocessing.Fill("MKpiChi2WideTest12", mass);
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
    TFile* outputFileHist = TFile::Open(outfileName.c_str(), "recreate");
    outsideprocessing.SaveToFile(outputFileHist);
    outputFileHist->Close();

    return 0;
}