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
bool almostAllChi2(map<string, double> chi2Map, string exception, double limit);
double getMass(StUPCTrack* assumed, StUPCTrack* calculated, double massAssumed);
void getTPCSector(StUPCEvent* event, StUPCTrack* track, double& eta, double& phi);

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
    outsideprocessing.AddHistogram(TH1D("pairInfo", "", 1, 0, 1));
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
    //chi2 tests
    outsideprocessing.AddHistogram(TH2D("chi2pipivsKpi1", ";#pi^{+}#pi^{-};K^{+}#pi^{-}", 100, 0, 100, 100, 0, 100));
    //mass histograms
    std::vector<std::string> pairTab = { "Kpi", "piK", "ppi", "pip", "KK", "pipi", "pp" };
    outsideprocessing.AddHistogram(TH1D("MKpiChi2", ";m_{K^{+}#pi^{-}} [GeV];Number of pairs", 200, 0.5, 2.0));
    outsideprocessing.AddHistogram(TH1D("MpiKChi2", ";m_{#pi^{+}K^{-}} [GeV];Number of pairs", 200, 0.5, 2.0));
    outsideprocessing.AddHistogram(TH1D("MppiChi2", ";m_{p^{+}#pi^{-}} [GeV];Number of pairs", 500, 1.0, 2.5));
    outsideprocessing.AddHistogram(TH1D("MpipChi2", ";m_{#pi^{+}p^{-}} [GeV];Number of pairs", 500, 1.0, 2.5));
    outsideprocessing.AddHistogram(TH1D("MKKChi2", ";m_{K^{+}K^{-}} [GeV];Number of pairs", 500, 0.9, 2.4));
    outsideprocessing.AddHistogram(TH1D("MpipiChi2", ";m_{#pi^{+}#pi^{-}} [GeV];Number of pairs", 600, 0.2, 1.4));
    outsideprocessing.AddHistogram(TH1D("MppChi2", ";m_{p^{+}p^{-}} [GeV];Number of pairs", 500, 1.5, 3.5));
    //closer histograms
    outsideprocessing.AddHistogram(TH1D("MKKChi2Close", ";m_{K^{+}K^{-}} [GeV];Number of pairs", 50, 0.99, 1.05));
    outsideprocessing.AddHistogram(TH1D("MKpiChi2Close", ";m_{K^{+}#pi^{-}} [GeV];Number of pairs", 50, 0.7, 1.1));
    outsideprocessing.AddHistogram(TH1D("MpiKChi2Close", ";m_{#pi^{+}K^{-}} [GeV];Number of pairs", 50, 0.7, 1.1));
    //adding mass histograms grouped by category
    getCategoryHistograms(outsideprocessing, pairTab);
    //eta-phi histograms
    outsideprocessing.AddHistogram(TH2D("EtaPhi_KpiChi2plus", ";#eta;#phi [rad]", 50, -0.7, 0.7, 50, -TMath::Pi(), TMath::Pi()));
    outsideprocessing.AddHistogram(TH2D("EtaPhi_KpiChi2minus", ";#eta;#phi [rad]", 50, -0.7, 0.7, 50, -TMath::Pi(), TMath::Pi()));
    outsideprocessing.AddHistogram(TH2D("EtaPhi_piKChi2plus", ";#eta;#phi [rad]", 50, -0.7, 0.7, 50, -TMath::Pi(), TMath::Pi()));
    outsideprocessing.AddHistogram(TH2D("EtaPhi_piKChi2minus", ";#eta;#phi [rad]", 50, -0.7, 0.7, 50, -TMath::Pi(), TMath::Pi()));
    outsideprocessing.AddHistogram(TH2D("EtaPhi_KKChi2plus", ";#eta;#phi [rad]", 50, -0.7, 0.7, 50, -TMath::Pi(), TMath::Pi()));
    outsideprocessing.AddHistogram(TH2D("EtaPhi_KKChi2minus", ";#eta;#phi [rad]", 50, -0.7, 0.7, 50, -TMath::Pi(), TMath::Pi()));
    //for tests, low resolution
    outsideprocessing.AddHistogram(TH2D("EtaPhi_KpiChi2plus_lowres", ";#eta;#phi [rad]", 14, -0.7, 0.7, 8, -TMath::Pi(), TMath::Pi()));
    outsideprocessing.AddHistogram(TH2D("EtaPhi_KpiChi2minus_lowres", ";#eta;#phi [rad]", 14, -0.7, 0.7, 8, -TMath::Pi(), TMath::Pi()));
    outsideprocessing.AddHistogram(TH2D("EtaPhi_piKChi2plus_lowres", ";#eta;#phi [rad]", 14, -0.7, 0.7, 8, -TMath::Pi(), TMath::Pi()));
    outsideprocessing.AddHistogram(TH2D("EtaPhi_piKChi2minus_lowres", ";#eta;#phi [rad]", 14, -0.7, 0.7, 8, -TMath::Pi(), TMath::Pi()));
    outsideprocessing.AddHistogram(TH2D("EtaPhi_KKChi2plus_lowres", ";#eta;#phi [rad]", 14, -0.7, 0.7, 8, -TMath::Pi(), TMath::Pi()));
    outsideprocessing.AddHistogram(TH2D("EtaPhi_KKChi2minus_lowres", ";#eta;#phi [rad]", 14, -0.7, 0.7, 8, -TMath::Pi(), TMath::Pi()));
    outsideprocessing.AddHistogram(TH2D("EtaPhi_pipiChi2plus_lowres", ";#eta;#phi [rad]", 14, -0.7, 0.7, 8, -TMath::Pi(), TMath::Pi()));
    outsideprocessing.AddHistogram(TH2D("EtaPhi_pipiChi2minus_lowres", ";#eta;#phi [rad]", 14, -0.7, 0.7, 8, -TMath::Pi(), TMath::Pi()));
    //tpc sector tests
    outsideprocessing.AddHistogram(TH3D("Sector_pipiChi2plus", ";fill;#eta;#phi [rad]", 1, 0, 1, 2, -0.7, 0.7, 12, -TMath::Pi(), TMath::Pi()));
    outsideprocessing.AddHistogram(TH3D("Sector_pipiChi2minus", ";fill;#eta;#phi [rad]", 1, 0, 1, 2, -0.7, 0.7, 12, -TMath::Pi(), TMath::Pi()));
    outsideprocessing.AddHistogram(TH1D("Sector_run_tracks_number", ";run;tracks", 1, 0, 1));

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
                if(tmptrk->getFlag(StUPCTrack::kTof)&&abs(tmptrk->getEta())<0.9&&tmptrk->getPt()>0.2){
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
                if(tempTrack->getPt()<=0.2 or abs(tempTrack->getEta())>=0.9){
                    continue;
                }
                if(tempTrack->getNhits()<=20){
                    continue;
                }
                if(tempTrack->getCharge()>0){
                    vector_Track_positive.push_back(tempTrack);
                    getTPCSector(tempUPCpointer, tempTrack, eta, phi);
                    insideprocessing.Fill("Sector_pipiChi2plus", to_string(tempUPCpointer->getRunNumber()).c_str(), eta, phi, 1.0);
                } else{
                    vector_Track_negative.push_back(tempTrack);
                    getTPCSector(tempUPCpointer, tempTrack, eta, phi);
                    insideprocessing.Fill("Sector_pipiChi2minus", to_string(tempUPCpointer->getRunNumber()).c_str(), eta, phi, 1.0);
                }
                nOfGoodTracks++;
            }
            insideprocessing.Fill("Sector_run_tracks_number", to_string(tempUPCpointer->getRunNumber()).c_str(), nOfGoodTracks);

            //filling a chi2 map with keys for all the possibilities
            for(auto const& imap:sigmaMap){
                chi2Map.insert({ imap.first, 0. });
            }

            //loop through identified particles
            for(long unsigned int i = 0; i<vector_Track_positive.size(); i++){
                for(long unsigned int j = 0; j<vector_Track_negative.size(); j++){
                    isdEdxOk = (vector_Track_positive[i]->getNhitsDEdx()>=15)&&(vector_Track_negative[j]->getNhitsDEdx()>=15);
                    isTOFOk = (vector_Track_positive[i]->getTofPathLength()>0)&&(vector_Track_positive[i]->getTofTime()>0)&&(vector_Track_negative[j]->getTofPathLength()>0)&&(vector_Track_negative[j]->getTofTime()>0);
                    if(isdEdxOk&&isTOFOk){
                        insideprocessing.Fill("pairInfo", "OK", 1.0);
                    } else if(isdEdxOk&&!isTOFOk){
                        insideprocessing.Fill("pairInfo", "TOF wrong", 1.0);
                        continue;
                    } else if(!isdEdxOk&&isTOFOk){
                        insideprocessing.Fill("pairInfo", "dEdx wrong", 1.0);
                        continue;
                    } else if(!isdEdxOk&&!isTOFOk){
                        insideprocessing.Fill("pairInfo", "Both wrong", 1.0);
                        continue;
                    }

                    //chi2
                    for(size_t pos = 0; pos<nParticlesExtended; pos++){
                        for(size_t neg = 0; neg<nParticlesExtended; neg++){
                            tempPairName = particleNicks[pos]+"_"+particleNicks[neg];
                            chi2Map[tempPairName] = getChi2(vector_Track_positive[i], vector_Track_negative[j], pos, neg, sigmaMap[tempPairName]);
                        }
                    }

                    //chi2 comparison
                    chi2pipi = chi2Map["pi_pi"];
                    chi2Kpi = chi2Map["K_pi"];
                    insideprocessing.Fill("chi2pipivsKpi1", chi2pipi, chi2Kpi);

                    //mass tests on different pairs
                    if(almostAllChi2(chi2Map, "K_pi", 9)){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Kaon]);
                        vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Pion]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        insideprocessing.Fill("MKpiChi2", mass);
                        insideprocessing.Fill("MKpiChi2Close", mass);
                        insideprocessing.Fill("MKpiChi2eta", mass, eta);
                        insideprocessing.Fill("MKpiChi2pT", mass, pT);
                        insideprocessing.Fill("EtaPhi_KpiChi2plus", vector_Track_positive[i]->getEta(), vector_Track_positive[i]->getPhi());
                        insideprocessing.Fill("EtaPhi_KpiChi2minus", vector_Track_negative[j]->getEta(), vector_Track_negative[j]->getPhi());
                        insideprocessing.Fill("EtaPhi_KpiChi2plus_lowres", vector_Track_positive[i]->getEta(), vector_Track_positive[i]->getPhi());
                        insideprocessing.Fill("EtaPhi_KpiChi2minus_lowres", vector_Track_negative[j]->getEta(), vector_Track_negative[j]->getPhi());
                    }
                    if(almostAllChi2(chi2Map, "pi_K", 9)){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                        vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Kaon]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        insideprocessing.Fill("MpiKChi2", mass);
                        insideprocessing.Fill("MpiKChi2Close", mass);
                        insideprocessing.Fill("MpiKChi2eta", mass, eta);
                        insideprocessing.Fill("MpiKChi2pT", mass, pT);
                        insideprocessing.Fill("EtaPhi_piKChi2plus", vector_Track_positive[i]->getEta(), vector_Track_positive[i]->getPhi());
                        insideprocessing.Fill("EtaPhi_piKChi2minus", vector_Track_negative[j]->getEta(), vector_Track_negative[j]->getPhi());
                        insideprocessing.Fill("EtaPhi_piKChi2plus_lowres", vector_Track_positive[i]->getEta(), vector_Track_positive[i]->getPhi());
                        insideprocessing.Fill("EtaPhi_piKChi2minus_lowres", vector_Track_negative[j]->getEta(), vector_Track_negative[j]->getPhi());
                    }
                    if(almostAllChi2(chi2Map, "p_pi", 9)){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Proton]);
                        vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Pion]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        insideprocessing.Fill("MppiChi2", mass);
                        insideprocessing.Fill("MppiChi2eta", mass, eta);
                        insideprocessing.Fill("MppiChi2pT", mass, pT);
                    }
                    if(almostAllChi2(chi2Map, "pi_p", 9)){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                        vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Proton]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        insideprocessing.Fill("MpipChi2", mass);
                        insideprocessing.Fill("MpipChi2eta", mass, eta);
                        insideprocessing.Fill("MpipChi2pT", mass, pT);
                    }
                    if(almostAllChi2(chi2Map, "K_K", 9)){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Kaon]);
                        vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Kaon]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        insideprocessing.Fill("MKKChi2", mass);
                        insideprocessing.Fill("MKKChi2Close", mass);
                        insideprocessing.Fill("MKKChi2eta", mass, eta);
                        insideprocessing.Fill("MKKChi2pT", mass, pT);
                        insideprocessing.Fill("EtaPhi_KKChi2plus", vector_Track_positive[i]->getEta(), vector_Track_positive[i]->getPhi());
                        insideprocessing.Fill("EtaPhi_KKChi2minus", vector_Track_negative[j]->getEta(), vector_Track_negative[j]->getPhi());
                        insideprocessing.Fill("EtaPhi_KKChi2plus_lowres", vector_Track_positive[i]->getEta(), vector_Track_positive[i]->getPhi());
                        insideprocessing.Fill("EtaPhi_KKChi2minus_lowres", vector_Track_negative[j]->getEta(), vector_Track_negative[j]->getPhi());
                    }
                    if(almostAllChi2(chi2Map, "pi_pi", 9)){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                        vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Pion]);
                        mass = (positive_track+negative_track).M();
                        eta = (positive_track+negative_track).Eta();
                        pT = (positive_track+negative_track).Pt();
                        insideprocessing.Fill("MpipiChi2", mass);
                        insideprocessing.Fill("MpipiChi2eta", mass, eta);
                        insideprocessing.Fill("MpipiChi2pT", mass, pT);
                        insideprocessing.Fill("EtaPhi_pipiChi2plus_lowres", vector_Track_positive[i]->getEta(), vector_Track_positive[i]->getPhi());
                        insideprocessing.Fill("EtaPhi_pipiChi2minus_lowres", vector_Track_negative[j]->getEta(), vector_Track_negative[j]->getPhi());
                    }
                    if(almostAllChi2(chi2Map, "p_p", 9)){
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

            //lambda finish
        }
        return 0;
    };

    TreeProc.Process(myFunction);

    //merging and tidying up
    outsideprocessing.Merge();
    outsideprocessing.GetPointerAfterMerge3D("Sector_pipiChi2plus")->LabelsDeflate();
    outsideprocessing.GetPointerAfterMerge3D("Sector_pipiChi2minus")->LabelsDeflate();
    outsideprocessing.GetPointerAfterMerge1D("Sector_run_tracks_number")->LabelsDeflate();

    //histogram of track numbers
    TH1D Sector_track_number_histogram("Sector_track_number_histogram", ";# of tracks;# of runs", 100, 0, 10000);
    for(int k = 1; k<=outsideprocessing.GetPointerAfterMerge1D("Sector_run_tracks_number")->GetXaxis()->GetNbins(); k++){
        double number_of_tracks = outsideprocessing.GetPointerAfterMerge1D("Sector_run_tracks_number")->GetBinContent(k);
        Sector_track_number_histogram.Fill(number_of_tracks);
    }

    //setting and filling the run test histogram
    TH1D Sector_filling_positive("Sector_filling_positive", ";fill;# of sectors", 100, 0, 1.0000000001);
    TH1D Sector_filling_negative("Sector_filling_negative", ";fill;# of sectors", 100, 0, 1.0000000001);
    TH1D Sector_filling_positive_bigger_runs("Sector_filling_positive_bigger_runs", ";fill;# of sectors", 100, 0, 1.0000000001);
    TH1D Sector_filling_negative_bigger_runs("Sector_filling_negative_bigger_runs", ";fill;# of sectors", 100, 0, 1.0000000001);
    TH1D Sector_max_tracks_with_empty_sector("Sector_max_tracks_with_empty_sector", ";max track number;runs", 100, 0, 100);
    int minimal_tracks_by_charge = 48; //so far chosen arbitrarily
    for(int k = 1; k<=outsideprocessing.GetPointerAfterMerge3D("Sector_pipiChi2plus")->GetXaxis()->GetNbins(); k++){
        outsideprocessing.GetPointerAfterMerge3D("Sector_pipiChi2plus")->GetXaxis()->SetRange(k, k);
        outsideprocessing.GetPointerAfterMerge3D("Sector_pipiChi2minus")->GetXaxis()->SetRange(k, k);
        TH1* sig_slice_positive = outsideprocessing.GetPointerAfterMerge3D("Sector_pipiChi2plus")->Project3D("yz");
        TH1* sig_slice_negative = outsideprocessing.GetPointerAfterMerge3D("Sector_pipiChi2minus")->Project3D("yz");
        double temp_max_positive = sig_slice_positive->GetMaximum();
        double temp_max_negative = sig_slice_negative->GetMaximum();
        if(temp_max_positive*temp_max_negative==0){
            continue;
        }
        sig_slice_positive->Scale(1./temp_max_positive);
        sig_slice_negative->Scale(1./temp_max_negative);
        bool runWasCheckedAndHasEmptySector = false;
        for(int xbin = 1; xbin<=sig_slice_positive->GetNbinsX(); xbin++){
            for(int ybin = 1; ybin<=sig_slice_positive->GetNbinsY(); ybin++){
                Sector_filling_positive.Fill(sig_slice_positive->GetBinContent(xbin, ybin));
                Sector_filling_negative.Fill(sig_slice_negative->GetBinContent(xbin, ybin));
                if(sig_slice_positive->Integral()*temp_max_positive>=minimal_tracks_by_charge)
                    Sector_filling_positive_bigger_runs.Fill(sig_slice_positive->GetBinContent(xbin, ybin));
                if(sig_slice_negative->Integral()*temp_max_negative>=minimal_tracks_by_charge)
                    Sector_filling_negative_bigger_runs.Fill(sig_slice_negative->GetBinContent(xbin, ybin));
                //checking if the run *has* empty sector
                if(sig_slice_positive->GetBinContent(xbin, ybin)==0.&&sig_slice_negative->GetBinContent(xbin, ybin)==0.&&!runWasCheckedAndHasEmptySector){
                    runWasCheckedAndHasEmptySector = true;
                    Sector_max_tracks_with_empty_sector.Fill(temp_max_positive+temp_max_negative);
                }
            }
        }
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
    Sector_track_number_histogram.Write();
    Sector_filling_positive.Write();
    Sector_filling_negative.Write();
    Sector_filling_positive_bigger_runs.Write();
    Sector_filling_negative_bigger_runs.Write();
    Sector_max_tracks_with_empty_sector.Write();

    outputFileHist->Close();

    return 0;
}

double getChi2(StUPCTrack* positive, StUPCTrack* negative, int positiveId, int negativeId, double sigmaT){
    double sigma1 = pow(positive->getNSigmasTPC(static_cast<StUPCTrack::Part>(positiveId)), 2);
    double sigma2 = pow(negative->getNSigmasTPC(static_cast<StUPCTrack::Part>(negativeId)), 2);
    double sigma3 = pow(DeltaT0(positive, negative, particleMassExtended[positiveId], particleMassExtended[negativeId])/sigmaT, 2);
    return sigma1+sigma2+sigma3;
}

bool almostAllChi2(map<string, double> chi2Map, string exception, double limit){
    //checks if almost all chi2 statistics stay above limit
    //except the exception, which needs to stay below
    //and some that just don't make sense
    for(auto&& i:chi2Map){
        if(i.first=="e_e"
            or i.first=="K_e"
            or i.first=="e_K"
            or i.first=="pi_e"
            or i.first=="e_pi"
            or i.first=="p_e"
            or i.first=="e_p"
            or i.first=="p_K"
            or i.first=="K_p"){
            continue;
        }
        if(i.first!=exception&&i.second<limit){
            return false;
        } else if(i.first==exception&&i.second>limit){
            return false;
        }
    }
    return true;
}

double getMass(StUPCTrack* assumed, StUPCTrack* calculated, double massAssumed){
    double intermediateSquareRoot = assumed->getTofPathLength()/calculated->getTofPathLength()-(assumed->getTofTime()-calculated->getTofTime())*100*0.299792458/calculated->getTofPathLength();
    TVector3 momentum;
    calculated->getMomentum(momentum);
    return momentum.Mag()*sqrt(intermediateSquareRoot*intermediateSquareRoot-1);
}

void getTPCSector(StUPCEvent* event, StUPCTrack* track, double& eta, double& phi){
    TVector3 momentum;
    track->getMomentum(momentum);
    StPicoPhysicalHelix helix(momentum, track->getOrigin(), event->getMagneticField()*kilogauss, track->getCharge());
    helix.moveOrigin(helix.pathLength(200.*centimeter).first);
    eta = helix.origin().Eta();
    phi = helix.origin().Phi();
}