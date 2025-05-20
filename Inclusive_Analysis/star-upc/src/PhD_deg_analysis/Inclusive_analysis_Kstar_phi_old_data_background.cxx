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
std::pair<int, int> extractExtendedParticlesNumbersFromPair(std::string pairName);

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
    //mixing TOF between events proved to be a failure
    //mass histograms with TOF first and mixing event pairs after
    outsideprocessing.AddHistogram(TH1D("MKpiChi2bcgTOF", ";m_{K^{+}#pi^{-}} [GeV];Number of pairs", 50, 0.5, 2.0));
    outsideprocessing.AddHistogram(TH1D("MpiKChi2bcgTOF", ";m_{#pi^{+}K^{-}} [GeV];Number of pairs", 50, 0.5, 2.0));
    outsideprocessing.AddHistogram(TH1D("MppiChi2bcgTOF", ";m_{p^{+}#pi^{-}} [GeV];Number of pairs", 50, 1.0, 2.5));
    outsideprocessing.AddHistogram(TH1D("MpipChi2bcgTOF", ";m_{#pi^{+}p^{-}} [GeV];Number of pairs", 50, 1.0, 2.5));
    outsideprocessing.AddHistogram(TH1D("MKKChi2bcgTOF", ";m_{K^{+}K^{-}} [GeV];Number of pairs", 50, 0.9, 2.4));
    outsideprocessing.AddHistogram(TH1D("MpipiChi2bcgTOF", ";m_{#pi^{+}#pi^{-}} [GeV];Number of pairs", 60, 0.2, 1.4));
    outsideprocessing.AddHistogram(TH1D("MppChi2bcgTOF", ";m_{p^{+}p^{-}} [GeV];Number of pairs", 50, 1.5, 3.5));
    //adding mass histograms grouped by category
    std::vector<std::string> pairTab = { "Kpi", "piK", "ppi", "pip", "KK", "pipi", "pp" };
    getCategoryHistograms(outsideprocessing, pairTab, "bcgTOF");
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

    //processing
    //defining TreeProcessor
    ROOT::TTreeProcessorMT TreeProc(*upcChain, nthreads);

    //other things
    int eventsProcessed = 0;
    int previousEventsMemorised = 15;

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
        std::vector<StUPCTrack*> vector_Track_positive = {};
        std::vector<StUPCTrack*> vector_Track_negative = {};
        map<string, std::vector<StUPCTrack*>> vector_Track_positive_TOF = {};
        map<string, std::vector<StUPCTrack*>> vector_Track_negative_TOF = {};
        std::vector<map<string, std::vector<StUPCTrack*>>> previous_vector_Tracks_positive_TOF = {};
        std::vector<map<string, std::vector<StUPCTrack*>>> previous_vector_Tracks_negative_TOF = {};
        for(int i = 0; i<previousEventsMemorised; i++){
            previous_vector_Tracks_positive_TOF.emplace_back();
            previous_vector_Tracks_negative_TOF.emplace_back();
        }

        StUPCTrack* tempTrack;
        TLorentzVector positive_track;
        TLorentzVector negative_track;
        double mass, pt, eta, phi;
        map<string, double> chi2Map;
        string tempPairName;
        TVector3 tempMomentum;
        string pairNames[] = { "K_pi", "pi_K", "p_pi", "pi_p", "K_K", "pi_pi", "p_p" };

        //actual loop
        while(myReader.Next()){
            //in a TTree, it *would* be constant, in TChain however not necessarily
            tempUPCpointer = StUPCEventInstance.Get();
            // tempRPpointer = StRPEventInstance.Get();

            //cleaning the loop
            vector_Track_positive.clear();
            vector_Track_negative.clear();
            for(auto&& name:pairNames){
                //removing previous pointers without removing their memory, because it is already invalidated
                vector_Track_positive_TOF[name].clear();
                vector_Track_negative_TOF[name].clear();
            }

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
                } else{
                    vector_Track_negative.push_back(tempTrack);
                }
            }

            //filling a chi2 map with keys for all the possibilities
            for(auto const& imap:sigmaMap){
                chi2Map.insert({ imap.first, 0. });
            }

            //loop for particles for TOF
            for(long unsigned int i = 0; i<vector_Track_positive.size(); i++){
                for(long unsigned int j = 0; j<vector_Track_negative.size(); j++){
                    //chi2
                    for(size_t pos = 0; pos<nParticlesExtended; pos++){
                        for(size_t neg = 0; neg<nParticlesExtended; neg++){
                            tempPairName = particleNicks[pos]+"_"+particleNicks[neg];
                            chi2Map[tempPairName] = getChi2(vector_Track_positive[i], vector_Track_negative[j], pos, neg, sigmaMap[tempPairName]);
                        }
                    }
                    //filling 
                    for(auto const& name:pairNames){
                        if(almostAllChi2(chi2Map, name, 9)){
                            vector_Track_positive_TOF[name].push_back(vector_Track_positive[i]);
                            vector_Track_negative_TOF[name].push_back(vector_Track_negative[j]);
                        }
                    }
                }
            }
            //mixing old and new tracks
            for(int prev_event = 0; prev_event<previousEventsMemorised; prev_event++){
                //from oldest to newest
                for(auto const& name:pairNames){
                    //previous + , current -
                    for(long unsigned int i = 0; i<previous_vector_Tracks_positive_TOF[prev_event][name].size(); i++){
                        for(long unsigned int j = 0; j<vector_Track_negative_TOF[name].size(); j++){
                            auto pair = extractExtendedParticlesNumbersFromPair(name);
                            //replaces $ with pair name without "_"
                            string tempHistName = "M$Chi2bcgTOF";
                            tempPairName = name;
                            tempPairName = tempPairName.erase(find(tempPairName.begin(), tempPairName.end(), '_')-tempPairName.begin(), 1);
                            tempHistName.replace(find(tempHistName.begin(), tempHistName.end(), '$')-tempHistName.begin(), 1, tempPairName);
                            //calculates mass to fill
                            previous_vector_Tracks_positive_TOF[prev_event][name][i]->getLorentzVector(positive_track, particleMassExtended[pair.first]);
                            vector_Track_negative_TOF[name][j]->getLorentzVector(negative_track, particleMassExtended[pair.second]);
                            mass = (positive_track+negative_track).M();
                            eta = (positive_track+negative_track).Eta();
                            pt = (positive_track+negative_track).Pt();
                            insideprocessing.Fill(tempHistName.c_str(), mass);
                            insideprocessing.Fill((tempHistName+"eta").c_str(), mass, eta);
                            insideprocessing.Fill((tempHistName+"pT").c_str(), mass, pt);
                        }
                    }
                    //current + , previous -
                    for(long unsigned int i = 0; i<vector_Track_positive_TOF[name].size(); i++){
                        for(long unsigned int j = 0; j<previous_vector_Tracks_negative_TOF[prev_event][name].size(); j++){
                            auto pair = extractExtendedParticlesNumbersFromPair(name);
                            //replaces $ with pair name without "_"
                            string tempHistName = "M$Chi2bcgTOF";
                            tempPairName = name;
                            tempPairName = tempPairName.erase(find(tempPairName.begin(), tempPairName.end(), '_')-tempPairName.begin(), 1);
                            tempHistName.replace(find(tempHistName.begin(), tempHistName.end(), '$')-tempHistName.begin(), 1, tempPairName);
                            //calculates mass to fill
                            vector_Track_positive_TOF[name][i]->getLorentzVector(positive_track, particleMassExtended[pair.first]);
                            previous_vector_Tracks_negative_TOF[prev_event][name][j]->getLorentzVector(negative_track, particleMassExtended[pair.second]);
                            mass = (positive_track+negative_track).M();
                            eta = (positive_track+negative_track).Eta();
                            pt = (positive_track+negative_track).Pt();
                            insideprocessing.Fill(tempHistName.c_str(), mass);
                            insideprocessing.Fill((tempHistName+"eta").c_str(), mass, eta);
                            insideprocessing.Fill((tempHistName+"pT").c_str(), mass, pt);
                        }
                    }
                }
            }

            //moving "current" pairs to "previous" storage if non-empty
            //and filling the current storage
            //positive
            for(auto&& name:pairNames){
                //if there is no good tracks in current event, there is no need to delete the old ones
                if(vector_Track_positive_TOF[name].size()==0){
                    continue;
                }
                //removing pointers from the oldest one
                for(auto&& i:previous_vector_Tracks_positive_TOF[0][name]){
                    delete i;
                }
                previous_vector_Tracks_positive_TOF[0][name].clear();
                //reverse bunny-hopping from the last one to the second-newest one
                //prev_event points to the one written to 
                for(int prev_event = 0; prev_event<previousEventsMemorised-1; prev_event++){
                    //copying empty track & setting needed values
                    for(long unsigned int i = 0; i<previous_vector_Tracks_positive_TOF[prev_event+1][name].size(); i++){
                        previous_vector_Tracks_positive_TOF[prev_event][name].push_back(new StUPCTrack());
                        previous_vector_Tracks_positive_TOF[prev_event+1][name][i]->getPtEtaPhi(pt, eta, phi);
                        previous_vector_Tracks_positive_TOF[prev_event][name].back()->setPtEtaPhi(pt, eta, phi);
                        previous_vector_Tracks_positive_TOF[prev_event][name].back()->setTofPathLength(previous_vector_Tracks_positive_TOF[prev_event+1][name][i]->getTofPathLength());
                        previous_vector_Tracks_positive_TOF[prev_event][name].back()->setTofTime(previous_vector_Tracks_positive_TOF[prev_event+1][name][i]->getTofTime());
                        for(size_t part = 0; part<nParticlesExtended; part++){
                            previous_vector_Tracks_positive_TOF[prev_event][name].back()->setNSigmasTPC(static_cast<StUPCTrack::Part>(part), previous_vector_Tracks_positive_TOF[prev_event+1][name][i]->getNSigmasTPC(static_cast<StUPCTrack::Part>(part)));
                        }
                    }
                    //clearing the currently moved one
                    for(auto&& i:previous_vector_Tracks_positive_TOF[prev_event+1][name]){
                        delete i;
                    }
                    previous_vector_Tracks_positive_TOF[prev_event+1][name].clear();
                }
                //moving the last one back with only the important parts filled
                for(long unsigned int i = 0; i<vector_Track_positive_TOF[name].size(); i++){
                    //copying empty track & setting needed values
                    previous_vector_Tracks_positive_TOF[previousEventsMemorised-1][name].push_back(new StUPCTrack());
                    vector_Track_positive_TOF[name][i]->getPtEtaPhi(pt, eta, phi);
                    previous_vector_Tracks_positive_TOF[previousEventsMemorised-1][name].back()->setPtEtaPhi(pt, eta, phi);
                    previous_vector_Tracks_positive_TOF[previousEventsMemorised-1][name].back()->setTofPathLength(vector_Track_positive_TOF[name][i]->getTofPathLength());
                    previous_vector_Tracks_positive_TOF[previousEventsMemorised-1][name].back()->setTofTime(vector_Track_positive_TOF[name][i]->getTofTime());
                    for(size_t part = 0; part<nParticlesExtended; part++){
                        previous_vector_Tracks_positive_TOF[previousEventsMemorised-1][name].back()->setNSigmasTPC(static_cast<StUPCTrack::Part>(part), vector_Track_positive_TOF[name][i]->getNSigmasTPC(static_cast<StUPCTrack::Part>(part)));
                    }
                }
            }
            //negative
            for(auto&& name:pairNames){
                //if there is no good tracks in current event, there is no need to delete the old ones
                if(vector_Track_negative_TOF[name].size()==0){
                    continue;
                }
                //removing pointers from the oldest one
                for(auto&& i:previous_vector_Tracks_negative_TOF[0][name]){
                    delete i;
                }
                previous_vector_Tracks_negative_TOF[0][name].clear();
                //reverse bunny-hopping from the last one to the second-newest one
                //prev_event points to the one written to 
                for(int prev_event = 0; prev_event<previousEventsMemorised-1; prev_event++){
                    for(long unsigned int i = 0; i<previous_vector_Tracks_negative_TOF[prev_event+1][name].size(); i++){
                        //copying empty track & setting needed values
                        previous_vector_Tracks_negative_TOF[prev_event][name].push_back(new StUPCTrack());
                        previous_vector_Tracks_negative_TOF[prev_event+1][name][i]->getPtEtaPhi(pt, eta, phi);
                        previous_vector_Tracks_negative_TOF[prev_event][name].back()->setPtEtaPhi(pt, eta, phi);
                        previous_vector_Tracks_negative_TOF[prev_event][name].back()->setTofPathLength(previous_vector_Tracks_negative_TOF[prev_event+1][name][i]->getTofPathLength());
                        previous_vector_Tracks_negative_TOF[prev_event][name].back()->setTofTime(previous_vector_Tracks_negative_TOF[prev_event+1][name][i]->getTofTime());
                        for(size_t part = 0; part<nParticlesExtended; part++){
                            previous_vector_Tracks_negative_TOF[prev_event][name].back()->setNSigmasTPC(static_cast<StUPCTrack::Part>(part), previous_vector_Tracks_negative_TOF[prev_event+1][name][i]->getNSigmasTPC(static_cast<StUPCTrack::Part>(part)));
                        }
                    }
                    //clearing the currently moved one
                    for(auto&& i:previous_vector_Tracks_negative_TOF[prev_event+1][name]){
                        delete i;
                    }
                    previous_vector_Tracks_negative_TOF[prev_event+1][name].clear();
                }
                //moving the last one back with only the important parts filled
                for(long unsigned int i = 0; i<vector_Track_negative_TOF[name].size(); i++){
                    //copying empty track
                    previous_vector_Tracks_negative_TOF[previousEventsMemorised-1][name].push_back(new StUPCTrack());
                    //setting needed values
                    vector_Track_negative_TOF[name][i]->getPtEtaPhi(pt, eta, phi);
                    previous_vector_Tracks_negative_TOF[previousEventsMemorised-1][name].back()->setPtEtaPhi(pt, eta, phi);
                    previous_vector_Tracks_negative_TOF[previousEventsMemorised-1][name].back()->setTofPathLength(vector_Track_negative_TOF[name][i]->getTofPathLength());
                    previous_vector_Tracks_negative_TOF[previousEventsMemorised-1][name].back()->setTofTime(vector_Track_negative_TOF[name][i]->getTofTime());
                    for(size_t part = 0; part<nParticlesExtended; part++){
                        previous_vector_Tracks_negative_TOF[previousEventsMemorised-1][name].back()->setNSigmasTPC(static_cast<StUPCTrack::Part>(part), vector_Track_negative_TOF[name][i]->getNSigmasTPC(static_cast<StUPCTrack::Part>(part)));
                    }
                }
            }


            //event loop finish
        }

        //lambda finish
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

std::pair<int, int> extractExtendedParticlesNumbersFromPair(std::string pairName){
    int first = find(particleNicks, particleNicks+nParticlesExtended, pairName.substr(0, find(pairName.begin(), pairName.end(), '_')-pairName.begin()))-particleNicks;
    int second = find(particleNicks, particleNicks+nParticlesExtended, pairName.substr(find(pairName.begin(), pairName.end(), '_')-pairName.begin()+1))-particleNicks;
    if(first>3 or second>3){
        printf("%s %d %d\n", pairName.c_str(), first, second);
        throw;
    }
    return std::pair<int, int>(first, second);
}