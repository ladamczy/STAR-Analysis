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
#include "StUPCV0.h"
#include "BeamPosition.h"

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
double RP_PV_z(StRPEvent* rpevent);
bool almostAllChi2(map<string, double> chi2Map, string exception, double limit);
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
    outsideprocessing.AddHistogram(TH1D("MKpiChi2", ";m_{K^{+}#pi^{-}} [GeV];Number of pairs", 50, 0.5, 2.0));
    outsideprocessing.AddHistogram(TH1D("MpiKChi2", ";m_{#pi^{+}K^{-}} [GeV];Number of pairs", 50, 0.5, 2.0));
    outsideprocessing.AddHistogram(TH1D("MppiChi2", ";m_{p^{+}#pi^{-}} [GeV];Number of pairs", 50, 1.0, 2.5));
    outsideprocessing.AddHistogram(TH1D("MpipChi2", ";m_{#pi^{+}p^{-}} [GeV];Number of pairs", 50, 1.0, 2.5));
    outsideprocessing.AddHistogram(TH1D("MKKChi2", ";m_{K^{+}K^{-}} [GeV];Number of pairs", 50, 0.9, 2.4));
    outsideprocessing.AddHistogram(TH1D("MpipiChi2", ";m_{#pi^{+}#pi^{-}} [GeV];Number of pairs", 60, 0.2, 1.4));
    outsideprocessing.AddHistogram(TH1D("MppChi2", ";m_{p^{+}p^{-}} [GeV];Number of pairs", 50, 1.5, 3.5));
    //adding mass histograms grouped by category
    getCategoryHistograms(outsideprocessing, pairTab);
    //adding TOF data test
    outsideprocessing.AddHistogram(TH2D("trackInfo", "TOF data;;", 0, 0, 0, 0, 0, 0));
    outsideprocessing.AddHistogram(TH1D("TOFlength", "TOF track length;TOF length [cm];Number of tracks", 100, 0, 1000));
    outsideprocessing.AddHistogram(TH1D("TOFtime", "TOF time;TOF time [ns];Number of tracks", 70, -1e4, 6e4));
    outsideprocessing.AddHistogram(TH1D("TOFtimePrecise", "TOF time;TOF time [ns];Number of tracks", 100, 0, 100));

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
        std::vector<string> TOFlength, TOFtime;
        std::vector<double> TOFlengthValue, TOFtimeValue;
        StUPCTrack* tempTrack;
        double beamValues[4];
        TVector3 vertexPrimary;
        StUPCV0* tempParticle;
        double mass, chi2pipi, chi2Kpi, eta, pT;
        map<string, double> chi2Map;
        bool isdEdxOk, isTOFOk, isWhicheverPrimary;
        string tempPairName;

        //filling pairInfo & trackInfo histograms in proper order
        insideprocessing.Fill("pairInfo", "OK", 0.0);
        insideprocessing.Fill("pairInfo", "TOF wrong, total", 0.0);
        insideprocessing.Fill("pairInfo", "TOF wrong, Primary", 0.0);
        insideprocessing.Fill("pairInfo", "TOF wrong, not Primary", 0.0);
        insideprocessing.Fill("pairInfo", "dEdx wrong", 0.0);
        insideprocessing.Fill("pairInfo", "Both wrong", 0.0);
        insideprocessing.Fill("trackInfo", "TOF length > 0", "TOF time > 0", 0.0);
        insideprocessing.Fill("trackInfo", "TOF length = 0", "TOF time = 0", 0.0);
        insideprocessing.Fill("trackInfo", "TOF length < 0", "TOF time < 0", 0.0);

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
            //exactly one vertex
            if(tempUPCpointer->getNumberOfVertices()!=1){
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

            //cuts & histogram filling
            //selecting tracks matching criteria:
            //new tracks (not doubled)
            //TOF
            //pt & eta 
            //Nhits

            //reusing good track counter
            nOfGoodTracks = 0;
            for(int i = 0; i<tempUPCpointer->getNumberOfTracks(); i++){
                tempTrack = tempUPCpointer->getTrack(i);
                //kV0 true
                if(!tempTrack->getFlag(StUPCTrack::kV0)){
                    continue;
                }
                if(!tempTrack->getFlag(StUPCTrack::kTof)){
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
                nOfGoodTracks++;
                //test if the track has TOF info
                //length
                if(tempTrack->getTofPathLength()>0){
                    TOFlength.push_back("TOF length > 0");
                } else if(tempTrack->getTofPathLength()==0){
                    TOFlength.push_back("TOF length = 0");
                } else if(tempTrack->getTofPathLength()<0){
                    TOFlength.push_back("TOF length < 0");
                }
                TOFlengthValue.push_back(tempTrack->getTofPathLength());
                //time
                if(tempTrack->getTofTime()>0){
                    TOFtime.push_back("TOF time > 0");
                } else if(tempTrack->getTofTime()==0){
                    TOFtime.push_back("TOF time = 0");
                } else if(tempTrack->getTofTime()<0){
                    TOFtime.push_back("TOF time < 0");
                }
                TOFtimeValue.push_back(tempTrack->getTofTime());
            }
            if(nOfGoodTracks>=2){
                //saving
                for(size_t i = 0; i<TOFlength.size(); i++){
                    insideprocessing.Fill("trackInfo", TOFlength[i].c_str(), TOFtime[i].c_str(), 1.0);
                }
                for(size_t i = 0; i<TOFlengthValue.size(); i++){
                    insideprocessing.Fill("TOFlength", TOFlengthValue[i]);
                    insideprocessing.Fill("TOFtime", TOFtimeValue[i]);
                    insideprocessing.Fill("TOFtimePrecise", TOFtimeValue[i]);
                }
            }
            TOFlength.clear();
            TOFtime.clear();
            TOFlengthValue.clear();
            TOFtimeValue.clear();

            //filling a chi2 map with keys for all the possibilities
            for(auto const& imap:sigmaMap){
                chi2Map.insert({ imap.first, 0. });
            }

            //stuff for proper StUPCV0 reconstruction
            beamValues[0] = tempUPCpointer->getBeamXPosition();
            beamValues[1] = tempUPCpointer->getBeamYPosition();
            beamValues[2] = tempUPCpointer->getBeamXSlope();
            beamValues[3] = tempUPCpointer->getBeamYSlope();
            // vertexPrimary = { beamValues[0]+RP_PV_z(tempRPpointer)*beamValues[2], beamValues[1]+RP_PV_z(tempRPpointer)*beamValues[3], RP_PV_z(tempRPpointer) };
            vertexPrimary = { tempUPCpointer->getVertex(0)->getPosX(), tempUPCpointer->getVertex(0)->getPosY(), tempUPCpointer->getVertex(0)->getPosZ() };

            //loop through identified particles
            for(long unsigned int i = 0; i<vector_Track_positive.size(); i++){
                for(long unsigned int j = 0; j<vector_Track_negative.size(); j++){
                    isdEdxOk = (vector_Track_positive[i]->getNhitsDEdx()>=15)&&(vector_Track_negative[j]->getNhitsDEdx()>=15);
                    isTOFOk = (vector_Track_positive[i]->getTofPathLength()>0)&&(vector_Track_positive[i]->getTofTime()>0)&&(vector_Track_negative[j]->getTofPathLength()>0)&&(vector_Track_negative[j]->getTofTime()>0);
                    isWhicheverPrimary = vector_Track_positive[i]->getFlag(StUPCTrack::kPrimary)||vector_Track_negative[j]->getFlag(StUPCTrack::kPrimary);
                    if(isdEdxOk&&isTOFOk){
                        insideprocessing.Fill("pairInfo", "OK", 1.0);
                    } else if(isdEdxOk&&!isTOFOk){
                        insideprocessing.Fill("pairInfo", "TOF wrong, total", 1.0);
                        if(isWhicheverPrimary){
                            insideprocessing.Fill("pairInfo", "TOF wrong, Primary", 1.0);
                        } else{
                            insideprocessing.Fill("pairInfo", "TOF wrong, not Primary", 1.0);
                        }
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
                        tempParticle = new StUPCV0(vector_Track_positive[i], vector_Track_negative[j], particleMass[Kaon], particleMass[Pion], 1, 1, vertexPrimary, beamValues, tempUPCpointer->getMagneticField(), false);
                        mass = tempParticle->m();
                        eta = tempParticle->eta();
                        pT = tempParticle->pt();
                        insideprocessing.Fill("MKpiChi2", mass);
                        insideprocessing.Fill("MKpiChi2eta", mass, eta);
                        insideprocessing.Fill("MKpiChi2pT", mass, pT);
                        delete tempParticle;
                    }
                    if(almostAllChi2(chi2Map, "pi_K", 9)){
                        tempParticle = new StUPCV0(vector_Track_positive[i], vector_Track_negative[j], particleMass[Pion], particleMass[Kaon], 1, 1, vertexPrimary, beamValues, tempUPCpointer->getMagneticField(), false);
                        mass = tempParticle->m();
                        eta = tempParticle->eta();
                        pT = tempParticle->pt();
                        insideprocessing.Fill("MpiKChi2", mass);
                        insideprocessing.Fill("MpiKChi2eta", mass, eta);
                        insideprocessing.Fill("MpiKChi2pT", mass, pT);
                        delete tempParticle;
                    }
                    if(almostAllChi2(chi2Map, "p_pi", 9)){
                        tempParticle = new StUPCV0(vector_Track_positive[i], vector_Track_negative[j], particleMass[Proton], particleMass[Pion], 1, 1, vertexPrimary, beamValues, tempUPCpointer->getMagneticField(), false);
                        mass = tempParticle->m();
                        eta = tempParticle->eta();
                        pT = tempParticle->pt();
                        insideprocessing.Fill("MppiChi2", mass);
                        insideprocessing.Fill("MppiChi2eta", mass, eta);
                        insideprocessing.Fill("MppiChi2pT", mass, pT);
                        delete tempParticle;
                    }
                    if(almostAllChi2(chi2Map, "pi_p", 9)){
                        tempParticle = new StUPCV0(vector_Track_positive[i], vector_Track_negative[j], particleMass[Pion], particleMass[Proton], 1, 1, vertexPrimary, beamValues, tempUPCpointer->getMagneticField(), false);
                        mass = tempParticle->m();
                        eta = tempParticle->eta();
                        pT = tempParticle->pt();
                        insideprocessing.Fill("MpipChi2", mass);
                        insideprocessing.Fill("MpipChi2eta", mass, eta);
                        insideprocessing.Fill("MpipChi2pT", mass, pT);
                        delete tempParticle;
                    }
                    if(almostAllChi2(chi2Map, "K_K", 9)){
                        tempParticle = new StUPCV0(vector_Track_positive[i], vector_Track_negative[j], particleMass[Kaon], particleMass[Kaon], 1, 1, vertexPrimary, beamValues, tempUPCpointer->getMagneticField(), false);
                        mass = tempParticle->m();
                        eta = tempParticle->eta();
                        pT = tempParticle->pt();
                        insideprocessing.Fill("MKKChi2", mass);
                        insideprocessing.Fill("MKKChi2eta", mass, eta);
                        insideprocessing.Fill("MKKChi2pT", mass, pT);
                        delete tempParticle;
                    }
                    if(almostAllChi2(chi2Map, "pi_pi", 9)){
                        tempParticle = new StUPCV0(vector_Track_positive[i], vector_Track_negative[j], particleMass[Pion], particleMass[Pion], 1, 1, vertexPrimary, beamValues, tempUPCpointer->getMagneticField(), false);
                        mass = tempParticle->m();
                        eta = tempParticle->eta();
                        pT = tempParticle->pt();
                        insideprocessing.Fill("MpipiChi2", mass);
                        insideprocessing.Fill("MpipiChi2eta", mass, eta);
                        insideprocessing.Fill("MpipiChi2pT", mass, pT);
                        delete tempParticle;
                    }
                    if(almostAllChi2(chi2Map, "p_p", 9)){
                        tempParticle = new StUPCV0(vector_Track_positive[i], vector_Track_negative[j], particleMass[Proton], particleMass[Proton], 1, 1, vertexPrimary, beamValues, tempUPCpointer->getMagneticField(), false);
                        mass = tempParticle->m();
                        eta = tempParticle->eta();
                        pT = tempParticle->pt();
                        insideprocessing.Fill("MppChi2", mass);
                        insideprocessing.Fill("MppChi2eta", mass, eta);
                        insideprocessing.Fill("MppChi2pT", mass, pT);
                        delete tempParticle;
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
    outsideprocessing.GetPointerAfterMerge1D("pairInfo")->LabelsDeflate();
    outsideprocessing.GetPointerAfterMerge2D("trackInfo")->LabelsDeflate("X");
    outsideprocessing.GetPointerAfterMerge2D("trackInfo")->LabelsDeflate("Y");

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

double RP_PV_z(StRPEvent* rpevent){
    StUPCRpsTrack* p_plus = rpevent->getTrack(0);
    p_plus->setEvent(rpevent);
    StUPCRpsTrack* p_minus;
    if(p_plus->getTrackPoint(0)->z()>0){
        p_plus = rpevent->getTrack(0);
        p_minus = rpevent->getTrack(1);
        p_plus->setEvent(rpevent);
        p_minus->setEvent(rpevent);
    } else{
        p_plus = rpevent->getTrack(1);
        p_minus = rpevent->getTrack(0);
        p_plus->setEvent(rpevent);
        p_minus->setEvent(rpevent);
    }
    double v1, v2, d1, d2, t1, t2;
    TLorentzVector temp;
    t1 = p_plus->time();
    t2 = p_minus->time();
    // printf("%lf %lf\n", t1*1e9, t2*1e9);
    d1 = 0.5*(p_plus->getTrackPoint(0)->z()+p_plus->getTrackPoint(1)->z());
    d2 = -0.5*(p_minus->getTrackPoint(0)->z()+p_minus->getTrackPoint(1)->z());
    temp.SetVectM(p_plus->pVec(), particleMass[Proton]);
    v1 = temp.BoostVector().Z()*299792458.0*100; //na cm/s
    temp.SetVectM(p_minus->pVec(), particleMass[Proton]);
    v2 = -temp.BoostVector().Z()*299792458.0*100; //na cm/s

    return (t1-t2-(d1/v1-d2/v2))/(1/v1+1/v2);
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

void getTPCSector(StUPCEvent* event, StUPCTrack* track, double& eta, double& phi){
    TVector3 momentum;
    track->getMomentum(momentum);
    StPicoPhysicalHelix helix(momentum, track->getOrigin(), event->getMagneticField()*kilogauss, track->getCharge());
    helix.moveOrigin(helix.pathLength(200.*centimeter).first);
    eta = helix.origin().Eta();
    phi = helix.origin().Phi();
}