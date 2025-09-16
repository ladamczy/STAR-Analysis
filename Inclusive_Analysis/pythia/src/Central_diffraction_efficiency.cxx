// Stdlib header file for input and output.
#include <iostream>

// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"

// ROOT, for histogramming.
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TF1.h"
#include "TEfficiency.h"
#include "THStack.h"
#include "TCanvas.h"

// My classes and stuff
#include "UsefulThingsPythia.h"
// #include "UsefulThings.h"

using namespace Pythia8;
using namespace std;

bool isInFiducial(TLorentzVector);
bool IsInXiElasticSpot(TLorentzVector, TLorentzVector);
bool IsInMomElasticSpot(TLorentzVector, TLorentzVector);
bool IsInXiFiducialNoBBCL(TLorentzVector, TLorentzVector, double);
bool IsInXiFiducialBBCL(TLorentzVector, TLorentzVector, double);

int main(int argc, char* argv[]){

    istringstream nEventsStream(argv[1]);
    int nEvents;
    if(!(nEventsStream>>nEvents)){
        std::cout<<"Invalid first argument!"<<endl;
        return 0;
    }

    istringstream mEnergyStream(argv[2]);
    int mEnergy;
    if(!(mEnergyStream>>mEnergy)){
        std::cout<<"Invalid second argument!"<<endl;
        return 0;
    }

    string folder = "";
    if(argc>3){
        istringstream mFolderStream(argv[3]);

        if(!(mFolderStream>>folder)){
            std::cout<<"Invalid third argument!"<<endl;
            return 0;
        }
    }

    std::cout<<"Number of events: "<<nEvents<<endl;
    std::cout<<"Energy: "<<mEnergy<<endl;
    string filename = to_string(nEvents)+"_"+to_string(mEnergy)+"GeV_"+string(argv[0]).substr(string(argv[0]).find_last_of("/\\")+1)+".root";
    // Create file on which histogram(s) can be saved.
    if(folder.length()!=0){
        filename = folder+"/"+filename;
    }
    TFile* outFile = new TFile(filename.c_str(), "RECREATE");
    std::cout<<"Output file: "<<filename<<endl;

    Pythia pythia;
    if(mEnergy==510){
        pythia.readString("Beams:eCM = 510"); //beam energy in GeV
    } else{
        std::cout<<"Input energy not equal to 510 GeV, if you wish to continue write \"yes\", otherwise write anything else"<<endl;
        string test;
        cin>>test;
        if(!strcmp(test.c_str(), ("yes"))){
            return 0;
        }
        pythia.readString("Beams:eCM = "+to_string(mEnergy));
    }

    //non-particle specific settings
    pythia.readString("SigmaTotal:zeroAXB = off");
    pythia.readString("SoftQCD:centralDiffractive = on");
    pythia.readString("Random:setSeed=on");
    pythia.readString("Random:seed=0");

    //For K0
    pythia.readString("310:onMode=0");
    pythia.readString("310:OnIfMatch=211 -211");
    pythia.readString("-310:onMode=0");
    pythia.readString("-310:OnIfMatch=-211 211");
    //For Lambda
    pythia.readString("3122:onMode=0");
    pythia.readString("3122:OnIfMatch=2212 -211");
    pythia.readString("-3122:onMode=0");
    pythia.readString("-3122:OnIfMatch=-2212 211");
    //For K*
    pythia.readString("313:onMode=0");
    pythia.readString("313:OnIfMatch=321 -211");
    pythia.readString("-313:onMode=0");
    pythia.readString("-313:OnIfMatch=-321 211");
    //For phi
    pythia.readString("333:onMode=0");
    pythia.readString("333:OnIfMatch=321 -321");
    pythia.readString("-333:onMode=0");
    pythia.readString("-333:OnIfMatch=-321 321");

    //setting counting to sparser
    pythia.readString("Next:numberCount=50000");
    //after the settings and all that
    pythia.init();

    //histograms
    std::vector<double> eta_efficiency_grid;
    for(double eta = -1.0; eta<1.11; eta += 0.05){
        eta_efficiency_grid.push_back(eta);
    }
    std::vector<double> pT_efficiency_grid;
    for(double pT = 0.0; pT<2.51; pT += 0.05){
        pT_efficiency_grid.push_back(pT);
    }
    TEfficiency efficiencyK0S("efficiencyK0S", "K^{0}_{S} efficiency;eta;p_{T}", eta_efficiency_grid.size()-1, eta_efficiency_grid.data(), pT_efficiency_grid.size()-1, pT_efficiency_grid.data());
    TEfficiency efficiencyLambda("efficiencyLambda", "#Lambda^{0} efficiency;eta;p_{T}", eta_efficiency_grid.size()-1, eta_efficiency_grid.data(), pT_efficiency_grid.size()-1, pT_efficiency_grid.data());
    TEfficiency efficiencyLambdabar("efficiencyLambdabar", "#bar{#Lambda}^{0} efficiency;eta;p_{T}", eta_efficiency_grid.size()-1, eta_efficiency_grid.data(), pT_efficiency_grid.size()-1, pT_efficiency_grid.data());
    TEfficiency efficiencyKstar("efficiencyKstar", "K^{*}(892) efficiency;eta;p_{T}", eta_efficiency_grid.size()-1, eta_efficiency_grid.data(), pT_efficiency_grid.size()-1, pT_efficiency_grid.data());
    TEfficiency efficiencyKstarbar("efficiencyKstarbar", "#bar{K}^{*}(892) efficiency;eta;p_{T}", eta_efficiency_grid.size()-1, eta_efficiency_grid.data(), pT_efficiency_grid.size()-1, pT_efficiency_grid.data());
    TEfficiency efficiencyphi("efficiencyphi", "#varphi(1020) efficiency;eta;p_{T}", eta_efficiency_grid.size()-1, eta_efficiency_grid.data(), pT_efficiency_grid.size()-1, pT_efficiency_grid.data());
    //Useful IDs
    const int K0sPDGid = 310;
    const int K0sbarPDGid = -310;
    const int LambdaPDGid = 3122;
    const int LambdabarPDGid = -3122;
    const int phiPDGid = 333;
    const int phibarPDGid = -333;
    const int KstarPDGid = 313;
    const int KstarbarPDGid = -313;
    const int piplusPDGid = 211;
    const int piminusPDGid = -211;
    const int KplusPDGid = 321;
    const int KminusPDGid = -321;
    const int pplusPDGid = 2212;
    const int pminusPDGid = -2212;

    //other things;
    int n_events = 0;
    std::vector<int> detected_particles_number;
    Rndm generator(0);
    bool AllDaughtersDetected;
    double momentumTPCRescale;

    // Begin event loop. Generate event; skip if generation aborted.
    for(int iEvent = 0; iEvent<nEvents; ++iEvent){
        if(!pythia.next()) continue;
        //checking if this event actually has particle I care about
        bool skipBecauseNoInterestingParticlesFound = true;
        for(int part_index = 0; part_index<pythia.event.size(); part_index++){
            switch(pythia.event[part_index].idAbs()){
            case K0sPDGid:
            case LambdaPDGid:
            case phiPDGid:
            case KstarPDGid:
                skipBecauseNoInterestingParticlesFound = false;
            default:
                break;
            }
            if(!skipBecauseNoInterestingParticlesFound){
                break;
            }
        }
        if(skipBecauseNoInterestingParticlesFound) continue;
        //cleaning up
        detected_particles_number.clear();

        //no proton cuts for now, so no proton safeguards needed

        //fishing for diffracted protons
        // int p1index = -1, p2index = -1;
        // for(int part_index = 0; part_index<pythia.event.size(); part_index++){
        //     if(pythia.event[part_index].id()!=pplusPDGid||pythia.event[pythia.event[part_index].mother1()].id()!=pplusPDGid){
        //         continue;
        //     }
        //     //choosing two with highest energies
        //     if(p1index<0){
        //         p1index = part_index;
        //     } else if(p2index<0){
        //         p2index = part_index;
        //         //switch so 1st proton has pz>0
        //         if(pythia.event[p1index].pz()<pythia.event[p2index].pz()){
        //             p2index = p1index;
        //             p1index = part_index;
        //         }
        //     } else if(pythia.event[p1index].pz()<pythia.event[part_index].pz()){
        //         p1index = part_index;
        //     } else if(pythia.event[p2index].pz()>pythia.event[part_index].pz()){
        //         p2index = part_index;
        //     }
        // }
        // //putting them into TLorentzVectors
        // TLorentzVector p1, p2;
        // p1.SetPxPyPzE(pythia.event[p1index].px(), pythia.event[p1index].py(), pythia.event[p1index].pz(), pythia.event[p1index].e());
        // p2.SetPxPyPzE(pythia.event[p2index].px(), pythia.event[p2index].py(), pythia.event[p2index].pz(), pythia.event[p2index].e());

        //CUTS
        // //proton fiducial
        // if(!(isInFiducial(p1)&&isInFiducial(p2))){
        //     continue;
        // }
        // //xi
        // if(!IsInXiFiducialNoBBCL(p1, p2, beamMomentum)){
        //     continue;
        // }
        // //non-elastic
        // if(IsInXiElasticSpot(p1, p2) or IsInMomElasticSpot(p1, p2)){
        //     continue;
        // }

        //making a list of detected particles
        //with simulated TPC acceptance
        //and simulated momentum reconstruction
        int TPCParticles = 0;
        for(int part_index = 0; part_index<pythia.event.size(); part_index++){
            if(!pythia.event[part_index].isFinal() or !pythia.event[part_index].isCharged()){
                continue;
            }
            momentumTPCRescale = 1+generator.gauss()*0.025;
            pythia.event[part_index].px(momentumTPCRescale*pythia.event[part_index].px());
            pythia.event[part_index].py(momentumTPCRescale*pythia.event[part_index].py());
            pythia.event[part_index].pz(momentumTPCRescale*pythia.event[part_index].pz());
            if(!isParticleInTPCAcceptance(pythia.event[part_index])){
                continue;
            }
            //reading acceptance from file ~/STAR-Analysis/share/etaPhiEfficiency_16_01_19_delta015_twoRuns.root
            if(isParticleDetected(&pythia.event[part_index])){
                TPCParticles++;
                detected_particles_number.push_back(part_index);
            }
        }

        // //production vertex closer to  primary vertex than 30mm
        // return abs(particle->eta())<0.7&&particle->vProd().pT()<30;

        //going through all the particles
        //just to  get to those we need for efficiency
        for(int part_index = 0; part_index<pythia.event.size(); part_index++){
            AllDaughtersDetected = true;
            //we are looking only for those particles
            switch(pythia.event[part_index].idAbs()){
            case K0sPDGid:
            case LambdaPDGid:
            case phiPDGid:
            case KstarPDGid:
                //check if all daugters of part_index particle are detected
                for(size_t daughter = 0; daughter<pythia.event[part_index].daughterList().size();daughter++){
                    if(std::find(detected_particles_number.begin(), detected_particles_number.end(), pythia.event[part_index].daughterList()[daughter])==detected_particles_number.end()){
                        AllDaughtersDetected = false;
                        break;
                    }
                }
                //the particle was created, but not detected because the whole event was thrown out
                if(TPCParticles<2){
                    AllDaughtersDetected = false;
                }
                //putting result (detection or not) into histogram
                switch(pythia.event[part_index].id()){
                case K0sPDGid:
                case K0sbarPDGid:
                    efficiencyK0S.Fill(AllDaughtersDetected, pythia.event[part_index].eta(), pythia.event[part_index].pT());
                    break;
                case LambdaPDGid:
                    efficiencyLambda.Fill(AllDaughtersDetected, pythia.event[part_index].eta(), pythia.event[part_index].pT());
                    break;
                case LambdabarPDGid:
                    efficiencyLambdabar.Fill(AllDaughtersDetected, pythia.event[part_index].eta(), pythia.event[part_index].pT());
                    break;
                case KstarPDGid:
                    efficiencyKstar.Fill(AllDaughtersDetected, pythia.event[part_index].eta(), pythia.event[part_index].pT());
                    break;
                case KstarbarPDGid:
                    efficiencyKstarbar.Fill(AllDaughtersDetected, pythia.event[part_index].eta(), pythia.event[part_index].pT());
                    break;
                case phiPDGid:
                case phibarPDGid:
                    efficiencyphi.Fill(AllDaughtersDetected, pythia.event[part_index].eta(), pythia.event[part_index].pT());
                    break;
                default:
                    break;
                }
            default:
                break;
            }

        }

        //total ratio
        n_events++;
    }

    //writing to file
    outFile->cd();
    outFile->Write();
    outFile->Close();

    // Statistics on event generation.
    pythia.stat();

    //other things
    std::cout<<"N_events: "<<n_events<<endl;

    return 0;
}

bool isInFiducial(TLorentzVector tested_particle){
    double px = tested_particle.X();
    double py = tested_particle.Y();
    bool f1 = (0.4<abs(py)&&abs(py)<0.8);
    bool f2 = (-0.27<px);
    bool f3 = (pow(px+0.6, 2)+pow(py, 2)<1.25);
    return f1&&f2&&f3;
}

bool IsInXiElasticSpot(TLorentzVector east, TLorentzVector west){
    //before Afterburner
    // double x_0 = 4.30588e-03;
    // double sigma_x = 2.02340e-03;
    // double y_0 = 1.72097e-03;
    // double sigma_y = 2.26638e-03;
    //after Afterburner
    double x_0 = -4.48170e-04;
    double sigma_x = 1.79095e-03;
    double y_0 = -8.04898e-04;
    double sigma_y = 2.12035e-03;
    return pow((xi(east)-x_0)/sigma_x, 2)+pow((xi(west)-y_0)/sigma_y, 2)<3*3;
}

bool IsInMomElasticSpot(TLorentzVector east, TLorentzVector west){
    //before Afterburner
    // double x_0 = -3.82151e-02;
    // double sigma_x = 3.67545e-02;
    // double y_0 = 1.98348e-03;
    // double sigma_y = 3.40440e-02;
    //after Afterburner
    double x_0 = 5.06472e-03;
    double sigma_x = 3.42004e-02;
    double y_0 = 5.98219e-04;
    double sigma_y = 3.15726e-02;
    double x = east.X()+west.X();
    double y = east.Y()+west.Y();
    return pow((x-x_0)/sigma_x, 2)+pow((y-y_0)/sigma_y, 2)<3*3;
}

bool IsInXiFiducialNoBBCL(TLorentzVector east, TLorentzVector west, double p = 254.867){
    if(xi(east, p)<0.005 or xi(east, p)>0.2){
        return false;
    }
    if(xi(west, p)<0.005 or xi(west, p)>0.2){
        return false;
    }
    return true;
}

bool IsInXiFiducialBBCL(TLorentzVector east, TLorentzVector west, double p = 254.867){
    if(xi(east, p)<0.005 or xi(east, p)>0.08){
        return false;
    }
    if(xi(west, p)<0.005 or xi(west, p)>0.08){
        return false;
    }
    return true;
}