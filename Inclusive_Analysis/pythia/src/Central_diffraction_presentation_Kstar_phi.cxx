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

// My classes and stuff
#include "UsefulThingsPythia.h"
// #include "UsefulThings.h"

using namespace Pythia8;
using namespace std;

bool isInFiducial(TLorentzVector);
bool IsInXiElasticSpot(TLorentzVector, TLorentzVector);
bool IsInMomElasticSpot(TLorentzVector, TLorentzVector);
bool IsInXiFiducialNoBBCL(TLorentzVector east, TLorentzVector west, double p = 254.867);

int main(int argc, char *argv[]){

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
    string filename = "Pythia_phi_simulation.root";
    // Create file on which histogram(s) can be saved.
    if(folder.length()!=0){
        filename = folder+filename;
    }
    TFile *outFile = new TFile(filename.c_str(), "RECREATE");
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

    pythia.readString("SigmaTotal:zeroAXB = off");
    pythia.readString("SoftQCD:centralDiffractive = on");

    //For phi
    pythia.readString("333:onMode=0");
    pythia.readString("333:OnIfMatch=321 -321");
    pythia.readString("-333:onMode=0");
    pythia.readString("-333:OnIfMatch=-321 321");
    pythia.init();

    //histogram variables
    int n_ptBins = 9;
    double ptBins[] = { 0,0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8 };
    int n_etaBins = 10;
    double etaBins[] = { -1., -0.8, -0.6, -0.4, -0.2, 0., 0.2, 0.4, 0.6, 0.8, 1.0 };

    //histograms
    TH1D PhiptHist("PhiptHist", "PhiptHist", n_ptBins, ptBins);
    TH1D PhietaHist("PhietaHist", "PhietaHist", n_etaBins, etaBins);


    //Useful IDs
    const int K0sPDGid = 310;
    const int K0sbarPDGid = -310;
    const int LambdaPDGid = 3122;
    const int LambdabarPDGid = -3122;
    const int piplusPDGid = 211;
    const int piminusPDGid = -211;
    const int pplusPDGid = 2212;
    const int pminusPDGid = -2212;

    //other things;
    int n_events = 0;
    int n_Phi = 0;

    // Begin event loop. Generate event; skip if generation aborted.
    for(int iEvent = 0; iEvent<nEvents; ++iEvent){
        if(!pythia.next()) continue;
        //fishing for diffracted protons
        int p1index = -1, p2index = -1;
        for(int part_index = 0; part_index<pythia.event.size(); part_index++){
            if(pythia.event[part_index].id()!=2212||pythia.event[pythia.event[part_index].mother1()].id()!=2212){
                continue;
            }
            //choosing two with highest energies
            if(p1index<0){
                p1index = part_index;
            } else if(p2index<0){
                p2index = part_index;
                //switch so 1st proton has pz>0
                if(pythia.event[p1index].pz()<pythia.event[p2index].pz()){
                    p2index = p1index;
                    p1index = part_index;
                }
            } else if(pythia.event[p1index].pz()<pythia.event[part_index].pz()){
                p1index = part_index;
            } else if(pythia.event[p2index].pz()>pythia.event[part_index].pz()){
                p2index = part_index;
            }
        }
        //putting them into TLorentzVectors
        TLorentzVector p1, p2;
        p1.SetPxPyPzE(pythia.event[p1index].px(), pythia.event[p1index].py(), pythia.event[p1index].pz(), pythia.event[p1index].e());
        p2.SetPxPyPzE(pythia.event[p2index].px(), pythia.event[p2index].py(), pythia.event[p2index].pz(), pythia.event[p2index].e());

        //CUTS
        //proton fiducial
        if(!(isInFiducial(p1)&&isInFiducial(p2))){
            continue;
        }
        //xi
        if(!IsInXiFiducialNoBBCL(p1, p2)){
            continue;
        }
        //non-elastic
        if(IsInXiElasticSpot(p1, p2) or IsInMomElasticSpot(p1, p2)){
            continue;
        }
        //at least 2 with TOF and opposite signs
        if(!hasAtLeastNParticlesInTPC(pythia.event, 2)){
            continue;
        }
        int charge = 0;
        int TPCParticles = 0;
        for(size_t i = 0; i<pythia.event.size(); i++){
            if(isParticleInTPCAcceptance(pythia.event[i])){
                TPCParticles++;
                //for some reason carge is in normal units now???
                charge += pythia.event[i].charge();
            }
        }
        if(TPCParticles==abs(charge)){
            continue;
        }

        //number of Phi
        int nPhi = 0;
        std::vector<int> Phiindices;
        int tempDaughter1, tempDaughter2;
        std::vector<Pythia8::Particle> acceptedParticlesPlusDaughters;
        //gathering info in main loop
        for(int part_index = 0; part_index<pythia.event.size(); part_index++){
            if(abs(pythia.event[part_index].id())==333){ Phiindices.push_back(part_index); }
        }
        //Phis
        for(long unsigned int Phicandidate = 0; Phicandidate<Phiindices.size(); Phicandidate++){
            tempDaughter1 = pythia.event[Phiindices[Phicandidate]].daughter1();
            tempDaughter2 = pythia.event[Phiindices[Phicandidate]].daughter2();
            if(isParticleInTPCAcceptance(pythia.event[tempDaughter1])&&isParticleInTPCAcceptance(pythia.event[tempDaughter2])&&(tempDaughter1+1==tempDaughter2)){
                //additional cause i got that |eta|<0.7 ugh
                if(abs(pythia.event[tempDaughter1].eta())<0.7&&abs(pythia.event[tempDaughter2].eta())<0.7){
                    nPhi++;
                    PhiptHist.Fill(pythia.event[Phiindices[Phicandidate]].pT());
                    PhietaHist.Fill(pythia.event[Phiindices[Phicandidate]].eta());
                }
            }
        }

        //total ratio
        n_events++;
        n_Phi += nPhi;
    }

    outFile->cd();
    outFile->Write();
    outFile->Close();

    // Statistics on event generation.
    pythia.stat();

    //other things
    std::cout<<"N_events: "<<n_events<<endl;
    std::cout<<"N_phi: "<<n_Phi<<endl;
    std::stringstream ratio_output;
    ratio_output<<std::fixed<<std::setprecision(4)<<double(n_Phi)/n_events;
    std::cout<<"N_phi/N_events: "<<ratio_output.str()<<endl;
    ratio_output.str(""); //clearing

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

bool IsInXiFiducialNoBBCL(TLorentzVector east, TLorentzVector west, double p){
    if(xi(east, p)<0.005 or xi(east, p)>0.2){
        return false;
    }
    if(xi(west, p)<0.005 or xi(west, p)>0.2){
        return false;
    }
    return true;
}