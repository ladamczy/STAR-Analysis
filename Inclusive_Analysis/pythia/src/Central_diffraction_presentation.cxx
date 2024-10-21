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

bool isInFiducial(TLorentzVector);
bool IsInXiElasticSpot(TLorentzVector, TLorentzVector);
bool IsInMomElasticSpot(TLorentzVector, TLorentzVector);
bool IsInXiFiducialNoBBCL(TLorentzVector, TLorentzVector, double);
bool IsInXiFiducialBBCL(TLorentzVector, TLorentzVector, double);

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
    string filename;
    filename = "Pythia_presentation_simulation.root";
    // Create file on which histogram(s) can be saved.
    if(folder.length()!=0){
        filename = folder+"/"+filename;
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

    //non-particle specific settings
    pythia.readString("SigmaTotal:zeroAXB = off");
    pythia.readString("SoftQCD:centralDiffractive = on");
    pythia.readString("Random:setSeed=on");
    pythia.readString("Random:seed=0");

    //For K0
    // pythia.readString("310:onMode=0");
    // pythia.readString("310:OnIfMatch=211 -211");
    // pythia.readString("-310:onMode=0");
    // pythia.readString("-310:OnIfMatch=-211 211");
    //For Lambda
    // pythia.readString("3122:onMode=0");
    // pythia.readString("3122:OnIfMatch=2212 -211");
    // pythia.readString("-3122:onMode=0");
    // pythia.readString("-3122:OnIfMatch=-2212 211");

    //after the settings and all that
    pythia.init();

    //histogram variables
    // int n_ptBins = 9;
    // double ptBins[] = { 0,0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8 };
    // int n_etaBins = 10;
    // double etaBins[] = { -1., -0.8, -0.6, -0.4, -0.2, 0., 0.2, 0.4, 0.6, 0.8, 1.0 };
    // int stateBins = 40;
    // double stateMassMax = 40;

    //histograms
    TH1D PythiaFigure3_6ab = TH1D("PythiaFigure3_6ab", ";#eta;number of tracks", 100, -2, 2);
    TH1D PythiaFigure3_6c = TH1D("PythiaFigure3_6c", ";p_{T} [GeV/c];number of tracks", 50, 0, 5);
    TH1D PythiaFigure3_6d = TH1D("PythiaFigure3_6d", ";#Phi [rad];number of tracks", 100, -TMath::Pi(), TMath::Pi());
    TEfficiency PythiaFigure3_7 = TEfficiency("PythiaFigure3_7", ";log_{10}(#xi_{1}#xi_{2});#varepsilon_{n_{sel}>=2}", 100, -3, 0);
    TEfficiency PythiaFigure3_7noAcceptance = TEfficiency("PythiaFigure3_7noAcceptance", ";log_{10}(#xi_{1}#xi_{2});#varepsilon_{n_{sel}>=2}", 100, -3, 0);
    TH1D PythiaFigure3_18g = TH1D("PythiaFigure3_18g", ";#xi_{1}#xi_{2};#frac{1}{N} #frac{dN}{d(#xi_{1}#xi_{2})}", 1000, 0, 0.04);
    TH1D PythiaFigure3_18gcloser = TH1D("PythiaFigure3_18gcloser", ";#xi_{1}#xi_{2};#frac{1}{N_{ev}} #frac{dN_{ev}}{d(#xi_{1}#xi_{2})}", 200, 0, 2e-4);
    TH1D PythiaFigure3_19a = TH1D("PythiaFigure3_19a", ";n_{sel};#frac{1}{N} #frac{dN}{dn_{sel}}", 10, 2, 12);
    TH1D PythiaFigure3_20 = TH1D("PythiaFigure3_20", ";p_{T} [GeV/c];#frac{1}{N} #frac{dN}{dp_{T}}", 50, 0, 5);
    TH1D PythiaFigure3_21 = TH1D("PythiaFigure3_21", ";#eta;#frac{1}{N} #frac{dN}{d#eta}", 100, -2, 2);

    //control histogram
    TH1D PythiaControl = TH1D("PythiaControl", ";cathegories;events", 1, 0, 1);


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
    Rndm generator(0);
    std::vector<int> detected_particles_number;
    Particle tempParticle;

    // Begin event loop. Generate event; skip if generation aborted.
    for(int iEvent = 0; iEvent<nEvents; ++iEvent){
        if(!pythia.next()) continue;
        //cleaning up
        detected_particles_number.clear();
        //fishing for diffracted protons
        int p1index = -1, p2index = -1;
        for(int part_index = 0; part_index<pythia.event.size(); part_index++){
            if(pythia.event[part_index].id()!=pplusPDGid||pythia.event[pythia.event[part_index].mother1()].id()!=pplusPDGid){
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
        if(!IsInXiFiducialNoBBCL(p1, p2, beamMomentum)){
            continue;
        }
        //non-elastic
        if(IsInXiElasticSpot(p1, p2) or IsInMomElasticSpot(p1, p2)){
            continue;
        }
        //making a list of detected particles
        //with TPC acceptance of 60%
        int TPCParticlesForNsel = 0;
        //and without
        int TPCParticlesForNch = 0;
        //loop itself
        for(size_t i = 0; i<pythia.event.size(); i++){
            if(!isParticleInTPCAcceptance(pythia.event[i])){
                continue;
            }
            //simulated acceptance
            if(generator.flat()<0.6){
                detected_particles_number.push_back(i);
                if(abs(pythia.event[i].eta())<0.7){
                    TPCParticlesForNsel++;
                }
            }
            //without simulated acceptance
            if(abs(pythia.event[i].eta())<0.7){
                TPCParticlesForNch++;
            }
        }

        //histograms just after RP cuts
        for(size_t i = 0; i<detected_particles_number.size(); i++){
            tempParticle = pythia.event[detected_particles_number[i]];
            PythiaFigure3_6ab.Fill(tempParticle.eta());
            PythiaFigure3_6c.Fill(tempParticle.pT());
            PythiaFigure3_6d.Fill(tempParticle.phi());
        }
        //histogram 3.7; a special case
        PythiaFigure3_7.Fill(TPCParticlesForNsel>=2, log10(xi(p1)*xi(p2)));
        PythiaFigure3_7noAcceptance.Fill(TPCParticlesForNsel>=2, log10(xi(p1)*xi(p2)));

        //n_sel filter
        PythiaControl.Fill("before n_sel", 1);
        if(TPCParticlesForNsel<2){
            continue;
        }
        PythiaControl.Fill("after n_sel", 1);

        //histograms after the n_sel>=2 filter
        PythiaFigure3_18g.Fill(xi(p1)*xi(p2));
        PythiaFigure3_18gcloser.Fill(xi(p1)*xi(p2));
        PythiaFigure3_19a.Fill(TPCParticlesForNsel);
        for(size_t i = 0; i<detected_particles_number.size(); i++){
            tempParticle = pythia.event[detected_particles_number[i]];
            PythiaFigure3_20.Fill(tempParticle.pT());
            PythiaFigure3_21.Fill(tempParticle.eta());
        }

        //total ratio
        n_events++;
    }

    //closing touches
    PythiaControl.LabelsDeflate();

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