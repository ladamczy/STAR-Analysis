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

int main(int argc, char *argv[]){

    istringstream nEventsStream(argv[1]);
    int nEvents;
    if(!(nEventsStream>>nEvents)){
        cout<<"Invalid first argument!"<<endl;
        return 0;
    }

    istringstream mEnergyStream(argv[2]);
    int mEnergy;
    if(!(mEnergyStream>>mEnergy)){
        cout<<"Invalid second argument!"<<endl;
        return 0;
    }

    string folder = "";
    if(argc>3){
        istringstream mFolderStream(argv[3]);

        if(!(mFolderStream>>folder)){
            cout<<"Invalid third argument!"<<endl;
            return 0;
        }
    }

    cout<<"Number of events: "<<nEvents<<endl;
    cout<<"Energy: "<<mEnergy<<endl;
    string filename = to_string(nEvents)+"_"+to_string(mEnergy)+"GeV_"+string(argv[0]).substr(string(argv[0]).find_last_of("/\\")+1)+".root";
    // Create file on which histogram(s) can be saved.
    if(folder.length()!=0){
        filename = folder+"/"+filename;
    }
    TFile *outFile = new TFile(filename.c_str(), "RECREATE");
    cout<<"Output file: "<<filename<<endl;

    Pythia pythia;
    if(mEnergy==510){
        pythia.readString("Beams:eCM = 510"); //beam energy in GeV
    } else{
        cout<<"Input energy not equal to 510 GeV, if you wish to continue write \"yes\", otherwise write anything else"<<endl;
        string test;
        cin>>test;
        if(!strcmp(test.c_str(), ("yes"))){
            return 0;
        }
        pythia.readString("Beams:eCM = "+to_string(mEnergy));
    }

    pythia.readString("SigmaTotal:zeroAXB = off");
    pythia.readString("SoftQCD:centralDiffractive = on");

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
    pythia.init();

    //histogram variables
    int stateBins = 40;
    double stateMassMax = 40;

    //histograms
    TH1D K0multiplicity("K0multiplicity", "K^{0}_{S} multiplicity", 10, 0, 10);
    TH1D K0AdditionalTracks("K0AdditionalTracks", "Additional tracks except K^{0}_{S} daughters", 20, 0, 20);
    TH2D K0multiplicityXmass("K0multiplicityXmass", "K^{0}_{S} multiplicity with respect to mass of a central state;m_{X} [GeV]", stateBins, 0, stateMassMax, 10, 0, 10);
    TEfficiency K0ratioXmass("K0ratioXmass", "Ratio of detected K^{0}_{S} to all measured events with respect to mass of a central state;m_{X} [GeV]", stateBins, 0, stateMassMax);
    TH2D K0AdditionalTracksXmass("K0AdditionalTracksXmass", "Additional tracks except K^{0}_{S} daughters;m_{X} [GeV]", stateBins, 0, stateMassMax, 20, 0, 20);
    TH1D K0AdditionalTracksWithAcceptance("K0AdditionalTracksWithAcceptance", "Additional tracks except K^{0}_{S} daughters, with TPC acceptance", 20, 0, 20);
    TH2D K0AdditionalTracksXmassWithAcceptance("K0AdditionalTracksXmassWithAcceptance", "Additional tracks except K^{0}_{S} daughters, with TPC acceptance;m_{X} [GeV]", stateBins, 0, stateMassMax, 20, 0, 20);


    TH1D Lambdamultiplicity("Lambdamultiplicity", "#Lambda^{0} multiplicity", 10, 0, 10);
    TH1D LambdaAdditionalTracks("LambdaAdditionalTracks", "Additional tracks except #Lambda^{0} daughters", 20, 0, 20);
    TH2D LambdamultiplicityXmass("LambdamultiplicityXmass", "#Lambda^{0} multiplicity with respect to mass of a central state;m_{X} [GeV]", stateBins, 0, stateMassMax, 10, 0, 10);
    TEfficiency LambdaratioXmass("LambdaratioXmass", "Ratio of detected #Lambda^{0} to all measured events with respect to mass of a central state;m_{X} [GeV]", stateBins, 0, stateMassMax);
    TH2D LambdaAdditionalTracksXmass("LambdaAdditionalTracksXmass", "Additional tracks except #Lambda^{0} daughters;m_{X} [GeV]", stateBins, 0, stateMassMax, 20, 0, 20);
    TH1D LambdaAdditionalTracksWithAcceptance("LambdaAdditionalTracksWithAcceptance", "Additional tracks except #Lambda^{0} daughters, with TPC acceptance", 20, 0, 20);
    TH2D LambdaAdditionalTracksXmassWithAcceptance("LambdaAdditionalTracksXmassWithAcceptance", "Additional tracks except #Lambda^{0} daughters, with TPC acceptance;m_{X} [GeV]", stateBins, 0, stateMassMax, 20, 0, 20);

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
    int n_K0 = 0;
    int n_Lambda = 0;
    Rndm generator(0);

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
        //no BBCL
        if(hasAtLeastNParticlesInBBCLarge(pythia.event, 1)){
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

        //number of K0 & Lambda
        int nTracksInTPCAcceptance = 0;
        int nK0 = 0;
        int nLambda = 0;
        std::vector<int> K0indices;
        std::vector<int> Lambdaindices;
        int tempDaughter1, tempDaughter2;
        std::vector<Pythia8::Particle> acceptedParticlesPlusDaughters;
        //gathering info in main loop
        for(int part_index = 0; part_index<pythia.event.size(); part_index++){
            if(isParticleInTPCAcceptance(pythia.event[part_index])){ nTracksInTPCAcceptance++; }
            if(abs(pythia.event[part_index].id())==K0sPDGid){ K0indices.push_back(part_index); }
            if(abs(pythia.event[part_index].id())==LambdaPDGid){ Lambdaindices.push_back(part_index); }
        }
        for(long unsigned int K0candidate = 0; K0candidate<K0indices.size(); K0candidate++){
            tempDaughter1 = pythia.event[K0indices[K0candidate]].daughter1();
            tempDaughter2 = pythia.event[K0indices[K0candidate]].daughter2();
            if(isParticleInTPCAcceptance(pythia.event[tempDaughter1])&&isParticleInTPCAcceptance(pythia.event[tempDaughter2])&&(tempDaughter1+1==tempDaughter2)){
                nK0++;
            }
        }
        for(long unsigned int Lambdacandidate = 0; Lambdacandidate<Lambdaindices.size(); Lambdacandidate++){
            tempDaughter1 = pythia.event[Lambdaindices[Lambdacandidate]].daughter1();
            tempDaughter2 = pythia.event[Lambdaindices[Lambdacandidate]].daughter2();
            if(isParticleInTPCAcceptance(pythia.event[tempDaughter1])&&isParticleInTPCAcceptance(pythia.event[tempDaughter2])&&(tempDaughter1+1==tempDaughter2)){
                nLambda++;
            }
        }
        TLorentzVector temp4Vector = { 0,0,0,0 };
        temp4Vector.SetE(mEnergy);
        double XStateMass = (temp4Vector-p1-p2).M();

        //total ratio
        n_events++;
        n_K0 += nK0;
        n_Lambda += nLambda;

        //for plots requiring simulated TPC acceptance
        for(int part_index = 0; part_index<pythia.event.size(); part_index++){
            if(!isParticleInTPCAcceptance(pythia.event[part_index])){
                continue;
            }
            if(abs(pythia.event[pythia.event[part_index].mother1()].id())==K0sPDGid&&nK0>0){
                acceptedParticlesPlusDaughters.push_back(pythia.event[part_index]);
            } else if(abs(pythia.event[pythia.event[part_index].mother1()].id())==LambdaPDGid&&nLambda>0){
                acceptedParticlesPlusDaughters.push_back(pythia.event[part_index]);
            } else if(generator.flat()<0.6){
                acceptedParticlesPlusDaughters.push_back(pythia.event[part_index]);
            }
        }

        //histogram filling
        K0multiplicity.Fill(nK0);
        K0multiplicityXmass.Fill(XStateMass, nK0);
        K0ratioXmass.Fill(nK0>0, XStateMass);
        if(nK0>0){
            K0AdditionalTracks.Fill(nTracksInTPCAcceptance-nK0*2);
            K0AdditionalTracksXmass.Fill(XStateMass, nTracksInTPCAcceptance-nK0*2);
            K0AdditionalTracksWithAcceptance.Fill(acceptedParticlesPlusDaughters.size()-nK0*2);
            K0AdditionalTracksXmassWithAcceptance.Fill(XStateMass, acceptedParticlesPlusDaughters.size()-nK0*2);
        }

        Lambdamultiplicity.Fill(nLambda);
        LambdamultiplicityXmass.Fill(XStateMass, nLambda);
        LambdaratioXmass.Fill(nLambda>0, XStateMass);
        if(nLambda>0){
            LambdaAdditionalTracks.Fill(nTracksInTPCAcceptance-nLambda*2);
            LambdaAdditionalTracksXmass.Fill(XStateMass, nTracksInTPCAcceptance-nLambda*2);
            LambdaAdditionalTracksWithAcceptance.Fill(acceptedParticlesPlusDaughters.size()-nLambda*2);
            LambdaAdditionalTracksXmassWithAcceptance.Fill(XStateMass, acceptedParticlesPlusDaughters.size()-nLambda*2);
        }
    }

    outFile->cd();
    outFile->Write();
    outFile->Close();

    // Statistics on event generation.
    pythia.stat();

    //other things
    cout<<"N_events: "<<n_events<<endl;
    cout<<"N_K0: "<<n_K0<<endl;
    cout<<"N_Lambda: "<<n_Lambda<<endl;
    std::stringstream ratio_output;
    ratio_output<<std::fixed<<std::setprecision(4)<<double(n_K0)/n_events;
    cout<<"N_K0/N_events: "<<ratio_output.str()<<endl;
    ratio_output.str(""); //clearing
    ratio_output<<std::fixed<<std::setprecision(4)<<double(n_Lambda)/n_events;
    cout<<"N_Lambda0/N_events: "<<ratio_output.str()<<endl;

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