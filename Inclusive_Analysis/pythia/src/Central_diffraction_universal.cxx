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
bool IsInMeasurementSpecificCuts(Particle*);

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
    string filename = to_string(nEvents)+"_"+to_string(mEnergy)+"GeV_"+string(argv[0]).substr(string(argv[0]).find_last_of("/\\")+1)+".root";
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

    //setting counting to sparser
    pythia.readString("Next:numberCount=50000");
    //after the settings and all that
    pythia.init();

    //histograms
    //general
    TH1D ParticlesDetected("ParticlesDetected", "ParticlesDetected", 20000, -10000, 10000);
    TH1D ParticlesReconstructedTwoDaughters("ParticlesReconstructedTwoDaughters", "ParticlesReconstructedTwoDaughters", 20000, -10000, 10000);
    TH1D ParticlesReconstructedMultipleDaughters("ParticlesReconstructedMultipleDaughters", "ParticlesReconstructedMultipleDaughters", 20000, -10000, 10000);
    TH1D MpipiSum("MpipiSum", ";m_{#pi#pi}", 500, 0, 5);
    TH1D MKKSum("MKKSum", ";m_{KK}", 500, 0, 5);
    TH1D MKpiSum("MKpiSum", ";m_{K#pi}", 500, 0, 5);
    TH1D MPIDSum("MPIDSum", ";m_{pair}", 500, 0, 5);
    TH1D MKKBackgroundPositiveSum("MKKBackgroundPositiveSum", ";m_{KK}", 500, 0, 5);
    TH1D MKKBackgroundNegativeSum("MKKBackgroundNegativeSum", ";m_{KK}", 500, 0, 5);
    TH1D t0Resonant("t0Resonant", ";t_{0} [ns]", 100, 0, 50);
    TH1D t0Nonresonant("t0Nonresonant", ";t_{0} [ns]", 100, 0, 50);
    //Detector simulation
    //everything as pion pairs
    TH1D MpipiNonresonant("MpipiNonresonant", ";m_{#pi#pi}", 500, 0, 5);
    TH1D MpipiResonant("MpipiResonant", ";m_{#pi#pi}", 500, 0, 5);
    TH2D MpipiParticles("MpipiParticles", "", 1, 0, 1, 500, 0, 5);
    //everything as kaon pairs
    TH1D MKKNonresonant("MKKNonresonant", ";m_{KK}", 500, 0, 5);
    TH1D MKKResonant("MKKResonant", ";m_{KK}", 500, 0, 5);
    TH2D MKKParticles("MKKParticles", "", 1, 0, 1, 500, 0, 5);
    //everything as kaon-pion pairs
    TH1D MKpiNonresonant("MKpiNonresonant", ";m_{K#pi}", 500, 0, 5);
    TH1D MKpiResonant("MKpiResonant", ";m_{K#pi}", 500, 0, 5);
    TH2D MKpiParticles("MKpiParticles", "", 1, 0, 1, 500, 0, 5);
    //the rest
    TH1D MPIDNonresonant("MPIDNonresonant", ";m_{pair}", 500, 0, 5);
    TH1D MPIDResonant("MPIDResonant", ";m_{pair}", 500, 0, 5);
    TH1D MKKBackgroundPositiveNonresonant("MKKBackgroundPositiveNonresonant", ";m_{KK}", 500, 0, 5);
    TH1D MKKBackgroundPositiveResonant("MKKBackgroundPositiveResonant", ";m_{KK}", 500, 0, 5);
    TH1D MKKBackgroundNegativeNonresonant("MKKBackgroundNegativeNonresonant", ";m_{KK}", 500, 0, 5);
    TH1D MKKBackgroundNegativeResonant("MKKBackgroundNegativeResonant", ";m_{KK}", 500, 0, 5);
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
    std::vector<int> detected_particles_number;
    std::vector<int> detected_particles_positive_number;
    std::vector<int> detected_particles_negative_number;
    std::vector<int> mother_particles_number;
    Rndm generator(0);
    TLorentzVector temp4VectorPositive;
    TLorentzVector temp4VectorNegative;
    bool AllDaughtersDetected;
    int DaughtersDetected;
    double momentumTPCRescale;

    // Begin event loop. Generate event; skip if generation aborted.
    for(int iEvent = 0; iEvent<nEvents; ++iEvent){
        if(!pythia.next()) continue;
        //cleaning up
        detected_particles_number.clear();
        detected_particles_positive_number.clear();
        detected_particles_negative_number.clear();
        mother_particles_number.clear();
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
            if(!IsInMeasurementSpecificCuts(&pythia.event[part_index])){
                continue;
            }
            if(isParticleDetected(&pythia.event[part_index])){
                TPCParticles++;
                detected_particles_number.push_back(part_index);
            }
        }
        if(TPCParticles<2){
            continue;
        }

        //gathering info in detection loop and assigning mothers
        for(size_t i = 0; i<detected_particles_number.size(); i++){
            ParticlesDetected.Fill(pythia.event[detected_particles_number[i]].name().c_str(), 1.);
            //getting mother: singular hadron
            //if hadron is not the mother, then it wont be saved
            //it works because all detected particles are sorted as positive or negative
            //and all particles without hadron, proper mothers, go to "Nonresonant" histogram
            if(pythia.event[detected_particles_number[i]].mother1()>0&&pythia.event[detected_particles_number[i]].mother2()==0){
                t0Resonant.Fill(pythia.event[detected_particles_number[i]].tProd()/299.792458); //300 mm/c = 1ns <-> 1 mm/c = 1/300 ns
                if(pythia.event[pythia.event[detected_particles_number[i]].mother1()].isHadron()){
                    mother_particles_number.push_back(pythia.event[detected_particles_number[i]].mother1());
                } else{
                    printf("Particle %d, a mother, is not a hadron\n", pythia.event[detected_particles_number[i]].mother1());
                    pythia.event.list(true, true);
                }
            } else{
                t0Nonresonant.Fill(pythia.event[detected_particles_number[i]].tProd()/299.792458); //300 mm/c = 1ns <-> 1 mm/c = 1/300 ns
            }
            //sorting into positive and negative
            if(pythia.event[detected_particles_number[i]].charge()>0){
                detected_particles_positive_number.push_back(detected_particles_number[i]);
            } else{
                detected_particles_negative_number.push_back(detected_particles_number[i]);
            }
        }

        //removing potential mother particle duplicates
        sort(mother_particles_number.begin(), mother_particles_number.end());
        mother_particles_number.erase(unique(mother_particles_number.begin(), mother_particles_number.end()), mother_particles_number.end());
        //mother loop
        for(size_t i = 0; i<mother_particles_number.size(); i++){
            DaughtersDetected = 0;
            AllDaughtersDetected = true;
            //check if all daugters of part_index particle are detected
            for(size_t daughter = 0; daughter<pythia.event[mother_particles_number[i]].daughterList().size();daughter++){
                if(std::find(detected_particles_number.begin(), detected_particles_number.end(), pythia.event[mother_particles_number[i]].daughterList()[daughter])!=detected_particles_number.end()){
                    DaughtersDetected++;
                } else{
                    AllDaughtersDetected = false;
                }
            }
            //doing stuff with their mothers
            if(AllDaughtersDetected&&DaughtersDetected==2){
                ParticlesReconstructedTwoDaughters.Fill(pythia.event[mother_particles_number[i]].name().c_str(), 1.);
            }
            if(AllDaughtersDetected&&DaughtersDetected>2){
                ParticlesReconstructedMultipleDaughters.Fill(pythia.event[mother_particles_number[i]].name().c_str(), 1.);
            }
        }

        // part for making a histogram of what is detected (opposite sign pairs)
        for(size_t i = 0; i<detected_particles_positive_number.size(); i++){
            for(size_t j = 0; j<detected_particles_negative_number.size(); j++){
                // assuming everything is a pion
                //pair 4momenta setting
                temp4VectorPositive.SetPtEtaPhiM(pythia.event[detected_particles_positive_number[i]].pT(),
                    pythia.event[detected_particles_positive_number[i]].eta(),
                    pythia.event[detected_particles_positive_number[i]].phi(),
                    particleMass[Pion]);
                temp4VectorNegative.SetPtEtaPhiM(pythia.event[detected_particles_negative_number[j]].pT(),
                    pythia.event[detected_particles_negative_number[j]].eta(),
                    pythia.event[detected_particles_negative_number[j]].phi(),
                    particleMass[Pion]);
                //pair 4momenta using and filling
                double mass = (temp4VectorPositive+temp4VectorNegative).M();
                int motherPositive = pythia.event[detected_particles_positive_number[i]].mother1();
                int motherNegative = pythia.event[detected_particles_negative_number[j]].mother1();
                bool t1 = std::find(mother_particles_number.begin(), mother_particles_number.end(), motherPositive)!=mother_particles_number.end();
                bool t2 = std::find(mother_particles_number.begin(), mother_particles_number.end(), motherNegative)!=mother_particles_number.end();
                bool t3 = motherNegative==motherPositive;
                bool t4 = pythia.event[detected_particles_positive_number[i]].mother2()==0&&pythia.event[detected_particles_negative_number[j]].mother2()==0;
                if(t1&&t2&&t3&&t4){
                    MpipiResonant.Fill(mass);
                    MpipiParticles.Fill(pythia.particleData.name(pythia.event[motherPositive].idAbs()).c_str(), mass, 1.);
                } else{
                    MpipiNonresonant.Fill(mass);
                }
                MpipiSum.Fill(mass);

                // assuming everything is a kaon
                //pair 4momenta setting
                temp4VectorPositive.SetPtEtaPhiM(pythia.event[detected_particles_positive_number[i]].pT(),
                    pythia.event[detected_particles_positive_number[i]].eta(),
                    pythia.event[detected_particles_positive_number[i]].phi(),
                    particleMass[Kaon]);
                temp4VectorNegative.SetPtEtaPhiM(pythia.event[detected_particles_negative_number[j]].pT(),
                    pythia.event[detected_particles_negative_number[j]].eta(),
                    pythia.event[detected_particles_negative_number[j]].phi(),
                    particleMass[Kaon]);
                //pair 4momenta using and filling
                mass = (temp4VectorPositive+temp4VectorNegative).M();
                motherPositive = pythia.event[detected_particles_positive_number[i]].mother1();
                motherNegative = pythia.event[detected_particles_negative_number[j]].mother1();
                t1 = std::find(mother_particles_number.begin(), mother_particles_number.end(), motherPositive)!=mother_particles_number.end();
                t2 = std::find(mother_particles_number.begin(), mother_particles_number.end(), motherNegative)!=mother_particles_number.end();
                t3 = motherNegative==motherPositive;
                t4 = pythia.event[detected_particles_positive_number[i]].mother2()==0&&pythia.event[detected_particles_negative_number[j]].mother2()==0;
                if(t1&&t2&&t3&&t4){
                    MKKResonant.Fill(mass);
                    MKKParticles.Fill(pythia.particleData.name(pythia.event[motherPositive].idAbs()).c_str(), mass, 1.);
                } else{
                    MKKNonresonant.Fill(mass);
                }
                MKKSum.Fill(mass);

                // assuming kaon-pion pairs, both directions
                //pair 4momenta setting
                temp4VectorPositive.SetPtEtaPhiM(pythia.event[detected_particles_positive_number[i]].pT(),
                    pythia.event[detected_particles_positive_number[i]].eta(),
                    pythia.event[detected_particles_positive_number[i]].phi(),
                    particleMass[Kaon]);
                temp4VectorNegative.SetPtEtaPhiM(pythia.event[detected_particles_negative_number[j]].pT(),
                    pythia.event[detected_particles_negative_number[j]].eta(),
                    pythia.event[detected_particles_negative_number[j]].phi(),
                    particleMass[Pion]);
                //pair 4momenta using and filling
                mass = (temp4VectorPositive+temp4VectorNegative).M();
                motherPositive = pythia.event[detected_particles_positive_number[i]].mother1();
                motherNegative = pythia.event[detected_particles_negative_number[j]].mother1();
                t1 = std::find(mother_particles_number.begin(), mother_particles_number.end(), motherPositive)!=mother_particles_number.end();
                t2 = std::find(mother_particles_number.begin(), mother_particles_number.end(), motherNegative)!=mother_particles_number.end();
                t3 = motherNegative==motherPositive;
                t4 = pythia.event[detected_particles_positive_number[i]].mother2()==0&&pythia.event[detected_particles_negative_number[j]].mother2()==0;
                if(t1&&t2&&t3&&t4){
                    MKpiResonant.Fill(mass);
                    MKpiParticles.Fill(pythia.particleData.name(pythia.event[motherPositive].idAbs()).c_str(), mass, 1.);
                } else{
                    MKpiNonresonant.Fill(mass);
                }
                MKpiSum.Fill(mass);
                //pair 4momenta setting
                temp4VectorPositive.SetPtEtaPhiM(pythia.event[detected_particles_positive_number[i]].pT(),
                    pythia.event[detected_particles_positive_number[i]].eta(),
                    pythia.event[detected_particles_positive_number[i]].phi(),
                    particleMass[Pion]);
                temp4VectorNegative.SetPtEtaPhiM(pythia.event[detected_particles_negative_number[j]].pT(),
                    pythia.event[detected_particles_negative_number[j]].eta(),
                    pythia.event[detected_particles_negative_number[j]].phi(),
                    particleMass[Kaon]);
                //pair 4momenta using and filling
                mass = (temp4VectorPositive+temp4VectorNegative).M();
                motherPositive = pythia.event[detected_particles_positive_number[i]].mother1();
                motherNegative = pythia.event[detected_particles_negative_number[j]].mother1();
                t1 = std::find(mother_particles_number.begin(), mother_particles_number.end(), motherPositive)!=mother_particles_number.end();
                t2 = std::find(mother_particles_number.begin(), mother_particles_number.end(), motherNegative)!=mother_particles_number.end();
                t3 = motherNegative==motherPositive;
                t4 = pythia.event[detected_particles_positive_number[i]].mother2()==0&&pythia.event[detected_particles_negative_number[j]].mother2()==0;
                if(t1&&t2&&t3&&t4){
                    MKpiResonant.Fill(mass);
                    MKpiParticles.Fill(pythia.particleData.name(pythia.event[motherPositive].idAbs()).c_str(), mass, 1.);
                } else{
                    MKpiNonresonant.Fill(mass);
                }
                MKpiSum.Fill(mass);

                // assuming proper cuts for kaons and protons (the rest is, you guessed it, assumed to be pions)
                //pair 4momenta setting
                int massTag = Pion;
                if(pythia.event[detected_particles_positive_number[i]].p().pAbs()<0.4&&pythia.event[detected_particles_positive_number[i]].idAbs()==321){
                    massTag = Kaon;
                } else if(pythia.event[detected_particles_positive_number[i]].p().pAbs()<0.8&&pythia.event[detected_particles_positive_number[i]].idAbs()==2212){
                    massTag = Proton;
                }
                temp4VectorPositive.SetPtEtaPhiM(pythia.event[detected_particles_positive_number[i]].pT(),
                    pythia.event[detected_particles_positive_number[i]].eta(),
                    pythia.event[detected_particles_positive_number[i]].phi(),
                    particleMass[massTag]);
                massTag = Pion;
                if(pythia.event[detected_particles_negative_number[j]].p().pAbs()<0.4&&pythia.event[detected_particles_negative_number[j]].idAbs()==321){
                    massTag = Kaon;
                } else if(pythia.event[detected_particles_negative_number[j]].p().pAbs()<0.8&&pythia.event[detected_particles_negative_number[j]].idAbs()==2212){
                    massTag = Proton;
                }
                temp4VectorNegative.SetPtEtaPhiM(pythia.event[detected_particles_negative_number[j]].pT(),
                    pythia.event[detected_particles_negative_number[j]].eta(),
                    pythia.event[detected_particles_negative_number[j]].phi(),
                    particleMass[massTag]);
                //pair 4momenta using and filling
                mass = (temp4VectorPositive+temp4VectorNegative).M();
                motherPositive = pythia.event[detected_particles_positive_number[i]].mother1();
                motherNegative = pythia.event[detected_particles_negative_number[j]].mother1();
                t1 = std::find(mother_particles_number.begin(), mother_particles_number.end(), motherPositive)!=mother_particles_number.end();
                t2 = std::find(mother_particles_number.begin(), mother_particles_number.end(), motherNegative)!=mother_particles_number.end();
                t3 = motherNegative==motherPositive;
                t4 = pythia.event[detected_particles_positive_number[i]].mother2()==0&&pythia.event[detected_particles_negative_number[j]].mother2()==0;
                if(t1&&t2&&t3&&t4){
                    MPIDResonant.Fill(mass);
                } else{
                    MPIDNonresonant.Fill(mass);
                }
                MPIDSum.Fill(mass);
            }
        }

        // part for making a histogram of what is detected (positive sign pairs)
        for(size_t i = 0; i<detected_particles_positive_number.size(); i++){
            for(size_t j = i+1; j<detected_particles_positive_number.size(); j++){
                // assuming everything is a kaon
                //here "negative" means nothing
                //pair 4momenta setting
                temp4VectorPositive.SetPtEtaPhiM(pythia.event[detected_particles_positive_number[i]].pT(),
                    pythia.event[detected_particles_positive_number[i]].eta(),
                    pythia.event[detected_particles_positive_number[i]].phi(),
                    particleMass[Kaon]);
                temp4VectorNegative.SetPtEtaPhiM(pythia.event[detected_particles_positive_number[j]].pT(),
                    pythia.event[detected_particles_positive_number[j]].eta(),
                    pythia.event[detected_particles_positive_number[j]].phi(),
                    particleMass[Kaon]);
                //pair 4momenta using and filling
                double mass = (temp4VectorPositive+temp4VectorNegative).M();
                int motherPositive = pythia.event[detected_particles_positive_number[i]].mother1();
                int motherNegative = pythia.event[detected_particles_positive_number[j]].mother1();
                bool t1 = std::find(mother_particles_number.begin(), mother_particles_number.end(), motherPositive)!=mother_particles_number.end();
                bool t2 = std::find(mother_particles_number.begin(), mother_particles_number.end(), motherNegative)!=mother_particles_number.end();
                bool t3 = motherNegative==motherPositive;
                bool t4 = pythia.event[detected_particles_positive_number[i]].mother2()==0&&pythia.event[detected_particles_positive_number[j]].mother2()==0;
                if(t1&&t2&&t3&&t4){
                    MKKBackgroundPositiveResonant.Fill(mass);
                } else{
                    MKKBackgroundPositiveNonresonant.Fill(mass);
                }
                MKKBackgroundPositiveSum.Fill(mass);
            }
        }

        // part for making a histogram of what is detected (negative sign pairs)
        for(size_t i = 0; i<detected_particles_negative_number.size(); i++){
            for(size_t j = i+1; j<detected_particles_negative_number.size(); j++){
                // assuming everything is a kaon
                //here "positive" means nothing
                //pair 4momenta setting
                temp4VectorPositive.SetPtEtaPhiM(pythia.event[detected_particles_negative_number[i]].pT(),
                    pythia.event[detected_particles_negative_number[i]].eta(),
                    pythia.event[detected_particles_negative_number[i]].phi(),
                    particleMass[Kaon]);
                temp4VectorNegative.SetPtEtaPhiM(pythia.event[detected_particles_negative_number[j]].pT(),
                    pythia.event[detected_particles_negative_number[j]].eta(),
                    pythia.event[detected_particles_negative_number[j]].phi(),
                    particleMass[Kaon]);
                //pair 4momenta using and filling
                double mass = (temp4VectorPositive+temp4VectorNegative).M();
                int motherPositive = pythia.event[detected_particles_negative_number[i]].mother1();
                int motherNegative = pythia.event[detected_particles_negative_number[j]].mother1();
                bool t1 = std::find(mother_particles_number.begin(), mother_particles_number.end(), motherPositive)!=mother_particles_number.end();
                bool t2 = std::find(mother_particles_number.begin(), mother_particles_number.end(), motherNegative)!=mother_particles_number.end();
                bool t3 = motherNegative==motherPositive;
                bool t4 = pythia.event[detected_particles_negative_number[i]].mother2()==0&&pythia.event[detected_particles_negative_number[j]].mother2()==0;
                if(t1&&t2&&t3&&t4){
                    MKKBackgroundNegativeResonant.Fill(mass);
                } else{
                    MKKBackgroundNegativeNonresonant.Fill(mass);
                }
                MKKBackgroundNegativeSum.Fill(mass);
            }
        }

        //total ratio
        n_events++;
    }

    //closing touches
    ParticlesDetected.LabelsDeflate();
    ParticlesDetected.SetMinimum(0);
    ParticlesDetected.LabelsOption("a", "X");
    ParticlesReconstructedTwoDaughters.LabelsDeflate();
    ParticlesReconstructedTwoDaughters.SetMinimum(0);
    ParticlesReconstructedTwoDaughters.LabelsOption("a", "X");
    ParticlesReconstructedMultipleDaughters.LabelsDeflate();
    ParticlesReconstructedMultipleDaughters.SetMinimum(0);
    ParticlesReconstructedMultipleDaughters.LabelsOption("a", "X");

    //writing to file
    outFile->cd();
    outFile->Write();
    outFile->Close();

    //pairs assumed as pions
    MpipiParticles.LabelsDeflate();
    THStack MpipiStack("MpipiStack", "#pi#pi pairs");
    MpipiNonresonant.SetName("Nonresonant or mismatched");
    MpipiNonresonant.SetFillColor(1);
    MpipiNonresonant.SetLineWidth(0);
    MpipiStack.Add(&MpipiNonresonant);
    std::vector<TH1D*> projections;
    for(Int_t i = 1; i<MpipiParticles.GetXaxis()->GetNbins()+1; i++){
        projections.push_back(new TH1D(*MpipiParticles.ProjectionY("", i, i)));
        projections[i-1]->SetName(MpipiParticles.GetXaxis()->GetBinLabel(i));
        projections[i-1]->SetFillColor(i+1+(i>=9)); //color 10 is white, cant have that
        projections[i-1]->SetLineWidth(0);
        MpipiStack.Add(projections[i-1]);
    }
    TCanvas canvas("Stack", "Stack", 4000, 2400);
    MpipiStack.Draw();
    canvas.BuildLegend(0.6, 0.5, 0.89, 0.89, "Particles reconstructed", "F");
    canvas.SaveAs(string(filename).replace(string(filename).find(".root"), 5, "Mpipi.pdf").c_str());
    //2nd version with interesting space zoomed into
    canvas.Clear();
    MpipiStack.GetXaxis()->SetRangeUser(0.2, 1.2);
    MpipiStack.Draw();
    canvas.BuildLegend(0.7, 0.6, 0.89, 0.89, "Particles reconstructed", "F");
    canvas.SaveAs(string(filename).replace(string(filename).find(".root"), 5, "Mpipiv2.pdf").c_str());
    for(size_t i = 0; i<projections.size(); i++){
        delete projections[i];
    }

    //pairs assumed as kaons
    MKKParticles.LabelsDeflate();
    THStack MKKStack("MKKStack", "KK pairs");
    MKKNonresonant.SetName("Nonresonant or mismatched");
    MKKNonresonant.SetFillColor(1);
    MKKNonresonant.SetLineWidth(0);
    MKKStack.Add(&MKKNonresonant);
    projections.clear();
    for(Int_t i = 1; i<MKKParticles.GetXaxis()->GetNbins()+1; i++){
        projections.push_back(new TH1D(*MKKParticles.ProjectionY("", i, i)));
        projections[i-1]->SetName(MKKParticles.GetXaxis()->GetBinLabel(i));
        projections[i-1]->SetFillColor(i+1+(i>=9)); //color 10 is white, cant have that
        projections[i-1]->SetLineWidth(0);
        MKKStack.Add(projections[i-1]);
    }
    canvas.Clear();
    MKKStack.Draw();
    canvas.BuildLegend(0.6, 0.5, 0.89, 0.89, "Particles reconstructed", "F");
    canvas.SaveAs(string(filename).replace(string(filename).find(".root"), 5, "MKK.pdf").c_str());
    //2nd version with interesting space zoomed into
    canvas.Clear();
    MKKStack.GetXaxis()->SetRangeUser(0.9, 1.7);
    MKKStack.Draw();
    canvas.BuildLegend(0.7, 0.6, 0.89, 0.89, "Particles reconstructed", "F");
    canvas.SaveAs(string(filename).replace(string(filename).find(".root"), 5, "MKKv2.pdf").c_str());
    for(size_t i = 0; i<projections.size(); i++){
        delete projections[i];
    }

    //pairs assumed as kaon-pion pairs
    MKpiParticles.LabelsDeflate();
    THStack MKpiStack("MKpiStack", "K#pi pairs");
    MKpiNonresonant.SetName("Nonresonant or mismatched");
    MKpiNonresonant.SetFillColor(1);
    MKpiNonresonant.SetLineWidth(0);
    MKpiStack.Add(&MKpiNonresonant);
    projections.clear();
    for(Int_t i = 1; i<MKpiParticles.GetXaxis()->GetNbins()+1; i++){
        projections.push_back(new TH1D(*MKpiParticles.ProjectionY("", i, i)));
        projections[i-1]->SetName(MKpiParticles.GetXaxis()->GetBinLabel(i));
        projections[i-1]->SetFillColor(i+1+(i>=9)); //color 10 is white, cant have that
        projections[i-1]->SetLineWidth(0);
        MKpiStack.Add(projections[i-1]);
    }
    canvas.Clear();
    MKpiStack.Draw();
    canvas.BuildLegend(0.6, 0.5, 0.89, 0.89, "Particles reconstructed", "F");
    canvas.SaveAs(string(filename).replace(string(filename).find(".root"), 5, "MKpi.pdf").c_str());
    //2nd version with interesting space zoomed into
    canvas.Clear();
    MKpiStack.GetXaxis()->SetRangeUser(0.6, 1.2);
    MKpiStack.Draw();
    canvas.BuildLegend(0.7, 0.6, 0.89, 0.89, "Particles reconstructed", "F");
    canvas.SaveAs(string(filename).replace(string(filename).find(".root"), 5, "MKpiv2.pdf").c_str());
    for(size_t i = 0; i<projections.size(); i++){
        delete projections[i];
    }

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

bool IsInMeasurementSpecificCuts(Particle* particle){
    //|eta|<0.7
    //production vertex closer to  primary vertex than 30mm
    return abs(particle->eta())<0.7&&particle->vProd().pT()<30;
}