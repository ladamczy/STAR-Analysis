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
bool IsParticleQuarkOrGluon(Particle);
void FindGoodMothers(Event *, int, std::vector<int> &);

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

    //histogram variables
    int n_ptBins = 9;
    double ptBins[] = { 0,0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8 };
    int n_etaBins = 10;
    double etaBins[] = { -1., -0.8, -0.6, -0.4, -0.2, 0., 0.2, 0.4, 0.6, 0.8, 1.0 };
    int stateBins = 40;
    double stateMassMax = 40;

    //histograms
    //general
    TH1D ParticlesDetected("ParticlesDetected", "ParticlesDetected", 20000, -10000, 10000);
    TH1D ParticlesReconstructed("ParticlesReconstructed", "ParticlesReconstructed", 20000, -10000, 10000);
    THStack MpipiStack("MpipiStack", "#pion mass assumed");
    //Deltas
    TH2D DeltappDaughters("DeltappDaughters", "DeltappDaughters", 1, 0, 1, 1, 0, 1);
    TH2D Delta0Daughters("Delta0Daughters", "Delta0Daughters", 1, 0, 1, 1, 0, 1);
    TH2D Deltabar0Daughters("Deltabar0Daughters", "Deltabar0Daughters", 1, 0, 1, 1, 0, 1);
    TH2D DeltabarmmDaughters("DeltabarmmDaughters", "DeltabarmmDaughters", 1, 0, 1, 1, 0, 1);
    TH2D DeltappMothers("DeltappMothers", "DeltappMothers", 1, 0, 1, 1, 0, 1);
    TH2D Delta0Mothers("Delta0Mothers", "Delta0Mothers", 1, 0, 1, 1, 0, 1);
    TH2D Deltabar0Mothers("Deltabar0Mothers", "Deltabar0Mothers", 1, 0, 1, 1, 0, 1);
    TH2D DeltabarmmMothers("DeltabarmmMothers", "DeltabarmmMothers", 1, 0, 1, 1, 0, 1);
    //K0s
    TH1D K0ptHist("K0ptHist", "K0ptHist", n_ptBins, ptBins);
    TH1D K0etaHist("K0etaHist", "K0etaHist", n_etaBins, etaBins);
    TH1D K0multiplicity("K0multiplicity", "K^{0}_{S} multiplicity", 10, 0, 10);
    TH1D K0AdditionalTracks("K0AdditionalTracks", "Additional tracks except K^{0}_{S} daughters", 20, 0, 20);
    TH2D K0multiplicityXmass("K0multiplicityXmass", "K^{0}_{S} multiplicity with respect to mass of a central state;m_{X} [GeV]", stateBins, 0, stateMassMax, 10, 0, 10);
    TEfficiency K0ratioXmass("K0ratioXmass", "Ratio of detected K^{0}_{S} to all measured events with respect to mass of a central state;m_{X} [GeV]", stateBins, 0, stateMassMax);
    TH2D K0AdditionalTracksXmass("K0AdditionalTracksXmass", "Additional tracks except K^{0}_{S} daughters;m_{X} [GeV]", stateBins, 0, stateMassMax, 20, 0, 20);
    //Lambda0
    TH1D LambdaptHist("LambdaptHist", "LambdaptHist", n_ptBins, ptBins);
    TH1D LambdaetaHist("LambdaetaHist", "LambdaetaHist", n_etaBins, etaBins);
    TH1D Lambdamultiplicity("Lambdamultiplicity", "#Lambda^{0} multiplicity", 10, 0, 10);
    TH1D LambdaAdditionalTracks("LambdaAdditionalTracks", "Additional tracks except #Lambda^{0} daughters", 20, 0, 20);
    TH2D LambdamultiplicityXmass("LambdamultiplicityXmass", "#Lambda^{0} multiplicity with respect to mass of a central state;m_{X} [GeV]", stateBins, 0, stateMassMax, 10, 0, 10);
    TEfficiency LambdaratioXmass("LambdaratioXmass", "Ratio of detected #Lambda^{0} to all measured events with respect to mass of a central state;m_{X} [GeV]", stateBins, 0, stateMassMax);
    TH2D LambdaAdditionalTracksXmass("LambdaAdditionalTracksXmass", "Additional tracks except #Lambda^{0} daughters;m_{X} [GeV]", stateBins, 0, stateMassMax, 20, 0, 20);
    //LambdaBar0
    TH1D LambdaBarptHist("LambdaBarptHist", "LambdaBarptHist", n_ptBins, ptBins);
    TH1D LambdaBaretaHist("LambdaBaretaHist", "LambdaBaretaHist", n_etaBins, etaBins);
    TH1D LambdaBarmultiplicity("LambdaBarmultiplicity", "#Lambda^{0} multiplicity", 10, 0, 10);
    TH1D LambdaBarAdditionalTracks("LambdaBarAdditionalTracks", "Additional tracks except #Lambda^{0} daughters", 20, 0, 20);
    TH2D LambdaBarmultiplicityXmass("LambdaBarmultiplicityXmass", "#Lambda^{0} multiplicity with respect to mass of a central state;m_{X} [GeV]", stateBins, 0, stateMassMax, 10, 0, 10);
    TEfficiency LambdaBarratioXmass("LambdaBarratioXmass", "Ratio of detected #Lambda^{0} to all measured events with respect to mass of a central state;m_{X} [GeV]", stateBins, 0, stateMassMax);
    TH2D LambdaBarAdditionalTracksXmass("LambdaBarAdditionalTracksXmass", "Additional tracks except #Lambda^{0} daughters;m_{X} [GeV]", stateBins, 0, stateMassMax, 20, 0, 20);
    //Technical
    TH1D MpipiNonresonant("MpipiNonresonant", ";m_{#pi#pi}", 500, 0, 5);
    // MpipiNonresonant.SetDirectory(nullptr);
    MpipiNonresonant.SetFillColor(kBlue);
    TH1D MpipiResonant("MpipiResonant", ";m_{#pi#pi}", 500, 0, 5);
    // MpipiResonant.SetDirectory(nullptr);
    MpipiResonant.SetFillColor(kRed);

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
    int n_LambdaBar = 0;
    std::vector<int> detected_particles_number;
    std::vector<int> detected_particles_positive_number;
    std::vector<int> detected_particles_negative_number;
    std::vector<int> all_daughters_detected_resonant_number;
    std::vector<int> all_daughters_detected_nonresonant_number;
    Rndm generator(0);
    TLorentzVector temp4Vector = { 0,0,0,0 };
    TLorentzVector temp4VectorPositive;
    TLorentzVector temp4VectorNegative;
    bool AllDaughtersDetected;
    int DaughtersDetected;
    std::vector<int> mother_particles_number;

    // Begin event loop. Generate event; skip if generation aborted.
    for(int iEvent = 0; iEvent<nEvents; ++iEvent){
        if(!pythia.next()) continue;
        //cleaning up
        detected_particles_number.clear();
        detected_particles_positive_number.clear();
        detected_particles_negative_number.clear();
        all_daughters_detected_resonant_number.clear();
        all_daughters_detected_nonresonant_number.clear();
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
        //with smulated TPC acceptance
        int TPCParticles = 0;
        for(int i = 0; i<pythia.event.size(); i++){
            // added current bound (0.7 instead of 0.9)
            if(isParticleInTPCAcceptance(pythia.event[i])&&isParticleDetected(&pythia.event[i])&&abs(pythia.event[i].eta())<0.7){
                TPCParticles++;
                detected_particles_number.push_back(i);
            }
        }
        if(TPCParticles<2){
            continue;
        }

        //number of K0 & Lambda
        int nK0 = 0;
        int nLambda = 0;
        int nLambdaBar = 0;
        std::vector<int> K0indices;
        std::vector<int> Lambdaindices;
        std::vector<int> LambdaBarindices;
        int tempDaughter1, tempDaughter2;
        std::vector<Pythia8::Particle> acceptedParticlesPlusDaughters;
        temp4Vector.SetE(mEnergy);
        double XStateMass = (temp4Vector-p1-p2).M();
        //gathering info in main loop
        for(int part_index = 0; part_index<pythia.event.size(); part_index++){
            if(abs(pythia.event[part_index].id())==K0sPDGid){ K0indices.push_back(part_index); }
            if(pythia.event[part_index].id()==LambdaPDGid){ Lambdaindices.push_back(part_index); }
            if(pythia.event[part_index].id()==LambdabarPDGid){ LambdaBarindices.push_back(part_index); }
            DaughtersDetected = 0;
            AllDaughtersDetected = true;
            //check if all daugters of part_index particle are detected
            for(size_t daughter = 0; daughter<pythia.event[part_index].daughterList().size();daughter++){
                if(std::find(detected_particles_number.begin(), detected_particles_number.end(), pythia.event[part_index].daughterList()[daughter])!=detected_particles_number.end()){
                    DaughtersDetected++;
                } else{
                    AllDaughtersDetected = false;
                }
            }
            //if all particles are detected, we check if we reconstruct a "viable" particle (not a gluon or quark)
            if(AllDaughtersDetected){
                mother_particles_number.clear();
                if(strcmp(pythia.event[part_index].name().c_str(), "g")==0){
                    printf("\nParticle checked is a gluon\n");
                    FindGoodMothers(&pythia.event, pythia.event[part_index].daughter1(), mother_particles_number);
                } else if(strcmp(pythia.event[part_index].name().c_str(), "d")==0){
                    printf("\nParticle checked is a d quark\n");
                    FindGoodMothers(&pythia.event, pythia.event[part_index].daughter1(), mother_particles_number);
                } else if(strcmp(pythia.event[part_index].name().c_str(), "dbar")==0){
                    printf("\nParticle checked is a dbar quark\n");
                    FindGoodMothers(&pythia.event, pythia.event[part_index].daughter1(), mother_particles_number);
                } else if(strcmp(pythia.event[part_index].name().c_str(), "u")==0){
                    printf("\nParticle checked is a u quark\n");
                    FindGoodMothers(&pythia.event, pythia.event[part_index].daughter1(), mother_particles_number);
                } else if(strcmp(pythia.event[part_index].name().c_str(), "ubar")==0){
                    printf("\nParticle checked is a ubar quark\n");
                    FindGoodMothers(&pythia.event, pythia.event[part_index].daughter1(), mother_particles_number);
                } else if(strcmp(pythia.event[part_index].name().c_str(), "s")==0){
                    printf("\nParticle checked is a s quark\n");
                    FindGoodMothers(&pythia.event, pythia.event[part_index].daughter1(), mother_particles_number);
                } else if(strcmp(pythia.event[part_index].name().c_str(), "sbar")==0){
                    printf("\nParticle checked is a sbar quark\n");
                    FindGoodMothers(&pythia.event, pythia.event[part_index].daughter1(), mother_particles_number);
                } else if(strcmp(pythia.event[part_index].name().c_str(), "Delta++")==0){
                    DeltappDaughters.Fill(pythia.event[pythia.event[part_index].daughter1()].name().c_str(), pythia.event[pythia.event[part_index].daughter2()].name().c_str(), 1.);
                    DeltappMothers.Fill(pythia.event[pythia.event[part_index].mother1()].name().c_str(), pythia.event[pythia.event[part_index].mother2()].name().c_str(), 1.);
                } else if(strcmp(pythia.event[part_index].name().c_str(), "Delta0")==0){
                    Delta0Daughters.Fill(pythia.event[pythia.event[part_index].daughter1()].name().c_str(), pythia.event[pythia.event[part_index].daughter2()].name().c_str(), 1.);
                    Delta0Mothers.Fill(pythia.event[pythia.event[part_index].mother1()].name().c_str(), pythia.event[pythia.event[part_index].mother2()].name().c_str(), 1.);
                } else if(strcmp(pythia.event[part_index].name().c_str(), "Deltabar0")==0){
                    Deltabar0Daughters.Fill(pythia.event[pythia.event[part_index].daughter1()].name().c_str(), pythia.event[pythia.event[part_index].daughter2()].name().c_str(), 1.);
                    Deltabar0Mothers.Fill(pythia.event[pythia.event[part_index].mother1()].name().c_str(), pythia.event[pythia.event[part_index].mother2()].name().c_str(), 1.);
                } else if(strcmp(pythia.event[part_index].name().c_str(), "Deltabar--")==0){
                    DeltabarmmDaughters.Fill(pythia.event[pythia.event[part_index].daughter1()].name().c_str(), pythia.event[pythia.event[part_index].daughter2()].name().c_str(), 1.);
                    DeltabarmmMothers.Fill(pythia.event[pythia.event[part_index].mother1()].name().c_str(), pythia.event[pythia.event[part_index].mother2()].name().c_str(), 1.);
                }
                //final checks to fill histograms
                //this way every particle gets filled exactly once
                //because only particles with TPC-detected decay products
                //get to this point of a loop
                //due to AllDaughtersDetected==true
                //so even if there was a copy earlier in the decay tree, it won't get registered
                if(!IsParticleQuarkOrGluon(pythia.event[part_index])){
                    ParticlesReconstructed.Fill(pythia.event[part_index].name().c_str(), 1.);
                    all_daughters_detected_resonant_number.push_back(part_index);
                } else{
                    //removing potential particle duplicates
                    sort(mother_particles_number.begin(), mother_particles_number.end());
                    mother_particles_number.erase(unique(mother_particles_number.begin(), mother_particles_number.end()), mother_particles_number.end());
                    //checking if there is something else than Pomerons
                    if(mother_particles_number.size()==1&&abs(pythia.event[mother_particles_number[0]].id())==990){
                        ParticlesReconstructed.Fill("#splitline{Single}{Pomeron}", 1.);
                        all_daughters_detected_nonresonant_number.push_back(part_index);
                    } else if(mother_particles_number.size()==2&&abs(pythia.event[mother_particles_number[0]].id())==990&&abs(pythia.event[mother_particles_number[1]].id())==990){
                        ParticlesReconstructed.Fill("#splitline{Double}{Pomeron}", 1.);
                        all_daughters_detected_nonresonant_number.push_back(part_index);
                    } else{
                        ParticlesReconstructed.Fill("other", 1.);
                        for(auto &&parr:mother_particles_number){
                            printf("%d: %s\n", parr, pythia.event[parr].name().c_str());
                        }
                        printf("\n");
                        pythia.event.list(true, true);
                    }
                }
            }
        }

        //gathering info in detection loop
        for(size_t part_index = 0; part_index<detected_particles_number.size(); part_index++){
            ParticlesDetected.Fill(pythia.event[detected_particles_number[part_index]].name().c_str(), 1.);
            if(pythia.event[detected_particles_number[part_index]].charge()>0){
                detected_particles_positive_number.push_back(part_index);
            } else{
                detected_particles_negative_number.push_back(part_index);
            }
        }

        //part for making a histogram of what is detected
        //assuming everything is a pion
        for(size_t i = 0; i<detected_particles_positive_number.size(); i++){
            for(size_t j = 0; j<detected_particles_negative_number.size(); j++){
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
                bool t1 = std::find(all_daughters_detected_resonant_number.begin(), all_daughters_detected_resonant_number.end(), motherPositive)!=all_daughters_detected_resonant_number.end();
                bool t2 = std::find(all_daughters_detected_resonant_number.begin(), all_daughters_detected_resonant_number.end(), motherNegative)!=all_daughters_detected_resonant_number.end();
                bool t3 = motherNegative==motherPositive;
                if(t1&&t2&&t3){
                    MpipiResonant.Fill(mass);
                } else{
                    MpipiNonresonant.Fill(mass);
                }
            }
        }
        //assuming everything is a kaon
        //assuming proper cuts for kaons and protons (the rest is, you guessed it, assumed to be pions)

        //tests for other particles
        //K0s
        for(long unsigned int K0candidate = 0; K0candidate<K0indices.size(); K0candidate++){
            tempDaughter1 = pythia.event[K0indices[K0candidate]].daughter1();
            tempDaughter2 = pythia.event[K0indices[K0candidate]].daughter2();
            if((std::find(detected_particles_number.begin(), detected_particles_number.end(), tempDaughter1)!=detected_particles_number.end())&&(std::find(detected_particles_number.begin(), detected_particles_number.end(), tempDaughter2)!=detected_particles_number.end())&&(tempDaughter1+1==tempDaughter2)){
                nK0++;
                K0ptHist.Fill(pythia.event[K0indices[K0candidate]].pT());
                K0etaHist.Fill(pythia.event[K0indices[K0candidate]].eta());
            }
        }
        //Lambda0
        for(long unsigned int Lambdacandidate = 0; Lambdacandidate<Lambdaindices.size(); Lambdacandidate++){
            tempDaughter1 = pythia.event[Lambdaindices[Lambdacandidate]].daughter1();
            tempDaughter2 = pythia.event[Lambdaindices[Lambdacandidate]].daughter2();
            if((std::find(detected_particles_number.begin(), detected_particles_number.end(), tempDaughter1)!=detected_particles_number.end())&&(std::find(detected_particles_number.begin(), detected_particles_number.end(), tempDaughter2)!=detected_particles_number.end())&&(tempDaughter1+1==tempDaughter2)){
                nLambda++;
                LambdaptHist.Fill(pythia.event[Lambdaindices[Lambdacandidate]].pT());
                LambdaetaHist.Fill(pythia.event[Lambdaindices[Lambdacandidate]].eta());
            }
        }
        //LambdaBar0
        for(long unsigned int LambdaBarcandidate = 0; LambdaBarcandidate<LambdaBarindices.size(); LambdaBarcandidate++){
            tempDaughter1 = pythia.event[LambdaBarindices[LambdaBarcandidate]].daughter1();
            tempDaughter2 = pythia.event[LambdaBarindices[LambdaBarcandidate]].daughter2();
            if((std::find(detected_particles_number.begin(), detected_particles_number.end(), tempDaughter1)!=detected_particles_number.end())&&(std::find(detected_particles_number.begin(), detected_particles_number.end(), tempDaughter2)!=detected_particles_number.end())&&(tempDaughter1+1==tempDaughter2)){
                nLambdaBar++;
                LambdaBarptHist.Fill(pythia.event[LambdaBarindices[LambdaBarcandidate]].pT());
                LambdaBaretaHist.Fill(pythia.event[LambdaBarindices[LambdaBarcandidate]].eta());
            }
        }

        //final histograms
        //K0s
        if(nK0>0){
            K0multiplicity.Fill(nK0);
            K0multiplicityXmass.Fill(XStateMass, nK0);
            K0ratioXmass.Fill(nK0>0, XStateMass);
            K0AdditionalTracks.Fill(detected_particles_number.size()-nK0*2);
            K0AdditionalTracksXmass.Fill(XStateMass, detected_particles_number.size()-nK0*2);
        }
        //Lambda0
        if(nLambda>0){
            Lambdamultiplicity.Fill(nLambda);
            LambdamultiplicityXmass.Fill(XStateMass, nLambda);
            LambdaratioXmass.Fill(nLambda>0, XStateMass);
            LambdaAdditionalTracks.Fill(detected_particles_number.size()-nLambda*2);
            LambdaAdditionalTracksXmass.Fill(XStateMass, detected_particles_number.size()-nLambda*2);
        }
        //LambdaBar0
        if(nLambdaBar>0){
            LambdaBarmultiplicity.Fill(nLambdaBar);
            LambdaBarmultiplicityXmass.Fill(XStateMass, nLambdaBar);
            LambdaBarratioXmass.Fill(nLambdaBar>0, XStateMass);
            LambdaBarAdditionalTracks.Fill(detected_particles_number.size()-nLambdaBar*2);
            LambdaBarAdditionalTracksXmass.Fill(XStateMass, detected_particles_number.size()-nLambdaBar*2);
        }

        //total ratio
        n_events++;
        n_K0 += nK0;
        n_Lambda += nLambda;
        n_LambdaBar += nLambdaBar;
    }

    //closing touches
    ParticlesDetected.LabelsDeflate();
    ParticlesDetected.SetMinimum(0);
    ParticlesDetected.LabelsOption("a", "X");
    ParticlesReconstructed.LabelsDeflate();
    ParticlesReconstructed.SetMinimum(0);
    ParticlesReconstructed.LabelsOption("a", "X");
    DeltappDaughters.LabelsDeflate();
    Delta0Daughters.LabelsDeflate();
    Deltabar0Daughters.LabelsDeflate();
    DeltabarmmDaughters.LabelsDeflate();
    DeltappMothers.LabelsDeflate();
    Delta0Mothers.LabelsDeflate();
    Deltabar0Mothers.LabelsDeflate();
    DeltabarmmMothers.LabelsDeflate();
    MpipiStack.Add(&MpipiNonresonant);
    MpipiStack.Add(&MpipiResonant);

    //writing to file
    outFile->cd();
    outFile->Write();
    outFile->Close();

    // Statistics on event generation.
    pythia.stat();

    //other things
    std::cout<<"N_events: "<<n_events<<endl;
    std::cout<<"N_K0: "<<n_K0<<endl;
    std::cout<<"N_Lambda: "<<n_Lambda<<endl;
    std::cout<<"N_LambdaBar: "<<n_LambdaBar<<endl;
    std::stringstream ratio_output;
    ratio_output<<std::fixed<<std::setprecision(6)<<double(n_K0)/n_events;
    std::cout<<"N_K0/N_events: "<<ratio_output.str()<<endl;
    ratio_output.str(""); //clearing
    ratio_output<<std::fixed<<std::setprecision(6)<<double(n_Lambda)/n_events;
    std::cout<<"N_Lambda0/N_events: "<<ratio_output.str()<<endl;
    ratio_output.str(""); //clearing
    ratio_output<<std::fixed<<std::setprecision(6)<<double(n_LambdaBar)/n_events;
    std::cout<<"N_LambdaBar0/N_events: "<<ratio_output.str()<<endl;

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

bool IsParticleQuarkOrGluon(Particle particle){
    if(abs(particle.id())<=9 or abs(particle.id())==21){
        return true;
    }
    return false;
}

void FindGoodMothers(Event *event, int particle, std::vector<int> &particle_set){
    //put all mothers of a particle on a list
    for(size_t i = 0; i<event->at(particle).motherList().size(); i++){
        particle_set.push_back(event->at(particle).motherList()[i]);
        //if particle is a quark or a gluon
        //recursively do the same
        if(IsParticleQuarkOrGluon(event->at(event->at(particle).motherList()[i]))){
            FindGoodMothers(event, event->at(particle).motherList()[i], particle_set);
        }
    }
    //and then, if the particle is on the list
    //remove it
    if(find(particle_set.begin(), particle_set.end(), particle)!=particle_set.end()&&IsParticleQuarkOrGluon(event->at(particle))){
        particle_set.erase(find(particle_set.begin(), particle_set.end(), particle));
    }
}