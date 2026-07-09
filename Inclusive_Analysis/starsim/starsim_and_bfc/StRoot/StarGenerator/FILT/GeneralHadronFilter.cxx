#include "GeneralHadronFilter.h"
#include "StarGenerator/EVENT/StarGenParticle.h"
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <memory>
#include "DllImport.h"

ClassImp(GeneralHadronFilter);

GeneralHadronFilter::GeneralHadronFilter():StarFilterMaker("generalhadronfilter"){
    this->SetNameTitle("generalhadronfilter", "Filter responsible for accepting generally all particles, with stipulation that at least some have to be in TPC fiducial region");
};
GeneralHadronFilter::~GeneralHadronFilter(){};

Int_t GeneralHadronFilter::Filter(StarGenEvent* event){
    //loop through all particles
    //either we find a particle that decays in Geant, and then there is no telling if it decays properly or not
    //or we find at least two particles roughly in TPC fiducial region
    bool hasAlreadyOneParticleInTPCFiducialRegion = false;

    StarGenParticle* particle = nullptr;
    for(size_t i = 0; i<event->GetNumberOfParticles(); i++){
        particle = (*event)[i];

        //scattered proton fiducial region special thing (currently off)
        //ID 3 and 4 because (excerpt from PYTHIA production):
        // [   0|   0|  -1] id=         0    Rootino stat=-201 p=(   0.000,   0.000,   0.000,   0.000;  510.000) v=(  0.0000,  0.0000,   0.000) [2212 2212] [0 15]
        // [   1|   1|  -1] id=      2212     proton stat=04 p=(  -0.241,   0.112, 254.998, 255.000;    0.938) v=( -0.1651, -0.1171,  46.621) [0 0] [3 5]
        // [   2|   2|  -1] id=      2212     proton stat=04 p=(   0.241,  -0.112,-254.998, 255.000;    0.938) v=( -0.1651, -0.1171,  46.621) [0 0] [3 5]
        // [   3|   3|  -1] id=      2212     proton stat=01 p=(  -0.360,   0.174, 254.969, 254.971;    0.938) v=( -0.1651, -0.1171,  46.621) [1 2] [0 0]
        // [   4|   4|  -1] id=      2212     proton stat=01 p=(  -0.151,  -0.201,-213.968, 213.970;    0.938) v=( -0.1651, -0.1171,  46.621) [1 2] [0 0]
        // [   5|   5|  -1] id=   9900110            stat=15 p=(   0.511,   0.028, -41.001,  41.058;  
        // if(i==3||i==4){
        //     double px = particle->GetPx();
        //     double py = particle->GetPy();
        //     bool f1 = (0.4<abs(py)&&abs(py)<0.8);
        //     bool f2 = (-0.27<px);
        //     bool f3 = (pow(px+0.6, 2)+pow(py, 2)<1.25);
        //     if(!(f1&&f2&&f3)){
        //         return StarGenEvent::kReject;
        //     }
        // }

        //"stable" particle check
        switch(abs(particle->GetId())){
        case 111:// pi0 stable to permit mother/daughter in star record
        case 211:// pi+/- stable
        case 221:// eta stable
        case 321:// K+/- stable
        case 310:// K short
        case 130:// K long
        case 3122:// Lambda 0 
        case 3112:// Sigma -
        case 3222:// Sigma +
        case 3212:// Sigma 0
        case 3312:// Xi -
        case 3322:// Xi 0
        case 3334:// Omega -
            //if there is a particle like that
            //event can be accepted
            return StarGenEvent::kAccept;
            break;
        default:
            break;
        }

        //checking only outgoing particles
        if(particle->GetStatus()!=1){
            continue;
        }
        //checking if the particle falls into (roughly) fiducial region for the detection
        //here |eta|<0.9+0.3, pT>0.2-0.05
        double eta = 0.9;
        double d_eta = 0.3;
        double pT = 0.2;
        double d_pT = 0.05;
        TLorentzVector particleFourmomentum = particle->momentum();
        if(fabs(particleFourmomentum.Eta())>=eta+d_eta||particleFourmomentum.Pt()<pT-d_pT){
            continue;
        }
        //decision
        if(hasAlreadyOneParticleInTPCFiducialRegion){
            return StarGenEvent::kAccept;
        } else{
            hasAlreadyOneParticleInTPCFiducialRegion = true;
        }
    }

    return StarGenEvent::kReject;
}