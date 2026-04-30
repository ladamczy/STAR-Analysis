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
        if(fabs(particle->Eta())>=eta+d_eta||particle->Pt()<pT-d_pT){
            continue;
        }
        //decision
        if(hasAlreadyOneParticleInTPCFiducialRegion){
            return StarGenEvent::kAccept;
        }else{
            hasAlreadyOneParticleInTPCFiducialRegion = true;
        }
    }

    return StarGenEvent::kReject;
}