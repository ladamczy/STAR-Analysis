#include "ParticleK0SFilter.h"
#include "StarGenerator/EVENT/StarGenParticle.h"
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <memory>
#include "DllImport.h"

ClassImp(ParticleK0SFilter);

ParticleK0SFilter::ParticleK0SFilter():StarFilterMaker("ParticleK0SFilter"){
    this->SetNameTitle("ParticleK0SFilter", "Filter responsible for accepting only events with K0S particles");
};
ParticleK0SFilter::~ParticleK0SFilter(){};

Int_t ParticleK0SFilter::Filter(StarGenEvent* event){
    //loop through all particles
    //if there is no K0S, event gets rejected
    //if there is even one - it gets accepted
    //also, there needs to be at least one particle with mStack>0
    StarGenParticle* particle = nullptr;
    for(size_t i = 0; i<event->GetNumberOfParticles(); i++){
        particle = (*event)[i];
        if(particle->GetStatus()!=1){
            continue;
        }
        //positive stack check
        if(particle->GetStack()<0){
            continue;
        }
        //K0S particle check
        if(abs(particle->GetId())!=310){
            continue;
        }

        // is in detectable region under the curve
        TLorentzVector mom = particle->momentum();
        if(fabs(mom.Eta())<=1.0){
            return StarGenEvent::kAccept;
        } else if(0.14/(fabs(mom.Eta())-1)+0.07>mom.Pt()){
            return StarGenEvent::kAccept;
        }

        return StarGenEvent::kAccept;
    }

    return StarGenEvent::kReject;
}
