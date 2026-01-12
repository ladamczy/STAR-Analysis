#include "ParticlePhiFilter.h"
#include "StarGenerator/EVENT/StarGenParticle.h"
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <memory>
#include "DllImport.h"

ClassImp(ParticlePhiFilter);

ParticlePhiFilter::ParticlePhiFilter():StarFilterMaker("ParticlePhiFilter"){
    this->SetNameTitle("ParticlePhiFilter", "Filter responsible for accepting only events with phi particles");
};
ParticlePhiFilter::~ParticlePhiFilter(){};

Int_t ParticlePhiFilter::Filter(StarGenEvent* event){
    //loop through all particles
    //if there is no phi, event gets rejected
    //if there is even one - it gets accepted
    //also, there needs to be at least one particle with mStack>0
    bool hasStrangeParticle = false;
    bool hasParticleWithPositiveStack = false;
    StarGenParticle* particle = nullptr;
    for(size_t i = 0; i<event->GetNumberOfParticles(); i++){
        particle = (*event)[i];
        if(particle->GetStatus()!=1){
            continue;
        }
        //strange particle check
        if(abs(particle->GetId())==333){
            //if there is a particle like that
            //event can be accepted
            hasStrangeParticle = true;
        }
        //positive stack check
        if(particle->GetStack()>0){
            hasParticleWithPositiveStack = true;
        }

        //decision
        if(hasStrangeParticle&&hasParticleWithPositiveStack){
            return StarGenEvent::kAccept;
        }
    }

    return StarGenEvent::kReject;
}