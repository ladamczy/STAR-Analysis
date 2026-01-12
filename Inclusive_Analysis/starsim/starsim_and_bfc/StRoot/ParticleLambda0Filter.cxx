#include "ParticleLambda0Filter.h"
#include "StarGenerator/EVENT/StarGenParticle.h"
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <memory>
#include "DllImport.h"

ClassImp(ParticleLambda0Filter);

ParticleLambda0Filter::ParticleLambda0Filter():StarFilterMaker("ParticleLambda0Filter"){
    this->SetNameTitle("ParticleLambda0Filter", "Filter responsible for accepting only events with Lambda0 particles");
};
ParticleLambda0Filter::~ParticleLambda0Filter(){};

Int_t ParticleLambda0Filter::Filter(StarGenEvent* event){
    //loop through all particles
    //if there is no Lambda0, event gets rejected
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
        if(abs(particle->GetId())==3122){
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