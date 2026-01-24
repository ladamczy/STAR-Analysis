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
    //due to short decay time, the decay is handled by pythia, and therefore it does not come out to the Geant
    //so checks like:

    //particle->GetStatus()!=1
    //particle->GetStack()>0

    //have been removed
    StarGenParticle* particle = nullptr;
    for(size_t i = 0; i<event->GetNumberOfParticles(); i++){
        particle = (*event)[i];
        //strange particle check
        if(abs(particle->GetId())==333){
            //if there is a particle like that
            //event can be accepted
            return StarGenEvent::kAccept;
        }
    }

    return StarGenEvent::kReject;
}
