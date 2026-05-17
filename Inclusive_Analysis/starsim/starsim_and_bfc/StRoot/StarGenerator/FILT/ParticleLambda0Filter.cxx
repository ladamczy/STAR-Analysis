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
        if(abs(particle->GetId())!=3122){
            continue;
        }

        //TODO: measure and fill it (for hogher efficiency of simulation, not for a higher physical goal)

        // //is in detectable region under the curve
        // TLorentzVector mom = particle->momentum();
        // printf("Eta:\t%lf\n", mom.Eta());
        // printf("pT:\t%lf\n", mom.Pt());
        // printf("Border:\t%lf", 0.14/(fabs(mom.Eta())-1)+0.07);
        // if(fabs(mom.Eta())<=1.0){
        //     return StarGenEvent::kAccept;
        // } else if(0.14/(fabs(mom.Eta())-1)+0.07>mom.Pt()){
        //     return StarGenEvent::kAccept;
        // }

        return StarGenEvent::kAccept;
    }

    return StarGenEvent::kReject;
}