#include "StrangeHadronFilter.h"
#include "StarGenerator/EVENT/StarGenParticle.h"
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <memory>
#include "DllImport.h"

ClassImp(StrangeHadronFilter);

StrangeHadronFilter::StrangeHadronFilter():StarFilterMaker("strangehadronfilter"){
    this->SetNameTitle("strangehadronfilter", "Filter responsible for accepting only events with K0S/Lambda0/K*/phi particles");
};
StrangeHadronFilter::~StrangeHadronFilter(){};

Int_t StrangeHadronFilter::Filter(StarGenEvent* event){
    //loop through all particles
    //if there is no K0S, Lambda0, phi or K*
    //event gets rejected
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
        switch(abs(particle->GetId())){
        case 310:       //K0S
        case 3122:      //Lambda0
        case 313:       //K*
        case 333:       //phi
            //if there is a particle like that
            //event can be accepted
            hasStrangeParticle = true;
            break;
        default:
            break;
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