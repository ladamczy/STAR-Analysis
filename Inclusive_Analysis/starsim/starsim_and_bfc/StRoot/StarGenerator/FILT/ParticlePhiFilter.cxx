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

    //have been removed and replaced
    StarGenParticle* particle = nullptr;
    for(size_t i = 0; i<event->GetNumberOfParticles(); i++){
        particle = (*event)[i];
        //strange particle check
        if(abs(particle->GetId())!=333){
            continue;
        }
        if(particle->GetStatus()!=2){
            continue;
        }
        //checking if the children fall into (roughly) fiducial region for the detection
        //here |eta|<0.9+0.1, pT>0.2-0.05
        double eta = 0.9;
        double d_eta = 0.3;
        double pT = 0.2;
        double d_pT = 0.05;
        TLorentzVector firstDaughter = (*event)[particle->GetFirstDaughter()]->momentum();
        TLorentzVector lastDaughter = (*event)[particle->GetLastDaughter()]->momentum();
        if(fabs(firstDaughter.Eta())>=eta+d_eta||firstDaughter.Pt()<pT-d_pT){
            continue;
        }
        if(fabs(lastDaughter.Eta())>=(eta+d_eta)||lastDaughter.Pt()<(pT-d_pT)){
            continue;
        }
        //checking the children (all combinations)
        if((*event)[particle->GetFirstDaughter()]->GetId()==321&&(*event)[particle->GetLastDaughter()]->GetId()==-321){
            return StarGenEvent::kAccept;
        }
        if((*event)[particle->GetFirstDaughter()]->GetId()==-321&&(*event)[particle->GetLastDaughter()]->GetId()==321){
            return StarGenEvent::kAccept;
        }
    }

    return StarGenEvent::kReject;
}
