#include "Pythia8/Pythia.h"

bool isParticleInTPCAcceptance(Pythia8::Particle particle){
    return abs(particle.eta())<0.9 && particle.pT()>0.2 && particle.isCharged() && particle.isFinal();
}

bool hasAtLeastNParticlesInTPC(Pythia8::Event event, int nParticles){
    int n = 0;
    for (int i = 0; i < event.size(); i++){
        if(isParticleInTPCAcceptance(event[i])){
            n++;
        }
    }

    return n>=nParticles;
}

bool hasAtLeastNParticlesInBBCLarge(Pythia8::Event event, int nParticles){
    int n = 0;
    for (int i = 0; i < event.size(); i++){
        if(abs(event[i].eta())>2.1 && abs(event[i].eta())<3.3 && event[i].pT()>0.5 && event[i].isFinal()){
            n++;
        }
    }

    return n>=nParticles;
}

double xi(TLorentzVector *fourvec, double s = 510){
    return (s/2-fourvec->E())/s*2;
}