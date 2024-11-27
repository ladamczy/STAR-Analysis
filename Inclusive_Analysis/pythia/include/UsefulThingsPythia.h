#include "Pythia8/Pythia.h"
#include <TLorentzVector.h>
#include <TH3F.h>
#include <TFile.h>

#include <string>

double beamMomentum = 254.867;

enum PARTICLES{ Pion = 0, Kaon = 1, Proton = 2, nParticles };
enum SIGN{ Minus, Plus, nSigns };

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

double xi(TLorentzVector fourvec, double p = 254.867){
    return (p-fourvec.P())/p;
}

bool readingAcceptanceFromFile(TH3F* Acceptance){
    TFile* file = TFile::Open("~/STAR-Analysis/share/etaPhiEfficiency_16_01_19_delta015_twoRuns.root");
    for(size_t i = 0; i<nSigns*nParticles; i++){
        Acceptance[i] = *dynamic_cast<TH3F*>(file->Get(("hTPCEffiCD"+std::to_string(i)+"121").c_str()));
    }
    file->Close();
    return true;
}

bool isParticleDetected(Pythia8::Particle* particle, double zVx){
    //filling the  histogram once at the beginning
    //and other static things
    static TH3F Acceptance[nSigns*nParticles];
    static bool isInitialised = readingAcceptanceFromFile(Acceptance);
    static Pythia8::Rndm generator(0);
    //getting the histogram
    int histID = -1;
    switch(particle->id()){
    case 211:
        histID = Pion+nParticles*Plus;
        break;
    case -211:
        histID = Pion+nParticles*Minus;
        break;
    case 321:
        histID = Kaon+nParticles*Plus;
        break;
    case -321:
        histID = Kaon+nParticles*Minus;
        break;
    case 2212:
        histID = Proton+nParticles*Plus;
        break;
    case -2212:
        histID = Proton+nParticles*Minus;
        break;
    default:
        return generator.flat()<0.6;
    }
    TH3F* temp = &Acceptance[histID];
    //reading from the histogram (with safeguards at the edges)
    double pT = particle->pT();
    double eta = particle->eta();
    if(zVx<temp->GetXaxis()->GetBinCenter(1) or zVx>temp->GetXaxis()->GetBinCenter(temp->GetNbinsX())){
        zVx += 1e-3*(1-2*(zVx>0));
    }
    if(pT<temp->GetYaxis()->GetBinCenter(1) or pT>temp->GetYaxis()->GetBinCenter(temp->GetNbinsY())){
        pT += 1e-3*(1-2*(pT>0));
    }
    if(eta<temp->GetZaxis()->GetBinCenter(1) or eta>temp->GetZaxis()->GetBinCenter(temp->GetNbinsZ())){
        eta += 1e-3*(1-2*(eta>0));
    }
    return temp->Interpolate(zVx, pT, eta);
}