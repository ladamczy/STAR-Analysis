#include "Pythia8/Pythia.h"
#include <TLorentzVector.h>
#include <TH3F.h>
#include <TFile.h>

#include <string>

double beamMomentum = 254.867;

enum PARTICLES{ Pion = 0, Kaon = 1, Proton = 2, nParticles };
enum SIGN{ Minus, Plus, nSigns };
const double particleMass[nParticles] = { 0.13957, 0.493677, 0.93827 }; // pion, kaon, proton in GeV /c^2 

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

bool isParticleDetected(Pythia8::Particle* particle){
    //filling the  histogram once at the beginning
    //and other static things
    static TH3F Acceptance[nSigns*nParticles];
    //solving initialization in a prettier way
    static bool isInitialised = false;
    if(!isInitialised){
        isInitialised = readingAcceptanceFromFile(Acceptance);
    }
    //random number generator
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
        if(particle->charge()>0){
            histID = Pion+nParticles*Plus;
        } else if(particle->charge()<0){
            histID = Pion+nParticles*Minus;
        }
    }
    TH3F* temp = &Acceptance[histID];
    //reading from the histogram (with safeguards at the edges)
    double zVx = particle->zProd()/10;  //mm->cm
    double pT = particle->pT();
    double eta = particle->eta();
    if(zVx<temp->GetXaxis()->GetBinCenter(1) or zVx>temp->GetXaxis()->GetBinCenter(temp->GetNbinsX())){
        zVx += 1e-3*(1-2*(zVx>0));
    }
    if(pT<temp->GetYaxis()->GetBinCenter(1) or pT>temp->GetYaxis()->GetBinCenter(temp->GetNbinsY())){
        pT += 1e-3*(1-2*(pT>1.0));
    }
    if(eta<temp->GetZaxis()->GetBinCenter(1) or eta>temp->GetZaxis()->GetBinCenter(temp->GetNbinsZ())){
        eta += 1e-3*(1-2*(eta>0));
    }
    try{
        return temp->Interpolate(zVx, pT, eta)>generator.flat();
    }
    catch(const std::exception& e){
        std::cerr<<e.what()<<'\n';
        printf("%lf, %lf, %lf\n", zVx, pT, eta);
    }
    return temp->Interpolate(zVx, pT, eta)>generator.flat();
}