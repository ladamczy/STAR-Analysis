#ifndef USEFUL_STUFF_H
#define USEFUL_STUFF_H

//normal stuff
#include <iostream> 
#include <string>
#include <fstream>

//ROOT & picoDST
#include "StUPCEvent.h"
#include "StUPCTrack.h"
#include "TChain.h"
#include "TFile.h"
#include "TParticle.h"
#include "TParticlePDG.h"
#include "StRPEvent.h"

using namespace std;

const int triggerID[] = { 570209, 570219, 570229, 570701, 570702, 570703, 570704, 570705, 
                  570709, 570711, 570712, 570719, 590701, 590703, 590705, 590708, 590709};
// 570702 RP_UPC // 570712 RP_UPC // 570703 RP_SDT // 570709 RP_ET // 570719 RP_ET // 570701 RP_CPT2 // 570711 RP_CPT2 // 570705 RP_CPT2noBBCL // 570704 RP_Zerobias // 590703 RP_SDT // 590709 RP_ET // 590701 RP_CPT2 // 590705 RP_CPT2noBBCL // 590708 RP_CPTnoBBCL // 570209 JPsi*HTTP // 570219 JPsi*HTTP // 570229 JPsi*HTTP
const int CEPtriggers[] = { 570701, 570705, 570711, 590701, 590705, 590708};
const int CEPCutDowntriggers[] = { 570701, 570705, 570711 };
const double beamMomentum = 254.867;

bool CheckTriggers(StUPCEvent* localupcEvt){
    bool CPTtrigger = false;
    const int nTriggers = 17;
    for(int var = 0; var < nTriggers; ++var){
        if(localupcEvt->isTrigger(triggerID[var])){
            //Checked if it is CPT trigger
            for (int i = 0; i < *(&CEPtriggers + 1) - CEPtriggers; ++i)
                if(triggerID[var] == CEPtriggers[i])
                    CPTtrigger=true;
        }
    }
    return CPTtrigger;
}

bool ConnectInput(int argc, char **argv, TChain *fileChain)
{
    int fileId = -1;
    string line;
    int lineId=0;
    //for testing if file exists
    TFile* infile;

    const string& input = argv[1];
    if(input.find(".list") != string::npos)
    {
        cout<<"Using list "<<input<<endl;
        ifstream instr(input.c_str());
        if (!instr.is_open())
        {
            cout<< "Couldn't open: "<<input.c_str()<<endl;
            return false;
        }


        while(getline(instr, line)) 
        {
            if(fileId==lineId || fileId== -1)
            {
                fileChain->AddFile(line.c_str());
                infile = TFile::Open(line.c_str(), "read");
                if(!infile)
                {
                    cout<< "Couldn't open: "<<line.c_str()<<endl;
                    return false;
                }
                infile->Close();
            }
            lineId++;
        }
        instr.close();
    }else if(input.find(".root") != string::npos){
        cout<<"Using root file "<<input<<endl;
        fileChain->AddFile(input.c_str());
        infile = TFile::Open(input.c_str(), "read");
        if(!infile){
            cout<< "Couldn't open: "<<input.c_str()<<endl;
            return false;
        }
        infile->Close();
    }
    
    cout<<"There are "<<fileChain->GetEntries()<<" entries"<<endl;
    return true;
}

double MetricCheck(TParticle* MCParticle, StUPCTrack* DETParticle){
    if(int(round(MCParticle->GetPDG()->Charge())) != DETParticle->getCharge()*3){
        return -1;
    }
    TLorentzVector tempMCvec;
    TVector3 tempDETvec;  
    MCParticle->Momentum(tempMCvec);
    DETParticle->getMomentum(tempDETvec);
    double phi_weight = 1.0;
    double eta_weight = 1.0;
    return sqrt(pow(phi_weight*(tempMCvec.Phi()-tempDETvec.Phi()),2)+pow(eta_weight*(tempMCvec.Eta()-tempDETvec.Eta()),2));
}

bool CutDownTriggers(StUPCEvent *localupcEvt){
    bool CPTtrigger = false;
    for(int var = 0; var<17; ++var){
        if(localupcEvt->isTrigger(triggerID[var])){
            //Checked if it is CPT trigger
            for(int i = 0; i<3; ++i){
                if(triggerID[var]==CEPCutDowntriggers[i]){
                    CPTtrigger = true;
                    return CPTtrigger;
                }
            }
        }
    }

    return CPTtrigger;
}

//cm->m->ns (0.3m=1ns), those are to remeber how it works
// deltaT1 = vector_Track_positive[0]->getTofPathLength()/100.0/0.299792458*sqrt(1+pow(particleMass[Kaon]/temp4Vector1.P(), 2));
// deltaT2 = vector_Track_negative[0]->getTofPathLength()/100.0/0.299792458*sqrt(1+pow(particleMass[Kaon]/temp4Vector1.P(), 2));

double t0(StUPCTrack* track, double mass){
    TVector3 mom;
    track->getMomentum(mom);
    //t1-(t1-t0), in ns
    return track->getTofTime()-track->getTofPathLength()/100.0/0.299792458*sqrt(1+pow(mass/mom.Mag(), 2));
}

double DeltaT0(StUPCTrack* track1, StUPCTrack* track2, double mass1, double mass2){
    return t0(track1, mass1)-t0(track2, mass2);
}

double M2TOF(StUPCTrack* track1, StUPCTrack* track2){
    double deltaT = (track1->getTofTime()-track2->getTofTime())*100.0*0.299792458;
    double L1 = track1->getTofPathLength();
    double L2 = track2->getTofPathLength();
    TVector3 mom1, mom2;
    track1->getMomentum(mom1);
    track2->getMomentum(mom2);
    double p1 = mom1.Mag();
    double p2 = mom2.Mag();
    double A = -2*pow(L1*L2/p1/p2, 2)+pow(L1/p1, 4)+pow(L2/p2, 4);
    double B = -2*pow(L1*L2, 2)*(pow(p1, -2)+pow(p2, -2))+2*pow(L1*L1/p1, 2)+2*pow(L2*L2/p2, 2)-2*pow(deltaT, 2)*(pow(L1/p1, 2)+pow(L2/p2, 2));
    double C = pow(deltaT, 4)-2*pow(deltaT, 2)*(L1*L1+L2*L2)+pow(L1*L1-L2*L2, 2);
    double m2TOF = (-B+sqrt(B*B-4*A*C))/2/A;
    return m2TOF;
}

#endif