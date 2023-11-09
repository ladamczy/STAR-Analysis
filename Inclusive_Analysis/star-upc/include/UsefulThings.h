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

using namespace std;

const int nTriggers = 17;
const int triggerID[] = { 570209, 570219, 570229, 570701, 570702, 570703, 570704, 570705, 
                  570709, 570711, 570712, 570719, 590701, 590703, 590705, 590708, 590709};
// 570702 RP_UPC // 570712 RP_UPC // 570703 RP_SDT // 570709 RP_ET // 570719 RP_ET // 570701 RP_CPT2 // 570711 RP_CPT2 // 570705 RP_CPT2noBBCL // 570704 RP_Zerobias // 590703 RP_SDT // 590709 RP_ET // 590701 RP_CPT2 // 590705 RP_CPT2noBBCL // 590708 RP_CPTnoBBCL // 570209 JPsi*HTTP // 570219 JPsi*HTTP // 570229 JPsi*HTTP
const int CEPtriggers[] = { 570701, 570705, 570711, 590701, 590705, 590708};

bool CheckTriggers(StUPCEvent* localupcEvt){
    bool CPTtrigger = false;
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

#endif