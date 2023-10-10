#include <iostream>
#include <string>

// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"

// ROOT, for histogramming.
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TF1.h"

using namespace Pythia8;

using namespace std;

int main(int argc, char const *argv[])
{
    istringstream nEventsStream(argv[1]);
    int nEvents;
    if (!(nEventsStream >> nEvents)){
        cout<<"Invalid first argument!"<<endl;
        return 0;
    }

    int mEnergy = 510;
    if(argc>1){
        istringstream mEnergyStream(argv[2]);
        
        if (!(mEnergyStream >> mEnergy)){
            cout<<"Invalid second argument!"<<endl;
            return 0;
        }
    }


    cout<<"Number of events: "<<nEvents<<endl;
    cout<<"Energy: "<<mEnergy<<endl;
    string filename = to_string(nEvents) + "_" + to_string(mEnergy) + "GeV_" + string(argv[0]).substr(string(argv[0]).find_last_of("/\\") + 1) + ".root";
    cout<<"Output file: "<<filename<<endl;

    Pythia pythia;

    if( mEnergy == 510){
        pythia.readString("Beams:eCM = 510"); //beam energy in GeV
    }else{
        cout<<"Input energy not equal to 510 GeV, if you wish to continue write \"yes\", otherwise write anything else"<<endl;
        string test;
        cin>>test;
        if(!strcmp(test.c_str(), ("yes"))){
            return 0;
        }
        pythia.readString("Beams:eCM = " + to_string(mEnergy));
    }

    pythia.readString("SoftQCD:centralDiffractive = on");
    pythia.readString("SigmaTotal:zeroAXB = off");
    //K0S exclusive
    pythia.readString("310:onMode=0");
    pythia.readString("310:OnIfMatch=211 -211");
    pythia.readString("-310:onMode=0");
    pythia.readString("-310:OnIfMatch=-211 211");

    pythia.init();

    //production loop
    for(int iEvent = 0; iEvent < nEvents; ++iEvent){
        if(!pythia.next()) continue;
    }

    // Statistics on event generation.
    pythia.stat();
                        
    return 0;
}
