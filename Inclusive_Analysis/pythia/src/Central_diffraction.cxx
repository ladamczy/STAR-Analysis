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

// My classes and stuff
#include "XiHistogramsClass.h"
#include "UsefulThingsPythia.h"

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
    if(argc>2){
        istringstream mEnergyStream(argv[2]);
        
        if (!(mEnergyStream >> mEnergy)){
            cout<<"Invalid second argument!"<<endl;
            return 0;
        }
    }

    string folder = "";
    if(argc>3){
        istringstream mFolderStream(argv[3]);
        
        if (!(mFolderStream >> folder)){
            cout<<"Invalid third argument!"<<endl;
            return 0;
        }
    }

    cout<<"Number of events: "<<nEvents<<endl;
    cout<<"Energy: "<<mEnergy<<endl;
    string filename = to_string(nEvents) + "_" + to_string(mEnergy) + "GeV_" + string(argv[0]).substr(string(argv[0]).find_last_of("/\\") + 1) + ".root";
    cout<<"Output file: "<<filename<<" in folder "<<folder<<endl;

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

    //is order important?
    pythia.readString("SoftQCD:centralDiffractive = on");
    pythia.readString("SigmaTotal:zeroAXB = off");
    //K0S exclusive
    pythia.readString("310:onMode=0");
    pythia.readString("310:OnIfMatch=211 -211");
    pythia.readString("-310:onMode=0");
    pythia.readString("-310:OnIfMatch=-211 211");

    pythia.init();

    //histograms & histogram class initialization
    vector<XiHistogramsClass> histvector;
    histvector.push_back(XiHistogramsClass(" before cuts"));
    histvector.push_back(XiHistogramsClass(" after TPC cut"));
    histvector.push_back(XiHistogramsClass(" after TPC & BBCL cuts"));
    histvector.push_back(XiHistogramsClass(" after TPC, BBCL & K^{0}_{S} cuts"));
    histvector.push_back(XiHistogramsClass(" after TPC, BBCL, K^{0}_{S} & extra particles cuts"));
    TH1D npart("npart", "Number of particles of not-K^{0}_{S} origin;particles;events", 30, 0, 30);

    //production loop
    for(int iEvent = 0; iEvent < nEvents; ++iEvent){
        if(!pythia.next()) continue;

        //fishing for diffracted protons
        int p1index=0, p2index=0;
        for (int part_index = 0; part_index < pythia.event.size(); part_index++){
            if (pythia.event[part_index].id() != 2212 || pythia.event[pythia.event[part_index].mother1()].id() != 2212){
                continue;
            }
            //choosing two with highest energies
            if(p1index==0){
                p1index = part_index;
            }else if(p2index==0){
                p2index = part_index;
                //switch so 1st proton has pz>0
                if(pythia.event[p1index].pz()<pythia.event[p2index].pz()){
                    p2index = p1index;
                    p1index = part_index;
                }
            }else if(pythia.event[p1index].pz()<pythia.event[part_index].pz()){
                p1index = part_index;
            }else if(pythia.event[p2index].pz()>pythia.event[part_index].pz()){
                p2index = part_index;
            }
        }

        //putting them into TLorentzVectors
        TLorentzVector p1, p2;
        p1.SetPxPyPzE(pythia.event[p1index].px(), pythia.event[p1index].py(), pythia.event[p1index].pz(), pythia.event[p1index].e());
        p2.SetPxPyPzE(pythia.event[p2index].px(), pythia.event[p2index].py(), pythia.event[p2index].pz(), pythia.event[p2index].e());

        histvector[0].AddProton(&p1, &p2);
        
        //TPC cut
        if(!hasAtLeastNParticlesInTPC(pythia.event, 2)){
            continue;
        }

        histvector[1].AddProton(&p1, &p2);

        //BBCL cut
        if(hasAtLeastNParticlesInBBCLarge(pythia.event, 1)){
            continue;
        }

        histvector[2].AddProton(&p1, &p2);

        //K0 decaying into 2pi0 cut
        bool doWeSeeK0 = false;
        for (int part_index = 0; part_index < pythia.event.size(); part_index++){
            if(pythia.event[part_index].id()==310 && isParticleInTPCAcceptance(pythia.event[pythia.event[part_index].daughter1()]) && isParticleInTPCAcceptance(pythia.event[pythia.event[part_index].daughter2()])){
                doWeSeeK0 = true;
                break;
            }
        }
        if(!doWeSeeK0){
            continue;
        }

        histvector[3].AddProton(&p1, &p2);

        //at least 2 non-K0 related particles
        int nonK0relatedParticles = 0;
        for (int part_index = 0; part_index < pythia.event.size(); part_index++){
            if(pythia.event[pythia.event[part_index].mother1()].id()!=310 && isParticleInTPCAcceptance(pythia.event[part_index])){
                nonK0relatedParticles++;
            }
        }

        npart.Fill(nonK0relatedParticles);

        if(nonK0relatedParticles<2){
            continue;
        }

        histvector[4].AddProton(&p1, &p2);
    }

    // Statistics on event generation.
    pythia.stat();

    //writing histograms to file
    if(folder.length()!=0){
        filename = folder + "/" + filename;
    }
    TFile* output1 = new TFile(filename.c_str(), "RECREATE");
    TFile* output2 = new TFile((filename.replace(filename.rfind("."), 5, "_ratios.root")).c_str(), "RECREATE");
    TH2D hist1;
    TH1D hist2;
    TH1D hist3;
    TH1D temphist2;
    TH1D temphist3;
    TH1D* temphist2towrite;
    TH1D* temphist3towrite;

    for (long unsigned int i = 0; i < histvector.size(); i++){
        histvector[i].RetrieveHistograms(hist1, hist2, hist3);
        output1->cd();
        hist1.Write();
        hist2.Write();
        hist3.Write();
        if(i==3){
            npart.Write();
        }
        if(i!=0){
            output2->cd();
            temphist2towrite = (TH1D*)hist2.Clone();
            temphist2towrite->Divide(&temphist2);
            temphist2towrite->SetTitle((string(temphist2towrite->GetTitle())+" (ratio);"+string(temphist2towrite->GetXaxis()->GetTitle())+";ratio").c_str());
            temphist2towrite->Write();
            temphist3towrite = (TH1D*)hist3.Clone();
            temphist3towrite->Divide(&temphist3);
            temphist3towrite->SetTitle((string(temphist3towrite->GetTitle())+" (ratio);"+string(temphist3towrite->GetXaxis()->GetTitle())+";ratio").c_str());
            temphist3towrite->Write();
        }
        temphist2 = *((TH1D*)hist2.Clone());
        temphist3 = *((TH1D*)hist3.Clone());
    }
    
    output1->Close();
    output2->Close();
                        
    return 0;
}
