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
    TH2D hist1tab[] = {TH2D(), TH2D(), TH2D(), TH2D(), TH2D()};
    TH1D hist2tab[] = {TH1D(), TH1D(), TH1D(), TH1D(), TH1D()};
    TH1D hist3tab[] = {TH1D(), TH1D(), TH1D(), TH1D(), TH1D()};
    XiHistogramsClass histgroup1(" before cuts", "1", &hist1tab[0], &hist2tab[0], &hist3tab[0]);
    XiHistogramsClass histgroup2(" after TPC cut", "2", &hist1tab[1], &hist2tab[1], &hist3tab[1]);
    XiHistogramsClass histgroup3(" after TPC & BBCL cuts", "3", &hist1tab[2], &hist2tab[2], &hist3tab[2]);
    XiHistogramsClass histgroup4(" after TPC, BBCL & K^{0}_{S} cuts", "4", &hist1tab[3], &hist2tab[3], &hist3tab[3]);
    XiHistogramsClass histgroup5(" after TPC, BBCL, K^{0}_{S} & extra particles cuts", "5", &hist1tab[4], &hist2tab[4], &hist3tab[4]);
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

        histgroup1.AddProton(&p1, &p2);
        
        //TPC cut
        if(!hasAtLeastNParticlesInTPC(pythia.event, 2)){
            continue;
        }

        histgroup2.AddProton(&p1, &p2);

        //BBCL cut
        if(hasAtLeastNParticlesInBBCLarge(pythia.event, 1)){
            continue;
        }

        histgroup3.AddProton(&p1, &p2);

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

        //at least 2 non-K0 related particles
        int nonK0relatedParticles = 0;
        for (int part_index = 0; part_index < pythia.event.size(); part_index++){
            if(pythia.event[part_index].mother1()!=310 && isParticleInTPCAcceptance(pythia.event[part_index])){
                nonK0relatedParticles++;
            }
        }

        npart.Fill(nonK0relatedParticles);
        histgroup4.AddProton(&p1, &p2);

        if(nonK0relatedParticles<2){
            continue;
        }

        histgroup5.AddProton(&p1, &p2);
    }

    // Statistics on event generation.
    pythia.stat();

    //writing histograms to file
    if(folder.length()!=0){
        filename = folder + "/" + filename;
    }
    TFile* output = new TFile(filename.c_str(), "RECREATE");
    output->cd();
    TH2D* hist1 = nullptr;
    TH1D* hist2 = nullptr;
    TH1D* hist3 = nullptr;

    histgroup1.RetrieveHistograms(hist1, hist2, hist3);
    hist1->Write();
    hist2->Write();
    hist3->Write();
    histgroup2.RetrieveHistograms(hist1, hist2, hist3);
    hist1->Write();
    hist2->Write();
    hist3->Write();
    histgroup3.RetrieveHistograms(hist1, hist2, hist3);
    hist1->Write();
    hist2->Write();
    hist3->Write();
    histgroup4.RetrieveHistograms(hist1, hist2, hist3);
    hist1->Write();
    hist2->Write();
    hist3->Write();
    histgroup5.RetrieveHistograms(hist1, hist2, hist3);
    hist1->Write();
    hist2->Write();
    hist3->Write();

    output->Close();
                        
    return 0;
}
