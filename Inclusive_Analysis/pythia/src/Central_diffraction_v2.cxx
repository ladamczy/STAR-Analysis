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
#include "TEfficiency.h"

// My classes and stuff
#include "UsefulThingsPythia.h"

using namespace Pythia8;

using namespace std;

int main(int argc, char const *argv[]){
    istringstream nEventsStream(argv[1]);
    int nEvents;
    if(!(nEventsStream>>nEvents)){
        cout<<"Invalid first argument!"<<endl;
        return 0;
    }

    int mEnergy = 510;
    if(argc>2){
        istringstream mEnergyStream(argv[2]);

        if(!(mEnergyStream>>mEnergy)){
            cout<<"Invalid second argument!"<<endl;
            return 0;
        }
    }

    string folder = "";
    if(argc>3){
        istringstream mFolderStream(argv[3]);

        if(!(mFolderStream>>folder)){
            cout<<"Invalid third argument!"<<endl;
            return 0;
        }
    }

    cout<<"Number of events: "<<nEvents<<endl;
    cout<<"Energy: "<<mEnergy<<endl;
    string filename = to_string(nEvents)+"_"+to_string(mEnergy)+"GeV_"+string(argv[0]).substr(string(argv[0]).find_last_of("/\\")+1)+".root";
    cout<<"Output file: "<<filename<<" in folder "<<folder<<endl;

    Pythia pythia;

    if(mEnergy==510){
        pythia.readString("Beams:eCM = 510"); //beam energy in GeV
    } else{
        cout<<"Input energy not equal to 510 GeV, if you wish to continue write \"yes\", otherwise write anything else"<<endl;
        string test;
        cin>>test;
        if(!strcmp(test.c_str(), ("yes"))){
            return 0;
        }
        pythia.readString("Beams:eCM = "+to_string(mEnergy));
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
    // xi_e_w = new TH2D(string("xi_e_w_"+identifier).c_str(), string("#xi_{E} vs #xi_{W}"+namepart+";#xi_{E};#xi_{W}").c_str(), 50, 0, 1, 50, 0, 1);
    // xi_multiplied = new TH1D(string("xi_multiplied_"+identifier).c_str(), string("#xi_{E}*#xi_{W}"+namepart+";#xi_{E}*#xi_{W};entries").c_str(), 60, 0, 0.3);
    // xi_log = new TH1D(string("xi_log_"+identifier).c_str(), string("ln(#xi_{E}/#xi_{W})"+namepart+";ln(#xi_{E}/#xi_{W});entries").c_str(), 100, -10, 10);
    // before
    TH1D LogCutEventsBefore("LogCutEventsBefore", "log(#xi_{E}*#xi_{W}) with |ln(#xi_{E}/#xi_{W})|<4 condition;log(#xi_{E}*#xi_{W});events", 60, -6, 0);
    TH1D MultiCutEventsBefore("MultiCutEventsBefore", "ln(#xi_{E}/#xi_{W}) with 0.007<#xi_{E,W}<0.08 conditions;ln(#xi_{E}/#xi_{W});events", 100, -10, 10);
    TH2D Log2DCutEventsBefore("Log2DCutEventsBefore", "log#xi_{E} vs log#xi_{W} with |ln(#xi_{E}/#xi_{W})|<4 condition;log#xi_{E};log#xi_{W}", 50, -4, 1, 50, -4, 1);
    //TPC
    TH1D LogTPCCutEventsAfter("LogTPCCutEventsAfter", "Effect of TPC cut on log(#xi_{E}*#xi_{W}) with |ln(#xi_{E}/#xi_{W})|<4 condition;log(#xi_{E}*#xi_{W});events", 60, -6, 0);
    TH1D MultiTPCCutEventsAfter("MultiTPCCutEventsAfter", "Effect of TPC cut on ln(#xi_{E}/#xi_{W}) with 0.007<#xi_{E,W}<0.08 conditions;ln(#xi_{E}/#xi_{W});events", 100, -10, 10);
    TH2D Log2DTPCCutEventsAfter("Log2DTPCCutEventsAfter", "Effect of TPC cut on log#xi_{E} vs log#xi_{W} with |ln(#xi_{E}/#xi_{W})|<4 condition;log#xi_{E};log#xi_{W}", 50, -4, 1, 50, -4, 1);
    //BBCL
    TH1D LogBBCLCutEventsAfter("LogBBCLCutEventsAfter", "Effect of BBCL cut on log(#xi_{E}*#xi_{W}) with |ln(#xi_{E}/#xi_{W})|<4 condition;log(#xi_{E}*#xi_{W});events", 60, -6, 0);
    TH1D MultiBBCLCutEventsAfter("MultiBBCLCutEventsAfter", "Effect of BBCL cut on ln(#xi_{E}/#xi_{W}) with 0.007<#xi_{E,W}<0.08 conditions;ln(#xi_{E}/#xi_{W});events", 100, -10, 10);
    TH2D Log2DBBCLCutEventsAfter("Log2DBBCLCutEventsAfter", "Effect of BBCL cut on log#xi_{E} vs log#xi_{W} with |ln(#xi_{E}/#xi_{W})|<4 condition;log#xi_{E};log#xi_{W}", 50, -4, 1, 50, -4, 1);
    //TPC+BBCL
    TH1D LogTPCBBCLCutEventsAfter("LogTPCBBCLCutEventsAfter", "Effect of TPC and BBCL cuts on log(#xi_{E}*#xi_{W}) with |ln(#xi_{E}/#xi_{W})|<4 condition;log(#xi_{E}*#xi_{W});events", 60, -6, 0);
    TH1D MultiTPCBBCLCutEventsAfter("MultiTPCBBCLCutEventsAfter", "Effect of TPC and BBCL cuts on ln(#xi_{E}/#xi_{W}) with 0.007<#xi_{E,W}<0.08 conditions;ln(#xi_{E}/#xi_{W});events", 100, -10, 10);
    TH2D Log2DTPCBBCLCutEventsAfter("Log2DTPCBBCLCutEventsAfter", "Effect of TPC and BBCL cut on log#xi_{E} vs log#xi_{W} with |ln(#xi_{E}/#xi_{W})|<4 condition;log#xi_{E};log#xi_{W}", 50, -4, 1, 50, -4, 1);
    //the rest
    TH2D nParticles("nParticles", "Number of tracks visible in TPC vs number of detected K_{S}^{0}, with both #xi conditions;tracks;# of K_{S}^{0}", 20, 0, 20, 4, 0, 4);
    TH2D nParticlesWithMoreCuts("nParticlesWithMoreCuts", "Number of tracks visible in TPC vs number of detected K_{S}^{0}, with TPC & BBCL cuts, and both #xi conditions;tracks;# of K_{S}^{0}", 20, 0, 20, 4, 0, 4);
    //more 
    TH2D pTEtaK0("pTEtaK0", "p_{T} vs #eta of K0 prticles, with both #xi conditions;p_{T} [GeV/c];#eta", 50, 0, 5, 20, -5, 5);
    TEfficiency K0detectionProbability("K0detectionProbability", "Probability of K_{S}^{0} detection vs log(#xi_{E}*#xi_{W});log(#xi_{E}*#xi_{W});probability", 60, -6, 0);
    //production loop
    for(int iEvent = 0; iEvent<nEvents; ++iEvent){
        if(!pythia.next()) continue;

        //fishing for diffracted protons
        int p1index = 0, p2index = 0;
        for(int part_index = 0; part_index<pythia.event.size(); part_index++){
            if(pythia.event[part_index].id()!=2212||pythia.event[pythia.event[part_index].mother1()].id()!=2212){
                continue;
            }
            //choosing two with highest energies
            if(p1index==0){
                p1index = part_index;
            } else if(p2index==0){
                p2index = part_index;
                //switch so 1st proton has pz>0
                if(pythia.event[p1index].pz()<pythia.event[p2index].pz()){
                    p2index = p1index;
                    p1index = part_index;
                }
            } else if(pythia.event[p1index].pz()<pythia.event[part_index].pz()){
                p1index = part_index;
            } else if(pythia.event[p2index].pz()>pythia.event[part_index].pz()){
                p2index = part_index;
            }
        }

        //putting them into TLorentzVectors
        TLorentzVector p1, p2;
        p1.SetPxPyPzE(pythia.event[p1index].px(), pythia.event[p1index].py(), pythia.event[p1index].pz(), pythia.event[p1index].e());
        p2.SetPxPyPzE(pythia.event[p2index].px(), pythia.event[p2index].py(), pythia.event[p2index].pz(), pythia.event[p2index].e());

        //CUTS
        //TPC cut
        bool TPCCut = hasAtLeastNParticlesInTPC(pythia.event, 2);
        //BBCL cut
        bool BBCLCut = !hasAtLeastNParticlesInBBCLarge(pythia.event, 1);
        //xi1>0.007 && xi2>0.007
        bool xi1Andxi2 = xi(&p1)>0.007&&xi(&p2)>0.007&&xi(&p1)<0.08&&xi(&p2)<0.08;
        //log cut
        bool LogCut = abs(log(xi(&p1)/xi(&p2)))<4;
        //2xK0 decaying into pi+- cut
        int nK0 = 0;
        int tempDaughter1, tempDaughter2;
        std::vector<int> K0indices;
        for(int part_index = 0; part_index<pythia.event.size(); part_index++){
            if(pythia.event[part_index].id()==310){ K0indices.push_back(part_index); }
        }
        for(long unsigned int K0candidate = 0; K0candidate<K0indices.size(); K0candidate++){
            tempDaughter1 = pythia.event[K0indices[K0candidate]].daughter1();
            tempDaughter2 = pythia.event[K0indices[K0candidate]].daughter2();
            if(isParticleInTPCAcceptance(pythia.event[tempDaughter1])&&isParticleInTPCAcceptance(pythia.event[tempDaughter2])&&(tempDaughter1+1==tempDaughter2)){
                nK0++;
                if(/*possible cuts*/true){
                    pTEtaK0.Fill(pythia.event[K0indices[K0candidate]].pT(), pythia.event[K0indices[K0candidate]].eta());
                }
            }
        }

        //HISTOGRAMS
        //point 3
        //before
        if(LogCut){ LogCutEventsBefore.Fill(log10(xi(&p1)*xi(&p2))); }
        if(xi1Andxi2){ MultiCutEventsBefore.Fill(log(xi(&p1)/xi(&p2))); }
        if(LogCut){ Log2DCutEventsBefore.Fill(log10(xi(&p1)), log10(xi(&p2))); }
        //TPC
        if(TPCCut&&LogCut){ LogTPCCutEventsAfter.Fill(log10(xi(&p1)*xi(&p2))); }
        if(TPCCut&&xi1Andxi2){ MultiTPCCutEventsAfter.Fill(log(xi(&p1)/xi(&p2))); }
        if(TPCCut&&LogCut){ Log2DTPCCutEventsAfter.Fill(log10(xi(&p1)), log10(xi(&p2))); }
        //BBCL
        if(BBCLCut&&LogCut){ LogBBCLCutEventsAfter.Fill(log10(xi(&p1)*xi(&p2))); }
        if(BBCLCut&&xi1Andxi2){ MultiBBCLCutEventsAfter.Fill(log(xi(&p1)/xi(&p2))); }
        if(BBCLCut&&LogCut){ Log2DBBCLCutEventsAfter.Fill(log10(xi(&p1)), log10(xi(&p2))); }
        //TOC & BBCL
        if(TPCCut&&BBCLCut&&LogCut){ LogTPCBBCLCutEventsAfter.Fill(log10(xi(&p1)*xi(&p2))); }
        if(TPCCut&&BBCLCut&&xi1Andxi2){ MultiTPCBBCLCutEventsAfter.Fill(log(xi(&p1)/xi(&p2))); }
        if(TPCCut&&BBCLCut&&LogCut){ Log2DTPCBBCLCutEventsAfter.Fill(log10(xi(&p1)), log10(xi(&p2))); }
        //point6
        int nPart = 0;
        if(xi1Andxi2&&LogCut){
            for(int i = 0; i<pythia.event.size(); i++){
                if(isParticleInTPCAcceptance(pythia.event[i])){
                    nPart++;
                }
            }
            nParticles.Fill(nPart, nK0);
            if(TPCCut&&BBCLCut){
                nParticlesWithMoreCuts.Fill(nPart, nK0);
            }
        }
        //even more
        K0detectionProbability.Fill(nK0>0&&xi1Andxi2&&LogCut&&TPCCut&&BBCLCut, log10(xi(&p1)*xi(&p2)));
    }

    // Statistics on event generation.
    pythia.stat();

    //HISTOGRAMS post-processing
    //using TEfficiency objects
    //point 3
    TEfficiency LogTPCCut = TEfficiency(LogTPCCutEventsAfter, LogCutEventsBefore);
    LogTPCCut.SetNameTitle("LogTPCCutEfficiency", "Effect of TPC cut on log(#xi_{E}*#xi_{W}) with |ln(#xi_{E}/#xi_{W})|<4 condition;log(#xi_{E}*#xi_{W});efficiency");
    TEfficiency MultiTPCCut = TEfficiency(MultiTPCCutEventsAfter, MultiCutEventsBefore);
    MultiTPCCut.SetNameTitle("MultiTPCCutEfficiency", "Effect of TPC cut on ln(#xi_{E}/#xi_{W}) with 0.007<#xi_{E,W}<0.08 conditions;ln(#xi_{E}/#xi_{W});efficiency");
    TEfficiency Log2DTPCCut = TEfficiency(Log2DTPCCutEventsAfter, Log2DCutEventsBefore);
    Log2DTPCCut.SetNameTitle("Log2DTPCCutEfficiency", "Effect of TPC cut on log#xi_{E} vs log#xi_{W} with |ln(#xi_{E}/#xi_{W})|<4 condition;log#xi_{E};log#xi_{W}");
    //point 4
    //division
    TEfficiency LogBBCLCut = TEfficiency(LogBBCLCutEventsAfter, LogCutEventsBefore);
    LogBBCLCut.SetNameTitle("LogBBCLCutEfficiency", "Effect of BBCL cut on log(#xi_{E}*#xi_{W}) with |ln(#xi_{E}/#xi_{W})|<4 condition;log(#xi_{E}*#xi_{W});efficiency");
    TEfficiency MultiBBCLCut = TEfficiency(MultiBBCLCutEventsAfter, MultiCutEventsBefore);
    MultiBBCLCut.SetNameTitle("MultiBBCLCutEfficiency", "Effect of BBCL cut on ln(#xi_{E}/#xi_{W}) with 0.007<#xi_{E,W}<0.08 conditions;ln(#xi_{E}/#xi_{W});efficiency");
    TEfficiency Log2DBBCLCut = TEfficiency(Log2DBBCLCutEventsAfter, Log2DCutEventsBefore);
    Log2DBBCLCut.SetNameTitle("Log2DBBCLCutEfficiency", "Effect of BBCL cut on log#xi_{E} vs log#xi_{W} with |ln(#xi_{E}/#xi_{W})|<4 condition;log#xi_{E};log#xi_{W}");
    //point 5
    //division
    TEfficiency LogTPCBBCLCut = TEfficiency(LogTPCBBCLCutEventsAfter, LogCutEventsBefore);
    LogTPCBBCLCut.SetNameTitle("LogTPCBBCLCutEfficiency", "Effect of TPC and BBCL cuts on log(#xi_{E}*#xi_{W}) with |ln(#xi_{E}/#xi_{W})|<4 condition;log(#xi_{E}*#xi_{W});efficiency");
    TEfficiency MultiTPCBBCLCut = TEfficiency(MultiTPCBBCLCutEventsAfter, MultiCutEventsBefore);
    MultiTPCBBCLCut.SetNameTitle("MultiTPCBBCLCutEfficiency", "Effect of TPC and BBCL cuts on ln(#xi_{E}/#xi_{W}) with 0.007<#xi_{E,W}<0.08 conditions;ln(#xi_{E}/#xi_{W});efficiency");
    TEfficiency Log2DTPCBBCLCut = TEfficiency(Log2DTPCBBCLCutEventsAfter, Log2DCutEventsBefore);
    Log2DTPCBBCLCut.SetNameTitle("Log2DTPCBBCLCutEfficiency", "Effect of TPC and BBCL cuts on log#xi_{E} vs log#xi_{W} with |ln(#xi_{E}/#xi_{W})|<4 condition;log#xi_{E};log#xi_{W}");

    //writing histograms to file
    if(folder.length()!=0){
        filename = folder+"/"+filename;
    }
    TFile *output1 = new TFile(filename.c_str(), "RECREATE");
    output1->cd();

    LogCutEventsBefore.Write();
    LogTPCCutEventsAfter.Write();
    LogTPCCut.Write();
    LogBBCLCutEventsAfter.Write();
    LogBBCLCut.Write();
    LogTPCBBCLCutEventsAfter.Write();
    LogTPCBBCLCut.Write();

    MultiCutEventsBefore.Write();
    MultiTPCCutEventsAfter.Write();
    MultiTPCCut.Write();
    MultiBBCLCutEventsAfter.Write();
    MultiBBCLCut.Write();
    MultiTPCBBCLCutEventsAfter.Write();
    MultiTPCBBCLCut.Write();

    Log2DCutEventsBefore.Write();
    Log2DTPCCutEventsAfter.Write();
    Log2DTPCCut.Write();
    Log2DBBCLCutEventsAfter.Write();
    Log2DBBCLCut.Write();
    Log2DTPCBBCLCutEventsAfter.Write();
    Log2DTPCBBCLCut.Write();

    //the other things
    nParticles.Write();
    nParticlesWithMoreCuts.Write();

    //the other other things
    pTEtaK0.Write();
    K0detectionProbability.Write();

    output1->Close();

    return 0;
}
