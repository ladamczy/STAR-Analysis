// main91.cc is a part of the PYTHIA event generator.
// Copyright (C) 2022 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Keywords: analysis; root;

// This is a simple test program.
// It studies the charged multiplicity distribution at the LHC.
// Modified by Rene Brun, Axel Naumann and Bernhard Meirose
// to use ROOT for histogramming.

// Stdlib header file for input and output.
#include <iostream>

// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"

// ROOT, for histogramming.
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TTree.h"
//#include "TChain.h"
#include "TLorentzVector.h"
#include "TF1.h"
//#include "TRandom3.h"
//#include "TString.h"
//#include "TMath.h"


// ROOT, for saving file.
//#include "TFile.h"

using namespace Pythia8;

using namespace std;

int main(int argc, char *argv[])
{

  istringstream nEventsStream(argv[1]);
  int nEvents;
  if (!(nEventsStream >> nEvents))
  {
    cout<<"Invalid first argument!"<<endl;
    return 0;
  }

  istringstream mEnergyStream(argv[2]);
  int mEnergy;
  if (!(mEnergyStream >> mEnergy))
  {
    cout<<"Invalid second argument!"<<endl;
    return 0;
  }
  
  cout<<"Number of events: "<<nEvents<<endl;
  cout<<"Energy: "<<mEnergy<<endl;
  cout<<"Output file: "<<argv[3]<<endl;

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

  // pythia.readString("SoftQCD:nonDiffractive = on");
  pythia.readString("SigmaTotal:zeroAXB = off");
  pythia.readString("SoftQCD:centralDiffractive = on");

                        

  //For K0
  pythia.readString("310:onMode=0");
  pythia.readString("310:OnIfMatch=211 -211");
  pythia.readString("-310:onMode=0");
  pythia.readString("-310:OnIfMatch=-211 211");
  //For Lambda
  pythia.readString("3122:onMode=0");
  pythia.readString("3122:OnIfMatch=2212 -211");
  pythia.readString("-3122:onMode=0");
  pythia.readString("-3122:OnIfMatch=-2212 211");
  pythia.init();

  // Create file on which histogram(s) can be saved.
  TFile* outFile = new TFile(argv[3], "RECREATE");

  //Useful IDs
  const int K0sPDGid = 310;
  const int K0sbarPDGid = -310;
  const int LambdaPDGid = 3122;
  const int LambdabarPDGid = -3122;
  const int piplusPDGid = 211;
  const int piminusPDGid = -211;
  const int pplusPDGid = 2212;
  const int pminusPDGid = 2212;

  //Histograms
  TH1D* K0S_Decay_R = new TH1D("K0S_Decay_R", "K^{0}_{S} decay distance;R [mm];N", 100, 0, 100);
  TH1D* K0S_Decay_R_filtered = new TH1D("K0S_Decay_R_filtered", "K^{0}_{S} decay distance, pis are filtered with min p_{T} and #eta;R [mm];N", 100, 0, 100);
  TH2D* K0S_pt_eta = new TH2D("K0S_pt_eta", "K^{0}_{S} #eta-p_{T} histogram;#eta;p_{T} [GeV]", 100, -3, 3, 100, 0, 2);
  TH2D* K0S_pt_eta_filtered = new TH2D("K0S_pt_eta_filtered", "K^{0}_{S} #eta-p_{T} histogram, pis are filtered with min p_{T} and #eta;#eta;p_{T} [GeV]", 100, -3, 3, 100, 0, 2);
  TH1D* Lambda_Decay_R = new TH1D("Lambda_Decay_R", "#Lambda decay distance;R [mm];N", 100, 0, 100);
  TH1D* Lambda_Decay_R_filtered = new TH1D("Lambda_Decay_R_filtered", "#Lambda decay distance, decay results are filtered with min p_{T} and #eta;R [mm];N", 100, 0, 100);
  TH2D* Lambda_pt_eta = new TH2D("Lambda_pt_eta", "#Lambda #eta-p_{T} histogram;#eta;p_{T} [GeV]", 100, -3, 3, 100, 0, 2);
  TH2D* Lambda_pt_eta_filtered = new TH2D("Lambda_pt_eta_filtered", "#Lambda #eta-p_{T} histogram, decay results are filtered with min p_{T} and #eta;#eta;p_{T} [GeV]", 100, -3, 3, 100, 0, 2);
  
  //Variables used in loops
  TVector3 K0S_Decay_Vertex;
  TVector3 Lambda_Decay_Vertex;
  int K0S_counter=0;
  int Lambda_counter=0;
  int K0S_filtered_counter=0;
  int Lambda_filtered_counter=0;
  int daughter_ID_1;
  int daughter_ID_2;
  TLorentzVector fourvector_1;
  TLorentzVector fourvector_2;

  // Begin event loop. Generate event; skip if generation aborted.
  for(int iEvent = 0; iEvent < nEvents; ++iEvent){
    if(!pythia.next()) continue;

    //loop over particles in event
    for(int i = 0; i < pythia.event.size(); ++i){
      if(fabs(pythia.event[i].id()) == K0sPDGid){
        K0S_counter++;
        K0S_Decay_Vertex.SetXYZ(pythia.event[i].xDec(), pythia.event[i].yDec(), pythia.event[i].zDec());
        K0S_Decay_R->Fill(K0S_Decay_Vertex.Mag());
        fourvector_1.SetPxPyPzE(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e());
        K0S_pt_eta->Fill(fourvector_1.Eta(), fourvector_1.Pt());
        //filtering
        daughter_ID_1 = pythia.event[i].daughter1();
        daughter_ID_2 = pythia.event[i].daughter2();
        fourvector_1.SetPxPyPzE(pythia.event[daughter_ID_1].px(), pythia.event[daughter_ID_1].py(), pythia.event[daughter_ID_1].pz(), pythia.event[daughter_ID_1].e());
        fourvector_2.SetPxPyPzE(pythia.event[daughter_ID_2].px(), pythia.event[daughter_ID_2].py(), pythia.event[daughter_ID_2].pz(), pythia.event[daughter_ID_2].e());
        if(fourvector_1.Pt()>0.2 && fabs(fourvector_1.Eta())<0.7 && fourvector_2.Pt()>0.2 && fabs(fourvector_2.Eta())<0.7){
          K0S_Decay_R_filtered->Fill(K0S_Decay_Vertex.Mag());
          K0S_pt_eta_filtered->Fill((fourvector_1+fourvector_2).Eta(), (fourvector_1+fourvector_2).Pt());
          K0S_filtered_counter++;
        }
      }//end if K0s
      if(fabs(pythia.event[i].id()) == LambdaPDGid){
        Lambda_counter++;
        Lambda_Decay_Vertex.SetXYZ(pythia.event[i].xDec(), pythia.event[i].yDec(), pythia.event[i].zDec());
        Lambda_Decay_R->Fill(K0S_Decay_Vertex.Mag());
        fourvector_1.SetPxPyPzE(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e());
        Lambda_pt_eta->Fill(fourvector_1.Eta(), fourvector_1.Pt());
        //filtering
        daughter_ID_1 = pythia.event[i].daughter1();
        daughter_ID_2 = pythia.event[i].daughter2();
        fourvector_1.SetPxPyPzE(pythia.event[daughter_ID_1].px(), pythia.event[daughter_ID_1].py(), pythia.event[daughter_ID_1].pz(), pythia.event[daughter_ID_1].e());
        fourvector_2.SetPxPyPzE(pythia.event[daughter_ID_2].px(), pythia.event[daughter_ID_2].py(), pythia.event[daughter_ID_2].pz(), pythia.event[daughter_ID_2].e());
        if(fourvector_1.Pt()>0.2 && fabs(fourvector_1.Eta())<0.7 && fourvector_2.Pt()>0.2 && fabs(fourvector_2.Eta())<0.7){
          Lambda_Decay_R_filtered->Fill(K0S_Decay_Vertex.Mag());
          Lambda_pt_eta_filtered->Fill((fourvector_1+fourvector_2).Eta(), (fourvector_1+fourvector_2).Pt());
          Lambda_filtered_counter++;
        }
      }//end if Lambda

    }//end particle loop
  }//end event loop
  
  outFile->cd();
  outFile->Write();
  outFile->Close();

  // Statistics on event generation.
  pythia.stat();

  //printout of K0S and Lambda efficiency
  cout<<"Produced "<<K0S_counter<<" K0S, filtered "<<K0S_filtered_counter<<" K0S, efficiency "<<double(K0S_filtered_counter)/K0S_counter*100<<"%."<<endl;
  cout<<"Produced "<<Lambda_counter<<" Lambda, filtered "<<Lambda_filtered_counter<<" Lambda, efficiency "<<double(Lambda_filtered_counter)/Lambda_counter*100<<"%."<<endl;

  // Done.
  return 0;
}
