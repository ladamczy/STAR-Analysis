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

bool can_be_detected(TLorentzVector);

int main(int argc, char *argv[]){

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
  const int pminusPDGid = -2212;

  //Histograms
  //K0
  TH1D* K0S_Decay_D = new TH1D("K0S_Decay_D", "K^{0}_{S} decay distance;R [mm];N", 100, 0, 100);
  TH1D* K0S_Decay_D_filtered = new TH1D("K0S_Decay_D_filtered", "K^{0}_{S} decay distance (3D), pis are filtered with min p_{T} and #eta;D [mm];N", 100, 0, 100);
  TH2D* K0S_pt_eta = new TH2D("K0S_pt_eta", "K^{0}_{S} #eta-p_{T} histogram;#eta;p_{T} [GeV]", 100, -3, 3, 100, 0, 2);
  TH2D* K0S_pt_eta_filtered = new TH2D("K0S_pt_eta_filtered", "K^{0}_{S} #eta-p_{T} histogram, pis are filtered with min p_{T} and #eta;#eta;p_{T} [GeV]", 100, -3, 3, 100, 0, 2);
  TH2D* K0S_detection = new TH2D();
  //Lambda
  TH1D* Lambda_Decay_D = new TH1D("Lambda_Decay_D", "#Lambda decay distance;R [mm];N", 100, 0, 100);
  TH1D* Lambda_Decay_D_filtered = new TH1D("Lambda_Decay_D_filtered", "#Lambda decay distance (3D), decay results are filtered with min p_{T} and #eta;D [mm];N", 100, 0, 100);
  TH2D* Lambda_pt_eta = new TH2D("Lambda_pt_eta", "#Lambda #eta-p_{T} histogram;#eta;p_{T} [GeV]", 100, -3, 3, 100, 0, 2);
  TH2D* Lambda_pt_eta_filtered = new TH2D("Lambda_pt_eta_filtered", "#Lambda #eta-p_{T} histogram, decay results are filtered with min p_{T} and #eta;#eta;p_{T} [GeV]", 100, -3, 3, 100, 0, 2);
  TH2D* Lambda_detection = new TH2D();
  //other
  TH1D* M_inv_pion_pairs = new TH1D("M_inv_pion_pairs", "Pion pairs;m_{inv} [GeV];Number of pairs", 100, 0.42, 0.56);

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
  vector<TLorentzVector> piplus_tab;
  vector<TLorentzVector> piminus_tab;

  // Begin event loop. Generate event; skip if generation aborted.
  for(int iEvent = 0; iEvent < nEvents; ++iEvent){
    if(!pythia.next()) continue;

    //loop over particles in event
    for(int i = 0; i < pythia.event.size(); ++i){
      if(fabs(pythia.event[i].id()) == K0sPDGid){
        K0S_counter++;
        K0S_Decay_Vertex.SetXYZ(pythia.event[i].xDec(), pythia.event[i].yDec(), pythia.event[i].zDec());
        K0S_Decay_D->Fill(K0S_Decay_Vertex.Mag());
        fourvector_1.SetPxPyPzE(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e());
        K0S_pt_eta->Fill(fourvector_1.Eta(), fourvector_1.Pt());
        //filtering
        daughter_ID_1 = pythia.event[i].daughter1();
        daughter_ID_2 = pythia.event[i].daughter2();
        fourvector_1.SetPxPyPzE(pythia.event[daughter_ID_1].px(), pythia.event[daughter_ID_1].py(), pythia.event[daughter_ID_1].pz(), pythia.event[daughter_ID_1].e());
        fourvector_2.SetPxPyPzE(pythia.event[daughter_ID_2].px(), pythia.event[daughter_ID_2].py(), pythia.event[daughter_ID_2].pz(), pythia.event[daughter_ID_2].e());
        if(can_be_detected(fourvector_1) && can_be_detected(fourvector_2)){
          K0S_Decay_D_filtered->Fill(K0S_Decay_Vertex.Mag());
          K0S_pt_eta_filtered->Fill((fourvector_1+fourvector_2).Eta(), (fourvector_1+fourvector_2).Pt());
          K0S_filtered_counter++;
        }
      }//end if K0s
      if(fabs(pythia.event[i].id()) == LambdaPDGid){
        Lambda_counter++;
        Lambda_Decay_Vertex.SetXYZ(pythia.event[i].xDec(), pythia.event[i].yDec(), pythia.event[i].zDec());
        Lambda_Decay_D->Fill(Lambda_Decay_Vertex.Mag());
        fourvector_1.SetPxPyPzE(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e());
        Lambda_pt_eta->Fill(fourvector_1.Eta(), fourvector_1.Pt());
        //filtering
        daughter_ID_1 = pythia.event[i].daughter1();
        daughter_ID_2 = pythia.event[i].daughter2();
        fourvector_1.SetPxPyPzE(pythia.event[daughter_ID_1].px(), pythia.event[daughter_ID_1].py(), pythia.event[daughter_ID_1].pz(), pythia.event[daughter_ID_1].e());
        fourvector_2.SetPxPyPzE(pythia.event[daughter_ID_2].px(), pythia.event[daughter_ID_2].py(), pythia.event[daughter_ID_2].pz(), pythia.event[daughter_ID_2].e());
        if(can_be_detected(fourvector_1) && can_be_detected(fourvector_2)){
          Lambda_Decay_D_filtered->Fill(Lambda_Decay_Vertex.Mag());
          Lambda_pt_eta_filtered->Fill((fourvector_1+fourvector_2).Eta(), (fourvector_1+fourvector_2).Pt());
          Lambda_filtered_counter++;
        }
      }//end if Lambda
      if(fabs(pythia.event[i].id()) == piplusPDGid){
        fourvector_1.SetPxPyPzE(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e());
        if(can_be_detected(fourvector_1)){
          if(pythia.event[i].id()>0){
            piplus_tab.push_back(fourvector_1);
          }else if(pythia.event[i].id()<0){
            piminus_tab.push_back(fourvector_1);
          } 
        }
      }//end if pi+-
    }//end particle loop

    //pion mixing
    for (size_t i = 0; i < piplus_tab.size(); i++){
      for (size_t j = 0; j < piminus_tab.size(); j++){
        M_inv_pion_pairs->Fill((piplus_tab[i]+piminus_tab[j]).M());
      }
    }
    
    //cleanup
    piplus_tab.clear();
    piminus_tab.clear();
  }//end event loop

  //Fixing last histograms
  K0S_detection = (TH2D*)K0S_pt_eta_filtered->Clone();
  K0S_detection->SetTitle("Probability of detecting K^{0}_{S} with given #eta and p_{T}");
  K0S_detection->SetName("K0S_detection");
  K0S_detection->GetXaxis()->SetTitle("#eta");
  K0S_detection->GetYaxis()->SetTitle("p_{T} [GeV]");
  K0S_detection->Divide(K0S_pt_eta);
  Lambda_detection = (TH2D*)Lambda_pt_eta_filtered->Clone();
  Lambda_detection->SetTitle("Probability of detecting #Lambda with given #eta and p_{T}");
  Lambda_detection->SetName("Lambda_detection");
  Lambda_detection->GetXaxis()->SetTitle("#eta");
  Lambda_detection->GetYaxis()->SetTitle("p_{T} [GeV]");
  Lambda_detection->Divide(Lambda_pt_eta);
  
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

bool can_be_detected(TLorentzVector tested_particle){
  if(fabs(tested_particle.Eta())<0.7 && tested_particle.Pt()>0.2){
    return true;
  }
  return false;
}