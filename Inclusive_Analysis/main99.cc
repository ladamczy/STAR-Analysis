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

  
  // Create Pythia instance and set it up to generate hard QCD processes
  // above pTHat = 20 GeV for pp collisions at 14 TeV.
  Pythia pythia;
  //pythia.readString("Beams:eCM = 510"); //beam energy in GeV
  if( mEnergy == 510) pythia.readString("Beams:eCM = 510"); //beam energy in GeV
  else if( mEnergy == 200) pythia.readString("Beams:eCM = 200");
  else
  {
    cout<<"Invalid input energy!"<<endl;
    return 0;
  
  }
  pythia.readString("SoftQCD:nonDiffractive = on"); //equivalent to MB
//  pythia.readString("310:onMode=0");
//  pythia.readString("310:OnIfMatch=211 -211");
//  pythia.readString("-310:onMode=0");
//  pythia.readString("-310:OnIfMatch=-211 211");
  pythia.readString("3122:onMode=0");
  pythia.readString("3122:OnIfMatch=2212 -211");
  pythia.readString("-3122:onMode=0");
  pythia.readString("-3122:OnIfMatch=-2212 211");
  pythia.init();

  // Create file on which histogram(s) can be saved.
  TFile* outFile = new TFile(argv[3], "RECREATE");

  // Book histogram.
  //TH1F *mult = new TH1F("mult","charged multiplicity", 100, -0.5, 799.5);
  
  const int nPtBins = 8;
  float const pT_bins[nPtBins+1] = { 0., 0.5, 1.,1.5, 2., 2.5, 3., 4., 5.};
  
  const int nPtBins_corr = 2;
  float const pT_bins_corr[nPtBins_corr+1] = { 0.5, 1.5, 5.};

  const int nEtaBins = 3;
  float const eta_bins[nEtaBins+1] = { -1, -0.4, 0.4, 1 };

//  const int K0sPDGid = 310;
//  const int K0sbarPDGid = -310;
  const int K0sPDGid = 3122;
  const int K0sbarPDGid = -3122;



  //variables for event and particle loop
  TLorentzVector K0_fourmom;
  TLorentzVector K0_fourmom_reverse; //for proton boost
  
  TVector3 K0s_production_vertex;
  TVector3 K0s_decay_vertex;

  TLorentzVector pi1_fourmom;
  TLorentzVector pi2_fourmom;


  //histograms
  TH1D* K0s_pt_hist = new TH1D("K0s_pt_hist","K0s_pt_hist",100,0,10);
  TH1D* K0s_pt_primary_hist = new TH1D("K0s_pt_primary_hist","K0s_pt_primary_hist",100,0,10);
  TH1D* K0s_p_hist = new TH1D("K0s_p_hist","K0s_p_hist",100,0,10);
  TH1D* K0s_eta_hist = new TH1D("K0s_eta_hist","K0s_eta_hist",100,-2,2);

//  TH1D* K0s_mass_hist = new TH1D("K0s_mass_hist","K0s_mass_hist",100,0.3,0.8);
  TH1D* K0s_mass_hist = new TH1D("K0s_mass_hist","K0s_mass_hist",100,1.0,1.4);
  TH1D* K0s_decayK0s_hist = new TH1D("K0s_decayK0s_hist", "K0s_decayK0s_hist", 100, 0, 500);
  TH1D* K0s_TdecayK0s_hist = new TH1D("K0s_TdecayK0s_hist", "K0s_TdecayK0s_hist", 100, 0, 500);  
	
  TH1D *K0s_thetaProdPlane[nPtBins+1][nEtaBins+1];
  TH1D *K0s_cosThetaProdPlane[nPtBins+1][nEtaBins+1];
  
  TH2D *K0s_y_vs_p_eta[nPtBins+1];
  TH2D *K0s_y_vs_pi_eta[nPtBins+1];
  
  TH1D *K0s_K0s_cosThetaProdPlane = new TH1D("K0s_K0s_cosThetaProdPlane", "K0s_K0s_cosThetaProdPlane", 10, -1, 1);
  
  TH1D *K0s_K0s_cosThetaProdPlane_pT_hist[nPtBins_corr][nPtBins_corr];

  TH1D *K0s_K0s_cosThetaProdPlane_eta_hist[nEtaBins][nEtaBins];
  
  
  //histograms after cuts
  TH1D *K0s_thetaProdPlane_cuts[nPtBins+1][nEtaBins+1];
  TH1D *K0s_cosThetaProdPlane_cuts[nPtBins+1][nEtaBins+1];
  
  TH2D *K0s_y_vs_p_eta_cuts[nPtBins+1];
  TH2D *K0s_y_vs_pi_eta_cuts[nPtBins+1];
  
  TH1D *K0s_K0s_cosThetaProdPlane_cuts = new TH1D("K0s_K0s_cosThetaProdPlane_cuts", "K0s_K0s_cosThetaProdPlane_cuts", 10, -1, 1);
  
  TH1D *K0s_K0s_cosThetaProdPlane_pT_cuts_hist[nPtBins_corr][nPtBins_corr];

  TH1D *K0s_K0s_cosThetaProdPlane_eta_cuts_hist[nEtaBins][nEtaBins];


  
  
  for(unsigned int pTbin = 0; pTbin < nPtBins+1; pTbin++)
  {
    K0s_y_vs_p_eta[pTbin] = new TH2D(Form("K0s_y_vs_p_eta_pT_%i", pTbin), Form("K0s_y_vs_p_eta_pT_%i", pTbin), 100, -1, 1, 500, -5, 5);
    K0s_y_vs_pi_eta[pTbin] = new TH2D(Form("K0s_y_vs_pi_eta_pT_%i", pTbin), Form("K0s_y_vs_pi_eta_pT_%i", pTbin), 100, -1, 1, 500, -5, 5);
    
    K0s_y_vs_p_eta_cuts[pTbin] = new TH2D(Form("K0s_y_vs_p_eta_cuts_pT_%i", pTbin), Form("K0s_y_vs_p_eta_cuts_pT_%i", pTbin), 100, -1, 1, 500, -5, 5);
    K0s_y_vs_pi_eta_cuts[pTbin] = new TH2D(Form("K0s_y_vs_pi_eta_cuts_pT_%i", pTbin), Form("K0s_y_vs_pi_eta_cuts_pT_%i", pTbin), 100, -1, 1, 500, -5, 5);
    
    
  
    for(unsigned int etaBin = 0; etaBin < nEtaBins+1; etaBin++)
    {

      K0s_thetaProdPlane[pTbin][etaBin] = new TH1D(Form("K0s_thetaProdPlane_pT_%i_eta_%i", pTbin, etaBin), Form("K0s_thetaProdPlane_pT_%i_eta_%i", pTbin, etaBin), 20, 0, TMath::Pi());
      K0s_cosThetaProdPlane[pTbin][etaBin] = new TH1D(Form("K0s_cosThetaProdPlane_pT_%i_eta_%i", pTbin, etaBin), Form("K0s_cosThetaProdPlane_pT_%i_eta_%i", pTbin, etaBin), 20, -1, 1);
      
      K0s_thetaProdPlane_cuts[pTbin][etaBin] = new TH1D(Form("K0s_thetaProdPlane_cuts_pT_%i_eta_%i", pTbin, etaBin), Form("K0s_thetaProdPlane_cuts_pT_%i_eta_%i", pTbin, etaBin), 20, 0, TMath::Pi());
      K0s_cosThetaProdPlane_cuts[pTbin][etaBin] = new TH1D(Form("K0s_cosThetaProdPlane_cuts_pT_%i_eta_%i", pTbin, etaBin), Form("K0s_cosThetaProdPlane_cuts_pT_%i_eta_%i", pTbin, etaBin), 20, -1, 1);
    
    }

  }
  
  
  //L-L correlation histograms in bins
  for(unsigned int pTbin1 = 0; pTbin1 < nPtBins_corr; pTbin1++)
  {
    for(unsigned int pTbin2 = 0; pTbin2 < nPtBins_corr; pTbin2++)
    {
      
      K0s_K0s_cosThetaProdPlane_pT_hist[pTbin1][pTbin2] = new TH1D(Form("K0s_K0s_cosThetaProdPlane_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("K0s_K0s_cosThetaProdPlane_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      K0s_K0s_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2] = new TH1D(Form("K0s_K0s_cosThetaProdPlane_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("K0s_K0s_cosThetaProdPlane_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);

    }
  }

  for(unsigned int etaBin1 = 0; etaBin1 < nEtaBins; etaBin1++)
  {
    for(unsigned int etaBin2 = 0; etaBin2 < nEtaBins; etaBin2++)
    {
      
      K0s_K0s_cosThetaProdPlane_eta_hist[etaBin1][etaBin2] = new TH1D(Form("K0s_K0s_cosThetaProdPlane_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("K0s_K0s_cosThetaProdPlane_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      K0s_K0s_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2] = new TH1D(Form("K0s_K0s_cosThetaProdPlane_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("K0s_K0s_cosThetaProdPlane_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);

    }
  }

  
  // Begin event loop. Generate event; skip if generation aborted.
  for (int iEvent = 0; iEvent < nEvents; ++iEvent)
  {
    if (!pythia.next()) continue;
    
    vector<int> K0s_pT_bin_vector;
    vector<int> K0s_pT_bin_cuts_vector;
    
    vector<int> K0s_eta_bin_vector;
    vector<int> K0s_eta_bin_cuts_vector;
    
    vector<TLorentzVector> pi_star_vector;
    vector<TLorentzVector> pi_star_cuts_vector;

    //loop over particles in event
    for (int i = 0; i < pythia.event.size(); ++i)
    {
      if( fabs(pythia.event[i].id()) == K0sPDGid)
      {
        //cout<<"Lambda found!"<<endl;
        
        //find index of decay daughters
        const int daughter1_Id = pythia.event[i].daughter1();
        const int daughter2_Id = pythia.event[i].daughter2();
        
        
        //Lambda fourmomentum
        K0_fourmom.SetPxPyPzE(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e());
        
        //Lambda MC decay veretex
        K0s_decay_vertex.SetXYZ(pythia.event[i].xDec(), pythia.event[i].yDec(), pythia.event[i].zDec());
        
        //Lambda MC production veretex (PV)
        K0s_production_vertex.SetXYZ(pythia.event[i].xProd(), pythia.event[i].yProd(), pythia.event[i].zProd());
        
        TVector3 K0s_decayK0s_vect = K0s_decay_vertex - K0s_production_vertex;
        float K0s_decayL = K0s_decayK0s_vect.Mag();
        float K0s_decayR = K0s_decayK0s_vect.Perp();        

//        K0s_decayK0s_hist->Fill(K0s_decayL);
        
//        if( fabs(K0_fourmom.Rapidity()) >= 1.) continue;
        
        
        K0_fourmom_reverse.SetPxPyPzE(-pythia.event[i].px(), -pythia.event[i].py(), -pythia.event[i].pz(), pythia.event[i].e());
        

        pi1_fourmom.SetPxPyPzE(pythia.event[daughter1_Id].px(), pythia.event[daughter1_Id].py(), pythia.event[daughter1_Id].pz(), pythia.event[daughter1_Id].e());
        pi2_fourmom.SetPxPyPzE(pythia.event[daughter2_Id].px(), pythia.event[daughter2_Id].py(), pythia.event[daughter2_Id].pz(), pythia.event[daughter2_Id].e());
        


        if( fabs(pi1_fourmom.Rapidity()) >= 1.) continue;
       	if( fabs(pi2_fourmom.Rapidity()) >= 1.) continue;
        K0s_eta_hist->Fill(K0_fourmom.Rapidity());
                
        if( fabs(K0_fourmom.Rapidity()) >= 1.) continue;

        K0s_decayK0s_hist->Fill(K0s_decayL);
        K0s_TdecayK0s_hist->Fill(K0s_decayR);
        K0s_pt_hist->Fill(K0_fourmom.Perp());
        K0s_p_hist->Fill(K0_fourmom.Beta()*K0_fourmom.Gamma());
        K0s_mass_hist->Fill(K0_fourmom.M());
        if(K0s_decayR<30) K0s_pt_primary_hist->Fill(K0_fourmom.Perp());



        TLorentzVector pi1_fourmom_star = pi1_fourmom;
        pi1_fourmom_star.Boost(K0_fourmom_reverse.BoostVector());  

        TVector3 beamVector(0.,0.,1.); //unity vector along the beam axis
        TVector3 mProdPlane = beamVector.Cross(K0_fourmom.Vect());
        mProdPlane = ( mProdPlane )*(1./mProdPlane.Mag() );

        float mThetaProdPlane = mProdPlane.Angle(pi1_fourmom_star.Vect());
        
        //fill all histograms for all pT and centrality bins
        int pT_bin = -1;

        //find pT bin of Lambda
        for(int j = 0; j < nPtBins; j++) //loop over pT bins
        {
          if(K0_fourmom.Pt() > pT_bins[j] && K0_fourmom.Pt() <= pT_bins[j+1])
          {
            pT_bin = j;
            break; //stop after pT bin is found
          }
        }

        if( pT_bin == -1 ) continue;
        
        
        int pT_bin_corr = -1;

        //find pT bin of Lambda
        for(int j = 0; j < nPtBins_corr; j++) //loop over pT bins
        {
          if(K0_fourmom.Pt() > pT_bins_corr[j] && K0_fourmom.Pt() <= pT_bins_corr[j+1])
          {
            pT_bin_corr = j;
            break; //stop after pT bin is found
          }
        }
        
        if( pT_bin_corr == -1 ) continue;


        //fill all histograms for all eta and centrality bins
        int eta_bin = -1;

        //find eta bin of Lambda
        for(int j = 0; j < nEtaBins; j++) //loop over eta bins
        {
          if(K0_fourmom.Eta() > eta_bins[j] && K0_fourmom.Eta() <= eta_bins[j+1])
          {
            eta_bin = j;
            break; //stop after eta bin is found
          }
        }

        if( eta_bin == -1 ) continue;
        //_____________________________________________________________________________________________
        
        
        K0s_pT_bin_vector.push_back(pT_bin_corr);
        K0s_eta_bin_vector.push_back(eta_bin);
        
        pi_star_vector.push_back(pi1_fourmom_star);
       
//        K0s_mass_hist->Fill(pythia.event[i].m());
//        K0s_pt_hist->Fill(K0_fourmom.Pt());

        K0s_thetaProdPlane[pT_bin][eta_bin]->Fill(mThetaProdPlane);
        K0s_thetaProdPlane[nPtBins][eta_bin]->Fill(mThetaProdPlane);  //pT integrated, eta bins
        K0s_thetaProdPlane[pT_bin][nEtaBins]->Fill(mThetaProdPlane);  //pT bins, -1 < eta < 1
        K0s_thetaProdPlane[nPtBins][nEtaBins]->Fill(mThetaProdPlane); //pT integrated and -1 < eta < 1


        K0s_cosThetaProdPlane[pT_bin][eta_bin]->Fill(cos(mThetaProdPlane));
        K0s_cosThetaProdPlane[nPtBins][eta_bin]->Fill(cos(mThetaProdPlane));
        K0s_cosThetaProdPlane[pT_bin][nEtaBins]->Fill(cos(mThetaProdPlane));
        K0s_cosThetaProdPlane[nPtBins][nEtaBins]->Fill(cos(mThetaProdPlane));
        
        K0s_y_vs_p_eta[pT_bin]->Fill(K0_fourmom.Rapidity(), pi1_fourmom.Eta());
        K0s_y_vs_pi_eta[pT_bin]->Fill(K0_fourmom.Rapidity(), pi2_fourmom.Eta());
        
        K0s_y_vs_p_eta[nPtBins]->Fill(K0_fourmom.Rapidity(), pi1_fourmom.Eta());
        K0s_y_vs_pi_eta[nPtBins]->Fill(K0_fourmom.Rapidity(), pi2_fourmom.Eta());
        //_____________________________________________________________________________________________
        
        
        //fill histograms after cuts
        
        //K0s cuts
        //decay length cuts in mm (default PYTHIA units)
        if(K0s_decayL < 5 || K0s_decayL > 250) continue;
        
        
        //daughter cuts
        if( fabs(pi1_fourmom.Eta()) >= 1. || fabs(pi2_fourmom.Eta()) >= 1. ) continue;       
        if( pi1_fourmom.Pt() < 0.15 || pi1_fourmom.Pt() > 20. ) continue;
        if( pi2_fourmom.Pt() < 0.15 || pi2_fourmom.Pt() > 20. ) continue;
        
        //add decay length cut on K0s
        
        pi_star_cuts_vector.push_back(pi1_fourmom_star);
        
        K0s_thetaProdPlane_cuts[pT_bin][eta_bin]->Fill(mThetaProdPlane);
        K0s_thetaProdPlane_cuts[nPtBins][eta_bin]->Fill(mThetaProdPlane);  //pT integrated, eta bins
        K0s_thetaProdPlane_cuts[pT_bin][nEtaBins]->Fill(mThetaProdPlane);  //pT bins, -1 < eta < 1
        K0s_thetaProdPlane_cuts[nPtBins][nEtaBins]->Fill(mThetaProdPlane); //pT integrated and -1 < eta < 1


        K0s_cosThetaProdPlane_cuts[pT_bin][eta_bin]->Fill(cos(mThetaProdPlane));
        K0s_cosThetaProdPlane_cuts[nPtBins][eta_bin]->Fill(cos(mThetaProdPlane));
        K0s_cosThetaProdPlane_cuts[pT_bin][nEtaBins]->Fill(cos(mThetaProdPlane));
        K0s_cosThetaProdPlane_cuts[nPtBins][nEtaBins]->Fill(cos(mThetaProdPlane));
        
        K0s_y_vs_p_eta_cuts[pT_bin]->Fill(K0_fourmom.Rapidity(), pi1_fourmom.Eta());
        K0s_y_vs_pi_eta_cuts[pT_bin]->Fill(K0_fourmom.Rapidity(), pi2_fourmom.Eta());
        
        K0s_y_vs_p_eta_cuts[nPtBins]->Fill(K0_fourmom.Rapidity(), pi1_fourmom.Eta());
        K0s_y_vs_pi_eta_cuts[nPtBins]->Fill(K0_fourmom.Rapidity(), pi2_fourmom.Eta());        
        

      }//end if K0s
    
    }//end particle loop
    
    //fill K0s-K0s correlation before cuts
    if(pi_star_vector.size() > 1)
    {
      for(unsigned int iK0s1 = 0; iK0s1 < pi_star_vector.size(); iK0s1++)
      {
        for(unsigned int iK0s2 = iK0s1+1; iK0s2 < pi_star_vector.size(); iK0s2++)
        {
          double theta_star = pi_star_vector.at(iK0s1).Angle(pi_star_vector.at(iK0s2).Vect());
          
          K0s_K0s_cosThetaProdPlane->Fill(TMath::Cos(theta_star));
          K0s_K0s_cosThetaProdPlane_pT_hist[K0s_pT_bin_vector.at(iK0s1)][K0s_pT_bin_vector.at(iK0s2)]->Fill(TMath::Cos(theta_star));
          K0s_K0s_cosThetaProdPlane_eta_hist[K0s_eta_bin_vector.at(iK0s1)][K0s_eta_bin_vector.at(iK0s2)]->Fill(TMath::Cos(theta_star));
        
        }
      
      }  
    
    }
    
    //fill K0s-K0s correlation after cuts
    if(pi_star_cuts_vector.size() > 1)
    {
      for(unsigned int iK0s1 = 0; iK0s1 < pi_star_cuts_vector.size(); iK0s1++)
      {
        for(unsigned int iK0s2 = iK0s1+1; iK0s2 < pi_star_cuts_vector.size(); iK0s2++)
        {
          double theta_star = pi_star_cuts_vector.at(iK0s1).Angle(pi_star_cuts_vector.at(iK0s2).Vect());
          
          K0s_K0s_cosThetaProdPlane_cuts->Fill(TMath::Cos(theta_star));
          K0s_K0s_cosThetaProdPlane_pT_cuts_hist[K0s_pT_bin_vector.at(iK0s1)][K0s_pT_bin_vector.at(iK0s2)]->Fill(TMath::Cos(theta_star));
          K0s_K0s_cosThetaProdPlane_eta_cuts_hist[K0s_eta_bin_vector.at(iK0s1)][K0s_eta_bin_vector.at(iK0s2)]->Fill(TMath::Cos(theta_star));
        
        }
      
      }  
    
    }
    
    pi_star_vector.clear();
    pi_star_cuts_vector.clear();

  }//end event loop
  
  
  
  outFile->cd();
  outFile->Write();
  outFile->Close();

  // Statistics on event generation.
  pythia.stat();
  
  //outFile->Close();

  // Done.
  return 0;
}
