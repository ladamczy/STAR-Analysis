#include "Includes.h"
#include <iostream> //
#include <string>    
#include <utility>
#include <sstream> 
#include <algorithm> 
#include <stdio.h> 
#include <stdlib.h> 
#include <vector> 
#include <fstream> 
#include <cmath> 
#include <cstdlib>
#include <sys/stat.h>
#include <iterator>
#include <ostream>
#include <TRandom.h>
#include <iomanip>
#include <stdexcept>
#include <limits>
#include "TROOT.h"
#include "TSystem.h"
#include "TThread.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include "TProfile.h"
#include <TH2.h> 
#include <TF1.h> 
#include <TF2.h> 
#include <THStack.h> 
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TStyle.h> 
#include <TGraph.h> 
#include <TGraph2D.h> 
#include <TGraphErrors.h> 
#include <TCanvas.h> 
#include <TLegend.h> 
#include <TGaxis.h> 
#include <TString.h> 
#include <TColor.h> 
#include <TLine.h> 
#include <TExec.h> 
#include <TFitResultPtr.h> 
#include <TFitResult.h> 
#include <TLatex.h> 
#include <TMath.h>
#include <TLorentzVector.h>
#include <ROOT/TThreadedObject.hxx>
#include <TTreeReader.h>
#include <ROOT/TTreeProcessorMT.hxx>
#include "StRPEvent.h"
#include "StUPCRpsTrack.h"
#include "StUPCRpsTrackPoint.h"
#include "StUPCEvent.h"
#include "StUPCTrack.h"
#include "StUPCBemcCluster.h"
#include "StUPCVertex.h"
#include "StUPCTofHit.h"
#include "StUPCV0.h"
#include "BeamPosition.h"
#include <Afterburner.h>
#include <StRPEvent.h>
#include <map>

using namespace std;

int main(int argc, char** argv)  
{

    ifstream inputFilePathList(argv[1]);
    if (!inputFilePathList) 
    {
        cerr << "Failed to open input file." << std::endl;
        return 1;
    }

    TChain *chain = new TChain("mUPCTree"); 

    string inputFileName;
    vector <string> rootFiles;
    while (std::getline(inputFilePathList, inputFileName))
    {
        chain->Add(inputFileName.c_str());
    }
    inputFilePathList.close();

    static StUPCEvent *upcEvt = 0x0;
    chain->SetBranchAddress("mUPCEvent", &upcEvt);

    TH1D* HistKaonE1 = new TH1D("HistKaonE1", "; E [GeV]; # events", 100 ,0, 3);
    TH1D* HistKaonE2 = new TH1D("HistKaonE2", "; E [GeV]; # events", 100 ,0, 3);
    TH1D* HistKaonM1 = new TH1D("HistKaonM1", "; m [GeV]; # events", 100 ,0.4, 0.6);
    TH1D* HistKaonM2 = new TH1D("HistKaonM2", "; m [GeV]; # events", 100 ,0.4, 0.6);
    TH1D* HistKaonPt1 = new TH1D("HistKaonPt1", "; pT [GeV]; # events", 100 ,0, 3);
    TH1D* HistKaonPt2 = new TH1D("HistKaonPt2", "; pT [GeV]; # events", 100 ,0, 3);
    TH1D* HistKaonEta1 = new TH1D("HistKaonEta1", "; #eta; # events", 100 ,-3.5, 3.5);
    TH1D* HistKaonEta2 = new TH1D("HistKaonEta2", "; #eta; # events", 100 ,-3.5, 3.5);

    TH1D* HistKaonTheta1 = new TH1D("HistKaonTheta1", "; #theta [rad]; # events", 51 ,0 , M_PI);
    TH1D* HistKaonTheta2 = new TH1D("HistKaonTheta2", "; #theta [rad]; # events", 51 ,0 , M_PI);

    TH1D* HistKaonPhi1 = new TH1D("HistKaonPhi1", "; #phi [rad]; # events", 51, -M_PI, M_PI);
    TH1D* HistKaonPhi2 = new TH1D("HistKaonPhi2", "; #phi [rad]; # events", 51, -M_PI, M_PI);

    TParticle* particle;
    TParticle* PosPion1; TParticle* NegPion1; TParticle*  PosPion2; TParticle*  NegPion2;
    vector <TParticle*> PosPions, NegPions;
    TLorentzVector lorentzVectorPosPion1, lorentzVectorNegPion1, lorentzVectorPosPion2, lorentzVectorNegPion2;
    TLorentzVector lorentzVectorNeutKaon1, lorentzVectorNeutKaon2;  


    for (Long64_t i = 0; i < chain->GetEntries(); ++i) 
    {
        chain->GetEntry(i);

        for (int i = 0; i < upcEvt->getNumberOfMCParticles(); i ++)
        {
            particle = upcEvt->getMCParticle(i);

            if (particle->GetPDG()->PdgCode() == 211)
            {
                PosPions.push_back(particle);
            }

            else if (particle->GetPDG()->PdgCode() == -211)
            {
                NegPions.push_back(particle);
            }
        }

        if (PosPions.size() == 2 and NegPions.size() == 2)
        {
            if (PosPions[0]->GetFirstMother() == NegPions[0]->GetFirstMother())
            {
                PosPion1 = PosPions[0];
                PosPion2 = PosPions[1];
                NegPion1 = NegPions[0];
                NegPion2 = NegPions[1];
            }

            else 
            {
                PosPion1 = PosPions[0];
                PosPion2 = PosPions[1];
                NegPion1 = NegPions[1];
                NegPion2 = NegPions[0];
            }

            lorentzVectorPosPion1.SetPxPyPzE(PosPion1->Px(), PosPion1->Py(), PosPion1->Pz(), PosPion1->Energy());
            lorentzVectorPosPion2.SetPxPyPzE(PosPion2->Px(), PosPion2->Py(), PosPion2->Pz(), PosPion2->Energy());
            lorentzVectorNegPion1.SetPxPyPzE(NegPion1->Px(), NegPion1->Py(), NegPion1->Pz(), NegPion1->Energy());
            lorentzVectorNegPion2.SetPxPyPzE(NegPion2->Px(), NegPion2->Py(), NegPion2->Pz(), NegPion2->Energy());

            lorentzVectorNeutKaon1 = lorentzVectorPosPion1 + lorentzVectorNegPion1;
            lorentzVectorNeutKaon2 = lorentzVectorPosPion2 + lorentzVectorNegPion2;

            HistKaonE1->Fill(lorentzVectorNeutKaon1.E());
            HistKaonM1->Fill(lorentzVectorNeutKaon1.M());
            HistKaonPt1->Fill(lorentzVectorNeutKaon1.Pt());
            HistKaonEta1->Fill(lorentzVectorNeutKaon1.Eta());
            HistKaonTheta1->Fill(lorentzVectorNeutKaon1.Theta());
            HistKaonPhi1->Fill(lorentzVectorNeutKaon1.Phi());

            HistKaonE2->Fill(lorentzVectorNeutKaon2.E());
            HistKaonM2->Fill(lorentzVectorNeutKaon2.M());
            HistKaonPt2->Fill(lorentzVectorNeutKaon2.Pt());
            HistKaonEta2->Fill(lorentzVectorNeutKaon2.Eta());
            HistKaonTheta2->Fill(lorentzVectorNeutKaon2.Theta());
            HistKaonPhi2->Fill(lorentzVectorNeutKaon2.Phi());

        }

        PosPions.clear();
        NegPions.clear();
    }
    
    TFile *outfile = TFile::Open(argv[2], "recreate"); 

    HistKaonE1->Write();
    HistKaonM1->Write();
    HistKaonPt1->Write();
    HistKaonEta1->Write();
    HistKaonTheta1->Write();
    HistKaonPhi1->Write();
    HistKaonE2->Write();
    HistKaonM2->Write();
    HistKaonPt2->Write();
    HistKaonEta2->Write();
    HistKaonTheta2->Write();
    HistKaonPhi2->Write();

    outfile->Close();

    return 0;
}
