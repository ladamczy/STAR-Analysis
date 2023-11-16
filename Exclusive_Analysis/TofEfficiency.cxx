#include <iostream>
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
#include "MatchFillPosition.h"
#include "ReadFillPositionFile.h"
using namespace std;

void FindProtons(bool isMC, StRPEvent *rpEvt, StUPCEvent *upcEvt, TVector3 & proton1, TVector3 & proton2);
void FillCutflow(int *i, TH1D *hist );

int main(int argc, char** argv)  
{

    float massProton =  938.272/1000.0;
    float massPion = 0.13957061;
    double massKaon =  497.611/1000.0;
    double massLambda = 1115.0/1000.0;
    TVector3 const tryVec(0,0,0);
    
    ifstream inputFilePathList(argv[1]);
    if (!inputFilePathList) 
    {
        cerr << "Failed to open input file." << std::endl;
        return 1;
    }

    TChain *chain = new TChain("mUPCTree"); 

    string inputFileName;
    while (std::getline(inputFilePathList, inputFileName))
    {
        chain->Add(inputFileName.c_str());
    }
    inputFilePathList.close();

    bool isMC = 0;
    static StUPCEvent *upcEvt = nullptr;
    static StRPEvent *rpEvt = nullptr;
    chain->SetBranchAddress("mUPCEvent", &upcEvt);
    chain->GetEntry(0);

    if (upcEvt->getRunNumber() == 1)
    {
        isMC = 1;
    }

    if (isMC == 0)
    {
        chain->SetBranchAddress("mRPEvent", &rpEvt);
    }

    cout << isMC << endl;
    vector <double> xyVec;
    vector <double> xVec;
    vector <double> yVec;
    vector <double> zVec;

    double beamPositionX, beamPositionY;
    double primVertexPosX, primVertexPosY, primVertexPosZ, primVertexPosErrX, primVertexPosErrY, primVertexPosErrZ;
    double distVertexBeamX, distVertexBeamY, distR, significance;
    int iCutFlow;

    const vector <int> triggerID = {570209, 570219, 570229, 570701, 570702, 570703, 570704, 570705, 570709, 570711, 570712, 570719, 590701, 590703, 590705, 590708, 590709};
    const vector <int> triggerCEP = {570701, 570705, 570711};

    // to be changed... there is no warning in case of the invalid input file
	vector <vector<double>> fillNumberWithPosition = ReadFillPositionData("../share/Run7PolarizationWithPosition.csv");
    TH1D* HistCutFlow = new TH1D("HistCutFlow", ";;count", 16, -0.5, 15.5);
    TH1D* HistKaonMassProbeWithoutTof = new TH1D("HistKaonMassProbeWithoutTof", "; m_{#pi^{+}#pi^{-}}^{tag} [GeV]; # events", 25 ,0.44, 0.54);  
    TH1D* HistKaonMassProbeWithTof = new TH1D("HistKaonMassProbeWithTof", "; m_{#pi^{+}#pi^{-}}^{probe} [GeV]; # events", 25 ,0.44, 0.54);

    for (Long64_t i = 0; i < chain->GetEntries(); ++i) 
    {
        chain->GetEntry(i);


        TLorentzVector trackVector;
        TVector3 proton1, proton2;

        if (i%100000 == 0)  
        {
            cout << i << "/" <<  chain->GetEntries() << endl;
        }

        iCutFlow = 0;
        FillCutflow(&iCutFlow, HistCutFlow); // cutflow: 0

        // select intact protons 
        FindProtons(isMC, rpEvt, upcEvt, proton1, proton2);
        double protonSignsY = proton1.Y() * proton2.Y();

        // select number of primary vertices per event
        Int_t numberOfPrimaryVertices = upcEvt->getNumberOfVertices();       	     
        if(numberOfPrimaryVertices!=1)
        {
            continue;
        }		
        FillCutflow(&iCutFlow, HistCutFlow); // cutflow: 1

        // extract T0F-matched tracks with DCA < 3.0 cm
		vector <StUPCTrack const*> tracksWithTofHit;
        vector <StUPCTrack const*> tracksWithSmallDca;
        for (Int_t i = 0; i<upcEvt->getNumberOfTracks(); i++)
		{
            if (sqrt( pow(upcEvt->getTrack(i)->getDcaXY(),2) + pow(upcEvt->getTrack(i)->getDcaZ(),2)) < 3.0)
            {
                tracksWithSmallDca.push_back(upcEvt->getTrack(i));
                if(upcEvt->getTrack(i)->getFlag(StUPCTrack::kTof));
                {
                    tracksWithTofHit.push_back(upcEvt->getTrack(i));
                } 
            }
		}       		

        // at least 3 TOF-matched tracks with DCA < 3.0 cm
		if (tracksWithTofHit.size() < 3)
		{
			continue;
		}
        FillCutflow(&iCutFlow, HistCutFlow); // cutflow: 2


        // loop over possible pairs
        StUPCTrack const * track1;
        StUPCTrack const * track2;
        for (int i = 0; i<tracksWithSmallDca.size(); i++)
		{
            for (int j = 1; j<tracksWithSmallDca.size(); j++)
		    {
                if (i != j)
                {
                    track1 = tracksWithSmallDca[i];
                    track2 = tracksWithSmallDca[j];
                    bool hasTofHitTrack1 = false;
                    bool hasTofHitTrack2 = false;

                    if (track1->getFlag(StUPCTrack::kTof))
                    {
                        hasTofHitTrack1 = true;
                    }

                    if (track2->getFlag(StUPCTrack::kTof))
                    {
                        hasTofHitTrack2 = true;
                    }

                    if (hasTofHitTrack1 == false and hasTofHitTrack2 == false)
                    {
                        continue;
                    }

                    StUPCV0 kaon(track1, track2, massPion, massPion, 1, 1, tryVec, upcEvt->getMagneticField(), bool(isMC), true);

                    if (kaon.dcaDaughters() < 1)
                    {   
                        if (int(hasTofHitTrack1) + int(hasTofHitTrack2) == 2)
                        {  
                            HistKaonMassProbeWithTof->Fill(kaon.m());
                        }
                        
                        else
                        {
                            HistKaonMassProbeWithoutTof->Fill(kaon.m());
                        }
                    }
                }
            }
        }
    }
    
    
    TFile *outfile = TFile::Open(argv[2], "recreate"); 

    HistCutFlow->Write();
    HistKaonMassProbeWithTof->Write();
    HistKaonMassProbeWithoutTof->Write();

    outfile->Close();

    return 0;
}

void FillCutflow(int *i, TH1D *hist )
{

    hist->Fill(*i);
    *i+=1;
}

void FindProtons(bool isMC, StRPEvent *rpEvt, StUPCEvent *upcEvt, TVector3 & proton1, TVector3 & proton2)
{
    if (isMC == 0)
    {
        for(unsigned int k = 0; k < rpEvt->getNumberOfTracks(); ++k)
        {		
            StUPCRpsTrack *trk = rpEvt->getTrack(k);
            trk->setEvent(rpEvt);
            if (k == 0)
            {
                proton1 = trk->pVec();
            }
            else if (k == 1)
            {
                proton2 = trk->pVec();
            }     
        }    
    }

    else if (isMC == 1)
    {
        vector <TParticle*> vProtons;
        TParticle* particle;
        
        for (int i = 0; i < upcEvt->getNumberOfMCParticles(); i++)
        {
            particle = upcEvt->getMCParticle(i);
            if (particle->GetPDG()->PdgCode() == 2212 and particle->GetFirstMother() == 1)
            {
                vProtons.push_back(particle);
            }
        }

        double sigma = 0.033;
        double px_reco1 = gRandom->Gaus(vProtons[0]->Px(), sigma);
        double py_reco1 = gRandom->Gaus(vProtons[0]->Py(), sigma);
        proton1.SetX(px_reco1);
        proton1.SetY(py_reco1);
        proton1.SetZ(vProtons[0]->Pz());

        double px_reco2 = gRandom->Gaus(vProtons[1]->Px(), sigma);
        double py_reco2 = gRandom->Gaus(vProtons[1]->Py(), sigma);
        proton2.SetX(px_reco2);
        proton2.SetY(py_reco2);
        proton2.SetZ(vProtons[1]->Pz());
        vProtons.clear();
    }
}

