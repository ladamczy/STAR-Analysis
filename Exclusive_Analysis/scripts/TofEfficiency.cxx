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
#include "BeamPosition.h"
#include <Afterburner.h>
#include <StRPEvent.h>
#include <map>

using namespace std;
void FindProtons(bool isMC, StRPEvent *rpEvt, StUPCEvent *upcEvt, TLorentzVector & proton1, TLorentzVector & proton2);

int main(int argc, char** argv)  
{
    float massPion = 0.13957061;
    double massKaon =  497.611/1000.0;
    double beamEnergy = 254.867;
    double massProton = 938.272/1000.0;

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
    static StUPCEvent *upcEvt = 0x0;
    static StRPEvent  *correctedRpEvent = 0x0;
    LoadOffsetFile("../share/OffSetsCorrectionsRun17.list", mCorrection);
    chain->SetBranchAddress("mUPCEvent", &upcEvt);

    chain->GetEntry(0);

    if (upcEvt->getRunNumber() == 1)
    {
        isMC = 1;
    }

    if (isMC == 0)
    {
        chain->SetBranchAddress("correctedRpEvent", &correctedRpEvent);
    }

    vector <int> triggerID = {570209, 570219, 570229, 570701, 570702, 570703, 570704, 570705, 570709, 570711, 570712, 570719, 590701, 590703, 590705, 590708, 590709};
    vector <int> triggerCEP = {570701, 570705, 570711};
    vector <int> fillMC = {20667,20704,20705,20668,20706,20707,20669,20708,20719,20720,20728,20733,20734,20735,20736,20670,20737,20738,20740,20741,20743,20745,20746,20747,20748,20751,20752,20756,20757,20758,20759,20760,20764,20765,20766,20768,20770,20771,20772,20773,20777,20671,20778,20780,20781,20783,20787,20788,20789,20790,20791,20799,20674,20800,20810,20811,20814,20815,20821,20823,20825,20827,20828,20835,20841,20845,20852,20675,20854,20856,20857,20859,20860,20861,20862,20863,20864,20865,20866,20867,20868,20870,20873,20874,20896,20897,20898,20900,20901,20902,20903,20680,20904,20906,20907,20910,20911,20920,20681,20924,20925,20926,20927,20928,20929,20682,20684,20687,20688,20692,20695,20697,20701,20702,20703};
    vector <int> fillRange;

    sort(fillMC.begin(), fillMC.end());
    for (int i = 0; i < fillMC.size(); i++)
    {
        if (i%23 == 0)
        {
          fillRange.push_back(i);
        }    
    }

    for (int i = 0; i < fillMC.size(); i++)
    {
        cout << fillMC[i] << ", " ;
    }
    cout << " \n==================\n";

    fillRange.push_back(fillMC.size()-1);
    for (int i = 0; i < fillRange.size(); i++)
    {
        cout << fillMC[fillRange[i]] << ", " ;
    }

    double pTmin = 0.0; double pTmax = 4.0; int nBinsPt =7;
    double etaMin = -3.0; double etaMax = 3.0; int nBinsEta = 7;
    double zMin = -80.0; double zMax = 80.0; int nBinsZ = 17;
    int nBinsFill = 7;
    
    double binsPt[7] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.2};
    double binsEta[7] = {-0.9, -0.6, -0.3, 0.0, 0.3, 0.6, 0.9};

    TH1D* HistKaonMassProbeWithoutTof = new TH1D("HistKaonMassProbeWithoutTof", "; m_{#pi^{+}#pi^{-}}^{tag} [GeV]; # events",60 ,0.44, 0.56);  
    TH1D* HistKaonMassProbeWithTof = new TH1D("HistKaonMassProbeWithTof", "; m_{#pi^{+}#pi^{-}}^{probe} [GeV]; # events", 60 ,0.44, 0.56);

    TH1D* HistKaonMassProbeWithoutTofPosProbe = new TH1D("HistKaonMassProbeWithoutTofPosProbe", "; m_{#pi^{+}#pi^{-}} [GeV]; events", 60, 0.44, 0.56);
    TH1D* HistKaonMassProbeWithoutTofNegProbe = new TH1D("HistKaonMassProbeWithoutTofNegProbe", "; m_{#pi^{+}#pi^{-}} [GeV]; events", 60, 0.44, 0.56);
    TH1D* HistKaonMassProbeWithTofPosProbe = new TH1D("HistKaonMassProbeWithTofPosProbe", "; m_{#pi^{+}#pi^{-}} [GeV]; events", 60, 0.44, 0.56);
    TH1D* HistKaonMassProbeWithTofNegProbe = new TH1D("HistKaonMassProbeWithTofNegProbe", "; m_{#pi^{+}#pi^{-}} [GeV]; events", 60, 0.44, 0.56);
 
    TH1D* HistKaonPtProbeWithoutTofPosProbe = new TH1D("HistKaonPtProbeWithoutTofPosProbe", "; p_{T} [GeV]; events",nBinsPt-1, binsPt);
    TH1D* HistKaonPtProbeWithoutTofNegProbe = new TH1D("HistKaonPtProbeWithoutTofNegProbe", "; P_{T} [GeV]; events", nBinsPt-1, binsPt);
    TH1D* HistKaonPtProbeWithTofPosProbe = new TH1D("HistKaonPtProbeWithTofPosProbe", "; p_{T} [GeV]; events", nBinsPt-1, binsPt);
    TH1D* HistKaonPtProbeWithTofNegProbe = new TH1D("HistKaonPtProbeWithTofNegProbe", "; p_{T} [GeV]; events", nBinsPt-1, binsPt);
 
    TH1D* HistKaonEtaProbeWithoutTofPosProbe = new TH1D("HistKaonEtaProbeWithoutTofPosProbe", "positive probe w/o TOF; #eta ; events",nBinsEta-1, binsEta);
    TH1D* HistKaonEtaProbeWithoutTofNegProbe = new TH1D("HistKaonEtaProbeWithoutTofNegProbe", "negative probe w/o TOF; #eta ; events",nBinsEta-1, binsEta);
    TH1D* HistKaonEtaProbeWithTofPosProbe = new TH1D("HistKaonEtaProbeWithTofPosProbe", "positive probe w/ TOF; #eta ; events",nBinsEta-1, binsEta);
    TH1D* HistKaonEtaProbeWithTofNegProbe = new TH1D("HistKaonEtaProbeWithTofNegProbe", "negative probe w/ TOF; #eta ; events",nBinsEta-1, binsEta);

    vector <TH1D*> HistKaonPtProbeWithoutTofPosProbeMass;
    vector <TH1D*> HistKaonPtProbeWithoutTofNegProbeMass;
    vector <TH1D*> HistKaonPtProbeWithTofPosProbeMass;
    vector <TH1D*> HistKaonPtProbeWithTofNegProbeMass;

    vector <TH1D*> HistKaonEtaProbeWithoutTofPosProbeMass;
    vector <TH1D*> HistKaonEtaProbeWithoutTofNegProbeMass;
    vector <TH1D*> HistKaonEtaProbeWithTofPosProbeMass;
    vector <TH1D*> HistKaonEtaProbeWithTofNegProbeMass;

    vector <TH1D*> HistKaonFillProbeWithoutTofPosProbeMass;
    vector <TH1D*> HistKaonFillProbeWithoutTofNegProbeMass;
    vector <TH1D*> HistKaonFillProbeWithTofPosProbeMass;
    vector <TH1D*> HistKaonFillProbeWithTofNegProbeMass;

    TH1D* HistKaonPtProbeWithoutTofPosProbeMassMC = new TH1D("HistKaonPtProbeWithoutTofPosProbeMassMC", "", 60, 0.44, 0.56);
    TH1D* HistKaonPtProbeWithoutTofNegProbeMassMC = new TH1D("HistKaonPtProbeWithoutTofNegProbeMassMC","",  60, 0.44, 0.56);
    TH1D* HistKaonPtProbeWithTofPosProbeMassMC = new TH1D("HistKaonPtProbeWithTofPosProbeMassMC","", 60, 0.44, 0.56);
    TH1D* HistKaonPtProbeWithTofNegProbeMassMC = new TH1D("HistKaonPtProbeWithTofNegProbeMassMC", "",60, 0.44, 0.56);

    TH1D* HistKaonEtaProbeWithoutTofPosProbeMassMC = new TH1D("HistKaonEtaProbeWithoutTofPosProbeMassMC", "", 60, 0.44, 0.56);
    TH1D* HistKaonEtaProbeWithoutTofNegProbeMassMC = new TH1D("HistKaonEtaProbeWithoutTofNegProbeMassMC","",  60, 0.44, 0.56);
    TH1D* HistKaonEtaProbeWithTofPosProbeMassMC = new TH1D("HistKaonEtaProbeWithTofPosProbeMassMC","", 60, 0.44, 0.56);
    TH1D* HistKaonEtaProbeWithTofNegProbeMassMC = new TH1D("HistKaonEtaProbeWithTofNegProbeMassMC", "",60, 0.44, 0.56);

    TH1D* HistPionPtWithTofMC = new TH1D("HistPionPtWithTofMC", "; p_{T} [GeV]; events",nBinsPt-1, binsPt);
    TH1D* HistPionPtWithoutTofMC = new TH1D("HistPionPtWithoutTofMC", "; P_{T} [GeV]; events", nBinsPt-1, binsPt);

    TH1D* HistPionEtaWithTofMC = new TH1D("HistPionEtaWithTofMC", "; #eta ; events",nBinsEta-1, binsEta);
    TH1D* HistPionEtaWithoutTofMC = new TH1D("HistPionEtaWithoutTofMC", "; #eta ; events",nBinsEta-1, binsEta);

    for (int i = 0; i < nBinsPt-1; ++i) 
    {
        TH1D* histPosWo = new TH1D(Form("HistKaonPtProbeWithoutTofPosProbeMass%d", i), Form("HistKaonPtProbeWithoutTofPosProbeMass %d", i), 60, 0.44, 0.56);
        TH1D* histNegWo = new TH1D(Form("HistKaonPtProbeWithoutTofNegProbeMass%d", i), Form("HistKaonPtProbeWithoutTofNegProbeMass %d", i), 60, 0.44, 0.56);
        TH1D* histPosW = new TH1D(Form("HistKaonPtProbeWithTofPosProbeMass%d", i), Form("HistKaonPtProbeWithTofPosProbeMass %d", i), 60, 0.44, 0.56);
        TH1D* histNegW = new TH1D(Form("HistKaonPtProbeWithTofNegProbeMass%d", i), Form("HistKaonPtProbeWithTofNegProbeMass %d", i), 60, 0.44, 0.56);
        HistKaonPtProbeWithoutTofPosProbeMass.push_back(histPosWo);
        HistKaonPtProbeWithoutTofNegProbeMass.push_back(histNegWo);
        HistKaonPtProbeWithTofPosProbeMass.push_back(histPosW);
        HistKaonPtProbeWithTofNegProbeMass.push_back(histNegW);
    }

    for (int i = 0; i < nBinsEta-1; ++i) 
    {
        TH1D* histPosWo = new TH1D(Form("HistKaonEtaProbeWithoutTofPosProbeMass%d", i), Form("HistKaonEtaProbeWithoutTofPosProbeMass %d", i), 60, 0.44, 0.56);
        TH1D* histNegWo = new TH1D(Form("HistKaonEtaProbeWithoutTofNegProbeMass%d", i), Form("HistKaonEtaProbeWithoutTofNegProbeMass %d", i), 60, 0.44, 0.56);
        TH1D* histPosW = new TH1D(Form("HistKaonEtaProbeWithTofPosProbeMass%d", i), Form("HistKaonEtaProbeWithTofPosProbeMass %d", i), 60, 0.44, 0.56);
        TH1D* histNegW = new TH1D(Form("HistKaonEtaProbeWithTofNegProbeMass%d", i), Form("HistKaonEtaProbeWithTofNegProbeMass %d", i), 60, 0.44, 0.56);
        HistKaonEtaProbeWithoutTofPosProbeMass.push_back(histPosWo);
        HistKaonEtaProbeWithoutTofNegProbeMass.push_back(histNegWo);
        HistKaonEtaProbeWithTofPosProbeMass.push_back(histPosW);
        HistKaonEtaProbeWithTofNegProbeMass.push_back(histNegW);
    }
    
    for (int i = 0; i < nBinsFill; ++i) 
    {
        TH1D* histPosWo = new TH1D(Form("HistKaonFillProbeWithoutTofPosProbeMass%d", i), Form("HistKaonFillProbeWithoutTofPosProbeMass %d", i), 60, 0.44, 0.56);
        TH1D* histNegWo = new TH1D(Form("HistKaonFillProbeWithoutTofNegProbeMass%d", i), Form("HistKaonFillProbeWithoutTofNegProbeMass %d", i), 60, 0.44, 0.56);
        TH1D* histPosW = new TH1D(Form("HistKaonFillProbeWithTofPosProbeMass%d", i), Form("HistKaonFillProbeWithTofPosProbeMass %d", i), 60, 0.44, 0.56);
        TH1D* histNegW = new TH1D(Form("HistKaonFillProbeWithTofNegProbeMass%d", i), Form("HistKaonFillProbeWithTofNegProbeMass %d", i), 60, 0.44, 0.56);
        HistKaonFillProbeWithoutTofPosProbeMass.push_back(histPosWo);
        HistKaonFillProbeWithoutTofNegProbeMass.push_back(histNegWo);
        HistKaonFillProbeWithTofPosProbeMass.push_back(histPosW);
        HistKaonFillProbeWithTofNegProbeMass.push_back(histNegW);
    }

    int fill;
    int nextFill = 0;
    int indxFill[7] = {0,1,2,3,4,5,6};

    for (Long64_t i = 0; i < chain->GetEntries(); ++i) 
    {    
        if (i%1000000 == 0)  
        {
            cout << i << "/" <<  chain->GetEntries() << endl;
        }

        chain->GetEntry(i);
        TLorentzVector lorentzVector,trackVector;
        TLorentzVector proton1, proton2;

        //grupowanie numerów fillów 
        fill = upcEvt->getFillNumber();
        int nextFill;

        if (fill < fillMC[fillRange[0]])
        {
            nextFill = 0;
        }
 
        for (int i = 0; i < fillRange.size()-1; i++)
        {
            if (i != fillRange.size()-2)
            {
            if (fill >= fillMC[fillRange[i]] and fill < fillMC[fillRange[i+1]])
            {
                nextFill = indxFill[i+1];

            }
            }

            else if (i == fillRange.size()-2)
            {
            if (fill >= fillMC[fillRange[i]] and fill <= fillMC[fillRange[i+1]])
            {
                nextFill = indxFill[i+1];

            }
            }
        }

        // vectors for tracks
        vector <StUPCTrack const*> tracksWithTofHit;
        vector <StUPCTrack const*> tracksWithSmallDca;
        vector <StUPCTrack const*> tracksForTrueLevel;
        // 1 step: SELECT tracks with kV0 flag, absEta < 0.9 and pT > 0.15 GeV and Nhit >= 20
        for (Int_t i = 0; i<upcEvt->getNumberOfTracks(); i++)
		{
            if (isMC == 0)
            {
                if ( (abs(upcEvt->getTrack(i)->getEta()) <= 0.9) and (upcEvt->getTrack(i)->getPt() >= 0.2) and (upcEvt->getTrack(i)->getNhitsFit() >= 20) and upcEvt->getTrack(i)->getFlag(StUPCTrack::kV0) )
                {
                    tracksWithSmallDca.push_back(upcEvt->getTrack(i));
                    if(upcEvt->getTrack(i)->getFlag(StUPCTrack::kTof));
                    {
                        tracksWithTofHit.push_back(upcEvt->getTrack(i));
                    } 
                }
            }
            else if (isMC == 1)
            {
                        tracksForTrueLevel.push_back(upcEvt->getTrack(i));
                 
                if ( (abs(upcEvt->getTrack(i)->getEta()) <= 0.9) and (upcEvt->getTrack(i)->getPt() >= 0.2) and (upcEvt->getTrack(i)->getNhitsFit() >= 20)) 
                {
                   tracksWithSmallDca.push_back(upcEvt->getTrack(i));
            
                    if(upcEvt->getTrack(i)->getFlag(StUPCTrack::kTof));
                    {
                        tracksWithTofHit.push_back(upcEvt->getTrack(i));
                    } 
                }
            }
		} 

        // analiza poziomu true - tylko dla MC = niezależna od tag and probe
        if (isMC==1) 
        {
            TParticle* particle;
            vector <TParticle *> vPions; // bierzemy pod uwagę wszystkie możliwe piony (niekoniecznie z rozpadu K0)
            for (int i = 0; i <  upcEvt->getNumberOfMCParticles(); i++)
            {
                particle = upcEvt->getMCParticle(i);
                if (particle->GetPDG()->PdgCode() == 211 or particle->GetPDG()->PdgCode() == -211)
                {
                    vPions.push_back(particle);
                }
            }

            TLorentzVector productionVertex; 
            TLorentzVector lorentzTrack;        
         
            for (int j = 0; j < vPions.size(); j++) // pętla po wszyskich pionach z poziomu true - szukanie matchingu
            {
                lorentzVector.SetPxPyPzE(vPions[j]->Px(), vPions[j]->Py(), vPions[j]->Pz(), vPions[j]->Energy()); 
                double pT = sqrt(pow(vPions[j]->Px(),2) + pow(vPions[j]->Py(),2)) ;
                double eta = lorentzVector.Eta();
                vPions[j]->ProductionVertex(productionVertex);

                double truthVertexR = sqrt(pow(productionVertex.X(),2) + pow(productionVertex.Y(),2)); //  poprzeczne położenie verteksu
                double truthVertexZ = productionVertex.Z(); // podłużne położenie werteksu
            
                if (abs(truthVertexZ) > 80 or abs(truthVertexR) > 3.0 )//or pT < 0.2 or abs(eta)> 0.9)// sprawdzenie czy true vertex i pt pionu jest w odpowiednim zakresie 
                {
                   continue; //ominięcie pionu poza zakresem vtxR,  vtxZ, pT true
               }
              
                vector <double> distr;
                vector <StUPCTrack const *> matchedTPCpions;
                if (tracksForTrueLevel.size() == 0) {break;}
                vector <int> indices;
                for (int i = 0; i < tracksForTrueLevel.size(); i++) // pętla po śladach globalnych dobrej jakości
                {
                    TLorentzVector pion, track;
                    pion.SetPxPyPzE(vPions[j]->Px(), vPions[j]->Py(), vPions[j]->Pz(), vPions[j]->Energy());
                    tracksForTrueLevel[i]->getLorentzVector(track, massPion);
                    double r = pion.DeltaR(track); // odległość w eta phi
                    matchedTPCpions.push_back(tracksForTrueLevel[i]); // tracki z tracksForTrueLevel są usuwane poniżej po znalezieniu najmniejszego r - dlatego zbieram jeszcze raz ślady do wektora
                    distr.push_back(r); // dla dnaego pionu z piomu true obliczenie wszsystkich możliwych odległości R dla tracków dobrej jakości tracksForTrueLevel
                    indices.push_back(i);
                }
          
                auto it = std::min_element(distr.begin(), distr.end()); //znajdź element o najmniejszym R
                int index = std::distance(distr.begin(), it); // odpowiedni indeks
                int inMinToRemove = indices[index];
                if (distr[index] < 0.15) // jeżeli odległość mniejsza niż .15 - akceptujemy przypadek
                {  
                    if (matchedTPCpions[index]->getPt()>=0.2 and abs(matchedTPCpions[index]->getEta()) <= 0.9 and matchedTPCpions[index]->getNhitsFit() >=20 )
                    {
                        if (matchedTPCpions[index]->getFlag(StUPCTrack::kTof)) // podział TOF nie TOF
                        {
                            HistPionPtWithTofMC->Fill(matchedTPCpions[index]->getPt());
                            HistPionEtaWithTofMC->Fill(matchedTPCpions[index]->getEta());
                        }

                        else
                        {
                            HistPionPtWithoutTofMC->Fill(matchedTPCpions[index]->getPt());
                            HistPionEtaWithoutTofMC->Fill(matchedTPCpions[index]->getEta());
                        } 
                        tracksForTrueLevel.erase(tracksForTrueLevel.begin()+inMinToRemove); // usunięcie z tracksForTrueLevel śladu
                    }
                }    
            }
        }
		if (tracksWithTofHit.size() < 3) // rozpatrujemy przypadki gdzie są przynajmniej trzy ślady zmatchowane z TOF
		{
			continue;
		}

   //             Int_t BBCLEast = upcEvt->getBBCLargeEast();
    //    Int_t BBCLWest = upcEvt->getBBCLargeWest();
//
  //      if (BBCLWest > 52 or BBCLEast > 52.0)
    //    {
      //      continue;
       // }

        StUPCTrack const * track1;
        StUPCTrack const * track2;
        for (int i = 0; i<tracksWithSmallDca.size(); i++) // pęta po trackach dobrej jakości
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

                    if (track1->getCharge() + track2->getCharge() !=0)
                    {
                        continue;
                    }

                    int nFillNumber = upcEvt->getFillNumber();      
                    vector <double> beamParVec;
                    double beamPar[4] = {};

                    if (isMC == 0)
                    {    
                        beamPar[0] = upcEvt->getBeamXPosition();
                        beamPar[1] = upcEvt->getBeamXSlope();
                        beamPar[2] = upcEvt->getBeamYPosition();
                        beamPar[3] = upcEvt->getBeamYSlope();
                    }

                    else
                    {
                        beamPar[0] = 0.0;
                        beamPar[1] = 0.0;
                        beamPar[2] = 0.0;
                        beamPar[3] = 0.0;
                    }

                    TVector3 const tryVec(0,0,0);
                    StUPCV0 kaon(track1, track2, massPion, massPion, 1, 1, tryVec,beamPar, upcEvt->getMagneticField(), true);
                    double dcaDau =  kaon.dcaDaughters();
                    double dcaBeam =  kaon.DCABeamLine();
                    double hypoAngle = kaon.pointingAngleHypo();
                    double d = kaon.decayLengthHypo();
                    double vtxZ = kaon.decayVertex().Z();
                    double vtxR = sqrt( pow(kaon.decayVertex().X() , 2)  +pow(kaon.decayVertex().Y() , 2)   );

                    if (  abs(vtxZ) < 80 and vtxR < 3)
                    {   
                     
                        if (int(hasTofHitTrack1) + int(hasTofHitTrack2) == 2)
                        {  
                            if (dcaDau < 1 and dcaBeam < 1 )
                            {

                            HistKaonMassProbeWithTof->Fill(kaon.m());
                            HistKaonMassProbeWithTof->Fill(kaon.m());

                            StUPCTrack const * tagWTof1 = track1;
                            StUPCTrack const * probeWTof1 = track2;  

                            StUPCTrack const * tagWTof2 = track2;
                            StUPCTrack const * probeWTof2 = track1;  

                            
                            for (int i = 0; i <nBinsPt-1; i++)
                            {
                                if ( probeWTof1->getPt() >= binsPt[i] and probeWTof1->getPt() < binsPt[i+1]  )
                               {
                                    HistKaonPtProbeWithTofPosProbeMass[i]->Fill(kaon.m());
                                }
                            }

                            for (int i = 0; i <nBinsEta-1; i++)
                            {
                                if ( probeWTof1->getEta() >= binsEta[i] and probeWTof1->getEta() < binsEta[i+1]  )
                                {
                                    HistKaonEtaProbeWithTofPosProbeMass[i]->Fill(kaon.m());
                                }
                            }
            
                            for (int i = 0; i <nBinsPt-1; i++)
                            {
                                if ( probeWTof2->getPt() >= binsPt[i] and probeWTof2->getPt() < binsPt[i+1]  )
                                {
                                    HistKaonPtProbeWithTofPosProbeMass[i]->Fill(kaon.m());
                                }
                            }

                            for (int i = 0; i <nBinsEta-1; i++)
                            {
                                if ( probeWTof2->getEta() >= binsEta[i] and probeWTof2->getEta() < binsEta[i+1]  )
                                {
                                    HistKaonEtaProbeWithTofPosProbeMass[i]->Fill(kaon.m());
                                }
                            }
                                HistKaonFillProbeWithTofPosProbeMass[nextFill]->Fill(kaon.m());   
                            }        
                        }
                        
                        else
                        {
                            if (dcaDau < 1  and dcaBeam <1)
                            {
                            HistKaonMassProbeWithoutTof->Fill(kaon.m());
                            StUPCTrack const * tagWoTof;
                            StUPCTrack const * probeWoTof;
 
                            if (hasTofHitTrack1)
                            {
                                tagWoTof = track1;
                                probeWoTof = track2;
                            }
                            else
                            {
                                tagWoTof = track2;
                                probeWoTof = track1;
                            }

                            for (int i = 0; i <nBinsPt-1; i++)
                            {
                                if ( probeWoTof->getPt() >= binsPt[i] and probeWoTof->getPt() < binsPt[i+1]  )
                                {
                                    HistKaonPtProbeWithoutTofPosProbeMass[i]->Fill(kaon.m());
                                }
                            }

                            for (int i = 0; i <nBinsEta-1; i++)
                            {
                                if ( probeWoTof->getEta() >= binsEta[i] and probeWoTof->getEta() < binsEta[i+1]  )
                                {
                                    HistKaonEtaProbeWithoutTofPosProbeMass[i]->Fill(kaon.m());
                                }
                            }
            
                            HistKaonFillProbeWithoutTofPosProbeMass[nextFill]->Fill(kaon.m());       
                            }
                        }
                        
                    }
                }
            }
        }
    }
            
    TFile *outfile = TFile::Open(argv[2], "recreate"); 

    HistKaonMassProbeWithTof->Write();
    HistKaonMassProbeWithoutTof->Write();

    HistKaonMassProbeWithoutTofPosProbe->Write();
    HistKaonMassProbeWithoutTofNegProbe->Write();
    HistKaonMassProbeWithTofPosProbe->Write();
    HistKaonMassProbeWithTofNegProbe->Write();

    HistKaonPtProbeWithoutTofPosProbe->Write();
    HistKaonPtProbeWithoutTofNegProbe->Write();
    HistKaonPtProbeWithTofPosProbe->Write();
    HistKaonPtProbeWithTofNegProbe->Write();
    HistKaonEtaProbeWithoutTofPosProbe->Write();
    HistKaonEtaProbeWithoutTofNegProbe->Write();
    HistKaonEtaProbeWithTofPosProbe->Write();
    HistKaonEtaProbeWithTofNegProbe->Write();

    HistKaonPtProbeWithoutTofNegProbeMassMC->Write();
    HistKaonPtProbeWithoutTofPosProbeMassMC->Write();
    HistKaonPtProbeWithTofNegProbeMassMC->Write();
    HistKaonPtProbeWithTofPosProbeMassMC->Write();

    HistKaonEtaProbeWithoutTofNegProbeMassMC->Write();
    HistKaonEtaProbeWithoutTofPosProbeMassMC->Write();
    HistKaonEtaProbeWithTofNegProbeMassMC->Write();
    HistKaonEtaProbeWithTofPosProbeMassMC->Write();

    HistPionPtWithTofMC->Write();
    HistPionEtaWithTofMC->Write();
    HistPionPtWithoutTofMC->Write();
    HistPionEtaWithoutTofMC->Write();
        

    for (int i = 0; i < HistKaonPtProbeWithoutTofPosProbeMass.size(); i++)
    {
        HistKaonPtProbeWithoutTofPosProbeMass[i]->Write();
        HistKaonPtProbeWithoutTofNegProbeMass[i]->Write();
        HistKaonPtProbeWithTofPosProbeMass[i]->Write();
        HistKaonPtProbeWithTofNegProbeMass[i]->Write();
    }
    for (int i = 0; i < HistKaonEtaProbeWithoutTofPosProbeMass.size(); i++)
    {
        HistKaonEtaProbeWithoutTofPosProbeMass[i]->Write();
        HistKaonEtaProbeWithoutTofNegProbeMass[i]->Write();
        HistKaonEtaProbeWithTofPosProbeMass[i]->Write();
        HistKaonEtaProbeWithTofNegProbeMass[i]->Write();
    }

    for (int i = 0; i < HistKaonFillProbeWithoutTofPosProbeMass.size(); i++)
    {
        HistKaonFillProbeWithoutTofPosProbeMass[i]->Write();
        HistKaonFillProbeWithoutTofNegProbeMass[i]->Write();
        HistKaonFillProbeWithTofPosProbeMass[i]->Write();
        HistKaonFillProbeWithTofNegProbeMass[i]->Write();
    }

    outfile->Close();

    return 0;
}


