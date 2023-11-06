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


void FillProtons(double protonSignsY, TH1D*HistProtonsSameSide, TH1D*HistProtonsOppositeSide, int iCutFlow)
{
    if (protonSignsY >= 0.0)
    {
        HistProtonsSameSide->Fill(iCutFlow);
    }
    else 
    {
       HistProtonsOppositeSide->Fill(iCutFlow);      
    }
}


int main(int argc, char** argv)  
{


    TH1D* HistCutFlow = new TH1D("HistCutFlow", ";;count", 16, -0.5, 15.5);
    TProfile* ProfileDistVertexBeamX = new TProfile("ProfileDistVertexBeamX","",639, 20511.5, 21150.5);
    TProfile* ProfileDistVertexBeamY = new TProfile("ProfileDistVertexBeamY","",639, 20511.5, 21150.5);

    TProfile* ProfileDistVertexBeamXAsAFunctionOfZ = new TProfile("ProfileDistVertexBeamXAsAFunctionOfZ","",100, -200, 200);
    TProfile* ProfileDistVertexBeamYAsAFunctionOfZ  = new TProfile("ProfileDistVertexBeamYAsAFunctionOfZ","",100, -200, 200);

    TH1D* HistDistVtxXLeadingKaon = new TH1D("HistDistVtxXLeadingKaon", "; x_{vtx}^{leading} - x_{beam} [cm]; # events", 100, -2, 2);
    TH1D* HistDistVtxYLeadingKaon = new TH1D("HistDistVtxYLeadingKaon"," ; y_{vtx}^{leading} - y_{beam} [cm]; # events", 100, -2, 2);
    TH1D* HistDistVtxZLeadingKaon = new TH1D("HistDistVtxZLeadingKaon", "; z_{vtx}^{leading} - z_{beam} [xm]; # events", 100, -200, 200);
    TH1D* HistRLeadingKaon = new TH1D("HistRLeadingKaon", ";R_{K0^{leading}} [cm]; # events", 30, 0, 6.0);

    TH1D* HistDistVtxXSubLeadingKaon = new TH1D("HistDistVtxXSubLeadingKaon", "; x_{vtx}^{SubLeading} - x_{beam} [cm]; # events", 100, -2, 2);
    TH1D* HistDistVtxYSubLeadingKaon = new TH1D("HistDistVtxYSubLeadingKaon"," ; y_{vtx}^{SubLeading} - y_{beam} [cm]; # events", 100, -2, 2);
    TH1D* HistDistVtxZSubLeadingKaon = new TH1D("HistDistVtxZSubLeadingKaon", "; z_{vtx}^{SubLeading} - z_{beam} [xm]; # events", 100, -200, 200);
    TH1D* HistRSubLeadingKaon = new TH1D("HistRSubLeadingKaon", ";R_{K0^{SubLeading}} [cm]; # events", 30, 0, 6.0);

    TH1D* HistNumPrimaryVertices = new TH1D("5_HistNumPrimaryVertices", " ;Number of primary vertices; count", 11, -0.5, 10.5); 
    TH1D* HistPrimaryVertexAbsPosZ = new TH1D("6_HistPrimaryVertexAbsPosZ", ";primary vertex z [cm]; count", 100, 0, 200);
    TH1D* HistNumTofMatchedTracks = new TH1D("7_HistNumTofMatchedTracks", ";number of TOF matched tracks; count", 16, -0.5, 15.5);    
    TH1D* HistTofMatchedTracksCharge = new TH1D("8_HistTofMatchedTracksCharge", ";charge [e]; count", 3, -1.5, 1.5);
    TH1D* HistTofMatchedTracksAbsEta = new TH1D("9_HistTofMatchedTracksAbsEta", ";#eta; count", 50, 0, 2);
    TH1D* HistTofMatchedTracksPt = new TH1D("9b_HistTofMatchedTracksPt", ";#eta; count", 50, 0, 2);
    TH1D* HistTofMatchedTracksNfit = new TH1D("12_HistTofMatchedTracksNfit", ";N_{fit}; count", 61, -0.5, 60.5);
    TH1D* HistTofMatchedTracksNdEdx = new TH1D("13_HistTofMatchedTracksNdEdx", ";N_{dE/dx}; count", 61, -0.5, 60.5);

    TH1D* HistNumPrimaryVerticesPostSelection = new TH1D("5_HistNumPrimaryVerticesPs", " ;Number of primary vertices; count", 11, -0.5, 10.5); 
    TH1D* HistPrimaryVertexAbsPosZPostSelection = new TH1D("6_HistPrimaryVertexAbsPosZPs", ";primary vertex z [cm]; count", 100, 0, 200);
    TH1D* HistNumTofMatchedTracksPostSelection = new TH1D("7_HistNumTofMatchedTracksPs", ";number of TOF matched tracks; count", 16, -0.5, 15.5);    
    TH1D* HistTofMatchedTracksChargePostSelection = new TH1D("8_HistTofMatchedTracksChargePs", ";charge [e]; count", 3, -1.5, 1.5);
    TH1D* HistTofMatchedTracksAbsEtaPostSelection = new TH1D("9_HistTofMatchedTracksAbsEtaPs", ";#eta; count", 50, 0, 2);
    TH1D* HistTofMatchedTracksPtPostSelection = new TH1D("9b_HistTofMatchedTracksPtPostSelection", ";#eta; count", 50, 0, 2);
    TH1D* HistTofMatchedTracksNfitPostSelection = new TH1D("12_HistTofMatchedTracksNfitPs", ";N_{fit}; count", 61, -0.5, 60.5);
    TH1D* HistTofMatchedTracksNdEdxPostSelection = new TH1D("13_HistTofMatchedTracksNdEdxPs", ";N_{dE/dx}; count", 61, -0.5, 60.5);

    TH2D* HistPtProtonPtKaonsX = new TH2D("HistPtProtonPtKaonsX", " ;[p_{K^{0}}p_{K^{0}}]_{x}; [p_{p'}^{W}p_{p'}^{E}]_{x} ",  100, -2.0, 2.0, 100, -2.0, 2.0);
    TH2D* HistPtProtonPtKaonsY = new TH2D("HistPtProtonPtKaonsY", ";[p_{K^{0}}p_{K^{0}}]_{y}; [p_{p'}^{W}p_{p'}^{E}]_{y} ",  100, -2.0, 2.0, 100, -2.0, 2.0);
    TH1D* HistPtProtonPtKaonsX1d = new TH1D("HistPtProtonPtKaonsX1d", " ;[p_{K^{0}}p_{K^{0}} + p_{p'}^{W}p_{p'}^{E}]_{x}; events ",  200, -5.0, 5.0);
    TH1D* HistPtProtonPtKaonsY1d = new TH1D("HistPtProtonPtKaonsY1d", ";[p_{K^{0}}p_{K^{0}} +  p_{p'}^{W}p_{p'}^{E}]_{y}; events ",  200, -5.0, 5.0);
    TH2D* HistPtProtonPtKaonsXOppY = new TH2D("HistPtProtonPtKaonsXOppY", "the opposite py signs;[p_{K^{0}}p_{K^{0}}]_{x}; [p_{p'}^{W}p_{p'}^{E}]_{x} ",  100, -2.0, 2.0, 100, -2.0, 2.0);
    TH2D* HistPtProtonPtKaonsYOppY = new TH2D("HistPtProtonPtKaonsYOppY", "the opposite py signs;[p_{K^{0}}p_{K^{0}}]_{y}; [p_{p'}^{W}p_{p'}^{E}]_{y} ",  100, -2.0, 2.0, 100, -2.0, 2.0);
    TH1D* HistPtProtonPtKaonsX1dOppY = new TH1D("HistPtProtonPtKaonsX1dOppY", "the opposite py signs;[p_{K^{0}}p_{K^{0}} + p_{p'}^{W}p_{p'}^{E}]_{x}; events ",  25, -3.0, 3.0);
    TH1D* HistPtProtonPtKaonsY1dOppY = new TH1D("HistPtProtonPtKaonsY1dOppY", "the opposite py signs;[p_{K^{0}}p_{K^{0}} +  p_{p'}^{W}p_{p'}^{E}]_{y}; events ",  25, -3.0, 3.0);
    TH2D* HistPtProtonPtKaonsXSameY = new TH2D("HistPtProtonPtKaonsXSameY", "the same py signs;[p_{K^{0}}p_{K^{0}}]_{x}; [p_{p'}^{W}p_{p'}^{E}]_{x} ",  100, -2.0, 2.0, 100, -2.0, 2.0);
    TH2D* HistPtProtonPtKaonsYSameY = new TH2D("HistPtProtonPtKaonsYSameY", "the same py signs;[p_{K^{0}}p_{K^{0}}]_{y}; [p_{p'}^{W}p_{p'}^{E}]_{y} ",  100, -2.0, 2.0, 100, -2.0, 2.0);
    TH1D* HistPtProtonPtKaonsX1dSameY = new TH1D("HistPtProtonPtKaonsX1dSameY", "the same py signs;[p_{K^{0}}p_{K^{0}} + p_{p'}^{W}p_{p'}^{E}]_{x}; events ",  25, -3,3);
    TH1D* HistPtProtonPtKaonsY1dSameY = new TH1D("HistPtProtonPtKaonsY1dSameY", "the same py signs;[p_{K^{0}}p_{K^{0}} + p_{p'}^{W}p_{p'}^{E}]_{y}; events ",  25, -3,3);

    TH2D* HistPtProtonPtKaonsXSameYMomentaGreaterThan0 = new TH2D("HistPtProtonPtKaonsXSameYMomentaGreaterThan0", ";[p_{K^{0}}p_{K^{0}}]_{x}; [p_{p'}^{W}p_{p'}^{E}]_{x} ",  100, -2.0, 2.0, 100, -2.0, 2.0);
    TH2D* HistPtProtonPtKaonsYSameYMomentaGreaterThan0 = new TH2D("HistPtProtonPtKaonsYSameYMomentaGreaterThan0", ";[p_{K^{0}}p_{K^{0}}]_{y}; [p_{p'}^{W}p_{p'}^{E}]_{y} ",  100, -2.0, 2.0, 100, -2.0, 2.0);
    TH1D* HistPtProtonPtKaonsX1dSameYMomentaGreaterThan0 = new TH1D("HistPtProtonPtKaonsX1dSameYMomentaGreaterThan0", "the same py signs, py > 0;[p_{K^{0}}p_{K^{0}} + p_{p'}^{W}p_{p'}^{E}]_{x}; events ",  25, -3.0, 3.0);
    TH1D* HistPtProtonPtKaonsY1dSameYMomentaGreaterThan0 = new TH1D("HistPtProtonPtKaonsY1dSameYMomentaGreaterThan0", "the same py signs, py > 0;[p_{K^{0}}p_{K^{0}} + p_{p'}^{W}p_{p'}^{E}]_{y}; events ",  25, -3.0, 3.0);

    TH2D* HistPtProtonPtKaonsXSameYMomentaLowerThan0 = new TH2D("HistPtProtonPtKaonsXSameYMomentaLowerThan0", ";[p_{K^{0}}p_{K^{0}}]_{x}; [p_{p'}^{W}p_{p'}^{E}]_{x} ",  100, -2.0, 2.0, 100, -2.0, 2.0);
    TH2D* HistPtProtonPtKaonsYSameYMomentaLowerThan0 = new TH2D("HistPtProtonPtKaonsYSameYMomentaLowerThan0", ";[p_{K^{0}}p_{K^{0}}]_{y}; [p_{p'}^{W}p_{p'}^{E}]_{y} ",  100, -2.0, 2.0, 100, -2.0, 2.0);
    TH1D* HistPtProtonPtKaonsX1dSameYMomentaLowerThan0 = new TH1D("HistPtProtonPtKaonsX1dSameYMomentaLowerThan0", "the same py signs, py < 0;[p_{K^{0}}p_{K^{0}} + p_{p'}^{W}p_{p'}^{E}]_{x}; events ",  25, -3.0, 3.0);
    TH1D* HistPtProtonPtKaonsY1dSameYMomentaLowerThan0 = new TH1D("HistPtProtonPtKaonsY1dSameYMomentaLowerThan0", "the same py signs, py < 0;[p_{K^{0}}p_{K^{0}} + p_{p'}^{W}p_{p'}^{E}]_{y}; events ",  25, -3.0, 3.0);

    TH1D* HistInvMassPiPiPeak = new TH1D("HistInvMassPiPiPeak", "; m_{#pi^{+}#pi^{-}} [GeV]; # events", 25 ,0.44, 0.54);
    TH1D* HistInvMassPiPiPeak2 = new TH1D("HistInvMassPiPiPeak2", "; m_{#pi^{+}#pi^{-}} [GeV]; # events", 25 ,0.44, 0.54);
    TH2D* HistInvMassPiPi2D =  new TH2D("HistInvMassPiPi2D", "; m_{#pi^{+}#pi^{-}}^{leading} [GeV];  m_{#pi^{+}#pi^{-}}^{sub-leading} [GeV]", 25 ,0.44, 0.54,  25 ,0.44, 0.54);
    TH2D* HistInvMassPiPi2DTofPtMiss =  new TH2D("HistInvMassPiPi2DTofPtMiss", "; m_{#pi^{+}#pi^{-}}^{leading} [GeV];  m_{#pi^{+}#pi^{-}}^{sub-leading} [GeV]", 25 , 0.44, 0.54,  25 ,0.44, 0.54);
    TH1D* HistNumOfClusters = new TH1D("HistNumOfClusters", "; # TOF clusters; # events", 40, -0.5, 40.5);
    TH1D* HistNumOfClustersPostSelection = new TH1D("HistNumOfClustersPs", "; # TOF clusters; # events", 40, -0.5, 40.5);

    TH1D* HistPtMiss = new TH1D("HistPtMiss", ";pT^{miss} [GeV]; # events", 60, 0.0, 3.0);
    TH1D* HistPtMissSameY = new TH1D("HistPtMissSameY", ";pT^{miss} [GeV]; # events", 100, 0.0, 3.0);
    TH1D* HistPtMissOppY = new TH1D("HistPtMissOppY", ";pT^{miss} [GeV]; # events", 100, 0.0, 3.0);
    TH1D* HistPtMissPostSelection = new TH1D("HistPtMissPs", ";pT^{miss} [GeV]; # events", 100, 0.0, 3.0);

    TH1D* HistKaonPtTruth = new TH1D("HistKaonPtTruth", "; pT_{K^{0}} [GeV]; # events", 25 ,0, 2.5);
    TH1D* HistKaonEtaTruth = new TH1D("HistKaonEtaTruth", "; #eta_{K^{0}}; # events", 25 ,-4.5, 4.5);
    TH1D* HistKaonVtxZTruth = new TH1D("HistKaonVtxZTruth", ";z_{vtx}^{#pi^{+}#pi^{-}} [cm]; # events", 10, -150, 150);
    TH1D* HistKaonVtxRTruth = new TH1D("HistKaonVtxRTruth", ";R_{vtx}^{#pi^{+}#pi^{-}} [cm]; # events", 10, 0, 20);
    TH1D* HistKaonPtDet = new TH1D("HistKaonPtDet", "; pT_{K^{0}} [GeV]; # events", 25 ,0, 2.5);
    TH1D* HistKaonEtaDet = new TH1D("HistKaonEtaDet", "; #eta_{K^{0}}; # events", 25 ,-4.5, 4.5);
    TH1D* HistKaonVtxZDet = new TH1D("HistKaonVtxZDet", "; z_{vtx}^{#pi^{+}#pi^{-}} [cm]; # events", 10, -150, 150);
    TH1D* HistKaonVtxRDet = new TH1D("HistKaonVtxRDet", "; R_{vtx}^{#pi^{+}#pi^{-}} [cm]; # events", 10, 0, 20);
    TH1D* HistProtonsSameSide = new TH1D("HistProtonsSameSide", " ; cutflow ;  # events where protons have the same sign of py",  16, -0.5, 15.5);
    TH1D* HistProtonsOppositeSide = new TH1D("HistProtonsOppositeSide", "; ; # events where protons have the opposite sign of py",  16, -0.5, 15.5);
    TH1D* HistDcaDaughtersLeadingKaon = new TH1D("HistDcaDaughtersLeadingKaon", "; DCA_{K0^{leading}} [cm] # events", 50, 0, 6);
    TH1D* HistDcaDaughtersSubLeadingKaon = new TH1D("HistDcaDaughtersSubLeadingKaon", "; DCA_{K0^subleading}} [cm]; #events", 50, 0, 6);

    TH1D* HistMassK0K0 = new TH1D("HistMassK0K0", " ; m_{K^{0}K^{0}} [GeV] ; events", 42, 0.9, 3);

    TH1D* HistPtK0K0 = new TH1D("HistPtK0K0", " ; p_{T, K^{0}K^{0}} [GeV] ; events", 20, 0.0, 2);
    TH1D* HistEtaK0K0 = new TH1D("HistEtaK0K0", " ; #eta_{K^{0}K^{0}} [GeV] ; events", 21, -2, 2);    
    TH1D* HistPtLeadingK0 = new TH1D("HistPtLeadingK0", " ; p_{T, K^{0}}^{leading} [GeV] ; events", 20, 0.0, 2);
    TH1D* HistEtaLeadingK0 = new TH1D("HistEtaLeadingK0", " ; #eta+{K^{0}}^{leading} [GeV] ; events", 50, -3, 3);
    TH1D* HistPtSubLeadingK0 = new TH1D("HistPtSubLeadingK0", " ; p_{T, K^{0}}^{leading} [GeV] ; events", 20, 0.0, 2);
    TH1D* HistEtaSubLeadingK0 = new TH1D("HistEtaSubLeadingK0", " ; #eta_{K^{0}}^{leading} [GeV] ; events", 50, -3, 3);
    TH1D* HistEtaLeadingSubLeadingDiff = new TH1D("HistEtaLeadingSubLeadingDiff", " ; #eta_{K^{0}}^{leading} [GeV] ; events", 50, -3, 3);

    TH1D * HistClosestDistanceBetweenKaons = new TH1D("HistClosestDistanceBetweenKaons", "; dist_{K^{0}K^{0}} ; events", 200, 0, 4);
    TH1D* HistZK0K0 = new TH1D("HistZK0K0", " ;vtx^{z}_{K0K0} [cm] ; events", 200, -200, 200);

    float massPion = 0.13957061;
    double massKaon =  497.611/1000.0;

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

    for (Long64_t i = 0; i < chain->GetEntries(); ++i) 
    {

        vector <TParticle*> vPosPionsMC, vNegPionsMC, vPosPionsMCPaired, vNegPionsMCPaired;
        TLorentzVector trackVector;
        TVector3 proton1, proton2;
        TParticle* particle;
        
        if (i%1000000 == 0)  
        {
            cout << i << "/" <<  chain->GetEntries() << endl;
        }

        chain->GetEntry(i);

        iCutFlow = 0;
        HistCutFlow->Fill(iCutFlow);  // cutflow: all data

   
        double protonSignsY = 1;
        vector <TParticle*> vProtons;
        if (isMC == 0)
        {
            for(unsigned int k = 0; k < rpEvt->getNumberOfTracks(); ++k)
            {		
                StUPCRpsTrack *trk = rpEvt->getTrack(k);
                trk->setEvent(rpEvt);
                if (k == 0)
                {
                    proton1 = trk->pVec();
                    protonSignsY *= proton1.Y();
                }
                else if (k == 1)
                {
                    proton2 = trk->pVec();
                    protonSignsY *= proton2.Y();
                }     
            }    
        }

        else if (isMC == 1)
        {
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

            protonSignsY = protonSignsY*vProtons[0]->Py()*vProtons[1]->Py();
        }

        FillProtons(protonSignsY, HistProtonsSameSide, HistProtonsOppositeSide,  iCutFlow); //all data

        // SELECTION: number of primary vertices
        Int_t numberOfPrimaryVertices = upcEvt->getNumberOfVertices();       	
        HistNumPrimaryVertices->Fill(numberOfPrimaryVertices);	      
        
        if(numberOfPrimaryVertices!=1)
        {
            vProtons.clear();   
            continue;
        }		
        iCutFlow+=1;
		HistCutFlow->Fill(iCutFlow); // cutflow: one vertex
        HistNumPrimaryVerticesPostSelection->Fill(numberOfPrimaryVertices);
        FillProtons(protonSignsY, HistProtonsSameSide, HistProtonsOppositeSide,  iCutFlow); // protons - 1 vertex


		// SELECTION: primary vertex is placed within |zvtx| < 80 cm lVecKaonComb2b.pt() << ", " << endl;
        HistPrimaryVertexAbsPosZ->Fill(abs(primVertexPosZ));
		if(abs(upcEvt->getVertex(0)->getPosZ())>=80)
		{
            vProtons.clear();
			continue;
        }   
        iCutFlow+=1;
		HistCutFlow->Fill(iCutFlow); // cutflow: primary vertex |zvtx| < 80 cm
        HistPrimaryVertexAbsPosZPostSelection->Fill(abs(upcEvt->getVertex(0)->getPosZ()));      
        FillProtons(protonSignsY, HistProtonsSameSide, HistProtonsOppositeSide,  iCutFlow); // protons - |zvtx| < 80 cm

		vector <StUPCTrack const*> tracksWithTofHit;
		vector <StUPCTrack const*> tracksWithoutTofHit;
        for (Int_t i = 0; i<upcEvt->getNumberOfTracks(); i++)
		{
			if(upcEvt->getTrack(i)->getFlag(StUPCTrack::kTof))
			{
                tracksWithTofHit.push_back(upcEvt->getTrack(i));
			} 
            else
            {
                tracksWithoutTofHit.push_back(upcEvt->getTrack(i));
            }
		}       		



		// SELECTION: four TOF-matched primary tracks
        HistNumTofMatchedTracks->Fill(tracksWithTofHit.size());
		bool isValidNumberOfTofMatchedTracks = 0;

    
		if (tracksWithTofHit.size() >= 3 and  tracksWithTofHit.size() <=4)
		{
			isValidNumberOfTofMatchedTracks = 1;
		}

		if (isValidNumberOfTofMatchedTracks == 0) 
        {
            vProtons.clear();
            tracksWithTofHit.clear();
            tracksWithoutTofHit.clear();
            continue;
        } 
        iCutFlow+=1;
		HistCutFlow->Fill(iCutFlow); // cutflow: four TOF-matched tracks        
        HistNumTofMatchedTracksPostSelection->Fill(tracksWithTofHit.size());
        FillProtons(protonSignsY, HistProtonsSameSide, HistProtonsOppositeSide,  iCutFlow); // protons - 4 tof


        
        vector <StUPCTrack const*> tracksWithoutTofHitGoodChargeThreeTof;
        vector <StUPCTrack const*> tracksWithoutTofHitGoodChargeTwoTof;
        int sumChargeTwoTracksTof = 0;
        int sumChargeThreeTracksTof = 0;
        if (tracksWithTofHit.size() == 4)
        {
            // SELECTION: TOF-matched tracks are of opposied signs		
            int sumCharge = 0;
            for (unsigned long int i = 0; i < tracksWithTofHit.size(); i++)
            {
                HistTofMatchedTracksCharge->Fill(tracksWithTofHit[i]->getCharge());
                sumCharge+=tracksWithTofHit[i]->getCharge();
            }

            if (sumCharge != 0)  
            {
                vProtons.clear();
                tracksWithTofHit.clear();
                tracksWithTofHit.clear();
                continue;
            }

            iCutFlow+=1;
            HistCutFlow->Fill(iCutFlow);  // cutflow: sum of the charge is zero
            FillProtons(protonSignsY, HistProtonsSameSide, HistProtonsOppositeSide,  iCutFlow); // charge 0

            sumCharge = 0;
            for (unsigned long int i = 0; i < tracksWithTofHit.size(); i++)
            {
                HistTofMatchedTracksChargePostSelection->Fill(tracksWithTofHit[i]->getCharge());       
            }
        }

        else if (tracksWithTofHit.size() == 3)
        {
            sumChargeThreeTracksTof = tracksWithTofHit[0]->getCharge() + tracksWithTofHit[1]->getCharge() + tracksWithTofHit[2]->getCharge();    

            if (sumChargeThreeTracksTof == 1)
            {
                for (int i = 0; i < tracksWithoutTofHit.size(); i++)
                {
                    if (tracksWithoutTofHit[i]->getCharge() == -1)
                    {
                        tracksWithoutTofHitGoodChargeThreeTof.push_back(tracksWithoutTofHit[i]);     
                    }
                }
            }

            else if (sumChargeThreeTracksTof == -1)
            {
                for (int i = 0; i < tracksWithoutTofHit.size(); i++)
                {
                    if (tracksWithoutTofHit[i]->getCharge() == 1)
                    {
                        tracksWithoutTofHitGoodChargeThreeTof.push_back(tracksWithoutTofHit[i]);     
                    }
                }
            }

            else 
            {
                vProtons.clear();
                tracksWithTofHit.clear();
                tracksWithoutTofHit.clear();
                tracksWithoutTofHitGoodChargeThreeTof.clear();
                continue;
            }

            if (tracksWithoutTofHitGoodChargeThreeTof.size() == 0) // delete event if there is no complementary trkcs of good sign without TOF
            {
                vProtons.clear();
                tracksWithTofHit.clear();
                tracksWithoutTofHit.clear();
                tracksWithoutTofHitGoodChargeThreeTof.clear();
                continue;
            }
        }

        else if (tracksWithTofHit.size() == 2)
        {
            sumChargeTwoTracksTof = tracksWithTofHit[0]->getCharge() + tracksWithTofHit[1]->getCharge();   

            if (sumChargeTwoTracksTof == 2)
            {
                for (int i = 0; i < tracksWithoutTofHit.size(); i++)
                {
                    if (tracksWithoutTofHit[i]->getCharge() == -1)
                    {
                        tracksWithoutTofHitGoodChargeTwoTof.push_back(tracksWithoutTofHit[i]);     
                    }
                }
            }

            else if (sumChargeTwoTracksTof == -2)
            {
                for (int i = 0; i < tracksWithoutTofHit.size(); i++)
                {
                    if (tracksWithoutTofHit[i]->getCharge() == 1)
                    {
                        tracksWithoutTofHitGoodChargeTwoTof.push_back(tracksWithoutTofHit[i]);     
                    }
                }
            }

            else if (sumChargeTwoTracksTof == 0)
            {
                for (int i = 0; i < tracksWithoutTofHit.size(); i++)
                {
                    //if (tracksWithoutTofHit[i]->getCharge() == 1)
                   // {
                        tracksWithoutTofHitGoodChargeTwoTof.push_back(tracksWithoutTofHit[i]);     
                    //}
                }
            }
            else 
            {
                vProtons.clear();
                tracksWithTofHit.clear();
                tracksWithoutTofHit.clear();
                tracksWithoutTofHitGoodChargeThreeTof.clear();
                tracksWithoutTofHitGoodChargeTwoTof.clear();
                continue;
            }

            if (tracksWithoutTofHitGoodChargeTwoTof.size() == 0) // delete event if there is no complementary trkcs of good sign without TOF
            {
                vProtons.clear();
                tracksWithTofHit.clear();
                tracksWithoutTofHit.clear();
                tracksWithoutTofHitGoodChargeThreeTof.clear();
                tracksWithoutTofHitGoodChargeTwoTof.clear();
                continue;
            }
        }

    
        // select good quality tracks in fiducial region from tracks without TOF hit.
        vector <StUPCTrack const *> goodComplementaryTracks3Tof;
        for (int i = 0; i < tracksWithoutTofHitGoodChargeThreeTof.size(); i++)
        {
            StUPCTrack  const*  tempTrack = tracksWithoutTofHitGoodChargeThreeTof[i];
            if (abs(tempTrack->getEta()) < 0.9 and  tempTrack->getPt() > 0.2 and tempTrack->getNhitsFit() >= 25 and tempTrack->getNhitsDEdx() >= 15)
            {
                goodComplementaryTracks3Tof.push_back(tempTrack);
            }
        }

        if (tracksWithTofHit.size() == 3)
        {
            if (goodComplementaryTracks3Tof.size()==0)
            {
                vProtons.clear();
                tracksWithTofHit.clear();
                tracksWithoutTofHit.clear();
                tracksWithoutTofHitGoodChargeThreeTof.clear();
                tracksWithoutTofHitGoodChargeTwoTof.clear();
                continue;
            }
        }

        // select good quality tracks in fiducial region from tracks without TOF hit.
        vector <StUPCTrack const *> goodComplementaryTracks2Tof;
        for (int i = 0; i < tracksWithoutTofHitGoodChargeTwoTof.size(); i++)
        {
            StUPCTrack  const*  tempTrack = tracksWithoutTofHitGoodChargeTwoTof[i];
            if (abs(tempTrack->getEta()) < 0.9 and  tempTrack->getPt() > 0.2 and tempTrack->getNhitsFit() >= 25 and tempTrack->getNhitsDEdx() >= 15)
            {
                goodComplementaryTracks2Tof.push_back(tempTrack);
            }
        }

        if (tracksWithTofHit.size() == 2)
        {
            if (goodComplementaryTracks2Tof.size() <= 1)
            {
                vProtons.clear();
                tracksWithTofHit.clear();
                tracksWithoutTofHit.clear();
                tracksWithoutTofHitGoodChargeThreeTof.clear();
                tracksWithoutTofHitGoodChargeTwoTof.clear();
                continue;
            }
        }


        bool isEtaValid = 1;
        for (unsigned long int i = 0; i < tracksWithTofHit.size(); i++)
        {	
            HistTofMatchedTracksAbsEta->Fill(abs(tracksWithTofHit[i]->getEta()));
            if (abs(tracksWithTofHit[i]->getEta()) >= 0.9)
            {
                isEtaValid = 0;
            }
        }

        if (isEtaValid == 0) 
        {
            vProtons.clear();
            tracksWithTofHit.clear();
            continue;
        }
        iCutFlow+=1;
		HistCutFlow->Fill(iCutFlow);  // cutflow: TOF-matched tracks within appropraite pseudorapidity range |eta| < 0.9
        FillProtons(protonSignsY, HistProtonsSameSide, HistProtonsOppositeSide,  iCutFlow); // protons:  TOF-matched tracks within appropraite pseudorapidity range |eta| < 0.9
        for (unsigned long int i = 0; i < tracksWithTofHit.size(); i++)
        {	
            HistTofMatchedTracksAbsEtaPostSelection->Fill(abs(tracksWithTofHit[i]->getEta()));
        }

        // SELECTION pT > 0.2 GeV
        bool isPtValid = 1;
        for (unsigned long int i = 0; i < tracksWithTofHit.size(); i++)
        {	
            HistTofMatchedTracksPt->Fill((tracksWithTofHit[i]->getPt()));
            if ((tracksWithTofHit[i]->getPt()) < 0.2)
            {
                isPtValid = 0;
            }
        }
        if (isPtValid == 0) 
        {
            vProtons.clear();
            tracksWithTofHit.clear();
            continue;
        }
        iCutFlow+=1;
		HistCutFlow->Fill(iCutFlow);  // cutflow: TOF-matched tracks within appropraite pseudorapidity range |eta| < 0.9
        FillProtons(protonSignsY, HistProtonsSameSide, HistProtonsOppositeSide,  iCutFlow); // protons:  TOF-matched tracks within appropraite pseudorapidity range |eta| < 0.9
        for (unsigned long int i = 0; i < tracksWithTofHit.size(); i++)
        {	
            HistTofMatchedTracksPtPostSelection->Fill((tracksWithTofHit[i]->getPt()));
        }


		// SELECTION: associated global tracks satisfy quality criteria N_{fit} >= 25
        bool isNfitAtLeast25 = 1;
        for (unsigned long int i = 0; i < tracksWithTofHit.size(); i++)
        {		
            HistTofMatchedTracksNfit->Fill(tracksWithTofHit[i]->getNhitsFit());
            if (tracksWithTofHit[i]->getNhitsFit() < 25)
            {
                isNfitAtLeast25 = 0;
            }
        }
        if (isNfitAtLeast25 == 0)
        {
            vProtons.clear();
            tracksWithTofHit.clear();
            continue;
        }
        iCutFlow+=1;
		HistCutFlow->Fill(iCutFlow);    // cutflow: associated global tracks quality: N_{fit} >= 25
        FillProtons(protonSignsY, HistProtonsSameSide, HistProtonsOppositeSide,  iCutFlow); // protons:  associated global tracks quality: N_{fit} >= 25

        for (unsigned long int i = 0; i < tracksWithTofHit.size(); i++)
        {		
            HistTofMatchedTracksNfitPostSelection->Fill(tracksWithTofHit[i]->getNhitsFit());
        }

		// SELECTION: associated global tracks satisfy quality criteria N_{dE/dx} >= 15 
        bool isNdEdxAtLeast15 = 1;
        for (unsigned long int i = 0; i < tracksWithTofHit.size(); i++)
        {		
            HistTofMatchedTracksNdEdx->Fill(tracksWithTofHit[i]->getNhitsDEdx());
            if (tracksWithTofHit[i]->getNhitsDEdx() < 15)
            {
                isNdEdxAtLeast15 = 0;
            }
        }
        if (isNdEdxAtLeast15 == 0) 
        {
            vProtons.clear();
            tracksWithTofHit.clear();
            continue;
        }
        iCutFlow+=1;
		HistCutFlow->Fill(iCutFlow);  // cutflow: associated global tracks quality N_{dE/dx} >= 15 
        FillProtons(protonSignsY, HistProtonsSameSide, HistProtonsOppositeSide,  iCutFlow); // protons:  associated global tracks quality: N_{dE/dx} >= 15
        for (unsigned long int i = 0; i < tracksWithTofHit.size(); i++)
        {		
            HistTofMatchedTracksNdEdxPostSelection->Fill(tracksWithTofHit[i]->getNhitsDEdx());
        }



        // SELECTION: both kaons in mass rangeisMC
        double kaonMassWindowNarrowLow = 0.48;
        double kaonMassWindowNarrowHigh = 0.51;

        double kaonMassWindowWideLow = 0.44;
        double kaonMassWindowWideHigh = 0.54;
        double distSumSquaredComb1, distSumSquaredComb2;

        vector <StUPCTrack const *> posPion;
        vector <StUPCTrack const *> negPion;
        bool fourComb1a = false;
        bool fourComb1b = false;
        bool fourComb2a = false;
        bool fourComb2b = false;
        vector <StUPCTrack const *> vPosNegPionLeadingKaon;
        vector <StUPCTrack const *> vPosNegPionSubLeadingKaon;

        if (tracksWithTofHit.size() == 4)
        {
            for (unsigned long int i = 0; i < tracksWithTofHit.size(); i++)
            {
        
                if( tracksWithTofHit[i]->getCharge() == 1)
                {
                    posPion.push_back(tracksWithTofHit[i]);
                }
                else if (tracksWithTofHit[i]->getCharge() == -1)
                {
                    negPion.push_back(tracksWithTofHit[i]);
                }
            }
            
            TVector3 const tryVec(0,0,0);
            StUPCV0 *lVecKaonComb1a = new StUPCV0(posPion[0], negPion[0], massPion, massPion, 1, 1, tryVec, upcEvt->getMagneticField(), bool(isMC), true);
            StUPCV0 *lVecKaonComb1b = new StUPCV0(posPion[1], negPion[1], massPion, massPion, 1, 1, tryVec, upcEvt->getMagneticField(), bool(isMC), true);
            StUPCV0 *lVecKaonComb2a = new StUPCV0(posPion[0], negPion[1], massPion, massPion, 1, 1, tryVec, upcEvt->getMagneticField(), bool(isMC), true);
            StUPCV0 *lVecKaonComb2b = new StUPCV0(posPion[1], negPion[0], massPion, massPion, 1, 1, tryVec, upcEvt->getMagneticField(), bool(isMC), true);
        
            distSumSquaredComb1 = sqrt(pow((lVecKaonComb1a->m()-massKaon), 2) + pow((lVecKaonComb1b->m()-massKaon), 2));
            distSumSquaredComb2 = sqrt(pow((lVecKaonComb2a->m()-massKaon), 2) + pow((lVecKaonComb2b->m()-massKaon), 2));
    
            if (distSumSquaredComb1 < distSumSquaredComb2)
            {
                if (lVecKaonComb1a->pt() > lVecKaonComb1b->pt())
                {
                    vPosNegPionLeadingKaon.push_back(posPion[0]);
                    vPosNegPionLeadingKaon.push_back(negPion[0]);
                    vPosNegPionSubLeadingKaon.push_back(posPion[1]);
                    vPosNegPionSubLeadingKaon.push_back(negPion[1]);
                }
                else
                {
                    vPosNegPionLeadingKaon.push_back(posPion[1]);
                    vPosNegPionLeadingKaon.push_back(negPion[1]);
                    vPosNegPionSubLeadingKaon.push_back(posPion[0]);
                    vPosNegPionSubLeadingKaon.push_back(negPion[0]);
                }
            }

            else
            {
                if (lVecKaonComb2a->pt() > lVecKaonComb2b->pt())
                {
                    vPosNegPionLeadingKaon.push_back(posPion[0]);
                    vPosNegPionLeadingKaon.push_back(negPion[1]);
                    vPosNegPionSubLeadingKaon.push_back(posPion[1]);
                    vPosNegPionSubLeadingKaon.push_back(negPion[0]);
                }
                else
                {
                    vPosNegPionLeadingKaon.push_back(posPion[1]);
                    vPosNegPionLeadingKaon.push_back(negPion[0]);
                    vPosNegPionSubLeadingKaon.push_back(posPion[0]);
                    vPosNegPionSubLeadingKaon.push_back(negPion[1]);
                }
            }

            delete lVecKaonComb1a;
            delete lVecKaonComb1b;
            delete lVecKaonComb2a;
            delete lVecKaonComb2b;
        }



        
    
        else if (tracksWithTofHit.size() == 3)
        {
            int complementaryCharge = goodComplementaryTracks3Tof[0]->getCharge();
            
            for (unsigned long int i = 0; i < tracksWithTofHit.size(); i++)
            {
        
                if( tracksWithTofHit[i]->getCharge() == 1)
                {
                    posPion.push_back(tracksWithTofHit[i]);
                }
                else if (tracksWithTofHit[i]->getCharge() == -1)
                {
                    negPion.push_back(tracksWithTofHit[i]);
                }
            }

            vector <double> dist1;
            vector <double> dist2;  
            vector <vector <StUPCTrack const *>> tracks1;
            vector <vector <StUPCTrack const *>> tracks2;

            for (int i = 0; i < goodComplementaryTracks3Tof.size(); i ++)
            {              
                if (complementaryCharge == 1)
                {
                    posPion[1] = goodComplementaryTracks3Tof[i];
                }
                
                else if (complementaryCharge == -1)
                {
                   negPion[1] = goodComplementaryTracks3Tof[i]; 
                }

                TVector3 const tryVec(0,0,0);
                StUPCV0 lVecKaonComb1aa(posPion[0], negPion[0], massPion, massPion, 1, 1, tryVec, upcEvt->getMagneticField(), bool(isMC), true);
                StUPCV0 lVecKaonComb1bb(posPion[1], negPion[1], massPion, massPion, 1, 1, tryVec, upcEvt->getMagneticField(), bool(isMC), true);
                vector <StUPCTrack const*> temp;
                temp.push_back(posPion[0]);
                temp.push_back(negPion[0]);
                temp.push_back(posPion[1]);
                temp.push_back(negPion[1]);
                tracks1.push_back(temp);

                distSumSquaredComb1 = sqrt(pow((lVecKaonComb1aa.m()-massKaon), 2) + pow((lVecKaonComb1bb.m()-massKaon), 2));
                dist1.push_back(distSumSquaredComb1);
            }
        
            for (int i = 0; i < goodComplementaryTracks3Tof.size(); i ++)
            {    
                if (complementaryCharge == 1)
                {
                    posPion[1] = goodComplementaryTracks3Tof[i];
                }
                
                else if (complementaryCharge == -1)
                {
                   negPion[1] = goodComplementaryTracks3Tof[i]; 
                }        

                TVector3 const tryVec(0,0,0);
                StUPCV0 lVecKaonComb2aa(posPion[0], negPion[1], massPion, massPion, 1, 1, tryVec, upcEvt->getMagneticField(), bool(isMC), true);
                StUPCV0 lVecKaonComb2bb(posPion[1], negPion[0], massPion, massPion, 1, 1, tryVec, upcEvt->getMagneticField(), bool(isMC), true);
                vector <StUPCTrack const*> temp;
                temp.push_back(posPion[0]);
                temp.push_back(negPion[1]);
                temp.push_back(posPion[1]);
                temp.push_back(negPion[0]);
                tracks2.push_back(temp);
                distSumSquaredComb2 = sqrt(pow((lVecKaonComb2aa.m()-massKaon), 2) + pow((lVecKaonComb2bb.m()-massKaon), 2));
                dist2.push_back(distSumSquaredComb2);
            }

            auto minElement1 = std::min_element(dist1.begin(), dist1.end());
            int minIndex1= std::distance(dist1.begin(), minElement1);
            
            auto minElement2 = std::min_element(dist2.begin(), dist2.end());
            int minIndex2= std::distance(dist2.begin(), minElement2);

            if (dist1[minIndex1] < dist2[minIndex2])
            {
                TVector3 const tryVec(0,0,0);
                StUPCV0 l1a(tracks1[minIndex1][0], tracks1[minIndex1][1],massPion, massPion, 1, 1, tryVec, upcEvt->getMagneticField(), bool(isMC), true);
                StUPCV0 l1b(tracks1[minIndex1][2], tracks1[minIndex1][3],massPion, massPion, 1, 1, tryVec, upcEvt->getMagneticField(), bool(isMC), true);
                
                if (l1a.pt() > l1b.pt())
                {
                    vPosNegPionLeadingKaon.push_back(tracks1[minIndex1][0]);
                    vPosNegPionLeadingKaon.push_back(tracks1[minIndex1][1]);
                    vPosNegPionSubLeadingKaon.push_back(tracks1[minIndex1][2]);
                    vPosNegPionSubLeadingKaon.push_back(tracks1[minIndex1][3]);
                }

                else
                {
                    vPosNegPionLeadingKaon.push_back(tracks1[minIndex1][2]);
                    vPosNegPionLeadingKaon.push_back(tracks1[minIndex1][3]);
                    vPosNegPionSubLeadingKaon.push_back(tracks1[minIndex1][0]);
                    vPosNegPionSubLeadingKaon.push_back(tracks1[minIndex1][1]);
                }
            }

            else
            {
                TVector3 const tryVec(0,0,0);
                StUPCV0 l2a(tracks2[minIndex2][0], tracks2[minIndex2][1],massPion, massPion, 1, 1, tryVec, upcEvt->getMagneticField(), bool(isMC), true);
                StUPCV0 l2b(tracks2[minIndex2][2], tracks2[minIndex2][3],massPion, massPion, 1, 1, tryVec, upcEvt->getMagneticField(), bool(isMC), true);
                
                if (l2a.pt() > l2b.pt())
                {
                    vPosNegPionLeadingKaon.push_back(tracks2[minIndex2][0]);
                    vPosNegPionLeadingKaon.push_back(tracks2[minIndex2][1]);
                    vPosNegPionSubLeadingKaon.push_back(tracks2[minIndex2][2]);
                    vPosNegPionSubLeadingKaon.push_back(tracks2[minIndex2][3]);
                }

                else
                {
                    vPosNegPionLeadingKaon.push_back(tracks2[minIndex2][2]);
                    vPosNegPionLeadingKaon.push_back(tracks2[minIndex2][3]);
                    vPosNegPionSubLeadingKaon.push_back(tracks2[minIndex2][0]);
                    vPosNegPionSubLeadingKaon.push_back(tracks2[minIndex2][1]);
                }
            }
        }

        else if (tracksWithTofHit.size() == 2)
        {
            for (unsigned long int i = 0; i < tracksWithTofHit.size(); i++)
            {
        
                if( tracksWithTofHit[i]->getCharge() == 1)
                {
                    posPion.push_back(tracksWithTofHit[i]);
                }
                else if (tracksWithTofHit[i]->getCharge() == -1)
                {
                    negPion.push_back(tracksWithTofHit[i]);
                }
            }

            if (posPion.size() == 2)
            {
                int sumMinus2 = 0;
                for (int i = 0; i < goodComplementaryTracks2Tof.size(); i++)
                {
                    if (goodComplementaryTracks2Tof[i]->getCharge() == -1 )
                    {
                        sumMinus2+=1;
                    }
                }
                
                if (sumMinus2 <= 1)
                {
                    continue;
                }
                double distance1;



            }

            else if (negPion.size() == 2)
            {
                int sumPlus2 = 0;
                for (int i = 0; i < goodComplementaryTracks2Tof.size(); i++)
                {
                    if (goodComplementaryTracks2Tof[i]->getCharge() == 1 )
                    {
                        sumPlus2+=1;
                    }
                }
                
                if (sumPlus2 <= 1)
                {
                    continue;
                }
            }

            else if (negPion.size() == 1 and posPion.size() == 1)
            {
                
            }

            else 
            {
                cout << "HEREEE we should not be xd" << endl;
            }
        }


        else
        {
            continue;
        }


        if (tracksWithTofHit.size() == 3)// or tracksWithTofHit.size() == 3)
        {
            
        bool areKaonsInNarrowMassWindow = 0;
        bool areKaonsInWideMassWindow = 0;
        bool isPtMissingSmall = 0;
        bool isNumberOfTofClusterSmall = 0;
        TVector3 const tryVec(0,0,0);
        StUPCV0 leadingKaon(vPosNegPionLeadingKaon[0], vPosNegPionLeadingKaon[1], massPion, massPion, 1, 1, tryVec, upcEvt->getMagneticField(), bool(isMC), true);
        StUPCV0 subLeadingKaon(vPosNegPionSubLeadingKaon[0], vPosNegPionSubLeadingKaon[1], massPion, massPion, 1, 1, tryVec, upcEvt->getMagneticField(), bool(isMC), true);

        if ( (leadingKaon.m() >= kaonMassWindowWideLow and  leadingKaon.m() < kaonMassWindowWideHigh) and (subLeadingKaon.m() >= kaonMassWindowWideLow and  subLeadingKaon.m() < kaonMassWindowWideHigh) )
        {
              areKaonsInWideMassWindow = 1;      
        }

        if ( (leadingKaon.m() >= kaonMassWindowNarrowLow and  leadingKaon.m() < kaonMassWindowNarrowHigh) and (subLeadingKaon.m() >= kaonMassWindowNarrowLow and  subLeadingKaon.m() < kaonMassWindowNarrowHigh) )
        {
              areKaonsInNarrowMassWindow = 1;      
        }

        if (areKaonsInWideMassWindow != 1) // delete event if kaon masses are not in mass window
        {
            vProtons.clear();
            tracksWithTofHit.clear();
            posPion.clear();
            negPion.clear();
            continue; 
        }

        // SELECTION: pT miss
        double pXmiss = leadingKaon.px() + subLeadingKaon.px() + proton1.X() + proton2.X();
        double pYmiss = leadingKaon.py() + subLeadingKaon.py()  + proton1.Y() + proton2.Y();
        double pTmiss = sqrt(pow(pXmiss,2) + pow(pYmiss,2));


        if (pTmiss <= 0.15) 
        {
                isPtMissingSmall = 1;
        }

        // SELECTION: number of cluster 
        Int_t nTofHits = upcEvt->getNumberOfHits();           
        vector <Int_t> vTray;
        vector <Int_t> vTrayUniqueVector;
        vector <Int_t> vTrayUniqueSet;
        vector <Int_t> vMmodule;
        
        for (Int_t i = 0; i < nTofHits; i++)
        {
            vTray.push_back(Int_t(upcEvt->getHit(i)->getTray()));
            vTrayUniqueVector.push_back(Int_t(upcEvt->getHit(i)->getTray()));
            vMmodule.push_back(Int_t(upcEvt->getHit(i)->getModule()));
        }

        sort(vTrayUniqueVector.begin(), vTrayUniqueVector.end());   
        auto last = unique(vTrayUniqueVector.begin(), vTrayUniqueVector.end());   
        vTrayUniqueVector.erase(last, vTrayUniqueVector.end()); 

        vector <vector <Int_t>> vModuleUniqueTray;
        for (int i = 0; i < vTrayUniqueVector.size(); i++)
        {
            vector <Int_t> vModuleUnique;
            for (int j = 0; j < nTofHits; j++)
            {
                if (vTrayUniqueVector[i] == vTray[j] )
                {
                    vModuleUnique.push_back(vMmodule[j]);
                }
            }
            vModuleUniqueTray.push_back(vModuleUnique);
            vModuleUnique.clear();
        }

        int totalCluster = 0;
        for (int i  = 0; i < vModuleUniqueTray.size(); i++)
        {
            vector <Int_t> vec =  vModuleUniqueTray[i] ;  
            sort(vec.begin(), vec.end());   
            auto last = unique(vec.begin(), vec.end());   
            vec.erase(last, vec.end()); 
            
            if (vec.size() == 1)
            {
                totalCluster+=1;
            }

            for (int j = 0; j < vec.size()-1; j++)
            {
                Int_t modNum = vec[j];
                int diff = 1;
                int num = 0;
                for (int z = j+1; z < vec.size(); z++)
                {
                    if (modNum+diff == vec[z])
                    {
                        num+=0;
                        diff+=1;

                        if (z == j+1)
                        {
                            num+=1;
                        }
                    }

                    else if ( j == (vec.size()-2) and  vec[vec.size()-2]+1 != vec[vec.size()-1])
                    {
                        num+=2; 
                        continue;
                    }

                    else
                    {
                        num+=1;
                    }

                    j+=(diff);
                }
                totalCluster+=num;
            }
        }


        if (totalCluster <= 9)
        {
            isNumberOfTofClusterSmall = 1;
        }

        //DCA Daughters
        double dcaDaughtersLeadingKaon = leadingKaon.dcaDaughters();
        double dcaDaughtersSubLeadingKaon = subLeadingKaon.dcaDaughters();
        

       // Beam Position based on fill number
        int nFillNumber = upcEvt->getFillNumber();
        if (isMC == 0)
        {
		    beamPositionX = FindPosition(nFillNumber, primVertexPosZ, fillNumberWithPosition[0], fillNumberWithPosition[1], fillNumberWithPosition[2], fillNumberWithPosition[3], fillNumberWithPosition[4], fillNumberWithPosition[5], fillNumberWithPosition[6], fillNumberWithPosition[7], fillNumberWithPosition[8])[0];
		    beamPositionY = FindPosition(nFillNumber, primVertexPosZ, fillNumberWithPosition[0], fillNumberWithPosition[1], fillNumberWithPosition[2], fillNumberWithPosition[3], fillNumberWithPosition[4], fillNumberWithPosition[5], fillNumberWithPosition[6], fillNumberWithPosition[7],fillNumberWithPosition[8] )[1];			  
        }

        else if (isMC == 1)
        {
            beamPositionX = 0.0;
            beamPositionY = 0.0;
        }

        // vertices
        TVector3 vertexLeadingKaon = leadingKaon.decayVertex();
        double leadingKaonVtxX = vertexLeadingKaon.X();
        double leadingKaonVtxY = vertexLeadingKaon.Y();
        double leadingKaonVtxZ = vertexLeadingKaon.Z();

        TVector3 vertexSubLeadingKaon = subLeadingKaon.decayVertex();
        double subLeadingKaonVtxX = vertexSubLeadingKaon.X();
        double subLeadingKaonVtxY = vertexSubLeadingKaon.Y();
        double subLeadingKaonVtxZ = vertexSubLeadingKaon.Z();
  

        double leadingKaonPx = leadingKaon.px();
        double leadingKaonPy = leadingKaon.py();
        double leadingKaonPz = leadingKaon.pz();

        double subLeadingKaonPx = subLeadingKaon.px();
        double subLeadingKaonPy = subLeadingKaon.py();
        double subLeadingKaonPz = subLeadingKaon.pz();

        double nx = leadingKaonPy*subLeadingKaonPz - leadingKaonPz*subLeadingKaonPy;
        double ny = leadingKaonPz*subLeadingKaonPx - leadingKaonPx*subLeadingKaonPz;
        double nz = leadingKaonPx*subLeadingKaonPy - leadingKaonPy*subLeadingKaonPx;
                    
        double dx = leadingKaonVtxX - subLeadingKaonVtxX;
        double dy = leadingKaonVtxY - subLeadingKaonVtxY;
        double dz = leadingKaonVtxZ - subLeadingKaonVtxZ;
             
        double ds = sqrt(pow(nx,2)+ pow(ny,2) + pow(nz,2)) ;
            
        double distance = (abs(nx*dx+dy*ny+dz*nz))/(ds);

        double distVertexBeamxLeadingKaon = leadingKaonVtxX - beamPositionY;
        double distVerteyBeamyLeadingKaon = leadingKaonVtxY - beamPositionX;
        double distRLeadingKaon = sqrt(pow(distVertexBeamxLeadingKaon,2) + pow(distVerteyBeamyLeadingKaon,2));

        double distVertexBeamxSubLeadingKaon = subLeadingKaonVtxX - beamPositionY;
        double distVerteyBeamySubLeadingKaon = subLeadingKaonVtxY - beamPositionX;
        double distRSubLeadingKaon = sqrt(pow(distVertexBeamxSubLeadingKaon,2) + pow(distVerteyBeamySubLeadingKaon,2));  


        double t1x = subLeadingKaonPy*nz - subLeadingKaonPz*ny;
        double t1y = subLeadingKaonPz*nx - subLeadingKaonPx*nz;
        double t1z = subLeadingKaonPx*ny - subLeadingKaonPy*nx;

        double t2x = leadingKaonPy*nz - leadingKaonPz*ny;
        double t2y = leadingKaonPz*nx - leadingKaonPx*nz;
        double t2z = leadingKaonPx*ny - leadingKaonPy*nx;


        double t1 = (t1x*dx + t1y*dy + t1z*dz)/ds/ds;
        double t2 = (t2x*dx + t2y*dy + t2z*dz)/ds/ds;

        double x1 = leadingKaonVtxX+t1*leadingKaonPx;
        double y1 = leadingKaonVtxY+t1*leadingKaonPy;
        double z1 = leadingKaonVtxZ+t1*leadingKaonPz;             

        double x11 = subLeadingKaonVtxX+t2*subLeadingKaonPx;
        double y11 = subLeadingKaonVtxY+t2*subLeadingKaonPy;
        double z11 = subLeadingKaonVtxZ+t2*subLeadingKaonPz;             

        double xs = (x1+x11)/2.0;
        double ys = (y1+y11)/2.0;   
        double zs = (z1+z11)/2.0;
        double xy = sqrt(pow(xs,2)+pow(ys,2));
        xyVec.push_back(xy);

        xVec.push_back(xs);
        yVec.push_back(ys);
        zVec.push_back(zs);
        HistZK0K0->Fill(zs);


        HistInvMassPiPi2D->Fill(leadingKaon.m(), subLeadingKaon.m());
        // TH1D* HistEtaLeadingSubLeadingDiff = new TH1D("HistEtaLeadingSubLeadingDiff", " ; #eta_{K^{0}}^{leading} [GeV] ; events", 100, -1, 1);
        if (areKaonsInNarrowMassWindow)
        {

            // DCA  
            HistClosestDistanceBetweenKaons->Fill(distance);
            HistDcaDaughtersLeadingKaon->Fill(dcaDaughtersLeadingKaon);
            HistDcaDaughtersLeadingKaon->Fill(dcaDaughtersSubLeadingKaon);

            if (dcaDaughtersLeadingKaon > 1.5 or dcaDaughtersSubLeadingKaon > 1.5 or distance >= 1.5)
            {
                    continue;
            }


            HistPtMiss->Fill(pTmiss);
            // VERTEX
            HistDistVtxXLeadingKaon->Fill(distVertexBeamxLeadingKaon);
            HistDistVtxYLeadingKaon->Fill(distVerteyBeamyLeadingKaon);
            HistDistVtxZLeadingKaon->Fill(leadingKaonVtxZ);
            HistRLeadingKaon->Fill(distRLeadingKaon);
            HistDistVtxXLeadingKaon->Fill(distVertexBeamxSubLeadingKaon);
            HistDistVtxYLeadingKaon->Fill(distVerteyBeamySubLeadingKaon);
            HistDistVtxZLeadingKaon->Fill(subLeadingKaonVtxZ);
            HistRLeadingKaon->Fill(distRSubLeadingKaon);
            
            HistNumOfClusters->Fill(totalCluster);

            // px, py plots for 4 pions + 2 protons - MISSING PT
            HistPtProtonPtKaonsX->Fill(leadingKaon.px()+subLeadingKaon.px(), proton1.X() +  proton2.X());
            HistPtProtonPtKaonsY->Fill(leadingKaon.py()+subLeadingKaon.py(), proton1.Y() +  proton2.Y());
            HistPtProtonPtKaonsX1d->Fill(leadingKaon.px()+subLeadingKaon.px() + proton1.X() +  proton2.X());
            HistPtProtonPtKaonsY1d->Fill(leadingKaon.py()+subLeadingKaon.py() +proton1.Y() +  proton2.Y());

            if (protonSignsY < 0.0)
            {
                HistPtProtonPtKaonsXOppY->Fill(leadingKaon.px()+subLeadingKaon.px(), proton1.X() +  proton2.X());
                HistPtProtonPtKaonsYOppY->Fill(leadingKaon.py()+subLeadingKaon.py(), proton1.Y() +  proton2.Y());
                HistPtProtonPtKaonsX1dOppY->Fill(leadingKaon.px()+subLeadingKaon.px() + proton1.X() +  proton2.X());
                HistPtProtonPtKaonsY1dOppY->Fill(leadingKaon.py()+subLeadingKaon.py() +proton1.Y() +  proton2.Y());
            }

            else
            {

                HistPtProtonPtKaonsXSameY->Fill(leadingKaon.px()+subLeadingKaon.px(), proton1.X() +  proton2.X());
                HistPtProtonPtKaonsYSameY->Fill(leadingKaon.py()+subLeadingKaon.py(), proton1.Y() +  proton2.Y());
                HistPtProtonPtKaonsX1dSameY->Fill(leadingKaon.px()+subLeadingKaon.px() + proton1.X() +  proton2.X());
                HistPtProtonPtKaonsY1dSameY->Fill(leadingKaon.py()+subLeadingKaon.py() +proton1.Y() +  proton2.Y());

                if (proton1.Y() >= 0.0 and proton2.Y() >= 0.0)
                {
                    HistPtProtonPtKaonsXSameYMomentaGreaterThan0->Fill(leadingKaon.px()+subLeadingKaon.px(), proton1.X() +  proton2.X());
                    HistPtProtonPtKaonsYSameYMomentaGreaterThan0->Fill(leadingKaon.py()+subLeadingKaon.py(), proton1.Y() +  proton2.Y());
                    HistPtProtonPtKaonsX1dSameYMomentaGreaterThan0->Fill(leadingKaon.px()+subLeadingKaon.px() + proton1.X() +  proton2.X());
                    HistPtProtonPtKaonsY1dSameYMomentaGreaterThan0->Fill(leadingKaon.py()+subLeadingKaon.py() +proton1.Y() +  proton2.Y());
                }

                else if (proton1.Y() < 0.0 and proton2.Y() < 0.0)
                {
                    HistPtProtonPtKaonsXSameYMomentaLowerThan0->Fill(leadingKaon.px()+subLeadingKaon.px(), proton1.X() +  proton2.X());
                    HistPtProtonPtKaonsYSameYMomentaLowerThan0->Fill(leadingKaon.py()+subLeadingKaon.py(), proton1.Y() +  proton2.Y());
                    HistPtProtonPtKaonsX1dSameYMomentaLowerThan0->Fill(leadingKaon.px()+subLeadingKaon.px() + proton1.X() +  proton2.X());
                    HistPtProtonPtKaonsY1dSameYMomentaLowerThan0->Fill(leadingKaon.py()+subLeadingKaon.py() +proton1.Y() +  proton2.Y());

                }
            }   

            if  (isPtMissingSmall and isNumberOfTofClusterSmall)
            {
                HistInvMassPiPi2DTofPtMiss->Fill(leadingKaon.m(), subLeadingKaon.m());


                TLorentzVector K01 = leadingKaon.lorentzVector();
                TLorentzVector K02 = subLeadingKaon.lorentzVector();
                TLorentzVector K0K0 = K01 + K02;
                HistMassK0K0->Fill(K0K0.M());
                HistPtK0K0->Fill(K0K0.Pt());
                HistEtaK0K0->Fill(K0K0.Eta());
                    
                HistPtLeadingK0->Fill(leadingKaon.pt());
                HistEtaLeadingK0->Fill(leadingKaon.eta());

                HistPtLeadingK0->Fill(subLeadingKaon.pt());
                HistEtaLeadingK0->Fill(subLeadingKaon.eta());
                HistEtaLeadingSubLeadingDiff->Fill(leadingKaon.eta()-subLeadingKaon.eta());
                
            } 
        }


        }
        posPion.clear();
        negPion.clear();  
     }

    int n = xVec.size();

    double xArr[n];
    double yArr[n];  
    double zArr[n]; 
    double xyArr[n]; 

    for (int i = 0; i< n; i++)
    {
        xArr[i] = xVec[i];
        yArr[i] = yVec[i];
        zArr[i] = zVec[i];
        xyArr[i] = xyVec[i];

    }

    TGraph * graph = new TGraph(xVec.size(), xArr, yArr);
    graph->SetLineWidth(0);
    graph->SetMarkerStyle(1);
    graph->GetXaxis()->SetTitle("x");
    graph->GetYaxis()->SetTitle("y");   
    TGraph * graphXY = new TGraph(xVec.size(), xyArr, zArr);
    graphXY->SetLineWidth(0);
    graphXY->SetMarkerStyle(1);
    graph->GetXaxis()->SetTitle("sqrt(x^2+y^2)");
    graph->GetYaxis()->SetTitle("z");         

    TFile *outfile = TFile::Open(argv[2], "recreate"); 
 
    HistCutFlow->Write(); 
    graph->Write();
    graphXY->Write();
    HistNumPrimaryVertices->Write(); 
    HistPrimaryVertexAbsPosZ->Write();
    HistNumTofMatchedTracks->Write();
    HistTofMatchedTracksCharge->Write();
    HistTofMatchedTracksAbsEta->Write();
    HistTofMatchedTracksPt->Write();
    HistTofMatchedTracksNfit->Write();
    HistTofMatchedTracksNdEdx->Write();

    HistNumPrimaryVerticesPostSelection->Write();
    HistPrimaryVertexAbsPosZPostSelection->Write();
    HistNumTofMatchedTracksPostSelection->Write();
    HistTofMatchedTracksChargePostSelection->Write();
    HistTofMatchedTracksAbsEtaPostSelection->Write();
    HistTofMatchedTracksPtPostSelection->Write();
    HistTofMatchedTracksNfitPostSelection->Write();
    HistTofMatchedTracksNdEdxPostSelection->Write();

    HistPtProtonPtKaonsX->Write();
    HistPtProtonPtKaonsY->Write();
    HistPtProtonPtKaonsX1d->Write();
    HistPtProtonPtKaonsY1d->Write();
    HistPtProtonPtKaonsXOppY->Write();
    HistPtProtonPtKaonsYOppY->Write();
    HistPtProtonPtKaonsX1dOppY->Write();
    HistPtProtonPtKaonsY1dOppY->Write();
    HistPtProtonPtKaonsXSameY->Write();
    HistPtProtonPtKaonsYSameY->Write();
    HistPtProtonPtKaonsX1dSameY->Write();
    HistPtProtonPtKaonsY1dSameY->Write();

    HistPtProtonPtKaonsXSameYMomentaGreaterThan0->Write();
    HistPtProtonPtKaonsYSameYMomentaGreaterThan0->Write();
    HistPtProtonPtKaonsX1dSameYMomentaGreaterThan0->Write();
    HistPtProtonPtKaonsY1dSameYMomentaGreaterThan0->Write();

    HistPtProtonPtKaonsXSameYMomentaLowerThan0->Write();
    HistPtProtonPtKaonsYSameYMomentaLowerThan0->Write();
    HistPtProtonPtKaonsX1dSameYMomentaLowerThan0->Write();
    HistPtProtonPtKaonsY1dSameYMomentaLowerThan0->Write();


    HistInvMassPiPiPeak->Write();
    HistInvMassPiPiPeak2->Write();
    HistInvMassPiPi2D->Write();
    HistNumOfClusters->Write();
    HistPtMiss->Write();

    if (isMC == 1)
    {
        HistKaonPtTruth->Write();
        HistKaonEtaTruth->Write();
        HistKaonVtxRTruth->Write();
        HistKaonVtxZTruth->Write();
        
        HistKaonPtDet->Write();
        HistKaonEtaDet->Write();
        HistKaonVtxRDet->Write();
        HistKaonVtxZDet->Write();
    }

    HistProtonsOppositeSide->Write();
    HistProtonsSameSide->Write();

    HistDcaDaughtersLeadingKaon->Write();
    HistDcaDaughtersSubLeadingKaon->Write();

    HistDistVtxXLeadingKaon->Write();
    HistDistVtxYLeadingKaon->Write();
    HistDistVtxZLeadingKaon->Write();
    HistRLeadingKaon->Write();
    HistClosestDistanceBetweenKaons->Write();
    HistDistVtxXSubLeadingKaon->Write();
    HistDistVtxYSubLeadingKaon->Write();
    HistDistVtxZSubLeadingKaon->Write();
    HistRSubLeadingKaon->Write();

    HistInvMassPiPi2DTofPtMiss->Write();
    HistPtLeadingK0->Write();
    HistEtaLeadingK0->Write();

    HistPtSubLeadingK0->Write();
    HistEtaSubLeadingK0->Write();
    HistEtaLeadingSubLeadingDiff->Write();

    HistMassK0K0->Write();
    HistPtK0K0->Write();
    HistEtaK0K0->Write();
    
    HistZK0K0->Write();
    outfile->Close();

    return 0;
}
