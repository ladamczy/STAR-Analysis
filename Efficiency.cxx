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
#include <TGaxis.h>

#include <TF1.h> 
#include <TF2.h> 
#include <THStack.h> 
#include <TArrow.h> 
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
void FillCutflow(int *i, TH1D *hist );

bool SelectNumberOfVertices(vector <int> vNumber, StUPCEvent *upcEvt, TH1D * HistNumPrimaryVertices, TH1D * HistNumPrimaryVerticesPostSelection);
bool VertexPositionZ(double limitAbsZ, StUPCEvent *upcEvt, TH1D * HistPrimaryVertexAbsPosZ, TH1D * HistPrimaryVertexAbsPosZPostSelection);
void SeparateTracks(StUPCEvent *upcEvt, vector <StUPCTrack const*> &tracksWithTofHit, vector <StUPCTrack const*> &tracksWithoutTofHit, bool isMC);
bool ValidNumberOfTofTracks(vector <int> vNums, vector <StUPCTrack const*> tracksWithTofHit, TH1D *HistNumTofMatchedTracks, TH1D *HistNumTofMatchedTracksPostSelection);
bool AreTofTracksGood(vector <StUPCTrack const*>  tracksWithTofHit, TH1D* HistTofMatchedTracksAbsEta,TH1D* HistTofMatchedTracksAbsEtaPostSelection, TH1D* HistTofMatchedTracksPt,TH1D*HistTofMatchedTracksPtPostSelection,TH1D* HistTofMatchedTracksNfit,TH1D* HistTofMatchedTracksNfitPostSelection,TH1D* HistTofMatchedTracksNdEdx, TH1D*HistTofMatchedTracksNdEdxPostSelection );
bool FindTracks( vector <StUPCTrack const*> &tracksWithTofHit, vector <StUPCTrack const*> &tracksWithoutTofHit, vector <StUPCTrack const*> &tracksWithoutTofHitGoodChargeThreeTof,  TH1D * HistTofMatchedTracksCharge, TH1D * HistTofMatchedTracksChargePostSelection);

void plotHistogram2(const char* inputFileData, const char* histogramName,  double val1, double val2,  double latexX, double latexY,  double xlmin, double ylmin ,string l1, string l2, int b ) ;
void plotHistogram(const char* inputFileData, const char* histogramName,  const char* histogramName2,  const char* histogramName3, double val1, double val2,  double latexX, double latexY,  double xlmin, double ylmin ,string l1, string l2,int , double XMIN, double XMAX) ;



double etaCutVal = 0.9;
double ptCutVal = 0.2;
double NhistFitCut = 20;

double duncErr(double T, double N,double t, double n)
{
    double a = sqrt(pow(N/pow((T+N),2)*t ,2) + pow(T/pow((T+N),2)*n,2) );
    return a;
}

void WriteCSV_file(string filename, vector <vector<double>> data)
{
 
    ofstream file(filename);

if (file.is_open())
{
  
    for (int i = 0; i < data.size(); i++)
    {
 

        for (int j = 0; j < 3; j++)
        {
            file << data[i][j];
            if (j!=2)
            {
                file << ",";
            }
        }
        file << endl;
    }



    file.close();
}



}

void SaveData( StUPCEvent *upcEvt, Long64_t &num, TLorentzVector &proton1, TLorentzVector &proton2);

    TH1D* HistPhi = new TH1D("HistPhi", ";;events", 101, -M_PI, M_PI);
    TH1D* HistPhi4 = new TH1D("HistPhi4", ";;events",101, -M_PI, M_PI);

    TH1D* HistEta = new TH1D("HistEta", ";;events",  101,  -3,3);
    TH1D* HistEta4 = new TH1D("HistEta4", ";;events", 101,  -3,3);


    TH2D* HistPhiEta = new TH2D("HistPhiEta", ";;events",  101, -M_PI, M_PI, 101,-3,3);
    TH2D* HistPhiEta4 = new TH2D("HistPhiEta4", ";;events", 101, -M_PI, M_PI, 101,-3,3);
    


int main(int argc, char** argv)  
{

    int d4 = 0;
    int d3 = 0;
    int d2 =0;

    TH2D* HistPtPion = new TH2D("HistPtPion", ";;", 40, 0.2, 1, 40, 0.2,1);
    TH2D* HistEtaPion = new TH2D("HistEtaPion", ";;",  28, -0.7, 0.7, 28, -0.7, 0.7);
    TH2D* HistMassK0K0 =  new TH2D("HistMassK0K0", ";;",  30, 0.44, 0.56, 30, 0.44,0.56);




    TH1D* HistCutFlow = new TH1D("HistCutFlow", ";;count", 16, -0.5, 15.5);
    TH1D* HistDistR = new TH1D("HistDistR", ";R [cm];events", 55, 0, 5.5);
    TH1D* HistDistRSum = new TH1D("HistDistRSum", "; #Sigma_{#pi} R [cm];events", 55, 0, 5.5);


    TH1D* HistNumTofTracks4 =  new TH1D("HistNumTofTracks_4", ";N_{TOF}^{tracks};events",16,-0.5,15.5);
    TH1D* HistMassLeadingK04 =  new TH1D("HistMassLeadingK0_4", " ;m_{#pi^{+}#pi^{-}} [GeV]; N_{K}/N_{m}",30 ,0.44, 0.56);
    TH1D* HistMassSubLeadingK04 =  new TH1D("HistMassSubLeadingK0_4", " ;m_{#pi^{+}#pi^{-}} [GeV]; events/N_{m}",30 ,0.44, 0.56);
    TH1D* HistPtMiss4 =  new TH1D("HistPtMiss_4", ";p_{T}^{miss} [GeV]; events/N_{m}",10, 0, 0.5);
    TH1D* HistNumOfClusters4 =  new TH1D("HistNumOfClusters_4", " ; N_{TOF}^{clusers}; events/N_{m}", 11, -0.5, 10.5);
    TH1D* HistDcaDaughtersLeadingKaon4 =  new TH1D("HistDcaDaughtersLeadingKaon_4", "; DCA_{daughters} [cm];  N_{K}/N_{m}",10, 0, 5);
    TH1D* HistDcaDaughtersSubLeadingKaon4 =  new TH1D("HistDcaDaughtersSubLeadingKaon_4", "; DCA_{K0} [cm];  N_{K}/N_{m}",10, 0, 5);
    TH1D* HistDcaBeamlineLeadingKaon4 =  new TH1D("HistDcaBeamlineLeadingKaon_4", " ;DCA_{beamline} [cm] ;  N_{K}/N_{m}",10, 0, 5);
    TH1D* HistDcaBeamlineSubLeadingKaon4 =  new TH1D("HistDcaBeamlineSubLeadingKaon_4", " ;DCA^{beamline} [cm] ;  N_{K}/N_{m}",10, 0, 5);
    TH1D* HistPointingAngleLeadingKaon4 =  new TH1D("HistPointingAngleLeadingKaon_4", " ;cos(#alpha_{p}) ;  N_{K}/N_{m}",20,0.9,1);
    TH1D* HistPointingAngleSubLeadingKaon4 =  new TH1D("HistPointingAngleSubLeadingKaon_4", " ;DCA_{beamline} [cm] ;  N_{K}/N_{m}",100,0.9,1);
    TH1D* HistHypoLengthLeadingKaon4 =  new TH1D("HistHypoLengthLeadingKaon_4", " ;l_{decay} [cm] ;  N_{K}/N_{m}",9,0,9);
    TH1D* HistHypoLengthSubLeadingKaon4 =  new TH1D("HistHypoLengthSubLeadingKaon_4", " ;DCA^{beamline} [cm] ;  N_{K}/N_{m}", 50,0, 10);
    TH1D* HistKsi14 =  new TH1D("HistKsi1_4", " ; #xi^{E}; events/N_{m}",30, -0.015, 0.015);
    TH1D* HistKsi24 =  new TH1D("HistKsi2_4", " ; #xi^{W}; events/N_{m}",30, -0.015, 0.015);
    TH1D* HistDiffProtonX4 =  new  TH1D("HistDiffProtonX_4", " ; p_{x}^{E} + p_{x}^{W} [GeV] ; events/N_{m}",20,-1,1);
    TH1D* HistDiffProtonY4 =  new  TH1D("HistDiffProtonY_4", " ; p_{y}^{E} + p_{y}^{W} [GeV] ; events/N_{m}", 40, -2,2);
    TH1D* HistKsiCorrelation4 =  new TH1D("HistKsiCorrelation_4", " ;m_{KK}/#sqrt{s} - #sqrt{#xi^{E}#xi^{W}} ; events/N_{m}", 10, -0.005,0.005);
    TH1D* HistEtaCorrelation4 =  new TH1D("HistEtaCorrelation_4", " ;y_{KK}-1/2ln(#xi^{E}/#xi^{W}) ; events/N_{m}",30, -1.5,1.5);
    TH1D* HistCorr14= new TH1D("HistCorr1_4", " ;m_{KK}e^{-y_{KK}} - #xi^{E} ; events/N_{m}",20, -0.01, 0.01);
    TH1D* HistCorr24= new TH1D("HistCorr2_4", " ;m_{KK}e^{y_{KK}} - #xi^{W} ; events/N_{m}",20, -0.01, 0.01);
    
    TH1D* HistNumTofTracks3 =  new TH1D("HistNumTofTracks_3", ";N_{tracks}^{TOF};events/N_{m}",11, -0.5, 10.5);
    TH1D* HistMassLeadingK03 =  new TH1D("HistMassLeadingK0_3", " ;m_{#pi^{+}#pi^{-}} [GeV]; events/N_{m}",30 ,0.44, 0.56);
    TH1D* HistMassSubLeadingK03 =  new TH1D("HistMassSubLeadingK0_3", " ;m_{#pi^{+}#pi^{-}} [GeV]; events/N_{m}",30 ,0.44, 0.56);
    TH1D* HistPtMiss3 =  new TH1D("HistPtMiss_3", ";p_{T}^{miss} [GeV]; events/N_{m}", 10, 0, 0.5);
    TH1D* HistNumOfClusters3 =  new TH1D("HistNumOfClusters_3", "; # TOF clusters; events/N_{m}",13, -0.5,12.5);
    TH1D* HistDcaDaughtersLeadingKaon3 =  new TH1D("HistDcaDaughtersLeadingKaon_3", "; DCA_{daughters} [cm]; events/N_{m}",10, 0, 5);
    TH1D* HistDcaDaughtersSubLeadingKaon3 =  new TH1D("HistDcaDaughtersSubLeadingKaon_3", "; DCA_{K0} [cm]; events/N_{m}",10, 0, 5);
    TH1D* HistDcaBeamlineLeadingKaon3 =  new TH1D("HistDcaBeamlineLeadingKaon_3", " ;DCA_{beamline} [cm] ; events/N_{m}",10, 0, 5);
    TH1D* HistDcaBeamlineSubLeadingKaon3 =  new TH1D("HistDcaBeamlineSubLeadingKaon_3", " ;DCA^{beamline} [cm] ; events/N_{m}",10, 0, 5);
    TH1D* HistPointingAngleLeadingKaon3 =  new TH1D("HistPointingAngleLeadingKaon_3", " ;cos(#alpha_{p}) ; events/N_{m}", 20, 0.9, 1.0);
    TH1D* HistPointingAngleSubLeadingKaon3 =  new TH1D("HistPointingAngleSubLeadingKaon_3", " ;DCA_{beamline} [cm] ; events/N_{m}", 100, -1, 1.0);
    TH1D* HistHypoLengthLeadingKaon3 =  new TH1D("HistHypoLengthLeadingKaon_3", " ;l_{decay} [cm] ; events/N_{m}",9,0,9);
    TH1D* HistHypoLengthSubLeadingKaon3 =  new TH1D("HistHypoLengthSubLeadingKaon_3", " ;DCA^{beamline} [cm] ; events/N_{m}", 50,0, 10);
    TH1D* HistKsi13 =  new TH1D("HistKsi1_3", " ; #xi; events/N_{m}",30, -0.015, 0.015);
    TH1D* HistKsi23 =  new TH1D("HistKsi2_3", " ; #xi; events/N_{m}",30, -0.015, 0.015);
    TH1D* HistDiffProtonX3 =  new  TH1D("HistDiffProtonX_3", " ; diff p_{x}^{proton} [GeV] ; events/N_{m}",20, -1,1);
    TH1D* HistDiffProtonY3 =  new  TH1D("HistDiffProtonY_3", " ; diff p_{y}^{proton} [GeV] ; events/N_{m}",40, -2,2);
    TH1D* HistKsiCorrelation3 =  new TH1D("HistKsiCorrelation_3", " ;m_{KK}/#sqrt{s} - #sqrt{#xi^{E}#xi^{W}} ; events/N_{m}", 10, -0.005,0.005);
    TH1D* HistEtaCorrelation3 =  new TH1D("HistEtaCorrelation_3", " ;#eta_{K0K0} - 1/2ln(#xi^{E}/#xi^{W}} ; events/N_{m}",  30, -1.5,1.5);
    TH1D* HistCorr13= new TH1D("HistCorr1_3", " ;m_{KK}e^{#pm y} - #xi_{W/E} ; events/N_{m}",20, -0.01, 0.01);
    TH1D* HistCorr23= new TH1D("HistCorr2_3", " ;m_{KK}/#sqrt{s} - #sqrt{#xi_1#xi_2} ; events/N_{m}",20, -0.01, 0.01);
    
    TH1D* HistNumTofTracks2 =  new TH1D("HistNumTofTracks_2", ";N_{tracks}^{TOF};events/N_{m}",16,-0.5,15.5);
    TH1D* HistMassLeadingK02 =  new TH1D("HistMassLeadingK0_2", " ;m_{#pi^{+}#pi^{-}} [GeV]; events/N_{m}",30 ,0.44, 0.56);
    TH1D* HistMassSubLeadingK02 =  new TH1D("HistMassSubLeadingK0_2", " ;m_{#pi^{+}#pi^{-}} [GeV]; events/N_{m}",30 ,0.44, 0.56);
    TH1D* HistPtMiss2 =  new TH1D("HistPtMiss_2", ";p_{T}^{miss} [GeV]; events/N_{m}",10, 0, 0.5);
    TH1D* HistNumOfClusters2 =  new TH1D("HistNumOfClusters_2", "; # TOF clusters; events/N_{m}",11, -0.5, 10.5);
    TH1D* HistDcaDaughtersLeadingKaon2 =  new TH1D("HistDcaDaughtersLeadingKaon_2", "; DCA_{daughters} [cm]; events/N_{m}",10, 0, 5);
    TH1D* HistDcaDaughtersSubLeadingKaon2 =  new TH1D("HistDcaDaughtersSubLeadingKaon_2", "; DCA_{K0} [cm]; events/N_{m}",10, 0, 5);
    TH1D* HistDcaBeamlineLeadingKaon2 =  new TH1D("HistDcaBeamlineLeadingKaon_2", " ;DCA_{beamline} [cm] ; events/N_{m}",10, 0, 5);
    TH1D* HistDcaBeamlineSubLeadingKaon2 =  new TH1D("HistDcaBeamlineSubLeadingKaon_2", " ;DCA^{beamline} [cm] ; events/N_{m}",10, 0, 5);
    TH1D* HistPointingAngleLeadingKaon2 =  new TH1D("HistPointingAngleLeadingKaon_2", " ;cos(#alpha_{p}) ; events/N_{m}", 20,0.9,1);
    TH1D* HistPointingAngleSubLeadingKaon2 =  new TH1D("HistPointingAngleSubLeadingKaon_2", " ;DCA_{beamline} [cm] ; events/N_{m}", 100, -1, 1.0);
    TH1D* HistHypoLengthLeadingKaon2 =  new TH1D("HistHypoLengthLeadingKaon_2", " ;l_{decay} [cm] ; events/N_{m}", 9,0,9);;
    TH1D* HistHypoLengthSubLeadingKaon2 =  new TH1D("HistHypoLengthSubLeadingKaon_2", " ;DCA^{beamline} [cm] ; events/N_{m}", 50,0, 10);
    TH1D* HistKsi12 =  new TH1D("HistKsi1_2", " ; #xi^{E}; events/N_{m}",20, -0.01, 0.01);
    TH1D* HistKsi22 =  new TH1D("HistKsi2_2", " ; #xi^{W}; events/N_{m}", 20, -0.01, 0.01);

    TH1D* HistDiffProtonX2 =  new  TH1D("HistDiffProtonX_2", " ; diff p_{x}^{proton} [GeV] ; events/N_{m}",20, 1,1);
    TH1D* HistDiffProtonY2 =  new  TH1D("HistDiffProtonY_2", " ; diff p_{y}^{proton} [GeV] ; events/N_{m}",40, -2,2);
    TH1D* HistKsiCorrelation2 =  new TH1D("HistKsiCorrelation_2", " ;m_{KK}/#sqrt{s} - #sqrt{#xi^{E}#xi^{W}} ; events/N_{m}", 10, -0.005,0.005);
    TH1D* HistEtaCorrelation2 =  new TH1D("HistEtaCorrelation_2", " ;#eta_{K0K0} - 1/2ln(#xi^{E}/#xi^{W}} ; events/N_{m}",  30, -1.5,1.5);
    TH1D* HistCorr12= new TH1D("HistCorr1_2", " ;m_{KK}e^{#pm y} - #xi_{W/E} ; events/N_{m}", 30, -0.015, 0.015);
    TH1D* HistCorr22= new TH1D("HistCorr2_2", " ;m_{KK}/#sqrt{s} - #sqrt{#xi_1#xi_2} ; events/N_{m}",  30, -0.015, 0.015);
    TH1D* histZDiff4 = new TH1D("histZDiff4", " ;vtx_{z}^{leading}-vtx_{z}^{subleading} [cm] ;events/N_{m}", 12, -30, 30);
    TH1D* histZMean4 = new TH1D("histZMean4", "  ;(vtx_{z}^{leading}+vtx_{z}^{subleading})/2 [cm]; events/N_{m}",  20, -100,100);
    TH1D* histZDiff3 = new TH1D("histZDiff3", " ;vtx_{z}^{leading}-vtx_{z}^{subleading} [cm] ;events/N_{m}", 12 , -30,30);
    TH1D* histZMean3 = new TH1D("histZMean3", "  ;(vtx_{z}^{leading}+vtx_{z}^{subleading} [cm])/2; events/N_{m}",  20, -100,100);
    TH1D* histZDiff2 = new TH1D("histZDiff2", " ;vtx_{z}^{leading}-vtx_{z}^{subleading} [cm] ;events/N_{m}", 12, -30,30);
    TH1D* histZMean2 = new TH1D("histZMean2", "  ;(vtx_{z}^{leading}+vtx_{z}^{subleading} [cm])/2; events/N_{m}", 20, -100,100);
    TH1D* histCosThetaStar4= new TH1D("histCosThetaStar4", "  ;cos(#theta*); N_{K}/N_{m}",  20, -1.0, 1.0);
    TH1D* histCosThetaStar3= new TH1D("histCosThetaStar3", "  ;cos(#theta*); N_{K}/N_{m}",  20, -1.0, 1.0);
    TH1D* histCosThetaStar2= new TH1D("histCosThetaStar2", "  ;cos(#theta*); N_{K}/N_{m}",  20, -1.0, 1.0);
    vector <TH1D *> HistNfitPionWithTot;
    vector <TH1D *> HistNfitPionWithoutTot;
      
    for (int i = 0; i < 3; i++)
    {
        TH1D* histNfitPionWithTof = new TH1D(Form("histNfitPionWithTof%d", i), Form("histNfitPionWithTof _%d", i), 101, -0.5, 100.5);
        TH1D* histNfitPionWithoutTot = new TH1D(Form("HistNfitPionWithoutTot%d", i), Form("HistNfitPionWithoutTot _%d", i), 101, -0.5, 100.5);
        HistNfitPionWithTot.push_back(histNfitPionWithTof);
        HistNfitPionWithoutTot.push_back(histNfitPionWithoutTot);
    }

    long int sumTofTracks[5] = {0,0,0,0,0};
    long int sumMassWindow[5] = {0,0,0,0,0};
    long int sumSmallPt[5] = {0,0,0,0,0};
    long int sumNtofCluster[5] = {0,0,0,0,0};
    long int sumDecayCos[5] = {0,0,0,0,0};
    long int sumDcaBeamline[5] = {0,0,0,0,0};
    long int sumDcaDaughters[5] = {0,0,0,0,0};
    long int sumAntiEllastic[5] = {0,0,0,0,0};
    long int sumKsiCorr[5] = {0,0,0,0,0};
    long int sumEtaCorr[5] = {0,0,0,0,0};
    long int sumSingleCorr[5] = {0,0,0,0,0};

    long int effZmean[5] = {0,0,0,0,0};
    long int effZdiff[5] = {0,0,0,0,0};
    long int effCosThetaStar[5] = {0,0,0,0,0};


    long int sumEvt = 0;

    float massPion = 0.13957061;
    double massKaon =  497.611/1000.0;
    double beamEnergy = 254.867;
    double massProton = 938.272/1000.0;
    int evtIn;
    Long64_t numsa = 0;

    ifstream inputFilePathList(argv[1]);

    TChain *chain = new TChain("mUPCTree"); 
    string inputFileName;

    while (std::getline(inputFilePathList, inputFileName))
    {
        chain->Add(inputFileName.c_str());     
    }
    inputFilePathList.close();


    bool isMC = 0;
    static StUPCEvent *upcEvt = 0x0;
 
    static StRPEvent *rpEvt = 0x0;
    static StRPEvent  *correctedRpEvent = 0x0;
    chain->SetBranchAddress("mRPEvent", &rpEvt);

    LoadOffsetFile("../share/OffSetsCorrectionsRun17.list", mCorrection);

    TTree* mUPCTree = new TTree("mUPCTree", "mUPCTree");
    mUPCTree->Branch("mUPCEvent", &upcEvt);
    mUPCTree->Branch("correctedRpEvent", &correctedRpEvent);



    int sumPassed;
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

    const vector <int> triggerID = {570209, 570219, 570229, 570701, 570702, 570703, 570704, 570705, 570709, 570711, 570712, 570719, 590701, 590703, 590705, 590708, 590709};
    const vector <int> triggerCEP = {570701, 570705, 570711};

    int flags = 0;
    for (Long64_t i = 0; i < chain->GetEntries(); ++i) 
    {
  
        if (i%1000000 == 0)  
        {
            cout << i << "/" <<  chain->GetEntries() << endl;
        }

        chain->GetEntry(i);
      
        TLorentzVector proton1, proton2;
        FindProtons(isMC, correctedRpEvent, upcEvt, proton1, proton2);  
        double protonSignsY = 1.0;
        TParticle* particle;   
        vector <TParticle *> vPions;       
        vector <TParticle *>  vPionsAll;
        vector <TParticle *>  vAntiPionsAll;
        vector <TParticle *>  vProtonsAll;
        
        vector <TParticle *>  vAntiProtonsAll;


  
            // extract all pions
            for (int i = 0; i <  upcEvt->getNumberOfMCParticles(); i++)
            {
                particle = upcEvt->getMCParticle(i);
                if (particle->GetPDG()->PdgCode() == 211 or particle->GetPDG()->PdgCode() == -211)
                {
                    vPionsAll.push_back(particle);
                }
            }

            TLorentzVector pioni, pionj;
            TLorentzVector productionVertexi, productionVertexj;
            
            // extract pion pairs coming from the same mother 
            for (int i = 0; i < vPionsAll.size(); i++)
            {
                for (int j = 0; j < vPionsAll.size(); j++)
                {
                    if ( (i!=j) and (vPionsAll[i]->GetFirstMother() ==  vPionsAll[j]->GetFirstMother())  and (vPionsAll[i]->GetPDG()->PdgCode() != vPionsAll[j]->GetPDG()->PdgCode()) )
                    {
                        auto it =find(vPions.begin(), vPions.end(),vPionsAll[i] );

                        if ( it == vPions.end())
                        {
                            pioni.SetPxPyPzE(vPionsAll[i]->Px(), vPionsAll[i]->Py(), vPionsAll[i]->Pz(), vPionsAll[i]->Energy() );
                            pionj.SetPxPyPzE(vPionsAll[j]->Px(), vPionsAll[j]->Py(), vPionsAll[j]->Pz(), vPionsAll[j]->Energy() );
                            vPionsAll[i]->ProductionVertex(productionVertexi);
                            vPionsAll[j]->ProductionVertex(productionVertexi);    
                            double truthVertexZ1 = productionVertexi.Z();
                            double truthVertexZ2 = productionVertexj.Z();

                            if(  pioni.Pt() > 0.2 and pionj.Pt() > 0.2 and abs(pioni.Eta())<0.9  and abs(pionj.Eta())<0.9)// and abs(truthVertexZ1) < 80.0)// and abs(truthVertexZ2) < 80.0 )
                            {
                                vPions.push_back(vPionsAll[i]);
                                vPions.push_back(vPionsAll[j]);
                            }
                        }
                    }
                }
                if (vPions.size() == 4)
                {
                    break;
                }
            }
        

        // CONTINUE
        if (vPions.size() != 4) // ok - uwzględniam przypadki z >= 5 kaonami
        {
            continue;
        }

        vector <TParticle*> vPionsNeg;
        vector <TParticle*> vPionsPos;



        for (int i = 0; i < 4; i++)
        {
            if ((vPions[i]->GetPDG()->PdgCode()==211))
            {
                vPionsPos.push_back(vPions[i]);
            }

            else
            {
                vPionsNeg.push_back(vPions[i]);  
            }
        }

        // move all TPC tracks to the vector
        vector <StUPCTrack const *> all; 
        for (int i = 0; i < upcEvt->getNumberOfTracks(); i++)
        {
            StUPCTrack *track = upcEvt->getTrack(i);
            all.push_back(track);
        }

        // CONTINUE event with less than 4 tpc tracks is removed
        if (all.size() < 4)
        {
            all.clear();
            continue;
        }

        double r;
        vector <StUPCTrack const *> matchedTPCpions;
        vector <double> smallestR;

        for (int j = 0; j < 4; j++)
        {
            vector <double> distr;
            for (int i = 0; i < all.size(); i++)
            {
                TLorentzVector pion, track;
                pion.SetPxPyPzE(vPions[j]->Px(), vPions[j]->Py(), vPions[j]->Pz(), vPions[j]->Energy());
                all[i]->getLorentzVector(track, massPion);
                r = pion.DeltaR(track);
                distr.push_back(r);
            }

            vector <double >sorted_distr = distr;
            sort(sorted_distr.begin(), sorted_distr.end() );

            for (int i = 0; i < sorted_distr.size(); i++)
            {

                auto it = std::find(distr.begin(), distr.end(), sorted_distr[i]);
                int index = std::distance(distr.begin(), it);
                if (all[index]->getCharge() * vPions[j]->GetPDG()->PdgCode() > 0.0   and all[index]->getPt()>0.2 and abs(all[index]->getEta())< 0.9 and all[index]->getNhitsFit() >=20  )
                {       
                    matchedTPCpions.push_back(all[index]);
                    smallestR.push_back(distr[index]);
                    all.erase(all.begin()+index);
                    distr.clear();
                    break;
                }

                else
                {
                    continue;
                }
            }
        }

        if ( matchedTPCpions.size() != 4)
        {
            continue;
        }
        
        for (int i = 0; i < matchedTPCpions.size(); i++)
        {
 
                if (vPions[0]->GetFirstMother() != vPions[1]->GetFirstMother()) {cout << "no" << endl;}
                TLorentzVector pion, track;
                pion.SetPxPyPzE(vPions[i]->Px(), vPions[i]->Py(), vPions[i]->Pz(), vPions[i]->Energy());
                double ptTrue = pion.Pt();
                double ptDet = matchedTPCpions[i]->getPt();
                double etaTrue = pion.Eta();
                double etaDet = matchedTPCpions[i]->getEta();

                HistPtPion->Fill(ptTrue,ptDet);
                HistEtaPion->Fill(etaTrue,etaDet);
        }

                double beamPar[4] = {0.0, 0.0, 0.0, 0.0};
     
        TVector3 const tryVec(0,0,0);
        StUPCV0 effKaon1(matchedTPCpions[0], matchedTPCpions[1], massPion, massPion, 1, 1, tryVec, beamPar,upcEvt->getMagneticField(), true);
        StUPCV0 effKaon2(matchedTPCpions[2], matchedTPCpions[3], massPion, massPion, 1, 1, tryVec,beamPar, upcEvt->getMagneticField(), true);

        TLorentzVector pion0, pion1, pion2, pion3;
        pion0.SetPxPyPzE(vPions[0]->Px(), vPions[0]->Py(), vPions[0]->Pz(), vPions[0]->Energy());
        pion1.SetPxPyPzE(vPions[1]->Px(), vPions[1]->Py(), vPions[1]->Pz(), vPions[1]->Energy());       
        pion2.SetPxPyPzE(vPions[2]->Px(), vPions[2]->Py(), vPions[2]->Pz(), vPions[2]->Energy());
        pion3.SetPxPyPzE(vPions[3]->Px(), vPions[3]->Py(), vPions[3]->Pz(), vPions[3]->Energy());    
        TLorentzVector kaon1 = pion0 + pion1;
        TLorentzVector kaon2 =pion2 + pion3;

        HistMassK0K0->Fill(kaon1.M(),effKaon1.m() );
        HistMassK0K0->Fill(kaon2.M(),effKaon2.m() );

        bool smallR = true;

        for (int i = 0; i < smallestR.size(); i++)
        {
            if (smallestR[i] >= 0.15)
            {
                smallR = false;
            }
        }

        // CONTINUE not matched
        if (smallR == false)
        {
            continue;
        }  

        Int_t BBCLEast = upcEvt->getBBCLargeEast();
        Int_t BBCLWest = upcEvt->getBBCLargeWest();

        if (BBCLWest > 52 or BBCLEast > 52.0)
        {
            continue;
        }

        sumEvt+=1;
    
        //N_{m}umber of TOF trakcs
        int Ntof = 0; 
        for (int i = 0; i < upcEvt->getNumberOfTracks(); i++)
        {   
            if (upcEvt->getTrack(i)->getFlag(StUPCTrack::kTof))
            {
                Ntof+=1;
            }          
        }

        //selection cut - mass - selecting kaon by pT
        vector <StUPCTrack const *> posPion;
        vector <StUPCTrack const *> negPion;
        vector <StUPCTrack const*> vPosNegPionLeadingKaon;
        vector <StUPCTrack const*> vPosNegPionSubLeadingKaon;


        double mass1,mass2,mass3,mass4;

            for (unsigned long int i = 0; i < matchedTPCpions.size(); i++)
            {
                if( matchedTPCpions[i]->getCharge() == 1)
                {
                    posPion.push_back(matchedTPCpions[i]);
                }
                else if (matchedTPCpions[i]->getCharge() == -1)
                {
                    negPion.push_back(matchedTPCpions[i]);
                }
            }


            StUPCV0 *lVecKaonComb1a = new StUPCV0(posPion[0], negPion[0], massPion, massPion, 1, 1, tryVec, beamPar, upcEvt->getMagneticField(), true);
            StUPCV0 *lVecKaonComb1b = new StUPCV0(posPion[1], negPion[1], massPion, massPion, 1, 1, tryVec, beamPar,upcEvt->getMagneticField(), true);
            StUPCV0 *lVecKaonComb2a = new StUPCV0(posPion[0], negPion[1], massPion, massPion, 1, 1, tryVec, beamPar,upcEvt->getMagneticField(), true);
            StUPCV0 *lVecKaonComb2b = new StUPCV0(posPion[1], negPion[0], massPion, massPion, 1, 1, tryVec, beamPar,upcEvt->getMagneticField(), true);
            
        
            double distSumSquaredComb1 = sqrt(pow((lVecKaonComb1a->m()-massKaon), 2) + pow((lVecKaonComb1b->m()-massKaon), 2));
            double distSumSquaredComb2 = sqrt(pow((lVecKaonComb2a->m()-massKaon), 2) + pow((lVecKaonComb2b->m()-massKaon), 2));

    
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
             mass1 = mass2 = mass3 =mass4 = massPion;
        


        StUPCV0 leadingKaon(vPosNegPionLeadingKaon[0], vPosNegPionLeadingKaon[1], mass1, mass2, 1, 1, tryVec, beamPar,upcEvt->getMagneticField(), true);
        StUPCV0 subLeadingKaon(vPosNegPionSubLeadingKaon[0], vPosNegPionSubLeadingKaon[1], mass3, mass4, 1, 1, tryVec,beamPar, upcEvt->getMagneticField(), true);

        // MISSING PTM
        double pXmiss = leadingKaon.px() + subLeadingKaon.px() + proton1.X() + proton2.X();
        double pYmiss = leadingKaon.py() + subLeadingKaon.py()  + proton1.Y() + proton2.Y();
        double pTmiss = sqrt(pow(pXmiss,2) + pow(pYmiss,2));



        // NTOF CLUSTER
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



        // protons ksi for MC
        double ksi1 =(255 - proton1.E())/255.;
        double ksi2 =(255 - proton2.E())/255.;
    
        // protons transverse momenta
        double sumProtonMomentumX = proton1.X() + proton2.X();
        double sumProtonMomentumY = proton1.Y() + proton2.Y();
        
 

        TLorentzVector K01 = leadingKaon.lorentzVector();
        TLorentzVector K02 = subLeadingKaon.lorentzVector();
        TLorentzVector K0K0 = K01 + K02;

        double ksi1New = ksi1;
        double ksi2New = ksi2;
                           
        if (ksi1 < 0.006)
        {
            ksi1New = 0.003;
        }

        if (ksi2 < 0.006)
        {
            ksi2New = 0.003;
        }
        
        double ksiCorrelation = K0K0.M()/510.0 - sqrt(ksi1New*ksi2New);
        double etaCorrelation = K0K0.Rapidity()  - 0.5*log(ksi1New/ksi2New);
 
        double newCorr1 = K0K0.M()/510.0*exp(-1*K0K0.Rapidity() )-ksi1;
        double newCorr2 = K0K0.M()/510.0*exp(K0K0.Rapidity() )-ksi2;

        int NtofMatched = 0;
        TVector3 vertexLeadingKaon = leadingKaon.decayVertex();
        TVector3 vertexSubLeadingKaon = subLeadingKaon.decayVertex();

        double zdiff = vertexLeadingKaon.Z()-vertexSubLeadingKaon.Z();
        double zmean = (vertexLeadingKaon.Z()+vertexSubLeadingKaon.Z())/2.0;

        double cosThetaStarLeading = leadingKaon.cosThetaStar();
        double cosThetaStarSubLeading = subLeadingKaon.cosThetaStar();

        vector <StUPCTrack const *> tof;
        vector <StUPCTrack const *> notof;
              
        for (int i = 0; i < matchedTPCpions.size(); i++)
        {
            if (matchedTPCpions[i]->getFlag(StUPCTrack::kTof))
            {
                NtofMatched+=1;
                tof.push_back(matchedTPCpions[i]);
            }
            else
            {
                notof.push_back(matchedTPCpions[i]);
            }
        }

        if (NtofMatched==4)
        {
            d4+=1;
            for (int i = 0; i < tof.size(); i++)
            {
                HistNfitPionWithTot[2]->Fill(tof[i]->getNhitsFit());
            }            
        }


        if (NtofMatched==3)
        
        
        {
             d3+=1;
            for (int i = 0; i < tof.size(); i++)
            {
                HistNfitPionWithTot[1]->Fill(tof[i]->getNhitsFit());
            }            
            for (int i = 0; i < notof.size(); i++)
            {
                HistNfitPionWithoutTot[1]->Fill(notof[i]->getNhitsFit());
            }               
        }


        if (NtofMatched==2)
        { d2+=1;
            for (int i = 0; i < tof.size(); i++)
            {
                HistNfitPionWithTot[0]->Fill(tof[i]->getNhitsFit());
            }            
            for (int i = 0; i < notof.size(); i++)
            {
                HistNfitPionWithoutTot[0]->Fill(notof[i]->getNhitsFit());
            }               
        }




        if (NtofMatched==4)
        {
      
            HistNumTofTracks4->Fill(Ntof);
            HistMassLeadingK04->Fill(leadingKaon.m());
            HistMassLeadingK04->Fill(subLeadingKaon.m());
            HistPtMiss4->Fill(pTmiss);
            HistNumOfClusters4->Fill(totalCluster);
            HistDcaBeamlineLeadingKaon4->Fill(leadingKaon.DCABeamLine());
            HistDcaBeamlineLeadingKaon4->Fill(subLeadingKaon.DCABeamLine());
            HistDcaDaughtersLeadingKaon4->Fill(leadingKaon.dcaDaughters());
            HistDcaDaughtersLeadingKaon4->Fill(subLeadingKaon.DCABeamLine());
            HistPointingAngleLeadingKaon4->Fill(leadingKaon.pointingAngleHypo());
            HistPointingAngleLeadingKaon4->Fill(subLeadingKaon.pointingAngleHypo());
            HistHypoLengthLeadingKaon4->Fill(leadingKaon.decayLengthHypo());
            HistHypoLengthLeadingKaon4->Fill(subLeadingKaon.decayLengthHypo());
            HistKsi14->Fill(ksi1);
            HistKsi24->Fill(ksi2);
            HistDiffProtonX4->Fill(sumProtonMomentumX);
            HistDiffProtonY4->Fill(sumProtonMomentumY);
            HistKsiCorrelation4->Fill(ksiCorrelation);
            HistEtaCorrelation4->Fill(etaCorrelation);
            HistCorr14->Fill(newCorr1);
            HistCorr24->Fill(newCorr2);
            histZDiff4->Fill(zdiff);
            histZMean4->Fill(zmean);

            histCosThetaStar4->Fill(cosThetaStarLeading);
            histCosThetaStar4->Fill(cosThetaStarSubLeading);    
        }


        if (NtofMatched==3)
        {
            HistNumTofTracks3->Fill(Ntof);
            HistMassLeadingK03->Fill(leadingKaon.m());
            HistMassLeadingK03->Fill(subLeadingKaon.m());
            HistPtMiss3->Fill(pTmiss);
            HistNumOfClusters3->Fill(totalCluster);
            HistDcaBeamlineLeadingKaon3->Fill(leadingKaon.DCABeamLine());
            HistDcaBeamlineLeadingKaon3->Fill(subLeadingKaon.DCABeamLine());
            // DCA daughters
            HistDcaDaughtersLeadingKaon3->Fill(leadingKaon.dcaDaughters());
            HistDcaDaughtersLeadingKaon3->Fill(subLeadingKaon.DCABeamLine());
            //Decay length and cosine of the pointing angle
            HistPointingAngleLeadingKaon3->Fill(leadingKaon.pointingAngleHypo());
            HistPointingAngleLeadingKaon3->Fill(subLeadingKaon.pointingAngleHypo());
            HistHypoLengthLeadingKaon3->Fill(leadingKaon.decayLengthHypo());
            HistHypoLengthLeadingKaon3->Fill(subLeadingKaon.decayLengthHypo());
            HistKsi13->Fill(ksi1);
            HistKsi23->Fill(ksi2);
            HistDiffProtonX3->Fill(sumProtonMomentumX);
            HistDiffProtonY3->Fill(sumProtonMomentumY);
            HistKsiCorrelation3->Fill(ksiCorrelation);
            HistEtaCorrelation3->Fill(etaCorrelation);
            HistCorr13->Fill(newCorr1);
            HistCorr23->Fill(newCorr2);
            histZDiff3->Fill(zdiff);
            histZMean3->Fill(zmean);

            histCosThetaStar3->Fill(cosThetaStarLeading);
            histCosThetaStar3->Fill(cosThetaStarSubLeading);   
        }


        if (NtofMatched==2)
        {
          
            HistNumTofTracks2->Fill(Ntof);
            HistMassLeadingK02->Fill(leadingKaon.m());
            HistMassLeadingK02->Fill(subLeadingKaon.m());
            HistPtMiss2->Fill(pTmiss);
            HistNumOfClusters2->Fill(totalCluster);
            HistDcaBeamlineLeadingKaon2->Fill(leadingKaon.DCABeamLine());
            HistDcaBeamlineLeadingKaon2->Fill(subLeadingKaon.DCABeamLine());

            HistDcaDaughtersLeadingKaon2->Fill(leadingKaon.dcaDaughters());
            HistDcaDaughtersLeadingKaon2->Fill(subLeadingKaon.DCABeamLine());

            HistPointingAngleLeadingKaon2->Fill(leadingKaon.pointingAngleHypo());
            HistPointingAngleLeadingKaon2->Fill(subLeadingKaon.pointingAngleHypo());
            HistHypoLengthLeadingKaon2->Fill(leadingKaon.decayLengthHypo());
            HistHypoLengthLeadingKaon2->Fill(subLeadingKaon.decayLengthHypo());
            HistKsi12->Fill(ksi1);
            HistKsi22->Fill(ksi2);
            HistDiffProtonX2->Fill(sumProtonMomentumX);
            HistDiffProtonY2->Fill(sumProtonMomentumY);
            HistKsiCorrelation2->Fill(ksiCorrelation);
            HistEtaCorrelation2->Fill(etaCorrelation);
            HistCorr12->Fill(newCorr1);
            HistCorr22->Fill(newCorr2);
            histZDiff2->Fill(zdiff);
            histZMean2->Fill(zmean);  
            
            histCosThetaStar2->Fill(cosThetaStarLeading);
            histCosThetaStar2->Fill(cosThetaStarSubLeading);   
        }



        if (NtofMatched == 4)
        {
            
            sumTofTracks[NtofMatched]+=1;

            if ((leadingKaon.m() >= 0.48 and leadingKaon.m() < 0.52) and (subLeadingKaon.m() >= 0.48 and subLeadingKaon.m() < 0.52))
            {
                sumMassWindow[NtofMatched]+=1;
            }

            if (pTmiss <= 0.15)
            {
                sumSmallPt[NtofMatched]+=1;
            }

            if (totalCluster  <= (2*Ntof+1))
            {
                sumNtofCluster[NtofMatched]+=1;
            }

            if (leadingKaon.DCABeamLine() < 2.5 and subLeadingKaon.DCABeamLine() < 2.5)
            {
                sumDcaBeamline[NtofMatched]+=1;
            }

            if (leadingKaon.dcaDaughters() < 2.5 and subLeadingKaon.dcaDaughters() < 2.5)
            {
                sumDcaDaughters[NtofMatched]+=1;            
            }
            
            if ((leadingKaon.decayLengthHypo() < 3.0 or leadingKaon.pointingAngleHypo() > 0.925) and (subLeadingKaon.decayLengthHypo() < 3.0 or subLeadingKaon.pointingAngleHypo() > 0.925))
            {
                sumDecayCos[NtofMatched]+=1;
            }

            if (abs(ksi1) >= 0.007 or abs(ksi2) >= 0.007 or abs(sumProtonMomentumX) >= 0.1 or abs(sumProtonMomentumY) >= 0.1 )
            {
                sumAntiEllastic[NtofMatched]+=1;
            }

            if (abs(ksiCorrelation )< 0.004)
            {
                sumKsiCorr[NtofMatched]+=1;
            }

            if (abs(etaCorrelation) < 0.9)
            {
                sumEtaCorr[NtofMatched]+=1;
            }


            if ( abs(newCorr1) < 0.004 and abs(newCorr2) < 0.004)
            {
                sumSingleCorr[NtofMatched]+=1;
            }

            if (abs(zmean)<80)
            {
                effZmean[NtofMatched]+=1;
            }

            if (abs(zdiff)<15)
            {
                effZdiff[NtofMatched]+=1;
            }
            if (abs(cosThetaStarLeading)<0.8 and abs(cosThetaStarSubLeading) < 0.8)
            {
                effCosThetaStar[NtofMatched]+=1;
            }
    }


       else if (NtofMatched == 3)
        {

   
                sumTofTracks[NtofMatched]+=1;


            if ((leadingKaon.m() >= 0.48 and leadingKaon.m() < 0.52) and (subLeadingKaon.m() >= 0.48 and subLeadingKaon.m() < 0.52))
            {
                sumMassWindow[NtofMatched]+=1;
            }

            if (pTmiss <= 0.15)
            {
                sumSmallPt[NtofMatched]+=1;
            }

            if (totalCluster  <= (2*Ntof+1))
            {
                sumNtofCluster[NtofMatched]+=1;
            }

            if (leadingKaon.DCABeamLine() < 2.5 and subLeadingKaon.DCABeamLine() < 1.5)
            {
                sumDcaBeamline[NtofMatched]+=1;
            }

            if (leadingKaon.dcaDaughters() < 2.5 and subLeadingKaon.dcaDaughters() < 1.5)
            {
                sumDcaDaughters[NtofMatched]+=1;            
            }
            
            if ((leadingKaon.decayLengthHypo() < 3.0 or leadingKaon.pointingAngleHypo() > 0.925) and (subLeadingKaon.pointingAngleHypo() > 0.95))
            {
                sumDecayCos[NtofMatched]+=1;
            }

            if (abs(ksi1) >= 0.007 or abs(ksi2) >= 0.007 or abs(sumProtonMomentumX) >= 0.1 or abs(sumProtonMomentumY) >= 0.1 )
            {
                sumAntiEllastic[NtofMatched]+=1;
            }

            if (abs(ksiCorrelation )< 0.004)
            {
                sumKsiCorr[NtofMatched]+=1;
            }

            if (abs(etaCorrelation) < 0.9)
            {
                sumEtaCorr[NtofMatched]+=1;
            }


            if ( abs(newCorr1) < 0.004 and abs(newCorr2) < 0.004)
            {
                sumSingleCorr[NtofMatched]+=1;
            }

                if (abs(zmean)<80)
            {
                effZmean[NtofMatched]+=1;
            }

            if (abs(zdiff)<15)
            {
                effZdiff[NtofMatched]+=1;
            }
            if (abs(cosThetaStarLeading)<0.8 and abs(cosThetaStarSubLeading) < 0.8)
            {
                effCosThetaStar[NtofMatched]+=1;
            }
    }


        if (NtofMatched == 2)
        {
   

                sumTofTracks[NtofMatched]+=1;
        

            if ((leadingKaon.m() >= 0.48 and leadingKaon.m() < 0.52) and (subLeadingKaon.m() >= 0.48 and subLeadingKaon.m() < 0.52))
            {
                sumMassWindow[NtofMatched]+=1;
            }

            if (pTmiss <= 0.15)
            {
                sumSmallPt[NtofMatched]+=1;
            }

            if (totalCluster  <= (2*Ntof+1))
            {
                sumNtofCluster[NtofMatched]+=1;
            }

            if (leadingKaon.DCABeamLine() < 1.5 and subLeadingKaon.DCABeamLine() < 1.5)
            {
                sumDcaBeamline[NtofMatched]+=1;
            }

            if (leadingKaon.dcaDaughters() < 1.5 and subLeadingKaon.dcaDaughters() < 1.5) 
            {
                sumDcaDaughters[NtofMatched]+=1;            
            }
            
            if ((leadingKaon.pointingAngleHypo() > 0.95) and (subLeadingKaon.pointingAngleHypo() > 0.95))
            {
                sumDecayCos[NtofMatched]+=1;
            }

            if (abs(ksi1) >= 0.007 or abs(ksi2) >= 0.007 or abs(sumProtonMomentumX) >= 0.1 or abs(sumProtonMomentumY) >= 0.1 )
            {
                sumAntiEllastic[NtofMatched]+=1;
            }

            if (abs(ksiCorrelation )< 0.004)
            {
                sumKsiCorr[NtofMatched]+=1;
            }

            if (abs(etaCorrelation) < 0.9)
            {
                sumEtaCorr[NtofMatched]+=1;
            }


            if ( abs(newCorr1) < 0.004 and abs(newCorr2) < 0.004)
            {
                sumSingleCorr[NtofMatched]+=1;
            }

            if (abs(zmean)<80)
            {
                effZmean[NtofMatched]+=1;
            }

            if (abs(zdiff)<15)
            {
                effZdiff[NtofMatched]+=1;
            }
            if (abs(cosThetaStarLeading)<0.8 and abs(cosThetaStarSubLeading) < 0.8)
            {
                effCosThetaStar[NtofMatched]+=1;
            }
    }


      if (flags < 1000)
        {     
            vector <vector <double>> dataGoodTof;
            vector <vector <double>> dataTof;
            vector <vector <double>> data;     

            for (int i = 0; i < upcEvt->getNumberOfTracks(); i++)
            {
          
                StUPCTrack *track = upcEvt->getTrack(i);
                vector <double> dGoodTof;
                vector <double> dTof;
                vector <double> d;

                if (track->getFlag(StUPCTrack::kTof)  and  track->getPt() >= ptCutVal and  (abs(track->getEta()) <= etaCutVal) and (track->getNhitsFit() >= NhistFitCut))
                {
                 
                    dGoodTof.push_back(track->getEta());
                    dGoodTof.push_back(track->getPhi());
                    dGoodTof.push_back(track->getPt());
    
                    dataGoodTof.push_back(dGoodTof);

                }

                else if (track->getFlag(StUPCTrack::kTof) )
                {
                    dTof.push_back(track->getEta());
                    dTof.push_back(track->getPhi());
                    dTof.push_back(track->getPt());
                                      
            
                      dataTof.push_back(dTof);

                
                }

                else
                {
                    d.push_back(track->getEta());
                    d.push_back(track->getPhi());
                    d.push_back(track->getPt());
                    data.push_back(d);
                }

   

            }
            

            if (data.size() < 50 and dataTof.size()!=0)
            {


        WriteCSV_file("goodTof"+to_string(flags)+".csv" ,dataGoodTof);
        WriteCSV_file("tof"+to_string(flags)+".csv" ,dataTof);
        WriteCSV_file("otherTpc"+to_string(flags)+".csv" ,data);         
               

            TLorentzVector posPion, negPion;
            vector <vector <double>> da;
            for (int i = 0; i < 2; i++)
            {
                vector <double> d1, d2;
                negPion.SetPxPyPzE(vPionsNeg[i]->Px(), vPionsNeg[i]->Py(), vPionsNeg[i]->Pz(), vPionsNeg[i]->Energy()); 
                posPion.SetPxPyPzE(vPionsPos[i]->Px(), vPionsPos[i]->Py(), vPionsPos[i]->Pz(), vPionsPos[i]->Energy()); 
                d1.push_back(negPion.Eta());
                d1.push_back(negPion.Phi());
                d1.push_back(negPion.Pt());            
                
                d2.push_back(posPion.Eta());
                d2.push_back(posPion.Phi());
                d2.push_back(posPion.Pt());   
                da.push_back(d1);
                da.push_back(d2);           
            }


            WriteCSV_file("true"+to_string(flags)+".csv" ,da);



                vector <vector<double>> prot;
                vector <double > p1;
                vector <double > p2;
                p1.push_back(proton1.Eta());
                p1.push_back(proton1.Phi());
                p1.push_back(proton1.Pt());            
              
                p2.push_back(proton2.Eta());
                p2.push_back(proton2.Phi());
                p2.push_back(proton2.Pt());   
 
                prot.push_back(p1);
                prot.push_back(p2);
            
                WriteCSV_file("pro"+to_string(flags)+".csv" ,prot);                    flags+=1;
            }



        }



        




    }
        

    TFile *outfile = TFile::Open(argv[2], "recreate"); 

    cout << "SUM " << endl; 
    cout << sumEvt << endl;
    HistNumTofTracks4->Scale(1.0/sumEvt);
    HistMassLeadingK04->Scale(1.0/sumEvt);
    HistPtMiss4->Scale(1.0/sumEvt);
    HistNumOfClusters4->Scale(1.0/sumEvt);
    HistDcaBeamlineLeadingKaon4->Scale(1.0/sumEvt);
    HistDcaDaughtersLeadingKaon4->Scale(1.0/sumEvt);
    HistPointingAngleLeadingKaon4->Scale(1.0/sumEvt);
    HistHypoLengthLeadingKaon4->Scale(1.0/sumEvt);
    HistKsi14->Scale(1.0/sumEvt);
    HistKsi24->Scale(1.0/sumEvt);
    HistDiffProtonX4->Scale(1.0/sumEvt);
    HistDiffProtonY4->Scale(1.0/sumEvt);
    HistKsiCorrelation4->Scale(1.0/sumEvt);
    HistEtaCorrelation4->Scale(1.0/sumEvt);
    HistCorr14->Scale(1.0/sumEvt);
    HistCorr24->Scale(1.0/sumEvt);

    histCosThetaStar4->Scale(1.0/sumEvt);
    histCosThetaStar3->Scale(1.0/sumEvt); 
    histCosThetaStar2->Scale(1.0/sumEvt);

    HistNumTofTracks3->Scale(1.0/sumEvt);
    HistMassLeadingK03->Scale(1.0/sumEvt);
    HistPtMiss3->Scale(1.0/sumEvt);
    HistNumOfClusters3->Scale(1.0/sumEvt);
    HistDcaBeamlineLeadingKaon3->Scale(1.0/sumEvt);
    HistDcaDaughtersLeadingKaon3->Scale(1.0/sumEvt);
    HistPointingAngleLeadingKaon3->Scale(1.0/sumEvt);

    HistHypoLengthLeadingKaon3->Scale(1.0/sumEvt);
    HistKsi13->Scale(1.0/sumEvt);
    HistKsi23->Scale(1.0/sumEvt);
    HistDiffProtonX3->Scale(1.0/sumEvt);
    HistDiffProtonY3->Scale(1.0/sumEvt);
    HistKsiCorrelation3->Scale(1.0/sumEvt);
    HistEtaCorrelation3->Scale(1.0/sumEvt);
    HistCorr13->Scale(1.0/sumEvt);
    HistCorr23->Scale(1.0/sumEvt);

    HistNumTofTracks2->Scale(1.0/sumEvt);
    HistMassLeadingK02->Scale(1.0/sumEvt);
    HistPtMiss2->Scale(1.0/sumEvt);
    HistNumOfClusters2->Scale(1.0/sumEvt);
    HistDcaBeamlineLeadingKaon2->Scale(1.0/sumEvt);
    HistDcaDaughtersLeadingKaon2->Scale(1.0/sumEvt);
    HistPointingAngleLeadingKaon2->Scale(1.0/sumEvt);

    HistHypoLengthLeadingKaon2->Scale(1.0/sumEvt);
    HistKsi12->Scale(1.0/sumEvt);
    HistKsi22->Scale(1.0/sumEvt);
    HistDiffProtonX2->Scale(1.0/sumEvt);
    HistDiffProtonY2->Scale(1.0/sumEvt);
    HistKsiCorrelation2->Scale(1.0/sumEvt);
    HistEtaCorrelation2->Scale(1.0/sumEvt);
    HistCorr12->Scale(1.0/sumEvt);
    HistCorr22->Scale(1.0/sumEvt);

    histZMean4->Scale(1.0/sumEvt);
    histZDiff4->Scale(1.0/sumEvt);
    histZMean3->Scale(1.0/sumEvt);
    histZDiff3->Scale(1.0/sumEvt);
    histZMean2->Scale(1.0/sumEvt);
    histZDiff2->Scale(1.0/sumEvt);
    //////////

    for (int i = 0; i < HistNfitPionWithoutTot.size(); i++)
    {
        HistNfitPionWithoutTot[i]->Write();
        HistNfitPionWithTot[i]->Write();
    }

    
    HistNumTofTracks4->Write();
    HistMassLeadingK04->Write();
    HistPtMiss4->Write();
    HistNumOfClusters4->Write();
    HistDcaBeamlineLeadingKaon4->Write();
    HistDcaDaughtersLeadingKaon4->Write();
    HistPointingAngleLeadingKaon4->Write();
    HistHypoLengthLeadingKaon4->Write();
    HistKsi14->Write();
    HistKsi24->Write();
    HistDiffProtonX4->Write();
    HistDiffProtonY4->Write();
    HistKsiCorrelation4->Write();
    HistEtaCorrelation4->Write();
    HistCorr14->Write();
    HistCorr24->Write();

    HistNumTofTracks3->Write();
    HistMassLeadingK03->Write();
    HistPtMiss3->Write();
    HistNumOfClusters3->Write();
    HistDcaBeamlineLeadingKaon3->Write();
    HistDcaDaughtersLeadingKaon3->Write();
    HistPointingAngleLeadingKaon3->Write();

    HistHypoLengthLeadingKaon3->Write();
    HistKsi13->Write();
    HistKsi23->Write();
    HistDiffProtonX3->Write();
    HistDiffProtonY3->Write();
    HistKsiCorrelation3->Write();
    HistEtaCorrelation3->Write();
    HistCorr13->Write();
    HistCorr23->Write();

    HistNumTofTracks2->Write();
    HistMassLeadingK02->Write();
    HistPtMiss2->Write();
    HistNumOfClusters2->Write();
    HistDcaBeamlineLeadingKaon2->Write();
    HistDcaDaughtersLeadingKaon2->Write();
    HistPointingAngleLeadingKaon2->Write();

    HistHypoLengthLeadingKaon2->Write();
    HistKsi12->Write();
    HistKsi22->Write();
    HistDiffProtonX2->Write();
    HistDiffProtonY2->Write();
    HistKsiCorrelation2->Write();
    HistEtaCorrelation2->Write();
    HistCorr12->Write();
    HistCorr22->Write();

    HistEtaPion->Write();
    HistPtPion->Write();
    HistMassK0K0->Write();
    histZMean4->Write();
    histZMean3->Write();
    histZMean2->Write();
    histZDiff4->Write();
    histZDiff3->Write();
    histZDiff2->Write();

    histCosThetaStar4->Write();
    histCosThetaStar3->Write();
    histCosThetaStar2->Write();
    
     
    
    
    outfile->Close();

      
    TCanvas* canvas = new TCanvas("canvas", "Histogram Canvas", 800, 600);
    canvas->SetMargin(0.11, 0.05, 0.14, 0.07);
    TH1F *histogram = new TH1F("histogram", "",13, 0.5, 13 + 0.5);


    histogram->SetBinContent(1,sumMassWindow[4]/double(sumEvt));
    histogram->SetBinContent(2,sumSmallPt[4]/double(sumEvt));
    histogram->SetBinContent(3,sumNtofCluster[4]/double(sumEvt));
    histogram->SetBinContent(4,sumDcaBeamline[4]/double(sumEvt));
    histogram->SetBinContent(5,sumDcaDaughters[4]/double(sumEvt));
    histogram->SetBinContent(6,sumDecayCos[4]/double(sumEvt));
    histogram->SetBinContent(7,sumAntiEllastic[4]/double(sumEvt));
    histogram->SetBinContent(8,sumKsiCorr[4]/double(sumEvt));
    histogram->SetBinContent(9,sumEtaCorr[4]/double(sumEvt));
    histogram->SetBinContent(10,sumSingleCorr[4]/double(sumEvt));
    histogram->SetBinContent(11,effZmean[4]/double(sumEvt));
    histogram->SetBinContent(12,effZdiff[4]/double(sumEvt));
    histogram->SetBinContent(13,effCosThetaStar[4]/double(sumEvt));



   /* histogram->GetXaxis()->SetBinLabel(1, "m_{#pi^{+}#pi^{-}}");
    histogram->GetXaxis()->SetBinLabel(2, "p_{T}^{miss}");
    histogram->GetXaxis()->SetBinLabel(3, "N_{TOF}^{cluster}");
    histogram->GetXaxis()->SetBinLabel(4, "DCA_{beamline}");
    histogram->GetXaxis()->SetBinLabel(5, "DCA_{daughters}");
    histogram->GetXaxis()->SetBinLabel(6, "l_{decay} and cos(#alpha_{p})");
    histogram->GetXaxis()->SetBinLabel(7, "anti-elastic");
    histogram->GetXaxis()->SetBinLabel(8, "m_{KK}/#sqrt{s}#minus#sqrt{#xi^{E}#xi^{W}}");
    histogram->GetXaxis()->SetBinLabel(9, "y_{KK}#minus1/2ln(#xi^{E}/#xi^{W})");
    histogram->GetXaxis()->SetBinLabel(10, "m_{KK}e^{#pmy}#minus#xi_{E/W}");
    histogram->GetXaxis()->SetBinLabel(11, "vtx_{z}^{avg} [cm]");
    histogram->GetXaxis()->SetBinLabel(12, "vtx_{z}^{diff} [cm]");
    histogram->GetXaxis()->SetBinLabel(13, "cos(#theta*)");*/
    histogram->SetLineColor(kBlack);
    histogram->SetLineWidth(2.5);      
    histogram->GetXaxis()->SetLabelSize(0.04);
    histogram->GetYaxis()->SetLabelSize(0.07); 
    histogram->GetYaxis()->SetTitle("#varepsilon");
    histogram->SetMinimum(0);
histogram->SetLabelOffset(0.01);
    histogram->GetYaxis()->SetTitleSize(0.07);    
    histogram->GetYaxis()->SetTitleOffset(1.7); 

    TH1F *histogram3 = new TH1F("histogram3", "",13, 0.5, 13 + 0.5);
    histogram3->SetBinContent(1,sumMassWindow[3]/double(sumEvt));
    histogram3->SetBinContent(2,sumSmallPt[3]/double(sumEvt));
    histogram3->SetBinContent(3,sumNtofCluster[3]/double(sumEvt));
    histogram3->SetBinContent(4,sumDcaBeamline[3]/double(sumEvt));
    histogram3->SetBinContent(5,sumDcaDaughters[3]/double(sumEvt));
    histogram3->SetBinContent(6,sumDecayCos[3]/double(sumEvt));
    histogram3->SetBinContent(7,sumAntiEllastic[3]/double(sumEvt));
    histogram3->SetBinContent(8,sumKsiCorr[3]/double(sumEvt));
    histogram3->SetBinContent(9,sumEtaCorr[3]/double(sumEvt));
    histogram3->SetBinContent(10,sumSingleCorr[3]/double(sumEvt));
    histogram3->SetBinContent(11,effZmean[3]/double(sumEvt));
    histogram3->SetBinContent(12,effZdiff[3]/double(sumEvt));
    histogram3->SetBinContent(13,effCosThetaStar[3]/double(sumEvt));

    histogram3->SetLineColor(kRed);
    histogram3->SetLineWidth(2.5);   

    TH1F *histogram2 = new TH1F("histogram2", "",13, 0.5, 13 + 0.5);
    histogram2->SetBinContent(1,sumMassWindow[2]/double(sumEvt));
    histogram2->SetBinContent(2,sumSmallPt[2]/double(sumEvt));
    histogram2->SetBinContent(3,sumNtofCluster[2]/double(sumEvt));
    histogram2->SetBinContent(4,sumDcaBeamline[2]/double(sumEvt));
    histogram2->SetBinContent(5,sumDcaDaughters[2]/double(sumEvt));
    histogram2->SetBinContent(6,sumDecayCos[2]/double(sumEvt));
    histogram2->SetBinContent(7,sumAntiEllastic[2]/double(sumEvt));
    histogram2->SetBinContent(8,sumKsiCorr[2]/double(sumEvt));
    histogram2->SetBinContent(9,sumEtaCorr[2]/double(sumEvt));
    histogram2->SetBinContent(10,sumSingleCorr[2]/double(sumEvt));
    histogram2->SetBinContent(11,effZmean[2]/double(sumEvt));
    histogram2->SetBinContent(12,effZdiff[2]/double(sumEvt));
    histogram2->SetBinContent(13,effCosThetaStar[2]/double(sumEvt));
    histogram2->SetLineColor(1);
    histogram2->SetLineWidth(2.5);   

    histogram->Draw("hist");
    histogram3->Draw("hist same");
    histogram2->Draw("hist same");
    histogram->SetStats(0);
  
    histogram3->SetStats(0);
    histogram2->SetStats(0);
    double maxRange = (std::max(histogram3->GetMaximum(), histogram2->GetMaximum() ))*1.15 ;
    histogram->GetYaxis()->SetRangeUser(0, maxRange);
    histogram2->GetYaxis()->SetRangeUser(0, maxRange);
    histogram3->GetYaxis()->SetRangeUser(0, maxRange);

    TLegend *legend = new TLegend(0.15, 0.2, .5, 0.42);
    legend->SetMargin(0.15);
    legend->AddEntry(histogram, "4 TOF trakcs", "f"); 
    legend->AddEntry(histogram3, "3 TOF trakcs", "f"); 
    legend->AddEntry(histogram2, "2 TOF trakcs", "f"); 
    legend->SetTextSize(0.16);
    legend->SetBorderSize(0);
    legend->SetFillColorAlpha(0,0.0);
    legend->Draw();
    canvas->SaveAs("DataAnalysis/perf.pdf");
    delete canvas;



    // Get the number of bins
    int nBins = histogram->GetNbinsX();

    // Create arrays to hold the x and y values
    double *x1 = new double[nBins];
    double *y1 = new double[nBins];
    double *x2 = new double[nBins];
    double *y2 = new double[nBins];
    double *x3 = new double[nBins];
    double *y3 = new double[nBins];

 
    double *err1 = new double[nBins];   
    double *err2 = new double[nBins];
    double *err3 = new double[nBins];
    double *e0 = new double[nBins];

  std::cout << std::fixed << std::setprecision(3)  << endl;
    // Extract data points from the histograms and slightly shift the x positions
    for (int i = 0; i < nBins; ++i) 
    {
        double binCenter = histogram->GetBinCenter(i+1);
        x1[i] = binCenter - 0.1; // Slightly shift left

        y1[i] = histogram->GetBinContent(i+1);

        x2[i] = binCenter; // No shift
        y2[i] = histogram2->GetBinContent(i+1);

        x3[i] = binCenter + 0.1; // Slightly shift right
        y3[i] = histogram3->GetBinContent(i+1);
        e0[i] = 0.0;
        
        err1[i] = duncErr(y1[i]*604, 604.0-y1[i]*604,  sqrt(y1[i]*604), sqrt(604.0-y1[i]*604 ) );
        err2[i] = duncErr(y2[i]*604, 604.0-y2[i]*604,  sqrt(y2[i]*604), sqrt(604.0-y2[i]*604) );
        err3[i] = duncErr(y3[i]*604, 604.0-y3[i]*604,       sqrt(y3[i]*604), sqrt(604.0-y3[i]* 604) );

        cout << "--------------------  " << i+1 << " -------------" << endl;
        cout << "&"<<  4 << " & " <<  y1[i]*604   <<" & " <<  sqrt(y1[i]*604)   <<" & " <<      604.0-y1[i]*604 <<" & " <<  sqrt( 604.0-y1[i]*604) <<" & " <<  y1[i]  <<" & " <<  err1[i] << "\\\\" << endl;
         cout  << "&" <<3 << " & " << y3[i]*604   <<" & " <<  sqrt(y3[i]*604)   <<" & " <<      604.0-y3[i]*604 <<" & " <<  sqrt( 604.0-y3[i]*604) <<" & " <<  y3[i]  <<" & " <<  err3[i] << "\\\\" << endl;

       cout  << "&"<< 2 << " & " << y2[i]*604   <<" & " <<  sqrt(y2[i]*604)   <<" & " <<      604.0-y2[i]*604 <<" & " <<  sqrt( 604.0-y2[i]*604) <<" & " <<  y2[i]  <<" & " <<  err2[i] << "\\\\" << endl;

      
        cout << endl << endl;

    }





    // Create TGraphs from the extracted data points
    TGraphErrors *graph1 = new TGraphErrors(nBins, x1, y1,e0, err1);
    TGraphErrors *graph2 = new TGraphErrors(nBins, x2, y2, e0, err2);
    TGraphErrors *graph3 = new TGraphErrors(nBins, x3, y3,e0,err3);

    // Set graph titles and colors
    graph1->SetTitle("Shifted Histograms to TGraphs;X axis;Y axis");
    graph1->SetLineColor(1);
    graph3->SetLineColor(1);
    graph2->SetLineColor(1);
    graph1->SetMarkerColor(1);
    graph3->SetMarkerColor(1);
    graph2->SetMarkerColor(1);
   graph1->SetMarkerStyle(21);
    graph3->SetMarkerStyle(20);
    graph2->SetMarkerStyle(22);
   graph1->SetMarkerSize(2);
    graph3->SetMarkerSize(2);
    graph2->SetMarkerSize(2);
    graph1->SetLineWidth(2);
    graph3->SetLineWidth(2);
    graph2->SetLineWidth(2);

    // Draw the histograms and graphs
    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
    c1->SetMargin(0.13, 0.05, 0.32, 0.03);
    histogram->SetLineColorAlpha(kAzure+6, 1);
    histogram3->SetLineColorAlpha(kViolet+6, 1);
    histogram2->SetLineColorAlpha(kMagenta-10,1);
histogram3->SetLineWidth(3);
histogram3->SetLineWidth(3);
histogram->SetLineWidth(3);

    histogram->Draw();
    histogram->GetYaxis()->SetNdivisions(5);
    histogram3->Draw("SAME");
    histogram2->Draw("SAME");
    histogram->GetXaxis()->SetBinLabel(1, "m_{#pi^{+}#pi^{-}}");
    histogram->GetXaxis()->SetBinLabel(2, "p_{T}^{miss}");
    histogram->GetXaxis()->SetBinLabel(3, "N_{TOF}^{cluster}");
    histogram->GetXaxis()->SetBinLabel(4, "DCA_{beamline}");
    histogram->GetXaxis()->SetBinLabel(5, "DCA_{daughters}");
    histogram->GetXaxis()->SetBinLabel(6, "l_{decay}, cos(#alpha_{p})");
    histogram->GetXaxis()->SetBinLabel(7, "anti-elastic");
    histogram->GetXaxis()->SetBinLabel(8, "m_{KK}/#sqrt{s}#minus#sqrt{#xi^{E}#xi^{W}}");
    histogram->GetXaxis()->SetBinLabel(9, "y_{KK}#minusln(#xi^{E}/#xi^{W})/2");
    histogram->GetXaxis()->SetBinLabel(10, "m_{KK}e^{#pmy}#minus#xi_{E/W}");
    histogram->GetXaxis()->SetBinLabel(11, "vtx_{z}^{mean}");
    histogram->GetXaxis()->SetBinLabel(12, "vtx_{z}^{diff}");
    histogram->GetXaxis()->SetBinLabel(13, "cos(#theta*)");
    histogram->LabelsOption("v"); histogram->GetXaxis()->SetLabelSize(0.07);
    graph1->Draw("PE SAME"); 
    graph3->Draw("PE SAME");
    graph2->Draw("PE SAME");histogram->GetYaxis()->SetTitle("#varepsilon");
    histogram->GetYaxis()->SetTitleOffset(0.76);
    histogram->SetLineStyle(2);
    histogram3->SetLineStyle(4);
    histogram2->SetLineStyle(7);

histogram->GetYaxis()->SetLabelSize(0.07);
histogram->GetXaxis()->SetLabelSize(0.07);
histogram->GetYaxis()->SetTitleSize(0.07);
histogram->GetXaxis()->SetTitleSize(0.07);

    TLegend *legend2 = new TLegend(0.18, 0.35, .62, 0.5);
    legend2->SetMargin(0.15);
    legend2->AddEntry(graph1, "4 TOF ", "p"); 
    legend2->AddEntry(graph3, "3 TOF ", "p"); 
    legend2->AddEntry(graph2, "2 TOF ", "p"); 
    legend2->SetTextSize(0.06);
    legend2->SetBorderSize(0);
    legend2->SetFillColorAlpha(0,0.0);
    legend2->Draw();
    TLatex* latex = new TLatex();
    latex->SetTextSize(0.05);
    latex->DrawLatexNDC(0.18, 0.92, "STAR pp #sqrt{s} = 510 GeV");

    c1->SaveAs("a.pdf");


    TFile* fileData = TFile::Open(argv[2]);


    TH2D* HistMassK0K0PLot = dynamic_cast<TH2D*>(fileData->Get("HistMassK0K0"));  
    TCanvas* CanvasHistMassK0K0PLot = new TCanvas("CanvasHistMassK0K0PLot", "Histogram Canvas", 800, 600);
    gStyle->SetPalette(kBlueRedYellow);
    HistMassK0K0PLot->Draw("COLZ");
    CanvasHistMassK0K0PLot->SetMargin(0.14, 0.15, 0.14, 0.04);
    HistMassK0K0PLot->GetXaxis()->SetLabelSize(0.07);
    HistMassK0K0PLot->GetYaxis()->SetLabelSize(0.07); 
    HistMassK0K0PLot->GetXaxis()->SetTitleSize(0.07);
    HistMassK0K0PLot->GetYaxis()->SetTitleSize(0.07); 
    HistMassK0K0PLot->GetZaxis()->SetTitleSize(0.07); 
    HistMassK0K0PLot->GetZaxis()->SetLabelSize(0.07);
    HistMassK0K0PLot->GetXaxis()->SetTitle("m_{p^{+}#pi^{-}}^{det}[GeV]");
    HistMassK0K0PLot->GetYaxis()->SetTitle("m_{p^{-}#pi^{+}}^{true} [GeV]");
    HistMassK0K0PLot->SetTitle("");
    HistMassK0K0PLot->GetYaxis()->SetTitleOffset(1.36);
    HistMassK0K0PLot->GetXaxis()->SetTitleOffset(1.25);
    HistMassK0K0PLot->SetStats(0);
    CanvasHistMassK0K0PLot->SaveAs("HistMassK0K0PLot.pdf");



    TH2D* HistPtPionPlot = dynamic_cast<TH2D*>(fileData->Get("HistPtPion"));  
    TCanvas* CanvasHistPtPionPlot = new TCanvas("CanvasHistPtPionPlot", "Histogram Canvas", 800, 600);
    gStyle->SetPalette(kBlueRedYellow);
    HistPtPionPlot->Draw("COLZ");
    CanvasHistPtPionPlot->SetMargin(0.14, 0.15, 0.14, 0.04);
    HistPtPionPlot->GetXaxis()->SetLabelSize(0.07);
    HistPtPionPlot->GetYaxis()->SetLabelSize(0.07); 
    HistPtPionPlot->GetXaxis()->SetTitleSize(0.07);
    HistPtPionPlot->GetYaxis()->SetTitleSize(0.07); 
    HistPtPionPlot->GetZaxis()->SetTitleSize(0.07); 
    HistPtPionPlot->GetZaxis()->SetLabelSize(0.07);
    HistPtPionPlot->GetXaxis()->SetTitle("m_{p^{+}#pi^{-}}^{det}[GeV]");
    HistPtPionPlot->GetYaxis()->SetTitle("m_{p^{-}#pi^{+}}^{true} [GeV]");
    HistPtPionPlot->SetTitle("");
    HistPtPionPlot->GetYaxis()->SetTitleOffset(1.36);
    HistPtPionPlot->GetXaxis()->SetTitleOffset(1.25);
    HistPtPionPlot->SetStats(0);
    CanvasHistPtPionPlot->SaveAs("HistPtPionPlot.pdf");



    TH2D* HistEtaPionPlot = dynamic_cast<TH2D*>(fileData->Get("HistEtaPion"));  
    TCanvas* CanvasHistEtaPionPlot = new TCanvas("CanvasHistEtaPionPlot", "Histogram Canvas", 800, 600);
    gStyle->SetPalette(kBlueRedYellow);
    HistEtaPionPlot->Draw("COLZ");
    CanvasHistEtaPionPlot->SetMargin(0.14, 0.15, 0.14, 0.04);
    HistEtaPionPlot->GetXaxis()->SetLabelSize(0.07);
    HistEtaPionPlot->GetYaxis()->SetLabelSize(0.07); 
    HistEtaPionPlot->GetXaxis()->SetTitleSize(0.07);
    HistEtaPionPlot->GetYaxis()->SetTitleSize(0.07); 
    HistEtaPionPlot->GetZaxis()->SetTitleSize(0.07); 
    HistEtaPionPlot->GetZaxis()->SetLabelSize(0.07);
    HistEtaPionPlot->GetXaxis()->SetTitle("m_{p^{+}#pi^{-}}^{det}[GeV]");
    HistEtaPionPlot->GetYaxis()->SetTitle("m_{p^{-}#pi^{+}}^{true} [GeV]");
    HistEtaPionPlot->SetTitle("");
    HistEtaPionPlot->GetYaxis()->SetTitleOffset(1.36);
    HistEtaPionPlot->GetXaxis()->SetTitleOffset(1.25);
    HistEtaPionPlot->SetStats(0);
    CanvasHistEtaPionPlot->SaveAs("HistEtaPionPlot.pdf");



    plotHistogram(argv[2],"HistMassLeadingK0_4", "HistMassLeadingK0_3", "HistMassLeadingK0_2", 0.52, 0.48,  0.2,0.88,  0.2 , 0.6 ,"leading kaon", "sub-leading kaon" ,2,0,0);
    plotHistogram(argv[2],"HistDcaDaughtersLeadingKaon_4", "HistDcaDaughtersLeadingKaon_3", "HistDcaDaughtersLeadingKaon_2", 2.5,2.5,  0.5,0.88,  0.75 , 0.6 ,"leading kaon", "sub-leading kaon" ,1,0,0);
    plotHistogram(argv[2],"HistDcaBeamlineLeadingKaon_4", "HistDcaBeamlineLeadingKaon_3", "HistDcaBeamlineLeadingKaon_2", 2.5, 2.5,  0.5, 0.88,  0.75 , 0.6 ,"leading kaon", "sub-leading kaon" ,1,0,0);
    plotHistogram(argv[2],"HistPointingAngleLeadingKaon_4", "HistPointingAngleLeadingKaon_3", "HistPointingAngleLeadingKaon_2",0.925,0.95,  0.2,0.88,  0.2 , 0.6 ,"leading kaon", "sub-leading kaon" ,0,0,0);
    plotHistogram(argv[2],"HistHypoLengthLeadingKaon_4", "HistHypoLengthLeadingKaon_3", "HistHypoLengthLeadingKaon_2", 3.0,3.0,  0.5,0.88,  0.75 , 0.6 ,"leading kaon", "sub-leading kaon" ,1,0,0);
    plotHistogram(argv[2],"HistKsi1_4", "HistKsi1_3", "HistKsi1_2", 0.007, -0.007,  0.2,0.87,  0.2 , 0.6 ,"p^{E}", "p^{W}" ,4,0,0);
    plotHistogram(argv[2],"HistKsi2_4", "HistKsi2_3", "HistKsi2_2", 0.007, -0.007,  0.2,0.87,  0.2 , 0.6 ,"p^{E}", "p^{W}" ,4,0,0);   

    plotHistogram(argv[2],"HistDiffProtonX_4", "HistDiffProtonX_3", "HistDiffProtonX_2", 0.1, -0.1,  0.2,0.89,  0.68 , 0.7 ,"p^{E}_{x}+p^{W}_{x}","p^{E}_{y}+p^{W}_{y}" ,4,0,0);
    plotHistogram(argv[2],"HistDiffProtonY_4", "HistDiffProtonY_3", "HistDiffProtonY_2", 0.1, -0.1,  0.2,0.89,  0.8 , 0.72 ,"p^{E}_{x}+p^{W}_{x}","p^{E}_{y}+p^{W}_{y}" ,4,0,0);

    plotHistogram(argv[2],"HistCorr1_4", "HistCorr1_3", "HistCorr1_2", 0.004, -0.004,  0.2,0.89,  0.72 , 0.7 ,"p^{E}", "p^{W}" ,2,-0.1, 0.1);
    plotHistogram(argv[2],"HistCorr2_4", "HistCorr2_3", "HistCorr2_2", 0.004, -0.004,  0.2,0.89,  0.72 , 0.7 ,"p^{E}", "p^{W}" ,2,-0.1,0.1);

    plotHistogram(argv[2], "HistPtMiss_4", "HistPtMiss_3", "HistPtMiss_2", 0.15, 0.15,  0.5, 0.88,  0.75,0.6  ,"a","a", 1,0,0);
    plotHistogram(argv[2], "HistNumOfClusters_4", "HistNumOfClusters_3", "HistNumOfClusters_2", 9, 9,  0.2, 0.88,  0.75,0.7  ,"a","a", 1,0,0);
    plotHistogram(argv[2], "HistKsiCorrelation_4", "HistKsiCorrelation_3", "HistKsiCorrelation_2", 0.004, -0.004,  0.5, 0.89,  0.2,0.7  ,"a","a", 2,0,0);
    plotHistogram(argv[2], "HistEtaCorrelation_4", "HistEtaCorrelation_3", "HistEtaCorrelation_2", 0.9, -0.9,  0.5, 0.89,  0.2,0.7  ,"a","a", 2,0,0);

    plotHistogram(argv[2], "HistNumTofTracks_4", "HistNumTofTracks_3", "HistNumTofTracks_2", 1.5, 4.5,  0.6, 0.88,  0.6,0.6  ,"a","a", 2,0,0);
    plotHistogram(argv[2], "histZMean4", "histZMean3", "histZMean2", 80, -80,  0.2, 0.88,  0.75,0.7  ,"a","a", 2,0,0);
    plotHistogram(argv[2], "histZDiff4", "histZDiff3", "histZDiff2", 15, -15,  0.2, 0.88,  0.75,0.7  ,"a","a", 2,0,0);

    plotHistogram(argv[2], "histCosThetaStar4", "histCosThetaStar3", "histCosThetaStar2", 0.8, -0.8,  0.2, 0.88,  0.75,0.68  ,"a","a", 2,0,0);
      
  
cout << d4 << endl << d3 << endl << d2 << endl;



    return 0;
}



void FindProtons(bool isMC, StRPEvent *rpEvt, StUPCEvent *upcEvt, TLorentzVector & proton1, TLorentzVector & proton2)
{
    double massProton = 938.272/1000.0;
    double beamEnergy = 254.867;
    if (isMC == 0)
    {
		
            StUPCRpsTrack *trk = rpEvt->getTrack(0);
            StUPCRpsTrack *trk2 = rpEvt->getTrack(1);
           trk->setEvent(rpEvt);
            trk2->setEvent(rpEvt);
            if (  trk->branch() < 2)
            {
                proton1.SetPxPyPzE(beamEnergy*trk->thetaRp(0), beamEnergy*trk->thetaRp(1), -(beamEnergy- beamEnergy*trk->thetaRp(1)-beamEnergy*trk->thetaRp(0)), beamEnergy);
                proton2.SetPxPyPzE(beamEnergy*trk2->thetaRp(0), beamEnergy*trk2->thetaRp(1),beamEnergy- beamEnergy*trk2->thetaRp(1)-beamEnergy*trk2->thetaRp(0), beamEnergy);    
            }
            else
            {
                proton2.SetPxPyPzE(beamEnergy*trk->thetaRp(0), beamEnergy*trk->thetaRp(1), -(beamEnergy- beamEnergy*trk->thetaRp(1)-beamEnergy*trk->thetaRp(0)), beamEnergy);
                proton1.SetPxPyPzE(beamEnergy*trk2->thetaRp(0), beamEnergy*trk2->thetaRp(1),beamEnergy- beamEnergy*trk2->thetaRp(1)-beamEnergy*trk2->thetaRp(0), beamEnergy);    
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
        double pz_reco1 = gRandom->Gaus(vProtons[0]->Pz(), 0.51);
        proton1.SetXYZM(px_reco1,py_reco1, pz_reco1, massProton);

        double px_reco2 = gRandom->Gaus(vProtons[1]->Px(), sigma);
        double py_reco2 = gRandom->Gaus(vProtons[1]->Py(), sigma);
        double pz_reco2 = gRandom->Gaus(vProtons[1]->Pz(), 0.51);
        proton2.SetXYZM(px_reco2,py_reco2, pz_reco2, massProton);
        

        if (pz_reco1 < 0.0)
        {
            proton1.SetXYZM(px_reco1,py_reco1, pz_reco1, massProton);
            proton2.SetXYZM(px_reco2,py_reco2, pz_reco2, massProton);
        }

        else
        {
            proton2.SetXYZM(px_reco1,py_reco1, pz_reco1, massProton);
            proton1.SetXYZM(px_reco2,py_reco2, pz_reco2, massProton);
        }



        vProtons.clear();
    }

   
}


///////////////

void SetYAxisRange(TH1* hist1, TH1* hist2)
{
    double maxRange = std::max(hist1->GetMaximum(), hist2->GetMaximum());

    hist1->GetYaxis()->SetRangeUser(0, maxRange);
    hist2->GetYaxis()->SetRangeUser(0, maxRange);
}

void plotHistogram(const char* inputFileData, const char* histogramName,  const char* histogramName2, const char* histogramName3, double val1, double val2,  double latexX, double latexY,  double xlmin, double ylmin ,string l1, string l2, int b , double XMIN , double XMAX) 
{

    TFile* fileData = TFile::Open(inputFileData);

    TH1D* histData = dynamic_cast<TH1D*>(fileData->Get(histogramName));
    TH1D* histData2 = dynamic_cast<TH1D*>(fileData->Get(histogramName2));
    TH1D* histData3 = dynamic_cast<TH1D*>(fileData->Get(histogramName3));


    std::string str = std::string(histogramName);
    
      
    TCanvas* canvas = new TCanvas("canvas", "Histogram Canvas", 800, 600);
    canvas->SetMargin(0.16, 0.04, 0.18, 0.04);

    
     
    histData->Draw("hist");
    histData2->Draw("hist same");
    histData3->Draw("hist same");
    histData->Draw("hist same");
    histData->SetStats(0);
    histData->SetLineWidth(2);

    histData->SetLineColor(kAzure+6);
    histData->SetMarkerColor(kAzure+6);
    histData->SetFillColorAlpha(kAzure+6,1);
    histData->SetMarkerColor(kAzure+6);
    histData->SetTitle("");
    histData->GetXaxis()->SetLabelSize(0.07);
    histData->GetYaxis()->SetLabelSize(0.07); 
    histData->GetXaxis()->SetTitleSize(0.07);
    histData->GetYaxis()->SetTitleSize(0.07); 
    histData->GetXaxis()->SetTitleOffset(1.16); 
    histData->GetYaxis()->SetTitleOffset(1.17); 
    histData->SetMarkerStyle(1);
    histData->SetMinimum(0);
    histData->SetStats(0);

    histData2->SetLineWidth(2);
    histData2->SetLineColor(kViolet+6);
    histData2->SetFillColorAlpha(kViolet+6,1);
    histData2->SetLineWidth(2);
    histData2->SetMarkerColor(kViolet+6);



    histData3->SetStats(0);
    histData3->SetLineWidth(2);
    histData3->SetFillColorAlpha(kMagenta-10,1);
    histData3->SetLineColor(kMagenta-10);
    histData3->SetMarkerColor(kMagenta-10);
   


    double maxRange = (std::max(histData2->GetMaximum(), histData3->GetMaximum() ))*1.15 ;
    histData->GetYaxis()->SetRangeUser(0, maxRange);
    histData2->GetYaxis()->SetRangeUser(0, maxRange);
    histData3->GetYaxis()->SetRangeUser(0, maxRange);


    SetYAxisRange(histData, histData2);



    TLegend* legend = new TLegend(xlmin,ylmin, xlmin+0.34, ylmin+0.22); 
    const char *cstr1 = l1.c_str();

    const char *cstr2 = l2.c_str();

    legend->AddEntry(histData, "4 TOF ", "f"); 
    legend->AddEntry(histData2, "3 TOF ", "f"); 
    legend->AddEntry(histData3, "2 TOF ", "f"); 
    legend->SetTextSize(0.06);
    legend->SetMargin(0.2);
    legend->SetBorderSize(0);
    legend->SetFillColorAlpha(0,0.0);
    legend->Draw();


    if (b == 2)
    {
        TBox *box = new TBox(val1, histData2->GetMinimum(), val2, histData->GetMaximum());
        box->SetFillColorAlpha(11, 0.27); // Set fill color to blue with transparency
        box->SetLineColor(1); // Set line color to blue
        box->SetLineWidth(0); // Set line width to 0 to hide the outline
        //box->Draw();

        auto ar2 = new TArrow(val1, 0.00, val1,histData2->GetMaximum()/2,1,"|->");
        ar2->SetLineWidth(2);
        ar2->Draw("same");

        auto ar22 = new TArrow(val1,  histData2->GetMaximum()/2, val1-histData->GetBinWidth(1), histData2->GetMaximum()/2,0.03,"|>");
        ar22->SetLineWidth(2);
         ar22->SetFillColor(1);
        ar22->Draw();
        auto ar1 = new TArrow(val2, 0.00, val2, histData2->GetMaximum()/2,1,"|->");
        ar1->SetLineWidth(2);
        ar1->Draw("same");
   
        auto ar12 = new TArrow(val2,histData2->GetMaximum()/2,  val2+histData->GetBinWidth(1), histData2->GetMaximum()/2,0.03,"|>");
        ar12->SetLineWidth(2);
        ar12->SetFillColor(1);
        ar12->Draw();

        ar12->SetLineColor(1);
        ar22->SetLineColor(1);
        ar2->SetLineColor(1);
        ar1->SetLineColor(1);

         if ("HistNumTofTracks_4" == histogramName)
        {
                         box->SetFillColorAlpha(kBlack, 0.07);
               ////box->Draw();
            TBox *box2 = new TBox(2.5, histData2->GetMinimum(), 3.5, histData->GetMaximum());
            box2->SetFillColorAlpha(kRed-7, 0.07); // Set fill color to blue with transparency
            box2->SetLineColor(1); // Set line color to blue
            box2->SetLineWidth(0); // Set line width to 0 to hide the outline
            ////////box2->Draw();
            
            TBox *box3 = new TBox(1.5, histData2->GetMinimum(), 2.5, histData->GetMaximum());
            box3->SetFillColorAlpha(1, 0.07); // Set fill color to blue with transparency
            box3->SetLineColor(1); // Set line color to blue
            box3->SetLineWidth(0); // Set line width to 0 to hide the outline
            //box3->Draw();
        }
    }

    if (b == 1)
    {


        auto ar2 = new TArrow(val1, 0.00, val1,histData2->GetMaximum()/2,1,"|->");
        ar2->SetLineWidth(2);
        ar2->Draw("same");
        ar2->SetLineColor(1);
        auto ar22 = new TArrow(val1,  histData2->GetMaximum()/2, val1-histData->GetBinWidth(1)/2, histData2->GetMaximum()/2,0.03,"|>");
        ar22->SetLineWidth(2);

        ar22->SetFillColor(1);
        ar22->SetLineColor(1);
        ar22->Draw();
        cout << histogramName << endl;


        if (histogramName == "HistDcaDaughtersLeadingKaon_4" or histogramName == "HistDcaBeamlineLeadingKaon_4")
        {
        auto ar2 = new TArrow(1.5, 0.00, 1.5,histData2->GetMaximum()/2,1,"|->");
        ar2->SetLineWidth(2);
        ar2->Draw("same");
        ar2->SetLineColor(1);
        auto ar22 = new TArrow(1.5,  histData2->GetMaximum()/2, 1.5-histData->GetBinWidth(1)/2, histData2->GetMaximum()/2,0.03,"|>");
        ar22->SetLineWidth(2);
        ar22->SetFillColor(1);
        ar22->SetLineColor(1);
        ar22->Draw();
        
        }
        if (histogramName == "HistPointingAngleLeadingKaon_4")
        {
            cout << "HEEE" << endl;
        auto ar2 = new TArrow(1.5, 0.00, 1.5,histData2->GetMaximum()/2,1,"|->");
        ar2->SetLineWidth(2);
        ar2->Draw("same");
        ar2->SetLineColor(1);
        auto ar22 = new TArrow(1.5,  histData2->GetMaximum()/2, 1.5-histData->GetBinWidth(1)/4, histData2->GetMaximum()/2,0.03,"|>");
        ar22->SetLineWidth(2);
        ar22->SetFillColor(1);
        ar22->SetLineColor(1);
        ar22->Draw();
        
        }
        if (histogramName =="HistNumOfClusters_4")
        {

            ar22->SetLineColor(kAzure);
            ar22->SetFillColor(kAzure);
            ar2->SetLineColor(kAzure);
        
            TBox *box3 = new TBox(histData2->GetXaxis()->GetXmin(), histData2->GetMinimum(), 5, histData->GetMaximum());
            box3->SetFillColorAlpha(1, 0.07); // Set fill color to blue with transparency
            box3->SetLineColor(1); // Set line color to blue
            box3->SetLineWidth(0); // Set line width to 0 to hide the outline


            auto ar1 = new TArrow(7, 0.00, 7,histData2->GetMaximum()/2,1,"|->");
            ar1->SetLineWidth(2);
            ar1->Draw("same");
            ar1->SetLineColor(kViolet+6);
            auto ar12 = new TArrow(7,  histData2->GetMaximum()/2, 7-histData->GetBinWidth(1)/2, histData2->GetMaximum()/2,0.03,"|>");
            ar12->SetLineWidth(2);
            ar12->SetFillColor(1);
            ar12->SetFillColor(kViolet+6);
            ar12->SetLineColor(kViolet+6);
            ar12->Draw();

            auto ar0 = new TArrow(5, 0.00, 5,histData2->GetMaximum()/2,1,"|->");
            ar0->SetLineWidth(2);
            ar0->Draw("same");
            ar0->SetLineColor(kMagenta-10);
            auto ar02 = new TArrow(5,  histData2->GetMaximum()/2, 5-histData->GetBinWidth(1)/2, histData2->GetMaximum()/2,0.03,"|>");
            ar02->SetLineWidth(2);
            ar02->SetFillColor(1);
            ar02->SetFillColor(kMagenta-10);
            ar02->SetLineColor(kMagenta-10);
            ar02->Draw();


            //box3->Draw();
        }


    }

    if (b == 0)
    {

        auto ar2 = new TArrow(val1, 0.00, val1,histData2->GetMaximum()/2,1,"|->");
        ar2->SetLineWidth(2);
        ar2->Draw("same");
ar2->SetLineColor(1);
    auto ar22 = new TArrow(val1,  histData2->GetMaximum()/2, val1+histData->GetBinWidth(1), histData2->GetMaximum()/2,0.03,"|>");
        ar22->SetLineWidth(2);
           ar22->SetFillColor(1);
        ar22->Draw();


        auto ar1 = new TArrow(val2, 0.00, val2,histData2->GetMaximum()/2,1,"|->");
        ar1->SetLineWidth(2);
        ar1->Draw("same");
ar1->SetLineColor(1);
    auto ar12 = new TArrow(val2,  histData2->GetMaximum()/2, val2+histData->GetBinWidth(1), histData2->GetMaximum()/2,0.03,"|>");
        ar12->SetLineWidth(2);
           ar12->SetFillColor(1);
           ar12->SetLineColor(1);
            ar22->SetLineColor(1);
        ar12->Draw();




        //box->Draw();
    }

        

    if (b == 4)
    {
        TBox *box = new TBox(val1, histData2->GetMinimum(), histData2->GetXaxis()->GetXmax(), histData->GetMaximum());
        box->SetFillColorAlpha(11, 0.27); // Set fill color to blue with transparency
        box->SetLineColor(1); // Set line color to blue
        box->SetLineWidth(0); // Set line width to 0 to hide the outline

        TBox *box3 = new TBox( histData2->GetXaxis()->GetXmin(), histData2->GetMinimum(),-1*val1, histData->GetMaximum());
        box3->SetFillColorAlpha(11, 0.27); // Set fill color to blue with transparency
        box3->SetLineColor(1); // Set line color to blue
        box3->SetLineWidth(0); // Set line width to 0 to hide the outline
        TLine line4 = TLine(-1*val1, histData->GetMinimum(), -1*val1, histData->GetMaximum());
        line4.SetLineColor(kRed);
        line4.SetLineWidth(4);
        line4.SetLineStyle(2);
        line4.Draw();


        auto ar2 = new TArrow(val1, 0.00, val1,histData2->GetMaximum()/2,1,"|->");
        ar2->SetLineWidth(2);
        ar2->Draw("same");

        auto ar22 = new TArrow(val1,  histData2->GetMaximum()/2, val1-histData->GetBinWidth(1), histData2->GetMaximum()/2,0.03,"|>");
        ar22->SetLineWidth(2);
         ar22->SetFillColor(1);
        ar22->Draw();
        auto ar1 = new TArrow(val2, 0.00, val2, histData2->GetMaximum()/2,1,"|->");
        ar1->SetLineWidth(2);
        ar1->Draw("same");
   
        auto ar12 = new TArrow(val2,histData2->GetMaximum()/2,  val2+histData->GetBinWidth(1), histData2->GetMaximum()/2,0.03,"|>");
        ar12->SetLineWidth(2);
        ar12->SetFillColor(1);
        ar12->Draw();

        ar12->SetLineColor(1);
        ar22->SetLineColor(1);
        ar2->SetLineColor(1);
        ar1->SetLineColor(1);

            line4.Draw();
    }



    TLine line = TLine(val1, histData->GetMinimum(), val1, histData->GetMaximum());
    line.SetLineColor(11);
    line.SetLineWidth(2);
    line.SetLineStyle(2);
       if (histogramName =="HistNumOfClusters_4" or (histogramName =="HistNumTofTracks_4") )
           line.SetLineWidth(0);
    ////line.Draw();

    TLine line2 = TLine(val2, histData->GetMinimum(), val2, histData->GetMaximum());
    
    line2.SetLineColor(11);
    line2.SetLineWidth(2);
           if (histogramName =="HistNumOfClusters_4" or "HistNumTofTracks_4" ==histogramName )
           line2.SetLineWidth(0);
    line2.SetLineStyle(2);
    //////line2.Draw();
    TLatex* latex = new TLatex();
    latex->SetTextSize(0.06);
    latex->DrawLatexNDC(latexX, latexY, "STAR pp #sqrt{s} = 510 GeV");
    if (XMIN !=0 or XMAX != 0)
    {   
        cout << "heeeeeeeee" << endl;
        histData->GetXaxis()->SetRangeUser(XMIN, XMAX);
         histData2->GetXaxis()->SetRangeUser(XMIN, XMAX);
            histData3->GetXaxis()->SetRangeUser(XMIN, XMAX);
    }

        histData->GetXaxis()->SetNdivisions(5);
     histData->GetYaxis()->SetNdivisions(6);

 TGaxis::SetMaxDigits(5);

    TString pdfFilename(histogramName);
    pdfFilename += ".pdf";

    canvas->SaveAs(pdfFilename);
    delete canvas;
    fileData->Close();
}





void SaveData( StUPCEvent *upcEvt, Long64_t &numsa, TLorentzVector &proton1, TLorentzVector &proton2)
{     

    if (numsa > 1300 and numsa < 1500)

    {
           vector <vector <double>> dataGoodTof;
            vector <vector <double>> dataTof;
            vector <vector <double>> data;     

            for (int i = 0; i < upcEvt->getNumberOfTracks(); i++)
            {
          
                StUPCTrack *track = upcEvt->getTrack(i);
                vector <double> dGoodTof;
                vector <double> dTof;
                vector <double> d;

                if (track->getFlag(StUPCTrack::kTof)  and  track->getPt() >= ptCutVal and  (abs(track->getEta()) <= etaCutVal) and (track->getNhitsFit() >= NhistFitCut))
                {
                 
                    dGoodTof.push_back(track->getEta());
                    dGoodTof.push_back(track->getPhi());
                    dGoodTof.push_back(track->getPt());
    
                    dataGoodTof.push_back(dGoodTof);

                }

                else if (track->getFlag(StUPCTrack::kTof) )
                {
                    dTof.push_back(track->getEta());
                    dTof.push_back(track->getPhi());
                    dTof.push_back(track->getPt());
                    dataTof.push_back(dTof);                
                }
                else
                {
                    d.push_back(track->getEta());
                    d.push_back(track->getPhi());
                    d.push_back(track->getPt());
                    data.push_back(d);
                }

            }

 
        
         //   WriteCSV_file("goodTof"+to_string(numsa)+".csv" ,dataGoodTof);
       //     WriteCSV_file("tof"+to_string(numsa)+".csv" ,dataTof);
          //  WriteCSV_file("otherTpc"+to_string(numsa)+".csv" ,data);         

            vector <vector<double>> prot;
            vector <double > p1;
            vector <double > p2;
            p1.push_back(proton1.Eta());
            p1.push_back(proton1.Phi());
            p1.push_back(proton1.Pt());            
              
            p2.push_back(proton2.Eta());
            p2.push_back(proton2.Phi());
            p2.push_back(proton2.Pt());   
 
            prot.push_back(p1);
            prot.push_back(p2);
            
           //WriteCSV_file("pro"+to_string(numsa)+".csv" ,prot);
           
    }
}

    
