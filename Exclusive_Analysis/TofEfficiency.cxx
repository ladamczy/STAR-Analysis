#include <iostream>  //tag and probe for true
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

// Function prototype for finding protons
void FindProtons(bool isMC, StRPEvent *rpEvt, StUPCEvent *upcEvt, TLorentzVector & proton1, TLorentzVector & proton2);

int main(int argc, char** argv)  
{
    // Particle masses and beam energy definitions
    float massPion = 0.13957061; // Mass of pion in GeV
    double massKaon = 497.611 / 1000.0; // Mass of kaon in GeV
    double beamEnergy = 254.867; // Beam energy in GeV
    double massProton = 938.272 / 1000.0; // Mass of proton in GeV

    //cut constants
    const double kDeltaRCut = 0.15;  //0.15

    // Open input file containing paths to data files
    ifstream inputFilePathList(argv[1]);
    if (!inputFilePathList) 
    {
        cerr << "Failed to open input file." << std::endl;
        return 1;
    }

    // Create TChain to read UPC tree from multiple files
    TChain *chain = new TChain("mUPCTree"); 
    string inputFileName;
    while (std::getline(inputFilePathList, inputFileName))
    {
        chain->Add(inputFileName.c_str());
    }
    inputFilePathList.close();

    
    // Initialize variables for event processing
    bool isMC = 0; // Flag to check if data is Monte Carlo (MC)
    static StUPCEvent *upcEvt = nullptr; // Pointer to UPC event
    static StRPEvent *correctedRpEvent = nullptr; // Pointer to corrected RP event

    // Load offset file for corrections (if applicable)
    //LoadOffsetFile("../data/OffSetsCorrectionsRun17.list", mCorrection);

    // Set branch addresses for UPC event and check if it's MC
    chain->SetBranchAddress("mUPCEvent", &upcEvt);
    chain->GetEntry(0);

    if (upcEvt->getRunNumber() == 1)
    {
        isMC = 1;
    }
    // Set branch address for corrected RP event if not MC
    if (isMC == 0)
    {
        chain->SetBranchAddress("correctedRpEvent", &correctedRpEvent);
    }

    // Define trigger IDs and fill numbers for analysis
    vector <int> triggerID = {570209, 570219, 570229, 570701, 570702, 570703, 570704, 570705, 570709, 570711, 570712, 570719, 590701, 590703, 590705, 590708, 590709};
    vector <int> triggerCEP = {570701, 570705, 570711};
    vector <int> fillMC = {20667,20704,20705,20668,20706,20707,20669,20708,20719,20720,20728,20733,20734,20735,20736,20670,20737,20738,20740,20741,20743,20745,20746,20747,20748,20751,20752,20756,20757,20758,20759,20760,20764,20765,20766,20768,20770,20771,20772,20773,20777,20671,20778,20780,20781,20783,20787,20788,20789,20790,20791,20799,20674,20800,20810,20811,20814,20815,20821,20823,20825,20827,20828,20835,20841,20845,20852,20675,20854,20856,20857,20859,20860,20861,20862,20863,20864,20865,20866,20867,20868,20870,20873,20874,20896,20897,20898,20900,20901,20902,20903,20680,20904,20906,20907,20910,20911,20920,20681,20924,20925,20926,20927,20928,20929,20682,20684,20687,20688,20692,20695,20697,20701,20702,20703};
    vector <int> fillRange;
    
     // Sort fill numbers and create a range for grouping fills
    sort(fillMC.begin(), fillMC.end());

    //for (int i = 0; i < fillMC.size(); i++)
    //{
        //cout << fillMC[i] << ", " ;
    //}
    //cout << " \n==================\n";

    for (int i = 0; i < fillMC.size(); i++)
    {
        if (i%23 == 0)
        {
          fillRange.push_back(i);
        }
        
    }
    fillRange.push_back(fillMC.size()-1);

     //for (int i = 0; i < fillRange.size(); i++)
    //{
        //cout << fillMC[fillRange[i]] << ", " ;
    //}

    // Define binning for histograms
    double pTmin = 0.0; double pTmax = 4.0; int nBinsPt =7;
    double etaMin = -3.0; double etaMax = 3.0; int nBinsEta = 7;
    double zMin = -80.0; double zMax = 80.0; int nBinsZ = 17;
    int nBinsFill = 7;
    
    double binsPt[7] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.2};
    double binsEta[7] = {-0.9, -0.6, -0.3, 0.0, 0.3, 0.6, 0.9};
     
    //the code explicitly pairs tracks with opposite charges, hence change PosProbe to BothProbe: just a suggestion

    TH1D* HistKaonMassProbeWithoutTof = new TH1D("HistKaonMassProbeWithoutTof", "; m_{#pi^{+}#pi^{-}}^{tag} [GeV]; # events",60 ,0.44, 0.56);  
    TH1D* HistKaonMassProbeWithTof = new TH1D("HistKaonMassProbeWithTof", "; m_{#pi^{+}#pi^{-}}^{probe} [GeV]; # events", 60 ,0.44, 0.56);

    TH1D* HistKaonMassProbeWithoutTofPosProbe = new TH1D("HistKaonMassProbeWithoutTofPosProbe", "; m_{#pi^{+}#pi^{-}} [GeV]; events", 60, 0.44, 0.56);
    //TH1D* HistKaonMassProbeWithoutTofNegProbe = new TH1D("HistKaonMassProbeWithoutTofNegProbe", "; m_{#pi^{+}#pi^{-}} [GeV]; events", 60, 0.44, 0.56);
    TH1D* HistKaonMassProbeWithTofPosProbe = new TH1D("HistKaonMassProbeWithTofPosProbe", "; m_{#pi^{+}#pi^{-}} [GeV]; events", 60, 0.44, 0.56);
    //TH1D* HistKaonMassProbeWithTofNegProbe = new TH1D("HistKaonMassProbeWithTofNegProbe", "; m_{#pi^{+}#pi^{-}} [GeV]; events", 60, 0.44, 0.56);
 
    TH1D* HistKaonPtProbeWithoutTofPosProbe = new TH1D("HistKaonPtProbeWithoutTofPosProbe", "; p_{T} [GeV]; events",nBinsPt-1, binsPt);
    //TH1D* HistKaonPtProbeWithoutTofNegProbe = new TH1D("HistKaonPtProbeWithoutTofNegProbe", "; P_{T} [GeV]; events", nBinsPt-1, binsPt);
    TH1D* HistKaonPtProbeWithTofPosProbe = new TH1D("HistKaonPtProbeWithTofPosProbe", "; p_{T} [GeV]; events", nBinsPt-1, binsPt);
    //TH1D* HistKaonPtProbeWithTofNegProbe = new TH1D("HistKaonPtProbeWithTofNegProbe", "; p_{T} [GeV]; events", nBinsPt-1, binsPt);
 
    TH1D* HistKaonEtaProbeWithoutTofPosProbe = new TH1D("HistKaonEtaProbeWithoutTofPosProbe", "positive probe w/o TOF; #eta ; events",nBinsEta-1, binsEta);
    //TH1D* HistKaonEtaProbeWithoutTofNegProbe = new TH1D("HistKaonEtaProbeWithoutTofNegProbe", "negative probe w/o TOF; #eta ; events",nBinsEta-1, binsEta);
    TH1D* HistKaonEtaProbeWithTofPosProbe = new TH1D("HistKaonEtaProbeWithTofPosProbe", "positive probe w/ TOF; #eta ; events",nBinsEta-1, binsEta);
    //TH1D* HistKaonEtaProbeWithTofNegProbe = new TH1D("HistKaonEtaProbeWithTofNegProbe", "negative probe w/ TOF; #eta ; events",nBinsEta-1, binsEta);


    vector <TH1D*> HistKaonPtProbeWithoutTofPosProbeMass;
    //vector <TH1D*> HistKaonPtProbeWithoutTofNegProbeMass;
    vector <TH1D*> HistKaonPtProbeWithTofPosProbeMass;
    //vector <TH1D*> HistKaonPtProbeWithTofNegProbeMass;

    vector <TH1D*> HistKaonEtaProbeWithoutTofPosProbeMass;
    //vector <TH1D*> HistKaonEtaProbeWithoutTofNegProbeMass;
    vector <TH1D*> HistKaonEtaProbeWithTofPosProbeMass;
    //vector <TH1D*> HistKaonEtaProbeWithTofNegProbeMass;


    vector <TH1D*> HistKaonFillProbeWithoutTofPosProbeMass;
    //vector <TH1D*> HistKaonFillProbeWithoutTofNegProbeMass;
    vector <TH1D*> HistKaonFillProbeWithTofPosProbeMass;
    //vector <TH1D*> HistKaonFillProbeWithTofNegProbeMass;

    TH1D* HistKaonPtProbeWithoutTofPosProbeMassMC = new TH1D("HistKaonPtProbeWithoutTofPosProbeMassMC", "", 60, 0.44, 0.56);
    //TH1D* HistKaonPtProbeWithoutTofNegProbeMassMC = new TH1D("HistKaonPtProbeWithoutTofNegProbeMassMC","",  60, 0.44, 0.56);
    TH1D* HistKaonPtProbeWithTofPosProbeMassMC = new TH1D("HistKaonPtProbeWithTofPosProbeMassMC","", 60, 0.44, 0.56);
    //TH1D* HistKaonPtProbeWithTofNegProbeMassMC = new TH1D("HistKaonPtProbeWithTofNegProbeMassMC", "",60, 0.44, 0.56);

    TH1D* HistKaonEtaProbeWithoutTofPosProbeMassMC = new TH1D("HistKaonEtaProbeWithoutTofPosProbeMassMC", "", 60, 0.44, 0.56);
    //TH1D* HistKaonEtaProbeWithoutTofNegProbeMassMC = new TH1D("HistKaonEtaProbeWithoutTofNegProbeMassMC","",  60, 0.44, 0.56);
    TH1D* HistKaonEtaProbeWithTofPosProbeMassMC = new TH1D("HistKaonEtaProbeWithTofPosProbeMassMC","", 60, 0.44, 0.56);
    //TH1D* HistKaonEtaProbeWithTofNegProbeMassMC = new TH1D("HistKaonEtaProbeWithTofNegProbeMassMC", "",60, 0.44, 0.56);

    TH1D* HistPionPtWithTofMC = new TH1D("HistPionPtWithTofMC", "; p_{T} [GeV]; events",nBinsPt-1, binsPt);
    TH1D* HistPionPtWithoutTofMC = new TH1D("HistPionPtWithoutTofMC", "; P_{T} [GeV]; events", nBinsPt-1, binsPt);

    TH1D* HistPionEtaWithTofMC = new TH1D("HistPionEtaWithTofMC", "; #eta ; events",nBinsEta-1, binsEta);
    TH1D* HistPionEtaWithoutTofMC = new TH1D("HistPionEtaWithoutTofMC", "; #eta ; events",nBinsEta-1, binsEta);

    // Create histograms for direct efficiency measurement
    TH1D* HistTruePionPt = new TH1D("HistTruePionPt", "; p_{T} [GeV]; All True Pions", nBinsPt-1, binsPt);
    TH1D* HistTruePionEta = new TH1D("HistTruePionEta", "; #eta ; All True Pions", nBinsEta-1, binsEta);
    TH1D* HistTruePionWithTofPt = new TH1D("HistTruePionWithTofPt", "; p_{T} [GeV]; True Pions with TOF", nBinsPt-1, binsPt);
    TH1D* HistTruePionWithTofEta = new TH1D("HistTruePionWithTofEta", "; #eta ; True Pions with TOF", nBinsEta-1, binsEta);

    
    for (int i = 0; i < nBinsPt-1; ++i) 
    {
        TH1D* histPosWo = new TH1D(Form("HistKaonPtProbeWithoutTofPosProbeMass%d", i), Form("HistKaonPtProbeWithoutTofPosProbeMass %d", i), 60, 0.44, 0.56);
        //TH1D* histNegWo = new TH1D(Form("HistKaonPtProbeWithoutTofNegProbeMass%d", i), Form("HistKaonPtProbeWithoutTofNegProbeMass %d", i), 60, 0.44, 0.56);
        TH1D* histPosW = new TH1D(Form("HistKaonPtProbeWithTofPosProbeMass%d", i), Form("HistKaonPtProbeWithTofPosProbeMass %d", i), 60, 0.44, 0.56);
        //TH1D* histNegW = new TH1D(Form("HistKaonPtProbeWithTofNegProbeMass%d", i), Form("HistKaonPtProbeWithTofNegProbeMass %d", i), 60, 0.44, 0.56);
        HistKaonPtProbeWithoutTofPosProbeMass.push_back(histPosWo);
        //HistKaonPtProbeWithoutTofNegProbeMass.push_back(histNegWo);
        HistKaonPtProbeWithTofPosProbeMass.push_back(histPosW);
        //HistKaonPtProbeWithTofNegProbeMass.push_back(histNegW);
    }

    for (int i = 0; i < nBinsEta-1; ++i) 
    {
        TH1D* histPosWo = new TH1D(Form("HistKaonEtaProbeWithoutTofPosProbeMass%d", i), Form("HistKaonEtaProbeWithoutTofPosProbeMass %d", i), 60, 0.44, 0.56);
        //TH1D* histNegWo = new TH1D(Form("HistKaonEtaProbeWithoutTofNegProbeMass%d", i), Form("HistKaonEtaProbeWithoutTofNegProbeMass %d", i), 60, 0.44, 0.56);
        TH1D* histPosW = new TH1D(Form("HistKaonEtaProbeWithTofPosProbeMass%d", i), Form("HistKaonEtaProbeWithTofPosProbeMass %d", i), 60, 0.44, 0.56);
        //TH1D* histNegW = new TH1D(Form("HistKaonEtaProbeWithTofNegProbeMass%d", i), Form("HistKaonEtaProbeWithTofNegProbeMass %d", i), 60, 0.44, 0.56);
        HistKaonEtaProbeWithoutTofPosProbeMass.push_back(histPosWo);
        //HistKaonEtaProbeWithoutTofNegProbeMass.push_back(histNegWo);
        HistKaonEtaProbeWithTofPosProbeMass.push_back(histPosW);
        //HistKaonEtaProbeWithTofNegProbeMass.push_back(histNegW);
    }

    
    for (int i = 0; i < nBinsFill; ++i) 
    {
        TH1D* histPosWo = new TH1D(Form("HistKaonFillProbeWithoutTofPosProbeMass%d", i), Form("HistKaonFillProbeWithoutTofPosProbeMass %d", i), 60, 0.44, 0.56);
        //TH1D* histNegWo = new TH1D(Form("HistKaonFillProbeWithoutTofNegProbeMass%d", i), Form("HistKaonFillProbeWithoutTofNegProbeMass %d", i), 60, 0.44, 0.56);
        TH1D* histPosW = new TH1D(Form("HistKaonFillProbeWithTofPosProbeMass%d", i), Form("HistKaonFillProbeWithTofPosProbeMass %d", i), 60, 0.44, 0.56);
        //TH1D* histNegW = new TH1D(Form("HistKaonFillProbeWithTofNegProbeMass%d", i), Form("HistKaonFillProbeWithTofNegProbeMass %d", i), 60, 0.44, 0.56);
        HistKaonFillProbeWithoutTofPosProbeMass.push_back(histPosWo);
        //HistKaonFillProbeWithoutTofNegProbeMass.push_back(histNegWo);
        HistKaonFillProbeWithTofPosProbeMass.push_back(histPosW);
        //HistKaonFillProbeWithTofNegProbeMass.push_back(histNegW);
    }
   
    // Histograms for truth-matched tag and probe analysis
    TH1D* HistKaonMassProbeWithTof_TruePions = new TH1D("HistKaonMassProbeWithTof_TruePions", "; m_{#pi^{+}#pi^{-}}^{probe} [GeV]; # events", 60, 0.44, 0.56);
    TH1D* HistKaonMassProbeWithoutTof_TruePions = new TH1D("HistKaonMassProbeWithoutTof_TruePions", "; m_{#pi^{+}#pi^{-}}^{tag} [GeV]; # events", 60, 0.44, 0.56);

    TH1D* HistPionPtWithTof_TruePions = new TH1D("HistPionPtWithTof_TruePions", "; p_{T} [GeV]; true pions with TOF", nBinsPt-1, binsPt);
    TH1D* HistPionPtWithoutTof_TruePions = new TH1D("HistPionPtWithoutTof_TruePions", "; p_{T} [GeV]; true pions without TOF", nBinsPt-1, binsPt);

    TH1D* HistPionEtaWithTof_TruePions = new TH1D("HistPionEtaWithTof_TruePions", "; #eta; true pions with TOF", nBinsEta-1, binsEta);
    TH1D* HistPionEtaWithoutTof_TruePions = new TH1D("HistPionEtaWithoutTof_TruePions", "; #eta; true pions without TOF", nBinsEta-1, binsEta);

    // True pion probe mass histograms (binned in pT and η)
    vector<TH1D*> HistKaonPtProbeWithTof_TruePionsMass;
    vector<TH1D*> HistKaonPtProbeWithoutTof_TruePionsMass;
    vector<TH1D*> HistKaonEtaProbeWithTof_TruePionsMass;
    vector<TH1D*> HistKaonEtaProbeWithoutTof_TruePionsMass;

    //histograms for true pions analysis
    // Histograms to store kaon mass for true pions in different pT bins
    vector <TH1D*> HistKaonPtTruePionWithTofMass;
    vector <TH1D*> HistKaonPtTruePionWithoutTofMass;

    // Histograms to store kaon mass for true pions in different eta bins
    vector <TH1D*> HistKaonEtaTruePionWithTofMass;
    vector <TH1D*> HistKaonEtaTruePionWithoutTofMass;

    // Initialize the histograms
    for (int i = 0; i < nBinsPt-1; ++i) {
        TH1D* histWithTof = new TH1D(Form("HistKaonPtTruePionWithTofMass%d", i), 
                                Form("Kaon Mass (True Pions) with TOF, %.2f < p_{T} < %.2f GeV", binsPt[i], binsPt[i+1]), 
                                60, 0.44, 0.56);
        TH1D* histWithoutTof = new TH1D(Form("HistKaonPtTruePionWithoutTofMass%d", i), 
                                    Form("Kaon Mass (True Pions) without TOF, %.2f < p_{T} < %.2f GeV", binsPt[i], binsPt[i+1]), 
                                    60, 0.44, 0.56);
        HistKaonPtTruePionWithTofMass.push_back(histWithTof);
        HistKaonPtTruePionWithoutTofMass.push_back(histWithoutTof);
    }

    for (int i = 0; i < nBinsEta-1; ++i) {
        TH1D* histWithTof = new TH1D(Form("HistKaonEtaTruePionWithTofMass%d", i), 
                                Form("Kaon Mass (True Pions) with TOF, %.2f < #eta < %.2f", binsEta[i], binsEta[i+1]), 
                                60, 0.44, 0.56);
        TH1D* histWithoutTof = new TH1D(Form("HistKaonEtaTruePionWithoutTofMass%d", i), 
                                    Form("Kaon Mass (True Pions) without TOF, %.2f < #eta < %.2f", binsEta[i], binsEta[i+1]), 
                                    60, 0.44, 0.56);
        HistKaonEtaTruePionWithTofMass.push_back(histWithTof);
        HistKaonEtaTruePionWithoutTofMass.push_back(histWithoutTof);
    }

    
    //cut flow histograms
    TH1D* HistEventCutFlow = new TH1D("HistEventCutFlow", ";Cut;Events", 5, 0, 5);
    HistEventCutFlow->GetXaxis()->SetBinLabel(1, "All Events");
    HistEventCutFlow->GetXaxis()->SetBinLabel(2, "Trigger Cuts");
    HistEventCutFlow->GetXaxis()->SetBinLabel(3, "Vertex Cuts");
    HistEventCutFlow->GetXaxis()->SetBinLabel(4, ">= 3 TOF Tracks");
    HistEventCutFlow->GetXaxis()->SetBinLabel(5, "Final Selection");

    TH1D* HistTrackCutFlow = new TH1D("HistTrackCutFlow", ";Cut;Tracks", 6, 0, 6);
    HistTrackCutFlow->GetXaxis()->SetBinLabel(1, "All Tracks");
    HistTrackCutFlow->GetXaxis()->SetBinLabel(2, "|#eta| < 0.9");
    HistTrackCutFlow->GetXaxis()->SetBinLabel(3, "p_{T} > 0.2 GeV");
    HistTrackCutFlow->GetXaxis()->SetBinLabel(4, "NhitsFit >= 20");
    HistTrackCutFlow->GetXaxis()->SetBinLabel(5, "kV0 Flag");
    //HistTrackCutFlow->GetXaxis()->SetBinLabel(6, "TOF Matched");

    // Loop over all entries in the TChain
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

        // Fill event cut flow: All events
        HistEventCutFlow->Fill(0);

        TLorentzVector lorentzVector,trackVector;
        TLorentzVector proton1, proton2;

        // Get fill number and determine the corresponding group index        
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
 
        // Select tracks based on quality criteria
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
                    if(upcEvt->getTrack(i)->getFlag(StUPCTrack::kTof))
                    {
                        tracksWithTofHit.push_back(upcEvt->getTrack(i));
                    } 
                }
            }
            else if (isMC == 1)
            {
                if ( (abs(upcEvt->getTrack(i)->getEta()) <= 0.9) and (upcEvt->getTrack(i)->getPt() >= 0.2) and (upcEvt->getTrack(i)->getNhitsFit() >= 20)) 
                {
                    tracksWithSmallDca.push_back(upcEvt->getTrack(i));
                    tracksForTrueLevel.push_back(upcEvt->getTrack(i));
                    if(upcEvt->getTrack(i)->getFlag(StUPCTrack::kTof))
                    {
                        tracksWithTofHit.push_back(upcEvt->getTrack(i));
                    } 
                }
            }
		} 
        // ============================
        // True Level Analysis (MC Only)
        // ============================
        if (isMC==1) 
        {
            //cout << "Processing MC truth section" << endl;
            //cout << "Number of MC particles: " << upcEvt->getNumberOfMCParticles() << endl;

           
            // Extract true pions from MC truth information
            TParticle* particle;
            vector <TParticle *> vPions; // Vector to store true pions
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
            // Loop over true pions and match them to reconstructed tracks
            for (int j = 0; j < vPions.size(); j++) 
            {
                 // Apply vertex and kinematic cuts 
                lorentzVector.SetPxPyPzE(vPions[j]->Px(), vPions[j]->Py(), vPions[j]->Pz(), vPions[j]->Energy()); 
                double pT = sqrt(pow(vPions[j]->Px(),2) + pow(vPions[j]->Py(),2)) ;
                double eta = lorentzVector.Eta();
                vPions[j]->ProductionVertex(productionVertex);

                double truthVertexR = sqrt(pow(productionVertex.X(),2) + pow(productionVertex.Y(),2)); // Transverse vertex position
                double truthVertexZ = productionVertex.Z(); // Longitudinal vertex position

                //cout << "Pion " << j << ": pT=" << pT << ", eta=" << eta 
                //<< ", vtxR=" << truthVertexR << ", vtxZ=" << truthVertexZ << endl;

                // Apply cuts on the production vertex
                if (abs(truthVertexZ) > 80 || truthVertexR > 3.0 || pT < 0.2 || abs(eta) > 0.9)
                {
                    continue; // Skip pions outside the acceptance
                }
                // Fill denominators
                HistTruePionPt->Fill(pT);
                HistTruePionEta->Fill(eta);
                
                // Match reconstructed tracks to true pions
                vector <double> distr;
                vector <StUPCTrack const *> matchedTPCpions;
                if (tracksForTrueLevel.size() == 0) {break;}
                for (int i = 0; i < tracksForTrueLevel.size(); i++) // pętla po śladach globalnych dobrej jakości
                {
                    TLorentzVector pion, track;
                    pion.SetPxPyPzE(vPions[j]->Px(), vPions[j]->Py(), vPions[j]->Pz(), vPions[j]->Energy());
                    tracksForTrueLevel[i]->getLorentzVector(track, massPion);
                    double r = pion.DeltaR(track);  // Spatial distance in η-φ space
                    matchedTPCpions.push_back(tracksForTrueLevel[i]); // tracki z tracksForTrueLevel są usuwane poniżej po znalezieniu najmniejszego r - dlatego zbieram jeszcze raz ślady do wektora
                    distr.push_back(r); // dla dnaego pionu z piomu true obliczenie wszsystkich możliwych odległości R dla tracków dobrej jakości tracksForTrueLevel
                }
   
                // Find the minimum distance and corresponding index
                auto it = std::min_element(distr.begin(), distr.end()); //znajdź element o najmniejszym R
                int index = std::distance(distr.begin(), it); // odpowiedni indeks
                //cout << "Best match ΔR: " << distr[index] << endl;                

                if (distr[index] < kDeltaRCut) // jeżeli odległość mniejsza niż .15 - akceptujemy przypadek
                {  
                    if (matchedTPCpions[index]->getFlag(StUPCTrack::kTof)) // podział TOF nie TOF
                    {
                        HistPionPtWithTofMC->Fill(pT);
                        HistPionEtaWithTofMC->Fill(eta);
                        // Fill numerators
                        HistTruePionWithTofPt->Fill(pT);
                        HistTruePionWithTofEta->Fill(eta);
                    }

                    else
                    {
                        HistPionPtWithoutTofMC->Fill(pT);
                        HistPionEtaWithoutTofMC->Fill(eta);
                    } 
                    // Remove the matched track to avoid double-counting
                    tracksForTrueLevel.erase(tracksForTrueLevel.begin()+index); // usunięcie z tracksForTrueLevel śladu
                }    
            }            
              
        }

        // Skip events with fewer than 3 TOF-matched tracks
		if (tracksWithTofHit.size() < 3) // rozpatrujemy przypadki gdzie są przynajmniej trzy ślady zmatchowane z TOF
		{
			continue;
		}
        
        // Fill event cut flow: >= 3 TOF-matched tracks
        HistEventCutFlow->Fill(3);

        //BBC cuts
        /*Int_t BBCLEast = upcEvt->getBBCLargeEast();
        Int_t BBCLWest = upcEvt->getBBCLargeWest();

        if (BBCLWest > 52 or BBCLEast > 52.0)
        {
            continue;
        }*/
        // ============================
        // Tag-and-Probe Analysis + matching tracks with true MC
        // ============================

        StUPCTrack const * track1;
        StUPCTrack const * track2;

        // For MC data, extract true pions
        vector<TParticle*> vTruePions;
        vector<TLorentzVector> vTruePionVectors;
        
        // For MC data, check if tracks are true pions
        if (isMC == 1) {
            // Extract and filter true pions
            for (int k = 0; k < upcEvt->getNumberOfMCParticles(); k++) {
                TParticle* particle = upcEvt->getMCParticle(k);
                if (particle->GetPDG()->PdgCode() == 211 || particle->GetPDG()->PdgCode() == -211) {
                    // Apply the same vertex and kinematic cuts as in the truth analysis
                    TLorentzVector pion;
                    pion.SetPxPyPzE(particle->Px(), particle->Py(), particle->Pz(), particle->Energy());
                    double pT = sqrt(pow(particle->Px(), 2) + pow(particle->Py(), 2));
                    double eta = pion.Eta();
                    
                    TLorentzVector productionVertex;
                    particle->ProductionVertex(productionVertex);
                    double truthVertexR = sqrt(pow(productionVertex.X(), 2) + pow(productionVertex.Y(), 2));
                    double truthVertexZ = productionVertex.Z();
                    
                    // Apply cuts on the production vertex
                    if (abs(truthVertexZ) <= 80 && truthVertexR <= 3.0 && pT >= 0.2 && abs(eta) <= 0.9) {
                        vTruePions.push_back(particle);
                        vTruePionVectors.push_back(pion);
                    }
                }
            }
        }

        // Pair tracks and apply additional cuts
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

                    // Check if the tracks are TOF-matched with true pions
                    bool isTrack1TruePion = false;
                    bool isTrack2TruePion = false;

                    // For MC data, check if tracks are true pions using the same matching approach as the truth analysis
                    if (isMC == 1 && !vTruePions.empty()) {
                        // Match track1 to true pions
                        vector<double> distr1;
                        TLorentzVector track1Vector;
                        track1->getLorentzVector(track1Vector, massPion);
                        
                        for (const auto& truePion : vTruePionVectors) {
                            double r = track1Vector.DeltaR(truePion);
                            distr1.push_back(r);
                        }
                        
                        if (!distr1.empty()) {
                            auto it1 = std::min_element(distr1.begin(), distr1.end());
                            int index1 = std::distance(distr1.begin(), it1);
                            if (distr1[index1] < kDeltaRCut) {
                                isTrack1TruePion = true;
                            }
                        }
                        
                        // Match track2 to true pions
                        vector<double> distr2;
                        TLorentzVector track2Vector;
                        track2->getLorentzVector(track2Vector, massPion);
                        
                        for (const auto& truePion : vTruePionVectors) {
                            double r = track2Vector.DeltaR(truePion);
                            distr2.push_back(r);
                        }
                        
                        if (!distr2.empty()) {
                            auto it2 = std::min_element(distr2.begin(), distr2.end());
                            int index2 = std::distance(distr2.begin(), it2);
                            if (distr2[index2] < kDeltaRCut) {
                                isTrack2TruePion = true;
                            }
                        }
                    }
                    
                    // Check if tracks have TOF hits
                    if (track1->getFlag(StUPCTrack::kTof))
                    {
                        hasTofHitTrack1 = true;
                    }

                    if (track2->getFlag(StUPCTrack::kTof))
                    {
                        hasTofHitTrack2 = true;
                    }
                     // Skip pairs without at least one TOF hit 
                    if (hasTofHitTrack1 == false and hasTofHitTrack2 == false)
                    {
                        continue;
                    }
                    // Skip pairs without opposite charges
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

                    // Apply V0 reconstruction and decay vertex cuts                    
                    TVector3 const tryVec(0,0,0);
                    StUPCV0 kaon(track1, track2, massPion, massPion, 1, 1, tryVec,beamPar, upcEvt->getMagneticField(), true);
                    double dcaDau = kaon.dcaDaughters();
                    double dcaBeam = kaon.DCABeamLine();
                    double hypoAngle = kaon.pointingAngleHypo();
                    double d = kaon.decayLengthHypo();
                    double vtxZ = kaon.decayVertex().Z();
                    double vtxR = sqrt(pow(kaon.decayVertex().X(), 2) + pow(kaon.decayVertex().Y(), 2));
                    // Cut on decay vertex position
                    if (abs(vtxZ) < 80 && vtxR < 3)
                    {
                        // Further cuts on DCA and TOF matching 
                        if (dcaDau < 1 and dcaBeam < 1 )// and  (kaon.pointingAngleHypo()>0.925 )  )
                        {   
                            //TOF matching 
                            if (int(hasTofHitTrack1) + int(hasTofHitTrack2) == 2)
                            {  
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
                                
                                // If both tracks are true pions, fill the truth-matched histograms
                                if (isMC == 1 && isTrack1TruePion && isTrack2TruePion) {
                                    HistKaonMassProbeWithTof_TruePions->Fill(kaon.m());
                                    
                                    // Fill pT and eta histograms for true pions with TOF
                                    HistPionPtWithTof_TruePions->Fill(track1->getPt());
                                    HistPionPtWithTof_TruePions->Fill(track2->getPt());
                                    HistPionEtaWithTof_TruePions->Fill(track1->getEta());
                                    HistPionEtaWithTof_TruePions->Fill(track2->getEta());

                                    // Fill pT-binned histograms
                                    for (int i = 0; i < nBinsPt-1; i++) {
                                        if (track1->getPt() >= binsPt[i] && track1->getPt() < binsPt[i+1]) {
                                            HistKaonPtTruePionWithTofMass[i]->Fill(kaon.m());
                                        }
                                        if (track2->getPt() >= binsPt[i] && track2->getPt() < binsPt[i+1]) {
                                            HistKaonPtTruePionWithTofMass[i]->Fill(kaon.m());
                                        }
                                    }
                                    
                                    // Fill eta-binned histograms
                                    for (int i = 0; i < nBinsEta-1; i++) {
                                        if (track1->getEta() >= binsEta[i] && track1->getEta() < binsEta[i+1]) {
                                            HistKaonEtaTruePionWithTofMass[i]->Fill(kaon.m());
                                        }
                                        if (track2->getEta() >= binsEta[i] && track2->getEta() < binsEta[i+1]) {
                                            HistKaonEtaTruePionWithTofMass[i]->Fill(kaon.m());
                                        }
                                    }
                                }   
                            }
                           
                            else     // without ToF
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

                                // Check if tag and probe are true pions (for MC only)
                                if (isMC == 1 && int(hasTofHitTrack1) + int(hasTofHitTrack2) == 1) {
                                    StUPCTrack const * tagTrack;
                                    StUPCTrack const * probeTrack;
                                    bool isTagTruePion = false;
                                    bool isProbeTruePion = false;
                                    
                                    if (hasTofHitTrack1) {
                                        tagTrack = track1;
                                        probeTrack = track2;
                                        isTagTruePion = isTrack1TruePion;
                                        isProbeTruePion = isTrack2TruePion;
                                    } else {
                                        tagTrack = track2;
                                        probeTrack = track1;
                                        isTagTruePion = isTrack2TruePion;
                                        isProbeTruePion = isTrack1TruePion;
                                    }
                                    if (isTagTruePion && isProbeTruePion) {
                                        HistKaonMassProbeWithoutTof_TruePions->Fill(kaon.m());
                                        
                                        // Fill pT and eta histograms for the probe track (without TOF) if it's a true pion
                                        HistPionPtWithoutTof_TruePions->Fill(probeWoTof->getPt());
                                        HistPionEtaWithoutTof_TruePions->Fill(probeWoTof->getEta());
                                    
                                        // Fill pT-binned histograms for the probe track
                                        for (int i = 0; i < nBinsPt-1; i++) {
                                            if (probeTrack->getPt() >= binsPt[i] && probeTrack->getPt() < binsPt[i+1]) {
                                                HistKaonPtTruePionWithoutTofMass[i]->Fill(kaon.m());
                                            }
                                        }
                                        
                                        // Fill eta-binned histograms for the probe track
                                        for (int i = 0; i < nBinsEta-1; i++) {
                                            if (probeTrack->getEta() >= binsEta[i] && probeTrack->getEta() < binsEta[i+1]) {
                                                HistKaonEtaTruePionWithoutTofMass[i]->Fill(kaon.m());
                                            }
                                        }
                                    }
                                }                           
                            }
                        }
                    }
                }
            }
        }

         // Fill event cut flow: Final selection
         HistEventCutFlow->Fill(4);

         //just for trackflow
        for (Int_t j = 0; j < upcEvt->getNumberOfTracks(); j++)
        {
            StUPCTrack const* track = upcEvt->getTrack(j);

            // Track cut flow: All tracks
            HistTrackCutFlow->Fill(0);

            // Track cut flow: |eta| < 0.9
            if (abs(track->getEta()) < 0.9)
            {
                HistTrackCutFlow->Fill(1);
                // Track cut flow: pT > 0.2 GeV
                if (track->getPt() < 0.2)
                {
                    HistTrackCutFlow->Fill(2);
                    // Track cut flow: NhitsFit >= 20
                    if (track->getNhitsFit() < 20)
                    {
                        HistTrackCutFlow->Fill(3);
                        // Track cut flow: kV0 flag
                        if (!track->getFlag(StUPCTrack::kV0))
                        {
                            HistTrackCutFlow->Fill(4);
                            // Track cut flow: TOF matched
                            if (track->getFlag(StUPCTrack::kTof))
                            {
                                tracksWithTofHit.push_back(track);
                                HistTrackCutFlow->Fill(5);
                            }
                        }
                    }  
                }
            } 
          
          
        }
    }
         
    // Write histograms to output file
    TFile *outfile = TFile::Open(argv[2], "recreate");

    HistEventCutFlow->Write();
    HistTrackCutFlow->Write();

    HistKaonMassProbeWithTof->Write();
    HistKaonMassProbeWithoutTof->Write();
    
    HistKaonMassProbeWithoutTofPosProbe->Write();
    //HistKaonMassProbeWithoutTofNegProbe->Write();
    HistKaonMassProbeWithTofPosProbe->Write();
    //HistKaonMassProbeWithTofNegProbe->Write();

    HistKaonPtProbeWithoutTofPosProbe->Write();
    //HistKaonPtProbeWithoutTofNegProbe->Write();
    HistKaonPtProbeWithTofPosProbe->Write();
    //HistKaonPtProbeWithTofNegProbe->Write();
    HistKaonEtaProbeWithoutTofPosProbe->Write();
    //HistKaonEtaProbeWithoutTofNegProbe->Write();
    HistKaonEtaProbeWithTofPosProbe->Write();
    //HistKaonEtaProbeWithTofNegProbe->Write();

    //HistKaonPtProbeWithoutTofNegProbeMassMC->Write();
    HistKaonPtProbeWithoutTofPosProbeMassMC->Write();
    //HistKaonPtProbeWithTofNegProbeMassMC->Write();
    HistKaonPtProbeWithTofPosProbeMassMC->Write();

    //HistKaonEtaProbeWithoutTofNegProbeMassMC->Write();
    HistKaonEtaProbeWithoutTofPosProbeMassMC->Write();
    //HistKaonEtaProbeWithTofNegProbeMassMC->Write();
    HistKaonEtaProbeWithTofPosProbeMassMC->Write();/**/

    HistPionPtWithTofMC->Write();
    HistPionEtaWithTofMC->Write();
    HistPionPtWithoutTofMC->Write();
    HistPionEtaWithoutTofMC->Write();
        
    // calculate and write efficiency histograms for true pions in tag and probe
    TH1D* HistTagProbeEfficiencyVsPt = new TH1D("HistTagProbeEfficiencyVsPt", "; p_{T} [GeV]; Tag & Probe TOF Efficiency", nBinsPt-1, binsPt);
    HistTagProbeEfficiencyVsPt->Divide(HistPionPtWithTof_TruePions, HistPionPtWithoutTof_TruePions, 1, 1, "B");
    //HistTagProbeEfficiencyVsPt->Write();

    TH1D* HistTagProbeEfficiencyVsEta = new TH1D("HistTagProbeEfficiencyVsEta", "; #eta; Tag & Probe TOF Efficiency", nBinsEta-1, binsEta);
    HistTagProbeEfficiencyVsEta->Divide(HistPionEtaWithTof_TruePions, HistPionEtaWithoutTof_TruePions, 1, 1, "B");
    //HistTagProbeEfficiencyVsEta->Write();    

    for (int i = 0; i < HistKaonPtProbeWithoutTofPosProbeMass.size(); i++)
    {
        HistKaonPtProbeWithoutTofPosProbeMass[i]->Write();
        //HistKaonPtProbeWithoutTofNegProbeMass[i]->Write();
        HistKaonPtProbeWithTofPosProbeMass[i]->Write();
        //HistKaonPtProbeWithTofNegProbeMass[i]->Write();
    }
    for (int i = 0; i < HistKaonEtaProbeWithoutTofPosProbeMass.size(); i++)
    {
        HistKaonEtaProbeWithoutTofPosProbeMass[i]->Write();
        //HistKaonEtaProbeWithoutTofNegProbeMass[i]->Write();
        HistKaonEtaProbeWithTofPosProbeMass[i]->Write();
        //HistKaonEtaProbeWithTofNegProbeMass[i]->Write();
    }

    for (int i = 0; i < HistKaonFillProbeWithoutTofPosProbeMass.size(); i++)
    {
        HistKaonFillProbeWithoutTofPosProbeMass[i]->Write();
        //HistKaonFillProbeWithoutTofNegProbeMass[i]->Write();
        HistKaonFillProbeWithTofPosProbeMass[i]->Write();
        //HistKaonFillProbeWithTofNegProbeMass[i]->Write();
    }/**/

    HistKaonMassProbeWithTof_TruePions->Write();
    HistKaonMassProbeWithoutTof_TruePions->Write();
    HistPionPtWithTof_TruePions->Write();
    HistPionPtWithoutTof_TruePions->Write();
    HistPionEtaWithTof_TruePions->Write();
    HistPionEtaWithoutTof_TruePions->Write();

    HistTruePionPt->Write();
    HistTruePionEta->Write();
    HistTruePionWithTofPt->Write();
    HistTruePionWithTofEta->Write();
    // Write true pions histograms for different pT bins
    for (int i = 0; i < HistKaonPtTruePionWithTofMass.size(); i++) {
        HistKaonPtTruePionWithTofMass[i]->Write();
        HistKaonPtTruePionWithoutTofMass[i]->Write();
    }
    // Write true pions histograms for different eta bins
    for (int i = 0; i < HistKaonEtaTruePionWithTofMass.size(); i++) {
        HistKaonEtaTruePionWithTofMass[i]->Write();
        HistKaonEtaTruePionWithoutTofMass[i]->Write();
    }
    outfile->Close();

    return 0;
}


