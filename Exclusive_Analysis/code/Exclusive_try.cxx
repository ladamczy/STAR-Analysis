#include "ExclusiveCode.h"
using namespace std;



int main(int argc, char** argv)  
{
    TH1D* HistNumWithTofTrakcs = new TH1D("HistNumWithTofTrakcs", ";N^{tracks};events", 53, -0.5, 52.5);
    TH1D* HistNumWithoutTofTrakcs = new TH1D("HistNumWithoutTofTrakcs", ";N^{tracks};events", 154, -0.5 , 153.5);  

    vector <TH1D*> HistPtPionWithTof;
    vector <TH1D*> HistPtPionWithoutTof;
    vector <TH1D*> HistEtaPionWithTof;
    vector <TH1D*> HistEtaPionWithoutTof;
    vector <TH1D*> HistNfitPionWithTof;
    vector <TH1D*> HistNfitPionWithoutTof;


    TH1D* histEvtMult = new TH1D("histEvtMult", ";N^{tracks};events", 6, -0.5, 5.5);
        TH1D* histEvtMult2 = new TH1D("histEvtMult2", ";N^{tracks};events", 6, -0.5, 5.5);
    for (int i = 2; i < 6; ++i) 
    {
        TH1D* histPtWith = new TH1D(Form("histPtWith_%d", i), Form("histPtWith _%d", i), 50, 0, 5);
        TH1D* histPtWithout = new TH1D(Form("histPtWithout_%d", i), Form("histPtWithout _%d", i), 50, 0, 5);
        TH1D* histEtaWith = new TH1D(Form("histEtaWith_%d", i), Form("histEtaWith _%d", i), 101, -5, 5);
        TH1D* histEtaWithout = new TH1D(Form("histEtaWithout_%d", i), Form("histEtaWithout _%d", i), 101, -5, 5);
        TH1D* histNfitPionWithTof = new TH1D(Form("histNfitPionWithTof_%d", i), Form("histNfitPionWithTof _%d", i), 101, -0.5, 100.5);
        TH1D* histNfitPionWithoutTof = new TH1D(Form("histNfitPionWithoutTof%d", i), Form("histNfitPionWithoutTof _%d", i), 101, -0.5, 100.5);
        HistPtPionWithTof.push_back(histPtWith);
        HistPtPionWithoutTof.push_back(histPtWithout);
        HistEtaPionWithTof.push_back(histEtaWith);
        HistEtaPionWithoutTof.push_back(histEtaWithout);
        HistNfitPionWithTof.push_back(histNfitPionWithTof);
        HistNfitPionWithoutTof.push_back(histNfitPionWithoutTof);
    }

    vector <TH1D *> HistPtMissBefore;
    vector <TH1D*> HistNTOFClusterBefore;
    vector <TH1D*> HistDCABeamlineBefore;
    vector <TH1D*> HistDCADaughtersBefore;
    vector <TH1D*> HistCosBefore;   
    vector <TH1D*> HistDecayBefore;            
    vector <TH1D*> HistCorrEBefore;
    vector <TH1D*> HistCorrWBefore;
    vector <TH1D*> HistCorrKsiBefore;
    vector <TH1D*> HistCorrEtaBefore;
    vector <TH1D*> HistKsiEBefore;
    vector <TH1D*> HistKsiWBefore;
    vector <TH1D*> HistSumProtonMomentaXBefore;
    vector <TH1D*> HistSumProtonMomentaYBefore;
    vector <TH1D*> HistZDiffBefore;
    vector <TH1D*> HistZMeanBefore;
    vector <TH2D*> HistInvMassPiPi2D;
    vector <TH1D*> HistPhiBefore;
    vector <TH1D*> HistCosThetaStarBefore;

    //BEFORE
    for (int i = 2; i < 6; i++)
    {        
        TH2D* histInvMassPiPi2D = new TH2D(Form("histInvMassPiPi2D%d", i), Form("histInvMassPiPi2D _%d", i), 30 ,0.44, 0.56, 30 , 0.44, 0.56);
        HistInvMassPiPi2D.push_back(histInvMassPiPi2D);

        TH1D* histPtMissBefore = new TH1D(Form("histPtMissBefore%d", i), Form("histPtMissBefore _%d", i),10,0,0.5);
        HistPtMissBefore.push_back(histPtMissBefore);
    
        TH1D* histNTOFClusterBefore = new TH1D(Form("histNTOFClusterBefore%d", i), Form("histNTOFClusterBefore _%d", i),  16, -0.5, 15.5);
        HistNTOFClusterBefore.push_back(histNTOFClusterBefore);

        TH1D* histDCABeamlineBefore = new TH1D(Form("histDCABeamlineBefore%d", i), Form("histDCABeamlineBefore _%d", i), 10, 0, 5.0);
        HistDCABeamlineBefore.push_back(histDCABeamlineBefore);

        TH1D* histDCADaughtersBefore = new TH1D(Form("histDCADaughtersBefore%d", i), Form("histDCADaughtersBefore _%d", i), 40, 0,20.0);
        HistDCADaughtersBefore.push_back(histDCADaughtersBefore);

        TH1D* histCosBefore = new TH1D(Form("histCosBefore%d", i), Form("histCosBefore _%d", i), 20,0.9,1);
        HistCosBefore.push_back(histCosBefore);

        TH1D* histDecayBefore = new TH1D(Form("histDecayBefore%d", i), Form("histDecayBefore _%d", i),9,0, 9);
        HistDecayBefore.push_back(histDecayBefore);

        TH1D* histCorrEBefore = new TH1D(Form("histCorrEBefore%d", i), Form("histCorrEBefore _%d", i),30, -0.015, 0.015);
        HistCorrEBefore.push_back(histCorrEBefore);

        TH1D* histCorrWBefore = new TH1D(Form("histCorrWBefore%d", i), Form("histCorrWBefore _%d", i),30, -0.015, 0.015);
        HistCorrWBefore.push_back(histCorrWBefore);

        TH1D* histCorrKsiBefore = new TH1D(Form("histCorrKsiBefore%d", i), Form("histCorrKsiBefore _%d", i),   10, -0.005,0.005);
        HistCorrKsiBefore.push_back(histCorrKsiBefore);

        TH1D* histCorrEtaBefore = new TH1D(Form("histCorrEtaBefore%d", i), Form("histCorrEtaBefore _%d", i),30, -1.5,1.5);
        HistCorrEtaBefore.push_back(histCorrEtaBefore);

        TH1D* histKsiEBefore = new TH1D(Form("histKsiEBefore%d", i), Form("histKsiEBefore _%d", i),30, -0.015, 0.015);
        HistKsiEBefore.push_back(histKsiEBefore);

        TH1D* histKsiWBefore = new TH1D(Form("histKsiWBefore%d", i), Form("histKsiWBefore _%d", i), 30, -0.015, 0.015);
        HistKsiWBefore.push_back(histKsiWBefore);

        TH1D* histSumProtonMomentaXBefore = new TH1D(Form("histSumProtonMomentaXBefore%d", i), Form("histSumProtonMomentaXBefore _%d", i),  20, -1,1);
        HistSumProtonMomentaXBefore.push_back(histSumProtonMomentaXBefore);

        TH1D* histSumProtonMomentaYBefore = new TH1D(Form("histSumProtonMomentaYBefore%d", i), Form("histSumProtonMomentaYBefore _%d", i),40, -2, 2);
        HistSumProtonMomentaYBefore.push_back(histSumProtonMomentaYBefore);

        TH1D* histZDiffBefore = new TH1D(Form("histZDiffBefore%d", i), Form("histZDiffBefore _%d", i),12,-30,30);
        HistZDiffBefore.push_back(histZDiffBefore);

        TH1D* histZMeanBefore = new TH1D(Form("histZMeanBefore%d", i), Form("histZMeanBefore _%d", i),20,-100,100);
        HistZMeanBefore.push_back(histZMeanBefore);

        TH1D* histCosThetaStarBefore = new TH1D(Form("histCosThetaStarBefore%d", i), Form("histCosThetaStarBefore _%d", i), 20, -1.0, 1.0);
        HistCosThetaStarBefore.push_back(histCosThetaStarBefore);
    }
    
    vector <TH1D*> HistInvMassPiPiN1;
    vector <TH1D *> HistPtMissN1;
    vector <TH1D*> HistNTOFClusterN1;
    vector <TH1D*> HistDCABeamlineN1;
    vector <TH1D*> HistDCADaughtersN1;
    vector <TH1D*> HistCosN1;   
    vector <TH1D*> HistDecayN1;            
    vector <TH1D*> HistCorrEN1;
    vector <TH1D*> HistCorrWN1;
    vector <TH1D*> HistCorrKsiN1;
    vector <TH1D*> HistCorrEtaN1;
    vector <TH1D*> HistKsiEN1;
    vector <TH1D*> HistKsiWN1;
    vector <TH1D*> HistSumProtonMomentaXN1;
    vector <TH1D*> HistSumProtonMomentaYN1;
    vector <TH1D*> HistZDiffN1;
    vector <TH1D*> HistZMeanN1;
    vector <TH1D*> HistPhiN1;
    vector <TH1D*> HistCosThetaStarN1;

    for (int i = 2; i < 6; i++)
    {       
        TH1D* histInvMassPiPiN1 = new TH1D(Form("histInvMassPiPiN1%d", i), Form("histInvMassPiPiN1 _%d", i), 30 ,0.44, 0.56);
        HistInvMassPiPiN1.push_back(histInvMassPiPiN1);
   
       TH1D* histPtMissN1 = new TH1D(Form("histPtMissN1%d", i), Form("histPtMissN1 _%d", i),10, 0,0.5);
       HistPtMissN1.push_back(histPtMissN1);

        TH1D* histNTOFClusterN1 = new TH1D(Form("histNTOFClusterN1%d", i), Form("histNTOFClusterN1 _%d", i), 16, -0.5, 15.5);
        HistNTOFClusterN1.push_back(histNTOFClusterN1);

        TH1D* histDCABeamlineN1 = new TH1D(Form("histDCABeamlineN1%d", i), Form("histDCABeamlineN1 _%d", i), 10, 0, 5.0);
        HistDCABeamlineN1.push_back(histDCABeamlineN1);

        TH1D* histDCADaughtersN1 = new TH1D(Form("histDCADaughtersN1%d", i), Form("histDCADaughtersN1 _%d", i), 10, 0, 5.0);
        HistDCADaughtersN1.push_back(histDCADaughtersN1);

        TH1D* histCosN1 = new TH1D(Form("histCosN1%d", i), Form("histCosN1 _%d", i),20,0.9,1);
        HistCosN1.push_back(histCosN1);

        TH1D* histDecayN1 = new TH1D(Form("histDecayN1%d", i), Form("histDecayN1 _%d", i), 18,0, 9);
        HistDecayN1.push_back(histDecayN1);

        TH1D* histCorrEN1 = new TH1D(Form("histCorrEN1%d", i), Form("histCorrEN1 _%d", i),30, -0.015, 0.015);
        HistCorrEN1.push_back(histCorrEN1);

        TH1D* histCorrWN1 = new TH1D(Form("histCorrWN1%d", i), Form("histCorrWN1 _%d", i),30, -0.015, 0.015);
        HistCorrWN1.push_back(histCorrWN1);

        TH1D* histCorrKsiN1 = new TH1D(Form("histCorrKsiN1%d", i), Form("histCorrKsiN1 _%d", i), 10, -0.005,0.005);
        HistCorrKsiN1.push_back(histCorrKsiN1);

        TH1D* histCorrEtaN1 = new TH1D(Form("histCorrEtaN1%d", i), Form("histCorrEtaN1 _%d", i),30, -1.5,1.5);
        HistCorrEtaN1.push_back(histCorrEtaN1);

        TH1D* histKsiEN1 = new TH1D(Form("histKsiEN1%d", i), Form("histKsiEN1 _%d", i), 30, -0.015, 0.015);
        HistKsiEN1.push_back(histKsiEN1);

        TH1D* histKsiWN1 = new TH1D(Form("histKsiWN1%d", i), Form("histKsiWN1 _%d", i),30, -0.015, 0.015);
        HistKsiWN1.push_back(histKsiWN1);

        TH1D* histSumProtonMomentaXN1 = new TH1D(Form("histSumProtonMomentaXN1%d", i), Form("histSumProtonMomentaXN1 _%d", i),  20, -1, 1);
        HistSumProtonMomentaXN1.push_back(histSumProtonMomentaXN1);

        TH1D* histSumProtonMomentaYN1 = new TH1D(Form("histSumProtonMomentaYN1%d", i), Form("histSumProtonMomentaYN1 _%d", i),40, -2, 2);
        HistSumProtonMomentaYN1.push_back(histSumProtonMomentaYN1);

        TH1D* histZDiffN1 = new TH1D(Form("histZDiffN1%d", i), Form("histZDiffN1 _%d", i),12 ,-30,30);
        HistZDiffN1.push_back(histZDiffN1);

        TH1D* histZMeanN1 = new TH1D(Form("histZMeanN1%d", i), Form("histZMeanN1 _%d", i),20,-100,100);
        HistZMeanN1.push_back(histZMeanN1);

        TH1D* histCosThetaStarN1 = new TH1D(Form("histCosThetaStarN1%d", i), Form("histCosThetaStarN1 _%d", i), 20, -1.0, 1.0);
        HistCosThetaStarN1.push_back(histCosThetaStarN1);
    }
    
    TH1D* HistMassK0K0 = new TH1D("HistMassK0K0", " ; m_{K^{0}K^{0}} [GeV] ; events", 21, 0.9, 3);
    TH1D* HistEtaK0K0 = new TH1D("HistEtaK0K0", " ; #eta_{K^{0}K^{0}} [GeV] ; events", 50, -5,5);
    TH1D* HistPtK0K0 = new TH1D("HistPtK0K0", " ; #p_{T,K^{0}K^{0}} [GeV] ; events",50, 0, 2);

    TH1D* HistPxK0K0 = new TH1D("HistPxK0K0", " ; #p_{T,K^{0}K^{0}} [GeV] ; events",50, -1, 1);
    TH1D* HistPyK0K0 = new TH1D("HistPyK0K0", " ; #p_{T,K^{0}K^{0}} [GeV] ; events",50, -2, 2);

    TH1D* HistKsiEK0K0 = new TH1D("HistKsiEK0K0", " ; #xi^{E} ; events",50, -0.005, 0.01);
    TH1D* HistKsiWK0K0 = new TH1D("HistKsiWK0K0", " ; #xi^{W} ; events",50, -0.005, 0.01);
    
    TH1D* HistPhiK0K0 = new TH1D("HistPhiK0K0", " ; #xi^{W} ; events",31, -M_PI, M_PI);
    
    TH2D* Hist2DEtaPhiK0K0 = new TH2D("Hist2DEtaPhiK0K0", " ; #Eta ; #Phi", 50, -5,5, 31, -M_PI, M_PI);
    
    TH2D* Hist2DKsi = new TH2D("Hist2DKsi", " ; #Eta ; #Phi", 20,0.001,0.006, 20,0.001,0.006);

    TH2D* Hist2DEta = new TH2D("Hist2DEta", " ; #Eta ; #Phi", 20, -0.8, 0.8, 20,-0.8,0.8);
    TH2D* Hist2DKsiE = new TH2D("Hist2DKsiE", " ; #Eta ; #Phi",20, -0.004,0.01,20, -0.004 ,0.01);
    TH2D* Hist2DKsiW = new TH2D("Hist2DKsiW", " ; #Eta ; #Phi", 20, -0.004,0.01,20, -0.004 ,0.01);



    TH1D* HistMassK0 = new TH1D("HistMassK0", " ;m [GeV] ; N_{K}", 30 ,0.44, 0.56);
    TH1D* HistPtK0 = new TH1D("HistPtK0", " ;m [GeV] ; N_{K}", 50, 0, 2);
    TH1D* HistEtaK0 = new TH1D("HistEtaK0", " ;m [GeV] ; N_{K}",  50, -5,5);
    TH1D* HistPhiK0 = new TH1D("HistPhiK0", " ;m [GeV] ; N_{K}",31, -M_PI, M_PI);


  
    TH1D* HistMassK0K0Opp = new TH1D("HistMassK0K0Opp", " ; m_{K^{0}K^{0}} [GeV] ; events", 21, 0.9, 3);
      TH1D* HistMassK0K0Same = new TH1D("HistMassK0K0Same", " ; m_{K^{0}K^{0}} [GeV] ; events", 21, 0.9, 3);
  




    // N-1

    vector <TH1D*> HistPtPionWithTofN1;
    vector <TH1D*> HistPtPionWithoutTofN1;
    vector <TH1D*> HistEtaPionWithTofN1;
    vector <TH1D*> HistEtaPionWithoutTofN1;
    vector <TH1D*> HistNfitPionWithTofN1;
    vector <TH1D*> HistNfitPionWithoutTofN1;




    for (int i = 2; i < 6; ++i) 
    {
        TH1D* histPtWithN1 = new TH1D(Form("histPtWithN1%d", i), Form("histPtWithN1%d", i), 50, 0, 5);
        TH1D* histPtWithoutN1 = new TH1D(Form("histPtWithoutN1%d", i), Form("histPtWithoutN1%d", i), 50, 0, 5);
        TH1D* histEtaWithN1 = new TH1D(Form("histEtaWithN1%d", i), Form("histEtaWithN1%d", i), 101, -5, 5);
        TH1D* histEtaWithoutN1 = new TH1D(Form("histEtaWithoutN1%d", i), Form("histEtaWithoutN1%d", i), 101, -5, 5);
        TH1D* histNfitPionWithTofN1 = new TH1D(Form("histNfitPionWithTofN1%d", i), Form("histNfitPionWithTofN1%d", i), 101, -0.5, 100.5);
        TH1D* histNfitPionWithoutTofN1 = new TH1D(Form("histNfitPionWithoutTofN1%d", i), Form("histNfitPionWithoutTofN1%d", i), 101, -0.5, 100.5);
        HistPtPionWithTofN1.push_back(histPtWithN1);
        HistPtPionWithoutTofN1.push_back(histPtWithoutN1);
        HistEtaPionWithTofN1.push_back(histEtaWithN1);
        HistEtaPionWithoutTofN1.push_back(histEtaWithoutN1);
        HistNfitPionWithTofN1.push_back(histNfitPionWithTofN1);
        HistNfitPionWithoutTofN1.push_back(histNfitPionWithoutTofN1);
    }
       vector <vector <TH1D*>> HistPtMissCF;
        vector <vector <TH1D*>> HistInvMassPiPiCF;
        vector <vector <TH1D*>> HistInvMassKKCF;
   
        for (int i = 2; i < 6; i++)
        {
            vector <TH1D*> vHistPtMissCF;
            vector <TH1D*> vHistInvMassPiPiCF;
            vector <TH1D*> vHistInvMassKKCF;

            for (int j = 0; j < 9; j++)
            {
                TH1D* histPtMissCF = new TH1D(Form("histPtMissCF%d%d", i,j), Form("histPtMissCF _%d%d", i,j),20, 0,1);
                vHistPtMissCF.push_back(histPtMissCF);

                TH1D* histInvMassPiPiCF = new TH1D(Form("histInvMassPiPiCF%d%d", i,j), Form("histInvMassPiPiCF _%d%d", i,j),  30 ,0.44, 0.56);
                vHistInvMassPiPiCF.push_back(histInvMassPiPiCF);        

                TH1D* histInvMassKKCF = new TH1D(Form("histInvMassKKCF%d%d", i,j), Form("histInvMassKKCF _%d%d", i,j), 42, 0.9, 3);
                vHistInvMassKKCF.push_back(histInvMassKKCF);
            }
            HistPtMissCF.push_back(vHistPtMissCF);
            HistInvMassPiPiCF.push_back(vHistInvMassPiPiCF);
            HistInvMassKKCF.push_back(vHistInvMassKKCF);
        }




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

    bool isMC = 0;
    static StUPCEvent *upcEvt = 0x0;
    static StRPEvent  *correctedRpEvent = 0x0;


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

    vector <double> xyVec;
    vector <double> xVec;
    vector <double> yVec;
    vector <double> zVec;

    double beamPositionX, beamPositionY;
    double primVertexPosX, primVertexPosY, primVertexPosZ, primVertexPosErrX, primVertexPosErrY, primVertexPosErrZ;
    double distVertexBeamX, distVertexBeamY, distR, significance;

    const vector <int> triggerID = {570209, 570219, 570229, 570701, 570702, 570703, 570704, 570705, 570709, 570711, 570712, 570719, 590701, 590703, 590705, 590708, 590709};
    const vector <int> triggerCEP = {570701, 570705, 570711};

    for (Long64_t i = 0; i < chain->GetEntries(); ++i) 
    {
        chain->GetEntry(i);
        if ( i%100000 == 0 )  
        {
           cout << i << "/" <<  chain->GetEntries() << endl;
        }

        TLorentzVector protonE, protonW;
        double protonSignsY = 1.0;
        FindProtons(isMC, correctedRpEvent, upcEvt, protonE, protonW);  
        protonSignsY = protonE.Y() * protonW.Y();
        
        // SEPARATE TRACKS
		vector <StUPCTrack const*> tracksWithTofHit;
		vector <StUPCTrack const*> tracksWithoutTofHit;
        SeparateTracks(upcEvt, tracksWithTofHit, tracksWithoutTofHit, isMC, HistNumWithTofTrakcs, HistNumWithoutTofTrakcs);
 
        // VALID NUMBER OF TRACKS
        vector <int> vNumberOfTofMatched = {2,3,4,5};
        bool isValidNumberOfTracks =  ValidNumberOfTofTracks(vNumberOfTofMatched, tracksWithTofHit,tracksWithoutTofHit, HistPtPionWithTof,HistPtPionWithoutTof,  HistEtaPionWithTof,HistEtaPionWithoutTof, HistNfitPionWithTof, HistNfitPionWithoutTof  );
        histEvtMult->Fill(tracksWithTofHit.size());

        if (isValidNumberOfTracks == false) {continue;}

        // SELECTION: are tof tracks good
        bool areTofTracksGood = AreTofTracksGood(tracksWithTofHit);
        if (areTofTracksGood == false) {continue;}

       
		// SELECTION: are complementary tracks 
        vector <StUPCTrack const*> goodTracksWithoutTofHit;    
        bool areComplementaryTracks= FindTracks(tracksWithTofHit, tracksWithoutTofHit, goodTracksWithoutTofHit);
        if (areComplementaryTracks == false) {continue;}
        histEvtMult2->Fill(tracksWithTofHit.size());
        double beamPar[4] = {};
        GetBeamPar(upcEvt,beamPar,isMC);

        double kaonMassWindowNarrowLow = 0.48;
        double kaonMassWindowNarrowHigh = 0.52;

        vector <StUPCTrack const *> vPosNegPionLeadingKaon;
        vector <StUPCTrack const *> vPosNegPionSubLeadingKaon;

        bool isFound = FindPions(upcEvt,tracksWithTofHit, tracksWithoutTofHit, goodTracksWithoutTofHit, vPosNegPionLeadingKaon ,vPosNegPionSubLeadingKaon,beamPar);
        if (isFound == false) continue;
        int iH = tracksWithTofHit.size()-2;

        TVector3 const tryVec(0,0,0);
        StUPCV0 leadingKaon(vPosNegPionLeadingKaon[0], vPosNegPionLeadingKaon[1], ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PION, 1, 1, tryVec, beamPar,upcEvt->getMagneticField(), true);
        StUPCV0 subLeadingKaon(vPosNegPionSubLeadingKaon[0], vPosNegPionSubLeadingKaon[1], ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PION, 1, 1, tryVec,beamPar, upcEvt->getMagneticField(), true);

        // SC3
        vector <bool *> vCuts;
        bool areKaonsInNarrowMassWindow = 0;
        bool isNumberOfTofClusterSmall = 0;
        bool isPtMissingSmall = 0;

        double leadingKaonMass = leadingKaon.m();
        double subleadingKaonMass = subLeadingKaon.m();


        double pTmiss;
        int totalCluster = 0;
          
        if ( (leadingKaonMass > kaonMassWindowNarrowLow and  leadingKaonMass < kaonMassWindowNarrowHigh) and (subleadingKaonMass > kaonMassWindowNarrowLow and  subleadingKaonMass < kaonMassWindowNarrowHigh) )
        {
              areKaonsInNarrowMassWindow = 1;      
        }

        isPtMissingSmall = CheckPtMiss(leadingKaon, subLeadingKaon,protonE, protonW, pTmiss);
        isNumberOfTofClusterSmall = CheckNumberOfClusters(upcEvt, tracksWithTofHit, totalCluster);


        vCuts.push_back(&areKaonsInNarrowMassWindow);
        vCuts.push_back(&isPtMissingSmall);
        vCuts.push_back(&isNumberOfTofClusterSmall);
        

        // SC4
        bool areBothDcaDaughtersSmall = false;
        bool areBothDcaBeamlineSmall = false;
        bool areBothPointingAngleOrDecayLength = false;

        double dcaDaughtersLeadingKaon = leadingKaon.dcaDaughters();
        double dcaDaughtersSubLeadingKaon = subLeadingKaon.dcaDaughters();

        double DCABeamlineLeadingKaon = leadingKaon.DCABeamLine();
        double DCABeamlineSubLeadingKaon = subLeadingKaon.DCABeamLine();

        double decayLengthLeadingKaon = leadingKaon.decayLengthHypo();
        double decayLengthSubLeadingKaon = subLeadingKaon.decayLengthHypo();

        double cosLeadingKaon = leadingKaon.pointingAngleHypo();
        double cosSubLeadingKaon = subLeadingKaon.pointingAngleHypo();


        double dcaBeamlineLEAD;
        double dcaDaughtersLEAD;

        double dcaBeamlineSUB;
        double dcaDaughtersSUB;
            
        if (tracksWithTofHit.size() == 4 or tracksWithTofHit.size() == 5)
        {
            dcaDaughtersLEAD = ExclusiveK0K0::DCA_BEAMLINE;
            dcaBeamlineLEAD =  ExclusiveK0K0::DCA_BEAMLINE;
            dcaDaughtersSUB = ExclusiveK0K0::DCA_BEAMLINE;
            dcaBeamlineSUB =  ExclusiveK0K0::DCA_BEAMLINE;
        }

        else if (tracksWithTofHit.size() == 3)
        {
            dcaDaughtersLEAD = ExclusiveK0K0::DCA_BEAMLINE;
            dcaBeamlineLEAD =  ExclusiveK0K0::DCA_BEAMLINE;
            dcaDaughtersSUB = ExclusiveK0K0::DCA_BEAMLINE_23;
            dcaDaughtersSUB =  ExclusiveK0K0::DCA_BEAMLINE_23;
        }

        else
        {
            dcaDaughtersLEAD =ExclusiveK0K0::DCA_BEAMLINE_23;
            dcaBeamlineLEAD =  ExclusiveK0K0::DCA_BEAMLINE_23;
            dcaDaughtersSUB = ExclusiveK0K0::DCA_BEAMLINE_23;
            dcaDaughtersSUB =  ExclusiveK0K0::DCA_BEAMLINE_23;
        }



        double cosPtAngle = ExclusiveK0K0::COS_PNT_ANG;
    


        if (dcaDaughtersLeadingKaon <= dcaDaughtersLEAD and dcaDaughtersSubLeadingKaon <= dcaDaughtersSUB)
        {
            areBothDcaDaughtersSmall = true;
        }

        if (DCABeamlineLeadingKaon <= dcaBeamlineLEAD and DCABeamlineSubLeadingKaon <= dcaBeamlineSUB)
        {
            areBothDcaBeamlineSmall = true;
        }

        if (tracksWithTofHit.size() == 5)
        {
            if ((decayLengthLeadingKaon <= ExclusiveK0K0::DECAY_LENGTH or cosLeadingKaon >= cosPtAngle) and (decayLengthSubLeadingKaon <= ExclusiveK0K0::DECAY_LENGTH or cosSubLeadingKaon >= cosPtAngle))
            {
                areBothPointingAngleOrDecayLength = true;
            }
        }

        else if (tracksWithTofHit.size() == 4)
        {
            if ((decayLengthLeadingKaon <= ExclusiveK0K0::DECAY_LENGTH or cosLeadingKaon >= cosPtAngle) and (decayLengthSubLeadingKaon <= ExclusiveK0K0::DECAY_LENGTH or cosSubLeadingKaon >= cosPtAngle))
            {
                areBothPointingAngleOrDecayLength = true;
            }
        }

        else if ((tracksWithTofHit.size() == 3)) 
        {
            if ((decayLengthLeadingKaon <= ExclusiveK0K0::DECAY_LENGTH or cosLeadingKaon >= cosPtAngle) and cosSubLeadingKaon >= 0.95)
            {
                areBothPointingAngleOrDecayLength = true;
            }
        }

        else if (tracksWithTofHit.size() == 2)
        {
            if (cosLeadingKaon >= 0.95 and cosSubLeadingKaon >= 0.95 )
            {
                areBothPointingAngleOrDecayLength = true;
            }
        }

        vCuts.push_back(&areBothDcaDaughtersSmall);
        vCuts.push_back(&areBothDcaBeamlineSmall);
        vCuts.push_back(&areBothPointingAngleOrDecayLength);
        

        // SC5 
        bool areProtonsKsiGood = false;
        bool areProtonsMomentaGood = false;
        bool antiElastic = false;


     
        double ksiE = 0;
        double ksiW = 0;
        double sumProtonMomentumX = protonE.X() + protonW.X();
        double sumProtonMomentumY = protonE.Y() + protonW.Y();
        
        if (isMC==0)
        {
            ksiE = correctedRpEvent->getTrack(0)->xi(ExclusiveK0K0::BEAM_ENERGY);
            ksiW = correctedRpEvent->getTrack(1)->xi(ExclusiveK0K0::BEAM_ENERGY);
        }

        else if (isMC==1)
        {
            ksiE =(255. - protonE.E())/255.;
            ksiW =(255.- protonW.E())/255.;
        }

        if (abs(ksiE) >= 0.007 or abs(ksiW) >= 0.007)
        {
            areProtonsKsiGood = true;
        }

        if (abs(sumProtonMomentumX) >= 0.1 or abs(sumProtonMomentumY) >= 0.1 )
        {
            areProtonsMomentaGood = true;
        }

        if (areProtonsMomentaGood or areProtonsKsiGood)
        {
            antiElastic = true;
        }


        vCuts.push_back(&antiElastic);

        //SC 6 
        bool isKsiCorrSmall = false;  
        bool isEtaCorrSmall = false;  
        bool areSeparateCorr = false; 



        TLorentzVector K01 = leadingKaon.lorentzVector();
        TLorentzVector K02 = subLeadingKaon.lorentzVector();
        TLorentzVector K0K0 = K01 + K02;
        
        double ksiENew = ksiE;
        double ksiWNew = ksiW;
        double massK0K0 = K0K0.M();

        if (ksiE < 0.006)
        {
            ksiENew = 0.003;
        }

        if (ksiW < 0.006)
        {
            ksiWNew = 0.003;
        }
       
        double ksiCorrelation = massK0K0/510.0-sqrt(ksiENew*ksiWNew);
        double etaCorrelation = K0K0.Rapidity()-0.5*log(ksiENew/ksiWNew);
        double corrE = massK0K0/(510.0)*exp(-1*K0K0.Rapidity())-ksiE;
        double corrW = massK0K0/(510.0)*exp(+1*K0K0.Rapidity())-ksiW;




        double corrEW;
        double corrEta;
        double corrKsi;

        if (tracksWithTofHit.size() == 2 or tracksWithTofHit.size() == 3)
        {
            corrEW = 0.004;
            corrEta = 0.9;
            corrKsi = 0.004;
        }

        else
        {
            corrEW = 0.004;
            corrEta = 0.9;
            corrKsi = 0.004;
        }

        if ( abs(corrE) < corrEW and abs(corrW) < corrEW)
        {
            areSeparateCorr = true;
        }

        if (abs(ksiCorrelation ) < corrKsi)
        {
            isKsiCorrSmall = true;
        }

        if (abs(etaCorrelation) < corrEta)
        {
            isEtaCorrSmall = true;
        }

        vCuts.push_back(&isKsiCorrSmall);
        vCuts.push_back(&isEtaCorrSmall);
        vCuts.push_back(&areSeparateCorr);

        // SC7
        TVector3 vertexLeadingKaon = leadingKaon.decayVertex();
        TVector3 vertexSubLeadingKaon = subLeadingKaon.decayVertex();
        double zdiff = vertexLeadingKaon.Z()-vertexSubLeadingKaon.Z();
        double zmean = (vertexLeadingKaon.Z()+vertexSubLeadingKaon.Z())/2.0;

        double cosThetaStarLeading = leadingKaon.cosThetaStar();
        double cosThetaStarSubLeading = subLeadingKaon.cosThetaStar();

    
        bool isZdiffSmall = false;
        if (abs(zdiff) <= 15.0)
        {
            isZdiffSmall = true;
        }
        vCuts.push_back(&isZdiffSmall);

        bool isZmeanSmall = false;
        if (abs(zmean) <= 80.0)
        {
            isZmeanSmall = true;
        }
        vCuts.push_back(&isZmeanSmall);
        // SC8


        bool areCosThetaStarSmall = false;
        if ( abs(cosThetaStarSubLeading) <= 0.8 and abs(cosThetaStarLeading) <= 0.8  ) 
        {
            areCosThetaStarSmall = true;
        }
    
        vCuts.push_back(&areCosThetaStarSmall);

      //  cout << vCuts.size() << endl;

        // N - 1

        //SC4
        bool N1_areKaonsInNarrowMassWindow = true;
        CheckN1(vCuts,  areKaonsInNarrowMassWindow,N1_areKaonsInNarrowMassWindow);
        if (N1_areKaonsInNarrowMassWindow==1) { HistInvMassPiPiN1[iH]->Fill(leadingKaonMass); HistInvMassPiPiN1[iH]->Fill(subleadingKaonMass);}

        bool N1_isPtMissingSmall = true;
        CheckN1(vCuts,  isPtMissingSmall,N1_isPtMissingSmall);
        if (N1_isPtMissingSmall==1) 
        { 
            HistPtMissN1[iH]->Fill(pTmiss);
            if (isPtMissingSmall)
            {
                     
                HistEtaK0K0->Fill(K0K0.Eta()); 
                if (protonSignsY > 0.0)
                {
                     HistMassK0K0Same->Fill(massK0K0); 
                       HistMassK0K0->Fill(massK0K0);  
                        HistPtK0K0->Fill(K0K0.Pt()); 
                }
                else
                {
                     HistMassK0K0Opp->Fill(massK0K0); 
                }
               
                HistPxK0K0->Fill(K0K0.Px()); 
                HistPyK0K0->Fill(K0K0.Py()); 
                HistKsiEK0K0->Fill(ksiE);  
                HistKsiWK0K0->Fill(ksiW);  
                HistPhiK0K0->Fill(K0K0.Phi());            
                Hist2DEtaPhiK0K0->Fill(K0K0.Eta(), K0K0.Phi() );

                Hist2DKsi->Fill(massK0K0/510.0,sqrt(ksiENew*ksiWNew));
                Hist2DEta->Fill( K0K0.Rapidity(), 0.5*log(ksiENew/ksiWNew));
                Hist2DKsiE->Fill(   massK0K0/(510.0)*exp(-1*K0K0.Rapidity()), ksiE);   
                Hist2DKsiW->Fill(massK0K0/(510.0)*exp(+1*K0K0.Rapidity()),ksiW);   


                HistMassK0->Fill(leadingKaon.m()); 
                HistMassK0->Fill(subLeadingKaon.m()); 

                HistPtK0->Fill(leadingKaon.pt()); 
                HistPtK0->Fill(subLeadingKaon.pt());             
        
                HistEtaK0->Fill(leadingKaon.eta()); 
                HistEtaK0->Fill(subLeadingKaon.eta());             
        
                HistPhiK0->Fill(leadingKaon.phi()); 
                HistPhiK0->Fill(subLeadingKaon.phi());             
        
  
  

            }
        }

        bool N1_isNumberOfTofClusterSmall = true;
        CheckN1(vCuts,  isNumberOfTofClusterSmall, N1_isNumberOfTofClusterSmall);
        if (N1_isNumberOfTofClusterSmall==1) { HistNTOFClusterN1[iH]->Fill(totalCluster);}

        //SC5
        bool N1_areBothDcaBeamlineSmall = true;
        CheckN1(vCuts,  areBothDcaBeamlineSmall, N1_areBothDcaBeamlineSmall);
        if (N1_areBothDcaBeamlineSmall==1) { HistDCABeamlineN1[iH]->Fill(DCABeamlineLeadingKaon);   HistDCABeamlineN1[iH]->Fill(DCABeamlineSubLeadingKaon); }

        bool N1_areBothDcaDaughtersSmall = true;
        CheckN1(vCuts,  areBothDcaDaughtersSmall, N1_areBothDcaDaughtersSmall);
        if (N1_areBothDcaDaughtersSmall==1) { HistDCADaughtersN1[iH]->Fill(dcaDaughtersLeadingKaon);   HistDCADaughtersN1[iH]->Fill(dcaDaughtersSubLeadingKaon); }

        bool N1_areBothPointingAngleOrDecayLength = true;
        CheckN1(vCuts, areBothPointingAngleOrDecayLength, N1_areBothPointingAngleOrDecayLength);
        if (N1_areBothPointingAngleOrDecayLength==1)
        { 
            HistCosN1[iH]->Fill(cosLeadingKaon);   HistCosN1[iH]->Fill(cosSubLeadingKaon);
            HistDecayN1[iH]->Fill(decayLengthLeadingKaon);   HistDecayN1[iH]->Fill(decayLengthSubLeadingKaon);
        }

        //SC6
        bool N1_antiElastic = true;
        CheckN1(vCuts, antiElastic, N1_antiElastic);
        if (N1_antiElastic==1)
        { 
            HistSumProtonMomentaXN1[iH]->Fill(sumProtonMomentumX);  
            HistSumProtonMomentaYN1[iH]->Fill(sumProtonMomentumY);  
            HistKsiEN1[iH]->Fill(ksiE);  
            HistKsiWN1[iH]->Fill(ksiW);             
        }
        
        //SC71
        bool N1_isKsiCorrSmall = true;
        CheckN1(vCuts, isKsiCorrSmall, N1_isKsiCorrSmall);
        if (N1_isKsiCorrSmall==1)
        { 
            HistCorrKsiN1[iH]->Fill(ksiCorrelation);             
        }
        
        //SC72
        bool N1_isEtaCorrSmall = true;
        CheckN1(vCuts, isEtaCorrSmall, N1_isEtaCorrSmall);
        if (N1_isEtaCorrSmall==1)
        { 
            HistCorrEtaN1[iH]->Fill(etaCorrelation);             
        }
        
        //SC73
        bool N1_areSeparateCorr = true;
        CheckN1(vCuts, areSeparateCorr, N1_areSeparateCorr);
        if (N1_areSeparateCorr==1)
        { 
            HistCorrEN1[iH]->Fill(corrE); 
            HistCorrWN1[iH]->Fill(corrW);                 
        }

        //SC81
        bool N1_isZdiffSmall = true;
        CheckN1(vCuts, isZdiffSmall, N1_isZdiffSmall);
        if (N1_isZdiffSmall==1)
        { 
            HistZDiffN1[iH]->Fill(zdiff);      
        }

        //SC82
        bool N1_isZmeanSmall = true;
        CheckN1(vCuts, isZmeanSmall, N1_isZmeanSmall);
        if (N1_isZmeanSmall==1)
        { 
           HistZMeanN1[iH]->Fill(zmean);  
        }

        //SC9
        bool N1_areCosThetaStarSmall = true;
        CheckN1(vCuts, areCosThetaStarSmall, N1_areCosThetaStarSmall);
        if (N1_areCosThetaStarSmall==1)
        { 
            HistCosThetaStarN1[iH]->Fill(cosThetaStarLeading);
            HistCosThetaStarN1[iH]->Fill(cosThetaStarSubLeading);
        }


        //Filling histograms BEFORE
        HistPtMissBefore[iH]->Fill(pTmiss); 
        HistNTOFClusterBefore[iH]->Fill(totalCluster);
        HistDCABeamlineBefore[iH]->Fill(DCABeamlineLeadingKaon);
        HistDCABeamlineBefore[iH]->Fill(DCABeamlineSubLeadingKaon);
        HistDCADaughtersBefore[iH]->Fill(dcaDaughtersLeadingKaon);
        HistDCADaughtersBefore[iH]->Fill(dcaDaughtersSubLeadingKaon);
        HistCosBefore[iH]->Fill(cosLeadingKaon);  
        HistCosBefore[iH]->Fill(cosSubLeadingKaon);   
        HistDecayBefore[iH]->Fill(decayLengthLeadingKaon);  
        HistDecayBefore[iH]->Fill(decayLengthSubLeadingKaon);       
        HistCorrEBefore[iH]->Fill(corrE);
        HistCorrWBefore[iH]->Fill(corrW);
        HistCorrKsiBefore[iH]->Fill(ksiCorrelation);
        HistCorrEtaBefore[iH]->Fill(etaCorrelation);
        HistKsiEBefore[iH]->Fill(ksiE);
        HistKsiWBefore[iH]->Fill(ksiW);
        HistSumProtonMomentaXBefore[iH]->Fill(sumProtonMomentumX);
        HistSumProtonMomentaYBefore[iH]->Fill(sumProtonMomentumY);
        HistZMeanBefore[iH]->Fill(zmean);
        HistZDiffBefore[iH]->Fill(zdiff);
        HistInvMassPiPi2D[iH]->Fill(leadingKaonMass, subleadingKaonMass);
        HistCosThetaStarBefore[iH]->Fill(cosThetaStarLeading);
        HistCosThetaStarBefore[iH]->Fill(cosThetaStarSubLeading);


        if (areKaonsInNarrowMassWindow)
        {  
            HistPtMissCF[iH][0]->Fill(pTmiss);  HistInvMassPiPiCF[iH][0]->Fill(leadingKaonMass); HistInvMassPiPiCF[iH][0]->Fill(subleadingKaonMass);  HistInvMassKKCF[iH][0]->Fill(massK0K0);


            if (isNumberOfTofClusterSmall)
            {
                
                HistPtMissCF[iH][1]->Fill(pTmiss);  HistInvMassPiPiCF[iH][1]->Fill(leadingKaonMass); HistInvMassPiPiCF[iH][1]->Fill(subleadingKaonMass);  HistInvMassKKCF[iH][1]->Fill(massK0K0);
                
                if (areBothDcaBeamlineSmall)
                {
                    HistPtMissCF[iH][2]->Fill(pTmiss);  HistInvMassPiPiCF[iH][2]->Fill(leadingKaonMass); HistInvMassPiPiCF[iH][2]->Fill(subleadingKaonMass);  HistInvMassKKCF[iH][2]->Fill(massK0K0);
                
                    if (areBothDcaDaughtersSmall)
                    {
                        HistPtMissCF[iH][3]->Fill(pTmiss);  HistInvMassPiPiCF[iH][3]->Fill(leadingKaonMass); HistInvMassPiPiCF[iH][3]->Fill(subleadingKaonMass);  HistInvMassKKCF[iH][3]->Fill(massK0K0);
                
                        if (areBothPointingAngleOrDecayLength)
                        {

                            HistPtMissCF[iH][4]->Fill(pTmiss);  HistInvMassPiPiCF[iH][4]->Fill(leadingKaonMass); HistInvMassPiPiCF[iH][4]->Fill(subleadingKaonMass);  HistInvMassKKCF[iH][4]->Fill(massK0K0);
                
                            if (antiElastic)
                            {
                                HistPtMissCF[iH][5]->Fill(pTmiss);  HistInvMassPiPiCF[iH][5]->Fill(leadingKaonMass); HistInvMassPiPiCF[iH][5]->Fill(subleadingKaonMass);  HistInvMassKKCF[iH][5]->Fill(massK0K0);
                
                                if (isKsiCorrSmall)
                                {
                                    HistPtMissCF[iH][6]->Fill(pTmiss);  HistInvMassPiPiCF[iH][6]->Fill(leadingKaonMass); HistInvMassPiPiCF[iH][6]->Fill(subleadingKaonMass);  HistInvMassKKCF[iH][6]->Fill(massK0K0);
                
                                    if (isEtaCorrSmall)
                                    {
                                        HistPtMissCF[iH][7]->Fill(pTmiss);  HistInvMassPiPiCF[iH][7]->Fill(leadingKaonMass); HistInvMassPiPiCF[iH][7]->Fill(subleadingKaonMass);  HistInvMassKKCF[iH][7]->Fill(massK0K0);
                
                                        if (areSeparateCorr)
                                        {               
                                            HistPtMissCF[iH][8]->Fill(pTmiss);  HistInvMassPiPiCF[iH][8]->Fill(leadingKaonMass); HistInvMassPiPiCF[iH][8]->Fill(subleadingKaonMass);  HistInvMassKKCF[iH][8]->Fill(massK0K0);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

    // end lop events
    }

    TFile *outfile = TFile::Open(argv[2], "recreate"); 

    HistNumWithTofTrakcs->Write();
    HistNumWithoutTofTrakcs->Write();

    for (int i = 0; i < 4; ++i) 
    {
        histEvtMult->Write();
        histEvtMult2->Write();
        // tracks distributions
        HistInvMassPiPi2D[i]->Write();
        HistPtPionWithTof[i]->Write();
        HistEtaPionWithTof[i]->Write();
        HistPtPionWithoutTof[i]->Write();
        HistEtaPionWithoutTof[i]->Write();
        HistNfitPionWithTof[i]->Write();
        HistNfitPionWithoutTof[i]->Write();
        
        // N1

        HistPtMissN1[i]->Write();
        HistNTOFClusterN1[i]->Write();
        HistDCABeamlineN1[i]->Write();
        HistDCABeamlineN1[i]->Write();
        HistDCADaughtersN1[i]->Write();
        HistDCADaughtersN1[i]->Write();
        HistCosN1[i]->Write();
        HistCosN1[i]->Write();  
        HistDecayN1[i]->Write();
        HistDecayN1[i]->Write();
        HistSumProtonMomentaXN1[i]->Write();
        HistSumProtonMomentaYN1[i]->Write();
        HistCorrEN1[i]->Write();
        HistCorrWN1[i]->Write();
        HistCorrKsiN1[i]->Write();
        HistCorrEtaN1[i]->Write();
        HistKsiEN1[i]->Write();
        HistKsiWN1[i]->Write();

        HistZMeanN1[i]->Write();
        HistZDiffN1[i]->Write();
        HistCosThetaStarN1[i]->Write();
        //before 
        HistPtMissBefore[i]->Write();
        HistNTOFClusterBefore[i]->Write();
        HistDCABeamlineBefore[i]->Write();
        HistDCABeamlineBefore[i]->Write();
        HistDCADaughtersBefore[i]->Write();
        HistDCADaughtersBefore[i]->Write();
        HistCosBefore[i]->Write();
        HistCosBefore[i]->Write();  
        HistDecayBefore[i]->Write();
        HistDecayBefore[i]->Write();
        HistCorrEBefore[i]->Write();
        HistCorrWBefore[i]->Write();
        HistCorrKsiBefore[i]->Write();
        HistCorrEtaBefore[i]->Write();
        HistKsiEBefore[i]->Write();
        HistKsiWBefore[i]->Write();
        HistSumProtonMomentaXBefore[i]->Write();
        HistSumProtonMomentaYBefore[i]->Write();
        HistZDiffBefore[i]->Write();
        HistZMeanBefore[i]->Write();
        HistCosThetaStarBefore[i]->Write();

HistMassK0K0Same->Write();
HistMassK0K0Opp->Write();
        //cutflow
        for (int j = 0; j < 9; j++)
        {
            HistPtMissCF[i][j]->Write();
            HistInvMassPiPiCF[i][j]->Write();
            HistInvMassKKCF[i][j]->Write();
        }
        //N-1

        HistPtPionWithTofN1[i]->Write();
        HistPtPionWithoutTofN1[i]->Write();
        HistEtaPionWithTofN1[i]->Write();
        HistEtaPionWithoutTofN1[i]->Write();
        HistNfitPionWithTofN1[i]->Write();
        HistNfitPionWithoutTofN1[i]->Write();

    }
    HistMassK0K0->Write();      
    HistEtaK0K0->Write();
    HistPtK0K0->Write();
    HistPxK0K0->Write();
    HistPyK0K0->Write(); 
    HistKsiEK0K0->Write();
    HistKsiWK0K0->Write();  
    HistPhiK0K0->Write();    
    Hist2DEtaPhiK0K0->Write();
    Hist2DKsi->Write();
    Hist2DEta->Write();
    Hist2DKsiE->Write();
    Hist2DKsiW->Write();

    HistMassK0->Write();
    HistPtK0->Write();
    HistEtaK0->Write();
    HistPhiK0->Write();



    outfile->Close();

    return 0;
}


