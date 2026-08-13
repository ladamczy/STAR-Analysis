#include "ExclusiveCode.h"
#include "SystematicCuts.h"
using namespace std;

// Helper: append systematic suffix to histogram base name
string MakeHistName(const char* base, const string& suffix) {
    return string(base) + "_" + suffix;
}

int main(int argc, char** argv)  
{
    // argv[1] = input file list
    // argv[2] = output file
    // argv[3] = (optional) 0 = nominal only (default), 1 = all systematics

    bool runAllSystematics = false;
    if (argc > 3) {
        runAllSystematics = (atoi(argv[3]) == 1);
    }

    vector<CutConfig> configs;
    if (runAllSystematics) {
        configs = CreateSystematicConfigs();
        cout << "ℹ️ Running " << configs.size() << " systematic variations" << endl;
    } else {
        configs.push_back(CutConfig("nominal"));
        cout << "ℹ️ Running NOMINAL only" << endl;
    }

    // Open chain once — shared across all systematics
    ifstream inputFilePathList(argv[1]);
    if (!inputFilePathList) {
        cerr << "Failed to open input file." << endl;
        return 1;
    }
    TChain *chain = new TChain("mUPCTree");
    string inputFileName;
    while (getline(inputFilePathList, inputFileName))
        chain->Add(inputFileName.c_str());
    inputFilePathList.close();

    bool isMC = 0;
    static StUPCEvent *upcEvt = 0x0;
    static StRPEvent  *correctedRpEvent = 0x0;
    chain->SetBranchAddress("mUPCEvent", &upcEvt);
    chain->GetEntry(0);
    if (upcEvt->getRunNumber() == 1) isMC = 1;
    if (isMC == 0) chain->SetBranchAddress("correctedRpEvent", &correctedRpEvent);
    //chain->SetBranchAddress("mUPCEvent", &upcEvt);
    //chain->SetBranchAddress("correctedRpEvent", &correctedRpEvent);
    //chain->GetEntry(0);
    //if (upcEvt->getRunNumber() == 1) isMC = 1;
    cout << "Running on " << (isMC ? "MC" : "Real Data") << endl;
    cout << "Total entries: " << chain->GetEntries() << endl;
    int nEntries = chain->GetEntries();
    //int nEntries = 100000; // for testing
    if (nEntries < chain->GetEntries()){
         cout << " ⚠️ Running on low statistics, testing for :" << nEntries << " entries" << endl;
    }
    else{
         cout << " ✅ Running on full statistics : " << nEntries << " entries" << endl;
    }
   

    // =========================================================
    // N_TOF CONFIGURATION
    // Change this vector to select which TOF track multiplicities to analyse.
    // Options: any combination of {2, 3, 4, 5}
    // e.g. {4}      = 4-TOF only
    //      {3, 4}   = 3 or 4 TOF (current nominal)
    //      {2,3,4,5} = all combinations
    // The analysis loops over all selected multiplicities simultaneously.
    // =========================================================
    const vector<int> N_TOF_SELECTED = {2,3,4,5};   // <-- change here

    TFile *outfile = TFile::Open(argv[2], "recreate");

    // =========================================================
    // LOOP OVER SYSTEMATIC CONFIGURATIONS
    // =========================================================
    for (const auto& config : configs) {

        cout << "\n========================================" << endl;
        cout << "Processing: " << config.name << endl;
        cout << "========================================" << endl;
        config.Print();

        string sf = config.name;  // suffix for histogram names

        // ---- Histogram booking ----

        TH1D* HistNumWithTofTrakcs    = new TH1D(MakeHistName("HistNumWithTofTrakcs",    sf).c_str(), ";N^{tracks};events", 52, 1.5, 53.5);
        TH1D* HistNumWithoutTofTrakcs = new TH1D(MakeHistName("HistNumWithoutTofTrakcs", sf).c_str(), ";N^{tracks};events", 150, 1.5, 151.5);

        vector<TH1D*> HistPtPionWithTof, HistPtPionWithoutTof;
        vector<TH1D*> HistEtaPionWithTof, HistEtaPionWithoutTof;
        vector<TH1D*> HistNfitPionWithTof, HistNfitPionWithoutTof;
        for (int i = 2; i < 6; ++i) {
            HistPtPionWithTof.push_back(    new TH1D(MakeHistName(Form("histPtWith_%d",i),    sf).c_str(), Form("histPtWith_%d",i),    50,0,5));
            HistPtPionWithoutTof.push_back( new TH1D(MakeHistName(Form("histPtWithout_%d",i), sf).c_str(), Form("histPtWithout_%d",i), 50,0,5));
            HistEtaPionWithTof.push_back(   new TH1D(MakeHistName(Form("histEtaWith_%d",i),   sf).c_str(), Form("histEtaWith_%d",i),   101,-5,5));
            HistEtaPionWithoutTof.push_back(new TH1D(MakeHistName(Form("histEtaWithout_%d",i),sf).c_str(), Form("histEtaWithout_%d",i),101,-5,5));
            HistNfitPionWithTof.push_back(  new TH1D(MakeHistName(Form("histNfitWith_%d",i),  sf).c_str(), Form("histNfitWith_%d",i),  101,-0.5,100.5));
            HistNfitPionWithoutTof.push_back(new TH1D(MakeHistName(Form("histNfitWithout_%d",i),sf).c_str(),Form("histNfitWithout_%d",i),101,-0.5,100.5));
        }

        // BEFORE histograms
        vector<TH1D*> HistPtMissBefore, HistNTOFClusterBefore;
        vector<TH1D*> HistDCABeamlineBefore, HistDCADaughtersBefore;
        vector<TH1D*> HistCosBefore, HistDecayBefore;
        vector<TH1D*> HistCorrEBefore, HistCorrWBefore, HistCorrKsiBefore, HistCorrEtaBefore;
        vector<TH1D*> HistKsiEBefore, HistKsiWBefore;
        vector<TH1D*> HistSumProtonMomentaXBefore, HistSumProtonMomentaYBefore;
        vector<TH1D*> HistZDiffBefore, HistZMeanBefore, HistCosThetaStarBefore;
        vector<TH2D*> HistInvMassPiPi2D;

        for (int i = 2; i < 6; i++) {
            HistInvMassPiPi2D.push_back(         new TH2D(MakeHistName(Form("histInvMassPiPi2D%d",i),sf).c_str(),"",30,0.44,0.56,30,0.44,0.56));
            HistPtMissBefore.push_back(          new TH1D(MakeHistName(Form("histPtMissBefore%d",i),sf).c_str(),"",20,0,1));
            HistNTOFClusterBefore.push_back(     new TH1D(MakeHistName(Form("histNTOFClusterBefore%d",i),sf).c_str(),"",40,-0.5,40.5));
            HistDCABeamlineBefore.push_back(     new TH1D(MakeHistName(Form("histDCABeamlineBefore%d",i),sf).c_str(),"",25,0,5.0));
            HistDCADaughtersBefore.push_back(    new TH1D(MakeHistName(Form("histDCADaughtersBefore%d",i),sf).c_str(),"",25,0,5.0));
            HistCosBefore.push_back(             new TH1D(MakeHistName(Form("histCosBefore%d",i),sf).c_str(),"",100,-1,1.0));
            HistDecayBefore.push_back(           new TH1D(MakeHistName(Form("histDecayBefore%d",i),sf).c_str(),"",50,0,10));
            HistCorrEBefore.push_back(           new TH1D(MakeHistName(Form("histCorrEBefore%d",i),sf).c_str(),"",51,-0.05,0.05));
            HistCorrWBefore.push_back(           new TH1D(MakeHistName(Form("histCorrWBefore%d",i),sf).c_str(),"",51,-0.05,0.05));
            HistCorrKsiBefore.push_back(         new TH1D(MakeHistName(Form("histCorrKsiBefore%d",i),sf).c_str(),"",51,-0.05,0.05));
            HistCorrEtaBefore.push_back(         new TH1D(MakeHistName(Form("histCorrEtaBefore%d",i),sf).c_str(),"",51,-3,3));
            HistKsiEBefore.push_back(            new TH1D(MakeHistName(Form("histKsiEBefore%d",i),sf).c_str(),"",60,-0.01,0.05));
            HistKsiWBefore.push_back(            new TH1D(MakeHistName(Form("histKsiWBefore%d",i),sf).c_str(),"",60,-0.01,0.05));
            HistSumProtonMomentaXBefore.push_back(new TH1D(MakeHistName(Form("histSumProtonMomentaXBefore%d",i),sf).c_str(),"",100,-5,5));
            HistSumProtonMomentaYBefore.push_back(new TH1D(MakeHistName(Form("histSumProtonMomentaYBefore%d",i),sf).c_str(),"",100,-5,5));
            HistZDiffBefore.push_back(           new TH1D(MakeHistName(Form("histZDiffBefore%d",i),sf).c_str(),"",101,-200,200));
            HistZMeanBefore.push_back(           new TH1D(MakeHistName(Form("histZMeanBefore%d",i),sf).c_str(),"",51,-200,200));
            HistCosThetaStarBefore.push_back(    new TH1D(MakeHistName(Form("histCosThetaStarBefore%d",i),sf).c_str(),"",50,-1.0,1.0));
        }

        // N-1 histograms
        vector<TH1D*> HistInvMassPiPiN1, HistPtMissN1, HistNTOFClusterN1;
        vector<TH1D*> HistDCABeamlineN1, HistDCADaughtersN1;
        vector<TH1D*> HistCosN1, HistDecayN1;
        vector<TH1D*> HistCorrEN1, HistCorrWN1, HistCorrKsiN1, HistCorrEtaN1;
        vector<TH1D*> HistKsiEN1, HistKsiWN1;
        vector<TH1D*> HistSumProtonMomentaXN1, HistSumProtonMomentaYN1;
        vector<TH1D*> HistZDiffN1, HistZMeanN1, HistCosThetaStarN1;

        for (int i = 2; i < 6; i++) {
            HistInvMassPiPiN1.push_back(         new TH1D(MakeHistName(Form("histInvMassPiPiN1%d",i),sf).c_str(),"",30,0.44,0.56));
            HistPtMissN1.push_back(              new TH1D(MakeHistName(Form("histPtMissN1%d",i),sf).c_str(),"",20,0,1));
            HistNTOFClusterN1.push_back(         new TH1D(MakeHistName(Form("histNTOFClusterN1%d",i),sf).c_str(),"",40,-0.5,40.5));
            HistDCABeamlineN1.push_back(         new TH1D(MakeHistName(Form("histDCABeamlineN1%d",i),sf).c_str(),"",25,0,5.0));
            HistDCADaughtersN1.push_back(        new TH1D(MakeHistName(Form("histDCADaughtersN1%d",i),sf).c_str(),"",25,0,5.0));
            HistCosN1.push_back(                 new TH1D(MakeHistName(Form("histCosN1%d",i),sf).c_str(),"",100,-1,1.0));
            HistDecayN1.push_back(               new TH1D(MakeHistName(Form("histDecayN1%d",i),sf).c_str(),"",50,0,10));
            HistCorrEN1.push_back(               new TH1D(MakeHistName(Form("histCorrEN1%d",i),sf).c_str(),"",51,-0.05,0.05));
            HistCorrWN1.push_back(               new TH1D(MakeHistName(Form("histCorrWN1%d",i),sf).c_str(),"",51,-0.05,0.05));
            HistCorrKsiN1.push_back(             new TH1D(MakeHistName(Form("histCorrKsiN1%d",i),sf).c_str(),"",51,-0.05,0.05));
            HistCorrEtaN1.push_back(             new TH1D(MakeHistName(Form("histCorrEtaN1%d",i),sf).c_str(),"",51,-3,3));
            HistKsiEN1.push_back(                new TH1D(MakeHistName(Form("histKsiEN1%d",i),sf).c_str(),"",60,-0.01,0.05));
            HistKsiWN1.push_back(                new TH1D(MakeHistName(Form("histKsiWN1%d",i),sf).c_str(),"",60,-0.01,0.05));
            HistSumProtonMomentaXN1.push_back(   new TH1D(MakeHistName(Form("histSumProtonMomentaXN1%d",i),sf).c_str(),"",100,-5,5));
            HistSumProtonMomentaYN1.push_back(   new TH1D(MakeHistName(Form("histSumProtonMomentaYN1%d",i),sf).c_str(),"",100,-5,5));
            HistZDiffN1.push_back(               new TH1D(MakeHistName(Form("histZDiffN1%d",i),sf).c_str(),"",101,-200,200));
            HistZMeanN1.push_back(               new TH1D(MakeHistName(Form("histZMeanN1%d",i),sf).c_str(),"",51,-200,200));
            HistCosThetaStarN1.push_back(        new TH1D(MakeHistName(Form("histCosThetaStarN1%d",i),sf).c_str(),"",50,-1.0,1.0));
        }

        TH1D* HistMassK0K0 = new TH1D(MakeHistName("HistMassK0K0",sf).c_str(), ";m_{K^{0}K^{0}} [GeV];events", 21, 0.9, 3);

        // N-1 track quality plots
        vector<TH1D*> HistPtPionWithTofN1, HistPtPionWithoutTofN1;
        vector<TH1D*> HistEtaPionWithTofN1, HistEtaPionWithoutTofN1;
        vector<TH1D*> HistNfitPionWithTofN1, HistNfitPionWithoutTofN1;
        for (int i = 2; i < 6; ++i) {
            HistPtPionWithTofN1.push_back(    new TH1D(MakeHistName(Form("histPtWithN1%d",i),sf).c_str(),"",50,0,5));
            HistPtPionWithoutTofN1.push_back( new TH1D(MakeHistName(Form("histPtWithoutN1%d",i),sf).c_str(),"",50,0,5));
            HistEtaPionWithTofN1.push_back(   new TH1D(MakeHistName(Form("histEtaWithN1%d",i),sf).c_str(),"",101,-5,5));
            HistEtaPionWithoutTofN1.push_back(new TH1D(MakeHistName(Form("histEtaWithoutN1%d",i),sf).c_str(),"",101,-5,5));
            HistNfitPionWithTofN1.push_back(  new TH1D(MakeHistName(Form("histNfitWithN1%d",i),sf).c_str(),"",101,-0.5,100.5));
            HistNfitPionWithoutTofN1.push_back(new TH1D(MakeHistName(Form("histNfitWithoutN1%d",i),sf).c_str(),"",101,-0.5,100.5));
        }

        // Cutflow: 11 steps
        vector<vector<TH1D*>> HistPtMissCF, HistInvMassPiPiCF, HistInvMassKKCF;
        for (int i = 2; i < 6; i++) {
            vector<TH1D*> vPt, vMass, vKK;
            for (int j = 0; j < 11; j++) {
                vPt.push_back(  new TH1D(MakeHistName(Form("histPtMissCF%d%d",i,j),   sf).c_str(),"",20,0,1));
                vMass.push_back(new TH1D(MakeHistName(Form("histInvMassPiPiCF%d%d",i,j),sf).c_str(),"",30,0.44,0.56));
                vKK.push_back(  new TH1D(MakeHistName(Form("histInvMassKKCF%d%d",i,j),sf).c_str(),"",42,0.9,3));
            }
            HistPtMissCF.push_back(vPt);
            HistInvMassPiPiCF.push_back(vMass);
            HistInvMassKKCF.push_back(vKK);
        }

        // =========================================================
        // EVENT COUNTERS
        // NOTE: SC 2 (CEP trigger, 2 RP tracks opposite sides,
        //             >=3/4 planes, fiducial region, >=2 TOF class-A)
        //       is applied in Preselection.cxx — the input files
        //       to this analyzer are ALREADY preselected on SC2.
        //       We count those events as "input after preselection".
        //
        // SC 3.1-3.4 (pT>0.2, |eta|<0.9, Nfit>=20, charge balance)
        //       are enforced inside AreTofTracksGood() and FindTracks().
        //       We count them explicitly below.
        // =========================================================
        long long cInput         = 0;  // all events entering this analyzer (= after SC2)
        long long cVertex        = 0;  // SC3.0: exactly 1 primary vertex, |z| < 80 cm
        long long cValidTracks   = 0;  // SC3: correct N_TOF multiplicity
        long long cTofGood       = 0;  // SC3.1-3.3: pT>0.2, |eta|<0.9, Nfit>=20
        long long cComplementary = 0;  // SC3.4: charge balance + nHitsDedx>=15
        long long cPionsFound    = 0;  // pion pair reconstruction succeeded
        long long cMassWindow    = 0;  // SC4.1: K0S mass window
        long long cPtMiss        = 0;  // SC4.2: pT miss
        long long cTofClusters   = 0;  // SC4.3: TOF clusters <= 2*N_TOF+1
        long long cDcaDaughters  = 0;  // SC5: DCA daughters
        long long cDcaBeamline   = 0;  // SC5: DCA beamline
        long long cPointing      = 0;  // SC5: pointing angle / decay length
        long long cAntiElastic   = 0;  // SC6: anti-elastic
        long long cKsiCorr       = 0;  // SC7.1: ksi correlation
        long long cEtaCorr       = 0;  // SC7.2: eta correlation
        long long cSepCorr       = 0;  // SC7.3: separate corr (corrE, corrW)
        long long cZmean         = 0;  // SC8.1: zmean
        long long cZdiff         = 0;  // SC8.2: zdiff
        long long cCosThetaStar  = 0;  // SC9: cos(theta*) — FINAL EXCLUSIVE SAMPLE

        // =========================================================
        // EVENT LOOP
        // =========================================================
        for (Long64_t i = 0; i < nEntries; ++i)
        {
            chain->GetEntry(i);
            if (i%100000 == 0)
                cout << i << "/" << nEntries << "  [" << config.name << "]" << endl;

            cInput++;

            // SC3.0: exactly 1 primary vertex AND |z_vtx| < 80 cm
            if (isMC == 0)
            {
                if (upcEvt->getNumberOfVertices() != 1) continue;
                if (abs(upcEvt->getVertex(0)->getPosZ()) >= 80.0) continue;
            }
            cVertex++;

            TLorentzVector protonE, protonW;
            FindProtons(isMC, correctedRpEvent, upcEvt, protonE, protonW);

            vector<StUPCTrack const*> tracksWithTofHit;
            vector<StUPCTrack const*> tracksWithoutTofHit;
            SeparateTracks(upcEvt, tracksWithTofHit, tracksWithoutTofHit, isMC,
                           HistNumWithTofTrakcs, HistNumWithoutTofTrakcs);

            vector<int> vNumberOfTofMatched = N_TOF_SELECTED;
            bool isValidNumberOfTracks = ValidNumberOfTofTracks(
                vNumberOfTofMatched, tracksWithTofHit, tracksWithoutTofHit,
                HistPtPionWithTof, HistPtPionWithoutTof,
                HistEtaPionWithTof, HistEtaPionWithoutTof,
                HistNfitPionWithTof, HistNfitPionWithoutTof);
            if (!isValidNumberOfTracks) continue;
            cValidTracks++;

            if (!AreTofTracksGood(tracksWithTofHit)) continue;
            cTofGood++;

            vector<StUPCTrack const*> goodTracksWithoutTofHit;
            if (!FindTracks(tracksWithTofHit, tracksWithoutTofHit, goodTracksWithoutTofHit)) continue;
            cComplementary++;

            double beamPar[4] = {};
            GetBeamPar(upcEvt, beamPar, isMC);

            vector<StUPCTrack const*> vPosNegPionLeadingKaon;
            vector<StUPCTrack const*> vPosNegPionSubLeadingKaon;
            if (!FindPions(upcEvt, tracksWithTofHit, tracksWithoutTofHit, goodTracksWithoutTofHit,
                           vPosNegPionLeadingKaon, vPosNegPionSubLeadingKaon, beamPar)) continue;
            cPionsFound++;

            int iH = tracksWithTofHit.size() - 2;

            TVector3 const tryVec(0,0,0);
            StUPCV0 leadingKaon(vPosNegPionLeadingKaon[0], vPosNegPionLeadingKaon[1],
                                ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PION,
                                1, 1, tryVec, beamPar, upcEvt->getMagneticField(), true);
            StUPCV0 subLeadingKaon(vPosNegPionSubLeadingKaon[0], vPosNegPionSubLeadingKaon[1],
                                   ExclusiveK0K0::MASS_PION, ExclusiveK0K0::MASS_PION,
                                   1, 1, tryVec, beamPar, upcEvt->getMagneticField(), true);

            vector<bool*> vCuts;

            // SC4 ---------------------------------------------------------
            double leadingKaonMass    = leadingKaon.m();
            double subleadingKaonMass = subLeadingKaon.m();
            double pTmiss;
            int    totalCluster = 0;

            bool areKaonsInNarrowMassWindow = false;
            if (leadingKaonMass    > config.massWinLow && leadingKaonMass    < config.massWinHigh &&
                subleadingKaonMass > config.massWinLow && subleadingKaonMass < config.massWinHigh)
                areKaonsInNarrowMassWindow = true;

            bool isPtMissingSmall          = CheckPtMiss(leadingKaon, subLeadingKaon, protonE, protonW, pTmiss);
            // Override with config threshold (for systematics)
            isPtMissingSmall = (pTmiss <= config.ptMissMax);

            // SC4.3: TOF clusters — call helper to compute totalCluster,
            // then apply threshold = 2*N_TOF + 1 (matches thesis exactly).
            // For systematics the config can tighten/loosen the threshold.
            CheckNumberOfClusters(upcEvt, tracksWithTofHit, totalCluster);
            int clusterThreshold = (int)(2 * tracksWithTofHit.size() + 1);
            bool isNumberOfTofClusterSmall = (totalCluster <= clusterThreshold);

            vCuts.push_back(&areKaonsInNarrowMassWindow);
            vCuts.push_back(&isPtMissingSmall);
            vCuts.push_back(&isNumberOfTofClusterSmall);

            if (areKaonsInNarrowMassWindow) cMassWindow++;
            if (areKaonsInNarrowMassWindow && isPtMissingSmall) cPtMiss++;
            if (areKaonsInNarrowMassWindow && isPtMissingSmall && isNumberOfTofClusterSmall) cTofClusters++;

            // SC5 ---------------------------------------------------------
            bool areBothDcaDaughtersSmall        = false;
            bool areBothDcaBeamlineSmall         = false;
            bool areBothPointingAngleOrDecayLength = false;

            double dcaDaughtersLeadingKaon    = leadingKaon.dcaDaughters();
            double dcaDaughtersSubLeadingKaon = subLeadingKaon.dcaDaughters();
            double DCABeamlineLeadingKaon     = leadingKaon.DCABeamLine();
            double DCABeamlineSubLeadingKaon  = subLeadingKaon.DCABeamLine();
            double decayLengthLeadingKaon     = leadingKaon.decayLengthHypo();
            double decayLengthSubLeadingKaon  = subLeadingKaon.decayLengthHypo();
            double cosLeadingKaon             = leadingKaon.pointingAngleHypo();
            double cosSubLeadingKaon          = subLeadingKaon.pointingAngleHypo();

            // Use config values for thresholds
            /*double dcaBeamlineCut  = (tracksWithTofHit.size() == 4) ? config.dcaBeamline4  : config.dcaBeamline23;
            double dcaDaughtersCut = (tracksWithTofHit.size() == 4) ? config.dcaDaughters4 : config.dcaDaughters23;
            double cosPA           = (tracksWithTofHit.size() <= 3) ? config.cosPA23       : config.cosPA;

            if (dcaDaughtersLeadingKaon <= dcaDaughtersCut && dcaDaughtersSubLeadingKaon <= dcaDaughtersCut)
                areBothDcaDaughtersSmall = true;

            if (DCABeamlineLeadingKaon <= dcaBeamlineCut && DCABeamlineSubLeadingKaon <= dcaBeamlineCut)
                areBothDcaBeamlineSmall = true;*/

            // 1. Determine Leading Kaon Cuts
            // If the event has 3, 4, or 5 TOF tracks, the Leading Kaon is guaranteed to have 2 TOF hits.
            // If the event has 2 TOF tracks, the Leading Kaon has 1 TOF hit.
            double dcaBeamLeadingCut = (tracksWithTofHit.size() >= 3) ? config.dcaBeamline4 : config.dcaBeamline23;
            double dcaDauLeadingCut  = (tracksWithTofHit.size() >= 3) ? config.dcaDaughters4 : config.dcaDaughters23;

            // 2. Determine SubLeading Kaon Cuts
            // If the event has 4 or 5 TOF tracks, the SubLeading Kaon has 2 TOF hits.
            // If the event has 2 or 3 TOF tracks, the SubLeading Kaon has 1 TOF hit.
            double dcaBeamSubLeadingCut = (tracksWithTofHit.size() >= 4) ? config.dcaBeamline4 : config.dcaBeamline23;
            double dcaDauSubLeadingCut  = (tracksWithTofHit.size() >= 4) ? config.dcaDaughters4 : config.dcaDaughters23;

            // 3. Apply the pointing angle cut (This is already correct in your code, keeping for context)
            double cosPA = (tracksWithTofHit.size() <= 3) ? config.cosPA23 : config.cosPA;

            // 4. Evaluate the conditions using the independent per-kaon cuts
            if (dcaDaughtersLeadingKaon <= dcaDauLeadingCut && dcaDaughtersSubLeadingKaon <= dcaDauSubLeadingCut)
                areBothDcaDaughtersSmall = true;

            if (DCABeamlineLeadingKaon <= dcaBeamLeadingCut && DCABeamlineSubLeadingKaon <= dcaBeamSubLeadingCut)
                areBothDcaBeamlineSmall = true;

            if (tracksWithTofHit.size() == 5 || tracksWithTofHit.size() == 4)
            {
                if ((decayLengthLeadingKaon    <= config.decayLength || cosLeadingKaon    >= config.cosPA) &&
                    (decayLengthSubLeadingKaon <= config.decayLength || cosSubLeadingKaon >= config.cosPA))
                    areBothPointingAngleOrDecayLength = true;
            }
            else if (tracksWithTofHit.size() == 3)
            {
                if ((decayLengthLeadingKaon <= config.decayLength || cosLeadingKaon >= config.cosPA) &&
                    cosSubLeadingKaon >= config.cosPA23)
                    areBothPointingAngleOrDecayLength = true;
            }
            else if (tracksWithTofHit.size() == 2)
            {
                if (cosLeadingKaon >= config.cosPA23 && cosSubLeadingKaon >= config.cosPA23)
                    areBothPointingAngleOrDecayLength = true;
            }

            vCuts.push_back(&areBothDcaDaughtersSmall);
            vCuts.push_back(&areBothDcaBeamlineSmall);
            vCuts.push_back(&areBothPointingAngleOrDecayLength);

            // running AND of all cuts so far for sequential counter
            bool passedSC4 = areKaonsInNarrowMassWindow && isPtMissingSmall && isNumberOfTofClusterSmall;
            if (passedSC4 && areBothDcaDaughtersSmall)                                       cDcaDaughters++;
            if (passedSC4 && areBothDcaDaughtersSmall && areBothDcaBeamlineSmall)            cDcaBeamline++;
            if (passedSC4 && areBothDcaDaughtersSmall && areBothDcaBeamlineSmall
                          && areBothPointingAngleOrDecayLength)                              cPointing++;

            // SC6 (anti-elastic) -----------------------------------------
            bool areProtonsKsiGood   = false;
            bool areProtonsMomentaGood = false;
            bool antiElastic         = false;
            vCuts.push_back(&antiElastic);

            double ksiE = 0, ksiW = 0;
            double sumProtonMomentumX = protonE.X() + protonW.X();
            double sumProtonMomentumY = protonE.Y() + protonW.Y();

            if (isMC == 0) {
                ksiE = correctedRpEvent->getTrack(0)->xi(ExclusiveK0K0::BEAM_ENERGY);
                ksiW = correctedRpEvent->getTrack(1)->xi(ExclusiveK0K0::BEAM_ENERGY);
            } else {
                ksiE = (255. - protonE.E()) / 255.;
                ksiW = (255. - protonW.E()) / 255.;
            }

            if (abs(ksiE) >= 0.007 || abs(ksiW) >= 0.007)         areProtonsKsiGood    = true;
            if (abs(sumProtonMomentumX) >= 0.1 || abs(sumProtonMomentumY) >= 0.1) areProtonsMomentaGood = true;
            if (areProtonsMomentaGood || areProtonsKsiGood)        antiElastic          = true;

            bool passedSC5 = passedSC4 && areBothDcaDaughtersSmall && areBothDcaBeamlineSmall && areBothPointingAngleOrDecayLength;
            if (passedSC5 && antiElastic) cAntiElastic++;

            // SC7 (kinematic correlations) --------------------------------
            bool isKsiCorrSmall   = false;
            bool isEtaCorrSmall   = false;
            bool areSeparateCorr  = false;

            TLorentzVector K0K0 = leadingKaon.lorentzVector() + subLeadingKaon.lorentzVector();
            double massK0K0 = K0K0.M();

            double ksiENew = (ksiE < 0.006) ? 0.003 : ksiE;
            double ksiWNew = (ksiW < 0.006) ? 0.003 : ksiW;

            double ksiCorrelation = massK0K0/510.0 - sqrt(ksiENew*ksiWNew);
            double etaCorrelation = K0K0.Rapidity() - 0.5*log(ksiENew/ksiWNew);
            double corrE = massK0K0/510.0 * exp(-K0K0.Rapidity()) - ksiE;
            double corrW = massK0K0/510.0 * exp(+K0K0.Rapidity()) - ksiW;

            if (abs(corrE)          < config.separateCorrMax && abs(corrW)          < config.separateCorrMax) areSeparateCorr = true;
            if (abs(ksiCorrelation) < config.ksiCorrMax)                                                       isKsiCorrSmall  = true;
            if (abs(etaCorrelation) < config.etaCorrMax)                                                       isEtaCorrSmall  = true;

            vCuts.push_back(&isKsiCorrSmall);
            vCuts.push_back(&isEtaCorrSmall);
            vCuts.push_back(&areSeparateCorr);

            bool passedSC6 = passedSC5 && antiElastic;
            if (passedSC6 && isKsiCorrSmall)                             cKsiCorr++;
            if (passedSC6 && isKsiCorrSmall && isEtaCorrSmall)           cEtaCorr++;
            if (passedSC6 && isKsiCorrSmall && isEtaCorrSmall
                          && areSeparateCorr)                            cSepCorr++;

            // SC8 (vertex positions) --------------------------------------
            TVector3 vertexLeadingKaon    = leadingKaon.decayVertex();
            TVector3 vertexSubLeadingKaon = subLeadingKaon.decayVertex();
            double zdiff = vertexLeadingKaon.Z() - vertexSubLeadingKaon.Z();
            double zmean = (vertexLeadingKaon.Z() + vertexSubLeadingKaon.Z()) / 2.0;

            bool isZmeanSmall = (abs(zmean) <= config.zmeanMax);
            bool isZdiffSmall = (abs(zdiff) <= config.zdiffMax);
            vCuts.push_back(&isZmeanSmall);
            vCuts.push_back(&isZdiffSmall);

            bool passedSC7 = passedSC6 && isKsiCorrSmall && isEtaCorrSmall && areSeparateCorr;
            if (passedSC7 && isZmeanSmall)                  cZmean++;
            if (passedSC7 && isZmeanSmall && isZdiffSmall)  cZdiff++;

            // SC9 (cos theta*) --------------------------------------------
            double cosThetaStarLeading    = leadingKaon.cosThetaStar();
            double cosThetaStarSubLeading = subLeadingKaon.cosThetaStar();

            bool areCosThetaStarSmall = (abs(cosThetaStarLeading)    <= config.cosThetaStarMax &&
                                         abs(cosThetaStarSubLeading) <= config.cosThetaStarMax);
            vCuts.push_back(&areCosThetaStarSmall);

            bool passedSC8 = passedSC7 && isZmeanSmall && isZdiffSmall;
            if (passedSC8 && areCosThetaStarSmall) cCosThetaStar++;

            // ---- N-1 histograms ----------------------------------------
            bool N1_mass = true; CheckN1(vCuts, areKaonsInNarrowMassWindow, N1_mass);
            if (N1_mass) { HistInvMassPiPiN1[iH]->Fill(leadingKaonMass); HistInvMassPiPiN1[iH]->Fill(subleadingKaonMass); }

            bool N1_ptmiss = true; CheckN1(vCuts, isPtMissingSmall, N1_ptmiss);
            if (N1_ptmiss) {
                HistPtMissN1[iH]->Fill(pTmiss);
                if (isPtMissingSmall) HistMassK0K0->Fill(massK0K0);
            }

            bool N1_cluster = true; CheckN1(vCuts, isNumberOfTofClusterSmall, N1_cluster);
            if (N1_cluster) HistNTOFClusterN1[iH]->Fill(totalCluster);

            bool N1_dcaBL = true; CheckN1(vCuts, areBothDcaBeamlineSmall, N1_dcaBL);
            if (N1_dcaBL) { HistDCABeamlineN1[iH]->Fill(DCABeamlineLeadingKaon); HistDCABeamlineN1[iH]->Fill(DCABeamlineSubLeadingKaon); }

            bool N1_dcaDau = true; CheckN1(vCuts, areBothDcaDaughtersSmall, N1_dcaDau);
            if (N1_dcaDau) { HistDCADaughtersN1[iH]->Fill(dcaDaughtersLeadingKaon); HistDCADaughtersN1[iH]->Fill(dcaDaughtersSubLeadingKaon); }

            bool N1_pa = true; CheckN1(vCuts, areBothPointingAngleOrDecayLength, N1_pa);
            if (N1_pa) {
                HistCosN1[iH]->Fill(cosLeadingKaon); HistCosN1[iH]->Fill(cosSubLeadingKaon);
                HistDecayN1[iH]->Fill(decayLengthLeadingKaon); HistDecayN1[iH]->Fill(decayLengthSubLeadingKaon);
            }

            bool N1_ae = true; CheckN1(vCuts, antiElastic, N1_ae);
            if (N1_ae) {
                HistSumProtonMomentaXN1[iH]->Fill(sumProtonMomentumX);
                HistSumProtonMomentaYN1[iH]->Fill(sumProtonMomentumY);
                HistKsiEN1[iH]->Fill(ksiE); HistKsiWN1[iH]->Fill(ksiW);
            }

            bool N1_ksi = true; CheckN1(vCuts, isKsiCorrSmall, N1_ksi);
            if (N1_ksi) HistCorrKsiN1[iH]->Fill(ksiCorrelation);

            bool N1_eta = true; CheckN1(vCuts, isEtaCorrSmall, N1_eta);
            if (N1_eta) HistCorrEtaN1[iH]->Fill(etaCorrelation);

            bool N1_sep = true; CheckN1(vCuts, areSeparateCorr, N1_sep);
            if (N1_sep) { HistCorrEN1[iH]->Fill(corrE); HistCorrWN1[iH]->Fill(corrW); }

            bool N1_zmean = true; CheckN1(vCuts, isZmeanSmall, N1_zmean);
            if (N1_zmean) HistZMeanN1[iH]->Fill(zmean);

            bool N1_zdiff = true; CheckN1(vCuts, isZdiffSmall, N1_zdiff);
            if (N1_zdiff) HistZDiffN1[iH]->Fill(zdiff);

            bool N1_costh = true; CheckN1(vCuts, areCosThetaStarSmall, N1_costh);
            if (N1_costh) { HistCosThetaStarN1[iH]->Fill(cosThetaStarLeading); HistCosThetaStarN1[iH]->Fill(cosThetaStarSubLeading); }

            // ---- BEFORE histograms -------------------------------------
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

            // ---- Cutflow (11 steps) ------------------------------------
            // Step 0: mass window
            // Step 1: + TOF clusters
            // Step 2: + DCA beamline
            // Step 3: + DCA daughters
            // Step 4: + pointing angle/decay length
            // Step 5: + anti-elastic
            // Step 6: + ksi correlation
            // Step 7: + eta correlation
            // Step 8: + separate corr
            // Step 9: + zmean
            // Step 10: + zdiff
            // (SC9 cosThetaStar is enforced via vCuts in the N-1 above)

            if (areKaonsInNarrowMassWindow) {
                HistPtMissCF[iH][0]->Fill(pTmiss); HistInvMassPiPiCF[iH][0]->Fill(leadingKaonMass); HistInvMassPiPiCF[iH][0]->Fill(subleadingKaonMass); HistInvMassKKCF[iH][0]->Fill(massK0K0);
                if (isNumberOfTofClusterSmall) {
                    HistPtMissCF[iH][1]->Fill(pTmiss); HistInvMassPiPiCF[iH][1]->Fill(leadingKaonMass); HistInvMassPiPiCF[iH][1]->Fill(subleadingKaonMass); HistInvMassKKCF[iH][1]->Fill(massK0K0);
                    if (areBothDcaBeamlineSmall) {
                        HistPtMissCF[iH][2]->Fill(pTmiss); HistInvMassPiPiCF[iH][2]->Fill(leadingKaonMass); HistInvMassPiPiCF[iH][2]->Fill(subleadingKaonMass); HistInvMassKKCF[iH][2]->Fill(massK0K0);
                        if (areBothDcaDaughtersSmall) {
                            HistPtMissCF[iH][3]->Fill(pTmiss); HistInvMassPiPiCF[iH][3]->Fill(leadingKaonMass); HistInvMassPiPiCF[iH][3]->Fill(subleadingKaonMass); HistInvMassKKCF[iH][3]->Fill(massK0K0);
                            if (areBothPointingAngleOrDecayLength) {
                                HistPtMissCF[iH][4]->Fill(pTmiss); HistInvMassPiPiCF[iH][4]->Fill(leadingKaonMass); HistInvMassPiPiCF[iH][4]->Fill(subleadingKaonMass); HistInvMassKKCF[iH][4]->Fill(massK0K0);
                                if (antiElastic) {
                                    HistPtMissCF[iH][5]->Fill(pTmiss); HistInvMassPiPiCF[iH][5]->Fill(leadingKaonMass); HistInvMassPiPiCF[iH][5]->Fill(subleadingKaonMass); HistInvMassKKCF[iH][5]->Fill(massK0K0);
                                    if (isKsiCorrSmall) {
                                        HistPtMissCF[iH][6]->Fill(pTmiss); HistInvMassPiPiCF[iH][6]->Fill(leadingKaonMass); HistInvMassPiPiCF[iH][6]->Fill(subleadingKaonMass); HistInvMassKKCF[iH][6]->Fill(massK0K0);
                                        if (isEtaCorrSmall) {
                                            HistPtMissCF[iH][7]->Fill(pTmiss); HistInvMassPiPiCF[iH][7]->Fill(leadingKaonMass); HistInvMassPiPiCF[iH][7]->Fill(subleadingKaonMass); HistInvMassKKCF[iH][7]->Fill(massK0K0);
                                            if (areSeparateCorr) {
                                                HistPtMissCF[iH][8]->Fill(pTmiss); HistInvMassPiPiCF[iH][8]->Fill(leadingKaonMass); HistInvMassPiPiCF[iH][8]->Fill(subleadingKaonMass); HistInvMassKKCF[iH][8]->Fill(massK0K0);
                                                if (isZmeanSmall) {
                                                    HistPtMissCF[iH][9]->Fill(pTmiss); HistInvMassPiPiCF[iH][9]->Fill(leadingKaonMass); HistInvMassPiPiCF[iH][9]->Fill(subleadingKaonMass); HistInvMassKKCF[iH][9]->Fill(massK0K0);
                                                    if (isZdiffSmall) {
                                                        HistPtMissCF[iH][10]->Fill(pTmiss); HistInvMassPiPiCF[iH][10]->Fill(leadingKaonMass); HistInvMassPiPiCF[iH][10]->Fill(subleadingKaonMass); HistInvMassKKCF[iH][10]->Fill(massK0K0);
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
            }

        } // end event loop

        // =========================================================
        // ANALYSIS SUMMARY
        // =========================================================
        long long total = chain->GetEntries();  // raw chain entries (before SC2 note)
        auto pct = [&](long long n) -> double {
            return cInput > 0 ? 100.0 * n / cInput : 0.0;
        };
        auto rel = [](long long n, long long prev) -> double {
            return prev > 0 ? 100.0 * n / prev : 0.0;
        };

        // Print N_TOF selection used
        cout << "\n ➡️  N_TOF selected: {";
        for (int k = 0; k < (int)N_TOF_SELECTED.size(); k++) {
            cout << N_TOF_SELECTED[k];
            if (k < (int)N_TOF_SELECTED.size()-1) cout << ", ";
        }
        cout << "}" << endl;

        cout << "\n";
        cout << "╔══════════════════════════════════════════════════════════════════════════╗" << endl;
        cout << "║ 👀 ANALYSIS SUMMARY 👀  —  " << config.name << endl;
        cout << "╠══════════════════════════════════════════════════════════════════════════╣" << endl;
        cout << Form("║  %-42s  %10lld             ║", "Raw chain entries",              total)                                          << endl;
        cout << "╠══════════════════════════════════════════════════════════════════════════╣" << endl;
        cout << "║  NOTE: SC 2 already applied in Preselection.cxx:                         ║" << endl;
        cout << "║    SC2.1  CEP trigger (570701, 570705, 570711)                           ║" << endl;
        cout << "║    SC2.2  Exactly 1 RP track per side (E+W)                              ║" << endl;
        cout << "║    SC2.3  >=3/4 RP planes used per trackpoint                            ║" << endl;
        cout << "║    SC2.4  Proton in fiducial region                                      ║" << endl;
        cout << "║    SC2.5  >=2 TOF class-A tracks                                         ║" << endl;
        cout << "║    SC3.0 (1 primary vertex, |z|<80 cm) applied here in EventAnalyzer.    ║" << endl;
        cout << "╠══════════════════════════════════════════════════════════════════════════╣" << endl;
        cout << Form("║  %-42s  %10lld  (100.00%%)  ║", "Input after preselection (SC2)", cInput)                                       << endl;
        cout << "╠══════════════════════════════════════════════════════════════════════════╣" << endl;
        cout << "║  -- SC 3: Track quality --                                               ║" << endl;
        cout << Form("║  %-42s  %10lld  (%6.2f%%)  rel:%5.1f%%  ║", "SC3.0: 1 vtx, |z_vtx|<80 cm",            cVertex,       pct(cVertex),       rel(cVertex,        cInput))        << endl;
        cout << Form("║  %-42s  %10lld  (%6.2f%%)  rel:%5.1f%%  ║", "SC3:   N_TOF in selected set",           cValidTracks,  pct(cValidTracks),  rel(cValidTracks,   cVertex))       << endl;
        cout << Form("║  %-42s  %10lld  (%6.2f%%)  rel:%5.1f%%  ║", "SC3.1-3.3: pT>0.2,|eta|<0.9,Nfit>=20",   cTofGood,      pct(cTofGood),      rel(cTofGood,       cValidTracks))  << endl;
        cout << Form("║  %-42s  %10lld  (%6.2f%%)  rel:%5.1f%%  ║", "SC3.4: charge balance + nHitsDedx>=15",  cComplementary,pct(cComplementary),rel(cComplementary, cTofGood))     << endl;
        cout << Form("║  %-42s  %10lld  (%6.2f%%)  rel:%5.1f%%  ║", "       Pion pairs reconstructed",        cPionsFound,   pct(cPionsFound),   rel(cPionsFound,    cComplementary)) << endl;
        cout << "╠══════════════════════════════════════════════════════════════════════════╣" << endl;
        cout << "║  -- SC 4: K0S exclusivity --                                             ║" << endl;
        cout << Form("║  %-42s  %10lld  (%6.2f%%)  rel:%5.1f%%  ║", "SC4.1: K0S mass window",                  cMassWindow,  pct(cMassWindow),  rel(cMassWindow,  cPionsFound))   << endl;
        cout << Form("║  %-42s  %10lld  (%6.2f%%)  rel:%5.1f%%  ║", "SC4.2: + pT miss",                        cPtMiss,      pct(cPtMiss),      rel(cPtMiss,      cMassWindow))   << endl;
        cout << Form("║  %-42s  %10lld  (%6.2f%%)  rel:%5.1f%%  ║", "SC4.3: + TOF clusters",                   cTofClusters, pct(cTofClusters), rel(cTofClusters, cPtMiss))       << endl;
        cout << "╠══════════════════════════════════════════════════════════════════════════╣" << endl;
        cout << "║  -- SC 5: K0S topology --                                                ║" << endl;
        cout << Form("║  %-42s  %10lld  (%6.2f%%)  rel:%5.1f%%  ║", "SC5:   + DCA daughters",                  cDcaDaughters,pct(cDcaDaughters),rel(cDcaDaughters, cTofClusters)) << endl;
        cout << Form("║  %-42s  %10lld  (%6.2f%%)  rel:%5.1f%%  ║", "SC5:   + DCA beamline",                   cDcaBeamline, pct(cDcaBeamline), rel(cDcaBeamline, cDcaDaughters)) << endl;
        cout << Form("║  %-42s  %10lld  (%6.2f%%)  rel:%5.1f%%  ║", "SC5:   + pointing angle / decay length",  cPointing,    pct(cPointing),    rel(cPointing,    cDcaBeamline))  << endl;
        cout << "╠══════════════════════════════════════════════════════════════════════════╣" << endl;
        cout << "║  -- SC 6: Anti-elastic --                                                ║" << endl;
        cout << Form("║  %-42s  %10lld  (%6.2f%%)  rel:%5.1f%%  ║", "SC6:   + anti-elastic (ksi/momenta)",     cAntiElastic, pct(cAntiElastic), rel(cAntiElastic, cPointing))     << endl;
        cout << "╠══════════════════════════════════════════════════════════════════════════╣" << endl;
        cout << "║  -- SC 7: Kinematic correlations --                                      ║" << endl;
        cout << Form("║  %-42s  %10lld  (%6.2f%%)  rel:%5.1f%%  ║", "SC7.1: + ksi correlation <0.004",         cKsiCorr,     pct(cKsiCorr),     rel(cKsiCorr,     cAntiElastic))  << endl;
        cout << Form("║  %-42s  %10lld  (%6.2f%%)  rel:%5.1f%%  ║", "SC7.2: + eta correlation <0.9",           cEtaCorr,     pct(cEtaCorr),     rel(cEtaCorr,     cKsiCorr))      << endl;
        cout << Form("║  %-42s  %10lld  (%6.2f%%)  rel:%5.1f%%  ║", "SC7.3: + corrE, corrW <0.004",            cSepCorr,     pct(cSepCorr),     rel(cSepCorr,     cEtaCorr))      << endl;
        cout << "╠══════════════════════════════════════════════════════════════════════════╣" << endl;
        cout << "║  -- SC 8: Decay vertex positions --                                      ║" << endl;
        cout << Form("║  %-42s  %10lld  (%6.2f%%)  rel:%5.1f%%  ║", "SC8.1: + |zmean| <=80 cm",                cZmean,       pct(cZmean),       rel(cZmean,       cSepCorr))      << endl;
        cout << Form("║  %-42s  %10lld  (%6.2f%%)  rel:%5.1f%%  ║", "SC8.2: + |zdiff| <=15 cm",                cZdiff,       pct(cZdiff),       rel(cZdiff,       cZmean))        << endl;
        cout << "╠══════════════════════════════════════════════════════════════════════════╣" << endl;
        cout << "║  -- SC 9: Angular distribution --                                        ║" << endl;
        cout << Form("║  %-42s  %10lld  (%6.2f%%)  rel:%5.1f%%  ║", "SC9:   + |cos(theta*)| <=0.8",            cCosThetaStar,pct(cCosThetaStar),rel(cCosThetaStar,cZdiff))        << endl;
        cout << "╠══════════════════════════════════════════════════════════════════════════╣" << endl;
        cout << Form("║  %-42s  %10lld                       ║", ">>> FINAL EXCLUSIVE SAMPLE <<<",              cCosThetaStar)                                                      << endl;
        cout << "╚══════════════════════════════════════════════════════════════════════════╝" << endl;
        cout << "\n";
        cout << " ✍ Writing histograms for: " << config.name << endl;

        HistNumWithTofTrakcs->Write();
        HistNumWithoutTofTrakcs->Write();

        for (int i = 0; i < 4; ++i) {
            HistInvMassPiPi2D[i]->Write();
            HistPtPionWithTof[i]->Write();   HistPtPionWithoutTof[i]->Write();
            HistEtaPionWithTof[i]->Write();  HistEtaPionWithoutTof[i]->Write();
            HistNfitPionWithTof[i]->Write(); HistNfitPionWithoutTof[i]->Write();

            HistInvMassPiPiN1[i]->Write();
            HistPtMissN1[i]->Write();        HistNTOFClusterN1[i]->Write();
            HistDCABeamlineN1[i]->Write();   HistDCADaughtersN1[i]->Write();
            HistCosN1[i]->Write();           HistDecayN1[i]->Write();
            HistSumProtonMomentaXN1[i]->Write(); HistSumProtonMomentaYN1[i]->Write();
            HistCorrEN1[i]->Write();  HistCorrWN1[i]->Write();
            HistCorrKsiN1[i]->Write(); HistCorrEtaN1[i]->Write();
            HistKsiEN1[i]->Write();   HistKsiWN1[i]->Write();
            HistZMeanN1[i]->Write();  HistZDiffN1[i]->Write();
            HistCosThetaStarN1[i]->Write();

            HistPtMissBefore[i]->Write();       HistNTOFClusterBefore[i]->Write();
            HistDCABeamlineBefore[i]->Write();  HistDCADaughtersBefore[i]->Write();
            HistCosBefore[i]->Write();          HistDecayBefore[i]->Write();
            HistCorrEBefore[i]->Write();        HistCorrWBefore[i]->Write();
            HistCorrKsiBefore[i]->Write();      HistCorrEtaBefore[i]->Write();
            HistKsiEBefore[i]->Write();         HistKsiWBefore[i]->Write();
            HistSumProtonMomentaXBefore[i]->Write(); HistSumProtonMomentaYBefore[i]->Write();
            HistZDiffBefore[i]->Write();        HistZMeanBefore[i]->Write();
            HistCosThetaStarBefore[i]->Write();

            for (int j = 0; j < 11; j++) {
                HistPtMissCF[i][j]->Write();
                HistInvMassPiPiCF[i][j]->Write();
                HistInvMassKKCF[i][j]->Write();
            }

            HistPtPionWithTofN1[i]->Write();    HistPtPionWithoutTofN1[i]->Write();
            HistEtaPionWithTofN1[i]->Write();   HistEtaPionWithoutTofN1[i]->Write();
            HistNfitPionWithTofN1[i]->Write();  HistNfitPionWithoutTofN1[i]->Write();
        }
        HistMassK0K0->Write();

    } // end systematic loop

    outfile->Close();
    cout << "\n ✅ All done. Output: " << argv[2] << endl;
    return 0;
}
