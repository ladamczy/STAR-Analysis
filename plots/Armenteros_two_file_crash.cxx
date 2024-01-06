#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TMath.h>
#include "StUPCEvent.h"
#include "StUPCTrack.h"
#include "StUPCV0.h"

using namespace std;

void AnalyzeData(TTree* tree, TH2D** hArmenteros, float mass1, float mass2, float lowerLimitOfInvMassV0, float upperLimitOfInvMassV0);

int main(int argc, char** argv) {

    if (argc != 4) {
        cerr << "Usage: " << argv[0] << " [input file K0 data] [input file Lambda data] [output file]" << endl;
        return 1;
    }

    TFile* inputFileK0 = new TFile(argv[1], "READ");
    TFile* inputFileLambda = new TFile(argv[2], "READ");
    if (!inputFileK0->IsOpen() || !inputFileLambda->IsOpen()) {
        cerr << "Failed to open one of the input files." << endl;
        inputFileK0->Close();
        inputFileLambda->Close();
        return 1;
    }

    TTree* treeK0 = static_cast<TTree*>(inputFileK0->Get("mUPCTree"));
    TTree* treeLambda = static_cast<TTree*>(inputFileLambda->Get("mUPCTree"));
    if (!treeK0 || !treeLambda) {
        cerr << "Failed to retrieve mUPCTree from one of the input files." << endl;
        inputFileK0->Close();
        inputFileLambda->Close();
        return 1;
    }

    TH2D* hArmenteros = new TH2D("hArmenteros", "Armenteros Plot;alpha;pT (GeV/c)", 200, -1, 1, 200, 0, 1);

    const float pionMass = 0.13957; // Mass of the pion in GeV/c^2
    const float protonMass = 0.938; // Mass of the proton in GeV/c^2
    const double centralOfInvMassLambda = 1.115; // Central value of Lambda invariant mass
    const double deltaOfInvMass = 0.035; // Delta of invariant mass for Lambda
    const double centralOfInvMassK0 = 0.495; // Central value of K0 invariant mass
    const double deltaOfInvMassK0 = 0.035; // Delta of invariant mass for K0

    double lowerLimitOfInvMassLambda = centralOfInvMassLambda - deltaOfInvMass;
    double upperLimitOfInvMassLambda = centralOfInvMassLambda + deltaOfInvMass;
    double lowerLimitOfInvMassK0 = centralOfInvMassK0 - deltaOfInvMassK0;
    double upperLimitOfInvMassK0 = centralOfInvMassK0 + deltaOfInvMassK0;

    AnalyzeData(treeK0, &hArmenteros, pionMass, pionMass, lowerLimitOfInvMassK0, upperLimitOfInvMassK0);
    AnalyzeData(treeLambda, &hArmenteros, protonMass, pionMass, lowerLimitOfInvMassLambda, upperLimitOfInvMassLambda);

    TFile outputFile(argv[3], "RECREATE");
    hArmenteros->Write();
    outputFile.Close();

    inputFileK0->Close();
    inputFileLambda->Close();

    return 0;
}



void AnalyzeData(TTree* tree, TH2D** hArmenteros, float mass1, float mass2, float lowerLimitOfInvMassV0, float upperLimitOfInvMassV0) {
    StUPCEvent* event = nullptr;
    tree->SetBranchAddress("mUPCEvent", &event);

    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);

        for (int j = 0; j < event->getNumberOfTracks(); ++j) {
            StUPCTrack* track1 = event->getTrack(j);
            for (int k = j + 1; k < event->getNumberOfTracks(); ++k) {
                StUPCTrack* track2 = event->getTrack(k);

                if ((track1->getCharge() > 0 && track2->getCharge() < 0) &&
                    track1->getNhits() > 15 && track2->getNhits() > 15 &&
                    track1->getPt() > 0.15 && track2->getPt() > 0.15 &&
                    abs(track1->getEta()) < 1.1 && abs(track2->getEta()) < 1.1 &&
                    (track1->getFlag(StUPCTrack::kTof) || track2->getFlag(StUPCTrack::kTof))) {
                    
                    StUPCV0 v0(track1, track2, mass1, mass2, j, k, TVector3(0,0,0), nullptr, event->getMagneticField(), true);

                    if (v0.m() > lowerLimitOfInvMassV0 && v0.m() < upperLimitOfInvMassV0) {
                        TVector3 pVec(v0.px(), v0.py(), v0.pz());
                        float pT = pVec.Perp();
                        float alpha = (2 * pVec.Dot(v0.decayVertex())) / (v0.lorentzVector().M2() - mass1*mass1 - mass2*mass2);
                        (*hArmenteros)->Fill(alpha, pT);
                    }
                }
            }
        }
    }
}