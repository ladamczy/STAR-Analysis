#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TMath.h>
#include <TVector3.h>
#include "StUPCEvent.h"
#include "StUPCTrack.h"
#include "StUPCV0.h"

using namespace std;

int main(int argc, char** argv) {
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " [input file data] [output file]" << endl;
        return 1;
    }

    TFile inputFile(argv[1], "READ");
    if (!inputFile.IsOpen()) {
        cerr << "Failed to open input file." << endl;
        return 1;
    }

    TTree* tree = static_cast<TTree*>(inputFile.Get("mUPCTree"));
    if (!tree) {
        cerr << "Failed to retrieve mUPCTree." << endl;
        inputFile.Close();
        return 1;
    }

    StUPCEvent* event = nullptr;
    tree->SetBranchAddress("mUPCEvent", &event);

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

    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);

        for (int j = 0; j < event->getNumberOfTracks(); ++j) {
            StUPCTrack* track1 = event->getTrack(j);
            for (int k = j + 1; k < event->getNumberOfTracks(); ++k) {
                StUPCTrack* track2 = event->getTrack(k);

                // Selekcja dla Lambda
                if ((track1->getCharge() > 0 && track2->getCharge() < 0) &&
                    track1->getNhits() > 15 && track2->getNhits() > 15 &&
                    track1->getPt() > 0.15 && track2->getPt() > 0.15 &&
                    abs(track1->getEta()) < 1.1 && abs(track2->getEta()) < 1.1 &&
                    (track1->getFlag(StUPCTrack::kTof) || track2->getFlag(StUPCTrack::kTof))) {

                    TVector3 const tryVec(0,0,0);
                    double beamLine[] = {0,0,0,0};

                    StUPCV0 v0(track1, track2, protonMass, pionMass, 1, 1, tryVec, beamLine, event->getMagneticField(), true);
                    if (v0.m() > lowerLimitOfInvMassLambda && v0.m() < upperLimitOfInvMassLambda) {
                        TVector3 pVec(v0.px(), v0.py(), v0.pz());
                        float pT = pVec.Perp();
                        TVector3 decayVec = v0.decayVertex() - tryVec;
                        float alpha = (decayVec.Dot(pVec)) / (decayVec.Mag() * pVec.Mag());
                        hArmenteros->Fill(alpha, pT);
                    }
                }

                // Selekcja dla LambdaBar
                if ((track1->getCharge() < 0 && track2->getCharge() > 0) &&
                    track1->getNhits() > 15 && track2->getNhits() > 15 &&
                    track1->getPt() > 0.15 && track2->getPt() > 0.15 &&
                    abs(track1->getEta()) < 1.1 && abs(track2->getEta()) < 1.1 &&
                    (track1->getFlag(StUPCTrack::kTof) || track2->getFlag(StUPCTrack::kTof))) {

                    TVector3 const tryVec(0,0,0);
                    double beamLine[] = {0,0,0,0};

                    StUPCV0 v0(track1, track2, protonMass, pionMass, 1, 1, tryVec, beamLine, event->getMagneticField(), true);
                    if (v0.m() > lowerLimitOfInvMassLambda && v0.m() < upperLimitOfInvMassLambda) {
                        TVector3 pVec(v0.px(), v0.py(), v0.pz());
                        float pT = pVec.Perp();
                        TVector3 decayVec = v0.decayVertex() - tryVec;
                        float alpha = (decayVec.Dot(pVec)) / (decayVec.Mag() * pVec.Mag());
                        hArmenteros->Fill(alpha, pT);
                    }
                }


                // // Selekcja dla K0
                // if (track1->getCharge() != track2->getCharge() &&
                //     track1->getNhits() > 15 && track2->getNhits() > 15 &&
                //     track1->getPt() > 0.15 && track2->getPt() > 0.15 &&
                //     abs(track1->getEta()) < 1.1 && abs(track2->getEta()) < 1.1 &&
                //     (track1->getFlag(StUPCTrack::kTof) || track2->getFlag(StUPCTrack::kTof))) {
                    
                //     TVector3 const tryVec(0,0,0);
                //     double beamLine[] = {0,0,0,0};

                //     StUPCV0 v0(track1, track2, pionMass, pionMass, 1, 1, tryVec, beamLine, event->getMagneticField(), true);
                //     if (v0.m() > lowerLimitOfInvMassK0 && v0.m() < upperLimitOfInvMassK0) {
                //         TVector3 pVec(v0.px(), v0.py(), v0.pz());
                //         float pT = pVec.Perp();
                //         TVector3 decayVec = v0.decayVertex() - TVector3(0,0,0);
                //         float alpha = (decayVec.Dot(pVec)) / (decayVec.Mag() * pVec.Mag());
                //         hArmenteros->Fill(alpha, pT);
                //     }
                // }
            }
        }
    }

    TFile outputFile(argv[2], "RECREATE");
    hArmenteros->Write();
    outputFile.Close();

    inputFile.Close();

    return 0;
}