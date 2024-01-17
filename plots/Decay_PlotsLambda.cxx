#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TMath.h>
#include "StUPCEvent.h"
#include "StUPCTrack.h"
#include "StUPCV0.h"

using namespace std;

int main(int argc, char** argv) {

    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " [input file] [output file]" << endl;
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

    const float pionMass = 0.13957;                  // Mass of the pion in GeV/c^2
    const float protonMass = 0.938;                  // Mass of the proton in GeV/c^2

    const double centralOfInvMassLambda = 1.115;     // Central value of the invariant mass of Lambda (LambdaBar)
    const double deltaOfInvMass = 0.035;             // Delta value of the invariant mass Lambda (LambdaBar)
    double lowerLimitOfInvMassLambda = centralOfInvMassLambda - deltaOfInvMass;
    double upperLimitOfInvMassLambda = centralOfInvMassLambda + deltaOfInvMass;

    TH1D* histDecayLengthHypo = new TH1D("DecayLengthHypo", ";Length (cm);Events", 100, 0, 9);

    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);

        for (int j = 0; j < event->getNumberOfTracks(); ++j) {
            StUPCTrack* track1 = event->getTrack(j);
            for (int k = j + 1; k < event->getNumberOfTracks(); ++k) {
                StUPCTrack* track2 = event->getTrack(k);

                if ((track1->getCharge() > 0 && track2->getCharge() < 0) &&  // track1 is proton and track2 is negative pion
                    track1->getNhits() > 15 && track2->getNhits() > 15 &&
                    track1->getPt() > 0.15 && track2->getPt() > 0.15 &&
                    abs(track1->getEta()) < 1.1 && abs(track2->getEta()) < 1.1 &&
                    (track1->getFlag(StUPCTrack::kTof) && track2->getFlag(StUPCTrack::kTof)) &&
                    abs(track1->getNSigmasTPCProton()) < 2.5 && abs(track2->getNSigmasTPCPion()) < 2.5) {
                    
                    TVector3 const tryVec(0,0,0);
                    double beamLine[] = {0,0,0,0};
                    StUPCV0 v0(track1, track2, protonMass, pionMass, 1, 1, tryVec, beamLine, event->getMagneticField(), true);

                    if (v0.dcaDaughters() < 1.5) {
                        if (v0.DCABeamLine() < 1.5) {
                            if (v0.pointingAngleHypo() < 0.925) {
                                if (v0.m() > lowerLimitOfInvMassLambda && v0.m() < upperLimitOfInvMassLambda) {
                                    histDecayLengthHypo->Fill(v0.decayLengthHypo());
                                }
                            }

                        }
                    }

                }
            }
        }
    }

    histDecayLengthHypo->SetStats(0);

    TFile outputFile(argv[2], "RECREATE");
    histDecayLengthHypo->Write();

    outputFile.Close();
    inputFile.Close();

    return 0;
}
