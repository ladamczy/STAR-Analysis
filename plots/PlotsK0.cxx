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

    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        for (int j = 0; j < event->getNumberOfTracks(); ++j) {
            StUPCTrack* track1 = event->getTrack(j);
            for (int k = j + 1; k < event->getNumberOfTracks(); ++k) {
                StUPCTrack* track2 = event->getTrack(k);

                if (track1->getCharge() != track2->getCharge() &&
                    track1->getNhits() > 15 && track2->getNhits() > 15 &&
                    track1->getPt() > 0.15 && track2->getPt() > 0.15 &&
                    abs(track1->getEta()) < 1.1 && abs(track2->getEta()) < 1.1 &&
                    (track1->getFlag(StUPCTrack::kTof) || track2->getFlag(StUPCTrack::kTof))) {
                    
                    TVector3 const tryVec(0,0,0);
                    double beamLine[] = {0,0,0,0};
                    StUPCV0 v0(track1, track2, pionMass, pionMass, 1, 1, tryVec, beamLine, event->getMagneticField(), true);
                   
                    double invMass = v0.m();

                    // Cuts to object V0

                    if (v0.DCABeamLine() < 1.5) {
                        if (v0.pointingAngleHypo() > 0.925) {
                            if (invMass > lowerLimitOfInvMassK0 && invMass < upperLimitOfInvMassK0) {
                                histDcaDaughters->Fill(v0.dcaDaughters());
                            }
                        }
                    }

                    if (v0.dcaDaughters() < 1.5) {
                        if (v0.DCABeamLine() < 1.5) {
                            if (v0.pointingAngleHypo() > 0.925) {
                                if (invMass > lowerLimitOfInvMassK0 && invMass < upperLimitOfInvMassK0) {
                                    histDecayLengthHypo->Fill(v0.decayLengthHypo());
                                }
                                histInvariantMassK0->Fill(invMass);
                            }
                            if (invMass > lowerLimitOfInvMassK0 && invMass < upperLimitOfInvMassK0) {
                                histPointingAngleHypo->Fill(v0.pointingAngleHypo());
                            }
                        }
                        if (v0.pointingAngleHypo() > 0.925) {
                            if (invMass > lowerLimitOfInvMassK0 && invMass < upperLimitOfInvMassK0) {
                                histDcaBeamLine->Fill(v0.DCABeamLine());
                            }
                        }
                    }
                }
            }
        }
    }

    TFile outputFile(argv[2], "RECREATE");
    histDcaDaughters->Write();
    histDcaBeamLine->Write();
    histPointingAngleHypo->Write();
    histInvariantMassK0->Write();
    histDecayLengthHypo->Write();

    outputFile.Close();
    inputFile.Close();

    return 0;
}
