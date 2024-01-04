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
    const double deltaOfInvMass = 0.035;             // Delta value of the invariant mass of Lambda (LambdaBar)

    double lowerLimitOfInvMassLambda = centralOfInvMassLambda - deltaOfInvMass;
    double upperLimitOfInvMassLambda = centralOfInvMassLambda + deltaOfInvMass;
    
    TH1D* histDcaDaughters = new TH1D("DcaDaughters", ";DCA (cm);Counts", 100, 0, 5);
    TH1D* histDcaBeamLine = new TH1D("DcaBeamLine", ";DCA (cm);Counts", 100, 0, 2.5);
    TH1D* histPointingAngleHypo = new TH1D("PointingAngleHypo", ";Cos(#theta);Counts", 100, -1, 1);
    TH1D* histInvariantMassLambdaBar = new TH1D("InvariantMassLambdaBar", ";Mass (GeV/c^2);Counts", 100, 1.08, 1.15);
    TH1D* histDecayLengthHypo = new TH1D("DecayLengthHypo", ";Length (cm);Counts", 100, 0, 20);

    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);

        for (int j = 0; j < event->getNumberOfTracks(); ++j) {
            StUPCTrack* track1 = event->getTrack(j);
            for (int k = j + 1; k < event->getNumberOfTracks(); ++k) {
                StUPCTrack* track2 = event->getTrack(k);

                if (track1->getCharge() < 0 && track2->getCharge() > 0 &&   // track1 is antiproton and track2 is positive pion
                    track1->getNhits() > 15 && track2->getNhits() > 15 &&
                    track1->getPt() > 0.15 && track2->getPt() > 0.15 &&
                    abs(track1->getEta()) < 1.1 && abs(track2->getEta()) < 1.1 &&
                    (track1->getFlag(StUPCTrack::kTof) || track2->getFlag(StUPCTrack::kTof))) {
                    
                    TVector3 const tryVec(0,0,0);
                    double beamLine[] = {0,0,0,0};
                    StUPCV0 v0(track1, track2, protonMass, pionMass, 1, 1, tryVec, beamLine, event->getMagneticField(), true);

                    double invMass = v0.m();

                    if (v0.DCABeamLine() < 1.5) {
                        if (v0.pointingAngleHypo() > 0.925) {
                            if (invMass > lowerLimitOfInvMassLambda && invMass < upperLimitOfInvMassLambda) {
                                histDcaDaughters->Fill(v0.dcaDaughters());
                            }
                        }
                    }

                    if (v0.dcaDaughters() < 1.5) {
                        if (v0.DCABeamLine() < 1.5) {
                            if (v0.pointingAngleHypo() > 0.925) {
                                if (invMass > lowerLimitOfInvMassLambda && invMass < upperLimitOfInvMassLambda) {
                                    histDecayLengthHypo->Fill(v0.decayLengthHypo());
                                }
                                histInvariantMassLambdaBar->Fill(invMass);
                            }
                            if (invMass > lowerLimitOfInvMassLambda && invMass < upperLimitOfInvMassLambda) {
                                histPointingAngleHypo->Fill(v0.pointingAngleHypo());
                            }
                        }
                        if (v0.pointingAngleHypo() > 0.925) {
                            if (invMass > lowerLimitOfInvMassLambda && invMass < upperLimitOfInvMassLambda) {
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
    histInvariantMassLambdaBar->Write();
    histDecayLengthHypo->Write();

    outputFile.Close();
    inputFile.Close();

    return 0;
}
