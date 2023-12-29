#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
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

    TH1D* histPointingAngle = new TH1D("histPointingAngle", "Pointing Angle;Angle;Counts", 100, 0, TMath::Pi());
    TH1D* histDecayLength = new TH1D("histDecayLength", "Decay Length;Length (cm);Counts", 100, 0, 10);
    TH1D* histDcaDaughters = new TH1D("histDcaDaughters", "DCA between Daughters;DCA (cm);Counts", 100, 0, 1);

    const float pionMass = 0.13957; // Pion mass in GeV/c^2

    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);

        for (int j = 0; j < event->getNumberOfTracks(); ++j) {
            StUPCTrack* track1 = event->getTrack(j);
            for (int k = j + 1; k < event->getNumberOfTracks(); ++k) {
                StUPCTrack* track2 = event->getTrack(k);
                TVector3 const tryVec(0,0,0);
                double beamLine[] = {0,0,0,0};
                StUPCV0 v0(track1, track2, pionMass, pionMass, 1, 1, tryVec, beamLine , event->getMagneticField(), true);

                histPointingAngle->Fill(v0.pointingAngle());
                histDecayLength->Fill(v0.decayLength());
                histDcaDaughters->Fill(v0.dcaDaughters());
            }
        }
    }

    TFile outputFile(argv[2], "RECREATE");
    histPointingAngle->Write();
    histDecayLength->Write();
    histDcaDaughters->Write();
    outputFile.Close();
    inputFile.Close();

    return 0;
}
