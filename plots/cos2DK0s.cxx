#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TMath.h>
#include "StUPCEvent.h"
#include "StUPCTrack.h"
#include "StUPCV0.h"
#include <TROOT.h> 
#include <TCanvas.h>

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
    
    const float pionMass = 0.13957;
    const float protonMass = 0.938;

    const double centralOfInvMassK0 = 0.495;
    const double deltaOfInvMass = 0.035;
    double lowerLimitOfInvMassK0 = centralOfInvMassK0 - deltaOfInvMass;
    double upperLimitOfInvMassK0 = centralOfInvMassK0 + deltaOfInvMass;

    // histogram 2D pokazujacy zależność cosPointingAngle od DecayLength.

    TH2D* cosPointingAngle_DecayLength = new TH2D("cosPointingAngle_DecayLength", ";cos(#theta);Length (cm)", 100, -1, 1, 100, 0, 20);

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
                    (track1->getFlag(StUPCTrack::kTof) && track2->getFlag(StUPCTrack::kTof)) &&
                    abs(track1->getNSigmasTPCPion()) < 2.5 && abs(track2->getNSigmasTPCPion()) < 2.5){
                    
                    TVector3 const tryVec(0,0,0);
                    double beamLine[] = {0,0,0,0};
                    StUPCV0 v0(track1, track2, pionMass, pionMass, 1, 1, tryVec, beamLine, event->getMagneticField(), true);
                   
                    if (v0.m() > lowerLimitOfInvMassK0 && v0.m() < upperLimitOfInvMassK0) {
                        float cospointingAngleHypo = v0.pointingAngleHypo();
                        float decayLengthHypo = v0.decayLengthHypo();
                        cosPointingAngle_DecayLength->Fill(cospointingAngleHypo, decayLengthHypo);
                    }
                }
            }
        }
    }



    TCanvas* canvas = new TCanvas("canvas", "Canvas", 2560, 1440);
    cosPointingAngle_DecayLength->SetStats(0);
    cosPointingAngle_DecayLength->Draw("colz");
    canvas->SetLogz();
    canvas->SaveAs("cosPointingAngle_DecayLength_K0_log.jpg");

    TFile outputFile(argv[2], "RECREATE");
    cosPointingAngle_DecayLength->Write();


    outputFile.Close();
    inputFile.Close();

    return 0;
}
