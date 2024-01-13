#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TMath.h>
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
        cerr << "Failed to retrieve mUPCTree from the input file." << endl;
        inputFile.Close();
        return 1;
    }

    StUPCEvent* event = nullptr;
    tree->SetBranchAddress("mUPCEvent", &event);


    const float pionMass = 0.13957;
    const float protonMass = 0.938;

    double lowerLimitOfInvMassLambda = 1.08;
    double upperLimitOfInvMassLambda = 1.15;
    double lowerLimitOfInvMassK0 = 0.46;
    double upperLimitOfInvMassK0 = 0.53;

    TH2D* hArmenteros = new TH2D("hArmenteros", "Armenteros Plot;alpha;pT (GeV/c)", 200, -1, 1, 100, 0, 0.3);

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

                    // Tworzenie obiektu StUPCV0
                    TVector3 const tryVec(0,0,0);
                    double beamLine[] = {0,0,0,0};
                    StUPCV0 v0(track1, track2, pionMass, pionMass, 1, 1, tryVec, beamLine, event->getMagneticField(), true);


                    if (v0.m() > lowerLimitOfInvMassK0 && v0.m() < upperLimitOfInvMassK0) {
                        TVector3 p1Vec, p2Vec;
                        track1->getMomentum(p1Vec);
                        track2->getMomentum(p2Vec);

                        TVector3 v0Momentum = TVector3(v0.px(), v0.py(), v0.pz()).Unit();
                        TVector3 zAxis(0, 0, 1);
                        TVector3 rotationAxis = zAxis.Cross(v0Momentum).Unit();
                        double rotationAngle = acos(zAxis.Dot(v0Momentum));
                        TRotation rotation;
                        rotation.Rotate(-rotationAngle, rotationAxis);

                        TVector3 p1MomentumRotated = rotation * p1Vec;
                        TVector3 p2MomentumRotated = rotation * p2Vec;

                        float pL1 = p1MomentumRotated.Z();
                        float pL2 = p2MomentumRotated.Z();

                        // Calculate transverse momentum components in the rotated system
                        float pT1 = p1MomentumRotated.Perp();
                        float pT2 = p2MomentumRotated.Perp();

                        // float alpha = (pL1 - pL2) / (pL1 + pL2);
                        // float pT = (pT1 + pT2) / 2;
                        // hArmenteros->Fill(alpha, pT);
                        
                        // Below is the code limiting the angle theta*

                        // Calculate cos(θ*) for each daughter particle in the V0 rest frame
                        TVector3 p1VecRest = p1Vec - v0Momentum * (v0Momentum.Dot(p1Vec) / v0Momentum.Mag2());
                        TVector3 p2VecRest = p2Vec - v0Momentum * (v0Momentum.Dot(p2Vec) / v0Momentum.Mag2());

                        float cosThetaStar1 = cos(p1VecRest.Angle(v0Momentum));
                        float cosThetaStar2 = cos(p2VecRest.Angle(v0Momentum));

                        if (cosThetaStar1 > -0.95 && cosThetaStar1 < 0.8 && cosThetaStar2 > -0.95 && cosThetaStar2 < 0.8) {
                            float alpha = (pL1 - pL2) / (pL1 + pL2);
                            float pT = (pT1 + pT2) / 2;

                            hArmenteros->Fill(alpha, pT);
                        }
                    }
                }
            }
        }
    }

    
    // TFile outputFile(argv[2], "RECREATE");
    TFile outputFile(argv[2], "UPDATE");
    hArmenteros->Write();
    outputFile.Close();
    inputFile.Close();

    return 0;
}