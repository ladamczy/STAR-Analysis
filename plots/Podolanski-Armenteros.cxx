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
    if (argc != 4) {
        cerr << "Usage: " << argv[0] << " [input file K0's data] [input file Lambda's data] [output file]" << endl;
        return 1;
    }

    TFile inputFile1(argv[1], "READ");
    if (!inputFile1.IsOpen()) {
        cerr << "Failed to open input file 1." << endl;
        return 1;
    }

    TFile inputFile2(argv[2], "READ");
    if (!inputFile2.IsOpen()) {
        cerr << "Failed to open input file 2." << endl;
        return 1;
    }

    TTree* tree1 = static_cast<TTree*>(inputFile1.Get("mUPCTree"));
    if (!tree1) {
        cerr << "Failed to retrieve mUPCTree from the input file 1." << endl;
        inputFile1.Close();
        return 1;
    }

    TTree* tree2 = static_cast<TTree*>(inputFile2.Get("mUPCTree"));
    if (!tree2) {
        cerr << "Failed to retrieve mUPCTree from the input file 2." << endl;
        inputFile2.Close();
        return 1;
    }

    StUPCEvent* event1 = nullptr;
    StUPCEvent* event2 = nullptr;
    tree1->SetBranchAddress("mUPCEvent", &event1);
    tree2->SetBranchAddress("mUPCEvent", &event2);

    const float pionMass = 0.13957;
    const float protonMass = 0.938;

    double lowerLimitOfInvMassLambda = 1.08;
    double upperLimitOfInvMassLambda = 1.15;
    double lowerLimitOfInvMassK0 = 0.46;
    double upperLimitOfInvMassK0 = 0.53;

    TH2D* hArmenteros = new TH2D("hArmenteros", "Armenteros Plot;alpha;pT (GeV/c)", 200, -1, 1, 100, 0, 0.3);

    // petla po pliku danych K0s
    Long64_t nEntries1 = tree1->GetEntries();
    for (Long64_t i = 0; i < nEntries1; ++i) {
        tree1->GetEntry(i);

        for (int j = 0; j < event1->getNumberOfTracks(); ++j) {
            StUPCTrack* track1 = event1->getTrack(j);
            for (int k = j + 1; k < event1->getNumberOfTracks(); ++k) {
                StUPCTrack* track2 = event1->getTrack(k);
                // K0s
                if ((track1->getCharge() > 0 && track2->getCharge() < 0) &&
                    track1->getNhits() > 15 && track2->getNhits() > 15 &&
                    track1->getPt() > 0.15 && track2->getPt() > 0.15 &&
                    abs(track1->getEta()) < 1.1 && abs(track2->getEta()) < 1.1 &&
                    (track1->getFlag(StUPCTrack::kTof) || track2->getFlag(StUPCTrack::kTof)) &&
                    abs(track1->getNSigmasTPCPion()) < 2.5 && abs(track2->getNSigmasTPCPion()) < 2.5) {

                    // Tworzenie obiektu StUPCV0
                    TVector3 const tryVec(0,0,0);
                    double beamLine[] = {0,0,0,0};
                    StUPCV0 v0(track1, track2, pionMass, pionMass, 1, 1, tryVec, beamLine, event1->getMagneticField(), true);


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

                        float alpha = (pL1 - pL2) / (pL1 + pL2);
                        float pT = (pT1 + pT2) / 2;
                        hArmenteros->Fill(alpha, pT);
                    }
                }
            }
        }
    }


    Long64_t nEntries2 = tree2->GetEntries();
    for (Long64_t i = 0; i < nEntries2; ++i) {
        tree2->GetEntry(i);

        for (int j = 0; j < event2->getNumberOfTracks(); ++j) {
            StUPCTrack* track1 = event2->getTrack(j);
            for (int k = j + 1; k < event2->getNumberOfTracks(); ++k) {
                StUPCTrack* track2 = event2->getTrack(k);
                // Lambda
                if ((track1->getCharge() > 0 && track2->getCharge() < 0) &&
                    track1->getNhits() > 15 && track2->getNhits() > 15 &&
                    track1->getPt() > 0.15 && track2->getPt() > 0.15 &&
                    abs(track1->getEta()) < 1.1 && abs(track2->getEta()) < 1.1 &&
                    (track1->getFlag(StUPCTrack::kTof) || track2->getFlag(StUPCTrack::kTof)) &&
                    abs(track1->getNSigmasTPCProton()) < 2.5 && abs(track2->getNSigmasTPCPion()) < 2.5) {

                    // Tworzenie obiektu StUPCV0
                    TVector3 const tryVec(0,0,0);
                    double beamLine[] = {0,0,0,0};
                    StUPCV0 v0(track1, track2, protonMass, pionMass, 1, 1, tryVec, beamLine, event2->getMagneticField(), true);

                    if (v0.m() > lowerLimitOfInvMassLambda && v0.m() < upperLimitOfInvMassLambda) {
                        TVector3 p1Vec, p2Vec;
                        track1->getMomentum(p1Vec);
                        track2->getMomentum(p2Vec);

                        TVector3 V0_ped = TVector3(v0.px(), v0.py(), v0.pz());
                        TVector3 v0Momentum = V0_ped.Unit();
                        TVector3 zAxis(0, 0, 1);
                        TVector3 rotationAxis = zAxis.Cross(v0Momentum).Unit();
                        double rotationAngle = acos(zAxis.Dot(v0Momentum));
                        TRotation rotation;
                        rotation.Rotate(-rotationAngle, rotationAxis);

                        TVector3 p1MomentumRotated = rotation * p1Vec;
                        TVector3 p2MomentumRotated = rotation * p2Vec;
                        // TVector3 ptMomentumRotated = rotation * V0_ped;

                        float pL1 = p1MomentumRotated.Z();
                        float pL2 = p2MomentumRotated.Z();

                        float pT1 = p1MomentumRotated.Perp();
                        float pT2 = p2MomentumRotated.Perp();

                        float alpha = (pL1 - pL2) / (pL1 + pL2);
                        float pT = (pT1 + pT2) / 2;
                        hArmenteros->Fill(alpha, pT);
                    }
                }
                // LambdaBar
                if ((track1->getCharge() < 0 && track2->getCharge() > 0) &&
                    track1->getNhits() > 15 && track2->getNhits() > 15 &&
                    track1->getPt() > 0.15 && track2->getPt() > 0.15 &&
                    abs(track1->getEta()) < 1.1 && abs(track2->getEta()) < 1.1 &&
                    (track1->getFlag(StUPCTrack::kTof) || track2->getFlag(StUPCTrack::kTof)) &&
                    abs(track1->getNSigmasTPCProton()) < 2.5 && abs(track2->getNSigmasTPCPion()) < 2.5) {

                    // Tworzenie obiektu StUPCV0
                    TVector3 const tryVec(0,0,0);
                    double beamLine[] = {0,0,0,0};
                    StUPCV0 v0(track1, track2, protonMass, pionMass, 1, 1, tryVec, beamLine, event2->getMagneticField(), true);


                    if (v0.m() > lowerLimitOfInvMassLambda && v0.m() < upperLimitOfInvMassLambda) {
                        TVector3 p1Vec, p2Vec;
                        track1->getMomentum(p1Vec);
                        track2->getMomentum(p2Vec);

                        TVector3 V0_ped = TVector3(v0.px(), v0.py(), v0.pz());
                        TVector3 v0Momentum = V0_ped.Unit();
                        TVector3 zAxis(0, 0, 1);
                        TVector3 rotationAxis = zAxis.Cross(v0Momentum).Unit();
                        double rotationAngle = acos(zAxis.Dot(v0Momentum));
                        TRotation rotation;
                        rotation.Rotate(-rotationAngle, rotationAxis);

                        TVector3 p1MomentumRotated = rotation * p1Vec;
                        TVector3 p2MomentumRotated = rotation * p2Vec;
                        // TVector3 ptMomentumRotated = rotation * V0_ped;

                        float pL1 = p1MomentumRotated.Z();
                        float pL2 = p2MomentumRotated.Z();

                        float pT1 = p1MomentumRotated.Perp();
                        float pT2 = p2MomentumRotated.Perp();

                        float alpha = (-pL1 + pL2) / (pL1 + pL2);
                        float pT = (pT1 + pT2) / 2;
                        hArmenteros->Fill(alpha, pT);
                    }
                }
            }
        }
    }

    
    TFile outputFile(argv[3], "RECREATE");
    hArmenteros->Write();
    outputFile.Close();
    inputFile1.Close();
    inputFile2.Close();

    return 0;
}