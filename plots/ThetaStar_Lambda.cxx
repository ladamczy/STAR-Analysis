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

    TH1D* hThetaStar = new TH1D("#theta Angle in CMS", ";#theta*;Counts", 100, -TMath::Pi(), TMath::Pi());
    
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);

        for (int j = 0; j < event->getNumberOfTracks(); ++j) {
            StUPCTrack* track1 = event->getTrack(j);
            for (int k = j + 1; k < event->getNumberOfTracks(); ++k) {
                StUPCTrack* track2 = event->getTrack(k);

                if ((track1->getCharge() < 0 && track2->getCharge() > 0) &&
                    track1->getNhits() > 15 && track2->getNhits() > 15 &&
                    track1->getPt() > 0.15 && track2->getPt() > 0.15 &&
                    abs(track1->getEta()) < 1.1 && abs(track2->getEta()) < 1.1 &&
                    (track1->getFlag(StUPCTrack::kTof) || track2->getFlag(StUPCTrack::kTof))) {

                    // Tworzenie obiektu StUPCV0
                    TVector3 const tryVec(0,0,0);
                    double beamLine[] = {0,0,0,0};
                    StUPCV0 v0(track1, track2, protonMass, pionMass, 1, 1, tryVec, beamLine, event->getMagneticField(), true);

                    // podejscie z TLorentzVector
                    if (v0.m() > lowerLimitOfInvMassLambda && v0.m() < upperLimitOfInvMassLambda) {
                        TLorentzVector p1, p2, v0;
                        TVector3 p1Vec, p2Vec;

                        // Get momenta from tracks
                        track1->getMomentum(p1Vec);
                        track2->getMomentum(p2Vec);

                        // Set TLorentzVector for each particle (assuming you have mass for each particle)
                        p1.SetVectM(p1Vec, protonMass); // mass1 is the mass of the first particle
                        p2.SetVectM(p2Vec, pionMass); // mass2 is the mass of the second particle

                        // Calculate the V0 Lorentz vector
                        v0 = p1 + p2;

                        // Boost to the V0 rest frame
                        TVector3 boostVector = -v0.BoostVector();
                        p1.Boost(boostVector);
                        p2.Boost(boostVector);

                        // Calculate cos(θ*) for each daughter particle in the V0 rest frame
                        float cosThetaStar1 = cos(p1.Vect().Angle(v0.Vect()));
                        float cosThetaStar2 = cos(p2.Vect().Angle(v0.Vect()));

                        // Calculate θ* angle for each daughter particle
                        float thetaStar1 = acos(cosThetaStar1);
                        float thetaStar2 = acos(cosThetaStar2);

                        // Fill the histogram with θ* angles for both daughters
                        hThetaStar->Fill(thetaStar1);
                        hThetaStar->Fill(thetaStar2);
                    }

                }
            }
        }
    }

    
    // TFile outputFile(argv[2], "RECREATE");
    TFile outputFile(argv[2], "UPDATE");
    hThetaStar->Write();
    outputFile.Close();
    inputFile.Close();

    return 0;
}