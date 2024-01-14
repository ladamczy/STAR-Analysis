#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TMath.h>
#include "StUPCEvent.h"
#include "StUPCTrack.h"
#include "StUPCV0.h"
#include <TLorentzVector.h>
#include <TVector3.h>

using namespace std;

void RotateToV0RestFrame(TLorentzVector& p, const TLorentzVector& v0) {
    // Określenie kierunku oś-Z w układzie V0
    TVector3 zAxis = v0.Vect().Unit();

    // Określenie kierunku oś-X w układzie V0
    TVector3 xAxis = v0.Vect().Cross(TVector3(0, 0, 1)).Unit();

    // Określenie kierunku oś-Y w układzie V0
    TVector3 yAxis = zAxis.Cross(xAxis).Unit();

    // Określenie kąta obrotu do płaszczyzny yz
    double theta_yz = TMath::ATan2(xAxis.Y(), xAxis.X());

    // Obrót do płaszczyzny yz
    p.RotateUz(zAxis);
    p.Rotate(theta_yz, TVector3(0, 1, 0));
}

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

    TH1D* hThetaStar = new TH1D("#theta Angle in CMS", ";cos(#theta*);Counts", 100, -TMath::Pi(), TMath::Pi());
    
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

                    TVector3 const tryVec(0,0,0);
                    double beamLine[] = {0,0,0,0};
                    StUPCV0 v0(track1, track2, pionMass, pionMass, 1, 1, tryVec, beamLine, event->getMagneticField(), true);

                    if (v0.m() > lowerLimitOfInvMassK0 && v0.m() < upperLimitOfInvMassK0) {
                        TLorentzVector p1, p2, v0;
                        TVector3 p1Vec, p2Vec, L;

                        track1->getMomentum(p1Vec);
                        track2->getMomentum(p2Vec);

                        p1.SetVectM(p1Vec, pionMass);
                        p2.SetVectM(p2Vec, pionMass);

                        v0 = p1 + p2;

                        TVector3 boostVector = -v0.BoostVector();
                        // p1.RotateUz(boostVector);
                        // p2.RotateUz(boostVector);
                        TVector3 zAxis = v0.Vect().Unit();

                        p1.RotateUz(zAxis);
                        p2.RotateUz(zAxis);
                        p1.Boost(boostVector);
                        p2.Boost(boostVector);;

                        float cosThetaStar1 = p1.Vect().Z() / p1.Vect().Mag();
                        hThetaStar->Fill(cosThetaStar1);

                    }

                }
            }
        }
    }

    
    TFile outputFile(argv[2], "RECREATE");
    // TFile outputFile(argv[2], "UPDATE");
    hThetaStar->Write();
    outputFile.Close();
    inputFile.Close();

    return 0;
}