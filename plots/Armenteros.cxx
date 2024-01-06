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
                if ((track1->getCharge() > 0 && track2->getCharge() < 0) && /* pozostałe warunki */) {
                    StUPCV0 v0(track1, track2, protonMass, pionMass, /* pozostałe parametry */);
                    if (v0.m() > lowerLimitOfInvMassLambda && v0.m() < upperLimitOfInvMassLambda) {
                        TVector3 p1, p2;
                        track1->getMomentum(p1); // Proton
                        track2->getMomentum(p2); // Pion

                        // Przeliczanie na układ środka masy (CMS) - przykładowe, wymaga dostosowania
                        TVector3 pCMS = p1 + p2;
                        TVector3 beta = -pCMS.BoostVector();
                        TLorentzVector p1Lorentz(p1, track1->getEnergy(protonMass));
                        TLorentzVector p2Lorentz(p2, track2->getEnergy(pionMass));
                        p1Lorentz.Boost(beta);
                        p2Lorentz.Boost(beta);

                        float alpha = (p1Lorentz.Pz() - p2Lorentz.Pz()) / (p1Lorentz.Pz() + p2Lorentz.Pz());
                        hArmenteros->Fill(alpha, v0.pt());
                    }
                }

                // Selekcja dla LambdaBar
                if ((track1->getCharge() < 0 && track2->getCharge() > 0) && /* pozostałe warunki */) {
                    StUPCV0 v0(track1, track2, protonMass, pionMass, /* pozostałe parametry */);
                    if (v0.m() > lowerLimitOfInvMassLambda && v0.m() < upperLimitOfInvMassLambda) {
                        TVector3 p1, p2;
                        track1->getMomentum(p1); // Antyproton
                        track2->getMomentum(p2); // Pion

                        // Przeliczanie na układ środka masy (CMS) - przykładowe, wymaga dostosowania
                        TVector3 pCMS = p1 + p2;
                        TVector3 beta = -pCMS.BoostVector();
                        TLorentzVector p1Lorentz(p1, track1->getEnergy(protonMass)); // Używamy masy protonu
                        TLorentzVector p2Lorentz(p2, track2->getEnergy(pionMass));
                        p1Lorentz.Boost(beta);
                        p2Lorentz.Boost(beta);

                        float alpha = (p2Lorentz.Pz() - p1Lorentz.Pz()) / (p2Lorentz.Pz() + p1Lorentz.Pz());
                        hArmenteros->Fill(alpha, v0.pt());
                    }
                }

                // Selekcja dla K0
                if (track1->getCharge() != track2->getCharge() && /* pozostałe warunki */) {
                    StUPCV0 v0(track1, track2, pionMass, pionMass, /* pozostałe parametry */);
                    if (v0.m() > lowerLimitOfInvMassK0 && v0.m() < upperLimitOfInvMassK0) {
                        TVector3 p1, p2;
                        track1->getMomentum(p1); // Pion
                        track2->getMomentum(p2); // Pion

                        // Przeliczanie na układ środka masy (CMS) - przykładowe, wymaga dostosowania
                        TVector3 pCMS = p1 + p2;
                        TVector3 beta = -pCMS.BoostVector();
                        TLorentzVector p1Lorentz(p1, track1->getEnergy(pionMass));
                        TLorentzVector p2Lorentz(p2, track2->getEnergy(pionMass));
                        p1Lorentz.Boost(beta);
                        p2Lorentz.Boost(beta);

                        float alpha = (p1Lorentz.Pz() - p2Lorentz.Pz()) / (p1Lorentz.Pz() + p2Lorentz.Pz());
                        hArmenteros->Fill(alpha, v0.pt());
                    }
                }
            }
        }
    }

    TFile outputFile(argv[2], "RECREATE");
    hArmenteros->Write();
    outputFile.Close();

    inputFile.Close();

    return 0;
}