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

    // Histograms
    TH1D* histPointingAngleHypo = new TH1D("histPointingAngleHypo", "Pointing Angle;Angle;Counts", 100, 0, TMath::Pi());
    TH1D* histDcaBeamLine = new TH1D("histDcaBeamLine", "DCA from Beam Line;DCA (cm);Counts", 100, 0, 5);
    TH1D* histInvariantMass = new TH1D("histInvariantMass", "Invariant Mass;Mass (GeV/c^2);Counts", 100, 0, 5);
    TH1D* histTransverseDecayVertex = new TH1D("histTransverseDecayVertex", "Transverse Decay Vertex;Distance (cm);Counts", 100, 0, 10);
    TH1D* histDecayLength = new TH1D("histDecayLength", "Decay Length;Length (cm);Counts", 100, 0, 10);
    TH1D* histDcaDaughters = new TH1D("histDcaDaughters", "DCA between Daughters;DCA (cm);Counts", 100, 0, 5);
    TH1D* histPt = new TH1D("histPt", "Transverse Momentum;Pt (GeV/c);Counts", 100, 0, 10);
    TH1D* histEta = new TH1D("histEta", "Pseudorapidity;Eta;Counts", 100, -5, 5);
    TH1D* histProdVertexHypo = new TH1D("histProdVertexHypo", "Production Vertex;Distance (cm);Counts", 100, 0, 10);
    TH1D* histParticle1Dca = new TH1D("histParticle1Dca", "DCA of Particle 1;DCA (cm);Counts", 100, 0, 5);
    TH1D* histParticle2Dca = new TH1D("histParticle2Dca", "DCA of Particle 2;DCA (cm);Counts", 100, 0, 5);    


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
                StUPCV0 v0(track1, track2, pionMass, pionMass, 1, 1, tryVec, beamLine, event->getMagneticField(), true);

                // Wypełnianie histogramów z odpowiednimi cięciami
                if (v0.pointingAngleHypo() > 0.925) {
                    histPointingAngleHypo->Fill(v0.pointingAngleHypo());
                }
                if (v0.DCABeamLine() < 1.5) {
                    histDcaBeamLine->Fill(v0.DCABeamLine());
                }
                // invariant mass of K0 candidate 0.495 +- 0.035 
                // invariant mass of Lambda0 candidate 1.115 +- 0.035
                // invariant mass of Lambda0Bar candidate 1.115 +- 0.035 
                if (v0.m() > 0.495 - 0.035 && v0.m() < 0.495 + 0.035) {
                    histInvariantMass->Fill(v0.m());
                }
                // tu nie wiem jaka wartość cięcia
                if (v0.decayVertex().Perp()) {
                    histTransverseDecayVertex->Fill(v0.decayVertex().Perp());
                }
                // decay length of K0 candidate: [1] 2.7 cm or [2] 2.0 cm ? czy zakres 1.9 - 2.8 jest ok?
                // decay length of Lambda0 candidate: [1] 7.9 cm or [2] 7.0 cm ? czy zakres 6.9 - 8.0 jest ok?
                // decay length of Lambda0Bar candidate ?
                // czy moze powinno decayLengthHypo zamias decayLength tak jak w przypadku kata?
                if (v0.decayLength()) {
                    histDecayLength->Fill(v0.decayLength());
                }
                if (v0.dcaDaughters() < 1.5) {
                    histDcaDaughters->Fill(v0.dcaDaughters());
                }
                if (v0.pt() > 0.15) {
                    histPt->Fill(v0.pt());
                }
                if (TMath::Abs(v0.eta()) < 1.1) {
                    histEta->Fill(v0.eta());
                }
                // tu poniżej również nie znam wartości cięcia
                if (v0.prodVertexHypo()) {
                    histProdVertexHypo->Fill(v0.prodVertexHypo());
                }
                if (v0.particle1Dca()) {
                    histParticle1Dca->Fill(v0.particle1Dca());
                }
                if (v0.particle2Dca()) {
                    histParticle2Dca->Fill(v0.particle2Dca());
                }
            }
        }
    }

    TFile outputFile(argv[2], "RECREATE");
    histPointingAngleHypo->Write();
    histDecayLength->Write();
    histDcaBeamLine->Write();
    histInvariantMass->Write();
    histTransverseDecayVertex->Write();
    histDcaDaughters->Write();
    histPt->Write();
    histEta->Write();
    histProdVertexHypo->Write();
    histParticle1Dca->Write();
    histParticle2Dca->Write();

    outputFile.Close();
    inputFile.Close();

    return 0;
}
