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

    TH2D* hArmenteros = new TH2D("hArmenteros", ";#alpha;pT (GeV/c)", 200, -1, 1, 100, 0, 0.3);

    Long64_t nEntries1 = tree1->GetEntries();
    for (Long64_t i = 0; i < nEntries1; ++i) {
        tree1->GetEntry(i);

        for (int j = 0; j < event1->getNumberOfTracks(); ++j) {
            StUPCTrack* track1 = event1->getTrack(j);
            for (int k = j + 1; k < event1->getNumberOfTracks(); ++k) {
                StUPCTrack* track2 = event1->getTrack(k);
                // K0s
                if ((track1->getCharge() != track2->getCharge()) &&
                    track1->getNhits() > 15 && track2->getNhits() > 15 &&
                    track1->getPt() > 0.15 && track2->getPt() > 0.15 &&
                    abs(track1->getEta()) < 1.1 && abs(track2->getEta()) < 1.1 &&
                    (track1->getFlag(StUPCTrack::kTof) && track2->getFlag(StUPCTrack::kTof)) &&
                    abs(track1->getNSigmasTPCPion()) < 2.5 && abs(track2->getNSigmasTPCPion()) < 2.5) {

                    TVector3 const tryVec(0,0,0);
                    double beamLine[] = {0,0,0,0};
                    StUPCV0 v0(track1, track2, pionMass, pionMass, 1, 1, tryVec, beamLine, event1->getMagneticField(), true);


                    if (v0.m() > lowerLimitOfInvMassK0 && v0.m() < upperLimitOfInvMassK0 && v0.decayLengthHypo() > 3.0) {
                        float alpha = v0.alphaAP();
                        float pT = v0.ptAP();
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
                    (track1->getFlag(StUPCTrack::kTof) && track2->getFlag(StUPCTrack::kTof)) &&
                    abs(track1->getNSigmasTPCProton()) < 2.5 && abs(track2->getNSigmasTPCPion()) < 2.5) {

                    TVector3 const tryVec(0,0,0);
                    double beamLine[] = {0,0,0,0};
                    StUPCV0 v0(track1, track2, protonMass, pionMass, 1, 1, tryVec, beamLine, event2->getMagneticField(), true);

                    if (v0.m() > lowerLimitOfInvMassLambda && v0.m() < upperLimitOfInvMassLambda && v0.decayLengthHypo() > 3.0) {
                        float alpha = v0.alphaAP();
                        float pT = v0.ptAP();
                        hArmenteros->Fill(alpha, pT);
                    }
                }
                // LambdaBar
                if ((track1->getCharge() < 0 && track2->getCharge() > 0) &&
                    track1->getNhits() > 15 && track2->getNhits() > 15 &&
                    track1->getPt() > 0.15 && track2->getPt() > 0.15 &&
                    abs(track1->getEta()) < 1.1 && abs(track2->getEta()) < 1.1 &&
                    (track1->getFlag(StUPCTrack::kTof) && track2->getFlag(StUPCTrack::kTof)) &&
                    abs(track1->getNSigmasTPCProton()) < 2.5 && abs(track2->getNSigmasTPCPion()) < 2.5) {

                    TVector3 const tryVec(0,0,0);
                    double beamLine[] = {0,0,0,0};
                    StUPCV0 v0(track1, track2, protonMass, pionMass, 1, 1, tryVec, beamLine, event2->getMagneticField(), true);


                    if (v0.m() > lowerLimitOfInvMassLambda && v0.m() < upperLimitOfInvMassLambda && v0.decayLengthHypo() > 3.0) {
                        float alpha = v0.alphaAP();
                        float pT = v0.ptAP();
                        hArmenteros->Fill(alpha, pT);
                    }
                }
            }
        }
    }

    TCanvas* canvas = new TCanvas("canvas", "Armenteros Plot Canvas", 2560, 1440);
    hArmenteros->SetStats(0);
    hArmenteros->Draw("colz");
    canvas->SetLogz();
    canvas->SaveAs("ARMENTEROS_BEZ.jpg");


    TFile outputFile(argv[3], "RECREATE");
    hArmenteros->Write();
    outputFile.Close();
    inputFile1.Close();
    inputFile2.Close();

    return 0;
}