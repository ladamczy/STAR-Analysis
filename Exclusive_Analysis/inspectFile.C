#include <TFile.h>
#include <TTree.h>
#include <iostream>

void inspectFile(const char* filename) {
    TFile* f = TFile::Open(filename);
    if (!f || f->IsZombie()) {
        std::cerr << "Error opening file " << filename << std::endl;
        return;
    }

    // List all keys in the file
    std::cout << "\nKeys in " << filename << ":\n";
    f->ls();

    // Get the tree (assuming it's named "mUPCTree")
    TTree* tree = dynamic_cast<TTree*>(f->Get("mUPCTree"));
    if (!tree) {
        std::cerr << "No tree named mUPCTree found.\n";
        f->Close();
        return;
    }

    // Print tree structure (branches)
    std::cout << "\nBranches in mUPCTree:\n";
    tree->Print();

    // Optional: read first event and dump some info
    // (requires linking against the STAR class libraries)
    // If you have the libraries loaded, you can do:

    /*
    StUPCEvent* upcEvt = nullptr;
    StRPEvent* rpEvt = nullptr;
    tree->SetBranchAddress("mUPCEvent", &upcEvt);
    tree->SetBranchAddress("mRPEvent", &rpEvt);

    tree->GetEntry(0);
    if (upcEvt) {
        std::cout << "Run number: " << upcEvt->getRunNumber() << std::endl;
        std::cout << "Number of tracks: " << upcEvt->getNumberOfTracks() << std::endl;
        // Check flags of first track
        if (upcEvt->getNumberOfTracks() > 0) {
            StUPCTrack* trk = upcEvt->getTrack(0);
            std::cout << "Track 0 flags: "
                      << " kV0=" << trk->getFlag(StUPCTrack::kV0)
                      << " kPrimary=" << trk->getFlag(StUPCTrack::kPrimary)
                      << " kCEP=" << trk->getFlag(StUPCTrack::kCEP)
                      << std::endl;
        }
    }
    if (rpEvt) {
        std::cout << "Number of RP tracks: " << rpEvt->getNumberOfTracks() << std::endl;
        // Optionally print first RP track theta to see if corrections are applied
        if (rpEvt->getNumberOfTracks() > 0) {
            StUPCRpsTrack* rptrk = rpEvt->getTrack(0);
            std::cout << "RP track 0 thetaX: " << rptrk->thetaRp(0) << std::endl;
        }
    }
    */

    f->Close();
}
