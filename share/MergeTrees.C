#include <TFile.h>
#include <TTree.h>

void MergeTrees() {
    const char* inputFile1 = "treepart1.root";
    const char* inputFile2 = "treepart2.root";
    const char* outputFile = "mergedTree.root";

    TFile* file1 = TFile::Open(inputFile1);
    TFile* file2 = TFile::Open(inputFile2);

    TTree* T1 = static_cast<TTree*>(file1->Get("T1"));
    TTree* T2 = static_cast<TTree*>(file2->Get("T2"));

    if (!T1 || !T2) {
        std::cerr << "Nie można odczytać drzewa T1 lub T2 z pliku." << std::endl;
        file1->Close();
        file2->Close();
        return;
    }

    Int_t Run, Event;
    Float_t x, y;

    T1->SetBranchAddress("Run", &Run);
    T1->SetBranchAddress("Event", &Event);
    T1->SetBranchAddress("x", &x);

    T2->SetBranchAddress("Run", &Run);
    T2->SetBranchAddress("Event", &Event);
    T2->SetBranchAddress("y", &y);

    TFile* outputFilePtr = TFile::Open(outputFile, "recreate");
    TTree* T3 = new TTree("T3", "");

    T3->Branch("Run", &Run, "Run/I");
    T3->Branch("Event", &Event, "Event/I");
    T3->Branch("x", &x, "x/F");
    T3->Branch("y", &y, "y/F");

    Long64_t entries1 = T1->GetEntries();
    Long64_t entries2 = T2->GetEntries();
    Long64_t maxEntries = (entries1 < entries2) ? entries1 : entries2;

    for (Long64_t i = 0; i < maxEntries; ++i) {
        T1->GetEntry(i);
        T2->GetEntry(i);
        T3->Fill();
    }

    T3->Write();
    outputFilePtr->Close();
    file1->Close();
    file2->Close();

    delete outputFilePtr;
    delete file1;
    delete file2;
}
