#include <TFile.h>
#include <TTree.h>

void MergeTrees() {
    const char* inputFile = "treepart2.root";
    const char* outputFile = "mergedTree.root";
    
    TFile* file = TFile::Open(inputFile);
    
    TTree* T1 = static_cast<TTree*>(file->Get("T1"));
    TTree* T2 = static_cast<TTree*>(file->Get("T2"));
    
    TFile* outputFilePtr = TFile::Open(outputFile, "recreate");
    TTree* T3 = new TTree("T3", "");
    
    Int_t Run, Event;
    Float_t x, y;
    
    T3->Branch("Run", &Run, "Run/I");
    T3->Branch("Event", &Event, "Event/I");
    T3->Branch("x", &x, "x/F");
    T3->Branch("y", &y, "y/F");
    
    T1->SetBranchAddress("Run", &Run);
    T1->SetBranchAddress("Event", &Event);
    T1->SetBranchAddress("x", &x);
    
    T2->SetBranchAddress("Run", &Run);
    T2->SetBranchAddress("Event", &Event);
    T2->SetBranchAddress("y", &y);
    
    for (Long64_t i = 0; i < T1->GetEntries(); ++i) {
        T1->GetEntry(i);
        T3->Fill();
    }
    
    for (Long64_t i = 0; i < T2->GetEntries(); ++i) {
        T2->GetEntry(i);
        T3->Fill();
    }
    
    T3->Write();
    outputFilePtr->Close();
    file->Close();
    
    delete outputFilePtr;
    delete file;
}