#include <TFile.h>
#include <TTree.h>
#include <map>

void MergeTrees_new() {
    const char* inputFile1 = "treepart1.root";
    const char* inputFile2 = "treepart2.root";
    const char* outputFile = "mergedTree.root";
    Int_t skipped = 0;

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

    Int_t Run1, Event1;
    Float_t x1;
    T1->SetBranchAddress("Run", &Run1);
    T1->SetBranchAddress("Event", &Event1);
    T1->SetBranchAddress("x", &x1);

    Int_t Run2, Event2;
    Float_t y2;
    T2->SetBranchAddress("Run", &Run2);
    T2->SetBranchAddress("Event", &Event2);
    T2->SetBranchAddress("y", &y2);

    std::map<std::pair<Int_t, Int_t>, Float_t> yMap;
    Long64_t entries2 = T2->GetEntries();
    for (Long64_t i = 0; i < entries2; ++i) {
        T2->GetEntry(i);
        yMap[std::make_pair(Run2, Event2)] = y2;
    }

    TFile* outputFilePtr = TFile::Open(outputFile, "recreate");
    TTree* T3 = new TTree("T3", "");

    T3->Branch("Run", &Run1, "Run/I");
    T3->Branch("Event", &Event1, "Event/I");
    T3->Branch("x", &x1, "x/F");
    T3->Branch("y", &y2, "y/F");

    Long64_t entries1 = T1->GetEntries();
    for (Long64_t i = 0; i < entries1; ++i) {
        T1->GetEntry(i);
        auto it = yMap.find(std::make_pair(Run1, Event1));
        if (it != yMap.end()) {
            y2 = it->second;
            T3->Fill();
        }
        else {
        	++skipped;
        }
    }

    T3->Write();
    outputFilePtr->Close();
    file1->Close();
    file2->Close();

    delete outputFilePtr;
    delete file1;
    delete file2;
    
    std::cout << "Polaczono drzewa" << inputFile1 << "oraz" << inputFile2 << std::endl;
    std::cout << "Pominiete przypadki: " << skipped << std::endl;
}

