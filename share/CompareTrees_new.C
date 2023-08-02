#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <iostream>
#include <map>
#include <cmath>

void CompareTrees_new() {
    TFile* file1 = TFile::Open("treeparrent.root");
    if (!file1 || file1->IsZombie()) {
        std::cerr << "Nie można otworzyć pliku treeparrent.root." << std::endl;
        return;
    }

    TFile* file2 = TFile::Open("mergedTree.root");
    if (!file2 || file2->IsZombie()) {
        std::cerr << "Nie można otworzyć pliku mergedTree.root." << std::endl;
        return;
    }

    TTree* treeT = static_cast<TTree*>(file1->Get("T"));
    if (!treeT) {
        std::cerr << "Nie znaleziono drzewa o nazwie T w pliku treeparrent.root." << std::endl;
        return;
    }

    TTree* treeT3 = static_cast<TTree*>(file2->Get("T3"));
    if (!treeT3) {
        std::cerr << "Nie znaleziono drzewa o nazwie T3 w pliku mergedTree.root." << std::endl;
        return;
    }

    Long64_t numEntriesT = treeT->GetEntries();
    Long64_t numEntriesT3 = treeT3->GetEntries();

    TObjArray* branchesT = treeT->GetListOfBranches();
    int numBranchesT = branchesT->GetEntries();

    TObjArray* branchesT3 = treeT3->GetListOfBranches();
    int numBranchesT3 = branchesT3->GetEntries();

    if (numBranchesT != numBranchesT3) {
        std::cout << "Liczba gałęzi w drzewach T i T3 jest różna." << std::endl;
        return;
    }

    std::map<std::pair<Float_t, Float_t>, std::pair<Float_t, Float_t>> leafMapT;

    for (int i = 0; i < numEntriesT; ++i) {
        treeT->GetEntry(i);
        Float_t x = *static_cast<Float_t*>(treeT->GetLeaf("x")->GetValuePointer());
        Float_t y = *static_cast<Float_t*>(treeT->GetLeaf("y")->GetValuePointer());
        Float_t run = *static_cast<Float_t*>(treeT->GetLeaf("Run")->GetValuePointer());
        Float_t event = *static_cast<Float_t*>(treeT->GetLeaf("Event")->GetValuePointer());
        leafMapT[std::make_pair(run, event)] = std::make_pair(x, y);
    }

    bool areEqual = true;
    for (int i = 0; i < numEntriesT3; ++i) {
        treeT3->GetEntry(i);
        Float_t x = *static_cast<Float_t*>(treeT3->GetLeaf("x")->GetValuePointer());
        Float_t y = *static_cast<Float_t*>(treeT3->GetLeaf("y")->GetValuePointer());
        Float_t run = *static_cast<Float_t*>(treeT3->GetLeaf("Run")->GetValuePointer());
        Float_t event = *static_cast<Float_t*>(treeT3->GetLeaf("Event")->GetValuePointer());

        auto it = leafMapT.find(std::make_pair(run, event));
        if (it != leafMapT.end()) {
            Float_t xT = it->second.first;
            Float_t yT = it->second.second;
            if (std::abs(x - xT) > std::numeric_limits<Float_t>::epsilon() ||
                std::abs(y - yT) > std::numeric_limits<Float_t>::epsilon()) {
                std::cout << "Różnica dla przypadku Run: " << run << ", Event: " << event << std::endl;
                areEqual = false;
            }
        } else {
            std::cout << "Nie znaleziono odpowiadającego przypadku w drzewie T dla Run: " << run << ", Event: " << event << std::endl;
            areEqual = false;
        }
    }

    if (areEqual) {
        std::cout << "Dane w drzewach T i T3 są identyczne dla wszystkich przypadków Run, Event." << std::endl;
    }

    file1->Close();
    file2->Close();

    delete file1;
    delete file2;
}
