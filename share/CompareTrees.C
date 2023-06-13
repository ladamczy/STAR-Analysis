#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <iostream>

void CompareTrees() {
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

    // Pobieranie liczby wpisów w drzewach
    Long64_t numEntriesT = treeT->GetEntries();
    Long64_t numEntriesT3 = treeT3->GetEntries();

    // Pobieranie liczby gałęzi w drzewach
    TObjArray* branchesT = treeT->GetListOfBranches();
    int numBranchesT = branchesT->GetEntries();
    TObjArray* branchesT3 = treeT3->GetListOfBranches();
    int numBranchesT3 = branchesT3->GetEntries();


    if (numBranchesT != numBranchesT3) {
        std::cout << "Liczba gałęzi w drzewach T i T3 jest różna." << std::endl;
        return;
    }

    for (int i = 0; i < numBranchesT; ++i) {
        TBranch* branchT = static_cast<TBranch*>(branchesT->At(i));
        TBranch* branchT3 = static_cast<TBranch*>(branchesT3->At(i));

        const char* branchName = branchT->GetName();
        int numEntriesBranchT = branchT->GetEntries();
        int numEntriesBranchT3 = branchT3->GetEntries();

        if (numEntriesBranchT != numEntriesBranchT3) {
            std::cout << "Liczba wpisów w gałęzi \"" << branchName << "\" jest różna." << std::endl;
            continue;
        }

        TLeaf* leafT = branchT->GetLeaf(branchName);
        TLeaf* leafT3 = branchT3->GetLeaf(branchName);

        if (!leafT || !leafT3) {
            std::cout << "Nie można pobrać liścia dla gałęzi \"" << branchName << "\"." << std::endl;
            continue;
        }

        // Pobieranie typu danych liścia
        const char* leafType = leafT->GetTypeName();

        // Porównywanie danych w liściach
        bool areEqual = true;
        for (Long64_t j = 0; j < numEntriesBranchT; ++j) {
            leafT->GetBranch()->GetEntry(j);
            leafT3->GetBranch()->GetEntry(j);

            // Porównywanie wartości liścia w zależności od typu danych
            if (strcmp(leafType, "Float_t") == 0) {
                float valueT = *(static_cast<float*>(leafT->GetValuePointer()));
                float valueT3 = *(static_cast<float*>(leafT3->GetValuePointer()));
                if (valueT != valueT3) {
                    std::cout << "Różnica w gałęzi \"" << branchName << "\" dla wpisu " << j << std::endl;
                    areEqual = false;
                }
            } else if (strcmp(leafType, "Double_t") == 0) {
                double valueT = *(static_cast<double*>(leafT->GetValuePointer()));
                double valueT3 = *(static_cast<double*>(leafT3->GetValuePointer()));
                if (valueT != valueT3) {
                    std::cout << "Różnica w gałęzi \"" << branchName << "\" dla wpisu " << j << std::endl;
                    areEqual = false;
                }
            } else if (strcmp(leafType, "Int_t") == 0) {
                int valueT = *(static_cast<int*>(leafT->GetValuePointer()));
                int valueT3 = *(static_cast<int*>(leafT3->GetValuePointer()));
                if (valueT != valueT3) {
                    std::cout << "Różnica w gałęzi \"" << branchName << "\" dla wpisu " << j << std::endl;
                    areEqual = false;
                }
            }
        }

        if (areEqual) {
            std::cout << "Gałąź \"" << branchName << "\" jest taka sama w obu drzewach." << std::endl;
        }
    }

    file1->Close();
    file2->Close();

    delete file1;
    delete file2;
}
