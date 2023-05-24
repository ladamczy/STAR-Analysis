#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TKey.h>

//	Kod czytający dowolne drzewo z pliku

void PrintBranches(const char* fileName) {
   TFile file(fileName, "READ");
   if (!file.IsOpen()) {
      std::cerr << "Nie można otworzyć pliku ROOT." << std::endl;
      return;
   }

   TTree* tree = nullptr;
   TIter nextKey(file.GetListOfKeys());
   TKey* key;
   while ((key = (TKey*)nextKey())) {
      if (key->GetClassName() == TString("TTree")) {
         tree = (TTree*)key->ReadObj();
         break;
      }
   }

   if (!tree) {
      std::cerr << "Nie znaleziono drzewa w pliku ROOT." << std::endl;
      file.Close();
      return;
   }

   TObjArray* branches = tree->GetListOfBranches();
   if (!branches) {
      std::cerr << "Nie znaleziono gałęzi w drzewie." << std::endl;
      file.Close();
      return;
   }

   std::cout << "Lista gałęzi w drzewie " << tree->GetName() << ":" << std::endl;
   for (int i = 0; i < branches->GetEntries(); ++i) {
      TObjString* branchName = (TObjString*)branches->At(i);
      std::cout << branchName->GetString().Data() << std::endl;
   }

   file.Close();
}

int main(void) {
   const char* fileName = "treepart1.root";
   PrintBranches(fileName);
   return 0;
}
