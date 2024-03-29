//default headers
#include <vector>
#include <string>

//ROOT headers
#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "TKey.h"
#include "TFile.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPaveStats.h"
#include "TEfficiency.h"

#include "Styles.h"

void loopdir(TDirectory*, std::vector<TKey*>*);

int ROOTtoPDF(const char* filename, std::string outputfile, int mainconfig, int oneconfig, int twoconfig)
{
    //new file
    TFile* fileToPDF = TFile::Open(filename, "READ");

    //vectors for histograms
    std::vector<TKey*> Hists;

    //filling the vectors
    loopdir(fileToPDF, &Hists);

    //drawing them in a PDF file
    //setting style
    MainStyle(gStyle, mainconfig);
    
    gROOT->ForceStyle();
    gROOT->SetBatch(kTRUE);
    TCanvas* c1 = new TCanvas("c1","c1",1600,1200);
    TPaveStats *st;
    TH1D* temp1D;
    TH2D* temp2D;
    TEfficiency *tempEff;
    if(outputfile.length()==0){
        outputfile = std::string(filename).substr(0, std::string(filename).find("."));
    }
    if(outputfile.find(".pdf")==std::string::npos){
        outputfile = outputfile + ".pdf";
    }
    c1->Print((outputfile+"[").c_str());

    //looping through both vectors
    for (int i = 0; i < Hists.size(); i++){
        //drawing histograms
        if(strstr(Hists[i]->GetClassName(), "TH1")){
            temp1D = (TH1D*)fileToPDF->Get(Hists[i]->GetName());
            TH1Style(c1, temp1D, oneconfig);

        }
        if(strstr(Hists[i]->GetClassName(), "TH2")){
            temp2D = (TH2D*)fileToPDF->Get(Hists[i]->GetName());
            TH2Style(c1, temp2D, twoconfig);
        }
        //TODO
        if(strstr(Hists[i]->GetClassName(), "TEfficiency")){
            tempEff = (TEfficiency *)fileToPDF->Get(Hists[i]->GetName());
            c1->SetLogz(0);
            tempEff->Draw();
        }
        //actual drawing to file
        gPad->RedrawAxis();
        c1->Print(outputfile.c_str());
    }
    c1->Print((outputfile+"]").c_str());

    fileToPDF->Close();
    Hists.clear();

    return 0;
}

void loopdir(TDirectory *dir, std::vector<TKey*>* hists){
//loop on all keys of dir including possible subdirs
//print a message for all keys with a class TH1 or TH2
//code taken from https://root-forum.cern.ch/t/finding-all-objects-in-a-root-file-regardless-of-directries/2351
TIter next (dir->GetListOfKeys());
TKey* key;
    while ((key = (TKey*)next())) {
        if(strstr(key->GetClassName(), "TH1")||strstr(key->GetClassName(), "TH2")||strstr(key->GetClassName(), "TEfficiency")){
            printf(" key : %s is a %s in %s\n", key->GetName(), key->GetClassName(), dir->GetPath());
            hists->push_back(key);
        }
        if(!strcmp(key->GetClassName(),"TDirectory")){
            dir->cd(key->GetName());
            TDirectory *subdir = gDirectory;
            loopdir(subdir, hists);
            dir->cd();
        }
    }
}

//to use maybe later, during code refactoring
    // //otwieranie nowego pliku
    // TGFileInfo* fileinfo = new TGFileInfo();
    // fileinfo->fFileTypeIdx = 2;
    // TGFileDialog* filedialog = new TGFileDialog(gClient->GetRoot(), nullptr, EFileDialogMode::kFDOpen, fileinfo);
    // cout<<fileinfo->fFilename<<endl;
    // TFile* anaoutput = TFile::Open(fileinfo->fFilename);
    // //zapisuje nazwę bez kropki
    // string filename = string(fileinfo->fFilename).substr(0, string(fileinfo->fFilename).find("."));
    // delete fileinfo;
    // delete filedialog;