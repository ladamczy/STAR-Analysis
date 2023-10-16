//default headers
#include <vector>

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

void loopdir(TDirectory*, std::vector<TKey*>*, std::vector<TKey*>*);

int ROOTtoPDF(const char* filename)
{
    //new file
    TFile* fileToPDF = TFile::Open(filename, "READ");

    //vectors for TH1 & TH2 histograms
    std::vector<TKey*> OneDimHists;
    std::vector<TKey*> TwoDimHists;

    //filling the vectors
    loopdir(fileToPDF, &OneDimHists, &TwoDimHists);

    //drawing them in a PDF file
    //setting style
    gStyle->SetHistLineWidth(1);
    gStyle->SetFrameLineWidth(3);
    gStyle->SetOptStat(111111);
    gStyle->SetStatY(0.89);
    gStyle->SetStatX(0.89);
    gStyle->SetStatW(0.1);
    gStyle->SetStatH(0.1); 
    gROOT->ForceStyle();
    gROOT->SetBatch(kTRUE);
    TCanvas* c1 = new TCanvas("c1","c1",1600,1200);
    TPaveStats *st;
    TH1D* temp1D;
    TH2D* temp2D;
    c1->Print("result.pdf[");

    //looping through both vectors
    for (int i = 0; i < OneDimHists.size(); i++){
        temp1D = (TH1D*)fileToPDF->Get(OneDimHists[i]->GetName());
        // temp1D->SetMarkerStyle(kFullCircle);
        // temp1D->Draw("E");
        c1->SetLogy(1);
        c1->SetLogz(0);
        temp1D->Draw("HIST");
        gPad->RedrawAxis();
        c1->Print("result.pdf");
    }
    for (int i = 0; i < TwoDimHists.size(); i++){
        temp2D = (TH2D*)fileToPDF->Get(TwoDimHists[i]->GetName());
        c1->SetLogy(0);
        c1->SetLogz(1);
        temp2D->Draw("COLZ");
        gPad->RedrawAxis();
        c1->Print("result.pdf");
    }
    c1->Print("result.pdf]");

    fileToPDF->Close();
    OneDimHists.clear();
    TwoDimHists.clear();

    return 0;
}

void loopdir(TDirectory *dir, std::vector<TKey*>* onedim, std::vector<TKey*>* twodim){
//loop on all keys of dir including possible subdirs
//print a message for all keys with a class TH1
//code taken from https://root-forum.cern.ch/t/finding-all-objects-in-a-root-file-regardless-of-directries/2351
TIter next (dir->GetListOfKeys());
TKey* key;
    while ((key = (TKey*)next())) {
        if(strstr(key->GetClassName(), "TH1")){
            printf(" key : %s is a %s in %s\n", key->GetName(), key->GetClassName(), dir->GetPath());
            onedim->push_back(key);
        }
        if(strstr(key->GetClassName(), "TH2")){
            printf(" key : %s is a %s in %s\n", key->GetName(), key->GetClassName(), dir->GetPath());
            twodim->push_back(key);
        }
        if(!strcmp(key->GetClassName(),"TDirectory")){
            dir->cd(key->GetName());
            TDirectory *subdir = gDirectory;
            loopdir(subdir, onedim, twodim);
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