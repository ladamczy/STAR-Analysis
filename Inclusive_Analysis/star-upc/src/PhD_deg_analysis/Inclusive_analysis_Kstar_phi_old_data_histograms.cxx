// Stdlib header file for input and output.
#include <iostream>
#include <cstring>

// ROOT, for histogramming.
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TEfficiency.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TStyle.h"

void draw_and_save(TH1D* data, std::string folderWithDiagonal, std::string name);

int main(int argc, char const* argv[]){
    TFile* input = TFile::Open(argv[1]);

    //getting histograms out
    TH1D* pairInfo = (TH1D*)input->Get("pairInfo");
    TH1D* MKpiChi2Wide = (TH1D*)input->Get("MKpiChi2Wide");
    TH1D* MpiKChi2Wide = (TH1D*)input->Get("MpiKChi2Wide");
    TH1D* MppiChi2Wide = (TH1D*)input->Get("MppiChi2Wide");
    TH1D* MpipChi2Wide = (TH1D*)input->Get("MpipChi2Wide");
    TH1D* MKKChi2Wide = (TH1D*)input->Get("MKKChi2Wide");
    TH1D* MpipiChi2Wide = (TH1D*)input->Get("MpipiChi2Wide");
    TH1D* MppChi2Wide = (TH1D*)input->Get("MppChi2Wide");
    /*
    MKpiChi2Wide
    MpiKChi2Wide
    MppiChi2Wide
    MpipChi2Wide
    MKKChi2Wide
    MpipiChi2Wide
    MppChi2Wide
    */

    std::string folderWithDiagonal = std::string(argv[2]);
    if(folderWithDiagonal[folderWithDiagonal.size()-1]!='/'){
        folderWithDiagonal += "/";
    }

    gStyle->SetOptStat(0);
    gStyle->SetFrameLineWidth(2);

    //drawing
    //pairInfo
    TCanvas* resultCanvas = new TCanvas("resultCanvas", "resultCanvas", 4000, 2400);
    pairInfo->SetMinimum(0);
    pairInfo->SetLineColor(kBlue+2);
    pairInfo->SetMarkerSize(2);
    //TODO: add fitting
    // pairInfo->Draw("E");
    pairInfo->Draw("hist");
    pairInfo->Draw("same text0");
    pairInfo->SetLineWidth(2);
    pairInfo->GetXaxis()->SetLabelSize(0.06);
    pairInfo->GetXaxis()->SetTitleSize(0.06);
    resultCanvas->SetTopMargin(0.05);
    gPad->Update();
    resultCanvas->SaveAs((folderWithDiagonal+"pairInfo.pdf").c_str());

    //the rest
    draw_and_save(MKpiChi2Wide, folderWithDiagonal, "MKpi");
    draw_and_save(MpiKChi2Wide, folderWithDiagonal, "MpiK");
    draw_and_save(MppiChi2Wide, folderWithDiagonal, "Mppi");
    draw_and_save(MpipChi2Wide, folderWithDiagonal, "Mpip");
    draw_and_save(MKKChi2Wide, folderWithDiagonal, "MKK");
    draw_and_save(MpipiChi2Wide, folderWithDiagonal, "Mpipi");
    draw_and_save(MppChi2Wide, folderWithDiagonal, "Mpp");

    return 0;
}

void draw_and_save(TH1D* data, std::string folderWithDiagonal, std::string name){
    TCanvas* resultCanvas = new TCanvas("resultCanvas", "resultCanvas", 4000, 2400);
    data->SetMinimum(0);
    data->SetMarkerStyle(kFullCircle);
    data->SetMarkerSize(2);
    data->SetMarkerColor(kBlue);
    data->SetLineColor(kBlue+2);
    //TODO: add fitting
    data->Draw("hist");
    data->SetLineWidth(2);
    data->GetXaxis()->SetLabelSize(0.05);
    data->GetXaxis()->SetTitleSize(0.05);
    data->GetYaxis()->SetLabelSize(0.05);
    data->GetYaxis()->SetTitleSize(0.05);
    resultCanvas->SetMargin(0.15, 0.05, 0.15, 0.05);
    gPad->Update();
    resultCanvas->SaveAs((folderWithDiagonal+name+".pdf").c_str());
}