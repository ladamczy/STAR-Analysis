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

#include "MyStyles.h"

void draw_and_save(TH1D* data, std::string folderWithDiagonal, std::string name, std::string title);
void draw_and_save_minus_background(TH1D* data, TH1D* bcg, std::string folderWithDiagonal, std::string name, std::string title, double bcg_region);

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

    TFile* inputbcg = TFile::Open(argv[2]);

    //getting background histograms out
    //MIXEDTOF, NORMALTOF, NOTOF
    std::string background = "NORMALTOF";
    TH1D* MKpiChi2bcg = (TH1D*)inputbcg->Get(("MKpiChi2bcg"+background).c_str());
    TH1D* MpiKChi2bcg = (TH1D*)inputbcg->Get(("MpiKChi2bcg"+background).c_str());
    TH1D* MppiChi2bcg = (TH1D*)inputbcg->Get(("MppiChi2bcg"+background).c_str());
    TH1D* MpipChi2bcg = (TH1D*)inputbcg->Get(("MpipChi2bcg"+background).c_str());
    TH1D* MKKChi2bcg = (TH1D*)inputbcg->Get(("MKKChi2bcg"+background).c_str());
    TH1D* MpipiChi2bcg = (TH1D*)inputbcg->Get(("MpipiChi2bcg"+background).c_str());
    TH1D* MppChi2bcg = (TH1D*)inputbcg->Get(("MppChi2bcg"+background).c_str());

    std::string folderWithDiagonal = std::string(argv[3]);
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
    draw_and_save(MKpiChi2Wide, folderWithDiagonal, "MKpi", "K^{+}#pi^{-}");
    draw_and_save(MpiKChi2Wide, folderWithDiagonal, "MpiK", "#pi^{+}K^{-}");
    draw_and_save(MppiChi2Wide, folderWithDiagonal, "Mppi", "p^{+}#pi^{-}");
    draw_and_save(MpipChi2Wide, folderWithDiagonal, "Mpip", "#pi^{+}p^{-}");
    draw_and_save(MKKChi2Wide, folderWithDiagonal, "MKK", "K^{+}K^{-}");
    draw_and_save(MpipiChi2Wide, folderWithDiagonal, "Mpipi", "#pi^{+}#pi^{-}");
    draw_and_save(MppChi2Wide, folderWithDiagonal, "Mpp", "p^{+}p^{-}");
    draw_and_save_minus_background(MKpiChi2Wide, MKpiChi2bcg, folderWithDiagonal, "MKpiFit", "K^{+}#pi^{-}", 1.0);
    draw_and_save_minus_background(MpiKChi2Wide, MpiKChi2bcg, folderWithDiagonal, "MpiKFit", "#pi^{+}K^{-}", 1.0);
    draw_and_save_minus_background(MppiChi2Wide, MppiChi2bcg, folderWithDiagonal, "MppiFit", "p^{+}#pi^{-}", 1.1);
    draw_and_save_minus_background(MpipChi2Wide, MpipChi2bcg, folderWithDiagonal, "MpipFit", "#pi^{+}p^{-}", 1.1);
    draw_and_save_minus_background(MKKChi2Wide, MKKChi2bcg, folderWithDiagonal, "MKKFit", "K^{+}K^{-}", 1.1);
    draw_and_save_minus_background(MpipiChi2Wide, MpipiChi2bcg, folderWithDiagonal, "MpipiFit", "#pi^{+}#pi^{-}", 0.8);
    draw_and_save_minus_background(MppChi2Wide, MppChi2bcg, folderWithDiagonal, "MppFit", "p^{+}p^{-}", 2.4);
    draw_and_save(MKpiChi2bcg, folderWithDiagonal, "MKpibcg", "K^{+}#pi^{-} - background");
    draw_and_save(MpiKChi2bcg, folderWithDiagonal, "MpiKbcg", "#pi^{+}K^{-} - background");
    draw_and_save(MppiChi2bcg, folderWithDiagonal, "Mppibcg", "p^{+}#pi^{-} - background");
    draw_and_save(MpipChi2bcg, folderWithDiagonal, "Mpipbcg", "#pi^{+}p^{-} - background");
    draw_and_save(MKKChi2bcg, folderWithDiagonal, "MKKbcg", "K^{+}K^{-} - background");
    draw_and_save(MpipiChi2bcg, folderWithDiagonal, "Mpipibcg", "#pi^{+}#pi^{-} - background");
    draw_and_save(MppChi2bcg, folderWithDiagonal, "Mppbcg", "p^{+}p^{-} - background");

    return 0;
}

void draw_and_save(TH1D* data, std::string folderWithDiagonal, std::string name, std::string title){
    TCanvas* resultCanvas = new TCanvas("resultCanvas", "resultCanvas", 4000, 2400);
    data->SetMinimum(0);
    data->SetMarkerStyle(kFullCircle);
    data->SetMarkerSize(2);
    data->SetMarkerColor(kBlue);
    data->SetLineColor(kBlue+2);
    data->SetTitle(title.c_str());
    //TODO: add fitting
    data->Draw("hist");
    data->SetLineWidth(2);
    data->GetXaxis()->SetLabelSize(0.05);
    data->GetXaxis()->SetTitleSize(0.05);
    data->GetYaxis()->SetLabelSize(0.05);
    data->GetYaxis()->SetTitleSize(0.05);
    resultCanvas->SetMargin(0.15, 0.05, 0.15, 0.10);
    data->SetTitle(title.c_str());
    // gPad->Update();
    resultCanvas->SaveAs((folderWithDiagonal+name+".pdf").c_str());
}

void draw_and_save_minus_background(TH1D* data, TH1D* bcg, std::string folderWithDiagonal, std::string name, std::string title, double bcg_region){
    MyStyles styleLibrary;
    TStyle tempStyle = styleLibrary.Hist2DQuarterSize(true);
    tempStyle.cd();
    gROOT->ForceStyle();
    TCanvas* resultCanvas = new TCanvas("resultCanvas", "resultCanvas", 4000, 2400);
    int bcg_bin = data->FindBin(bcg_region);
    bcg->Scale(data->Integral(bcg_bin, -1)/bcg->Integral(bcg_bin, -1));
    data->Add(bcg, -1.);
    data->SetMarkerStyle(kFullCircle);
    data->SetMarkerSize(2);
    data->SetMarkerColor(kBlue);
    data->SetTitle(title.c_str());
    //TODO: add fitting
    data->Draw("e");

    resultCanvas->SaveAs((folderWithDiagonal+name+".pdf").c_str());
}