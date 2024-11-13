#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1D.h"
#include "TStyle.h"

#include "string.h"

int main(int argc, char const *argv[]){
    TFile *data = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/AnaOutput_Inclusive_analysis_Kstar_phi_old_data.root");

    TH1D *signal = (TH1D *)data->Get("MKKWideNoVeto");
    TH1D *bcg = (TH1D *)data->Get("MKKWideBackground");
    signal->SetName("signal");
    bcg->SetName("bcg");

    //1st idea
    // int left_bin = bcg->FindBin(1.5);
    // int right_bin = bcg->FindBin(4.9);
    // printf("Signal: %lf\nBackground: %lf\nRatio S/B: %lf\n", signal->Integral(left_bin, right_bin), bcg->Integral(left_bin, right_bin), signal->Integral(left_bin, right_bin)/bcg->Integral(left_bin, right_bin));
    // bcg->Scale(signal->Integral(left_bin, right_bin)/bcg->Integral(left_bin, right_bin));
    //2nd idea
    int highest_bin = bcg->GetMaximumBin();
    bcg->Scale(signal->GetBinContent(highest_bin)/bcg->GetBinContent(highest_bin));

    TH1D *difference = (TH1D *)signal->Clone();
    difference->Add(bcg, -1.0);

    //signal & background
    TCanvas *resultCanvas = new TCanvas("resultCanvas", "resultCanvas", 1800, 1600);
    gStyle->SetOptStat(0);
    signal->SetMinimum(0);
    signal->SetLineColor(kBlue);
    signal->SetLineWidth(2);
    signal->Draw("hist");
    bcg->SetLineColor(kRed);
    bcg->SetLineWidth(2);
    bcg->Draw("Hist same");
    TLegend *legend = new TLegend(0.7, 0.7, 0.89, 0.89);
    legend->SetTextSize(0.025);
    legend->SetHeader("#bf{pp, #sqrt{s} = 510 GeV}", "C");
    legend->AddEntry(signal->GetName(), "signal");
    legend->AddEntry(bcg->GetName(), "background");
    legend->SetBorderSize(0);
    legend->Draw("SAME");
    gPad->Update();
    resultCanvas->SaveAs("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/AnaOutput_Inclusive_analysis_Kstar_phi_old_data.pdf");

    //difference
    resultCanvas->Clear();
    difference->SetLineColor(kBlue);
    difference->SetLineWidth(2);
    difference->Draw("hist");
    // TLegend *legend = new TLegend(0.7, 0.7, 0.89, 0.89);
    // legend->SetTextSize(0.025);
    // legend->SetHeader("#bf{pp, #sqrt{s} = 510 GeV}", "C");
    // legend->AddEntry(signal->GetName(), "signal");
    // legend->AddEntry(bcg->GetName(), "background");
    // legend->SetBorderSize(0);
    // legend->Draw("SAME");
    gPad->Update();
    resultCanvas->SaveAs("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/AnaOutput_Inclusive_analysis_Kstar_phi_old_data_diff.pdf");
    return 0;
}
