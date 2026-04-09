#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TDirectory.h>
#include <iostream>
#include <string>
#include <TGraphAsymmErrors.h>

#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TF1.h>
#include <vector>
#include <string>
#include <iostream>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TMath.h>
#include <TLatex.h>
#include <iomanip>

void sethist(TH1* h1, TH1* h2) {
    h1->SetLineColor(kRed);
    h1->SetMarkerStyle(20);
    h1->SetMarkerSize(1.2); // Reduced size for side-by-side visibility
    h1->SetMarkerColor(kRed);
    h2->SetLineColor(kBlue);
    h2->SetMarkerStyle(21);
    h2->SetMarkerSize(1.2);
    h2->SetMarkerColor(kBlue);
}

void Plot_EffCorrections()
{
    string rootfile = "PhyEffCorr.root";
    TFile* fcorr = TFile::Open(rootfile.c_str(), "READ");

    // Load Histograms
    TH1D* hPt_P_C = (TH1D*)fcorr->Get("h1D_CorrectedPt_P");
    TH1D* hPt_P_R = (TH1D*)fcorr->Get("h1D_RawPt_P");
    TH1D* hPt_N_C = (TH1D*)fcorr->Get("h1D_CorrectedPt_N");
    TH1D* hPt_N_R = (TH1D*)fcorr->Get("h1D_RawPt_N");

    TH1D* hEta_P_C = (TH1D*)fcorr->Get("h1D_CorrectedEta_P");
    TH1D* hEta_P_R = (TH1D*)fcorr->Get("h1D_RawEta_P");
    TH1D* hEta_N_C = (TH1D*)fcorr->Get("h1D_CorrectedEta_N");
    TH1D* hEta_N_R = (TH1D*)fcorr->Get("h1D_RawEta_N");

    TH1D* hVerZ_P_C = (TH1D*)fcorr->Get("h1D_CorrectedVz_P");
    TH1D* hVerZ_P_R = (TH1D*)fcorr->Get("h1D_RawVz_P");
    TH1D* hVerZ_N_C = (TH1D*)fcorr->Get("h1D_CorrectedVz_N");
    TH1D* hVerZ_N_R = (TH1D*)fcorr->Get("h1D_RawVz_N");

    TH1F *hFinalMass = (TH1F*)fcorr->Get("h1D_Reco_InvMass_Corrected_4pi_4TOF")->Clone("hFinalMass");
    hFinalMass->Add( (TH1F*)fcorr->Get("h1D_Reco_InvMass_Corrected_4pi_3TOF") );
    hFinalMass->Add( (TH1F*)fcorr->Get("h1D_Reco_InvMass_Corrected_4pi_2TOF") );

    sethist(hPt_P_C, hPt_P_R);
    sethist(hPt_N_C, hPt_N_R);
    sethist(hEta_P_C, hEta_P_R);
    sethist(hEta_N_C, hEta_N_R);
    sethist(hVerZ_P_C, hVerZ_P_R);
    sethist(hVerZ_N_C, hVerZ_N_R);

    //Gstyle
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    // --- Plotting PT (Positive vs Negative) ---
    TCanvas* cPt = new TCanvas("cPt", "pT Comparison", 1600, 800);
    cPt->Divide(2, 1);

    cPt->cd(1); // Left side: Positive
    auto rp1 = new TRatioPlot(hPt_P_C, hPt_P_R, "divsym");
    rp1->Draw();
    rp1->GetLowerRefGraph()->SetMinimum(-1.0);
    rp1->GetLowerRefGraph()->SetMaximum(5.0);
    rp1->GetUpperPad()->SetLogy(); // Often useful for pT distributions

    cPt->cd(2); // Right side: Negative
    auto rp2 = new TRatioPlot(hPt_N_C, hPt_N_R, "divsym");
    rp2->Draw();
    rp2->GetLowerRefGraph()->SetMinimum(-1.0);
    rp2->GetLowerRefGraph()->SetMaximum(5.0);
    rp2->GetUpperPad()->SetLogy();

    cPt->SaveAs("plots/test_pT_PosNeg.png");

    // --- Plotting ETA (Positive vs Negative) ---
    TCanvas* cEta = new TCanvas("cEta", "Eta Comparison", 1600, 800);
    cEta->Divide(2, 1);

    cEta->cd(1);
    auto rp3 = new TRatioPlot(hEta_P_C, hEta_P_R, "divsym");
    rp3->Draw();
    rp3->GetLowerRefGraph()->SetMinimum(-1.0);
    rp3->GetLowerRefGraph()->SetMaximum(5.0);

    cEta->cd(2);
    auto rp4 = new TRatioPlot(hEta_N_C, hEta_N_R, "divsym");
    rp4->Draw();
    rp4->GetLowerRefGraph()->SetMinimum(-1.0);
    rp4->GetLowerRefGraph()->SetMaximum(5.0);

    cEta->SaveAs("plots/test_Eta_PosNeg.png");

    TCanvas* cVerZ = new TCanvas("cVerZ", "VerZ Comparison", 1600, 800);
    cVerZ->Divide(2, 1);

    cVerZ->cd(1);
    auto rp5 = new TRatioPlot(hVerZ_P_C, hVerZ_P_R, "divsym");
    rp5->Draw();
    rp5->GetLowerRefGraph()->SetMinimum(-1.0);
    rp5->GetLowerRefGraph()->SetMaximum(5.0);

    cVerZ->cd(2);
    auto rp6 = new TRatioPlot(hVerZ_N_C, hVerZ_N_R, "divsym");
    rp6->Draw();
    rp6->GetLowerRefGraph()->SetMinimum(-1.0);
    rp6->GetLowerRefGraph()->SetMaximum(5.0);

    cVerZ->SaveAs("plots/test_VerZ_PosNeg.png");

    // --- Plotting Final Mass ---
    TCanvas* cMass = new TCanvas("cMass", "Final Mass Distribution", 800, 600);
    hFinalMass->SetMarkerStyle(20);
    hFinalMass->SetMarkerColor(kBlack);
    hFinalMass->SetLineColor(kBlack);
    hFinalMass->GetXaxis()->SetRangeUser(0.9, 3.0);
    hFinalMass->GetXaxis()->SetTitle("m_{K_{S}^{0}K_{S}^{0}} [GeV/c^{2}]");
    hFinalMass->GetYaxis()->SetTitle("events");
    hFinalMass->GetYaxis()->SetTitleOffset(0.9);
    hFinalMass->GetXaxis()->SetTitleOffset(0.6);
    hFinalMass->GetXaxis()->SetTitleSize(0.06); //axis title size
    hFinalMass->GetYaxis()->SetTitleSize(0.06);

    hFinalMass->GetXaxis()->SetLabelSize(0.05);
    hFinalMass->GetYaxis()->SetLabelSize(0.05);
    
    hFinalMass->Draw("PE");

    TF1* fitFunc = new TF1("fitFunc", "gaus(0) + pol1(3)", 1.0, 2.5);
    fitFunc->SetParameters(15, 1.68, 0.08, 10, -2); 
    fitFunc->SetLineColor(kRed);
    hFinalMass->Fit(fitFunc, "R"); 

    TLegend* leg = new TLegend(0.55, 0.75, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.045);
    leg->AddEntry(hFinalMass, "Data", "pe");
    leg->AddEntry(fitFunc, "Fit: gaus + pol1", "l");
    leg->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.045);
    latex.DrawLatex(0.55, 0.68, Form("A = %.2f", fitFunc->GetParameter(0)));
    latex.DrawLatex(0.55, 0.62, Form("#mu = %.2f GeV/c^{2}", fitFunc->GetParameter(1)));
    latex.DrawLatex(0.55, 0.56, Form("#sigma = %.2f GeV/c^{2}", fitFunc->GetParameter(2)));

    cMass->SaveAs("plots/test_FinalMass.png");
}

