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

void For1atatime (TH1D* corrected, TH1D* raw, const char* xtitle, double pullMinY1, double pullMaxY2) {
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0.012);
    pad1->SetLeftMargin(0.15);
    pad1->SetRightMargin(0.05);
    pad1->Draw();
    pad1->cd();
    corrected->SetLineColor(kBlue);
    corrected->SetMarkerColor(kBlue);
    corrected->SetMarkerStyle(22);
    corrected->SetLineWidth(2);
    corrected->SetMarkerSize(1.5);
    raw->SetLineColor(kRed);
    raw->SetMarkerColor(kRed);
    raw->SetMarkerStyle(23);
    raw->SetLineWidth(2);
    raw->SetMarkerSize(1.5);
    corrected->GetYaxis()->SetMaxDigits(3);
    corrected->GetXaxis()->SetTitleOffset(1.5);

    TLegend *leg4 = new TLegend(0.69,0.69,0.88,0.88);
    leg4->SetBorderSize(0);
    leg4->SetFillStyle(0);
    leg4->SetTextSize(0.07);
    leg4->AddEntry(corrected, "Corrected", "lep");
    leg4->AddEntry(raw, "Raw", "lep");
    corrected->GetYaxis()->SetLabelSize(0.06);
    corrected->GetYaxis()->SetTitleSize(0.05);
    //corrected->GetYaxis()->SetRangeUser(ymin, ymax);
    corrected->Draw("E");
    raw->Draw("E SAME");
    leg4->Draw();
    gPad->GetMother()->cd();
    // Create lower pad for ratio
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
    pad2->SetTopMargin(0.018);
    pad2->SetBottomMargin(0.35);
    pad2->SetLeftMargin(0.15);
    pad2->SetRightMargin(0.05);
    pad2->Draw();
    pad2->cd();
    //draw ratio of corrected to raw
    TH1D *ratio = (TH1D*)corrected->Clone("ratio");
    ratio->Divide(raw);
    ratio->SetLineColor(kBlue);
    ratio->SetMarkerColor(kBlue);
    ratio->SetMarkerStyle(22);
    ratio->SetLineWidth(2);
    ratio->SetMarkerSize(1.5);
    ratio->GetYaxis()->SetTitle("Corrected/Raw");
    ratio->GetYaxis()->SetTitleSize(0.12);
    ratio->GetYaxis()->SetLabelSize(0.10);
    ratio->GetYaxis()->SetTitleOffset(0.5);
    ratio->GetYaxis()->SetNdivisions(505);
    ratio->GetYaxis()->SetRangeUser(pullMinY1, pullMaxY2);
    //ratio->GetXaxis()->SetTitle("m_{#pi^{+} #pi^{-}} (GeV/c^{2})");
    ratio->GetXaxis()->SetTitle(xtitle);
    ratio->GetXaxis()->SetTitleSize(0.15);
    ratio->GetXaxis()->SetLabelSize(0.15);
    ratio->GetXaxis()->SetTitleOffset(1.0);
    ratio->SetTitle("");
    ratio->Draw("E");

    // Draw reference line at 1.0
    TLine *line = new TLine(ratio->GetXaxis()->GetXmin(), 0.0, 
                            ratio->GetXaxis()->GetXmax(), 0.0);
    line->SetLineColor(kRed);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw("SAME");

    // Draw ±1σ bands
    TLine *line1sigma_up = new TLine(ratio->GetXaxis()->GetXmin(), 3.0, 
                                    ratio->GetXaxis()->GetXmax(), 3.0);
    line1sigma_up->SetLineColor(kGray);
    line1sigma_up->SetLineStyle(2);
    line1sigma_up->Draw("SAME");

    TLine *line1sigma_down = new TLine(ratio->GetXaxis()->GetXmin(), -3.0, 
                                      ratio->GetXaxis()->GetXmax(), -3.0);
    line1sigma_down->SetLineColor(kGray);
    line1sigma_down->SetLineStyle(2);
    line1sigma_down->Draw("SAME");
}

void PlotInvMass(TH1D* corr,TH1D* raw, const char* xtitle, double fitmin, double fitmax, double pullMinY1, double pullMaxY2) {
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0.012);
    pad1->SetLeftMargin(0.15);
    pad1->SetRightMargin(0.05);
    pad1->Draw();
    pad1->cd();
    
    corr->SetMarkerStyle(20);
    corr->SetMarkerColor(kBlack);
    corr->SetLineColor(kBlack);
    corr->GetXaxis()->SetRangeUser(0.9, 3.0);
    corr->GetXaxis()->SetTitle("m_{K_{S}^{0}K_{S}^{0}} [GeV/c^{2}]");
    corr->GetYaxis()->SetTitle("events");
    raw->GetYaxis()->SetTitleOffset(0.9);
    raw->GetXaxis()->SetTitleOffset(0.6);
    raw->GetXaxis()->SetTitleSize(0.06); //axis title size
    raw->GetYaxis()->SetTitleSize(0.06);
    raw->SetTitle("");
    raw->SetLineColor(kGreen);
    raw->SetLineWidth(3);
    raw->SetMarkerStyle(21);
    raw->SetMarkerColor(kGreen);
    corr->GetXaxis()->SetLabelSize(0.05);
    corr->GetYaxis()->SetLabelSize(0.05);
    //set log scale for better visibility of the background
    //cMass->SetLogy();
    corr->Draw("PE");
    TLegend *leg4 = new TLegend(0.69,0.69,0.88,0.88);
    leg4->SetBorderSize(0);
    leg4->SetFillStyle(0);
    leg4->SetTextSize(0.07);
    leg4->AddEntry(corr, "Corrected", "lep");
    leg4->AddEntry(raw, "Raw", "lep");
    leg4->Draw(); 

    TF1* fitFunc = new TF1("fitFunc", "gaus(0) + pol1(3)", fitmin, fitmax);
    fitFunc->SetParameters(15, 1.68, 0.08, 10, -2); 
    fitFunc->SetLineColor(kRed);
    corr->Fit(fitFunc, "R"); 

    //The background part separately as a dashed line to show what it's doing
    TF1 *bgFunc = new TF1("bgFunc", "pol1(0)", 1.0, 2.4);
    bgFunc->SetParameter(0, fitFunc->GetParameter(3));
    bgFunc->SetParameter(1, fitFunc->GetParameter(4));
    bgFunc->SetLineStyle(2); // Dashed
    bgFunc->SetLineColor(kBlue);
    bgFunc->Draw("SAME");

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.045);
    latex.DrawLatex(0.55, 0.68, Form("A = %.2f", fitFunc->GetParameter(0)));
    latex.DrawLatex(0.55, 0.62, Form("#mu = %.2f GeV/c^{2}", fitFunc->GetParameter(1)));
    latex.DrawLatex(0.55, 0.56, Form("#sigma = %.2f GeV/c^{2}", fitFunc->GetParameter(2)));

    raw->Draw("PE same");
    corr->Draw("PE same");
    gPad->GetMother()->cd();
    // Create lower pad for ratio
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
    pad2->SetTopMargin(0.018);
    pad2->SetBottomMargin(0.35);
    pad2->SetLeftMargin(0.15);
    pad2->SetRightMargin(0.05);
    pad2->Draw();
    pad2->cd();
    //draw ratio of corrected to raw
    TH1D *ratio = (TH1D*)corr->Clone("ratio");
    ratio->GetListOfFunctions()->Clear(); 
    ratio->Divide(raw);
    ratio->SetLineColor(kBlue);
    ratio->SetMarkerColor(kBlue);
    ratio->SetMarkerStyle(22);
    ratio->SetLineWidth(2);
    ratio->SetMarkerSize(1.5);
    ratio->GetYaxis()->SetTitle("Corrected/Raw");
    ratio->GetYaxis()->SetTitleSize(0.12);
    ratio->GetYaxis()->SetLabelSize(0.10);
    ratio->GetYaxis()->SetTitleOffset(0.5);
    ratio->GetYaxis()->SetNdivisions(505);
    ratio->GetYaxis()->SetRangeUser(pullMinY1, pullMaxY2);
    //ratio->GetXaxis()->SetTitle("m_{#pi^{+} #pi^{-}} (GeV/c^{2})");
    ratio->GetXaxis()->SetTitleSize(0.15);
    ratio->GetXaxis()->SetLabelSize(0.15);
    ratio->GetXaxis()->SetTitleOffset(1.0);
    ratio->SetTitle("");
    ratio->Draw("PE");

    // Draw reference line at 1.0
    TLine *line = new TLine(ratio->GetXaxis()->GetXmin(), 0.0, 
                            ratio->GetXaxis()->GetXmax(), 0.0);
    line->SetLineColor(kRed);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw("SAME");

    // Draw ±1σ bands
    TLine *line1sigma_up = new TLine(ratio->GetXaxis()->GetXmin(), 3.0, 
                                    ratio->GetXaxis()->GetXmax(), 3.0);
    line1sigma_up->SetLineColor(kGray);
    line1sigma_up->SetLineStyle(2);
    line1sigma_up->Draw("SAME");

    TLine *line1sigma_down = new TLine(ratio->GetXaxis()->GetXmin(), 5.0, 
                                      ratio->GetXaxis()->GetXmax(), 5.0);
    line1sigma_down->SetLineColor(kGray);
    line1sigma_down->SetLineStyle(2);
    line1sigma_down->Draw("SAME");   
}

void Plot_EffCorrections()
{
    string DataFile = "EffCorrWithNewEffFileWithNewCap.root"; //EffCorrWithOldEffSameAsResultsApril14 EffCorrWithNewEffFile
    TFile* fcorr = TFile::Open(DataFile.c_str(), "READ");

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

    TH1D* h1D_Reco_InvMass_Corrected_4pi_4TOF = (TH1D*)fcorr->Get("h1D_Reco_InvMass_Corrected_4pi_4TOF");
    TH1D* h1D_Reco_InvMass_Corrected_4pi_3TOF = (TH1D*)fcorr->Get("h1D_Reco_InvMass_Corrected_4pi_3TOF");
    TH1D* h1D_Reco_InvMass_Corrected_4pi_2TOF = (TH1D*)fcorr->Get("h1D_Reco_InvMass_Corrected_4pi_2TOF");

    TH1D* h1D_Reco_InvMass_Raw_4pi_4TOF = (TH1D*)fcorr->Get("h1D_Reco_InvMass_Raw_4pi_4TOF");
    TH1D* h1D_Reco_InvMass_Raw_4pi_3TOF = (TH1D*)fcorr->Get("h1D_Reco_InvMass_Raw_4pi_3TOF");
    TH1D* h1D_Reco_InvMass_Raw_4pi_2TOF = (TH1D*)fcorr->Get("h1D_Reco_InvMass_Raw_4pi_2TOF");

    TH1D *hInvMassCorr = (TH1D*)h1D_Reco_InvMass_Corrected_4pi_4TOF->Clone("hInvMassCorr");
    hInvMassCorr->Add(h1D_Reco_InvMass_Corrected_4pi_3TOF);
    hInvMassCorr->Add(h1D_Reco_InvMass_Corrected_4pi_2TOF);
    hInvMassCorr->Scale(1.0 / 3.0); // Average over the three TOF categories

    TH1D *hInvMassRaw = (TH1D*)h1D_Reco_InvMass_Raw_4pi_4TOF->Clone("hInvMassRaw");
    hInvMassRaw->Add(h1D_Reco_InvMass_Raw_4pi_3TOF);
    hInvMassRaw->Add(h1D_Reco_InvMass_Raw_4pi_2TOF);

    sethist(hPt_P_C, hPt_P_R);
    sethist(hPt_N_C, hPt_N_R);
    sethist(hEta_P_C, hEta_P_R);
    sethist(hEta_N_C, hEta_N_R);
    sethist(hVerZ_P_C, hVerZ_P_R);
    sethist(hVerZ_N_C, hVerZ_N_R);

    gStyle->SetOptStat(0); // Disable stats box for cleaner plots

    // --- Plotting PT (Positive vs Negative) ---
    TCanvas* cPt = new TCanvas("cPt", "pT Comparison", 2000, 1000);
    cPt->Divide(2);
    cPt->cd(1);
    //set x axis title
    For1atatime(hPt_P_C, hPt_P_R, "p_{T} (GeV/c)", -1.0, 9.0);
    cPt->cd(2);
    For1atatime(hPt_N_C, hPt_N_R, "p_{T} (GeV/c)", -1.0, 9.0);
    cPt->SaveAs("plots/EffCorrData_pT_PosNeg3.png");

    // --- Plotting ETA (Positive vs Negative) ---
    TCanvas* cEta = new TCanvas("cEta", "Eta Comparison", 2000, 1000);
    cEta->Divide(2);
    cEta->cd(1);
    For1atatime(hEta_P_C, hEta_P_R, "#eta", -1.0, 9.0);
    cEta->cd(2);
    For1atatime(hEta_N_C, hEta_N_R, "#eta", -1.0, 9.0);
    cEta->SaveAs("plots/EffCorrData_Eta_PosNeg3.png");

    TCanvas* cVerZ = new TCanvas("cVerZ", "VerZ Comparison", 2000, 1000);
    cVerZ->Divide(2);
    cVerZ->cd(1);
    For1atatime(hVerZ_P_C, hVerZ_P_R, "Vz (cm)", -1.0, 9.0);
    cVerZ->cd(2);
    For1atatime(hVerZ_N_C, hVerZ_N_R, "Vz (cm)", -1.0, 9.0);
    cVerZ->SaveAs("plots/EffCorrData_VerZ_PosNeg3.png");

    //std::cout << "Raw entries: " << hInvMassRaw->GetEntries() << ", Max bin: " << hInvMassRaw->GetMaximum() << std::endl;
    //std::cout << "Corr entries: " << hInvMassCorr->GetEntries() << ", Max bin: " << hInvMassCorr->GetMaximum() << std::endl;
    // --- Plotting Final Mass ---
    TCanvas* cMass = new TCanvas("cMass", "Final Mass Distribution", 1000, 1000);
    PlotInvMass(hInvMassCorr,hInvMassRaw, "m_{K_{S}^{0}K_{S}^{0}} (GeV/c^{2})", 1.0, 2.5, -2.0, 15.0);
    //hInvMassCorr->SetTitle("Corrected");
    /*TH1D *hInvMassRaw_scaled = (TH1D*)hInvMassRaw->Clone("hInvMassRaw_scaled");
    //scaling the raw histogram to match the corrected one for better visual comparison
    double scaleFactor = 5.0; // Adjust to make it visible (204/4 ≈ 50)
    hInvMassRaw_scaled->Scale(scaleFactor);
    hInvMassRaw_scaled->SetLineColor(kGreen);
    hInvMassRaw_scaled->SetLineWidth(3);
    //hInvMassRaw_scaled->Draw("HIST same");
    hInvMassRaw->SetLineColor(kGreen);
    hInvMassRaw->SetLineWidth(3);
    hInvMassRaw->SetMarkerStyle(21);
    hInvMassRaw->SetMarkerColor(kGreen);
    hInvMassRaw->Draw("PE same");
    hInvMassCorr->Draw("PE same");
     
    TLegend* leg = new TLegend(0.55, 0.75, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.045);
    leg->AddEntry(hInvMassCorr, "Eff Corrected Data", "pe");
    //leg->AddEntry(fitFunc, "Fit: gaus + pol1", "l");
    //leg->AddEntry(hInvMassRaw_scaled, Form("Scaled Raw (%.1f)", scaleFactor), "l");
    leg->AddEntry(hInvMassRaw, "Raw Data", "l");*/
    //leg->Draw();
    cMass->SaveAs("plots/EffCorrData_K0K0Mass3.png");

    TCanvas* cMassRaw = new TCanvas("cMassRaw", "Raw Mass Distribution", 1000, 1000);
    hInvMassRaw->SetLineColor(kGreen);
    hInvMassRaw->SetLineWidth(3);
    hInvMassRaw->SetMarkerStyle(21);
    hInvMassRaw->SetMarkerColor(kGreen);
    hInvMassRaw->Draw("PE");
    cMassRaw->SaveAs("plots/EffCorrData_Raw_K0K0Mass3.png");

    TCanvas* c4TOF = new TCanvas("c4TOF", "4 TOF Hits Comparison", 1000, 1000);
    h1D_Reco_InvMass_Corrected_4pi_4TOF->SetTitle("4 TOF Hits");
    PlotInvMass(h1D_Reco_InvMass_Corrected_4pi_4TOF, h1D_Reco_InvMass_Raw_4pi_4TOF, "m_{K_{S}^{0}K_{S}^{0}} (GeV/c^{2})", 1.1, 2.5, -2.0, 10.0);
    c4TOF->SaveAs("plots/EffCorrData_4TOF.png");

     TCanvas* c3TOF = new TCanvas("c3TOF", "3 TOF Hits Comparison", 1000, 1000);
    h1D_Reco_InvMass_Corrected_4pi_3TOF->SetTitle("3 TOF Hits");
    PlotInvMass(h1D_Reco_InvMass_Corrected_4pi_3TOF, h1D_Reco_InvMass_Raw_4pi_3TOF, "m_{K_{S}^{0}K_{S}^{0}} (GeV/c^{2})", 1.1, 2.5, -2.0, 10.0);
    c3TOF->SaveAs("plots/EffCorrData_3TOF.png");

    TCanvas* c2TOF = new TCanvas("c2TOF", "2 TOF Hits Comparison", 1000, 1000);
    h1D_Reco_InvMass_Corrected_4pi_2TOF->SetTitle("2 TOF Hits");
    PlotInvMass(h1D_Reco_InvMass_Corrected_4pi_2TOF, h1D_Reco_InvMass_Raw_4pi_2TOF, "m_{K_{S}^{0}K_{S}^{0}} (GeV/c^{2})", 1.1, 2.5, -2.0, 100.0);
    c2TOF->SaveAs("plots/EffCorrData_2TOF.png");

}


