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


void SetHistogramStyle(TH1* hist, const char* title, 
    int markerStyle, int markerColor, float markerSize) {
    hist->GetXaxis()->SetTitleOffset(0.9);
    hist->GetYaxis()->SetTitleOffset(0.93);
    hist->GetXaxis()->SetTitleSize(0.07);
    hist->GetYaxis()->SetTitleSize(0.07);
    hist->GetXaxis()->SetLabelSize(0.05);
    hist->GetYaxis()->SetLabelSize(0.06);
    hist->GetXaxis()->SetNdivisions(6);
    hist->GetYaxis()->SetNdivisions(6);   
    //hist->GetYaxis()->SetRangeUser(10, 2000);
    int b_max = hist->GetMaximumBin();
    Double_t x_max = hist->GetBinCenter(b_max);
    Double_t y_max = hist->GetBinContent(b_max); 
    hist->GetYaxis()->SetRangeUser(x_max-200, y_max*1.2);
    
    hist->SetMarkerStyle(markerStyle);
    hist->SetMarkerColor(markerColor);
    hist->SetMarkerSize(markerSize);
    hist->SetLineColor(markerColor);

    if(title) hist->SetTitle(title);
}

void AxisSetting(TH1* hist, double labelsize, double titlesize, int div) {
    hist->GetXaxis()->SetTitleOffset(0.9);
    hist->GetYaxis()->SetTitleOffset(0.9);
    hist->GetXaxis()->SetTitleSize(titlesize);
    hist->GetYaxis()->SetTitleSize(titlesize);
    hist->GetXaxis()->SetLabelSize(labelsize);
    hist->GetYaxis()->SetLabelSize(labelsize);
    hist->GetXaxis()->SetNdivisions(div);
    hist->GetYaxis()->SetNdivisions(div);
}

void Combined_efficiency() {
    TFile* rootFile = TFile::Open("InputData/tofeff_withRP_0.3_v3.root");//tofeff_thesis
    if (!rootFile || rootFile->IsZombie()) {
        std::cerr << "Error opening input file" << std::endl;
        return;
    }

    // Get histograms for pT analysis
    TH1* HistPionPtWithTofMC = (TH1D*)rootFile->Get("HistPionPtWithTofMC;1");
    TH1* HistPionPtWithoutTofMC = (TH1D*)rootFile->Get("HistPionPtWithoutTofMC;1");
    TH1* HistTruePionPt = (TH1D*)rootFile->Get("HistTruePionPt;1");

    // Get histograms for eta analysis
    TH1* HistPionEtaWithTofMC = (TH1D*)rootFile->Get("HistPionEtaWithTofMC;1");
    TH1* HistPionEtaWithoutTofMC = (TH1D*)rootFile->Get("HistPionEtaWithoutTofMC;1"); 
    TH1* HistTruePionEta = (TH1D*)rootFile->Get("HistTruePionEta;1"); 

    // Get other histograms that were in the original code
    TH1* HistKaonMassProbeWithTof_TruePions = (TH1D*)rootFile->Get("HistKaonMassProbeWithTof_TruePions;1");
    TH1* HistKaonMassProbeWithoutTof_TruePions = (TH1D*)rootFile->Get("HistKaonMassProbeWithoutTof_TruePions;1");
    TH1* HistPionPtWithTof_TruePions = (TH1D*)rootFile->Get("HistPionPtWithTof_TruePions;1");
    TH1* HistPionPtWithoutTof_TruePions = (TH1D*)rootFile->Get("HistPionPtWithoutTof_TruePions;1");
    TH1* HistPionEtaWithTof_TruePions = (TH1D*)rootFile->Get("HistPionEtaWithTof_TruePions;1");
    TH1* HistPionEtaWithoutTof_TruePions = (TH1D*)rootFile->Get("HistPionEtaWithoutTof_TruePions;1");

    // --- pT ANALYSIS ---
    
    // Calculate efficiency for pT
    TH1D* HistTrueTofEfficiencyPt = (TH1D*)HistPionPtWithTofMC->Clone("HistTrueTofEfficiencyPt");
    //HistTrueTofEfficiencyPt->Divide(HistTruePionPt);  // Efficiency = (TOF-matched / All true pions)
    HistTrueTofEfficiencyPt->Divide(HistPionPtWithTofMC, HistTruePionPt, 1.0, 1.0, "B");                 
    // Create a histogram with (1 - efficiency) for pT
    TH1D* HistTrueTofInefficiencyPt = (TH1D*)HistTruePionPt->Clone("HistTrueTofInefficiencyPt");
    HistTrueTofInefficiencyPt->Add(HistPionPtWithTofMC, -1.0);
    HistTrueTofInefficiencyPt->Divide(HistTruePionPt);

    // Correct non-TOF pions using inefficiency for pT
    TH1D* HistCorrectedNonTOF_Pt = (TH1D*)HistPionPtWithoutTofMC->Clone("HistCorrectedNonTOF_Pt");
    HistCorrectedNonTOF_Pt->Divide(HistTrueTofInefficiencyPt);

    // Create a histogram with sum of TOF and non-TOF matched pions for pT
    TH1D* HistObservedTotal_Pt = (TH1D*)HistPionPtWithTofMC->Clone("HistObservedTotal_Pt");
    HistObservedTotal_Pt->Add(HistPionPtWithoutTofMC);

    // Keep TOF-matched pions as they are for pT (they're already detected)
    TH1D* HistCorrectedTOF_Pt = (TH1D*)HistPionPtWithTofMC->Clone("HistCorrectedTOF_Pt");

    // Create corrected total for pT (should match true total if everything works)
    TH1D* HistCorrectedTotal_Pt = (TH1D*)HistCorrectedTOF_Pt->Clone("HistCorrectedTotal_Pt");
    HistCorrectedTotal_Pt->Add(HistCorrectedNonTOF_Pt);

    // --- ETA ANALYSIS  ---
    
    // Calculate efficiency for eta
    TH1D* HistTrueTofEfficiencyEta = (TH1D*)HistPionEtaWithTofMC->Clone("HistTrueTofEfficiencyEta");
    //HistTrueTofEfficiencyEta->Divide(HistTruePionEta);  // Efficiency = (TOF-matched / All true pions)
    HistTrueTofEfficiencyEta->Divide(HistPionEtaWithTofMC, HistTruePionEta, 1.0, 1.0, "B");                 
    // Create a histogram with (1 - efficiency) for eta
    TH1D* HistTrueTofInefficiencyEta = (TH1D*)HistTruePionEta->Clone("HistTrueTofInefficiencyEta");
    HistTrueTofInefficiencyEta->Add(HistPionEtaWithTofMC, -1.0);
    HistTrueTofInefficiencyEta->Divide(HistTruePionEta);

    // Correct non-TOF pions using inefficiency for eta
    TH1D* HistCorrectedNonTOF_Eta = (TH1D*)HistPionEtaWithoutTofMC->Clone("HistCorrectedNonTOF_Eta");
    HistCorrectedNonTOF_Eta->Divide(HistTrueTofInefficiencyEta);

    // Create a histogram with sum of TOF and non-TOF matched pions for eta
    TH1D* HistObservedTotal_Eta = (TH1D*)HistPionEtaWithTofMC->Clone("HistObservedTotal_Eta");
    HistObservedTotal_Eta->Add(HistPionEtaWithoutTofMC);

    // Keep TOF-matched pions as they are for eta (they're already detected)
    TH1D* HistCorrectedTOF_Eta = (TH1D*)HistPionEtaWithTofMC->Clone("HistCorrectedTOF_Eta");

    // Create corrected total for eta (should match true total if everything works)
    TH1D* HistCorrectedTotal_Eta = (TH1D*)HistCorrectedTOF_Eta->Clone("HistCorrectedTotal_Eta");
    HistCorrectedTotal_Eta->Add(HistCorrectedNonTOF_Eta);

    // Set style options
    gStyle->SetStatY(0.990);
    gStyle->SetStatX(0.990);
    gStyle->SetStatW(0.15);
    gStyle->SetStatH(0.18);
    gStyle->SetPalette(kBird);
    gStyle->SetOptStat(0);

    // --- DRAW PT PLOTS ---
    TCanvas* c1 = new TCanvas("c1", "TOF Efficiency vs pT", 1800, 1300);

    // Create upper pad for main histograms
    TPad* pad1_pt = new TPad("pad1_pt", "pad1_pt", 0, 0.3, 1, 1.0);
    pad1_pt->SetLeftMargin(0.15);
    pad1_pt->SetRightMargin(0.10);
    pad1_pt->SetTopMargin(0.10);
    pad1_pt->SetBottomMargin(0.0); // No bottom margin for upper pad
    pad1_pt->Draw();
    
    // Create lower pad for efficiency
    TPad* pad2_pt = new TPad("pad2_pt", "pad2_pt", 0, 0.05, 1, 0.3);
    pad2_pt->SetLeftMargin(0.15);
    pad2_pt->SetRightMargin(0.10);
    pad2_pt->SetTopMargin(0.0); // No top margin for lower pad
    pad2_pt->SetBottomMargin(0.38);
    pad2_pt->Draw();
    
    // Draw main histograms in upper pad
    pad1_pt->cd();
    SetHistogramStyle(HistTruePionPt, nullptr, 20, kGreen+3, 3.0);    
    SetHistogramStyle(HistPionPtWithTofMC, nullptr, 21, kRed-3, 3.0);  
    SetHistogramStyle(HistCorrectedTotal_Pt, nullptr, 22, kBlue-3, 4.0);
    
    HistTruePionPt->GetXaxis()->SetRangeUser(0.2, 1.4);
    HistTruePionPt->GetXaxis()->SetLabelSize(0); // Hide x-axis labels
    HistTruePionPt->GetXaxis()->SetTitle("");    // Hide x-axis title
    HistTruePionPt->GetYaxis()->SetTitle("counts"); 
    HistTruePionPt->Draw("PE"); 
    HistPionPtWithTofMC->Draw("PE same"); 
    HistCorrectedTotal_Pt->Draw("PE same");
    
    TLegend* legend_pt = new TLegend(0.72, 0.65, 0.9, 0.85);
    legend_pt->SetBorderSize(0);
    legend_pt->SetFillStyle(0);
    legend_pt->AddEntry(HistTruePionPt, "True", "p");
    legend_pt->AddEntry(HistPionPtWithTofMC, "ToF matched", "p");
    legend_pt->AddEntry(HistCorrectedTotal_Pt, "Corrected Total", "p");
    legend_pt->Draw();
    
    // Add a title to the plot
    TLatex* title_pt = new TLatex(0.15, 0.95, "TOF Efficiency vs p_{T}");
    title_pt->SetNDC();
    title_pt->SetTextSize(0.045);
    title_pt->Draw();
    
    // Draw efficiency in lower pad
    pad2_pt->cd();
    SetHistogramStyle(HistTrueTofEfficiencyPt, nullptr, 21, kBlack-3, 3.0);  
    
    HistTrueTofEfficiencyPt->GetXaxis()->SetRangeUser(0.2, 1.4);
    HistTrueTofEfficiencyPt->GetXaxis()->SetTitle("p_{T} [GeV]");
    HistTrueTofEfficiencyPt->GetYaxis()->SetTitle("Efficiency");
    HistTrueTofEfficiencyPt->GetYaxis()->SetRangeUser(0, 1.1); // Adjust for percentage
    
    AxisSetting(HistTrueTofEfficiencyPt, 0.15, 0.19, 3);
    // Draw efficiency and reference line
    HistTrueTofEfficiencyPt->Draw("PE"); 
    TLine* line_pt = new TLine(0.2, 1.0, 1.4, 1.0); // 100% reference line
    line_pt->SetLineColor(kRed);
    line_pt->SetLineStyle(2);
    //line_pt->Draw();
    
    c1->Print("plots/TofEff_pT_analysis.png");

    // --- DRAW ETA PLOTS ---
    TCanvas* c2 = new TCanvas("c2", "TOF Efficiency vs Eta", 1800, 1300);

    // Create upper pad for main histograms
    TPad* pad1_eta = new TPad("pad1_eta", "pad1_eta", 0, 0.3, 1, 1.0);
    pad1_eta->SetLeftMargin(0.15);
    pad1_eta->SetRightMargin(0.10);
    pad1_eta->SetTopMargin(0.10);
    pad1_eta->SetBottomMargin(0.0); // No bottom margin for upper pad
    pad1_eta->Draw();
    
    // Create lower pad for efficiency
    TPad* pad2_eta = new TPad("pad2_eta", "pad2_eta", 0, 0.05, 1, 0.3);
    pad2_eta->SetLeftMargin(0.15);
    pad2_eta->SetRightMargin(0.10);
    pad2_eta->SetTopMargin(0.0); // No top margin for lower pad
    pad2_eta->SetBottomMargin(0.38);
    pad2_eta->Draw();
    
    // Draw main histograms in upper pad
    pad1_eta->cd();
    SetHistogramStyle(HistTruePionEta, nullptr, 20, kGreen+3, 3.0);    
    SetHistogramStyle(HistPionEtaWithTofMC, nullptr, 21, kRed-3, 3.0);  
    SetHistogramStyle(HistCorrectedTotal_Eta, nullptr, 22, kBlue-3, 4.0);
    
    // Set reasonable eta range (typical range is -2.5 to 2.5)
    HistTruePionEta->GetXaxis()->SetRangeUser(-1.5, 1.5); 
    HistTruePionEta->GetYaxis()->SetRangeUser(-100, 2500); // Adjust for counts
    HistTruePionEta->GetXaxis()->SetLabelSize(0); // Hide x-axis labels
    HistTruePionEta->GetXaxis()->SetTitle("");    // Hide x-axis title
    HistTruePionEta->GetYaxis()->SetTitle("counts"); 
    HistTruePionEta->Draw("PE"); 
    HistPionEtaWithTofMC->Draw("PE same"); 
    HistCorrectedTotal_Eta->Draw("PE same");
    
    TLegend* legend_eta = new TLegend(0.72, 0.65, 0.9, 0.85);
    legend_eta->SetBorderSize(0);
    legend_eta->SetFillStyle(0);
    legend_eta->AddEntry(HistTruePionEta, "True", "p");
    legend_eta->AddEntry(HistPionEtaWithTofMC, "ToF matched", "p");
    legend_eta->AddEntry(HistCorrectedTotal_Eta, "Corrected Total", "p");
    legend_eta->Draw();
    
    // Add a title to the plot
    TLatex* title_eta = new TLatex(0.15, 0.95, "TOF Efficiency vs #eta");
    title_eta->SetNDC();
    title_eta->SetTextSize(0.045);
    title_eta->Draw();
    
    // Draw efficiency in lower pad
    pad2_eta->cd();
    SetHistogramStyle(HistTrueTofEfficiencyEta, nullptr, 21, kBlack-3, 3.0);  
    
    HistTrueTofEfficiencyEta->GetXaxis()->SetRangeUser(-1.5, 1.5);
    HistTrueTofEfficiencyEta->GetXaxis()->SetTitle("#eta");
    HistTrueTofEfficiencyEta->GetYaxis()->SetTitle("Efficiency");
    HistTrueTofEfficiencyEta->GetYaxis()->SetRangeUser(0, 1.1); // Adjust for percentage
    
    AxisSetting(HistTrueTofEfficiencyEta, 0.15, 0.19, 3);
    // Draw efficiency and reference line
    HistTrueTofEfficiencyEta->Draw("PE"); 
    TLine* line_eta = new TLine(-1.5, 1.0, 1.5, 1.0); // 100% reference line
    line_eta->SetLineColor(kRed);
    line_eta->SetLineStyle(2);
   // line_eta->Draw();
    
    c2->Print("plots/TofEff_eta_analysis.png");

    // --- ADDITIONAL PLOT: RESIDUAL ANALYSIS ---
    // Create residual plots (how well correction recovers true distribution)
    
    // For pT
    TCanvas* c3 = new TCanvas("c3", "pT Residual Analysis", 1800, 600);
    TH1D* HistPtResidual = (TH1D*)HistCorrectedTotal_Pt->Clone("HistPtResidual");
    HistPtResidual->Add(HistTruePionPt, -1.0); // Subtract true from corrected
    HistPtResidual->Divide(HistTruePionPt);    // Divide by true to get relative residual
    HistPtResidual->Scale(100.0);              // Convert to percentage
    
    HistPtResidual->GetXaxis()->SetTitle("p_{T} [GeV]");
    HistPtResidual->GetYaxis()->SetTitle("Residual (Corrected-True)/True [%]");
    HistPtResidual->GetXaxis()->SetRangeUser(0.2, 1.4);
    HistPtResidual->GetYaxis()->SetRangeUser(-20, 20); // +/- 20% range
    HistPtResidual->SetMarkerStyle(21);
    HistPtResidual->SetMarkerColor(kBlue);
    HistPtResidual->SetMarkerSize(2.0);
    HistPtResidual->SetLineColor(kBlue);
    HistPtResidual->Draw("PE");
    
    // Add a zero line for reference
    TLine* zero_line_pt = new TLine(0.2, 0, 1.4, 0);
    zero_line_pt->SetLineColor(kRed);
    zero_line_pt->SetLineStyle(2);
    zero_line_pt->Draw();
    
    //c3->Print("plots/TofEff_pT_residual.png");
    
    // For eta
    TCanvas* c4 = new TCanvas("c4", "Eta Residual Analysis", 1800, 600);
    TH1D* HistEtaResidual = (TH1D*)HistCorrectedTotal_Eta->Clone("HistEtaResidual");
    HistEtaResidual->Add(HistTruePionEta, -1.0); // Subtract true from corrected
    HistEtaResidual->Divide(HistTruePionEta);    // Divide by true to get relative residual
    HistEtaResidual->Scale(100.0);                 // Convert to percentage
    
    HistEtaResidual->GetXaxis()->SetTitle("#eta");
    HistEtaResidual->GetYaxis()->SetTitle("Residual (Corrected-True)/True [%]");
    HistEtaResidual->GetXaxis()->SetRangeUser(-1.5, 1.5);
    HistEtaResidual->GetYaxis()->SetRangeUser(-20, 20); // +/- 20% range
    HistEtaResidual->SetMarkerStyle(21);
    HistEtaResidual->SetMarkerColor(kBlue);
    HistEtaResidual->SetMarkerSize(2.0);
    HistEtaResidual->SetLineColor(kBlue);
    HistEtaResidual->Draw("PE");
    
    // Add a zero line for reference
    TLine* zero_line_eta = new TLine(-1.5, 0, 1.5, 0);
    zero_line_eta->SetLineColor(kRed);
    zero_line_eta->SetLineStyle(2);
    zero_line_eta->Draw();
    
   // c4->Print("plots/TofEff_eta_residual.png");

    // --- ADDITIONAL PLOT: 2D EFFICIENCY MAP ---
    // Create a 2D histogram showing efficiency vs. pT and eta
    // This would require additional data or rebinning existing histograms
    
    // Close the file
    rootFile->Close();
}