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
#include <iostream>
#include <cmath>
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"


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

void EvaluatePtMissSystematics(TH1D* hPtMiss) {
    // ---------------------------------------------------------
    // 1. Define Regions
    // ---------------------------------------------------------
    double sigMin = 0.0;
    double sigMax = 0.15;
    double binWidth = hPtMiss->GetBinWidth(1);

    // Calculate total events in the signal region
    int binMin = hPtMiss->FindBin(sigMin);
    int binMax = hPtMiss->FindBin(sigMax - 0.0001); // Avoid double-counting the boundary
    double totalEvents = hPtMiss->Integral(binMin, binMax);

    std::cout << "======================================================\n";
    std::cout << "Total Events in Signal Region (0 - 0.15): " << totalEvents << "\n";
    std::cout << "======================================================\n";

    // Create a canvas to draw all the fits
    TCanvas* cSys = new TCanvas("cSys", "pTmiss Systematics", 1000, 800);
    hPtMiss->SetMarkerStyle(20);
    hPtMiss->SetLineColor(kBlack);
    hPtMiss->GetXaxis()->SetLabelSize(0.05);
    hPtMiss->Draw("E1");

    TLegend* leg = new TLegend(0.11, 0.71, 0.48, 0.88);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.045);

    // ---------------------------------------------------------
    // 2. NOMINAL FIT: Rayleigh-like function, Range [0.20, 0.80]
    // ---------------------------------------------------------
    TF1* fitNominal = new TF1("fitNominal", "[0]*x*exp(-[1]*x)", 0.20, 0.80);
    fitNominal->SetParameters(10000, 5); // Starting guesses
    fitNominal->SetLineColor(kRed);
    fitNominal->SetLineWidth(3);
    
    hPtMiss->Fit(fitNominal, "RQ"); // 'R' = use range, 'Q' = quiet mode
    
    double bgNominal = fitNominal->Integral(sigMin, sigMax) / binWidth;
    double sigNominal = totalEvents - bgNominal;
    
    leg->AddEntry(fitNominal, Form("Nominal: Sig = %.1f", sigNominal), "l");

    // ---------------------------------------------------------
    // 3. SYS 1: Rayleigh-like function, Shifted Range High [0.25, 0.90]
    // ---------------------------------------------------------
    TF1* fitSys1 = new TF1("fitSys1", "[0]*x*exp(-[1]*x)", 0.25, 0.90);
    fitSys1->SetParameters(10000, 5);
    fitSys1->SetLineColor(kBlue);
    fitSys1->SetLineStyle(2); // Dashed line
    fitSys1->SetLineWidth(2);
    
    hPtMiss->Fit(fitSys1, "RQ+"); // '+' = add to existing drawing
    
    double bgSys1 = fitSys1->Integral(sigMin, sigMax) / binWidth;
    double sigSys1 = totalEvents - bgSys1;
    
    leg->AddEntry(fitSys1, Form("Sys1: Sig = %.1f", sigSys1), "l");

    // ---------------------------------------------------------
    // 4. SYS 2: Rayleigh-like function, Shifted Range Low [0.18, 0.70]
    // ---------------------------------------------------------
    TF1* fitSys2 = new TF1("fitSys2", "[0]*x*exp(-[1]*x)", 0.18, 0.70);
    fitSys2->SetParameters(10000, 5);
    fitSys2->SetLineColor(kGreen+2);
    fitSys2->SetLineStyle(3); // Dotted line
    fitSys2->SetLineWidth(2);
    
    hPtMiss->Fit(fitSys2, "RQ+");
    
    double bgSys2 = fitSys2->Integral(sigMin, sigMax) / binWidth;
    double sigSys2 = totalEvents - bgSys2;
    
    leg->AddEntry(fitSys2, Form("Sys2: Sig = %.1f", sigSys2), "l");

    // ---------------------------------------------------------
    // 5. SYS 3: Change Mathematical Function to a Modified Rayleigh
    // ---------------------------------------------------------
    // Adds a quadratic term to the exponent: x * exp(-B*x - C*x^2)
    TF1* fitSys3 = new TF1("fitSys3", "[0]*x*exp(-[1]*x - [2]*x*x)", 0.20, 0.80);
    fitSys3->SetParameters(10000, 5, 0.1); // Starting guesses
    fitSys3->SetLineColor(kMagenta);
    fitSys3->SetLineStyle(4); 
    fitSys3->SetLineWidth(2);
    
    hPtMiss->Fit(fitSys3, "RQ+");
    
    double bgSys3 = fitSys3->Integral(sigMin, sigMax) / binWidth;
    double sigSys3 = totalEvents - bgSys3;
    
    leg->AddEntry(fitSys3, Form("Sys3 (Mod. Rayleigh): Sig = %.1f", sigSys3), "l");

    // ---------------------------------------------------------
    // 6. Calculate Final Systematic Error
    // ---------------------------------------------------------
    double diff1 = std::abs(sigNominal - sigSys1);
    double diff2 = std::abs(sigNominal - sigSys2);
    double diff3 = std::abs(sigNominal - sigSys3);

    // Find the maximum deviation
    double maxSysError = std::max({diff1, diff2, diff3});

    leg->Draw();
    cSys->Update();
    cSys->SaveAs("plots/pTmiss_Systematics.png");

    std::cout << "\n=== BACKGROUND SUBTRACTION RESULTS ===\n";
    std::cout << "Nominal Yield:      " << sigNominal << " events\n";
    std::cout << "Sys 1 Yield (Diff): " << sigSys1 << " (" << diff1 << ")\n";
    std::cout << "Sys 2 Yield (Diff): " << sigSys2 << " (" << diff2 << ")\n";
    std::cout << "Sys 3 Yield (Diff): " << sigSys3 << " (" << diff3 << ")\n";
    std::cout << "------------------------------------------------------\n";
    std::cout << "FINAL REPORTED YIELD: " << sigNominal << " +/- " << maxSysError << " (syst.)\n";
    std::cout << "======================================================\n";
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

void ExtractDifferentialCrossSection(TH2D* h2) {

    int nBins = h2->GetNbinsX();

    // 2. Create the final output histograms
    // One for the signal yield (statistical errors), one for the systematic boxes
    TH1D* hFinalYield = (TH1D*)h2->ProjectionX("hFinalYield");
    hFinalYield->Reset(); // Clear the contents, keep the binning
    hFinalYield->SetTitle("Pure K_{S}^{0}K_{S}^{0} Yield;M(K_{S}^{0}K_{S}^{0}) [GeV/c^{2}];Events");

    // We use a TGraphErrors to draw the grey systematic boxes later
    TGraphErrors* grSystematics = new TGraphErrors(nBins);

    // 3. Define your fit functions (Modified Rayleigh)
    TF1* fitNominal = new TF1("fitNominal", "[0]*x*exp(-[1]*x - [2]*x*x)", 0.20, 0.80);
    TF1* fitSys1    = new TF1("fitSys1",    "[0]*x*exp(-[1]*x - [2]*x*x)", 0.25, 0.90);
    TF1* fitSys2    = new TF1("fitSys2",    "[0]*x*exp(-[1]*x - [2]*x*x)", 0.18, 0.70);

    // 4. THE LOOP
    for (int i = 1; i <= nBins; i++) {
        // Slice the 2D histogram to get the pTmiss for THIS mass bin only
        TH1D* hPtMiss = h2->ProjectionY(Form("ptMiss_bin_%d", i), i, i);
        double binWidth = hPtMiss->GetBinWidth(1);

        int binMin = hPtMiss->FindBin(0.0);
        int binMax = hPtMiss->FindBin(0.1499);
        double totalEvents = hPtMiss->Integral(binMin, binMax);

        // -> ADD THIS: Check how many events are actually in the tail!
        int tailMin = hPtMiss->FindBin(0.20);
        int tailMax = hPtMiss->FindBin(0.80);
        double tailEvents = hPtMiss->Integral(tailMin, tailMax);

        // --- NEW SAFETY CHECK ---
        // If there are very few events in the signal region, OR the tail is too empty to fit:
        if (totalEvents < 5 || tailEvents < 10) {
            // Set yield to 0 (or totalEvents if you prefer an upper limit)
            hFinalYield->SetBinContent(i, 0); 
            hFinalYield->SetBinError(i, 0);                        
            // Set systematic error to 0 to prevent giant grey boxes
            grSystematics->SetPoint(i-1, hFinalYield->GetBinCenter(i), 0);
            grSystematics->SetPointError(i-1, hFinalYield->GetBinWidth(i)/2.0, 0);
            
            std::cout << "Bin " << i << " skipped due to low stats.\n";
            continue; // Skip the fitting entirely!
        }

        // --- Fit Nominal ---
        fitNominal->SetParameters(1000, 5, 0.1); 
        // Optional: Stop the fitter from guessing insane negative slopes
        fitNominal->SetParLimits(1, 0.0, 50.0); 
        hPtMiss->Fit(fitNominal, "RQ0"); 
        
        double bgNominal = fitNominal->Integral(0.0, 0.15) / binWidth;
        if (bgNominal > totalEvents) bgNominal = totalEvents; // THE REALITY CAP
        if (bgNominal < 0) bgNominal = 0; // Prevent negative background guesses
        double sigNominal = totalEvents - bgNominal;

        // --- Fit Sys 1 ---
        fitSys1->SetParameters(1000, 5, 0.1);
        fitSys1->SetParLimits(1, 0.0, 50.0);
        hPtMiss->Fit(fitSys1, "RQ0");
        
        double bgSys1 = fitSys1->Integral(0.0, 0.15) / binWidth;
        if (bgSys1 > totalEvents) bgSys1 = totalEvents; // THE REALITY CAP
        if (bgSys1 < 0) bgSys1 = 0;
        double sigSys1 = totalEvents - bgSys1;

        // --- Fit Sys 2 ---
        fitSys2->SetParameters(1000, 5, 0.1);
        fitSys2->SetParLimits(1, 0.0, 50.0);
        hPtMiss->Fit(fitSys2, "RQ0");
        
        double bgSys2 = fitSys2->Integral(0.0, 0.15) / binWidth;
        if (bgSys2 > totalEvents) bgSys2 = totalEvents; // THE REALITY CAP
        if (bgSys2 < 0) bgSys2 = 0;
        double sigSys2 = totalEvents - bgSys2;

        // --- Calculate Systematics ---
        double diff1 = std::abs(sigNominal - sigSys1);
        double diff2 = std::abs(sigNominal - sigSys2);
        double maxSysError = std::max(diff1, diff2);

        // Prevent negative yields in pure noise bins
        if (sigNominal < 0) sigNominal = 0; 

        // 5. STORE THE RESULTS
        // Set the statistical data point
        hFinalYield->SetBinContent(i, sigNominal);
        
        // The statistical error of the background subtraction is roughly sqrt(Total + Bg)
        double statError = sqrt(totalEvents + bgNominal); 
        hFinalYield->SetBinError(i, statError);

        // Set the systematic grey box
        double massCenter = hFinalYield->GetBinCenter(i);
        double massWidth  = hFinalYield->GetBinWidth(i) / 2.0;
        
        grSystematics->SetPoint(i-1, massCenter, sigNominal);
        // X-error is the bin width, Y-error is your systematic error
        grSystematics->SetPointError(i-1, massWidth, maxSysError); 

        std::cout << "Bin " << i << " (Mass " << massCenter << "): Yield = " 
                  << sigNominal << " +/- " << statError << " (stat) +/- " 
                  << maxSysError << " (syst)\n";
    }

    gStyle->SetOptTitle(0); // Removes the default "Graph" title box at the top

    // 6. DRAW THE FINAL PLOT
    TCanvas* cFinal = new TCanvas("cFinal", "Differential Cross Section", 800, 600);
    
    // Draw the grey systematic boxes first
    grSystematics->SetFillColor(kGray);
    //set axis titles and ranges for the graph (since it's just boxes, it won't auto-set them)
    grSystematics->GetXaxis()->SetTitle("m_{K_{S}^{0}K_{S}^{0}} [GeV/c^{2}]");
    grSystematics->GetYaxis()->SetTitle("Corrected Yield / Bin");
    //grSystematics->GetXaxis()->SetLimits(0.9, 3.0);
    //grSystematics->GetYaxis()->SetRangeUser(0, hFinalYield->GetMaximum() * 1.5);
    grSystematics->Draw("A2"); // 'A2' draws boxes
    
    // Draw the statistical points on top
    hFinalYield->SetMarkerStyle(20);
    hFinalYield->SetLineColor(kBlack);
    hFinalYield->SetMinimum(0.0);
    hFinalYield->SetTitle("Central Exclusive Production: K_{S}^{0}K_{S}^{0}; m_{K_{S}^{0}K_{S}^{0}} [GeV/c^{2}]; Corrected Yield / Bin");

    hFinalYield->Draw("E1 SAME");
    
    cFinal->SaveAs("plots/Final_K0K0_Spectrum.png");
}


void Plot_EffCorrections()
{
    string DataFile = "EffCorr_nHits20_AorB_pTmiss3.root"; //EffCorrWithOldEffSameAsResultsApril14 EffCorrWithNewEffFile
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

    TH1D *h1D_PtMiss_Raw = (TH1D*)fcorr->Get("hPtMiss_Raw");
    TH1D *h1D_PtMiss_Corrected = (TH1D*)fcorr->Get("hPtMiss_Corrected");
    TH2D *h2_PtMiss_Vs_Mass = (TH2D*)fcorr->Get("h2_PtMiss_Vs_Mass");
    TH2D *h2_PtMiss_Vs_Mass_Corrected = (TH2D*)fcorr->Get("h2_PtMiss_Vs_Mass_Corrected");

    ExtractDifferentialCrossSection(h2_PtMiss_Vs_Mass_Corrected);

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

    

    TCanvas* cPtMiss = new TCanvas("cPtMiss", "pTmiss Comparison", 1000, 1000);
    For1atatime(h1D_PtMiss_Corrected, h1D_PtMiss_Raw, "p_{T}^{miss} (GeV/c)", -1.0, 100.0);
    cPtMiss->SaveAs("plots/EffCorrData_pTmiss3.png");

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

    EvaluatePtMissSystematics(h1D_PtMiss_Corrected);

    gStyle->SetPaintTextFormat(".0f"); // Set text format for 2 decimal places in the 2D histogram
   
    TCanvas* c2D = new TCanvas("c2D", "pTmiss vs Mass", 2000, 1000);
    c2D->Divide(2);
    c2D->cd(1);
    h2_PtMiss_Vs_Mass->SetTitle("Raw p_{T}^{miss} vs Mass");
    h2_PtMiss_Vs_Mass->Draw("TEXT COLZ");
    c2D->cd(2);
    h2_PtMiss_Vs_Mass_Corrected->SetTitle("Corrected p_{T}^{miss} vs Mass");
    h2_PtMiss_Vs_Mass_Corrected->Draw("TEXT COLZ");
    c2D->SaveAs("plots/pTmiss_vs_Mass.png");

}


