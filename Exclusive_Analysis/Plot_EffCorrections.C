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
#include <algorithm>
#include <TGraphAsymmErrors.h>

void CalculateScaleFactorSystematic(double rawTailCounts, double estimatedBg, double netSignal, const char* label) {
    // 1. Calculate the Relative Poisson Error from the raw tail counts
    double poissonErr = (rawTailCounts > 0) ? std::sqrt(rawTailCounts) : 0.0;
    double relError = (rawTailCounts > 0) ? (poissonErr / rawTailCounts) : 0.0;
    
    // 2. Propagate to the Absolute Background Count
    double absBgSys = estimatedBg * relError;
    
    // 3. Compare to the Signal Statistical Error
    double sigStatErr = (netSignal > 0) ? std::sqrt(netSignal) : 0.0;
    double totalErr = std::sqrt((sigStatErr * sigStatErr) + (absBgSys * absBgSys));
    
    // Calculate the percentage increase on the final error bar
    double impactPercent = (sigStatErr > 0) ? ((totalErr - sigStatErr) / sigStatErr) * 100.0 : 0.0;

    // Print Formatted Table Entry
    std::cout << "--- Systematics Table Entry (" << label << ") ---\n";
    std::cout << "Systematic Source    : " << label << " Hybrid Scale Factor\n";
    std::cout << "Assigned Variation   : +/- " << std::fixed << std::setprecision(1) << (relError * 100.0) << "%\n";
    std::cout << "Absolute Yield Shift : +/- " << std::setprecision(2) << absBgSys << " events\n";
    
    // Round impact up to the nearest 0.1% for conservative reporting
    double displayImpact = (impactPercent < 0.1 && impactPercent > 0) ? 0.1 : std::ceil(impactPercent * 10.0) / 10.0;
    std::cout << "Impact on Total Error: < " << std::setprecision(1) << displayImpact << "%\n";
    std::cout << "======================================================\n\n";
}

TGraphAsymmErrors* CreateSysBand(TH1D* nom, TH1D* sys1, TH1D* sys2, Color_t color = kGray+1) {
    int nBins = nom->GetNbinsX();
    TGraphAsymmErrors* gr = new TGraphAsymmErrors(nBins);
    
    for (int i = 1; i <= nBins; ++i) {
        double x = nom->GetBinCenter(i);
        double y = nom->GetBinContent(i);
        double ex = nom->GetBinWidth(i) / 2.0;
        
        double val1 = sys1->GetBinContent(i);
        double val2 = sys2->GetBinContent(i);
        
        // Calculate maximum positive and negative deviations from nominal
        double errUp = std::max({0.0, val1 - y, val2 - y});
        double errDown = std::max({0.0, y - val1, y - val2});
        
        gr->SetPoint(i - 1, x, y);
        gr->SetPointError(i - 1, ex, ex, errDown, errUp);
    }
    
    // Set box visual styling
    gr->SetFillColorAlpha(color, 0.4); // 40% opaque for visibility
    gr->SetFillStyle(1001);            // Solid fill
    gr->SetLineWidth(0);               // No borders around the boxes
    
    return gr;
}

void sethist(TH1* h1, TH1* h2) {
    h1->SetLineColor(kRed);
    h1->SetMarkerStyle(20);
    h1->SetMarkerSize(1.2); 
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
    corrected->SetMarkerSize(3.0);
    raw->SetLineColor(kRed);
    raw->SetMarkerColor(kRed);
    raw->SetMarkerStyle(23);
    raw->SetLineWidth(2);
    raw->SetMarkerSize(3.0);
    corrected->GetYaxis()->SetMaxDigits(3);
    corrected->GetXaxis()->SetTitleOffset(1.5);

    TLegend *leg4 = new TLegend(0.69,0.69,0.88,0.88);
    leg4->SetBorderSize(0);
    leg4->SetFillStyle(0);
    leg4->SetTextSize(0.07);
    leg4->AddEntry(corrected, "Corrected", "lep");
    leg4->AddEntry(raw, "Raw", "lep");

    corrected->GetYaxis()->SetLabelSize(0.07);
    corrected->GetYaxis()->SetTitleSize(0.07);
    corrected->GetYaxis()->SetTitleOffset(0.8);

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
    TH1D *ratio = (TH1D*)corrected->Clone(Form("ratio_%s", corrected->GetName()));
    ratio->Divide(raw);
    ratio->SetLineColor(kBlue);
    ratio->SetMarkerColor(kBlue);
    ratio->SetMarkerStyle(22);
    ratio->SetLineWidth(2);
    ratio->SetMarkerSize(3.0);
    ratio->GetYaxis()->SetTitle("#frac{Corrected}{Raw}");
    ratio->GetYaxis()->SetLabelSize(0.15);
    ratio->GetYaxis()->SetTitleSize(0.19);
    ratio->GetYaxis()->SetTitleOffset(0.3);
    ratio->GetYaxis()->SetNdivisions(505);
    ratio->GetYaxis()->SetRangeUser(pullMinY1, pullMaxY2);
    ratio->GetXaxis()->SetTitle(xtitle);
    ratio->GetXaxis()->SetTitleSize(0.19);
    ratio->GetXaxis()->SetLabelSize(0.19);
    ratio->GetXaxis()->SetTitleOffset(0.9);
    ratio->SetTitle("");
    ratio->Draw("E");

    TLine *line = new TLine(ratio->GetXaxis()->GetXmin(), 0.0, ratio->GetXaxis()->GetXmax(), 0.0);
    line->SetLineColor(kRed);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw("SAME");

    TLine *line1sigma_up = new TLine(ratio->GetXaxis()->GetXmin(), 3.0, ratio->GetXaxis()->GetXmax(), 3.0);
    line1sigma_up->SetLineColor(kGray);
    line1sigma_up->SetLineStyle(2);
    line1sigma_up->Draw("SAME");

    TLine *line1sigma_down = new TLine(ratio->GetXaxis()->GetXmin(), -3.0, ratio->GetXaxis()->GetXmax(), -3.0);
    line1sigma_down->SetLineColor(kGray);
    line1sigma_down->SetLineStyle(2);
    //line1sigma_down->Draw("SAME");
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
    raw->GetXaxis()->SetTitleSize(0.06); 
    raw->GetYaxis()->SetTitleSize(0.06);
    raw->SetTitle("");
    raw->SetLineColor(kGreen);
    raw->SetLineWidth(3);
    raw->SetMarkerStyle(21);
    raw->SetMarkerColor(kGreen);
    corr->GetXaxis()->SetLabelSize(0.05);
    corr->GetYaxis()->SetLabelSize(0.05);
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

    TF1 *bgFunc = new TF1("bgFunc", "pol1(0)", 1.0, 2.4);
    bgFunc->SetParameter(0, fitFunc->GetParameter(3));
    bgFunc->SetParameter(1, fitFunc->GetParameter(4));
    bgFunc->SetLineStyle(2);
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
    
    TH1D *ratio = (TH1D*)corr->Clone(Form("ratioMass_%s", corr->GetName()));
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
    ratio->GetXaxis()->SetTitleSize(0.15);
    ratio->GetXaxis()->SetLabelSize(0.15);
    ratio->GetXaxis()->SetTitleOffset(1.0);
    ratio->SetTitle("");
    ratio->Draw("PE");

    TLine *line = new TLine(ratio->GetXaxis()->GetXmin(), 0.0, ratio->GetXaxis()->GetXmax(), 0.0);
    line->SetLineColor(kRed);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw("SAME");

    TLine *line1sigma_up = new TLine(ratio->GetXaxis()->GetXmin(), 3.0, ratio->GetXaxis()->GetXmax(), 3.0);
    line1sigma_up->SetLineColor(kGray);
    line1sigma_up->SetLineStyle(2);
    line1sigma_up->Draw("SAME");

    TLine *line1sigma_down = new TLine(ratio->GetXaxis()->GetXmin(), 5.0, ratio->GetXaxis()->GetXmax(), 5.0);
    line1sigma_down->SetLineColor(kGray);
    line1sigma_down->SetLineStyle(2);
    line1sigma_down->Draw("SAME");   
}

void PlotSys(TH1D* nom, TH1D* sys1, TH1D* sys2, const char* xtitle, const char* sys1Name, const char* sys2Name, double pullMin, double pullMax) {
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0.012);
    pad1->SetLeftMargin(0.15);
    pad1->SetRightMargin(0.05);
    pad1->Draw();
    pad1->cd();

    // Style the Nominal points
    nom->SetLineColor(kBlack);
    nom->SetMarkerColor(kBlack);
    nom->SetMarkerStyle(20);
    nom->SetLineWidth(2);
    nom->SetMarkerSize(2.0);

    // Style the variations (optional to draw them as points, but good for debugging)
    sys1->SetLineColor(kBlue);
    sys1->SetMarkerColor(kBlue);
    sys1->SetMarkerStyle(22);
    sys1->SetMarkerSize(2.0);
    sys2->SetLineColor(kRed);
    sys2->SetMarkerColor(kRed);
    sys2->SetMarkerStyle(23);
    sys2->SetMarkerSize(2.0);

    nom->GetListOfFunctions()->Clear();

    // Auto-scale max Y to fit all plots
    double maxY = std::max({nom->GetMaximum(), sys1->GetMaximum(), sys2->GetMaximum()}) * 1.2;
    nom->SetMaximum(maxY);
    nom->SetTitle("");
    nom->GetYaxis()->SetMaxDigits(3);
    nom->GetYaxis()->SetLabelSize(0.07);
    nom->GetYaxis()->SetTitleSize(0.06);
    double Axymin = std::min({nom->GetMinimum(), sys1->GetMinimum(), sys2->GetMinimum()});
    nom->GetYaxis()->SetRangeUser(Axymin * 1.5, nom->GetMaximum() * 1.2); //150);// 
    //nom->GetYaxis()->SetRangeUser(-2, nom->GetMaximum() * 1.2); //150);//

    // Draw an empty frame first using nominal
    nom->Draw("AXIS");

    // 1. Create and draw the Systematic Band
    TGraphAsymmErrors* sysBand = CreateSysBand(nom, sys1, sys2, kGray+2);
    sysBand->Draw("2 SAME"); // '2' draws filled rectangles

    // 2. Draw variations as markers on top (optional, can comment out if too cluttered)
    sys1->Draw("P SAME");
    sys2->Draw("P SAME");
    
    // 3. Draw Nominal data on top of the band
    nom->Draw("E SAME");

    double fitmin = 1.0;
    double fitmax = 2.6;
    /*TF1* fitFunc = new TF1("fitFunc", "gaus(0) + pol1(3)", fitmin, fitmax);
    fitFunc->SetParameters(15, 1.68, 0.08, 10, -2); 
    fitFunc->SetLineColor(kRed);
    nom->Fit(fitFunc, "R"); 

    TF1 *bgFunc = new TF1("bgFunc", "pol1(0)", 1.0, 2.4);
    bgFunc->SetParameter(0, fitFunc->GetParameter(3));
    bgFunc->SetParameter(1, fitFunc->GetParameter(4));
    bgFunc->SetLineStyle(2);
    bgFunc->SetLineColor(kBlue);
    bgFunc->Draw("SAME");

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.045);
    latex.DrawLatex(0.55, 0.58, Form("A = %.2f", fitFunc->GetParameter(0)));
    latex.DrawLatex(0.55, 0.52, Form("#mu = %.2f GeV/c^{2}", fitFunc->GetParameter(1)));
    latex.DrawLatex(0.55, 0.46, Form("#sigma = %.2f GeV/c^{2}", fitFunc->GetParameter(2)));*/
    

    TLegend *leg = new TLegend(0.55, 0.65, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.06);
    leg->AddEntry(nom, "Nominal (nHits 20)", "lep");
    leg->AddEntry(sys1, sys1Name, "p");
    leg->AddEntry(sys2, sys2Name, "p");
    leg->AddEntry(sysBand, "Syst. Uncertainty", "f");
    leg->Draw();

    gPad->GetMother()->cd();

    // LOWER PAD: Ratio
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
    pad2->SetTopMargin(0.018);
    pad2->SetBottomMargin(0.455);
    pad2->SetLeftMargin(0.15);
    pad2->SetRightMargin(0.05);
    pad2->Draw();
    pad2->cd();

    // Create ratios
    TH1D *ratioNom = (TH1D*)nom->Clone(Form("ratioNom_%s", nom->GetName()));
    ratioNom->Divide(nom); // Centered at 1.0

    TH1D *ratio1 = (TH1D*)sys1->Clone(Form("ratio1_%s", sys1->GetName()));
    ratio1->GetListOfFunctions()->Clear();
    ratio1->Divide(nom);
    ratio1->SetLineColor(kBlue);
    ratio1->SetMarkerColor(kBlue);
    ratio1->SetMarkerStyle(22);

    TH1D *ratio2 = (TH1D*)sys2->Clone(Form("ratio2_%s", sys2->GetName()));
    ratio2->GetListOfFunctions()->Clear();
    ratio2->Divide(nom);
    ratio2->SetLineColor(kRed);
    ratio2->SetMarkerColor(kRed);
    ratio2->SetMarkerStyle(23);

    // Format axes using ratio1 as the frame
    ratio1->GetYaxis()->SetTitle("#frac{Sys}{Nom}");
    ratio1->GetYaxis()->SetTitleSize(0.18);
    ratio1->GetYaxis()->SetLabelSize(0.12);
    ratio1->GetYaxis()->SetTitleOffset(0.38);
    ratio1->GetYaxis()->SetNdivisions(505);
    ratio1->GetYaxis()->SetRangeUser(pullMin, pullMax);
    ratio1->GetXaxis()->SetTitle(xtitle);
    ratio1->GetXaxis()->SetTitleSize(0.23);
    ratio1->GetXaxis()->SetLabelSize(0.18);
    ratio1->GetXaxis()->SetTitleOffset(0.78);
    ratio1->SetTitle("");
    
    // Draw empty ratio frame
    ratio1->Draw("AXIS");

    // Create and draw systematic ratio band
    TGraphAsymmErrors* ratioBand = CreateSysBand(ratioNom, ratio1, ratio2, kGray+2);
    ratioBand->Draw("2 SAME");

    // Draw reference line
    TLine *line = new TLine(ratio1->GetXaxis()->GetXmin(), 1.0, ratio1->GetXaxis()->GetXmax(), 1.0);
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);
    line->Draw("SAME");

    // Draw the actual ratio markers
    ratio1->Draw("P SAME");
    ratio2->Draw("P SAME");
}

void EvaluatePtMissSystematics2(TH1D* h1D, TH2D* h2D, const char* name) {
    TH1D* hPtMiss = nullptr;

    if (h1D != nullptr) {
        hPtMiss = h1D;
    }
    else if (h2D != nullptr) {
        hPtMiss = h2D->ProjectionX(Form("%s_projX", name));
    }
    else {
        std::cerr << "Error: No valid histogram provided for " << name << std::endl;
        return;
    }

    if (hPtMiss == nullptr || hPtMiss->GetEntries() == 0) {
        std::cerr << "Warning: hPtMiss is empty or invalid for " << name << std::endl;
        return;
    }

    // ---------------------------------------------------------
    // 1. Define Regions (Synchronized with 2D Shape Method)
    // ---------------------------------------------------------
    double sigMin  = 0.0;
    double sigMax  = 0.15; // Your K0K0 signal cut
    double tailMin = 0.20; // Start of background tail
    double tailMax = 0.80; // End of background tail
    double binWidth = hPtMiss->GetBinWidth(1);

    int binMin = hPtMiss->FindBin(sigMin);
    int binMax = hPtMiss->FindBin(sigMax - 0.0001);
    double totalEvents = hPtMiss->Integral(binMin, binMax);

    std::cout << "\n======================================================\n";
    std::cout << " GLOBAL 1D BACKGROUND FIT (" << name << ")\n";
    std::cout << "======================================================\n";
    std::cout << "Total Inclusive Events in Signal Region (" << sigMin << " - " << sigMax << "): " << totalEvents << "\n";

    TCanvas* cSys = new TCanvas(Form("cSys_%s", name), "Global pTmiss Systematics", 1000, 800);
    hPtMiss->SetMarkerStyle(20);
    hPtMiss->SetLineColor(kBlack);
    hPtMiss->GetXaxis()->SetLabelSize(0.05);
    hPtMiss->GetXaxis()->SetTitle("p_{T}^{miss} [GeV/c]");
    hPtMiss->GetYaxis()->SetTitle("Events");
    //draw a darker frame box
    gPad->SetFrameLineWidth(3);
    hPtMiss->Draw("E1");

    TLegend* leg = new TLegend(0.11, 0.71, 0.48, 0.88);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.045);

    // ---------------------------------------------------------
    // 2. NOMINAL FIT: Anchored pol2, Range [0.20, 0.80]
    // ---------------------------------------------------------
    TF1* fitNominal = new TF1("fitNominal_1D", "pol2", tailMin, tailMax);
    fitNominal->SetParameter(1, 10);
    fitNominal->SetParameter(2, 1);
    fitNominal->FixParameter(0, 0.0); // THE ANCHOR
    fitNominal->SetLineColor(kRed);
    fitNominal->SetLineWidth(3);
    
    hPtMiss->Fit(fitNominal, "RLQ"); 
    
    double bgNominal = fitNominal->Integral(sigMin, sigMax) / binWidth;
    if (bgNominal > totalEvents) bgNominal = totalEvents;
    if (bgNominal < 0) bgNominal = 0;
    double sigNominal = totalEvents - bgNominal;
    
    leg->AddEntry(fitNominal, Form("Nominal Fit: Sig = %.1f", sigNominal), "l");

    // ---------------------------------------------------------
    // 3. SYS 1: Anchored pol2, Shifted Range High [0.25, 0.85]
    // ---------------------------------------------------------
    TF1* fitSys1 = new TF1("fitSys1_1D", "pol2", tailMin + 0.05, tailMax + 0.05);
    fitSys1->SetParameter(1, 10);
    fitSys1->SetParameter(2, 1);
    fitSys1->FixParameter(0, 0.0);
    fitSys1->SetLineColor(kBlue);
    fitSys1->SetLineStyle(2); 
    fitSys1->SetLineWidth(2);
    
    hPtMiss->Fit(fitSys1, "RLQ+"); 
    
    double bgSys1 = fitSys1->Integral(sigMin, sigMax) / binWidth;
    if (bgSys1 > totalEvents) bgSys1 = totalEvents;
    if (bgSys1 < 0) bgSys1 = 0;
    double sigSys1 = totalEvents - bgSys1;
    
    leg->AddEntry(fitSys1, Form("Sys1 (Range High): Sig = %.1f", sigSys1), "l");

    // ---------------------------------------------------------
    // 4. SYS 2: Anchored pol2, Shifted Range Low [0.18, 0.75]
    // ---------------------------------------------------------
    TF1* fitSys2 = new TF1("fitSys2_1D", "pol2", tailMin - 0.02, tailMax - 0.05);
    fitSys2->SetParameter(1, 10);
    fitSys2->SetParameter(2, 1);
    fitSys2->FixParameter(0, 0.0);
    fitSys2->SetLineColor(kGreen+2);
    fitSys2->SetLineStyle(3); 
    fitSys2->SetLineWidth(2);
    
    hPtMiss->Fit(fitSys2, "RLQ+");
    
    double bgSys2 = fitSys2->Integral(sigMin, sigMax) / binWidth;
    if (bgSys2 > totalEvents) bgSys2 = totalEvents;
    if (bgSys2 < 0) bgSys2 = 0;
    double sigSys2 = totalEvents - bgSys2;
    
    leg->AddEntry(fitSys2, Form("Sys2 (Range Low): Sig = %.1f", sigSys2), "l");

    // ---------------------------------------------------------
    // 5. Calculate Final Systematic Error & Report
    // ---------------------------------------------------------
    double diff1 = std::abs(sigNominal - sigSys1);
    double diff2 = std::abs(sigNominal - sigSys2);
    double maxSysError = std::max(diff1, diff2);

    leg->Draw();
    cSys->Update();
    cSys->SaveAs(Form("plots/pTmiss_Systematics_%s_Global.png", name));

    std::cout << "\n--- GLOBAL BACKGROUND EXTRACTION --- \n";
    std::cout << "Nominal Global Background: " << bgNominal << " events\n";
    std::cout << "Nominal Global Signal:     " << sigNominal << " events\n";
    std::cout << "Max Systematic Error:      " << maxSysError << " events\n";
    std::cout << "======================================================\n";
}

TH1D* ExtractDifferentialCrossSection2(TH2D* h2, const char* nameSuffix) {
    std::cout << "\n======================================================\n";
    std::cout << " EXTRACTING DIFFERENTIAL CROSS SECTION (SHAPE METHOD) \n";
    std::cout << "======================================================\n";

    // ---------------------------------------------------------
    // 1. Define Kinematic Regions
    // ---------------------------------------------------------
    double sigMax = 0.15;  // Signal region limit (0.0 to 0.15 GeV/c)
    
    // The tail region used to extract the background mass shape and fit
    // Tomáš suggests this should be close to the signal region but outside it.
    double tailMin = 0.20; 
    double tailMax = 0.80; 

    // ---------------------------------------------------------
    // 2. Global Fit: Find Total Background in Signal Region
    // ---------------------------------------------------------
    // Project all masses to get the global 1D pTmiss distribution
    TH1D* hPtMissGlobal = h2->ProjectionY("hPtMiss_Global");
    double binWidthY = hPtMissGlobal->GetBinWidth(1);

    // Fit Nominal
    TF1* fitNominal = new TF1("fitNominal", "pol2", tailMin, tailMax);
    fitNominal->SetParameter(1, 10);
    fitNominal->SetParameter(2, 1);
    fitNominal->FixParameter(0, 0.0); // Anchored to 0
    hPtMissGlobal->Fit(fitNominal, "RLQ0"); // Q0 = Quiet, don't draw
    
    double totalBgNominal = fitNominal->Integral(0.0, sigMax) / binWidthY;
    if (totalBgNominal < 0) totalBgNominal = 0;

    // Fit Sys 1 (Shifted High)
    TF1* fitSys1 = new TF1("fitSys1", "pol2", tailMin + 0.05, tailMax + 0.05);
    fitSys1->SetParameter(1, 10);
    fitSys1->SetParameter(2, 1);
    fitSys1->FixParameter(0, 0.0);
    hPtMissGlobal->Fit(fitSys1, "RLQ0");
    double totalBgSys1 = std::max(0.0, fitSys1->Integral(0.0, sigMax) / binWidthY);

    // Fit Sys 2 (Shifted Low)
    TF1* fitSys2 = new TF1("fitSys2", "pol2", tailMin - 0.02, tailMax - 0.05);
    fitSys2->SetParameter(1, 10);
    fitSys2->SetParameter(2, 1);
    fitSys2->FixParameter(0, 0.0);
    hPtMissGlobal->Fit(fitSys2, "RLQ0");
    double totalBgSys2 = std::max(0.0, fitSys2->Integral(0.0, sigMax) / binWidthY);

    std::cout << "Global Background Yield in Signal Region: " << totalBgNominal << "\n";

    // ---------------------------------------------------------
    // 3. Extract 1D Mass Shapes
    // ---------------------------------------------------------
    // Signal + Background Shape (project X from Y bins [0, sigMax])
    int binSigMin = h2->GetYaxis()->FindBin(0.0);
    int binSigMax = h2->GetYaxis()->FindBin(sigMax - 0.0001);
    TH1D* hSigShape = h2->ProjectionX("hSigShape", binSigMin, binSigMax);

    // Background Shape (project X from Y bins [tailMin, tailMax])
    int binTailMin = h2->GetYaxis()->FindBin(tailMin);
    int binTailMax = h2->GetYaxis()->FindBin(tailMax - 0.0001);
    TH1D* hBgShape = h2->ProjectionX("hBgShape", binTailMin, binTailMax);

    // ---------------------------------------------------------
    // 4. Normalize Background Shape
    // ---------------------------------------------------------
    double tailShapeIntegral = hBgShape->Integral();

    TH1D* hBgNominal = (TH1D*)hBgShape->Clone("hBgNominal");
    TH1D* hBgSys1    = (TH1D*)hBgShape->Clone("hBgSys1");
    TH1D* hBgSys2    = (TH1D*)hBgShape->Clone("hBgSys2");

    if (tailShapeIntegral > 0) {
        hBgNominal->Scale(totalBgNominal / tailShapeIntegral);
        hBgSys1->Scale(totalBgSys1 / tailShapeIntegral);
        hBgSys2->Scale(totalBgSys2 / tailShapeIntegral);
    }
    /*
    // ---------------------------------------------------------
    // 5. Subtract Background & Calculate Systematics
    // ---------------------------------------------------------
    TH1D* hFinalYield = (TH1D*)hSigShape->Clone("hFinalYield");
    hFinalYield->Add(hBgNominal, -1.0); // Subtract Nominal Background
    hFinalYield->SetTitle("Central Exclusive Production: K_{S}^{0}K_{S}^{0}; m_{K_{S}^{0}K_{S}^{0}} [GeV/c^{2}]; Corrected Yield / Bin");

    int nBins = hFinalYield->GetNbinsX();
    TGraphErrors* grSystematics = new TGraphErrors(nBins);

    for (int i = 1; i <= nBins; i++) {
        double yieldNominal = hFinalYield->GetBinContent(i);
        if (yieldNominal < 0) {
            hFinalYield->SetBinContent(i, 0);
            yieldNominal = 0;
        }

        // Calculate yields for systematic shifts
        double yieldSys1 = hSigShape->GetBinContent(i) - hBgSys1->GetBinContent(i);
        double yieldSys2 = hSigShape->GetBinContent(i) - hBgSys2->GetBinContent(i);
        if (yieldSys1 < 0) yieldSys1 = 0;
        if (yieldSys2 < 0) yieldSys2 = 0;

        // Max deviation defines the systematic error
        double diff1 = std::abs(yieldNominal - yieldSys1);
        double diff2 = std::abs(yieldNominal - yieldSys2);
        double maxSysError = std::max(diff1, diff2);

        // Store Statistical Error (sqrt(Signal + Background))
        double statError = sqrt(hSigShape->GetBinContent(i)); 
        hFinalYield->SetBinError(i, statError);

        // Store Systematic Error (Grey Boxes)
        double massCenter = hFinalYield->GetBinCenter(i);
        double massWidth  = hFinalYield->GetBinWidth(i) / 2.0;
        grSystematics->SetPoint(i-1, massCenter, yieldNominal);
        grSystematics->SetPointError(i-1, massWidth, maxSysError);
    }

    // ---------------------------------------------------------
    // 6. Draw Final Plot
    // ---------------------------------------------------------
    
    TCanvas* cFinal = new TCanvas("cFinal", "Differential Cross Section", 800, 600);
    
    grSystematics->SetFillColor(kGray);
    grSystematics->GetXaxis()->SetTitle("m_{K_{S}^{0}K_{S}^{0}} [GeV/c^{2}]");
    grSystematics->GetYaxis()->SetTitle("Corrected Yield / Bin");
    grSystematics->GetYaxis()->SetRangeUser(-20, hFinalYield->GetMaximum() * 1.5);
    grSystematics->Draw("A2"); 
    
    hFinalYield->SetMarkerStyle(20);
    hFinalYield->SetLineColor(kBlack);
    hFinalYield->Draw("E1 SAME");

    TLegend* leg = new TLegend(0.45, 0.75, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.045);
    leg->AddEntry(hFinalYield, "Statistical Uncertainty", "lep");
    leg->AddEntry(grSystematics, "Systematic Uncertainty", "f");
    leg->Draw();
    
    cFinal->SaveAs("plots/Final_K0K0_Spectrum_ShapeMethod.png");*/
    // ---------------------------------------------------------
    // 5. Subtract Background & Calculate Systematics
    // ---------------------------------------------------------
    TH1D* hFinalYield = (TH1D*)hSigShape->Clone("hFinalYield");
    hFinalYield->Add(hBgNominal, -1.0); // Subtract Nominal Background
    hFinalYield->SetTitle("Central Exclusive Production: K_{S}^{0}K_{S}^{0}; m_{K_{S}^{0}K_{S}^{0}} [GeV/c^{2}]; Corrected Yield / Bin");

    int nBins = hFinalYield->GetNbinsX();
    TGraphErrors* grSystematics = new TGraphErrors(nBins);

    std::cout << "\n--- BIN-BY-BIN BACKGROUND QUALITY CHECK ---\n";
    std::cout << std::setw(15) << "Mass Bin (GeV)" 
              << std::setw(15) << "Total Events" 
              << std::setw(15) << "Background" 
              << std::setw(15) << "Subtracted" << "\n";
    std::cout << "--------------------------------------------------------------\n";

    for (int i = 1; i <= nBins; i++) {
        double yieldNominal = hFinalYield->GetBinContent(i);
        if (yieldNominal < 0) {
            hFinalYield->SetBinContent(i, 0);
            yieldNominal = 0;
        }

        // Calculate yields for systematic shifts
        double yieldSys1 = hSigShape->GetBinContent(i) - hBgSys1->GetBinContent(i);
        double yieldSys2 = hSigShape->GetBinContent(i) - hBgSys2->GetBinContent(i);
        if (yieldSys1 < 0) yieldSys1 = 0;
        if (yieldSys2 < 0) yieldSys2 = 0;

        // Max deviation defines the systematic error
        double diff1 = std::abs(yieldNominal - yieldSys1);
        double diff2 = std::abs(yieldNominal - yieldSys2);
        double maxSysError = std::max(diff1, diff2);

        // Store Statistical Error (sqrt(Total events in signal region))
        //double statError = sqrt(hSigShape->GetBinContent(i)); 
        //hFinalYield->SetBinError(i, statError);
        double statError = hFinalYield->GetBinError(i); 
        hFinalYield->SetBinError(i, statError);

        double massCenter = hFinalYield->GetBinCenter(i);
        double massWidth  = hFinalYield->GetBinWidth(i) / 2.0;
        grSystematics->SetPoint(i-1, massCenter, yieldNominal);
        grSystematics->SetPointError(i-1, massWidth, maxSysError);

        // --- PRINT TO TERMINAL FOR SENIOR ---
        // Only print for bins inside your actual mass plotting range (e.g., 0.9 to 3.0)
        if (massCenter >= 0.9 && massCenter <= 3.0) {
            std::cout << std::setw(15) << massCenter 
                      << std::setw(15) << hSigShape->GetBinContent(i) 
                      << std::setw(15) << hBgNominal->GetBinContent(i) 
                      << std::setw(15) << yieldNominal << "\n";
        }
    }
    std::cout << "--------------------------------------------------------------\n";
    

    // ---------------------------------------------------------
    // 6. Draw Final Plot 
    // ---------------------------------------------------------
   
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    TCanvas* cFinal = new TCanvas("cFinal", "Differential Cross Section", 800, 600);
    
    grSystematics->SetFillColor(kGray);
    grSystematics->GetXaxis()->SetTitle("m_{K_{S}^{0}K_{S}^{0}} [GeV/c^{2}]");
    grSystematics->GetYaxis()->SetTitle("Corrected Yield / Bin");
    grSystematics->GetYaxis()->SetRangeUser(-20, hFinalYield->GetMaximum() * 1.5);
    grSystematics->SetTitle("Final K_{S}^{0}K_{S}^{0} Yield with Systematics (Shape Method)");
    grSystematics->GetXaxis()->SetTitleSize(0.05);
    grSystematics->GetYaxis()->SetTitleSize(0.05);
    grSystematics->GetXaxis()->SetLabelSize(0.04);
    grSystematics->GetYaxis()->SetLabelSize(0.04);
    grSystematics->GetYaxis()->SetTitleOffset(0.8);
    grSystematics->GetXaxis()->SetTitleOffset(0.8);
    //draw a darker frame box
    gPad->SetFrameLineWidth(3);
    grSystematics->Draw("A2"); 
    
    hFinalYield->SetMarkerStyle(20);
    hFinalYield->SetLineColor(kBlack);
    hFinalYield->Draw("E1 SAME");

    TLegend* leg = new TLegend(0.45, 0.75, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.045);
    leg->AddEntry(hFinalYield, "Statistical Uncertainty", "lep");
    leg->AddEntry(grSystematics, "Systematic Uncertainty", "f");
    leg->Draw();
    
    //cFinal->SaveAs("plots/Final_K0K0_Spectrum_ShapeMethod.png");
    //cFinal->SaveAs(Form("plots/Final_K0K0_Spectrum_ShapeMethod_%s.png", nameSuffix));

    // ---------------------------------------------------------
    // 7. NEW: Quality Check Overlay Plot for your Senior
    // ---------------------------------------------------------
    TCanvas* cQuality = new TCanvas("cQuality", "Background Quality Check", 800, 600);
    
    // Total Events in Signal Region (Before Subtraction)
    hSigShape->SetLineColor(kBlue);
    hSigShape->SetMarkerColor(kBlue);
    hSigShape->SetMarkerStyle(21);
    hSigShape->GetXaxis()->SetRangeUser(0.9, 3.0); // Focus on interesting mass range
    hSigShape->GetYaxis()->SetRangeUser(0, hSigShape->GetMaximum() * 1.3);
    hSigShape->GetYaxis()->SetTitle("Events / Bin");
    hSigShape->GetXaxis()->SetTitle("m_{K_{S}^{0}K_{S}^{0}} [GeV/c^{2}]");
    hSigShape->SetTitle("Quality Check: Background Subtraction per Mass Bin");
    hSigShape->Draw("PE"); 

    // Estimated Background (Filled area to show what gets removed)
    hBgNominal->SetLineColor(kRed);
    hBgNominal->SetFillColorAlpha(kRed, 0.3); // Semi-transparent red
    hBgNominal->Draw("HIST SAME"); 

    // Subtracted Final Yield
    hFinalYield->Draw("PE SAME"); 

    TLegend* legQ = new TLegend(0.45, 0.70, 0.88, 0.88);
    legQ->SetBorderSize(0);
    legQ->AddEntry(hSigShape, Form("Total Events (pT^{miss} < %.2f)", sigMax), "pe");
    legQ->AddEntry(hBgNominal, "Estimated Background", "f");
    legQ->AddEntry(hFinalYield, "Final Subtracted Signal", "pe");
    legQ->Draw();

    //cQuality->SaveAs("plots/QualityCheck_BackgroundPerBin.png");

    hFinalYield->SetName(Form("hFinalYield_%s", nameSuffix));
    return hFinalYield;
}
/*
void CheckBinByBinPtMissFits(TH2D* h2) {
    std::cout << "\n======================================================\n";
    std::cout << " GENERATING BIN-BY-BIN pTmiss FITS (PDF REPORT) \n";
    std::cout << "======================================================\n";

    // Define the mass range we care about (0.9 to 3.0 GeV)
    int binMin = h2->GetXaxis()->FindBin(0.9);
    int binMax = h2->GetXaxis()->FindBin(3.0);
    
    // Limits for pTmiss
    double sigMax = 0.15;
    double tailMin = 0.20;
    double tailMax = 0.80;


    TCanvas* cFit = new TCanvas("cFit", "Bin-by-Bin Fits", 800, 600);
    
    // Open a multi-page PDF
    TString pdfName = "plots/BinByBin_PtMiss_Fits_Report.pdf";
    cFit->Print(pdfName + "["); // The "[" tells ROOT to open the file but not close it

    gStyle->SetOptStat(0);
    //gStyle->SetOptFit(1111); // Show fit parameters box

    for (int i = binMin; i <= binMax; i++) {
        double massCenter = h2->GetXaxis()->GetBinCenter(i);
        
        // Slice the 2D histogram for this specific X (Mass) bin
        TH1D* hSlice = h2->ProjectionY(Form("hSlice_%d", i), i, i);
        
        // Skip bins with very low statistics to avoid broken fits
        if (hSlice->GetEntries() < 5) continue;

        hSlice->SetTitle(Form("p_{T}^{miss} for Mass = %.2f GeV/c^{2}; p_{T}^{miss} [GeV/c]; Events", massCenter));
        hSlice->SetMarkerStyle(20);
        hSlice->SetLineColor(kBlack);
        hSlice->GetXaxis()->SetRangeUser(0.0, 1.0); // Zoom in on the relevant pTmiss region

        // Define and apply the background fit to the tail
        TF1* fitBg = new TF1(Form("fitBg_%d", i), "pol2", tailMin, tailMax);
        fitBg->SetParameter(1, 10);
        fitBg->SetParameter(2, 1);
        fitBg->FixParameter(0, 0.0); // Anchored to 0 just like your global fit
        fitBg->SetLineColor(kRed);
        fitBg->SetLineWidth(2);

        // Fit quietly (Q) and store results
        TFitResultPtr r = hSlice->Fit(fitBg, "RLSQ"); 
        //draw a darker frame box
        gPad->SetFrameLineWidth(3);
        hSlice->Draw("E1");

        // Draw a line to show the signal region cutoff
        TLine* lineSig = new TLine(sigMax, 0, sigMax, hSlice->GetMaximum());
        lineSig->SetLineColor(kBlue);
        lineSig->SetLineStyle(2);
        lineSig->Draw("SAME");

        // Calculate integrals for the legend
        double binWidthY = hSlice->GetBinWidth(1);
        double bgYield = fitBg->Integral(0.0, sigMax) / binWidthY;
        int sigBinMax = hSlice->FindBin(sigMax - 0.0001);
        double totalYield = hSlice->Integral(1, sigBinMax);
        double subYield = totalYield - bgYield;

        //latex box for slice info
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.04);
        latex.SetTextColor(kBlue);
        latex.DrawLatex(0.15, 0.65, Form("Mass Bin Center: %.2f GeV/c^{2}", massCenter));
        //instead of mass bin center, write the mass range for this bin
        double massBinLow = h2->GetXaxis()->GetBinLowEdge(i);
        double massBinHigh = h2->GetXaxis()->GetBinUpEdge(i);
        //latex.DrawLatex(0.15, 0.65, Form("Mass Bin: [%.2f, %.2f] GeV/c^{2}", massBinLow, massBinHigh));
       
        TLegend* leg = new TLegend(0.1, 0.75, 0.38, 0.88);
        leg->SetBorderSize(0);
        leg->SetTextSize(0.05);
        leg->AddEntry(hSlice, Form("Total in Signal Region: %.0f", totalYield), "lep");
        leg->AddEntry(fitBg, Form("Est. Background: %.0f", bgYield), "l");
        if (subYield < 0) subYield = 0; // Avoid negative signal yield in legend
        leg->AddEntry((TObject*)0, Form("Final Subtracted Signal: %.0f", subYield), "");
        leg->Draw();

        cFit->Update();
        cFit->Print(pdfName); // Add this canvas as a page in the PDF
        
        delete hSlice; // Clean up memory
    }

    // Close the PDF
    cFit->Print(pdfName + "]");
    std::cout << "Report saved to: " << pdfName << "\n";
    std::cout << "======================================================\n";
}
*/
void CheckBinByBinPtMissFits(TH2D* h2) {
    std::cout << "\n======================================================\n";
    std::cout << " GENERATING BIN-BY-BIN pTmiss FITS (PDF REPORT) \n";
    std::cout << "======================================================\n";

    // Limits for pTmiss
    double sigMax = 0.15;
    double tailMin = 0.20;
    double tailMax = 0.80;
    
    // --- STEP 1: GET THE GLOBAL SHAPE ---
    // Project the entire 2D plot to get high-statistics global pTmiss
    TH1D* hGlobal = h2->ProjectionY("hGlobal_for_shape");
    
    TF1* globalFit = new TF1("globalFit", "pol2", tailMin, tailMax);
    globalFit->FixParameter(0, 0.0);
    hGlobal->Fit(globalFit, "RQ0"); // Fit quietly
    
    double p1_glob = globalFit->GetParameter(1);
    double p2_glob = globalFit->GetParameter(2);
    
    // Calculate the 'b' parameter (shape) as requested by the senior
    double shape_b = p2_glob / p1_glob; 
    
    std::cout << "Global Shape Parameter 'b' fixed at: " << shape_b << "\n\n";

    // --- STEP 2: BIN-BY-BIN LOOP ---
    int binMin = h2->GetXaxis()->FindBin(0.9);
    int binMax = h2->GetXaxis()->FindBin(3.0);

    TCanvas* cFit = new TCanvas("cFit", "Bin-by-Bin Fits", 800, 600);
    TString pdfName = "plots/BinByBin_PtMiss_Fits_Report.pdf";
    cFit->Print(pdfName + "["); 

    gStyle->SetOptStat(0);

    for (int i = binMin; i <= binMax; i++) {
        double massCenter = h2->GetXaxis()->GetBinCenter(i);
        
        TH1D* hSlice = h2->ProjectionY(Form("hSlice_%d", i), i, i);
        
        // REBINNING: Merge 4 bins into 1 to eliminate empty bins in the tail.
        // change this to 2 or 5 
        hSlice->Rebin(10); 
        
        if (hSlice->GetEntries() < 3) continue;

         //pad margins for better fit display
        //gPad->SetLeftMargin(0.12);
        //gPad->SetBottomMargin(0.13); 
        hSlice->GetYaxis()->SetTitleOffset(0.8);
        hSlice->GetYaxis()->SetTitleSize(0.06);
        hSlice->GetYaxis()->SetLabelSize(0.06);
        /*
        hFinal->GetXaxis()->SetTitleOffset(0.9);
        hFinal->GetXaxis()->SetTitleSize(0.06);
        hFinal->GetXaxis()->SetLabelSize(0.05);*/

        hSlice->SetTitle(Form("p_{T}^{miss} for Mass = %.2f GeV/c^{2}; p_{T}^{miss} [GeV/c]; Events", massCenter));
        hSlice->SetMarkerStyle(20);
        hSlice->SetLineColor(kBlack);
        hSlice->GetYaxis()->SetRangeUser(0, hSlice->GetMaximum() * 1.99);
        hSlice->GetXaxis()->SetRangeUser(0.0, 1.0); 

        // NEW FITTING: a * (x + b*x^2)
        // [0] is 'a' (Normalization - floating)
        // [1] is 'b' (Shape - fixed)
        TF1* fitBg = new TF1(Form("fitBg_%d", i), "[0] * (x + [1]*x*x)", tailMin, tailMax);
        fitBg->FixParameter(1, shape_b); // Fix the shape from the global fit
        
        // Give 'a' a reasonable starting guess so the fit converges easily
        fitBg->SetParameter(0, hSlice->GetMaximum()); 
        
        fitBg->SetLineColor(kRed);
        fitBg->SetLineWidth(2);

        // Fit quietly (Q) and store results
        TFitResultPtr r = hSlice->Fit(fitBg, "RLSQ"); 
        
        gPad->SetFrameLineWidth(3);
        hSlice->Draw("E1");

        TLine* lineSig = new TLine(sigMax, 0, sigMax, hSlice->GetMaximum());
        lineSig->SetLineColor(kBlue);
        lineSig->SetLineStyle(2);
        lineSig->Draw("SAME");

        // Calculate integrals 
        double binWidthY = hSlice->GetBinWidth(1); // Get new rebinned width
        double bgYield = fitBg->Integral(0.0, sigMax) / binWidthY;
        int sigBinMax = hSlice->FindBin(sigMax - 0.0001);
        double totalYield = hSlice->Integral(1, sigBinMax);
        double subYield = totalYield - bgYield;

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.04);
        latex.SetTextColor(kBlue);
        latex.DrawLatex(0.15, 0.65, Form("Mass Bin Center: %.2f GeV/c^{2}", massCenter));
       
        TLegend* leg = new TLegend(0.1, 0.75, 0.38, 0.88);
        leg->SetBorderSize(0);
        leg->SetTextSize(0.05);
        leg->AddEntry(hSlice, Form("Total in Signal Region: %.0f", totalYield), "lep");
        leg->AddEntry(fitBg, Form("Est. Background: %.0f", bgYield), "l");
        if (subYield < 0) subYield = 0; 
        leg->AddEntry((TObject*)0, Form("Final Subtracted Signal: %.0f", subYield), "");
        leg->Draw();

        cFit->Update();
        cFit->Print(pdfName); 
        
        delete hSlice; 
    }

    cFit->Print(pdfName + "]");
    std::cout << "Report saved to: " << pdfName << "\n";
    std::cout << "======================================================\n";
}

void CheckTofBackgroundShapes(TH1D* h2TOF, TH1D* h3TOF, TH1D* h4TOF) {
    TCanvas* cShape = new TCanvas("cShape", "TOF Background Shapes", 800, 600);
    
    // Clone so we don't mess up your original histograms
    TH1D* h2 = (TH1D*)h2TOF->Clone("h2_shape");
    TH1D* h3 = (TH1D*)h3TOF->Clone("h3_shape");
    TH1D* h4 = (TH1D*)h4TOF->Clone("h4_shape");

    // Normalize all of them to an area of 1 so we can compare the pure SHAPE
    if (h2->Integral() > 0) h2->Scale(1.0 / h2->Integral());
    if (h3->Integral() > 0) h3->Scale(1.0 / h3->Integral());
    if (h4->Integral() > 0) h4->Scale(1.0 / h4->Integral());

    // Styling
    h2->SetLineColor(kBlue);
    h2->SetLineWidth(2);
    h3->SetLineColor(kGreen+2);
    h3->SetLineWidth(2);
    h4->SetLineColor(kRed);
    h4->SetLineWidth(2);

    // Zoom in on the relevant tail region (e.g., 0 to 1.0)
    h2->GetXaxis()->SetRangeUser(0.0, 1.0);
    // Find the max Y so they all fit on the plot
    double maxY = std::max({h2->GetMaximum(), h3->GetMaximum(), h4->GetMaximum()});
    h2->GetYaxis()->SetRangeUser(0, maxY * 1.2);

    h2->SetTitle("Normalized p_{T}^{miss} Shape Comparison; p_{T}^{miss} [GeV/c]; Normalized Yield");
    
    // Draw them overlaid
    h2->Draw("HIST");
    h3->Draw("HIST SAME");
    h4->Draw("HIST SAME");

    // Add a legend
    TLegend* leg = new TLegend(0.55, 0.70, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->AddEntry(h2, "2 TOF Hits", "l");
    leg->AddEntry(h3, "3 TOF Hits", "l");
    leg->AddEntry(h4, "4 TOF Hits", "l");
    leg->Draw();

    cShape->SaveAs("plots/pTmiss_Shape_Comparison.png");
}

void EstimateTofBackgrounds(TH1D* h2, TH1D* h3, TH1D* h4, const char* nameSuffix) {
    std::cout << "\n======================================================\n";
    std::cout << " TOF MULTIPLICITY BACKGROUND ESTIMATES \n";
    std::cout << "======================================================\n";

    TH1D* hists[3] = {h2, h3, h4};
    int tofCounts[3] = {2, 3, 4};
    
    // Create a wide canvas to show all 3 fits side-by-side
    TCanvas* cBg = new TCanvas("cBg", "TOF Backgrounds", 3000, 900);
    cBg->Divide(3, 1);

    double sigMax = 0.15;
    double tailMin = 0.20;
    double tailMax = 0.80;

    for (int i = 0; i < 3; i++) {
        cBg->cd(i + 1);
        TH1D* h = hists[i];

        //rebin the hist
        h->Rebin(10); // Adjust the rebinning factor as needed to ensure smooth fits

        // Format histogram
        //h->SetTitle(Form("%d TOF Hits p_{T}^{miss}; p_{T}^{miss} [GeV/c]; Events", tofCounts[i]));
        h->SetTitle("");
        h->SetMarkerStyle(20);
        h->SetMarkerSize(3.0);
        h->SetLineColor(kBlack);
        h->GetXaxis()->SetRangeUser(0.0, 1.0);
        //set y axis range
        h->GetYaxis()->SetRangeUser(0, h->GetMaximum() * 1.9);

        /*h->GetYaxis()->SetTitleOffset(1.95);
        h->GetYaxis()->SetTitleSize(0.08);
        h->GetYaxis()->SetLabelSize(0.06);
        
        h->GetXaxis()->SetTitleOffset(1.9);
        h->GetXaxis()->SetTitleSize(0.06);
        h->GetXaxis()->SetLabelSize(0.05);*/
        
        // Define standard pol2 background fit anchored at 0
        TF1* fit = new TF1(Form("fit_%dTOF", tofCounts[i]), "pol2", tailMin, tailMax);
        fit->FixParameter(0, 0.0);
        fit->SetLineColor(kRed);
        fit->SetLineWidth(2);
        
        h->Fit(fit, "RQ"); // R = range, Q = quiet
        h->Draw("E1");
        
        // Calculate yields
        int binMin = h->FindBin(0.0);
        int binMax = h->FindBin(sigMax - 0.0001);
        double totalYield = h->Integral(binMin, binMax);
        
        // Estimate background in signal region by integrating the fit
        double bgYield = fit->Integral(0.0, sigMax) / h->GetBinWidth(1);
        if (bgYield < 0) bgYield = 0; // Prevent negative background
        
        double sigYield = totalYield - bgYield;
        if (sigYield < 0) sigYield = 0;
        
        // Calculate what percentage of the signal region is just background
        double bgFraction = (totalYield > 0) ? (bgYield / totalYield) * 100.0 : 0.0;
       
        // --- Extract raw tail counts for systematic calculation ---
        int binTailMin = h->FindBin(tailMin);
        int binTailMax = h->FindBin(tailMax - 0.0001);
        double rawTailCounts = h->Integral(binTailMin, binTailMax);

        // Draw signal line cutoff
        TLine* line = new TLine(sigMax, 0, sigMax, h->GetMaximum());
        line->SetLineColor(kBlue);
        line->SetLineStyle(2);
        line->Draw("SAME");
        
        // Add info to plot
        TLegend* leg = new TLegend(0.35, 0.65, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->SetTextSize(0.045);
        leg->AddEntry(h, Form("Total: %.0f", totalYield), "lep");
        leg->AddEntry(fit, Form("Est. BG: %.0f (%.1f%%)", bgYield, bgFraction), "l");
        leg->AddEntry((TObject*)0, Form("Net Signal: %.0f", sigYield), "");
        leg->Draw();
        
        // Print clean summary to terminal for your senior
        std::cout << "--- " << tofCounts[i] << " TOF Hits ---\n";
        std::cout << "Total Events in Signal Region: " << totalYield << "\n";
        std::cout << "Estimated Background:          " << bgYield << " (" << bgFraction << "% background)\n";
        std::cout << "Extracted Net Signal:          " << sigYield << "\n\n";

        CalculateScaleFactorSystematic(rawTailCounts, bgYield, sigYield, Form("%d-TOF", tofCounts[i]));
    }
    
    //cBg->SaveAs("plots/TOF_Background_Estimates.png");
    cBg->SaveAs(Form("plots/TOF_Background_Estimates_%s.png", nameSuffix));
    std::cout << "======================================================\n";
}

TH1D* CombineTofHistograms(TH1D* h2, TH1D* h3, TH1D* h4) {
    // Clone one to get the same binning structure
    TH1D* hFinal = (TH1D*)h3->Clone("hFinal_CombinedMass");
    hFinal->Reset(); // Clear old contents
    hFinal->SetTitle("Final Combined K_{S}^{0}K_{S}^{0} Corrected Yield; m_{K_{S}^{0}K_{S}^{0}} [GeV/c^{2}]; Corrected Yield");
    
    int nBins = hFinal->GetNbinsX();
    
    for (int i = 1; i <= nBins; i++) {
        double y2 = h2->GetBinContent(i);
        double e2 = h2->GetBinError(i);
        double y3 = h3->GetBinContent(i);
        double e3 = h3->GetBinError(i);
        double y4 = h4->GetBinContent(i);
        double e4 = h4->GetBinError(i);
        
        double sumWeights = 0.0;
        double sumYields = 0.0;
        
        // Only include bins that actually have valid errors (prevent division by zero)
        if (e2 > 0) { double w2 = 1.0 / (e2 * e2); sumWeights += w2; sumYields += w2 * y2; }
        if (e3 > 0) { double w3 = 1.0 / (e3 * e3); sumWeights += w3; sumYields += w3 * y3; }
        if (e4 > 0) { double w4 = 1.0 / (e4 * e4); sumWeights += w4; sumYields += w4 * y4; }
        
        if (sumWeights > 0) {
            double finalYield = sumYields / sumWeights;
            double finalError = std::sqrt(1.0 / sumWeights);
            
            hFinal->SetBinContent(i, finalYield);
            hFinal->SetBinError(i, finalError);
        } else {
            hFinal->SetBinContent(i, 0);
            hFinal->SetBinError(i, 0);
        }
    }
    
    // Style the final plot
    hFinal->SetMarkerStyle(20);
    hFinal->SetMarkerColor(kBlack);
    hFinal->SetLineColor(kBlack);
    
    return hFinal;
}

TH1D* CalculateWeightedAverage(TH1D* h2, TH1D* h3, TH1D* h4) {
    std::cout << "\n======================================================\n";
    std::cout << " CALCULATING INVERSE-VARIANCE WEIGHTED AVERAGE \n";
    std::cout << "======================================================\n";

    // Clone one histogram to create the final framework
    TH1D* hFinal = (TH1D*)h4->Clone("hFinal_WeightedAverage");
    hFinal->Reset(); // Clear the contents, keep the binning
    hFinal->SetTitle("Final K_{S}^{0}K_{S}^{0} Corrected Yield (Weighted Avg); m_{K_{S}^{0}K_{S}^{0}} [GeV/c^{2}]; Corrected Yield");

    int nBins = hFinal->GetNbinsX();

    for (int i = 1; i <= nBins; i++) {
        double y2 = h2->GetBinContent(i); double e2 = h2->GetBinError(i);
        double y3 = h3->GetBinContent(i); double e3 = h3->GetBinError(i);
        double y4 = h4->GetBinContent(i); double e4 = h4->GetBinError(i);

        // If error is 0, the weight would be infinity (divide by zero). Skip it.
        double w2 = (e2 > 0) ? 1.0 / (e2 * e2) : 0;
        double w3 = (e3 > 0) ? 1.0 / (e3 * e3) : 0;
        double w4 = (e4 > 0) ? 1.0 / (e4 * e4) : 0;

        double sumWeights = w2 + w3 + w4;

        if (sumWeights > 0) {
            // The Weighted Average Formula
            double yAvg = (y2 * w2 + y3 * w3 + y4 * w4) / sumWeights;
            
            // The New Propagated Error Formula
            double eAvg = std::sqrt(1.0 / sumWeights);

            hFinal->SetBinContent(i, yAvg);
            hFinal->SetBinError(i, eAvg);
        } else {
            hFinal->SetBinContent(i, 0);
            hFinal->SetBinError(i, 0);
        }
    }

    // --- DRAW AND FIT THE FINAL RESULT ---
    TCanvas* cFinalAvg = new TCanvas("cFinalAvg", "Final Weighted Avg", 1100, 900);
    gPad->SetFrameLineWidth(3);
    //pad margins for better fit display
    gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.13); 
    hFinal->SetMarkerStyle(20);
    hFinal->SetMarkerSize(1.5);
    hFinal->SetMarkerColor(kBlack);
    hFinal->SetLineColor(kBlack);
    hFinal->GetXaxis()->SetRangeUser(1.0, 2.6); // Zoom in to your fit region
    hFinal->GetYaxis()->SetTitleOffset(0.8);
    hFinal->GetYaxis()->SetTitleSize(0.06);
    hFinal->GetYaxis()->SetLabelSize(0.06);
    
    hFinal->GetXaxis()->SetTitleOffset(0.9);
    hFinal->GetXaxis()->SetTitleSize(0.06);
    hFinal->GetXaxis()->SetLabelSize(0.05);
    // Fit the final, clean, background-subtracted data!
    TF1* finalFit = new TF1("finalFit", "gaus(0) + pol1(3)", 1.0, 2.6);
    // Give it reasonable starting guesses
    finalFit->SetParameters(hFinal->GetMaximum(), 1.68, 0.08, 0, 0); 
    finalFit->SetLineColor(kRed);
    finalFit->SetLineWidth(2);
    
    
    hFinal->Fit(finalFit, "R");
    hFinal->Draw("PE1");

    TLegend* leg = new TLegend(0.45, 0.70, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.04);
    leg->AddEntry(hFinal, "Weighted Avg (2,3,4 TOF)", "pe");
    leg->AddEntry(finalFit, "Signal Fit", "l");
    leg->Draw();

    // Print parameters on plot
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.60, 0.60, Form("Mass = %.3f GeV/c^{2}", finalFit->GetParameter(1)));
    latex.DrawLatex(0.60, 0.55, Form("Width = %.3f GeV/c^{2}", finalFit->GetParameter(2)));

    cFinalAvg->SaveAs("plots/Final_InvariantMass_WeightedAverage.png");
    
    return hFinal;
}

void CheckCorrectedYields(TH1D* h2, TH1D* h3, TH1D* h4) {
    std::cout << "\n======================================================\n";
    std::cout << " TOF MULTIPLICITY CORRECTED YIELDS & ERRORS \n";
    std::cout << "======================================================\n";

    TH1D* hists[3] = {h2, h3, h4};
    int tof[3] = {2, 3, 4};
    double sigMax = 0.15;
    double tailMin = 0.20;
    double tailMax = 0.80;

    for (int i = 0; i < 3; i++) {
        TH1D* h = hists[i];
        //h->Rebin(10); // Ensure same rebinning as fits for consistency
        // 1. Get the Statistical Error of the Total Region correctly using SumW2
        int binMin = h->FindBin(0.0);
        int binMax = h->FindBin(sigMax - 0.0001);
        double errTot = 0;
        
        // IntegralAndError automatically uses the weighted SumW2 errors!
        double yTot = h->IntegralAndError(binMin, binMax, errTot);

        // 2. Fit the background (using 'S' to save the fit result matrix)
        TF1* fit = new TF1(Form("fit_%d", tof[i]), "pol2", tailMin, tailMax);
        fit->FixParameter(0, 0.0);
        TFitResultPtr r = h->Fit(fit, "RQS"); 

        // 3. Background Yield
        double binWidth = h->GetBinWidth(1);
        double yBg = fit->Integral(0.0, sigMax) / binWidth;

        // 4. Background Error (Propagated from the fit's covariance matrix)
        double errBg = 0;
        if(r->IsValid()) {
            double params[3] = {fit->GetParameter(0), fit->GetParameter(1), fit->GetParameter(2)};
            TMatrixDSym covMatrix = r->GetCovarianceMatrix();
            errBg = fit->IntegralError(0.0, sigMax, params, covMatrix.GetMatrixArray()) / binWidth;
        }

        // 5. Final Corrected Signal
        double sigYield = yTot - yBg;
        
        // Errors add in quadrature
        double sigErr = std::sqrt(errTot * errTot + errBg * errBg);

        std::cout << "--- " << tof[i] << " TOF Hits ---\n";
        std::cout << "Total Corrected in Signal: " << yTot << " +/- " << errTot << "\n";
        std::cout << "Est. Corrected Background: " << yBg << " +/- " << errBg << "\n";
        std::cout << "NET CORRECTED SIGNAL:      " << sigYield << " +/- " << sigErr << "\n\n";
    }
    std::cout << "======================================================\n";
}

TH1D* ProcessFinalYield(TH2D* h2_2TOF, TH2D* h2_3TOF, TH2D* h2_4TOF, TH1D* hRaw1D, TH1D* hCorr1D, const char* suffix) {
    std::cout << "\n======================================================\n";
    std::cout << " PROCESSING DATASET: " << suffix << "\n";
    std::cout << "======================================================\n";

    // 1. Extract the Raw Signal using your shape method
    // We add the suffix to the names so ROOT doesn't overwrite them in memory
    TH1D* hCleanRaw_2TOF = ExtractDifferentialCrossSection2(h2_2TOF, Form("2TOF_Raw_%s", suffix));
    TH1D* hCleanRaw_3TOF = ExtractDifferentialCrossSection2(h2_3TOF, Form("3TOF_Raw_%s", suffix));
    TH1D* hCleanRaw_4TOF = ExtractDifferentialCrossSection2(h2_4TOF, Form("4TOF_Raw_%s", suffix));

    // 2. SUM the raw slices together
    TH1D* hFinalCleanRaw = (TH1D*)hCleanRaw_4TOF->Clone(Form("hFinalCleanRaw_%s", suffix));
    hFinalCleanRaw->Add(hCleanRaw_3TOF);
    hFinalCleanRaw->Add(hCleanRaw_2TOF);

    // 3. Derive the specific 1D Efficiency Curve for this dataset
    TH1D* hEff1D = (TH1D*)hRaw1D->Clone(Form("hEff1D_%s", suffix));
    hEff1D->Divide(hCorr1D); 

    // 4. Clone the clean raw signal and apply the efficiency correction
    TH1D* hFinalCorrected = (TH1D*)hFinalCleanRaw->Clone(Form("hFinalCorrected_%s", suffix));
    hFinalCorrected->SetTitle(Form("Final Corrected Signal (%s); m_{K_{S}^{0}K_{S}^{0}} [GeV/c^{2}]; Corrected Yield", suffix));
    hFinalCorrected->Divide(hEff1D); 

    // 5. Draw, Fit, and Save
    TCanvas* cFinal = new TCanvas(Form("cFinal_%s", suffix), Form("Final Yield %s", suffix), 1100, 900);
    gPad->SetFrameLineWidth(3);
    gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.13); 

    hFinalCorrected->SetMarkerStyle(20);
    hFinalCorrected->SetMarkerSize(1.5);
    hFinalCorrected->SetMarkerColor(kBlack);
    hFinalCorrected->SetLineColor(kBlack);
    hFinalCorrected->GetXaxis()->SetRangeUser(1.0, 2.6);
    hFinalCorrected->GetYaxis()->SetTitleOffset(0.97);
    hFinalCorrected->GetYaxis()->SetTitleSize(0.06);
    hFinalCorrected->GetYaxis()->SetLabelSize(0.06);
    hFinalCorrected->GetXaxis()->SetTitleOffset(0.9);
    hFinalCorrected->GetXaxis()->SetTitleSize(0.06);
    hFinalCorrected->GetXaxis()->SetLabelSize(0.05);

    TF1* finalFit = new TF1(Form("fit_%s", suffix), "gaus(0) + pol1(3)", 1.0, 2.6);
    finalFit->SetParameters(hFinalCorrected->GetMaximum(), 1.68, 0.08, 0, 0); 
    finalFit->SetLineColor(kRed);
    finalFit->SetLineWidth(2);
    finalFit->SetNpx(1000); 

    hFinalCorrected->Fit(finalFit, "R");
    hFinalCorrected->Draw("PE1");

    TLegend* legSum = new TLegend(0.55, 0.70, 0.88, 0.88);
    legSum->SetBorderSize(0);
    legSum->SetTextSize(0.04);
    legSum->AddEntry(hFinalCorrected, Form("Corrected Yield"), "pe");
    legSum->AddEntry(finalFit, "Signal Fit", "l");
    legSum->Draw();

    TLatex latexSum;
    latexSum.SetNDC();
    latexSum.SetTextSize(0.04);
    latexSum.DrawLatex(0.60, 0.60, Form("Mass = %.3f GeV/c^{2}", finalFit->GetParameter(1)));
    latexSum.DrawLatex(0.60, 0.55, Form("Width = %.3f GeV/c^{2}", finalFit->GetParameter(2)));

    cFinal->Update();
    cFinal->SaveAs(Form("plots/Final_InvariantMass_Corrected_%s.png", suffix));

    // Return the histogram so it can be used for the Systematics Plot!
    return hFinalCorrected;
}

TH1D* ExtractHybridCorrectedYield(TH2D* h2Raw, TH2D* h2Corr, const char* nameSuffix) {
    double sigMax = 0.15;
    double tailMin = 0.20;
    double tailMax = 0.80;

    // -------------------------------------------------------------
    // STEP 1: Use RAW data to find the Scale Factor (S)
    // -------------------------------------------------------------
    TH1D* hPtMissRaw = h2Raw->ProjectionY(Form("tempRawY_%s", nameSuffix));
    TF1* fitRaw = new TF1(Form("fitRaw_%s", nameSuffix), "pol2", tailMin, tailMax);
    fitRaw->FixParameter(0, 0.0);
    hPtMissRaw->Fit(fitRaw, "RQ0");
    
    double binWidth = hPtMissRaw->GetBinWidth(1);
    double bgInSigRegionRaw = fitRaw->Integral(0.0, sigMax) / binWidth;
    
    int binTailMinRaw = hPtMissRaw->FindBin(tailMin);
    int binTailMaxRaw = hPtMissRaw->FindBin(tailMax - 0.0001);
    double bgInTailRegionRaw = hPtMissRaw->Integral(binTailMinRaw, binTailMaxRaw);
    
    // The ratio of background under the peak to background in the tail
    double scaleFactor = (bgInTailRegionRaw > 0) ? (bgInSigRegionRaw / bgInTailRegionRaw) : 0;

    // -------------------------------------------------------------
    // STEP 2: Extract Mass Shapes from CORRECTED (Weighted) data
    // -------------------------------------------------------------
    int binSigMin = h2Corr->GetYaxis()->FindBin(0.0);
    int binSigMax = h2Corr->GetYaxis()->FindBin(sigMax - 0.0001);
    TH1D* hSigRegionCorr = h2Corr->ProjectionX(Form("hSigRegionCorr_%s", nameSuffix), binSigMin, binSigMax);

    int binTailMin = h2Corr->GetYaxis()->FindBin(tailMin);
    int binTailMax = h2Corr->GetYaxis()->FindBin(tailMax - 0.0001);
    TH1D* hTailRegionCorr = h2Corr->ProjectionX(Form("hTailRegionCorr_%s", nameSuffix), binTailMin, binTailMax);

    // -------------------------------------------------------------
    // STEP 3: Scale the Corrected Tail & Subtract
    // -------------------------------------------------------------
    TH1D* hBgEstimateCorr = (TH1D*)hTailRegionCorr->Clone(Form("hBgEstimateCorr_%s", nameSuffix));
    hBgEstimateCorr->Scale(scaleFactor);

    TH1D* hFinalCorrectedYield = (TH1D*)hSigRegionCorr->Clone(Form("hFinalCorrectedYield_%s", nameSuffix));
    hFinalCorrectedYield->Add(hBgEstimateCorr, -1.0);
    
    return hFinalCorrectedYield;
}
/*

void ClosureTest_TailStability(TH2D* h2Corr, const char* suffix) {
    std::cout << "\n======================================================\n";
    std::cout << " CLOSURE TEST 3: pTmiss Tail Mass Shape Stability (" << suffix << ")\n";
    std::cout << "======================================================\n";

    // Define two sub-regions in the pTmiss tail
    // Region A: 0.20 to 0.40 GeV/c
    // Region B: 0.40 to 0.80 GeV/c
    int binA_min = h2Corr->GetYaxis()->FindBin(0.20);
    int binA_max = h2Corr->GetYaxis()->FindBin(0.40 - 0.0001);
    int binB_min = h2Corr->GetYaxis()->FindBin(0.40);
    int binB_max = h2Corr->GetYaxis()->FindBin(0.80 - 0.0001);

    TH1D* hMassA = h2Corr->ProjectionX(Form("hMassA_%s", suffix), binA_min, binA_max);
    TH1D* hMassB = h2Corr->ProjectionX(Form("hMassB_%s", suffix), binB_min, binB_max);

    // Rebin if necessary to ensure enough stats for a valid KS test
    hMassA->Rebin(2);
    hMassB->Rebin(2);

    // Normalize both to an area of 1 so we are strictly comparing SHAPE
    if (hMassA->Integral() > 0) hMassA->Scale(1.0 / hMassA->Integral());
    if (hMassB->Integral() > 0) hMassB->Scale(1.0 / hMassB->Integral());

    // Perform the Kolmogorov-Smirnov Test
    double ksProb = hMassA->KolmogorovTest(hMassB);
    std::cout << "KS Test Probability (Region A vs Region B): " << ksProb << "\n";
    if (ksProb > 0.05) {
        std::cout << "-> SUCCESS: Shapes are statistically consistent.\n";
    } else {
        std::cout << "-> WARNING: Shapes show statistical divergence.\n";
    }

    // --- Plotting ---
    TCanvas* cClosure3 = new TCanvas(Form("cClosure3_%s", suffix), "Tail Stability", 800, 600);
    gPad->SetFrameLineWidth(3);
    gPad->SetLeftMargin(0.12);

    hMassA->SetLineColor(kBlue);
    hMassA->SetMarkerColor(kBlue);
    hMassA->SetMarkerStyle(20);
    hMassA->SetTitle(Form("Background Mass Shape Stability (%s); m_{K_{S}^{0}K_{S}^{0}} [GeV/c^{2}]; Normalized Yield", suffix));
    hMassA->GetXaxis()->SetRangeUser(0.9, 3.0);
    
    hMassB->SetLineColor(kRed);
    hMassB->SetMarkerColor(kRed);
    hMassB->SetMarkerStyle(21);

    double maxY = std::max(hMassA->GetMaximum(), hMassB->GetMaximum());
    hMassA->GetYaxis()->SetRangeUser(0, maxY * 1.4);

    hMassA->Draw("E1");
    hMassB->Draw("E1 SAME");

    TLegend* leg = new TLegend(0.35, 0.75, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->AddEntry(hMassA, "Region A: 0.20 < p_{T}^{miss} < 0.40", "lep");
    leg->AddEntry(hMassB, "Region B: 0.40 < p_{T}^{miss} < 0.80", "lep");
    leg->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    if (ksProb > 0.05) latex.SetTextColor(kGreen+2);
    else latex.SetTextColor(kRed);
    latex.DrawLatex(0.40, 0.68, Form("KS Test Prob: %.3f", ksProb));

    cClosure3->SaveAs(Form("plots/ClosureTest_TailStability_%s.png", suffix));
}
*/
void ClosureTest_TailStability(TH2D* h2Corr, const char* suffix) {
    std::cout << "\n======================================================\n";
    std::cout << " CLOSURE TEST 3: pTmiss Tail Mass Shape (Chi-Square) (" << suffix << ")\n";
    std::cout << "======================================================\n";

    // Define two sub-regions in the pTmiss tail
    // Region A: 0.20 to 0.40 GeV/c
    // Region B: 0.40 to 0.80 GeV/c
    int binA_min = h2Corr->GetYaxis()->FindBin(0.20);
    int binA_max = h2Corr->GetYaxis()->FindBin(0.40 - 0.0001);
    int binB_min = h2Corr->GetYaxis()->FindBin(0.40);
    int binB_max = h2Corr->GetYaxis()->FindBin(0.80 - 0.0001);

    TH1D* hMassA = h2Corr->ProjectionX(Form("hMassA_%s", suffix), binA_min, binA_max);
    TH1D* hMassB = h2Corr->ProjectionX(Form("hMassB_%s", suffix), binB_min, binB_max);

    // Rebin if necessary to ensure enough stats per bin for a valid Chi2 test
    // You may need to increase this rebin factor if ROOT warns about empty bins!
    hMassA->Rebin(2);
    hMassB->Rebin(2);

    // ------------------------------------------------------------------
    // PERFORM CHI-SQUARE TEST *BEFORE* SCALING
    // Option "WW N": 
    // "WW" = Both histograms are weighted (efficiency corrected)
    // "N"  = Normalize them to the same area before comparing shapes
    // "P"  = Print the detailed Chi2 result to the terminal
    // ------------------------------------------------------------------
    double chi2Prob = hMassA->Chi2Test(hMassB, "WW N P");
    
    std::cout << "Chi2 Test Probability (Region A vs Region B): " << chi2Prob << "\n";
    if (chi2Prob > 0.05) {
        std::cout << "-> SUCCESS: Shapes are statistically consistent.\n";
    } else {
        std::cout << "-> WARNING: Shapes show statistical divergence.\n";
    }

    // Now normalize both to an area of 1 so we can overlay them visually on the plot
    if (hMassA->Integral() > 0) hMassA->Scale(1.0 / hMassA->Integral());
    if (hMassB->Integral() > 0) hMassB->Scale(1.0 / hMassB->Integral());

    // --- Plotting ---
    TCanvas* cClosure3 = new TCanvas(Form("cClosure3_%s", suffix), "Tail Stability", 800, 600);
    gPad->SetFrameLineWidth(3);
    gPad->SetLeftMargin(0.12);

    hMassA->SetLineColor(kBlue);
    hMassA->SetMarkerColor(kBlue);
    hMassA->SetMarkerStyle(20);
    hMassA->SetTitle(Form("Background Mass Shape Stability (%s); m_{K_{S}^{0}K_{S}^{0}} [GeV/c^{2}]; Normalized Yield", suffix));
    hMassA->GetXaxis()->SetRangeUser(0.9, 3.0);
    
    hMassB->SetLineColor(kRed);
    hMassB->SetMarkerColor(kRed);
    hMassB->SetMarkerStyle(21);

    double maxY = std::max(hMassA->GetMaximum(), hMassB->GetMaximum());
    hMassA->GetYaxis()->SetRangeUser(0, maxY * 1.4);

    hMassA->Draw("E1");
    hMassB->Draw("E1 SAME");

    TLegend* leg = new TLegend(0.35, 0.75, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->AddEntry(hMassA, "Region A: 0.20 < p_{T}^{miss} < 0.40", "lep");
    leg->AddEntry(hMassB, "Region B: 0.40 < p_{T}^{miss} < 0.80", "lep");
    leg->Draw();

    // Print the Chi2 Probability on the plot
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    if (chi2Prob > 0.05) latex.SetTextColor(kGreen+2);
    else latex.SetTextColor(kRed);
    latex.DrawLatex(0.40, 0.68, Form("#chi^{2} Test Prob: %.3f", chi2Prob));

    cClosure3->SaveAs(Form("plots/ClosureTest_TailStability_%s_Chi2.png", suffix));
}

void ClosureTest_EfficiencyBias(TH2D* h2Raw, TH2D* h2Corr, const char* suffix) {
    std::cout << "\n======================================================\n";
    std::cout << " CLOSURE TEST 4: Efficiency Bias on pTmiss Shape (" << suffix << ")\n";
    std::cout << "======================================================\n";

    double sigMax = 0.15;
    double tailMin = 0.20;
    double tailMax = 0.80;

    TH1D* hPtRaw = h2Raw->ProjectionY(Form("hPtRawClosure_%s", suffix));
    TH1D* hPtCorr = h2Corr->ProjectionY(Form("hPtCorrClosure_%s", suffix));

    // Rebin exactly as you did in EstimateTofBackgrounds for stability
    hPtRaw->Rebin(10);
    hPtCorr->Rebin(10);

    // --- 1. FIT RAW DATA ---
    TF1* fitRaw = new TF1(Form("fitRawClosure_%s", suffix), "pol2", tailMin, tailMax);
    fitRaw->FixParameter(0, 0.0);
    hPtRaw->Fit(fitRaw, "RQ0"); // RQ0 = Range, Quiet, Do not draw
    
    double bgSigRaw = fitRaw->Integral(0.0, sigMax) / hPtRaw->GetBinWidth(1);
    int binTMinRaw = hPtRaw->FindBin(tailMin);
    int binTMaxRaw = hPtRaw->FindBin(tailMax - 0.0001);
    double bgTailRaw = hPtRaw->Integral(binTMinRaw, binTMaxRaw);
    double scaleRaw = (bgTailRaw > 0) ? (bgSigRaw / bgTailRaw) : 0;

    // --- 2. FIT CORRECTED DATA ---
    TF1* fitCorr = new TF1(Form("fitCorrClosure_%s", suffix), "pol2", tailMin, tailMax);
    fitCorr->FixParameter(0, 0.0);
    hPtCorr->Fit(fitCorr, "RQ0");
    
    double bgSigCorr = fitCorr->Integral(0.0, sigMax) / hPtCorr->GetBinWidth(1);
    int binTMinCorr = hPtCorr->FindBin(tailMin);
    int binTMaxCorr = hPtCorr->FindBin(tailMax - 0.0001);
    double bgTailCorr = hPtCorr->Integral(binTMinCorr, binTMaxCorr);
    double scaleCorr = (bgTailCorr > 0) ? (bgSigCorr / bgTailCorr) : 0;

    // --- 3. COMPARE ---
    double diffPercent = 0.0;
    if (scaleRaw > 0) diffPercent = std::abs(scaleRaw - scaleCorr) / scaleRaw * 100.0;

    std::cout << "Calculated Scale Factor (Raw Data):       " << scaleRaw << "\n";
    std::cout << "Calculated Scale Factor (Corrected Data): " << scaleCorr << "\n";
    std::cout << "Percentage Difference:                    " << diffPercent << "%\n";
    if (diffPercent < 5.0) {
        std::cout << "-> SUCCESS: Efficiency weights do not bias the scaling.\n";
    } else {
        std::cout << "-> WARNING: Efficiency weights significantly alter the background shape.\n";
    }

    // --- Plotting Normalized Shapes to visually prove it ---
    TCanvas* cClosure4 = new TCanvas(Form("cClosure4_%s", suffix), "Efficiency Bias", 800, 600);
    gPad->SetFrameLineWidth(3);
    
    TH1D* hPtRawNorm = (TH1D*)hPtRaw->Clone(Form("hPtRawNorm_%s", suffix));
    TH1D* hPtCorrNorm = (TH1D*)hPtCorr->Clone(Form("hPtCorrNorm_%s", suffix));
    
    // Normalize to 1
    if (hPtRawNorm->Integral() > 0) hPtRawNorm->Scale(1.0 / hPtRawNorm->Integral());
    if (hPtCorrNorm->Integral() > 0) hPtCorrNorm->Scale(1.0 / hPtCorrNorm->Integral());

    hPtRawNorm->SetLineColor(kBlack);
    hPtRawNorm->SetMarkerColor(kBlack);
    hPtRawNorm->SetMarkerStyle(20);
    hPtRawNorm->SetTitle(Form("p_{T}^{miss} Shape Comparison (%s); p_{T}^{miss} [GeV/c]; Normalized Yield", suffix));
    hPtRawNorm->GetXaxis()->SetRangeUser(0.0, 1.0);
    
    hPtCorrNorm->SetLineColor(kBlue);
    hPtCorrNorm->SetMarkerColor(kBlue);
    hPtCorrNorm->SetMarkerStyle(21);

    double maxY = std::max(hPtRawNorm->GetMaximum(), hPtCorrNorm->GetMaximum());
    hPtRawNorm->GetYaxis()->SetRangeUser(0, maxY * 1.3);

    hPtRawNorm->Draw("E1");
    hPtCorrNorm->Draw("E1 SAME");

    TLegend* leg = new TLegend(0.45, 0.75, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->AddEntry(hPtRawNorm, "Raw Shape", "lep");
    leg->AddEntry(hPtCorrNorm, "Corrected Shape", "lep");
    leg->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    latex.DrawLatex(0.50, 0.68, Form("Scale Factor Raw: %.3f", scaleRaw));
    latex.DrawLatex(0.50, 0.63, Form("Scale Factor Corr: %.3f", scaleCorr));
    latex.DrawLatex(0.50, 0.58, Form("Difference: %.1f%%", diffPercent));

    cClosure4->SaveAs(Form("plots/ClosureTest_EfficiencyBias_%s.png", suffix));
}

void Plot_EffCorrections()
{
    // === LOAD FILES ===
    string DataFileNom = "Sep17Data_MW4753.root";//EffCorr_June22//~/Downloads/EffCorr_nHits20_AorB_pTmiss_May21.root";  //EffCorr_nHits20_AorB_2  EffCorr_nHits20_AorB_pTmiss_May21
    string DataFile17  = "~/Downloads/EffCorr_nHits17_AorB_pTmiss_May25.root";
    string DataFile22  = "~/Downloads/EffCorr_nHits22_AorB_pTmiss_May25.root";

    TFile* fcorr    = TFile::Open(DataFileNom.c_str(), "READ");
    TFile* fcorr_17 = TFile::Open(DataFile17.c_str(), "READ");
    TFile* fcorr_22 = TFile::Open(DataFile22.c_str(), "READ");

    // === NOMINAL HISTOGRAMS ===
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
    hInvMassCorr->Scale(1.0/3.0);

    TH1D *hInvMassRaw = (TH1D*)h1D_Reco_InvMass_Raw_4pi_4TOF->Clone("hInvMassRaw");
    hInvMassRaw->Add(h1D_Reco_InvMass_Raw_4pi_3TOF);
    hInvMassRaw->Add(h1D_Reco_InvMass_Raw_4pi_2TOF);
    hInvMassRaw->Scale(1.0/3.0);

    TH1D *h1D_PtMiss_Raw = (TH1D*)fcorr->Get("hPtMiss_Raw");
    TH1D *h1D_PtMiss_Corrected = (TH1D*)fcorr->Get("hPtMiss_Corrected");

    TH2D *h2_PtMiss_Vs_Mass = (TH2D*)fcorr->Get("h2_PtMiss_Vs_Mass");
    TH2D *h2_PtMiss_Vs_Mass_Corrected = (TH2D*)fcorr->Get("h2_PtMiss_Vs_Mass_Corrected");
    
    TH2D *h2_PtMiss_Vs_Mass_4TOF = (TH2D*)fcorr->Get("h2_PtMiss_Vs_Mass_4TOF");
    TH2D *h2_PtMiss_Vs_Mass_4TOF_Corrected = (TH2D*)fcorr->Get("h2_PtMiss_Vs_Mass_4TOF_Corrected");

    TH2D *h2_PtMiss_Vs_Mass_3TOF = (TH2D*)fcorr->Get("h2_PtMiss_Vs_Mass_3TOF");
    TH2D *h2_PtMiss_Vs_Mass_3TOF_Corrected = (TH2D*)fcorr->Get("h2_PtMiss_Vs_Mass_3TOF_Corrected");
    
    TH2D *h2_PtMiss_Vs_Mass_2TOF = (TH2D*)fcorr->Get("h2_PtMiss_Vs_Mass_2TOF");
    TH2D *h2_PtMiss_Vs_Mass_2TOF_Corrected = (TH2D*)fcorr->Get("h2_PtMiss_Vs_Mass_2TOF_Corrected");

    TH1D* h1D_Reco_PtMiss_Raw_4pi_4TOF = (TH1D*)fcorr->Get("h1D_Reco_PtMiss_Raw_4pi_4TOF");
    TH1D* h1D_Reco_PtMiss_Raw_4pi_3TOF = (TH1D*)fcorr->Get("h1D_Reco_PtMiss_Raw_4pi_3TOF");
    TH1D* h1D_Reco_PtMiss_Raw_4pi_2TOF = (TH1D*)fcorr->Get("h1D_Reco_PtMiss_Raw_4pi_2TOF");
    TH1D* h1D_Reco_PtMiss_Corrected_4pi_4TOF = (TH1D*)fcorr->Get("h1D_Reco_PtMiss_Corrected_4pi_4TOF");
    TH1D* h1D_Reco_PtMiss_Corrected_4pi_3TOF = (TH1D*)fcorr->Get("h1D_Reco_PtMiss_Corrected_4pi_3TOF");
    TH1D* h1D_Reco_PtMiss_Corrected_4pi_2TOF = (TH1D*)fcorr->Get("h1D_Reco_PtMiss_Corrected_4pi_2TOF");


    // === SYSTEMATICS (nHits = 17) HISTOGRAMS ===
    TH1D* hPt_P_C_17   = (TH1D*)fcorr_17->Get("h1D_CorrectedPt_P");
    TH1D* hPt_N_C_17   = (TH1D*)fcorr_17->Get("h1D_CorrectedPt_N");
    TH1D* hEta_P_C_17  = (TH1D*)fcorr_17->Get("h1D_CorrectedEta_P");
    TH1D* hEta_N_C_17  = (TH1D*)fcorr_17->Get("h1D_CorrectedEta_N");
    TH1D* hVerZ_P_C_17 = (TH1D*)fcorr_17->Get("h1D_CorrectedVz_P");
    TH1D* hVerZ_N_C_17 = (TH1D*)fcorr_17->Get("h1D_CorrectedVz_N");

   
    TH1D* h1D_Reco_InvMass_Raw_4pi_4TOF_17 = (TH1D*)fcorr_17->Get("h1D_Reco_InvMass_Raw_4pi_4TOF");
    TH1D* h1D_Reco_InvMass_Raw_4pi_3TOF_17 = (TH1D*)fcorr_17->Get("h1D_Reco_InvMass_Raw_4pi_3TOF");
    TH1D* h1D_Reco_InvMass_Raw_4pi_2TOF_17 = (TH1D*)fcorr_17->Get("h1D_Reco_InvMass_Raw_4pi_2TOF");

    TH1D *hInvMassRaw_17 = (TH1D*)h1D_Reco_InvMass_Raw_4pi_4TOF_17->Clone("hInvMassRaw");
    hInvMassRaw_17->Add(h1D_Reco_InvMass_Raw_4pi_3TOF_17);
    hInvMassRaw_17->Add(h1D_Reco_InvMass_Raw_4pi_2TOF_17);
    hInvMassRaw_17->Scale(1.0/3.0);

    TH1D* h1D_Mass_17_4TOF = (TH1D*)fcorr_17->Get("h1D_Reco_InvMass_Corrected_4pi_4TOF");
    TH1D* h1D_Mass_17_3TOF = (TH1D*)fcorr_17->Get("h1D_Reco_InvMass_Corrected_4pi_3TOF");
    TH1D* h1D_Mass_17_2TOF = (TH1D*)fcorr_17->Get("h1D_Reco_InvMass_Corrected_4pi_2TOF");
    
    TH1D* hInvMassCorr_17 = (TH1D*)h1D_Mass_17_4TOF->Clone("hInvMassCorr_17");
    hInvMassCorr_17->Add(h1D_Mass_17_3TOF);
    hInvMassCorr_17->Add(h1D_Mass_17_2TOF);
    hInvMassCorr_17->Scale(1.0/3.0);

    TH2D *h2_PtMiss_Vs_Mass_17 = (TH2D*)fcorr_17->Get("h2_PtMiss_Vs_Mass");
    TH2D *h2_PtMiss_Vs_Mass_Corrected_17 = (TH2D*)fcorr_17->Get("h2_PtMiss_Vs_Mass_Corrected");

    TH2D *h2_PtMiss_Vs_Mass_4TOF_17 = (TH2D*)fcorr_17->Get("h2_PtMiss_Vs_Mass_4TOF");
    TH2D *h2_PtMiss_Vs_Mass_4TOF_Corrected_17 = (TH2D*)fcorr_17->Get("h2_PtMiss_Vs_Mass_4TOF_Corrected");

    TH2D *h2_PtMiss_Vs_Mass_3TOF_17 = (TH2D*)fcorr_17->Get("h2_PtMiss_Vs_Mass_3TOF");
    TH2D *h2_PtMiss_Vs_Mass_3TOF_Corrected_17 = (TH2D*)fcorr_17->Get("h2_PtMiss_Vs_Mass_3TOF_Corrected");

    TH2D *h2_PtMiss_Vs_Mass_2TOF_17 = (TH2D*)fcorr_17->Get("h2_PtMiss_Vs_Mass_2TOF");
    TH2D *h2_PtMiss_Vs_Mass_2TOF_Corrected_17 = (TH2D*)fcorr_17->Get("h2_PtMiss_Vs_Mass_2TOF_Corrected");
    

    // === SYSTEMATICS (nHits = 22) HISTOGRAMS ===
    TH1D* hPt_P_C_22   = (TH1D*)fcorr_22->Get("h1D_CorrectedPt_P");
    TH1D* hPt_N_C_22   = (TH1D*)fcorr_22->Get("h1D_CorrectedPt_N");
    TH1D* hEta_P_C_22  = (TH1D*)fcorr_22->Get("h1D_CorrectedEta_P");
    TH1D* hEta_N_C_22  = (TH1D*)fcorr_22->Get("h1D_CorrectedEta_N");
    TH1D* hVerZ_P_C_22 = (TH1D*)fcorr_22->Get("h1D_CorrectedVz_P");
    TH1D* hVerZ_N_C_22 = (TH1D*)fcorr_22->Get("h1D_CorrectedVz_N");
    
    TH1D* h1D_Reco_InvMass_Raw_4pi_4TOF_22 = (TH1D*)fcorr_22->Get("h1D_Reco_InvMass_Raw_4pi_4TOF");
    TH1D* h1D_Reco_InvMass_Raw_4pi_3TOF_22 = (TH1D*)fcorr_22->Get("h1D_Reco_InvMass_Raw_4pi_3TOF");
    TH1D* h1D_Reco_InvMass_Raw_4pi_2TOF_22 = (TH1D*)fcorr_22->Get("h1D_Reco_InvMass_Raw_4pi_2TOF");

    TH1D *hInvMassRaw_22 = (TH1D*)h1D_Reco_InvMass_Raw_4pi_4TOF_22->Clone("hInvMassRaw");
    hInvMassRaw_22->Add(h1D_Reco_InvMass_Raw_4pi_3TOF_22);
    hInvMassRaw_22->Add(h1D_Reco_InvMass_Raw_4pi_2TOF_22);
    hInvMassRaw_22->Scale(1.0/3.0);

    TH1D* h1D_Mass_22_4TOF = (TH1D*)fcorr_22->Get("h1D_Reco_InvMass_Corrected_4pi_4TOF");
    TH1D* h1D_Mass_22_3TOF = (TH1D*)fcorr_22->Get("h1D_Reco_InvMass_Corrected_4pi_3TOF");
    TH1D* h1D_Mass_22_2TOF = (TH1D*)fcorr_22->Get("h1D_Reco_InvMass_Corrected_4pi_2TOF");

    TH1D *hInvMassCorr_22 = (TH1D*)h1D_Mass_22_4TOF->Clone("hInvMassCorr_22");
    hInvMassCorr_22->Add(h1D_Mass_22_3TOF);
    hInvMassCorr_22->Add(h1D_Mass_22_2TOF);
    hInvMassCorr_22->Scale(1.0/3.0);

    TH2D *h2_PtMiss_Vs_Mass_22 = (TH2D*)fcorr_22->Get("h2_PtMiss_Vs_Mass");
    TH2D *h2_PtMiss_Vs_Mass_Corrected_22 = (TH2D*)fcorr_22->Get("h2_PtMiss_Vs_Mass_Corrected");

    TH2D *h2_PtMiss_Vs_Mass_4TOF_22 = (TH2D*)fcorr_22->Get("h2_PtMiss_Vs_Mass_4TOF");
    TH2D *h2_PtMiss_Vs_Mass_4TOF_Corrected_22 = (TH2D*)fcorr_22->Get("h2_PtMiss_Vs_Mass_4TOF_Corrected");

    TH2D *h2_PtMiss_Vs_Mass_3TOF_22 = (TH2D*)fcorr_22->Get("h2_PtMiss_Vs_Mass_3TOF");
    TH2D *h2_PtMiss_Vs_Mass_3TOF_Corrected_22 = (TH2D*)fcorr_22->Get("h2_PtMiss_Vs_Mass_3TOF_Corrected");

    TH2D *h2_PtMiss_Vs_Mass_2TOF_22 = (TH2D*)fcorr_22->Get("h2_PtMiss_Vs_Mass_2TOF");
    TH2D *h2_PtMiss_Vs_Mass_2TOF_Corrected_22 = (TH2D*)fcorr_22->Get("h2_PtMiss_Vs_Mass_2TOF_Corrected");

 
    sethist(hPt_P_C, hPt_P_R);
    sethist(hPt_N_C, hPt_N_R);
    sethist(hEta_P_C, hEta_P_R);
    sethist(hEta_N_C, hEta_N_R);
    sethist(hVerZ_P_C, hVerZ_P_R);
    sethist(hVerZ_N_C, hVerZ_N_R);

    gStyle->SetOptStat(0); // Disable stats box for cleaner plots

    // ==========================================
    // === EXISTING NOMINAL VS RAW PLOTS ========
    // ==========================================
    
    TCanvas* cPt = new TCanvas("cPt", "pT Comparison", 2000, 1000);
    cPt->Divide(2);
    cPt->cd(1);
    For1atatime(hPt_P_C, hPt_P_R, "p_{T} (GeV/c)", -1.0, 9.0);
    cPt->cd(2);
    For1atatime(hPt_N_C, hPt_N_R, "p_{T} (GeV/c)", -1.0, 9.0);
    //cPt->SaveAs("plots/EffCorrData_pT_PosNeg3.png");

    TCanvas* cEta = new TCanvas("cEta", "Eta Comparison", 2000, 1000);
    cEta->Divide(2);
    cEta->cd(1);
    For1atatime(hEta_P_C, hEta_P_R, "#eta", -1.0, 9.0);
    cEta->cd(2);
    For1atatime(hEta_N_C, hEta_N_R, "#eta", -1.0, 9.0);
    //cEta->SaveAs("plots/EffCorrData_Eta_PosNeg3.png");

    TCanvas* cVerZ = new TCanvas("cVerZ", "VerZ Comparison", 2000, 1000);
    cVerZ->Divide(2);
    cVerZ->cd(1);
    For1atatime(hVerZ_P_C, hVerZ_P_R, "Vz (cm)", -1.0, 9.0);
    cVerZ->cd(2);
    For1atatime(hVerZ_N_C, hVerZ_N_R, "Vz (cm)", -1.0, 9.0);
    //cVerZ->SaveAs("plots/EffCorrData_VerZ_PosNeg3.png");

    TCanvas* cMass = new TCanvas("cMass", "Final Mass Distribution", 1000, 1000);
    PlotInvMass(hInvMassCorr,hInvMassRaw, "m_{K_{S}^{0}K_{S}^{0}} (GeV/c^{2})", 1.0, 2.6, -2.0, 50.0);
    //cMass->SaveAs("plots/EffCorrData_K0K0Mass3.png");

    TCanvas* c4TOF = new TCanvas("c4TOF", "4 TOF Hits Comparison", 1000, 1000);
    h1D_Reco_InvMass_Corrected_4pi_4TOF->SetTitle("4 TOF Hits");
    PlotInvMass(h1D_Reco_InvMass_Corrected_4pi_4TOF, h1D_Reco_InvMass_Raw_4pi_4TOF, "m_{K_{S}^{0}K_{S}^{0}} (GeV/c^{2})", 1.0, 2.6, -2.0, 30.0);
    //c4TOF->SaveAs("plots/EffCorrData_4TOF.png");

    TCanvas* c3TOF = new TCanvas("c3TOF", "3 TOF Hits Comparison", 1000, 1000);
    h1D_Reco_InvMass_Corrected_4pi_3TOF->SetTitle("3 TOF Hits");
    PlotInvMass(h1D_Reco_InvMass_Corrected_4pi_3TOF, h1D_Reco_InvMass_Raw_4pi_3TOF, "m_{K_{S}^{0}K_{S}^{0}} (GeV/c^{2})", 1.0, 2.6, -2.0, 20.0);
    //c3TOF->SaveAs("plots/EffCorrData_3TOF.png");

    TCanvas* c2TOF = new TCanvas("c2TOF", "2 TOF Hits Comparison", 1000, 1000);
    h1D_Reco_InvMass_Corrected_4pi_2TOF->SetTitle("2 TOF Hits");
    PlotInvMass(h1D_Reco_InvMass_Corrected_4pi_2TOF, h1D_Reco_InvMass_Raw_4pi_2TOF, "m_{K_{S}^{0}K_{S}^{0}} (GeV/c^{2})", 1.0, 2.6, -2.0, 50.0);
    //c2TOF->SaveAs("plots/EffCorrData_2TOF.png");

    // ==========================================
    // === SYSTEMATICS PLOTS (Sys vs Nom) =======
    // ==========================================
    // Plotting standard range around 1.0 for the systematics ratio (e.g., 0.5 to 1.5)
    double sysMinY = 0.9;
    double sysMaxY = 1.1;

    TCanvas* cPt_Sys = new TCanvas("cPt_Sys", "pT Systematics", 2000, 1000);
    cPt_Sys->Divide(2);
    cPt_Sys->cd(1);
    PlotSys(hPt_P_C, hPt_P_C_17, hPt_P_C_22, "p_{T} (GeV/c)", "nHits=17", "nHits=22", sysMinY, sysMaxY);
    cPt_Sys->cd(2);
    PlotSys(hPt_N_C, hPt_N_C_17, hPt_N_C_22, "p_{T} (GeV/c)", "nHits=17", "nHits=22", sysMinY, sysMaxY);
    //cPt_Sys->SaveAs("plots/EffCorrSys_pT_PosNeg.png");

    TCanvas* cEta_Sys = new TCanvas("cEta_Sys", "Eta Systematics", 2000, 1000);
    cEta_Sys->Divide(2);
    cEta_Sys->cd(1);
    PlotSys(hEta_P_C, hEta_P_C_17, hEta_P_C_22, "#eta", "nHits=17", "nHits=22", sysMinY, sysMaxY);
    cEta_Sys->cd(2);
    PlotSys(hEta_N_C, hEta_N_C_17, hEta_N_C_22, "#eta", "nHits=17", "nHits=22", sysMinY, sysMaxY);
    //cEta_Sys->SaveAs("plots/EffCorrSys_Eta_PosNeg.png");

    TCanvas* cVerZ_Sys = new TCanvas("cVerZ_Sys", "VerZ Systematics", 2000, 1000);
    cVerZ_Sys->Divide(2);
    cVerZ_Sys->cd(1);
    PlotSys(hVerZ_P_C, hVerZ_P_C_17, hVerZ_P_C_22, "Vz (cm)", "nHits=17", "nHits=22", sysMinY, sysMaxY);
    cVerZ_Sys->cd(2);
    PlotSys(hVerZ_N_C, hVerZ_N_C_17, hVerZ_N_C_22, "Vz (cm)", "nHits=17", "nHits=22", sysMinY, sysMaxY);
    //cVerZ_Sys->SaveAs("plots/EffCorrSys_VerZ_PosNeg.png");*/

    TCanvas* cMass_Sys = new TCanvas("cMass_Sys", "Mass Systematics", 1000, 1000);
    PlotSys(hInvMassCorr, hInvMassCorr_17, hInvMassCorr_22, "m_{K_{S}^{0}K_{S}^{0}} (GeV/c^{2})", "nHits=17", "nHits=22", -1.0, 2.5);
    //cMass_Sys->SaveAs("plots/EffCorrSys_K0K0Mass.png");
    

    
    //EvaluatePtMissSystematics2(h1D_PtMiss_Corrected, nullptr, "1D");

    //CheckBinByBinPtMissFits(h2_PtMiss_Vs_Mass_Corrected);
    CheckBinByBinPtMissFits(h2_PtMiss_Vs_Mass);

    CheckTofBackgroundShapes(h1D_Reco_PtMiss_Raw_4pi_2TOF, h1D_Reco_PtMiss_Raw_4pi_3TOF, h1D_Reco_PtMiss_Raw_4pi_4TOF);

    //CheckTofBackgroundShapes(hPtMiss_2TOF_Raw, hPtMiss_3TOF_Raw, hPtMiss_4TOF_Raw);

    //CheckTofBackgroundShapes(hPtMiss_2TOF_Corr, hPtMiss_3TOF_Corr, hPtMiss_4TOF_Corr);

    EstimateTofBackgrounds(h1D_Reco_PtMiss_Raw_4pi_2TOF, h1D_Reco_PtMiss_Raw_4pi_3TOF, h1D_Reco_PtMiss_Raw_4pi_4TOF, "2D_Raw");
    // Inputs: rawTailCounts (10), estimatedBg (1.0), netSignal (16.0), label ("4-TOF")
    //CalculateScaleFactorSystematic(10.0, 1.0, 16.0, "4-TOF");
    // Run Closure Tests for 4 tof
    ClosureTest_TailStability(h2_PtMiss_Vs_Mass_4TOF_Corrected, "4TOF");
    ClosureTest_EfficiencyBias(h2_PtMiss_Vs_Mass_4TOF, h2_PtMiss_Vs_Mass_4TOF_Corrected, "4TOF");
    // Run Closure Tests for 3 tof
    ClosureTest_TailStability(h2_PtMiss_Vs_Mass_3TOF_Corrected, "3TOF");
    ClosureTest_EfficiencyBias(h2_PtMiss_Vs_Mass_3TOF, h2_PtMiss_Vs_Mass_3TOF_Corrected, "3TOF");
    // Run Closure Tests for 2 tof
    ClosureTest_TailStability(h2_PtMiss_Vs_Mass_2TOF_Corrected, "2TOF");
    ClosureTest_EfficiencyBias(h2_PtMiss_Vs_Mass_2TOF, h2_PtMiss_Vs_Mass_2TOF_Corrected, "2TOF");

    TH1D* hCleanCorr_2TOF = ExtractHybridCorrectedYield(h2_PtMiss_Vs_Mass_2TOF, h2_PtMiss_Vs_Mass_2TOF_Corrected, "2TOF");
    TH1D* hCleanCorr_3TOF = ExtractHybridCorrectedYield(h2_PtMiss_Vs_Mass_3TOF, h2_PtMiss_Vs_Mass_3TOF_Corrected, "3TOF");
    TH1D* hCleanCorr_4TOF = ExtractHybridCorrectedYield(h2_PtMiss_Vs_Mass_4TOF, h2_PtMiss_Vs_Mass_4TOF_Corrected, "4TOF");

    // Sum the properly weighted and subtracted yields!
    TH1D* hFinalCorrected = (TH1D*)hCleanCorr_4TOF->Clone("hFinalCorrected");
    hFinalCorrected->Add(hCleanCorr_3TOF);
    hFinalCorrected->Add(hCleanCorr_2TOF);
    hFinalCorrected->Scale(1.0/3.0); // Average over the 3 TOF categories to get a single corrected yield

    // 3. DRAWING AND SAVING
    TCanvas* cFinalSum = new TCanvas("cFinalSum", "Final Summed Yield", 1100, 900);
    gPad->SetFrameLineWidth(3);
    gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.13); 

    hFinalCorrected->SetMarkerStyle(20);
    hFinalCorrected->SetMarkerSize(1.5);
    hFinalCorrected->SetMarkerColor(kBlack);
    hFinalCorrected->SetLineColor(kBlack);
    hFinalCorrected->GetXaxis()->SetRangeUser(1.0, 2.6);
    hFinalCorrected->GetYaxis()->SetTitleOffset(0.97);
    hFinalCorrected->GetYaxis()->SetTitleSize(0.06);
    hFinalCorrected->GetYaxis()->SetLabelSize(0.06);
    hFinalCorrected->GetXaxis()->SetTitleOffset(0.9);
    hFinalCorrected->GetXaxis()->SetTitleSize(0.06);
    hFinalCorrected->GetXaxis()->SetLabelSize(0.05);


    // After your Add() commands:
    hFinalCorrected->SetTitle("; m_{K_{S}^{0}K_{S}^{0}} [GeV/c^{2}]; Corrected Yield");

    // Fit the final, CORRECTED dat
    TF1* finalFit = new TF1("finalFit", "gaus(0) + pol1(3)", 1.0, 2.6);
    finalFit->SetParameters(hFinalCorrected->GetMaximum(), 1.68, 0.08, 0, 0); 
    finalFit->SetLineColor(kRed);
    finalFit->SetLineWidth(2);
    finalFit->SetNpx(1000); // Fixes the gap in the red line!

    hFinalCorrected->Fit(finalFit, "R");
    hFinalCorrected->Draw("PE1");

    // Legend and text...
    TLegend* legSum = new TLegend(0.55, 0.70, 0.88, 0.88);
    legSum->SetBorderSize(0);
    legSum->SetTextSize(0.04);
    legSum->AddEntry(hFinalCorrected, "Corrected Yield", "pe");
    legSum->AddEntry(finalFit, "Signal Fit", "l");
    legSum->Draw();

    TLatex latexSum;
    latexSum.SetNDC();
    latexSum.SetTextSize(0.04);
    latexSum.DrawLatex(0.60, 0.60, Form("Mass = %.3f GeV/c^{2}", finalFit->GetParameter(1)));
    latexSum.DrawLatex(0.60, 0.55, Form("Width = %.3f GeV/c^{2}", finalFit->GetParameter(2)));

    cFinalSum->Update();
    cFinalSum->SaveAs("plots/Final_InvariantMass_Corrected.png");

    TH1D* hCleanCorr_2TOF_17 = ExtractHybridCorrectedYield(h2_PtMiss_Vs_Mass_2TOF_17, h2_PtMiss_Vs_Mass_2TOF_Corrected_17, "2TOF");
    TH1D* hCleanCorr_3TOF_17 = ExtractHybridCorrectedYield(h2_PtMiss_Vs_Mass_3TOF_17, h2_PtMiss_Vs_Mass_3TOF_Corrected_17, "3TOF");
    TH1D* hCleanCorr_4TOF_17 = ExtractHybridCorrectedYield(h2_PtMiss_Vs_Mass_4TOF_17, h2_PtMiss_Vs_Mass_4TOF_Corrected_17, "4TOF");

    // Sum the properly weighted and subtracted yields!
    TH1D* hFinalCorrected_17 = (TH1D*)hCleanCorr_4TOF_17->Clone("hFinalCorrected_17");
    hFinalCorrected_17->Add(hCleanCorr_3TOF_17);
    hFinalCorrected_17->Add(hCleanCorr_2TOF_17);
    hFinalCorrected_17->Scale(1.0/3.0); // Average over the 3 TOF categories to get a single corrected yield

    TH1D* hCleanCorr_2TOF_22 = ExtractHybridCorrectedYield(h2_PtMiss_Vs_Mass_2TOF_22, h2_PtMiss_Vs_Mass_2TOF_Corrected_22, "2TOF");
    TH1D* hCleanCorr_3TOF_22 = ExtractHybridCorrectedYield(h2_PtMiss_Vs_Mass_3TOF_22, h2_PtMiss_Vs_Mass_3TOF_Corrected_22, "3TOF");
    TH1D* hCleanCorr_4TOF_22 = ExtractHybridCorrectedYield(h2_PtMiss_Vs_Mass_4TOF_22, h2_PtMiss_Vs_Mass_4TOF_Corrected_22, "4TOF");

    // Sum the properly weighted and subtracted yields!
    TH1D* hFinalCorrected_22 = (TH1D*)hCleanCorr_4TOF_22->Clone("hFinalCorrected_22");
    hFinalCorrected_22->Add(hCleanCorr_3TOF_22);
    hFinalCorrected_22->Add(hCleanCorr_2TOF_22);
    hFinalCorrected_22->Scale(1.0/3.0); // Average over the 3 TOF categories to get a single corrected yield

    TCanvas* cMass_Sys_corr = new TCanvas("cMass_Sys_corr", "Mass Systematics corrected", 1000, 1000);
    PlotSys(hFinalCorrected, hFinalCorrected_17, hFinalCorrected_22, "m_{K_{S}^{0}K_{S}^{0}} (GeV/c^{2})", "nHits=17", "nHits=22", -1.0, 2.5);
    //cMass_Sys_corr->SaveAs("plots/Final_InvariantMass_Corrected_Sys.png");

    TCanvas* cMass_Raw_Sys = new TCanvas("cMass_Raw_Sys", "Mass Systematics raw", 1000, 1000);
    PlotSys(hInvMassRaw, hInvMassRaw_17, hInvMassRaw_22, "m_{K_{S}^{0}K_{S}^{0}} (GeV/c^{2})", "nHits=17", "nHits=22", -1.0, 2.5);
    //cMass_Raw_Sys->SaveAs("plots/Final_InvariantMass_Raw_Sys.png");

    // 1. Process Nominal Dataset
    //TH1D* hYield_Nominal = ProcessFinalYield(h2_PtMiss_Vs_Mass_2TOF, h2_PtMiss_Vs_Mass_3TOF, h2_PtMiss_Vs_Mass_4TOF, hInvMassRaw, hInvMassCorr, "Nominal");
    
    /*gStyle->SetPaintTextFormat(".0f"); // Set text format for 2 decimal places in the 2D histogram
   
    TCanvas* c2D = new TCanvas("c2D", "pTmiss vs Mass", 2000, 1000);
    c2D->Divide(2);
    c2D->cd(1);
    h2_PtMiss_Vs_Mass->SetTitle("Raw p_{T}^{miss} vs Mass");
    h2_PtMiss_Vs_Mass->Draw("TEXT COLZ");
    c2D->cd(2);
    h2_PtMiss_Vs_Mass_Corrected->SetTitle("Corrected p_{T}^{miss} vs Mass");
    h2_PtMiss_Vs_Mass_Corrected->Draw("TEXT COLZ");
    c2D->SaveAs("plots/pTmiss_vs_Mass.png");

    TCanvas* cMassRaw = new TCanvas("cMassRaw", "Raw Mass Distribution", 1100, 1000);
    hInvMassRaw->SetLineColor(kGreen);
    hInvMassRaw->SetLineWidth(3);
    hInvMassRaw->SetMarkerStyle(21);
    hInvMassRaw->SetMarkerColor(kGreen);
    hInvMassRaw->Draw("PE");
    hInvMassRaw->SetMarkerStyle(20);
    hInvMassRaw->SetMarkerColor(kBlack);
    hInvMassRaw->SetLineColor(kBlack);
    hInvMassRaw->GetXaxis()->SetRangeUser(0.9, 3.0);
    hInvMassRaw->GetXaxis()->SetTitle("m_{K_{S}^{0}K_{S}^{0}} [GeV/c^{2}]");
    hInvMassRaw->GetYaxis()->SetTitle("events");
    hInvMassRaw->GetYaxis()->SetTitleOffset(0.8);
    hInvMassRaw->GetXaxis()->SetTitleOffset(0.8);
    hInvMassRaw->GetXaxis()->SetTitleSize(0.05); //axis title size
    hInvMassRaw->GetYaxis()->SetTitleSize(0.06);

    hInvMassRaw->GetXaxis()->SetLabelSize(0.05);
    hInvMassRaw->GetYaxis()->SetLabelSize(0.05);
    
    hInvMassRaw->Draw("PE");

    TF1* fitFunc = new TF1("fitFunc", "gaus(0) + pol1(3)", 1.0, 2.6);
    fitFunc->SetParameters(15, 1.68, 0.08, 10, -2); 
    fitFunc->SetLineColor(kRed);
    hInvMassRaw->Fit(fitFunc, "R"); 

    TLegend* leg = new TLegend(0.55, 0.75, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.045);
    leg->AddEntry(hInvMassRaw, "Data", "pe");
    leg->AddEntry(fitFunc, "Fit: gaus + pol1", "l");
    leg->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.045);
    latex.DrawLatex(0.55, 0.68, Form("A = %.2f", fitFunc->GetParameter(0)));
    latex.DrawLatex(0.55, 0.62, Form("#mu = %.2f GeV/c^{2}", fitFunc->GetParameter(1)));
    latex.DrawLatex(0.55, 0.56, Form("#sigma = %.2f GeV/c^{2}", fitFunc->GetParameter(2)));

    cMassRaw->SaveAs("plots/EffCorrData_Raw_K0K0Mass3.png");*/

}