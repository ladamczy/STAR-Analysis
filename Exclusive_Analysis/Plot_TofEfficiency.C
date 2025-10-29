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
//
// Add these as global variables at the top of the file
double gFixedMean = 0.0;
double gFixedSigma = 0.0;

// Gaussian function definition
Double_t Gaussian(Double_t *x, Double_t *par) {
    return par[0] * TMath::Gaus(x[0], par[1], par[2]);
}

// Combined function with Gaussian + polynomial background
Double_t CombinedFunction(Double_t *x, Double_t *par) {
    Double_t signalValue = Gaussian(x, par);
    Double_t polyValue = par[3] + par[4] * x[0];
    return signalValue + polyValue;
}

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
    hist->GetXaxis()->SetRangeUser(0.45, 0.55);    
    int b_max = hist->FindBin(0.495);//GetMaximumBin();
    Double_t x_max = hist->GetBinCenter(b_max);
    Double_t y_max = hist->GetBinContent(b_max); 
    hist->GetYaxis()->SetRangeUser(x_max, y_max*2.9); //(0,100000);//
    hist->GetXaxis()->SetTitle("m_{#pi^{+}#pi^{-}} [GeV/c^{2}]");
    hist->GetYaxis()->SetTitle("N_{K}");

    hist->SetMarkerStyle(markerStyle);
    hist->SetMarkerColor(markerColor);
    hist->SetMarkerSize(markerSize);
    hist->SetLineColor(markerColor);

    if(title) hist->SetTitle(title);
}

void SetHistogramStyle2(TH1* hist, int markerColor) {
    hist->SetMinimum(0);
    hist->SetMaximum(1.1);
    hist->SetMarkerStyle(20);
    hist->SetMarkerSize(3);
    hist->SetLineWidth(3);
    hist->SetMarkerColor(markerColor);
    hist->SetLineColor(markerColor);
    hist->GetXaxis()->SetTitleOffset(0.9);
    hist->GetYaxis()->SetTitleOffset(0.93);
    hist->GetXaxis()->SetTitleSize(0.07);
    hist->GetYaxis()->SetTitleSize(0.07);
    hist->GetXaxis()->SetLabelSize(0.06);
    hist->GetYaxis()->SetLabelSize(0.06);
    hist->GetXaxis()->SetNdivisions(6);
    hist->GetYaxis()->SetNdivisions(6);
}

double CalculateYield(TF1* func, double xmin, double xmax, double binWidth) {
    const int nSteps = 1000;
    double step = (xmax - xmin) / nSteps;
    double sum = 0.0;
    
    for (int i = 0; i < nSteps; i++) {
        double x = xmin + (i + 0.5) * step; // Evaluate at bin center
        sum += func->Eval(x);
    }
    
    return sum * step / binWidth;
}

double FitHistogramAndGetYield(TH1* hist, TF1* &fitFunc, 
    const std::string& suffix, 
    bool draw = true, bool useTightRange = false,
    bool useFixedParams = false) {
    
    // Set fit range based on flag
    double f1 = useTightRange ? 0.485 : 0.48;
    double f2 = useTightRange ? 0.505 : 0.52;

    if (!fitFunc) {
        TString funcName = Form("fitFunc_%s", suffix.c_str());
        
        // Gaussian has 5 parameters (3 for Gaussian + 2 for linear background)
        int nParams = 5;
        
        fitFunc = new TF1(funcName, CombinedFunction, f1, f2, nParams);
        
        // Set parameter names and initial values for Gaussian
        fitFunc->SetParNames("Norm", "Mean", "Sigma", "Poly_a", "Poly_b");
        
        if (useFixedParams && gFixedMean > 0 && gFixedSigma > 0) {
            // Use the global parameters if fixed params are requested
            fitFunc->SetParameters(hist->GetMaximum(), gFixedMean, gFixedSigma, 0.0, 0.0);
            // Fix the parameters
            fitFunc->FixParameter(1, gFixedMean);
            fitFunc->FixParameter(2, gFixedSigma);
        } else {
            // Use initial guesses as before
            fitFunc->SetParameters(hist->GetMaximum(), 0.498, 0.005, 0.0, 0.0);
            fitFunc->SetParLimits(0, 0, hist->GetMaximum() * 5);
            fitFunc->SetParLimits(1, 0.49, 0.506);
            fitFunc->SetParLimits(2, 0.001, 0.01);
        }
    }

    TFitResultPtr fitResult = hist->Fit(fitFunc, draw ? "RS" : "RSQ+S");

    // Create signal component
    TString gaussName = Form("gaussian_%s", suffix.c_str());
    TF1* signalFunc = new TF1(gaussName, Gaussian, f1, f2, 3);
    for (int i = 0; i < 3; i++) {
        signalFunc->SetParameter(i, fitFunc->GetParameter(i));
    }
    
    signalFunc->SetLineColor(kBlue);
    signalFunc->SetLineStyle(2);

    // Create background component
    TString polyName = Form("poly_%s", suffix.c_str());
    TF1* poly = new TF1(polyName, "pol1", f1, f2);
    poly->SetParameter(0, fitFunc->GetParameter(3));
    poly->SetParameter(1, fitFunc->GetParameter(4));
    
    poly->SetLineColor(kOrange+4);
    poly->SetLineStyle(3);

    if (draw) {
        poly->Draw("same");
    }   

    double yield = signalFunc->Integral(f1, f2);

    return yield;
}

void DrawParameterBlock(double xStart, double yStart, 
    TF1* fitFunc, const char* title) {
    TLatex *text = new TLatex();
    text->SetNDC();
    text->SetTextSize(0.03);
    text->SetTextFont(42);
    text->SetTextColor(kBlue);

    // Format parameter strings
    auto formatParam = [](const char* name, double val, double err) {
        return Form("%s = %.4f #pm %.4f", name, val, err);
    };

    // Title
    text->DrawLatex(xStart, yStart, title);

    // Parameters
    double yPos = yStart - 0.03; // Start below title
    text->DrawLatex(xStart, yPos, formatParam("#mu", 
        fitFunc->GetParameter(1), fitFunc->GetParError(1)));
    yPos -= 0.03;
    
    text->DrawLatex(xStart, yPos, formatParam("#sigma", 
        fitFunc->GetParameter(2), fitFunc->GetParError(2)));
    yPos -= 0.03;

    text->DrawLatex(xStart, yPos, Form("#chi^{2}/NDF = %.3f", 
        fitFunc->GetChisquare()/fitFunc->GetNDF()));
}

double FitHistogramPairAndGetEfficiency(TH1* histWith, TH1* histWithout,
    TF1* &fitWith, TF1* &fitWithout,
    const std::string& suffix, 
    bool draw = true, bool useTightRange = false,
    bool useFixedParams = false) {
    
    // Set fit range based on flag
    double f1 = useTightRange ? 0.485 : 0.46;
    double f2 = useTightRange ? 0.505 : 0.53;  

    // First fit the "with TOF" histogram
    if (!fitWith) {
        TString funcNameWith = Form("fitFuncWith_%s", suffix.c_str());
        
        // Gaussian has 5 parameters
        int nParams = 5;
        
        fitWith = new TF1(funcNameWith, CombinedFunction, f1, f2, nParams);
        
        // Set parameter names and initial values for Gaussian
        fitWith->SetParNames("Norm", "Mean", "Sigma", "Poly_a", "Poly_b");
        
        if (useFixedParams && gFixedMean > 0 && gFixedSigma > 0) {
            // Use the global parameters if fixed params are requested
            fitWith->SetParameters(histWith->GetMaximum(), gFixedMean, gFixedSigma, 0.0, 0.0);
            // Fix the parameters
            fitWith->FixParameter(1, gFixedMean);
            fitWith->FixParameter(2, gFixedSigma);
        } else {
            // Use initial guesses as before
            fitWith->SetParameters(histWith->GetMaximum(), 0.498, 0.005, 0.0, 0.0);
            fitWith->SetParLimits(0, 0, histWith->GetMaximum() * 2);
            fitWith->SetParLimits(1, 0.49, 0.506);
            fitWith->SetParLimits(2, 0.001, 0.01);
        }
    }
    
    fitWith->SetLineColor(kRed);
    fitWith->SetLineWidth(3);
    TFitResultPtr fitResultWith = histWith->Fit(fitWith, draw ? "RS" : "RSQ");

    // Now fit the "without TOF" histogram, using the mean value from "with TOF" fit
    if (!fitWithout) {
        TString funcNameWithout = Form("fitFuncWithout_%s", suffix.c_str());
        
        // Gaussian has 5 parameters
        int nParams = 5;
        
        fitWithout = new TF1(funcNameWithout, CombinedFunction, f1, f2, nParams);
        
        // Set parameters
        fitWithout->SetParNames("Norm", "Mean", "Sigma", "Poly_a", "Poly_b");
        fitWithout->SetParameters(
            histWithout->GetMaximum(),       // Adjust normalization for this histogram
            fitWith->GetParameter(1),        // Use mean from "with TOF" fit
            fitWith->GetParameter(2),        // Initial sigma from "with TOF" fit
            0.0,                             // Reset background parameters
            0.0
        );

        fitWithout->SetParLimits(0, 0, histWithout->GetMaximum() * 2); // norm
        
        // Fix mean and sigma to the values from the "with TOF" fit
        fitWithout->FixParameter(1, fitWith->GetParameter(1));
        fitWithout->FixParameter(2, fitWith->GetParameter(2)); 
    }
    
    fitWithout->SetLineColor(kRed+2);
    fitWithout->SetLineWidth(3);
    TFitResultPtr fitResultWithout = histWithout->Fit(fitWithout, draw ? "RS" : "RSQ");

    // Create signal functions for yield calculations
    TString gaussNameWith = Form("gaussianWith_%s", suffix.c_str());
    TF1* signalWith = new TF1(gaussNameWith, Gaussian, f1, f2, 3);
    for (int i = 0; i < 3; i++) {
        signalWith->SetParameter(i, fitWith->GetParameter(i));
    }
    
    TString gaussNameWithout = Form("gaussianWithout_%s", suffix.c_str());
    TF1* signalWithout = new TF1(gaussNameWithout, Gaussian, f1, f2, 3);
    for (int i = 0; i < 3; i++) {
        signalWithout->SetParameter(i, fitWithout->GetParameter(i));
    }
    
    signalWith->SetLineColor(kBlue);
    signalWith->SetLineStyle(2);
    signalWithout->SetLineColor(kGreen+2);
    signalWithout->SetLineStyle(2);

    // Create background components
    TString polyNameWith = Form("polyWith_%s", suffix.c_str());
    TF1* polyWith = new TF1(polyNameWith, "pol1", f1, f2);
    
    TString polyNameWithout = Form("polyWithout_%s", suffix.c_str());
    TF1* polyWithout = new TF1(polyNameWithout, "pol1", f1, f2);
    
    polyWith->SetParameter(0, fitWith->GetParameter(3));
    polyWith->SetParameter(1, fitWith->GetParameter(4));
    polyWithout->SetParameter(0, fitWithout->GetParameter(3));
    polyWithout->SetParameter(1, fitWithout->GetParameter(4));
    
    polyWith->SetLineColor(kBlack);
    polyWith->SetLineStyle(3);
    polyWithout->SetLineColor(kBlack);
    polyWithout->SetLineStyle(3);

    if (draw) {
        polyWith->Draw("same");
        polyWithout->Draw("same");
    }

    double yieldWith = signalWith->Integral(f1, f2);
    double yieldWithout = signalWithout->Integral(f1, f2);

    return yieldWith / (yieldWith + yieldWithout); // Return the efficiency
}


// Create a function to extract and store the parameters from the tag-and-probe fit
void ExtractTagAndProbeParameters(TH1* histWith, TH1* histWithout) {
    TF1* fitWith = nullptr;
    TF1* fitWithout = nullptr;
    
    // Perform the fit (use true for draw, false for useTightRange)
    FitHistogramPairAndGetEfficiency(
        histWith, histWithout, 
        fitWith, fitWithout,
        "tagprobe", 
        true, false
    );
    
    // Store the mean and sigma globally
    gFixedMean = fitWith->GetParameter(1);
    gFixedSigma = fitWith->GetParameter(2);
    
    std::cout << "Extracted parameters from tag-and-probe fit:" << std::endl;
    std::cout << "Mean (mu): " << gFixedMean << std::endl;
    std::cout << "Sigma: " << gFixedSigma << std::endl;
}

void CreateSinglePlot(TCanvas* canvas, TH1* histWithout, TH1* histWith, 
    const char* title, bool useTightRange = false, 
    bool useFixedParams = false) {
    canvas->Clear();
    canvas->cd();
    canvas->SetFrameLineWidth(3);

    SetHistogramStyle(histWithout, title, 20, kGreen+3, 2.0);    
    SetHistogramStyle(histWith, nullptr, 21, kGreen-3, 2.0);

    histWithout->Draw("PE");
    histWith->Draw("PE same");

    // Use the fit function that fits "with TOF" first and then "without TOF"
    TF1 *fitWith = nullptr;
    TF1 *fitWithout = nullptr;
    double efficiency = FitHistogramPairAndGetEfficiency(
    histWith, histWithout, 
    fitWith, fitWithout,
    Form("bin_%d", canvas->GetNumber()), 
    true, useTightRange, useFixedParams
    );
    
    // Draw the fitted functions
    DrawParameterBlock(0.82, 0.82, fitWithout, "#bf{Without TOF}");
    DrawParameterBlock(0.82, 0.32, fitWith, "#bf{With TOF}");

    TLegend* legend = new TLegend(0.18, 0.75, 0.45, 0.89);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(histWithout, "Probe without TOF", "p");
    legend->AddEntry(histWith, "Probe with TOF", "p");
    legend->AddEntry(fitWith, "Fit (with TOF)", "l");
    legend->AddEntry(fitWithout, "Fit (without TOF)", "l");

    // Retrieve components with unique names
    TF1* polyWith = static_cast<TF1*>(gROOT->GetListOfFunctions()->
    FindObject(Form("polyWith_bin_%d", canvas->GetNumber())));
    TF1* polyWithout = static_cast<TF1*>(gROOT->GetListOfFunctions()->
    FindObject(Form("polyWithout_bin_%d", canvas->GetNumber())));

    if (polyWith && polyWithout) {
        legend->AddEntry(polyWith, "Poly1", "l");
    }

    legend->Draw();

    TLatex *effText = new TLatex();
    effText->SetNDC();
    effText->SetTextSize(0.05);
    effText->SetTextFont(42);
    effText->DrawLatex(0.38, 0.17, Form("#color[6]{Efficiency = %.1f%%}", efficiency * 100));
}


    
const char* etaTitles[] = {
    "#eta [-0.9 to -0.6]",
    "#eta [-0.6 to -0.3]",
    "#eta [-0.3 to 0.0]",
    "#eta [0.0 to 0.3]",
    "#eta [0.3 to 0.6]",
    "#eta [0.6 to 0.9]"
};

const char* ptTitles[] = {
    "p_{T} [0.2 to 0.3 GeV]",
    "p_{T} [0.3 to 0.4 GeV]",
    "p_{T} [0.4 to 0.5 GeV]",
    "p_{T} [0.5 to 0.6 GeV]",
    "p_{T} [0.6 to 0.8 GeV]",
    "p_{T} [0.8 to 1.2 GeV]"
};

const char* some[] = {"True pion matching",
    "Tag and Probe without true pion matching",
    "Tag and Probe for the data"};

const char* fitfunction[] = {"Gaussian"}; 
const char* fixvar[] = {"#sigma & #mu"};
double thresholdIndex = 0.3;  

// Latex for text
TLatex Tl;
Tl.SetTextAlign(10);
Tl.SetTextSize(0.05);

TLatex *fixedParamText = new TLatex();
fixedParamText->SetNDC();
fixedParamText->SetTextSize(0.035);
fixedParamText->SetTextFont(42);
fixedParamText->SetTextColor(kViolet+3);
 
double etaBins[7] = {-0.9, -0.6, -0.3, 0.0, 0.3, 0.6, 0.9};
double ptBins[7] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.2};

// Summary plots
TH1D* effVsEta = new TH1D("effVsEta", "TOF Efficiency vs #eta;#eta;Efficiency", 6, etaBins);
TH1D* effVsPt = new TH1D("effVsPt", "TOF Efficiency vs p_{T};p_{T} [GeV/c];Efficiency", 6, ptBins);

void Plot_MC(){   
    
    // Root style settings
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);
    
    // Open input file
    TFile* rootFile = TFile::Open("InputData/tofeff_withRP_0.3_v3.root");
    if (!rootFile || rootFile->IsZombie()) {
        std::cerr << "Error opening MC input file" << std::endl;
        return;
    }
    
    // Get histograms
    TH1* HistKaonMassProbeWithoutTof = (TH1D*)rootFile->Get("HistKaonMassProbeWithoutTof;1");
    TH1* HistKaonMassProbeWithTof = (TH1D*)rootFile->Get("HistKaonMassProbeWithTof;1");
     
    HistKaonMassProbeWithoutTof->Scale(1.0, "width");
    HistKaonMassProbeWithTof->Scale(1.0, "width");
    std::vector<TH1*> etaWithTof(6), etaWithoutTof(6);
    std::vector<TH1*> ptWithTof(6), ptWithoutTof(6);

    for (int i = 0; i < 6; i++) {
        etaWithTof[i] = static_cast<TH1*>(rootFile->Get(Form("HistKaonEtaProbeWithTofPosProbeMass%d;1", i)));
        etaWithoutTof[i] = static_cast<TH1*>(rootFile->Get(Form("HistKaonEtaProbeWithoutTofPosProbeMass%d;1", i)));
        ptWithTof[i] = static_cast<TH1*>(rootFile->Get(Form("HistKaonPtProbeWithTofPosProbeMass%d;1", i)));
        ptWithoutTof[i] = static_cast<TH1*>(rootFile->Get(Form("HistKaonPtProbeWithoutTofPosProbeMass%d;1", i)));
        
        if (!etaWithTof[i] || !etaWithoutTof[i] || !ptWithTof[i] || !ptWithoutTof[i]) {
            std::cerr << "Missing histogram for index " << i << std::endl;
            rootFile->Close();
            return;
        }

        etaWithTof[i]->Scale(1.0, "width");
        etaWithoutTof[i]->Scale(1.0, "width");
        ptWithTof[i]->Scale(1.0, "width");
        ptWithoutTof[i]->Scale(1.0, "width");

    }

    TCanvas* canvas = new TCanvas("canvas", "TOF Efficiency", 1500, 1048);
    canvas->SetLeftMargin(0.13);
    canvas->SetRightMargin(0.19);
    canvas->SetTopMargin(0.10);
    canvas->SetBottomMargin(0.15);


    // First fit the tag and probe histograms to extract the parameters
    //CreateSinglePlot(canvas, HistKaonMassProbeWithoutTof, HistKaonMassProbeWithTof, some[0], false);
    // tag and probe plots  
    CreateSinglePlot(canvas, HistKaonMassProbeWithoutTof, HistKaonMassProbeWithTof, some[1], false);
    //HistKaonMassProbeWithTof->GetYaxis()->SetRangeUser(-20, 150); 
    int hist_max2 = HistKaonMassProbeWithTof->GetMaximumBin();
    Double_t x_max2 = HistKaonMassProbeWithTof->GetBinCenter(hist_max2);
    Double_t y_max2 = HistKaonMassProbeWithTof->GetBinContent(hist_max2); 
    Tl.DrawLatex(0.51, y_max2 + 50, Form("#scale[0.6]{matching threshold < %.2f}", thresholdIndex));   
    fixedParamText->DrawLatex(0.18, 0.7, Form("#splitline{%s fixed}{Fit with: %s distribution}",fixvar[0],fitfunction[0]));
    canvas->Print("plots/TofEff_tagNprobe_thesis_MC_GS_0.3_v3.png");
    
    // Extract and store the parameters
    ExtractTagAndProbeParameters(HistKaonMassProbeWithTof, HistKaonMassProbeWithoutTof);

    // Eta plots
    canvas->Print("plots/TofEff_eta_thesis_MC_GS_0.3_v3.pdf[");
    for (int i = 0; i < 6; i++) {
        CreateSinglePlot(canvas, etaWithoutTof[i], etaWithTof[i], etaTitles[i], false, true);
        fixedParamText->DrawLatex(0.18, 0.7, Form("#splitline{%s fixed}{Fit with: %s distribution}",fixvar[0],fitfunction[0]));
        int hist_max = etaWithoutTof[i]->GetMaximumBin();
        Double_t x_max = etaWithoutTof[i]->GetBinCenter(hist_max);
        Double_t y_max = etaWithoutTof[i]->GetBinContent(hist_max);   
        Tl.DrawLatex(0.51, y_max + 10, Form("#scale[0.6]{matching threshold < %.2f}", thresholdIndex));
        canvas->Print("plots/TofEff_eta_thesis_MC_GS_0.3_v3.pdf");
    }
    canvas->Print("plots/TofEff_eta_thesis_MC_GS_0.3_v3.pdf]");

    // pT plots
    canvas->Print("plots/TofEff_pt_thesis_MC_GS_0.3_v3.pdf[");
    for (int i = 0; i < 6; i++) {
        CreateSinglePlot(canvas, ptWithoutTof[i], ptWithTof[i], ptTitles[i], false, true);
        fixedParamText->DrawLatex(0.18, 0.7, Form("#splitline{%s fixed}{Fit with: %s distribution}",fixvar[0],fitfunction[0]));
        int hist_max = ptWithoutTof[i]->GetMaximumBin();
        Double_t x_max = ptWithoutTof[i]->GetBinCenter(hist_max);
        Double_t y_max = ptWithoutTof[i]->GetBinContent(hist_max);   
        Tl.DrawLatex(0.51, y_max + 10, Form("#scale[0.6]{matching threshold < %.2f}", thresholdIndex));
        canvas->Print("plots/TofEff_pt_thesis_MC_GS_0.3_v3.pdf");
    }
    canvas->Print("plots/TofEff_pt_thesis_MC_GS_0.3_v3.pdf]");
   
    for (int i = 0; i < 6; i++) {
        TF1 *fitEtaWith = nullptr, *fitEtaWithout = nullptr;
        TF1 *fitPtWith = nullptr, *fitPtWithout = nullptr;

        // Calculate efficiency for eta bins with fixed parameters
        double effEta = FitHistogramPairAndGetEfficiency(
            etaWithTof[i], etaWithoutTof[i], 
            fitEtaWith, fitEtaWithout,
            Form("eta_%d", i), false, false, true // Pass true for useFixedParams
        );

        // Calculate efficiency for pt bins with fixed parameters
        double effPt = FitHistogramPairAndGetEfficiency(
            ptWithTof[i], ptWithoutTof[i], 
            fitPtWith, fitPtWithout,
            Form("pt_%d", i), false, false, true // Pass true for useFixedParams
        );
        
        double yieldEtaWith = FitHistogramAndGetYield(etaWithTof[i], fitEtaWith, Form("eta_w_%d", i), false);
        double yieldEtaWithout = FitHistogramAndGetYield(etaWithoutTof[i], fitEtaWithout, Form("eta_wo_%d", i), false);
        double yieldPtWith = FitHistogramAndGetYield(ptWithTof[i], fitPtWith, Form("pt_w_%d", i), false);
        double yieldPtWithout = FitHistogramAndGetYield(ptWithoutTof[i], fitPtWithout, Form("pt_wo_%d", i), false);

        // Set bin content and errors for efficiency histograms
        effVsEta->SetBinContent(i+1, effEta);
        effVsEta->SetBinError(i+1, sqrt(effEta*(1-effEta)/(yieldEtaWith + yieldEtaWithout)));
        effVsPt->SetBinContent(i+1, effPt);
        effVsPt->SetBinError(i+1, sqrt(effPt*(1-effPt)/(yieldPtWith + yieldPtWithout)));

        // Clean up
        delete fitEtaWith;
        delete fitEtaWithout;
        delete fitPtWith;
        delete fitPtWithout;
    }
        
   
    delete canvas;
    //delete effVsEta;
    //delete effVsPt;
    rootFile->Close();
}

// Summary plots
TH1D* effVsEta_data = new TH1D("effVsEta_data", "TOF Efficiency vs #eta;#eta;Efficiency", 6, etaBins);
TH1D* effVsPt_data = new TH1D("effVsPt_data", "TOF Efficiency vs p_{T};p_{T} [GeV/c];Efficiency", 6, ptBins);

void Plot_data(){
        
    // Root style settings
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

   
    // Open input file
    TFile* rootFile = TFile::Open("InputData/tofeff_data_0.3_v3.root");
    if (!rootFile || rootFile->IsZombie()) {
        std::cerr << "Error opening MC input file" << std::endl;
        return;
    }
    
    // Get histograms
    TH1* HistKaonMassProbeWithoutTof_data = (TH1D*)rootFile->Get("HistKaonMassProbeWithoutTof;1");
    TH1* HistKaonMassProbeWithTof_data = (TH1D*)rootFile->Get("HistKaonMassProbeWithTof;1");
     
    HistKaonMassProbeWithoutTof_data->Scale(1.0, "width");
    HistKaonMassProbeWithTof_data->Scale(1.0, "width");

    std::vector<TH1*> etaWithTof_data(6), etaWithoutTof_data(6);
    std::vector<TH1*> ptWithTof_data(6), ptWithoutTof_data(6);

    for (int i = 0; i < 6; i++) {
        etaWithTof_data[i] = static_cast<TH1*>(rootFile->Get(Form("HistKaonEtaProbeWithTofPosProbeMass%d;1", i)));
        etaWithoutTof_data[i] = static_cast<TH1*>(rootFile->Get(Form("HistKaonEtaProbeWithoutTofPosProbeMass%d;1", i)));
        ptWithTof_data[i] = static_cast<TH1*>(rootFile->Get(Form("HistKaonPtProbeWithTofPosProbeMass%d;1", i)));
        ptWithoutTof_data[i] = static_cast<TH1*>(rootFile->Get(Form("HistKaonPtProbeWithoutTofPosProbeMass%d;1", i)));
        
        if (!etaWithTof_data[i] || !etaWithoutTof_data[i] || !ptWithTof_data[i] || !ptWithoutTof_data[i]) {
            std::cerr << "Missing histogram for index " << i << std::endl;
            rootFile->Close();
            return;
        }

        etaWithTof_data[i]->Scale(1.0, "width");
        etaWithoutTof_data[i]->Scale(1.0, "width");
        ptWithTof_data[i]->Scale(1.0, "width");
        ptWithoutTof_data[i]->Scale(1.0, "width");

    }

    TCanvas* canvas = new TCanvas("canvas", "TOF Efficiency", 1500, 1048);
    canvas->SetLeftMargin(0.13);
    canvas->SetRightMargin(0.19);
    canvas->SetTopMargin(0.10);
    canvas->SetBottomMargin(0.15);


    // tag and probe plots  
    CreateSinglePlot(canvas, HistKaonMassProbeWithoutTof_data, HistKaonMassProbeWithTof_data, some[2], false);
    //HistKaonMassProbeWithTof_data->GetYaxis()->SetRangeUser(-20, 150); 
    int hist_max2 = HistKaonMassProbeWithTof_data->GetMaximumBin();
    Double_t x_max2 = HistKaonMassProbeWithTof_data->GetBinCenter(hist_max2);
    Double_t y_max2 = HistKaonMassProbeWithTof_data->GetBinContent(hist_max2); 
    Tl.DrawLatex(0.51, y_max2 + 50, Form("#scale[0.6]{matching threshold < %.2f}", thresholdIndex));   
    fixedParamText->DrawLatex(0.18, 0.7, Form("#splitline{%s fixed}{Fit with: %s distribution}",fixvar[0],fitfunction[0]));
    canvas->Print("plots/TofEff_tagNprobe_thesis_data_GS_0.3_v3.png");
    
    // Extract and store the parameters
    ExtractTagAndProbeParameters(HistKaonMassProbeWithTof_data, HistKaonMassProbeWithoutTof_data);

    // Eta plots
    canvas->Print("plots/TofEff_eta_thesis_data_GS_0.3_v3.pdf[");
    for (int i = 0; i < 6; i++) {
        CreateSinglePlot(canvas, etaWithoutTof_data[i], etaWithTof_data[i], etaTitles[i], false, true);
        fixedParamText->DrawLatex(0.18, 0.7, Form("#splitline{%s fixed}{Fit with: %s distribution}",fixvar[0],fitfunction[0]));
        int hist_max = etaWithoutTof_data[i]->GetMaximumBin();
        Double_t x_max = etaWithoutTof_data[i]->GetBinCenter(hist_max);
        Double_t y_max = etaWithoutTof_data[i]->GetBinContent(hist_max);   
        Tl.DrawLatex(0.51, y_max + 10, Form("#scale[0.6]{matching threshold < %.2f}", thresholdIndex));
        canvas->Print("plots/TofEff_eta_thesis_data_GS_0.3_v3.pdf");
    }
    canvas->Print("plots/TofEff_eta_thesis_data_GS_0.3_v3.pdf]");

    // pT plots
    canvas->Print("plots/TofEff_pt_thesis_data_GS_0.3_v3.pdf[");
    for (int i = 0; i < 6; i++) {
        CreateSinglePlot(canvas, ptWithoutTof_data[i], ptWithTof_data[i], ptTitles[i], false, true);
        fixedParamText->DrawLatex(0.18, 0.7, Form("#splitline{%s fixed}{Fit with: %s distribution}",fixvar[0],fitfunction[0]));
        int hist_max = ptWithoutTof_data[i]->GetMaximumBin();
        Double_t x_max = ptWithoutTof_data[i]->GetBinCenter(hist_max);
        Double_t y_max = ptWithoutTof_data[i]->GetBinContent(hist_max);   
        Tl.DrawLatex(0.51, y_max + 10, Form("#scale[0.6]{matching threshold < %.2f}", thresholdIndex));
        canvas->Print("plots/TofEff_pt_thesis_data_GS_0.3_v3.pdf");
    }
    canvas->Print("plots/TofEff_pt_thesis_data_GS_0.3_v3.pdf]");
   
    for (int i = 0; i < 6; i++) {
        TF1 *fitEtaWith = nullptr, *fitEtaWithout = nullptr;
        TF1 *fitPtWith = nullptr, *fitPtWithout = nullptr;

        // Calculate efficiency for eta bins with fixed parameters
        double effEta = FitHistogramPairAndGetEfficiency(
            etaWithTof_data[i], etaWithoutTof_data[i], 
            fitEtaWith, fitEtaWithout,
            Form("eta_%d", i), false, false, true // Pass true for useFixedParams
        );

        // Calculate efficiency for pt bins with fixed parameters
        double effPt = FitHistogramPairAndGetEfficiency(
            ptWithTof_data[i], ptWithoutTof_data[i], 
            fitPtWith, fitPtWithout,
            Form("pt_%d", i), false, false, true // Pass true for useFixedParams
        );
        
        double yieldEtaWith = FitHistogramAndGetYield(etaWithTof_data[i], fitEtaWith, Form("eta_w_%d", i), false);
        double yieldEtaWithout = FitHistogramAndGetYield(etaWithoutTof_data[i], fitEtaWithout, Form("eta_wo_%d", i), false);
        double yieldPtWith = FitHistogramAndGetYield(ptWithTof_data[i], fitPtWith, Form("pt_w_%d", i), false);
        double yieldPtWithout = FitHistogramAndGetYield(ptWithoutTof_data[i], fitPtWithout, Form("pt_wo_%d", i), false);

        // Set bin content and errors for efficiency histograms
        effVsEta_data->SetBinContent(i+1, effEta);
        effVsEta_data->SetBinError(i+1, sqrt(effEta*(1-effEta)/(yieldEtaWith + yieldEtaWithout)));
        effVsPt_data->SetBinContent(i+1, effPt);
        effVsPt_data->SetBinError(i+1, sqrt(effPt*(1-effPt)/(yieldPtWith + yieldPtWithout)));

        // Clean up
        delete fitEtaWith;
        delete fitEtaWithout;
        delete fitPtWith;
        delete fitPtWithout;
    }
        
   
    delete canvas;
    //delete effVsEta_data;
    //delete effVsPt_data;
    rootFile->Close();
}

  
TH1D* effVsEta_true = new TH1D("effVsEta_true", "TOF Efficiency vs #eta;#eta;Efficiency", 6, etaBins);
TH1D* effVsPt_true = new TH1D("effVsPt_true", "TOF Efficiency vs p_{T};p_{T} [GeV/c];Efficiency", 6, ptBins);

void Plot_true(){
       
    // Root style settings
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);
     
    // Open input file
    TFile* rootFile = TFile::Open("InputData/tofeff_withRP_0.3_v3.root");
    if (!rootFile || rootFile->IsZombie()) {
        std::cerr << "Error opening MC input file" << std::endl;
        return;
    }
    
    // Get histograms
    TH1* HistKaonMassProbeWithoutTof_true = (TH1D*)rootFile->Get("HistKaonMassProbeWithoutTof_TruePions;1");
    TH1* HistKaonMassProbeWithTof_true = (TH1D*)rootFile->Get("HistKaonMassProbeWithTof_TruePions;1");
     
    HistKaonMassProbeWithoutTof_true->Scale(1.0, "width");
    HistKaonMassProbeWithTof_true->Scale(1.0, "width");

    std::vector<TH1*> etaWithTof_true(6), etaWithoutTof_true(6);
    std::vector<TH1*> ptWithTof_true(6), ptWithoutTof_true(6);

    for (int i = 0; i < 6; i++) {
        etaWithTof_true[i] = static_cast<TH1*>(rootFile->Get(Form("HistKaonEtaTruePionWithTofMass%d;1", i)));
        etaWithoutTof_true[i] = static_cast<TH1*>(rootFile->Get(Form("HistKaonEtaTruePionWithoutTofMass%d;1", i)));
        ptWithTof_true[i] = static_cast<TH1*>(rootFile->Get(Form("HistKaonPtTruePionWithTofMass%d;1", i)));
        ptWithoutTof_true[i] = static_cast<TH1*>(rootFile->Get(Form("HistKaonPtTruePionWithoutTofMass%d;1", i)));
        
        if (!etaWithTof_true[i] || !etaWithoutTof_true[i] || !ptWithTof_true[i] || !ptWithoutTof_true[i]) {
            std::cerr << "Missing histogram for index " << i << std::endl;
            rootFile->Close();
            return;
        }

        etaWithTof_true[i]->Scale(1.0, "width");
        etaWithoutTof_true[i]->Scale(1.0, "width");
        ptWithTof_true[i]->Scale(1.0, "width");
        ptWithoutTof_true[i]->Scale(1.0, "width");

    }

    TCanvas* canvas = new TCanvas("canvas", "TOF Efficiency", 1500, 1048);
    canvas->SetLeftMargin(0.13);
    canvas->SetRightMargin(0.19);
    canvas->SetTopMargin(0.10);
    canvas->SetBottomMargin(0.15);


    // First fit the tag and probe histograms to extract the parameters
    //CreateSinglePlot(canvas, HistKaonMassProbeWithoutTof_true, HistKaonMassProbeWithTof_true, some[0], false);
    // tag and probe plots  
    CreateSinglePlot(canvas, HistKaonMassProbeWithoutTof_true, HistKaonMassProbeWithTof_true, some[0], true);
    //HistKaonMassProbeWithTof_true->GetYaxis()->SetRangeUser(-20, 150); 
    int hist_max2 = HistKaonMassProbeWithTof_true->GetMaximumBin();
    Double_t x_max2 = HistKaonMassProbeWithTof_true->GetBinCenter(hist_max2);
    Double_t y_max2 = HistKaonMassProbeWithTof_true->GetBinContent(hist_max2); 
    Tl.DrawLatex(0.51, y_max2 + 50, Form("#scale[0.6]{matching threshold < %.2f}", thresholdIndex));   
    fixedParamText->DrawLatex(0.18, 0.7, Form("#splitline{%s fixed}{Fit with: %s distribution}",fixvar[0],fitfunction[0]));
    canvas->Print("plots/TofEff_tagNprobe_thesis_true_GS_0.3_v3.png");
    
    // Extract and store the parameters
    ExtractTagAndProbeParameters(HistKaonMassProbeWithTof_true, HistKaonMassProbeWithoutTof_true);

    // Eta plots
    canvas->Print("plots/TofEff_eta_thesis_true_GS_0.3_v3.pdf[");
    for (int i = 0; i < 6; i++) {
        CreateSinglePlot(canvas, etaWithoutTof_true[i], etaWithTof_true[i], etaTitles[i], true, true);
        fixedParamText->DrawLatex(0.18, 0.7, Form("#splitline{%s fixed}{Fit with: %s distribution}",fixvar[0],fitfunction[0]));
        int hist_max = etaWithoutTof_true[i]->GetMaximumBin();
        Double_t x_max = etaWithoutTof_true[i]->GetBinCenter(hist_max);
        Double_t y_max = etaWithoutTof_true[i]->GetBinContent(hist_max);   
        Tl.DrawLatex(0.51, y_max + 10, Form("#scale[0.6]{matching threshold < %.2f}", thresholdIndex));
        canvas->Print("plots/TofEff_eta_thesis_true_GS_0.3_v3.pdf");
    }
    canvas->Print("plots/TofEff_eta_thesis_true_GS_0.3_v3.pdf]");

    // pT plots
    canvas->Print("plots/TofEff_pt_thesis_true_GS_0.3_v3.pdf[");
    for (int i = 0; i < 6; i++) {
        CreateSinglePlot(canvas, ptWithoutTof_true[i], ptWithTof_true[i], ptTitles[i], true, true);
        fixedParamText->DrawLatex(0.18, 0.7, Form("#splitline{%s fixed}{Fit with: %s distribution}",fixvar[0],fitfunction[0]));
        int hist_max = ptWithoutTof_true[i]->GetMaximumBin();
        Double_t x_max = ptWithoutTof_true[i]->GetBinCenter(hist_max);
        Double_t y_max = ptWithoutTof_true[i]->GetBinContent(hist_max);   
        Tl.DrawLatex(0.51, y_max + 10, Form("#scale[0.6]{matching threshold < %.2f}", thresholdIndex));
        canvas->Print("plots/TofEff_pt_thesis_true_GS_0.3_v3.pdf");
    }
    canvas->Print("plots/TofEff_pt_thesis_true_GS_0.3_v3.pdf]");
   
    for (int i = 0; i < 6; i++) {
        TF1 *fitEtaWith = nullptr, *fitEtaWithout = nullptr;
        TF1 *fitPtWith = nullptr, *fitPtWithout = nullptr;

        // Calculate efficiency for eta bins with fixed parameters
        double effEta = FitHistogramPairAndGetEfficiency(
            etaWithTof_true[i], etaWithoutTof_true[i], 
            fitEtaWith, fitEtaWithout,
            Form("eta_%d", i), false, false, true // Pass true for useFixedParams
        );

        // Calculate efficiency for pt bins with fixed parameters
        double effPt = FitHistogramPairAndGetEfficiency(
            ptWithTof_true[i], ptWithoutTof_true[i], 
            fitPtWith, fitPtWithout,
            Form("pt_%d", i), false, false, true // Pass true for useFixedParams
        );
        
        double yieldEtaWith = FitHistogramAndGetYield(etaWithTof_true[i], fitEtaWith, Form("eta_w_%d", i), false);
        double yieldEtaWithout = FitHistogramAndGetYield(etaWithoutTof_true[i], fitEtaWithout, Form("eta_wo_%d", i), false);
        double yieldPtWith = FitHistogramAndGetYield(ptWithTof_true[i], fitPtWith, Form("pt_w_%d", i), false);
        double yieldPtWithout = FitHistogramAndGetYield(ptWithoutTof_true[i], fitPtWithout, Form("pt_wo_%d", i), false);

        // Set bin content and errors for efficiency histograms
        effVsEta_true->SetBinContent(i+1, effEta);
        effVsEta_true->SetBinError(i+1, sqrt(effEta*(1-effEta)/(yieldEtaWith + yieldEtaWithout)));
        effVsPt_true->SetBinContent(i+1, effPt);
        effVsPt_true->SetBinError(i+1, sqrt(effPt*(1-effPt)/(yieldPtWith + yieldPtWithout)));

        // Clean up
        delete fitEtaWith;
        delete fitEtaWithout;
        delete fitPtWith;
        delete fitPtWithout;
    }
        
   
    delete canvas;
    //delete effVsEta_true;
    //delete effVsPt_true;
    rootFile->Close();
}


void Plot_TofEfficiency() {

    Plot_MC();

    Plot_data();

    Plot_true();

    TCanvas* canvas = new TCanvas("canvas", "TOF Efficiency", 1500, 1048);
    // Draw summary plots
    canvas->Clear();
    canvas->SetLeftMargin(0.15);
    canvas->SetRightMargin(0.10);
    canvas->SetTopMargin(0.10);
    canvas->SetBottomMargin(0.15);
    canvas->SetFrameLineWidth(6);
    
    SetHistogramStyle2(effVsEta, kBlue);
    SetHistogramStyle2(effVsEta_true, kRed);
    SetHistogramStyle2(effVsEta_data, kGreen);
    SetHistogramStyle2(effVsPt, kBlue);
    SetHistogramStyle2(effVsPt_true, kRed);
    SetHistogramStyle2(effVsPt_data, kGreen);

    TLegend* legend = new TLegend(0.72, 0.25, 0.9, 0.45);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(effVsEta, "MC", "p");
    legend->AddEntry(effVsEta_true, "True", "p");
    legend->AddEntry(effVsEta_data, "Data", "p");

    effVsEta->Draw("PE");
    effVsEta_data->Draw("PE same");
    effVsEta_true->Draw("PE same");
    legend->Draw();
    Tl.DrawLatex(-0.8,0.1, Form("#splitline{%s fixed}{Fit with: %s distribution}",fixvar[0],fitfunction[0]));
    canvas->Print("plots/TofEff_summary_eta_thesis_MC_GS_0.3_v3.png");

    canvas->Clear();
    SetHistogramStyle2(effVsPt, kBlue);
    effVsPt->Draw("PE");
    effVsPt_data->Draw("PE same");
    effVsPt_true->Draw("PE same");
    legend->Draw();
    
    Tl.DrawLatex(0.25,0.1, Form("#splitline{%s fixed}{Fit with: %s distribution}",fixvar[0],fitfunction[0]));
    canvas->Print("plots/TofEff_summary_pt_thesis_MC_GS_0.3_v3.png");

    
}