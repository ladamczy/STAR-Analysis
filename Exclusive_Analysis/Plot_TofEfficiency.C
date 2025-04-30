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
// Crystal Ball function definition
Double_t CrystalBall(Double_t *x, Double_t *par) {
    Double_t t = (x[0] - par[1]) / par[2];
    Double_t absAlpha = fabs(par[3]);
    if (par[3] < 0) t = -t;
    
    if (t >= -absAlpha) {
        return par[0] * exp(-0.5 * t * t);
    } else {
        Double_t A = pow(par[4]/absAlpha, par[4]) * exp(-0.5 * absAlpha*absAlpha);
        Double_t B = par[4]/absAlpha - absAlpha;
        return par[0] * A * pow(B - t, -par[4]);
    }
}

Double_t Gaussian(Double_t *x, Double_t *par) {
    return par[0] * TMath::Gaus(x[0], par[1], par[2]);
}

/*
Double_t VoigtFunction(Double_t *x, Double_t *par) {
    Double_t sigma = par[2]; // Gaussian width (detector resolution)
    Double_t lg = par[3];    // Lorentzian width (natural width)
    return par[0] * TMath::Voigt(x[0]-par[1], sigma, lg);
}*/

// 
Double_t VoigtFunction(Double_t *x, Double_t *par) {
    Double_t sigma = par[2]; // Gaussian width (detector resolution)
    Double_t lg = par[3];    // Lorentzian width (natural width)
    
    // Avoid numerical issues with tiny values
    if (sigma < 1e-6) sigma = 1e-6;
    if (lg < 1e-6) lg = 1e-6;
    
    return par[0] * TMath::Voigt(x[0]-par[1], sigma, lg);
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
    //hist->GetYaxis()->SetRangeUser(10, 2000);
    int b_max = hist->GetMaximumBin();
    Double_t x_max = hist->GetBinCenter(b_max);
    Double_t y_max = hist->GetBinContent(b_max); 
    hist->GetYaxis()->SetRangeUser(x_max-2, y_max*1.8);
    //hist->Scale(1.0/1000.0);
    hist->GetXaxis()->SetTitle("m_{#pi^{+}#pi^{-}} [GeV/c^{2}]");
    //hist->GetYaxis()->SetTitle("N_{K}*10^{3}");//
    hist->GetYaxis()->SetTitle("N_{K}");//
    

    hist->SetMarkerStyle(markerStyle);
    hist->SetMarkerColor(markerColor);
    hist->SetMarkerSize(markerSize);
    hist->SetLineColor(markerColor);

    if(title) hist->SetTitle(title);
}

void SetHistogramStyle2(TH1* hist, int markerColor) {
    hist->SetMinimum(0);
    hist->SetMaximum(1);
    hist->SetMarkerStyle(20);
    hist->SetMarkerSize(3);hist->SetLineWidth(3);
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

// Add an enum for fit function types
enum FitFunctionType {
    kCrystalBall,
    kGaussian,
    kVoigt
};

// Modified combined function that uses a global variable to determine which function to use
FitFunctionType gSelectedFitFunc = kCrystalBall; // Default to Crystal Ball

Double_t CombinedFunction(Double_t *x, Double_t *par) {
    Double_t signalValue = 0.0;
    Double_t polyValue = 0.0;
    
    // Select the appropriate signal function based on the global setting
    switch(gSelectedFitFunc) {
        case kCrystalBall:
            signalValue = CrystalBall(x, par);
            polyValue = par[5] + par[6] * x[0];    
            break;
        case kGaussian:
            signalValue = Gaussian(x, par);
            polyValue = par[3] + par[4] * x[0];
            break;
        case kVoigt:
            signalValue = VoigtFunction(x, par);
            polyValue = par[4] + par[5] * x[0];
            break;
    }
    
    return signalValue + polyValue;
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
// Modified FitHistogramAndGetYield function to handle different fit functions
double FitHistogramAndGetYield(TH1* hist, TF1* &fitFunc, 
    const std::string& suffix, 
    bool draw = true, bool useTightRange = false) {
    
    // Set fit range based on flag
    double f1 = useTightRange ? 0.485 : 0.48;
    double f2 = useTightRange ? 0.505 : 0.52;

    if (!fitFunc) {
        TString funcName = Form("fitFunc_%s", suffix.c_str());
        
        // Number of parameters depends on the fit function type
        int nParams = 0;
        switch(gSelectedFitFunc) {
            case kCrystalBall: nParams = 7; break;
            case kGaussian: nParams = 5; break;
            case kVoigt: nParams = 6; break;
        }
        
        fitFunc = new TF1(funcName, CombinedFunction,  f1, f2, nParams);
        
        // Set parameter names and initial values based on fit function type
        switch(gSelectedFitFunc) {
            case kCrystalBall:
                fitFunc->SetParNames("Norm", "Mean", "Sigma", "Alpha", "n", "Poly_a", "Poly_b");
                fitFunc->SetParameters(hist->GetMaximum(), 0.485, 0.0005, 1.0, 5.0, 0.0, 0.0);
                
                fitFunc->SetParLimits(0, 0, hist->GetMaximum() * 5);
                fitFunc->SetParLimits(1, 0.485, 0.515);
                fitFunc->SetParLimits(2, 0.0005, 0.015);
                fitFunc->SetParLimits(3, 0.05, 20.0);
                fitFunc->SetParLimits(4, 0.50, 20.0);
                break;
                
            case kGaussian:
                fitFunc->SetParNames("Norm", "Mean", "Sigma", "Poly_a", "Poly_b");
                fitFunc->SetParameters(hist->GetMaximum(), 0.498, 0.005, 0.0, 0.0);
                
                fitFunc->SetParLimits(0, 0, hist->GetMaximum() * 5);
                fitFunc->SetParLimits(1, 0.49, 0.506);
                fitFunc->SetParLimits(2, 0.001, 0.01);
                break;
                
            case kVoigt:
                fitFunc->SetParNames("Norm", "Mean", "Sigma", "Gamma", "Poly_a", "Poly_b");
                fitFunc->SetParameters(hist->GetMaximum(), 0.498, 0.005, 0.001, 0.0, 0.0);
                
                fitFunc->SetParLimits(0, 0, hist->GetMaximum() * 5);
                fitFunc->SetParLimits(1, 0.49, 0.506);
                fitFunc->SetParLimits(2, 0.001, 0.01);
                fitFunc->SetParLimits(3, 0.0001, 0.005);
                break;
        }
    }

    TFitResultPtr fitResult = hist->Fit(fitFunc, draw ? "RS" : "RSQ+S");

    // Create signal component
    TF1* signalFunc = nullptr;
    switch(gSelectedFitFunc) {
        case kCrystalBall: {
            TString cbName = Form("crystalBall_%s", suffix.c_str());
            signalFunc = new TF1(cbName, CrystalBall,  f1, f2, 5);
            for (int i = 0; i < 5; i++) {
                signalFunc->SetParameter(i, fitFunc->GetParameter(i));
            }
            break;
        }
        case kGaussian: {
            TString gaussName = Form("gaussian_%s", suffix.c_str());
            signalFunc = new TF1(gaussName, Gaussian,  f1, f2, 3);
            for (int i = 0; i < 3; i++) {
                signalFunc->SetParameter(i, fitFunc->GetParameter(i));
            }
            break;
        }
        case kVoigt: {
            TString voigtName = Form("voigt_%s", suffix.c_str());
            signalFunc = new TF1(voigtName, VoigtFunction,  f1, f2, 4);
            for (int i = 0; i < 4; i++) {
                signalFunc->SetParameter(i, fitFunc->GetParameter(i));
            }
            break;
        }
    }
    
    if (signalFunc) {
        signalFunc->SetLineColor(kBlue);
        signalFunc->SetLineStyle(2);
    }

    // Create background component
    TString polyName = Form("poly_%s", suffix.c_str());
    TF1* poly = new TF1(polyName, "pol1",  f1, f2);
    
    switch(gSelectedFitFunc) {
        case kCrystalBall:
            poly->SetParameter(0, fitFunc->GetParameter(5));
            poly->SetParameter(1, fitFunc->GetParameter(6));
            break;
        case kGaussian:
            poly->SetParameter(0, fitFunc->GetParameter(3));
            poly->SetParameter(1, fitFunc->GetParameter(4));
            break;
        case kVoigt:
            poly->SetParameter(0, fitFunc->GetParameter(4));
            poly->SetParameter(1, fitFunc->GetParameter(5));
            break;
    }
    
    poly->SetLineColor(kOrange+4);
    poly->SetLineStyle(3);

    if (draw) {
        //signalFunc->Draw("same");
        poly->Draw("same");
    }   

    double yield = signalFunc->Integral(f1, f2);// / hist->GetBinWidth(1);

    return yield;
}

// Modified function to draw parameter blocks based on fit function type
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

    // Draw different parameters based on the fit function type
    switch(gSelectedFitFunc) {
        case kCrystalBall:
            text->DrawLatex(xStart, yPos, formatParam("#sigma", 
                fitFunc->GetParameter(2), fitFunc->GetParError(2)));
            yPos -= 0.03;
            
            text->DrawLatex(xStart, yPos, formatParam("#alpha", 
                fitFunc->GetParameter(3), fitFunc->GetParError(3)));
            yPos -= 0.03;
            
            text->DrawLatex(xStart, yPos, formatParam("n", 
                fitFunc->GetParameter(4), fitFunc->GetParError(4)));
            yPos -= 0.03;
            break;
            
        case kGaussian:
            text->DrawLatex(xStart, yPos, formatParam("#sigma", 
                fitFunc->GetParameter(2), fitFunc->GetParError(2)));
            yPos -= 0.03;
            break;
            
        case kVoigt:
            text->DrawLatex(xStart, yPos, formatParam("#sigma", 
                fitFunc->GetParameter(2), fitFunc->GetParError(2)));
            yPos -= 0.03;
            
            text->DrawLatex(xStart, yPos, formatParam("#gamma", 
                fitFunc->GetParameter(3), fitFunc->GetParError(3)));
            yPos -= 0.03;
            break;
    }

    text->DrawLatex(xStart, yPos, Form("#chi^{2}/NDF = %.3f", 
        fitFunc->GetChisquare()/fitFunc->GetNDF()));
}

// Modify FitHistogramPairAndGetEfficiency similarly to handle different fit functions
double FitHistogramPairAndGetEfficiency(TH1* histWith, TH1* histWithout,
    TF1* &fitWith, TF1* &fitWithout,
    const std::string& suffix, 
    bool draw = true, bool useTightRange = false) {

    //histWithout->Scale(1.0, "width");
    //histWith->Scale(1.0, "width");
    
    // Set fit range based on flag
    double f1 = useTightRange ? 0.485 : 0.48;
    double f2 = useTightRange ? 0.505 : 0.52;  

    // First fit the "with TOF" histogram
    if (!fitWith) {
        TString funcNameWith = Form("fitFuncWith_%s", suffix.c_str());
        
        // Number of parameters depends on the fit function type
        int nParams = 0;
        switch(gSelectedFitFunc) {
            case kCrystalBall: nParams = 7; break;
            case kGaussian: nParams = 5; break;
            case kVoigt: nParams = 6; break;
        }
        
        fitWith = new TF1(funcNameWith, CombinedFunction,  f1, f2, nParams);
        
        // Set parameter names and initial values based on fit function type
        switch(gSelectedFitFunc) {
            case kCrystalBall:
                fitWith->SetParNames("Norm", "Mean", "Sigma", "Alpha", "n", "Poly_a", "Poly_b");
                fitWith->SetParameters(histWith->GetMaximum(), 0.498, 0.005, 1.0, 5.0, 0.0, 0.0);
                
                fitWith->SetParLimits(0, 0, histWith->GetMaximum() * 2);
                fitWith->SetParLimits(1, 0.49, 0.506);
                fitWith->SetParLimits(2, 0.001, 0.01);
                //fitWith->FixParameter(2, 0.0030);// specific value
                fitWith->SetParLimits(3, 0.1, 10.0);
                fitWith->SetParLimits(4, 1.0, 10.0);
                break;
                
            case kGaussian:
                fitWith->SetParNames("Norm", "Mean", "Sigma", "Poly_a", "Poly_b");
                fitWith->SetParameters(histWith->GetMaximum(), 0.498, 0.005, 0.0, 0.0);
                
                fitWith->SetParLimits(0, 0, histWith->GetMaximum() * 2);
                fitWith->SetParLimits(1, 0.49, 0.506);
                fitWith->SetParLimits(2, 0.001, 0.01);
                //fitWith->FixParameter(2, 0.0030);// specific value
                break;
                
            case kVoigt:
                fitWith->SetParNames("Norm", "Mean", "Sigma", "Gamma", "Poly_a", "Poly_b");
                fitWith->SetParameters(histWith->GetMaximum(), 0.498, 0.005, 0.001, 0.0, 0.0);
                
                fitWith->SetParLimits(0, 0, histWith->GetMaximum() * 2);
                fitWith->SetParLimits(1, 0.49, 0.506);
                fitWith->SetParLimits(2, 0.001, 0.01);
                //fitWith->FixParameter(2, 0.0030);// specific value
                fitWith->SetParLimits(3, 0.0001, 0.005);
                break;
        }
    }
    fitWith->SetLineColor(kRed);fitWith->SetLineWidth(3);
    TFitResultPtr fitResultWith = histWith->Fit(fitWith, draw ? "RS" : "RSQ");//RSW //RSQ

    // Now fit the "without TOF" histogram, using the mean value from "with TOF" fit
    if (!fitWithout) {
        TString funcNameWithout = Form("fitFuncWithout_%s", suffix.c_str());
        
        // Number of parameters depends on the fit function type
        int nParams = 0;
        switch(gSelectedFitFunc) {
            case kCrystalBall: nParams = 7; break;
            case kGaussian: nParams = 5; break;
            case kVoigt: nParams = 6; break;
        }
        
        fitWithout = new TF1(funcNameWithout, CombinedFunction,  f1, f2, nParams);
        
        // Set parameters based on fit function type
        switch(gSelectedFitFunc) {
            case kCrystalBall:
                fitWithout->SetParNames("Norm", "Mean", "Sigma", "Alpha", "n", "Poly_a", "Poly_b");
                fitWithout->SetParameters(
                    histWithout->GetMaximum(),           // Adjust normalization for this histogram
                    fitWith->GetParameter(1),            // Use mean from "with TOF" fit
                    fitWith->GetParameter(2),            // Initial sigma from "with TOF" fit
                    fitWith->GetParameter(3),            // Initial alpha from "with TOF" fit
                    fitWith->GetParameter(4),            // Initial n from "with TOF" fit
                    0.0,                                 // Reset background parameters
                    0.0
                );
                break;
                
            case kGaussian:
                fitWithout->SetParNames("Norm", "Mean", "Sigma", "Poly_a", "Poly_b");
                fitWithout->SetParameters(
                    histWithout->GetMaximum(),           // Adjust normalization for this histogram
                    fitWith->GetParameter(1),            // Use mean from "with TOF" fit
                    fitWith->GetParameter(2),            // Initial sigma from "with TOF" fit
                    0.0,                                 // Reset background parameters
                    0.0
                );
                break;
                
            case kVoigt:
                fitWithout->SetParNames("Norm", "Mean", "Sigma", "Gamma", "Poly_a", "Poly_b");
                fitWithout->SetParameters(
                    histWithout->GetMaximum(),           // Adjust normalization for this histogram
                    fitWith->GetParameter(1),            // Use mean from "with TOF" fit
                    fitWith->GetParameter(2),            // Initial sigma from "with TOF" fit
                    fitWith->GetParameter(3),            // Initial gamma from "with TOF" fit
                    0.0,                                 // Reset background parameters
                    0.0
                );
                break;
        }

        fitWithout->SetParLimits(0, 0, histWithout->GetMaximum() * 2); //norm
        //if no fixing anything 
        //fitWithout->SetParLimits(1, 0.49, 0.506);  
        //fitWithout->SetParLimits(2, 0.001, 0.01); 
        
        fitWithout->FixParameter(1, fitWith->GetParameter(1));
        fitWithout->FixParameter(2, fitWith->GetParameter(2)); 

        //fitWithout->FixParameter(1, 0.498);// specific value
        //fitWithout->FixParameter(2, 0.0030);// specific value

        /*bool fixMean = false; // true or false, depending on requirement

        if (fixMean) {
            // Fix the mean and set sigma limits
            fitWithout->FixParameter(1, fitWith->GetParameter(1)); // Fix mean (par 1)
            fitWithout->SetParLimits(2, 0.001, 0.01);             // Set sigma (par 2) limits
        } else {
            // Fix sigma and set mean limits
            fitWithout->FixParameter(2, fitWith->GetParameter(2)); // Fix sigma (par 2)
            //fitWithout->FixParameter(2, 0.0030);// specific value
            fitWithout->SetParLimits(1, 0.49, 0.506);             // Set mean (par 1) limits
        }*/

        
        // Set other parameter limits based on fit function type
        switch(gSelectedFitFunc) {
            case kCrystalBall:                
                fitWithout->SetParLimits(3, 0.1, 10.0);
                fitWithout->SetParLimits(4, 1.0, 10.0);
                break;
                
            case kGaussian:                
                break;
                
            case kVoigt:               
                fitWithout->SetParLimits(3, 0.0001, 0.005);
                break;
        }
    }
    
    fitWithout->SetLineColor(kRed+2);fitWithout->SetLineWidth(3);
    TFitResultPtr fitResultWithout = histWithout->Fit(fitWithout, draw ? "RS" : "RSQ");

    // Create signal functions for yield calculations
    TF1* signalWith = nullptr;
    TF1* signalWithout = nullptr;
    
    switch(gSelectedFitFunc) {
        case kCrystalBall: {
            TString cbNameWith = Form("crystalBallWith_%s", suffix.c_str());
            signalWith = new TF1(cbNameWith, CrystalBall,  f1, f2, 5);
            for (int i = 0; i < 5; i++) {
                signalWith->SetParameter(i, fitWith->GetParameter(i));
            }
            
            TString cbNameWithout = Form("crystalBallWithout_%s", suffix.c_str());
            signalWithout = new TF1(cbNameWithout, CrystalBall,  f1, f2, 5);
            for (int i = 0; i < 5; i++) {
                signalWithout->SetParameter(i, fitWithout->GetParameter(i));
            }
            break;
        }
        
        case kGaussian: {
            TString gaussNameWith = Form("gaussianWith_%s", suffix.c_str());
            signalWith = new TF1(gaussNameWith, Gaussian,  f1, f2, 3);
            for (int i = 0; i < 3; i++) {
                signalWith->SetParameter(i, fitWith->GetParameter(i));
            }
            
            TString gaussNameWithout = Form("gaussianWithout_%s", suffix.c_str());
            signalWithout = new TF1(gaussNameWithout, Gaussian,  f1, f2, 3);
            for (int i = 0; i < 3; i++) {
                signalWithout->SetParameter(i, fitWithout->GetParameter(i));
            }
            break;
        }
        
        case kVoigt: {
            TString voigtNameWith = Form("voigtWith_%s", suffix.c_str());
            signalWith = new TF1(voigtNameWith, VoigtFunction,  f1, f2, 4);
            for (int i = 0; i < 4; i++) {
                signalWith->SetParameter(i, fitWith->GetParameter(i));
            }
            
            TString voigtNameWithout = Form("voigtWithout_%s", suffix.c_str());
            signalWithout = new TF1(voigtNameWithout, VoigtFunction,  f1, f2, 4);
            for (int i = 0; i < 4; i++) {
                signalWithout->SetParameter(i, fitWithout->GetParameter(i));
            }
            break;
        }
    }
    
    signalWith->SetLineColor(kBlue);
    signalWith->SetLineStyle(2);
    signalWithout->SetLineColor(kGreen+2);
    signalWithout->SetLineStyle(2);

    // Create background components
    TString polyNameWith = Form("polyWith_%s", suffix.c_str());
    TF1* polyWith = new TF1(polyNameWith, "pol1",  f1, f2);
    
    TString polyNameWithout = Form("polyWithout_%s", suffix.c_str());
    TF1* polyWithout = new TF1(polyNameWithout, "pol1",  f1, f2);
    
    switch(gSelectedFitFunc) {
        case kCrystalBall:
            polyWith->SetParameter(0, fitWith->GetParameter(5));
            polyWith->SetParameter(1, fitWith->GetParameter(6));
            polyWithout->SetParameter(0, fitWithout->GetParameter(5));
            polyWithout->SetParameter(1, fitWithout->GetParameter(6));
            break;
            
        case kGaussian:
            polyWith->SetParameter(0, fitWith->GetParameter(3));
            polyWith->SetParameter(1, fitWith->GetParameter(4));
            polyWithout->SetParameter(0, fitWithout->GetParameter(3));
            polyWithout->SetParameter(1, fitWithout->GetParameter(4));
            break;
            
        case kVoigt:
            polyWith->SetParameter(0, fitWith->GetParameter(4));
            polyWith->SetParameter(1, fitWith->GetParameter(5));
            polyWithout->SetParameter(0, fitWithout->GetParameter(4));
            polyWithout->SetParameter(1, fitWithout->GetParameter(5));
            break;
    }
    
    polyWith->SetLineColor(kBlack);
    polyWith->SetLineStyle(3);
    polyWithout->SetLineColor(kBlack);
    polyWithout->SetLineStyle(3);

    if (draw) {
        polyWith->Draw("same");
        polyWithout->Draw("same");
    }

    double yieldWith = signalWith->Integral(f1, f2);// / histWith->GetBinWidth(1);
    double yieldWithout = signalWithout->Integral(f1, f2);// / histWithout->GetBinWidth(1);

    return yieldWith / (yieldWith + yieldWithout); // Return the efficiency
}

// Modified function to create a single plot with the new fitting approach
void CreateSinglePlot(TCanvas* canvas, TH1* histWithout, TH1* histWith, const char* title,bool useTightRange = false) {
    canvas->Clear();
    canvas->cd();
    canvas->SetFrameLineWidth(3);

    SetHistogramStyle(histWithout, title, 20, kGreen+3, 2.0);    
    SetHistogramStyle(histWith, nullptr, 21, kGreen-3, 2.0);
   
    histWithout->Draw("PE");
    histWith->Draw("PE same");

    // Use the new fit function that fits "with TOF" first and then "without TOF"
    TF1 *fitWith = nullptr;
    TF1 *fitWithout = nullptr;
    double efficiency = FitHistogramPairAndGetEfficiency(
        histWith, histWithout, 
        fitWith, fitWithout,
        Form("bin_%d", canvas->GetNumber()), 
        true, useTightRange
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
        //legend->AddEntry(polyWithout, "Poly1 (without TOF)", "l");
    }

    legend->Draw();

    TLatex *effText = new TLatex();
    effText->SetNDC();
    effText->SetTextSize(0.05);
    effText->SetTextFont(42);
    effText->DrawLatex(0.38, 0.17, Form("#color[6]{Efficiency = %.1f%%}", efficiency * 100));

}

// Example of how to modify the main function to select fit type
void Plot_TofEfficiency(const char* fitType = "Gaussian") {  //Gaussian//empty for crystalball

    const char* fitfunction[] = {"Gaussian"}; //Crystalball//check above what 
    const char* fixvar[] = {"#sigma & #mu"};//#mean #sigma//check above what 
    double thresholdIndex = 0.3;

    // Set the fit function type based on input parameter
    if (strcmp(fitType, "Gaussian") == 0) {
        gSelectedFitFunc = kGaussian;
        std::cout << "Using Gaussian fit" << std::endl;
    } else if (strcmp(fitType, "Voigt") == 0) {
        gSelectedFitFunc = kVoigt;
        std::cout << "Using Voigt fit" << std::endl;
    } else {
        gSelectedFitFunc = kCrystalBall;
        std::cout << "Using Crystal Ball fit" << std::endl;
    }
    
    // Rest of your code remains the same
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    //latex
    TLatex Tl;
    Tl.SetTextAlign(10);
    Tl.SetTextSize(0.05);

    TLatex *fixedParamText = new TLatex();
    fixedParamText->SetNDC();
    fixedParamText->SetTextSize(0.035);
    fixedParamText->SetTextFont(42);
    fixedParamText->SetTextColor(kViolet+3);
    
    TFile* rootFile = TFile::Open("InputData/tofeff_withRP_0.3_v3.root");//tofeff_thesis
    if (!rootFile || rootFile->IsZombie()) {
        std::cerr << "Error opening MC input file" << std::endl;
        return;
    }
    TFile* rootFile2 = TFile::Open("InputData/tofeff_data_0.3_v3.root");//tofeff_thesis
    if (!rootFile2 || rootFile2->IsZombie()) {
        std::cerr << "Error opening data input file" << std::endl;
        return;
    }
    
    TH1* HistKaonMassProbeWithoutTof_TruePions = (TH1D*)rootFile->Get("HistKaonMassProbeWithoutTof_TruePions;1");
    TH1* HistKaonMassProbeWithTof_TruePions = (TH1D*)rootFile->Get("HistKaonMassProbeWithTof_TruePions;1");
    
    TH1* HistKaonMassProbeWithoutTof_MC = (TH1D*)rootFile->Get("HistKaonMassProbeWithoutTof;1");
    TH1* HistKaonMassProbeWithTof_MC = (TH1D*)rootFile->Get("HistKaonMassProbeWithTof;1");
         
    TH1* HistKaonMassProbeWithoutTof_data = (TH1D*)rootFile2->Get("HistKaonMassProbeWithoutTof;1");
    TH1* HistKaonMassProbeWithTof_data = (TH1D*)rootFile2->Get("HistKaonMassProbeWithTof;1");
   
    std::vector<TH1*> etaWithTof(6), etaWithoutTof(6);
    std::vector<TH1*> ptWithTof(6), ptWithoutTof(6);
    std::vector<TH1*> etaWithTof_true(6), etaWithoutTof_true(6);
    std::vector<TH1*> ptWithTof_true(6), ptWithoutTof_true(6);

    std::vector<TH1*> etaWithTof_data(6), etaWithoutTof_data(6);
    std::vector<TH1*> ptWithTof_data(6), ptWithoutTof_data(6);

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

        etaWithTof_true[i] = static_cast<TH1*>(rootFile->Get(Form("HistKaonPtTruePionWithTofMass%d;1", i)));
        etaWithoutTof_true[i] = static_cast<TH1*>(rootFile->Get(Form("HistKaonEtaTruePionWithTofMass%d;1", i)));
        ptWithTof_true[i] = static_cast<TH1*>(rootFile->Get(Form("HistKaonPtTruePionWithoutTofMass%d;1", i)));
        ptWithoutTof_true[i] = static_cast<TH1*>(rootFile->Get(Form("HistKaonEtaTruePionWithoutTofMass%d;1", i)));
        
        if (!etaWithTof_true[i] || !etaWithoutTof_true[i] || !ptWithTof_true[i] || !ptWithoutTof_true[i]) {
            std::cerr << "Missing true histogram for index " << i << std::endl;
            rootFile->Close();
            return;
        }

        etaWithTof_data[i] = static_cast<TH1*>(rootFile2->Get(Form("HistKaonEtaProbeWithTofPosProbeMass%d;1", i)));
        etaWithoutTof_data[i] = static_cast<TH1*>(rootFile2->Get(Form("HistKaonEtaProbeWithoutTofPosProbeMass%d;1", i)));
        ptWithTof_data[i] = static_cast<TH1*>(rootFile2->Get(Form("HistKaonPtProbeWithTofPosProbeMass%d;1", i)));
        ptWithoutTof_data[i] = static_cast<TH1*>(rootFile2->Get(Form("HistKaonPtProbeWithoutTofPosProbeMass%d;1", i)));

        if (!etaWithTof_data[i] || !etaWithoutTof_data[i] || !ptWithTof_data[i] || !ptWithoutTof_data[i]) {
            std::cerr << "Missing data histogram for index " << i << std::endl;
            rootFile2->Close();
            return;
        }

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
                          "Tag and Probe without true pion matching"};

    TCanvas* canvas = new TCanvas("canvas", "TOF Efficiency", 1500, 1048);//1000, 700 //2860, 2000
    canvas->SetLeftMargin(0.13);
    canvas->SetRightMargin(0.19);
    canvas->SetTopMargin(0.10);
    canvas->SetBottomMargin(0.15);

    // Eta plots
    canvas->Print("plots/TofEff_eta_thesis_MC_GS_0.3_v3.pdf[");
    for (int i = 0; i < 6; i++) {
        CreateSinglePlot(canvas, etaWithoutTof[i], etaWithTof[i], etaTitles[i],  false);
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
        CreateSinglePlot(canvas, ptWithoutTof[i], ptWithTof[i], ptTitles[i],  false);
        fixedParamText->DrawLatex(0.18, 0.7, Form("#splitline{%s fixed}{Fit with: %s distribution}",fixvar[0],fitfunction[0]));
        int hist_max = ptWithoutTof[i]->GetMaximumBin();
        Double_t x_max = ptWithoutTof[i]->GetBinCenter(hist_max);
        Double_t y_max = ptWithoutTof[i]->GetBinContent(hist_max);   
        Tl.DrawLatex(0.51, y_max + 10, Form("#scale[0.6]{matching threshold < %.2f}", thresholdIndex));
        canvas->Print("plots/TofEff_pt_thesis_MC_GS_0.3_v3.pdf");
    }
    canvas->Print("plots/TofEff_pt_thesis_MC_GS_0.3_v3.pdf]");

    // tag and probe plots  
    CreateSinglePlot(canvas, HistKaonMassProbeWithoutTof_MC, HistKaonMassProbeWithTof_MC, some[1],  false);
    HistKaonMassProbeWithTof_MC->GetYaxis()->SetRangeUser(-20, 150); 
    int hist_max2 = HistKaonMassProbeWithTof_MC->GetMaximumBin();
    Double_t x_max2 = HistKaonMassProbeWithTof_MC->GetBinCenter(hist_max2);
    Double_t y_max2 = HistKaonMassProbeWithTof_MC->GetBinContent(hist_max2); 
    Tl.DrawLatex(0.51, y_max2 + 50, Form("#scale[0.6]{matching threshold < %.2f}", thresholdIndex));   
    fixedParamText->DrawLatex(0.18, 0.7, Form("#splitline{%s fixed}{Fit with: %s distribution}",fixvar[0],fitfunction[0]));
    canvas->Print("plots/TofEff_tagNprobe_thesis_MC_GS_0.3_v3.png");

    
    // true plots  
    CreateSinglePlot(canvas, HistKaonMassProbeWithoutTof_TruePions, HistKaonMassProbeWithTof_TruePions, some[0], true);
    HistKaonMassProbeWithoutTof_TruePions->GetYaxis()->SetRangeUser(-20, 150);
    int hist_max = HistKaonMassProbeWithoutTof_TruePions->GetMaximumBin();
    Double_t x_max = HistKaonMassProbeWithoutTof_TruePions->GetBinCenter(hist_max);
    Double_t y_max = HistKaonMassProbeWithoutTof_TruePions->GetBinContent(hist_max);   
    Tl.DrawLatex(0.51, y_max + 10, Form("#scale[0.6]{matching threshold < %.2f}", thresholdIndex));
    fixedParamText->DrawLatex(0.18, 0.7, Form("#splitline{%s fixed}{Fit with: %s distribution}",fixvar[0],fitfunction[0]));
    canvas->Print("plots/TofEff_true_thesis_MC_GS_0.3_v3.png"); 
    
    // True Eta plots
    canvas->Print("plots/TofEff_true_eta_thesis_MC_GS_0.3_v3.pdf[");
    for (int i = 0; i < 6; i++) {
        CreateSinglePlot(canvas, etaWithoutTof_true[i], etaWithTof_true[i], etaTitles[i], true);
        fixedParamText->DrawLatex(0.18, 0.7, Form("#splitline{%s fixed}{Fit with: %s distribution}",fixvar[0],fitfunction[0]));
        int hist_max = etaWithoutTof_true[i]->GetMaximumBin();
        Double_t x_max = etaWithoutTof_true[i]->GetBinCenter(hist_max);
        Double_t y_max = etaWithoutTof_true[i]->GetBinContent(hist_max);   
        Tl.DrawLatex(0.51, y_max + 10, Form("#scale[0.6]{matching threshold < %.2f}", thresholdIndex));
        canvas->Print("plots/TofEff_true_eta_thesis_MC_GS_0.3_v3.pdf");
    }
    canvas->Print("plots/TofEff_true_eta_thesis_MC_GS_0.3_v3.pdf]");

    // True pT plots
    canvas->Print("plots/TofEff_true_pt_thesis_MC_GS_0.3_v3.pdf[");
    for (int i = 0; i < 6; i++) {
        CreateSinglePlot(canvas, ptWithoutTof_true[i], ptWithTof_true[i], ptTitles[i], true);
        fixedParamText->DrawLatex(0.18, 0.7, Form("#splitline{%s fixed}{Fit with: %s distribution}",fixvar[0],fitfunction[0]));
        int hist_max = ptWithoutTof_true[i]->GetMaximumBin();
        Double_t x_max = ptWithoutTof_true[i]->GetBinCenter(hist_max);
        Double_t y_max = ptWithoutTof_true[i]->GetBinContent(hist_max);   
        Tl.DrawLatex(0.51, y_max + 10, Form("#scale[0.6]{matching threshold < %.2f}", thresholdIndex));
        canvas->Print("plots/TofEff_true_pt_thesis_MC_GS_0.3_v3.pdf");
    }
    canvas->Print("plots/TofEff_true_pt_thesis_MC_GS_0.3_v3.pdf]");

    //Data
    // Eta plots
    canvas->Print("plots/TofEff_eta_thesis_data_GS_0.3_v3.pdf[");
    for (int i = 0; i < 6; i++) {
        CreateSinglePlot(canvas, etaWithoutTof_data[i], etaWithTof_data[i], etaTitles[i],  false);
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
        CreateSinglePlot(canvas, ptWithoutTof_data[i], ptWithTof_data[i], ptTitles[i], false);
        fixedParamText->DrawLatex(0.18, 0.7, Form("#splitline{%s fixed}{Fit with: %s distribution}",fixvar[0],fitfunction[0]));
        int hist_max = ptWithoutTof_data[i]->GetMaximumBin();
        Double_t x_max = ptWithoutTof_data[i]->GetBinCenter(hist_max);
        Double_t y_max = ptWithoutTof_data[i]->GetBinContent(hist_max);   
        Tl.DrawLatex(0.51, y_max + 10, Form("#scale[0.6]{matching threshold < %.2f}", thresholdIndex));
        canvas->Print("plots/TofEff_pt_thesis_data_GS_0.3_v3.pdf");
    }
    canvas->Print("plots/TofEff_pt_thesis_data_GS_0.3_v3.pdf]");
    
    // Global 
    CreateSinglePlot(canvas, HistKaonMassProbeWithoutTof_data, HistKaonMassProbeWithTof_data, some[1], false);
    HistKaonMassProbeWithTof_data->GetYaxis()->SetRangeUser(-20, 150); 
    int hist_max3 = HistKaonMassProbeWithTof_data->GetMaximumBin();
    Double_t x_max3 = HistKaonMassProbeWithTof_data->GetBinCenter(hist_max3);
    Double_t y_max3 = HistKaonMassProbeWithTof_data->GetBinContent(hist_max3); 
    Tl.DrawLatex(0.51, y_max2 + 50, Form("#scale[0.6]{matching threshold < %.2f}", thresholdIndex));   
    fixedParamText->DrawLatex(0.18, 0.7, Form("#splitline{%s fixed}{Fit with: %s distribution}",fixvar[0],fitfunction[0]));
    canvas->Print("plots/TofEff_tagNprobe_thesis_data_GS_0.3_v3.png");
    

    double etaBins[7] = {-0.9, -0.6, -0.3, 0.0, 0.3, 0.6, 0.9}; // Use your actual bin edges
    double ptBins[7] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.2};    // Use your actual bin edges

    const int nPoints = 6;
    double etaValues[nPoints] = {-0.75, -0.45, -0.15, 0.15, 0.45, 0.75};
    double ptValues[nPoints] = {0.25, 0.35, 0.45, 0.55, 0.7, 1.0};
    double effEtaValues[nPoints] = {0}; //  fill these in the loop
    double effPtValues[nPoints] = {0};  //  fill these in the loop
    double effEtaErrors[nPoints] = {0}; //  fill these in the loop
    double effPtErrors[nPoints] = {0};  //  fill these in the loop

    // Define horizontal errors (half bin width)
    double etaXErrors[nPoints] = {0.15, 0.15, 0.15, 0.15, 0.15, 0.15}; // Based on bin width of 0.3
    double ptXErrors[nPoints] = {0.05, 0.05, 0.05, 0.05, 0.1, 0.2};    // Based on variable bin widths
    
    // Summary plots
    TH1D* effVsEta = new TH1D("effVsEta", "TOF Efficiency vs #eta;#eta;Efficiency", 6, etaBins);
    TH1D* effVsPt = new TH1D("effVsPt", "TOF Efficiency vs p_{T};p_{T} [GeV/c];Efficiency", 6, ptBins);
    
    TH1D* effVsEta_true = new TH1D("effVsEta_true", "TOF Efficiency vs #eta;#eta;Efficiency", 6, etaBins);
    TH1D* effVsPt_true = new TH1D("effVsPt_true", "TOF Efficiency vs p_{T};p_{T} [GeV/c];Efficiency", 6, ptBins);

    TH1D* effVsEta_data = new TH1D("effVsEta_data", "TOF Efficiency vs #eta;#eta;Efficiency", 6, etaBins);
    TH1D* effVsPt_data = new TH1D("effVsPt_data", "TOF Efficiency vs p_{T};p_{T} [GeV/c];Efficiency", 6, ptBins);

    for (int i = 0; i < 6; i++) {

        TF1 *fitEtaWith = nullptr, *fitEtaWithout = nullptr;
        TF1 *fitPtWith = nullptr, *fitPtWithout = nullptr;

        // Use the paired function for summary calculation
        double effEta = FitHistogramPairAndGetEfficiency(
            etaWithTof[i], etaWithoutTof[i], 
            fitEtaWith, fitEtaWithout,
            Form("eta_%d", i), false
        );

        double effPt = FitHistogramPairAndGetEfficiency(
            ptWithTof[i], ptWithoutTof[i], 
            fitPtWith, fitPtWithout,
            Form("pt_%d", i), false
        );
        double yieldEtaWith = FitHistogramAndGetYield(etaWithTof[i], fitEtaWith, Form("eta_w_%d", i), false);
        double yieldEtaWithout = FitHistogramAndGetYield(etaWithoutTof[i], fitEtaWithout, Form("eta_wo_%d", i), false);
        double yieldPtWith = FitHistogramAndGetYield(ptWithTof[i], fitPtWith, Form("pt_w_%d", i), false);
        double yieldPtWithout = FitHistogramAndGetYield(ptWithoutTof[i], fitPtWithout, Form("pt_wo_%d", i), false);

        effEtaValues[i] = effEta;
        effPtValues[i] = effPt;
        effEtaErrors[i] = sqrt(effEta*(1-effEta)/(yieldEtaWith + yieldEtaWithout));
        effPtErrors[i] = sqrt(effPt*(1-effPt)/(yieldPtWith + yieldPtWithout));

        
        effVsEta->SetBinContent(i+1, effEta);
        effVsEta->SetBinError(i+1, sqrt(effEta*(1-effEta)/(yieldEtaWith + yieldEtaWithout)));
        effVsPt->SetBinContent(i+1, effPt);
        effVsPt->SetBinError(i+1, sqrt(effPt*(1-effPt)/(yieldPtWith + yieldPtWithout)));
        

        TF1 *fitEtaWith_true = nullptr, *fitEtaWithout_true = nullptr;
        TF1 *fitPtWith_true = nullptr, *fitPtWithout_true = nullptr;

        // Use the paired function for summary calculation
        double effEta_true= FitHistogramPairAndGetEfficiency(
            etaWithTof_true[i], etaWithoutTof_true[i], 
            fitEtaWith_true, fitEtaWithout_true,
            Form("eta_true_%d", i), false
        );

        double effPt_true = FitHistogramPairAndGetEfficiency(
            ptWithTof_true[i], ptWithoutTof_true[i], 
            fitPtWith_true, fitPtWithout_true,
            Form("pt_true_%d", i), false
        );

        
        double yieldEtaWith_true = FitHistogramAndGetYield(etaWithTof_true[i],       fitEtaWith_true,    Form("eta_w_true_%d",  i), false);
        double yieldEtaWithout_true = FitHistogramAndGetYield(etaWithoutTof_true[i], fitEtaWithout_true, Form("eta_wo_true_%d", i), false);
        double yieldPtWith_true = FitHistogramAndGetYield(ptWithTof_true[i],         fitPtWith_true,     Form("pt_w_true_%d",   i), false);
        double yieldPtWithout_true = FitHistogramAndGetYield(ptWithoutTof_true[i],   fitPtWithout_true,  Form("pt_wo_true_%d",  i), false);

        effVsEta_true->SetBinContent(i+1, effEta_true);
        effVsEta_true->SetBinError(i+1, sqrt(effEta_true*(1-effEta_true)/(yieldEtaWith_true + yieldEtaWithout_true)));
        effVsPt_true->SetBinContent(i+1, effPt_true);
        effVsPt_true->SetBinError(i+1, sqrt(effPt_true*(1-effPt_true)/(yieldPtWith_true + yieldPtWithout_true)));
        
        TF1 *fitEtaWith_data = nullptr, *fitEtaWithout_data = nullptr;
        TF1 *fitPtWith_data = nullptr, *fitPtWithout_data = nullptr;

        // Use the paired function for summary calculation
        double effEta_data= FitHistogramPairAndGetEfficiency(
            etaWithTof_data[i], etaWithoutTof_data[i], 
            fitEtaWith_data, fitEtaWithout_data,
            Form("eta_data_%d", i), false
        );

        double effPt_data = FitHistogramPairAndGetEfficiency(
            ptWithTof_data[i], ptWithoutTof_data[i], 
            fitPtWith_data, fitPtWithout_data,
            Form("pt_data_%d", i), false
        );
        double yieldEtaWith_data = FitHistogramAndGetYield(etaWithTof_data[i],       fitEtaWith_data,    Form("eta_w_data_%d",  i), false);
        double yieldEtaWithout_data = FitHistogramAndGetYield(etaWithoutTof_data[i], fitEtaWithout_data, Form("eta_wo_data_%d", i), false);
        double yieldPtWith_data = FitHistogramAndGetYield(ptWithTof_data[i],         fitPtWith_data,     Form("pt_w_data_%d",   i), false);
        double yieldPtWithout_data = FitHistogramAndGetYield(ptWithoutTof_data[i],   fitPtWithout_data,  Form("pt_wo_data_%d",  i), false);

        effVsEta_data->SetBinContent(i+1, effEta_data);
        effVsEta_data->SetBinError(i+1, sqrt(effEta_data*(1-effEta_data)/(yieldEtaWith_data + yieldEtaWithout_data)));
        effVsPt_data->SetBinContent(i+1, effPt_data);
        effVsPt_data->SetBinError(i+1, sqrt(effPt_data*(1-effPt_data)/(yieldPtWith_data + yieldPtWithout_data)));

        // Clean up
        delete fitEtaWith;
        delete fitEtaWithout;
        delete fitPtWith;
        delete fitPtWithout;
        
        delete fitEtaWith_true;
        delete fitEtaWithout_true;
        delete fitPtWith_true;
        delete fitPtWithout_true;

        delete fitEtaWith_data;
        delete fitEtaWithout_data;
        delete fitPtWith_data;
        delete fitPtWithout_data;
    }
        
    canvas->Clear();
    canvas->SetLeftMargin(0.15);
    canvas->SetRightMargin(0.10);
    canvas->SetTopMargin(0.10);
    canvas->SetBottomMargin(0.15);
    canvas->SetFrameLineWidth(6);
    
    SetHistogramStyle2(effVsEta, kBlue);
    SetHistogramStyle2(effVsEta_true, kRed);
    SetHistogramStyle2(effVsEta_data, kGreen);

    TLegend* legend = new TLegend(0.72, 0.25, 0.9, 0.45);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(effVsEta, "MC", "p");
    legend->AddEntry(effVsEta_true, "True", "p");
    legend->AddEntry(effVsEta_data, "Data", "p");
    effVsEta->Draw("PE");
    effVsEta_true->Draw("PE same");
    effVsEta_data->Draw("PE same");
    legend->Draw();
    //TGraphErrors* graphEta = new TGraphErrors(nPoints, etaValues, effEtaValues, etaXErrors, effEtaErrors);
    //graphEta->SetTitle("TOF Efficiency vs #eta");
    //graphEta->GetXaxis()->SetTitle("#eta");
    //graphEta->GetYaxis()->SetTitle("Efficiency");
    //graphEta->Draw("APE");  legend->Draw();
    Tl.DrawLatex(-0.8,0.1, Form("#splitline{%s fixed}{Fit with: %s distribution}",fixvar[0],fitfunction[0]));
    canvas->Print("plots/TofEff_summary_eta_thesis_MC_GS_0.3_v3.png");

    canvas->Clear();
    SetHistogramStyle2(effVsPt, kBlue);
    SetHistogramStyle2(effVsPt_true, kRed);
    SetHistogramStyle2(effVsPt_data, kGreen);
    effVsPt->Draw("PE");
    effVsPt_true->Draw("PE same");
    effVsPt_data->Draw("PE same");
    legend->Draw();
    //TGraphErrors* graphPt = new TGraphErrors(nPoints, ptValues, effPtValues, ptXErrors, effPtErrors);
    //graphPt->SetTitle("TOF Efficiency vs p_{T}");
    //graphPt->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    //graphPt->GetYaxis()->SetTitle("Efficiency");
    //graphPt->Draw("APE"); 
    Tl.DrawLatex(0.25,0.1, Form("#splitline{%s fixed}{Fit with: %s distribution}",fixvar[0],fitfunction[0]));
    canvas->Print("plots/TofEff_summary_pt_thesis_MC_GS_0.3_v3.png");

    delete canvas;
    delete effVsEta;
    delete effVsPt;
    //delete graphEta;
    //delete graphPt;
    rootFile->Close();
}