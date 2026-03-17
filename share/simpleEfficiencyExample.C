//==============================================================================
// EXAMPLE: How to Use StEfficiencyCorrector3D Class
// Author: Sneha Bhosale (AGH Krakow)
// Purpose: Learn efficiency correction with proper 3D treatment 
// The example is correcting to itself (closure test)
// shown for one particle type (positive pions) and Combined TPC+TOF efficiency
//==============================================================================

#include "StEfficiencyCorrector3D.h" // Efficiency corrector class
#include <TFile.h>
#include <TH3F.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <iostream>

void simpleEfficiencyExample() {
    
    // SILENCE INTERPOLATION WARNINGS (they're non-fatal)
    gErrorIgnoreLevel = kError + 1;  // Suppress all Error messages, show Fatal only
    
    std::cout << "=== Simple Efficiency Correction Tutorial ===" << std::endl;
    
    //==========================================================================
    // STEP 1: Open your efficiency file
    //==========================================================================
    TFile* effFile = TFile::Open("../Downloads/PionsAug12_0.root", "READ");
    if (!effFile || effFile->IsZombie()) {
        std::cout << "ERROR: Cannot open efficiency file!" << std::endl;
        return;
    }
    std::cout << "✓ Step 1: Opened efficiency file" << std::endl;
    
    //==========================================================================
    // STEP 2: Load numerator and denominator histograms
    // Format: h3D_[detector]_[num/den]_[charge]
    // Axes: X=pT, Y=eta, Z=Vz
    //==========================================================================
    TH3F* h3D_TPC_num_plus = (TH3F*)effFile->Get("h3D_TPC_num_P");
    TH3F* h3D_TPC_den_plus = (TH3F*)effFile->Get("h3D_TPC_den_P");
    TH3F* h3D_TOF_num_plus = (TH3F*)effFile->Get("h3D_TOF_num_P");
    TH3F* h3D_TOF_den_plus = (TH3F*)effFile->Get("h3D_TOF_den_P");
    
    if (!h3D_TPC_num_plus || !h3D_TPC_den_plus) {
        std::cout << "ERROR: Cannot load histograms!" << std::endl;
        return;
    }
    std::cout << "✓ Step 2: Loaded TPC and TOF histograms" << std::endl;
    
    //==========================================================================
    // STEP 3: Calculate efficiency = numerator / denominator
    //==========================================================================
    TH3F* tpcEfficiency = (TH3F*)h3D_TPC_num_plus->Clone("tpcEfficiency");
    TH3F* tofEfficiency = (TH3F*)h3D_TOF_num_plus->Clone("tofEfficiency");
    
    // Binomial division for proper error handling
    tpcEfficiency->Divide(h3D_TPC_num_plus, h3D_TPC_den_plus, 1, 1, "B");
    tofEfficiency->Divide(h3D_TOF_num_plus, h3D_TOF_den_plus, 1, 1, "B");
    
    std::cout << "✓ Step 3: Calculated efficiencies" << std::endl;
    
    //==========================================================================
    // STEP 4: Create and initialize the efficiency corrector
    //==========================================================================
    StEfficiencyCorrector3D corrector;
    
    int charge = +1;        // Positive charge
    int particlePid = 0;    // Pion=0, Kaon=1, Proton=2
    bool cloneHist = true;  // Make internal copy
    
    // Set TPC efficiency
    corrector.setTpcEfficiency(tpcEfficiency, charge, particlePid, cloneHist);
    
    // Set TOF efficiency
    corrector.setTofEfficiency(tofEfficiency, charge, particlePid, cloneHist);
    
    std::cout << "✓ Step 4: Initialized corrector with efficiencies" << std::endl;
    
    //==========================================================================
    // STEP 5: Apply efficiency corrections in 3D (CORRECT METHOD)
    //==========================================================================
    std::cout << "✓ Step 5: Applying 3D efficiency corrections..." << std::endl;
    
    // Create corrected histogram by cloning the reconstructed one
    TH3F* h3D_corrected = (TH3F*)h3D_TPC_num_plus->Clone("h3D_corrected");
    h3D_corrected->Reset(); // Clear contents but keep structure
    
    // Statistics counters
    int totalBins = 0;
    int correctedBins = 0;
    int skippedBins = 0;
    
    // Loop over all bins in 3D space
    for (int binX = 1; binX <= h3D_TPC_num_plus->GetNbinsX(); binX++) {
        for (int binY = 1; binY <= h3D_TPC_num_plus->GetNbinsY(); binY++) {
            for (int binZ = 1; binZ <= h3D_TPC_num_plus->GetNbinsZ(); binZ++) {
                
                totalBins++;
                
                // Get bin center coordinates
                Double_t pT = h3D_TPC_num_plus->GetXaxis()->GetBinCenter(binX);
                Double_t eta = h3D_TPC_num_plus->GetYaxis()->GetBinCenter(binY);
                Double_t vz = h3D_TPC_num_plus->GetZaxis()->GetBinCenter(binZ);
                
                // Get original content and error
                Double_t content = h3D_TPC_num_plus->GetBinContent(binX, binY, binZ);
                Double_t error = h3D_TPC_num_plus->GetBinError(binX, binY, binZ);
                
                // Skip empty bins
                if (content < 1e-10) {
                    h3D_corrected->SetBinContent(binX, binY, binZ, 0);
                    h3D_corrected->SetBinError(binX, binY, binZ, 0);
                    continue;
                }
                
                // Get efficiency at this specific point in phase space
                Double_t eff = corrector.getCombinedEfficiency(eta, pT, vz, charge, particlePid);

                //Optinal: // Check if within analysis ranges
                
                // Apply correction factor = 1 / efficiency
                Double_t corrFactor = 1.0;
                if (eff > 1e-4) {  // Avoid division by zero
                    corrFactor = 1.0 / eff;
                    correctedBins++;
                } else {
                    skippedBins++;
                }
                
                // Set corrected values
                h3D_corrected->SetBinContent(binX, binY, binZ, content * corrFactor);
                h3D_corrected->SetBinError(binX, binY, binZ, error * corrFactor);
            }
        }
    }
    
    std::cout << "   Total bins processed: " << totalBins << std::endl;
    std::cout << "   Bins corrected: " << correctedBins << std::endl;
    std::cout << "   Bins skipped (low efficiency): " << skippedBins << std::endl;
    
    //==========================================================================
    // STEP 6: Project to 1D for visualization
    //==========================================================================
    Double_t refEta = 0.0;
    Double_t refVz = 0.0;
    
    // Define projection ranges
    Int_t etaBinMin = TMath::Max(1, h3D_TPC_den_plus->GetYaxis()->FindBin(refEta - 0.1));
    Int_t etaBinMax = TMath::Min(h3D_TPC_den_plus->GetNbinsY(), 
                                  h3D_TPC_den_plus->GetYaxis()->FindBin(refEta + 0.1));
    Int_t vzBinMin = TMath::Max(1, h3D_TPC_den_plus->GetZaxis()->FindBin(refVz - 5.0));
    Int_t vzBinMax = TMath::Min(h3D_TPC_den_plus->GetNbinsZ(), 
                                 h3D_TPC_den_plus->GetZaxis()->FindBin(refVz + 5.0));
    
    // Project True (denominator)
    TH1D* hTrue = h3D_TPC_den_plus->ProjectionX("hTrue", etaBinMin, etaBinMax, 
                                                 vzBinMin, vzBinMax);
    hTrue->SetTitle("True MC;p_{T} (GeV/c);Counts");
    
    // Project Reconstructed (numerator)
    TH1D* hReco = h3D_TPC_num_plus->ProjectionX("hReco", etaBinMin, etaBinMax, 
                                                 vzBinMin, vzBinMax);
    hReco->SetTitle("Reconstructed;p_{T} (GeV/c);Counts");
    
    // Project Corrected
    TH1D* hCorrected = h3D_corrected->ProjectionX("hCorrected", etaBinMin, etaBinMax, 
                                                   vzBinMin, vzBinMax);
    hCorrected->SetTitle("Corrected Spectrum;p_{T} (GeV/c);Counts");
    
    std::cout << "✓ Step 6: Created 1D projections" << std::endl;
    
    // Print integral statistics
    std::cout << "\n=== Integral Check ===" << std::endl;
    std::cout << "Total True: " << hTrue->Integral() << std::endl;
    std::cout << "Total Reco: " << hReco->Integral() << std::endl;
    std::cout << "Total Corrected: " << hCorrected->Integral() << std::endl;
    std::cout << "Global efficiency: " << hReco->Integral() / hTrue->Integral() << std::endl;
    std::cout << "Closure (Corrected/True): " << hCorrected->Integral() / hTrue->Integral() << std::endl;
    
    //==========================================================================
    // STEP 7: Visualize the results
    //==========================================================================
    TCanvas* canvas = new TCanvas("canvas", "Efficiency Correction Example", 1400, 1000);
    canvas->Divide(2, 2);
    
    // Plot 1: True vs Reco vs Corrected (Normalized)
    canvas->cd(1);
    gPad->SetLogy(0);
    
    hTrue->SetLineColor(kBlack);
    hTrue->SetLineWidth(2);
    hTrue->SetTitle("Comparison (Normalized);p_{T} (GeV/c);Normalized Counts");
    hTrue->DrawNormalized("HIST");
    
    hReco->SetLineColor(kRed);
    hReco->SetLineWidth(2);
    hReco->SetLineStyle(2);
    hReco->DrawNormalized("HIST SAME");
    
    hCorrected->SetLineColor(kBlue);
    hCorrected->SetLineWidth(2);
    hCorrected->DrawNormalized("HIST SAME");
    
    TLegend* leg1 = new TLegend(0.55, 0.65, 0.88, 0.88);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->AddEntry(hTrue, "MC Truth", "l");
    leg1->AddEntry(hReco, "Reconstructed", "l");
    leg1->AddEntry(hCorrected, "Corrected", "l");
    leg1->Draw();
    
    // Plot 2: Efficiency vs pT
    canvas->cd(2);
    TH1F* hEffVsPt = new TH1F("hEffVsPt", "Efficiency vs p_{T};p_{T} (GeV/c);Efficiency", 
                               50, 0.25, 2.0);
    for (int i = 1; i <= hEffVsPt->GetNbinsX(); i++) {
        Double_t pT = hEffVsPt->GetXaxis()->GetBinCenter(i);
        Double_t eff = corrector.getCombinedEfficiency(refEta, pT, refVz, charge, particlePid);
        hEffVsPt->SetBinContent(i, eff);
    }
    hEffVsPt->SetLineColor(kBlack);
    hEffVsPt->SetLineWidth(2);
    hEffVsPt->SetMinimum(0.0);
    hEffVsPt->SetMaximum(1.0);
    hEffVsPt->Draw("HIST");
    
    // Plot 3: Correction Factor vs pT
    canvas->cd(3);
    TH1F* hCorrFactorVsPt = new TH1F("hCorrFactorVsPt", 
                                      "Correction Factor vs p_{T};p_{T} (GeV/c);Correction Factor", 
                                      50, 0.25, 2.0);
    for (int i = 1; i <= hCorrFactorVsPt->GetNbinsX(); i++) {
        Double_t pT = hCorrFactorVsPt->GetXaxis()->GetBinCenter(i);
        Double_t eff = corrector.getCombinedEfficiency(refEta, pT, refVz, charge, particlePid);
        Double_t corrFactor = (eff > 1e-4) ? 1.0/eff : 1.0;
        hCorrFactorVsPt->SetBinContent(i, corrFactor);
    }
    hCorrFactorVsPt->SetLineColor(kGreen+2);
    hCorrFactorVsPt->SetLineWidth(2);
    hCorrFactorVsPt->SetMinimum(1.0);
    hCorrFactorVsPt->SetMaximum(5.0);
    hCorrFactorVsPt->Draw("HIST");
    
    // Plot 4: Closure Test (Corrected / True ratio)
    canvas->cd(4);
    //TH1D* hClosure = (TH1D*)hCorrected->Clone("hClosure");
    //hClosure->Divide(hTrue);
    // Slghly different approach - only divide bins with sufficient statistics
    TH1D* hClosure = (TH1D*)hCorrected->Clone("hClosure");
    for (int bin = 1; bin <= hClosure->GetNbinsX(); bin++) {
        Double_t trueVal = hTrue->GetBinContent(bin);
        Double_t corrVal = hCorrected->GetBinContent(bin);
        
        if (trueVal > 10.0) {  // Require at least 10 counts in truth
            Double_t ratio = corrVal / trueVal;
            Double_t error = ratio * TMath::Sqrt(
                TMath::Power(hCorrected->GetBinError(bin)/corrVal, 2) +
                TMath::Power(hTrue->GetBinError(bin)/trueVal, 2)
            );
            hClosure->SetBinContent(bin, ratio);
            hClosure->SetBinError(bin, error);
        } else {
            hClosure->SetBinContent(bin, 0);
            hClosure->SetBinError(bin, 0);
        }
    }

    hClosure->SetTitle("Closure Test: Corrected/True;p_{T} (GeV/c);Ratio");
    hClosure->SetLineColor(kBlue);
    hClosure->SetLineWidth(2);
    hClosure->SetMarkerStyle(20);
    hClosure->SetMarkerSize(0.8);
    hClosure->SetMarkerColor(kBlue);
    hClosure->SetMinimum(0.5);
    hClosure->SetMaximum(1.5);
    hClosure->Draw("E");
    
    // Add reference lines
    TLine* unity = new TLine(hClosure->GetXaxis()->GetXmin(), 1.0,
                             hClosure->GetXaxis()->GetXmax(), 1.0);
    unity->SetLineStyle(2);
    unity->SetLineColor(kRed);
    unity->SetLineWidth(2);
    unity->Draw();
    
    TLine* upper = new TLine(hClosure->GetXaxis()->GetXmin(), 1.1,
                             hClosure->GetXaxis()->GetXmax(), 1.1);
    upper->SetLineStyle(3);
    upper->SetLineColor(kGray+1);
    upper->Draw();
    
    TLine* lower = new TLine(hClosure->GetXaxis()->GetXmin(), 0.9,
                             hClosure->GetXaxis()->GetXmax(), 0.9);
    lower->SetLineStyle(3);
    lower->SetLineColor(kGray+1);
    lower->Draw();
    
    // Add statistics box
    TLegend* leg4 = new TLegend(0.15, 0.75, 0.45, 0.88);
    leg4->SetBorderSize(1);
    leg4->SetFillStyle(1001);
    leg4->SetFillColor(kWhite);
    leg4->AddEntry((TObject*)0, Form("Mean: %.3f", hClosure->GetMean(2)), "");
    leg4->AddEntry((TObject*)0, Form("RMS: %.3f", hClosure->GetRMS(2)), "");
    leg4->Draw();
    
    canvas->SaveAs("simple_efficiency_correction_example.png");
    std::cout << "✓ Step 7: Created visualization plots" << std::endl;
    
    //==========================================================================
    // KEY CONCEPTS SUMMARY
    //==========================================================================
    std::cout << "\n=== Key Concepts ===" << std::endl;
    std::cout << "1. Efficiency = Reconstructed / True (in 3D: eta, pT, Vz)" << std::endl;
    std::cout << "2. Correction Factor = 1 / Efficiency" << std::endl;
    std::cout << "3. Apply correction to EACH bin in 3D BEFORE projection" << std::endl;
    std::cout << "4. Good closure: Corrected ≈ True (ratio ≈ 1.0 ± 0.1)" << std::endl;
    
    std::cout << "\n=== Important Methods ===" << std::endl;
    std::cout << "- setTpcEfficiency(): Load TPC efficiency histogram" << std::endl;
    std::cout << "- setTofEfficiency(): Load TOF efficiency histogram" << std::endl;
    std::cout << "- getCombinedEfficiency(eta, pT, vz, charge, pid): Get efficiency" << std::endl;
    
    std::cout << "\n=== Common Mistakes to Avoid ===" << std::endl;
    std::cout << "✗ DON'T project first, then correct" << std::endl;
    std::cout << "✓ DO correct in 3D, then project" << std::endl;
    std::cout << "✗ DON'T use fixed efficiency for all bins" << std::endl;
    std::cout << "✓ DO use bin-specific efficiency" << std::endl;
    
    // Cleanup
    effFile->Close();
    delete effFile;
    
    std::cout << "\n=== Analysis Complete ===" << std::endl;

}

/*
// STEP 1: Calculate efficiency from MC simulation
TFile* mcFile = TFile::Open("MC_simulation.root");
TH3F* mc_true = (TH3F*)mcFile->Get("h3D_true");
TH3F* mc_reco = (TH3F*)mcFile->Get("h3D_reco");
TH3F* efficiency = mc_reco / mc_true;

corrector.setTpcEfficiency(efficiency, charge, pid, true);

// STEP 2: Apply to REAL EXPERIMENTAL DATA (different file!)
TFile* dataFile = TFile::Open("RealCollisionData_2024.root");
TH3F* real_data_raw = (TH3F*)dataFile->Get("h3D_pion_raw");

// Correct the real data
TH3F* real_data_corrected = (TH3F*)real_data_raw->Clone();
for (all bins) {
    efficiency = corrector.getCombinedEfficiency(eta, pT, vz, charge, pid);
    real_data_corrected->SetBinContent(bin, raw_content / efficiency);
}
*/
