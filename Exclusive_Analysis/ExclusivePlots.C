#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TString.h>
#include <TF1.h>
#include <TLatex.h>
#include <TBox.h>
#include <iostream>

using namespace std;

// Standard colors matching the presentation (5 TOF=Blue, 4 TOF=Black, 3 TOF=Red, 2 TOF=Green)
const int colors[4] = {kGreen+1, kRed, kBlack, kBlue};
const TString tofLabels[4] = {"2 TOF", "3 TOF", "4 TOF", "5 TOF"};

void SetMyStyle() {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    
    // Increased margins to make room for larger text
    gStyle->SetPadLeftMargin(0.14);
    gStyle->SetPadBottomMargin(0.14);
    
    // Increased label and title sizes
    gStyle->SetTitleSize(0.06, "XYZ");
    gStyle->SetLabelSize(0.05, "XYZ");
}

// Helper function to plot overlaid MC with Data for N=2,3,4,5
void PlotStackedMultiplicity(TFile* fData, TFile* fMC, TString baseName, TString xTitle, 
                             double xMin, double xMax, TString outName, bool rebin = false) 
{
    TCanvas* c = new TCanvas("c_"+outName, "c", 800, 600);
    THStack* hsMC = new THStack("hsMC_"+outName, "");
    
    // Adjusted legend to not clash with larger axes
    TLegend* leg = new TLegend(0.52, 0.65, 0.88, 0.88);
    leg->SetNColumns(2);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->SetHeader("#font[22]{STAR pp #sqrt{s} = 510 GeV}");

    double totalData = 0;
    double totalMC = 0;

    TH1D* hData[4];
    TH1D* hMC[4];

    // First loop: Get histograms and calculate scale factor
    for (int i = 0; i < 4; i++) {
        int nTof = i + 2; // 2, 3, 4, 5
        TString histName = Form(baseName.Data(), nTof) + TString("_nominal");
        
        hData[i] = (TH1D*)fData->Get(histName);
        hMC[i]   = (TH1D*)fMC->Get(histName);

        if (!hData[i] || !hMC[i]) return; 

        if (rebin) {
            hData[i]->Rebin(2);
            hMC[i]->Rebin(2);
        }

        totalData += hData[i]->Integral();
        totalMC += hMC[i]->Integral();
    }

    double scale = (totalMC > 0) ? totalData / totalMC : 1.0;
    double maxY = 0;

    // Second loop: Style, scale, and add to stack
    for (int i = 0; i < 4; i++) {
        hMC[i]->Scale(scale);
        
        if (hData[i]->GetMaximum() > maxY) maxY = hData[i]->GetMaximum();
        if (hMC[i]->GetMaximum() > maxY) maxY = hMC[i]->GetMaximum();
        
        // MC Styling 
        hMC[i]->SetFillColorAlpha(colors[i], 0.3); 
        hMC[i]->SetLineColor(colors[i]);
        hMC[i]->SetLineWidth(2);
        //hsMC->GetYaxis()->SetMaxDigits(3);
        hsMC->Add(hMC[i]);

        // Data Styling
        hData[i]->SetMarkerStyle(20);
        hData[i]->SetMarkerSize(0.8);
        hData[i]->GetYaxis()->SetMaxDigits(3);
        hData[i]->SetLineColor(colors[i]);

        leg->AddEntry(hData[i], Form("Data %s", tofLabels[i].Data()), "pe");
        leg->AddEntry(hMC[i], Form("MC %s", tofLabels[i].Data()), "f");
    }
    hsMC->Draw("HIST NOSTACK");
    
    hsMC->GetXaxis()->SetTitle(xTitle);
    hsMC->GetYaxis()->SetTitle("events");
    hsMC->GetXaxis()->SetRangeUser(xMin, xMax);
    hsMC->SetMaximum(maxY * 1.3); 
    
    // Offsets to push the titles away from the larger labels
    hsMC->GetYaxis()->SetTitleOffset(1.1);
    hsMC->GetXaxis()->SetTitleOffset(1.1);
    
    for (int i = 0; i < 4; i++) {
        hData[i]->Draw("PE SAME");
    }
    
    leg->Draw();
    c->SaveAs(outName + ".pdf");
    delete c;
}

// Helper function for 2D Invariant Mass Plots (Slide 12)
void Plot2DMassWindow(TFile* file, TString dataType, TString outName) {
    TCanvas* c = new TCanvas("c_2d_"+outName, "2D Mass", 1000, 800);
    c->Divide(2, 2);

    for (int i = 0; i < 4; i++) {
        c->cd(i + 1);
        
        gPad->SetRightMargin(0.15); 
        gPad->SetLeftMargin(0.18);
        gPad->SetBottomMargin(0.15);

        int nTof = i + 2; 
        TString histName = Form("histInvMassPiPi2D%d_nominal", nTof);
        TH2D* h2 = (TH2D*)file->Get(histName);

        if (!h2) continue;

        h2->GetXaxis()->SetTitle("m_{#pi^{+}#pi^{-}}^{Leading} [GeV/c^{2}]");
        h2->GetYaxis()->SetTitle("m_{#pi^{+}#pi^{-}}^{Subleading} [GeV/c^{2}]");
        //h2->GetXaxis()->SetTitleOffset(1.9);
        h2->GetYaxis()->SetTitleOffset(0.0);
        h2->GetXaxis()->SetTitleSize(0.06); //axis title size
        h2->GetYaxis()->SetTitleSize(0.06);

        h2->GetXaxis()->SetLabelSize(0.05);
        h2->GetYaxis()->SetLabelSize(0.05);

        h2->SetStats(0); 
        
        h2->Draw("COLZ");

        TBox* cutBox = new TBox(0.48, 0.48, 0.52, 0.52);
        cutBox->SetFillStyle(0);      
        cutBox->SetLineColor(kRed);   
        cutBox->SetLineWidth(2);
        cutBox->SetLineStyle(2);      
        cutBox->Draw("SAME");

        // NEW: Draw "TOF = X" in the upper right corner
        TLatex latex;
        latex.SetNDC();
        latex.SetTextFont(62); // Bold font
        latex.SetTextSize(0.06);
        latex.SetTextColor(kRed+1);
        latex.SetTextAlign(33); // Align Top-Right
        
        // Coordinates: x=0.82 (just left of COLZ palette), y=0.88 (just below top frame)
        latex.DrawLatex(0.82, 0.88, Form("TOF = %d", nTof));
    }

    c->SaveAs(outName + ".pdf");
    delete c;
}

// Helper function for the invariant mass plot from Slide 20
void PlotFinalMass(TFile* fData, TFile* fMC) {
    TCanvas* c = new TCanvas("c_mass", "Mass", 800, 600);
    TH1D* hMass = (TH1D*)fData->Get("HistMassK0K0_nominal");
    TH1D* hMassMC = (TH1D*)fMC->Get("HistMassK0K0_nominal");
    
    if (!hMass) return;

    hMass->SetMarkerStyle(20);
    hMass->SetMarkerColor(kBlack);
    hMass->SetLineColor(kBlack);
    hMass->GetXaxis()->SetRangeUser(0.9, 3.0);
    hMass->GetXaxis()->SetTitle("m_{K_{S}^{0}K_{S}^{0}} [GeV/c^{2}]");
    hMass->GetYaxis()->SetTitle("events");
    //hMass->GetYaxis()->SetTitleOffset(1.1);
    //hMass->GetXaxis()->SetTitleOffset(1.1);
    hMass->GetXaxis()->SetTitleSize(0.06); //axis title size
    hMass->GetYaxis()->SetTitleSize(0.06);

    hMass->GetXaxis()->SetLabelSize(0.05);
    hMass->GetYaxis()->SetLabelSize(0.05);
    
    hMass->Draw("PE");
    //hMassMC->Draw("HIST SAME");

    TF1* fitFunc = new TF1("fitFunc", "gaus(0) + pol1(3)", 1.0, 2.5);
    fitFunc->SetParameters(15, 1.68, 0.08, 10, -2); 
    fitFunc->SetLineColor(kRed);
    hMass->Fit(fitFunc, "R"); 

    TLegend* leg = new TLegend(0.55, 0.75, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.045);
    leg->AddEntry(hMass, "Data", "pe");
    leg->AddEntry(fitFunc, "Fit: gaus + pol1", "l");
    leg->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.045);
    latex.DrawLatex(0.55, 0.68, Form("A = %.2f", fitFunc->GetParameter(0)));
    latex.DrawLatex(0.55, 0.62, Form("#mu = %.2f GeV/c^{2}", fitFunc->GetParameter(1)));
    latex.DrawLatex(0.55, 0.56, Form("#sigma = %.2f GeV/c^{2}", fitFunc->GetParameter(2)));

    c->SaveAs("Plot_K0K0_InvariantMass.pdf");
    delete c;
}


void ExclusivePlots() {
    SetMyStyle();

    TFile* fData = TFile::Open("DataExclusiveClassA.root");
    TFile* fMC = TFile::Open("MCExclusiveClassA.root");

    if (!fData || !fMC) {
        cout << "Error opening ROOT files!" << endl;
        return;
    }

    PlotStackedMultiplicity(fData, fMC, "histPtWith_%d", "p_{T} [GeV/c]", 0, 2.0, "Plot_pT_WithTOF", true);
    PlotStackedMultiplicity(fData, fMC, "histPtWithout_%d", "p_{T} [GeV/c]", 0, 2.0, "Plot_pT_WithoutTOF", true);
    PlotStackedMultiplicity(fData, fMC, "histEtaWith_%d", "#eta", -2.5, 2.5, "Plot_Eta_WithTOF");
    PlotStackedMultiplicity(fData, fMC, "histEtaWithout_%d", "#eta", -2.5, 2.5, "Plot_Eta_WithoutTOF");
    PlotStackedMultiplicity(fData, fMC, "histNfitWith_%d", "N_{fit}^{hit}", 10, 50, "Plot_Nfit_WithTOF");
    PlotStackedMultiplicity(fData, fMC, "histNfitWithout_%d", "N_{fit}^{hit}", 10, 50, "Plot_Nfit_WithoutTOF");
    
    // Slide 12: 2D Narrow Mass Window on mpipi
    Plot2DMassWindow(fData, "Data", "Plot_2DMassWindow_Data");
    Plot2DMassWindow(fMC, "MC", "Plot_2DMassWindow_MC");

    PlotStackedMultiplicity(fData, fMC, "histPtMissBefore%d", "p_{T}^{miss} [GeV/c]", 0, 0.5, "Plot_PtMiss");
    PlotStackedMultiplicity(fData, fMC, "histNTOFClusterBefore%d", "N_{TOF}^{clusters}", 0, 15, "Plot_NTOFClusters");
    PlotStackedMultiplicity(fData, fMC, "histDCADaughtersBefore%d", "DCA_{daughters} [cm]", 0, 5.0, "Plot_DCADaughters");
    PlotStackedMultiplicity(fData, fMC, "histDCABeamlineBefore%d", "DCA_{beamline} [cm]", 0, 5.0, "Plot_DCABeamline");
    PlotStackedMultiplicity(fData, fMC, "histCosBefore%d", "cos(#alpha_{p})", -1.0, 1.0, "Plot_CosPointingAngle");
    PlotStackedMultiplicity(fData, fMC, "histDecayBefore%d", "l_{decay} [cm]", 0, 10.0, "Plot_DecayLength");
    PlotStackedMultiplicity(fData, fMC, "histKsiEBefore%d", "#xi_{E}", -0.005, 0.02, "Plot_KsiE");
    PlotStackedMultiplicity(fData, fMC, "histKsiWBefore%d", "#xi_{W}", -0.005, 0.02, "Plot_KsiW");
    PlotStackedMultiplicity(fData, fMC, "histSumProtonMomentaXBefore%d", "p_{x}^{E} + p_{x}^{W} [GeV/c]", -1.0, 1.0, "Plot_SumPx");
    PlotStackedMultiplicity(fData, fMC, "histSumProtonMomentaYBefore%d", "p_{y}^{E} + p_{y}^{W} [GeV/c]", -1.0, 1.0, "Plot_SumPy");
    PlotStackedMultiplicity(fData, fMC, "histCorrKsiBefore%d", "m_{K_{S}^{0}K_{S}^{0}} / #sqrt{s} - #sqrt{#xi_{E}#xi_{W}}", -0.006, 0.006, "Plot_CorrKsi");
    PlotStackedMultiplicity(fData, fMC, "histZDiffBefore%d", "vtx_{z}^{leading} - vtx_{z}^{subleading} [cm]", -50, 50, "Plot_ZDiff");
    PlotStackedMultiplicity(fData, fMC, "histCosThetaStarBefore%d", "cos(#theta^{*})", -1.0, 1.0, "Plot_CosThetaStar");

    PlotFinalMass(fData, fMC);


    cout << "All plots successfully generated and saved as PDFs!" << endl;
}