#include <TF1.h>
#include <TF2.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TEllipse.h>
#include <TStyle.h>
#include <TLegend.h>

#include <stdio.h>
#include <string>

int main(int argc, char *argv[]){
    //something used so that the histograms would draw
    //https://stackoverflow.com/questions/30932725/painting-a-tcanvas-to-the-screen-in-a-compiled-root-cern-application
    TApplication theApp("App", &argc, argv);

    TFile *HistFile = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/AnaOutput_Inclusive_analysis_with_STUPCV0.root");

    TH2D *XiHist;
    HistFile->GetObject("Xi2DProtons", XiHist);
    XiHist->SetDirectory(0);
    TH2D *pHist;
    HistFile->GetObject("sump2DProtonsExact", pHist);
    pHist->SetDirectory(0);
    TH1D *mHist;
    HistFile->GetObject("Mpipi", mHist);
    mHist->SetDirectory(0);
    TH1D *mHistAfter;
    HistFile->GetObject("MpipiAfterElasticCut", mHistAfter);
    mHistAfter->SetDirectory(0);
    TH1D mHistDifference = *(TH1D *)mHist->Clone()-*(TH1D *)mHistAfter->Clone();
    TH2D *thetaHist;
    HistFile->GetObject("sumTheta2DProtonsExact", thetaHist);
    thetaHist->SetDirectory(0);

    HistFile->Close();

    //fitting
    TF2 *XiGauss = new TF2("XiGauss", "[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])", -0.01, 0.015, -0.01, 0.015);
    XiGauss->SetParameters(40, 0.005, 0.001, 0.005, 0.001);
    XiGauss->SetParNames("const", "x_0", "sigma_x", "y_0", "sigma_y");
    XiHist->Fit(XiGauss, "0");

    TF2 *pGauss = new TF2("pGauss", "[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4]) + [5]*TMath::Gaus(x,[6],[7])*TMath::Gaus(y,[8],[9])", -0.6, 1., -0.5, 0.5);
    pGauss->SetParameters(40, -0.04, 0.04, 0.0, 0.04, 5, -0.04, 0.2, 0.0, 0.2);
    pGauss->SetParNames("const_n", "x_0_n", "sigma_x_n", "y_0_n", "sigma_y_n", "const_w", "x_0_w", "sigma_x_w", "y_0_w", "sigma_y_w");
    pHist->Fit(pGauss, "0");

    TEllipse *XiEllipse = new TEllipse(XiGauss->GetParameter(1), XiGauss->GetParameter(3), 3*XiGauss->GetParameter(2), 3*XiGauss->GetParameter(4));
    XiEllipse->SetLineColor(2);
    XiEllipse->SetLineWidth(2);
    XiEllipse->SetFillStyle(0);

    TEllipse *pEllipse = new TEllipse(pGauss->GetParameter(1), pGauss->GetParameter(3), 3*pGauss->GetParameter(2), 3*pGauss->GetParameter(4));
    pEllipse->SetLineColor(2);
    pEllipse->SetLineWidth(2);
    pEllipse->SetFillStyle(0);

    //drawing new things
    TCanvas *c1 = new TCanvas("c1", "c1", 1600, 1200);
    // gStyle->SetOptStat(10);
    // gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    XiHist->SetTitle("");
    XiHist->Draw("colz");
    // XiEllipse->Draw("same");
    c1->SetLeftMargin(0.15);
    c1->SetRightMargin(0.15);
    c1->SetRealAspectRatio();
    c1->Draw();
    c1->Print("Xi.pdf");

    TCanvas *c2 = new TCanvas("c2", "c2", 1600, 1200);
    // gStyle->SetOptStat(10);
    // gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    pHist->SetTitle("");
    pHist->Draw("colz");
    // pEllipse->Draw("same");
    c2->SetRealAspectRatio();
    c2->Draw();
    c2->Print("p.pdf");

    TCanvas *c3 = new TCanvas("c3", "c3", 1600, 1200);
    gStyle->SetOptStat(0);
    mHist->SetMinimum(0);
    mHist->SetLineWidth(2);
    mHistAfter->SetLineWidth(2);
    mHistAfter->SetLineColor(2);
    mHistDifference.SetLineWidth(2);
    mHistDifference.SetLineColor(3);
    mHist->SetTitle("");
    mHist->Draw();
    mHistAfter->Draw("same");
    mHistDifference.Draw("same");
    TLegend legend = TLegend(0.65, 0.7, 0.89, 0.89);
    legend.AddEntry(mHist, "Before elastic cuts", "l");
    legend.AddEntry(mHistAfter, "After elastic cuts", "l");
    legend.AddEntry(&mHistDifference, "Difference", "l");
    legend.SetBorderSize(0);
    legend.Draw("same");
    c3->Draw();
    c3->Print("K0.pdf");

    TCanvas *c4 = new TCanvas("c4", "c4", 1600, 1200);
    // gStyle->SetOptStat(10);
    // gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    thetaHist->SetTitle(";#Sigma#theta_{x};#Sigma#theta_{y}");
    thetaHist->Draw("colz");
    // pEllipse->Draw("same");
    c4->SetLeftMargin(0.15);
    c4->SetRightMargin(0.15);
    c4->SetRealAspectRatio();
    c4->Draw();
    c4->Print("theta.pdf");

    theApp.Run();
    return 0;
}
