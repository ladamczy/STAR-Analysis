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
#include <TPad.h>
#include <TPaveStats.h>

#include <stdio.h>
#include <string>

int main(int argc, char *argv[]){
    //something used so that the histograms would draw
    //https://stackoverflow.com/questions/30932725/painting-a-tcanvas-to-the-screen-in-a-compiled-root-cern-application
    TApplication theApp("App", &argc, argv);

    TFile *HistFile = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/AnaOutput_Inclusive_zerobias_background_analysis_with_STUPCV0.root");

    TH2D *XiHist;
    HistFile->GetObject("Xi2DProtons", XiHist);
    XiHist->SetDirectory(0);
    TH2D *pHist;
    HistFile->GetObject("sump2DProtonsExact", pHist);
    pHist->SetDirectory(0);
    TH2D *thetaHist;
    HistFile->GetObject("sumTheta2DProtons", thetaHist);
    thetaHist->SetDirectory(0);

    HistFile->Close();

    //fitting
    TF2 *XiGauss = new TF2("XiGauss", "[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])", -0.01, 0.015, -0.01, 0.015);
    XiGauss->SetParameters(40, 0.005, 0.001, 0.005, 0.001);
    XiGauss->SetParNames("const", "x_{0}", "#sigma_{x}", "y_{0}", "#sigma_{y}");
    XiHist->Fit(XiGauss, "0R");

    TF2 *pGauss = new TF2("pGauss", "[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])", -0.2, 0.2, -0.2, 0.2);
    pGauss->SetParameters(40, -0.04, 0.04, 0.0, 0.04);
    pGauss->SetParNames("const", "x_{0}", "#sigma_{x}", "y_{0}", "#sigma_{y}");
    pHist->Fit(pGauss, "0R");

    TF2 *thetaGauss = new TF2("thetaGauss", "[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])", -0.0006, 0.0006, -0.0006, 0.0006);
    thetaGauss->SetParameters(40, -0.0, 0.0002, 0.0, 0.0002);
    thetaGauss->SetParNames("const", "x_{0}", "#sigma_{x}", "y_{0}", "#sigma_{y}");
    thetaHist->Fit(thetaGauss, "0");

    TEllipse *XiEllipse = new TEllipse(XiGauss->GetParameter(1), XiGauss->GetParameter(3), 3*XiGauss->GetParameter(2), 3*XiGauss->GetParameter(4));
    XiEllipse->SetLineColor(2);
    XiEllipse->SetLineWidth(2);
    XiEllipse->SetFillStyle(0);

    TEllipse *pEllipse = new TEllipse(pGauss->GetParameter(1), pGauss->GetParameter(3), 3*pGauss->GetParameter(2), 3*pGauss->GetParameter(4));
    pEllipse->SetLineColor(2);
    pEllipse->SetLineWidth(2);
    pEllipse->SetFillStyle(0);

    TEllipse *thetaEllipse = new TEllipse(thetaGauss->GetParameter(1), thetaGauss->GetParameter(3), 3*thetaGauss->GetParameter(2), 3*thetaGauss->GetParameter(4));
    thetaEllipse->SetLineColor(2);
    thetaEllipse->SetLineWidth(2);
    thetaEllipse->SetFillStyle(0);

    //drawing new things
    TCanvas *c1 = new TCanvas("c1", "c1", 1600, 1200);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1);
    XiHist->SetTitle("");
    XiHist->Draw("colz");
    XiEllipse->Draw("same");
    c1->SetLeftMargin(0.15);
    c1->SetRightMargin(0.15);
    c1->SetRealAspectRatio();
    gPad->Update();
    TPaveStats *st = (TPaveStats *)XiHist->FindObject("stats");
    st->SetX1NDC(0.52);
    st->SetX2NDC(0.84);
    st->SetY1NDC(0.55);
    st->SetY2NDC(0.89);
    st->SetBorderSize(0);
    c1->Draw();
    c1->Print("Xi_zerobias.pdf");


    TCanvas *c2 = new TCanvas("c2", "c2", 1600, 1200);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1);
    pHist->SetTitle("");
    pHist->Draw("colz");
    pEllipse->Draw("same");
    c2->SetRealAspectRatio();
    gPad->Update();
    TPaveStats *st2 = (TPaveStats *)pHist->FindObject("stats");
    st2->SetX1NDC(0.6);
    st2->SetX2NDC(0.89);
    st2->SetY1NDC(0.7);
    st2->SetY2NDC(0.89);
    st2->SetBorderSize(0);
    c2->Draw();
    c2->Print("p_zerobias.pdf");

    TCanvas *c3 = new TCanvas("c3", "c3", 1600, 1200);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1);
    thetaHist->SetTitle(";#Sigma#theta_{x};#Sigma#theta_{y}");
    thetaHist->Draw("colz");
    thetaEllipse->Draw("same");
    c3->SetLeftMargin(0.15);
    c3->SetRightMargin(0.15);
    c3->SetRealAspectRatio();
    gPad->Update();
    TPaveStats *st3 = (TPaveStats *)thetaHist->FindObject("stats");
    st3->SetX1NDC(0.58);
    st3->SetX2NDC(0.84);
    st3->SetY1NDC(0.68);
    st3->SetY2NDC(0.89);
    st3->SetBorderSize(0);
    c3->Draw();
    c3->Print("theta_zerobias.pdf");

    theApp.Run();
    return 0;
}
