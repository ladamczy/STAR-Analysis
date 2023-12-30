#include <TF1.h>
#include <TF2.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TApplication.h>

#include <stdio.h>
#include <string>

int main(int argc, char *argv[]){
    //something used so that the histograms would draw
    //https://stackoverflow.com/questions/30932725/painting-a-tcanvas-to-the-screen-in-a-compiled-root-cern-application
    TApplication theApp("App", &argc, argv);

    TFile *HistFile = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/AnaOutput_Inclusive_analysis_with_STUPCV0.root");

    TH2D *XiHist;
    TH2D *pHist;
    HistFile->GetObject("Xi2DProtons", XiHist);
    HistFile->GetObject("sump2DProtonsExact", pHist);
    XiHist->SetDirectory(0);
    pHist->SetDirectory(0);

    HistFile->Close();

    //fitting
    TF2 *XiGauss = new TF2("XiGauss", "[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])", -0.01, 0.01, -0.01, 0.01);
    XiGauss->SetParameters(40, 0.005, 0.001, 0.005, 0.001);
    XiGauss->SetParNames("const", "x_0", "sigma_x", "y_0", "sigma_y");
    // XiHist->Fit(XiGauss, "0");
    XiHist->Fit(XiGauss, "", "colz");

    TF2 *pGauss = new TF2("pGauss", "[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4]) + [5]*TMath::Gaus(x,[6],[7])*TMath::Gaus(y,[8],[9])", -0.6, 1., -0.5, 0.5);
    pGauss->SetParameters(40, -0.04, 0.04, 0.0, 0.04, 5, -0.04, 0.2, 0.0, 0.2);
    pGauss->SetParNames("const_n", "x_0_n", "sigma_x_n", "y_0_n", "sigma_y_n", "const_w", "x_0_w", "sigma_x_w", "y_0_w", "sigma_y_w");
    pHist->Fit(pGauss, "0");

    //drawing new things
    // TH1D *XiHistXprofile = XiHist->ProjectionX("XiHistXprofile");
    // TF1 *XiGaussXprofile = new TF1("XiGaussXprofile", "gaus", -0.01, 0.02);
    // XiGaussXprofile->SetParameter(0, XiGauss->GetParameter(0));
    // XiGaussXprofile->SetParameter(1, XiGauss->GetParameter(1));
    // XiGaussXprofile->SetParameter(2, XiGauss->GetParameter(2));
    // XiHistXprofile->Draw();
    // XiGaussXprofile->Draw("same");

    theApp.Run();
    return 0;
}
