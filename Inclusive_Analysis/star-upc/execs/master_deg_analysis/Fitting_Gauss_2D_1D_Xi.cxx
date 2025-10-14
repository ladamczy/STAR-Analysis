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
#include <TROOT.h>

#include <stdio.h>
#include <string>

#include "MyStyles.h"

int main(int argc, char *argv[]){
    //something used so that the histograms would draw
    //https://stackoverflow.com/questions/30932725/painting-a-tcanvas-to-the-screen-in-a-compiled-root-cern-application
    // TApplication theApp("App", &argc, argv);

    MyStyles styleLibrary;
    TStyle mystyle = styleLibrary.Hist2DNormalSize(false);
    mystyle.cd();
    gROOT->ForceStyle();

    TFile *HistDataFile = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/AnaOutput_Inclusive_analysis_with_STUPCV0.root");
    TFile *HistZerobiasFile = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/AnaOutput_Inclusive_zerobias_background_analysis_with_STUPCV0.root");

    TH1D *XiEDataHist;
    TH1D *XiEDataHistTemp;
    TH1D *XiWDataHist;
    TH1D *XiWDataHistTemp;
    TH1D *XiEZerobiasHist;
    TH1D *XiWZerobiasHist;
    HistDataFile->GetObject("XiEprotoncloser", XiEDataHistTemp);
    HistDataFile->GetObject("XiWprotoncloser", XiWDataHistTemp);
    HistZerobiasFile->GetObject("XiEprotoncloser", XiEZerobiasHist);
    HistZerobiasFile->GetObject("XiWprotoncloser", XiWZerobiasHist);
    XiEDataHistTemp->SetDirectory(0);
    XiWDataHistTemp->SetDirectory(0);
    XiEZerobiasHist->SetDirectory(0);
    XiWZerobiasHist->SetDirectory(0);

    HistDataFile->Close();
    HistZerobiasFile->Close();

    double maxXiEData, maxXiWData, maxXiEZerobias, maxXiWZerobias;
    XiEDataHist = (TH1D *)XiEDataHistTemp->Clone();
    XiWDataHist = (TH1D *)XiWDataHistTemp->Clone();
    XiEDataHistTemp->GetXaxis()->SetRange(0, 112);
    XiWDataHistTemp->GetXaxis()->SetRange(0, 112);
    maxXiEData = XiEDataHistTemp->GetMaximum();
    maxXiWData = XiWDataHistTemp->GetMaximum();
    maxXiEZerobias = XiEZerobiasHist->GetMaximum();
    maxXiWZerobias = XiWZerobiasHist->GetMaximum();

    TF1 *functionToMultiply = new TF1("functionToMultiply", "[0]", -100, 100);

    functionToMultiply->SetParameter(0, maxXiEData/maxXiEZerobias);
    XiEZerobiasHist->Multiply(functionToMultiply);
    functionToMultiply->SetParameter(0, maxXiWData/maxXiWZerobias);
    XiWZerobiasHist->Multiply(functionToMultiply);

    TCanvas* c1 = new TCanvas("c1", "c1", 1600, 1200);
    XiEDataHist->SetTitle("");
    XiEDataHist->Draw();
    XiEZerobiasHist->SetLineColor(2);
    XiEZerobiasHist->Draw("same");
    //the legend
    TLegend legend1 = TLegend(0.5, 0.75, 0.89, 0.94);
    legend1.AddEntry(XiEDataHist, "pp collision data", "l");
    legend1.AddEntry(XiEZerobiasHist, "Normalised zerobias trigger", "l");
    legend1.Draw("same");
    c1->Draw();
    c1->Print("~/XiE.pdf");

    TCanvas* c2 = new TCanvas("c2", "c2", 1600, 1200);
    XiWDataHist->SetTitle("");
    XiWDataHist->Draw();
    XiWZerobiasHist->SetLineColor(2);
    XiWZerobiasHist->Draw("same");
    //the legend
    TLegend legend2 = TLegend(0.5, 0.75, 0.89, 0.94);
    legend2.AddEntry(XiWDataHist, "pp collision data", "l");
    legend2.AddEntry(XiWZerobiasHist, "Normalised zerobias trigger", "l");
    legend2.Draw("same");
    c2->Draw();
    c2->Print("~/XiW.pdf");

    // theApp.Run();
    return 0;
}
