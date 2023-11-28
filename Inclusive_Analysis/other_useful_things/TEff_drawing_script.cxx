#include "TEfficiency.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"

#include <string>

int TEff_drawing_script(TEfficiency *TEff, const char *beforeLegendName, const char *afterLegendName){
    if(TEff->GetTotalHistogram()->GetDimension()==1){
        TH1D *beforeHist = (TH1D *)TEff->GetTotalHistogram();
        TH1D *afterHist = (TH1D *)TEff->GetPassedHistogram();
        beforeHist->SetLineWidth(2);
        afterHist->SetLineWidth(2);
        afterHist->SetLineColor(2);
        gStyle->SetOptStat(0);

        TCanvas c1("c1", "", 1600, 900);
        c1.SetGrid();

        //data
        TLegend legend = TLegend(0.65, 0.7, 0.89, 0.89);
        legend.AddEntry(beforeHist, beforeLegendName, "l");
        legend.AddEntry(afterHist, afterLegendName, "l");
        legend.SetBorderSize(0);

        beforeHist->Draw();
        afterHist->Draw("same");
        legend.Draw("same");
        c1.Print((std::string(TEff->GetName())+"data.pdf").c_str());
        // c1.SaveAs((std::string(TEff->GetName())+"data.root").c_str());

        //ratio
        c1.Clear();
        beforeHist->GetYaxis()->SetTitle("ratio");
        afterHist->GetYaxis()->SetTitle("ratio");
        TEfficiency *tempEff = new TEfficiency(*afterHist, *beforeHist);
        TEff->SetMarkerColor(kBlue);
        TEff->SetMarkerStyle(7);
        TEff->Draw("EP");
        c1.Update();
        TEff->GetPaintedGraph()->SetMinimum(0.0);
        TEff->GetPaintedGraph()->GetYaxis()->SetTitle("ratio");
        c1.Modified();
        c1.Update();
        c1.Print((std::string(TEff->GetName())+"ratio.pdf").c_str());
        // c1.SaveAs((std::string(TEff->GetName())+"ratio.root").c_str());
        delete tempEff;
    } else if(TEff->GetTotalHistogram()->GetDimension()==2){
        TH2D *beforeHist = (TH2D *)TEff->GetTotalHistogram();
        TH2D *afterHist = (TH2D *)TEff->GetPassedHistogram();
        gStyle->SetOptStat(0);

        TCanvas c1("c1", "", 1600, 900);
        c1.SetGrid();

        beforeHist->SetMinimum(0.0);
        beforeHist->Draw("colz");
        c1.Print((std::string(TEff->GetName())+"databefore.pdf").c_str());
        c1.Clear();
        afterHist->SetMinimum(0.0);
        afterHist->Draw("colz");
        c1.Print((std::string(TEff->GetName())+"dataafter.pdf").c_str());
        // c1.SaveAs((std::string(TEff->GetName())+"data.root").c_str());

        //ratio
        c1.Clear();
        TEff->Draw("colz");
        c1.Update();
        TEff->GetPaintedHistogram()->SetMinimum(0.0);
        TEff->GetPaintedHistogram()->SetMaximum(1.0);
        c1.Print((std::string(TEff->GetName())+"ratio.pdf").c_str());
        // c1.SaveAs((std::string(TEff->GetName())+"ratio.root").c_str());
    }

    return 0;
}