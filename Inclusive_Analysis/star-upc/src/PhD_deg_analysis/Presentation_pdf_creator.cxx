#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStyle.h"

#include "string.h"

void draw_and_save(TH1D *data, TH1D *sim, std::string fileTitle, double x1, double y1, double x2, double y2);
void simpler_draw_and_save(TH1 *hist, bool is3D = false, bool isYlogarithmic = false, double normalisation = 1);

int main(int argc, char const *argv[]){
    // TFile *sim = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/Pythia_simulation.root");
    TFile *data = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/AnaOutput_Histograms_after_RP_cuts.root");

    TH1D *DataFigure3_4a = (TH1D *)data->Get("DataFigure3_4a");
    TH1D *DataFigure3_4b = (TH1D *)data->Get("DataFigure3_4b");
    TH1D *DataFigure3_5a = (TH1D *)data->Get("DataFigure3_5a");
    TH1D *DataFigure3_5b = (TH1D *)data->Get("DataFigure3_5b");
    TH1D *DataFigure3_5c = (TH1D *)data->Get("DataFigure3_5c");
    TH1D *DataFigure3_5d = (TH1D *)data->Get("DataFigure3_5d");
    TH1D *DataFigure3_6ab = (TH1D *)data->Get("DataFigure3_6ab");
    TH1D *DataFigure3_6c = (TH1D *)data->Get("DataFigure3_6c");
    TH1D *DataFigure3_6d = (TH1D *)data->Get("DataFigure3_6d");
    TH1D *DataFigure3_18g = (TH1D *)data->Get("DataFigure3_18g");
    TH1D *DataFigure3_19a = (TH1D *)data->Get("DataFigure3_19a");
    TH1D *DataFigure3_20 = (TH1D *)data->Get("DataFigure3_20");
    TH1D *DataFigure3_21 = (TH1D *)data->Get("DataFigure3_21");
    TH2D *DataFigure3_60 = (TH2D *)data->Get("DataFigure3_60");
    TH1D *DataFigure3_62a = (TH1D *)data->Get("DataFigure3_62a");
    TH1D *DataFigure3_62b = (TH1D *)data->Get("DataFigure3_62b");
    TH1D *DataFigure3_62c = (TH1D *)data->Get("DataFigure3_62c");
    TH1D *DataFigure3_62d = (TH1D *)data->Get("DataFigure3_62d");

    simpler_draw_and_save(DataFigure3_4a);
    simpler_draw_and_save(DataFigure3_4b);
    simpler_draw_and_save(DataFigure3_5a);
    simpler_draw_and_save(DataFigure3_5b);
    simpler_draw_and_save(DataFigure3_5c);
    simpler_draw_and_save(DataFigure3_5d);
    simpler_draw_and_save(DataFigure3_6ab);
    simpler_draw_and_save(DataFigure3_6c);
    simpler_draw_and_save(DataFigure3_6d);
    simpler_draw_and_save(DataFigure3_18g, false, false, 1/DataFigure3_18g->Integral());
    simpler_draw_and_save(DataFigure3_19a, false, false, 1/DataFigure3_19a->Integral());
    simpler_draw_and_save(DataFigure3_20, false, false, 1/DataFigure3_20->Integral());
    simpler_draw_and_save(DataFigure3_21, false, false, 1/DataFigure3_21->Integral());
    simpler_draw_and_save(DataFigure3_60, true);
    simpler_draw_and_save(DataFigure3_62a, false, true);
    simpler_draw_and_save(DataFigure3_62b, false, true);
    simpler_draw_and_save(DataFigure3_62c, false, true);
    simpler_draw_and_save(DataFigure3_62d, false, true);

    return 0;
}

void draw_and_save(TH1D *data, TH1D *sim, std::string fileTitle, double x1, double y1, double x2, double y2){
    TCanvas *resultCanvas = new TCanvas("resultCanvas", "resultCanvas", 1800, 1600);
    sim->Scale(data->Integral()/sim->GetEntries());
    data->SetMinimum(0);
    data->SetMaximum(std::max(data->GetMaximum(), sim->GetMaximum())*1.1);
    data->SetMarkerStyle(kFullCircle);
    data->SetMarkerSize(2);
    data->SetMarkerColor(kBlue);
    data->SetLineColor(kBlue+2);
    data->Draw("E");
    sim->SetMarkerColor(kRed);
    sim->SetLineColor(kRed+2);
    // sim->SetLineStyle(k); do lepszego końca linii
    sim->Draw("Hist same");
    TLegend *legend = new TLegend(x1, y1, x2, y2);
    legend->SetTextSize(0.025);
    legend->SetHeader("#bf{pp, #sqrt{s} = 510 GeV}", "C");
    legend->AddEntry(data->GetName(), "#bf{extended upcDST}");
    legend->AddEntry(sim->GetName(), "#bf{PYTHIA8 true level}");
    legend->SetBorderSize(0);
    legend->Draw("SAME");
    resultCanvas->SetLeftMargin(0.15);
    resultCanvas->SetRightMargin(0.05);
    gPad->Update();
    resultCanvas->SaveAs(("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/"+fileTitle).c_str());
}

void simpler_draw_and_save(TH1 *hist, bool is3D, bool isYlogarithmic, double normalisation){
    TCanvas *resultCanvas = new TCanvas("resultCanvas", "resultCanvas", 1800, 1600);
    gStyle->SetFrameLineWidth(1);
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);
    resultCanvas->SetLogy(isYlogarithmic);
    resultCanvas->SetLeftMargin(0.15);
    std::string fileTitle = hist->GetName();
    hist->SetMinimum(0);
    if(isYlogarithmic){
        hist->SetMinimum(0.1);
    }
    // hist->SetMarkerStyle(kFullCircle);
    // hist->SetMarkerSize(2);
    // hist->SetMarkerColor(kBlue);
    hist->SetLineColor(kBlue+2);
    if(is3D){
        resultCanvas->SetLogz();
        hist->Draw("colz");
    } else{
        hist->Draw("hist");
    }
    resultCanvas->SaveAs(("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/presentation_pdfs/"+fileTitle+".pdf").c_str());

}
