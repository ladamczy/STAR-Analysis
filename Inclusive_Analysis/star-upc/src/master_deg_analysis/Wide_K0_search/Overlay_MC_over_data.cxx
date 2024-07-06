#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1D.h"

#include "string.h"

void draw_and_save(TH1D *data, TH1D *sim, std::string fileTitle, double x1, double y1, double x2, double y2);

int main(int argc, char const *argv[]){
    TFile *K0ptFile = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/Differential_crossection_K0_pt_new_data.root");
    TFile *K0etaFile = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/Differential_crossection_K0_eta_new_data.root");
    TFile *LambdaptFile = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/Differential_crossection_Lambda_pt_new_data.root");
    TFile *LambdaetaFile = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/Differential_crossection_Lambda_eta_new_data.root");
    TFile *LambdaBarptFile = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/Differential_crossection_LambdaBar_pt_new_data.root");
    TFile *LambdaBaretaFile = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/Differential_crossection_LambdaBar_eta_new_data.root");
    TFile *sim = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/Pythia_simulation.root");

    TH1D *K0ptHistData = (TH1D *)K0ptFile->Get("newdataK01D");
    TH1D *K0ptHistSim = (TH1D *)sim->Get("K0ptHist");
    K0ptHistData->SetName("K0ptHistData");
    K0ptHistSim->SetName("K0ptHistSim");

    TH1D *K0etaHistData = (TH1D *)K0etaFile->Get("newdataK01D");
    TH1D *K0etaHistSim = (TH1D *)sim->Get("K0etaHist");
    K0etaHistData->SetName("K0etaHistData");
    K0etaHistSim->SetName("K0etaHistSim");

    TH1D *LambdaptHistData = (TH1D *)LambdaptFile->Get("newdataLambda1D");
    TH1D *LambdaptHistSim = (TH1D *)sim->Get("LambdaptHist");
    LambdaptHistData->SetName("LambdaptHistData");
    LambdaptHistSim->SetName("LambdaptHistSim");

    TH1D *LambdaetaHistData = (TH1D *)LambdaetaFile->Get("newdataLambda1D");
    TH1D *LambdaetaHistSim = (TH1D *)sim->Get("LambdaetaHist");
    LambdaetaHistData->SetName("LambdaetaHistData");
    LambdaetaHistSim->SetName("LambdaetaHistSim");

    TH1D *LambdaBarptHistData = (TH1D *)LambdaBarptFile->Get("newdataLambda1D");
    TH1D *LambdaBarptHistSim = (TH1D *)sim->Get("LambdaBarptHist");
    LambdaBarptHistData->SetName("LambdaBarptHistData");
    LambdaBarptHistSim->SetName("LambdaBarptHistSim");

    TH1D *LambdaBaretaHistData = (TH1D *)LambdaBaretaFile->Get("newdataLambda1D");
    TH1D *LambdaBaretaHistSim = (TH1D *)sim->Get("LambdaBaretaHist");
    LambdaBaretaHistData->SetName("LambdaBaretaHistData");
    LambdaBaretaHistSim->SetName("LambdaBaretaHistSim");

    draw_and_save(K0ptHistData, K0ptHistSim, "K0ptHist.pdf", 0.65, 0.65, 0.85, 0.85);
    draw_and_save(K0etaHistData, K0etaHistSim, "K0etaHist.pdf", 0.4, 0.21, 0.6, 0.41);
    draw_and_save(LambdaptHistData, LambdaptHistSim, "LambdaptHist.pdf", 0.65, 0.65, 0.85, 0.85);
    draw_and_save(LambdaetaHistData, LambdaetaHistSim, "LambdaetaHist.pdf", 0.4, 0.21, 0.6, 0.41);
    draw_and_save(LambdaBarptHistData, LambdaBarptHistSim, "LambdaBarptHist.pdf", 0.65, 0.65, 0.85, 0.85);
    draw_and_save(LambdaBaretaHistData, LambdaBaretaHistSim, "LambdaBaretaHist.pdf", 0.6, 0.21, 0.8, 0.41);

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
    legend->AddEntry(sim->GetName(), "#bf{PYTHIA simuation}");
    legend->SetBorderSize(0);
    legend->Draw("SAME");
    resultCanvas->SetLeftMargin(0.15);
    resultCanvas->SetRightMargin(0.05);
    gPad->Update();
    resultCanvas->SaveAs(("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/"+fileTitle).c_str());
}
