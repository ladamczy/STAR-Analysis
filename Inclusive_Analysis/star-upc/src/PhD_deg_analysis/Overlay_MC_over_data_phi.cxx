#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1D.h"

#include "string.h"

void draw_and_save(TH1D *data, TH1D *sim, std::string fileTitle, double x1, double y1, double x2, double y2);

int main(int argc, char const *argv[]){
    TFile *K0ptFile = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/Differential_crossection_phi_pt_new_data.root");
    TFile *K0etaFile = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/Differential_crossection_phi_eta_new_data.root");
    TFile *sim = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/Pythia_phi_simulation.root");

    TH1D *K0ptHistData = (TH1D *)K0ptFile->Get("newdataK01D");
    TH1D *K0ptHistSim = (TH1D *)sim->Get("PhiptHist");
    K0ptHistData->SetName("K0ptHistData");
    K0ptHistSim->SetName("K0ptHistSim");

    TH1D *K0etaHistData = (TH1D *)K0etaFile->Get("newdataK01D");
    TH1D *K0etaHistSim = (TH1D *)sim->Get("PhietaHist");
    K0etaHistData->SetName("K0etaHistData");
    K0etaHistSim->SetName("K0etaHistSim");

    draw_and_save(K0ptHistData, K0ptHistSim, "phiptHist.pdf", 0.65, 0.65, 0.85, 0.85);
    draw_and_save(K0etaHistData, K0etaHistSim, "phietaHist.pdf", 0.4, 0.21, 0.6, 0.41);

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
    legend->AddEntry(data->GetName(), "#bf{upcDST}");
    legend->AddEntry(sim->GetName(), "#bf{PYTHIA8 true level}");
    legend->SetBorderSize(0);
    legend->Draw("SAME");
    resultCanvas->SetLeftMargin(0.15);
    resultCanvas->SetRightMargin(0.05);
    gPad->Update();
    resultCanvas->SaveAs(("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/presentation_pdfs/"+fileTitle).c_str());
}
