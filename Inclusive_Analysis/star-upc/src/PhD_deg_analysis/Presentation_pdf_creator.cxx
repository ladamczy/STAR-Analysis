#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TEfficiency.h"
#include "TLatex.h"
#include "TStyle.h"

#include "string.h"

void draw_and_save(TH1D *data, TH1D *sim, std::string fileTitle, double x1, double y1, double x2, double y2);
void simpler_draw_and_save(TH1 *hist, std::string filename = "", bool is3D = false, bool isYlogarithmic = false);
void simpler_draw_and_save(TEfficiency *hist, std::string filename = "");

int main(int argc, char const *argv[]){
    TFile *pythia = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/Pythia_presentation_simulation.root");
    TFile *data = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/AnaOutput_Histograms_after_RP_cuts.root");
    TFile *pairs = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/AnaOutput_Inclusive_analysis_Kstar_phi_old_data.root");

    //data
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
    TH1D *DataFigure3_18gcloser = (TH1D *)data->Get("DataFigure3_18gcloser");
    TH1D *DataFigure3_19a = (TH1D *)data->Get("DataFigure3_19a");
    TH1D *DataFigure3_20 = (TH1D *)data->Get("DataFigure3_20");
    TH1D *DataFigure3_21 = (TH1D *)data->Get("DataFigure3_21");
    TH2D *DataFigure3_60 = (TH2D *)data->Get("DataFigure3_60");
    TH1D *DataFigure3_62a = (TH1D *)data->Get("DataFigure3_62a");
    TH1D *DataFigure3_62b = (TH1D *)data->Get("DataFigure3_62b");
    TH1D *DataFigure3_62c = (TH1D *)data->Get("DataFigure3_62c");
    TH1D *DataFigure3_62d = (TH1D *)data->Get("DataFigure3_62d");
    TH1D *DataControl = (TH1D *)data->Get("DataControl");
    //Pythia
    TH1D *PythiaFigure3_6ab = (TH1D *)pythia->Get("PythiaFigure3_6ab");
    TH1D *PythiaFigure3_6c = (TH1D *)pythia->Get("PythiaFigure3_6c");
    TH1D *PythiaFigure3_6d = (TH1D *)pythia->Get("PythiaFigure3_6d");
    TEfficiency *PythiaFigure3_7 = (TEfficiency *)pythia->Get("PythiaFigure3_7");
    TH1D *PythiaFigure3_18g = (TH1D *)pythia->Get("PythiaFigure3_18g");
    TH1D *PythiaFigure3_18gcloser = (TH1D *)pythia->Get("PythiaFigure3_18gcloser");
    TH1D *PythiaFigure3_19a = (TH1D *)pythia->Get("PythiaFigure3_19a");
    TH1D *PythiaFigure3_20 = (TH1D *)pythia->Get("PythiaFigure3_20");
    TH1D *PythiaFigure3_21 = (TH1D *)pythia->Get("PythiaFigure3_21");
    TH1D *PythiaControl = (TH1D *)pythia->Get("PythiaControl");
    //pairs
    TH1D *MKKWide = (TH1D *)pairs->Get("MKKWide");
    TH1D *MKpiWide = (TH1D *)pairs->Get("MKpiWide");
    TH1D *MpipiWide = (TH1D *)pairs->Get("MpipiWide");
    TH1D *MKKWidedEdx = (TH1D *)pairs->Get("MKKWidedEdx");
    TH1D *MKpiWidedEdx = (TH1D *)pairs->Get("MKpiWidedEdx");
    TH1D *MpipiWidedEdx = (TH1D *)pairs->Get("MpipiWidedEdx");

    //preprocessing of some histograms
    int Data_N_events_after_nsel = DataControl->GetBinContent(2);
    int Pythia_N_events_after_nsel = PythiaControl->GetBinContent(2);
    //Data
    DataFigure3_18g->Scale(1/DataFigure3_18g->GetBinWidth(1)/Data_N_events_after_nsel);
    DataFigure3_18gcloser->Scale(1/DataFigure3_18gcloser->GetBinWidth(1)/Data_N_events_after_nsel);
    DataFigure3_20->Scale(1/DataFigure3_20->GetBinWidth(1)/Data_N_events_after_nsel);
    DataFigure3_21->Scale(1/DataFigure3_21->GetBinWidth(1)/Data_N_events_after_nsel);
    DataFigure3_62a->Scale(1/DataFigure3_62a->GetBinWidth(1)/Data_N_events_after_nsel);
    DataFigure3_62b->Scale(1/DataFigure3_62b->GetBinWidth(1)/Data_N_events_after_nsel);
    DataFigure3_62c->Scale(1/DataFigure3_62c->GetBinWidth(1)/Data_N_events_after_nsel);
    DataFigure3_62d->Scale(1/DataFigure3_62d->GetBinWidth(1)/Data_N_events_after_nsel);
    //Pythia
    PythiaFigure3_6ab->Scale(DataFigure3_6ab->Integral()/PythiaFigure3_6ab->Integral());
    PythiaFigure3_6c->Scale(DataFigure3_6c->Integral()/PythiaFigure3_6c->Integral());
    PythiaFigure3_6d->Scale(DataFigure3_6d->Integral()/PythiaFigure3_6d->Integral());
    PythiaFigure3_18g->Scale(1/PythiaFigure3_18g->GetBinWidth(1)/Pythia_N_events_after_nsel);
    PythiaFigure3_18gcloser->Scale(1/PythiaFigure3_18gcloser->GetBinWidth(1)/Pythia_N_events_after_nsel);
    PythiaFigure3_20->Scale(1/PythiaFigure3_20->GetBinWidth(1)/Pythia_N_events_after_nsel);
    PythiaFigure3_21->Scale(1/PythiaFigure3_21->GetBinWidth(1)/Pythia_N_events_after_nsel);

    //pairs drawing and labeling

    //drawing and saving - data &PYTHIA
    simpler_draw_and_save(DataFigure3_4a, "Figure3_4a");
    simpler_draw_and_save(DataFigure3_4b, "Figure3_4b");
    simpler_draw_and_save(DataFigure3_5a, "Figure3_5a");
    simpler_draw_and_save(DataFigure3_5b, "Figure3_5b");
    simpler_draw_and_save(DataFigure3_5c, "Figure3_5c");
    simpler_draw_and_save(DataFigure3_5d, "Figure3_5d");
    draw_and_save(DataFigure3_6ab, PythiaFigure3_6ab, "Figure3_6ab", 0.4, 0.2, 0.6, 0.4);
    draw_and_save(DataFigure3_6c, PythiaFigure3_6c, "Figure3_6c", 0.64, 0.69, 0.84, 0.89);
    draw_and_save(DataFigure3_6d, PythiaFigure3_6d, "Figure3_6d", 0.4, 0.2, 0.6, 0.4);
    simpler_draw_and_save(PythiaFigure3_7, "Figure3_7");
    draw_and_save(DataFigure3_18g, PythiaFigure3_18g, "Figure3_18g", 0.64, 0.69, 0.84, 0.89);
    draw_and_save(DataFigure3_18gcloser, PythiaFigure3_18gcloser, "Figure3_18gcloser", 0.5, 0.2, 0.7, 0.4);
    draw_and_save(DataFigure3_19a, PythiaFigure3_19a, "Figure3_19a", 0.64, 0.69, 0.84, 0.89);
    draw_and_save(DataFigure3_20, PythiaFigure3_20, "Figure3_20", 0.64, 0.69, 0.84, 0.89);
    draw_and_save(DataFigure3_21, PythiaFigure3_21, "Figure3_21", 0.42, 0.2, 0.62, 0.4);
    simpler_draw_and_save(DataFigure3_60, "Figure3_60", true);
    simpler_draw_and_save(DataFigure3_62a, "Figure3_62a", false, true);
    simpler_draw_and_save(DataFigure3_62b, "Figure3_62b", false, true);
    simpler_draw_and_save(DataFigure3_62c, "Figure3_62c", false, true);
    simpler_draw_and_save(DataFigure3_62d, "Figure3_62d", false, true);

    //pairs drawing
    //setting
    TCanvas *resultCanvas = new TCanvas("resultCanvas", "resultCanvas", 1800, 1600);
    gStyle->SetFrameLineWidth(1);
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);
    resultCanvas->SetLeftMargin(0.15);
    //KK without PID
    MKKWide->SetTitle("");
    MKKWide->SetMinimum(0);
    MKKWide->SetLineColor(kBlue+2);
    MKKWide->Draw("hist");
    TLatex newText(1.16, 23500, "Unidentified resonance (~1070)");
    newText.SetTextColor(kRed);
    newText.SetTextSize(0.03);
    newText.Draw("same");
    TLegend *legend = new TLegend(0.64, 0.59, 0.84, 0.79);
    legend->SetTextSize(0.025);
    legend->SetHeader("#bf{pp, #sqrt{s} = 510 GeV}", "C");
    legend->AddEntry(MKKWide->GetName(), "#bf{upcDST data}");
    legend->SetBorderSize(0);
    legend->DrawClone("SAME");
    resultCanvas->SaveAs("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/presentation_pdfs/MKKWide.pdf");
    resultCanvas->Clear();
    //KK with PID
    MKKWidedEdx->SetTitle("");
    MKKWidedEdx->SetMinimum(0);
    MKKWidedEdx->SetLineColor(kBlue+2);
    MKKWidedEdx->Draw("hist");
    newText.SetText(0.38, 840, "#phi(1020)");
    newText.DrawClone("same");
    newText.SetText(1.15, 880, "Unidentified resonance (~1070)");
    newText.DrawClone("same");
    newText.SetText(1.3, 740, "Unidentified resonance (~a)");
    newText.DrawClone("same");
    newText.SetText(1.7, 660, "Unidentified resonance (~f)");
    newText.DrawClone("same");
    delete legend;
    legend = new TLegend(0.64, 0.39, 0.84, 0.59);
    legend->SetTextSize(0.025);
    legend->SetHeader("#bf{pp, #sqrt{s} = 510 GeV}", "C");
    legend->AddEntry(MKKWidedEdx->GetName(), "#bf{upcDST data}");
    legend->SetBorderSize(0);
    legend->DrawClone("SAME");
    resultCanvas->SaveAs("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/presentation_pdfs/MKKWidedEdx.pdf");
    resultCanvas->Clear();
    //Kpi without PID
    MKpiWide->SetTitle("");
    MKpiWide->SetMinimum(0);
    MKpiWide->SetLineColor(kBlue+2);
    MKpiWide->Draw("hist");
    newText.SetText(0.9, 21000, "#splitline{Unidentified}{resonance (~730)}");
    newText.DrawClone("same");
    newText.SetText(1., 29000, "K^{*}(892)");
    newText.DrawClone("same");
    delete legend;
    legend = new TLegend(0.64, 0.69, 0.84, 0.89);
    legend->SetTextSize(0.025);
    legend->SetHeader("#bf{pp, #sqrt{s} = 510 GeV}", "C");
    legend->AddEntry(MKpiWide->GetName(), "#bf{upcDST data}");
    legend->SetBorderSize(0);
    legend->DrawClone("SAME");
    resultCanvas->SaveAs("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/presentation_pdfs/MKpiWide.pdf");
    resultCanvas->Clear();
    //Kpi with PID
    MKpiWidedEdx->SetTitle("");
    MKpiWidedEdx->SetMinimum(0);
    MKpiWidedEdx->SetLineColor(kBlue+2);
    MKpiWidedEdx->Draw("hist");
    newText.SetText(0.8, 3400, "#splitline{Unidentified}{resonance (~700)}");
    newText.DrawClone("same");
    newText.SetText(1., 5600, "K^{*}(892)");
    newText.DrawClone("same");
    delete legend;
    legend = new TLegend(0.64, 0.69, 0.84, 0.89);
    legend->SetTextSize(0.025);
    legend->SetHeader("#bf{pp, #sqrt{s} = 510 GeV}", "C");
    legend->AddEntry(MKpiWidedEdx->GetName(), "#bf{upcDST data}");
    legend->SetBorderSize(0);
    legend->DrawClone("SAME");
    resultCanvas->SaveAs("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/presentation_pdfs/MKpiWidedEdx.pdf");
    resultCanvas->Clear();
    //pipi without PID
    MpipiWide->SetTitle("");
    MpipiWide->SetMinimum(0);
    MpipiWide->SetLineColor(kBlue+2);
    MpipiWide->Draw("hist");
    newText.SetText(0.6, 14000, "K^{0}_{S}");
    newText.DrawClone("same");
    delete legend;
    legend = new TLegend(0.64, 0.69, 0.84, 0.89);
    legend->SetTextSize(0.025);
    legend->SetHeader("#bf{pp, #sqrt{s} = 510 GeV}", "C");
    legend->AddEntry(MpipiWide->GetName(), "#bf{upcDST data}");
    legend->SetBorderSize(0);
    legend->DrawClone("SAME");
    resultCanvas->SaveAs("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/presentation_pdfs/MpipiWide.pdf");
    resultCanvas->Clear();
    //pipi with PID
    MpipiWidedEdx->SetTitle("");
    MpipiWidedEdx->SetMinimum(0);
    MpipiWidedEdx->SetLineColor(kBlue+2);
    MpipiWidedEdx->Draw("hist");
    newText.SetText(0.6, 13000, "K^{0}_{S}");
    newText.DrawClone("same");
    delete legend;
    legend = new TLegend(0.64, 0.69, 0.84, 0.89);
    legend->SetTextSize(0.025);
    legend->SetHeader("#bf{pp, #sqrt{s} = 510 GeV}", "C");
    legend->AddEntry(MpipiWidedEdx->GetName(), "#bf{upcDST data}");
    legend->SetBorderSize(0);
    legend->DrawClone("SAME");
    resultCanvas->SaveAs("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/presentation_pdfs/MpipiWidedEdx.pdf");
    resultCanvas->Clear();

    return 0;
}

void draw_and_save(TH1D *data, TH1D *sim, std::string fileTitle, double x1, double y1, double x2, double y2){
    TCanvas *resultCanvas = new TCanvas("resultCanvas", "resultCanvas", 1800, 1600);
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
    if(x1*y1*x2*y2!=0){
        TLegend *legend = new TLegend(x1, y1, x2, y2);
        legend->SetTextSize(0.025);
        legend->SetHeader("#bf{pp, #sqrt{s} = 510 GeV}", "C");
        legend->AddEntry(data->GetName(), "#bf{upcDST data}");
        legend->AddEntry(sim->GetName(), "#bf{PYTHIA8 simulation}");
        legend->SetBorderSize(0);
        legend->Draw("SAME");
    }
    resultCanvas->SetLeftMargin(0.15);
    resultCanvas->SetRightMargin(0.07);
    gPad->Update();
    resultCanvas->SaveAs(("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/presentation_pdfs/"+fileTitle+".pdf").c_str());
}

void simpler_draw_and_save(TH1 *hist, std::string filename, bool is3D, bool isYlogarithmic){
    TCanvas *resultCanvas = new TCanvas("resultCanvas", "resultCanvas", 1800, 1600);
    std::string fileTitle;
    gStyle->SetFrameLineWidth(1);
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);
    resultCanvas->SetLogy(isYlogarithmic);
    resultCanvas->SetLeftMargin(0.15);
    if(isYlogarithmic){
        // hist->SetMinimum(0.1);
    } else{
        hist->SetMinimum(0);
    }
    hist->SetLineColor(kBlue+2);
    if(is3D){
        resultCanvas->SetLogz();
        hist->Draw("colz");
    } else{
        hist->Draw("hist");
    }
    if(filename.length()==0){
        fileTitle = hist->GetName();
    } else{
        fileTitle = filename;
    }
    resultCanvas->SaveAs(("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/presentation_pdfs/"+fileTitle+".pdf").c_str());
}

void simpler_draw_and_save(TEfficiency *hist, std::string filename){
    TCanvas *resultCanvas = new TCanvas("resultCanvas", "resultCanvas", 1800, 1600);
    std::string fileTitle;
    gStyle->SetFrameLineWidth(1);
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);
    resultCanvas->SetLeftMargin(0.15);
    hist->SetLineColor(kBlue+2);
    hist->SetMarkerStyle(kFullCircle);
    hist->SetMarkerSize(2);
    hist->SetMarkerColor(kBlue);
    hist->Draw("e");
    if(filename.length()==0){
        fileTitle = hist->GetName();
    } else{
        fileTitle = filename;
    }
    resultCanvas->SaveAs(("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/presentation_pdfs/"+fileTitle+".pdf").c_str());
}
