#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TApplication.h"

#include <iomanip>
#include <sstream>

using namespace std;

void extract_histograms(TH2D &newdata, std::string histName);
void draw_and_save(TH1D &newdata, string histTitle, string fileTitle, double x1 = .65, double y1 = .79, double x2 = .86, double y2 = .89);
void phi_differential_crossection_fit(TH1D &newdataK01D, string histName);

int main(int argc, char *argv[]){
    //something used so that the histograms would draw
    //https://stackoverflow.com/questions/30932725/painting-a-tcanvas-to-the-screen-in-a-compiled-root-cern-application
    TApplication theApp("App", &argc, argv);

    gStyle->SetFrameLineWidth(1);
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);

    //data for master thesis
    //phi
    TH1D newdataphipt;
    phi_differential_crossection_fit(newdataphipt, "Mphipt2DHist");
    draw_and_save(newdataphipt, "Uncorrected #varphi(1020) yields vs. p_{T};p_{T} [GeV];n_{#varphi(1020)}", "Differential_crossection_phi_pt_new_data.pdf");
    TH1D newdataphieta;
    phi_differential_crossection_fit(newdataphieta, "Mphieta2DHist");
    draw_and_save(newdataphieta, "Uncorrected #varphi(1020) yields vs.  #eta;#eta;n_{#varphi(1020)}", "Differential_crossection_phi_eta_new_data.pdf", 0.37, 0.2, 0.63, 0.3);
    theApp.Run();
    return 0;
}

void extract_histograms(TH2D &newdata, string histName){
    TFile *anaoutputnewdata1 = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/AnaOutput_Inclusive_analysis_Kstar_phi_old_data.root");
    newdata = *(TH2D *)anaoutputnewdata1->Get(histName.c_str());
}

void draw_and_save(TH1D &newdata, string histTitle, string fileTitle, double x1, double y1, double x2, double y2){
    TCanvas *resultCanvas = new TCanvas("resultCanvas", "resultCanvas", 1800, 1600);
    //main ones, to start
    newdata.SetMinimum(0);
    newdata.SetTitle(histTitle.c_str());
    //the rest of settings
    newdata.SetMarkerStyle(kFullCircle);
    newdata.SetMarkerColor(kBlue);
    newdata.SetLineColor(kBlue+2);
    newdata.Draw("E");
    TLegend *legend = new TLegend(x1, y1, x2, y2);
    legend->SetTextSize(0.02);
    legend->SetHeader("#bf{pp, #sqrt{s} = 510 GeV}", "C");
    legend->AddEntry(newdata.GetName(), "#bf{new data, new tracks}");
    legend->SetBorderSize(0);
    legend->Draw("SAME");
    resultCanvas->SetLeftMargin(0.15);
    resultCanvas->SetRightMargin(0.05);
    gPad->Update();
    printf("Confirm by pressing \"Enter\"\n");
    getchar();
    resultCanvas->SaveAs(("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/"+fileTitle).c_str());
    fileTitle.replace(fileTitle.find(".pdf"), 4, ".root");
    newdata.SaveAs(("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/"+fileTitle).c_str());
    resultCanvas->Close();
    delete resultCanvas;
}

void phi_differential_crossection_fit(TH1D &newdataK01D, string histName){

    TH2D newdataK0;
    extract_histograms(newdataK0, histName);

    std::vector<double> binEdgesNewData;
    for(size_t i = 0; i<newdataK0.GetNbinsY(); i++)
        binEdgesNewData.push_back(newdataK0.GetYaxis()->GetBinLowEdge(i+1));
    binEdgesNewData.push_back(newdataK0.GetYaxis()->GetBinLowEdge(newdataK0.GetNbinsY())+newdataK0.GetYaxis()->GetBinWidth(newdataK0.GetNbinsY()));

    newdataK01D = TH1D("newdataK01D", "newdataK01D", newdataK0.GetNbinsY(), binEdgesNewData.data());

    for(int i = 0; i<newdataK0.GetNbinsY(); i++){
        TH1D resultKHistNewDataNewTracks = *newdataK0.ProjectionX("K01DHist", i+1, i+1);

        //K0
        //setting up the functions and data
        double lowrangefit = 0.98, upperrangefit = 1.08;
        TCanvas *resultK = new TCanvas("resultK", "resultK", 1600, 800);
        resultK->SetLeftMargin(0.15);
        TF1 *GfitK = new TF1("GfitK", "gausn(0) + pol0(3)", lowrangefit, upperrangefit);
        TF1 *GfitKBcg = new TF1("GfitKBcg", "pol0", lowrangefit, upperrangefit);
        TF1 *GfitK1Sig = new TF1("GfitK1Sig", "gausn", lowrangefit, upperrangefit);
        GfitK->SetParNames("Constant", "Mean", "Sigma", "a");
        GfitK->SetParameters(100, 1.020, 5e-3, 0);
        GfitK->FixParameter(3, 0);
        //old data, the base of drawing
        resultKHistNewDataNewTracks.SetMinimum(0);
        resultKHistNewDataNewTracks.SetMarkerStyle(kFullCircle);
        resultKHistNewDataNewTracks.SetMarkerColor(kBlue);
        resultKHistNewDataNewTracks.SetLineColor(kBlue+2);
        resultKHistNewDataNewTracks.Fit(GfitK, "0BR");
        // TH1D *tempHist = (TH1D *)resultKHistNewDataNewTracks.DrawCopy("E", "NewDataNewTracks");
        resultKHistNewDataNewTracks.DrawCopy("E", "NewDataNewTracks");
        gPad->Update();
        // TPaveStats *st = (TPaveStats *)tempHist->FindObject("stats");
        // st->SetX1NDC(0.15);
        // st->SetX2NDC(0.35);
        // st->SetY1NDC(0.65);
        // st->SetY2NDC(0.85);
        // st->SetBorderSize(0);
        Double_t paramsK[5];
        GfitK->GetParameters(paramsK);
        GfitKBcg->SetParameters(paramsK+3);
        GfitK1Sig->SetParameters(paramsK);
        GfitKBcg->SetLineStyle(3);
        GfitKBcg->SetLineWidth(3);
        GfitK1Sig->SetLineColor(kBlue);
        GfitK->SetNpx(1000);
        TF1 *newDataNewTracksFit = GfitK->DrawCopy("CSAME");
        TF1 *newDataNewTracksBcgFit = GfitKBcg->DrawCopy("CSAME");
        TF1 *newDataNewTracksSigFit = GfitK1Sig->DrawCopy("CSAME");
        newDataNewTracksSigFit->SetParErrors(GfitK->GetParErrors());

        TLegend *legendK = new TLegend(.18, .65, .45, .89);
        legendK->SetTextSize(0.02);
        legendK->AddEntry("K01DHistNewDataNewTracks", "pp new data, new tracks, #sqrt{s} = 510 GeV");
        legendK->SetBorderSize(0);
        legendK->Draw("SAME");
        gPad->Update();

        //interactive part
        std::cout<<"Enter which parameter (0-4) and how much change"<<std::endl;
        std::cout<<"Writing \"abs\" before just sets the parameter"<<std::endl;
        std::cout<<"Or write \"fit\" to fit"<<std::endl;
        std::cout<<"Or write \"zero\" if there is no data suitable to fit"<<std::endl;
        std::cout<<"Or press \"Enter\" if everything ok"<<std::endl;
        string response;
        //needed to read whole line, not just until first space
        std::getline(std::cin, response);
        while(response.size()!=0&&response.find("zero")==string::npos){
            //fitting/changing part
            if(response.find("fit")!=string::npos){
                GfitK->SetParameters(newDataNewTracksFit->GetParameters());
                resultKHistNewDataNewTracks.Fit(GfitK, "0BR");
                printf("Chi2 = %f, ndof = %d\n", GfitK->GetChisquare(), GfitK->GetNDF());
                newDataNewTracksFit->SetParameters(GfitK->GetParameters());
                newDataNewTracksBcgFit->SetParameters(GfitK->GetParameters()+3);
                newDataNewTracksSigFit->SetParameters(GfitK->GetParameters());
                newDataNewTracksSigFit->SetParErrors(GfitK->GetParErrors());
            } else if(response.find("zero")!=string::npos){
                newDataNewTracksSigFit->SetParameter(0, 0);
                newDataNewTracksSigFit->SetParError(0, 0);
            } else{
                int coeff_number;
                double change, relative;
                if(response.find("abs")==string::npos){
                    relative = 1.;
                    sscanf(response.c_str(), "%d%lf", &coeff_number, &change);
                } else{
                    relative = 0.;
                    sscanf(response.c_str(), "abs %d%lf", &coeff_number, &change);
                }
                newDataNewTracksFit->SetParameter(coeff_number, newDataNewTracksFit->GetParameter(coeff_number)*relative+change);
                if(coeff_number>2)
                    newDataNewTracksBcgFit->SetParameter(coeff_number-3, newDataNewTracksBcgFit->GetParameter(coeff_number-3)*relative+change);
                else
                    newDataNewTracksSigFit->SetParameter(coeff_number, newDataNewTracksSigFit->GetParameter(coeff_number)*relative+change);
            }
            //drawing part
            newDataNewTracksFit->Draw("CSAME");
            newDataNewTracksBcgFit->Draw("CSAME");
            newDataNewTracksSigFit->Draw("CSAME");
            gPad->Update();
            //interactive part
            std::cout<<"Enter which parameter (0-4) and how much change"<<std::endl;
            std::cout<<"Writing \"abs\" before just sets the parameter"<<std::endl;
            std::cout<<"Or write \"fit\" to fit"<<std::endl;
            std::cout<<"Or write \"zero\" if there is no data suitable to fit"<<std::endl;
            std::cout<<"Or press \"Enter\" if everything ok"<<std::endl;
            //needed to read whole line, not just until first space
            std::getline(std::cin, response);
        }
        //adding to the result histogram
        if(response.find("zero")!=string::npos){
            newdataK01D.SetBinContent(i+1, 0.);
            newdataK01D.SetBinError(i+1, 0.);
        } else{
            newdataK01D.SetBinContent(i+1, newDataNewTracksSigFit->GetParameter(0)/newdataK0.GetXaxis()->GetBinWidth(0)*newdataK0.GetYaxis()->GetBinWidth(i+1));
            newdataK01D.SetBinError(i+1, newDataNewTracksSigFit->GetParError(0)/newdataK0.GetXaxis()->GetBinWidth(0)*newdataK0.GetYaxis()->GetBinWidth(i+1));
        }
        printf("Accepted current parameters\n");
        //on last loop close the TCanvas
        if(i==newdataK0.GetNbinsY()-1)
            resultK->Close();
    }
}