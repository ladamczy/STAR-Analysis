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
void K0_differential_crossection_fit(TH1D &newdataK01D, string histName);
void Lambda_differential_crossection_fit(TH1D &newdataLambda1D, string histName);

int main(int argc, char *argv[]){
    //something used so that the histograms would draw
    //https://stackoverflow.com/questions/30932725/painting-a-tcanvas-to-the-screen-in-a-compiled-root-cern-application
    TApplication theApp("App", &argc, argv);

    gStyle->SetFrameLineWidth(1);
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);

    //K0
    TH1D newdataK0pt;
    K0_differential_crossection_fit(newdataK0pt, "K0pt2DHist");
    draw_and_save(newdataK0pt, "Uncorrected K^{ 0}_{ S} yields vs. p_{T};p_{T} [GeV];n_{K^{0}_{S}}", "Differential_crossection_K0_pt_new_data.pdf");
    TH1D newdataK0eta;
    K0_differential_crossection_fit(newdataK0eta, "K0eta2DHist");
    draw_and_save(newdataK0eta, "Uncorrected K^{ 0}_{ S} yields vs.  #eta;#eta;n_{K^{0}_{S}}", "Differential_crossection_K0_eta_new_data.pdf", 0.37, 0.2, 0.63, 0.3);
    TH1D newdataK0XiMulti;
    K0_differential_crossection_fit(newdataK0XiMulti, "K0XiMulti2DHist");
    draw_and_save(newdataK0XiMulti, "Uncorrected K^{ 0}_{ S} yields vs. log(#xi_{E}#cdot#xi_{W});log(#xi_{E}#cdot#xi_{W});n_{K^{0}_{S}}", "Differential_crossection_K0_XiMulti_new_data.pdf");
    TH1D newdataK0XiSum;
    K0_differential_crossection_fit(newdataK0XiSum, "K0XiSum2DHist");
    draw_and_save(newdataK0XiSum, "Uncorrected K^{ 0}_{ S} yields vs. #xi_{E}+#xi_{W};#xi_{E}+#xi_{W};n_{K^{0}_{S}}", "Differential_crossection_K0_XiSum_new_data.pdf");
    //Lambda
    TH1D newdataLambdapt;
    Lambda_differential_crossection_fit(newdataLambdapt, "Lambdapt2DHist");
    draw_and_save(newdataLambdapt, "Uncorrected #Lambda^{ 0} yields vs. p_{T};p_{T} [GeV];n_{#Lambda^{0}}", "Differential_crossection_Lambda_pt_new_data.pdf");
    TH1D newdataLambdaeta;
    Lambda_differential_crossection_fit(newdataLambdaeta, "Lambdaeta2DHist");
    draw_and_save(newdataLambdaeta, "Uncorrected #Lambda^{ 0} yields vs.  #eta;#eta;n_{#Lambda^{0}}", "Differential_crossection_Lambda_eta_new_data.pdf", 0.37, 0.2, 0.63, 0.3);
    TH1D newdataLambdaXiMulti;
    Lambda_differential_crossection_fit(newdataLambdaXiMulti, "LambdaXiMulti2DHist");
    draw_and_save(newdataLambdaXiMulti, "Uncorrected #Lambda^{ 0} yields vs. log(#xi_{E}#cdot#xi_{W});log(#xi_{E}#cdot#xi_{W});n_{#Lambda^{0}}", "Differential_crossection_Lambda_XiMulti_new_data.pdf");
    TH1D newdataLambdaXiSum;
    Lambda_differential_crossection_fit(newdataLambdaXiSum, "LambdaXiSum2DHist");
    draw_and_save(newdataLambdaXiSum, "Uncorrected #Lambda^{ 0} yields vs. #xi_{E}+#xi_{W};#xi_{E}+#xi_{W};n_{#Lambda^{0}}", "Differential_crossection_Lambda_XiSum_new_data.pdf");

    theApp.Run();
    return 0;
}

void extract_histograms(TH2D &newdata, string histName){
    TFile *anaoutputnewdata1 = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/CPT/AnaOutput_Inclusive_analysis_with_STUPCV0_with_extended_range_noAfterburner_new_data_new_tracks.root");
    TFile *anaoutputnewdata2 = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/CPTnoBBCL/AnaOutput_Inclusive_analysis_with_STUPCV0_with_extended_range_noAfterburner_new_data_new_tracks.root");

    newdata = *(TH2D *)anaoutputnewdata1->Get(histName.c_str());
    newdata.Add((TH2D *)anaoutputnewdata2->Get(histName.c_str()));
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
    resultCanvas->Close();
    delete resultCanvas;
}

void K0_differential_crossection_fit(TH1D &newdataK01D, string histName){

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
        double lowrangefit = 0.48, upperrangefit = 0.52;
        TCanvas *resultK = new TCanvas("resultK", "resultK", 1600, 800);
        resultK->SetLeftMargin(0.15);
        TF1 *GfitK = new TF1("GfitK", "gausn(0) + pol1(3)", lowrangefit, upperrangefit);
        TF1 *GfitKBcg = new TF1("GfitKBcg", "pol1", lowrangefit, upperrangefit);
        TF1 *GfitK1Sig = new TF1("GfitK1Sig", "gausn", lowrangefit, upperrangefit);
        GfitK->SetParNames("Constant", "Mean", "Sigma", "b", "a");
        GfitK->SetParameters(100, 0.497, 5e-3);
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

void Lambda_differential_crossection_fit(TH1D &newdataLambda1D, string histName){

    TH2D newdataLambda;
    extract_histograms(newdataLambda, histName);

    std::vector<double> binEdgesNewData;
    for(size_t i = 0; i<newdataLambda.GetNbinsY(); i++)
        binEdgesNewData.push_back(newdataLambda.GetYaxis()->GetBinLowEdge(i+1));
    binEdgesNewData.push_back(newdataLambda.GetYaxis()->GetBinLowEdge(newdataLambda.GetNbinsY())+newdataLambda.GetYaxis()->GetBinWidth(newdataLambda.GetNbinsY()));

    newdataLambda1D = TH1D("newdataLambda1D", "newdataLambda1D", newdataLambda.GetNbinsY(), binEdgesNewData.data());

    for(int i = 0; i<newdataLambda.GetNbinsY(); i++){
        TH1D resultLambdaHistNewDataNewTracks = *newdataLambda.ProjectionX("Lambda1DHist", i+1, i+1);

        //Lambda
        //setting up the functions and data
        double fitmin = 1.1;
        double fitmax = 1.14;
        TCanvas *resultLambda = new TCanvas("resultLambda", "resultLambda", 1600, 800);
        resultLambda->SetLeftMargin(0.15);
        TF1 *GfitLambda = new TF1("GfitLambda", "gausn(0) + pol2(3)", fitmin, fitmax);
        TF1 *GfitLambdaBcg = new TF1("GfitLambdaBcg", "pol2", fitmin, fitmax);
        TF1 *GfitLambda1Sig = new TF1("GfitLambda1Sig", "gausn", fitmin, fitmax);
        GfitLambda->SetParNames("Constant", "Mean", "Sigma", "c", "b", "a");
        // GfitLambda->SetParameters(100, 1.115, 5e-3, 2500000, -4200000, 1200000);
        // GfitLambda->SetParameters(3, 1.115, 5e-3, 900000, -1600000, 780000);
        GfitLambda->SetParameters(3, 1.115, 5e-3, 100, 0, 0);
        GfitLambda->SetParLimits(0, 0, 200);
        GfitLambda->SetParLimits(1, 1.11, 1.12);
        GfitLambda->SetParLimits(2, 5e-5, 0.01);
        //old data, the base of drawing
        resultLambdaHistNewDataNewTracks.SetMinimum(0);
        resultLambdaHistNewDataNewTracks.SetMarkerStyle(kFullCircle);
        resultLambdaHistNewDataNewTracks.SetMarkerColor(kBlue);
        resultLambdaHistNewDataNewTracks.SetLineColor(kBlue+2);
        resultLambdaHistNewDataNewTracks.Fit(GfitLambda, "0BR");
        // TH1D *tempHist = (TH1D *)resultLambdaHistNewDataNewTracks.DrawCopy("E", "NewDataNewTracks");
        resultLambdaHistNewDataNewTracks.DrawCopy("E", "NewDataNewTracks");
        gPad->Update();
        // TPaveStats *st = (TPaveStats *)tempHist->FindObject("stats");
        // st->SetX1NDC(0.15);
        // st->SetX2NDC(0.35);
        // st->SetY1NDC(0.65);
        // st->SetY2NDC(0.85);
        // st->SetBorderSize(0);
        Double_t paramsLambda[6];
        GfitLambda->GetParameters(paramsLambda);
        // cout<<"PARAMETERS:"<<endl;
        // for(int i = 0; i<5; i++){
        //     cout<<paramsLambda[i]<<endl;
        // }
        GfitLambdaBcg->SetParameters(paramsLambda+3);
        GfitLambda1Sig->SetParameters(paramsLambda);
        GfitLambdaBcg->SetLineStyle(3);
        GfitLambdaBcg->SetLineWidth(3);
        GfitLambda1Sig->SetLineColor(kBlue);
        GfitLambda->SetNpx(1000);
        TF1 *newDataNewTracksFit = GfitLambda->DrawCopy("CSAME");
        TF1 *newDataNewTracksBcgFit = GfitLambdaBcg->DrawCopy("CSAME");
        TF1 *newDataNewTracksSigFit = GfitLambda1Sig->DrawCopy("CSAME");
        newDataNewTracksSigFit->SetParErrors(GfitLambda->GetParErrors());

        TLegend *legendLambda = new TLegend(.18, .65, .45, .89);
        legendLambda->SetTextSize(0.02);
        legendLambda->AddEntry("Lambda1DHistNewDataNewTracks", "pp new data, new tracks, #sqrt{s} = 510 GeV");
        legendLambda->SetBorderSize(0);
        legendLambda->Draw("SAME");
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
                GfitLambda->SetParameters(newDataNewTracksFit->GetParameters());
                resultLambdaHistNewDataNewTracks.Fit(GfitLambda, "0BR");
                newDataNewTracksFit->SetParameters(GfitLambda->GetParameters());
                newDataNewTracksBcgFit->SetParameters(GfitLambda->GetParameters()+3);
                newDataNewTracksSigFit->SetParameters(GfitLambda->GetParameters());
                newDataNewTracksSigFit->SetParErrors(GfitLambda->GetParErrors());
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
            newdataLambda1D.SetBinContent(i+1, 0.);
            newdataLambda1D.SetBinError(i+1, 0.);
        } else{
            newdataLambda1D.SetBinContent(i+1, newDataNewTracksSigFit->GetParameter(0)/newdataLambda.GetXaxis()->GetBinWidth(0)*newdataLambda.GetYaxis()->GetBinWidth(i+1));
            newdataLambda1D.SetBinError(i+1, newDataNewTracksSigFit->GetParError(0)/newdataLambda.GetXaxis()->GetBinWidth(0)*newdataLambda.GetYaxis()->GetBinWidth(i+1));
        }
        printf("Accepted current parameters\n");
        if(i==newdataLambda.GetNbinsY()-1)
            resultLambda->Close();
    }
}