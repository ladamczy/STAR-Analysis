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

void extract_histograms(TH2D &olddata, TH2D &oldtracks, TH2D &newdata, std::string histName);
void draw_and_save(TH1D &olddata, TH1D &oldtracks, TH1D &newdata, string histTitle, string fileTitle);
void K0_differential_crossection_fit(TH1D &olddataK01D, TH1D &oldtracksK01D, TH1D &newdataK01D, string histName);
void Lambda_differential_crossection_fit(TH1D &olddataLambda1D, TH1D &oldtracksLambda1D, TH1D &newdataLambda1D, string histName);

int main(int argc, char *argv[]){
    //something used so that the histograms would draw
    //https://stackoverflow.com/questions/30932725/painting-a-tcanvas-to-the-screen-in-a-compiled-root-cern-application
    TApplication theApp("App", &argc, argv);

    gStyle->SetFrameLineWidth(1);
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);

    //K0
    TH1D olddataK0pt, oldtracksK0pt, newdataK0pt;
    K0_differential_crossection_fit(olddataK0pt, oldtracksK0pt, newdataK0pt, "K0pt2DHist");
    draw_and_save(olddataK0pt, oldtracksK0pt, newdataK0pt, "Differential crossection for K^{ 0}_{ S} with respect to p_{T};p_{T} [GeV];n_{K^{0}_{S}}", "Differential_crossection_K0_pt.pdf");
    TH1D olddataK0eta, oldtracksK0eta, newdataK0eta;
    K0_differential_crossection_fit(olddataK0eta, oldtracksK0eta, newdataK0eta, "K0eta2DHist");
    draw_and_save(olddataK0eta, oldtracksK0eta, newdataK0eta, "Differential crossection for K^{ 0}_{ S} with respect to #eta;#eta;n_{K^{0}_{S}}", "Differential_crossection_K0_eta.pdf");
    TH1D olddataK0XiMulti, oldtracksK0XiMulti, newdataK0XiMulti;
    K0_differential_crossection_fit(olddataK0XiMulti, oldtracksK0XiMulti, newdataK0XiMulti, "K0XiMulti2DHist");
    draw_and_save(olddataK0XiMulti, oldtracksK0XiMulti, newdataK0XiMulti, "Differential crossection for K^{ 0}_{ S} with respect to log(#xi_{E}#cdot#xi_{W});log(#xi_{E}#cdot#xi_{W});n_{K^{0}_{S}}", "Differential_crossection_K0_XiMulti.pdf");
    TH1D olddataK0XiSum, oldtracksK0XiSum, newdataK0XiSum;
    K0_differential_crossection_fit(olddataK0XiSum, oldtracksK0XiSum, newdataK0XiSum, "K0XiSum2DHist");
    draw_and_save(olddataK0XiSum, oldtracksK0XiSum, newdataK0XiSum, "Differential crossection for K^{ 0}_{ S} with respect to #xi_{E}+#xi_{W};#xi_{E}+#xi_{W};n_{K^{0}_{S}}", "Differential_crossection_K0_XiSum.pdf");
    //Lambda
    TH1D olddataLambdapt, oldtracksLambdapt, newdataLambdapt;
    Lambda_differential_crossection_fit(olddataLambdapt, oldtracksLambdapt, newdataLambdapt, "Lambdapt2DHist");
    draw_and_save(olddataLambdapt, oldtracksLambdapt, newdataLambdapt, "Differential crossection for #Lambda^{ 0} with respect to p_{T};p_{T} [GeV];n_{#Lambda^{0}}", "Differential_crossection_Lambda_pt.pdf");
    TH1D olddataLambdaeta, oldtracksLambdaeta, newdataLambdaeta;
    Lambda_differential_crossection_fit(olddataLambdaeta, oldtracksLambdaeta, newdataLambdaeta, "Lambdaeta2DHist");
    draw_and_save(olddataLambdaeta, oldtracksLambdaeta, newdataLambdaeta, "Differential crossection for #Lambda^{ 0} with respect to #eta;#eta;n_{#Lambda^{0}}", "Differential_crossection_Lambda_eta.pdf");
    TH1D olddataLambdaXiMulti, oldtracksLambdaXiMulti, newdataLambdaXiMulti;
    Lambda_differential_crossection_fit(olddataLambdaXiMulti, oldtracksLambdaXiMulti, newdataLambdaXiMulti, "LambdaXiMulti2DHist");
    draw_and_save(olddataLambdaXiMulti, oldtracksLambdaXiMulti, newdataLambdaXiMulti, "Differential crossection for #Lambda^{ 0} with respect to log(#xi_{E}#cdot#xi_{W});log(#xi_{E}#cdot#xi_{W});n_{#Lambda^{0}}", "Differential_crossection_Lambda_XiMulti.pdf");
    TH1D olddataLambdaXiSum, oldtracksLambdaXiSum, newdataLambdaXiSum;
    Lambda_differential_crossection_fit(olddataLambdaXiSum, oldtracksLambdaXiSum, newdataLambdaXiSum, "LambdaXiSum2DHist");
    draw_and_save(olddataLambdaXiSum, oldtracksLambdaXiSum, newdataLambdaXiSum, "Differential crossection for #Lambda^{ 0} with respect to #xi_{E}+#xi_{W};#xi_{E}+#xi_{W};n_{#Lambda^{0}}", "Differential_crossection_Lambda_XiSum.pdf");

    theApp.Run();
    return 0;
}

void extract_histograms(TH2D &olddata, TH2D &oldtracks, TH2D &newdata, string histName){
    TFile *anaoutputolddata1 = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/CPT/AnaOutput_Inclusive_analysis_with_STUPCV0_with_extended_range_noAfterburner.root");
    TFile *anaoutputolddata2 = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/CPTnoBBCL/AnaOutput_Inclusive_analysis_with_STUPCV0_with_extended_range_noAfterburner.root");
    TFile *anaoutputoldtracks1 = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/CPT/AnaOutput_Inclusive_analysis_with_STUPCV0_with_extended_range_noAfterburner_new_data_old_tracks.root");
    TFile *anaoutputoldtracks2 = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/CPTnoBBCL/AnaOutput_Inclusive_analysis_with_STUPCV0_with_extended_range_noAfterburner_new_data_old_tracks.root");
    TFile *anaoutputnewdata1 = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/CPT/AnaOutput_Inclusive_analysis_with_STUPCV0_with_extended_range_noAfterburner_new_data_new_tracks.root");
    TFile *anaoutputnewdata2 = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/CPTnoBBCL/AnaOutput_Inclusive_analysis_with_STUPCV0_with_extended_range_noAfterburner_new_data_new_tracks.root");

    olddata = *(TH2D *)anaoutputolddata1->Get(histName.c_str());
    olddata.Add((TH2D *)anaoutputolddata2->Get(histName.c_str()));
    oldtracks = *(TH2D *)anaoutputoldtracks1->Get(histName.c_str());
    oldtracks.Add((TH2D *)anaoutputoldtracks2->Get(histName.c_str()));
    newdata = *(TH2D *)anaoutputnewdata1->Get(histName.c_str());
    newdata.Add((TH2D *)anaoutputnewdata2->Get(histName.c_str()));
}

void draw_and_save(TH1D &olddata, TH1D &oldtracks, TH1D &newdata, string histTitle, string fileTitle){
    TCanvas *resultCanvas = new TCanvas("resultCanvas", "resultCanvas", 1800, 1600);
    //main ones, to start
    olddata.SetMinimum(0);
    olddata.SetMaximum(max(olddata.GetMaximum(), max(oldtracks.GetMaximum(), newdata.GetMaximum()))*1.1);
    olddata.SetTitle(histTitle.c_str());
    //the rest of settings
    olddata.SetMarkerStyle(kFullCircle);
    olddata.SetMarkerColor(kRed);
    olddata.SetLineColor(kRed+2);
    olddata.Draw("E");
    oldtracks.SetMarkerStyle(kFullCircle);
    oldtracks.SetMarkerColor(kGreen);
    oldtracks.SetLineColor(kGreen+2);
    oldtracks.Draw("ESAME");
    newdata.SetMarkerStyle(kFullCircle);
    newdata.SetMarkerColor(kBlue);
    newdata.SetLineColor(kBlue+2);
    newdata.Draw("ESAME");
    TLegend *legend = new TLegend(.65, .65, .86, .89);
    legend->SetTextSize(0.02);
    legend->SetHeader("#bf{pp, #sqrt{s} = 510 GeV}", "C");
    legend->AddEntry(olddata.GetName(), "#bf{old data}");
    legend->AddEntry(oldtracks.GetName(), "#bf{new data, old tracks}");
    legend->AddEntry(newdata.GetName(), "#bf{new data, new tracks}");
    legend->SetBorderSize(0);
    legend->Draw("SAME");
    gPad->Update();
    printf("Confirm by pressing \"Enter\"\n");
    getchar();
    resultCanvas->SaveAs(("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/"+fileTitle).c_str());
    resultCanvas->Close();
    delete resultCanvas;
}

void K0_differential_crossection_fit(TH1D &olddataK01D, TH1D &oldtracksK01D, TH1D &newdataK01D, string histName){

    TH2D olddataK0, oldtracksK0, newdataK0;
    extract_histograms(olddataK0, oldtracksK0, newdataK0, histName);

    std::vector<double> binEdgesOldData, binEdgesOldTracks, binEdgesNewData;
    for(size_t i = 0; i<olddataK0.GetNbinsY(); i++)
        binEdgesOldData.push_back(olddataK0.GetYaxis()->GetBinLowEdge(i+1));
    binEdgesOldData.push_back(olddataK0.GetYaxis()->GetBinLowEdge(olddataK0.GetNbinsY())+olddataK0.GetYaxis()->GetBinWidth(olddataK0.GetNbinsY()));
    for(size_t i = 0; i<oldtracksK0.GetNbinsY(); i++)
        binEdgesOldTracks.push_back(oldtracksK0.GetYaxis()->GetBinLowEdge(i+1));
    binEdgesOldTracks.push_back(oldtracksK0.GetYaxis()->GetBinLowEdge(oldtracksK0.GetNbinsY())+oldtracksK0.GetYaxis()->GetBinWidth(oldtracksK0.GetNbinsY()));
    for(size_t i = 0; i<newdataK0.GetNbinsY(); i++)
        binEdgesNewData.push_back(newdataK0.GetYaxis()->GetBinLowEdge(i+1));
    binEdgesNewData.push_back(newdataK0.GetYaxis()->GetBinLowEdge(newdataK0.GetNbinsY())+newdataK0.GetYaxis()->GetBinWidth(newdataK0.GetNbinsY()));

    olddataK01D = TH1D("olddataK01D", "olddataK01D", olddataK0.GetNbinsY(), binEdgesOldData.data());
    oldtracksK01D = TH1D("oldtracksK01D", "oldtracksK01D", oldtracksK0.GetNbinsY(), binEdgesOldTracks.data());
    newdataK01D = TH1D("newdataK01D", "newdataK01D", newdataK0.GetNbinsY(), binEdgesNewData.data());

    for(int i = 0; i<olddataK0.GetNbinsY(); i++){
        TH1D resultKHistOldData = *olddataK0.ProjectionX("K01DHist", i+1, i+1);
        TH1D resultKHistNewDataOldTracks = *oldtracksK0.ProjectionX("K01DHist", i+1, i+1);
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
        resultKHistOldData.SetMinimum(0);
        resultKHistOldData.SetMaximum(max(resultKHistOldData.GetMaximum(), max(resultKHistNewDataOldTracks.GetMaximum(), resultKHistNewDataNewTracks.GetMaximum()))*1.1);
        resultKHistOldData.SetMarkerStyle(kFullCircle);
        resultKHistOldData.SetMarkerColor(kRed);
        resultKHistOldData.SetLineColor(kRed+2);
        resultKHistOldData.Fit(GfitK, "0BR");
        // TH1D *tempHist = (TH1D *)resultKHistOldData.DrawCopy("E", "OldData");
        resultKHistOldData.DrawCopy("E", "OldData");
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
        GfitK1Sig->SetLineColor(kRed);
        GfitK->SetNpx(1000);
        TF1 *oldDataFit = GfitK->DrawCopy("CSAME");
        TF1 *oldDataBcgFit = GfitKBcg->DrawCopy("CSAME");
        TF1 *oldDataSigFit = GfitK1Sig->DrawCopy("CSAME");
        oldDataSigFit->SetParErrors(GfitK->GetParErrors());
        //new data old tracks
        resultKHistNewDataOldTracks.SetMinimum(0);
        resultKHistNewDataOldTracks.SetMarkerStyle(kFullCircle);
        resultKHistNewDataOldTracks.SetMarkerColor(kGreen);
        resultKHistNewDataOldTracks.SetLineColor(kGreen+2);
        resultKHistNewDataOldTracks.Fit(GfitK, "0BR");
        resultKHistNewDataOldTracks.DrawCopy("E SAME", "NewDataOldTracks");
        gPad->Update();
        GfitK->GetParameters(paramsK);
        GfitKBcg->SetParameters(paramsK+3);
        GfitK1Sig->SetParameters(paramsK);
        GfitKBcg->SetLineStyle(3);
        GfitKBcg->SetLineWidth(3);
        GfitK->SetLineColor(kGreen);
        GfitKBcg->SetLineColor(kGreen);
        GfitK1Sig->SetLineColor(kGreen);
        GfitK->SetNpx(1000);
        TF1 *newDataOldTracksFit = GfitK->DrawCopy("CSAME");
        TF1 *newDataOldTracksBcgFit = GfitKBcg->DrawCopy("CSAME");
        TF1 *newDataOldTracksSigFit = GfitK1Sig->DrawCopy("CSAME");
        newDataOldTracksSigFit->SetParErrors(GfitK->GetParErrors());
        //new data new tracks
        resultKHistNewDataNewTracks.SetMinimum(0);
        resultKHistNewDataNewTracks.SetMarkerStyle(kFullCircle);
        resultKHistNewDataNewTracks.SetMarkerColor(kBlue);
        resultKHistNewDataNewTracks.Fit(GfitK, "0BR");
        resultKHistNewDataNewTracks.DrawCopy("E SAME", "NewDataNewTracks");
        gPad->Update();
        GfitK->GetParameters(paramsK);
        GfitKBcg->SetParameters(paramsK+3);
        GfitK1Sig->SetParameters(paramsK);
        GfitKBcg->SetLineStyle(3);
        GfitKBcg->SetLineWidth(3);
        GfitK->SetLineColor(kBlue);
        GfitKBcg->SetLineColor(kBlue);
        GfitK1Sig->SetLineColor(kBlue);
        GfitK->SetNpx(1000);
        TF1 *newDataNewTracksFit = GfitK->DrawCopy("CSAME");
        TF1 *newDataNewTracksBcgFit = GfitKBcg->DrawCopy("CSAME");
        TF1 *newDataNewTracksSigFit = GfitK1Sig->DrawCopy("CSAME");
        newDataNewTracksSigFit->SetParErrors(GfitK->GetParErrors());

        TLegend *legendK = new TLegend(.18, .65, .45, .89);
        legendK->SetTextSize(0.02);
        legendK->AddEntry("K01DHistOldData", "pp old data, #sqrt{s} = 510 GeV");
        legendK->AddEntry("K01DHistNewDataOldTracks", "pp new data, old tracks, #sqrt{s} = 510 GeV");
        legendK->AddEntry("K01DHistNewDataNewTracks", "pp new data, new tracks, #sqrt{s} = 510 GeV");
        legendK->SetBorderSize(0);
        legendK->Draw("SAME");
        gPad->Update();

        //interactive part
        std::cout<<"Enter which data (0,1,2), which parameter (0-4) and how much change"<<std::endl;
        std::cout<<"Writing \"abs\" before just sets the parameter"<<std::endl;
        std::cout<<"Or write \"fit <data_number>\" to fit (3 for fitting all)"<<std::endl;
        std::cout<<"Or write \"copy <data_number_from> <data_number_to>\" to copy"<<std::endl;
        std::cout<<"Or write \"zero\" if there is no data suitable to fit"<<std::endl;
        std::cout<<"Or press \"Enter\" if everything ok"<<std::endl;
        string response;
        //needed to read whole line, not just until first space
        std::getline(std::cin, response);
        while(response.size()!=0&&response.find("zero")==string::npos){
            //fitting/changing part
            if(response.find("fit")!=string::npos){
                switch(response[4]-'0'){ //presents response as int by substracting ASCII code for 0
                case 0:
                    GfitK->SetParameters(oldDataFit->GetParameters());
                    resultKHistOldData.Fit(GfitK, "0BR");
                    oldDataFit->SetParameters(GfitK->GetParameters());
                    oldDataBcgFit->SetParameters(GfitK->GetParameters()+3);
                    oldDataSigFit->SetParameters(GfitK->GetParameters());
                    oldDataSigFit->SetParErrors(GfitK->GetParErrors());
                    break;
                case 1:
                    GfitK->SetParameters(newDataOldTracksFit->GetParameters());
                    resultKHistNewDataOldTracks.Fit(GfitK, "0BR");
                    newDataOldTracksFit->SetParameters(GfitK->GetParameters());
                    newDataOldTracksBcgFit->SetParameters(GfitK->GetParameters()+3);
                    newDataOldTracksSigFit->SetParameters(GfitK->GetParameters());
                    newDataOldTracksSigFit->SetParErrors(GfitK->GetParErrors());
                    break;
                case 2:
                    GfitK->SetParameters(newDataNewTracksFit->GetParameters());
                    resultKHistNewDataNewTracks.Fit(GfitK, "0BR");
                    newDataNewTracksFit->SetParameters(GfitK->GetParameters());
                    newDataNewTracksBcgFit->SetParameters(GfitK->GetParameters()+3);
                    newDataNewTracksSigFit->SetParameters(GfitK->GetParameters());
                    newDataNewTracksSigFit->SetParErrors(GfitK->GetParErrors());
                    break;
                case 3:
                    GfitK->SetParameters(oldDataFit->GetParameters());
                    resultKHistOldData.Fit(GfitK, "0BR");
                    oldDataFit->SetParameters(GfitK->GetParameters());
                    oldDataBcgFit->SetParameters(GfitK->GetParameters()+3);
                    oldDataSigFit->SetParameters(GfitK->GetParameters());
                    oldDataSigFit->SetParErrors(GfitK->GetParErrors());
                    GfitK->SetParameters(newDataOldTracksFit->GetParameters());
                    resultKHistNewDataOldTracks.Fit(GfitK, "0BR");
                    newDataOldTracksFit->SetParameters(GfitK->GetParameters());
                    newDataOldTracksBcgFit->SetParameters(GfitK->GetParameters()+3);
                    newDataOldTracksSigFit->SetParameters(GfitK->GetParameters());
                    newDataOldTracksSigFit->SetParErrors(GfitK->GetParErrors());
                    GfitK->SetParameters(newDataNewTracksFit->GetParameters());
                    resultKHistNewDataNewTracks.Fit(GfitK, "0BR");
                    newDataNewTracksFit->SetParameters(GfitK->GetParameters());
                    newDataNewTracksBcgFit->SetParameters(GfitK->GetParameters()+3);
                    newDataNewTracksSigFit->SetParameters(GfitK->GetParameters());
                    newDataNewTracksSigFit->SetParErrors(GfitK->GetParErrors());
                    break;
                default:
                    std::cout<<"Whoops, something went wrong with fit"<<std::endl;
                    break;
                }
            } else if(response.find("copy")!=string::npos){
                int copy_from, copy_to;
                TF1 *fit_tab[] = { oldDataFit, newDataOldTracksFit, newDataNewTracksFit };
                sscanf(response.c_str(), "copy %d%d", &copy_from, &copy_to);
                GfitK->SetParameters(fit_tab[copy_from]->GetParameters());
                fit_tab[copy_to]->SetParameters(GfitK->GetParameters());
            } else if(response.find("zero")!=string::npos){
                oldDataSigFit->SetParameter(0, 0);
                oldDataSigFit->SetParError(0, 0);
                newDataOldTracksSigFit->SetParameter(0, 0);
                newDataOldTracksSigFit->SetParError(0, 0);
                newDataNewTracksSigFit->SetParameter(0, 0);
                newDataNewTracksSigFit->SetParError(0, 0);
            } else{
                int data_number, coeff_number;
                double change, relative;
                if(response.find("abs")==string::npos){
                    relative = 1.;
                    sscanf(response.c_str(), "%d%d%lf", &data_number, &coeff_number, &change);
                } else{
                    relative = 0.;
                    sscanf(response.c_str(), "abs %d%d%lf", &data_number, &coeff_number, &change);
                }
                switch(data_number){
                case 0:
                    oldDataFit->SetParameter(coeff_number, oldDataFit->GetParameter(coeff_number)*relative+change);
                    if(coeff_number>2)
                        oldDataBcgFit->SetParameter(coeff_number-3, oldDataBcgFit->GetParameter(coeff_number-3)*relative+change);
                    else
                        oldDataSigFit->SetParameter(coeff_number, oldDataSigFit->GetParameter(coeff_number)*relative+change);
                    break;
                case 1:
                    newDataOldTracksFit->SetParameter(coeff_number, newDataOldTracksFit->GetParameter(coeff_number)*relative+change);
                    if(coeff_number>2)
                        newDataOldTracksBcgFit->SetParameter(coeff_number-3, newDataOldTracksBcgFit->GetParameter(coeff_number-3)*relative+change);
                    else
                        newDataOldTracksSigFit->SetParameter(coeff_number, newDataOldTracksSigFit->GetParameter(coeff_number)*relative+change);
                    break;
                case 2:
                    newDataNewTracksFit->SetParameter(coeff_number, newDataNewTracksFit->GetParameter(coeff_number)*relative+change);
                    if(coeff_number>2)
                        newDataNewTracksBcgFit->SetParameter(coeff_number-3, newDataNewTracksBcgFit->GetParameter(coeff_number-3)*relative+change);
                    else
                        newDataNewTracksSigFit->SetParameter(coeff_number, newDataNewTracksSigFit->GetParameter(coeff_number)*relative+change);
                    break;
                default:
                    std::cout<<"Whoops, something went wrong with param change"<<std::endl;
                    break;
                }
            }
            //drawing part
            oldDataFit->Draw("CSAME");
            oldDataBcgFit->Draw("CSAME");
            oldDataSigFit->Draw("CSAME");
            newDataOldTracksFit->Draw("CSAME");
            newDataOldTracksBcgFit->Draw("CSAME");
            newDataOldTracksSigFit->Draw("CSAME");
            newDataNewTracksFit->Draw("CSAME");
            newDataNewTracksBcgFit->Draw("CSAME");
            newDataNewTracksSigFit->Draw("CSAME");
            gPad->Update();
            //interactive part
            std::cout<<"Enter which data (0,1,2), which parameter (0-4) and how much change"<<std::endl;
            std::cout<<"Writing \"abs\" before just sets the parameter"<<std::endl;
            std::cout<<"Or write \"fit <data_number>\" to fit (3 for fitting all)"<<std::endl;
            std::cout<<"Or write \"copy <data_number_from> <data_number_to>\" to copy"<<std::endl;
            std::cout<<"Or write \"zero\" if there is no data suitable to fit"<<std::endl;
            std::cout<<"Or press \"Enter\" if everything ok"<<std::endl;
            //needed to read whole line, not just until first space
            std::getline(std::cin, response);
        }
        //adding to the result histogram
        if(response.find("zero")!=string::npos){
            olddataK01D.SetBinContent(i+1, 0.);
            olddataK01D.SetBinError(i+1, 0.);
            oldtracksK01D.SetBinContent(i+1, 0.);
            oldtracksK01D.SetBinError(i+1, 0.);
            newdataK01D.SetBinContent(i+1, 0.);
            newdataK01D.SetBinError(i+1, 0.);
        } else{
            olddataK01D.SetBinContent(i+1, oldDataSigFit->GetParameter(0)/olddataK0.GetXaxis()->GetBinWidth(0)*olddataK0.GetYaxis()->GetBinWidth(i+1));
            olddataK01D.SetBinError(i+1, oldDataSigFit->GetParError(0)/olddataK0.GetXaxis()->GetBinWidth(0)*olddataK0.GetYaxis()->GetBinWidth(i+1));
            oldtracksK01D.SetBinContent(i+1, newDataOldTracksSigFit->GetParameter(0)/oldtracksK0.GetXaxis()->GetBinWidth(0)*oldtracksK0.GetYaxis()->GetBinWidth(i+1));
            oldtracksK01D.SetBinError(i+1, newDataOldTracksSigFit->GetParError(0)/oldtracksK0.GetXaxis()->GetBinWidth(0)*oldtracksK0.GetYaxis()->GetBinWidth(i+1));
            newdataK01D.SetBinContent(i+1, newDataNewTracksSigFit->GetParameter(0)/newdataK0.GetXaxis()->GetBinWidth(0)*newdataK0.GetYaxis()->GetBinWidth(i+1));
            newdataK01D.SetBinError(i+1, newDataNewTracksSigFit->GetParError(0)/newdataK0.GetXaxis()->GetBinWidth(0)*newdataK0.GetYaxis()->GetBinWidth(i+1));
        }
        printf("Accepted current parameters\n");
        //on last loop close the TCanvas
        if(i==olddataK0.GetNbinsY()-1)
            resultK->Close();
    }
}

void Lambda_differential_crossection_fit(TH1D &olddataLambda1D, TH1D &oldtracksLambda1D, TH1D &newdataLambda1D, string histName){

    TH2D olddataLambda, oldtracksLambda, newdataLambda;
    extract_histograms(olddataLambda, oldtracksLambda, newdataLambda, histName);

    std::vector<double> binEdgesOldData, binEdgesOldTracks, binEdgesNewData;
    for(size_t i = 0; i<olddataLambda.GetNbinsY(); i++)
        binEdgesOldData.push_back(olddataLambda.GetYaxis()->GetBinLowEdge(i+1));
    binEdgesOldData.push_back(olddataLambda.GetYaxis()->GetBinLowEdge(olddataLambda.GetNbinsY())+olddataLambda.GetYaxis()->GetBinWidth(olddataLambda.GetNbinsY()));
    for(size_t i = 0; i<oldtracksLambda.GetNbinsY(); i++)
        binEdgesOldTracks.push_back(oldtracksLambda.GetYaxis()->GetBinLowEdge(i+1));
    binEdgesOldTracks.push_back(oldtracksLambda.GetYaxis()->GetBinLowEdge(oldtracksLambda.GetNbinsY())+oldtracksLambda.GetYaxis()->GetBinWidth(oldtracksLambda.GetNbinsY()));
    for(size_t i = 0; i<newdataLambda.GetNbinsY(); i++)
        binEdgesNewData.push_back(newdataLambda.GetYaxis()->GetBinLowEdge(i+1));
    binEdgesNewData.push_back(newdataLambda.GetYaxis()->GetBinLowEdge(newdataLambda.GetNbinsY())+newdataLambda.GetYaxis()->GetBinWidth(newdataLambda.GetNbinsY()));

    olddataLambda1D = TH1D("olddataLambda1D", "olddataLambda1D", olddataLambda.GetNbinsY(), binEdgesOldData.data());
    oldtracksLambda1D = TH1D("oldtracksLambda1D", "oldtracksLambda1D", oldtracksLambda.GetNbinsY(), binEdgesOldTracks.data());
    newdataLambda1D = TH1D("newdataLambda1D", "newdataLambda1D", newdataLambda.GetNbinsY(), binEdgesNewData.data());

    for(int i = 0; i<olddataLambda.GetNbinsY(); i++){
        TH1D resultLambdaHistOldData = *olddataLambda.ProjectionX("Lambda1DHist", i+1, i+1);
        TH1D resultLambdaHistNewDataOldTracks = *oldtracksLambda.ProjectionX("Lambda1DHist", i+1, i+1);
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
        resultLambdaHistOldData.SetMinimum(0);
        resultLambdaHistOldData.SetMaximum(max(resultLambdaHistOldData.GetMaximum(), max(resultLambdaHistNewDataOldTracks.GetMaximum(), resultLambdaHistNewDataNewTracks.GetMaximum()))*1.1);
        resultLambdaHistOldData.SetMarkerStyle(kFullCircle);
        resultLambdaHistOldData.SetMarkerColor(kRed);
        resultLambdaHistOldData.SetLineColor(kRed+2);
        resultLambdaHistOldData.Fit(GfitLambda, "0BR");
        // TH1D *tempHist = (TH1D *)resultLambdaHistOldData.DrawCopy("E", "OldData");
        resultLambdaHistOldData.DrawCopy("E", "OldData");
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
        GfitLambda1Sig->SetLineColor(kRed);
        GfitLambda->SetNpx(1000);
        TF1 *oldDataFit = GfitLambda->DrawCopy("CSAME");
        TF1 *oldDataBcgFit = GfitLambdaBcg->DrawCopy("CSAME");
        TF1 *oldDataSigFit = GfitLambda1Sig->DrawCopy("CSAME");
        oldDataSigFit->SetParErrors(GfitLambda->GetParErrors());
        //new data old tracks
        resultLambdaHistNewDataOldTracks.SetMinimum(0);
        resultLambdaHistNewDataOldTracks.SetMarkerStyle(kFullCircle);
        resultLambdaHistNewDataOldTracks.SetMarkerColor(kGreen);
        resultLambdaHistNewDataOldTracks.SetLineColor(kGreen+2);
        // GfitLambda->SetParameters(100, 1.115, 5e-3, 2500000, -4200000, 1000000);
        GfitLambda->SetParameters(3, 1.115, 5e-3, 100, 0, 0);
        resultLambdaHistNewDataOldTracks.Fit(GfitLambda, "0BR");
        resultLambdaHistNewDataOldTracks.DrawCopy("E SAME", "NewDataOldTracks");
        gPad->Update();
        GfitLambda->GetParameters(paramsLambda);
        GfitLambdaBcg->SetParameters(paramsLambda+3);
        GfitLambda1Sig->SetParameters(paramsLambda);
        GfitLambdaBcg->SetLineStyle(3);
        GfitLambdaBcg->SetLineWidth(3);
        GfitLambda->SetLineColor(kGreen);
        GfitLambdaBcg->SetLineColor(kGreen);
        GfitLambda1Sig->SetLineColor(kGreen);
        GfitLambda->SetNpx(1000);
        TF1 *newDataOldTracksFit = GfitLambda->DrawCopy("CSAME");
        TF1 *newDataOldTracksBcgFit = GfitLambdaBcg->DrawCopy("CSAME");
        TF1 *newDataOldTracksSigFit = GfitLambda1Sig->DrawCopy("CSAME");
        newDataOldTracksSigFit->SetParErrors(GfitLambda->GetParErrors());
        //new data new tracks
        resultLambdaHistNewDataNewTracks.SetMinimum(0);
        resultLambdaHistNewDataNewTracks.SetMarkerStyle(kFullCircle);
        resultLambdaHistNewDataNewTracks.SetMarkerColor(kBlue);
        // GfitLambda->SetParameters(100, 1.115, 5e-3, 2000000, -4860000, 1500000);
        // GfitLambda->SetParameters(6, 1.115, 5e-3, 900000, -1600000, 780000);
        GfitLambda->SetParameters(3, 1.115, 5e-3, 100, 0, 0);
        resultLambdaHistNewDataNewTracks.Fit(GfitLambda, "0BR");
        resultLambdaHistNewDataNewTracks.DrawCopy("E SAME", "NewDataNewTracks");
        gPad->Update();
        GfitLambda->GetParameters(paramsLambda);
        GfitLambdaBcg->SetParameters(paramsLambda+3);
        GfitLambda1Sig->SetParameters(paramsLambda);
        GfitLambdaBcg->SetLineStyle(3);
        GfitLambdaBcg->SetLineWidth(3);
        GfitLambda->SetLineColor(kBlue);
        GfitLambdaBcg->SetLineColor(kBlue);
        GfitLambda1Sig->SetLineColor(kBlue);
        GfitLambda->SetNpx(1000);
        TF1 *newDataNewTracksFit = GfitLambda->DrawCopy("CSAME");
        TF1 *newDataNewTracksBcgFit = GfitLambdaBcg->DrawCopy("CSAME");
        TF1 *newDataNewTracksSigFit = GfitLambda1Sig->DrawCopy("CSAME");
        newDataNewTracksSigFit->SetParErrors(GfitLambda->GetParErrors());

        TLegend *legendLambda = new TLegend(.18, .65, .45, .89);
        legendLambda->SetTextSize(0.02);
        legendLambda->AddEntry("Lambda1DHistOldData", "pp old data, #sqrt{s} = 510 GeV");
        legendLambda->AddEntry("Lambda1DHistNewDataOldTracks", "pp new data, old tracks, #sqrt{s} = 510 GeV");
        legendLambda->AddEntry("Lambda1DHistNewDataNewTracks", "pp new data, new tracks, #sqrt{s} = 510 GeV");
        legendLambda->SetBorderSize(0);
        legendLambda->Draw("SAME");
        gPad->Update();

        //interactive part
        std::cout<<"Enter which data (0,1,2), which parameter (0-4) and how much change"<<std::endl;
        std::cout<<"Writing \"abs\" before just sets the parameter"<<std::endl;
        std::cout<<"Or write \"fit <data_number>\" to fit (3 for fitting all)"<<std::endl;
        std::cout<<"Or write \"copy <data_number_from> <data_number_to>\" to copy"<<std::endl;
        std::cout<<"Or write \"zero\" if there is no data suitable to fit"<<std::endl;
        std::cout<<"Or press \"Enter\" if everything ok"<<std::endl;
        string response;
        //needed to read whole line, not just until first space
        std::getline(std::cin, response);
        while(response.size()!=0&&response.find("zero")==string::npos){
            //fitting/changing part
            if(response.find("fit")!=string::npos){
                switch(response[4]-'0'){ //presents response as int by substracting ASCII code for 0
                case 0:
                    GfitLambda->SetParameters(oldDataFit->GetParameters());
                    resultLambdaHistOldData.Fit(GfitLambda, "0BR");
                    oldDataFit->SetParameters(GfitLambda->GetParameters());
                    oldDataBcgFit->SetParameters(GfitLambda->GetParameters()+3);
                    oldDataSigFit->SetParameters(GfitLambda->GetParameters());
                    oldDataSigFit->SetParErrors(GfitLambda->GetParErrors());
                    break;
                case 1:
                    GfitLambda->SetParameters(newDataOldTracksFit->GetParameters());
                    resultLambdaHistNewDataOldTracks.Fit(GfitLambda, "0BR");
                    newDataOldTracksFit->SetParameters(GfitLambda->GetParameters());
                    newDataOldTracksBcgFit->SetParameters(GfitLambda->GetParameters()+3);
                    newDataOldTracksSigFit->SetParameters(GfitLambda->GetParameters());
                    newDataOldTracksSigFit->SetParErrors(GfitLambda->GetParErrors());
                    break;
                case 2:
                    GfitLambda->SetParameters(newDataNewTracksFit->GetParameters());
                    resultLambdaHistNewDataNewTracks.Fit(GfitLambda, "0BR");
                    newDataNewTracksFit->SetParameters(GfitLambda->GetParameters());
                    newDataNewTracksBcgFit->SetParameters(GfitLambda->GetParameters()+3);
                    newDataNewTracksSigFit->SetParameters(GfitLambda->GetParameters());
                    newDataNewTracksSigFit->SetParErrors(GfitLambda->GetParErrors());
                    break;
                case 3:
                    GfitLambda->SetParameters(oldDataFit->GetParameters());
                    resultLambdaHistOldData.Fit(GfitLambda, "0BR");
                    oldDataFit->SetParameters(GfitLambda->GetParameters());
                    oldDataBcgFit->SetParameters(GfitLambda->GetParameters()+3);
                    oldDataSigFit->SetParameters(GfitLambda->GetParameters());
                    oldDataSigFit->SetParErrors(GfitLambda->GetParErrors());
                    GfitLambda->SetParameters(newDataOldTracksFit->GetParameters());
                    resultLambdaHistNewDataOldTracks.Fit(GfitLambda, "0BR");
                    newDataOldTracksFit->SetParameters(GfitLambda->GetParameters());
                    newDataOldTracksBcgFit->SetParameters(GfitLambda->GetParameters()+3);
                    newDataOldTracksSigFit->SetParameters(GfitLambda->GetParameters());
                    newDataOldTracksSigFit->SetParErrors(GfitLambda->GetParErrors());
                    GfitLambda->SetParameters(newDataNewTracksFit->GetParameters());
                    resultLambdaHistNewDataNewTracks.Fit(GfitLambda, "0BR");
                    newDataNewTracksFit->SetParameters(GfitLambda->GetParameters());
                    newDataNewTracksBcgFit->SetParameters(GfitLambda->GetParameters()+3);
                    newDataNewTracksSigFit->SetParameters(GfitLambda->GetParameters());
                    newDataNewTracksSigFit->SetParErrors(GfitLambda->GetParErrors());
                    break;
                default:
                    std::cout<<"Whoops, something went wrong with fit"<<std::endl;
                    break;
                }
            } else if(response.find("copy")!=string::npos){
                int copy_from, copy_to;
                TF1 *fit_tab[] = { oldDataFit, newDataOldTracksFit, newDataNewTracksFit };
                sscanf(response.c_str(), "copy %d%d", &copy_from, &copy_to);
                GfitLambda->SetParameters(fit_tab[copy_from]->GetParameters());
                fit_tab[copy_to]->SetParameters(GfitLambda->GetParameters());
            } else if(response.find("zero")!=string::npos){
                oldDataSigFit->SetParameter(0, 0);
                oldDataSigFit->SetParError(0, 0);
                newDataOldTracksSigFit->SetParameter(0, 0);
                newDataOldTracksSigFit->SetParError(0, 0);
                newDataNewTracksSigFit->SetParameter(0, 0);
                newDataNewTracksSigFit->SetParError(0, 0);
            } else{
                int data_number, coeff_number;
                double change, relative;
                if(response.find("abs")==string::npos){
                    relative = 1.;
                    sscanf(response.c_str(), "%d%d%lf", &data_number, &coeff_number, &change);
                } else{
                    relative = 0.;
                    sscanf(response.c_str(), "abs %d%d%lf", &data_number, &coeff_number, &change);
                }
                switch(data_number){
                case 0:
                    oldDataFit->SetParameter(coeff_number, oldDataFit->GetParameter(coeff_number)*relative+change);
                    if(coeff_number>2)
                        oldDataBcgFit->SetParameter(coeff_number-3, oldDataBcgFit->GetParameter(coeff_number-3)*relative+change);
                    else
                        oldDataSigFit->SetParameter(coeff_number, oldDataSigFit->GetParameter(coeff_number)*relative+change);
                    break;
                case 1:
                    newDataOldTracksFit->SetParameter(coeff_number, newDataOldTracksFit->GetParameter(coeff_number)*relative+change);
                    if(coeff_number>2)
                        newDataOldTracksBcgFit->SetParameter(coeff_number-3, newDataOldTracksBcgFit->GetParameter(coeff_number-3)*relative+change);
                    else
                        newDataOldTracksSigFit->SetParameter(coeff_number, newDataOldTracksSigFit->GetParameter(coeff_number)*relative+change);
                    break;
                case 2:
                    newDataNewTracksFit->SetParameter(coeff_number, newDataNewTracksFit->GetParameter(coeff_number)*relative+change);
                    if(coeff_number>2)
                        newDataNewTracksBcgFit->SetParameter(coeff_number-3, newDataNewTracksBcgFit->GetParameter(coeff_number-3)*relative+change);
                    else
                        newDataNewTracksSigFit->SetParameter(coeff_number, newDataNewTracksSigFit->GetParameter(coeff_number)*relative+change);
                    break;
                default:
                    std::cout<<"Whoops, something went wrong with param change"<<std::endl;
                    break;
                }
            }
            //drawing part
            oldDataFit->Draw("CSAME");
            oldDataBcgFit->Draw("CSAME");
            oldDataSigFit->Draw("CSAME");
            newDataOldTracksFit->Draw("CSAME");
            newDataOldTracksBcgFit->Draw("CSAME");
            newDataOldTracksSigFit->Draw("CSAME");
            newDataNewTracksFit->Draw("CSAME");
            newDataNewTracksBcgFit->Draw("CSAME");
            newDataNewTracksSigFit->Draw("CSAME");
            gPad->Update();
            //interactive part
            std::cout<<"Enter which data (0,1,2), which parameter (0-4) and how much change"<<std::endl;
            std::cout<<"Writing \"abs\" before just sets the parameter"<<std::endl;
            std::cout<<"Or write \"fit <data_number>\" to fit (3 for fitting all)"<<std::endl;
            std::cout<<"Or write \"copy <data_number_from> <data_number_to>\" to copy"<<std::endl;
            std::cout<<"Or write \"zero\" if there is no data suitable to fit"<<std::endl;
            std::cout<<"Or press \"Enter\" if everything ok"<<std::endl;
            //needed to read whole line, not just until first space
            std::getline(std::cin, response);
        }
        //adding to the result histogram
        if(response.find("zero")!=string::npos){
            olddataLambda1D.SetBinContent(i+1, 0.);
            olddataLambda1D.SetBinError(i+1, 0.);
            oldtracksLambda1D.SetBinContent(i+1, 0.);
            oldtracksLambda1D.SetBinError(i+1, 0.);
            newdataLambda1D.SetBinContent(i+1, 0.);
            newdataLambda1D.SetBinError(i+1, 0.);
        } else{
            olddataLambda1D.SetBinContent(i+1, oldDataSigFit->GetParameter(0)/olddataLambda.GetXaxis()->GetBinWidth(0)*olddataLambda.GetYaxis()->GetBinWidth(i+1));
            olddataLambda1D.SetBinError(i+1, oldDataSigFit->GetParError(0)/olddataLambda.GetXaxis()->GetBinWidth(0)*olddataLambda.GetYaxis()->GetBinWidth(i+1));
            oldtracksLambda1D.SetBinContent(i+1, newDataOldTracksSigFit->GetParameter(0)/oldtracksLambda.GetXaxis()->GetBinWidth(0)*oldtracksLambda.GetYaxis()->GetBinWidth(i+1));
            oldtracksLambda1D.SetBinError(i+1, newDataOldTracksSigFit->GetParError(0)/oldtracksLambda.GetXaxis()->GetBinWidth(0)*oldtracksLambda.GetYaxis()->GetBinWidth(i+1));
            newdataLambda1D.SetBinContent(i+1, newDataNewTracksSigFit->GetParameter(0)/newdataLambda.GetXaxis()->GetBinWidth(0)*newdataLambda.GetYaxis()->GetBinWidth(i+1));
            newdataLambda1D.SetBinError(i+1, newDataNewTracksSigFit->GetParError(0)/newdataLambda.GetXaxis()->GetBinWidth(0)*newdataLambda.GetYaxis()->GetBinWidth(i+1));
        }
        printf("Accepted current parameters\n");
        if(i==olddataLambda.GetNbinsY()-1)
            resultLambda->Close();
    }
}