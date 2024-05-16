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

int main(int argc, char *argv[]){
    //something used so that the histograms would draw
    //https://stackoverflow.com/questions/30932725/painting-a-tcanvas-to-the-screen-in-a-compiled-root-cern-application
    TApplication theApp("App", &argc, argv);



    // delete fileinfo;
    // delete filedialog;
    TFile *anaoutputolddata1 = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/CPT/AnaOutput_Inclusive_analysis_with_STUPCV0_with_extended_range_noAfterburner.root");
    TFile *anaoutputolddata2 = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/CPTnoBBCL/AnaOutput_Inclusive_analysis_with_STUPCV0_with_extended_range_noAfterburner.root");
    TFile *anaoutputoldtracks1 = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/CPT/AnaOutput_Inclusive_analysis_with_STUPCV0_with_extended_range_noAfterburner_new_data_old_tracks.root");
    TFile *anaoutputoldtracks2 = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/CPTnoBBCL/AnaOutput_Inclusive_analysis_with_STUPCV0_with_extended_range_noAfterburner_new_data_old_tracks.root");
    TFile *anaoutputnewtracks1 = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/CPT/AnaOutput_Inclusive_analysis_with_STUPCV0_with_extended_range_noAfterburner_new_data_new_tracks.root");
    TFile *anaoutputnewtracks2 = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/CPTnoBBCL/AnaOutput_Inclusive_analysis_with_STUPCV0_with_extended_range_noAfterburner_new_data_new_tracks.root");

    gStyle->SetFrameLineWidth(1);
    // gStyle->SetOptFit(111);
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);

    TH2D *anaoutputolddata1_2dpt = (TH2D *)anaoutputolddata1->Get("K0pt2DHist");
    TH2D *anaoutputolddata2_2dpt = (TH2D *)anaoutputolddata2->Get("K0pt2DHist");
    TH2D *anaoutputoldtracks1_2dpt = (TH2D *)anaoutputoldtracks1->Get("K0pt2DHist");
    TH2D *anaoutputoldtracks2_2dpt = (TH2D *)anaoutputoldtracks2->Get("K0pt2DHist");
    TH2D *anaoutputnewtracks1_2dpt = (TH2D *)anaoutputnewtracks1->Get("K0pt2DHist");
    TH2D *anaoutputnewtracks2_2dpt = (TH2D *)anaoutputnewtracks2->Get("K0pt2DHist");

    std::vector<double> binEdgesOldData, binEdgesOldTracks, binEdgesNewData;
    for(size_t i = 0; i<anaoutputolddata1_2dpt->GetNbinsY(); i++)
        binEdgesOldData.push_back(anaoutputolddata1_2dpt->GetYaxis()->GetBinLowEdge(i+1));
    binEdgesOldData.push_back(anaoutputolddata1_2dpt->GetYaxis()->GetBinLowEdge(anaoutputolddata1_2dpt->GetNbinsY())+anaoutputolddata1_2dpt->GetYaxis()->GetBinWidth(anaoutputolddata1_2dpt->GetNbinsY()));
    for(size_t i = 0; i<anaoutputoldtracks1_2dpt->GetNbinsY(); i++)
        binEdgesOldTracks.push_back(anaoutputoldtracks1_2dpt->GetYaxis()->GetBinLowEdge(i+1));
    binEdgesOldTracks.push_back(anaoutputoldtracks1_2dpt->GetYaxis()->GetBinLowEdge(anaoutputoldtracks1_2dpt->GetNbinsY())+anaoutputoldtracks1_2dpt->GetYaxis()->GetBinWidth(anaoutputoldtracks1_2dpt->GetNbinsY()));
    for(size_t i = 0; i<anaoutputnewtracks1_2dpt->GetNbinsY(); i++)
        binEdgesNewData.push_back(anaoutputnewtracks1_2dpt->GetYaxis()->GetBinLowEdge(i+1));
    binEdgesNewData.push_back(anaoutputnewtracks1_2dpt->GetYaxis()->GetBinLowEdge(anaoutputnewtracks1_2dpt->GetNbinsY())+anaoutputnewtracks1_2dpt->GetYaxis()->GetBinWidth(anaoutputnewtracks1_2dpt->GetNbinsY()));

    TH1D *olddataK0pt = new TH1D("olddataK0pt", "olddataK0pt", anaoutputolddata1_2dpt->GetNbinsY(), binEdgesOldData.data());
    TH1D *oldtracksK0pt = new TH1D("oldtracksK0pt", "oldtracksK0pt", anaoutputoldtracks1_2dpt->GetNbinsY(), binEdgesOldTracks.data());
    TH1D *newdataK0pt = new TH1D("newdataK0pt", "newdataK0pt", anaoutputnewtracks1_2dpt->GetNbinsY(), binEdgesNewData.data());

    for(int i = 0; i<anaoutputolddata1_2dpt->GetNbinsY(); i++){
        //old data
        TH1D *MKinvClose1 = anaoutputolddata1_2dpt->ProjectionX("K0pt1DHist", i+1, i+1);
        TH1D *MKinvClose2 = anaoutputolddata2_2dpt->ProjectionX("K0pt1DHist", i+1, i+1);
        TH1D resultKHistOldData = *MKinvClose1;
        resultKHistOldData.Add(MKinvClose2);
        // resultKHistOldData.SetStats(1);
        //new data, old tracks
        MKinvClose1 = anaoutputoldtracks1_2dpt->ProjectionX("K0pt1DHist", i+1, i+1);
        MKinvClose2 = anaoutputoldtracks2_2dpt->ProjectionX("K0pt1DHist", i+1, i+1);
        TH1D resultKHistNewDataOldTracks = *MKinvClose1;
        resultKHistNewDataOldTracks.Add(MKinvClose2);
        // resultKHistNewDataOldTracks.SetStats(1);
        //new data, new tracks
        MKinvClose1 = anaoutputnewtracks1_2dpt->ProjectionX("K0pt1DHist", i+1, i+1);
        MKinvClose2 = anaoutputnewtracks2_2dpt->ProjectionX("K0pt1DHist", i+1, i+1);
        TH1D resultKHistNewDataNewTracks = *MKinvClose1;
        resultKHistNewDataNewTracks.Add(MKinvClose2);
        // resultKHistNewDataNewTracks.SetStats(1);

        //K0
        //setting up the functions and data
        double lowrangefit = 0.48, upperrangefit = 0.52;
        TCanvas *resultK = new TCanvas("resultK", "resultK", 1600, 800);
        resultK->SetLeftMargin(0.15);
        TF1 *GfitK = new TF1("GfitK", "gausn(0) + pol1(3)", lowrangefit, upperrangefit);
        TF1 *GfitKBcg = new TF1("GfitKBcg", "pol1", lowrangefit, upperrangefit);
        TF1 *GfitK1Sig = new TF1("GfitK1Sig", "gausn", lowrangefit, upperrangefit);
        GfitK->SetParNames("Constant", "Mean", "Sigma", "b", "a");
        GfitK->SetParameters(100, 0.5, 5e-3, 100, 0.48, 5e-3);
        double nK0OldData, nK0NewDataOldTracks, nK0NewDataNewTracks;
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
        // nK0OldData = paramsK[0];
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
        // nK0NewDataOldTracks = paramsK[0];
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
        // nK0NewDataNewTracks = paramsK[0];
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
        legendK->AddEntry("K0pt1DHistOldData", "pp old data, #sqrt{s} = 510 GeV");
        legendK->AddEntry("K0pt1DHistNewDataOldTracks", "pp new data, old tracks, #sqrt{s} = 510 GeV");
        legendK->AddEntry("K0pt1DHistNewDataNewTracks", "pp new data, new tracks, #sqrt{s} = 510 GeV");
        legendK->SetBorderSize(0);
        legendK->Draw("SAME");
        gPad->Update();

        //interactive part
        std::cout<<"Enter which data (0,1,2), which parameter (0-4) and how much change"<<std::endl;
        std::cout<<"Writing \"abs\" before just sets the parameter"<<std::endl;
        std::cout<<"Or write \"fit <data_number>\" to fit"<<std::endl;
        std::cout<<"Or write \"ok\" if everything ok"<<std::endl;
        string response;
        //needed to read whole line, not just until first space
        std::getline(std::cin, response);
        while(response.find("ok")==string::npos){
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
                default:
                    std::cout<<"Whoops, something went wrong with fit"<<std::endl;
                    break;
                }
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
            std::cout<<"Or write \"fit <data_number>\" to fit"<<std::endl;
            std::cout<<"Or write \"ok\" if everything ok"<<std::endl;
            //needed to read whole line, not just until first space
            std::getline(std::cin, response);
        }
        //adding to the result histogram
        olddataK0pt->SetBinContent(i+1, oldDataSigFit->GetParameter(0));
        olddataK0pt->SetBinError(i+1, oldDataSigFit->GetParError(0));
        oldtracksK0pt->SetBinContent(i+1, newDataOldTracksSigFit->GetParameter(0));
        oldtracksK0pt->SetBinError(i+1, newDataOldTracksSigFit->GetParError(0));
        newdataK0pt->SetBinContent(i+1, newDataNewTracksSigFit->GetParameter(0));
        newdataK0pt->SetBinError(i+1, newDataNewTracksSigFit->GetParError(0));
        printf("Accepted current parameters\n");
    }

    TFile *outputFileHist = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/Differential_crossection_results.root", "recreate");
    outputFileHist->cd();
    olddataK0pt->Write();
    oldtracksK0pt->Write();
    newdataK0pt->Write();
    outputFileHist->Close();

    // gPad->RedrawAxis();
    TCanvas *resultKpt = new TCanvas("resultKpt", "resultKpt", 1600, 800);
    olddataK0pt->SetTitle("Differential crossection for K^{  0}_{  S} with respect to p_{T}");
    olddataK0pt->SetMarkerStyle(kFullCircle);
    olddataK0pt->SetMarkerColor(kRed);
    olddataK0pt->SetLineColor(kRed+2);
    olddataK0pt->Draw("E");
    oldtracksK0pt->SetMarkerStyle(kFullCircle);
    oldtracksK0pt->SetMarkerColor(kGreen);
    oldtracksK0pt->SetLineColor(kGreen+2);
    oldtracksK0pt->Draw("ESAME");
    newdataK0pt->SetMarkerStyle(kFullCircle);
    newdataK0pt->SetMarkerColor(kBlue);
    newdataK0pt->SetLineColor(kBlue+2);
    newdataK0pt->Draw("ESAME");
    TLegend *legendKpt = new TLegend(.62, .65, .89, .89);
    legendKpt->SetTextSize(0.02);
    legendKpt->AddEntry("olddataK0pt", "pp old data, #sqrt{s} = 510 GeV");
    legendKpt->AddEntry("oldtracksK0pt", "pp new data, old tracks, #sqrt{s} = 510 GeV");
    legendKpt->AddEntry("newdataK0pt", "pp new data, new tracks, #sqrt{s} = 510 GeV");
    legendKpt->SetBorderSize(0);
    legendKpt->Draw("SAME");
    resultKpt->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/Differential_crossection_K0_pt.pdf").c_str());
    theApp.Run();
    return 0;


    //Lambda
    //old data
    TH1D *MLambdainvClose1 = (TH1D *)anaoutputolddata1->Get("MppiWide");
    TH1D *MLambdainvClose2 = (TH1D *)anaoutputolddata2->Get("MppiWide");
    TH1D resultLambdaHistOldData = *MLambdainvClose1;
    resultLambdaHistOldData.Add(MLambdainvClose2);
    // resultLambdaHistOldData.SetStats(1);
    //new data, old tracks
    MLambdainvClose1 = (TH1D *)anaoutputoldtracks1->Get("MppiWide");
    MLambdainvClose2 = (TH1D *)anaoutputoldtracks2->Get("MppiWide");
    TH1D resultLambdaHistNewDataOldTracks = *MLambdainvClose1;
    resultLambdaHistNewDataOldTracks.Add(MLambdainvClose2);
    // resultLambdaHistNewDataOldTracks.SetStats(1);
    //new data, new tracks
    MLambdainvClose1 = (TH1D *)anaoutputnewtracks1->Get("MppiWide");
    MLambdainvClose2 = (TH1D *)anaoutputnewtracks2->Get("MppiWide");
    TH1D resultLambdaHistNewDataNewTracks = *MLambdainvClose1;
    resultLambdaHistNewDataNewTracks.Add(MLambdainvClose2);
    // resultLambdaHistNewDataNewTracks.SetStats(1);

    //Lambda
    //setting up the functions and data
    double fitmin = 1.1;
    double fitmax = 1.14;
    TCanvas *resultLambda = new TCanvas("resultLambda", "resultLambda", 1800, 1600);
    resultLambda->SetLeftMargin(0.15);
    TF1 *GfitLambda = new TF1("GfitLambda", "gausn(0) + pol2(3)", fitmin, fitmax);
    TF1 *GfitLambdaBcg = new TF1("GfitLambdaBcg", "pol2", fitmin, fitmax);
    TF1 *GfitLambda1Sig = new TF1("GfitLambda1Sig", "gausn", fitmin, fitmax);
    GfitLambda->SetParNames("Constant", "Mean", "Sigma", "c", "b", "a");
    // GfitLambda->SetParameters(100, 1.115, 5e-3, 2500000, -4200000, 1200000);
    GfitLambda->SetParameters(3, 1.115, 5e-3, 900000, -1600000, 780000);
    GfitLambda->SetParLimits(0, 0, 200);
    GfitLambda->SetParLimits(1, 1.11, 1.12);
    GfitLambda->SetParLimits(2, 5e-5, 0.01);
    double nLambda0OldData, nLambda0NewDataOldTracks, nLambda0NewDataNewTracks;
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
    nLambda0OldData = paramsLambda[0];
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
    GfitLambda->DrawCopy("CSAME");
    GfitLambdaBcg->DrawCopy("CSAME");
    GfitLambda1Sig->DrawCopy("CSAME");
    //new data old tracks
    resultLambdaHistNewDataOldTracks.SetMinimum(0);
    resultLambdaHistNewDataOldTracks.SetMarkerStyle(kFullCircle);
    resultLambdaHistNewDataOldTracks.SetMarkerColor(kGreen);
    resultLambdaHistNewDataOldTracks.SetLineColor(kGreen+2);
    GfitLambda->SetParameters(100, 1.115, 5e-3, 2500000, -4200000, 1000000);
    resultLambdaHistNewDataOldTracks.Fit(GfitLambda, "0BR");
    resultLambdaHistNewDataOldTracks.DrawCopy("E SAME", "NewDataOldTracks");
    gPad->Update();
    GfitLambda->GetParameters(paramsLambda);
    nLambda0NewDataOldTracks = paramsLambda[0];
    GfitLambdaBcg->SetParameters(paramsLambda+3);
    GfitLambda1Sig->SetParameters(paramsLambda);
    GfitLambdaBcg->SetLineStyle(3);
    GfitLambdaBcg->SetLineWidth(3);
    GfitLambda->SetLineColor(kGreen);
    GfitLambdaBcg->SetLineColor(kGreen);
    GfitLambda1Sig->SetLineColor(kGreen);
    GfitLambda->SetNpx(1000);
    GfitLambda->DrawCopy("CSAME");
    GfitLambdaBcg->DrawCopy("CSAME");
    GfitLambda1Sig->DrawCopy("CSAME");
    //new data new tracks
    resultLambdaHistNewDataNewTracks.SetMinimum(0);
    resultLambdaHistNewDataNewTracks.SetMarkerStyle(kFullCircle);
    resultLambdaHistNewDataNewTracks.SetMarkerColor(kBlue);
    // GfitLambda->SetParameters(100, 1.115, 5e-3, 2000000, -4860000, 1500000);
    GfitLambda->SetParameters(6, 1.115, 5e-3, 900000, -1600000, 780000);
    resultLambdaHistNewDataNewTracks.Fit(GfitLambda, "0BR");
    resultLambdaHistNewDataNewTracks.DrawCopy("E SAME", "NewDataNewTracks");
    gPad->Update();
    GfitLambda->GetParameters(paramsLambda);
    nLambda0NewDataNewTracks = paramsLambda[0];
    GfitLambdaBcg->SetParameters(paramsLambda+3);
    GfitLambda1Sig->SetParameters(paramsLambda);
    GfitLambdaBcg->SetLineStyle(3);
    GfitLambdaBcg->SetLineWidth(3);
    GfitLambda->SetLineColor(kBlue);
    GfitLambdaBcg->SetLineColor(kBlue);
    GfitLambda1Sig->SetLineColor(kBlue);
    GfitLambda->SetNpx(1000);
    GfitLambda->DrawCopy("CSAME");
    GfitLambdaBcg->DrawCopy("CSAME");
    GfitLambda1Sig->DrawCopy("CSAME");

    TLegend *legendLambda = new TLegend(.18, .65, .45, .89);
    legendLambda->SetTextSize(0.02);
    legendLambda->AddEntry("MppiWideOldData", "pp old data, #sqrt{s} = 510 GeV");
    legendLambda->AddEntry("MppiWideNewDataOldTracks", "pp new data, old tracks, #sqrt{s} = 510 GeV");
    legendLambda->AddEntry("MppiWideNewDataNewTracks", "pp new data, new tracks, #sqrt{s} = 510 GeV");
    // legendLambda->AddEntry("GfitLambda", "Data fit", "l");
    // legendLambda->AddEntry("GfitLambdaBcg", "Background", "l");
    // legendLambda->AddEntry("GfitLambda1Sig", "K^{0} fit", "l");
    legendLambda->SetBorderSize(0);
    legendLambda->Draw("SAME");

    nLambda0OldData /= resultLambdaHistOldData.GetBinWidth(1);
    nLambda0NewDataOldTracks /= resultLambdaHistOldData.GetBinWidth(1);
    nLambda0NewDataNewTracks /= resultLambdaHistOldData.GetBinWidth(1);
    std::stringstream nLambda0OldDataStream, nLambda0NewDataOldTracksStream, nLambda0NewDataNewTracksStream;
    nLambda0OldDataStream<<std::fixed<<std::setprecision(2)<<nLambda0OldData;
    nLambda0NewDataOldTracksStream<<std::fixed<<std::setprecision(2)<<nLambda0NewDataOldTracks;
    nLambda0NewDataNewTracksStream<<std::fixed<<std::setprecision(2)<<nLambda0NewDataNewTracks;
    TPaveStats *pt = new TPaveStats(.18, .45, .45, .65, "NB NDC");
    pt->SetFillColor(kWhite);
    TText *tempText = pt->AddText(0, 0, "Number of #Lambda^{0} detected in:");
    tempText->SetTextSize(3);
    pt->AddText(0, 0.33, ("Old data = "+nLambda0OldDataStream.str()).c_str());
    pt->AddText(0, 0.66, ("New data, old tracks = "+nLambda0NewDataOldTracksStream.str()).c_str());
    pt->AddText(0, 1, ("New data, new tracks = "+nLambda0NewDataNewTracksStream.str()).c_str());
    // pt->AddText("They are added to the pave using the AddText method.");
    // pt->AddLine(.0, .5, 1., .5);
    // pt->AddText("Even complex TLatex formulas can be added:");
    // TText *t1 = pt->AddText("F(t) = #sum_{i=-#infty}^{#infty}A(i)cos#[]{#frac{i}{t+i}}");
    // t1->SetTextColor(kBlue);
    pt->Draw("SAME");

    // gPad->RedrawAxis();
    // resultLambda->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/Lambda0afterV0finder.pdf").c_str());

    theApp.Run();
    return 0;
}

// void extract_histograms(TH2D *olddata, TH2D *oldtracks, TH2D *newdata, string histName){

// }

// void K0_differential_crossection_fit(string histName){

// }