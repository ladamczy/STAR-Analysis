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
#include "TEllipse.h"

#include <iomanip>
#include <sstream>

using namespace std;

int main(){
    // TGFileInfo* fileinfo = new TGFileInfo();
    // fileinfo->fFileTypeIdx = 2;
    // TGFileDialog* filedialog = new TGFileDialog(gClient->GetRoot(), nullptr, EFileDialogMode::kFDOpen, fileinfo);

    // cout<<fileinfo->fFilename<<endl;
    // TFile* anaoutput = TFile::Open(fileinfo->fFilename);
    // cout<<"File opened"<<endl;

    // delete fileinfo;
    // delete filedialog;
    TFile *anaoutputnewtracks1 = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/CPT/AnaOutput_Inclusive_analysis_with_STUPCV0_with_extended_range_noAfterburner_new_data_new_tracks.root");
    TFile *anaoutputnewtracks2 = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/CPTnoBBCL/AnaOutput_Inclusive_analysis_with_STUPCV0_with_extended_range_noAfterburner_new_data_new_tracks.root");

    gStyle->SetFrameLineWidth(1);
    // gStyle->SetOptFit(111);
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);

    // resultKHistNewDataOldTracks.SetStats(1);
    //new data, new tracks
    TH1D *MKinvClose1 = (TH1D *)anaoutputnewtracks1->Get("MpipiWide");
    TH1D *MKinvClose2 = (TH1D *)anaoutputnewtracks2->Get("MpipiWide");
    TH1D resultKHistNewDataNewTracks = *MKinvClose1;
    resultKHistNewDataNewTracks.Add(MKinvClose2);
    // resultKHistNewDataNewTracks.SetStats(1);

    //K0
    //setting up the functions and data
    TCanvas *resultK = new TCanvas("resultK", "resultK", 1800, 1600);
    resultK->SetLeftMargin(0.15);
    TF1 *GfitK = new TF1("GfitK", "gausn(0) + pol1(3)", resultKHistNewDataNewTracks.GetXaxis()->GetXmin(), resultKHistNewDataNewTracks.GetXaxis()->GetXmax());
    TF1 *GfitKBcg = new TF1("GfitKBcg", "pol1", resultKHistNewDataNewTracks.GetXaxis()->GetXmin(), resultKHistNewDataNewTracks.GetXaxis()->GetXmax());
    TF1 *GfitK1Sig = new TF1("GfitK1Sig", "gausn", resultKHistNewDataNewTracks.GetXaxis()->GetXmin(), resultKHistNewDataNewTracks.GetXaxis()->GetXmax());
    GfitK->SetParNames("Constant", "Mean", "Sigma", "b", "a");
    GfitK->SetParameters(100, 0.5, 5e-3, 100, 0.48, 5e-3);
    double nK0NewDataNewTracks;
    //old data, the base of drawing
    resultKHistNewDataNewTracks.SetMinimum(0);
    resultKHistNewDataNewTracks.SetMarkerStyle(kFullCircle);
    resultKHistNewDataNewTracks.SetMarkerColor(kBlue);
    resultKHistNewDataNewTracks.SetLineColor(kBlue+2);
    resultKHistNewDataNewTracks.Fit(GfitK, "0");
    // TH1D *tempHist = (TH1D *)resultKHistNewDataNewTracks.DrawCopy("E", "OldData");
    resultKHistNewDataNewTracks.DrawCopy("E", "OldData");
    gPad->Update();
    // TPaveStats *st = (TPaveStats *)tempHist->FindObject("stats");
    // st->SetX1NDC(0.15);
    // st->SetX2NDC(0.35);
    // st->SetY1NDC(0.65);
    // st->SetY2NDC(0.85);
    // st->SetBorderSize(0);
    Double_t paramsK[5];
    GfitK->GetParameters(paramsK);
    nK0NewDataNewTracks = paramsK[0];
    // cout<<"PARAMETERS:"<<endl;
    // for(int i = 0; i<5; i++){
    //     cout<<paramsK[i]<<endl;
    // }
    GfitKBcg->SetParameters(paramsK+3);
    GfitK1Sig->SetParameters(paramsK);
    GfitKBcg->SetLineStyle(3);
    GfitKBcg->SetLineWidth(3);
    GfitK->SetLineColor(kBlue);
    GfitKBcg->SetLineColor(kBlue);
    GfitK1Sig->SetLineColor(kBlue);
    GfitK->SetNpx(1000);
    GfitK->DrawCopy("CSAME");
    GfitKBcg->DrawCopy("CSAME");
    GfitK1Sig->DrawCopy("CSAME");

    TLegend *legendK = new TLegend(.18, .65, .45, .89);
    legendK->SetTextSize(0.02);
    legendK->AddEntry("MpipiWideNewDataNewTracks", "pp new data, new tracks, #sqrt{s} = 510 GeV");
    legendK->SetBorderSize(0);
    legendK->Draw("SAME");

    nK0NewDataNewTracks /= resultKHistNewDataNewTracks.GetBinWidth(1);
    std::stringstream nK0NewDataNewTracksStream;
    nK0NewDataNewTracksStream<<std::fixed<<std::setprecision(2)<<nK0NewDataNewTracks;
    TPaveStats *pt = new TPaveStats(0.62, 0.7, 0.89, 0.89, "NB NDC");
    TText *tempText;
    pt->SetFillColor(kWhite);
    tempText = pt->AddText(0, 0, "Number of K^{0}_{S} detected in:");
    tempText->SetTextSize(3);
    tempText = pt->AddText(0, 1, ("New data = "+nK0NewDataNewTracksStream.str()).c_str());
    tempText->SetTextSize(0.029);
    // pt->AddText("They are added to the pave using the AddText method.");
    // pt->AddLine(.0, .5, 1., .5);
    // pt->AddText("Even complex TLatex formulas can be added:");
    // TText *t1 = pt->AddText("F(t) = #sum_{i=-#infty}^{#infty}A(i)cos#[]{#frac{i}{t+i}}");
    // t1->SetTextColor(kBlue);
    pt->Draw("SAME");

    // gPad->RedrawAxis();
    resultK->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/K0afterV0finder_new_data.pdf").c_str());



    //Lambda
    //new data, new tracks
    TH1D *MLambdainvClose1 = (TH1D *)anaoutputnewtracks1->Get("MppiWide");
    TH1D *MLambdainvClose2 = (TH1D *)anaoutputnewtracks2->Get("MppiWide");
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
    double nLambda0NewDataNewTracks;
    //old data, the base of drawing
    resultLambdaHistNewDataNewTracks.SetMinimum(0);
    resultLambdaHistNewDataNewTracks.SetMarkerStyle(kFullCircle);
    resultLambdaHistNewDataNewTracks.SetMarkerColor(kBlue);
    resultLambdaHistNewDataNewTracks.SetLineColor(kBlue+2);
    resultLambdaHistNewDataNewTracks.Fit(GfitLambda, "0BR");
    // TH1D *tempHist = (TH1D *)resultLambdaHistNewDataNewTracks.DrawCopy("E", "OldData");
    resultLambdaHistNewDataNewTracks.DrawCopy("E", "OldData");
    gPad->Update();
    // TPaveStats *st = (TPaveStats *)tempHist->FindObject("stats");
    // st->SetX1NDC(0.15);
    // st->SetX2NDC(0.35);
    // st->SetY1NDC(0.65);
    // st->SetY2NDC(0.85);
    // st->SetBorderSize(0);
    Double_t paramsLambda[6];
    GfitLambda->GetParameters(paramsLambda);
    nLambda0NewDataNewTracks = paramsLambda[0];
    // cout<<"PARAMETERS:"<<endl;
    // for(int i = 0; i<5; i++){
    //     cout<<paramsLambda[i]<<endl;
    // }
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
    legendLambda->AddEntry("MppiWideNewDataNewTracks", "pp new data, new tracks, #sqrt{s} = 510 GeV");
    // legendLambda->AddEntry("GfitLambda", "Data fit", "l");
    // legendLambda->AddEntry("GfitLambdaBcg", "Background", "l");
    // legendLambda->AddEntry("GfitLambda1Sig", "K^{0} fit", "l");
    legendLambda->SetBorderSize(0);
    legendLambda->Draw("SAME");

    nLambda0NewDataNewTracks /= resultLambdaHistNewDataNewTracks.GetBinWidth(1);
    std::stringstream nLambda0NewDataNewTracksStream;
    nLambda0NewDataNewTracksStream<<std::fixed<<std::setprecision(2)<<nLambda0NewDataNewTracks;
    pt = new TPaveStats(.18, .45, .45, .65, "NB NDC");
    pt->SetFillColor(kWhite);
    tempText = pt->AddText(0, 0, "Number of #Lambda^{0} detected in:");
    tempText->SetTextSize(3);
    tempText = pt->AddText(0, 1, ("New data = "+nLambda0NewDataNewTracksStream.str()).c_str());
    tempText->SetTextSize(0.029);
    // pt->AddText("They are added to the pave using the AddText method.");
    // pt->AddLine(.0, .5, 1., .5);
    // pt->AddText("Even complex TLatex formulas can be added:");
    // TText *t1 = pt->AddText("F(t) = #sum_{i=-#infty}^{#infty}A(i)cos#[]{#frac{i}{t+i}}");
    // t1->SetTextColor(kBlue);
    pt->Draw("SAME");

    // gPad->RedrawAxis();
    resultLambda->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/Lambda0afterV0finder_new_data.pdf").c_str());



    //rest of the histograms



    //NotherTracks
    TCanvas *c1 = new TCanvas("c1", "c1", 1800, 1600);
    c1->SetGrid();
    gPad->SetLogy();
    TH1D *NotherTracksTab[6];
    NotherTracksTab[4] = (TH1D *)anaoutputnewtracks1->Get("NotherTracks");
    NotherTracksTab[5] = (TH1D *)anaoutputnewtracks2->Get("NotherTracks");
    NotherTracksTab[4]->SetName("NotherTracks4");
    NotherTracksTab[5]->SetName("NotherTracks5");
    for(size_t i = 4; i<6; i++){
        NotherTracksTab[i]->Scale(1/NotherTracksTab[i]->Integral());
        NotherTracksTab[i]->SetLineWidth(2);
        if(i%2==1){
            NotherTracksTab[i]->SetLineStyle(3);
        }
        if(i<2){
            NotherTracksTab[i]->SetLineColor(kRed);
        } else if(i<4){
            NotherTracksTab[i]->SetLineColor(kGreen);
        } else{
            NotherTracksTab[i]->SetLineColor(kBlue);
        }
        NotherTracksTab[i]->SetStats(0);
    }
    NotherTracksTab[4]->Draw("HIST");
    for(size_t i = 5; i<6; i++){
        NotherTracksTab[i]->Draw("HIST SAME");
    }
    TLegend *legendNotherTracks = new TLegend(0.5, 0.6, 0.89, 0.89);
    legendNotherTracks->SetTextSize(0.02);
    legendNotherTracks->AddEntry("NotherTracks4", "CPT2, new data, new tracks", "l");
    legendNotherTracks->AddEntry("NotherTracks5", "CPT2noBBCL, new data, new tracks", "l");
    legendNotherTracks->SetBorderSize(1);
    legendNotherTracks->Draw("SAME");
    c1->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/NotherTracksAll_new_data.pdf").c_str());
    delete c1;

    //K0multiplicity
    c1 = new TCanvas("c1", "c1", 1800, 1600);
    c1->SetGrid();
    gPad->SetLogy();
    TH1D *K0multiplicityTab[6];
    K0multiplicityTab[4] = (TH1D *)anaoutputnewtracks1->Get("K0multiplicity");
    K0multiplicityTab[5] = (TH1D *)anaoutputnewtracks2->Get("K0multiplicity");
    K0multiplicityTab[4]->SetName("K0multiplicity4");
    K0multiplicityTab[5]->SetName("K0multiplicity5");
    for(size_t i = 4; i<6; i++){
        K0multiplicityTab[i]->Scale(1/K0multiplicityTab[i]->Integral());
        K0multiplicityTab[i]->SetLineWidth(2);
        if(i%2==1){
            K0multiplicityTab[i]->SetLineStyle(3);
        }
        if(i<2){
            K0multiplicityTab[i]->SetLineColor(kRed);
        } else if(i<4){
            K0multiplicityTab[i]->SetLineColor(kGreen);
        } else{
            K0multiplicityTab[i]->SetLineColor(kBlue);
        }
        K0multiplicityTab[i]->SetStats(0);
    }
    K0multiplicityTab[4]->Draw("HIST");
    for(size_t i = 5; i<6; i++){
        K0multiplicityTab[i]->Draw("HIST SAME");
    }
    TLegend *legendK0multiplicity = new TLegend(0.5, 0.6, 0.89, 0.89);
    legendK0multiplicity->AddEntry("K0multiplicity4", "CPT2, new data, new tracks", "l");
    legendK0multiplicity->AddEntry("K0multiplicity5", "CPT2noBBCL, new data, new tracks", "l");
    legendK0multiplicity->SetBorderSize(1);
    legendK0multiplicity->Draw("SAME");
    c1->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/K0multiplicityAll_new_data.pdf").c_str());
    delete c1;





    //K0 backgrounds
    TCanvas *resultKBkg = new TCanvas("resultKBkg", "resultKBkg", 1800, 1600);
    resultKBkg->SetLeftMargin(0.15);
    //new data, new tracks
    TH1D *MKinvClose1Bkg = (TH1D *)anaoutputnewtracks1->Get("MpipiWideMissing");
    TH1D *MKinvClose2Bkg = (TH1D *)anaoutputnewtracks2->Get("MpipiWideMissing");
    TH1D resultKHistNewDataNewTracksBkg = *MKinvClose1Bkg;
    resultKHistNewDataNewTracks.Add(MKinvClose2Bkg);
    resultKHistNewDataNewTracksBkg.SetTitle("#pi^{+}#pi^{-} pairs rejected by basic cuts");
    resultKHistNewDataNewTracksBkg.SetMinimum(600);
    resultKHistNewDataNewTracksBkg.SetMaximum(1100);
    //new data new tracks
    resultKHistNewDataNewTracksBkg.SetMarkerStyle(kFullCircle);
    resultKHistNewDataNewTracksBkg.SetMarkerColor(kBlue);
    resultKHistNewDataNewTracksBkg.DrawCopy("E", "NewDataNewTracksBkg");

    TLegend *legendKBkg = new TLegend(.18, .65, .45, .89);
    legendKBkg->SetTextSize(0.02);
    legendKBkg->AddEntry("MpipiWideMissingNewDataNewTracksBkg", "pp new data, #sqrt{s} = 510 GeV");
    legendKBkg->SetBorderSize(0);
    legendKBkg->Draw("SAME");

    resultKBkg->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/K0afterV0finderBkg_new_data.pdf").c_str());


    //K0 after PID
    //new data, new tracks
    TH1D *MKinvPIDClose1 = (TH1D *)anaoutputnewtracks1->Get("MpipiWideWithPidEcut");
    TH1D *MKinvPIDClose2 = (TH1D *)anaoutputnewtracks2->Get("MpipiWideWithPidEcut");
    TH1D resultKPIDHistNewDataNewTracks = *MKinvPIDClose1;
    resultKPIDHistNewDataNewTracks.Add(MKinvPIDClose2);
    //setting up the functions and data
    TCanvas *resultKPID = new TCanvas("resultKPID", "resultKPID", 1800, 1600);
    resultKPID->SetLeftMargin(0.15);
    TF1 *GfitKPID = new TF1("GfitKPID", "gausn(0) + pol1(3)", resultKPIDHistNewDataNewTracks.GetXaxis()->GetXmin(), resultKPIDHistNewDataNewTracks.GetXaxis()->GetXmax());
    TF1 *GfitKPIDBcg = new TF1("GfitKPIDBcg", "pol1", resultKPIDHistNewDataNewTracks.GetXaxis()->GetXmin(), resultKPIDHistNewDataNewTracks.GetXaxis()->GetXmax());
    TF1 *GfitKPID1Sig = new TF1("GfitKPID1Sig", "gausn", resultKPIDHistNewDataNewTracks.GetXaxis()->GetXmin(), resultKPIDHistNewDataNewTracks.GetXaxis()->GetXmax());
    GfitKPID->SetParNames("Constant", "Mean", "Sigma", "b", "a");
    GfitKPID->SetParameters(100, 0.5, 5e-3, 100, 0.48, 5e-3);
    // double nK0OldData, nK0NewDataOldTracks, nK0NewDataNewTracks;
    //old data, the base of drawing
    resultKPIDHistNewDataNewTracks.SetTitle("K^{0}_{S} mass in wide range with PID cuts");
    resultKPIDHistNewDataNewTracks.SetMinimum(0);
    resultKPIDHistNewDataNewTracks.SetMarkerStyle(kFullCircle);
    resultKPIDHistNewDataNewTracks.SetMarkerColor(kBlue);
    resultKPIDHistNewDataNewTracks.Fit(GfitKPID, "0");
    resultKPIDHistNewDataNewTracks.DrawCopy("E", "NewDataNewTracks");
    gPad->Update();
    GfitKPID->GetParameters(paramsK);
    nK0NewDataNewTracks = paramsK[0];
    GfitKPIDBcg->SetParameters(paramsK+3);
    GfitKPID1Sig->SetParameters(paramsK);
    GfitKPIDBcg->SetLineStyle(3);
    GfitKPIDBcg->SetLineWidth(3);
    GfitKPID->SetLineColor(kBlue);
    GfitKPIDBcg->SetLineColor(kBlue);
    GfitKPID1Sig->SetLineColor(kBlue);
    GfitKPID->SetNpx(1000);
    GfitKPID->DrawCopy("CSAME");
    GfitKPIDBcg->DrawCopy("CSAME");
    GfitKPID1Sig->DrawCopy("CSAME");

    TLegend *legendKPID = new TLegend(.18, .74, .45, .84);
    legendKPID->SetTextSize(0.02);
    legendKPID->SetHeader("#bf{pp, #sqrt{s} = 510 GeV}", "C");
    legendKPID->AddEntry("MpipiWideWithPidEcutNewDataNewTracks", "#bf{new data, new tracks}");
    legendKPID->SetBorderSize(0);
    legendKPID->Draw("SAME");

    nK0NewDataNewTracks /= resultKPIDHistNewDataNewTracks.GetBinWidth(1);
    std::stringstream nK0PIDNewDataNewTracksStream;
    nK0PIDNewDataNewTracksStream<<std::fixed<<std::setprecision(2)<<nK0NewDataNewTracks;
    // TPaveStats *pt = new TPaveStats(0.62, 0.7, 0.89, 0.89, "NB NDC");
    pt = new TPaveStats(0.62, 0.7, 0.89, 0.89, "NB NDC");
    // TText *tempText;
    pt->SetFillColor(kWhite);
    tempText = pt->AddText(0, 0, "Number of K^{0}_{S} detected in:");
    tempText->SetTextSize(3);
    tempText = pt->AddText(0, 1, ("New data = "+nK0PIDNewDataNewTracksStream.str()).c_str());
    tempText->SetTextSize(0.029);
    // pt->AddText("They are added to the pave using the AddText method.");
    // pt->AddLine(.0, .5, 1., .5);
    // pt->AddText("Even complex TLatex formulas can be added:");
    // TText *t1 = pt->AddText("F(t) = #sum_{i=-#infty}^{#infty}A(i)cos#[]{#frac{i}{t+i}}");
    // t1->SetTextColor(kBlue);
    pt->Draw("SAME");

    // gPad->RedrawAxis();
    resultKPID->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/K0afterV0finderPID_new_data.pdf").c_str());





    //K0 pointAngleHypo
    TCanvas *pointAngleHypoKCanvas = new TCanvas("pointAngleHypoKCanvas", "pointAngleHypoKCanvas", 1800, 1600);
    pointAngleHypoKCanvas->SetGrid();
    // gPad->SetLogy();
    TH1D *pointAngleHypoKTab[6];
    pointAngleHypoKTab[4] = (TH1D *)anaoutputnewtracks1->Get("K0pointingAngleHypoSignalDetailed");
    pointAngleHypoKTab[4]->Add((TH1D *)anaoutputnewtracks2->Get("K0pointingAngleHypoSignalDetailed"));
    pointAngleHypoKTab[5] = (TH1D *)anaoutputnewtracks1->Get("K0pointingAngleHypoBackgroundDetailed");
    pointAngleHypoKTab[5]->Add((TH1D *)anaoutputnewtracks2->Get("K0pointingAngleHypoBackgroundDetailed"));
    pointAngleHypoKTab[4]->SetName("K0pointingAngleHypoSignalDetailed2");
    pointAngleHypoKTab[5]->SetName("K0pointingAngleHypoBackgroundDetailed2");
    for(size_t i = 4; i<6; i++){
        pointAngleHypoKTab[i]->Scale(1/pointAngleHypoKTab[i]->GetBinContent(50));
        pointAngleHypoKTab[i]->SetLineWidth(2);
        if(i%2==1){
            pointAngleHypoKTab[i]->SetLineStyle(3);
        }
        if(i<2){
            pointAngleHypoKTab[i]->SetLineColor(kRed);
        } else if(i<4){
            pointAngleHypoKTab[i]->SetLineColor(kGreen);
        } else{
            pointAngleHypoKTab[i]->SetLineColor(kBlue);
        }
        pointAngleHypoKTab[i]->SetStats(0);
    }
    pointAngleHypoKTab[4]->SetTitle("pointingAngleHypo(), with normalization to the last bin");
    pointAngleHypoKTab[4]->Draw("HIST");
    for(size_t i = 5; i<6; i++){
        pointAngleHypoKTab[i]->Draw("HIST SAME");
    }
    TLegend *legendpointAngleHypoK = new TLegend(0.15, 0.55, 0.6, 0.85);
    legendpointAngleHypoK->SetTextSize(0.02);
    legendpointAngleHypoK->AddEntry("K0pointingAngleHypoSignalDetailed2", "K^{0}_{S} signal, new data, new tracks", "l");
    legendpointAngleHypoK->AddEntry("K0pointingAngleHypoBackgroundDetailed2", "K^{0}_{S} background, new data, new tracks", "l");
    legendpointAngleHypoK->SetBorderSize(1);
    legendpointAngleHypoK->Draw("SAME");
    pointAngleHypoKCanvas->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/K0pointAngleHypoAll_new_data.pdf").c_str());
    delete pointAngleHypoKCanvas;





    //K0 after final cuts
    //new data, new tracks
    TH1D *MKinvFinalCutClose1 = (TH1D *)anaoutputnewtracks1->Get("MpipiWideCheck");
    TH1D *MKinvFinalCutClose2 = (TH1D *)anaoutputnewtracks2->Get("MpipiWideCheck");
    TH1D resultKFinalCutHistNewDataNewTracks = *MKinvFinalCutClose1;
    resultKFinalCutHistNewDataNewTracks.Add(MKinvFinalCutClose2);
    //setting up the functions and data
    TCanvas *resultKFinalCut = new TCanvas("resultKFinalCut", "resultKFinalCut", 1800, 1600);
    resultKFinalCut->SetLeftMargin(0.15);
    TF1 *GfitKFinalCut = new TF1("GfitKFinalCut", "gausn(0) + pol1(3)", resultKFinalCutHistNewDataNewTracks.GetXaxis()->GetXmin(), resultKFinalCutHistNewDataNewTracks.GetXaxis()->GetXmax());
    TF1 *GfitKFinalCutBcg = new TF1("GfitKFinalCutBcg", "pol1", resultKFinalCutHistNewDataNewTracks.GetXaxis()->GetXmin(), resultKFinalCutHistNewDataNewTracks.GetXaxis()->GetXmax());
    TF1 *GfitKFinalCut1Sig = new TF1("GfitKFinalCut1Sig", "gausn", resultKFinalCutHistNewDataNewTracks.GetXaxis()->GetXmin(), resultKFinalCutHistNewDataNewTracks.GetXaxis()->GetXmax());
    GfitKFinalCut->SetParNames("Constant", "Mean", "Sigma", "b", "a");
    GfitKFinalCut->SetParameters(100, 0.5, 5e-3, 100, 0.48, 5e-3);
    // double nK0OldData, nK0NewDataOldTracks, nK0NewDataNewTracks;
    //new data new tracks
    resultKFinalCutHistNewDataNewTracks.SetMinimum(0);
    resultKFinalCutHistNewDataNewTracks.SetMarkerStyle(kFullCircle);
    resultKFinalCutHistNewDataNewTracks.SetMarkerColor(kBlue);
    resultKFinalCutHistNewDataNewTracks.Fit(GfitKFinalCut, "0");
    resultKFinalCutHistNewDataNewTracks.DrawCopy("E", "NewDataNewTracks");
    gPad->Update();
    GfitKFinalCut->GetParameters(paramsK);
    nK0NewDataNewTracks = paramsK[0];
    GfitKFinalCutBcg->SetParameters(paramsK+3);
    GfitKFinalCut1Sig->SetParameters(paramsK);
    GfitKFinalCutBcg->SetLineStyle(3);
    GfitKFinalCutBcg->SetLineWidth(3);
    GfitKFinalCut->SetLineColor(kBlue);
    GfitKFinalCutBcg->SetLineColor(kBlue);
    GfitKFinalCut1Sig->SetLineColor(kBlue);
    GfitKFinalCut->SetNpx(1000);
    GfitKFinalCut->DrawCopy("CSAME");
    GfitKFinalCutBcg->DrawCopy("CSAME");
    GfitKFinalCut1Sig->DrawCopy("CSAME");

    TLegend *legendKFinalCut = new TLegend(.18, .65, .45, .89);
    legendKFinalCut->SetTextSize(0.02);
    legendKFinalCut->AddEntry("MpipiWideCheckNewDataNewTracks", "pp new data, new tracks, #sqrt{s} = 510 GeV");
    legendKFinalCut->SetBorderSize(0);
    legendKFinalCut->Draw("SAME");

    nK0NewDataNewTracks /= resultKFinalCutHistNewDataNewTracks.GetBinWidth(1);
    std::stringstream nK0FinalCutNewDataNewTracksStream;
    nK0FinalCutNewDataNewTracksStream<<std::fixed<<std::setprecision(2)<<nK0NewDataNewTracks;
    // TPaveStats *pt = new TPaveStats(0.62, 0.7, 0.89, 0.89, "NB NDC");
    pt = new TPaveStats(0.62, 0.7, 0.89, 0.89, "NB NDC");
    // TText *tempText;
    pt->SetFillColor(kWhite);
    tempText = pt->AddText(0, 0, "Number of K^{0}_{S} detected in:");
    tempText->SetTextSize(3);
    tempText = pt->AddText(0, 1, ("New data = "+nK0FinalCutNewDataNewTracksStream.str()).c_str());
    tempText->SetTextSize(0.029);
    // pt->AddText("They are added to the pave using the AddText method.");
    // pt->AddLine(.0, .5, 1., .5);
    // pt->AddText("Even complex TLatex formulas can be added:");
    // TText *t1 = pt->AddText("F(t) = #sum_{i=-#infty}^{#infty}A(i)cos#[]{#frac{i}{t+i}}");
    // t1->SetTextColor(kBlue);
    pt->Draw("SAME");

    // gPad->RedrawAxis();
    resultKFinalCut->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/K0afterV0finderAfterCuts_new_data.pdf").c_str());













    ////////////////////////////////////////////////////////

    // LAMBDA

    ////////////////////////////////////////////////////////













    //Lambda0 backgrounds
    TCanvas *resultLambdaBkg = new TCanvas("resultLambdaBkg", "resultLambdaBkg", 1800, 1600);
    resultLambdaBkg->SetLeftMargin(0.15);
    //new data, new tracks
    TH1D *MLambdainvClose1Bkg = (TH1D *)anaoutputnewtracks1->Get("MppiWideMissing");
    TH1D *MLambdainvClose2Bkg = (TH1D *)anaoutputnewtracks2->Get("MppiWideMissing");
    TH1D resultLambdaHistNewDataNewTracksBkg = *MLambdainvClose1Bkg;
    resultLambdaHistNewDataNewTracks.Add(MLambdainvClose2Bkg);
    //new data new tracks
    resultLambdaHistNewDataNewTracksBkg.SetTitle("p^{#pm}#pi^{#mp} pairs rejected by basic cuts");
    resultLambdaHistNewDataNewTracksBkg.SetMinimum(0);
    resultLambdaHistNewDataNewTracksBkg.SetMarkerStyle(kFullCircle);
    resultLambdaHistNewDataNewTracksBkg.SetMarkerColor(kBlue);
    resultLambdaHistNewDataNewTracksBkg.DrawCopy("E", "NewDataNewTracksBkg");

    TLegend *legendLambdaBkg = new TLegend(.18, .65, .45, .89);
    legendLambdaBkg->SetTextSize(0.02);
    legendLambdaBkg->AddEntry("MppiWideMissingNewDataNewTracksBkg", "pp new data, #sqrt{s} = 510 GeV");
    legendLambdaBkg->SetBorderSize(0);
    legendLambdaBkg->Draw("SAME");

    resultLambdaBkg->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/Lambda0afterV0finderBkg_new_data.pdf").c_str());


    //Lambda0 after PID
    //new data, new tracks
    TH1D *MLambdainvPIDClose1 = (TH1D *)anaoutputnewtracks1->Get("MppiWideWithPidEcut");
    TH1D *MLambdainvPIDClose2 = (TH1D *)anaoutputnewtracks2->Get("MppiWideWithPidEcut");
    TH1D resultLambdaPIDHistNewDataNewTracks = *MLambdainvPIDClose1;
    resultLambdaPIDHistNewDataNewTracks.Add(MLambdainvPIDClose2);
    //setting up the functions and data
    TCanvas *resultLambdaPID = new TCanvas("resultLambdaPID", "resultLambdaPID", 1800, 1600);
    resultLambdaPID->SetLeftMargin(0.15);
    TF1 *GfitLambdaPID = new TF1("GfitLambdaPID", "gausn(0) + pol2(3)", resultLambdaPIDHistNewDataNewTracks.GetXaxis()->GetXmin(), resultLambdaPIDHistNewDataNewTracks.GetXaxis()->GetXmax());
    TF1 *GfitLambdaPIDBcg = new TF1("GfitLambdaPIDBcg", "pol2", resultLambdaPIDHistNewDataNewTracks.GetXaxis()->GetXmin(), resultLambdaPIDHistNewDataNewTracks.GetXaxis()->GetXmax());
    TF1 *GfitLambdaPID1Sig = new TF1("GfitLambdaPID1Sig", "gausn", resultLambdaPIDHistNewDataNewTracks.GetXaxis()->GetXmin(), resultLambdaPIDHistNewDataNewTracks.GetXaxis()->GetXmax());
    GfitLambdaPID->SetParNames("Constant", "Mean", "Sigma", "c", "b", "a");
    GfitLambdaPID->SetParameters(100, 1.115, 5e-3, -35000, 82000, -48000);
    GfitLambdaPID->SetParLimits(0, 0, 200);
    GfitLambdaPID->SetParLimits(1, 1.11, 1.12);
    GfitLambdaPID->SetParLimits(2, 0, 0.01);
    //new data new tracks
    resultLambdaPIDHistNewDataNewTracks.SetTitle("#Lambda^{0} mass in wide range with PID cuts");
    resultLambdaPIDHistNewDataNewTracks.SetMinimum(0);
    resultLambdaPIDHistNewDataNewTracks.SetMarkerStyle(kFullCircle);
    resultLambdaPIDHistNewDataNewTracks.SetMarkerColor(kBlue);
    resultLambdaPIDHistNewDataNewTracks.Fit(GfitLambdaPID, "0BR");
    resultLambdaPIDHistNewDataNewTracks.DrawCopy("E", "NewDataNewTracks");
    gPad->Update();
    GfitLambdaPID->GetParameters(paramsLambda);
    nLambda0NewDataNewTracks = paramsLambda[0];
    GfitLambdaPIDBcg->SetParameters(paramsLambda+3);
    GfitLambdaPID1Sig->SetParameters(paramsLambda);
    GfitLambdaPIDBcg->SetLineStyle(3);
    GfitLambdaPIDBcg->SetLineWidth(3);
    GfitLambdaPID->SetLineColor(kBlue);
    GfitLambdaPIDBcg->SetLineColor(kBlue);
    GfitLambdaPID1Sig->SetLineColor(kBlue);
    GfitLambdaPID->SetNpx(1000);
    GfitLambdaPID->DrawCopy("CSAME");
    GfitLambdaPIDBcg->DrawCopy("CSAME");
    GfitLambdaPID1Sig->DrawCopy("CSAME");

    TLegend *legendLambdaPID = new TLegend(.18, .74, .45, .84);
    legendLambdaPID->SetTextSize(0.02);
    legendLambdaPID->SetHeader("#bf{pp, #sqrt{s} = 510 GeV}", "C");
    legendLambdaPID->AddEntry("MppiWideWithPidEcutNewDataNewTracks", "#bf{new data,  new tracks}");
    legendLambdaPID->SetBorderSize(0);
    legendLambdaPID->Draw("SAME");

    nLambda0NewDataNewTracks /= resultLambdaPIDHistNewDataNewTracks.GetBinWidth(1);
    std::stringstream nLambda0PIDNewDataNewTracksStream;
    nLambda0PIDNewDataNewTracksStream<<std::fixed<<std::setprecision(2)<<nLambda0NewDataNewTracks;
    // TPaveStats *pt = new TPaveStats(0.62, 0.7, 0.89, 0.89, "NB NDC");
    pt = new TPaveStats(0.62, 0.7, 0.89, 0.89, "NB NDC");
    // TText *tempText;
    pt->SetFillColor(kWhite);
    tempText = pt->AddText(0, 0, "Number of #Lambda^{0} detected in:");
    tempText->SetTextSize(0.02);
    tempText = pt->AddText(0, 1, ("New data = "+nLambda0PIDNewDataNewTracksStream.str()).c_str());
    tempText->SetTextSize(0.029);
    // pt->AddText("They are added to the pave using the AddText method.");
    // pt->AddLine(.0, .5, 1., .5);
    // pt->AddText("Even complex TLatex formulas can be added:");
    // TText *t1 = pt->AddText("F(t) = #sum_{i=-#infty}^{#infty}A(i)cos#[]{#frac{i}{t+i}}");
    // t1->SetTextColor(kBlue);
    pt->Draw("SAME");

    // gPad->RedrawAxis();
    resultLambdaPID->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/Lambda0afterV0finderPID_new_data.pdf").c_str());





    //Lambda0 pointAngleHypo
    TCanvas *pointAngleHypoLambdaCanvas = new TCanvas("pointAngleHypoLambdaCanvas", "pointAngleHypoLambdaCanvas", 1800, 1600);
    pointAngleHypoLambdaCanvas->SetGrid();
    // gPad->SetLogy();
    TH1D *pointAngleHypoLambdaTab[6];
    pointAngleHypoLambdaTab[4] = (TH1D *)anaoutputnewtracks1->Get("LambdapointingAngleHypoSignalDetailed");
    pointAngleHypoLambdaTab[4]->Add((TH1D *)anaoutputnewtracks2->Get("LambdapointingAngleHypoSignalDetailed"));
    pointAngleHypoLambdaTab[5] = (TH1D *)anaoutputnewtracks1->Get("LambdapointingAngleHypoBackgroundDetailed");
    pointAngleHypoLambdaTab[5]->Add((TH1D *)anaoutputnewtracks2->Get("LambdapointingAngleHypoBackgroundDetailed"));
    pointAngleHypoLambdaTab[4]->SetName("Lambda0pointingAngleHypoSignalDetailed2");
    pointAngleHypoLambdaTab[5]->SetName("Lambda0pointingAngleHypoBackgroundDetailed2");
    for(size_t i = 4; i<6; i++){
        pointAngleHypoLambdaTab[i]->Scale(1/pointAngleHypoLambdaTab[i]->GetBinContent(50));
        pointAngleHypoLambdaTab[i]->SetLineWidth(2);
        if(i%2==1){
            pointAngleHypoLambdaTab[i]->SetLineStyle(3);
        }
        if(i<2){
            pointAngleHypoLambdaTab[i]->SetLineColor(kRed);
        } else if(i<4){
            pointAngleHypoLambdaTab[i]->SetLineColor(kGreen);
        } else{
            pointAngleHypoLambdaTab[i]->SetLineColor(kBlue);
        }
        pointAngleHypoLambdaTab[i]->SetStats(0);
    }
    pointAngleHypoLambdaTab[4]->SetTitle("pointingAngleHypo(), with normalization to the last bin");
    pointAngleHypoLambdaTab[4]->Draw("HIST");
    for(size_t i = 5; i<6; i++){
        pointAngleHypoLambdaTab[i]->Draw("HIST SAME");
    }
    TLegend *legendpointAngleHypoLambda = new TLegend(0.15, 0.55, 0.6, 0.85);
    legendpointAngleHypoLambda->SetTextSize(0.02);
    legendpointAngleHypoLambda->AddEntry("Lambda0pointingAngleHypoSignalDetailed2", "#Lambda^{0} signal, new data, new tracks", "l");
    legendpointAngleHypoLambda->AddEntry("Lambda0pointingAngleHypoBackgroundDetailed2", "#Lambda^{0} background, new data, new tracks", "l");
    legendpointAngleHypoLambda->SetBorderSize(1);
    legendpointAngleHypoLambda->Draw("SAME");
    pointAngleHypoLambdaCanvas->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/Lambda0pointAngleHypoAll_new_data.pdf").c_str());
    delete pointAngleHypoLambdaCanvas;





    //Lambda0 after final cuts
    //new data, new tracks
    TH1D *MLambdainvFinalCutClose1 = (TH1D *)anaoutputnewtracks1->Get("MppiWideCheck");
    TH1D *MLambdainvFinalCutClose2 = (TH1D *)anaoutputnewtracks2->Get("MppiWideCheck");
    TH1D resultLambdaFinalCutHistNewDataNewTracks = *MLambdainvFinalCutClose1;
    resultLambdaFinalCutHistNewDataNewTracks.Add(MLambdainvFinalCutClose2);
    //setting up the functions and data
    TCanvas *resultLambdaFinalCut = new TCanvas("resultLambdaFinalCut", "resultLambdaFinalCut", 1800, 1600);
    resultLambdaFinalCut->SetLeftMargin(0.15);
    TF1 *GfitLambdaFinalCut = new TF1("GfitLambdaFinalCut", "gausn(0) + pol2(3)", resultLambdaFinalCutHistNewDataNewTracks.GetXaxis()->GetXmin(), resultLambdaFinalCutHistNewDataNewTracks.GetXaxis()->GetXmax());
    TF1 *GfitLambdaFinalCutBcg = new TF1("GfitLambdaFinalCutBcg", "pol2", resultLambdaFinalCutHistNewDataNewTracks.GetXaxis()->GetXmin(), resultLambdaFinalCutHistNewDataNewTracks.GetXaxis()->GetXmax());
    TF1 *GfitLambdaFinalCut1Sig = new TF1("GfitLambdaFinalCut1Sig", "gausn", resultLambdaFinalCutHistNewDataNewTracks.GetXaxis()->GetXmin(), resultLambdaFinalCutHistNewDataNewTracks.GetXaxis()->GetXmax());
    GfitLambdaFinalCut->SetParNames("Constant", "Mean", "Sigma", "c", "b", "a");
    // GfitLambdaFinalCut->SetParameters(100, 1.115, 5e-3, 2000000, -4200000, 1500000);
    GfitLambdaFinalCut->SetParameters(100, 1.115, 5e-3, -33000, 58000, -25000);
    GfitLambdaFinalCut->SetParLimits(0, 0, 200);
    GfitLambdaFinalCut->SetParLimits(1, 1.11, 1.12);
    GfitLambdaFinalCut->SetParLimits(2, 0, 0.01);
    //new data new tracks
    resultLambdaFinalCutHistNewDataNewTracks.SetMinimum(0);
    resultLambdaFinalCutHistNewDataNewTracks.SetMarkerStyle(kFullCircle);
    resultLambdaFinalCutHistNewDataNewTracks.SetMarkerColor(kBlue);
    resultLambdaFinalCutHistNewDataNewTracks.Fit(GfitLambdaFinalCut, "0BR");
    resultLambdaFinalCutHistNewDataNewTracks.DrawCopy("E", "NewDataNewTracks");
    gPad->Update();
    GfitLambdaFinalCut->GetParameters(paramsLambda);
    nLambda0NewDataNewTracks = paramsLambda[0];
    GfitLambdaFinalCutBcg->SetParameters(paramsLambda+3);
    GfitLambdaFinalCut1Sig->SetParameters(paramsLambda);
    GfitLambdaFinalCutBcg->SetLineStyle(3);
    GfitLambdaFinalCutBcg->SetLineWidth(3);
    GfitLambdaFinalCut->SetLineColor(kBlue);
    GfitLambdaFinalCutBcg->SetLineColor(kBlue);
    GfitLambdaFinalCut1Sig->SetLineColor(kBlue);
    GfitLambdaFinalCut->SetNpx(1000);
    GfitLambdaFinalCut->DrawCopy("CSAME");
    GfitLambdaFinalCutBcg->DrawCopy("CSAME");
    GfitLambdaFinalCut1Sig->DrawCopy("CSAME");

    TLegend *legendLambdaFinalCut = new TLegend(.18, .65, .45, .89);
    legendLambdaFinalCut->SetTextSize(0.02);
    legendLambdaFinalCut->AddEntry("MppiWideCheckNewDataNewTracks", "pp new data, new tracks, #sqrt{s} = 510 GeV");
    legendLambdaFinalCut->SetBorderSize(0);
    legendLambdaFinalCut->Draw("SAME");

    nLambda0NewDataNewTracks /= resultLambdaFinalCutHistNewDataNewTracks.GetBinWidth(1);
    std::stringstream nLambda0FinalCutNewDataNewTracksStream;
    nLambda0FinalCutNewDataNewTracksStream<<std::fixed<<std::setprecision(2)<<nLambda0NewDataNewTracks;
    // TPaveStats *pt = new TPaveStats(0.62, 0.7, 0.89, 0.89, "NB NDC");
    pt = new TPaveStats(0.62, 0.7, 0.89, 0.89, "NB NDC");
    // TText *tempText;
    pt->SetFillColor(kWhite);
    tempText = pt->AddText(0, 0, "Number of #Lambda^{0} detected in:");
    tempText->SetTextSize(3);
    tempText = pt->AddText(0, 1, ("New data = "+nLambda0FinalCutNewDataNewTracksStream.str()).c_str());
    tempText->SetTextSize(0.029);
    // pt->AddText("They are added to the pave using the AddText method.");
    // pt->AddLine(.0, .5, 1., .5);
    // pt->AddText("Even complex TLatex formulas can be added:");
    // TText *t1 = pt->AddText("F(t) = #sum_{i=-#infty}^{#infty}A(i)cos#[]{#frac{i}{t+i}}");
    // t1->SetTextColor(kBlue);
    pt->Draw("SAME");

    // gPad->RedrawAxis();
    resultLambdaFinalCut->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/Lambda0afterV0finderAfterCuts_new_data.pdf").c_str());



    //dEdx
    TH2D *dEdx1 = (TH2D *)anaoutputnewtracks1->Get("dEdxMomentumEnergyLoss");
    TH2D *dEdx2 = (TH2D *)anaoutputnewtracks2->Get("dEdxMomentumEnergyLoss");
    TH2D dEdx = *dEdx1;
    dEdx.Add(dEdx2);
    dEdx.SetTitle(";pq [GeV/c];dE/dx [keV/cm]");
    c1 = new TCanvas("c1", "c1", 1600, 1200);
    c1->SetLogz();
    c1->SetMargin(0.10, 0.13, 0.1, 0.05);
    dEdx.Draw("colz");
    c1->Update();
    TLatex text1(-0.1, 2.5, "#bf{#pi^{#pm}}");
    text1.DrawClone("same");
    TLatex text2(-0.15, 11, "K^{#pm}");
    text2.DrawClone("same");
    TLatex text3(-0.1, 22, "p^{#pm}");
    text3.DrawClone("same");
    TLatex text4(0.9, 18, "deuterons");
    text4.DrawClone("same");
    c1->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/master_thesis/dEdxData.pdf").c_str());
    delete c1;

    //2d plot
    TH2D *Lambdapt2DHist1 = (TH2D *)anaoutputnewtracks1->Get("Lambdapt2DHist");
    TH2D *Lambdapt2DHist2 = (TH2D *)anaoutputnewtracks2->Get("Lambdapt2DHist");
    TH2D Lambdapt2DHist = *Lambdapt2DHist1;
    Lambdapt2DHist.Add(Lambdapt2DHist2);
    Lambdapt2DHist.SetTitle(";m_{p#pi} [GeV];p_T [GeV/c]");
    c1 = new TCanvas("c1", "c1", 1600, 1000);
    c1->SetMargin(0.10, 0.13, 0.1, 0.05);
    Lambdapt2DHist.Draw("colz");
    c1->Update();
    TEllipse newEllipse(1.115, 1.1, 0.005, 0.6);
    newEllipse.SetLineWidth(2);
    newEllipse.SetLineStyle(kDashed);
    newEllipse.SetLineColor(kRed);
    newEllipse.SetFillColorAlpha(0, 0);
    newEllipse.Draw("same");
    c1->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/master_thesis/Lambdapt2DHist.pdf").c_str());
    delete c1;

    //1d example plot
    gStyle->SetOptFit(1);
    int n = 5;
    TH1D *Lambdapt1DHist = Lambdapt2DHist.ProjectionX("Lambdapt1DHist", n, n);
    c1 = new TCanvas("c1", "c1", 1600, 1200);
    c1->SetMargin(0.05, 0.05, 0.1, 0.05);
    //Lambda
    GfitLambda->SetParNames("Constant", "Mean", "Sigma", "c", "b", "a");
    GfitLambda->SetParameters(3, 1.115, 5e-3, 0, 0., -10000000);
    GfitLambda->SetParLimits(0, 0, 200);
    GfitLambda->SetParLimits(1, 1.11, 1.12);
    GfitLambda->SetParLimits(2, 5e-5, 0.01);
    //old data, the base of drawing
    Lambdapt1DHist->SetMinimum(0);
    Lambdapt1DHist->SetMarkerStyle(kFullCircle);
    Lambdapt1DHist->SetMarkerColor(kBlue);
    Lambdapt1DHist->SetLineColor(kBlue+2);
    Lambdapt1DHist->Fit(GfitLambda, "0BR");
    Lambdapt1DHist->Fit(GfitLambda, "0BR");
    Lambdapt1DHist->Fit(GfitLambda, "0BR");
    Lambdapt1DHist->Fit(GfitLambda, "0BR");
    Lambdapt1DHist->Fit(GfitLambda, "0BR");
    Lambdapt1DHist->DrawCopy("E", "NewDataNewTracks");
    gPad->Update();
    GfitLambda->GetParameters(paramsLambda);
    GfitLambdaBcg->SetParameters(paramsLambda+3);
    GfitLambda1Sig->SetParameters(paramsLambda);
    GfitLambdaBcg->SetLineStyle(3);
    GfitLambdaBcg->SetLineWidth(3);
    GfitLambda1Sig->SetLineColor(kBlue);
    GfitLambda->SetNpx(1000);
    GfitLambda->Draw("CSAME");
    GfitLambdaBcg->Draw("CSAME");
    GfitLambda1Sig->Draw("CSAME");

    TPaveStats *ps = (TPaveStats *)c1->GetPrimitive("stats");
    ps->SetX1(1.065);
    ps->SetX2(1.1);
    ps->SetY1(600);
    ps->SetY2(900);
    ps->SetBorderSize(0);

    legendLambda = new TLegend(.1, .35, .45, .59);
    legendLambda->SetTextSize(0.035);
    legendLambda->AddEntry("Lambdapt1DHistNewDataNewTracks", "#Lambda^{0}, 0.8 GeV/c < p_{T} < 1.0 GeV/c");
    legendLambda->SetBorderSize(0);
    legendLambda->Draw("SAME");
    gPad->Update();
    c1->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/master_thesis/Lambdapt1DHist.pdf").c_str());
    delete c1;

    return 0;
}