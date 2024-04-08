#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPaveText.h"

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
    TFile *anaoutputolddata1 = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/CPT/AnaOutput_Inclusive_analysis_with_STUPCV0_with_extended_range_noAfterburner.root");
    TFile *anaoutputolddata2 = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/CPTnoBBCL/AnaOutput_Inclusive_analysis_with_STUPCV0_with_extended_range_noAfterburner.root");
    TFile *anaoutputoldtracks1 = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/CPT/AnaOutput_Inclusive_analysis_with_STUPCV0_with_extended_range_noAfterburner_new_data_old_tracks.root");
    TFile *anaoutputoldtracks2 = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/CPTnoBBCL/AnaOutput_Inclusive_analysis_with_STUPCV0_with_extended_range_noAfterburner_new_data_old_tracks.root");
    TFile *anaoutputnewtracks1 = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/CPT/AnaOutput_Inclusive_analysis_with_STUPCV0_with_extended_range_noAfterburner_new_data_new_tracks.root");
    TFile *anaoutputnewtracks2 = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/CPTnoBBCL/AnaOutput_Inclusive_analysis_with_STUPCV0_with_extended_range_noAfterburner_new_data_new_tracks.root");

    gStyle->SetFrameLineWidth(3);
    // gStyle->SetOptFit(111);
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);

    //old data
    TH1D *MKinvClose1 = (TH1D *)anaoutputolddata1->Get("MpipiWide");
    TH1D *MKinvClose2 = (TH1D *)anaoutputolddata2->Get("MpipiWide");
    TH1D resultKHistOldData = *MKinvClose1;
    resultKHistOldData.Add(MKinvClose2);
    // resultKHistOldData.SetStats(1);
    //new data, old tracks
    MKinvClose1 = (TH1D *)anaoutputoldtracks1->Get("MpipiWide");
    MKinvClose2 = (TH1D *)anaoutputoldtracks2->Get("MpipiWide");
    TH1D resultKHistNewDataOldTracks = *MKinvClose1;
    resultKHistNewDataOldTracks.Add(MKinvClose2);
    // resultKHistNewDataOldTracks.SetStats(1);
    //new data, new tracks
    MKinvClose1 = (TH1D *)anaoutputnewtracks1->Get("MpipiWide");
    MKinvClose2 = (TH1D *)anaoutputnewtracks2->Get("MpipiWide");
    TH1D resultKHistNewDataNewTracks = *MKinvClose1;
    resultKHistNewDataNewTracks.Add(MKinvClose2);
    // resultKHistNewDataNewTracks.SetStats(1);

    //K0
    //setting up the functions and data
    TCanvas *resultK = new TCanvas("resultK", "resultK", 1800, 1600);
    resultK->SetLeftMargin(0.15);
    TF1 *GfitK = new TF1("GfitK", "gausn(0) + pol1(3)", resultKHistOldData.GetXaxis()->GetXmin(), resultKHistOldData.GetXaxis()->GetXmax());
    TF1 *GfitKBcg = new TF1("GfitKBcg", "pol1", resultKHistOldData.GetXaxis()->GetXmin(), resultKHistOldData.GetXaxis()->GetXmax());
    TF1 *GfitK1Sig = new TF1("GfitK1Sig", "gausn", resultKHistOldData.GetXaxis()->GetXmin(), resultKHistOldData.GetXaxis()->GetXmax());
    GfitK->SetParNames("Constant", "Mean", "Sigma", "b", "a");
    GfitK->SetParameters(100, 0.5, 5e-3, 100, 0.48, 5e-3);
    double nK0OldData, nK0NewDataOldTracks, nK0NewDataNewTracks;
    //old data, the base of drawing
    resultKHistOldData.SetMinimum(0);
    resultKHistOldData.SetMaximum(max(resultKHistOldData.GetMaximum(), max(resultKHistNewDataOldTracks.GetMaximum(), resultKHistNewDataNewTracks.GetMaximum()))*1.1);
    resultKHistOldData.SetMarkerStyle(kFullCircle);
    resultKHistOldData.SetMarkerColor(kRed);
    resultKHistOldData.SetLineColor(kRed+2);
    resultKHistOldData.Fit(GfitK, "0");
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
    nK0OldData = paramsK[0];
    // cout<<"PARAMETERS:"<<endl;
    // for(int i = 0; i<5; i++){
    //     cout<<paramsK[i]<<endl;
    // }
    GfitKBcg->SetParameters(paramsK+3);
    GfitK1Sig->SetParameters(paramsK);
    GfitKBcg->SetLineStyle(3);
    GfitKBcg->SetLineWidth(3);
    GfitK1Sig->SetLineColor(kRed);
    GfitK->SetNpx(1000);
    GfitK->DrawCopy("CSAME");
    GfitKBcg->DrawCopy("CSAME");
    GfitK1Sig->DrawCopy("CSAME");
    //new data old tracks
    resultKHistNewDataOldTracks.SetMinimum(0);
    resultKHistNewDataOldTracks.SetMarkerStyle(kFullCircle);
    resultKHistNewDataOldTracks.SetMarkerColor(kGreen);
    resultKHistNewDataOldTracks.SetLineColor(kGreen+2);
    resultKHistNewDataOldTracks.Fit(GfitK, "0");
    resultKHistNewDataOldTracks.DrawCopy("E SAME", "NewDataOldTracks");
    gPad->Update();
    GfitK->GetParameters(paramsK);
    nK0NewDataOldTracks = paramsK[0];
    GfitKBcg->SetParameters(paramsK+3);
    GfitK1Sig->SetParameters(paramsK);
    GfitKBcg->SetLineStyle(3);
    GfitKBcg->SetLineWidth(3);
    GfitK->SetLineColor(kGreen);
    GfitKBcg->SetLineColor(kGreen);
    GfitK1Sig->SetLineColor(kGreen);
    GfitK->SetNpx(1000);
    GfitK->DrawCopy("CSAME");
    GfitKBcg->DrawCopy("CSAME");
    GfitK1Sig->DrawCopy("CSAME");
    //new data new tracks
    resultKHistNewDataNewTracks.SetMinimum(0);
    resultKHistNewDataNewTracks.SetMarkerStyle(kFullCircle);
    resultKHistNewDataNewTracks.SetMarkerColor(kBlue);
    resultKHistNewDataNewTracks.Fit(GfitK, "0");
    resultKHistNewDataNewTracks.DrawCopy("E SAME", "NewDataNewTracks");
    gPad->Update();
    GfitK->GetParameters(paramsK);
    nK0NewDataNewTracks = paramsK[0];
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
    legendK->AddEntry("MpipiWideOldData", "pp old data, #sqrt{s} = 510 GeV");
    legendK->AddEntry("MpipiWideNewDataOldTracks", "pp new data, old tracks, #sqrt{s} = 510 GeV");
    legendK->AddEntry("MpipiWideNewDataNewTracks", "pp new data, new tracks, #sqrt{s} = 510 GeV");
    legendK->SetBorderSize(0);
    legendK->Draw("SAME");

    nK0OldData /= resultKHistOldData.GetBinWidth(1);
    nK0NewDataOldTracks /= resultKHistOldData.GetBinWidth(1);
    nK0NewDataNewTracks /= resultKHistOldData.GetBinWidth(1);
    std::stringstream nK0OldDataStream, nK0NewDataOldTracksStream, nK0NewDataNewTracksStream;
    nK0OldDataStream<<std::fixed<<std::setprecision(2)<<nK0OldData;
    nK0NewDataOldTracksStream<<std::fixed<<std::setprecision(2)<<nK0NewDataOldTracks;
    nK0NewDataNewTracksStream<<std::fixed<<std::setprecision(2)<<nK0NewDataNewTracks;
    TPaveStats *pt = new TPaveStats(0.62, 0.7, 0.89, 0.89, "NB NDC");
    TText *tempText;
    pt->SetFillColor(kWhite);
    tempText = pt->AddText(0, 0, "Number of K^{0}_{S} detected in:");
    tempText->SetTextSize(3);
    pt->AddText(0, 0.33, ("Old data = "+nK0OldDataStream.str()).c_str());
    pt->AddText(0, 0.66, ("New data, old tracks = "+nK0NewDataOldTracksStream.str()).c_str());
    pt->AddText(0, 1, ("New data, new tracks = "+nK0NewDataNewTracksStream.str()).c_str());
    // pt->AddText("They are added to the pave using the AddText method.");
    // pt->AddLine(.0, .5, 1., .5);
    // pt->AddText("Even complex TLatex formulas can be added:");
    // TText *t1 = pt->AddText("F(t) = #sum_{i=-#infty}^{#infty}A(i)cos#[]{#frac{i}{t+i}}");
    // t1->SetTextColor(kBlue);
    pt->Draw("SAME");

    // gPad->RedrawAxis();
    resultK->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/K0afterV0finder.pdf").c_str());



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
    GfitLambda->SetParameters(100, 1.115, 5e-3, 2000000, -4200000, 1500000);
    GfitLambda->SetParLimits(0, 0, 200);
    GfitLambda->SetParLimits(1, 1.11, 1.12);
    GfitLambda->SetParLimits(2, 0, 0.01);
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
    GfitLambda->SetParameters(100, 1.115, 5e-3, 2000000, -4860000, 1500000);
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
    GfitLambda->SetParameters(100, 1.115, 5e-3, 2000000, -4860000, 1500000);
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
    pt = new TPaveStats(.18, .45, .45, .65, "NB NDC");
    pt->SetFillColor(kWhite);
    tempText = pt->AddText(0, 0, "Number of #Lambda^{0} detected in:");
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
    resultLambda->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/Lambda0afterV0finder.pdf").c_str());



    //rest of the histograms



    //NotherTracks
    TCanvas *c1 = new TCanvas("c1", "c1", 1800, 1600);
    c1->SetGrid();
    gPad->SetLogy();
    TH1D *NotherTracksTab[6];
    NotherTracksTab[0] = (TH1D *)anaoutputolddata1->Get("NotherTracks");
    NotherTracksTab[1] = (TH1D *)anaoutputolddata2->Get("NotherTracks");
    NotherTracksTab[2] = (TH1D *)anaoutputoldtracks1->Get("NotherTracks");
    NotherTracksTab[3] = (TH1D *)anaoutputoldtracks2->Get("NotherTracks");
    NotherTracksTab[4] = (TH1D *)anaoutputnewtracks1->Get("NotherTracks");
    NotherTracksTab[5] = (TH1D *)anaoutputnewtracks2->Get("NotherTracks");
    NotherTracksTab[0]->SetName("NotherTracks0");
    NotherTracksTab[1]->SetName("NotherTracks1");
    NotherTracksTab[2]->SetName("NotherTracks2");
    NotherTracksTab[3]->SetName("NotherTracks3");
    NotherTracksTab[4]->SetName("NotherTracks4");
    NotherTracksTab[5]->SetName("NotherTracks5");
    for(size_t i = 0; i<6; i++){
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
    NotherTracksTab[0]->Draw("HIST");
    for(size_t i = 1; i<6; i++){
        NotherTracksTab[i]->Draw("HIST SAME");
    }
    TLegend *legendNotherTracks = new TLegend(0.5, 0.6, 0.89, 0.89);
    legendNotherTracks->SetTextSize(0.02);
    legendNotherTracks->AddEntry("NotherTracks0", "CPT2, old data", "l");
    legendNotherTracks->AddEntry("NotherTracks1", "CPT2noBBCL, old data", "l");
    legendNotherTracks->AddEntry("NotherTracks2", "CPT2, new data, old tracks", "l");
    legendNotherTracks->AddEntry("NotherTracks3", "CPT2noBBCL, new data, old tracks", "l");
    legendNotherTracks->AddEntry("NotherTracks4", "CPT2, new data, new tracks", "l");
    legendNotherTracks->AddEntry("NotherTracks5", "CPT2noBBCL, new data, new tracks", "l");
    legendNotherTracks->SetBorderSize(1);
    legendNotherTracks->Draw("SAME");
    c1->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/NotherTracksAll.pdf").c_str());
    delete c1;

    //K0multiplicity
    c1 = new TCanvas("c1", "c1", 1800, 1600);
    c1->SetGrid();
    gPad->SetLogy();
    TH1D *K0multiplicityTab[6];
    K0multiplicityTab[0] = (TH1D *)anaoutputolddata1->Get("K0multiplicity");
    K0multiplicityTab[1] = (TH1D *)anaoutputolddata2->Get("K0multiplicity");
    K0multiplicityTab[2] = (TH1D *)anaoutputoldtracks1->Get("K0multiplicity");
    K0multiplicityTab[3] = (TH1D *)anaoutputoldtracks2->Get("K0multiplicity");
    K0multiplicityTab[4] = (TH1D *)anaoutputnewtracks1->Get("K0multiplicity");
    K0multiplicityTab[5] = (TH1D *)anaoutputnewtracks2->Get("K0multiplicity");
    K0multiplicityTab[0]->SetName("K0multiplicity0");
    K0multiplicityTab[1]->SetName("K0multiplicity1");
    K0multiplicityTab[2]->SetName("K0multiplicity2");
    K0multiplicityTab[3]->SetName("K0multiplicity3");
    K0multiplicityTab[4]->SetName("K0multiplicity4");
    K0multiplicityTab[5]->SetName("K0multiplicity5");
    for(size_t i = 0; i<6; i++){
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
    K0multiplicityTab[0]->Draw("HIST");
    for(size_t i = 1; i<6; i++){
        K0multiplicityTab[i]->Draw("HIST SAME");
    }
    TLegend *legendK0multiplicity = new TLegend(0.5, 0.6, 0.89, 0.89);
    legendK0multiplicity->SetTextSize(0.02);
    legendK0multiplicity->AddEntry("K0multiplicity0", "CPT2, old data", "l");
    legendK0multiplicity->AddEntry("K0multiplicity1", "CPT2noBBCL, old data", "l");
    legendK0multiplicity->AddEntry("K0multiplicity2", "CPT2, new data, old tracks", "l");
    legendK0multiplicity->AddEntry("K0multiplicity3", "CPT2noBBCL, new data, old tracks", "l");
    legendK0multiplicity->AddEntry("K0multiplicity4", "CPT2, new data, new tracks", "l");
    legendK0multiplicity->AddEntry("K0multiplicity5", "CPT2noBBCL, new data, new tracks", "l");
    legendK0multiplicity->SetBorderSize(1);
    legendK0multiplicity->Draw("SAME");
    c1->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/K0multiplicityAll.pdf").c_str());
    delete c1;

    return 0;
}