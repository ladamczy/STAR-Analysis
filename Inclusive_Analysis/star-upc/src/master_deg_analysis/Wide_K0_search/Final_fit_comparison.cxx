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

    gStyle->SetFrameLineWidth(1);
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





    //K0 backgrounds
    TCanvas *resultKBkg = new TCanvas("resultKBkg", "resultKBkg", 1800, 1600);
    resultKBkg->SetLeftMargin(0.15);
    //old data
    TH1D *MKinvClose1Bkg = (TH1D *)anaoutputolddata1->Get("MpipiWideMissing");
    TH1D *MKinvClose2Bkg = (TH1D *)anaoutputolddata2->Get("MpipiWideMissing");
    TH1D resultKHistOldDataBkg = *MKinvClose1Bkg;
    resultKHistOldData.Add(MKinvClose2Bkg);
    // resultKHistOldData.SetStats(1);
    //new data, old tracks
    MKinvClose1Bkg = (TH1D *)anaoutputoldtracks1->Get("MpipiWideMissing");
    MKinvClose2Bkg = (TH1D *)anaoutputoldtracks2->Get("MpipiWideMissing");
    TH1D resultKHistNewDataOldTracksBkg = *MKinvClose1Bkg;
    resultKHistNewDataOldTracks.Add(MKinvClose2Bkg);
    // resultKHistNewDataOldTracks.SetStats(1);
    //new data, new tracks
    MKinvClose1Bkg = (TH1D *)anaoutputnewtracks1->Get("MpipiWideMissing");
    MKinvClose2Bkg = (TH1D *)anaoutputnewtracks2->Get("MpipiWideMissing");
    TH1D resultKHistNewDataNewTracksBkg = *MKinvClose1Bkg;
    resultKHistNewDataNewTracks.Add(MKinvClose2Bkg);
    //old data, the base of drawing
    resultKHistOldDataBkg.SetTitle("#pi^{+}#pi^{-} pairs rejected by basic cuts");
    resultKHistOldDataBkg.SetMinimum(0);
    resultKHistOldDataBkg.SetMaximum(max(resultKHistOldDataBkg.GetMaximum(), max(resultKHistNewDataOldTracksBkg.GetMaximum(), resultKHistNewDataNewTracksBkg.GetMaximum()))*1.1);
    resultKHistOldDataBkg.SetMarkerStyle(kFullCircle);
    resultKHistOldDataBkg.SetMarkerColor(kRed);
    resultKHistOldDataBkg.SetLineColor(kRed+2);
    resultKHistOldDataBkg.DrawCopy("E", "OldDataBkg");
    //new data old tracks
    resultKHistNewDataOldTracksBkg.SetMinimum(0);
    resultKHistNewDataOldTracksBkg.SetMarkerStyle(kFullCircle);
    resultKHistNewDataOldTracksBkg.SetMarkerColor(kGreen);
    resultKHistNewDataOldTracksBkg.SetLineColor(kGreen+2);
    resultKHistNewDataOldTracksBkg.DrawCopy("E SAME", "NewDataOldTracksBkg");
    //new data new tracks
    resultKHistNewDataNewTracksBkg.SetMinimum(0);
    resultKHistNewDataNewTracksBkg.SetMarkerStyle(kFullCircle);
    resultKHistNewDataNewTracksBkg.SetMarkerColor(kBlue);
    resultKHistNewDataNewTracksBkg.DrawCopy("E SAME", "NewDataNewTracksBkg");
    resultKBkg->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/K0afterV0finderBkg.pdf").c_str());


    //K0 after PID
    //old data
    TH1D *MKinvPIDClose1 = (TH1D *)anaoutputolddata1->Get("MpipiWideWithPidEcut");
    TH1D *MKinvPIDClose2 = (TH1D *)anaoutputolddata2->Get("MpipiWideWithPidEcut");
    TH1D resultKPIDHistOldData = *MKinvPIDClose1;
    resultKPIDHistOldData.Add(MKinvPIDClose2);
    // resultKPIDHistOldData.SetStats(1);
    //new data, old tracks
    MKinvPIDClose1 = (TH1D *)anaoutputoldtracks1->Get("MpipiWideWithPidEcut");
    MKinvPIDClose2 = (TH1D *)anaoutputoldtracks2->Get("MpipiWideWithPidEcut");
    TH1D resultKPIDHistNewDataOldTracks = *MKinvPIDClose1;
    resultKPIDHistNewDataOldTracks.Add(MKinvPIDClose2);
    // resultKPIDHistNewDataOldTracks.SetStats(1);
    //new data, new tracks
    MKinvPIDClose1 = (TH1D *)anaoutputnewtracks1->Get("MpipiWideWithPidEcut");
    MKinvPIDClose2 = (TH1D *)anaoutputnewtracks2->Get("MpipiWideWithPidEcut");
    TH1D resultKPIDHistNewDataNewTracks = *MKinvPIDClose1;
    resultKPIDHistNewDataNewTracks.Add(MKinvPIDClose2);
    //setting up the functions and data
    TCanvas *resultKPID = new TCanvas("resultKPID", "resultKPID", 1800, 1600);
    resultKPID->SetLeftMargin(0.15);
    TF1 *GfitKPID = new TF1("GfitKPID", "gausn(0) + pol1(3)", resultKHistOldData.GetXaxis()->GetXmin(), resultKHistOldData.GetXaxis()->GetXmax());
    TF1 *GfitKPIDBcg = new TF1("GfitKPIDBcg", "pol1", resultKHistOldData.GetXaxis()->GetXmin(), resultKHistOldData.GetXaxis()->GetXmax());
    TF1 *GfitKPID1Sig = new TF1("GfitKPID1Sig", "gausn", resultKHistOldData.GetXaxis()->GetXmin(), resultKHistOldData.GetXaxis()->GetXmax());
    GfitKPID->SetParNames("Constant", "Mean", "Sigma", "b", "a");
    GfitKPID->SetParameters(100, 0.5, 5e-3, 100, 0.48, 5e-3);
    // double nK0OldData, nK0NewDataOldTracks, nK0NewDataNewTracks;
    //old data, the base of drawing
    resultKPIDHistOldData.SetTitle("K^{0}_{S} mass in wide range with PID cuts");
    resultKPIDHistOldData.SetMinimum(0);
    resultKPIDHistOldData.SetMaximum(max(resultKPIDHistOldData.GetMaximum(), max(resultKPIDHistNewDataOldTracks.GetMaximum(), resultKPIDHistNewDataNewTracks.GetMaximum()))*1.1);
    resultKPIDHistOldData.SetMarkerStyle(kFullCircle);
    resultKPIDHistOldData.SetMarkerColor(kRed);
    resultKPIDHistOldData.SetLineColor(kRed+2);
    resultKPIDHistOldData.Fit(GfitKPID, "0");
    // TH1D *tempHist = (TH1D *)resultKPIDHistOldData.DrawCopy("E", "OldData");
    resultKPIDHistOldData.DrawCopy("E", "OldData");
    gPad->Update();
    // TPaveStats *st = (TPaveStats *)tempHist->FindObject("stats");
    // st->SetX1NDC(0.15);
    // st->SetX2NDC(0.35);
    // st->SetY1NDC(0.65);
    // st->SetY2NDC(0.85);
    // st->SetBorderSize(0);
    // Double_t paramsK[5];
    GfitKPID->GetParameters(paramsK);
    nK0OldData = paramsK[0];
    // cout<<"PARAMETERS:"<<endl;
    // for(int i = 0; i<5; i++){
    //     cout<<paramsK[i]<<endl;
    // }
    GfitKPIDBcg->SetParameters(paramsK+3);
    GfitKPID1Sig->SetParameters(paramsK);
    GfitKPIDBcg->SetLineStyle(3);
    GfitKPIDBcg->SetLineWidth(3);
    GfitKPID1Sig->SetLineColor(kRed);
    GfitKPID->SetNpx(1000);
    GfitKPID->DrawCopy("CSAME");
    GfitKPIDBcg->DrawCopy("CSAME");
    GfitKPID1Sig->DrawCopy("CSAME");
    //new data old tracks
    resultKPIDHistNewDataOldTracks.SetMinimum(0);
    resultKPIDHistNewDataOldTracks.SetMarkerStyle(kFullCircle);
    resultKPIDHistNewDataOldTracks.SetMarkerColor(kGreen);
    resultKPIDHistNewDataOldTracks.SetLineColor(kGreen+2);
    resultKPIDHistNewDataOldTracks.Fit(GfitKPID, "0");
    resultKPIDHistNewDataOldTracks.DrawCopy("E SAME", "NewDataOldTracks");
    gPad->Update();
    GfitKPID->GetParameters(paramsK);
    nK0NewDataOldTracks = paramsK[0];
    GfitKPIDBcg->SetParameters(paramsK+3);
    GfitKPID1Sig->SetParameters(paramsK);
    GfitKPIDBcg->SetLineStyle(3);
    GfitKPIDBcg->SetLineWidth(3);
    GfitKPID->SetLineColor(kGreen);
    GfitKPIDBcg->SetLineColor(kGreen);
    GfitKPID1Sig->SetLineColor(kGreen);
    GfitKPID->SetNpx(1000);
    GfitKPID->DrawCopy("CSAME");
    GfitKPIDBcg->DrawCopy("CSAME");
    GfitKPID1Sig->DrawCopy("CSAME");
    //new data new tracks
    resultKPIDHistNewDataNewTracks.SetMinimum(0);
    resultKPIDHistNewDataNewTracks.SetMarkerStyle(kFullCircle);
    resultKPIDHistNewDataNewTracks.SetMarkerColor(kBlue);
    resultKPIDHistNewDataNewTracks.Fit(GfitKPID, "0");
    resultKPIDHistNewDataNewTracks.DrawCopy("E SAME", "NewDataNewTracks");
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

    TLegend *legendKPID = new TLegend(.18, .65, .45, .89);
    legendKPID->SetTextSize(0.02);
    legendKPID->AddEntry("MpipiWideWithPidEcutOldData", "pp old data, #sqrt{s} = 510 GeV");
    legendKPID->AddEntry("MpipiWideWithPidEcutNewDataOldTracks", "pp new data, old tracks, #sqrt{s} = 510 GeV");
    legendKPID->AddEntry("MpipiWideWithPidEcutNewDataNewTracks", "pp new data, new tracks, #sqrt{s} = 510 GeV");
    legendKPID->SetBorderSize(0);
    legendKPID->Draw("SAME");

    nK0OldData /= resultKPIDHistOldData.GetBinWidth(1);
    nK0NewDataOldTracks /= resultKPIDHistOldData.GetBinWidth(1);
    nK0NewDataNewTracks /= resultKPIDHistOldData.GetBinWidth(1);
    std::stringstream nK0PIDOldDataStream, nK0PIDNewDataOldTracksStream, nK0PIDNewDataNewTracksStream;
    nK0PIDOldDataStream<<std::fixed<<std::setprecision(2)<<nK0OldData;
    nK0PIDNewDataOldTracksStream<<std::fixed<<std::setprecision(2)<<nK0NewDataOldTracks;
    nK0PIDNewDataNewTracksStream<<std::fixed<<std::setprecision(2)<<nK0NewDataNewTracks;
    // TPaveStats *pt = new TPaveStats(0.62, 0.7, 0.89, 0.89, "NB NDC");
    pt = new TPaveStats(0.62, 0.7, 0.89, 0.89, "NB NDC");
    // TText *tempText;
    pt->SetFillColor(kWhite);
    tempText = pt->AddText(0, 0, "Number of K^{0}_{S} detected in:");
    tempText->SetTextSize(3);
    pt->AddText(0, 0.33, ("Old data = "+nK0PIDOldDataStream.str()).c_str());
    pt->AddText(0, 0.66, ("New data, old tracks = "+nK0PIDNewDataOldTracksStream.str()).c_str());
    pt->AddText(0, 1, ("New data, new tracks = "+nK0PIDNewDataNewTracksStream.str()).c_str());
    // pt->AddText("They are added to the pave using the AddText method.");
    // pt->AddLine(.0, .5, 1., .5);
    // pt->AddText("Even complex TLatex formulas can be added:");
    // TText *t1 = pt->AddText("F(t) = #sum_{i=-#infty}^{#infty}A(i)cos#[]{#frac{i}{t+i}}");
    // t1->SetTextColor(kBlue);
    pt->Draw("SAME");

    // gPad->RedrawAxis();
    resultKPID->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/K0afterV0finderPID.pdf").c_str());





    //K0 pointAngleHypo
    TCanvas *pointAngleHypoKCanvas = new TCanvas("pointAngleHypoKCanvas", "pointAngleHypoKCanvas", 1800, 1600);
    pointAngleHypoKCanvas->SetGrid();
    // gPad->SetLogy();
    TH1D *pointAngleHypoKTab[6];
    pointAngleHypoKTab[0] = (TH1D *)anaoutputolddata1->Get("K0pointingAngleHypoSignalDetailed");
    pointAngleHypoKTab[0]->Add((TH1D *)anaoutputolddata2->Get("K0pointingAngleHypoSignalDetailed"));
    pointAngleHypoKTab[1] = (TH1D *)anaoutputolddata1->Get("K0pointingAngleHypoBackgroundDetailed");
    pointAngleHypoKTab[1]->Add((TH1D *)anaoutputolddata2->Get("K0pointingAngleHypoBackgroundDetailed"));
    pointAngleHypoKTab[2] = (TH1D *)anaoutputoldtracks1->Get("K0pointingAngleHypoSignalDetailed");
    pointAngleHypoKTab[2]->Add((TH1D *)anaoutputoldtracks2->Get("K0pointingAngleHypoSignalDetailed"));
    pointAngleHypoKTab[3] = (TH1D *)anaoutputoldtracks1->Get("K0pointingAngleHypoBackgroundDetailed");
    pointAngleHypoKTab[3]->Add((TH1D *)anaoutputoldtracks2->Get("K0pointingAngleHypoBackgroundDetailed"));
    pointAngleHypoKTab[4] = (TH1D *)anaoutputnewtracks1->Get("K0pointingAngleHypoSignalDetailed");
    pointAngleHypoKTab[4]->Add((TH1D *)anaoutputnewtracks2->Get("K0pointingAngleHypoSignalDetailed"));
    pointAngleHypoKTab[5] = (TH1D *)anaoutputnewtracks1->Get("K0pointingAngleHypoBackgroundDetailed");
    pointAngleHypoKTab[5]->Add((TH1D *)anaoutputnewtracks2->Get("K0pointingAngleHypoBackgroundDetailed"));
    pointAngleHypoKTab[0]->SetName("K0pointingAngleHypoSignalDetailed0");
    pointAngleHypoKTab[1]->SetName("K0pointingAngleHypoBackgroundDetailed0");
    pointAngleHypoKTab[2]->SetName("K0pointingAngleHypoSignalDetailed1");
    pointAngleHypoKTab[3]->SetName("K0pointingAngleHypoBackgroundDetailed1");
    pointAngleHypoKTab[4]->SetName("K0pointingAngleHypoSignalDetailed2");
    pointAngleHypoKTab[5]->SetName("K0pointingAngleHypoBackgroundDetailed2");
    for(size_t i = 0; i<6; i++){
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
    pointAngleHypoKTab[0]->SetTitle("pointingAngleHypo(), with normalization to the last bin");
    pointAngleHypoKTab[0]->Draw("HIST");
    for(size_t i = 1; i<6; i++){
        pointAngleHypoKTab[i]->Draw("HIST SAME");
    }
    TLegend *legendpointAngleHypoK = new TLegend(0.15, 0.55, 0.6, 0.85);
    legendpointAngleHypoK->SetTextSize(0.02);
    legendpointAngleHypoK->AddEntry("K0pointingAngleHypoSignalDetailed0", "K^{0}_{S} signal, old data", "l");
    legendpointAngleHypoK->AddEntry("K0pointingAngleHypoBackgroundDetailed0", "K^{0}_{S} background, old data", "l");
    legendpointAngleHypoK->AddEntry("K0pointingAngleHypoSignalDetailed1", "K^{0}_{S} signal, new data, old tracks", "l");
    legendpointAngleHypoK->AddEntry("K0pointingAngleHypoBackgroundDetailed1", "K^{0}_{S} background, new data, old tracks", "l");
    legendpointAngleHypoK->AddEntry("K0pointingAngleHypoSignalDetailed2", "K^{0}_{S} signal, new data, new tracks", "l");
    legendpointAngleHypoK->AddEntry("K0pointingAngleHypoBackgroundDetailed2", "K^{0}_{S} background, new data, new tracks", "l");
    legendpointAngleHypoK->SetBorderSize(1);
    legendpointAngleHypoK->Draw("SAME");
    pointAngleHypoKCanvas->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/K0pointAngleHypoAll.pdf").c_str());
    delete pointAngleHypoKCanvas;





    //K0 after final cuts
    //old data
    TH1D *MKinvFinalCutClose1 = (TH1D *)anaoutputolddata1->Get("MpipiWideCheck");
    TH1D *MKinvFinalCutClose2 = (TH1D *)anaoutputolddata2->Get("MpipiWideCheck");
    TH1D resultKFinalCutHistOldData = *MKinvFinalCutClose1;
    resultKFinalCutHistOldData.Add(MKinvFinalCutClose2);
    // resultKFinalCutHistOldData.SetStats(1);
    //new data, old tracks
    MKinvFinalCutClose1 = (TH1D *)anaoutputoldtracks1->Get("MpipiWideCheck");
    MKinvFinalCutClose2 = (TH1D *)anaoutputoldtracks2->Get("MpipiWideCheck");
    TH1D resultKFinalCutHistNewDataOldTracks = *MKinvFinalCutClose1;
    resultKFinalCutHistNewDataOldTracks.Add(MKinvFinalCutClose2);
    // resultKFinalCutHistNewDataOldTracks.SetStats(1);
    //new data, new tracks
    MKinvFinalCutClose1 = (TH1D *)anaoutputnewtracks1->Get("MpipiWideCheck");
    MKinvFinalCutClose2 = (TH1D *)anaoutputnewtracks2->Get("MpipiWideCheck");
    TH1D resultKFinalCutHistNewDataNewTracks = *MKinvFinalCutClose1;
    resultKFinalCutHistNewDataNewTracks.Add(MKinvFinalCutClose2);
    //setting up the functions and data
    TCanvas *resultKFinalCut = new TCanvas("resultKFinalCut", "resultKFinalCut", 1800, 1600);
    resultKFinalCut->SetLeftMargin(0.15);
    TF1 *GfitKFinalCut = new TF1("GfitKFinalCut", "gausn(0) + pol1(3)", resultKHistOldData.GetXaxis()->GetXmin(), resultKHistOldData.GetXaxis()->GetXmax());
    TF1 *GfitKFinalCutBcg = new TF1("GfitKFinalCutBcg", "pol1", resultKHistOldData.GetXaxis()->GetXmin(), resultKHistOldData.GetXaxis()->GetXmax());
    TF1 *GfitKFinalCut1Sig = new TF1("GfitKFinalCut1Sig", "gausn", resultKHistOldData.GetXaxis()->GetXmin(), resultKHistOldData.GetXaxis()->GetXmax());
    GfitKFinalCut->SetParNames("Constant", "Mean", "Sigma", "b", "a");
    GfitKFinalCut->SetParameters(100, 0.5, 5e-3, 100, 0.48, 5e-3);
    // double nK0OldData, nK0NewDataOldTracks, nK0NewDataNewTracks;
    //old data, the base of drawing
    resultKFinalCutHistOldData.SetMinimum(0);
    resultKFinalCutHistOldData.SetMaximum(max(resultKFinalCutHistOldData.GetMaximum(), max(resultKFinalCutHistNewDataOldTracks.GetMaximum(), resultKFinalCutHistNewDataNewTracks.GetMaximum()))*1.1);
    resultKFinalCutHistOldData.SetMarkerStyle(kFullCircle);
    resultKFinalCutHistOldData.SetMarkerColor(kRed);
    resultKFinalCutHistOldData.SetLineColor(kRed+2);
    resultKFinalCutHistOldData.Fit(GfitKFinalCut, "0");
    // TH1D *tempHist = (TH1D *)resultKFinalCutHistOldData.DrawCopy("E", "OldData");
    resultKFinalCutHistOldData.DrawCopy("E", "OldData");
    gPad->Update();
    // TPaveStats *st = (TPaveStats *)tempHist->FindObject("stats");
    // st->SetX1NDC(0.15);
    // st->SetX2NDC(0.35);
    // st->SetY1NDC(0.65);
    // st->SetY2NDC(0.85);
    // st->SetBorderSize(0);
    // Double_t paramsK[5];
    GfitKFinalCut->GetParameters(paramsK);
    nK0OldData = paramsK[0];
    // cout<<"PARAMETERS:"<<endl;
    // for(int i = 0; i<5; i++){
    //     cout<<paramsK[i]<<endl;
    // }
    GfitKFinalCutBcg->SetParameters(paramsK+3);
    GfitKFinalCut1Sig->SetParameters(paramsK);
    GfitKFinalCutBcg->SetLineStyle(3);
    GfitKFinalCutBcg->SetLineWidth(3);
    GfitKFinalCut1Sig->SetLineColor(kRed);
    GfitKFinalCut->SetNpx(1000);
    GfitKFinalCut->DrawCopy("CSAME");
    GfitKFinalCutBcg->DrawCopy("CSAME");
    GfitKFinalCut1Sig->DrawCopy("CSAME");
    //new data old tracks
    resultKFinalCutHistNewDataOldTracks.SetMinimum(0);
    resultKFinalCutHistNewDataOldTracks.SetMarkerStyle(kFullCircle);
    resultKFinalCutHistNewDataOldTracks.SetMarkerColor(kGreen);
    resultKFinalCutHistNewDataOldTracks.SetLineColor(kGreen+2);
    resultKFinalCutHistNewDataOldTracks.Fit(GfitKFinalCut, "0");
    resultKFinalCutHistNewDataOldTracks.DrawCopy("E SAME", "NewDataOldTracks");
    gPad->Update();
    GfitKFinalCut->GetParameters(paramsK);
    nK0NewDataOldTracks = paramsK[0];
    GfitKFinalCutBcg->SetParameters(paramsK+3);
    GfitKFinalCut1Sig->SetParameters(paramsK);
    GfitKFinalCutBcg->SetLineStyle(3);
    GfitKFinalCutBcg->SetLineWidth(3);
    GfitKFinalCut->SetLineColor(kGreen);
    GfitKFinalCutBcg->SetLineColor(kGreen);
    GfitKFinalCut1Sig->SetLineColor(kGreen);
    GfitKFinalCut->SetNpx(1000);
    GfitKFinalCut->DrawCopy("CSAME");
    GfitKFinalCutBcg->DrawCopy("CSAME");
    GfitKFinalCut1Sig->DrawCopy("CSAME");
    //new data new tracks
    resultKFinalCutHistNewDataNewTracks.SetMinimum(0);
    resultKFinalCutHistNewDataNewTracks.SetMarkerStyle(kFullCircle);
    resultKFinalCutHistNewDataNewTracks.SetMarkerColor(kBlue);
    resultKFinalCutHistNewDataNewTracks.Fit(GfitKFinalCut, "0");
    resultKFinalCutHistNewDataNewTracks.DrawCopy("E SAME", "NewDataNewTracks");
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
    legendKFinalCut->AddEntry("MpipiWideCheckOldData", "pp old data, #sqrt{s} = 510 GeV");
    legendKFinalCut->AddEntry("MpipiWideCheckNewDataOldTracks", "pp new data, old tracks, #sqrt{s} = 510 GeV");
    legendKFinalCut->AddEntry("MpipiWideCheckNewDataNewTracks", "pp new data, new tracks, #sqrt{s} = 510 GeV");
    legendKFinalCut->SetBorderSize(0);
    legendKFinalCut->Draw("SAME");

    nK0OldData /= resultKFinalCutHistOldData.GetBinWidth(1);
    nK0NewDataOldTracks /= resultKFinalCutHistOldData.GetBinWidth(1);
    nK0NewDataNewTracks /= resultKFinalCutHistOldData.GetBinWidth(1);
    std::stringstream nK0FinalCutOldDataStream, nK0FinalCutNewDataOldTracksStream, nK0FinalCutNewDataNewTracksStream;
    nK0FinalCutOldDataStream<<std::fixed<<std::setprecision(2)<<nK0OldData;
    nK0FinalCutNewDataOldTracksStream<<std::fixed<<std::setprecision(2)<<nK0NewDataOldTracks;
    nK0FinalCutNewDataNewTracksStream<<std::fixed<<std::setprecision(2)<<nK0NewDataNewTracks;
    // TPaveStats *pt = new TPaveStats(0.62, 0.7, 0.89, 0.89, "NB NDC");
    pt = new TPaveStats(0.62, 0.7, 0.89, 0.89, "NB NDC");
    // TText *tempText;
    pt->SetFillColor(kWhite);
    tempText = pt->AddText(0, 0, "Number of K^{0}_{S} detected in:");
    tempText->SetTextSize(3);
    pt->AddText(0, 0.33, ("Old data = "+nK0FinalCutOldDataStream.str()).c_str());
    pt->AddText(0, 0.66, ("New data, old tracks = "+nK0FinalCutNewDataOldTracksStream.str()).c_str());
    pt->AddText(0, 1, ("New data, new tracks = "+nK0FinalCutNewDataNewTracksStream.str()).c_str());
    // pt->AddText("They are added to the pave using the AddText method.");
    // pt->AddLine(.0, .5, 1., .5);
    // pt->AddText("Even complex TLatex formulas can be added:");
    // TText *t1 = pt->AddText("F(t) = #sum_{i=-#infty}^{#infty}A(i)cos#[]{#frac{i}{t+i}}");
    // t1->SetTextColor(kBlue);
    pt->Draw("SAME");

    // gPad->RedrawAxis();
    resultKFinalCut->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/K0afterV0finderAfterCuts.pdf").c_str());













    ////////////////////////////////////////////////////////

    // LAMBDA

    ////////////////////////////////////////////////////////













    //Lambda0 backgrounds
    TCanvas *resultLambdaBkg = new TCanvas("resultLambdaBkg", "resultLambdaBkg", 1800, 1600);
    resultLambdaBkg->SetLeftMargin(0.15);
    //old data
    TH1D *MLambdainvClose1Bkg = (TH1D *)anaoutputolddata1->Get("MppiWideMissing");
    TH1D *MLambdainvClose2Bkg = (TH1D *)anaoutputolddata2->Get("MppiWideMissing");
    TH1D resultLambdaHistOldDataBkg = *MLambdainvClose1Bkg;
    resultLambdaHistOldData.Add(MLambdainvClose2Bkg);
    // resultLambdaHistOldData.SetStats(1);
    //new data, old tracks
    MLambdainvClose1Bkg = (TH1D *)anaoutputoldtracks1->Get("MppiWideMissing");
    MLambdainvClose2Bkg = (TH1D *)anaoutputoldtracks2->Get("MppiWideMissing");
    TH1D resultLambdaHistNewDataOldTracksBkg = *MLambdainvClose1Bkg;
    resultLambdaHistNewDataOldTracks.Add(MLambdainvClose2Bkg);
    // resultLambdaHistNewDataOldTracks.SetStats(1);
    //new data, new tracks
    MLambdainvClose1Bkg = (TH1D *)anaoutputnewtracks1->Get("MppiWideMissing");
    MLambdainvClose2Bkg = (TH1D *)anaoutputnewtracks2->Get("MppiWideMissing");
    TH1D resultLambdaHistNewDataNewTracksBkg = *MLambdainvClose1Bkg;
    resultLambdaHistNewDataNewTracks.Add(MLambdainvClose2Bkg);
    //old data, the base of drawing
    resultLambdaHistOldDataBkg.SetTitle("p^{#pm}#pi^{#mp} pairs rejected by basic cuts");
    resultLambdaHistOldDataBkg.SetMinimum(0);
    resultLambdaHistOldDataBkg.SetMaximum(max(resultLambdaHistOldDataBkg.GetMaximum(), max(resultLambdaHistNewDataOldTracksBkg.GetMaximum(), resultLambdaHistNewDataNewTracksBkg.GetMaximum()))*1.1);
    resultLambdaHistOldDataBkg.SetMarkerStyle(kFullCircle);
    resultLambdaHistOldDataBkg.SetMarkerColor(kRed);
    resultLambdaHistOldDataBkg.SetLineColor(kRed+2);
    resultLambdaHistOldDataBkg.DrawCopy("E", "OldDataBkg");
    //new data old tracks
    resultLambdaHistNewDataOldTracksBkg.SetMinimum(0);
    resultLambdaHistNewDataOldTracksBkg.SetMarkerStyle(kFullCircle);
    resultLambdaHistNewDataOldTracksBkg.SetMarkerColor(kGreen);
    resultLambdaHistNewDataOldTracksBkg.SetLineColor(kGreen+2);
    resultLambdaHistNewDataOldTracksBkg.DrawCopy("E SAME", "NewDataOldTracksBkg");
    //new data new tracks
    resultLambdaHistNewDataNewTracksBkg.SetMinimum(0);
    resultLambdaHistNewDataNewTracksBkg.SetMarkerStyle(kFullCircle);
    resultLambdaHistNewDataNewTracksBkg.SetMarkerColor(kBlue);
    resultLambdaHistNewDataNewTracksBkg.DrawCopy("E SAME", "NewDataNewTracksBkg");
    resultLambdaBkg->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/Lambda0afterV0finderBkg.pdf").c_str());


    //Lambda0 after PID
    //old data
    TH1D *MLambdainvPIDClose1 = (TH1D *)anaoutputolddata1->Get("MppiWideWithPidEcut");
    TH1D *MLambdainvPIDClose2 = (TH1D *)anaoutputolddata2->Get("MppiWideWithPidEcut");
    TH1D resultLambdaPIDHistOldData = *MLambdainvPIDClose1;
    resultLambdaPIDHistOldData.Add(MLambdainvPIDClose2);
    // resultLambdaPIDHistOldData.SetStats(1);
    //new data, old tracks
    MLambdainvPIDClose1 = (TH1D *)anaoutputoldtracks1->Get("MppiWideWithPidEcut");
    MLambdainvPIDClose2 = (TH1D *)anaoutputoldtracks2->Get("MppiWideWithPidEcut");
    TH1D resultLambdaPIDHistNewDataOldTracks = *MLambdainvPIDClose1;
    resultLambdaPIDHistNewDataOldTracks.Add(MLambdainvPIDClose2);
    // resultLambdaPIDHistNewDataOldTracks.SetStats(1);
    //new data, new tracks
    MLambdainvPIDClose1 = (TH1D *)anaoutputnewtracks1->Get("MppiWideWithPidEcut");
    MLambdainvPIDClose2 = (TH1D *)anaoutputnewtracks2->Get("MppiWideWithPidEcut");
    TH1D resultLambdaPIDHistNewDataNewTracks = *MLambdainvPIDClose1;
    resultLambdaPIDHistNewDataNewTracks.Add(MLambdainvPIDClose2);
    //setting up the functions and data
    TCanvas *resultLambdaPID = new TCanvas("resultLambdaPID", "resultLambdaPID", 1800, 1600);
    resultLambdaPID->SetLeftMargin(0.15);
    TF1 *GfitLambdaPID = new TF1("GfitLambdaPID", "gausn(0) + pol2(3)", resultLambdaHistOldData.GetXaxis()->GetXmin(), resultLambdaHistOldData.GetXaxis()->GetXmax());
    TF1 *GfitLambdaPIDBcg = new TF1("GfitLambdaPIDBcg", "pol2", resultLambdaHistOldData.GetXaxis()->GetXmin(), resultLambdaHistOldData.GetXaxis()->GetXmax());
    TF1 *GfitLambdaPID1Sig = new TF1("GfitLambdaPID1Sig", "gausn", resultLambdaHistOldData.GetXaxis()->GetXmin(), resultLambdaHistOldData.GetXaxis()->GetXmax());
    GfitLambdaPID->SetParNames("Constant", "Mean", "Sigma", "c", "b", "a");
    // GfitLambdaPID->SetParameters(100, 0.5, 5e-3, 100, 0.48, 5e-3);
    GfitLambdaPID->SetParameters(100, 1.115, 5e-3, -35000, 82000, -48000);
    GfitLambdaPID->SetParLimits(0, 0, 200);
    GfitLambdaPID->SetParLimits(1, 1.11, 1.12);
    GfitLambdaPID->SetParLimits(2, 0, 0.01);
    // double nLambda0OldData, nLambda0NewDataOldTracks, nLambda0NewDataNewTracks;
    //old data, the base of drawing
    resultLambdaPIDHistOldData.SetTitle("#Lambda^{0} mass in wide range with PID cuts");
    resultLambdaPIDHistOldData.SetMinimum(0);
    resultLambdaPIDHistOldData.SetMaximum(max(resultLambdaPIDHistOldData.GetMaximum(), max(resultLambdaPIDHistNewDataOldTracks.GetMaximum(), resultLambdaPIDHistNewDataNewTracks.GetMaximum()))*1.1);
    resultLambdaPIDHistOldData.SetMarkerStyle(kFullCircle);
    resultLambdaPIDHistOldData.SetMarkerColor(kRed);
    resultLambdaPIDHistOldData.SetLineColor(kRed+2);
    resultLambdaPIDHistOldData.Fit(GfitLambdaPID, "0BR");
    // TH1D *tempHist = (TH1D *)resultLambdaPIDHistOldData.DrawCopy("E", "OldData");
    resultLambdaPIDHistOldData.DrawCopy("E", "OldData");
    gPad->Update();
    // TPaveStats *st = (TPaveStats *)tempHist->FindObject("stats");
    // st->SetX1NDC(0.15);
    // st->SetX2NDC(0.35);
    // st->SetY1NDC(0.65);
    // st->SetY2NDC(0.85);
    // st->SetBorderSize(0);
    // Double_t paramsK[5];
    GfitLambdaPID->GetParameters(paramsLambda);
    nLambda0OldData = paramsLambda[0];
    // cout<<"PARAMETERS:"<<endl;
    // for(int i = 0; i<5; i++){
    //     cout<<paramsLambda[i]<<endl;
    // }
    GfitLambdaPIDBcg->SetParameters(paramsLambda+3);
    GfitLambdaPID1Sig->SetParameters(paramsLambda);
    GfitLambdaPIDBcg->SetLineStyle(3);
    GfitLambdaPIDBcg->SetLineWidth(3);
    GfitLambdaPID1Sig->SetLineColor(kRed);
    GfitLambdaPID->SetNpx(1000);
    GfitLambdaPID->DrawCopy("CSAME");
    GfitLambdaPIDBcg->DrawCopy("CSAME");
    GfitLambdaPID1Sig->DrawCopy("CSAME");
    //new data old tracks
    resultLambdaPIDHistNewDataOldTracks.SetMinimum(0);
    resultLambdaPIDHistNewDataOldTracks.SetMarkerStyle(kFullCircle);
    resultLambdaPIDHistNewDataOldTracks.SetMarkerColor(kGreen);
    resultLambdaPIDHistNewDataOldTracks.SetLineColor(kGreen+2);
    resultLambdaPIDHistNewDataOldTracks.Fit(GfitLambdaPID, "0BR");
    resultLambdaPIDHistNewDataOldTracks.DrawCopy("E SAME", "NewDataOldTracks");
    gPad->Update();
    GfitLambdaPID->GetParameters(paramsLambda);
    nLambda0NewDataOldTracks = paramsLambda[0];
    GfitLambdaPIDBcg->SetParameters(paramsLambda+3);
    GfitLambdaPID1Sig->SetParameters(paramsLambda);
    GfitLambdaPIDBcg->SetLineStyle(3);
    GfitLambdaPIDBcg->SetLineWidth(3);
    GfitLambdaPID->SetLineColor(kGreen);
    GfitLambdaPIDBcg->SetLineColor(kGreen);
    GfitLambdaPID1Sig->SetLineColor(kGreen);
    GfitLambdaPID->SetNpx(1000);
    GfitLambdaPID->DrawCopy("CSAME");
    GfitLambdaPIDBcg->DrawCopy("CSAME");
    GfitLambdaPID1Sig->DrawCopy("CSAME");
    //new data new tracks
    GfitLambdaPID->SetParameters(100, 1.115, 5e-3, -53000, 124000, -72000);
    resultLambdaPIDHistNewDataNewTracks.SetMinimum(0);
    resultLambdaPIDHistNewDataNewTracks.SetMarkerStyle(kFullCircle);
    resultLambdaPIDHistNewDataNewTracks.SetMarkerColor(kBlue);
    resultLambdaPIDHistNewDataNewTracks.Fit(GfitLambdaPID, "0BR");
    resultLambdaPIDHistNewDataNewTracks.DrawCopy("E SAME", "NewDataNewTracks");
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

    TLegend *legendLambdaPID = new TLegend(.18, .65, .45, .89);
    legendLambdaPID->SetTextSize(0.02);
    legendLambdaPID->AddEntry("MppiWideWithPidEcutOldData", "pp old data, #sqrt{s} = 510 GeV");
    legendLambdaPID->AddEntry("MppiWideWithPidEcutNewDataOldTracks", "pp new data, old tracks, #sqrt{s} = 510 GeV");
    legendLambdaPID->AddEntry("MppiWideWithPidEcutNewDataNewTracks", "pp new data, new tracks, #sqrt{s} = 510 GeV");
    legendLambdaPID->SetBorderSize(0);
    legendLambdaPID->Draw("SAME");

    nLambda0OldData /= resultLambdaPIDHistOldData.GetBinWidth(1);
    nLambda0NewDataOldTracks /= resultLambdaPIDHistOldData.GetBinWidth(1);
    nLambda0NewDataNewTracks /= resultLambdaPIDHistOldData.GetBinWidth(1);
    std::stringstream nLambda0PIDOldDataStream, nLambda0PIDNewDataOldTracksStream, nLambda0PIDNewDataNewTracksStream;
    nLambda0PIDOldDataStream<<std::fixed<<std::setprecision(2)<<nLambda0OldData;
    nLambda0PIDNewDataOldTracksStream<<std::fixed<<std::setprecision(2)<<nLambda0NewDataOldTracks;
    nLambda0PIDNewDataNewTracksStream<<std::fixed<<std::setprecision(2)<<nLambda0NewDataNewTracks;
    // TPaveStats *pt = new TPaveStats(0.62, 0.7, 0.89, 0.89, "NB NDC");
    pt = new TPaveStats(0.62, 0.7, 0.89, 0.89, "NB NDC");
    // TText *tempText;
    pt->SetFillColor(kWhite);
    tempText = pt->AddText(0, 0, "Number of #Lambda^{0} detected in:");
    tempText->SetTextSize(3);
    pt->AddText(0, 0.33, ("Old data = "+nLambda0PIDOldDataStream.str()).c_str());
    pt->AddText(0, 0.66, ("New data, old tracks = "+nLambda0PIDNewDataOldTracksStream.str()).c_str());
    pt->AddText(0, 1, ("New data, new tracks = "+nLambda0PIDNewDataNewTracksStream.str()).c_str());
    // pt->AddText("They are added to the pave using the AddText method.");
    // pt->AddLine(.0, .5, 1., .5);
    // pt->AddText("Even complex TLatex formulas can be added:");
    // TText *t1 = pt->AddText("F(t) = #sum_{i=-#infty}^{#infty}A(i)cos#[]{#frac{i}{t+i}}");
    // t1->SetTextColor(kBlue);
    pt->Draw("SAME");

    // gPad->RedrawAxis();
    resultLambdaPID->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/Lambda0afterV0finderPID.pdf").c_str());





    //Lambda0 pointAngleHypo
    TCanvas *pointAngleHypoLambdaCanvas = new TCanvas("pointAngleHypoLambdaCanvas", "pointAngleHypoLambdaCanvas", 1800, 1600);
    pointAngleHypoLambdaCanvas->SetGrid();
    // gPad->SetLogy();
    TH1D *pointAngleHypoLambdaTab[6];
    pointAngleHypoLambdaTab[0] = (TH1D *)anaoutputolddata1->Get("LambdapointingAngleHypoSignalDetailed");
    pointAngleHypoLambdaTab[0]->Add((TH1D *)anaoutputolddata2->Get("LambdapointingAngleHypoSignalDetailed"));
    pointAngleHypoLambdaTab[1] = (TH1D *)anaoutputolddata1->Get("LambdapointingAngleHypoBackgroundDetailed");
    pointAngleHypoLambdaTab[1]->Add((TH1D *)anaoutputolddata2->Get("LambdapointingAngleHypoBackgroundDetailed"));
    pointAngleHypoLambdaTab[2] = (TH1D *)anaoutputoldtracks1->Get("LambdapointingAngleHypoSignalDetailed");
    pointAngleHypoLambdaTab[2]->Add((TH1D *)anaoutputoldtracks2->Get("LambdapointingAngleHypoSignalDetailed"));
    pointAngleHypoLambdaTab[3] = (TH1D *)anaoutputoldtracks1->Get("LambdapointingAngleHypoBackgroundDetailed");
    pointAngleHypoLambdaTab[3]->Add((TH1D *)anaoutputoldtracks2->Get("LambdapointingAngleHypoBackgroundDetailed"));
    pointAngleHypoLambdaTab[4] = (TH1D *)anaoutputnewtracks1->Get("LambdapointingAngleHypoSignalDetailed");
    pointAngleHypoLambdaTab[4]->Add((TH1D *)anaoutputnewtracks2->Get("LambdapointingAngleHypoSignalDetailed"));
    pointAngleHypoLambdaTab[5] = (TH1D *)anaoutputnewtracks1->Get("LambdapointingAngleHypoBackgroundDetailed");
    pointAngleHypoLambdaTab[5]->Add((TH1D *)anaoutputnewtracks2->Get("LambdapointingAngleHypoBackgroundDetailed"));
    pointAngleHypoLambdaTab[0]->SetName("Lambda0pointingAngleHypoSignalDetailed0");
    pointAngleHypoLambdaTab[1]->SetName("Lambda0pointingAngleHypoBackgroundDetailed0");
    pointAngleHypoLambdaTab[2]->SetName("Lambda0pointingAngleHypoSignalDetailed1");
    pointAngleHypoLambdaTab[3]->SetName("Lambda0pointingAngleHypoBackgroundDetailed1");
    pointAngleHypoLambdaTab[4]->SetName("Lambda0pointingAngleHypoSignalDetailed2");
    pointAngleHypoLambdaTab[5]->SetName("Lambda0pointingAngleHypoBackgroundDetailed2");
    for(size_t i = 0; i<6; i++){
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
    pointAngleHypoLambdaTab[0]->SetTitle("pointingAngleHypo(), with normalization to the last bin");
    pointAngleHypoLambdaTab[0]->Draw("HIST");
    for(size_t i = 1; i<6; i++){
        pointAngleHypoLambdaTab[i]->Draw("HIST SAME");
    }
    TLegend *legendpointAngleHypoLambda = new TLegend(0.15, 0.55, 0.6, 0.85);
    legendpointAngleHypoLambda->SetTextSize(0.02);
    legendpointAngleHypoLambda->AddEntry("Lambda0pointingAngleHypoSignalDetailed0", "#Lambda^{0} signal, old data", "l");
    legendpointAngleHypoLambda->AddEntry("Lambda0pointingAngleHypoBackgroundDetailed0", "#Lambda^{0} background, old data", "l");
    legendpointAngleHypoLambda->AddEntry("Lambda0pointingAngleHypoSignalDetailed1", "#Lambda^{0} signal, new data, old tracks", "l");
    legendpointAngleHypoLambda->AddEntry("Lambda0pointingAngleHypoBackgroundDetailed1", "#Lambda^{0} background, new data, old tracks", "l");
    legendpointAngleHypoLambda->AddEntry("Lambda0pointingAngleHypoSignalDetailed2", "#Lambda^{0} signal, new data, new tracks", "l");
    legendpointAngleHypoLambda->AddEntry("Lambda0pointingAngleHypoBackgroundDetailed2", "#Lambda^{0} background, new data, new tracks", "l");
    legendpointAngleHypoLambda->SetBorderSize(1);
    legendpointAngleHypoLambda->Draw("SAME");
    pointAngleHypoLambdaCanvas->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/Lambda0pointAngleHypoAll.pdf").c_str());
    delete pointAngleHypoLambdaCanvas;





    //Lambda0 after final cuts
    //old data
    TH1D *MLambdainvFinalCutClose1 = (TH1D *)anaoutputolddata1->Get("MppiWideCheck");
    TH1D *MLambdainvFinalCutClose2 = (TH1D *)anaoutputolddata2->Get("MppiWideCheck");
    TH1D resultLambdaFinalCutHistOldData = *MLambdainvFinalCutClose1;
    resultLambdaFinalCutHistOldData.Add(MLambdainvFinalCutClose2);
    // resultLambdaFinalCutHistOldData.SetStats(1);
    //new data, old tracks
    MLambdainvFinalCutClose1 = (TH1D *)anaoutputoldtracks1->Get("MppiWideCheck");
    MLambdainvFinalCutClose2 = (TH1D *)anaoutputoldtracks2->Get("MppiWideCheck");
    TH1D resultLambdaFinalCutHistNewDataOldTracks = *MLambdainvFinalCutClose1;
    resultLambdaFinalCutHistNewDataOldTracks.Add(MLambdainvFinalCutClose2);
    // resultLambdaFinalCutHistNewDataOldTracks.SetStats(1);
    //new data, new tracks
    MLambdainvFinalCutClose1 = (TH1D *)anaoutputnewtracks1->Get("MppiWideCheck");
    MLambdainvFinalCutClose2 = (TH1D *)anaoutputnewtracks2->Get("MppiWideCheck");
    TH1D resultLambdaFinalCutHistNewDataNewTracks = *MLambdainvFinalCutClose1;
    resultLambdaFinalCutHistNewDataNewTracks.Add(MLambdainvFinalCutClose2);
    //setting up the functions and data
    TCanvas *resultLambdaFinalCut = new TCanvas("resultLambdaFinalCut", "resultLambdaFinalCut", 1800, 1600);
    resultLambdaFinalCut->SetLeftMargin(0.15);
    TF1 *GfitLambdaFinalCut = new TF1("GfitLambdaFinalCut", "gausn(0) + pol2(3)", resultLambdaHistOldData.GetXaxis()->GetXmin(), resultLambdaHistOldData.GetXaxis()->GetXmax());
    TF1 *GfitLambdaFinalCutBcg = new TF1("GfitLambdaFinalCutBcg", "pol2", resultLambdaHistOldData.GetXaxis()->GetXmin(), resultLambdaHistOldData.GetXaxis()->GetXmax());
    TF1 *GfitLambdaFinalCut1Sig = new TF1("GfitLambdaFinalCut1Sig", "gausn", resultLambdaHistOldData.GetXaxis()->GetXmin(), resultLambdaHistOldData.GetXaxis()->GetXmax());
    GfitLambdaFinalCut->SetParNames("Constant", "Mean", "Sigma", "c", "b", "a");
    // GfitLambdaFinalCut->SetParameters(100, 1.115, 5e-3, 2000000, -4200000, 1500000);
    GfitLambdaFinalCut->SetParameters(100, 1.115, 5e-3, -33000, 58000, -25000);
    GfitLambdaFinalCut->SetParLimits(0, 0, 200);
    GfitLambdaFinalCut->SetParLimits(1, 1.11, 1.12);
    GfitLambdaFinalCut->SetParLimits(2, 0, 0.01);
    // double nLambda0OldData, nLambda0NewDataOldTracks, nLambda0NewDataNewTracks;
    //old data, the base of drawing
    resultLambdaFinalCutHistOldData.SetMinimum(0);
    resultLambdaFinalCutHistOldData.SetMaximum(max(resultLambdaFinalCutHistOldData.GetMaximum(), max(resultLambdaFinalCutHistNewDataOldTracks.GetMaximum(), resultLambdaFinalCutHistNewDataNewTracks.GetMaximum()))*1.1);
    resultLambdaFinalCutHistOldData.SetMarkerStyle(kFullCircle);
    resultLambdaFinalCutHistOldData.SetMarkerColor(kRed);
    resultLambdaFinalCutHistOldData.SetLineColor(kRed+2);
    resultLambdaFinalCutHistOldData.Fit(GfitLambdaFinalCut, "0BR");
    // TH1D *tempHist = (TH1D *)resultLambdaFinalCutHistOldData.DrawCopy("E", "OldData");
    resultLambdaFinalCutHistOldData.DrawCopy("E", "OldData");
    gPad->Update();
    // TPaveStats *st = (TPaveStats *)tempHist->FindObject("stats");
    // st->SetX1NDC(0.15);
    // st->SetX2NDC(0.35);
    // st->SetY1NDC(0.65);
    // st->SetY2NDC(0.85);
    // st->SetBorderSize(0);
    // Double_t paramsLambda[5];
    GfitLambdaFinalCut->GetParameters(paramsLambda);
    nLambda0OldData = paramsLambda[0];
    // cout<<"PARAMETERS:"<<endl;
    // for(int i = 0; i<5; i++){
    //     cout<<paramsLambda[i]<<endl;
    // }
    GfitLambdaFinalCutBcg->SetParameters(paramsLambda+3);
    GfitLambdaFinalCut1Sig->SetParameters(paramsLambda);
    GfitLambdaFinalCutBcg->SetLineStyle(3);
    GfitLambdaFinalCutBcg->SetLineWidth(3);
    GfitLambdaFinalCut1Sig->SetLineColor(kRed);
    GfitLambdaFinalCut->SetNpx(1000);
    GfitLambdaFinalCut->DrawCopy("CSAME");
    GfitLambdaFinalCutBcg->DrawCopy("CSAME");
    GfitLambdaFinalCut1Sig->DrawCopy("CSAME");
    //new data old tracks
    GfitLambdaFinalCut->SetParameters(100, 1.115, 5e-3, -33000, 58000, -25000);
    resultLambdaFinalCutHistNewDataOldTracks.SetMinimum(0);
    resultLambdaFinalCutHistNewDataOldTracks.SetMarkerStyle(kFullCircle);
    resultLambdaFinalCutHistNewDataOldTracks.SetMarkerColor(kGreen);
    resultLambdaFinalCutHistNewDataOldTracks.SetLineColor(kGreen+2);
    resultLambdaFinalCutHistNewDataOldTracks.Fit(GfitLambdaFinalCut, "0BR");
    resultLambdaFinalCutHistNewDataOldTracks.DrawCopy("E SAME", "NewDataOldTracks");
    gPad->Update();
    GfitLambdaFinalCut->GetParameters(paramsLambda);
    nLambda0NewDataOldTracks = paramsLambda[0];
    GfitLambdaFinalCutBcg->SetParameters(paramsLambda+3);
    GfitLambdaFinalCut1Sig->SetParameters(paramsLambda);
    GfitLambdaFinalCutBcg->SetLineStyle(3);
    GfitLambdaFinalCutBcg->SetLineWidth(3);
    GfitLambdaFinalCut->SetLineColor(kGreen);
    GfitLambdaFinalCutBcg->SetLineColor(kGreen);
    GfitLambdaFinalCut1Sig->SetLineColor(kGreen);
    GfitLambdaFinalCut->SetNpx(1000);
    GfitLambdaFinalCut->DrawCopy("CSAME");
    GfitLambdaFinalCutBcg->DrawCopy("CSAME");
    GfitLambdaFinalCut1Sig->DrawCopy("CSAME");
    //new data new tracks
    GfitLambdaFinalCut->SetParameters(100, 1.115, 5e-3, -33000, 58000, -25000);
    resultLambdaFinalCutHistNewDataNewTracks.SetMinimum(0);
    resultLambdaFinalCutHistNewDataNewTracks.SetMarkerStyle(kFullCircle);
    resultLambdaFinalCutHistNewDataNewTracks.SetMarkerColor(kBlue);
    resultLambdaFinalCutHistNewDataNewTracks.Fit(GfitLambdaFinalCut, "0BR");
    resultLambdaFinalCutHistNewDataNewTracks.DrawCopy("E SAME", "NewDataNewTracks");
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
    legendLambdaFinalCut->AddEntry("MppiWideCheckOldData", "pp old data, #sqrt{s} = 510 GeV");
    legendLambdaFinalCut->AddEntry("MppiWideCheckNewDataOldTracks", "pp new data, old tracks, #sqrt{s} = 510 GeV");
    legendLambdaFinalCut->AddEntry("MppiWideCheckNewDataNewTracks", "pp new data, new tracks, #sqrt{s} = 510 GeV");
    legendLambdaFinalCut->SetBorderSize(0);
    legendLambdaFinalCut->Draw("SAME");

    nLambda0OldData /= resultLambdaFinalCutHistOldData.GetBinWidth(1);
    nLambda0NewDataOldTracks /= resultLambdaFinalCutHistOldData.GetBinWidth(1);
    nLambda0NewDataNewTracks /= resultLambdaFinalCutHistOldData.GetBinWidth(1);
    std::stringstream nLambda0FinalCutOldDataStream, nLambda0FinalCutNewDataOldTracksStream, nLambda0FinalCutNewDataNewTracksStream;
    nLambda0FinalCutOldDataStream<<std::fixed<<std::setprecision(2)<<nLambda0OldData;
    nLambda0FinalCutNewDataOldTracksStream<<std::fixed<<std::setprecision(2)<<nLambda0NewDataOldTracks;
    nLambda0FinalCutNewDataNewTracksStream<<std::fixed<<std::setprecision(2)<<nLambda0NewDataNewTracks;
    // TPaveStats *pt = new TPaveStats(0.62, 0.7, 0.89, 0.89, "NB NDC");
    pt = new TPaveStats(0.62, 0.7, 0.89, 0.89, "NB NDC");
    // TText *tempText;
    pt->SetFillColor(kWhite);
    tempText = pt->AddText(0, 0, "Number of #Lambda^{0} detected in:");
    tempText->SetTextSize(3);
    pt->AddText(0, 0.33, ("Old data = "+nLambda0FinalCutOldDataStream.str()).c_str());
    pt->AddText(0, 0.66, ("New data, old tracks = "+nLambda0FinalCutNewDataOldTracksStream.str()).c_str());
    pt->AddText(0, 1, ("New data, new tracks = "+nLambda0FinalCutNewDataNewTracksStream.str()).c_str());
    // pt->AddText("They are added to the pave using the AddText method.");
    // pt->AddLine(.0, .5, 1., .5);
    // pt->AddText("Even complex TLatex formulas can be added:");
    // TText *t1 = pt->AddText("F(t) = #sum_{i=-#infty}^{#infty}A(i)cos#[]{#frac{i}{t+i}}");
    // t1->SetTextColor(kBlue);
    pt->Draw("SAME");

    // gPad->RedrawAxis();
    resultLambdaFinalCut->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/Lambda0afterV0finderAfterCuts.pdf").c_str());


    return 0;
}