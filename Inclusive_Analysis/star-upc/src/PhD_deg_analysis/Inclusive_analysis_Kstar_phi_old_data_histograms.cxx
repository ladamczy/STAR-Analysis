// Stdlib header file for input and output.
#include <iostream>
#include <cstring>
#include <fstream>
#include <sstream>

// ROOT, for histogramming.
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TEfficiency.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"

#include "MyStyles.h"

void draw_and_save(TH1D* data, std::string folderWithDiagonal, std::string name, std::string title);
void draw_and_save_minus_background(TH1D* data, TH1D* bcg, std::string folderWithDiagonal, std::string name, std::string title, double bcg_region);

int main(int argc, char const* argv[]){
    TFile* input = TFile::Open(argv[1]);

    //getting histograms out
    TH1D* pairInfo = (TH1D*)input->Get("pairInfo");
    TH1D* MKpiChi2 = (TH1D*)input->Get("MKpiChi2");
    TH1D* MpiKChi2 = (TH1D*)input->Get("MpiKChi2");
    TH1D* MppiChi2 = (TH1D*)input->Get("MppiChi2");
    TH1D* MpipChi2 = (TH1D*)input->Get("MpipChi2");
    TH1D* MKKChi2 = (TH1D*)input->Get("MKKChi2");
    TH1D* MpipiChi2 = (TH1D*)input->Get("MpipiChi2");
    TH1D* MppChi2 = (TH1D*)input->Get("MppChi2");
    /*
    MKpiChi2
    MpiKChi2
    MppiChi2
    MpipChi2
    MKKChi2
    MpipiChi2
    MppChi2
    */

    TFile* inputbcg = TFile::Open(argv[2]);

    //getting background histograms out
    //MIXEDTOF, NORMALTOF, NOTOF
    std::string background = "NORMALTOF";
    TH1D* MKpiChi2bcg = (TH1D*)inputbcg->Get(("MKpiChi2bcg"+background).c_str());
    TH1D* MpiKChi2bcg = (TH1D*)inputbcg->Get(("MpiKChi2bcg"+background).c_str());
    TH1D* MppiChi2bcg = (TH1D*)inputbcg->Get(("MppiChi2bcg"+background).c_str());
    TH1D* MpipChi2bcg = (TH1D*)inputbcg->Get(("MpipChi2bcg"+background).c_str());
    TH1D* MKKChi2bcg = (TH1D*)inputbcg->Get(("MKKChi2bcg"+background).c_str());
    TH1D* MpipiChi2bcg = (TH1D*)inputbcg->Get(("MpipiChi2bcg"+background).c_str());
    TH1D* MppChi2bcg = (TH1D*)inputbcg->Get(("MppChi2bcg"+background).c_str());
    //getting NORMALTOF diffractive background and signal
    // std::vector<std::string> pairTab = { "Kpi", "piK", "ppi", "pip", "KK", "pipi", "pp" };
    // std::vector<double> bcgRegion = { 1.0, 1.0, 1.1, 1.1, 1.1, 0.2, 2.4 };
    std::vector<std::string> pairTab = { "Kpi", "piK", "KK" };
    std::vector<double> bcgRegion = { 1.0, 1.0, 1.1 };
    //creating list of categories
    std::vector<std::string> allCategories;
    std::ifstream infile("STAR-Analysis/Inclusive_Analysis/star-upc/src/PhD_deg_analysis/Differential_crossection_values.txt");
    std::string line, buf;
    while(std::getline(infile, line)){
        std::stringstream ss(line);
        ss>>buf;
        //erasing ":" after the category
        if(buf.find(':')!=std::string::npos){
            buf.erase(buf.end()-1);
        }
        allCategories.push_back(buf);
    }
    //filling vectors of background and signal
    std::vector<TH2D*> background_vector, signal_vector;
    std::vector<TH1D*> result_vector;
    std::string tempSignalName = "M$Chi2";
    std::string tempBackgroundName = "M$Chi2bcgNORMALTOF";
    std::string tempHistName;
    for(auto&& pair:pairTab){
        for(auto&& category:allCategories){
            //signal
            tempHistName = tempSignalName;
            tempHistName.replace(find(tempHistName.begin(), tempHistName.end(), '$')-tempHistName.begin(), 1, pair);
            signal_vector.push_back((TH2D*)input->Get((tempHistName+category).c_str()));
            //result, but uses signal template
            result_vector.push_back(new TH1D((tempHistName+"Result"+category).c_str(), (pair+" "+category).c_str(), signal_vector.back()->GetNbinsY(), signal_vector.back()->GetYaxis()->GetXmin(), signal_vector.back()->GetYaxis()->GetXmax()));
            //background
            tempHistName = tempBackgroundName;
            tempHistName.replace(find(tempHistName.begin(), tempHistName.end(), '$')-tempHistName.begin(), 1, pair);
            background_vector.push_back((TH2D*)inputbcg->Get((tempHistName+category).c_str()));
        }
    }
    //fitting function
    TF1* fitFunc = new TF1("fitFunc", "breitwigner+pol1(3)", 0.8, 1.0);
    //substracting one from another
    //fitting the difference
    //and filling the result
    for(size_t i = 0; i<pairTab.size(); i++){
        for(size_t j = 0; j<allCategories.size(); j++){
            TH2D* bcg_pointer = background_vector[i*allCategories.size()+j];
            TH2D* sig_pointer = signal_vector[i*allCategories.size()+j];
            for(size_t k = 0; k<sig_pointer->GetNbinsY(); k++){
                //getting a slice and removing background
                TH1D* bcg_slice = bcg_pointer->ProjectionX("_px", k+1, k+1, "e");
                TH1D* sig_slice = sig_pointer->ProjectionX("_px", k+1, k+1, "e");
                int bcg_bin = bcg_slice->FindBin(bcgRegion[i]);
                if(sig_slice->Integral(bcg_bin, -1)==0||bcg_slice->Integral(bcg_bin, -1)==0){
                    result_vector[i*allCategories.size()+j]->SetBinContent(k+1, TMath::QuietNaN());
                    result_vector[i*allCategories.size()+j]->SetBinError(k+1, TMath::QuietNaN());
                    continue;
                }
                bcg_slice->Scale(sig_slice->Integral(bcg_bin, -1)/bcg_slice->Integral(bcg_bin, -1));
                sig_slice->Add(bcg_slice, -1.);
                //fitting and filling result
                double max_center = sig_slice->GetBinCenter(sig_slice->GetMaximumBin());
                fitFunc->SetParameters(100., 0.1, max_center, 0., 0.);
                fitFunc->SetRange(max_center-0.1, max_center+0.1);
                sig_slice->Fit(fitFunc, "R");
                sig_slice->Draw();
                result_vector[i*allCategories.size()+j]->SetBinContent(k+1, fitFunc->GetParameter(0));
                result_vector[i*allCategories.size()+j]->SetBinError(k+1, fitFunc->GetParError(0));
            }
        }
    }

    std::string folderWithDiagonal = std::string(argv[3]);
    if(folderWithDiagonal[folderWithDiagonal.size()-1]!='/'){
        folderWithDiagonal += "/";
    }

    gStyle->SetOptStat(0);
    gStyle->SetFrameLineWidth(2);

    //drawing
    //pairInfo
    TCanvas* resultCanvas = new TCanvas("resultCanvas", "resultCanvas", 4000, 2400);
    pairInfo->SetMinimum(0);
    pairInfo->SetLineColor(kBlue+2);
    pairInfo->SetMarkerSize(2);
    //TODO: add fitting
    // pairInfo->Draw("E");
    pairInfo->Draw("hist");
    pairInfo->Draw("same text0");
    pairInfo->SetLineWidth(2);
    pairInfo->GetXaxis()->SetLabelSize(0.06);
    pairInfo->GetXaxis()->SetTitleSize(0.06);
    resultCanvas->SetTopMargin(0.05);
    gPad->Update();
    resultCanvas->SaveAs((folderWithDiagonal+"pairInfo.pdf").c_str());

    //the rest
    draw_and_save(MKpiChi2, folderWithDiagonal, "MKpi", "K^{+}#pi^{-}");
    draw_and_save(MpiKChi2, folderWithDiagonal, "MpiK", "#pi^{+}K^{-}");
    draw_and_save(MppiChi2, folderWithDiagonal, "Mppi", "p^{+}#pi^{-}");
    draw_and_save(MpipChi2, folderWithDiagonal, "Mpip", "#pi^{+}p^{-}");
    draw_and_save(MKKChi2, folderWithDiagonal, "MKK", "K^{+}K^{-}");
    draw_and_save(MpipiChi2, folderWithDiagonal, "Mpipi", "#pi^{+}#pi^{-}");
    draw_and_save(MppChi2, folderWithDiagonal, "Mpp", "p^{+}p^{-}");
    draw_and_save_minus_background(MKpiChi2, MKpiChi2bcg, folderWithDiagonal, "MKpiFit", "K^{+}#pi^{-} background removed", 1.0);
    draw_and_save_minus_background(MpiKChi2, MpiKChi2bcg, folderWithDiagonal, "MpiKFit", "#pi^{+}K^{-} background removed", 1.0);
    draw_and_save_minus_background(MppiChi2, MppiChi2bcg, folderWithDiagonal, "MppiFit", "p^{+}#pi^{-} background removed", 1.1);
    draw_and_save_minus_background(MpipChi2, MpipChi2bcg, folderWithDiagonal, "MpipFit", "#pi^{+}p^{-} background removed", 1.1);
    draw_and_save_minus_background(MKKChi2, MKKChi2bcg, folderWithDiagonal, "MKKFit", "K^{+}K^{-} background removed", 1.1);
    draw_and_save_minus_background(MpipiChi2, MpipiChi2bcg, folderWithDiagonal, "MpipiFit", "#pi^{+}#pi^{-} background removed", 0.2); //originally 0.8
    draw_and_save_minus_background(MppChi2, MppChi2bcg, folderWithDiagonal, "MppFit", "p^{+}p^{-} background removed", 2.4);
    draw_and_save(MKpiChi2bcg, folderWithDiagonal, "MKpibcg", "K^{+}#pi^{-} background");
    draw_and_save(MpiKChi2bcg, folderWithDiagonal, "MpiKbcg", "#pi^{+}K^{-} background");
    draw_and_save(MppiChi2bcg, folderWithDiagonal, "Mppibcg", "p^{+}#pi^{-} background");
    draw_and_save(MpipChi2bcg, folderWithDiagonal, "Mpipbcg", "#pi^{+}p^{-} background");
    draw_and_save(MKKChi2bcg, folderWithDiagonal, "MKKbcg", "K^{+}K^{-} background");
    draw_and_save(MpipiChi2bcg, folderWithDiagonal, "Mpipibcg", "#pi^{+}#pi^{-} background");
    draw_and_save(MppChi2bcg, folderWithDiagonal, "Mppbcg", "p^{+}p^{-} background");
    //fits
    for(size_t i = 0; i<pairTab.size(); i++){
        for(size_t j = 0; j<allCategories.size(); j++){
            draw_and_save(result_vector[i*allCategories.size()+j], "~/star-upc/", "M"+pairTab[i]+allCategories[j], pairTab[i]+" "+allCategories[j]+";"+allCategories[j]+";entries");
        }
    }
    return 0;
}

void draw_and_save(TH1D* data, std::string folderWithDiagonal, std::string name, std::string title){
    MyStyles styleLibrary;
    TStyle tempStyle = styleLibrary.Hist2DQuarterSize(true);
    tempStyle.cd();
    gROOT->ForceStyle();
    TCanvas* resultCanvas = new TCanvas("resultCanvas", "resultCanvas", 4000, 2400);
    //TODO: add fitting
    data->Draw("hist");
    resultCanvas->UseCurrentStyle();
    data->SetMinimum(0);
    data->SetMarkerStyle(kFullCircle);
    data->SetMarkerSize(2);
    data->SetMarkerColor(kBlue);
    data->SetTitle(title.c_str());
    resultCanvas->SaveAs((folderWithDiagonal+name+".pdf").c_str());
}

void draw_and_save_minus_background(TH1D* data, TH1D* bcg, std::string folderWithDiagonal, std::string name, std::string title, double bcg_region){
    MyStyles styleLibrary;
    TStyle tempStyle = styleLibrary.Hist2DQuarterSize(true);
    tempStyle.cd();
    gROOT->ForceStyle();
    TCanvas* resultCanvas = new TCanvas("resultCanvas", "resultCanvas", 4000, 2400);
    int bcg_bin = data->FindBin(bcg_region);
    bcg->Scale(data->Integral(bcg_bin, -1)/bcg->Integral(bcg_bin, -1));
    data->Add(bcg, -1.);
    data->Draw("e");
    resultCanvas->UseCurrentStyle();
    data->SetMarkerStyle(kFullCircle);
    data->SetMarkerSize(2);
    data->SetMarkerColor(kBlue);
    data->SetTitle(title.c_str());

    resultCanvas->SaveAs((folderWithDiagonal+name+".pdf").c_str());
}