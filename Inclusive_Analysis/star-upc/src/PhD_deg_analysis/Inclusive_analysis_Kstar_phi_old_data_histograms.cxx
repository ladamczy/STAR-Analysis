// Stdlib header file for input and output.
#include <iostream>
#include <cstring>
#include <fstream>
#include <sstream>

// ROOT, for histogramming.
#include "TH1.h"
#include "TH2.h"
#include "TPaveStats.h"
#include "TFile.h"
#include "TEfficiency.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"
#include "TF1Convolution.h"
#include "TApplication.h"
#include "TFitResult.h"

#include "MyStyles.h"

int GetFirstNonzeroBinNumber(TH1* input);
void draw_and_save(TH1D* data, std::string folderWithDiagonal, std::string name, std::string title, std::string options = "");
void draw_and_save_minus_background(TH1D* data, TH1D* bcg, std::string folderWithDiagonal, std::string name, std::string title, double bcg_region);
TFitResult differential_crossection_fit(TPad* pad, TH1D* slice, TF1* fitting_function_signal, TF1* fitting_function_bcg, std::string draw_full_path = "");

int main(int argc, char* argv[]){
    //something used so that the histograms would draw
    //https://stackoverflow.com/questions/30932725/painting-a-tcanvas-to-the-screen-in-a-compiled-root-cern-application
    TApplication theApp("App", &argc, argv);
    argv = theApp.Argv();

    TFile* input = TFile::Open(static_cast<const char*>(argv[1]));

    //getting histograms out
    TH1D* pairInfo = (TH1D*)input->Get("pairInfo");
    TH1D* MKpiChi2 = (TH1D*)input->Get("MKpiChi2");
    TH1D* MpiKChi2 = (TH1D*)input->Get("MpiKChi2");
    TH1D* MppiChi2 = (TH1D*)input->Get("MppiChi2");
    TH1D* MpipChi2 = (TH1D*)input->Get("MpipChi2");
    TH1D* MKKChi2 = (TH1D*)input->Get("MKKChi2");
    TH1D* MpipiChi2 = (TH1D*)input->Get("MpipiChi2");
    TH1D* MppChi2 = (TH1D*)input->Get("MppChi2");
    TH1D* MKKChi2Close = (TH1D*)input->Get("MKKChi2Close");
    TH1D* MKpiChi2Close = (TH1D*)input->Get("MKpiChi2Close");
    TH1D* MpiKChi2Close = (TH1D*)input->Get("MpiKChi2Close");
    /*
    MKpiChi2
    MpiKChi2
    MppiChi2
    MpipChi2
    MKKChi2
    MpipiChi2
    MppChi2
    */

    TFile* inputbcg = TFile::Open(static_cast<const char*>(argv[2]));

    //getting background histograms out
    TH1D* MKpiChi2bcg = (TH1D*)inputbcg->Get("MKpiChi2bcgTOF");
    TH1D* MpiKChi2bcg = (TH1D*)inputbcg->Get("MpiKChi2bcgTOF");
    TH1D* MppiChi2bcg = (TH1D*)inputbcg->Get("MppiChi2bcgTOF");
    TH1D* MpipChi2bcg = (TH1D*)inputbcg->Get("MpipChi2bcgTOF");
    TH1D* MKKChi2bcg = (TH1D*)inputbcg->Get("MKKChi2bcgTOF");
    TH1D* MpipiChi2bcg = (TH1D*)inputbcg->Get("MpipiChi2bcgTOF");
    TH1D* MppChi2bcg = (TH1D*)inputbcg->Get("MppChi2bcgTOF");
    //getting diffractive background and signal
    // std::vector<std::string> pairTab = { "Kpi", "piK", "ppi", "pip", "KK", "pipi", "pp" };
    // std::vector<double> bcgRegion = { 1.0, 1.0, 1.1, 1.1, 1.1, 0.2, 2.4 };
    std::vector<std::string> pairTab = { "Kpi", "piK", "KK" };
    std::vector<double> bcgRegion = { 1.0, 1.0, 1.1 };
    std::vector<double> fitMaximum = { 0.892, 0.892, 1.02 };
    std::vector<double> fitWidth = { 0.0514, 0.0514, 0.004 };//51.4 MeV for K*(892), 4.43 MeV for phi(1020)
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

    //directory manipulation
    std::string folderWithDiagonal = std::string(static_cast<const char*>(argv[3]));
    if(folderWithDiagonal[folderWithDiagonal.size()-1]!='/'){
        folderWithDiagonal += "/";
    }

    //filling vectors of background and signal
    //and Chi2 with and without the background
    std::vector<TH2D*> background_vector, signal_vector;
    std::vector<TH1D*> result_vector, result_vector_nobcgfit, result_vector_nobcgremoval;
    std::vector<TH1D*> width_vector, width_vector_nobcgfit, width_vector_nobcgremoval;
    std::vector<TH1D*> Chi2withbcg_vector, Chi2withoutbcg_vector, Chi2withoutremovingbcg_vector;
    std::string tempSignalName = "M$Chi2";
    std::string tempBackgroundName = "M$Chi2bcgTOF";
    std::string tempHistName;
    for(auto&& pair:pairTab){
        for(auto&& category:allCategories){
            //signal
            tempHistName = tempSignalName;
            tempHistName.replace(find(tempHistName.begin(), tempHistName.end(), '$')-tempHistName.begin(), 1, pair);
            signal_vector.push_back((TH2D*)input->Get((tempHistName+category).c_str()));
            //background
            tempHistName = tempBackgroundName;
            tempHistName.replace(find(tempHistName.begin(), tempHistName.end(), '$')-tempHistName.begin(), 1, pair);
            background_vector.push_back((TH2D*)inputbcg->Get((tempHistName+category).c_str()));
            //result, but uses signal template
            result_vector.push_back(new TH1D((tempHistName+"Result"+category).c_str(), (pair+" "+category).c_str(), signal_vector.back()->GetNbinsY(), signal_vector.back()->GetYaxis()->GetXmin(), signal_vector.back()->GetYaxis()->GetXmax()));
            result_vector_nobcgfit.push_back(new TH1D((tempHistName+"Result"+category).c_str(), (pair+" "+category).c_str(), signal_vector.back()->GetNbinsY(), signal_vector.back()->GetYaxis()->GetXmin(), signal_vector.back()->GetYaxis()->GetXmax()));
            result_vector_nobcgremoval.push_back(new TH1D((tempHistName+"Result"+category).c_str(), (pair+" "+category).c_str(), signal_vector.back()->GetNbinsY(), signal_vector.back()->GetYaxis()->GetXmin(), signal_vector.back()->GetYaxis()->GetXmax()));
            //width of the fits
            width_vector.push_back(new TH1D((tempHistName+"Width"+category).c_str(), (pair+" "+category).c_str(), signal_vector.back()->GetNbinsY(), signal_vector.back()->GetYaxis()->GetXmin(), signal_vector.back()->GetYaxis()->GetXmax()));
            width_vector_nobcgfit.push_back(new TH1D((tempHistName+"Width"+category).c_str(), (pair+" "+category).c_str(), signal_vector.back()->GetNbinsY(), signal_vector.back()->GetYaxis()->GetXmin(), signal_vector.back()->GetYaxis()->GetXmax()));
            width_vector_nobcgremoval.push_back(new TH1D((tempHistName+"Width"+category).c_str(), (pair+" "+category).c_str(), signal_vector.back()->GetNbinsY(), signal_vector.back()->GetYaxis()->GetXmin(), signal_vector.back()->GetYaxis()->GetXmax()));
            //Chi2 comparison
            Chi2withbcg_vector.push_back(new TH1D((tempHistName+"Chi2withbcg"+category).c_str(), (pair+" "+category+" bcg").c_str(), signal_vector.back()->GetNbinsY(), signal_vector.back()->GetYaxis()->GetXmin(), signal_vector.back()->GetYaxis()->GetXmax()));
            Chi2withoutbcg_vector.push_back(new TH1D((tempHistName+"Chi2withoutbcg"+category).c_str(), (pair+" "+category+" no fitted bcg").c_str(), signal_vector.back()->GetNbinsY(), signal_vector.back()->GetYaxis()->GetXmin(), signal_vector.back()->GetYaxis()->GetXmax()));
            Chi2withoutremovingbcg_vector.push_back(new TH1D((tempHistName+"Chi2withoutremovingbcg"+category).c_str(), (pair+" "+category+" not removed bcg").c_str(), signal_vector.back()->GetNbinsY(), signal_vector.back()->GetYaxis()->GetXmin(), signal_vector.back()->GetYaxis()->GetXmax()));
        }
    }

    //for skipping
    std::string wannaSkip;
    bool skippedFitting = false;
    bool skippedFittingNoBackground = false;
    bool skippedFittingNotRemovedBackground = false;
    //###########################################################
    //                FITTING
    //###########################################################
    TCanvas* result = new TCanvas("result", "result", -1, 0, 1600, 900);
    //fitting ONCE the total m_KK
    TF1Convolution conv_sig("breitwigner", "gausn", 0.5, 1.5); //extra range for all your convolution needs
    conv_sig.SetNofPointsFFT(10000);
    TF1 fit_func_custom_sig("fit_func_custom_sig", conv_sig, 0.99, 1.05, 6);
    TF1 fit_func_custom_bcg("fit_func_custom_bcg", "pol2", 0.99, 1.05);
    printf("If you want to skip the fitting of one-time things, write \"yes\"\n");
    std::getline(std::cin, wannaSkip);
    if(wannaSkip.find("yes")!=std::string::npos){
        goto skipOneTimeFitting;
    }
    fit_func_custom_sig.SetParameter(0, 9);
    fit_func_custom_sig.SetParameter(1, 1.02);
    fit_func_custom_sig.FixParameter(2, 0.00443);
    fit_func_custom_sig.FixParameter(3, 1.);
    fit_func_custom_sig.FixParameter(4, 0.);
    fit_func_custom_sig.SetParameter(5, 0.0013);
    fit_func_custom_sig.SetParNames("N_{BW}", "m_{0}", "#Gamma_{0,BW}", "N_{Gauss}", "m_{0,Gauss}", "#sigma_{Gauss}");
    fit_func_custom_bcg.SetParNames("c", "b", "a");
    differential_crossection_fit(result, MKKChi2Close, &fit_func_custom_sig, &fit_func_custom_bcg, folderWithDiagonal+"MKKwhole.pdf");
    //fitting ONCE the total m_Kpi
    fit_func_custom_sig.SetRange(0.7, 1.1);
    fit_func_custom_bcg.SetRange(0.7, 1.1);
    fit_func_custom_sig.SetParameter(0, 9);
    fit_func_custom_sig.SetParameter(1, 0.892);
    fit_func_custom_sig.FixParameter(2, 0.0514);
    fit_func_custom_sig.FixParameter(3, 1.);
    fit_func_custom_sig.FixParameter(4, 0.);
    fit_func_custom_sig.SetParameter(5, 0.0013);
    differential_crossection_fit(result, MKpiChi2Close, &fit_func_custom_sig, &fit_func_custom_bcg, folderWithDiagonal+"MKpiwhole.pdf");
    //fitting ONCE the total m_piK
    fit_func_custom_sig.SetRange(0.7, 1.1);
    fit_func_custom_bcg.SetRange(0.7, 1.1);
    fit_func_custom_sig.SetParameter(0, 9);
    fit_func_custom_sig.SetParameter(1, 0.892);
    fit_func_custom_sig.FixParameter(2, 0.0514);
    fit_func_custom_sig.FixParameter(3, 1.);
    fit_func_custom_sig.FixParameter(4, 0.);
    fit_func_custom_sig.SetParameter(5, 0.0013);
    differential_crossection_fit(result, MpiKChi2Close, &fit_func_custom_sig, &fit_func_custom_bcg, folderWithDiagonal+"MpiKwhole.pdf");
skipOneTimeFitting:

    //fitting functions
    TF1* fit_func_sig = new TF1("fit_func_sig", "breitwigner", 0.8, 1.0);
    //backgrounds - one custom, with parameters p0, p1
    //being values at ends (p2, p3)
    auto custom_background = [&](double* x, double* p){
        //fitting behaves weird with demand of p0 & p1 >0
        //so p0 & p1 will actually be square roots of the end values
        // double a = (p[1]*p[1]-p[0]*p[0])/(p[3]-p[2]);
        // return a*(x[0]-p[2])+p[0]*p[0];
        // it didn't work very well, going back to linear background
        double a = (p[1]-p[0])/(p[3]-p[2]);
        return a*(x[0]-p[2])+p[0];
    };
    TF1* fit_func_bcg = new TF1("fit_func_bcg", custom_background, 0.8, 1.0, 4);
    TF1* fit_func_empty_bcg = new TF1("fit_func_bcg", "0.", 0.8, 1.0);
    //custom backgroundfrom the paper
    auto even_more_custom_background = [&](double* x, double* p){
        double A = p[0];
        double B = p[1];
        double C = p[2];
        double m = x[0];
        double mK = 0.493677;
        return (1-exp((2*mK-m)/C))*pow(m/(2*mK), A)+B*(m/(2*mK)-1);
    };
    TF1* fit_func_paper_bcg = new TF1("fit_func_paper_bcg", even_more_custom_background, 0.8, 1.0, 3);
    TF1* fit_func_for_Kstar = new TF1("fit_func_for_Kstar", "pol2", 0.8, 1.0);
    //substracting one from another
    //fitting the difference
    //and filling the result
    printf("If you want to skip the fitting, write \"yes\"\n");
    std::getline(std::cin, wannaSkip);
    if(wannaSkip.find("yes")!=std::string::npos){
        skippedFitting = true;
        goto fittingBackground;
    }
    for(size_t i = 0; i<pairTab.size(); i++){
        for(size_t j = 0; j<allCategories.size(); j++){
            TH2D* bcg_pointer = (TH2D*)background_vector[i*allCategories.size()+j]->Clone();
            TH2D* sig_pointer = (TH2D*)signal_vector[i*allCategories.size()+j]->Clone();
            //background fitting moved from binned to overall sum
            //background still removed like that
            int bcg_bin = bcg_pointer->GetXaxis()->FindBin(bcgRegion[i]);
            bcg_pointer->Scale(sig_pointer->Integral(bcg_bin, -1, 0, -1)/bcg_pointer->Integral(bcg_bin, -1, 0, -1));
            sig_pointer->Add(bcg_pointer, -1.);
            //for keeping title
            std::string baseOfTitle = std::string(sig_pointer->GetTitle())+" ";
            for(Int_t k = 0; k<sig_pointer->GetNbinsY(); k++){
                TH1D* sig_slice = sig_pointer->ProjectionX("_bcg", k+1, k+1, "e");
                //fitting and filling result
                //setting lower range for m0-3*gamma
                //if lower bound is lower than m0-4*gamma
                //and 2*width otherwise
                double lower_range = fitMaximum[i]-3*fitWidth[i];
                if(sig_slice->GetBinLowEdge(GetFirstNonzeroBinNumber(sig_slice))>fitMaximum[i]-4*fitWidth[i]){
                    lower_range = fitMaximum[i]-2*fitWidth[i];
                }
                fit_func_sig->SetParameters(10., fitMaximum[i], fitWidth[i]);
                fit_func_sig->SetRange(lower_range, 1.1);
                //custom setting background function - p2 & p3 are the ends of the range
                fit_func_bcg->SetRange(lower_range, 1.1);
                fit_func_bcg->SetParameters(0., 0., fit_func_bcg->GetXmin(), fit_func_bcg->GetXmax());
                fit_func_bcg->FixParameter(2, fit_func_bcg->GetXmin());
                fit_func_bcg->FixParameter(3, fit_func_bcg->GetXmax());
                double par_value, par_error;
                std::string newTitle = baseOfTitle;
                newTitle += std::to_string(sig_pointer->GetYaxis()->GetBinLowEdge(k+1))+" - ";
                newTitle += std::to_string(sig_pointer->GetYaxis()->GetBinUpEdge(k+1));
                sig_slice->SetTitle(newTitle.c_str());
                TFitResult tempResult = differential_crossection_fit(result, sig_slice, fit_func_sig, fit_func_bcg);
                if(tempResult.Chi2()==0){
                    result_vector[i*allCategories.size()+j]->SetBinContent(k+1, 0.);
                    result_vector[i*allCategories.size()+j]->SetBinError(k+1, 0.);
                    width_vector[i*allCategories.size()+j]->SetBinContent(k+1, 0.);
                    width_vector[i*allCategories.size()+j]->SetBinError(k+1, 0.);
                } else{
                    double bin_width = result_vector[i*allCategories.size()+j]->GetBinWidth(k+1);
                    par_value = abs(tempResult.Parameter(0)/sig_slice->GetXaxis()->GetBinWidth(1)*bin_width);
                    par_error = tempResult.ParError(0)/sig_slice->GetXaxis()->GetBinWidth(1)*bin_width;
                    result_vector[i*allCategories.size()+j]->SetBinContent(k+1, par_value);
                    result_vector[i*allCategories.size()+j]->SetBinError(k+1, par_error);
                    width_vector[i*allCategories.size()+j]->SetBinContent(k+1, tempResult.Parameter(2));
                    width_vector[i*allCategories.size()+j]->SetBinError(k+1, tempResult.Error(2));
                }
                Chi2withbcg_vector[i*allCategories.size()+j]->SetBinContent(k+1, tempResult.Chi2()/tempResult.Ndf());
            }
        }
    }
fittingBackground:
    //special fitting without the background
    printf("If you want to skip the fitting (excluding background), write \"yes\"\n");
    std::getline(std::cin, wannaSkip);
    if(wannaSkip.find("yes")!=std::string::npos){
        skippedFittingNoBackground = true;
        goto notFittingBackground;
    }
    for(size_t i = 0; i<pairTab.size(); i++){
        for(size_t j = 0; j<allCategories.size(); j++){
            TH2D* bcg_pointer = (TH2D*)background_vector[i*allCategories.size()+j]->Clone();
            TH2D* sig_pointer = (TH2D*)signal_vector[i*allCategories.size()+j]->Clone();
            //background fitting moved from binned to overall sum
            //background still removed like that
            int bcg_bin = bcg_pointer->GetXaxis()->FindBin(bcgRegion[i]);
            bcg_pointer->Scale(sig_pointer->Integral(bcg_bin, -1, 0, -1)/bcg_pointer->Integral(bcg_bin, -1, 0, -1));
            sig_pointer->Add(bcg_pointer, -1.);
            //for keeping title
            std::string baseOfTitle = std::string(sig_pointer->GetTitle())+" ";
            for(Int_t k = 0; k<sig_pointer->GetNbinsY(); k++){
                TH1D* sig_slice = sig_pointer->ProjectionX("_nobcg", k+1, k+1, "e");
                //fitting and filling result
                //setting lower range for m0-3*gamma
                //if lower bound is lower than m0-4*gamma
                //and 2*width otherwise
                double lower_range = fitMaximum[i]-3*fitWidth[i];
                if(sig_slice->GetBinLowEdge(GetFirstNonzeroBinNumber(sig_slice))>fitMaximum[i]-4*fitWidth[i]){
                    lower_range = fitMaximum[i]-2*fitWidth[i];
                }
                fit_func_sig->SetParameters(10., fitMaximum[i], fitWidth[i]);
                fit_func_sig->SetRange(lower_range, 1.1);
                fit_func_empty_bcg->SetRange(lower_range, 1.1);
                double par_value, par_error;
                std::string newTitle = baseOfTitle;
                newTitle += std::to_string(sig_pointer->GetYaxis()->GetBinLowEdge(k+1))+" - ";
                newTitle += std::to_string(sig_pointer->GetYaxis()->GetBinUpEdge(k+1));
                newTitle += " ZERO BACKGROUND";
                sig_slice->SetTitle(newTitle.c_str());
                TFitResult tempResult = differential_crossection_fit(result, sig_slice, fit_func_sig, fit_func_empty_bcg);
                if(tempResult.Chi2()==0){
                    result_vector_nobcgfit[i*allCategories.size()+j]->SetBinContent(k+1, 0.);
                    result_vector_nobcgfit[i*allCategories.size()+j]->SetBinError(k+1, 0.);
                    width_vector_nobcgfit[i*allCategories.size()+j]->SetBinContent(k+1, 0.);
                    width_vector_nobcgfit[i*allCategories.size()+j]->SetBinError(k+1, 0.);
                } else{
                    double bin_width = result_vector_nobcgfit[i*allCategories.size()+j]->GetBinWidth(k+1);
                    par_value = abs(tempResult.Parameter(0)/sig_slice->GetXaxis()->GetBinWidth(1)*bin_width);
                    par_error = tempResult.ParError(0)/sig_slice->GetXaxis()->GetBinWidth(1)*bin_width;
                    result_vector_nobcgfit[i*allCategories.size()+j]->SetBinContent(k+1, par_value);
                    result_vector_nobcgfit[i*allCategories.size()+j]->SetBinError(k+1, par_error);
                    width_vector_nobcgfit[i*allCategories.size()+j]->SetBinContent(k+1, tempResult.Parameter(2));
                    width_vector_nobcgfit[i*allCategories.size()+j]->SetBinError(k+1, tempResult.Error(2));
                }
                Chi2withoutbcg_vector[i*allCategories.size()+j]->SetBinContent(k+1, tempResult.Chi2()/tempResult.Ndf());
            }
        }
    }
notFittingBackground:
    //even more special fitting without removing the background
    printf("If you want to skip the fitting without removing background, write \"yes\"\n");
    std::getline(std::cin, wannaSkip);
    if(wannaSkip.find("yes")!=std::string::npos){
        skippedFittingNotRemovedBackground = true;
        goto noBackground;
    }
    for(size_t i = 0; i<pairTab.size(); i++){
        TF1* bcg_func;
        //pol2 background for Kstar (0,1)
        //custom for phi (2)
        if(i!=2){
            bcg_func = fit_func_for_Kstar;
            fit_func_sig->SetRange(0.75, 1.05);
            bcg_func->SetRange(0.75, 1.05);
        } else if(i==2){
            bcg_func = fit_func_paper_bcg;
            bcg_func->SetParameters(1, 1, 1);
            fit_func_sig->SetRange(1., 1.04);
            bcg_func->SetRange(1., 1.04);
        }
        //actual fitting
        for(size_t j = 0; j<allCategories.size(); j++){
            TH2D* sig_pointer = (TH2D*)signal_vector[i*allCategories.size()+j]->Clone();
            //for keeping title
            std::string baseOfTitle = std::string(sig_pointer->GetTitle())+" ";
            for(Int_t k = 0; k<sig_pointer->GetNbinsY(); k++){
                TH1D* sig_slice = sig_pointer->ProjectionX("_nobcg", k+1, k+1, "e");
                //fitting and filling result
                fit_func_sig->SetParameters(10., fitMaximum[i], fitWidth[i]);
                double par_value, par_error;
                std::string newTitle = baseOfTitle;
                newTitle += std::to_string(sig_pointer->GetYaxis()->GetBinLowEdge(k+1))+" - ";
                newTitle += std::to_string(sig_pointer->GetYaxis()->GetBinUpEdge(k+1));
                newTitle += " NOT REMOVED BACKGROUND";
                sig_slice->SetTitle(newTitle.c_str());
                TFitResult tempResult = differential_crossection_fit(result, sig_slice, fit_func_sig, bcg_func);
                if(tempResult.Chi2()==0){
                    result_vector_nobcgremoval[i*allCategories.size()+j]->SetBinContent(k+1, 0.);
                    result_vector_nobcgremoval[i*allCategories.size()+j]->SetBinError(k+1, 0.);
                    width_vector_nobcgremoval[i*allCategories.size()+j]->SetBinContent(k+1, 0.);
                    width_vector_nobcgremoval[i*allCategories.size()+j]->SetBinError(k+1, 0.);
                } else{
                    double bin_width = result_vector_nobcgremoval[i*allCategories.size()+j]->GetBinWidth(k+1);
                    par_value = abs(tempResult.Parameter(0)/sig_slice->GetXaxis()->GetBinWidth(1)*bin_width);
                    par_error = tempResult.ParError(0)/sig_slice->GetXaxis()->GetBinWidth(1)*bin_width;
                    result_vector_nobcgremoval[i*allCategories.size()+j]->SetBinContent(k+1, par_value);
                    result_vector_nobcgremoval[i*allCategories.size()+j]->SetBinError(k+1, par_error);
                    width_vector_nobcgremoval[i*allCategories.size()+j]->SetBinContent(k+1, tempResult.Parameter(2));
                    width_vector_nobcgremoval[i*allCategories.size()+j]->SetBinError(k+1, tempResult.Error(2));
                }
                Chi2withoutremovingbcg_vector[i*allCategories.size()+j]->SetBinContent(k+1, tempResult.Chi2()/tempResult.Ndf());
            }
        }
    }
noBackground:

    result->Close();
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
            if(!skippedFitting){
                draw_and_save(result_vector[i*allCategories.size()+j], folderWithDiagonal, "M"+pairTab[i]+allCategories[j], pairTab[i]+" "+allCategories[j]+";"+allCategories[j]+";entries", "e");
                draw_and_save(width_vector[i*allCategories.size()+j], folderWithDiagonal, "M"+pairTab[i]+allCategories[j]+"Width", pairTab[i]+" "+allCategories[j]+";"+allCategories[j]+";entries", "e");
            }
            if(!skippedFittingNoBackground){
                draw_and_save(result_vector_nobcgfit[i*allCategories.size()+j], folderWithDiagonal, "M"+pairTab[i]+allCategories[j]+"nobcgfit", pairTab[i]+" "+allCategories[j]+";"+allCategories[j]+";entries", "e");
                draw_and_save(width_vector_nobcgfit[i*allCategories.size()+j], folderWithDiagonal, "M"+pairTab[i]+allCategories[j]+"Width"+"nobcgfit", pairTab[i]+" "+allCategories[j]+";"+allCategories[j]+";entries", "e");
            }
            if(!skippedFittingNotRemovedBackground){
                draw_and_save(result_vector_nobcgremoval[i*allCategories.size()+j], folderWithDiagonal, "M"+pairTab[i]+allCategories[j]+"nobcgremoval", pairTab[i]+" "+allCategories[j]+";"+allCategories[j]+";entries", "e");
                draw_and_save(width_vector_nobcgremoval[i*allCategories.size()+j], folderWithDiagonal, "M"+pairTab[i]+allCategories[j]+"Width"+"nobcgremoval", pairTab[i]+" "+allCategories[j]+";"+allCategories[j]+";entries", "e");
            }
        }
    }
    //Chi2 drawing
    for(size_t i = 0; i<pairTab.size(); i++){
        for(size_t j = 0; j<allCategories.size(); j++){
            MyStyles styleLibrary;
            TStyle tempStyle = styleLibrary.Hist2DNormalSize(true);
            tempStyle.cd();
            gROOT->ForceStyle();
            TCanvas* resultCanvas = new TCanvas("resultCanvas", "resultCanvas", 4000, 2400);
            Chi2withbcg_vector[i*allCategories.size()+j]->SetMinimum(0.);
            Chi2withbcg_vector[i*allCategories.size()+j]->SetMaximum(1.1*std::max(std::max(Chi2withbcg_vector[i*allCategories.size()+j]->GetMaximum(), Chi2withoutbcg_vector[i*allCategories.size()+j]->GetMaximum()), Chi2withoutremovingbcg_vector[i*allCategories.size()+j]->GetMaximum()));
            Chi2withbcg_vector[i*allCategories.size()+j]->Draw("hist");
            Chi2withoutbcg_vector[i*allCategories.size()+j]->Draw("hist same");
            Chi2withoutremovingbcg_vector[i*allCategories.size()+j]->Draw("hist same");
            resultCanvas->UseCurrentStyle();
            Chi2withbcg_vector[i*allCategories.size()+j]->SetTitle((pairTab[i]+" "+allCategories[j]+";"+allCategories[j]+";#chi^{2}/ndf").c_str());
            Chi2withbcg_vector[i*allCategories.size()+j]->SetLineColor(kBlue);
            Chi2withoutbcg_vector[i*allCategories.size()+j]->SetLineColor(kRed);
            Chi2withoutremovingbcg_vector[i*allCategories.size()+j]->SetLineColor(kGreen+2);
            resultCanvas->BuildLegend();
            resultCanvas->Update();
            resultCanvas->SaveAs((folderWithDiagonal+"M"+pairTab[i]+allCategories[j]+"Chi2.pdf").c_str());
        }
    }

    //the end
    theApp.Run();
    return 0;
}

int GetFirstNonzeroBinNumber(TH1* input){
    for(int i = 1; i<=input->GetNbinsX(); i++){
        if(input->GetBinContent(i)!=0.0)
            return i;
    }
    return -1;
}

void draw_and_save(TH1D* data, std::string folderWithDiagonal, std::string name, std::string title, std::string options){
    MyStyles styleLibrary;
    TStyle tempStyle = styleLibrary.Hist2DQuarterSize(true);
    tempStyle.cd();
    gROOT->ForceStyle();
    TCanvas* resultCanvas = new TCanvas("resultCanvas", "resultCanvas", 4000, 2400);
    //TODO: add fitting
    data->Draw(options.c_str());
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
    if(data->Integral(bcg_bin, -1)==0 or bcg->Integral(bcg_bin, -1)){
        data->SetTitle((std::string(data->GetTitle())+" no background  removed").c_str());
    } else{
        bcg->Scale(data->Integral(bcg_bin, -1)/bcg->Integral(bcg_bin, -1));
        data->Add(bcg, -1.);
    }
    data->Draw("e");
    resultCanvas->UseCurrentStyle();
    data->SetMarkerStyle(kFullCircle);
    data->SetMarkerSize(2);
    data->SetMarkerColor(kBlue);
    data->SetTitle(title.c_str());

    resultCanvas->SaveAs((folderWithDiagonal+name+".pdf").c_str());
}

TFitResult differential_crossection_fit(TPad* pad, TH1D* slice, TF1* fitting_function_signal, TF1* fitting_function_bcg, std::string draw_full_path){
    //setting up the functions
    pad->Clear();
    gROOT->SetSelectedPad(pad);
    gStyle->SetOptStat(0);
    gStyle->SetHistMinimumZero();
    gStyle->SetOptFit();
    gStyle->SetStatBorderSize(0.);
    gStyle->SetStatX(0.94);
    gStyle->SetStatY(0.89);
    gStyle->SetStatW(0.18);
    gStyle->SetFitFormat("6.5g");
    gStyle->SetStatColor(0);
    gPad->SetMargin(0.1, 0.05, 0.1, 0.1);
    TFitResultPtr fitPointer;
    int signalparams, bcgparams, totalparams;
    double rangemin, rangemax;
    signalparams = fitting_function_signal->GetNpar();
    bcgparams = fitting_function_bcg->GetNpar();
    totalparams = signalparams+bcgparams;
    fitting_function_signal->GetRange(rangemin, rangemax);
    auto fitting_function_sum = [&](double* x, double* par){
        return fitting_function_signal->EvalPar(x, par)+fitting_function_bcg->EvalPar(x, par+signalparams);
    };
    TF1* fitting_function_total = new TF1("fitting_function_total", fitting_function_sum, rangemin, rangemax, totalparams);
    double lowlim, highlim;
    for(int i = 0; i<signalparams; i++){
        fitting_function_signal->GetParLimits(i, lowlim, highlim);
        if(lowlim==highlim&&lowlim*highlim!=0){
            fitting_function_total->FixParameter(i, lowlim);
        } else{
            fitting_function_total->SetParLimits(i, lowlim, highlim);
        }
        //if not default name (beginning with p), copy it
        if(std::string(fitting_function_signal->GetParName(i))[0]!='p')
            fitting_function_total->SetParName(i, fitting_function_signal->GetParName(i));
    }
    for(int i = 0; i<bcgparams; i++){
        fitting_function_bcg->GetParLimits(i, lowlim, highlim);
        if(lowlim==highlim&&lowlim*highlim!=0){
            fitting_function_total->FixParameter(i+signalparams, lowlim);
        } else{
            fitting_function_total->SetParLimits(i+signalparams, lowlim, highlim);
        }
        //if not default name (beginning with p), copy it
        if(std::string(fitting_function_bcg->GetParName(i))[0]!='p')
            fitting_function_total->SetParName(i+signalparams, fitting_function_bcg->GetParName(i));
    }
    //old data, the base of drawing
    slice->SetMarkerStyle(kFullCircle);
    slice->SetMarkerColor(kBlue);
    slice->SetLineColor(kBlue+2);
    //setting beginning params for function
    Double_t params[totalparams];
    fitting_function_signal->GetParameters(params);
    fitting_function_bcg->GetParameters(params+signalparams);
    fitting_function_total->SetParameters(params);
    fitting_function_bcg->SetLineStyle(3);
    fitting_function_bcg->SetLineWidth(3);
    fitting_function_signal->SetLineStyle(3);
    fitting_function_signal->SetLineWidth(3);
    fitting_function_signal->SetNpx(1000);
    fitting_function_total->SetNpx(1000);
    //fitting and drawing
    fitPointer = slice->Fit(fitting_function_total, "0BRS");
    slice->Draw("E");
    fitting_function_signal->SetParameters(fitting_function_total->GetParameters());
    fitting_function_signal->SetParErrors(fitting_function_total->GetParErrors());
    fitting_function_bcg->SetParameters(fitting_function_total->GetParameters()+signalparams);
    TF1* fitting_function_total_COPY = fitting_function_total->DrawCopy("CSAME");
    TF1* fitting_function_bcg_COPY = fitting_function_bcg->DrawCopy("CSAME");
    TF1* fitting_function_signal_COPY = fitting_function_signal->DrawCopy("CSAME");
    //drawing prameters
    gPad->Update();

    //interactive part
    std::cout<<"Enter which parameter and how much change"<<std::endl;
    std::cout<<"Writing \"abs\" before just sets the parameter"<<std::endl;
    std::cout<<"Or write \"fit\" to fit"<<std::endl;
    std::cout<<"Or write \"zero\" if there is no data suitable to fit"<<std::endl;
    std::cout<<"Or press \"Enter\" if everything ok"<<std::endl;
    std::string response;
    //needed to read whole line, not just until first space
    std::getline(std::cin, response);
    while(response.size()!=0&&response.find("zero")==std::string::npos){
        //fitting/changing part
        if(response.find("fit")!=std::string::npos){
            fitting_function_total->SetParameters(fitting_function_total_COPY->GetParameters());
            fitPointer = slice->Fit(fitting_function_total, "0BRS");
            printf("Chi2 from fit: %f, ndof = %d\n", fitting_function_total->GetChisquare(), fitting_function_total->GetNDF());
            fitting_function_total_COPY->SetParameters(fitting_function_total->GetParameters());
            fitting_function_bcg_COPY->SetParameters(fitting_function_total->GetParameters()+signalparams);
            fitting_function_signal_COPY->SetParameters(fitting_function_total->GetParameters());
            fitting_function_signal_COPY->SetParErrors(fitting_function_total->GetParErrors());
        } else if(response.find("zero")!=std::string::npos){
            fitting_function_signal_COPY->SetParameter(0, 0);
            fitting_function_signal_COPY->SetParError(0, 0);
        } else{
            int coeff_number;
            double change, relative;
            if(response.find("abs")==std::string::npos){
                relative = 1.;
                sscanf(response.c_str(), "%d%lf", &coeff_number, &change);
            } else{
                relative = 0.;
                sscanf(response.c_str(), "abs %d%lf", &coeff_number, &change);
            }
            fitting_function_total_COPY->SetParameter(coeff_number, fitting_function_total_COPY->GetParameter(coeff_number)*relative+change);
            if(coeff_number>=signalparams)
                fitting_function_bcg_COPY->SetParameter(coeff_number-signalparams, fitting_function_bcg_COPY->GetParameter(coeff_number-signalparams)*relative+change);
            else
                fitting_function_signal_COPY->SetParameter(coeff_number, fitting_function_signal_COPY->GetParameter(coeff_number)*relative+change);
        }
        //drawing part
        fitting_function_total_COPY->Draw("CSAME");
        fitting_function_bcg_COPY->Draw("CSAME");
        fitting_function_signal_COPY->Draw("CSAME");
        gPad->Update();
        //interactive part
        std::cout<<"Chi2 after edit/fit: "<<slice->Chisquare(fitting_function_total_COPY, "R")<<", ndof: "<<fitting_function_total_COPY->GetNDF()<<std::endl;
        std::cout<<"Enter which parameter and how much change"<<std::endl;
        std::cout<<"Writing \"abs\" before just sets the parameter"<<std::endl;
        std::cout<<"Or write \"fit\" to fit"<<std::endl;
        std::cout<<"Or write \"zero\" if there is no data suitable to fit"<<std::endl;
        std::cout<<"Or press \"Enter\" if everything ok"<<std::endl;
        //needed to read whole line, not just until first space
        std::getline(std::cin, response);
    }
    //finishing touches
    printf("Accepted current parameters\n");
    if(draw_full_path.size()!=0){
        pad->SaveAs(draw_full_path.c_str());
    }
    gStyle->SetOptStat(1);
    //giving result
    //if there is no particles to fit
    //then chi2=0 is set as a sign
    //there was a problem with zero data to fit ans TFitResultPtr not having... a result
    //so this is workaround
    if(response.find("zero")!=std::string::npos){
        TFitResult zeroResult;
        zeroResult.SetChi2AndNdf(0., 1.);
        return zeroResult;
    }
    return *(fitPointer.Get());
}