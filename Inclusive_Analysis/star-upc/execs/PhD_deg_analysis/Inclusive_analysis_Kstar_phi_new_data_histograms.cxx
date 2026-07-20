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
#include "THStack.h"
#include "TLine.h"
#include "TLegend.h"

#include "MyStyles.h"

std::string rountToNSignificantFigures(double input, int n = 2);
int GetFirstNonzeroBinNumber(TH1* input);
void set_background_fitting(TCanvas* canvas, TH1D* data, TH1D* bcg, double& bcgRegionStart, double& bcgRegionStop, std::string name, std::string title, std::string draw_full_path = "");
void draw_and_save(TH1D* data, std::string folderWithDiagonal, std::string name, std::string title, std::string options = "");
void draw_and_save_minus_background(TH1D* data, TH1D* bcg, std::string folderWithDiagonal, std::string name, std::string title, double bcg_region);
void draw_bulk(std::vector<TH1D*> data, std::string folderWithDiagonal, std::string name, std::string title, std::string options);
TFitResult differential_crossection_fit(TPad* pad, TH1D* slice, TF1* fitting_function_signal, TF1* fitting_function_bcg, std::string draw_full_path = "");

int main(int argc, char* argv[]){
    //something used so that the histograms would draw
    //https://stackoverflow.com/questions/30932725/painting-a-tcanvas-to-the-screen-in-a-compiled-root-cern-application
    TApplication theApp("App", &argc, argv);
    argv = theApp.Argv();
    argc = theApp.Argc();

    //input file
    TFile* input = TFile::Open(static_cast<const char*>(argv[1]));
    printf("Type of background (3rd argument):\n");
    printf("1 - same-sign\n");
    printf("2 - track rotation\n");
    printf("3 - random track rotation\n");
    printf("4 - mixed-event\n");
    std::string BcgType;
    if(argc<4){
        printf("Background not chosen\n");
        return 1;
    }
    switch(atoi(argv[3])){
    case 1:
        BcgType = "SameSign";
        break;
    case 2:
        BcgType = "TrackRotation";
        break;
    case 3:
        BcgType = "RandomTrackRotation";
        break;
    case 4:
        BcgType = "MixedEvent";
        break;
    default:
        printf("Background not chosen\n");
        return 1;
        break;
    }
    printf("Background option chosen: %d (%s)\n", atoi(argv[3]), BcgType.c_str());

    //getting histograms out
    TH1D* pairInfoSignal = (TH1D*)input->Get("pairInfoSignal");
    TH1D* pairInfoBackgroundSameSign = (TH1D*)input->Get("pairInfoBackgroundSameSign");
    TH1D* pairInfoBackgroundTrackRotation = (TH1D*)input->Get("pairInfoBackgroundTrackRotation");
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
    //getting background histograms out
    TH1D* MKpiChi2bcg = (TH1D*)input->Get(("MKpiChi2Bcg"+BcgType).c_str());
    TH1D* MpiKChi2bcg = (TH1D*)input->Get(("MpiKChi2Bcg"+BcgType).c_str());
    TH1D* MppiChi2bcg = (TH1D*)input->Get(("MppiChi2Bcg"+BcgType).c_str());
    TH1D* MpipChi2bcg = (TH1D*)input->Get(("MpipChi2Bcg"+BcgType).c_str());
    TH1D* MKKChi2bcg = (TH1D*)input->Get(("MKKChi2Bcg"+BcgType).c_str());
    TH1D* MpipiChi2bcg = (TH1D*)input->Get(("MpipiChi2Bcg"+BcgType).c_str());
    TH1D* MppChi2bcg = (TH1D*)input->Get(("MppChi2Bcg"+BcgType).c_str());
    //getting diffractive background and signal
    std::vector<std::string> pairTab = { "Kpi", "piK", "KK" };
    std::vector<double> bcgRegionStart = { 1.05, 1.05, 1.8 };
    std::vector<double> bcgRegionStop = { 1.25, 1.25, 2.4 };
    std::vector<double> fitMaximum = { 0.892, 0.892, 1.02 };
    std::vector<double> fitWidth = { 0.0514, 0.0514, 0.00443 };//51.4 MeV for K*(892), 4.43 MeV for phi(1020)
    //creating list of categories
    std::vector<std::string> allCategories;
    std::ifstream infile("STAR-Analysis/Inclusive_Analysis/star-upc/execs/PhD_deg_analysis/Differential_crossection_values.txt");
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
    std::string folderWithDiagonal = std::string(static_cast<const char*>(argv[2]));
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
    std::string tempBackgroundName = "M$Chi2Bcg"+BcgType;
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
            background_vector.push_back((TH2D*)input->Get((tempHistName+category).c_str()));
            //result, but uses signal template
            result_vector.push_back(new TH1D((tempHistName+"Result"+category).c_str(), (pair+" "+category+" bcg").c_str(), signal_vector.back()->GetNbinsY(), signal_vector.back()->GetYaxis()->GetXmin(), signal_vector.back()->GetYaxis()->GetXmax()));
            result_vector_nobcgfit.push_back(new TH1D((tempHistName+"Result"+category).c_str(), (pair+" "+category+" no fitted bcg").c_str(), signal_vector.back()->GetNbinsY(), signal_vector.back()->GetYaxis()->GetXmin(), signal_vector.back()->GetYaxis()->GetXmax()));
            result_vector_nobcgremoval.push_back(new TH1D((tempHistName+"Result"+category).c_str(), (pair+" "+category+" not removed bcg").c_str(), signal_vector.back()->GetNbinsY(), signal_vector.back()->GetYaxis()->GetXmin(), signal_vector.back()->GetYaxis()->GetXmax()));
            //width of the fits
            width_vector.push_back(new TH1D((tempHistName+"Width"+category).c_str(), (pair+" "+category+" bcg").c_str(), signal_vector.back()->GetNbinsY(), signal_vector.back()->GetYaxis()->GetXmin(), signal_vector.back()->GetYaxis()->GetXmax()));
            width_vector_nobcgfit.push_back(new TH1D((tempHistName+"Width"+category).c_str(), (pair+" "+category+" no fitted bcg").c_str(), signal_vector.back()->GetNbinsY(), signal_vector.back()->GetYaxis()->GetXmin(), signal_vector.back()->GetYaxis()->GetXmax()));
            width_vector_nobcgremoval.push_back(new TH1D((tempHistName+"Width"+category).c_str(), (pair+" "+category+" not removed bcg").c_str(), signal_vector.back()->GetNbinsY(), signal_vector.back()->GetYaxis()->GetXmin(), signal_vector.back()->GetYaxis()->GetXmax()));
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
    TCanvas* result = MyStyles::DefaultCanvas("result");
    MyStyles styleLibrary;
    TStyle mystyle = styleLibrary.Hist2DDisplay(false);
    mystyle.cd();
    result->UseCurrentStyle();
    // gROOT->ForceStyle();
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
    fit_func_custom_bcg.SetParNames("a_{0}", "a_{1}", "a_{2}");
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

    //checking fitting bounds
    TStyle mystyle2 = styleLibrary.Hist2DDisplay(true);
    mystyle2.cd();
    result->UseCurrentStyle();
    set_background_fitting(result, MKpiChi2, MKpiChi2bcg, bcgRegionStart[0], bcgRegionStop[0], "KpiRatio", "K^{+}#pi^{-} Background/Signal ratio with "+BcgType+" background", folderWithDiagonal+"MKpiRatio.pdf");
    mystyle2.cd();
    result->UseCurrentStyle();
    set_background_fitting(result, MpiKChi2, MpiKChi2bcg, bcgRegionStart[1], bcgRegionStop[1], "piKRatio", "#pi^{+}K^{-} Background/Signal ratio with "+BcgType+" background", folderWithDiagonal+"MpiKRatio.pdf");
    mystyle2.cd();
    result->UseCurrentStyle();
    set_background_fitting(result, MKKChi2, MKKChi2bcg, bcgRegionStart[2], bcgRegionStop[2], "KKRatio", "K^{+}K^{-} Background/Signal ratio with "+BcgType+" background", folderWithDiagonal+"MKKRatio.pdf");

    //fitting functions
    TF1* fit_func_sig = new TF1("fit_func_sig", "breitwigner", 0.8, 1.0);
    fit_func_sig->SetParNames("N", "m_{0}", "#Gamma");
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
    fit_func_bcg->SetParNames("bcg(x_{1})", "bcg(x_{2})");
    TF1* fit_func_empty_bcg = new TF1("fit_func_empty_bcg", "0.", 0.8, 1.0);
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
    fit_func_paper_bcg->SetParNames("A", "B", "C");
    // TF1* fit_func_for_Kstar = new TF1("fit_func_for_Kstar", "pol3", 0.8, 1.0);
    // fit_func_for_Kstar->SetParNames("a_{0}", "a_{1}", "a_{2}", "a_{3}");
    auto weibull = [&](double* x, double* p){
        double N = p[0];
        double k = p[1];
        double lambda = p[2];
        double x0 = p[3];
        double xi = x[0]-x0;
        if(xi<=0){
            return 0.;
        }
        return N*k/lambda*pow(xi/lambda, k-1)*exp(-pow(xi/lambda, k));
    };
    TF1* fit_func_for_Kstar = new TF1("fit_func_for_Kstar", weibull, 0.8, 1.0, 4);
    fit_func_for_Kstar->SetParNames("N_{bcg}", "k", "#lambda", "x_{0}");
    fit_func_for_Kstar->SetParLimits(0, 0, 100000);
    fit_func_for_Kstar->SetParLimits(1, 0, 100000);
    fit_func_for_Kstar->SetParLimits(2, 0, 100000);
    fit_func_for_Kstar->SetParLimits(3, 0.6, 1.0);
    TF1* fit_func_for_phi = new TF1("fit_func_for_phi", "pol2", 0.8, 1.0);
    fit_func_for_phi->SetParNames("a_{0}", "a_{1}", "a_{2}");
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
            int bcg_bin_start = bcg_pointer->GetXaxis()->FindBin(bcgRegionStart[i]);
            int bcg_bin_stop = bcg_pointer->GetXaxis()->FindBin(bcgRegionStop[i]);
            bcg_pointer->Scale(sig_pointer->Integral(bcg_bin_start, bcg_bin_stop, 0, -1)/bcg_pointer->Integral(bcg_bin_start, bcg_bin_stop, 0, -1));
            sig_pointer->Add(bcg_pointer, -1.);
            //for keeping title
            std::string baseOfTitle = std::string(sig_pointer->GetTitle())+" ";
            for(Int_t k = 0; k<sig_pointer->GetNbinsY(); k++){
                TH1D* sig_slice = sig_pointer->ProjectionX("_bcg", k+1, k+1, "e1");
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
                newTitle += "("+rountToNSignificantFigures(sig_pointer->GetYaxis()->GetBinLowEdge(k+1))+", ";
                newTitle += rountToNSignificantFigures(sig_pointer->GetYaxis()->GetBinUpEdge(k+1))+")";
                sig_slice->SetTitle(newTitle.c_str());
                TFitResult tempResult = differential_crossection_fit(result, sig_slice, fit_func_sig, fit_func_bcg);
                if(tempResult.Chi2()==0){
                    result_vector[i*allCategories.size()+j]->SetBinContent(k+1, 0.);
                    result_vector[i*allCategories.size()+j]->SetBinError(k+1, 0.);
                    width_vector[i*allCategories.size()+j]->SetBinContent(k+1, 0.);
                    width_vector[i*allCategories.size()+j]->SetBinError(k+1, 0.);
                } else{
                    double bin_width = result_vector[i*allCategories.size()+j]->GetBinWidth(k+1);
                    par_value = fabs(tempResult.Parameter(0)/sig_slice->GetXaxis()->GetBinWidth(1)*bin_width);
                    par_error = tempResult.ParError(0)/sig_slice->GetXaxis()->GetBinWidth(1)*bin_width;
                    result_vector[i*allCategories.size()+j]->SetBinContent(k+1, par_value);
                    result_vector[i*allCategories.size()+j]->SetBinError(k+1, par_error);
                    width_vector[i*allCategories.size()+j]->SetBinContent(k+1, fabs(tempResult.Parameter(2)));
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
            int bcg_bin_start = bcg_pointer->GetXaxis()->FindBin(bcgRegionStart[i]);
            int bcg_bin_stop = bcg_pointer->GetXaxis()->FindBin(bcgRegionStop[i]);
            bcg_pointer->Scale(sig_pointer->Integral(bcg_bin_start, bcg_bin_stop, 0, -1)/bcg_pointer->Integral(bcg_bin_start, bcg_bin_stop, 0, -1));
            sig_pointer->Add(bcg_pointer, -1.);
            //for keeping title
            std::string baseOfTitle = std::string(sig_pointer->GetTitle())+" ";
            for(Int_t k = 0; k<sig_pointer->GetNbinsY(); k++){
                TH1D* sig_slice = sig_pointer->ProjectionX("_nobcg", k+1, k+1, "e1");
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
                newTitle += "("+rountToNSignificantFigures(sig_pointer->GetYaxis()->GetBinLowEdge(k+1))+", ";
                newTitle += rountToNSignificantFigures(sig_pointer->GetYaxis()->GetBinUpEdge(k+1))+")";
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
                    par_value = fabs(tempResult.Parameter(0)/sig_slice->GetXaxis()->GetBinWidth(1)*bin_width);
                    par_error = tempResult.ParError(0)/sig_slice->GetXaxis()->GetBinWidth(1)*bin_width;
                    result_vector_nobcgfit[i*allCategories.size()+j]->SetBinContent(k+1, par_value);
                    result_vector_nobcgfit[i*allCategories.size()+j]->SetBinError(k+1, par_error);
                    width_vector_nobcgfit[i*allCategories.size()+j]->SetBinContent(k+1, fabs(tempResult.Parameter(2)));
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
            // bcg_func = fit_func_paper_bcg;
            // bcg_func->SetParameters(1, 1, 1);
            bcg_func = fit_func_for_phi;
            fit_func_sig->SetRange(1., 1.04);
            bcg_func->SetRange(1., 1.04);
        }
        //actual fitting
        for(size_t j = 0; j<allCategories.size(); j++){
            TH2D* sig_pointer = (TH2D*)signal_vector[i*allCategories.size()+j]->Clone();
            //for keeping title
            std::string baseOfTitle = std::string(sig_pointer->GetTitle())+" ";
            for(Int_t k = 0; k<sig_pointer->GetNbinsY(); k++){
                TH1D* sig_slice = sig_pointer->ProjectionX("_nobcg", k+1, k+1, "e1");
                //if for the case of small tail
                if(i!=2&&k==0&&j==0){
                    fit_func_sig->SetRange(0.8, 1.05);
                    bcg_func->SetRange(0.8, 1.05);
                } else if(i!=2){
                    fit_func_sig->SetRange(0.75, 1.05);
                    bcg_func->SetRange(0.75, 1.05);
                }
                //prefitting linear background seed for  the polynomial
                fit_func_sig->SetParameters(10., fitMaximum[i], fitWidth[i]);
                //TODO: fix
                // double y1, y2;
                // y1 = sig_slice->GetBinContent(sig_slice->GetXaxis()->FindBin(0.8));
                // y2 = sig_slice->GetBinContent(sig_slice->GetXaxis()->FindBin(0.8));
                // bcg_func->SetParameters((y1+y2)/2, (y2-y1)/0.2, 0., 0.);
                if(i!=2){
                    //"N_{bcg}", "k", "#lambda", "x_{0}"
                    bcg_func->SetParameters(110., 2.1, 0.26, 0.795);
                } else{
                    double y1, y2;
                    y1 = sig_slice->GetBinContent(sig_slice->GetXaxis()->FindBin(0.8));
                    y2 = sig_slice->GetBinContent(sig_slice->GetXaxis()->FindBin(0.8));
                    bcg_func->SetParameters((y1+y2)/2, (y2-y1)/0.2, 0., 0.);
                }
                //fitting and filling result
                double par_value, par_error;
                std::string newTitle = baseOfTitle;
                newTitle += "("+rountToNSignificantFigures(sig_pointer->GetYaxis()->GetBinLowEdge(k+1))+", ";
                newTitle += rountToNSignificantFigures(sig_pointer->GetYaxis()->GetBinUpEdge(k+1))+")";
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
                    par_value = fabs(tempResult.Parameter(0)/sig_slice->GetXaxis()->GetBinWidth(1)*bin_width);
                    par_error = tempResult.ParError(0)/sig_slice->GetXaxis()->GetBinWidth(1)*bin_width;
                    result_vector_nobcgremoval[i*allCategories.size()+j]->SetBinContent(k+1, par_value);
                    result_vector_nobcgremoval[i*allCategories.size()+j]->SetBinError(k+1, par_error);
                    width_vector_nobcgremoval[i*allCategories.size()+j]->SetBinContent(k+1, fabs(tempResult.Parameter(2)));
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

    printf("Fitting was done, now drawing begins\nFirst, pairInfo\n");

    //drawing
    //pairInfo
    TCanvas* resultCanvas = new TCanvas("resultCanvas", "resultCanvas", 4000, 2400);
    //signal
    pairInfoSignal->SetMinimum(0);
    pairInfoSignal->SetLineColor(kBlue+2);
    pairInfoSignal->SetMarkerSize(2);
    //TODO: add fitting
    // pairInfoSignal->Draw("E");
    pairInfoSignal->Draw("hist");
    pairInfoSignal->Draw("same text0");
    pairInfoSignal->SetLineWidth(2);
    pairInfoSignal->GetXaxis()->SetLabelSize(0.06);
    pairInfoSignal->GetXaxis()->SetTitleSize(0.06);
    resultCanvas->SetTopMargin(0.05);
    gPad->Update();
    resultCanvas->SaveAs((folderWithDiagonal+"pairInfoSignal.pdf").c_str());
    resultCanvas->Clear();
    //background (same-sign)
    pairInfoBackgroundSameSign->SetMinimum(0);
    pairInfoBackgroundSameSign->SetLineColor(kBlue+2);
    pairInfoBackgroundSameSign->SetMarkerSize(2);
    //TODO: add fitting
    // pairInfoBackgroundSameSign->Draw("E");
    pairInfoBackgroundSameSign->Draw("hist");
    pairInfoBackgroundSameSign->Draw("same text0");
    pairInfoBackgroundSameSign->SetLineWidth(2);
    pairInfoBackgroundSameSign->GetXaxis()->SetLabelSize(0.06);
    pairInfoBackgroundSameSign->GetXaxis()->SetTitleSize(0.06);
    resultCanvas->SetTopMargin(0.05);
    gPad->Update();
    resultCanvas->SaveAs((folderWithDiagonal+"pairInfoBackgroundSameSign.pdf").c_str());
    resultCanvas->Clear();
    //background (track rotation)
    pairInfoBackgroundTrackRotation->SetMinimum(0);
    pairInfoBackgroundTrackRotation->SetLineColor(kBlue+2);
    pairInfoBackgroundTrackRotation->SetMarkerSize(2);
    //TODO: add fitting
    // pairInfoBackgroundTrackRotation->Draw("E");
    pairInfoBackgroundTrackRotation->Draw("hist");
    pairInfoBackgroundTrackRotation->Draw("same text0");
    pairInfoBackgroundTrackRotation->SetLineWidth(2);
    pairInfoBackgroundTrackRotation->GetXaxis()->SetLabelSize(0.06);
    pairInfoBackgroundTrackRotation->GetXaxis()->SetTitleSize(0.06);
    resultCanvas->SetTopMargin(0.05);
    gPad->Update();
    resultCanvas->SaveAs((folderWithDiagonal+"pairInfoBackgroundTrackRotation.pdf").c_str());
    delete resultCanvas;

    printf("Second, the rest\n");

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

    //drawing results - chi2, number of detected decays and width of resonances
    for(size_t i = 0; i<pairTab.size(); i++){
        for(size_t j = 0; j<allCategories.size(); j++){
            //Chi2
            std::vector<TH1D*> AllvectorChi2;
            AllvectorChi2.push_back(Chi2withbcg_vector[i*allCategories.size()+j]);
            AllvectorChi2.push_back(Chi2withoutbcg_vector[i*allCategories.size()+j]);
            AllvectorChi2.push_back(Chi2withoutremovingbcg_vector[i*allCategories.size()+j]);
            draw_bulk(AllvectorChi2, folderWithDiagonal, "M"+pairTab[i]+allCategories[j]+"Chi2", pairTab[i]+" "+allCategories[j]+" #chi^{2}/ndf;"+allCategories[j]+";#chi^{2}/ndf", "hist");
            //results & width
            std::vector<TH1D*> AllvectorResults, AllvectorWidths;
            if(!skippedFitting){
                AllvectorResults.push_back(result_vector[i*allCategories.size()+j]);
                AllvectorWidths.push_back(width_vector[i*allCategories.size()+j]);
            }
            if(!skippedFittingNoBackground){
                AllvectorResults.push_back(result_vector_nobcgfit[i*allCategories.size()+j]);
                AllvectorWidths.push_back(width_vector_nobcgfit[i*allCategories.size()+j]);
            }
            if(!skippedFittingNotRemovedBackground){
                AllvectorResults.push_back(result_vector_nobcgremoval[i*allCategories.size()+j]);
                AllvectorWidths.push_back(width_vector_nobcgremoval[i*allCategories.size()+j]);
            }
            draw_bulk(AllvectorResults, folderWithDiagonal, "M"+pairTab[i]+allCategories[j]+"Results", pairTab[i]+" "+allCategories[j]+" result;"+allCategories[j]+";entries", "e1");
            draw_bulk(AllvectorWidths, folderWithDiagonal, "M"+pairTab[i]+allCategories[j]+"Widths", pairTab[i]+" "+allCategories[j]+" width;"+allCategories[j]+";#Gamma (GeV)", "e1");
        }
    }

    //the end
    theApp.Run();
    return 0;
}

std::string rountToNSignificantFigures(double input, int n){
    char buffer[100];
    sprintf(buffer, ("%."+std::to_string(n)+"f").c_str(), input);
    return std::string(buffer);
}

int GetFirstNonzeroBinNumber(TH1* input){
    for(int i = 1; i<=input->GetNbinsX(); i++){
        if(input->GetBinContent(i)!=0.0)
            return i;
    }
    return -1;
}

void set_background_fitting(TCanvas* canvas, TH1D* data, TH1D* bcg, double& bcgRegionStart, double& bcgRegionStop, std::string name, std::string title, std::string draw_full_path){
    //preparing for drawing
    canvas->Clear();
    printf("%d\n", gStyle->GetLegendFont());
    //preparing and drawing ratio plot
    TH1D SigBcgRatio(*bcg);
    SigBcgRatio.SetMarkerStyle(kFullCircle);
    SigBcgRatio.SetMarkerColor(kBlue);
    SigBcgRatio.SetLineColor(kBlue+2);
    data->Sumw2();
    SigBcgRatio.Divide(data);
    SigBcgRatio.SetNameTitle(name.c_str(), title.c_str());
    SigBcgRatio.Rebin(2);
    SigBcgRatio.Draw("e1");
    canvas->Update();
    //drawing two lines at proper coordinates
    TLine line1(bcgRegionStart, canvas->GetUymin(), bcgRegionStart, canvas->GetUymax());
    TLine line2(bcgRegionStop, canvas->GetUymin(), bcgRegionStop, canvas->GetUymax());
    line1.SetLineColor(kRed);
    line1.SetLineStyle(kDashed);
    line1.SetLineWidth(2);
    line1.Draw("same");
    line2.SetLineColor(kRed);
    line2.SetLineStyle(kDashed);
    line2.SetLineWidth(2);
    line2.Draw("same");
    //drawing the legend
    TLegend legend_for_background_fitting(0.56, 0.2, 0.89, 0.39);
    legend_for_background_fitting.SetTextSize(0.04);
    legend_for_background_fitting.AddEntry(&SigBcgRatio, "Background to total ratio");
    legend_for_background_fitting.AddEntry(&line1, "Tested background region", "l");
    legend_for_background_fitting.SetBorderSize(0);
    legend_for_background_fitting.DrawClone("SAME");
    std::string inputString;
    printf("Write lower & higher bound of the background fit and press \"Enter\"\n");
    printf("Or just press \"Enter\" to keep current: %lf, %lf\n", bcgRegionStart, bcgRegionStop);
    do{
        sscanf(inputString.c_str(), "%lf %lf", &bcgRegionStart, &bcgRegionStop);
        printf("Chosen bounds for \"%s\" histogram: %lf, %lf\n", name.c_str(), bcgRegionStart, bcgRegionStop);
        line1.SetX1(bcgRegionStart);
        line1.SetX2(bcgRegionStart);
        line2.SetX1(bcgRegionStop);
        line2.SetX2(bcgRegionStop);
        canvas->ModifiedUpdate();
        std::getline(std::cin, inputString);
    } while(inputString.size()!=0);
    //saving and clearing the canvas
    if(draw_full_path.size()!=0){
        //setting new style
        MyStyles styleLibrary;
        TStyle mystyle = styleLibrary.Hist2DQuarterSize();
        mystyle.cd();
        canvas->UseCurrentStyle();
        SigBcgRatio.SetDrawOption("e1");
        SigBcgRatio.UseCurrentStyle();
        canvas->ModifiedUpdate();
        canvas->RedrawAxis();
        //drawing
        canvas->SaveAs(draw_full_path.c_str());
    }
    //resetting style and canvas
    //needed because really WEIRD things were happening (somehow BELLE2 style set itself up)
    //TODO fix that
    canvas->Clear();
}

void draw_and_save(TH1D* data, std::string folderWithDiagonal, std::string name, std::string title, std::string options){
    //failsave in case nullptr was passed
    if(data==nullptr){
        printf("WARNING!!! A (data) nullptr has been passed to draw_and_save function!\n");
    }
    //normal proceeding
    MyStyles styleLibrary;
    TStyle tempStyle = styleLibrary.Hist2DQuarterSize(true);
    tempStyle.cd();
    gROOT->ForceStyle();
    TCanvas* resultCanvas = new TCanvas("resultCanvas", "resultCanvas", 4000, 2400);
    data->Draw(options.c_str());
    resultCanvas->UseCurrentStyle();
    data->SetMinimum(0);
    data->SetMarkerStyle(kFullCircle);
    data->SetMarkerSize(2);
    data->SetMarkerColor(kBlue);
    data->SetTitle(title.c_str());
    resultCanvas->SaveAs((folderWithDiagonal+name+".pdf").c_str());
    resultCanvas->Clear();
    delete resultCanvas;
}

void draw_and_save_minus_background(TH1D* data, TH1D* bcg, std::string folderWithDiagonal, std::string name, std::string title, double bcg_region){
    //failsave in case nullptr was passed
    if(data==nullptr){
        printf("WARNING!!! A (data) nullptr has been passed to draw_and_save_minus_background function!\n");
    } else if(bcg==nullptr){
        printf("WARNING!!! A (bcg) nullptr has been passed to draw_and_save_minus_background function!\n");
    }
    //normal proceeding
    MyStyles styleLibrary;
    TStyle tempStyle = styleLibrary.Hist2DQuarterSize(true);
    tempStyle.cd();
    gROOT->ForceStyle();
    TCanvas* resultCanvas = new TCanvas("resultCanvas", "resultCanvas", 4000, 2400);
    int bcg_bin = data->FindBin(bcg_region);
    if(data->Integral(bcg_bin, -1)==0 or bcg->Integral(bcg_bin, -1)==0){
        data->SetTitle((std::string(data->GetTitle())+" no background  removed").c_str());
    } else{
        bcg->Scale(data->Integral(bcg_bin, -1)/bcg->Integral(bcg_bin, -1));
        data->Add(bcg, -1.);
    }
    data->Draw("e1");
    resultCanvas->UseCurrentStyle();
    data->SetMarkerStyle(kFullCircle);
    data->SetMarkerSize(2);
    data->SetMarkerColor(kBlue);
    data->SetTitle(title.c_str());

    resultCanvas->SaveAs((folderWithDiagonal+name+".pdf").c_str());
    resultCanvas->Clear();
    delete resultCanvas;
}

void draw_bulk(std::vector<TH1D*> data, std::string folderWithDiagonal, std::string name, std::string title, std::string options){
    MyStyles styleLibrary;
    TStyle tempStyle = styleLibrary.Hist2DNormalSize(true);
    tempStyle.cd();
    gROOT->ForceStyle();
    TCanvas* resultCanvas = new TCanvas("resultCanvas", "resultCanvas", 4000, 2400);
    THStack drawingStack(name.c_str(), title.c_str());
    drawingStack.SetMinimum(0);
    for(size_t i = 0; i<data.size(); i++){
        drawingStack.Add(data[i], options.c_str());
    }
    drawingStack.Draw("nostack");
    TList* hists = drawingStack.GetHists();
    for(size_t i = 0; i<drawingStack.GetNhists(); i++){
        //0 & 1 are black & white, and 10 is white too
        ((TH1D*)hists->At(i))->UseCurrentStyle();
        ((TH1D*)hists->At(i))->SetLineColor(i+2+(i+2>=10));
        ((TH1D*)hists->At(i))->SetMarkerColor(TColor::GetColorDark(i+2+(i+2>=10)));
        //writing the output, to just copy the fitted values already
        ((TH1D*)hists->At(i))->Print("all");

    }
    resultCanvas->BuildLegend();
    resultCanvas->Update();
    resultCanvas->SaveAs((folderWithDiagonal+name+".pdf").c_str());
    resultCanvas->Clear();
    delete resultCanvas;
}

TFitResult differential_crossection_fit(TPad* pad, TH1D* slice, TF1* fitting_function_signal, TF1* fitting_function_bcg, std::string draw_full_path){
    //failsave in case nullptr was passed
    if(slice==nullptr){
        printf("WARNING!!! A nullptr has been passed to differential_crossection_fit function!\n");
        TFitResult zeroResult;
        zeroResult.SetChi2AndNdf(0., 1.);
        return zeroResult;
    }
    //normal proceeding
    pad->Clear();
    gROOT->SetSelectedPad(pad);
    //setting up the functions
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
    slice->Draw("E1");
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
                sscanf(response.c_str(), "%d %lf", &coeff_number, &change);
            } else{
                relative = 0.;
                sscanf(response.c_str(), "abs %d %lf", &coeff_number, &change);
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
        for(size_t i = 0; i<totalparams; i++){
            printf("%s\t\t%lf\n", fitting_function_total_COPY->GetParName(i), fitting_function_total_COPY->GetParameter(i));
        }
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
    pad->Clear();
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