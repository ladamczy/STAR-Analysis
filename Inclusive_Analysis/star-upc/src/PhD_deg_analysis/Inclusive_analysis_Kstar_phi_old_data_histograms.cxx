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
#include "TApplication.h"

#include "MyStyles.h"

void draw_and_save(TH1D* data, std::string folderWithDiagonal, std::string name, std::string title, std::string options = "");
void draw_and_save_minus_background(TH1D* data, TH1D* bcg, std::string folderWithDiagonal, std::string name, std::string title, double bcg_region);
std::pair<double, double> differential_crossection_fit(TPad* pad, TH1D* slice, TF1* fitting_function_signal, TF1* fitting_function_bcg, double& param_value, double& param_error);

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
    bool fixWidthsInPlace = true;
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
    //and Chi2 with and without the background
    std::vector<TH2D*> background_vector, signal_vector;
    std::vector<TH1D*> result_vector, Chi2withbcg_vector, Chi2withoutbcg_vector;
    std::string tempSignalName = "M$Chi2";
    std::string tempBackgroundName = "M$Chi2bcgTOF";
    std::string tempHistName;
    for(auto&& pair:pairTab){
        for(auto&& category:allCategories){
            //signal
            tempHistName = tempSignalName;
            tempHistName.replace(find(tempHistName.begin(), tempHistName.end(), '$')-tempHistName.begin(), 1, pair);
            signal_vector.push_back((TH2D*)input->Get((tempHistName+category).c_str()));
            //result, but uses signal template
            result_vector.push_back(new TH1D((tempHistName+"Result"+category).c_str(), (pair+" "+category).c_str(), signal_vector.back()->GetNbinsY(), signal_vector.back()->GetYaxis()->GetXmin(), signal_vector.back()->GetYaxis()->GetXmax()));
            //Chi2 comparison
            Chi2withbcg_vector.push_back(new TH1D((tempHistName+"Chi2withbcg"+category).c_str(), (pair+" "+category+" bcg").c_str(), signal_vector.back()->GetNbinsY(), signal_vector.back()->GetYaxis()->GetXmin(), signal_vector.back()->GetYaxis()->GetXmax()));
            Chi2withoutbcg_vector.push_back(new TH1D((tempHistName+"Chi2withoutbcg"+category).c_str(), (pair+" "+category+" no bcg").c_str(), signal_vector.back()->GetNbinsY(), signal_vector.back()->GetYaxis()->GetXmin(), signal_vector.back()->GetYaxis()->GetXmax()));
            //background
            tempHistName = tempBackgroundName;
            tempHistName.replace(find(tempHistName.begin(), tempHistName.end(), '$')-tempHistName.begin(), 1, pair);
            background_vector.push_back((TH2D*)inputbcg->Get((tempHistName+category).c_str()));
        }
    }
    //fitting functions
    TF1* fit_func_sig = new TF1("fit_func_sig", "breitwigner", 0.8, 1.0);
    //backgrounds - one custom, with parameters p0, p1
    //being values at ends (p2, p3)
    auto custom_background = [&](double* x, double* p){
        double a = (p[1]-p[0])/(p[3]-p[2]);
        return a*(x[0]-p[2])+p[0];
    };
    TF1* fit_func_bcg = new TF1("fit_func_bcg", custom_background, 0.8, 1.0, 4);
    TF1* fit_func_empty_bcg = new TF1("fit_func_bcg", "0.", 0.8, 1.0);
    //substracting one from another
    //fitting the difference
    //and filling the result
    TCanvas* result = new TCanvas("result", "result", 1600, 800);
    for(size_t i = 0; i<pairTab.size(); i++){
        for(size_t j = 0; j<allCategories.size(); j++){
            TH2D* bcg_pointer = background_vector[i*allCategories.size()+j];
            TH2D* sig_pointer = signal_vector[i*allCategories.size()+j];
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
                fit_func_sig->SetParameters(6., fitMaximum[i], fitWidth[i]);
                fit_func_sig->SetParLimits(0, 0, 1e9);
                fit_func_sig->SetParLimits(1, fitMaximum[i]-fitWidth[i]*5, fitMaximum[i]+fitWidth[i]*5);
                if(fixWidthsInPlace){
                    fit_func_sig->FixParameter(2, fitWidth[i]);
                } else{
                    fit_func_sig->SetParLimits(2, 0, fitWidth[i]*10);
                }
                fit_func_sig->SetRange(fitMaximum[i]-std::max(fitWidth[i]*5, 0.1), sig_slice->GetXaxis()->GetXmax());
                //custom setting background function - p2 & p3 are the ends of the range
                fit_func_bcg->SetRange(fitMaximum[i]-std::max(fitWidth[i]*5, 0.1), sig_slice->GetXaxis()->GetXmax());
                fit_func_bcg->SetParameters(0., 0., fit_func_bcg->GetXmin(), fit_func_bcg->GetXmax());
                fit_func_bcg->FixParameter(2, fit_func_bcg->GetXmin());
                fit_func_bcg->FixParameter(3, fit_func_bcg->GetXmax());
                fit_func_bcg->SetParLimits(0, 0., 1e9);
                fit_func_bcg->SetParLimits(1, 0., 1e9);
                double par_value, par_error;
                std::string newTitle = baseOfTitle;
                newTitle += std::to_string(sig_pointer->GetYaxis()->GetBinLowEdge(k+1))+" - ";
                newTitle += std::to_string(sig_pointer->GetYaxis()->GetBinUpEdge(k+1));
                sig_slice->SetTitle(newTitle.c_str());
                auto tempPair = differential_crossection_fit(result, sig_slice, fit_func_sig, fit_func_bcg, par_value, par_error);
                double chi2, ndf;
                chi2 = tempPair.first;
                ndf = tempPair.second;
                double bin_width = result_vector[i*allCategories.size()+j]->GetBinWidth(k+1);
                par_value *= bin_width;
                par_error *= bin_width;
                result_vector[i*allCategories.size()+j]->SetBinContent(k+1, par_value);
                result_vector[i*allCategories.size()+j]->SetBinError(k+1, par_error);
                Chi2withbcg_vector[i*allCategories.size()+j]->SetBinContent(k+1, chi2/ndf);
            }
        }
    }
    //special fitting without the background
    for(size_t i = 0; i<pairTab.size(); i++){
        for(size_t j = 0; j<allCategories.size(); j++){
            TH2D* bcg_pointer = background_vector[i*allCategories.size()+j];
            TH2D* sig_pointer = signal_vector[i*allCategories.size()+j];
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
                fit_func_sig->SetParameters(6., fitMaximum[i], fitWidth[i]);
                fit_func_sig->SetParLimits(0, 0, 1e9);
                fit_func_sig->SetParLimits(1, fitMaximum[i]-fitWidth[i]*5, fitMaximum[i]+fitWidth[i]*5);
                if(fixWidthsInPlace){
                    fit_func_sig->FixParameter(2, fitWidth[i]);
                } else{
                    fit_func_sig->SetParLimits(2, 0, fitWidth[i]*10);
                }
                fit_func_sig->SetRange(fitMaximum[i]-std::max(fitWidth[i]*5, 0.1), sig_slice->GetXaxis()->GetXmax());
                fit_func_empty_bcg->SetRange(fitMaximum[i]-std::max(fitWidth[i]*5, 0.1), sig_slice->GetXaxis()->GetXmax());
                double par_value, par_error;
                std::string newTitle = baseOfTitle;
                newTitle += std::to_string(sig_pointer->GetYaxis()->GetBinLowEdge(k+1))+" - ";
                newTitle += std::to_string(sig_pointer->GetYaxis()->GetBinUpEdge(k+1));
                newTitle += " ZERO BACKGROUND";
                sig_slice->SetTitle(newTitle.c_str());
                auto tempPair = differential_crossection_fit(result, sig_slice, fit_func_sig, fit_func_empty_bcg, par_value, par_error);
                double chi2, ndf;
                chi2 = tempPair.first;
                ndf = tempPair.second;
                Chi2withoutbcg_vector[i*allCategories.size()+j]->SetBinContent(k+1, chi2/ndf);
            }
        }
    }
    result->Close();

    std::string folderWithDiagonal = std::string(static_cast<const char*>(argv[3]));
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
            draw_and_save(result_vector[i*allCategories.size()+j], "~/star-upc/", "M"+pairTab[i]+allCategories[j], pairTab[i]+" "+allCategories[j]+";"+allCategories[j]+";entries", "e");
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
            Chi2withbcg_vector[i*allCategories.size()+j]->SetMinimum(0.5);
            Chi2withbcg_vector[i*allCategories.size()+j]->SetMaximum(1.1*std::max(Chi2withbcg_vector[i*allCategories.size()+j]->GetMaximum(), Chi2withoutbcg_vector[i*allCategories.size()+j]->GetMaximum()));
            Chi2withbcg_vector[i*allCategories.size()+j]->Draw("hist");
            Chi2withoutbcg_vector[i*allCategories.size()+j]->Draw("hist same");
            resultCanvas->UseCurrentStyle();
            Chi2withbcg_vector[i*allCategories.size()+j]->SetTitle((pairTab[i]+" "+allCategories[j]+";"+allCategories[j]+";#chi^{2}/ndf").c_str());
            Chi2withbcg_vector[i*allCategories.size()+j]->SetLineColor(kBlue);
            Chi2withoutbcg_vector[i*allCategories.size()+j]->SetLineColor(kRed);
            resultCanvas->BuildLegend();
            resultCanvas->Update();
            resultCanvas->SaveAs(("~/star-upc/M"+pairTab[i]+allCategories[j]+"Chi2.pdf").c_str());
        }
    }

    //the end
    theApp.Run();
    return 0;
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

std::pair<double, double> differential_crossection_fit(TPad* pad, TH1D* slice, TF1* fitting_function_signal, TF1* fitting_function_bcg, double& param_value, double& param_error){
    //setting up the functions
    gROOT->SetSelectedPad(pad);
    gStyle->SetOptStat(0);
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
    }
    for(int i = 0; i<bcgparams; i++){
        fitting_function_bcg->GetParLimits(i, lowlim, highlim);
        if(lowlim==highlim&&lowlim*highlim!=0){
            fitting_function_total->FixParameter(i+signalparams, lowlim);
        } else{
            fitting_function_total->SetParLimits(i+signalparams, lowlim, highlim);
        }
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
    slice->Fit(fitting_function_total, "0BR");
    slice->DrawCopy("E", "NewDataNewTracks");
    fitting_function_signal->SetParameters(fitting_function_total->GetParameters());
    fitting_function_signal->SetParErrors(fitting_function_total->GetParErrors());
    fitting_function_bcg->SetParameters(fitting_function_total->GetParameters()+signalparams);
    TF1* fitting_function_total_COPY = fitting_function_total->DrawCopy("CSAME");
    TF1* fitting_function_bcg_COPY = fitting_function_bcg->DrawCopy("CSAME");
    TF1* fitting_function_signal_COPY = fitting_function_signal->DrawCopy("CSAME");
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
            slice->Fit(fitting_function_total, "0BR");
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
    //adding to the result histogram
    if(response.find("zero")!=std::string::npos){
        param_value = 0.;
        param_error = 0.;
    } else{
        param_value = fitting_function_signal_COPY->GetParameter(0)/slice->GetXaxis()->GetBinWidth(0);
        param_error = fitting_function_signal_COPY->GetParError(0)/slice->GetXaxis()->GetBinWidth(0);
    }
    printf("Accepted current parameters\n");
    //on last loop close the TCanvas
    gStyle->SetOptStat(1);
    pad->Clear();
    return std::pair<double, double>(fitting_function_total->GetChisquare(), fitting_function_total->GetNDF());
}