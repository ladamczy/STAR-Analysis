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
void differential_crossection_fit(TH1D* slice, TF1* fitting_function_signal, TF1* fitting_function_bcg, double& param_value, double& param_error);

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
    TH1D* MKpiChi2bcg = (TH1D*)inputbcg->Get("MKpiChi2bcg");
    TH1D* MpiKChi2bcg = (TH1D*)inputbcg->Get("MpiKChi2bcg");
    TH1D* MppiChi2bcg = (TH1D*)inputbcg->Get("MppiChi2bcg");
    TH1D* MpipChi2bcg = (TH1D*)inputbcg->Get("MpipChi2bcg");
    TH1D* MKKChi2bcg = (TH1D*)inputbcg->Get("MKKChi2bcg");
    TH1D* MpipiChi2bcg = (TH1D*)inputbcg->Get("MpipiChi2bcg");
    TH1D* MppChi2bcg = (TH1D*)inputbcg->Get("MppChi2bcg");
    //getting diffractive background and signal
    // std::vector<std::string> pairTab = { "Kpi", "piK", "ppi", "pip", "KK", "pipi", "pp" };
    // std::vector<double> bcgRegion = { 1.0, 1.0, 1.1, 1.1, 1.1, 0.2, 2.4 };
    std::vector<std::string> pairTab = { "Kpi", "piK", "KK" };
    std::vector<double> bcgRegion = { 1.0, 1.0, 1.1 };
    std::vector<double> fitMaximum = { 0.892, 0.892, 1.02 };
    std::vector<double> fitWidth = { 0.05, 0.05, 0.004 };//50 MeV for K*(892), 4 MeV for phi(1020)
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
            //background
            tempHistName = tempBackgroundName;
            tempHistName.replace(find(tempHistName.begin(), tempHistName.end(), '$')-tempHistName.begin(), 1, pair);
            background_vector.push_back((TH2D*)inputbcg->Get((tempHistName+category).c_str()));
        }
    }
    //fitting functions
    TF1* fit_func_sig = new TF1("fit_func_sig", "breitwigner", 0.8, 1.0);
    TF1* fit_func_bcg = new TF1("fit_func_bcg", "pol1", 0.8, 1.0);
    //substracting one from another
    //fitting the difference
    //and filling the result
    for(size_t i = 0; i<pairTab.size(); i++){
        for(size_t j = 0; j<allCategories.size(); j++){
            TH2D* bcg_pointer = background_vector[i*allCategories.size()+j];
            TH2D* sig_pointer = signal_vector[i*allCategories.size()+j];
            for(Int_t k = 0; k<sig_pointer->GetNbinsY(); k++){
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
                fit_func_sig->SetParameters(2.3, fitMaximum[i], fitWidth[i]);
                fit_func_sig->SetParLimits(0, 0, 1e9);
                fit_func_sig->SetParLimits(1, fitMaximum[i]-fitWidth[i]*5, fitMaximum[i]+fitWidth[i]*5);
                fit_func_sig->SetParLimits(2, 0, fitWidth[i]*10);
                fit_func_sig->SetRange(fitMaximum[i]-fitWidth[i]*5, sig_slice->GetXaxis()->GetXmax());
                fit_func_bcg->SetParameters(0., 0.);
                fit_func_bcg->SetRange(fitMaximum[i]-fitWidth[i]*5, sig_slice->GetXaxis()->GetXmax());
                double par_value, par_error;
                std::string newTitle = std::string(sig_slice->GetTitle())+" ";
                newTitle += std::to_string(sig_pointer->GetYaxis()->GetBinLowEdge(k+1))+" - ";
                newTitle += std::to_string(sig_pointer->GetYaxis()->GetBinUpEdge(k+1));
                sig_slice->SetTitle(newTitle.c_str());
                differential_crossection_fit(sig_slice, fit_func_sig, fit_func_bcg, par_value, par_error);
                double bin_width = result_vector[i*allCategories.size()+j]->GetBinWidth(k+1);
                par_value *= bin_width;
                par_error *= bin_width;
                result_vector[i*allCategories.size()+j]->SetBinContent(k+1, par_value);
                result_vector[i*allCategories.size()+j]->SetBinError(k+1, par_error);
            }
        }
    }

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

void differential_crossection_fit(TH1D* slice, TF1* fitting_function_signal, TF1* fitting_function_bcg, double& param_value, double& param_error){
    //setting up the functions
    TCanvas* result = new TCanvas("result", "result", 1600, 800);
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
            printf("Chi2 = %f, ndof = %d\n", fitting_function_total->GetChisquare(), fitting_function_total->GetNDF());
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
    result->Close();
    gStyle->SetOptStat(1);
}