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
TFitResult fit_and_draw_and_save(TH1D* data, TF1* signal, TF1* background, std::string folderWithDiagonal, std::string name, std::string title, std::string options);
void custom_draw_and_save(TH1D* data, double expected_value, std::string folderWithDiagonal, std::string name, std::string title, std::string options);

int main(int argc, char* argv[]){
    //something used so that the histograms would draw
    //https://stackoverflow.com/questions/30932725/painting-a-tcanvas-to-the-screen-in-a-compiled-root-cern-application
    TApplication theApp("App", &argc, argv);
    argv = theApp.Argv();
    argc = theApp.Argc();

    //input file
    TFile* input = TFile::Open(static_cast<const char*>(argv[1]));

    // //getting signal histograms out
    // TH1D* MKpiChi2 = (TH1D*)input->Get("MKpiChi2");
    // TH1D* MpiKChi2 = (TH1D*)input->Get("MpiKChi2");

    // //getting background histograms out
    // TH1D* MKpiChi2bcg = (TH1D*)input->Get("MKpiChi2BcgMixedEvent");
    // TH1D* MpiKChi2bcg = (TH1D*)input->Get("MpiKChi2BcgMixedEvent");

    //getting diffractive background and signal
    std::vector<std::string> pairTab = { "Kpi", "piK" };
    double fitMaximum = 0.892;
    double fitWidth = 0.0514;//51.4 MeV for K*(892)
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
    std::vector<TH1D*> Chi2withbcg_vector,
        result_vector,
        mass_vector,
        width_vector,
        resolution_vector,
        background_normalisation_vector;
    std::string tempSignalName = "M$Chi2";
    std::string tempBackgroundName = "M$Chi2BcgMixedEvent";
    std::string tempHistName;
    for(auto&& pair:pairTab){
        for(auto&& category:allCategories){
            //signal
            tempHistName = tempSignalName;
            tempHistName.replace(find(tempHistName.begin(), tempHistName.end(), '$')-tempHistName.begin(), 1, pair);
            signal_vector.push_back((TH2D*)input->Get((tempHistName+category).c_str()));
            //results, width, chi2
            Chi2withbcg_vector.push_back(new TH1D((tempHistName+"_Chi2_"+category).c_str(), (pair+" "+category+" Chi2/NDF;"+category+";#chi^{2}/NDF").c_str(), signal_vector.back()->GetNbinsY(), signal_vector.back()->GetYaxis()->GetXmin(), signal_vector.back()->GetYaxis()->GetXmax()));
            result_vector.push_back(new TH1D((tempHistName+"_Result_"+category).c_str(), (pair+" "+category+" results;"+category+";Number of pairs").c_str(), signal_vector.back()->GetNbinsY(), signal_vector.back()->GetYaxis()->GetXmin(), signal_vector.back()->GetYaxis()->GetXmax()));
            mass_vector.push_back(new TH1D((tempHistName+"_Mass_"+category).c_str(), (pair+" "+category+" mass;"+category+";m_{0} [GeV/c^{2}]").c_str(), signal_vector.back()->GetNbinsY(), signal_vector.back()->GetYaxis()->GetXmin(), signal_vector.back()->GetYaxis()->GetXmax()));
            width_vector.push_back(new TH1D((tempHistName+"_Width_"+category).c_str(), (pair+" "+category+" width;"+category+";Width [GeV/c^{2}]").c_str(), signal_vector.back()->GetNbinsY(), signal_vector.back()->GetYaxis()->GetXmin(), signal_vector.back()->GetYaxis()->GetXmax()));
            resolution_vector.push_back(new TH1D((tempHistName+"_Resolution_"+category).c_str(), (pair+" "+category+" resolution;"+category+";Resolution [GeV/c^{2}]").c_str(), signal_vector.back()->GetNbinsY(), signal_vector.back()->GetYaxis()->GetXmin(), signal_vector.back()->GetYaxis()->GetXmax()));
            background_normalisation_vector.push_back(new TH1D((tempHistName+"_Bcgnorm_"+category).c_str(), (pair+" "+category+" results;"+category+";Bcg norm").c_str(), signal_vector.back()->GetNbinsY(), signal_vector.back()->GetYaxis()->GetXmin(), signal_vector.back()->GetYaxis()->GetXmax()));
            //background
            tempHistName = tempBackgroundName;
            tempHistName.replace(find(tempHistName.begin(), tempHistName.end(), '$')-tempHistName.begin(), 1, pair);
            background_vector.push_back((TH2D*)input->Get((tempHistName+category).c_str()));
        }
    }

    //###########################################################
    //                FITTING
    //###########################################################
    TF1Convolution conv_sig("breitwigner", "gausn", 0.6, 2.4); //extra range for all your convolution needs
    conv_sig.SetNofPointsFFT(10000);
    TF1 fit_func_sig("fit_func_sig", conv_sig, 0.99, 1.05, 6);

    fit_func_sig.SetParameter(0, 1.);           //N_BW
    fit_func_sig.SetParameter(1, fitMaximum);   //m0
    fit_func_sig.SetParameter(2, fitWidth);     //gamma0
    fit_func_sig.FixParameter(3, 1.);           //N_Gauss  (fixed to 1)
    fit_func_sig.FixParameter(4, 0.);           //mu_Gauss (fixed to 0)
    fit_func_sig.SetParameter(5, 0.0013);       //sigma_Gauss
    fit_func_sig.SetParNames("N_{BW}", "m_{0}", "#Gamma_{0,BW}", "N_{Gauss}", "m_{0,Gauss}", "#sigma_{Gauss}");

    for(size_t i = 0; i<pairTab.size(); i++){
        for(size_t j = 0; j<allCategories.size(); j++){
            //copying total histogram for given pair and cathegory
            TH2D* sig_pointer = (TH2D*)signal_vector[i*allCategories.size()+j]->Clone();
            TH2D* bcg_pointer = (TH2D*)background_vector[i*allCategories.size()+j]->Clone();
            //for keeping title
            std::string baseOfTitle = std::string(sig_pointer->GetTitle())+" ";
            for(Int_t k = 0; k<sig_pointer->GetNbinsY(); k++){
                TH1D* sig_slice = sig_pointer->ProjectionX("_sig", k+1, k+1, "e1");
                sig_slice->GetYaxis()->SetTitle("Number of pairs");
                TH1D* bcg_slice = bcg_pointer->ProjectionX("_bcg", k+1, k+1, "e1");
                //creating background function - a scaled slice of mixed event background
                auto mixed_event_background = [&](double* x, double* p){
                    return p[0]*bcg_slice->GetBinContent(bcg_slice->FindBin(x[0]));
                };
                TF1 fit_func_bcg("fit_func_bcg", mixed_event_background, 0.99, 1.05, 1, 1);

                //fitting and filling result
                //setting lower range for m0-3*gamma
                //if lower bound is lower than m0-4*gamma
                //and 2*width otherwise
                double lower_range = fitMaximum-3*fitWidth;
                if(sig_slice->GetBinLowEdge(GetFirstNonzeroBinNumber(sig_slice))>fitMaximum-4*fitWidth){
                    lower_range = fitMaximum-2*fitWidth;
                }
                fit_func_sig.SetRange(lower_range, 1.5);
                fit_func_bcg.SetRange(lower_range, 1.5);
                fit_func_bcg.SetParNames("N_{bcg}");

                //setting initial fitting values
                double initial_bcg_scale = sig_slice->GetBinContent(sig_slice->FindBin(fitMaximum+fitWidth))/bcg_slice->GetBinContent(bcg_slice->FindBin(fitMaximum+fitWidth));
                fit_func_bcg.SetParameters(initial_bcg_scale);
                double initial_sig_height = (sig_slice->GetBinContent(sig_slice->FindBin(fitMaximum))-sig_slice->GetBinContent(sig_slice->FindBin(fitMaximum+fitWidth)))*TMath::PiOver2()*fitWidth;
                fit_func_sig.SetParameters(initial_sig_height, fitMaximum, fitWidth);

                //fitting singular slice
                double par_value, par_error;
                std::string newTitle = baseOfTitle;
                newTitle += "("+rountToNSignificantFigures(sig_pointer->GetYaxis()->GetBinLowEdge(k+1))+", ";
                newTitle += rountToNSignificantFigures(sig_pointer->GetYaxis()->GetBinUpEdge(k+1))+")";
                sig_slice->SetTitle(newTitle.c_str());
                TFitResult tempResult = fit_and_draw_and_save(sig_slice, &fit_func_sig, &fit_func_bcg, folderWithDiagonal, newTitle, newTitle, "e1");

                //saving fit results for further analysys
                //Chi2
                Chi2withbcg_vector[i*allCategories.size()+j]->SetBinContent(k+1, tempResult.Chi2()/tempResult.Ndf());
                //result
                double bin_width = result_vector[i*allCategories.size()+j]->GetBinWidth(k+1);
                par_value = fabs(tempResult.Parameter(0)/sig_slice->GetXaxis()->GetBinWidth(1)*bin_width);
                par_error = tempResult.ParError(0)/sig_slice->GetXaxis()->GetBinWidth(1)*bin_width;
                result_vector[i*allCategories.size()+j]->SetBinContent(k+1, par_value);
                result_vector[i*allCategories.size()+j]->SetBinError(k+1, par_error);
                //mass
                mass_vector[i*allCategories.size()+j]->SetBinContent(k+1, fabs(tempResult.Parameter(1)));
                mass_vector[i*allCategories.size()+j]->SetBinError(k+1, tempResult.Error(1));
                //width
                width_vector[i*allCategories.size()+j]->SetBinContent(k+1, fabs(tempResult.Parameter(2)));
                width_vector[i*allCategories.size()+j]->SetBinError(k+1, tempResult.Error(2));
                //resolution
                resolution_vector[i*allCategories.size()+j]->SetBinContent(k+1, fabs(tempResult.Parameter(5)));
                resolution_vector[i*allCategories.size()+j]->SetBinError(k+1, tempResult.Error(5));
                //background normalisation
                background_normalisation_vector[i*allCategories.size()+j]->SetBinContent(k+1, fabs(tempResult.Parameter(6)));
                background_normalisation_vector[i*allCategories.size()+j]->SetBinError(k+1, tempResult.Error(6));
            }
        }
    }

    //drawing results - chi2, number of detected decays and width of resonances
    for(size_t i = 0; i<pairTab.size(); i++){
        for(size_t j = 0; j<allCategories.size(); j++){
            custom_draw_and_save(Chi2withbcg_vector[i*allCategories.size()+j], 1., folderWithDiagonal, "", "", "hist min0");
            custom_draw_and_save(result_vector[i*allCategories.size()+j], 0., folderWithDiagonal, "", "", "e1 min0");
            custom_draw_and_save(mass_vector[i*allCategories.size()+j], fitMaximum, folderWithDiagonal, "", "", "e1");
            custom_draw_and_save(width_vector[i*allCategories.size()+j], fitWidth, folderWithDiagonal, "", "", "e1 min0");
            custom_draw_and_save(resolution_vector[i*allCategories.size()+j], 0., folderWithDiagonal, "", "", "e1 min0");
            custom_draw_and_save(background_normalisation_vector[i*allCategories.size()+j], 0., folderWithDiagonal, "", "", "e1 min0");
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

TFitResult fit_and_draw_and_save(TH1D* data, TF1* signal, TF1* background, std::string folderWithDiagonal, std::string name, std::string title, std::string options){
    //failsave in case nullptr was passed
    if(data==nullptr){
        printf("WARNING!!! A (data) nullptr has been passed to fit_and_draw_and_save function!\n");
    }

    //creating total fitting function
    int signalparams, bcgparams, totalparams;
    double rangemin, rangemax;
    signalparams = signal->GetNpar();
    bcgparams = background->GetNpar();
    totalparams = signalparams+bcgparams;
    signal->GetRange(rangemin, rangemax);
    auto fitting_function_sum = [&](double* x, double* par){
        return signal->EvalPar(x, par)+background->EvalPar(x, par+signalparams);
    };
    TF1* fitting_function_total = new TF1("fitting_function_total", fitting_function_sum, rangemin, rangemax, totalparams);
    double lowlim, highlim;
    for(int i = 0; i<signalparams; i++){
        signal->GetParLimits(i, lowlim, highlim);
        if(lowlim==highlim&&lowlim*highlim!=0){
            fitting_function_total->FixParameter(i, lowlim);
        } else{
            fitting_function_total->SetParLimits(i, lowlim, highlim);
        }
        //if not default name (beginning with p), copy it
        if(std::string(signal->GetParName(i))[0]!='p')
            fitting_function_total->SetParName(i, signal->GetParName(i));
    }
    for(int i = 0; i<bcgparams; i++){
        background->GetParLimits(i, lowlim, highlim);
        if(lowlim==highlim&&lowlim*highlim!=0){
            fitting_function_total->FixParameter(i+signalparams, lowlim);
        } else{
            fitting_function_total->SetParLimits(i+signalparams, lowlim, highlim);
        }
        //if not default name (beginning with p), copy it
        if(std::string(background->GetParName(i))[0]!='p')
            fitting_function_total->SetParName(i+signalparams, background->GetParName(i));
    }
    //setting parameters for total function
    Double_t params[totalparams];
    signal->GetParameters(params);
    background->GetParameters(params+signalparams);
    fitting_function_total->SetParameters(params);
    //setting function look
    fitting_function_total->SetLineColor(kRed);
    fitting_function_total->SetLineWidth(2);

    //normal proceeding
    TCanvas* resultCanvas = MyStyles::DefaultCanvas("resultCanvas");
    MyStyles styleLibrary;
    TStyle tempStyle = styleLibrary.Hist2DNormalSize(true);
    tempStyle.SetMarkerSize(0.5);
    tempStyle.cd();
    tempStyle.SetOptFit();
    resultCanvas->UseCurrentStyle();
    //drawing data and function
    data->SetTitle(title.c_str());
    data->SetMarkerStyle(kFullCircle);
    data->SetMarkerColor(kBlue);
    TFitResultPtr fitPointer = data->Fit(fitting_function_total, "BRS");
    data->Draw(options.c_str());
    //drawing legend
    double width = 0.3;
    double height = 0.1;
    double leftedge = 1.0-resultCanvas->GetRightMargin()-0.01-width;
    double loweredge = 1.0-resultCanvas->GetTopMargin()-0.01-height;
    TLegend legend_for_background_fitting(leftedge, loweredge, leftedge+width, loweredge+height);
    legend_for_background_fitting.SetTextSize(0.03);
    legend_for_background_fitting.AddEntry(data, "Data");
    legend_for_background_fitting.AddEntry(fitting_function_total, "Mixed background + signal fit", "l");
    legend_for_background_fitting.SetBorderSize(0);
    legend_for_background_fitting.DrawClone("SAME");
    //setting proper look of fitting parameters
    double stats_height = 0.25;
    gPad->Update();
    TPaveStats* stats = static_cast<TPaveStats*>(data->GetListOfFunctions()->FindObject("stats"));
    stats->SetX1NDC(leftedge);
    stats->SetX2NDC(leftedge+width);
    stats->SetY1NDC(loweredge-stats_height);
    stats->SetY2NDC(loweredge);
    stats->SetTextSize(0.03);
    stats->SetBorderSize(0);
    //drawing before saving
    gPad->ModifiedUpdate();
    //saving
    resultCanvas->SaveAs((folderWithDiagonal+name+".pdf").c_str());
    resultCanvas->Clear();
    delete resultCanvas;
    delete fitting_function_total;
    //returning result for further analysis
    return *(fitPointer.Get());
}

void custom_draw_and_save(TH1D* data, double expected_value, std::string folderWithDiagonal, std::string name, std::string title, std::string options){
    //failsave in case nullptr was passed
    if(data==nullptr){
        printf("WARNING!!! A (data) nullptr has been passed to custom_draw_and_save function!\n");
    }
    //normal proceeding
    MyStyles styleLibrary;
    TStyle tempStyle = styleLibrary.Hist2DNormalSize(true);
    tempStyle.cd();
    gROOT->ForceStyle();
    TCanvas* resultCanvas = MyStyles::DefaultCanvas("resultCanvas");
    resultCanvas->UseCurrentStyle();
    //setting name and title
    if(name.length()!=0){
        data->SetName(name.c_str());
    }
    if(title.length()!=0){
        data->SetTitle(title.c_str());
    }
    //drawing data
    data->SetMarkerStyle(kFullCircle);
    data->SetMarkerColor(kBlue);
    data->Draw(options.c_str());
    //drawing default lines and legend
    resultCanvas->Update(); //needed to register the change in Ux values
    TLine line1(resultCanvas->GetUxmin(), expected_value, resultCanvas->GetUxmax(), expected_value);
    line1.SetLineColor(kRed);
    line1.SetLineStyle(kDashed);
    line1.SetLineWidth(2);
    if(expected_value>0){
        line1.Draw("same");
    }
    //saving
    resultCanvas->ModifiedUpdate();
    resultCanvas->RedrawAxis();
    resultCanvas->SaveAs((folderWithDiagonal+std::string(data->GetName())+".pdf").c_str());
    resultCanvas->Clear();
    delete resultCanvas;
}