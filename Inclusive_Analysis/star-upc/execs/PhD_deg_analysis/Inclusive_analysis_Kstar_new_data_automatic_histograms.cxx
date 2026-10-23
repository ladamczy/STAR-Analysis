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
#include "TApplication.h"
#include "TFitResult.h"
#include "THStack.h"
#include "TLine.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TRootCanvas.h"
#include "TVirtualX.h"
#include "Math/MinimizerOptions.h"

#include "MyStyles.h"

std::string rountToNSignificantFigures(double input, int n = 2);
int GetFirstNonzeroBinNumber(TH1* input);
double SignalNormalizationEstimation(TH1D* signal, double x0, double gamma);
double BackgroundNormalizationEstimation(TH1D* signal, TH1D* background, std::vector<std::pair<double, double>> peak_positions_widths);
TFitResult fit_and_draw_and_save(TH1D* data, TF1* signal, TF1* background, std::string folderWithDiagonal, std::string name, std::string title, std::string options, TH1D* only_mixed_background = nullptr, std::string total_pdf_file_title = "");
void custom_draw_and_save(TH1D* data, double expected_value, std::string folderWithDiagonal, std::string name, std::string title, std::string options, std::string total_pdf_file_title);
void custom_draw_and_save_status(TH2D* data, std::string folderWithDiagonal, std::string name, std::string title, std::string options, std::string total_pdf_file_title);
std::string remove_bad_characters(std::string input);
void writing_to_one_pdf(TCanvas* input, std::string full_path);

int main(int argc, char* argv[]){
    //argv:
    //1 - .root file or .list list
    //2 - output folder
    //3 - (optional) batch mode (default no, accepts case insensitive yes/y/true/batch/-b)

    bool batch = false;
    //checking the number of arguments and setting up batch mode if needed
    switch(argc){
    case 4:
        //needs to set up the batch mode
        //brackets needed to initialize the batch_mode
    {
        std::string batch_mode(argv[3]);

        std::transform(batch_mode.begin(), batch_mode.end(),
            batch_mode.begin(),
            [](unsigned char c){
            return static_cast<char>(std::tolower(c));
        });

        if(batch_mode=="yes"||batch_mode=="y"||
            batch_mode=="true"||batch_mode=="batch"||
            batch_mode=="-b"){
            batch = true;
        }
    }
    case 3:
        //nothing happens (batch mode false), but still valid number of arguments
        break;
    default:
        printf("Invalid number of arguments (%d). Proper argument usage:\n", argc-1);
        printf("argv[1] - .root file or .list list\n");
        printf("argv[2] - output folder\n");
        printf("argv[3] - (optional) batch mode (default no, accepts case insensitive yes/y/true/batch/-b))\n");
        return 1;
        break;
    }

    // saving arguments before ROOT processes -b flag for batch mode
    std::string input_file(argv[1]);
    std::string output_folder(argv[2]);

    //creating argv for TApplication
    //arrays are here because with direct initialization warnings arise
    char app_name[] = "App";
    char batch_option[] = "-b";
    std::vector<char*> new_argv;
    new_argv.push_back(app_name);
    if(batch)
        new_argv.push_back(batch_option);
    int new_argc = static_cast<int>(new_argv.size());

    //something used so that the histograms would draw
    //https://stackoverflow.com/questions/30932725/painting-a-tcanvas-to-the-screen-in-a-compiled-root-cern-application
    TApplication theApp("App", &new_argc, new_argv.data());
    //just in case, setting batch mode explicitly
    if(batch)
        gROOT->SetBatch(kTRUE);
    printf("Batch mode set to %s\n", gROOT->IsBatch() ? "true" : "false");
    printf("Display is set to %s\n", gVirtualX->GetName());

    //input file
    TFile* input = TFile::Open(input_file.c_str());

    //directory manipulation
    std::string folderWithDiagonal = output_folder;
    if(folderWithDiagonal[folderWithDiagonal.size()-1]!='/'){
        folderWithDiagonal += "/";
    }

    //getting diffractive background and signal
    std::vector<std::string> pairTab = { "Kpi", "piK" };
    //data taken from PDG for neutral only, as a good central point of fitting
    double K892fitMaximum = 0.89556;
    double K892fitWidth = 0.0471;

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

    //###########################################################
    //       SETTING UP FUNCTIONS FOR SIGNAL AND BACKGROUND
    //###########################################################

    //SIGNAL FUNCTIONS
    //convolution of BW and gauss (default range [0. 1.] will be changed later anyway)
    //!!! GAMMA AND GAUSS STD ARE SWAPPED TO PRESERVE m, gamma, std ORDER !!!
    TF1 fit_func_sig_template("fit_func_sig_template", "[0]*TMath::Voigt(x-[1], [3], [2])");
    //vector of potential signals, filled with functions
    std::vector<TF1> signal_function_vector;
    //vector of PDG values of mass and width
    std::vector<std::pair<double, double>> PDG_data;
    //map of starting values for fits (might be different from PDG ones)
    std::vector<std::pair<double, double>> start_values;

    //setting up signal particles, in order
    //K*(700)
    signal_function_vector.push_back(*static_cast<TF1*>(fit_func_sig_template.Clone("Kstar700")));
    signal_function_vector.back().SetParNames("N_{K^{*}(700)}", "m_{K^{*}(700)}", "#Gamma_{K^{*}(700)}", "#sigma_{K^{*}(700)}");
    PDG_data.emplace_back(0.838, 0.463);
    start_values.emplace_back(0.838, 0.463);
    //custom limits
    signal_function_vector.back().SetParLimits(1, 0.7, 0.89);
    signal_function_vector.back().SetParLimits(2, 0.2, 0.6);

    //K*(892)
    signal_function_vector.push_back(*static_cast<TF1*>(fit_func_sig_template.Clone("Kstar892")));
    signal_function_vector.back().SetParNames("N_{K^{*}(892)}", "m_{K^{*}(892)}", "#Gamma_{K^{*}(892)}", "#sigma_{K^{*}(892)}");
    PDG_data.emplace_back(0.89556, 0.0471);
    start_values.emplace_back(0.89556, 0.0471);
    //custom limits
    signal_function_vector.back().SetParLimits(1, 0.85, 0.95);
    signal_function_vector.back().SetParLimits(2, 0.02, 0.06);

    //K*(1430)
    signal_function_vector.push_back(*static_cast<TF1*>(fit_func_sig_template.Clone("Kstar1430")));
    signal_function_vector.back().SetParNames("N_{K^{*}(1430)}", "m_{K^{*}(1430)}", "#Gamma_{K^{*}(1430)}", "#sigma_{K^{*}(1430)}");
    PDG_data.emplace_back(1.425, 0.27);
    start_values.emplace_back(1.425, 0.15);
    //custom limits
    signal_function_vector.back().SetParLimits(1, 1.35, 1.5);
    signal_function_vector.back().SetParLimits(2, 0.08, 0.3);

    //setting default parameter limits
    for(auto&& function:signal_function_vector){
        //made especially for current template, will have to be changed upon signal template change
        function.SetParLimits(0, 0., 1e5);
        function.SetParLimits(3, 0., 0.1);
    }


    //###########################################################
    //                 SETTING UP HISTOGRAMS
    //###########################################################

    //creating folders for saving pdfs
    for(auto&& pair:pairTab){
        for(auto&& category:allCategories){
            gSystem->mkdir((folderWithDiagonal+pair+"/"+category).c_str(), true);
            gSystem->mkdir((folderWithDiagonal+pair+"/"+category+"/whole_background").c_str(), true);
        }
    }

    //filling vectors of background and signal and Chi2 with and without the background
    //this scary nested map exists so I can do map[pair][category] = histogram
    //data
    std::map<std::string, std::map<std::string, TH2D*>> background_vector, signal_vector;
    //total chi2 and fit status
    std::map<std::string, std::map<std::string, TH1D*>> Chi2withbcg_vector;
    std::map<std::string, std::map<std::string, TH2D*>> Status_vector;
    //separate histograms for each signal
    std::map<std::string, std::map<std::string, std::vector<TH1D*>>> result_vector;
    std::map<std::string, std::map<std::string, std::vector<TH1D*>>> mass_vector;
    std::map<std::string, std::map<std::string, std::vector<TH1D*>>> width_vector;
    std::map<std::string, std::map<std::string, std::vector<TH1D*>>> resolution_vector;
    //background
    std::map<std::string, std::map<std::string, TH1D*>> background_normalisation_vector;
    std::map<std::string, std::map<std::string, TH2D*>> background_additional_parameters_vector;

    //creating histograms
    std::string tempSignalName = "M$Chi2";
    std::string tempBackgroundName = "M$Chi2BcgMixedEvent";
    std::string tempHistName;
    for(auto&& pair:pairTab){
        for(auto&& category:allCategories){
            //signal
            tempHistName = tempSignalName;
            tempHistName.replace(find(tempHistName.begin(), tempHistName.end(), '$')-tempHistName.begin(), 1, pair);
            signal_vector[pair][category] = (TH2D*)input->Get((tempHistName+category).c_str());
            //getting the common part of the title
            std::string histogram_title = signal_vector[pair][category]->GetTitle();
            std::string axis_title = signal_vector[pair][category]->GetYaxis()->GetTitle();
            //creating histograms
            Chi2withbcg_vector[pair][category] = new TH1D((tempHistName+"_"+category+"_Chi2").c_str(), (histogram_title+"-dependent  #chi^{2}/NDF;"+axis_title+";#chi^{2}/NDF").c_str(), signal_vector[pair][category]->GetNbinsY(), signal_vector[pair][category]->GetYaxis()->GetXmin(), signal_vector[pair][category]->GetYaxis()->GetXmax());
            Status_vector[pair][category] = new TH2D((tempHistName+"_"+category+"_Fit_status").c_str(), (histogram_title+" fit status;"+axis_title+";Fit error code").c_str(), signal_vector[pair][category]->GetNbinsY(), signal_vector[pair][category]->GetYaxis()->GetXmin(), signal_vector[pair][category]->GetYaxis()->GetXmax(), 1, 0, 1);
            background_normalisation_vector[pair][category] = new TH1D((tempHistName+"_"+category+"_Bcgnorm").c_str(), (histogram_title+"-dependent background normalization;"+axis_title+";Bcg normalization").c_str(), signal_vector[pair][category]->GetNbinsY(), signal_vector[pair][category]->GetYaxis()->GetXmin(), signal_vector[pair][category]->GetYaxis()->GetXmax());
            for(size_t i = 0; i<signal_function_vector.size(); i++){
                //name of particle (takes N_{XYZ} and removes first 3 and last characters, leaving XYZ)
                std::string particle_name = std::string(signal_function_vector[i].GetParName(0)).substr(3, std::string(signal_function_vector[i].GetParName(0)).length()-4);
                //title
                std::string title_part = signal_function_vector[i].GetName();
                //histograms
                result_vector[pair][category].push_back(new TH1D((tempHistName+"_"+category+"_Result_"+title_part).c_str(), (histogram_title+"-dependent yield of "+particle_name+";"+axis_title+";Number of pairs").c_str(), signal_vector[pair][category]->GetNbinsY(), signal_vector[pair][category]->GetYaxis()->GetXmin(), signal_vector[pair][category]->GetYaxis()->GetXmax()));
                mass_vector[pair][category].push_back(new TH1D((tempHistName+"_"+category+"_Mass_"+title_part).c_str(), (histogram_title+"-dependent mass of "+particle_name+";"+axis_title+";m_{0} [GeV/c^{2}]").c_str(), signal_vector[pair][category]->GetNbinsY(), signal_vector[pair][category]->GetYaxis()->GetXmin(), signal_vector[pair][category]->GetYaxis()->GetXmax()));
                width_vector[pair][category].push_back(new TH1D((tempHistName+"_"+category+"_Width_"+title_part).c_str(), (histogram_title+"-dependent B-W width of "+particle_name+";"+axis_title+";#Gamma_{0} [GeV/c^{2}]").c_str(), signal_vector[pair][category]->GetNbinsY(), signal_vector[pair][category]->GetYaxis()->GetXmin(), signal_vector[pair][category]->GetYaxis()->GetXmax()));
                resolution_vector[pair][category].push_back(new TH1D((tempHistName+"_"+category+"_Resolution_"+title_part).c_str(), (histogram_title+"-dependent Gaussian width of "+particle_name+";"+axis_title+";#sigma [GeV/c^{2}]").c_str(), signal_vector[pair][category]->GetNbinsY(), signal_vector[pair][category]->GetYaxis()->GetXmin(), signal_vector[pair][category]->GetYaxis()->GetXmax()));
            }
            //background
            tempHistName = tempBackgroundName;
            tempHistName.replace(find(tempHistName.begin(), tempHistName.end(), '$')-tempHistName.begin(), 1, pair);
            background_vector[pair][category] = ((TH2D*)input->Get((tempHistName+category).c_str()));
            //creating histograms for additional background parameters
            //additional "  _{}" after ": " is because there were problems with lower index (and a space immediately after)
            background_additional_parameters_vector[pair][category] = new TH2D((tempHistName+"_"+category+"_").c_str(), (histogram_title+"-dependent background extra parameter:   _{};;"+axis_title).c_str(), 1, 0, 1, signal_vector[pair][category]->GetNbinsY(), signal_vector[pair][category]->GetYaxis()->GetXmin(), signal_vector[pair][category]->GetYaxis()->GetXmax());
        }
    }

    //###########################################################
    //                FITTING
    //###########################################################

    //parameters
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(500000);
    ROOT::Math::MinimizerOptions::SetDefaultMaxIterations(50000);
    ROOT::Math::MinimizerOptions::SetDefaultTolerance(1.0e-5);
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);

    //fitting loop
    for(auto&& pair:pairTab){
        for(auto&& category:allCategories){
            //opening pdf files
            std::string folderToSave = folderWithDiagonal+pair+"/"+category+"/";
            writing_to_one_pdf(nullptr, folderToSave+pair+"_"+category+"_TOTAL.pdf");
            writing_to_one_pdf(nullptr, folderToSave+pair+"_"+category+"_TOTAL_whole_background.pdf");
            //copying total histogram for given pair and cathegory
            TH2D* sig_pointer = (TH2D*)signal_vector[pair][category]->Clone();
            TH2D* bcg_pointer = (TH2D*)background_vector[pair][category]->Clone();
            //for keeping title
            std::string baseOfTitle = std::string(sig_pointer->GetTitle())+" #in ";
            //fitting loop
            for(Int_t k = 0; k<sig_pointer->GetNbinsY(); k++){
                //slices of signal and background
                TH1D* sig_slice = sig_pointer->ProjectionX("_sig", k+1, k+1, "e1");
                sig_slice->GetYaxis()->SetTitle("Number of pairs");
                TH1D* bcg_slice = bcg_pointer->ProjectionX("_bcg", k+1, k+1, "e1");
                //normalising background to the signal, so that fitting can be done in [0,1) range
                bcg_slice->Scale(sig_slice->Integral()/bcg_slice->Integral());

                //setting fitting range
                double lower_range = 0.65;
                double upper_range = 1.6;


                //creating background function - a scaled slice of mixed event background
                auto mixed_event_background = [&](double* x, double* p){
                    return p[0]*bcg_slice->GetBinContent(bcg_slice->FindBin(x[0]));
                };
                TF1 fit_func_bcg("fit_func_bcg", mixed_event_background, lower_range, upper_range, 1, 1);
                fit_func_bcg.SetParNames("N_{bcg}");
                //setting initial fitting values for background
                fit_func_bcg.SetParameters(BackgroundNormalizationEstimation(sig_slice, bcg_slice, start_values));
                fit_func_bcg.SetParLimits(0, 0., 1.);


                //creating signal function - a sum of all the signals
                auto signal_sum_function = [&](double* x, double* p){
                    double result = 0;
                    int starting_parameter = 0;
                    double lower_edge = bcg_slice->GetBinLowEdge(bcg_slice->FindBin(x[0]));
                    double upper_edge = bcg_slice->GetBinLowEdge(bcg_slice->FindBin(x[0])+1);
                    for(auto&& sig_func_partial:signal_function_vector){
                        sig_func_partial.SetParameters(p+starting_parameter);
                        //relaxed eps from 1e-12 to 1e-6, as it is more than enough for our purposes
                        //and tight constraints cause errors because of Voigt function implementation I guess
                        result += sig_func_partial.Integral(lower_edge, upper_edge, 1e-6);
                        starting_parameter += sig_func_partial.GetNpar();
                    }
                    return result;
                };
                int total_number_of_signal_parameters = 0;
                for(auto&& sig_func_partial:signal_function_vector){
                    total_number_of_signal_parameters += sig_func_partial.GetNpar();
                }
                TF1 fit_func_sig("fit_func_sig", signal_sum_function, lower_range, upper_range, total_number_of_signal_parameters, 1);
                //setting initial values for the function
                int temp_par_number = 0;
                for(size_t n_func = 0; n_func<signal_function_vector.size(); n_func++){
                    for(size_t n_par = 0; n_par<signal_function_vector[n_func].GetNpar(); n_par++){
                        //names
                        fit_func_sig.SetParName(temp_par_number, signal_function_vector[n_func].GetParName(n_par));
                        //limits
                        double min_limit, max_limit;
                        signal_function_vector[n_func].GetParLimits(n_par, min_limit, max_limit);
                        fit_func_sig.SetParLimits(temp_par_number, min_limit, max_limit);

                        //parameters
                        switch(n_par){
                        case 0:
                            //normalisation
                            fit_func_sig.SetParameter(temp_par_number, SignalNormalizationEstimation(sig_slice, start_values[n_func].first, start_values[n_func].second));
                            break;
                        case 1:
                            //mass
                            fit_func_sig.SetParameter(temp_par_number, start_values[n_func].first);
                            break;
                        case 2:
                            //gamma
                            fit_func_sig.SetParameter(temp_par_number, start_values[n_func].second);
                            break;
                        case 3:
                            //sigma
                            fit_func_sig.SetParameter(temp_par_number, 0.006);
                            break;
                        default:
                            break;
                        }
                        //the end - increment number of total parameter in full function
                        temp_par_number++;
                    }
                }


                //fitting singular slice
                double par_value, par_error;
                std::string newTitle = baseOfTitle;
                std::string newName = pair+" "+category;
                std::string boundaries = "("+rountToNSignificantFigures(sig_pointer->GetYaxis()->GetBinLowEdge(k+1))+", ";
                boundaries += rountToNSignificantFigures(sig_pointer->GetYaxis()->GetBinUpEdge(k+1))+")";
                sig_slice->SetTitle(newTitle.c_str());
                //replacing title fragments
                //title is in format:
                //"pair_name category #in (lower_bound, upper_bound)"
                //so I exchange first space to " pair mass, ", so the whole title is:
                //"pair_name pair mass, category #in (lower_bound, upper_bound)"
                newTitle.replace(newTitle.find(" "), 1, " pair mass, ");
                newTitle += boundaries;
                newName += " "+boundaries;
                TFitResult tempResult = fit_and_draw_and_save(sig_slice, &fit_func_sig, &fit_func_bcg, folderToSave, newName, newTitle, "e1", bcg_slice, pair+"_"+category+"_TOTAL");
                //saving fit results for further analysys
                //Chi2
                Chi2withbcg_vector[pair][category]->SetBinContent(k+1, tempResult.Chi2()/tempResult.Ndf());
                //Status - careful, this is 2D histogram
                Status_vector[pair][category]->SetBinContent(k+1, 1, tempResult.Status());
                for(size_t i = 0; i<signal_function_vector.size(); i++){
                    //result
                    double bin_width = result_vector[pair][category][i]->GetBinWidth(k+1);
                    //TODO: check later if taking constant bin width does not mess things up
                    //TODO: change i*fit_func_sig_template.GetNpar() into proper values
                    par_value = fabs(tempResult.Parameter(i*fit_func_sig_template.GetNpar()+0)/sig_slice->GetXaxis()->GetBinWidth(1)*bin_width);
                    par_error = tempResult.ParError(i*fit_func_sig_template.GetNpar()+0)/sig_slice->GetXaxis()->GetBinWidth(1)*bin_width;
                    result_vector[pair][category][i]->SetBinContent(k+1, par_value);
                    result_vector[pair][category][i]->SetBinError(k+1, par_error);
                    //mass
                    mass_vector[pair][category][i]->SetBinContent(k+1, tempResult.Parameter(i*fit_func_sig_template.GetNpar()+1));
                    mass_vector[pair][category][i]->SetBinError(k+1, tempResult.Error(i*fit_func_sig_template.GetNpar()+1));
                    //width
                    width_vector[pair][category][i]->SetBinContent(k+1, tempResult.Parameter(i*fit_func_sig_template.GetNpar()+2));
                    width_vector[pair][category][i]->SetBinError(k+1, tempResult.Error(i*fit_func_sig_template.GetNpar()+2));
                    //resolution
                    resolution_vector[pair][category][i]->SetBinContent(k+1, tempResult.Parameter(i*fit_func_sig_template.GetNpar()+3));
                    resolution_vector[pair][category][i]->SetBinError(k+1, tempResult.Error(i*fit_func_sig_template.GetNpar()+3));
                }
                //background normalisation
                //TODO: change that when introducing multiple backgrounds
                background_normalisation_vector[pair][category]->SetBinContent(k+1, tempResult.Parameter(total_number_of_signal_parameters));
                background_normalisation_vector[pair][category]->SetBinError(k+1, tempResult.Error(total_number_of_signal_parameters));
                //background additional parameters
                int total_signal_parameters = fit_func_sig.GetNpar();
                //TODO: change that when introducing multiple backgrounds
                int total_mixed_event_background_parameters = 1;
                for(size_t par = total_signal_parameters+total_mixed_event_background_parameters; par<tempResult.NTotalParameters(); par++){
                    //this way we start fit_func_bcg parameters from after signal and after mixed background
                    int fit_func_bcg_additional_param_number = par-total_signal_parameters;
                    std::string parameter_name = fit_func_bcg.GetParName(fit_func_bcg_additional_param_number);
                    double bin_center = background_additional_parameters_vector[pair][category]->GetYaxis()->GetBinCenter(k+1);
                    background_additional_parameters_vector[pair][category]->Fill(parameter_name.c_str(), bin_center, tempResult.Parameter(par));
                    int parameter_bin_number = background_additional_parameters_vector[pair][category]->GetXaxis()->FindFixBin(parameter_name.c_str());
                    background_additional_parameters_vector[pair][category]->SetBinError(parameter_bin_number, k+1, tempResult.Error(par));
                }
            }
        }
    }

    //drawing results - chi2, number of detected decays and width of resonances and more
    for(auto&& pair:pairTab){
        for(auto&& category:allCategories){
            std::string folderToSave = folderWithDiagonal+pair+"/"+category+"/";
            //always present
            custom_draw_and_save(Chi2withbcg_vector[pair][category], 1., folderToSave, "", "", "hist min0", pair+"_"+category+"_TOTAL");
            custom_draw_and_save_status(Status_vector[pair][category], folderToSave, "", "", "col text min0", pair+"_"+category+"_TOTAL");
            for(size_t i = 0; i<signal_function_vector.size(); i++){
                custom_draw_and_save(result_vector[pair][category][i], 0., folderToSave, "", "", "e1 min0", pair+"_"+category+"_TOTAL");
                custom_draw_and_save(mass_vector[pair][category][i], PDG_data[i].first, folderToSave, "", "", "e1", pair+"_"+category+"_TOTAL");
                custom_draw_and_save(width_vector[pair][category][i], PDG_data[i].second, folderToSave, "", "", "e1 min0", pair+"_"+category+"_TOTAL");
                custom_draw_and_save(resolution_vector[pair][category][i], 0., folderToSave, "", "", "e1 min0", pair+"_"+category+"_TOTAL");
            }
            custom_draw_and_save(background_normalisation_vector[pair][category], 0., folderToSave, "", "", "e1 min0", pair+"_"+category+"_TOTAL");
            //custom parameters (including those from linear background)
            TH2D* hist_pointer = background_additional_parameters_vector[pair][category];
            hist_pointer->LabelsDeflate("X");
            for(size_t bin_number = 0; bin_number<hist_pointer->GetXaxis()->GetNbins(); bin_number++){
                //name (if there is nothing added, there is no additional parameters and no histogram should be made)
                std::string temp_name = hist_pointer->GetName();
                temp_name += remove_bad_characters(hist_pointer->GetXaxis()->GetBinLabel(bin_number+1));
                if(std::string(hist_pointer->GetXaxis()->GetBinLabel(bin_number+1)).length()==0){
                    continue;
                }
                //title
                std::string temp_title = hist_pointer->GetTitle();
                temp_title += hist_pointer->GetXaxis()->GetBinLabel(bin_number+1);
                temp_title += ";";
                temp_title += hist_pointer->GetYaxis()->GetTitle();
                temp_title += ";Value";
                //drawing
                custom_draw_and_save(hist_pointer->ProjectionY("_parameter", bin_number+1, bin_number+1), 0., folderToSave, temp_name, temp_title, "e1 text0", pair+"_"+category+"_TOTAL");
            }
            //end of drawing
            writing_to_one_pdf(nullptr, folderToSave+pair+"_"+category+"_TOTAL.pdf");
            writing_to_one_pdf(nullptr, folderToSave+pair+"_"+category+"_TOTAL_whole_background.pdf");
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

double SignalNormalizationEstimation(TH1D* signal, double x0, double gamma){
    double maximum = signal->GetBinContent(signal->FindBin(x0));
    double half_maximum = (signal->GetBinContent(signal->FindBin(x0+gamma/2))+signal->GetBinContent(signal->FindBin(x0-gamma/2)))/2;
    //by definition of FWHM, the height gamma/2 from peak will be half of the height of the peak
    //and I'm taking the mean from both sides, because the ground may be uneven
    double height = 2*(maximum-half_maximum);
    if(height<0)
        return 100.;
    //height of normalised to 1 B-W function at peak is 2/(pi*gamma)
    //for a very thin bin of dx width, the integral of function normalized to N therefore would be N*2*dx/(pi*gamma)
    //and this is exactly what height is - number of counts, compared to integral. Therefore:
    //height = N*2*dx/(pi*gamma)
    //N = height*pi*gamma/(2*dx)
    double N = height*TMath::Pi()*gamma/2/signal->GetBinWidth(signal->FindBin(x0));
    //now, safeguards:
    //if N<0, then there is small amount of events - so lets just return 100
    if(N<0){
        return 100.;
    }
    //if N>1e5, then probably gamma is big - so lets just return 1000 and see how it goes
    if(N>1e5){
        return 1000.;
    }
    //if everything is fine, we just return the N
    return N;
}

double BackgroundNormalizationEstimation(TH1D* signal, TH1D* background, std::vector<std::pair<double, double>> peak_positions_widths){
    //making the list of signal peaks, extended with the position of the lowest bin and a point 1 gamma after the last peak
    std::vector<int> points_of_interest;
    points_of_interest.push_back(signal->GetBinLowEdge(GetFirstNonzeroBinNumber(signal)));
    for(auto&& sig_pair:peak_positions_widths){
        //masses
        points_of_interest.push_back(sig_pair.first);
    }
    points_of_interest.push_back(peak_positions_widths.back().first+peak_positions_widths.back().second);

    //finding the position of biggest difference
    double difference = 0, middle_of_difference;
    for(size_t i = 1; i<points_of_interest.size(); i++){
        if(points_of_interest[i]-points_of_interest[i-1]>difference){
            difference = points_of_interest[i]-points_of_interest[i-1];
            middle_of_difference = (points_of_interest[i]+points_of_interest[i-1])/2;
        }
    }

    //looking up a bin in the middle and taking a sum of a few adjacent bins
    int adjacent_bins_to_sum = 2; //can be changed
    double signal_sum = 0;
    double background_sum = 0;
    int central_bin = signal->FindBin(middle_of_difference);
    for(int i = -adjacent_bins_to_sum; i<=adjacent_bins_to_sum; i++){
        signal_sum += signal->GetBinContent(i+central_bin);
        background_sum += background->GetBinContent(i+central_bin);
    }
    //returning ratio by which scale the background
    return signal_sum/background_sum;
}

TFitResult fit_and_draw_and_save(TH1D* data, TF1* signal, TF1* background, std::string folderWithDiagonal, std::string name, std::string title, std::string options, TH1D* only_mixed_background, std::string total_pdf_file_title){
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
    fitting_function_total->SetNpx(1000);

    //normal proceeding
    TCanvas* resultCanvas = MyStyles::DefaultCanvas("resultCanvas");
    TStyle tempStyle = MyStyles::Hist2DNormalSize(true);
    tempStyle.SetMarkerSize(0.5);
    tempStyle.cd();
    tempStyle.SetOptFit();
    resultCanvas->UseCurrentStyle();
    data->UseCurrentStyle();
    data->SetTitle(title.c_str());
    //other bits
    data->SetMarkerStyle(kFullCircle);
    data->SetMarkerColor(kBlue);
    //drawing data and function
    TFitResultPtr fitPointer = data->Fit(fitting_function_total, "B R S");
    data->Draw(options.c_str());
    //drawing legend
    double width = 0.3;
    double height = 0.1;
    double leftedge = 1.0-resultCanvas->GetRightMargin()-0.01-width;
    double loweredge = 1.0-resultCanvas->GetTopMargin()-0.01-height;
    TLegend legend_for_background_fitting(leftedge, loweredge, leftedge+width, loweredge+height);
    legend_for_background_fitting.UseCurrentStyle();
    legend_for_background_fitting.AddEntry(data, "Data");
    legend_for_background_fitting.AddEntry(fitting_function_total, "Mixed background + signal fit", "l");
    legend_for_background_fitting.Draw("SAME");
    //setting proper look of fitting parameters
    double stats_height = 0.35;
    gPad->Update();
    TPaveStats* stats = static_cast<TPaveStats*>(data->GetListOfFunctions()->FindObject("stats"));
    stats->SetX1NDC(leftedge);
    stats->SetX2NDC(leftedge+width);
    stats->SetY1NDC(loweredge-stats_height);
    stats->SetY2NDC(loweredge);
    stats->SetTextSize(0.025);
    //drawing before saving
    gPad->ModifiedUpdate();
    //saving
    resultCanvas->SaveAs((folderWithDiagonal+name+".pdf").c_str());
    writing_to_one_pdf(resultCanvas, folderWithDiagonal+total_pdf_file_title+".pdf");
    if(only_mixed_background!=nullptr){
        //drawing two lines at proper coordinates
        double lower_limit, upper_limit;
        fitting_function_total->GetRange(lower_limit, upper_limit);
        TLine line1(lower_limit, resultCanvas->GetUymin(), lower_limit, resultCanvas->GetUymax());
        TLine line2(upper_limit, resultCanvas->GetUymin(), upper_limit, resultCanvas->GetUymax());
        line1.SetLineColor(kRed);
        line1.SetLineStyle(kDashed);
        line1.SetLineWidth(2);
        line1.Draw("same");
        line2.SetLineColor(kRed);
        line2.SetLineStyle(kDashed);
        line2.SetLineWidth(2);
        line2.Draw("same");

        //redrawing fitting functions in full range and adding drawing components
        //preparations
        TF1* fitting_function_total = static_cast<TF1*>(data->GetListOfFunctions()->FindObject("fitting_function_total"));
        fitting_function_total->SetRange(resultCanvas->GetUxmin(), resultCanvas->GetUxmax());
        //signals
        signal->SetParameters(fitting_function_total->GetParameters());
        signal->SetLineColor(kRed);
        signal->SetLineStyle(kDashed);
        signal->SetNpx(1000);
        TF1* signal_copy = signal->DrawCopy("same");//fix to draw signal in the whooooole region
        signal_copy->SetRange(resultCanvas->GetUxmin(), resultCanvas->GetUxmax());
        //mixed background
        auto only_mixed_event_background = [&](double* x, double* p){
            return p[0]*only_mixed_background->GetBinContent(only_mixed_background->FindBin(x[0]));
        };
        TF1 fitting_function_mixed_background("fitting_function_mixed_background", only_mixed_event_background, resultCanvas->GetUxmin(), resultCanvas->GetUxmax(), 1, 1);
        fitting_function_mixed_background.SetParameter(0, fitting_function_total->GetParameter(signalparams));//we start from 0 here, so signalparams-1 is the last signal parameter
        fitting_function_mixed_background.SetLineColor(kBlue);
        fitting_function_mixed_background.SetLineStyle(kDashed);
        fitting_function_mixed_background.SetNpx(1000);
        fitting_function_mixed_background.Draw("same");
        //remaining, non-mixed background
        auto remaining_background = [&](double* x, double* p){
            return background->EvalPar(x, fitting_function_total->GetParameters()+signalparams)-fitting_function_mixed_background.Eval(x[0]);
        };
        TF1 fitting_function_remaining_background("fitting_function_remaining_background", remaining_background, resultCanvas->GetUxmin(), resultCanvas->GetUxmax(), 0, 1);
        fitting_function_remaining_background.SetLineColor(kGreen);
        fitting_function_remaining_background.SetLineStyle(kDashed);
        fitting_function_remaining_background.SetNpx(1000);
        fitting_function_remaining_background.Draw("same");

        //adding fitting region and other function markers to legend
        legend_for_background_fitting.AddEntry(&fitting_function_mixed_background, "Mixed events background", "l");
        legend_for_background_fitting.AddEntry(&fitting_function_remaining_background, "Other parts of the background", "l");
        legend_for_background_fitting.AddEntry(&line1, "Fitting region", "l");
        //moving legend and statbox
        loweredge -= 0.1;
        legend_for_background_fitting.SetY1NDC(loweredge-0.05);
        stats->SetY1NDC(loweredge-stats_height-0.05);
        stats->SetY2NDC(loweredge-0.05);

        //saving in subfolder
        gPad->ModifiedUpdate();
        gSystem->mkdir((folderWithDiagonal+"whole_background").c_str(), true);
        resultCanvas->SaveAs((folderWithDiagonal+"whole_background/"+name+".pdf").c_str());
        writing_to_one_pdf(resultCanvas, folderWithDiagonal+total_pdf_file_title+"_whole_background.pdf");
    }
    resultCanvas->Clear();
    delete resultCanvas;
    delete fitting_function_total;
    //returning result for further analysis
    return *(fitPointer.Get());
}

void custom_draw_and_save(TH1D* data, double expected_value, std::string folderWithDiagonal, std::string name, std::string title, std::string options, std::string total_pdf_file_title){
    //failsave in case nullptr was passed
    if(data==nullptr){
        printf("WARNING!!! A (data) nullptr has been passed to custom_draw_and_save function!\n");
    }
    //normal proceeding
    TStyle tempStyle = MyStyles::Hist2DNormalSize(true);
    //setting to align middle of text over middle of bin in "text0 e1" histograms
    tempStyle.SetTextAlign(kHAlignCenter+kVAlignBottom);
    tempStyle.cd();
    gROOT->ForceStyle();
    TCanvas* resultCanvas = MyStyles::DefaultCanvas("resultCanvas");
    resultCanvas->UseCurrentStyle();
    data->UseCurrentStyle();
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
    writing_to_one_pdf(resultCanvas, folderWithDiagonal+total_pdf_file_title+".pdf");
    writing_to_one_pdf(resultCanvas, folderWithDiagonal+total_pdf_file_title+"_whole_background.pdf");
    resultCanvas->Clear();
    delete resultCanvas;
}

void custom_draw_and_save_status(TH2D* data, std::string folderWithDiagonal, std::string name, std::string title, std::string options, std::string total_pdf_file_title){
    //failsave in case nullptr was passed
    if(data==nullptr){
        printf("WARNING!!! A (data) nullptr has been passed to custom_draw_and_save_status function!\n");
    }
    //normal proceeding
    TStyle tempStyle = MyStyles::Hist2DNormalSize(true);
    //setting to align middle of text over middle of bin in "text0 e1" histograms
    tempStyle.SetTextAlign(kHAlignCenter+kVAlignBottom);
    tempStyle.cd();
    gROOT->ForceStyle();
    TCanvas* resultCanvas = MyStyles::DefaultCanvas("resultCanvas");
    resultCanvas->UseCurrentStyle();
    data->UseCurrentStyle();
    //setting name and title
    if(name.length()!=0){
        data->SetName(name.c_str());
    }
    if(title.length()!=0){
        data->SetTitle(title.c_str());
    }
    //creating custom red-green palete and applying it to gStyle
    int palette[] = { kGreen, kRed };
    gStyle->SetPalette(2, palette);
    resultCanvas->SetFrameFillColor(TColor::GetColorPalette(0));
    //drawing data (with removed Y axis, but left the description so that the width will visually correspond to the bins)
    data->GetYaxis()->CenterTitle(true);
    data->GetYaxis()->SetLabelSize(0);
    data->GetYaxis()->SetTickLength(0);
    data->SetStats(false);
    data->SetMinimum(0.);
    data->SetMaximum(1.);
    data->SetMarkerSize(4);
    data->SetMarkerColor(kBlack);
    resultCanvas->SetGrid(1, 0);
    data->Draw(options.c_str());
    //saving
    resultCanvas->SaveAs((folderWithDiagonal+std::string(data->GetName())+".pdf").c_str());
    writing_to_one_pdf(resultCanvas, folderWithDiagonal+total_pdf_file_title+".pdf");
    writing_to_one_pdf(resultCanvas, folderWithDiagonal+total_pdf_file_title+"_whole_background.pdf");
    resultCanvas->Clear();
    //restoring the default palette
    gStyle->SetPalette();
    resultCanvas->SetGrid(0, 0);
    delete resultCanvas;
}

std::string remove_bad_characters(std::string input){
    //removes {,},#,^
    for(std::string substring : { "{", "}", "#", "^" }){
        while(input.find(substring)!=std::string::npos){
            input.replace(input.find(substring), substring.size(), "");
        }
    }
    //replaces * with "star"
    for(std::string substring : { "*" }){
        while(input.find(substring)!=std::string::npos){
            input.replace(input.find(substring), substring.size(), "star");
        }
    }
    return input;
}

void writing_to_one_pdf(TCanvas* input, std::string full_path){
    static std::map<std::string, int> page_numbers;
    static TCanvas* temp = MyStyles::DefaultCanvas("temp");
    //to hide the temporary window when not in batch mode
    if(!gROOT->IsBatch()){
        ((TRootCanvas*)temp->GetCanvasImp())->UnmapWindow();
    }
    //emplacing successful means it is the first call of a given file
    if(page_numbers.emplace(full_path, 0).second&&input==nullptr){
        temp->SaveAs((full_path+"[").c_str());
    } else if(input==nullptr){
        temp->SaveAs((full_path+"]").c_str());
        printf("Saved multipage pdf file %s with %d pages\n", full_path.c_str(), page_numbers[full_path]);
        page_numbers.erase(full_path);
    } else{
        input->SaveAs(full_path.c_str());
        page_numbers[full_path]++;
    }
}