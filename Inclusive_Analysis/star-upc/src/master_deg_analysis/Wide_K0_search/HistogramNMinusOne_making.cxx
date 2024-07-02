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
#include "TLine.h"
#include "TEllipse.h"
#include "TLatex.h"

#include <iomanip>
#include <sstream>
#include <fstream>
#include <numeric>

using namespace std;

int main(){
    TFile *histograms = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/AnaOutput_Histograms_NMinusOne_new_data_new_tracks.root");

    gStyle->SetFrameLineWidth(1);
    // gStyle->SetOptFit(111);
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);

    //lines for fiducial
    TLine *fidVertOne, *fidVertTwo, *fidHorOne, *fidHorTwo, *fidHorThree, *fidHorFour;
    TEllipse *fidEllOne, *fidEllTwo;
    fidVertOne = new TLine(-0.27, -0.8, -0.27, -0.4);
    fidVertOne->SetLineWidth(2);
    fidVertOne->SetLineColor(kRed);
    fidVertTwo = new TLine(-0.27, 0.4, -0.27, 0.8);
    fidVertTwo->SetLineWidth(2);
    fidVertTwo->SetLineColor(kRed);
    fidHorOne = new TLine(-0.27, -0.8, 0.1824, -0.8);
    fidHorOne->SetLineWidth(2);
    fidHorOne->SetLineColor(kRed);
    fidHorTwo = new TLine(-0.27, -0.4, 0.4451, -0.4);
    fidHorTwo->SetLineWidth(2);
    fidHorTwo->SetLineColor(kRed);
    fidHorThree = new TLine(-0.27, 0.4, 0.4451, 0.4);
    fidHorThree->SetLineWidth(2);
    fidHorThree->SetLineColor(kRed);
    fidHorFour = new TLine(-0.27, 0.8, 0.1824, 0.8);
    fidHorFour->SetLineWidth(2);
    fidHorFour->SetLineColor(kRed);
    fidEllOne = new TEllipse(-0.6, 0, sqrt(1.25), sqrt(1.25), -20.9725, -45.6649);
    fidEllOne->SetLineWidth(2);
    fidEllOne->SetLineColor(kRed);
    fidEllOne->SetFillColorAlpha(0, 0);
    fidEllOne->SetNoEdges();
    fidEllTwo = new TEllipse(-0.6, 0, sqrt(1.25), sqrt(1.25), 20.9725, 45.6649);
    fidEllTwo->SetLineWidth(2);
    fidEllTwo->SetLineColor(kRed);
    fidEllTwo->SetFillColorAlpha(0, 0);
    fidEllTwo->SetNoEdges();
    //fiducial
    TH2D *RP_FIDUCIAL_east = (TH2D *)histograms->Get("RP_FIDUCIAL_east");
    TH2D *RP_FIDUCIAL_west = (TH2D *)histograms->Get("RP_FIDUCIAL_west");
    RP_FIDUCIAL_east->SetTitle(";p_{x} [GeV/c];p_{y} [GeV/c]");
    RP_FIDUCIAL_west->SetTitle(";p_{x} [GeV/c];p_{y} [GeV/c]");
    TCanvas *c1 = new TCanvas("c1", "c1", 1000, 1500);
    c1->SetMargin(0.15, 0.15, 0.1, 0.05);
    RP_FIDUCIAL_east->Draw("colz");
    c1->SetRealAspectRatio();
    fidVertOne->Draw("same");
    fidVertTwo->Draw("same");
    fidHorOne->Draw("same");
    fidHorTwo->Draw("same");
    fidHorThree->Draw("same");
    fidHorFour->Draw("same");
    fidEllOne->Draw("same");
    fidEllTwo->Draw("same");
    c1->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/master_thesis/RP_FIDUCIAL_east.pdf").c_str());
    delete c1;
    c1 = new TCanvas("c1", "c1", 1000, 1500);
    c1->SetMargin(0.15, 0.15, 0.1, 0.05);
    RP_FIDUCIAL_west->Draw("colz");
    c1->SetRealAspectRatio();
    fidVertOne->Draw("same");
    fidVertTwo->Draw("same");
    fidHorOne->Draw("same");
    fidHorTwo->Draw("same");
    fidHorThree->Draw("same");
    fidHorFour->Draw("same");
    fidEllOne->Draw("same");
    fidEllTwo->Draw("same");
    c1->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/master_thesis/RP_FIDUCIAL_west.pdf").c_str());
    delete c1;



    TH1D *RP_XI_east_CPT2 = (TH1D *)histograms->Get("RP_XI_east_CPT2");
    TH1D *RP_XI_west_CPT2 = (TH1D *)histograms->Get("RP_XI_west_CPT2");
    RP_XI_east_CPT2->SetTitle(";#xi;protons");
    RP_XI_west_CPT2->SetTitle(";#xi;protons");
    c1 = new TCanvas("c1", "c1", 1600, 1000);
    c1->SetRightMargin(0.05);
    RP_XI_east_CPT2->Draw();
    RP_XI_west_CPT2->SetLineColor(kBlack);
    RP_XI_west_CPT2->Draw("same");
    TLegend *legendXi = new TLegend(0.6, 0.7, 0.79, 0.89);
    legendXi->AddEntry("RP_XI_east_CPT2", "East side", "l");
    legendXi->AddEntry("RP_XI_west_CPT2", "West side", "l");
    legendXi->SetBorderSize(0);
    legendXi->Draw("SAME");
    c1->Update();
    TLine *lineXi1 = new TLine(0.005, 0, 0.005, c1->GetUymax());
    lineXi1->SetLineStyle(kDashed);
    lineXi1->SetLineColor(kRed);
    lineXi1->SetLineWidth(1);
    lineXi1->Draw("same");
    TLine *lineXi2 = new TLine(0.2, 0, 0.2, c1->GetUymax());
    lineXi2->SetLineStyle(kDashed);
    lineXi2->SetLineColor(kRed);
    lineXi2->SetLineWidth(1);
    lineXi2->Draw("same");
    c1->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/master_thesis/RP_XI_CPT2.pdf").c_str());
    delete c1;

    TH1D *RP_XI_east_CPT2noBBCL = (TH1D *)histograms->Get("RP_XI_east_CPT2noBBCL");
    TH1D *RP_XI_west_CPT2noBBCL = (TH1D *)histograms->Get("RP_XI_west_CPT2noBBCL");
    RP_XI_east_CPT2noBBCL->SetTitle(";#xi;protons");
    RP_XI_west_CPT2noBBCL->SetTitle(";#xi;protons");
    c1 = new TCanvas("c1", "c1", 1600, 1000);
    c1->SetRightMargin(0.05);
    RP_XI_east_CPT2noBBCL->Draw();
    RP_XI_west_CPT2noBBCL->SetLineColor(kBlack);
    RP_XI_west_CPT2noBBCL->Draw("same");
    legendXi = new TLegend(0.6, 0.7, 0.79, 0.89);
    legendXi->AddEntry("RP_XI_east_CPT2noBBCL", "East side", "l");
    legendXi->AddEntry("RP_XI_west_CPT2noBBCL", "West side", "l");
    legendXi->SetBorderSize(0);
    legendXi->Draw("SAME");
    c1->Update();
    lineXi1 = new TLine(0.005, 0, 0.005, c1->GetUymax());
    lineXi1->SetLineStyle(kDashed);
    lineXi1->SetLineColor(kRed);
    lineXi1->SetLineWidth(1);
    lineXi1->Draw("same");
    lineXi2 = new TLine(0.08, 0, 0.08, c1->GetUymax());
    lineXi2->SetLineStyle(kDashed);
    lineXi2->SetLineColor(kRed);
    lineXi2->SetLineWidth(1);
    lineXi2->Draw("same");
    c1->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/master_thesis/RP_XI_CPT2noBBCL.pdf").c_str());
    delete c1;



    TH2D *RP_ELASTIC_xi = (TH2D *)histograms->Get("RP_ELASTIC_xi");
    RP_ELASTIC_xi->SetTitle(";#xi_{E};#xi_{W}");
    c1 = new TCanvas("c1", "c1", 1800, 1400);
    c1->SetMargin(0.15, 0.15, 0.1, 0.05);
    double x_0 = -4.48170e-04;
    double sigma_x = 1.79095e-03;
    double y_0 = -8.04898e-04;
    double sigma_y = 2.12035e-03;
    TEllipse *elasticXi = new TEllipse(x_0, y_0, 3*sigma_x, 3*sigma_y);
    elasticXi->SetLineColor(kRed);
    RP_ELASTIC_xi->Draw("colz");
    c1->Update();
    elasticXi->Draw("same");
    c1->SetRealAspectRatio();
    c1->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/master_thesis/RP_ELASTIC_xi.pdf").c_str());
    delete c1;

    TH2D *RP_ELASTIC_theta = (TH2D *)histograms->Get("RP_ELASTIC_theta");
    RP_ELASTIC_theta->SetTitle(";#theta_{XZ};#theta_{YZ}");
    c1 = new TCanvas("c1", "c1", 2200, 1400);
    c1->SetMargin(0.15, 0.15, 0.1, 0.05);
    RP_ELASTIC_theta->Draw("colz");
    c1->SetRealAspectRatio();
    c1->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/master_thesis/RP_ELASTIC_theta.pdf").c_str());
    delete c1;

    TH2D *RP_ELASTIC_p = (TH2D *)histograms->Get("RP_ELASTIC_p");
    RP_ELASTIC_p->SetTitle(";p_{X};p_{Y}");
    c1 = new TCanvas("c1", "c1", 2200, 1400);
    c1->SetMargin(0.10, 0.10, 0.1, 0.05);
    x_0 = 5.06472e-03;
    sigma_x = 3.42004e-02;
    y_0 = 5.98219e-04;
    sigma_y = 3.15726e-02;
    TEllipse *elasticP = new TEllipse(x_0, y_0, 3*sigma_x, 3*sigma_y);
    elasticP->SetLineColor(kRed);
    elasticP->SetFillColorAlpha(0, 0);
    RP_ELASTIC_p->Draw("colz");
    c1->Update();
    elasticP->Draw("same");
    c1->SetRealAspectRatio();
    c1->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/master_thesis/RP_ELASTIC_p.pdf").c_str());
    delete c1;



    TH1D *TRACKS_TOF = (TH1D *)histograms->Get("TRACKS_TOF");
    TH1D *TRACKS_TOF_OK = (TH1D *)histograms->Get("TRACKS_TOF_OK");
    TRACKS_TOF->SetTitle(";Number of tracks per event;Number of events");
    TRACKS_TOF_OK->SetTitle(";Number of tracks per event;Number of events");
    c1 = new TCanvas("c1", "c1", 1600, 1000);
    c1->SetRightMargin(0.05);
    TRACKS_TOF->Draw();
    TRACKS_TOF_OK->SetLineColor(kRed);
    TRACKS_TOF_OK->Draw("same");
    TLegend *legendToF = new TLegend(0.6, 0.7, 0.89, 0.79);
    legendToF->AddEntry("TRACKS_TOF", "All tracks", "l");
    legendToF->AddEntry("TRACKS_TOF_OK", "Tracks with TOF flag", "l");
    legendToF->SetBorderSize(0);
    legendToF->Draw("SAME");
    c1->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/master_thesis/TRACKS_TOF.pdf").c_str());
    delete c1;

    TH1D *TRACKS_PT = (TH1D *)histograms->Get("TRACKS_PT");
    TRACKS_PT->SetTitle(";p_{T} [GeV/c];Number of tracks");
    c1 = new TCanvas("c1", "c1", 1600, 1000);
    c1->SetRightMargin(0.05);
    TRACKS_PT->Draw();
    c1->Update();
    TLine *lineP = new TLine(0.2, 0, 0.2, c1->GetUymax());
    lineP->SetLineStyle(kDashed);
    lineP->SetLineColor(kRed);
    lineP->SetLineWidth(1);
    lineP->Draw("same");
    c1->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/master_thesis/TRACKS_PT.pdf").c_str());
    delete c1;

    TH1D *TRACKS_ETA = (TH1D *)histograms->Get("TRACKS_ETA");
    TRACKS_ETA->SetTitle(";#eta;Number of tracks");
    c1 = new TCanvas("c1", "c1", 1600, 1000);
    c1->SetRightMargin(0.05);
    TRACKS_ETA->Draw();
    c1->Update();
    TLine *lineEta1 = new TLine(-0.9, 0, -0.9, c1->GetUymax());
    lineEta1->SetLineStyle(kDashed);
    lineEta1->SetLineColor(kRed);
    lineEta1->SetLineWidth(1);
    lineEta1->Draw("same");
    TLine *lineEta2 = new TLine(0.9, 0, 0.9, c1->GetUymax());
    lineEta2->SetLineStyle(kDashed);
    lineEta2->SetLineColor(kRed);
    lineEta2->SetLineWidth(1);
    lineEta2->Draw("same");
    c1->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/master_thesis/TRACKS_ETA.pdf").c_str());
    delete c1;

    TH1D *TRACKS_NHITS = (TH1D *)histograms->Get("TRACKS_NHITS");
    TRACKS_NHITS->SetTitle(";Number of hits;Number of tracks");
    c1 = new TCanvas("c1", "c1", 1600, 1000);
    c1->SetLogy();
    c1->SetRightMargin(0.05);
    TRACKS_NHITS->Draw();
    c1->Update();
    // TLine *lineNhits = new TLine(20, c1->GetUymin(), 20, c1->GetUymax());
    TLine *lineNhits = new TLine(20, 3, 20, 4.5e5);
    lineNhits->SetLineStyle(kDashed);
    lineNhits->SetLineColor(kRed);
    lineNhits->SetLineWidth(1);
    lineNhits->Draw("same");
    c1->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/master_thesis/TRACKS_NHITS.pdf").c_str());
    delete c1;



    TH1D *PAIRS_DCADAUGHTERS_K0 = (TH1D *)histograms->Get("PAIRS_DCADAUGHTERS_K0");
    TH1D *PAIRS_DCADAUGHTERS_Lambda0 = (TH1D *)histograms->Get("PAIRS_DCADAUGHTERS_Lambda0");
    PAIRS_DCADAUGHTERS_K0->SetTitle(";DCA of decay products [cm];particles");
    PAIRS_DCADAUGHTERS_Lambda0->SetTitle(";DCA of decay products [cm];particles");
    c1 = new TCanvas("c1", "c1", 1600, 1000);
    c1->SetRightMargin(0.05);
    PAIRS_DCADAUGHTERS_K0->Draw();
    PAIRS_DCADAUGHTERS_Lambda0->SetLineColor(kBlack);
    PAIRS_DCADAUGHTERS_Lambda0->Draw("same");
    TLegend *legendDcaDaughters = new TLegend(0.7, 0.7, 0.89, 0.89);
    legendDcaDaughters->AddEntry("PAIRS_DCADAUGHTERS_K0", "K^{0}_{S} particles", "l");
    legendDcaDaughters->AddEntry("PAIRS_DCADAUGHTERS_Lambda0", "#Lambda^{0} particles", "l");
    legendDcaDaughters->SetBorderSize(0);
    legendDcaDaughters->Draw("SAME");
    c1->Update();
    TLine *lineDcaDaughters = new TLine(2.5, 0, 2.5, c1->GetUymax());
    lineDcaDaughters->SetLineStyle(kDashed);
    lineDcaDaughters->SetLineColor(kRed);
    lineDcaDaughters->SetLineWidth(1);
    lineDcaDaughters->Draw("same");
    c1->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/master_thesis/PAIRS_DCADAUGHTERS.pdf").c_str());
    delete c1;

    TH1D *PAIRS_DCABEAMLINE_K0 = (TH1D *)histograms->Get("PAIRS_DCABEAMLINE_K0");
    TH1D *PAIRS_DCABEAMLINE_Lambda0 = (TH1D *)histograms->Get("PAIRS_DCABEAMLINE_Lambda0");
    PAIRS_DCABEAMLINE_K0->SetTitle(";DCA of reconstructions to the beamline [cm];particles");
    PAIRS_DCABEAMLINE_Lambda0->SetTitle(";DCA of reconstructions to the beamline [cm];particles");
    c1 = new TCanvas("c1", "c1", 1600, 1000);
    c1->SetRightMargin(0.05);
    PAIRS_DCABEAMLINE_K0->Draw();
    PAIRS_DCABEAMLINE_Lambda0->SetLineColor(kBlack);
    PAIRS_DCABEAMLINE_Lambda0->Draw("same");
    TLegend *legendDcaBeamline = new TLegend(0.7, 0.7, 0.89, 0.89);
    legendDcaBeamline->AddEntry("PAIRS_DCABEAMLINE_K0", "K^{0}_{S} particles", "l");
    legendDcaBeamline->AddEntry("PAIRS_DCABEAMLINE_Lambda0", "#Lambda^{0} particles", "l");
    legendDcaBeamline->SetBorderSize(0);
    legendDcaBeamline->Draw("SAME");
    c1->Update();
    TLine *lineDcaBeamline = new TLine(2.5, 0, 2.5, c1->GetUymax());
    lineDcaBeamline->SetLineStyle(kDashed);
    lineDcaBeamline->SetLineColor(kRed);
    lineDcaBeamline->SetLineWidth(1);
    lineDcaBeamline->Draw("same");
    c1->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/master_thesis/PAIRS_DCABEAMLINE.pdf").c_str());
    delete c1;

    TH1D *PAIRS_DECAYLENGTH_K0_LEN = (TH1D *)histograms->Get("PAIRS_DECAYLENGTH_K0_LEN");
    TH1D *PAIRS_DECAYLENGTH_Lambda0_LEN = (TH1D *)histograms->Get("PAIRS_DECAYLENGTH_Lambda0_LEN");
    PAIRS_DECAYLENGTH_K0_LEN->SetTitle(";Decay length of reconstructed particles [cm];particles");
    PAIRS_DECAYLENGTH_Lambda0_LEN->SetTitle(";Decay length of reconstructed particles [cm];particles");
    c1 = new TCanvas("c1", "c1", 1600, 1000);
    c1->SetRightMargin(0.05);
    PAIRS_DECAYLENGTH_K0_LEN->Draw();
    PAIRS_DECAYLENGTH_Lambda0_LEN->SetLineColor(kBlack);
    PAIRS_DECAYLENGTH_Lambda0_LEN->Draw("same");
    TLegend *legendDecayLength = new TLegend(0.7, 0.7, 0.89, 0.89);
    legendDecayLength->AddEntry("PAIRS_DECAYLENGTH_K0_LEN", "K^{0}_{S} particles", "l");
    legendDecayLength->AddEntry("PAIRS_DECAYLENGTH_Lambda0_LEN", "#Lambda^{0} particles", "l");
    legendDecayLength->SetBorderSize(0);
    legendDecayLength->Draw("SAME");
    c1->Update();
    TLine *lineDecayLength = new TLine(3.0, 0, 3.0, c1->GetUymax());
    lineDecayLength->SetLineStyle(kDashed);
    lineDecayLength->SetLineColor(kRed);
    lineDecayLength->SetLineWidth(1);
    lineDecayLength->Draw("same");
    c1->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/master_thesis/PAIRS_DECAYLENGTH_LEN.pdf").c_str());
    delete c1;

    TH1D *PAIRS_DECAYLENGTH_K0_ANGLE = (TH1D *)histograms->Get("PAIRS_DECAYLENGTH_K0_ANGLE");
    TH1D *PAIRS_DECAYLENGTH_Lambda0_ANGLE = (TH1D *)histograms->Get("PAIRS_DECAYLENGTH_Lambda0_ANGLE");
    PAIRS_DECAYLENGTH_K0_ANGLE->SetTitle(";Decay length of reconstructed particles [cm];particles");
    PAIRS_DECAYLENGTH_Lambda0_ANGLE->SetTitle(";Decay length of reconstructed particles [cm];particles");
    c1 = new TCanvas("c1", "c1", 1600, 1000);
    c1->SetRightMargin(0.05);
    PAIRS_DECAYLENGTH_K0_ANGLE->Draw();
    PAIRS_DECAYLENGTH_Lambda0_ANGLE->SetLineColor(kBlack);
    PAIRS_DECAYLENGTH_Lambda0_ANGLE->Draw("same");
    TLegend *legendDecayAngle = new TLegend(0.7, 0.7, 0.89, 0.89);
    legendDecayAngle->AddEntry("PAIRS_DECAYLENGTH_K0_ANGLE", "K^{0}_{S} particles", "l");
    legendDecayAngle->AddEntry("PAIRS_DECAYLENGTH_Lambda0_ANGLE", "#Lambda^{0} particles", "l");
    legendDecayAngle->SetBorderSize(0);
    legendDecayAngle->Draw("SAME");
    c1->Update();
    TLine *lineDecayAngle = new TLine(0.925, 0, 0.925, c1->GetUymax());
    lineDecayAngle->SetLineStyle(kDashed);
    lineDecayAngle->SetLineColor(kRed);
    lineDecayAngle->SetLineWidth(1);
    lineDecayAngle->Draw("same");
    c1->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/master_thesis/PAIRS_DECAYLENGTH_ANGLE.pdf").c_str());
    delete c1;

    TH2D *PAIRS_DECAYLENGTH_K0 = (TH2D *)histograms->Get("PAIRS_DECAYLENGTH_K0");
    TH2D *PAIRS_DECAYLENGTH_Lambda0 = (TH2D *)histograms->Get("PAIRS_DECAYLENGTH_Lambda0");
    PAIRS_DECAYLENGTH_K0->SetTitle(";p_{x} [GeV/c];p_{y} [GeV/c]");
    PAIRS_DECAYLENGTH_Lambda0->SetTitle(";p_{x} [GeV/c];p_{y} [GeV/c]");
    c1 = new TCanvas("c1", "c1", 1600, 1200);
    c1->SetMargin(0.15, 0.15, 0.1, 0.05);
    PAIRS_DECAYLENGTH_K0->Draw("colz");
    c1->Update();
    TLine *lineHorizontal = new TLine(3., 0.925, c1->GetUxmax(), 0.925);
    lineHorizontal->SetLineStyle(kDashed);
    lineHorizontal->SetLineColor(kRed);
    lineHorizontal->SetLineWidth(2);
    lineHorizontal->Draw("same");
    TLine *lineVertical = new TLine(3., c1->GetUymin(), 3., 0.925);
    lineVertical->SetLineStyle(kDashed);
    lineVertical->SetLineColor(kRed);
    lineVertical->SetLineWidth(2);
    lineVertical->Draw("same");
    c1->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/master_thesis/PAIRS_DECAYLENGTH_K0.pdf").c_str());
    c1->Clear();
    PAIRS_DECAYLENGTH_Lambda0->Draw("colz");
    c1->Update();
    lineHorizontal->Draw("same");
    lineVertical->Draw("same");
    c1->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/master_thesis/PAIRS_DECAYLENGTH_Lambda0.pdf").c_str());
    delete c1;



    //Other things
    //luminosity of RP_CPT2 trigger
    ifstream inputFile("/home/adam/STAR-Analysis/share/lum_perrun_RP_CPT2.txt");
    std::string line, col2, col3, col4, col5;
    int col1;
    double col6;
    std::vector<int> runNumber;
    std::vector<double> prescale;
    while(getline(inputFile, line)){
        std::istringstream ss(line);
        ss>>col1>>col2>>col3>>col4>>col5>>col6;
        runNumber.push_back(col1);
        prescale.push_back(col6);
    }
    inputFile.close();
    std::vector<int> numberOfRuns;
    numberOfRuns.push_back(0);
    std::vector<std::string> etiquette;
    etiquette.push_back("18057---");
    int tempRunNumber = 18057;
    for(size_t i = 0; i<runNumber.size(); i++){
        // if(runNumber[i]/1000==tempRunNumber){
        if(runNumber[i]/1000==tempRunNumber){
            numberOfRuns[numberOfRuns.size()-1]++;
        } else{
            numberOfRuns.push_back(1);
            etiquette.push_back(to_string(runNumber[i]/1000)+"---");
            tempRunNumber = runNumber[i]/1000;
        }
    }

    //filling hitogram and dealing with the axis
    TH1D prescaleHist("prescaleHist", "RP_CPT2 trigger prescaling;fill number;prescale", runNumber.size(), 0, runNumber.size());
    for(size_t i = 0; i<runNumber.size(); i++){
        prescaleHist.SetBinContent(i+1, prescale[i]);
    }
    double prescaleAxisBinLimits[numberOfRuns.size()+1];
    prescaleAxisBinLimits[0] = 0;
    for(size_t i = 0; i<numberOfRuns.size(); i++){
        prescaleAxisBinLimits[i+1] = accumulate(numberOfRuns.begin(), numberOfRuns.begin()+i, 0);
    }
    TH1D prescaleAxisHist("prescaleAxisHist", "RP_CPT2 trigger prescaling;fill number;prescale", numberOfRuns.size(), prescaleAxisBinLimits);
    //0.9 to get rid of the last label covering up axis description
    for(size_t i = 0; i<numberOfRuns.size()*0.9; i++){
        if(i%10==0){
            prescaleAxisHist.GetXaxis()->ChangeLabel(i+1, 60, -1, -1, -1, -1, etiquette[i].c_str());

        } else{
            prescaleAxisHist.GetXaxis()->SetBinLabel(i+1, "");
        }
    }

    //drawing
    c1 = new TCanvas("c1", "c1", 1600, 1000);
    c1->SetRightMargin(0.05);
    c1->SetLogy();
    prescaleAxisHist.SetMaximum(100);
    prescaleAxisHist.SetMinimum(0.5);
    prescaleAxisHist.Draw("axis");
    prescaleHist.Draw("ah same");
    c1->Update();
    TLine *lineOne = new TLine(c1->GetUxmin(), 1, c1->GetUxmax(), 1);
    lineOne->SetLineStyle(kDashed);
    lineOne->SetLineColor(kRed);
    lineOne->SetLineWidth(1);
    lineOne->Draw("same");
    TLine *lineVert = new TLine(603, c1->GetUymin(), 603, 100);
    lineVert->SetLineStyle(kDashed);
    lineVert->SetLineColor(kRed);
    lineVert->SetLineWidth(1);
    lineVert->Draw("same");
    TLatex *textOne = new TLatex(590, 2.6, "up to 18083025");
    textOne->SetTextAngle(90);
    textOne->SetTextColor(kRed);
    textOne->SetTextSize(0.03);
    textOne->Draw("same");
    TLatex *textTwo = new TLatex(650, 2.6, "over 18083025");
    textTwo->SetTextAngle(90);
    textTwo->SetTextColor(kRed);
    textTwo->SetTextSize(0.03);
    textTwo->Draw("same");
    c1->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/master_thesis/prescaleHist.pdf").c_str());
    delete c1;


    TH1D *REMOVING_DUPLICATES_K0_BEFORE = (TH1D *)histograms->Get("REMOVING_DUPLICATES_K0_BEFORE");
    TH1D *REMOVING_DUPLICATES_K0_AFTER = (TH1D *)histograms->Get("REMOVING_DUPLICATES_K0_AFTER");
    TH1D *REMOVING_DUPLICATES_K0_DIFFERENCE = new TH1D();
    REMOVING_DUPLICATES_K0_BEFORE->Copy(*REMOVING_DUPLICATES_K0_DIFFERENCE);
    REMOVING_DUPLICATES_K0_DIFFERENCE->Add(REMOVING_DUPLICATES_K0_AFTER, -1.);
    REMOVING_DUPLICATES_K0_DIFFERENCE->SetName("REMOVING_DUPLICATES_K0_DIFFERENCE");
    c1 = new TCanvas("c1", "c1", 1600, 1000);
    c1->SetRightMargin(0.05);
    REMOVING_DUPLICATES_K0_BEFORE->SetMinimum(0);
    REMOVING_DUPLICATES_K0_BEFORE->SetTitle(";m_{#pi^{+}#pi^{-}} [GeV];reconstructed K^{0}_{S}");
    REMOVING_DUPLICATES_K0_BEFORE->SetLineColor(kBlack);
    REMOVING_DUPLICATES_K0_BEFORE->Draw();
    REMOVING_DUPLICATES_K0_AFTER->SetTitle(";m_{#pi^{+}#pi^{-}} [GeV];reconstructed K^{0}_{S}");
    REMOVING_DUPLICATES_K0_AFTER->SetLineColor(kBlue);
    REMOVING_DUPLICATES_K0_AFTER->Draw("same");
    REMOVING_DUPLICATES_K0_DIFFERENCE->SetTitle(";m_{#pi^{+}#pi^{-}} [GeV];reconstructed K^{0}_{S}");
    REMOVING_DUPLICATES_K0_DIFFERENCE->SetLineColor(kRed);
    REMOVING_DUPLICATES_K0_DIFFERENCE->Draw("same");
    TLegend *legendRemovingDuplicateK0 = new TLegend(0.2, 0.65, 0.4, 0.85);
    legendRemovingDuplicateK0->AddEntry("REMOVING_DUPLICATES_K0_BEFORE", "Before removal", "l");
    legendRemovingDuplicateK0->AddEntry("REMOVING_DUPLICATES_K0_AFTER", "After removal", "l");
    legendRemovingDuplicateK0->AddEntry("REMOVING_DUPLICATES_K0_DIFFERENCE", "Difference", "l");
    legendRemovingDuplicateK0->SetBorderSize(0);
    legendRemovingDuplicateK0->Draw("SAME");
    c1->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/master_thesis/REMOVING_DUPLICATES_K0.pdf").c_str());
    delete c1;

    TH1D *REMOVING_DUPLICATES_Lambda0_BEFORE = (TH1D *)histograms->Get("REMOVING_DUPLICATES_Lambda0_BEFORE");
    TH1D *REMOVING_DUPLICATES_Lambda0_AFTER = (TH1D *)histograms->Get("REMOVING_DUPLICATES_Lambda0_AFTER");
    TH1D *REMOVING_DUPLICATES_Lambda0_DIFFERENCE = new TH1D();
    REMOVING_DUPLICATES_Lambda0_BEFORE->Copy(*REMOVING_DUPLICATES_Lambda0_DIFFERENCE);
    REMOVING_DUPLICATES_Lambda0_DIFFERENCE->Add(REMOVING_DUPLICATES_Lambda0_AFTER, -1.);
    REMOVING_DUPLICATES_Lambda0_DIFFERENCE->SetName("REMOVING_DUPLICATES_Lambda0_DIFFERENCE");
    c1 = new TCanvas("c1", "c1", 1600, 1000);
    c1->SetRightMargin(0.05);
    REMOVING_DUPLICATES_Lambda0_BEFORE->SetMinimum(0);
    REMOVING_DUPLICATES_Lambda0_BEFORE->SetTitle(";m_{#pi^{+}#pi^{-}} [GeV];reconstructed K^{0}_{S}");
    REMOVING_DUPLICATES_Lambda0_BEFORE->SetLineColor(kBlack);
    REMOVING_DUPLICATES_Lambda0_BEFORE->Draw();
    REMOVING_DUPLICATES_Lambda0_AFTER->SetTitle(";m_{#pi^{+}#pi^{-}} [GeV];reconstructed K^{0}_{S}");
    REMOVING_DUPLICATES_Lambda0_AFTER->SetLineColor(kBlue);
    REMOVING_DUPLICATES_Lambda0_AFTER->Draw("same");
    REMOVING_DUPLICATES_Lambda0_DIFFERENCE->SetTitle(";m_{#pi^{+}#pi^{-}} [GeV];reconstructed K^{0}_{S}");
    REMOVING_DUPLICATES_Lambda0_DIFFERENCE->SetLineColor(kRed);
    REMOVING_DUPLICATES_Lambda0_DIFFERENCE->Draw("same");
    TLegend *legendRemovingDuplicateLambda0 = new TLegend(0.2, 0.65, 0.4, 0.85);
    legendRemovingDuplicateLambda0->AddEntry("REMOVING_DUPLICATES_Lambda0_BEFORE", "Before removal", "l");
    legendRemovingDuplicateLambda0->AddEntry("REMOVING_DUPLICATES_Lambda0_AFTER", "After removal", "l");
    legendRemovingDuplicateLambda0->AddEntry("REMOVING_DUPLICATES_Lambda0_DIFFERENCE", "Difference", "l");
    legendRemovingDuplicateLambda0->SetBorderSize(0);
    legendRemovingDuplicateLambda0->Draw("SAME");
    c1->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/master_thesis/REMOVING_DUPLICATES_Lambda0.pdf").c_str());
    delete c1;

    return 0;
}