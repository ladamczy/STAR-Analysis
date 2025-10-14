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

int main(int argc, char** argv){
    TFile* histograms = TFile::Open(argv[1]);
    std::string folderWithDiagonal = std::string(argv[2]);
    if(folderWithDiagonal[folderWithDiagonal.size()-1]!='/'){
        folderWithDiagonal += "/";
    }

    gStyle->SetFrameLineWidth(1);
    // gStyle->SetOptFit(111);
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);

    //lines for fiducial
    TLine* fidVertOne, * fidVertTwo, * fidHorOne, * fidHorTwo, * fidHorThree, * fidHorFour;
    TEllipse* fidEllOne, * fidEllTwo;
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
    TH2D* RP_FIDUCIAL_east = (TH2D*)histograms->Get("RP_FIDUCIAL_east");
    TH2D* RP_FIDUCIAL_west = (TH2D*)histograms->Get("RP_FIDUCIAL_west");
    RP_FIDUCIAL_east->SetTitle(";p_{x} [GeV/c];p_{y} [GeV/c]");
    RP_FIDUCIAL_west->SetTitle(";p_{x} [GeV/c];p_{y} [GeV/c]");
    TCanvas* c1 = new TCanvas("c1", "c1", 1000, 1500);
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
    c1->SaveAs((folderWithDiagonal+"RP_FIDUCIAL_east.pdf").c_str());
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
    c1->SaveAs((folderWithDiagonal+"RP_FIDUCIAL_west.pdf").c_str());
    delete c1;
}