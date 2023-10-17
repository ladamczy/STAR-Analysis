#include "TStyle.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"

//here you can add your own after corresponding number
//so far it will be one number, which will turn on all styles

//style applied to gStyle
void MainStyle(TStyle* style, int styleNumber=0){
    switch (styleNumber){
    case 0:
        style->SetHistLineWidth(1);
        style->SetFrameLineWidth(3);
        style->SetOptStat(111111);
        style->SetStatY(0.89);
        style->SetStatX(0.89);
        style->SetStatW(0.1);
        style->SetStatH(0.1); 
        break;
    case 1:
        //your style here 
        break;
    
    default:
        break;
    }
}

//style applied to 1 dimensional histograms
void TH1Style(TCanvas* c, TH1D* onedimhist, int styleNumber=0){
    switch (styleNumber){
    case 0:
        // temp1D->SetMarkerStyle(kFullCircle);
        // temp1D->Draw("E");
        c->SetLogx(0);
        c->SetLogy(1);
        c->SetLogz(0);
        onedimhist->Draw("HIST");
        break;
    case 1:
        //your style here 
        break;
    
    default:
        break;
    }
}

//style applied to 2 dimensional histograms
void TH2Style(TCanvas* c, TH2D* twodimhist, int styleNumber=0){
    switch (styleNumber){
    case 0:
        c->SetLogx(0);
        c->SetLogy(0);
        c->SetLogz(1);
        twodimhist->Draw("COLZ");
        break;
    case 1:
        //your style here 
        break;
    
    default:
        break;
    }
}