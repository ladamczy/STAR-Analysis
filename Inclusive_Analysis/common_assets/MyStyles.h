#include "TStyle.h"

class MyStyles{
private:
    TStyle currentMyStyle;
    void SetDefault();
public:
    MyStyles(/* args */);
    ~MyStyles();
    TStyle* GetPointer();
    TStyle Hist2DNormalSize(bool containsTitle = true);
    TStyle Hist2DQuarterSize(bool containsTitle = true);
};

MyStyles::MyStyles(/* args */){
    this->SetDefault();
}

MyStyles::~MyStyles(){}

void MyStyles::SetDefault(){
    TStyle* modernStyle = gROOT->GetStyle("Modern");
    currentMyStyle = TStyle("My style", "A style with some of my modifications");
    modernStyle->Copy(currentMyStyle);
    //changes from the default ("Modern") style
    //pad
    currentMyStyle.SetPadLeftMargin(0.1);
    currentMyStyle.SetPadRightMargin(0.05);
    currentMyStyle.SetPadBottomMargin(0.12);
    currentMyStyle.SetPadTopMargin(0.10);
    //histograms
    currentMyStyle.SetOptStat(0);
    currentMyStyle.SetLegendBorderSize(0);
    currentMyStyle.SetLegendTextSize(0.04);
    currentMyStyle.SetTitleFontSize(0.06);
    currentMyStyle.SetTitleX(0.55);
    currentMyStyle.SetHistLineWidth(2);
    //axis
    currentMyStyle.SetAxisMaxDigits(3);
    currentMyStyle.SetLabelSize(0.045, "xyz");
    currentMyStyle.SetTitleSize(0.045, "xyz");
    currentMyStyle.SetTitleOffset(1.01, "Y");
}

TStyle* MyStyles::GetPointer(){
    return &currentMyStyle;
}

TStyle MyStyles::Hist2DNormalSize(bool containsTitle){
    this->SetDefault();
    if(!containsTitle){
        currentMyStyle.SetPadTopMargin(0.05);
    }
    return currentMyStyle;
}

TStyle MyStyles::Hist2DQuarterSize(bool containsTitle){
    this->SetDefault();
    if(!containsTitle){
        currentMyStyle.SetPadTopMargin(0.05);
    }
    //pad
    currentMyStyle.SetPadLeftMargin(0.15);
    currentMyStyle.SetPadBottomMargin(0.15);
    //histograms
    currentMyStyle.SetLegendTextSize(0.04);
    currentMyStyle.SetTitleFontSize(0.07);
    //axis
    currentMyStyle.SetLabelSize(0.06, "xyz");
    currentMyStyle.SetTitleSize(0.06, "xyz");
    return currentMyStyle;
}
