#include "TStyle.h"
#include "TROOT.h"

class MyStyles{
private:
    static int internalCanvasCounter;
    TStyle currentMyStyle;
    void SetDefault();
public:
    MyStyles(/* args */);
    ~MyStyles();
    //default canvas source, 16:9 ratio
    static TCanvas* DefaultCanvas(std::string name = "", std::string title = "");
    TStyle* GetPointer();
    TStyle Hist2DDisplay(bool containsTitle = true);
    TStyle Hist2DNormalSize(bool containsTitle = true);
    TStyle Hist2DQuarterSize(bool containsTitle = true);
};

//initialisation of canvas counter
int MyStyles::internalCanvasCounter{ 1 };

//constructor & destructor

MyStyles::MyStyles(/* args */){
    this->SetDefault();
}

MyStyles::~MyStyles(){}

//all the other functions

TCanvas* MyStyles::DefaultCanvas(std::string name, std::string title){
    //if name is empty, assign default one (c#)
    if(name.length()==0){
        name = "c"+std::to_string(internalCanvasCounter);
        internalCanvasCounter++;
    }
    //if title is empty, just copy the name of the canvas
    if(title.length()==0){
        title = name;
    }

    //return new canvas with default values (16:9, no menu)
    return new TCanvas(name.c_str(), title.c_str(), -1, 0, 1600, 900);
}

void MyStyles::SetDefault(){
    // TStyle* modernStyle = gROOT->GetStyle("Modern");
    // currentMyStyle = TStyle("My style", "A style with some of my modifications");
    // modernStyle->Copy(currentMyStyle);
    TStyle* modernStyle = gROOT->GetStyle("Modern");
    currentMyStyle = TStyle();
    modernStyle->Copy(currentMyStyle);
    currentMyStyle.SetNameTitle("MyBaseStyle", "A style with some of my modifications");
    //changes from the default ("Modern") style
    //pad
    currentMyStyle.SetPadLeftMargin(0.1);
    currentMyStyle.SetPadRightMargin(0.05);
    currentMyStyle.SetPadBottomMargin(0.1);
    currentMyStyle.SetPadTopMargin(0.1);
    currentMyStyle.SetFrameLineWidth(2);
    //histograms
    currentMyStyle.SetOptStat(0);
    currentMyStyle.SetLegendBorderSize(0);
    currentMyStyle.SetLegendTextSize(0.04);
    currentMyStyle.SetTitleFontSize(0.06);
    currentMyStyle.SetTitleX(0.55);
    currentMyStyle.SetHistLineWidth(1);
    currentMyStyle.SetHistLineColor(kBlue+2);
    currentMyStyle.SetMarkerStyle(kFullCircle);
    currentMyStyle.SetMarkerColor(kBlue);
    //axis
    currentMyStyle.SetAxisMaxDigits(3);
    currentMyStyle.SetLabelSize(0.045, "xyz");
    currentMyStyle.SetTitleSize(0.045, "xyz");
    currentMyStyle.SetTitleOffset(1.01, "Y");
}

TStyle* MyStyles::GetPointer(){
    return &currentMyStyle;
}

TStyle MyStyles::Hist2DDisplay(bool containsTitle){
    this->SetDefault();
    //additional, general settings
    currentMyStyle.SetNameTitle("2DDisplay", "A style to display (not save) things");
    currentMyStyle.SetHistMinimumZero();
    currentMyStyle.SetOptFit();
    currentMyStyle.SetStatBorderSize(0.);
    currentMyStyle.SetStatX(0.94);
    currentMyStyle.SetStatY(0.89);
    currentMyStyle.SetStatW(0.18);
    currentMyStyle.SetFitFormat("6.5g");
    currentMyStyle.SetStatColor(0);
    //fixing positions if slide has no title
    if(!containsTitle){
        currentMyStyle.SetPadTopMargin(0.05);
        currentMyStyle.SetStatY(0.94);
    }
    return currentMyStyle;
}

TStyle MyStyles::Hist2DNormalSize(bool containsTitle){
    this->SetDefault();
    currentMyStyle.SetNameTitle("2DNormalSize", "A style to save things to show at fullscreen");
    if(!containsTitle){
        currentMyStyle.SetPadTopMargin(0.05);
    }
    return currentMyStyle;
}

TStyle MyStyles::Hist2DQuarterSize(bool containsTitle){
    this->SetDefault();
    currentMyStyle.SetNameTitle("2DQuarterSize", "A style to save things to show at 1/4th of fullscreen");
    if(!containsTitle){
        currentMyStyle.SetPadTopMargin(0.05);
    }
    //pad
    currentMyStyle.SetPadLeftMargin(0.15);
    currentMyStyle.SetPadBottomMargin(0.15);
    //histograms
    currentMyStyle.SetLegendTextSize(0.04);
    currentMyStyle.SetTitleFontSize(0.07);
    currentMyStyle.SetHistLineWidth(2);
    currentMyStyle.SetEndErrorSize(4.);
    //axis
    currentMyStyle.SetLabelSize(0.06, "xyz");
    currentMyStyle.SetTitleSize(0.06, "xyz");
    return currentMyStyle;
}
