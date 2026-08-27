#include "TStyle.h"
#include "TROOT.h"

class MyStyles{
private:
    static int internalCanvasCounter;
    static TStyle currentMyStyle;
    static void SetDefault(TStyle& changedStyle);
public:
    MyStyles(/* args */);
    ~MyStyles();
    //default canvas source, 16:9 ratio
    static TCanvas* DefaultCanvas(std::string name = "", std::string title = "");
    TStyle* GetPointer();
    static TStyle Hist2DDisplay(bool containsTitle = true);
    static TStyle Hist2DNormalSize(bool containsTitle = true);
    static TStyle Hist2DQuarterSize(bool containsTitle = true);
};

//initialisation of canvas counter
int MyStyles::internalCanvasCounter{ 1 };
TStyle MyStyles::currentMyStyle{ TStyle() };

//constructor & destructor
MyStyles::MyStyles(/* args */){
    SetDefault(currentMyStyle);
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

void MyStyles::SetDefault(TStyle& changedStyle){
    TStyle* modernStyle = gROOT->GetStyle("Modern");
    modernStyle->Copy(changedStyle);
    changedStyle.SetNameTitle("MyBaseStyle", "A style with some of my modifications");
    //changes from the default ("Modern") style
    //pad
    changedStyle.SetPadLeftMargin(0.1);
    changedStyle.SetPadRightMargin(0.05);
    changedStyle.SetPadBottomMargin(0.1);
    changedStyle.SetPadTopMargin(0.1);
    changedStyle.SetFrameLineWidth(2);
    //histograms
    changedStyle.SetOptStat(0);
    changedStyle.SetTitleFontSize(0.06);
    changedStyle.SetTitleX(0.55);
    changedStyle.SetHistLineWidth(1);
    changedStyle.SetHistLineColor(kBlue+2);
    changedStyle.SetMarkerStyle(kFullCircle);
    changedStyle.SetMarkerColor(kBlue);
    //legend and stats
    currentMyStyle.SetLegendTextSize(0.03);
    currentMyStyle.SetLegendBorderSize(0);
    currentMyStyle.SetStatFontSize(0.03);
    currentMyStyle.SetStatBorderSize(0);
    //axis
    changedStyle.SetAxisMaxDigits(3);
    changedStyle.SetLabelSize(0.045, "xyz");
    changedStyle.SetTitleSize(0.045, "xyz");
    changedStyle.SetTitleOffset(1.01, "Y");
}

TStyle* MyStyles::GetPointer(){
    return &currentMyStyle;
}

TStyle MyStyles::Hist2DDisplay(bool containsTitle){
    SetDefault(currentMyStyle);
    //additional, general settings
    currentMyStyle.SetNameTitle("2DDisplay", "A style to display (not save) things");
    currentMyStyle.SetHistMinimumZero();
    currentMyStyle.SetOptFit();
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
    SetDefault(currentMyStyle);
    currentMyStyle.SetNameTitle("2DNormalSize", "A style to save things to show at fullscreen");
    if(!containsTitle){
        currentMyStyle.SetPadTopMargin(0.05);
    }
    return currentMyStyle;
}

TStyle MyStyles::Hist2DQuarterSize(bool containsTitle){
    SetDefault(currentMyStyle);
    currentMyStyle.SetNameTitle("2DQuarterSize", "A style to save things to show at 1/4th of fullscreen");
    if(!containsTitle){
        currentMyStyle.SetPadTopMargin(0.05);
    }
    //pad
    currentMyStyle.SetPadLeftMargin(0.15);
    currentMyStyle.SetPadBottomMargin(0.15);
    //histograms
    currentMyStyle.SetTitleFontSize(0.07);
    currentMyStyle.SetHistLineWidth(2);
    currentMyStyle.SetEndErrorSize(4.);
    //legend and stats
    currentMyStyle.SetLegendTextSize(0.04);
    currentMyStyle.SetStatFontSize(0.04);
    //axis
    currentMyStyle.SetLabelSize(0.06, "xyz");
    currentMyStyle.SetTitleSize(0.06, "xyz");
    return currentMyStyle;
}
