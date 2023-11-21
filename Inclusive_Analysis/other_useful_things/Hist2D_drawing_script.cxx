#include "TCanvas.h"
#include "TROOT.h"
#include "TH2.h"
#include "TEfficiency.h"
#include "TF1.h"
#include "TLine.h"

int Hist2D_drawing_script(TEfficiency *TEff){
    TCanvas c1("c1", "", 1600, 900);
    c1.SetLeftMargin(0.15);
    double upper = 0.08;
    double lower = 0.005;
    double K = log(lower)/log(upper);

    //1st cut
    TLine cut1a = TLine(log10(lower), log10(upper), log10(upper), log10(upper));
    cut1a.SetLineWidth(1);
    cut1a.SetLineColor(2);
    TLine cut1b = TLine(log10(upper), log10(upper), log10(upper), log10(lower));
    cut1b.SetLineWidth(1);
    cut1b.SetLineColor(2);
    TLine cut1c = TLine(log10(upper), log10(lower), log10(lower), log10(lower));
    cut1c.SetLineWidth(1);
    cut1c.SetLineColor(2);
    TLine cut1d = TLine(log10(lower), log10(lower), log10(lower), log10(upper));
    cut1d.SetLineWidth(1);
    cut1d.SetLineColor(2);

    //2nd cut
    TF1 cut2a = TF1("cut2a", "2*log10([0])-x", -1.75, -0.5);
    cut2a.SetParameter(0, upper);
    cut2a.SetLineWidth(1);
    cut2a.SetLineColor(3);
    TF1 cut2b = TF1("cut2b", "2*log10([0])-x", -2.9, -1.75);
    cut2b.SetParameter(0, lower);
    cut2b.SetLineWidth(1);
    cut2b.SetLineColor(3);
    TF1 cut2c = TF1("cut2c", "log10(exp(log(pow(10, x))+[0]))", -2.9, -1.3);
    cut2c.SetParameter(0, K);
    cut2c.SetLineWidth(1);
    cut2c.SetLineColor(3);
    TF1 cut2d = TF1("cut2d", "log10(exp(log(pow(10, x))-[0]))", -1.9, -0.5);
    cut2d.SetParameter(0, K);
    cut2d.SetLineWidth(1);
    cut2d.SetLineColor(3);

    TEff->Draw("colz");
    TEff->GetPaintedHistogram()->SetMinimum(0.);
    TEff->GetPaintedHistogram()->SetMaximum(0.7);
    cut1a.Draw("same");
    cut1b.Draw("same");
    cut1c.Draw("same");
    cut1d.Draw("same");
    cut2a.Draw("same");
    cut2b.Draw("same");
    cut2c.Draw("same");
    cut2d.Draw("same");
    c1.SetRealAspectRatio();
    c1.Update();
    c1.Print((std::string(TEff->GetName())+"ratio.pdf").c_str());
    return 0;
}