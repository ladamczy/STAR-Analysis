#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TEfficiency.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TROOT.h"

#include "string.h"
#include <cmath>

#include "MyStyles.h"

void draw_and_save(TH1D* data, TH1D* sim, std::string fileTitle, bool isYlogarithmic = false, double x1 = 0, double y1 = 0, double x2 = 0, double y2 = 0, double line1 = 0, double line2 = 0);
void simpler_draw_and_save(TH1* hist, std::string filename = "", bool is3D = false, bool isYlogarithmic = false, double line1 = 0, double line2 = 0);
void simpler_draw_and_save(TEfficiency* hist, std::string filename = "");

int main(int argc, char const* argv[]){
    //setting style
    // TStyle mystyle = MyStyles::Hist2DNormalSize(false);
    TStyle mystyle = MyStyles::Hist2DQuarterSize();
    mystyle.cd();
    gROOT->ForceStyle();

    //data taken using PhD_data_new folder
    TFile* newdata = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/AnaOutput_Inclusive_analysis_Kstar_phi_new_data_target_chi2_ONLY.root");
    //data taken using PhD_data_old folder, created by command:
    //time ~/star-upc/build/bin/Preselection_olddata_from_beginning ~/OneTouch/starlist.list ~/OneTouch/PhD_data_old/ 6
    TFile* olddata = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/OLDTEST.root");

    //simulation to compare new data to (and also source to +-1ns problem)
    TFile* KKsimudata = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/SIMUKK.root");
    TFile* Kpisimudata = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/SIMUKpi.root");

    //nsigma issue
    TFile* generalTimingData = TFile::Open("/home/adam/Pulpit/filelistGeneralPreselectedSampleMCIdentificationPhiHistograms.root");

    //pairs
    TH1D* MKKnew = (TH1D*)newdata->Get("MKKChi2");
    MKKnew->SetName("MKKChi2New");
    TH1D* MKKold = (TH1D*)olddata->Get("MKKChi2");
    MKKold->SetName("MKKChi2Old");
    TH1D* MKpinew = (TH1D*)newdata->Get("MKpiChi2");
    MKpinew->SetName("MKpiChi2New");
    TH1D* MKpiold = (TH1D*)olddata->Get("MKpiChi2");
    MKpiold->SetName("MKpiChi2Old");
    TH1D* MpiKnew = (TH1D*)newdata->Get("MpiKChi2");
    MpiKnew->SetName("MpiKChi2New");
    TH1D* MpiKold = (TH1D*)olddata->Get("MpiKChi2");
    MpiKold->SetName("MpiKChi2Old");

    //simulated numbers
    TH1D* MKKsimueta = ((TH2D*)KKsimudata->Get("MKKChi2Identifiedeta"))->ProjectionY();
    TH1D* MKKsimupT = ((TH2D*)KKsimudata->Get("MKKChi2IdentifiedpT"))->ProjectionY();
    //real numbers, from fitting, so there is no root file, just printout
    TH1D* MKKneweta = new TH1D("MKKneweta", "", 10, -1.0, 1.0);
    MKKneweta->SetBinContent(1, 26.21410);
    MKKneweta->SetBinContent(2, 169.1920);
    MKKneweta->SetBinContent(3, 243.4890);
    MKKneweta->SetBinContent(4, 246.4290);
    MKKneweta->SetBinContent(5, 215.3030);
    MKKneweta->SetBinContent(6, 234.4910);
    MKKneweta->SetBinContent(7, 270.5140);
    MKKneweta->SetBinContent(8, 245.070);
    MKKneweta->SetBinContent(9, 177.0690);
    MKKneweta->SetBinContent(10, 34.65790);
    MKKneweta->SetBinError(1, 5.46303);
    MKKneweta->SetBinError(2, 14.0211);
    MKKneweta->SetBinError(3, 17.4732);
    MKKneweta->SetBinError(4, 16.1021);
    MKKneweta->SetBinError(5, 14.1829);
    MKKneweta->SetBinError(6, 14.8291);
    MKKneweta->SetBinError(7, 18.4275);
    MKKneweta->SetBinError(8, 16.2559);
    MKKneweta->SetBinError(9, 15.9723);
    MKKneweta->SetBinError(10, 6.92361);
    TH1D* MKKnewpT = new TH1D("MKKnewpT", "", 12, 0.0, 2.4);
    MKKnewpT->SetBinContent(1, 0);
    MKKnewpT->SetBinContent(2, 0);
    MKKnewpT->SetBinContent(3, 37.2622);
    MKKnewpT->SetBinContent(4, 260.811);
    MKKnewpT->SetBinContent(5, 415.285);
    MKKnewpT->SetBinContent(6, 376.893);
    MKKnewpT->SetBinContent(7, 232.338);
    MKKnewpT->SetBinContent(8, 204.612);
    MKKnewpT->SetBinContent(9, 117.182);
    MKKnewpT->SetBinContent(10, 80.6733);
    MKKnewpT->SetBinContent(11, 60.2078);
    MKKnewpT->SetBinContent(12, 17.3023);
    MKKnewpT->SetBinError(1, 0);
    MKKnewpT->SetBinError(2, 0);
    MKKnewpT->SetBinError(3, 4.40634);
    MKKnewpT->SetBinError(4, 13.4975);
    MKKnewpT->SetBinError(5, 18.0108);
    MKKnewpT->SetBinError(6, 18.0827);
    MKKnewpT->SetBinError(7, 17.5927);
    MKKnewpT->SetBinError(8, 19.3959);
    MKKnewpT->SetBinError(9, 18.8536);
    MKKnewpT->SetBinError(10, 12.6103);
    MKKnewpT->SetBinError(11, 12.0681);
    MKKnewpT->SetBinError(12, 4.09201);
    //chi2 of real numbers
    TH1D* MKKnewetaChi2 = new TH1D("MKKnewetaChi2", "", 10, -1.0, 1.0);
    MKKnewetaChi2->SetBinContent(1, 1.15901);
    MKKnewetaChi2->SetBinContent(2, 0.765579);
    MKKnewetaChi2->SetBinContent(3, 1.70266);
    MKKnewetaChi2->SetBinContent(4, 1.34132);
    MKKnewetaChi2->SetBinContent(5, 0.733261);
    MKKnewetaChi2->SetBinContent(6, 1.16306);
    MKKnewetaChi2->SetBinContent(7, 1.27507);
    MKKnewetaChi2->SetBinContent(8, 1.22938);
    MKKnewetaChi2->SetBinContent(9, 0.836511);
    MKKnewetaChi2->SetBinContent(10, 0.964468);
    TH1D* MKKnewpTChi2 = new TH1D("MKKnewpTChi2", "", 12, 0.0, 2.4);
    MKKnewpTChi2->SetBinContent(1, 0);
    MKKnewpTChi2->SetBinContent(2, 0);
    MKKnewpTChi2->SetBinContent(3, 2.00508);
    MKKnewpTChi2->SetBinContent(4, 1.21756);
    MKKnewpTChi2->SetBinContent(5, 0.608935);
    MKKnewpTChi2->SetBinContent(6, 0.91743);
    MKKnewpTChi2->SetBinContent(7, 1.36787);
    MKKnewpTChi2->SetBinContent(8, 1.20226);
    MKKnewpTChi2->SetBinContent(9, 1.29679);
    MKKnewpTChi2->SetBinContent(10, 1.18959);
    MKKnewpTChi2->SetBinContent(11, 1.08141);
    MKKnewpTChi2->SetBinContent(12, 1.59401);

    //+-1ns issue histograms
    TH1D* deltaT0KK = (TH1D*)KKsimudata->Get("deltaT0KK");
    TH1D* deltaT0Kpi = (TH1D*)Kpisimudata->Get("deltaT0Kpi");
    TH1D* deltaT0piK = (TH1D*)Kpisimudata->Get("deltaT0piK");

    //nsigma issue histograms
    TH1D* nsigmaPionPositive = (TH1D*)generalTimingData->Get("NSigmaTestMCPionTPCPIDPionPositive");
    TH1D* nsigmaKaonPositive = (TH1D*)generalTimingData->Get("NSigmaTestMCKaonTPCPIDKaonPositive");
    TH1D* nsigmaProtonPositive = (TH1D*)generalTimingData->Get("NSigmaTestMCProtonTPCPIDProtonPositive");
    TH1D* nsigmaPionNegative = (TH1D*)generalTimingData->Get("NSigmaTestMCPionTPCPIDPionNegative");
    TH1D* nsigmaKaonNegative = (TH1D*)generalTimingData->Get("NSigmaTestMCKaonTPCPIDKaonNegative");
    TH1D* nsigmaProtonNegative = (TH1D*)generalTimingData->Get("NSigmaTestMCProtonTPCPIDProtonNegative");

    //preprocessing of some histograms
    //original had a bit different binning, so I'm fixing it: 0.015 GeV/bin, starting from 0.9 GeV, ending at 2.4
    MKKold->Rebin(5);
    MKpiold->Rebin(5);
    MpiKold->Rebin(5);
    //new has to also be rebinned this way, to look presentable
    MKKnew->Rebin(5);
    MKpinew->Rebin(5);
    MpiKnew->Rebin(5);
    //changing range of +-1ns issue histograms
    deltaT0KK->SetAxisRange(-2., 2.);
    deltaT0Kpi->SetAxisRange(-2., 2.);
    deltaT0piK->SetAxisRange(-2., 2.);
    //changing range of nsigma histograms
    nsigmaPionPositive->SetAxisRange(-5., 5.);
    nsigmaKaonPositive->SetAxisRange(-5., 5.);
    nsigmaProtonPositive->SetAxisRange(-5., 5.);
    nsigmaPionNegative->SetAxisRange(-5., 5.);
    nsigmaKaonNegative->SetAxisRange(-5., 5.);
    nsigmaProtonNegative->SetAxisRange(-5., 5.);

    //general settings
    TCanvas* resultCanvas = new TCanvas("resultCanvas", "resultCanvas", 1600, 900);
    resultCanvas->UseCurrentStyle();
    resultCanvas->SetLogy();
    resultCanvas->SetFrameLineWidth(2);
    // resultCanvas->SetMargin(0.1, 0.05, 0.1, 0.05);

    //################################### KK
    //new and setting the canvas
    MKKnew->SetTitle("K^{+}K^{-} total yield");
    MKKnew->SetMinimum(0.8);
    MKKnew->SetLineColor(kBlue+2);
    MKKnew->Draw("hist");
    //old
    MKKold->SetLineColor(kRed);
    MKKold->Draw("hist same");
    TLegend* legend = new TLegend(0.29, 0.24, 0.49, 0.44);
    legend->SetTextSize(0.025);
    legend->SetHeader("#bf{pp, #sqrt{s} = 510 GeV}", "C");
    legend->AddEntry(MKKold, "#bf{yield with previous cuts}");
    legend->AddEntry(MKKnew, "#bf{yield with current cuts}");
    legend->SetBorderSize(0);
    legend->DrawClone("SAME");
    resultCanvas->SaveAs("/home/adam/Pulpit/STAR_Presentation_custom/MKKComparison.pdf");
    resultCanvas->Clear();
    //################################### Kpi
    //new and setting the canvas
    MKpinew->SetTitle("K^{+}#pi^{-} total yield");
    MKpinew->SetMinimum(0.8);
    MKpinew->SetLineColor(kBlue+2);
    MKpinew->Draw("hist");
    //old
    MKpiold->SetLineColor(kRed);
    MKpiold->Draw("hist same");
    legend = new TLegend(0.34, 0.29, 0.54, 0.49);
    legend->SetTextSize(0.025);
    legend->SetHeader("#bf{pp, #sqrt{s} = 510 GeV}", "C");
    legend->AddEntry(MKpiold, "#bf{yield with previous cuts}");
    legend->AddEntry(MKpinew, "#bf{yield with current cuts}");
    legend->SetBorderSize(0);
    legend->DrawClone("SAME");
    resultCanvas->SaveAs("/home/adam/Pulpit/STAR_Presentation_custom/MKpiComparison.pdf");
    resultCanvas->Clear();
    delete legend;
    //################################### piK
    //new and setting the canvas
    MpiKnew->SetTitle("#pi^{+}K^{-} total yield");
    MpiKnew->SetMinimum(0.8);
    MpiKnew->SetLineColor(kBlue+2);
    MpiKnew->Draw("hist");
    //old
    MpiKold->SetLineColor(kRed);
    MpiKold->Draw("hist same");
    legend = new TLegend(0.34, 0.29, 0.54, 0.49);
    legend->SetTextSize(0.025);
    legend->SetHeader("#bf{pp, #sqrt{s} = 510 GeV}", "C");
    legend->AddEntry(MpiKold, "#bf{yield with previous cuts}");
    legend->AddEntry(MpiKnew, "#bf{yield with current cuts}");
    legend->SetBorderSize(0);
    legend->DrawClone("SAME");
    resultCanvas->SaveAs("/home/adam/Pulpit/STAR_Presentation_custom/MpiKComparison.pdf");
    resultCanvas->Clear();
    delete legend;


    resultCanvas->SetLogy(0);


    //################################### KK eta results + comparison with simulation
    //data
    MKKneweta->SetTitle("K^{+}K^{-} #eta-dependent yield;#eta;pairs");
    MKKneweta->SetMinimum(0.);
    MKKneweta->SetLineColor(kBlue+2);
    MKKneweta->Draw("e1");
    //simulation
    MKKsimueta->SetLineColor(kRed);
    //scaling to number of events before drawing
    MKKsimueta->Scale(MKKneweta->Integral()/MKKsimueta->Integral());
    MKKsimueta->Draw("hist same");
    legend = new TLegend(0.4, 0.29, 0.6, 0.49);
    legend->SetTextSize(0.025);
    legend->SetHeader("#bf{pp, #sqrt{s} = 510 GeV}", "C");
    legend->AddEntry(MKKsimueta, "#bf{yield from starsim (scaled)}");
    legend->AddEntry(MKKneweta, "#bf{yield from data}");
    legend->SetBorderSize(0);
    legend->DrawClone("SAME");
    resultCanvas->SaveAs("/home/adam/Pulpit/STAR_Presentation_custom/MKKetaResults.pdf");
    resultCanvas->Clear();
    delete legend;
    //################################### KK eta Chi2
    MKKnewetaChi2->SetTitle("K^{+}K^{-} #eta-dependent  #chi^{2}/ndf;#eta;#chi^{2}/ndf");
    MKKnewetaChi2->SetMinimum(0.);
    MKKnewetaChi2->SetLineColor(kBlue+2);
    MKKnewetaChi2->Draw("hist");
    TLine* chi2ndfOne = new TLine(resultCanvas->GetUxmin(), 1., resultCanvas->GetUxmax(), 1.);
    chi2ndfOne->SetLineColor(kRed);
    chi2ndfOne->SetLineStyle(kDashed);
    chi2ndfOne->SetLineWidth(2);
    chi2ndfOne->Draw("same");
    legend = new TLegend(0.6, 0.25, 0.8, 0.45);
    legend->SetTextSize(0.025);
    legend->SetHeader("#bf{pp, #sqrt{s} = 510 GeV}", "C");
    legend->AddEntry(chi2ndfOne, "#bf{#chi^{2}/ndf = 1}", "l");
    legend->AddEntry(MKKnewetaChi2, "#bf{#chi^{2}/ndf of data fit}");
    legend->SetBorderSize(0);
    legend->DrawClone("SAME");
    resultCanvas->SaveAs("/home/adam/Pulpit/STAR_Presentation_custom/MKKetaChi2.pdf");
    resultCanvas->Clear();
    delete legend;
    //################################### KK pT results + comparison with simulation
    double temp_top_margin = resultCanvas->GetTopMargin();
    resultCanvas->SetTopMargin(0.12);//because T of pT clips the upper edge
    //drawing scaled simulation forst cause it gets higher results
    //simulation
    MKKsimupT->SetTitle("K^{+}K^{-} p_{T}-dependent yield;p_{T} [GeV/c];pairs");
    MKKsimupT->SetMinimum(0.);
    //scaling to number of events before drawing
    MKKsimupT->Scale(MKKnewpT->Integral()/MKKsimupT->Integral());
    MKKsimupT->SetLineColor(kRed);
    MKKsimupT->Draw("hist");
    //data
    MKKnewpT->SetLineColor(kBlue+2);
    MKKnewpT->Draw("e1 same");
    legend = new TLegend(0.69, 0.67, 0.89, 0.87);
    legend->SetTextSize(0.025);
    legend->SetHeader("#bf{pp, #sqrt{s} = 510 GeV}", "C");
    legend->AddEntry(MKKsimupT, "#bf{yield from starsim (scaled)}");
    legend->AddEntry(MKKnewpT, "#bf{yield from data}");
    legend->SetBorderSize(0);
    legend->DrawClone("SAME");
    resultCanvas->SaveAs("/home/adam/Pulpit/STAR_Presentation_custom/MKKpTResults.pdf");
    resultCanvas->SetTopMargin(temp_top_margin); //bringing back the previous margin
    resultCanvas->Clear();
    delete legend;
    //################################### KK pT Chi2
    temp_top_margin = resultCanvas->GetTopMargin();
    resultCanvas->SetTopMargin(0.12);//because T of pT clips the upper edge
    MKKnewpTChi2->SetTitle("K^{+}K^{-} p_{T}-dependent  #chi^{2}/ndf;p_{T} [GeV/c];#chi^{2}/ndf");
    MKKnewpTChi2->SetMinimum(0.);
    MKKnewpTChi2->SetLineColor(kBlue+2);
    MKKnewpTChi2->Draw("hist");
    chi2ndfOne->SetX1(resultCanvas->GetUxmin());
    chi2ndfOne->SetX2(resultCanvas->GetUxmax());
    chi2ndfOne->SetLineColor(kRed);
    chi2ndfOne->SetLineStyle(kDashed);
    chi2ndfOne->SetLineWidth(2);
    chi2ndfOne->Draw("same");
    legend = new TLegend(0.6, 0.25, 0.8, 0.45);
    legend->SetTextSize(0.025);
    legend->SetHeader("#bf{pp, #sqrt{s} = 510 GeV}", "C");
    legend->AddEntry(chi2ndfOne, "#bf{#chi^{2}/ndf = 1}", "l");
    legend->AddEntry(MKKnewpTChi2, "#bf{#chi^{2}/ndf of data fit}");
    legend->SetBorderSize(0);
    legend->DrawClone("SAME");
    resultCanvas->SaveAs("/home/adam/Pulpit/STAR_Presentation_custom/MKKpTChi2.pdf");
    resultCanvas->SetTopMargin(temp_top_margin); //bringing back the previous margin
    resultCanvas->Clear();
    delete legend;



    //############################################ +-1ns issue, KK pair
    deltaT0KK->SetTitle("#Delta t_{0} of K^{+}K^{-} pair;#Delta t_{0} [ns];pairs");
    deltaT0KK->SetMinimum(0.);
    deltaT0KK->SetLineColor(kBlue+2);
    deltaT0KK->Draw("hist");
    resultCanvas->SaveAs("/home/adam/Pulpit/STAR_Presentation_custom/deltaT0KK.pdf");
    resultCanvas->Clear();
    //############################################ +-1ns issue, Kpi pair
    deltaT0Kpi->SetTitle("#Delta t_{0} of K^{+}#pi^{-} pair;#Delta t_{0} [ns];pairs");
    deltaT0Kpi->SetMinimum(0.);
    deltaT0Kpi->SetLineColor(kBlue+2);
    deltaT0Kpi->Draw("hist");
    resultCanvas->SaveAs("/home/adam/Pulpit/STAR_Presentation_custom/deltaT0Kpi.pdf");
    resultCanvas->Clear();
    //############################################ +-1ns issue, piK pair
    deltaT0piK->SetTitle("#Delta t_{0} of #pi^{+}K^{-} pair;#Delta t_{0} [ns];pairs");
    deltaT0piK->SetMinimum(0.);
    deltaT0piK->SetLineColor(kBlue+2);
    deltaT0piK->Draw("hist");
    resultCanvas->SaveAs("/home/adam/Pulpit/STAR_Presentation_custom/deltaT0piK.pdf");
    resultCanvas->Clear();



    //############################################ nsigma pi+
    nsigmaPionPositive->SetTitle("n#sigma_{#pi} for #pi^{+} tracks;n#sigma_{#pi};tracks");
    nsigmaPionPositive->SetMinimum(0.);
    nsigmaPionPositive->SetLineColor(kBlue+2);
    nsigmaPionPositive->Draw("hist");
    nsigmaPionPositive->SetNdivisions(520);
    resultCanvas->ModifiedUpdate(); //for Uymin and Uymax to update
    TLine* nsigmaZero = new TLine(0., resultCanvas->GetUymin(), 0., resultCanvas->GetUymax());
    nsigmaZero->SetLineColor(kRed);
    nsigmaZero->SetLineStyle(kDashed);
    nsigmaZero->SetLineWidth(2);
    nsigmaZero->Draw("same");
    TLine* nsigmaMean = new TLine(nsigmaPionPositive->GetMean(), resultCanvas->GetUymin(), nsigmaPionPositive->GetMean(), resultCanvas->GetUymax());
    nsigmaMean->SetLineColor(kBlue);
    nsigmaMean->SetLineStyle(kDashed);
    nsigmaMean->SetLineWidth(2);
    nsigmaMean->Draw("same");
    legend = new TLegend(0.2, 0.69, 0.4, 0.89);
    legend->SetTextSize(0.025);
    legend->SetHeader("#bf{pp, #sqrt{s} = 510 GeV}", "C");
    legend->AddEntry(nsigmaZero, "n#sigma = 0", "l");
    legend->AddEntry(nsigmaMean, "n#sigma mean", "l");
    legend->AddEntry(nsigmaPionPositive, "n#sigma of simulated tracks");
    legend->SetBorderSize(0);
    legend->DrawClone("SAME");
    resultCanvas->SaveAs("/home/adam/Pulpit/STAR_Presentation_custom/nsigmaPionPositive.pdf");
    resultCanvas->Clear();
    delete legend;
    //############################################ nsigma K+
    nsigmaKaonPositive->SetTitle("n#sigma_{K} for K^{+} tracks;n#sigma_{K};tracks");
    nsigmaKaonPositive->SetMinimum(0.);
    nsigmaKaonPositive->SetLineColor(kBlue+2);
    nsigmaKaonPositive->Draw("hist");
    nsigmaKaonPositive->SetNdivisions(520);
    resultCanvas->ModifiedUpdate(); //for Uymin and Uymax to update
    nsigmaZero = new TLine(0., resultCanvas->GetUymin(), 0., resultCanvas->GetUymax());
    nsigmaZero->SetLineColor(kRed);
    nsigmaZero->SetLineStyle(kDashed);
    nsigmaZero->SetLineWidth(2);
    nsigmaZero->Draw("same");
    nsigmaMean = new TLine(nsigmaKaonPositive->GetMean(), resultCanvas->GetUymin(), nsigmaKaonPositive->GetMean(), resultCanvas->GetUymax());
    nsigmaMean->SetLineColor(kBlue);
    nsigmaMean->SetLineStyle(kDashed);
    nsigmaMean->SetLineWidth(2);
    nsigmaMean->Draw("same");
    legend = new TLegend(0.2, 0.69, 0.4, 0.89);
    legend->SetTextSize(0.025);
    legend->SetHeader("#bf{pp, #sqrt{s} = 510 GeV}", "C");
    legend->AddEntry(nsigmaZero, "n#sigma = 0", "l");
    legend->AddEntry(nsigmaMean, "n#sigma mean", "l");
    legend->AddEntry(nsigmaKaonPositive, "n#sigma of simulated tracks");
    legend->SetBorderSize(0);
    legend->DrawClone("SAME");
    resultCanvas->SaveAs("/home/adam/Pulpit/STAR_Presentation_custom/nsigmaKaonPositive.pdf");
    resultCanvas->Clear();
    delete legend;
    //############################################ nsigma p+
    nsigmaProtonPositive->SetTitle("n#sigma_{p} for p^{+} tracks;n#sigma_{p};tracks");
    nsigmaProtonPositive->SetMinimum(0.);
    nsigmaProtonPositive->SetLineColor(kBlue+2);
    nsigmaProtonPositive->Draw("hist");
    nsigmaProtonPositive->SetNdivisions(520);
    resultCanvas->ModifiedUpdate(); //for Uymin and Uymax to update
    nsigmaZero = new TLine(0., resultCanvas->GetUymin(), 0., resultCanvas->GetUymax());
    nsigmaZero->SetLineColor(kRed);
    nsigmaZero->SetLineStyle(kDashed);
    nsigmaZero->SetLineWidth(2);
    nsigmaZero->Draw("same");
    nsigmaMean = new TLine(nsigmaProtonPositive->GetMean(), resultCanvas->GetUymin(), nsigmaProtonPositive->GetMean(), resultCanvas->GetUymax());
    nsigmaMean->SetLineColor(kBlue);
    nsigmaMean->SetLineStyle(kDashed);
    nsigmaMean->SetLineWidth(2);
    nsigmaMean->Draw("same");
    legend = new TLegend(0.2, 0.69, 0.4, 0.89);
    legend->SetTextSize(0.025);
    legend->SetHeader("#bf{pp, #sqrt{s} = 510 GeV}", "C");
    legend->AddEntry(nsigmaZero, "n#sigma = 0", "l");
    legend->AddEntry(nsigmaMean, "n#sigma mean", "l");
    legend->AddEntry(nsigmaProtonPositive, "n#sigma of simulated tracks");
    legend->SetBorderSize(0);
    legend->DrawClone("SAME");
    resultCanvas->SaveAs("/home/adam/Pulpit/STAR_Presentation_custom/nsigmaProtonPositive.pdf");
    resultCanvas->Clear();
    delete legend;
    //############################################ nsigma pi-
    nsigmaPionNegative->SetTitle("n#sigma_{#pi} for #pi^{-} tracks;n#sigma_{#pi};tracks");
    nsigmaPionNegative->SetMinimum(0.);
    nsigmaPionNegative->SetLineColor(kBlue+2);
    nsigmaPionNegative->Draw("hist");
    nsigmaPionNegative->SetNdivisions(520);
    resultCanvas->ModifiedUpdate(); //for Uymin and Uymax to update
    nsigmaZero = new TLine(0., resultCanvas->GetUymin(), 0., resultCanvas->GetUymax());
    nsigmaZero->SetLineColor(kRed);
    nsigmaZero->SetLineStyle(kDashed);
    nsigmaZero->SetLineWidth(2);
    nsigmaZero->Draw("same");
    nsigmaMean = new TLine(nsigmaPionNegative->GetMean(), resultCanvas->GetUymin(), nsigmaPionNegative->GetMean(), resultCanvas->GetUymax());
    nsigmaMean->SetLineColor(kBlue);
    nsigmaMean->SetLineStyle(kDashed);
    nsigmaMean->SetLineWidth(2);
    nsigmaMean->Draw("same");
    legend = new TLegend(0.2, 0.69, 0.4, 0.89);
    legend->SetTextSize(0.025);
    legend->SetHeader("#bf{pp, #sqrt{s} = 510 GeV}", "C");
    legend->AddEntry(nsigmaZero, "n#sigma = 0", "l");
    legend->AddEntry(nsigmaMean, "n#sigma mean", "l");
    legend->AddEntry(nsigmaPionNegative, "n#sigma of simulated tracks");
    legend->SetBorderSize(0);
    legend->DrawClone("SAME");
    resultCanvas->SaveAs("/home/adam/Pulpit/STAR_Presentation_custom/nsigmaPionNegative.pdf");
    resultCanvas->Clear();
    delete legend;
    //############################################ nsigma K-
    nsigmaKaonNegative->SetTitle("n#sigma_{K} for K^{-} tracks;n#sigma_{K};tracks");
    nsigmaKaonNegative->SetMinimum(0.);
    nsigmaKaonNegative->SetLineColor(kBlue+2);
    nsigmaKaonNegative->Draw("hist");
    nsigmaKaonNegative->SetNdivisions(520);
    resultCanvas->ModifiedUpdate(); //for Uymin and Uymax to update
    nsigmaZero = new TLine(0., resultCanvas->GetUymin(), 0., resultCanvas->GetUymax());
    nsigmaZero->SetLineColor(kRed);
    nsigmaZero->SetLineStyle(kDashed);
    nsigmaZero->SetLineWidth(2);
    nsigmaZero->Draw("same");
    nsigmaMean = new TLine(nsigmaKaonNegative->GetMean(), resultCanvas->GetUymin(), nsigmaKaonNegative->GetMean(), resultCanvas->GetUymax());
    nsigmaMean->SetLineColor(kBlue);
    nsigmaMean->SetLineStyle(kDashed);
    nsigmaMean->SetLineWidth(2);
    nsigmaMean->Draw("same");
    legend = new TLegend(0.2, 0.69, 0.4, 0.89);
    legend->SetTextSize(0.025);
    legend->SetHeader("#bf{pp, #sqrt{s} = 510 GeV}", "C");
    legend->AddEntry(nsigmaZero, "n#sigma = 0", "l");
    legend->AddEntry(nsigmaMean, "n#sigma mean", "l");
    legend->AddEntry(nsigmaKaonNegative, "n#sigma of simulated tracks");
    legend->SetBorderSize(0);
    legend->DrawClone("SAME");
    resultCanvas->SaveAs("/home/adam/Pulpit/STAR_Presentation_custom/nsigmaKaonNegative.pdf");
    resultCanvas->Clear();
    delete legend;
    //############################################ nsigma p-
    nsigmaProtonNegative->SetTitle("n#sigma_{p} for p^{-} tracks;n#sigma_{p};tracks");
    nsigmaProtonNegative->SetMinimum(0.);
    nsigmaProtonNegative->SetLineColor(kBlue+2);
    nsigmaProtonNegative->Draw("hist");
    nsigmaProtonNegative->SetNdivisions(520);
    resultCanvas->ModifiedUpdate(); //for Uymin and Uymax to update
    nsigmaZero = new TLine(0., resultCanvas->GetUymin(), 0., resultCanvas->GetUymax());
    nsigmaZero->SetLineColor(kRed);
    nsigmaZero->SetLineStyle(kDashed);
    nsigmaZero->SetLineWidth(2);
    nsigmaZero->Draw("same");
    nsigmaMean = new TLine(nsigmaProtonNegative->GetMean(), resultCanvas->GetUymin(), nsigmaProtonNegative->GetMean(), resultCanvas->GetUymax());
    nsigmaMean->SetLineColor(kBlue);
    nsigmaMean->SetLineStyle(kDashed);
    nsigmaMean->SetLineWidth(2);
    nsigmaMean->Draw("same");
    legend = new TLegend(0.2, 0.69, 0.4, 0.89);
    legend->SetTextSize(0.025);
    legend->SetHeader("#bf{pp, #sqrt{s} = 510 GeV}", "C");
    legend->AddEntry(nsigmaZero, "n#sigma = 0", "l");
    legend->AddEntry(nsigmaMean, "n#sigma mean", "l");
    legend->AddEntry(nsigmaProtonNegative, "n#sigma of simulated tracks");
    legend->SetBorderSize(0);
    legend->DrawClone("SAME");
    resultCanvas->SaveAs("/home/adam/Pulpit/STAR_Presentation_custom/nsigmaProtonNegative.pdf");
    resultCanvas->Clear();
    delete legend;

    return 0;
}

void draw_and_save(TH1D* data, TH1D* sim, std::string fileTitle, bool isYlogarithmic, double x1, double y1, double x2, double y2, double line1, double line2){
    TStyle mystyle = MyStyles::Hist2DNormalSize(false);
    mystyle.cd();
    gROOT->ForceStyle();
    TCanvas* resultCanvas = new TCanvas("resultCanvas", "resultCanvas", 1800, 1600);
    //it'll work
    resultCanvas->SetLogy(isYlogarithmic);
    if(!isYlogarithmic){
        data->SetMinimum(0);
        data->SetMaximum(std::max(data->GetMaximum(), sim->GetMaximum())*1.1);
    } else{
        data->SetMaximum(std::max(data->GetMaximum(), sim->GetMaximum())*2.1);
    }
    data->Draw("E");
    sim->Draw("Hist same");
    if(x1*y1*x2*y2!=0){
        TLegend* legend = new TLegend(x1, y1, x2, y2);
        legend->SetHeader("#bf{pp, #sqrt{s} = 510 GeV}", "C");
        legend->AddEntry(data->GetName(), "#bf{upcDST data}");
        legend->AddEntry(sim->GetName(), "#bf{PYTHIA8 simulation}");
        legend->Draw("SAME");
    }
    resultCanvas->UseCurrentStyle();
    data->SetMarkerStyle(kFullCircle);
    data->SetMarkerSize(2);
    data->SetMarkerColor(kBlue);
    data->SetLineColor(kBlue+2);
    sim->SetMarkerColor(kRed);
    sim->SetLineColor(kRed+2);
    gPad->Update();

    resultCanvas->SaveAs(("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/presentation_pdfs/"+fileTitle+".pdf").c_str());
}

void simpler_draw_and_save(TH1* hist, std::string filename, bool is3D, bool isYlogarithmic, double line1, double line2){
    TCanvas* resultCanvas = new TCanvas("resultCanvas", "resultCanvas", 1800, 1600);
    std::string fileTitle;
    gStyle->SetFrameLineWidth(1);
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);
    resultCanvas->SetLogy(isYlogarithmic);
    resultCanvas->SetLeftMargin(0.15);
    if(isYlogarithmic){
        // hist->SetMinimum(0.1);
    } else{
        hist->SetMinimum(0);
    }
    hist->SetLineColor(kBlue+2);
    if(is3D){
        resultCanvas->SetLogz();
        hist->Draw("colz");
    } else{
        hist->Draw("hist");
    }
    resultCanvas->Draw();
    if(line1!=0){
        double min = resultCanvas->GetUymin();
        double max = resultCanvas->GetUymax();
        if(isYlogarithmic){
            min = pow(10, resultCanvas->GetUymin());
            max = pow(10, resultCanvas->GetUymax());
        }
        TLine* lineDraw1 = new TLine(line1, min, line1, max);
        lineDraw1->SetLineColor(kRed);
        lineDraw1->SetLineStyle(kDashed);
        lineDraw1->SetLineWidth(2.);
        lineDraw1->Draw("same");
    }
    if(line2!=0){
        double min = resultCanvas->GetUymin();
        double max = resultCanvas->GetUymax();
        if(isYlogarithmic){
            min = pow(10, resultCanvas->GetUymin());
            max = pow(10, resultCanvas->GetUymax());
        }
        TLine* lineDraw2 = new TLine(line2, min, line2, max);
        lineDraw2->SetLineColor(kRed);
        lineDraw2->SetLineStyle(kDashed);
        lineDraw2->SetLineWidth(2.);
        lineDraw2->Draw("same");
    }
    if(filename.length()==0){
        fileTitle = hist->GetName();
    } else{
        fileTitle = filename;
    }
    resultCanvas->SaveAs(("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/presentation_pdfs/"+fileTitle+".pdf").c_str());
}

void simpler_draw_and_save(TEfficiency* hist, std::string filename){
    TCanvas* resultCanvas = new TCanvas("resultCanvas", "resultCanvas", 1800, 1600);
    std::string fileTitle;
    gStyle->SetFrameLineWidth(1);
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);
    resultCanvas->SetLeftMargin(0.15);
    resultCanvas->SetRightMargin(0.05);
    hist->SetLineColor(kBlue+2);
    hist->SetMarkerStyle(kFullCircle);
    hist->SetMarkerSize(2);
    hist->SetMarkerColor(kBlue);
    hist->Draw("ap");
    if(filename.length()==0){
        fileTitle = hist->GetName();
    } else{
        fileTitle = filename;
    }
    resultCanvas->SaveAs(("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/scripts/presentation_pdfs/"+fileTitle+".pdf").c_str());
}
