int Hist(){
    TFile *anaoutput = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/CPT/AnaOutput_Inclusive_analysis_with_STUPCV0_with_extended_range_noAfterburner.root");
    TFile *anaoutput2 = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/CPTnoBBCL/AnaOutput_Inclusive_analysis_with_STUPCV0_with_extended_range_noAfterburner.root");

    gStyle->SetHistLineWidth(3);
    gStyle->SetFrameLineWidth(3);
    gROOT->ForceStyle();

    //NotherTracks
    TCanvas *c1 = new TCanvas("c1", "c1", 1800, 1600);
    c1->SetGrid();
    gPad->SetLogy();
    TH1D *NotherTracks = (TH1D *)anaoutput->Get("NotherTracks");
    TH1D *NotherTracks2 = (TH1D *)anaoutput2->Get("NotherTracks");
    NotherTracks2->SetName("NotherTracks2");
    NotherTracks->Scale(1/NotherTracks->Integral());
    NotherTracks2->Scale(1/NotherTracks2->Integral());
    NotherTracks->SetLineWidth(2);
    NotherTracks2->SetLineWidth(2);
    NotherTracks->SetLineColor(kBlue);
    NotherTracks2->SetLineColor(kRed);
    NotherTracks->SetStats(0);
    NotherTracks2->SetStats(0);
    NotherTracks->Draw("HIST");
    NotherTracks2->Draw("HIST SAME");
    TLegend *legendK = new TLegend(0.7, 0.7, 0.89, 0.89);
    legendK->AddEntry("NotherTracks", "CTP2", "l");
    legendK->AddEntry("NotherTracks2", "CTP2noBBCL", "l");
    legendK->SetBorderSize(1);
    legendK->Draw("SAME");
    c1->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/CPT/NotherTracksAll.pdf").c_str());
    delete c1;

    //VertexMultiplicity
    c1 = new TCanvas("c1", "c1", 1800, 1600);
    TH2D *VertexMultiplicity = (TH2D *)anaoutput->Get("VertexMultiplicity");
    c1->SetLogz();
    c1->SetMargin(0.15, 0.15, 0.1, 0.1);
    VertexMultiplicity->SetStats(0);
    VertexMultiplicity->Draw("COLZ");
    gPad->RedrawAxis();
    c1->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/CPT/VertexMultiplicity.pdf").c_str());
    delete c1;

    //K0multiplicity
    c1 = new TCanvas("c1", "c1", 1800, 1600);
    c1->SetGrid();
    gPad->SetLogy();
    TH1D *K0multiplicity = (TH1D *)anaoutput->Get("K0multiplicity");
    TH1D *K0multiplicity2 = (TH1D *)anaoutput2->Get("K0multiplicity");
    K0multiplicity2->SetName("K0multiplicity2");
    K0multiplicity->Scale(1/K0multiplicity->Integral());
    K0multiplicity2->Scale(1/K0multiplicity2->Integral());
    K0multiplicity->SetLineWidth(2);
    K0multiplicity2->SetLineWidth(2);
    K0multiplicity->SetLineColor(kBlue);
    K0multiplicity2->SetLineColor(kRed);
    K0multiplicity->SetStats(0);
    K0multiplicity2->SetStats(0);
    K0multiplicity->Draw("HIST");
    K0multiplicity2->Draw("HIST SAME");
    TLegend *legendK0 = new TLegend(0.7, 0.7, 0.89, 0.89);
    legendK0->AddEntry("K0multiplicity", "CTP2", "l");
    legendK0->AddEntry("K0multiplicity2", "CTP2noBBCL", "l");
    legendK0->SetBorderSize(1);
    legendK0->Draw("SAME");
    c1->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/CPT/K0multiplicityAll.pdf").c_str());
    delete c1;

    //PVV0K0dist
    c1 = new TCanvas("c1", "c1", 1800, 1600);
    c1->SetGrid();
    // gPad->SetLogy();
    TH1D *PVV0K0dist = (TH1D *)anaoutput->Get("PVV0K0dist");
    PVV0K0dist->SetLineWidth(2);
    // PVV0K0dist->GetXaxis()->SetTitle("Number of primary vertices");
    PVV0K0dist->SetStats(0);
    PVV0K0dist->Draw();
    c1->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/CPT/PVV0K0dist.pdf").c_str());
    delete c1;

    //K0PVdistance
    c1 = new TCanvas("c1", "c1", 1800, 1600);
    c1->SetGrid();
    // gPad->SetLogy();
    TH1D *K0PVdistance = (TH1D *)anaoutput->Get("K0PVdistance");
    K0PVdistance->SetLineWidth(2);
    // K0PVdistance->GetXaxis()->SetTitle("Number of primary vertices");
    K0PVdistance->SetStats(0);
    K0PVdistance->Draw();
    c1->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/CPT/K0PVdistance.pdf").c_str());
    delete c1;



    return 0;
}
