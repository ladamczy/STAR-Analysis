int FinalFit(){
    // TGFileInfo* fileinfo = new TGFileInfo();
    // fileinfo->fFileTypeIdx = 2;
    // TGFileDialog* filedialog = new TGFileDialog(gClient->GetRoot(), nullptr, EFileDialogMode::kFDOpen, fileinfo);

    // cout<<fileinfo->fFilename<<endl;
    // TFile* anaoutput = TFile::Open(fileinfo->fFilename);
    // cout<<"File opened"<<endl;

    // delete fileinfo;
    // delete filedialog;
    TFile *anaoutput1 = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/CPT/AnaOutput_Inclusive_analysis_with_STUPCV0_with_extended_range_noAfterburner.root");
    TFile *anaoutput2 = TFile::Open("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/CPTnoBBCL/AnaOutput_Inclusive_analysis_with_STUPCV0_with_extended_range_noAfterburner.root");

    gStyle->SetFrameLineWidth(3);
    gStyle->SetOptFit(111);
    gStyle->SetOptStat(0);

    TH1D *MKinvClose1 = (TH1D *)anaoutput1->Get("MpipiWide");
    TH1D *MKinvClose2 = (TH1D *)anaoutput2->Get("MpipiWide");
    TH1D resultKHist = *MKinvClose1;
    resultKHist.Add(MKinvClose2);
    resultKHist.SetStats(1);

    //K0
    TCanvas *resultK = new TCanvas("resultK", "resultK", 1800, 1600);
    TF1 *GfitK = new TF1("GfitK", "gausn(0) + pol1(3)", resultKHist.GetXaxis()->GetXmin(), resultKHist.GetXaxis()->GetXmax());
    TF1 *GfitKBcg = new TF1("GfitKBcg", "pol1", resultKHist.GetXaxis()->GetXmin(), resultKHist.GetXaxis()->GetXmax());
    TF1 *GfitK1Sig = new TF1("GfitK1Sig", "gausn", resultKHist.GetXaxis()->GetXmin(), resultKHist.GetXaxis()->GetXmax());
    GfitK->SetParNames("Constant", "Mean", "Sigma", "b", "a");
    GfitK->SetParameters(100, 0.5, 5e-3, 100, 0.48, 5e-3);
    resultKHist.SetMinimum(0);
    resultKHist.SetMarkerStyle(kFullCircle);
    resultKHist.Fit(GfitK, "0");
    resultKHist.DrawClone("E");
    Double_t paramsK[5];
    GfitK->GetParameters(paramsK);
    cout<<"PARAMETERS:"<<endl;
    for(int i = 0; i<5; i++){
        cout<<paramsK[i]<<endl;
    }
    GfitKBcg->SetParameters(paramsK+3);
    GfitK1Sig->SetParameters(paramsK);
    GfitKBcg->SetLineStyle(3);
    GfitKBcg->SetLineWidth(3);
    GfitK1Sig->SetLineColor(3);
    GfitK->SetNpx(1000);
    GfitK->Draw("CSAME");
    GfitKBcg->Draw("CSAME");
    GfitK1Sig->Draw("CSAME");

    TLegend *legendK = new TLegend(0.7, 0.7, 0.89, 0.89);
    legendK->AddEntry("MpipiWide", "pp data, #sqrt{s} = 510 GeV");
    legendK->AddEntry("GfitK", "Data fit", "l");
    legendK->AddEntry("GfitKBcg", "Background", "l");
    legendK->AddEntry("GfitK1Sig", "K^{0} fit", "l");
    legendK->SetBorderSize(0);
    legendK->Draw("SAME");


    // gPad->RedrawAxis();
    resultK->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/CPT/K0afterV0finder.pdf").c_str());

    return 0;
}