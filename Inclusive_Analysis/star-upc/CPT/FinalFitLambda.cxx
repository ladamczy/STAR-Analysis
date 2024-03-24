int FinalFitLambda(){
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

    TH1D *MLambdainvClose1 = (TH1D *)anaoutput1->Get("MppiWide");
    TH1D *MLambdainvClose2 = (TH1D *)anaoutput2->Get("MppiWide");
    TH1D resultLambdaHist = *MLambdainvClose1;
    resultLambdaHist.Add(MLambdainvClose2);
    resultLambdaHist.SetStats(1);

    // double fitmin = resultLambdaHist.GetXaxis()->GetXmin();
    // double fitmax = resultLambdaHist.GetXaxis()->GetXmax();
    double fitmin = 1.1;
    double fitmax = 1.14;

    //Lambda
    TCanvas *resultLambda = new TCanvas("resultLambda", "resultLambda", 1800, 1600);
    TF1 *GfitLambda = new TF1("GfitLambda", "gausn(0) + pol2(3)", fitmin, fitmax);
    TF1 *GfitLambdaBcg = new TF1("GfitLambdaBcg", "pol2", fitmin, fitmax);
    TF1 *GfitLambda1Sig = new TF1("GfitLambda1Sig", "gausn", fitmin, fitmax);
    GfitLambda->SetParNames("Constant", "Mean", "Sigma", "c", "b", "a");
    GfitLambda->SetParameters(100, 1.115, 5e-3, 2376400, -4360000, 3000000);
    GfitLambda->SetParLimits(0, 0, 200);
    GfitLambda->SetParLimits(1, 1.11, 1.12);
    GfitLambda->SetParLimits(2, 0, 0.01);
    resultLambdaHist.SetMinimum(0);
    resultLambdaHist.SetMarkerStyle(kFullCircle);
    resultLambdaHist.Fit(GfitLambda, "0BR");
    resultLambdaHist.DrawClone("E");
    Double_t paramsK[6];
    GfitLambda->GetParameters(paramsK);
    cout<<"PARAMETERS:"<<endl;
    for(int i = 0; i<6; i++){
        cout<<paramsK[i]<<endl;
    }
    GfitLambdaBcg->SetParameters(paramsK+3);
    GfitLambda1Sig->SetParameters(paramsK);
    GfitLambdaBcg->SetLineStyle(3);
    GfitLambdaBcg->SetLineWidth(3);
    GfitLambda1Sig->SetLineColor(3);
    GfitLambda->SetNpx(1000);
    GfitLambda1Sig->SetNpx(1000);
    GfitLambda->Draw("CSAME");
    GfitLambdaBcg->Draw("CSAME");
    GfitLambda1Sig->Draw("CSAME");

    TLegend *legendLambda = new TLegend(0.7, 0.7, 0.89, 0.89);
    legendLambda->AddEntry("MpipiWide", "pp data, #sqrt{s} = 510 GeV");
    legendLambda->AddEntry("GfitLambda", "Data fit", "l");
    legendLambda->AddEntry("GfitLambdaBcg", "Background", "l");
    legendLambda->AddEntry("GfitLambda1Sig", "#Lambda^{0} fit", "l");
    legendLambda->SetBorderSize(0);
    legendLambda->Draw("SAME");


    // gPad->RedrawAxis();
    resultLambda->SaveAs(string("/home/adam/STAR-Analysis/Inclusive_Analysis/star-upc/CPT/Lambda0afterV0finder.pdf").c_str());

    return 0;
}