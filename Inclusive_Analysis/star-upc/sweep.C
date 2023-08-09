int sweep(){
    // TGFileInfo* fileinfo = new TGFileInfo();
    // fileinfo->fFileTypeIdx = 2;
    // TGFileDialog* filedialog = new TGFileDialog(gClient->GetRoot(), nullptr, EFileDialogMode::kFDOpen, fileinfo);

    // cout<<fileinfo->fFilename<<endl;
    // TFile* anaoutput = TFile::Open(fileinfo->fFilename);
    // cout<<"File opened"<<endl;

    // delete fileinfo;
    // delete filedialog;

    TFile* anaoutput = TFile::Open("AnaOutput_DCASweep.root");

    TH2D* Minv = (TH2D*)anaoutput->Get("MinvDCAsweepHist");
    cout<<Minv<<endl;

    TH1D* temp_hist = new TH1D();
    TF1* Gfit = new TF1("Gfit", "gausn(0) + gausn(3) + pol1(6)", Minv->GetXaxis()->GetXmin(), Minv->GetXaxis()->GetXmax());
    Gfit->SetParNames("Constant", "Mean", "Sigma", "Constant2", "Mean2", "Sigma2", "a", "b");
    Gfit->SetParameters(100, 0.5, 5e-3, 100, 0.48, 5e-3);
    Gfit->SetParLimits(0, 0, 100);
    Gfit->SetParLimits(1, 0.49, 0.51);
    Gfit->SetParLimits(3, 0, 100);
    Gfit->SetParLimits(4, 0.47, 0.49);
    double x[50] = {0};
    double y[50] = {0};
    Double_t paramsnobcg[8];

    for (int i = 0; i < 50; i++){
        cout << i << endl;
        temp_hist = Minv->ProjectionX("temp_hist", i, i+1);
        temp_hist->Fit(Gfit, "Q0");
        Gfit->GetParameters(paramsnobcg);
        cout<<"PARAMETERS:"<<endl;
        for (int i = 0; i < 8; i++){
            cout<<paramsnobcg[i]<<endl;
        }
        x[i] = i*0.2;
        y[i] = paramsnobcg[3];
    }

    TCanvas *resultnobcgCanvas = new TCanvas("resultnobcgCanvas","resultnobcgCanvas",1800,1600);
    TGraph* func1 = new TGraph(50, x, y);
    func1->Draw();
    

    return 0;
}