int ROOTtoPDF(){
    //otwieranie nowego pliku
    TGFileInfo* fileinfo = new TGFileInfo();
    fileinfo->fFileTypeIdx = 2;
    TGFileDialog* filedialog = new TGFileDialog(gClient->GetRoot(), nullptr, EFileDialogMode::kFDOpen, fileinfo);
    cout<<fileinfo->fFilename<<endl;
    TFile* anaoutput = TFile::Open(fileinfo->fFilename);
    //zapisuje nazwę bez kropki
    string filename = string(fileinfo->fFilename).substr(0, string(fileinfo->fFilename).find("."));
    delete fileinfo;
    delete filedialog;

    //ustalanie stylu wykresów
    gStyle->SetHistLineWidth(1);
    gStyle->SetFrameLineWidth(3);
    gStyle->SetOptStat(111111);
    gStyle->SetStatY(0.89);
    gStyle->SetStatX(0.89);
    gStyle->SetStatW(0.1);
    gStyle->SetStatH(0.1); 
    gROOT->ForceStyle();

    //zmienne do rysowania
    //ustawianie batch mode żeby nie otwierały się okna podczas rysowania
    gROOT->SetBatch(kTRUE);
    TCanvas* c1 = new TCanvas("c1","c1",1600,1200);
    TPaveStats *st;
    TH1D* temp1D;
    TH2D* temp2D;
    string name;
    c1->Print((filename + ".pdf[").c_str());
    

    //przetwarzanie wykresów, oryginał kodu z 
    //https://root-forum.cern.ch/t/finding-all-objects-in-a-root-file-regardless-of-directries/2351
    //przetwarza również podfoldery w pliku .root

    TDirectory *dir = anaoutput;
    TIter next(dir->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)next())) {
        if (strstr(key->GetClassName(),"TH1")){
            printf (" key : %s is a %s in %s\n", key->GetName(), key->GetClassName(), dir->GetPath());
            name = key->GetName();
            temp1D = (TH1D*)anaoutput->Get(name.c_str());
            temp1D->SetMarkerStyle(kFullCircle);
            temp1D->Draw("E");
            gPad->RedrawAxis();
            c1->Print((filename + ".pdf").c_str());
        }
        if (strstr(key->GetClassName(),"TH2")){
            printf (" key : %s is a %s in %s\n", key->GetName(), key->GetClassName(), dir->GetPath());
            name = key->GetName();
            temp2D = (TH2D*)anaoutput->Get(name.c_str());
            temp2D->SetMarkerStyle(kFullCircle);
            temp2D->Draw("COLZ");
            gPad->RedrawAxis();
            c1->Print((filename + ".pdf").c_str());
        }
        //TODO: loopdir wymaga rekurencji, zajmę się tym później
        //strcmp zwraca 0 gdy stringi są takie same
        // if (!strcmp(key->GetClassName(),"TDirectory")){
        //     dir->cd(key->GetName());
        //     TDirectory *subdir = gDirectory;
        //     loopdir(subdir);
        //     dir->cd();
        // }
    }
    c1->Print((filename + ".pdf]").c_str());

    return 0;
}