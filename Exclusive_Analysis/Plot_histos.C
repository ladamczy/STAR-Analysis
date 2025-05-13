#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TKey.h>
#include <TList.h>
#include <iostream>
#include <map>
#include <string>
#include <vector>

/*
void DrawSelectedHistograms(const char* filename) {
    // Open the ROOT file
    TFile* file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    // Create output PDF file
    TCanvas* c = new TCanvas("c", "c", 1800, 1300);
    c->Print("plots/histograms.pdf[");  // Open PDF file

    // Set style options
    gStyle->SetStatY(0.990);
    gStyle->SetStatX(0.990);
    gStyle->SetStatW(0.15);
    gStyle->SetStatH(0.18);
    gStyle->SetPalette(kBird);
    //gStyle->SetPadLeftMargin(0.4);
    //gStyle->SetPadRightMargin(0.4);
    //gStyle->SetPadTopMargin(0.1);
    //gStyle->SetPadBottomMargin(0.99);
    gStyle->SetOptStat(0);

    // Set pad margins directly on the canvas pad
    c->SetLeftMargin(0.15);   // Adjust left margin (default is 0.1)
    c->SetRightMargin(0.10);  // Adjust right margin (default is 0.1)
    c->SetTopMargin(0.10);    // Adjust top margin (default is 0.1)
    c->SetBottomMargin(0.15); // Adjust bottom margin (default is 0.1)

    

    // List of histograms to draw
    std::vector<std::string> histogramsToDraw = {
        "HistMassK0K0Both",
        "HistMassK0K0Opp",
        "HistMassK0K0Same",
        "HistEtaK0K0",
        "HistPtK0K0Same",
        "HistPtK0K0Opp",
        "HistPtK0K0Both",
        //"HistPxK0K0",
        //"HistPyK0K0"
        "HistKsiEK0K0",
        "HistKsiWK0K0",
        "HistMassK0",
        "HistPtK0",
        "HistEtaK0",
        "HistPhiK0"
    };

    // Loop through the selected histograms
    for (const auto& histName : histogramsToDraw) {
        TObject* obj = file->Get(histName.c_str());
        if (!obj) {
            std::cerr << "Warning: Histogram " << histName << " not found in the file." << std::endl;
            continue;
        }
        
        // Clone the object to prevent it from being deleted when getting the next one
        obj = obj->Clone();

        // Check if the object is a 1D histogram
        if (obj->InheritsFrom(TH1::Class()) && !obj->InheritsFrom(TH2::Class())) {
            TH1* h1 = (TH1*)obj;
            
            // Apply rebinning to specific histograms
            if (histName == "HistPtK0K0Both" || histName == "HistEtaK0K0") {
                std::cout << "Rebinning histogram: " << histName << " by factor of 2" << std::endl;
                h1->Rebin(2);
            }
            // Set specific axis ranges based on histogram name
            if (histName == "HistMassK0K0Both") {
                h1->GetXaxis()->SetRangeUser(0.9, 3.0); // X
                h1->GetYaxis()->SetRangeUser(0.0, 25.0); // Y
            }else if (histName == "HistMassK0K0Same") {
                h1->GetXaxis()->SetRangeUser(0.9, 3.0); // X
                h1->GetYaxis()->SetRangeUser(0.0, 22.0); // Y
            }else if (histName == "HistMassK0K0Opp") {
                h1->GetXaxis()->SetRangeUser(0.9, 3.0); // X
                h1->GetYaxis()->SetRangeUser(0.0, 13.0); // Y
            } else if (histName == "HistEtaK0K0") {
                h1->GetXaxis()->SetRangeUser(-5.0, 5.0); //X
                h1->GetYaxis()->SetRangeUser(0.0, 50.0); // Y
            } else if (histName == "HistPtK0K0Both") {
                h1->GetXaxis()->SetRangeUser(0, 2.0); // X
                h1->GetYaxis()->SetRangeUser(0.0, 30.0); // Y
            }
           //else if (histName == "HistPtK0K0Same" || histName == "HistPtK0K0Opp" || histName == "HistPtK0K0Both") {
           //  h1->GetXaxis()->SetRangeUser(0, 5.0); // X
           // } 

            c->cd();
            c->Clear();
            h1->GetXaxis()->SetTitleOffset(0.9);
            h1->GetYaxis()->SetTitleOffset(0.8);
            h1->GetXaxis()->SetTitleSize(0.06);
            h1->GetYaxis()->SetTitleSize(0.06);
            h1->GetXaxis()->SetLabelSize(0.05);
            h1->GetYaxis()->SetLabelSize(0.06);
            h1->GetXaxis()->SetNdivisions(6);
            h1->GetYaxis()->SetNdivisions(6);
            h1->SetMarkerStyle(4);
            h1->SetMarkerSize(0.9);
            h1->SetMarkerColor(kBlack);
            h1->SetLineColor(kBlack);
            h1->Draw("PE");
            c->Update();
            c->Print("plots/histograms.pdf");  // Add this histogram to PDF            
        }
        // Check if the object is a 2D histogram
        else if (obj->InheritsFrom(TH2::Class())) {
            TH2* h2 = (TH2*)obj;
            c->cd();
            c->Clear();
            h2->GetXaxis()->SetTitleOffset(0.6);
            h2->GetYaxis()->SetTitleOffset(0.8);
            h2->GetXaxis()->SetTitleSize(0.06);
            h2->GetYaxis()->SetTitleSize(0.05);
            h2->GetXaxis()->SetLabelSize(0.05);
            h2->GetYaxis()->SetLabelSize(0.06);
            h2->GetXaxis()->SetNdivisions(6);
            h2->GetYaxis()->SetNdivisions(6);
            h2->Draw("colz");
            c->Update();
            c->Print("plots/histograms.pdf");  // Add this histogram to PDF
        }
    }

    // Close PDF file
    c->Print("plots/histograms.pdf]");

    // Clean up
    delete c;
    file->Close();
    
    std::cout << "Successfully created PDF with " << histogramsToDraw.size() << " selected histograms." << std::endl;
}

void Plot_histos() {
    const char* filename = "InputData/ex_try.root";//"build/ex3.root";
    DrawSelectedHistograms(filename);
}
*/
void DrawAllHistograms(const char* filename) {
    // Open the ROOT file
    TFile* file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    // Create output PDF file
    TCanvas* c = new TCanvas("c", "c", 1800, 1300);
    c->Print("plots/histograms_tof_newcuts.pdf[");  // Open PDF file
    //TPad *pad1 = (TPad*)c->GetPad(1);
    // Set margins (left, right, bottom, top)
    // Set pad margins directly on the canvas pad
    c->SetLeftMargin(0.15);   // Adjust left margin (default is 0.1)
    c->SetRightMargin(0.10);  // Adjust right margin (default is 0.1)
    c->SetTopMargin(0.10);    // Adjust top margin (default is 0.1)
    c->SetBottomMargin(0.15); // Adjust bottom margin (default is 0.1)


    gStyle->SetStatY(0.990);
    gStyle->SetStatX(0.990);
    gStyle->SetStatW(0.15);
    gStyle->SetStatH(0.18);
    gStyle->SetPalette(kBird);
    //gStyle->SetPadLeftMargin(0.4);
    //gStyle->SetPadRightMargin(0.4);
    //gStyle->SetPadTopMargin(0.1);
    //gStyle->SetPadBottomMargin(0.99);
    //gStyle->SetOptStat(0);
      //gStyle->SetPadBorderSize(2);

    // Get the list of keys in the file
    TList* keys = file->GetListOfKeys();
    if (!keys) {
        std::cerr << "Error: No keys found in the file." << std::endl;
        file->Close();
        return;
    }

    // Loop over all keys in the file
    TIter next(keys);
    TKey* key;
    while ((key = (TKey*)next())) {
        // Get the object associated with the key
        TObject* obj = key->ReadObj();
        if (!obj) continue;

        // Check if the object is a 1D histogram
        if (obj->InheritsFrom(TH1::Class()) && !obj->InheritsFrom(TH2::Class())) {
            TH1* h1 = (TH1*)obj;
            c->cd();
            c->Clear();
            h1->GetXaxis()->SetTitleOffset(0.9);
            h1->GetYaxis()->SetTitleOffset(0.8);
            h1->GetXaxis()->SetTitleSize(0.06);
            h1->GetYaxis()->SetTitleSize(0.06);
            h1->GetXaxis()->SetLabelSize(0.05);
            h1->GetYaxis()->SetLabelSize(0.06);
            h1->GetXaxis()->SetNdivisions(6);
            h1->GetYaxis()->SetNdivisions(6);
            h1->SetMarkerStyle(4);
            h1->SetMarkerSize(0.9);
            h1->SetMarkerColor(kBlack);
            h1->SetLineColor(kBlack);
            //gPad->SetLogy(1);
            //h1->Draw();
            //h1->GetYaxis()->SetRangeUser(-10.0, 100.0); // Set Y-axis range
            h1->Draw("PE");
            c->Update();
            c->Print("plots/histograms_tof_newcuts.pdf");  // Add this histogram to PDF
        }
        // Check if the object is a 2D histogram
        else if (obj->InheritsFrom(TH2::Class())) {
            TH2* h2 = (TH2*)obj;
            c->cd();
            c->Clear();
            h2->GetXaxis()->SetTitleOffset(0.6);
            h2->GetYaxis()->SetTitleOffset(0.8);
            h2->GetXaxis()->SetTitleSize(0.06);
            h2->GetYaxis()->SetTitleSize(0.05);
            h2->GetXaxis()->SetLabelSize(0.05);
            h2->GetYaxis()->SetLabelSize(0.06);
            h2->GetXaxis()->SetNdivisions(6);
            h2->GetYaxis()->SetNdivisions(6);
            h2->Draw("colz");
            c->Update();
            c->Print("plots/histograms_tof_newcuts.pdf");  // Add this histogram to PDF
        }
    }

    // Close PDF file
    c->Print("plots/histograms_tof_newcuts.pdf]");

    // Clean up
    delete c;
    file->Close();

}

void Plot_histos() {
    const char* filename = "InputData/tofeff_withRP_0.3_v3.root"; //InputData/tofeff.root";  tofeffsim_star  tofeffsim_withRP
    DrawAllHistograms(filename);

}
