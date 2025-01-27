// Stdlib header file for input and output.
#include <iostream>
#include <cstring>

// ROOT, for histogramming.
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TEfficiency.h"
#include "THStack.h"
#include "TCanvas.h"

void drawStack(TH2D* particles, TH1D* background, std::string folderWithDiagonal, std::string name, std::string title, double x1, double x2);

int main(int argc, char const* argv[]){
    TFile* input = TFile::Open(argv[1]);

    //getting histograms out
    TH1D* MpipiNonresonant = (TH1D*)input->Get("MpipiNonresonant");
    TH2D* MpipiParticles = (TH2D*)input->Get("MpipiParticles");
    TH1D* MKpiNonresonant = (TH1D*)input->Get("MKpiNonresonant");
    TH2D* MKpiParticles = (TH2D*)input->Get("MKpiParticles");
    TH1D* MKKNonresonant = (TH1D*)input->Get("MKKNonresonant");
    TH2D* MKKParticles = (TH2D*)input->Get("MKKParticles");

    std::string folderWithDiagonal = std::string(argv[2]);
    if(folderWithDiagonal[folderWithDiagonal.size()-1]!='/'){
        folderWithDiagonal += "/";
    }

    //drawing them
    drawStack(MpipiParticles, nullptr, folderWithDiagonal, "MpipiNoBcg", "Assuming  #pi#pi masses, only resonant production shown", 0.2, 1.0);
    drawStack(MpipiParticles, MpipiNonresonant, folderWithDiagonal, "Mpipi", "Assuming  #pi#pi masses", 0.2, 1.0);
    drawStack(MKpiParticles, nullptr, folderWithDiagonal, "MKpiNoBcg", "Assuming  K#pi masses, only resonant production shown", 0.5, 1.3);
    drawStack(MKpiParticles, MKpiNonresonant, folderWithDiagonal, "MKpi", "Assuming  K#pi masses", 0.5, 1.3);
    drawStack(MKKParticles, nullptr, folderWithDiagonal, "MKKNoBcg", "Assuming  KK masses, only resonant production shown", 0.9, 1.7);
    drawStack(MKKParticles, MKKNonresonant, folderWithDiagonal, "MKK", "Assuming  KK masses", 0.9, 1.7);

    return 0;
}

void drawStack(TH2D* particles, TH1D* background, std::string folderWithDiagonal, std::string name, std::string title, double x1, double x2){
    THStack MStack(name.c_str(), title.c_str());
    std::vector<TH1D*> projections;
    for(Int_t i = 1; i<particles->GetXaxis()->GetNbins()+1; i++){
        projections.push_back(new TH1D(*particles->ProjectionY("", i, i)));
        projections[projections.size()-1]->SetName(particles->GetXaxis()->GetBinLabel(i));
    }
    //adding background if necessary
    if(background!=nullptr){
        background->SetName("Nonresonant or mismatched");
        background->SetFillColor(1);
        background->SetLineWidth(0);
        MStack.Add(background);
    }
    //sorting in terms of input to color only the x highest
    sort(projections.begin(), projections.end(), [](TH1D* x, TH1D* y){return x->GetEntries()>y->GetEntries();});
    int n_to_colour = 5;
    double opacity = 0.5;
    //setting apropriate colours and opacity, including small-contribution rest of particles
    for(Int_t i = 0; i<projections.size(); i++){
        if(i<n_to_colour){
            projections[i]->SetFillColorAlpha(i+2+(i+2>=10), 1.-(background==nullptr)*opacity); //color 10 is white, cant have that
            projections[i]->SetLineWidth(0);
            MStack.Add(projections[i]);
        } else if(i==n_to_colour){
            projections[i]->SetName("other");
            projections[i]->SetFillColorAlpha(kGray+2, 1.-(background==nullptr)*opacity);
            projections[i]->SetLineWidth(0);
        } else{
            projections[n_to_colour]->Add(projections[i]);
        }
        //adding the "other" projection last
        if(i==projections.size()-1){
            MStack.Add(projections[n_to_colour]);
        }
    }
    //actual drawing
    TCanvas canvas("Stack", "Stack", 4000, 2400);
    std::string stackOptions = "";
    if(background==nullptr){
        stackOptions = "nostack";
    }
    MStack.Draw(stackOptions.c_str());
    MStack.GetXaxis()->SetRangeUser(x1, x2);
    MStack.GetXaxis()->SetTitle("m_{inv} [GeV/c^{2}]");
    MStack.GetYaxis()->SetTitle("number of particles");
    canvas.Update();
    // canvas.BuildLegend(0.7, 0.6, 0.89, 0.89, "Particles reconstructed", "F");
    canvas.BuildLegend(0.7, 0.7, 0.89, 0.89, "Particles reconstructed", "F");
    canvas.SaveAs((folderWithDiagonal+name+".pdf").c_str());
    for(size_t i = 0; i<projections.size(); i++){
        delete projections[i];
    }
    canvas.Clear();
}