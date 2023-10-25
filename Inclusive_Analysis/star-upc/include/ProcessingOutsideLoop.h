#ifndef PROCESSING_OUTSIDE_LOOP_H
#define PROCESSING_OUTSIDE_LOOP_H

#include "TClass.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TObject.h"
#include <stdio.h>
#include <string>
#include <ROOT/TThreadedObject.hxx>
#include "TFile.h"

//to avoid circular reference we "dummy define" class here i guess? like with functions
class ProcessingInsideLoop;

class ProcessingOutsideLoop {
private:
    //so that Inside loop instances would have access to these things
    friend ProcessingInsideLoop;
    std::vector<ROOT::TThreadedObject<TH1D> *> hist1dtab;
    std::vector<ROOT::TThreadedObject<TH2D> *> hist2dtab;
    std::vector<ROOT::TThreadedObject<TH3D> *> hist3dtab;
    std::vector<std::shared_ptr<TH1D>> hist1dtabFinal;
    std::vector<std::shared_ptr<TH2D>> hist2dtabFinal;
    std::vector<std::shared_ptr<TH3D>> hist3dtabFinal;

public:
    ProcessingOutsideLoop(/* args */);
    ~ProcessingOutsideLoop();
    void AddHistogram(TH1D);
    void AddHistogram(TH2D);
    void AddHistogram(TH3D);
    void Merge();
    void SaveToFile(TFile *);
};

ProcessingOutsideLoop::ProcessingOutsideLoop(/* args */) {
}

ProcessingOutsideLoop::~ProcessingOutsideLoop() {
}

void ProcessingOutsideLoop::AddHistogram(TH1D hist){
    hist1dtab.push_back(new ROOT::TThreadedObject<TH1D>(hist));
    hist2dtab.push_back(nullptr);
    hist3dtab.push_back(nullptr);
}

void ProcessingOutsideLoop::AddHistogram(TH2D hist){
    hist1dtab.push_back(nullptr);
    hist2dtab.push_back(new ROOT::TThreadedObject<TH2D>(hist));
    hist3dtab.push_back(nullptr);
}

void ProcessingOutsideLoop::AddHistogram(TH3D hist){
    hist1dtab.push_back(nullptr);
    hist2dtab.push_back(nullptr);
    hist3dtab.push_back(new ROOT::TThreadedObject<TH3D>(hist));
}

void ProcessingOutsideLoop::Merge(){
    for(int i = 0; i<hist1dtab.size(); i++){
        hist1dtabFinal.push_back(nullptr);
        hist2dtabFinal.push_back(nullptr);
        hist3dtabFinal.push_back(nullptr);
        if(hist1dtab[i]!=nullptr){
            hist1dtabFinal[i] = hist1dtab[i]->Merge();
        } else if(hist2dtab[i]!=nullptr){
            hist2dtabFinal[i] = hist2dtab[i]->Merge();
        } else if(hist3dtab[i]!=nullptr){
            hist3dtabFinal[i] = hist3dtab[i]->Merge();
        }
    }
}

void ProcessingOutsideLoop::SaveToFile(TFile *file){
    file->cd();
    for(int i = 0; i<hist1dtabFinal.size(); i++){
        if(hist1dtabFinal[i]!=nullptr){
            hist1dtabFinal[i]->Write();
        } else if(hist2dtabFinal[i]!=nullptr){
            hist2dtabFinal[i]->Write();
        } else if(hist3dtabFinal[i]!=nullptr){
            hist3dtabFinal[i]->Write();
        }
    }
}

#endif