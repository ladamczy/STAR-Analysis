#ifndef PROCESSING_INSIDE_LOOP_H
#define PROCESSING_INSIDE_LOOP_H

#include "ProcessingOutsideLoop.h"

class ProcessingInsideLoop {
private:
    std::vector<std::shared_ptr<TH1D>> hist1dtabLocal;
    std::vector<std::shared_ptr<TH2D>> hist2dtabLocal;
    std::vector<std::shared_ptr<TH3D>> hist3dtabLocal;
public:
    ProcessingInsideLoop(/* args */);
    ~ProcessingInsideLoop();
    void GetLocalHistograms(ProcessingOutsideLoop *);
    void Fill(int, double);
    void Fill(int, double, double);
};

ProcessingInsideLoop::ProcessingInsideLoop(/* args */) {
}

ProcessingInsideLoop::~ProcessingInsideLoop() {
}

void ProcessingInsideLoop::GetLocalHistograms(ProcessingOutsideLoop *outsideloop){
    for(long unsigned int i = 0; i<outsideloop->hist1dtab.size(); i++){
        hist1dtabLocal.push_back(nullptr);
        hist2dtabLocal.push_back(nullptr);
        hist3dtabLocal.push_back(nullptr);
        if(outsideloop->hist1dtab[i]!=nullptr){
            hist1dtabLocal[i] = outsideloop->hist1dtab[i]->Get();
        } else if(outsideloop->hist2dtab[i]!=nullptr){
            hist2dtabLocal[i] = outsideloop->hist2dtab[i]->Get();
        } else if(outsideloop->hist3dtab[i]!=nullptr){
            hist3dtabLocal[i] = outsideloop->hist3dtab[i]->Get();
        }
    }

}

void ProcessingInsideLoop::Fill(int hist_number, double x){
    hist1dtabLocal[hist_number]->Fill(x);
}

void ProcessingInsideLoop::Fill(int hist_number, double x, double y_or_w){
    if(hist1dtabLocal[hist_number]!=nullptr){
        hist1dtabLocal[hist_number]->Fill(x, y_or_w);
    } else{
        hist2dtabLocal[hist_number]->Fill(x, y_or_w);
    }

}

#endif