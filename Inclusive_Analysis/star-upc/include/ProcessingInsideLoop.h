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
    void Fill(const char *, double);
    void Fill(int, double, double);
    void Fill(const char *, double, double);
    void Fill(int, const char*, double);
    void Fill(const char*, const char*, double);
    void Fill(int, const char*, double, double);
    void Fill(const char*, const char*, double, double);
    void Fill(int, const char*, const char*, double);
    void Fill(const char*, const char*, const char*, double);
    void Fill(int, const char*, double, double, double);
    void Fill(const char*, const char*, double, double, double);
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

//1d histogram, no weight
void ProcessingInsideLoop::Fill(int hist_number, double x){
    if(hist1dtabLocal[hist_number]!=nullptr){
        hist1dtabLocal[hist_number]->Fill(x);
    } else{
        throw std::invalid_argument("1D histogram with number \""+std::to_string(hist_number)+"\" does not exist.");
    }
}
void ProcessingInsideLoop::Fill(const char *hist_name, double x){
    for(long unsigned int i = 0; i<hist1dtabLocal.size(); i++){
        if(hist1dtabLocal[i]!=nullptr&&strcmp(hist_name, hist1dtabLocal[i]->GetName())==0){
            ProcessingInsideLoop::Fill(i, x);
            return;
        }
    }
    throw std::invalid_argument("Histogram with name \""+std::string(hist_name)+"\" could not be found.");
}

//1d histogram, weight, or 2d histogram, no weight
void ProcessingInsideLoop::Fill(int hist_number, double x, double y_or_w){
    if(hist1dtabLocal[hist_number]!=nullptr){
        hist1dtabLocal[hist_number]->Fill(x, y_or_w);
    } else if(hist2dtabLocal[hist_number]!=nullptr){
        hist2dtabLocal[hist_number]->Fill(x, y_or_w);
    } else{
        throw std::invalid_argument("1D nor 2D histogram with number \""+std::to_string(hist_number)+"\" does not exist.");
    }
}
void ProcessingInsideLoop::Fill(const char *hist_name, double x, double y_or_w){
    //looking through 1d histograms
    for(long unsigned int i = 0; i<hist1dtabLocal.size(); i++){
        if(hist1dtabLocal[i]!=nullptr&&strcmp(hist_name, hist1dtabLocal[i]->GetName())==0){
            ProcessingInsideLoop::Fill(i, x, y_or_w);
            return;
        }
    }
    //looking through 2d histograms
    for(long unsigned int i = 0; i<hist2dtabLocal.size(); i++){
        if(hist2dtabLocal[i]!=nullptr&&strcmp(hist_name, hist2dtabLocal[i]->GetName())==0){
            ProcessingInsideLoop::Fill(i, x, y_or_w);
            return;
        }
    }
    throw std::invalid_argument("Histogram with name \""+std::string(hist_name)+"\" could not be found.");
}

//1d histogram, alphanumeric, with weight
void ProcessingInsideLoop::Fill(int hist_number, const char* key, double w){
    if(hist1dtabLocal[hist_number]!=nullptr){
        hist1dtabLocal[hist_number]->Fill(key, w);
    } else{
        throw std::invalid_argument("1D histogram with number \""+std::to_string(hist_number)+"\" does not exist.");
    }
}
void ProcessingInsideLoop::Fill(const char* hist_name, const char* key, double w){
    //looking through 1d histograms
    for(long unsigned int i = 0; i<hist1dtabLocal.size(); i++){
        if(hist1dtabLocal[i]!=nullptr&&strcmp(hist_name, hist1dtabLocal[i]->GetName())==0){
            ProcessingInsideLoop::Fill(i, key, w);
            return;
        }
    }
    throw std::invalid_argument("Histogram with name \""+std::string(hist_name)+"\" could not be found.");
}

//2d histogram, alphanumeric, with weight
void ProcessingInsideLoop::Fill(int hist_number, const char* key, double y, double w){
    if(hist2dtabLocal[hist_number]!=nullptr){
        hist2dtabLocal[hist_number]->Fill(key, y, w);
    } else{
        throw std::invalid_argument("1D histogram with number \""+std::to_string(hist_number)+"\" does not exist.");
    }
}
void ProcessingInsideLoop::Fill(const char* hist_name, const char* key, double y, double w){
    //looking through 2d histograms
    for(long unsigned int i = 0; i<hist2dtabLocal.size(); i++){
        if(hist2dtabLocal[i]!=nullptr&&strcmp(hist_name, hist2dtabLocal[i]->GetName())==0){
            ProcessingInsideLoop::Fill(i, key, y, w);
            return;
        }
    }
    throw std::invalid_argument("Histogram with name \""+std::string(hist_name)+"\" could not be found.");
}

//2d histogram, 2xalphanumeric, with weight
void ProcessingInsideLoop::Fill(int hist_number, const char* key_x, const char* key_y, double w){
    if(hist2dtabLocal[hist_number]!=nullptr){
        hist2dtabLocal[hist_number]->Fill(key_x, key_y, w);
    } else{
        throw std::invalid_argument("1D histogram with number \""+std::to_string(hist_number)+"\" does not exist.");
    }
}
void ProcessingInsideLoop::Fill(const char* hist_name, const char* key_x, const char* key_y, double w){
    //looking throught2d histograms
    for(long unsigned int i = 0; i<hist2dtabLocal.size(); i++){
        if(hist2dtabLocal[i]!=nullptr&&strcmp(hist_name, hist2dtabLocal[i]->GetName())==0){
            ProcessingInsideLoop::Fill(i, key_x, key_y, w);
            return;
        }
    }
    throw std::invalid_argument("Histogram with name \""+std::string(hist_name)+"\" could not be found.");
}

//3d histogram, alphanumeric, with weight
//TODO: expand
void ProcessingInsideLoop::Fill(int hist_number, const char* key, double y, double z, double w){
    if(hist3dtabLocal[hist_number]!=nullptr){
        hist3dtabLocal[hist_number]->Fill(key, y, z, w);
    } else{
        throw std::invalid_argument("3D histogram with number \""+std::to_string(hist_number)+"\" does not exist.");
    }
}
void ProcessingInsideLoop::Fill(const char* hist_name, const char* key, double y, double z, double w){
    //looking through 1d histograms
    for(long unsigned int i = 0; i<hist3dtabLocal.size(); i++){
        if(hist3dtabLocal[i]!=nullptr&&strcmp(hist_name, hist3dtabLocal[i]->GetName())==0){
            ProcessingInsideLoop::Fill(i, key, y, z, w);
            return;
        }
    }
    throw std::invalid_argument("Histogram with name \""+std::string(hist_name)+"\" could not be found.");
}

#endif