#ifndef PROCESSING_OUTSIDE_LOOP_H
#define PROCESSING_OUTSIDE_LOOP_H

#include "TClass.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THnBase.h"
#include "TEfficiency.h"
#include "TObject.h"
#include <stdio.h>
#include <string>
#include <ROOT/TThreadedObject.hxx>
#include "TFile.h"

#include <memory>

//to avoid circular reference we "dummy define" class here i guess? like with functions
class ProcessingInsideLoop;

class ProcessingOutsideLoop {
private:
    //so that Inside loop instances would have access to these things
    friend ProcessingInsideLoop;
    //wrappers for abstract classes for TThreadedObject (needs to be a descendant of TObject)
    template<class T>
    class Wrapper : public TObject{
        //so that Inside loop instances would have access to these things
        T* inside_hist = nullptr;
    public:
        //constructors
        Wrapper(): inside_hist(nullptr){}
        Wrapper(T* h): inside_hist(h){}
        Wrapper(const TObject& h){
            this->inside_hist = static_cast<T*>(h.Clone());
        }

        //copy constructor - for TThreadedObject
        Wrapper(const Wrapper& other){
            if(other.inside_hist!=nullptr){
                this->inside_hist = static_cast<T*>(other.inside_hist->Clone());
            } else{
                this->inside_hist = nullptr;
            }
        }

        //move constructor
        Wrapper(Wrapper&& other) noexcept: inside_hist(other.inside_hist){
            other.inside_hist = nullptr;
        }

        //destructor
        ~Wrapper(){
            delete inside_hist;
        }

        //method needed for TList to be overwritten
        Bool_t IsFolder() const override{ return kFALSE; }

        //for getting pointer to the object this wrapper contains
        T* GetInternalHist(){ return inside_hist; }

        //Merge method
        Long64_t Merge(TCollection* list){
            if(!inside_hist||!list) return -1;

            //we need to strip elements from their Wrappers
            TList rawList;
            TIter next(list);
            while(TObject* obj = next()){
                Wrapper<T>* processed_wrapper = dynamic_cast<Wrapper<T>*>(obj);
                //the part that does the stripping
                if(processed_wrapper&&processed_wrapper->GetInternalHist()){
                    rawList.Add(processed_wrapper->GetInternalHist());
                }
            }

            //merging of actual histograms
            return inside_hist->Merge(&rawList);
        }
    };
    //internal memory of objects from main three classes of histograms
    std::vector<ROOT::TThreadedObject<Wrapper<TH1>>*> histTH1tab;
    std::vector<ROOT::TThreadedObject<Wrapper<THnBase>>*> histTHnBasetab;
    std::vector<ROOT::TThreadedObject<TEfficiency>*> histTEfficiencytab;
    std::vector<std::shared_ptr<TH1>> histTH1tabFinal;
    std::vector<std::shared_ptr<THnBase>> histTHnBasetabFinal;
    std::vector<std::shared_ptr<TEfficiency>> histTEfficiencytabFinal;
    //map where each histogram is, for faster lookup
    std::map<std::string, int> histogramMap;
public:
    //constructors
    ProcessingOutsideLoop(){};
    ~ProcessingOutsideLoop(){};

    //functions
    // void AddHistogram(const TObject&);
    void AddHistogram(const TH1&);
    void AddHistogram(const THnBase&);
    void AddHistogram(const TEfficiency&);
    void Merge();
    void SaveToFile(TFile*);
    int GetNumberOfHistograms();

    //members
    template<class T>
    T* GetPointer(int);
    template<class T>
    T* GetPointer(const char*);
    template<class T>
    std::shared_ptr<T> GetPointerAfterMerge(int);
    template<class T>
    std::shared_ptr<T> GetPointerAfterMerge(const char*);
};

void ProcessingOutsideLoop::AddHistogram(const TH1& hist){
    //filling vectors with dummy pointers
    histTH1tab.push_back(nullptr);
    histTHnBasetab.push_back(nullptr);
    histTEfficiencytab.push_back(nullptr);
    //replacing pointers and invoking Sumw2 for correct error calculation
    histTH1tab.back() = new ROOT::TThreadedObject<Wrapper<TH1>>(hist);
    //lone Get() needed to populate histograms first
    histTH1tab.back()->Get();
    histTH1tab.back()->Get()->GetInternalHist()->Sumw2();
    histogramMap.emplace(hist.GetName(), histTH1tab.size()-1);
}
void ProcessingOutsideLoop::AddHistogram(const THnBase& hist){
    //filling vectors with dummy pointers
    histTH1tab.push_back(nullptr);
    histTHnBasetab.push_back(nullptr);
    histTEfficiencytab.push_back(nullptr);
    //replacing pointers and invoking Sumw2 for correct error calculation
    histTHnBasetab.back() = new ROOT::TThreadedObject<Wrapper<THnBase>>(hist);
    //lone Get() needed to populate histograms first
    histTHnBasetab.back()->Get();
    histTHnBasetab.back()->Get()->GetInternalHist()->Sumw2();
    histogramMap.emplace(hist.GetName(), histTHnBasetab.size()-1);
}
void ProcessingOutsideLoop::AddHistogram(const TEfficiency& hist){
    //filling vectors with dummy pointers
    histTH1tab.push_back(nullptr);
    histTHnBasetab.push_back(nullptr);
    histTEfficiencytab.push_back(nullptr);
    //replacing pointers and invoking Sumw2 for correct error calculation
    histTEfficiencytab.back() = new ROOT::TThreadedObject<TEfficiency>(hist);
    //lone Get() needed to populate histograms first
    histTEfficiencytab.back()->Get();
    histogramMap.emplace(hist.GetName(), histTEfficiencytab.size()-1);
}

void ProcessingOutsideLoop::Merge(){
    //using the first one for counting, as they all are the same size
    for(long unsigned int i = 0; i<histTH1tab.size(); i++){
        //filling vectors with dummy pointers
        histTH1tabFinal.push_back(nullptr);
        histTHnBasetabFinal.push_back(nullptr);
        histTEfficiencytabFinal.push_back(nullptr);
        //we use so-called "aliasing constructor" to get TH1 out of its wrapper while still having it be a shared pointer
        if(histTH1tab[i]!=nullptr){
            std::shared_ptr<Wrapper<TH1>> temp_wrapper_TH1 = histTH1tab[i]->Merge();
            histTH1tabFinal[i] = std::shared_ptr<TH1>(temp_wrapper_TH1, temp_wrapper_TH1->GetInternalHist());
        }
        if(histTHnBasetab[i]!=nullptr){
            std::shared_ptr<Wrapper<THnBase>> temp_wrapper_THnBase = histTHnBasetab[i]->Merge();
            histTHnBasetabFinal[i] = std::shared_ptr<THnBase>(temp_wrapper_THnBase, temp_wrapper_THnBase->GetInternalHist());
        }
        histTEfficiencytabFinal[i] = histTEfficiencytab[i]!=nullptr ? histTEfficiencytab[i]->Merge() : nullptr;
    }
}

int ProcessingOutsideLoop::GetNumberOfHistograms(){
    return histTH1tab.size();
}

void ProcessingOutsideLoop::SaveToFile(TFile* file){
    file->cd();
    //using the first one for counting, as they all are the same size
    for(long unsigned int i = 0; i<histTH1tabFinal.size(); i++){
        if(histTH1tabFinal[i]!=nullptr)
            histTH1tabFinal[i]->Write();
        if(histTHnBasetabFinal[i]!=nullptr)
            histTHnBasetabFinal[i]->Write();
        if(histTEfficiencytabFinal[i]!=nullptr)
            histTEfficiencytabFinal[i]->Write();
    }
}

template<class T>
T* ProcessingOutsideLoop::GetPointer(int i){
    //dynamic_cast gives back nullptr if trying to cast uncastable type
    //hence why I'm using it everywhere despite inability to cast TEff. into TH1D
    if(histTH1tab[i]!=nullptr)
        return dynamic_cast<T*>(histTH1tab[i]->Get()->GetInternalHist());
    if(histTHnBasetab[i]!=nullptr)
        return dynamic_cast<T*>(histTHnBasetab[i]->Get()->GetInternalHist());
    if(histTEfficiencytab[i]!=nullptr)
        return dynamic_cast<T*>(histTEfficiencytab[i]->Get().get());

    //in case nothing worked
    return nullptr;
}
template<class T>
T* ProcessingOutsideLoop::GetPointer(const char* hist_name){
    if(histogramMap.find(hist_name)!=histogramMap.end())
        return GetPointer<T>(histogramMap[hist_name]);
    throw std::invalid_argument("Object with name \""+std::string(hist_name)+"\" could not be found.");
}

template<class T>
std::shared_ptr<T> ProcessingOutsideLoop::GetPointerAfterMerge(int i){
    if(histTH1tabFinal[i]!=nullptr)
        return dynamic_pointer_cast<T>(histTH1tabFinal[i]);
    if(histTHnBasetabFinal[i]!=nullptr)
        return dynamic_pointer_cast<T>(histTHnBasetabFinal[i]);
    if(histTEfficiencytabFinal[i]!=nullptr)
        return dynamic_pointer_cast<T>(histTEfficiencytabFinal[i]);

    //if something goes wrong
    return nullptr;
}
template<class T>
std::shared_ptr<T> ProcessingOutsideLoop::GetPointerAfterMerge(const char* hist_name){
    if(histogramMap.find(hist_name)!=histogramMap.end())
        return GetPointerAfterMerge<T>(histogramMap[hist_name]);
    throw std::invalid_argument("Object with name \""+std::string(hist_name)+"\" could not be found.");
}

#endif