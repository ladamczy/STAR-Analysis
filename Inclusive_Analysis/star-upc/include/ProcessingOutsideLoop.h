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
    void AddHistogram(const TH1&);
    void AddHistogram(const THnBase&);
    void AddHistogram(const TEfficiency&);
    void PrintHistograms();
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
    //returns whether insertion was successful (e.g. the key provided does not exist in the map)
    bool insertion_result = histogramMap.emplace(hist.GetName(), histTH1tab.size()-1).second;

    //making sure that name was correct
    if(std::string(hist.GetName()).length()==0){
        throw std::invalid_argument("Histogram added "+std::string(histTH1tab.size()==0 ? "as the first one" : "after \""+histogramMap.rbegin()->first+"\"")+" has an empty name!");
    } else if(!insertion_result){
        throw std::invalid_argument("Histogram of name \""+std::string(hist.GetName())+"\" already exists!");
    }
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
    //returns whether insertion was successful (e.g. the key provided does not exist in the map)
    bool insertion_result = histogramMap.emplace(hist.GetName(), histTHnBasetab.size()-1).second;

    //making sure that name was correct
    if(std::string(hist.GetName()).length()==0){
        throw std::invalid_argument("Histogram added "+std::string(histTH1tab.size()==0 ? "as the first one" : "after \""+histogramMap.rbegin()->first+"\"")+" has an empty name!");
    } else if(!insertion_result){
        throw std::invalid_argument("Histogram of name \""+std::string(hist.GetName())+"\" already exists!");
    }
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
    //returns whether insertion was successful (e.g. the key provided does not exist in the map)
    bool insertion_result = histogramMap.emplace(hist.GetName(), histTEfficiencytab.size()-1).second;

    //making sure that name was correct
    if(std::string(hist.GetName()).length()==0){
        throw std::invalid_argument("Histogram added "+std::string(histTH1tab.size()==0 ? "as the first one" : "after \""+histogramMap.rbegin()->first+"\"")+" has an empty name!");
    } else if(!insertion_result){
        throw std::invalid_argument("Histogram of name \""+std::string(hist.GetName())+"\" already exists!");
    }
}

void ProcessingOutsideLoop::PrintHistograms(){
    std::vector<std::string> types, names;
    for(long unsigned int i = 0; i<histTH1tab.size(); i++){
        if(histTH1tab[i]!=nullptr){
            types.push_back(histTH1tab[i]->Get()->GetInternalHist()->ClassName());
            names.push_back(histTH1tab[i]->Get()->GetInternalHist()->GetName());
        } else if(histTHnBasetab[i]!=nullptr){
            types.push_back(histTHnBasetab[i]->Get()->GetInternalHist()->ClassName());
            names.push_back(histTHnBasetab[i]->Get()->GetInternalHist()->GetName());
        } else if(histTEfficiencytab[i]!=nullptr){
            types.push_back(histTEfficiencytab[i]->Get().get()->ClassName());
            names.push_back(histTEfficiencytab[i]->Get().get()->GetName());
        } else{
            throw std::invalid_argument("Histogram #"+to_string(i)+" is somehow represented by the wrong amount of pointers (should be exactly one)!");
        }
    }
    int max_len_number = ceil(log10(histTH1tab.size()-1));
    int max_len_type = std::max_element(types.begin(), types.end(), [](const std::string& a, const std::string& b){return a.length()<b.length();})->length();
    int max_len_name = std::max_element(names.begin(), names.end(), [](const std::string& a, const std::string& b){return a.length()<b.length();})->length();

    //printf with the - option allows for padding to the left; | will provide nice column divider; example: "%-9s|%-9s|%-9s\n"
    //2,4,4 are lengths of names of columns - "ID", "Type", "Name"
    std::string printf_options = "%-"+to_string(std::max(max_len_number, 2)+3)+"s%-"+to_string(std::max(max_len_type, 4)+2)+"s%-"+to_string(std::max(max_len_name, 4)+2)+"s\n";
    int total_length = std::max(max_len_number, 2)+3+std::max(max_len_type, 4)+2+std::max(max_len_name, 4)+2;

    //actual printing
    printf(printf_options.c_str(), "ID", "Type", "Name");
    for(size_t i = 0; i<total_length; i++) putchar('-'); printf("\n");
    for(size_t i = 0; i<types.size(); i++){
        printf(printf_options.c_str(), to_string(i).c_str(), types[i].c_str(), names[i].c_str());
    };
}

void ProcessingOutsideLoop::Merge(){
    //using the first one for counting, as they all are the same size
    for(long unsigned int i = 0; i<histTH1tab.size(); i++){
        //checking that exactly one pointer is not a nullptr
        if((histTH1tab[i]!=nullptr)+(histTHnBasetab[i]!=nullptr)+(histTEfficiencytab[i]!=nullptr)!=1){
            throw std::invalid_argument("Histogram #"+to_string(i)+" is somehow represented by the wrong amount of pointers (should be exactly one)!");
        }
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
        //checking that exactly one pointer is not a nullptr
        if((histTH1tabFinal[i]!=nullptr)+(histTHnBasetabFinal[i]!=nullptr)+(histTEfficiencytabFinal[i]!=nullptr)!=1){
            throw std::invalid_argument("Histogram #"+to_string(i)+" after calling Merge() is somehow represented by the wrong amount of pointers (should be exactly one)!");
        }
    }
}

int ProcessingOutsideLoop::GetNumberOfHistograms(){
    return histTH1tab.size();
}

void ProcessingOutsideLoop::SaveToFile(TFile* file){
    //checking the file opening
    if(!file->IsOpen()||file->IsZombie()){
        throw std::invalid_argument("File \""+std::string(file->GetName())+"\" is wrongly opened or something!");
    }

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
    //checking that exactly one pointer is not a nullptr
    if((histTH1tab[i]!=nullptr)+(histTHnBasetab[i]!=nullptr)+(histTEfficiencytab[i]!=nullptr)!=1){
        throw std::invalid_argument("Histogram #"+to_string(i)+" after calling GetPointer() is somehow represented by the wrong amount of pointers (should be exactly one)!");
    }

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
    //checking that exactly one pointer is not a nullptr
    if((histTH1tabFinal[i]!=nullptr)+(histTHnBasetabFinal[i]!=nullptr)+(histTEfficiencytabFinal[i]!=nullptr)!=1){
        throw std::invalid_argument("Histogram #"+to_string(i)+" after calling GetPointerAfterMerge() is somehow represented by the wrong amount of pointers (should be exactly one)!");
    }

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