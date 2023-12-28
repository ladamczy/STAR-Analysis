#ifndef PARTIAL_HISTOGRAMS_H
#define PARTIAL_HISTOGRAMS_H

#include <StRPEvent.h>
#include <StUPCEvent.h>
#include <TH1.h>

class PartialHistograms{
private:
    StRPEvent *privateRPEventPointer;
    StUPCEvent *privateUPCEventPointer;
    std::vector<bool(*)(StRPEvent *, StUPCEvent *, TH1 *)> functionVector;
    std::vector<TH1 *> histVector;
    std::vector<bool> boolVector;
public:
    PartialHistograms(/* args */){}
    ~PartialHistograms(){}
    void GetEventPointers(StRPEvent *, StUPCEvent *);
    void AddCut(bool(StRPEvent *, StUPCEvent *, TH1 *), TH1 *);
    void ProcessEvent();
};

void PartialHistograms::AddCut(bool f(StRPEvent *, StUPCEvent *, TH1 *), TH1 *hist){
    functionVector.push_back(f);
    histVector.push_back(hist);
    boolVector.push_back(false);
}

void PartialHistograms::GetEventPointers(StRPEvent *RPEventPointer, StUPCEvent *UPCEventPointer){
    privateRPEventPointer = RPEventPointer;
    privateUPCEventPointer = UPCEventPointer;
}

void PartialHistograms::ProcessEvent(){
    for(long unsigned int i = 0; i<functionVector.size(); i++){
        boolVector[i] = functionVector[i](privateRPEventPointer, privateUPCEventPointer, nullptr);
    }
    if(std::count(boolVector.begin(), boolVector.end(), false)==1){
        int index_of_false = std::find(boolVector.begin(), boolVector.end(), false)-boolVector.begin();
        functionVector[index_of_false](privateRPEventPointer, privateUPCEventPointer, histVector[index_of_false]);
    }
}

#endif