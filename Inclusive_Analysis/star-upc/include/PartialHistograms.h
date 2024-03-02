#ifndef PARTIAL_HISTOGRAMS_H
#define PARTIAL_HISTOGRAMS_H

#include <StRPEvent.h>
#include <StUPCEvent.h>
#include <TH1.h>

class PartialHistograms{
private:
    StRPEvent *privateRPEventPointer;
    StUPCEvent *privateUPCEventPointer;
    std::vector<int> *privateAllowedRPTracks;
    std::vector<int> *privateAllowedUPCTracks;
    std::vector<bool(*)(StRPEvent *, StUPCEvent *, std::vector<int> *, std::vector<int> *)> necessaryCutVector;
    std::vector<bool(*)(StRPEvent *, StUPCEvent *, TH1 *, std::vector<int> *, std::vector<int> *)> histogramCutVector;
    //checks if whole event is deleted or just some tracks
    std::vector<bool> IsNecessaryCutFull;
    std::vector<bool> IsHistogramCutFull;
    std::vector<TH1 *> histVector;
public:
    PartialHistograms(/* args */){
        privateAllowedRPTracks = new std::vector<int>;
        privateAllowedUPCTracks = new std::vector<int>;
    }
    ~PartialHistograms(){}
    void SetEventPointers(StRPEvent *, StUPCEvent *);
    void AddNecessaryCut(bool (*f)(StRPEvent *, StUPCEvent *, std::vector<int> *, std::vector<int> *), bool);
    void AddHistogramCut(bool (*f)(StRPEvent *, StUPCEvent *, TH1 *, std::vector<int> *, std::vector<int> *), TH1 *, bool);
    void ProcessEvent();
};

void PartialHistograms::AddNecessaryCut(bool f(StRPEvent *, StUPCEvent *, std::vector<int> *, std::vector<int> *), bool DoesItCutOutWholeEvent){
    necessaryCutVector.push_back(f);
    IsNecessaryCutFull.push_back(DoesItCutOutWholeEvent);
}

void PartialHistograms::AddHistogramCut(bool f(StRPEvent *, StUPCEvent *, TH1 *, std::vector<int> *, std::vector<int> *), TH1 *hist, bool DoesItCutOutWholeEvent){
    histogramCutVector.push_back(f);
    histVector.push_back(hist);
    IsHistogramCutFull.push_back(DoesItCutOutWholeEvent);
}

void PartialHistograms::SetEventPointers(StRPEvent *RPEventPointer, StUPCEvent *UPCEventPointer){
    //cleaning
    privateAllowedRPTracks->clear();
    privateAllowedUPCTracks->clear();
    //filling again
    privateRPEventPointer = RPEventPointer;
    privateUPCEventPointer = UPCEventPointer;
    for(long unsigned int i = 0; i<privateRPEventPointer->getNumberOfTracks(); i++){
        privateAllowedRPTracks->push_back(i);
    }
    for(long unsigned int i = 0; i<privateUPCEventPointer->getNumberOfTracks(); i++){
        privateAllowedUPCTracks->push_back(i);
    }

}

void PartialHistograms::ProcessEvent(){
    //if there's something wrong, it just stops analysing that event, thanks to return

    //necessary cuts
    for(long unsigned int i = 0; i<necessaryCutVector.size(); i++){
        if(IsNecessaryCutFull[i]){
            if(!necessaryCutVector[i](privateRPEventPointer, privateUPCEventPointer, privateAllowedRPTracks, privateAllowedUPCTracks)){
                return;
            }
        } else{
            necessaryCutVector[i](privateRPEventPointer, privateUPCEventPointer, privateAllowedRPTracks, privateAllowedUPCTracks);
        }
    }
    //histogram cuts
    for(long unsigned int omittedCut = 0; omittedCut<histogramCutVector.size(); omittedCut++){
        //reset for identical starting conditions
        privateAllowedRPTracks->clear();
        privateAllowedUPCTracks->clear();
        for(int i = 0; i<privateRPEventPointer->getNumberOfTracks(); i++){
            privateAllowedRPTracks->push_back(i);
        }
        for(int i = 0; i<privateUPCEventPointer->getNumberOfTracks(); i++){
            privateAllowedUPCTracks->push_back(i);
        }
        //applying every cut except one
        for(long unsigned int i = 0; i<histogramCutVector.size(); i++){
            if(i==omittedCut){
                continue;
            }
            if(IsHistogramCutFull[i]){
                if(!histogramCutVector[i](privateRPEventPointer, privateUPCEventPointer, nullptr, privateAllowedRPTracks, privateAllowedUPCTracks)){
                    return;
                }
            } else{
                histogramCutVector[i](privateRPEventPointer, privateUPCEventPointer, nullptr, privateAllowedRPTracks, privateAllowedUPCTracks);
            }
        }
        //filling histogram associated with that one cut
        histogramCutVector[omittedCut](privateRPEventPointer, privateUPCEventPointer, histVector[omittedCut], privateAllowedRPTracks, privateAllowedUPCTracks);
    }
}

#endif