#ifndef PROCESSING_INSIDE_LOOP_H
#define PROCESSING_INSIDE_LOOP_H

#include "ProcessingOutsideLoop.h"

#include <utility>

#include <type_traits>

class ProcessingInsideLoop {
private:
    std::vector<std::shared_ptr<ProcessingOutsideLoop::Wrapper<TH1>>> histTH1tabLocal;
    std::vector<std::shared_ptr<ProcessingOutsideLoop::Wrapper<THnBase>>> histTHnBasetabLocal;
    std::vector<std::shared_ptr<TEfficiency>> histTEfficiencytabLocal;
    std::map<std::string, int> histogramMap;
public:
    ProcessingInsideLoop(/* args */);
    ~ProcessingInsideLoop();
    void GetLocalHistograms(ProcessingOutsideLoop *);
    template<typename... Args>
    void Fill(int, Args...);
    template<typename... Args>
    void Fill(const char*, Args...);
};

ProcessingInsideLoop::ProcessingInsideLoop(/* args */) {
}

ProcessingInsideLoop::~ProcessingInsideLoop() {
}

void ProcessingInsideLoop::GetLocalHistograms(ProcessingOutsideLoop *outsideloop){
    //using the first one for counting, as they all are the same size
    for(long unsigned int i = 0; i<outsideloop->histTH1tab.size(); i++){
        histTH1tabLocal.push_back(outsideloop->histTH1tab[i]!=nullptr ? outsideloop->histTH1tab[i]->Get() : nullptr);
        histTHnBasetabLocal.push_back(outsideloop->histTHnBasetab[i]!=nullptr ? outsideloop->histTHnBasetab[i]->Get() : nullptr);
        histTEfficiencytabLocal.push_back(outsideloop->histTEfficiencytab[i]!=nullptr ? outsideloop->histTEfficiencytab[i]->Get() : nullptr);
    }
    //getting the map too
    histogramMap = outsideloop->histogramMap;
}

//helpers to check if Fill function will work properly
//first one - default false
template <typename T, typename, typename... Args>
struct can_call_fill : std::false_type{};
//second one - the one that actually checks
template <typename T, typename... Args>
struct can_call_fill<T, std::void_t<decltype(std::declval<T*>()->Fill(std::declval<Args>()...))>, Args...>
    : std::true_type{};
//wrapper for those two above
template <typename T, typename... Args>
inline constexpr bool can_call_fill_v = can_call_fill<T, void, Args...>::value;

//helper to check if there is any text in the args
template <typename... Args>
inline constexpr bool has_string_arg_v = (
    std::is_convertible_v<Args, std::string>||...
    );

template<typename ...Args>
void ProcessingInsideLoop::Fill(int hist_number, Args ...args){
    //TH1 and derivative classes, with TH2 and TH3 the encompass all things stemming from TH1
    //TH3 was needed for some particular filling method with many arguments
    if(histTH1tabLocal[hist_number]!=nullptr){
        if constexpr(can_call_fill_v<TH1, Args...>){
            dynamic_cast<TH1*>(histTH1tabLocal[hist_number]->GetInternalHist())->Fill(std::forward<Args>(args)...);
        } else if constexpr(can_call_fill_v<TH2, Args...>){
            dynamic_cast<TH2*>(histTH1tabLocal[hist_number]->GetInternalHist())->Fill(std::forward<Args>(args)...);
        } else if constexpr(can_call_fill_v<TH3, Args...>){
            dynamic_cast<TH3*>(histTH1tabLocal[hist_number]->GetInternalHist())->Fill(std::forward<Args>(args)...);
        }
        return;
    }
    //THnBase - here check was needed whether it has words (const char*) - Fill method here does not
    //So it would try to convert an array of doubles (double*) and cause problems
    if(histTHnBasetabLocal[hist_number]!=nullptr){
        //compiler was complaining, so here is a safeguard against trying to fill THnBase with something weird
        if constexpr(can_call_fill_v<THnBase, Args...>&&!has_string_arg_v<Args...>)
            histTHnBasetabLocal[hist_number]->GetInternalHist()->Fill(std::forward<Args>(args)...);
        return;
    }
    //TEfficiency - everything is straightforward here
    if(histTEfficiencytabLocal[hist_number]!=nullptr){
        //compiler was complaining, so here is a safeguard against trying to fill TEff. with one parameter only
        if constexpr(can_call_fill_v<TEfficiency, Args...>)
            histTEfficiencytabLocal[hist_number]->Fill(std::forward<Args>(args)...);
        return;
    }
}

template<typename ...Args>
void ProcessingInsideLoop::Fill(const char* hist_name, Args ...args){
    if(histogramMap.find(hist_name)!=histogramMap.end()){
        Fill(histogramMap[hist_name], std::forward<Args>(args)...);
        return;
    }
    throw std::invalid_argument("Object with name \""+std::string(hist_name)+"\" could not be found.");
}

#endif