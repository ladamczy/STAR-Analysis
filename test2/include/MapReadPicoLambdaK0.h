#ifndef MAP_READPICO_LAMBDA_K0_H
#define MAP_READPICO_LAMBDA_K0_H

#include "Includes.h"

class MapReadPicoLambdaK0 {
public:
    MapReadPicoLambdaK0(TTree* chain2);
    void ProcessData(Long64_t i, StUPCEvent* upcEvt, TChain *chain, TTree* chain2);

    const std::map<Long64_t, std::vector<double>>& getEventDataMap() const;

private:
    Int_t eventId, p1_hasTOFinfo, p2_hasTOFinfo, pair_charge, p1_ch;
    Float_t lead_pt, lead_phi, lead_eta,
             sublead_pt, sublead_phi, sublead_eta,
             p1_pt, p1_phi, p1_eta,
             p2_pt, p2_phi, p2_eta,
             pair_phi, pair_eta, pair_pt, pair_mass;

    std::map<Long64_t, std::vector<double>> eventDataMap;
};

#endif  // MAP_READPICO_LAMBDA_K0_H
