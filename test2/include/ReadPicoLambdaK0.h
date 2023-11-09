#ifndef READ_PICO_LAMBDA_K0_H
#define READ_PICO_LAMBDA_K0_H

#include <Includes.h>

class ReadPicoLambdaK0 {
public:
    ReadPicoLambdaK0(TTree* chain2);
    void ProcessData(Long64_t i, StUPCEvent* upcEvt, TChain *chain, TTree* chain2);

    std::vector<Int_t> eventIdVectors;
    std::vector<Float_t> leadPtVectors;
    std::vector<Float_t> leadPhiVectors;
    std::vector<Float_t> leadEtaVectors;
    std::vector<Float_t> subleadPtVectors;
    std::vector<Float_t> subleadPhiVectors;
    std::vector<Float_t> subleadEtaVectors;
    std::vector<Float_t> p1PtVectors;
    std::vector<Float_t> p1PhiVectors;
    std::vector<Float_t> p1EtaVectors;
    std::vector<Int_t> p1ChVectors;
    std::vector<Int_t> p1HasTOFInfoVectors;
    std::vector<Float_t> p2PtVectors;
    std::vector<Float_t> p2PhiVectors;
    std::vector<Float_t> p2EtaVectors;
    std::vector<Int_t> p2HasTOFInfoVectors;
    std::vector<Int_t> pairChargeVectors;
    std::vector<Float_t> pairPhiVectors;
    std::vector<Float_t> pairEtaVectors;
    std::vector<Float_t> pairPtVectors;
    std::vector<Float_t> pairMassVectors;

    //std::vector<Long64_t> matchedEvents;

private:
    Int_t eventId, p1_hasTOFinfo, p2_hasTOFinfo, pair_charge, p1_ch;
    Float_t lead_pt, lead_phi, lead_eta,
            sublead_pt, sublead_phi, sublead_eta,
            p1_pt, p1_phi, p1_eta,
            p2_pt, p2_phi, p2_eta,
            pair_phi, pair_eta, pair_pt, pair_mass;
    int j=0;
};

#endif  // READ_PICO_LAMBDA_K0_H
