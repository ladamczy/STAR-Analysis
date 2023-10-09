#ifndef READ_PICO_LAMBDA_K0_H
#define READ_PICO_LAMBDA_K0_H

#include <Includes.h>

class ReadPicoLambdaK0 {
public:
    ReadPicoLambdaK0(TTree* chain2);
    void ProcessData(Long64_t i, StUPCEvent* upcEvt, Tchain *chain, TTree* chain2);

    std::vector<double> eventIdVectors;
    std::vector<double> leadPtVectors;
    std::vector<double> leadPhiVectors;
    std::vector<double> leadEtaVectors;
    std::vector<double> subleadPtVectors;
    std::vector<double> subleadPhiVectors;
    std::vector<double> subleadEtaVectors;
    std::vector<double> p1PtVectors;
    std::vector<double> p1PhiVectors;
    std::vector<double> p1EtaVectors;
    std::vector<double> p1ChVectors;
    std::vector<double> p1HasTOFInfoVectors;
    std::vector<double> p2PtVectors;
    std::vector<double> p2PhiVectors;
    std::vector<double> p2EtaVectors;
    std::vector<double> p2HasTOFInfoVectors;
    std::vector<double> pairChargeVectors;
    std::vector<double> pairPhiVectors;
    std::vector<double> pairEtaVectors;
    std::vector<double> pairPtVectors;
    std::vector<double> pairMassVectors;

private:
    Int_t eventId, p1_hasTOFinfo, p2_hasTOFinfo, pair_charge, p1_ch;
    Double_t lead_pt, lead_phi, lead_eta,
             sublead_pt, sublead_phi, sublead_eta,
             p1_pt, p1_phi, p1_eta,
             p2_pt, p2_phi, p2_eta,
             pair_phi, pair_eta, pair_pt, pair_mass;
    int j=0;
};

#endif  // READ_PICO_LAMBDA_K0_H