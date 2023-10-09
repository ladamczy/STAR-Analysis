#include "ReadPicoLambdaK0.h"

ReadPicoLambdaK0::ReadPicoLambdaK0(TTree* chain2) {

    chain2->SetBranchAddress("eventId", &eventId);
    chain2->SetBranchAddress("lead_pt", &lead_pt);
    chain2->SetBranchAddress("lead_phi", &lead_phi);
    chain2->SetBranchAddress("lead_eta", &lead_eta);
    chain2->SetBranchAddress("sublead_pt", &sublead_pt);
    chain2->SetBranchAddress("sublead_phi", &sublead_phi);
    chain2->SetBranchAddress("sublead_eta", &sublead_eta);
    chain2->SetBranchAddress("p1_pt", &p1_pt);
    chain2->SetBranchAddress("p1_phi", &p1_phi);
    chain2->SetBranchAddress("p1_eta", &p1_eta);
    chain2->SetBranchAddress("p1_ch", &p1_ch);
    chain2->SetBranchAddress("p1_hasTOFinfo", &p1_hasTOFinfo);
    chain2->SetBranchAddress("p2_pt", &p2_pt);
    chain2->SetBranchAddress("p2_phi", &p2_phi);
    chain2->SetBranchAddress("p2_eta", &p2_eta);
    chain2->SetBranchAddress("p2_hasTOFinfo", &p2_hasTOFinfo);
    chain2->SetBranchAddress("pair_charge", &pair_charge);
    chain2->SetBranchAddress("pair_phi", &pair_phi);
    chain2->SetBranchAddress("pair_eta", &pair_eta);
    chain2->SetBranchAddress("pair_pt", &pair_pt);
    chain2->SetBranchAddress("pair_mass", &pair_mass);
}

void ReadPicoLambdaK0::ProcessData(Long64_t i, StUPCEvent* upcEvt, Tchain *chain, TTree* chain2) {
    eventIdVectors.clear();
    leadPtVectors.clear();
    leadPhiVectors.clear();
    leadEtaVectors.clear();
    subleadPtVectors.clear();
    subleadPhiVectors.clear();
    subleadEtaVectors.clear();
    p1PtVectors.clear();
    p1PhiVectors.clear();
    p1EtaVectors.clear();
    p1ChVectors.clear();
    p1HasTOFInfoVectors.clear();
    p2PtVectors.clear();
    p2PhiVectors.clear();
    p2EtaVectors.clear();
    p2HasTOFInfoVectors.clear();
    pairChargeVectors.clear();
    pairPhiVectors.clear();
    pairEtaVectors.clear();
    pairPtVectors.clear();
    pairMassVectors.clear();

    chain->GetEntry(i);
    Long64_t Event1 = upcEvt->getEventNumber();
    for(int jj = j ; jj < chain2->GetEntries(); ++jj)
    {
        chain2->GetEntry(jj);
        Long64_t Event2 = eventId;
        if (Event1 == Event2)
        {
            ++j;
            eventIdVectors.push_back(eventId);
            leadPtVectors.push_back(lead_pt);
            leadPhiVectors.push_back(lead_phi);
            leadEtaVectors.push_back(lead_eta);
            subleadPtVectors.push_back(sublead_pt);
            subleadPhiVectors.push_back(sublead_phi);
            subleadEtaVectors.push_back(sublead_eta);
            p1PtVectors.push_back(p1_pt);
            p1PhiVectors.push_back(p1_phi);
            p1EtaVectors.push_back(p1_eta);
            p1ChVectors.push_back(p1_ch);
            p1HasTOFInfoVectors.push_back(p1_hasTOFinfo);
            p2PtVectors.push_back(p2_pt);
            p2PhiVectors.push_back(p2_phi);
            p2EtaVectors.push_back(p2_eta);
            p2HasTOFInfoVectors.push_back(p2_hasTOFinfo);
            pairChargeVectors.push_back(pair_charge);
            pairPhiVectors.push_back(pair_phi);
            pairEtaVectors.push_back(pair_eta);
            pairPtVectors.push_back(pair_pt);
            pairMassVectors.push_back(pair_mass);
        }
        else if (Event1 < Event2)
            break;
        else if (Event1 > Event2)
            continue;
         else 
            break;
    }
}