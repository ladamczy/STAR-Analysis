#include "MapReadPicoLambdaK0.h"

struct MyData {
    Long64_t eventId;
    Float_t lead_pt;
    Float_t lead_phi;
    Float_t lead_eta;
    Float_t sublead_pt;
    Float_t sublead_phi;
    Float_t sublead_eta;
    Float_t p1_pt;
    Float_t p1_phi;
    Float_t p1_eta;
    Int_t p1_ch;
    Int_t p1_hasTOFinfo;
    Float_t p2_pt;
    Float_t p2_phi;
    Float_t p2_eta;
    Int_t p2_hasTOFinfo;
    Int_t pair_charge;
    Float_t pair_phi;
    Float_t pair_eta;
    Float_t pair_pt;
    Float_t pair_mass;
};

ReadPicoLambdaK0::ReadPicoLambdaK0(TTree* chain2) {
    std::vector<MyData> myDataVector;

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

    for (Long64_t j = 0; j < GetEntriesFast(); ++j) {
        chain2->GetEntry(j);
        MyData data;
        data.eventId = eventId;
        data.lead_pt = lead_pt;
        data.lead_phi = lead_phi;
        data.lead_eta = lead_eta;
        data.sublead_pt = sublead_pt;
        data.sublead_phi = sublead_phi;
        data.sublead_eta = sublead_eta;
        data.p1_pt = p1_pt;
        data.p1_phi = p1_phi;
        data.p1_eta = p1_eta;
        data.p1_ch = p1_ch;
        data.p1_hasTOFinfo = p1_hasTOFinfo;
        data.p2_pt = p2_pt;
        data.p2_phi = p2_phi;
        data.p2_eta = p2_eta;
        data.p2_hasTOFinfo = p2_hasTOFinfo;
        data.pair_charge = pair_charge;
        data.pair_phi = pair_phi;
        data.pair_eta = pair_eta;
        data.pair_pt = pair_pt;
        data.pair_mass = pair_mass;
        
        myDataVector.push_back(data);
    }

    for (const MyData& data : myDataVector) {
        eventDataMap[data.eventId] = {
            data.lead_pt,
            data.lead_phi,
            data.lead_eta,
            data.sublead_pt,
            data.sublead_phi,
            data.sublead_eta,
            data.p1_pt,
            data.p1_phi,
            data.p1_eta,
            data.p1_ch,
            data.p1_hasTOFinfo,
            data.p2_pt,
            data.p2_phi,
            data.p2_eta,
            data.p2_hasTOFinfo,
            data.pair_charge,
            data.pair_phi,
            data.pair_eta,
            data.pair_pt,
            data.pair_mass
        };
    }
}

void ReadPicoLambdaK0::ProcessData(Long64_t i, StUPCEvent* upcEvt, TChain *chain, TTree* chain2) {
    
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
    for(int jj = j ; jj < chain2->GetEntriesFast(); ++jj)
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
            //2624 przyporzadkowane
            matchedEvents.push_back(eventId);
        }
        else if (Event1 < Event2){
            //unmatchedEventsTree1s.push_back(Event1);
            break;
        }
        else if (Event1 > Event2){
            unmatchedEventsTree2.push_back(Event2);
            continue;
        }
        else{
            //unmatchedEventsTree1s.push_back(Event1);
            break;
        }

    }
}