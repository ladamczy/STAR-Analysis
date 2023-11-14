#include "Includes.h"

using namespace std;

int main(int argc, char** argv)  {
    if (argc != 3) 
    {
        cerr << "two input files required" << std::endl;
        return 1;
    }


    TFile* file1 = TFile::Open(argv[1]);
    if (!file1) 
    {
        cerr << "Failed to open input file." << endl;
        return 1;
    }
    TTree* tree_1 = static_cast<TTree*>(file1->Get("mUPCTree")); 


    ////
    static StUPCEvent *upcEvt = 0x0;
    tree_1->SetBranchAddress("mUPCEvent", &upcEvt);
    ////

    std::vector<Int_t> eventIdVectors;
    std::vector<Float_t> leadPtVectors, leadPhiVectors, leadEtaVectors;
    std::vector<Float_t> subleadPtVectors, subleadPhiVectors, subleadEtaVectors;
    std::vector<Float_t> p1PtVectors, p1PhiVectors, p1EtaVectors;
    std::vector<Int_t> p1ChVectors, p1HasTOFInfoVectors;
    std::vector<Float_t> p2PtVectors, p2PhiVectors, p2EtaVectors;
    std::vector<Int_t> p2HasTOFInfoVectors;
    std::vector<Int_t> pairChargeVectors;
    std::vector<Float_t> pairPhiVectors, pairEtaVectors, pairPtVectors, pairMassVectors;

    TBranch *beventIdVectors = tree_1->Branch("eventIdVectors", &eventIdVectors);
    TBranch *bleadPtVectors = tree_1->Branch("leadPtVectors", &leadPtVectors);
    TBranch *bleadPhiVectors = tree_1->Branch("leadPhiVectors", &leadPhiVectors);
    TBranch *bleadEtaVectors = tree_1->Branch("leadEtaVectors", &leadEtaVectors);
    TBranch *bsubleadPtVectors = tree_1->Branch("subleadPtVectors", &subleadPtVectors);
    TBranch *bsubleadPhiVectors = tree_1->Branch("subleadPhiVectors", &subleadPhiVectors);
    TBranch *bsubleadEtaVectors = tree_1->Branch("subleadEtaVectors", &subleadEtaVectors);
    TBranch *bp1PtVectors = tree_1->Branch("p1PtVectors", &p1PtVectors);
    TBranch *bp1PhiVectors = tree_1->Branch("p1PhiVectors", &p1PhiVectors);
    TBranch *bp1EtaVectors = tree_1->Branch("p1EtaVectors", &p1EtaVectors);
    TBranch *bp1ChVectors = tree_1->Branch("p1ChVectors", &p1ChVectors);
    TBranch *bp1HasTOFInfoVectors = tree_1->Branch("p1HasTOFInfoVectors", &p1HasTOFInfoVectors);
    TBranch *bp2PtVectors = tree_1->Branch("p2PtVectors", &p2PtVectors);
    TBranch *bp2PhiVectors = tree_1->Branch("p2PhiVectors", &p2PhiVectors);
    TBranch *bp2EtaVectors = tree_1->Branch("p2EtaVectors", &p2EtaVectors);
    TBranch *bp2HasTOFInfoVectors = tree_1->Branch("p2HasTOFInfoVectors", &p2HasTOFInfoVectors);
    TBranch *bpairChargeVectors = tree_1->Branch("pairChargeVectors", &pairChargeVectors);
    TBranch *bpairPhiVectors = tree_1->Branch("pairPhiVectors", &pairPhiVectors);
    TBranch *bpairEtaVectors = tree_1->Branch("pairEtaVectors", &pairEtaVectors);
    TBranch *bpairPtVectors = tree_1->Branch("pairPtVectors", &pairPtVectors);
    TBranch *bpairMassVectors = tree_1->Branch("pairMassVectors", &pairMassVectors);


    const char* inputFile2 = argv[2];
    TFile* file2 = TFile::Open(inputFile2);
    if (!file2) 
    {
        cerr << "Failed to open second input file." << endl;
        return 1;
    }
    TTree* tree_2 = static_cast<TTree*>(file2->Get("ntp_K0s"));


    ReadPicoLambdaK0 Read_K0(tree_2);

    for (Long64_t i = 0; i < tree_1->GetEntries(); ++i) {

        Read_K0.ProcessData(i, upcEvt, tree_1, tree_2);

        eventIdVectors = Read_K0.eventIdVectors;
        leadPtVectors = Read_K0.leadPtVectors;
        leadPhiVectors = Read_K0.leadPhiVectors;
        leadEtaVectors = Read_K0.leadEtaVectors;
        subleadPtVectors = Read_K0.subleadPtVectors;
        subleadPhiVectors = Read_K0.subleadPhiVectors;
        subleadEtaVectors = Read_K0.subleadEtaVectors;
        p1PtVectors = Read_K0.p1PtVectors;
        p1PhiVectors = Read_K0.p1PhiVectors;
        p1EtaVectors = Read_K0.p1EtaVectors;
        p1ChVectors = Read_K0.p1ChVectors;
        p1HasTOFInfoVectors = Read_K0.p1HasTOFInfoVectors;
        p2PtVectors = Read_K0.p2PtVectors;
        p2PhiVectors = Read_K0.p2PhiVectors;
        p2EtaVectors = Read_K0.p2EtaVectors;
        p2HasTOFInfoVectors = Read_K0.p2HasTOFInfoVectors;
        pairChargeVectors = Read_K0.pairChargeVectors;
        pairPhiVectors = Read_K0.pairPhiVectors;
        pairEtaVectors = Read_K0.pairEtaVectors;
        pairPtVectors = Read_K0.pairPtVectors;
        pairMassVectors = Read_K0.pairMassVectors;

        beventIdVectors->Fill();
        bleadPtVectors->Fill();
        bleadPhiVectors->Fill();
        bleadEtaVectors->Fill();
        bsubleadPtVectors->Fill();
        bsubleadPhiVectors->Fill();
        bsubleadEtaVectors->Fill();
        bp1PtVectors->Fill();
        bp1PhiVectors->Fill();
        bp1EtaVectors->Fill();
        bp1ChVectors->Fill();
        bp1HasTOFInfoVectors->Fill();
        bp2PtVectors->Fill();
        bp2PhiVectors->Fill();
        bp2EtaVectors->Fill();
        bp2HasTOFInfoVectors->Fill();
        bpairChargeVectors->Fill();
        bpairPhiVectors->Fill();
        bpairEtaVectors->Fill();
        bpairPtVectors->Fill();
        bpairMassVectors->Fill();
    }   

    tree_1->Write();
    file1->Close();
    file2->Close();
   


    return 0;
}