#include "Includes.h"

using namespace std;

int main(int argc, char** argv) {

    if (argc != 3) {
        cerr << "two input files required" << endl;
        return 1;
    }


    TFile file1(argv[1], "update");
    if (!file1.IsOpen()) {
        cerr << "Failed to open input file1." << endl;
        return 1;
    }

    TTree *UPCTree = static_cast<TTree*>(file1.Get("mUPCTree"));
    if (!UPCTree) {
        cerr << "Failed to retrieve mUPCTree tree." << endl;
        file1.Close();
        return 1;
    }

    static StUPCEvent *upcEvt = 0x0;
    UPCTree->SetBranchAddress("mUPCEvent", &upcEvt);

    vector<Int_t> eventIdVectors;
    vector<Float_t> leadPtVectors, leadPhiVectors, leadEtaVectors;
    vector<Float_t> subleadPtVectors, subleadPhiVectors, subleadEtaVectors;
    vector<Float_t> p1PtVectors, p1PhiVectors, p1EtaVectors;
    vector<Int_t> p1ChVectors, p1HasTOFInfoVectors;
    vector<Float_t> p2PtVectors, p2PhiVectors, p2EtaVectors;
    vector<Int_t> p2HasTOFInfoVectors;
    vector<Int_t> pairChargeVectors;
    vector<Float_t> pairPhiVectors, pairEtaVectors, pairPtVectors, pairMassVectors;

    TBranch *beventIdVectors = UPCTree->Branch("eventIdVectors", &eventIdVectors);
    TBranch *bleadPtVectors = UPCTree->Branch("leadPtVectors", &leadPtVectors);
    TBranch *bleadPhiVectors = UPCTree->Branch("leadPhiVectors", &leadPhiVectors);
    TBranch *bleadEtaVectors = UPCTree->Branch("leadEtaVectors", &leadEtaVectors);
    TBranch *bsubleadPtVectors = UPCTree->Branch("subleadPtVectors", &subleadPtVectors);
    TBranch *bsubleadPhiVectors = UPCTree->Branch("subleadPhiVectors", &subleadPhiVectors);
    TBranch *bsubleadEtaVectors = UPCTree->Branch("subleadEtaVectors", &subleadEtaVectors);
    TBranch *bp1PtVectors = UPCTree->Branch("p1PtVectors", &p1PtVectors);
    TBranch *bp1PhiVectors = UPCTree->Branch("p1PhiVectors", &p1PhiVectors);
    TBranch *bp1EtaVectors = UPCTree->Branch("p1EtaVectors", &p1EtaVectors);
    TBranch *bp1ChVectors = UPCTree->Branch("p1ChVectors", &p1ChVectors);
    TBranch *bp1HasTOFInfoVectors = UPCTree->Branch("p1HasTOFInfoVectors", &p1HasTOFInfoVectors);
    TBranch *bp2PtVectors = UPCTree->Branch("p2PtVectors", &p2PtVectors);
    TBranch *bp2PhiVectors = UPCTree->Branch("p2PhiVectors", &p2PhiVectors);
    TBranch *bp2EtaVectors = UPCTree->Branch("p2EtaVectors", &p2EtaVectors);
    TBranch *bp2HasTOFInfoVectors = UPCTree->Branch("p2HasTOFInfoVectors", &p2HasTOFInfoVectors);
    TBranch *bpairChargeVectors = UPCTree->Branch("pairChargeVectors", &pairChargeVectors);
    TBranch *bpairPhiVectors = UPCTree->Branch("pairPhiVectors", &pairPhiVectors);
    TBranch *bpairEtaVectors = UPCTree->Branch("pairEtaVectors", &pairEtaVectors);
    TBranch *bpairPtVectors = UPCTree->Branch("pairPtVectors", &pairPtVectors);
    TBranch *bpairMassVectors = UPCTree->Branch("pairMassVectors", &pairMassVectors);

    const char* inputFile2 = argv[2];
    TFile* file2 = TFile::Open(inputFile2);
    if (!file2) {
        cerr << "Failed to open second input file2." << endl;
        file1.Close();
        return 1;
    }
    TTree* tree_2 = static_cast<TTree*>(file2->Get("ntp_K0s"));


    ReadPicoLambdaK0 Read_K0(tree_2);

    for (Long64_t i = 0; i < UPCTree->GetEntries(); ++i) {

        Read_K0.ProcessData(i, upcEvt, UPCTree, tree_2);

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
    file1.cd();
    // UPCTree->Write();
    UPCTree->Write("", TObject::kOverwrite);
    file1.Close();
    file2->Close();
    // delete upcEvt;

    return 0;
}