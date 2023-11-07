#include "Includes.h"

using namespace std;

int main(int argc, char** argv)  
{
    if (argc != 3 && argc != 4) 
    {
        cerr << "three input files required ([input] [input2] [ouput])" << std::endl;
        return 1;
    }

    double massPion = 0.13957061;
    double massKaon =  497.611;

    ifstream inputFilePathList(argv[1]);
    if (!inputFilePathList) 
    {
        cerr << "Failed to open input file." << std::endl;
        return 1;
    }

    TChain *chain = new TChain("mUPCTree"); 
    // string inputFileName;
    // vector <string> rootFiles;
    // while (std::getline(inputFilePathList, inputFileName))
    // {
    //    chain->Add(inputFileName.c_str());
    chain->AddFile(argv[1]);
    // }
    // inputFilePathList.close();

    static StUPCEvent *upcEvt = 0x0;
    chain->SetBranchAddress("mUPCEvent", &upcEvt);

    TH1D* HistKaonPtTruth = new TH1D("HistKaonPtTruth", "; pT_{K^{0}} [GeV]; # events", 25 ,0, 2.5);
    TH1D* HistKaonEtaTruth = new TH1D("HistKaonEtaTruth", "; #eta_{K^{0}}; # events", 25 ,-4.5, 4.5);
    TH1D* HistKaonVtxZTruth = new TH1D("HistKaonVtxZTruth", ";z_{vtx}^{#pi^{+}#pi^{-}} [cm]; # events", 10, -150, 150);
    TH1D* HistKaonVtxRTruth = new TH1D("HistKaonVtxRTruth", ";R_{vtx}^{#pi^{+}#pi^{-}} [cm]; # events", 10, 0, 20);

    TH1D* HistKaonPtDet = new TH1D("HistKaonPtDet", "; pT_{K^{0}} [GeV]; # events", 25 ,0, 2.5);
    TH1D* HistKaonEtaDet = new TH1D("HistKaonEtaDet", "; #eta_{K^{0}}; # events", 25 ,-4.5, 4.5);
    TH1D* HistKaonVtxZDet = new TH1D("HistKaonVtxZDet", "; z_{vtx}^{#pi^{+}#pi^{-}} [cm]; # events", 10, -150, 150);
    TH1D* HistKaonVtxRDet = new TH1D("HistKaonVtxRDet", "; R_{vtx}^{#pi^{+}#pi^{-}} [cm]; # events", 10, 0, 20);


    TParticle* particle;
    TParticle* PosPion1; TParticle* NegPion1; TParticle*  PosPion2; TParticle*  NegPion2;
    vector <TParticle*> PosPions, NegPions, Protons, PosPionsPaired, NegPionsPaired;
    vector <StUPCTrack*> tpcTrack;

    TLorentzVector lorentzVectorPosPion, lorentzVectorNegPion;
    TLorentzVector lorentzVectorNeutKaon;
    TLorentzVector trackVector;

    double truthVertexR, truthVertexZ, truthEta, truthPt;
    double detVertexR, detVertexZ, detEta, detPt;
             
    const char* inputFile2 = argv[2];
    TFile* file2 = TFile::Open(inputFile2);
    TTree* chain2 = static_cast<TTree*>(file2->Get("ntp_K0s"));

    ReadPicoLambdaK0 Read_K0(chain2);
    int xx=0;

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

    for (Long64_t i = 0; i < chain->GetEntries(); ++i) 
    {
        Read_K0.ProcessData(i, upcEvt, chain, chain2);

        // if(i <50) {
        //     xx = 0;
        //     for (size_t j = 0; j < Read_K0.eventIdVectors.size(); ++j){
        //         xx++;
        //         std::cout << Read_K0.eventIdVectors[j] << " ";
        //     }
        //     if(xx>0)
        //         std::cout << std::endl;
        // }    

        const std::vector<Int_t>& eventIdVectors = Read_K0.eventIdVectors;
        const std::vector<Float_t>& leadPtVectors = Read_K0.leadPtVectors;
        const std::vector<Float_t>& leadPhiVectors = Read_K0.leadPhiVectors;
        const std::vector<Float_t>& leadEtaVectors = Read_K0.leadEtaVectors;
        const std::vector<Float_t>& subleadPtVectors = Read_K0.subleadPtVectors;
        const std::vector<Float_t>& subleadPhiVectors = Read_K0.subleadPhiVectors;
        const std::vector<Float_t>& subleadEtaVectors = Read_K0.subleadEtaVectors;
        const std::vector<Float_t>& p1PtVectors = Read_K0.p1PtVectors;
        const std::vector<Float_t>& p1PhiVectors = Read_K0.p1PhiVectors;
        const std::vector<Float_t>& p1EtaVectors = Read_K0.p1EtaVectors;
        const std::vector<Int_t>& p1ChVectors = Read_K0.p1ChVectors;
        const std::vector<Int_t>& p1HasTOFInfoVectors = Read_K0.p1HasTOFInfoVectors;
        const std::vector<Float_t>& p2PtVectors = Read_K0.p2PtVectors;
        const std::vector<Float_t>& p2PhiVectors = Read_K0.p2PhiVectors;
        const std::vector<Float_t>& p2EtaVectors = Read_K0.p2EtaVectors;
        const std::vector<Int_t>& p2HasTOFInfoVectors = Read_K0.p2HasTOFInfoVectors;
        const std::vector<Int_t>& pairChargeVectors = Read_K0.pairChargeVectors;
        const std::vector<Float_t>& pairPhiVectors = Read_K0.pairPhiVectors;
        const std::vector<Float_t>& pairEtaVectors = Read_K0.pairEtaVectors;
        const std::vector<Float_t>& pairPtVectors = Read_K0.pairPtVectors;
        const std::vector<Float_t>& pairMassVectors = Read_K0.pairMassVectors;
    

        // extract all Pi+, Pi- and diffractive protons
        for (int i = 0; i < upcEvt->getNumberOfMCParticles(); i++)
        {
            particle = upcEvt->getMCParticle(i);
        
            if (particle->GetPDG()->PdgCode() == 211)
            {
                PosPions.push_back(particle);
            }

            else if (particle->GetPDG()->PdgCode() == -211)
            {
                NegPions.push_back(particle);
            }
            
            if (particle->GetPDG()->PdgCode() == 2212 and particle->GetFirstMother() == 1)
            {
                Protons.push_back(particle);
            }
        }

        //only paired pions - based on the GetFirstMother() method
        for (int i = 0; i < PosPions.size(); i++)
        {
            for (int j = 0; j < NegPions.size(); j++)
            {
                if (PosPions[i]->GetFirstMother() == NegPions[j]->GetFirstMother())
                {
                    PosPionsPaired.push_back(PosPions[i]);
                    NegPionsPaired.push_back(NegPions[j]);
                }
            }
        }

        if (PosPionsPaired.size() == 0 or NegPionsPaired.size() == 0 ) 
        {
            PosPions.clear();
            NegPions.clear();
            PosPionsPaired.clear();
            NegPionsPaired.clear();
            Protons.clear();
            tpcTrack.clear();
            continue;
        }
        
        for (int i = 0; i < PosPionsPaired.size(); i++)
        {   
            vector <StUPCTrack*> tpcTrack;
            int num = 0;

            lorentzVectorPosPion.SetPxPyPzE(PosPionsPaired[i]->Px(), PosPionsPaired[i]->Py(), PosPionsPaired[i]->Pz(), PosPionsPaired[i]->Energy());
            lorentzVectorNegPion.SetPxPyPzE(NegPionsPaired[i]->Px(), NegPionsPaired[i]->Py(), NegPionsPaired[i]->Pz(), NegPionsPaired[i]->Energy());

            // if at least one pion is outside of the acceptance then proceed to the next pair...
            if ( lorentzVectorPosPion.Pt() <= 0.1 or abs(lorentzVectorPosPion.Eta()) >= 1 )
            {
                PosPions.clear();
                NegPions.clear();
                PosPionsPaired.clear();
                NegPionsPaired.clear();
                Protons.clear();
                tpcTrack.clear();
                continue;
            }
            
            if ( lorentzVectorNegPion.Pt() <= 0.1 or abs(lorentzVectorNegPion.Eta()) >= 1 )
            {
                PosPions.clear();
                NegPions.clear();
                PosPionsPaired.clear();
                NegPionsPaired.clear();
                Protons.clear();
                tpcTrack.clear();
                continue;
            }

            lorentzVectorNeutKaon = lorentzVectorPosPion + lorentzVectorNegPion;

            TLorentzVector productionVertex;
            PosPionsPaired[i]->ProductionVertex(productionVertex);

            truthVertexR = sqrt(pow(productionVertex.X(),2) + pow(productionVertex.Y(),2));
            truthVertexZ = productionVertex.Z();
            truthPt = lorentzVectorNeutKaon.Pt();
            truthEta = lorentzVectorNeutKaon.Eta();

            if ( abs(truthEta) < 1.2 and truthVertexR < 3.0 and abs(truthVertexZ) < 200)
            {
                HistKaonPtTruth->Fill(truthPt);
            }

            if ( truthVertexR < 3.0 and abs(truthVertexZ) < 200)
            {
                HistKaonEtaTruth->Fill(truthEta);
            }

            if ( abs(truthEta) < 1.2 and truthVertexR < 3.0 )
            {
                HistKaonVtxZTruth->Fill(truthVertexZ);
            }

            if (( abs(truthEta) < 1.2 and abs(truthVertexZ) < 200) )
            {
                HistKaonVtxRTruth->Fill(truthVertexR);
            }

            // detector level - at least two tracks 
            if ( upcEvt->getNumberOfTracks() >= 2)
            {               
                
                for (int j = 0; j < upcEvt->getNumberOfTracks(); j++)
                {
                    if ( abs(lorentzVectorPosPion.Eta() - upcEvt->getTrack(j)->getEta()) < 0.1  and  abs(lorentzVectorPosPion.Phi() - upcEvt->getTrack(j)->getPhi()) < 0.1  and upcEvt->getTrack(i)->getFlag(StUPCTrack::kTof)) 
                    {
                        tpcTrack.push_back(upcEvt->getTrack(j));
                        num+=1;
                    }

                    if ( abs(lorentzVectorNegPion.Eta() - upcEvt->getTrack(j)->getEta()) < 0.1  and  abs(lorentzVectorNegPion.Phi() - upcEvt->getTrack(j)->getPhi()) < 0.1  and upcEvt->getTrack(i)->getFlag(StUPCTrack::kTof)) 
                    {
                        tpcTrack.push_back(upcEvt->getTrack(j));
                        num+=1;
                    }
                }         

                // phi - eta consistency with truth level for exacly two tracks                               
                if (num != 2) 
                {       
                    PosPions.clear();
                    NegPions.clear();
                    PosPionsPaired.clear();
                    NegPionsPaired.clear();
                    Protons.clear();
                    tpcTrack.clear();
                    continue;
                }
            }

            else 
            {       
                PosPions.clear();
                NegPions.clear();
                PosPionsPaired.clear();
                NegPionsPaired.clear();
                Protons.clear();
                tpcTrack.clear();
                continue;
            }
            
            TLorentzVector tpcKaon = {0,0,0,0};

            int c = 0;
            double vertexTpcR = 0;
            double vertexTpcZ = 0;

            for (int i = 0; i < tpcTrack.size(); i++)
            {
                tpcTrack[i]->getLorentzVector(trackVector, massPion);
                if ( tpcTrack[i]->getCharge() < 0.0)
                {
                    tpcKaon+=trackVector;
                    c+=1;
                }

                else
                {
                    tpcKaon+=trackVector;
                    c+=1;
                }
                vertexTpcR+=sqrt( pow(tpcTrack[i]->getVertex()->getPosX() ,2)+ pow(tpcTrack[i]->getVertex()->getPosY() ,2) )/2.0;
                vertexTpcZ+=tpcTrack[i]->getVertex()->getPosZ()/2.0;
            }    
   
            if (c!=2) 
            {
                PosPions.clear();
                NegPions.clear();
                PosPionsPaired.clear();
                NegPionsPaired.clear();
                tpcTrack.clear();
                Protons.clear();
                continue;
            }

            detVertexR = vertexTpcR;
            detVertexZ = vertexTpcZ;
            detPt = tpcKaon.Pt();
            detEta = tpcKaon.Eta();

            if ( abs(detEta) < 1.2 and detVertexR < 3.0 and abs(detVertexZ) < 200)
            {
                HistKaonPtDet->Fill(detPt);
            }

            if ( detVertexR < 3.0 and abs(detVertexZ) < 200)
            {
                HistKaonEtaDet->Fill(detEta);
            }

            if ( abs(detEta) < 1.2 and detVertexR < 3.0 )
            {
                HistKaonVtxZDet->Fill(detVertexZ);
            }

            if (( abs(detEta) < 1.2 and abs(detVertexZ) < 200) )
            {
                HistKaonVtxRDet->Fill(detVertexR);
            }
            
        }    

        PosPions.clear();
        NegPions.clear();
        PosPionsPaired.clear();
        NegPionsPaired.clear();
        tpcTrack.clear();
        Protons.clear();   
    
    }
    //********************************************
    // TFile *outfile = TFile::Open(argv[3], "recreate"); 

    // HistKaonPtTruth->Write();
    // HistKaonEtaTruth->Write();
    // HistKaonVtxRTruth->Write();
    // HistKaonVtxZTruth->Write();
    
    // HistKaonPtDet->Write();
    // HistKaonEtaDet->Write();
    // HistKaonVtxRDet->Write();
    // HistKaonVtxZDet->Write();
    // outfile->Close();
    //********************************************


    return 0;
}
