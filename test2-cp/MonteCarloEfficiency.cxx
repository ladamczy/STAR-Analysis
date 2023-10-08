#include "Includes.h"

using namespace std;

int main(int argc, char** argv)  
{
    if (argc != 3 && argc != 4) 
    {
        cerr << "two or three input files required (input ouput input2)" << std::endl;
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
    string inputFileName;
    vector <string> rootFiles;
    while (std::getline(inputFilePathList, inputFileName))
    {
        chain->Add(inputFileName.c_str());
    }
    inputFilePathList.close();

    const char* inputFile2 = argv[3];
    TFile* file2 = TFile::Open(inputFile2);
    TTree* chain2 = static_cast<TTree*>(file2->Get("ntp_K0s"));

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

    vector<double> eventIdVectors;
    vector<double> leadPtVectors;
    vector<double> leadPhiVectors;
    vector<double> leadEtaVectors;
    vector<double> subleadPtVectors;
    vector<double> subleadPhiVectors;
    vector<double> subleadEtaVectors;
    vector<double> p1PtVectors;
    vector<double> p1PhiVectors;
    vector<double> p1EtaVectors;
    vector<double> p1ChVectors;
    vector<double> p1HasTOFInfoVectors;
    vector<double> p2PtVectors;
    vector<double> p2PhiVectors;
    vector<double> p2EtaVectors;
    vector<double> p2HasTOFInfoVectors;
    vector<double> pairChargeVectors;
    vector<double> pairPhiVectors;
    vector<double> pairEtaVectors;
    vector<double> pairPtVectors;
    vector<double> pairMassVectors;

    Int_t eventId, p1_hasTOFinfo, p2_hasTOFinfo, pair_charge, p1_ch;
    Double_t lead_pt, lead_phi, lead_eta, 
             sublead_pt, sublead_phi, sublead_eta,
             p1_pt, p1_phi, p1_eta,
             p2_pt, p2_phi, p2_eta,
             pair_phi, pair_eta, pair_pt, pair_mass;
             
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

    int j = 0;

    for (Long64_t i = 0; i < chain->GetEntries(); ++i) 
    {
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
        Long64_t Event1 = upcEvt->GetEventNumber();

        //cout << Event1 << endl;                  debug 1
        for(int jj = j ; jj < chain2->GetEntries(); ++jj)
        {
            chain2->GetEntry(jj);
            Long64_t Event2 = eventId;
            //cout << Event2 << endl;              debug 2
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
        //cout << sizeof(leadPtVectors) << endl;                  debug 3


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
    
    TFile *outfile = TFile::Open(argv[2], "recreate"); 

    HistKaonPtTruth->Write();
    HistKaonEtaTruth->Write();
    HistKaonVtxRTruth->Write();
    HistKaonVtxZTruth->Write();
    
    HistKaonPtDet->Write();
    HistKaonEtaDet->Write();
    HistKaonVtxRDet->Write();
    HistKaonVtxZDet->Write();
    outfile->Close();

    return 0;
}
