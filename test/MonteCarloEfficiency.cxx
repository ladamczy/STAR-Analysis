#include "Includes.h"

using namespace std;

int main(int argc, char** argv)  
{
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

//    while (std::getline(inputFilePathList, inputFileName))
//    {
//        std::cout << inputFileName.c_str() << std::endl; 
        chain->AddFile(argv[1]);
//    }
//    inputFilePathList.close();

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
    TH1D* HistKaonMDet = new TH1D("HistKaonMDet", "; M^{#pi^{+}#pi^{-}} [GeV]; # events", 100, 0.4, 0.6);

    TH1D* HistXM4T = new TH1D("HistXM4T", "; M^{#pi^{+}#pi^{-}} [GeV]; # events",21, 0.9, 3.0);
    TH1D* HistXM2P = new TH1D("HistXM2P", "; M^{#pi^{+}#pi^{-}} [GeV]; # events",21, 0.9, 3.0);

    TH1D* HistV0MDet = new TH1D("HistV0MDet", "; M^{V0} [GeV]; # events", 100, 0.4, 0.6);
    TH1D* HistV0DCAD = new TH1D("HistV0DCAD", "; DCAD [cm]; # events", 100, 0.0, 5.0);
    TH1D* HistV0R = new TH1D("HistV0R", "; R [cm]; # events", 100, 0.0, 5.0);

    TH1D* HistPYPY = new TH1D("HistPYPY", "; Py*Py [GeV^2]; # events", 10, -5.0, 5.0);

    TParticle* particle;
    TParticle* PosPion1; TParticle* NegPion1; TParticle*  PosPion2; TParticle*  NegPion2;    
    vector <TParticle*> PosPions, NegPions, Protons, PosPionsPaired, NegPionsPaired;
    vector <StUPCTrack*> tpcTrack;
    vector <StUPCTrack*> tpcTrackUM;
    Int_t ToFPair = 0;
    Int_t ToFTrack= 0;
    TLorentzVector lorentzVectorPosPion, lorentzVectorNegPion;
    TLorentzVector lorentzVectorNeutKaon;
    TLorentzVector trackVector;

    double truthVertexR, truthVertexZ, truthEta, truthPt;
    double detVertexR, detVertexZ, detEta, detPt, detM, V0M;

//    std::cout << "entries " << chain->GetEntries() << std::endl;

    for (Long64_t i = 0; i < chain->GetEntries(); ++i) 
    {
       TLorentzVector lorentzVectorX;
       chain->GetEntry(i);

//     cout << "ntracks total  " << upcEvt->getNumberOfTracks() << endl;

     for (int j = 0; j < upcEvt->getNumberOfTracks(); j++) {
             for (int jj = j; jj < upcEvt->getNumberOfTracks(); jj++) {       
               if( upcEvt->getTrack(j)->getCharge() != upcEvt->getTrack(jj)->getCharge()) {
        TVector3 vertex(0,0,0);
        int id1,id2;
        StUPCV0 V0(upcEvt->getTrack(j),upcEvt->getTrack(jj), massPion, massPion,id1,id2, vertex, upcEvt->getMagneticField(), false, false);
        HistV0DCAD->Fill(V0.dcaDaughters());
        HistV0R->Fill(sqrt(V0.v0x()*V0.v0x()+V0.v0y()*V0.v0y()));
        HistV0MDet->Fill(V0.m()); 
               }
          }
     }

      Double_t PYPY=1;

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
                PYPY = PYPY*particle->Py();
            }
        }


        HistPYPY->Fill(PYPY);
        
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
        cout << " true level " << PosPionsPaired.size() << " " << NegPionsPaired.size()  << endl;
        if (PosPionsPaired.size()!=2 or NegPionsPaired.size() != 2 ) 
        {
            PosPions.clear();
            NegPions.clear();
            PosPionsPaired.clear();
            NegPionsPaired.clear();
            Protons.clear();
            ToFTrack=0;
            ToFPair=0;
            tpcTrackUM.clear();
            continue;
        }
        
        for (int i = 0; i < PosPionsPaired.size(); i++)
        {   
            vector <StUPCTrack*> tpcTrack;
            int num = 0;
            int numToF=0;
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
       	    ToFTrack=0;
       	    ToFPair=0;
                tpcTrackUM.clear();
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
       	    ToFTrack=0;
       	    ToFPair=0;
                tpcTrackUM.clear();
                continue;
            }

            lorentzVectorNeutKaon = lorentzVectorPosPion + lorentzVectorNegPion;
            lorentzVectorX +=lorentzVectorNeutKaon;
            TLorentzVector productionVertex;
            PosPionsPaired[i]->ProductionVertex(productionVertex);

            truthVertexR = sqrt(pow(productionVertex.X(),2) + pow(productionVertex.Y(),2));
//            cout << "vertex  " << productionVertex.X() << " " << productionwVertex.Y() << " " << productionVertex.Z() << endl; 
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
                    if ( abs(lorentzVectorPosPion.Eta() - upcEvt->getTrack(j)->getEta()) < 0.1  and  abs(lorentzVectorPosPion.Phi() - upcEvt->getTrack(j)->getPhi()) < 0.1) 
                    {
                        tpcTrack.push_back(upcEvt->getTrack(j));
                        tpcTrackUM.push_back(upcEvt->getTrack(j));
                        num+=1;
                        if(upcEvt->getTrack(j)->getFlag(StUPCTrack::kTof)) numToF+=1;
                    }

                    if ( abs(lorentzVectorNegPion.Eta() - upcEvt->getTrack(j)->getEta()) < 0.1  and  abs(lorentzVectorNegPion.Phi() - upcEvt->getTrack(j)->getPhi()) < 0.1) 
                    {
                        tpcTrack.push_back(upcEvt->getTrack(j));
                        tpcTrackUM.push_back(upcEvt->getTrack(j));
                        num+=1;
       	       	       	if(upcEvt->getTrack(j)->getFlag(StUPCTrack::kTof)) numToF+=1;
                    }
                }         

                // phi - eta consistency with truth level for exacly two tracks    
                if(numToF>0) ToFPair+=1;
                ToFTrack+=numToF;           
                if (num != 2) 
                {       
                    PosPions.clear();
                    NegPions.clear();
                    PosPionsPaired.clear();
                    NegPionsPaired.clear();
                    Protons.clear();
                    tpcTrack.clear();
       	    ToFTrack=0;
       	    ToFPair=0;
                    tpcTrackUM.clear();
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
       	    ToFTrack=0;
       	    ToFPair=0;
                tpcTrackUM.clear();
                continue;
            }

            TLorentzVector tpcKaon = {0,0,0,0};
            int c = 0;
            double vertexTpcR = 0;
            double vertexTpcZ = 0;


//            cout << "number of tracks " << tpcTrack.size() << endl;
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
//                cout << "tpc vert " << tpcTrack[i]->getVertex()->getPosX() << " " 
//                                    << tpcTrack[i]->getVertex()->getPosY() << "	"
//                                    << tpcTrack[i]->getVertex()->getPosZ() << endl;

                vertexTpcR+=sqrt( pow(tpcTrack[i]->getVertex()->getPosX() ,2)+ pow(tpcTrack[i]->getVertex()->getPosY() ,2) )/2.0;
                vertexTpcZ+=tpcTrack[i]->getVertex()->getPosZ()/2.0;
            }    


//	cout << "number of tracks, c " << tpcTrack.size() << " " << c << endl;;

   
            if (c!=2) 
            {
                PosPions.clear();
                NegPions.clear();
                PosPionsPaired.clear();
                NegPionsPaired.clear();
                tpcTrack.clear();
       	    ToFTrack=0;
       	    ToFPair=0;
                tpcTrackUM.clear();
                Protons.clear();
                continue;
            }


            TVector3 vertex(tpcTrack[0]->getVertex()->getPosX(),tpcTrack[0]->getVertex()->getPosY(),tpcTrack[0]->getVertex()->getPosZ());
            TVector3 vertex2(tpcTrack[1]->getVertex()->getPosX(),tpcTrack[1]->getVertex()->getPosY(),tpcTrack[1]->getVertex()->getPosZ());
//            cout << "paricle1 vertex " << vertex.X() << " " << vertex2.X() << endl;
            Int_t id1, id2;
//upcEvt->getMagneticField()
            StUPCV0 V0(tpcTrack[0],tpcTrack[1], massPion, massPion,id1,id2, vertex, upcEvt->getMagneticField(), false, false);
            cout << "dca in the pair " << V0.dcaDaughters() << "len " << V0.decayLength() << "dcs to PV " 
                 << V0.DcaToPrimaryVertex() << "mass " << V0.m() << "dca1 2 " << V0.particle1Dca() << " " << 
                 V0.particle2Dca() << endl;
 

            detVertexR = vertexTpcR;
            detVertexZ = vertexTpcZ;
            detPt = tpcKaon.Pt();
            detEta = tpcKaon.Eta();
            detM = tpcKaon.M();

            if ( abs(detEta) < 1.2 and detVertexR < 3.0 and abs(detVertexZ) < 200)
            {
                HistKaonPtDet->Fill(detPt);
                HistKaonMDet->Fill(detM);
                HistV0MDet->Fill(V0.m());
                HistV0DCAD->Fill(V0.dcaDaughters());
                HistV0R->Fill(sqrt(V0.v0x()*V0.v0x()+V0.v0y()*V0.v0y()));
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


        cout << "Event " << PosPions.size() << " " << NegPions.size() << " " << 
                            PosPionsPaired.size() << " " << NegPionsPaired.size() << " " <<
                            tpcTrack.size() << " " << tpcTrackUM.size() << "ToF Tracks " << ToFTrack << 
                            "ToF Pairs " << ToFPair << endl;

        if(ToFTrack==4&&tpcTrackUM.size()==4&&abs(truthVertexZ)<80) HistXM4T->Fill(lorentzVectorX.M());
        if(ToFPair==2&&tpcTrackUM.size()==4&&abs(truthVertexZ)<80) HistXM2P->Fill(lorentzVectorX.M());

        PosPions.clear();
        NegPions.clear();
        PosPionsPaired.clear();
        NegPionsPaired.clear();
        tpcTrack.clear();
       	    ToFTrack=0;
       	    ToFPair=0;
        tpcTrackUM.clear();
        Protons.clear();   
    
    }
    
    TFile *outfile = TFile::Open(argv[2], "recreate"); 

    HistXM4T->Write();
    HistXM2P->Write();
    HistKaonPtTruth->Write();
    HistKaonEtaTruth->Write();
    HistKaonVtxRTruth->Write();
    HistKaonVtxZTruth->Write();
    
    HistKaonPtDet->Write();
    HistKaonEtaDet->Write();
    HistKaonVtxRDet->Write();
    HistKaonVtxZDet->Write();
    HistKaonMDet->Write();
    HistV0MDet->Write();
    HistV0DCAD->Write();
    HistV0R->Write();
    HistPYPY->Write();
    outfile->Close();

    return 0;
}
