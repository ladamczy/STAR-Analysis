#include "Includes.h"

using namespace std;

int main(int argc, char** argv)  
{
    double massPion = 0.13957061;
    double massKaon =  0.497611;
    double massProton = 0.93827;     

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

//    TH1D* HistV0MDet = new TH1D("HistV0MDet", "; M^{V0} [GeV]; # events", 100, 0.4, 0.6);
//    TH1D* HistV0MDetK0 = new TH1D("HistV0MDetK0", "; M^{V0} [GeV]; # events", 100, 0.4,  0.6);
//    TH2D* HistV0MvsDecayLDet = new TH2D("HistV0MvsDecayLDet", "; M^{V0} [GeV]; Decay Length", 50, 0., 25., 100, 0.4, 0.6);
    TH1D* HistV0MDet = new TH1D("HistV0MDet", "; M^{V0} [GeV]; # events", 100, 1.0, 1.2);
    TH2D* HistV0MvsDecayLDet = new TH2D("HistV0MvsDecayLDet", "; M^{V0} [GeV]; Decay Length", 50, 0., 25., 100, 1.0, 1.2);
    TH1D* HistV0MDetK0 = new TH1D("HistV0MDetK0", "; M^{V0} [GeV]; # events", 100, 1.0,  1.2);

    TH1D* HistV0DCAD = new TH1D("HistV0DCAD", "; DCAD [cm]; # events", 100, 0.0, 5.0);
    TH1D* HistV0R = new TH1D("HistV0R", "; R [cm]; # events", 150, 0.0, 15.0);
    TH1D* HistV0DCABeamLine = new TH1D("HistV0DCABeamLine", "; DCABeamLine [cm]; # events", 150, 0.0, 15.0);
    TH1D* HistV0PointingAngle = new TH1D("HistV0PointingAngle", "; PointingAngle ; # events", 200, -1.0, 1.0);
    TH1D* HistV0DecayLength = new TH1D("HistV0DecayLength", "; DecayLength [cm]; # events", 150, 0.0, 15.0);

    TH1D* HistV0MDetRejected = new TH1D("HistV0MDetRejected", "; M^{V0} [GeV]; # events", 100, 0.4, 0.6);
    TH1D* HistV0MDetRejectedDCA = new TH1D("HistV0MDetRejectedDCA", "; M^{V0} [GeV]; # events", 100, 0.4, 0.6);
    TH1D* HistV0MDetRejectedPA = new TH1D("HistV0MDetRejectedPA", "; M^{V0} [GeV]; # events", 100, 0.4, 0.6);
    TH1D* HistV0MDetRejectedDL = new TH1D("HistV0MDetRejectedDL", "; M^{V0} [GeV]; # events", 100, 0.4, 0.6);
    TH1D* HistV0DCADK0 = new TH1D("HistV0DCADK0", "; DCAD [cm]; # events", 100, 0.0, 5.0);
    TH1D* HistV0RK0 = new TH1D("HistV0RK0", "; R [cm]; # events", 150, 0.0, 15.0);
    TH1D* HistV0DCABeamLineK0 = new TH1D("HistV0DCABeamLineK0", "; DCABeamLine [cm]; # events", 150, 0.0, 15.0);
    TH1D* HistV0PointingAngleK0 = new TH1D("HistV0PointingAngleK0", "; PointingAngle ; # events", 200, -1.0, 1.0);
    TH1D* HistV0DecayLengthK0 = new TH1D("HistV0DecayLengthK0", "; DecayLength [cm]; # events", 150, 0.0, 15.0);
    TH1D* HistMatch = new TH1D("HistMatch", "; Delta; # events", 100, 0.0, 1.0);

    TH1D* HistDeltaZ = new TH1D("HistDeltaZ", "; Delta Z; # events", 200, -5.0, 5.0);
    TH1D* HistDeltaX = new TH1D("HistDeltaX", "; Delta X; # events", 200, -0.5, 0.5);
    TH2D* HistDeltaL = new TH2D("HistDeltaL", "; Delta L; # events",100,0.0,20.0,100,0.0,20.0);
    TH1D* HistDeltaR = new TH1D("HistDeltaR", "; Delta R; # events", 200, -5.0, 5.0);


    TH1D* HistPYPY = new TH1D("HistPYPY", "; Py*Py [GeV^2]; # events", 10, -5.0, 5.0);

    TParticle* particle;
    TParticle* PosPion1; TParticle* NegPion1; TParticle*  PosPion2; TParticle*  NegPion2;    
    vector <TParticle*> PosPions, NegPions, Protons, PosPionsPaired, NegPionsPaired, TrueK0;
    vector <StUPCTrack*> tpcTrack;
    vector <StUPCTrack*> tpcTrackUM;
    Int_t ToFPair = 0;
    Int_t ToFTrack= 0;
    TLorentzVector lorentzVectorPosPion, lorentzVectorNegPion;
    TLorentzVector lorentzVectorNeutKaon;
    TLorentzVector trackVector;
    Double_t PYPY=1;

    int good_events = 0;
    double truthVertexR, truthVertexZ, truthEta, truthPt;
    double detVertexR, detVertexZ, detEta, detPt, detM, V0M;

    double truthDecayVertexR, truthProdVertxZ, truthDecayLength;

//    std::cout << "entries " << chain->GetEntries() << std::endl;

    vector <vector<double>> fillNumberWithPosition = ReadFillPositionData("../share/Run7PolarizationWithPosition.csv");
    
//    cout << "Size of Fill file " << fillNumberWithPosition[0].size() << endl;
    int fillNumberOld = -1;
    double x0=0,y0=0,xs=0,ys=0;

    for (Long64_t i = 0; i < chain->GetEntries(); ++i) 
//    for (Long64_t i = 0; i < 1000; ++i)
    {
       TLorentzVector lorentzVectorX;
       chain->GetEntry(i);

       if( fillNumberOld != upcEvt->getFillNumber() ) {
	fillNumberOld = upcEvt->getFillNumber();

	for (unsigned int i = 0; i < fillNumberWithPosition[0].size(); i++)
	{
//             cout << fillNumberWithPosition[0][i] << endl;
		if (upcEvt->getFillNumber()  == fillNumberWithPosition[0][i] )
		{
			x0 = fillNumberWithPosition[1][i];	
			y0 = fillNumberWithPosition[2][i];
                        xs = fillNumberWithPosition[5][i];
                        ys = fillNumberWithPosition[6][i];
                        break;
		}
	}
      } 

//       cout << "Fill " << upcEvt->getFillNumber() << " " << x0 << " " << y0 << " " << xs << " " << ys << endl;
      double beamline[4] = {x0, y0, xs, ys} ;
// double beamline[4] = {0, 0, 0, 0};
     cout << "ntracks total  " << upcEvt->getNumberOfTracks() << endl;

        TrueK0.clear();

        // extract all K0, Pi+, Pi- and diffractive protons
        for (int i = 0; i < upcEvt->getNumberOfMCParticles(); i++)
        {
            particle = upcEvt->getMCParticle(i);
//            cout << "MC " << i << "PDG Code " << particle->GetPDG()->PdgCode() << "Mother " << particle->GetFirstMother() << endl;

//            cout << "Code " << particle->GetPDG()->PdgCode() << endl; 
    
        if (particle->GetPDG()->PdgCode() == 310)
////           if (   abs(particle->GetPDG()->PdgCode()) == 3122)
           {
           TrueK0.push_back(particle);
           }

        if (particle->GetPDG()->PdgCode() == 211 && abs(particle->Eta())<1.0 && particle->Pt()>0.1 )
////            if (particle->GetPDG()->PdgCode() == 2212 and particle->GetFirstMother() != 1 or particle->GetPDG()->PdgCode() == 211)
            {
             	PosPions.push_back(particle);
            }

            else if (particle->GetPDG()->PdgCode() == -211 && abs(particle->Eta())<1.0 && particle->Pt()>0.1 )
////            else if (particle->GetPDG()->PdgCode() == -2212 and particle->GetFirstMother() != 1 or particle->GetPDG()->PdgCode() == -211)
            {
             	NegPions.push_back(particle);
            }

            if (particle->GetPDG()->PdgCode() == 2212 and particle->GetFirstMother() == 1)
            {
             	Protons.push_back(particle);
                PYPY = PYPY*particle->Py();
            }
	}

     if (TrueK0.size()==2 && PosPions.size()==2 && NegPions.size()==2)  good_events++;
     
     cout << "True KO, pos, neg " << TrueK0.size() << " " << PosPions.size() << " " << NegPions.size() 
          << " good events " << good_events << endl;

     for (int j = 0; j < upcEvt->getNumberOfTracks(); j++) {
             for (int jj = j; jj < upcEvt->getNumberOfTracks(); jj++) {       
               if( (upcEvt->getTrack(j)->getCharge() != upcEvt->getTrack(jj)->getCharge()) 
                 && upcEvt->getTrack(j)->getNhits()>25 && upcEvt->getTrack(jj)->getNhits()>25 
                 && (upcEvt->getTrack(j)->getFlag(StUPCTrack::kTof) 
                 && upcEvt->getTrack(jj)->getFlag(StUPCTrack::kTof)) 
                 && upcEvt->getTrack(jj)->getPt()>0.15 && upcEvt->getTrack(j)->getPt()>0.15 
                 && abs(upcEvt->getTrack(jj)->getEta())<1.0 &&  abs(upcEvt->getTrack(j)->getEta())<1.0 ) {

        TVector3 vertex(0,0,0);
        int id1,id2;
        TLorentzVector  v1,  v2;
        if(upcEvt->getTrack(j)->getCharge() > 0) {
          upcEvt->getTrack(j)->getLorentzVector(v1, massPion);
          upcEvt->getTrack(jj)->getLorentzVector(v2, massPion);
        } else {
          upcEvt->getTrack(j)->getLorentzVector(v2, massPion);
          upcEvt->getTrack(jj)->getLorentzVector(v1, massPion);
        }
//        continue;
//        if ( (v1+v2).M()>0.6 || (v1+v2).M()<0.4 ) continue;  
        StUPCV0 V0(upcEvt->getTrack(j),upcEvt->getTrack(jj), massProton, massPion,id1,id2, vertex, beamline, upcEvt->getMagneticField(), false);
////        continue;
        std::cout << "V0 " << V0.dcaDaughters() 
                  << " " << V0.DCABeamLine() 
  << " " <<sqrt( V0.prodVertexHypo().X()*V0.prodVertexHypo().X()+V0.prodVertexHypo().Y()*V0.prodVertexHypo().Y()) 
        <<  std::endl;
//        continue;
        TParticle V0Part;
        V0Part.SetMomentum(V0.px(),V0.py(),V0.pz(),sqrt(V0.m()*V0.m()+V0.px()*V0.px()+V0.py()*V0.py()+V0.pz()*V0.pz()));

        double Delta1 = 1000.;
        double Delta2 = 1000.;

        if ( TrueK0.size() == 2) {
        double pi=acos(-1.);
        double DeltaPhi1=V0Part.Phi()-TrueK0[0]->Phi();
       	double DeltaPhi2=V0Part.Phi()-TrueK0[1]->Phi();

        if ( DeltaPhi1 < -pi ) DeltaPhi1 =2*pi+DeltaPhi1;
       	if ( DeltaPhi1 > pi )  DeltaPhi1 =2*pi-DeltaPhi1;
       	if ( DeltaPhi2 < -pi ) DeltaPhi2 =2*pi+DeltaPhi2;
        if ( DeltaPhi2 > pi )  DeltaPhi2 =2*pi-DeltaPhi2;


        Delta1 = sqrt((V0Part.Eta()-TrueK0[0]->Eta())*(V0Part.Eta()-TrueK0[0]->Eta())+DeltaPhi1*DeltaPhi1);
        Delta2 = sqrt((V0Part.Eta()-TrueK0[1]->Eta())*(V0Part.Eta()-TrueK0[1]->Eta())+DeltaPhi2*DeltaPhi2);


        HistMatch->Fill(sqrt((V0Part.Eta()-TrueK0[0]->Eta())*(V0Part.Eta()-TrueK0[0]->Eta())+DeltaPhi1*DeltaPhi1));
        HistMatch->Fill(sqrt((V0Part.Eta()-TrueK0[1]->Eta())*(V0Part.Eta()-TrueK0[1]->Eta())+DeltaPhi2*DeltaPhi2));


        }
        bool MasCut = V0.m()<0.53 && V0.m()>0.46;
//          bool MasCut = V0.m()<1.15 && V0.m()>1.08;
//        MasCut = true;
        bool DCADCut = V0.dcaDaughters()<1.5;
        bool DCABLCut = V0.DCABeamLine()<1.5;
        bool RCut = sqrt(V0.v0x()*V0.v0x()+V0.v0y()*V0.v0y())>1.5;
        RCut = true; 
//        continue;
        bool PACut = V0.pointingAngleHypo()>0.925;
        TLorentzVector ProdVert0,  DecayVert0(0,0,0,0), ProdVert1,  DecayVert1(0,0,0,0) , RecK0, mom1, mom2;

        if( Delta1<0.05 || Delta2<0.05 ) { 
             if(Delta1<0.05) {
                    for (int ii=0; ii < PosPions.size(); ii++)
                         for (int jj=0; jj< NegPions.size(); jj++ ) {
                    PosPions[ii]->Momentum(mom1); 
                    NegPions[jj]->Momentum(mom2);
                    RecK0 = mom1 + mom2;
                       if ( abs(RecK0.Eta()-TrueK0[0]->Eta())<0.001 ) {
//                         cout << "Pos pions " << abs(mom1.Eta()-v1.Eta()) << " " << abs(mom2.Eta()-v2.Eta()) << endl; 
                        if( abs(mom1.Eta()-v1.Eta()) < 0.05 && abs(mom2.Eta()-v2.Eta()) < 0.05)  PosPions[ii]->ProductionVertex(DecayVert0);
                       }     
                    }
                   TrueK0[0]->ProductionVertex(ProdVert0);
                 ProdVert0.Print();
//                 DecayVert0.Print();
//                 std::cout << "0 Z,R,L : " << TrueK0[0]->Vz() << " " << sqrt(DecayVert0.X()*DecayVert0.X()+DecayVert0.Y()*DecayVert0.Y())
//                                                           << " " << (DecayVert0-ProdVert0).Vect().Mag() << std::endl;
             }   


             if(Delta2<0.05) {
                    for (int ii=0; ii < PosPions.size(); ii++)
                         for (int jj=0; jj< NegPions.size(); jj++ ) {
                    PosPions[ii]->Momentum(mom1);
                    NegPions[jj]->Momentum(mom2);
                    RecK0 = mom1 + mom2;
                       if ( abs(RecK0.Eta()-TrueK0[1]->Eta())<0.001 ) {
//                         cout << "Delta pions " << abs(mom1.Eta()-v1.Eta()) << " " << abs(mom2.Eta()-v2.Eta())<< endl;
                         if( abs(mom1.Eta()-v1.Eta()) < 0.05 && abs(mom2.Eta()-v2.Eta()) < 0.05) PosPions[ii]->ProductionVertex(DecayVert1);
                       }
                    }
                   TrueK0[1]->ProductionVertex(ProdVert1);
//                 ProdVert1.Print();
//                 DecayVert1.Print();
//                 std::cout << "1 Z,R,L : " << TrueK0[1]->Vz() << " " << sqrt(DecayVert1.X()*DecayVert1.X()+DecayVert1.Y()*DecayVert1.Y())
//                                                            << " " << (DecayVert1-ProdVert1).Vect().Mag() << std::endl;
             }


//        DecayVert0.Print();
//        DecayVert1.Print();
     
          if ( DecayVert0.Rho()!=0 || DecayVert1.Rho()!=0 ) {
            if( DCADCut&&DCABLCut&& RCut&&PACut ) HistV0MDetK0->Fill(V0.m());
            if( MasCut&&DCABLCut&& RCut&&PACut) HistV0DCADK0->Fill(V0.dcaDaughters());
            if( DCADCut&&MasCut&&DCABLCut&&PACut) HistV0RK0->Fill(sqrt(V0.v0x()*V0.v0x()+V0.v0y()*V0.v0y()));
            if( DCADCut&&MasCut&& RCut&&PACut) HistV0DCABeamLineK0->Fill(V0.DCABeamLine());
            if( DCADCut&&MasCut&&DCABLCut&& RCut) HistV0PointingAngleK0->Fill(V0.pointingAngleHypo());
            if( DCADCut&&MasCut&&DCABLCut&& RCut &&PACut) HistV0DecayLengthK0->Fill(V0.decayLengthHypo());
            if( DCADCut&&MasCut&&DCABLCut&& RCut &&PACut) {

              if( DecayVert0.Rho()!=0 ) {
              HistDeltaZ->Fill(ProdVert0.Vect().Z()-V0.prodVertexHypo().Z());
              HistDeltaX->Fill(ProdVert0.Vect().X()-V0.prodVertexHypo().X());
 
              HistDeltaR->Fill(sqrt( (V0.v0x()-DecayVert0.X())*(V0.v0x()-DecayVert0.X())+
                                     (V0.v0y()-DecayVert0.Y())*(V0.v0y()-DecayVert0.Y())+
                                     (V0.v0z()-DecayVert0.Z())*(V0.v0z()-DecayVert0.Z()))); 
              HistDeltaL->Fill( (DecayVert0-ProdVert0).Vect().Mag(),V0.decayLengthHypo() );
              }
       	      if( DecayVert1.Rho()!=0 )	{
              HistDeltaZ->Fill(ProdVert1.Vect().Z()-V0.prodVertexHypo().Z());
              HistDeltaX->Fill(ProdVert1.Vect().X()-V0.prodVertexHypo().X());

              HistDeltaR->Fill(sqrt( (V0.v0x()-DecayVert1.X())*(V0.v0x()-DecayVert1.X())+  
       	       	       	       	     (V0.v0y()-DecayVert1.Y())*(V0.v0y()-DecayVert1.Y())+
       	       	       	       	     (V0.v0z()-DecayVert1.Z())*(V0.v0z()-DecayVert1.Z())));
              HistDeltaL->Fill( (DecayVert1-ProdVert1).Vect().Mag(),V0.decayLengthHypo() );
              }
         }
        }
        } else {
         if( DCADCut&&DCABLCut&& RCut &&PACut) { 
              HistV0MDet->Fill(V0.m());
              HistV0MvsDecayLDet->Fill(V0.decayLengthHypo(),V0.m());
         }
         if( MasCut&&DCABLCut&& RCut &&PACut) HistV0DCAD->Fill(V0.dcaDaughters());
         if( DCADCut&&MasCut&&DCABLCut &&PACut) HistV0R->Fill(sqrt(V0.v0x()*V0.v0x()+V0.v0y()*V0.v0y()));
         if( DCADCut&&MasCut&& RCut&&PACut) HistV0DCABeamLine->Fill(V0.DCABeamLine());
         if( DCADCut&&MasCut&&DCABLCut&& RCut) HistV0PointingAngle->Fill(V0.pointingAngleHypo());
         if( DCADCut&&MasCut&&DCABLCut&& RCut&&PACut) HistV0DecayLength->Fill(V0.decayLengthHypo());
      
        }
        if ( V0.dcaDaughters()>1.5 ) {
         HistV0MDetRejected->Fill(V0.m());
        }
        if ( V0.DCABeamLine()>1.0 ) {
         HistV0MDetRejectedDCA->Fill(V0.m());
        }
       }
      }
     }
      std::cout << "event: " << i << std::endl;
//      continue;
    


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
            continue;
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
            StUPCV0 V0(tpcTrack[0],tpcTrack[1], massPion, massPion,id1,id2, vertex, beamline, upcEvt->getMagneticField(), false);
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
    HistV0MDetK0->Write();
    HistV0DCAD->Write();
    HistV0R->Write();
    HistV0MDetRejected->Write();
    HistV0MDetRejectedDCA->Write();
    HistV0DCADK0->Write();
    HistV0DCABeamLineK0->Write();
    HistV0DCABeamLine->Write();
    HistV0PointingAngle->Write();
    HistV0PointingAngleK0->Write();
    HistV0DecayLength->Write();
    HistV0DecayLengthK0->Write();
    HistMatch->Write();
    HistV0RK0->Write();
    HistPYPY->Write();
    HistDeltaZ->Write();
    HistDeltaX->Write();
    HistDeltaL->Write();
    HistDeltaR->Write();
    HistV0MvsDecayLDet->Write();

    outfile->Close();

    return 0;
}
