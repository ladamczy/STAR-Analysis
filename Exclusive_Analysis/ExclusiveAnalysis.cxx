#include "Includes.h"

using namespace std;

int main(int argc, char** argv)  
{

    TH1D* HistCutFlow = new TH1D("HistCutFlow", ";;count", 16, -0.5, 15.5);
    TProfile* ProfileDistVertexBeamX = new TProfile("ProfileDistVertexBeamX","",639, 20511.5, 21150.5);
    TProfile* ProfileDistVertexBeamY = new TProfile("ProfileDistVertexBeamY","",639, 20511.5, 21150.5);

    TH1D* HistVtxX = new TH1D("HistVtxX", "; vtx_{x} [cm]; # events", 100, -2, 2);
    TH1D* HistVtxY = new TH1D("HistVtxY"," ; vtx_{y} [cm]; # events", 100, -2, 2);
    TH1D* HistVtxZ = new TH1D("HistVtxZ", "; vtx_{z} [cm]; # events", 100, -200, 200);
    TH1D* HistVtxErrX = new TH1D("HistVtxErrX", ";err vtx_{x} [cm]; # events", 100, 0.0, 1);
    TH1D* HistVtxErrY = new TH1D("HistVtxErrY"," ;err vtx_{y} [cm]; # events", 100, 0.0, 1);
    TH1D* HistVtxErrZ = new TH1D("HistVtxErrZ", ";err vtx_{z} [cm]; # events", 100, 0.0, 0.5);
    TH1D* HistR = new TH1D("HistR", ";R [cm]; # events", 150, 0, 0.3);
    TH1D* HistSig = new TH1D("HistSig"," ; significance; # events", 250, 0.0, 10);

    TH1D* HistNumPrimaryVertices = new TH1D("5_HistNumPrimaryVertices", " ;Number of primary vertices; count", 11, -0.5, 10.5); 
    TH1D* HistPrimaryVertexAbsPosZ = new TH1D("6_HistPrimaryVertexAbsPosZ", ";primary vertex z [cm]; count", 100, 0, 200);
    TH1D* HistNumTofMatchedTracks = new TH1D("7_HistNumTofMatchedTracks", ";number of TOF matched tracks; count", 16, -0.5, 15.5);    
    TH1D* HistTofMatchedTracksCharge = new TH1D("8_HistTofMatchedTracksCharge", ";charge [e]; count", 3, -1.5, 1.5);
    TH1D* HistTofMatchedTracksAbsEta = new TH1D("9_HistTofMatchedTracksAbsEta", ";#eta; count", 50, 0, 2);
    TH1D* HistTofMatchedTracksAbsDcaZ = new TH1D("10_HistTofMatchedTracksAbsDcaZ", ";DCA_{Z} [cm]; count ", 60, 0, 4);
    TH1D* HistTofMatchedTracksDcaXY = new TH1D("11_HistTofMatchedTracksDcaXY", ";DCA_{XY} [cm]; count", 60, 0, 4);
    TH1D* HistTofMatchedTracksNfit = new TH1D("12_HistTofMatchedTracksNfit", ";N_{fit}; count", 61, -0.5, 60.5);
    TH1D* HistTofMatchedTracksNdEdx = new TH1D("13_HistTofMatchedTracksNdEdx", ";N_{dE/dx}; count", 61, -0.5, 60.5);

    TH1D* HistNumPrimaryVerticesPostSelection = new TH1D("5_HistNumPrimaryVerticesPs", " ;Number of primary vertices; count", 11, -0.5, 10.5); 
    TH1D* HistPrimaryVertexAbsPosZPostSelection = new TH1D("6_HistPrimaryVertexAbsPosZPs", ";primary vertex z [cm]; count", 100, 0, 200);
    TH1D* HistNumTofMatchedTracksPostSelection = new TH1D("7_HistNumTofMatchedTracksPs", ";number of TOF matched tracks; count", 16, -0.5, 15.5);    
    TH1D* HistTofMatchedTracksChargePostSelection = new TH1D("8_HistTofMatchedTracksChargePs", ";charge [e]; count", 3, -1.5, 1.5);
    TH1D* HistTofMatchedTracksAbsEtaPostSelection = new TH1D("9_HistTofMatchedTracksAbsEtaPs", ";#eta; count", 50, 0, 2);
    TH1D* HistTofMatchedTracksAbsDcaZPostSelection = new TH1D("10_HistTofMatchedTracksAbsDcaZPs", ";DCA_{Z} [cm]; count ", 60, 0, 4);
    TH1D* HistTofMatchedTracksDcaXYPostSelection = new TH1D("11_HistTofMatchedTracksDcaXYPs", ";DCA_{XY} [cm]; count", 60, 0, 4);
    TH1D* HistTofMatchedTracksNfitPostSelection = new TH1D("12_HistTofMatchedTracksNfitPs", ";N_{fit}; count", 61, -0.5, 60.5);
    TH1D* HistTofMatchedTracksNdEdxPostSelection = new TH1D("13_HistTofMatchedTracksNdEdxPs", ";N_{dE/dx}; count", 61, -0.5, 60.5);

    TH2D* HistPtProtonPtKaonsX = new TH2D("HistPtProtonPtKaonsX", ";[p_{K^{0}}p_{K^{0}}]_{x}; [p_{p'}^{W}p_{p'}^{E}]_{x} ",  100, -2.0, 2.0, 100, -2.0, 2.0);
    TH2D* HistPtProtonPtKaonsY = new TH2D("HistPtProtonPtKaonsY", ";[p_{K^{0}}p_{K^{0}}]_{y}; [p_{p'}^{W}p_{p'}^{E}]_{y} ",  100, -2.0, 2.0, 100, -2.0, 2.0);
    TH1D* HistPtProtonPtKaonsX1d = new TH1D("HistPtProtonPtKaonsX1d", ";[p_{K^{0}}p_{K^{0}} + [p_{p'}^{W}p_{p'}^{E}]_{x}; events ",  200, -5.0, 5.0);
    TH1D* HistPtProtonPtKaonsY1d = new TH1D("HistPtProtonPtKaonsY1d", ";[p_{K^{0}}p_{K^{0}} +  [p_{p'}^{W}p_{p'}^{E}]_{y}; events ",  200, -5.0, 5.0);
    TH2D* HistPtProtonPtKaonsXOppY = new TH2D("HistPtProtonPtKaonsXOppY", ";[p_{K^{0}}p_{K^{0}}]_{x}; [p_{p'}^{W}p_{p'}^{E}]_{x} ",  100, -2.0, 2.0, 100, -2.0, 2.0);
    TH2D* HistPtProtonPtKaonsYOppY = new TH2D("HistPtProtonPtKaonsYOppY", ";[p_{K^{0}}p_{K^{0}}]_{y}; [p_{p'}^{W}p_{p'}^{E}]_{y} ",  100, -2.0, 2.0, 100, -2.0, 2.0);
    TH1D* HistPtProtonPtKaonsX1dOppY = new TH1D("HistPtProtonPtKaonsX1dOppY", ";[p_{K^{0}}p_{K^{0}} + [p_{p'}^{W}p_{p'}^{E}]_{x}; events ",  200, -5.0, 5.0);
    TH1D* HistPtProtonPtKaonsY1dOppY = new TH1D("HistPtProtonPtKaonsY1dOppY", ";[p_{K^{0}}p_{K^{0}} +  [p_{p'}^{W}p_{p'}^{E}]_{y}; events ",  200, -5.0, 5.0);
    TH2D* HistPtProtonPtKaonsXSameY = new TH2D("HistPtProtonPtKaonsXSameY", ";[p_{K^{0}}p_{K^{0}}]_{x}; [p_{p'}^{W}p_{p'}^{E}]_{x} ",  100, -2.0, 2.0, 100, -2.0, 2.0);
    TH2D* HistPtProtonPtKaonsYSameY = new TH2D("HistPtProtonPtKaonsYSameY", ";[p_{K^{0}}p_{K^{0}}]_{y}; [p_{p'}^{W}p_{p'}^{E}]_{y} ",  100, -2.0, 2.0, 100, -2.0, 2.0);
    TH1D* HistPtProtonPtKaonsX1dSameY = new TH1D("HistPtProtonPtKaonsX1dSameY", ";[p_{K^{0}}p_{K^{0}} + [p_{p'}^{W}p_{p'}^{E}]_{x}; events ",  200, -5.0, 5.0);
    TH1D* HistPtProtonPtKaonsY1dSameY = new TH1D("HistPtProtonPtKaonsY1dSameY", ";[p_{K^{0}}p_{K^{0}} +  [p_{p'}^{W}p_{p'}^{E}]_{y}; events ",  200, -5.0, 5.0);

    TH1D* HistInvMassPiPi = new TH1D("HistInvMassPiPi", "; m_{#pi^{+}#pi^{-}} [GeV]; # events", 250 ,0, 3);
    TH1D* HistInvMassPiPiPeak = new TH1D("HistInvMassPiPiPeak", "; m_{#pi^{+}#pi^{-}} [GeV]; # events", 100 ,0.4, 0.6);
    TH1D* HistInvMassPiPi2 = new TH1D("HistInvMassPiPi2", "; m_{#pi^{+}#pi^{-}} [GeV]; # events", 250 ,0, 3);
    TH1D* HistInvMassPiPiPeak2 = new TH1D("HistInvMassPiPiPeak2", "; m_{#pi^{+}#pi^{-}} [GeV]; # events", 100 ,0.4, 0.6);

    TH1D* HistInvMassPiPiPostSelection = new TH1D("HistInvMassPiPiPs", "; m_{#pi^{+}#pi^{-}} [GeV]; # events", 250 ,0, 3);
    TH1D* HistInvMassPiPiPeakPostSelection = new TH1D("HistInvMassPiPiPeakPs", "; m_{#pi^{+}#pi^{-}} [GeV]; # events", 100 ,0.4, 0.6);
    TH1D* HistInvMassPiPi2PostSelection = new TH1D("HistInvMassPiPi2Ps", "; m_{#pi^{+}#pi^{-}} [GeV]; # events", 250 ,0, 3);
    TH1D* HistInvMassPiPiPeak2PostSelection = new TH1D("HistInvMassPiPiPeak2Ps", "; m_{#pi^{+}#pi^{-}} [GeV]; # events", 100 ,0.4, 0.6);

    TH1D* HistKaonPtTruth = new TH1D("HistKaonPtTruth", "; pT_{K^{0}} [GeV]; # events", 25 ,0, 2.5);
    TH1D* HistKaonEtaTruth = new TH1D("HistKaonEtaTruth", "; #eta_{K^{0}}; # events", 25 ,-4.5, 4.5);
    TH1D* HistKaonVtxZTruth = new TH1D("HistKaonVtxZTruth", ";z_{vtx}^{#pi^{+}#pi^{-}} [cm]; # events", 10, -150, 150);
    TH1D* HistKaonVtxRTruth = new TH1D("HistKaonVtxRTruth", ";R_{vtx}^{#pi^{+}#pi^{-}} [cm]; # events", 10, 0, 20);

    TH1D* HistKaonPtDet = new TH1D("HistKaonPtDet", "; pT_{K^{0}} [GeV]; # events", 25 ,0, 2.5);
    TH1D* HistKaonEtaDet = new TH1D("HistKaonEtaDet", "; #eta_{K^{0}}; # events", 25 ,-4.5, 4.5);
    TH1D* HistKaonVtxZDet = new TH1D("HistKaonVtxZDet", "; z_{vtx}^{#pi^{+}#pi^{-}} [cm]; # events", 10, -150, 150);
    TH1D* HistKaonVtxRDet = new TH1D("HistKaonVtxRDet", "; R_{vtx}^{#pi^{+}#pi^{-}} [cm]; # events", 10, 0, 20);


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
    while (std::getline(inputFilePathList, inputFileName))
    {
        chain->Add(inputFileName.c_str());
    }
    inputFilePathList.close();

    bool isMC = 0;
    static StUPCEvent *upcEvt = 0x0;
    static StRPEvent *rpEvt = 0x0;
    chain->SetBranchAddress("mUPCEvent", &upcEvt);
    chain->GetEntry(0);
    if (upcEvt->getRunNumber() == 1)
    {
        isMC = 1;
    }

    if (isMC == 0)
    {
        chain->SetBranchAddress("mRPEvent", &rpEvt);
    }

    vector <TLorentzVector> posPion;
    vector <TLorentzVector> negPion;
    TLorentzVector trackVector;
	TVector3 proton1, proton2;
    TParticle* particle;

    vector <TParticle*> vPosPionsMC, vNegPionsMC, vPosPionsMCPaired, vNegPionsMCPaired;

    double beamPositionX, beamPositionY;
    double primVertexPosX, primVertexPosY, primVertexPosZ, primVertexPosErrX, primVertexPosErrY, primVertexPosErrZ;
    double distVertexBeamX, distVertexBeamY, distR, significance;
    int iCutFlow;

    const vector <int> triggerID = {570209, 570219, 570229, 570701, 570702, 570703, 570704, 570705, 570709, 570711, 570712, 570719, 590701, 590703, 590705, 590708, 590709};
    const vector <int> triggerCEP = {570701, 570705, 570711, 590701, 590705, 590708};


    // to be changed... there is no warning in case of the invalid input file
	vector <vector<double>> fillNumberWithPosition = ReadFillPositionData("InputData/Run7PolarizationWithPosition.csv");

    for (Long64_t i = 0; i < chain->GetEntries(); ++i) 
    {
        if (i%1000000 == 0)  
        {
            cout << i << "/" <<  chain->GetEntries() << endl;
        }

        chain->GetEntry(i);

        iCutFlow = 0;
        HistCutFlow->Fill(iCutFlow);  // cutflow: all data

        double protonSignsY = 1;
        vector <TParticle*> vProtons;
        if (isMC == 0)
        {
            for(unsigned int k = 0; k < rpEvt->getNumberOfTracks(); ++k)
            {		
                StUPCRpsTrack *trk = rpEvt->getTrack(k);
                trk->setEvent(rpEvt);
                if (k == 0)
                {
                    proton1 = trk->pVec();
                    protonSignsY *= trk->pVec().Y();
                }
                else if (k == 1)
                {
                    proton2 = trk->pVec();
                    protonSignsY *= trk->pVec().Y();
                }
            
            }    
        }


        else if (isMC == 1)
        {
            for (int i = 0; i < upcEvt->getNumberOfMCParticles(); i++)
            {
                particle = upcEvt->getMCParticle(i);
                if (particle->GetPDG()->PdgCode() == 2212 and particle->GetFirstMother() == 1)
                {
                    vProtons.push_back(particle);
                }
            }

            proton1.SetX(vProtons[0]->Px());
            proton1.SetY(vProtons[0]->Py());
            proton1.SetZ(vProtons[0]->Pz());

            proton2.SetX(vProtons[1]->Px());
            proton2.SetY(vProtons[1]->Py());
            proton2.SetZ(vProtons[1]->Pz());
            protonSignsY = protonSignsY*vProtons[0]->Py()*vProtons[1]->Py();
        }

        // number of primary vertices - fos single Kaon we demand exacly one primary vertex
        Int_t numberOfPrimaryVertices = upcEvt->getNumberOfVertices();       	
        HistNumPrimaryVertices->Fill(numberOfPrimaryVertices);		
        if(numberOfPrimaryVertices!=1)
        {
            continue;
        }		
        iCutFlow+=1;
		HistCutFlow->Fill(iCutFlow); // cutflow: one vertex
        HistNumPrimaryVerticesPostSelection->Fill(numberOfPrimaryVertices);	

        int nFillNumber = upcEvt->getFillNumber();
        if (isMC == 0)
        {
		    beamPositionX = FindPosition(nFillNumber, fillNumberWithPosition[0], fillNumberWithPosition[1], fillNumberWithPosition[2], fillNumberWithPosition[3], fillNumberWithPosition[4])[0];
		    beamPositionY = FindPosition(nFillNumber, fillNumberWithPosition[0], fillNumberWithPosition[1], fillNumberWithPosition[2], fillNumberWithPosition[3], fillNumberWithPosition[4])[1];			  
        }

        else if (isMC == 1)
        {
            beamPositionX = 0.0;
            beamPositionY = 0.0;
        }

        // profile - distance between vertex and beam in x-y plane, parameter R and significance
        primVertexPosX = upcEvt->getVertex(0)->getPosX();
        primVertexPosY = upcEvt->getVertex(0)->getPosY();
        primVertexPosZ = upcEvt->getVertex(0)->getPosZ();

        primVertexPosErrX = upcEvt->getVertex(0)->getErrX();
        primVertexPosErrY = upcEvt->getVertex(0)->getErrY();
        primVertexPosErrZ = upcEvt->getVertex(0)->getErrZ();

        distVertexBeamX = primVertexPosX - beamPositionX;
        distVertexBeamY = primVertexPosY - beamPositionY;
        
        distR = sqrt(pow(distVertexBeamX,2) + pow(distVertexBeamY,2));
   		significance = sqrt(pow((distVertexBeamX/primVertexPosErrX),2) + pow((distVertexBeamY/primVertexPosErrY),2));       

        // fill number and beam position
        ProfileDistVertexBeamX->Fill(nFillNumber, distVertexBeamX); 
        ProfileDistVertexBeamY->Fill(nFillNumber, distVertexBeamY); 
        HistR->Fill(distR);
   		HistSig->Fill(significance);    

        HistVtxX->Fill(distVertexBeamX);
        HistVtxY->Fill(distVertexBeamY);
        HistVtxZ->Fill(primVertexPosZ);
        
        HistVtxErrX->Fill(primVertexPosErrX);
        HistVtxErrY->Fill(primVertexPosErrY);
        HistVtxErrZ->Fill(primVertexPosErrZ);


		//primary vertex is placed within |zvtx| < 80 cm
        HistPrimaryVertexAbsPosZ->Fill(abs(primVertexPosZ));
		if(abs(upcEvt->getVertex(0)->getPosZ())>=80)
		{
			continue;
        }   
        iCutFlow+=1;
		HistCutFlow->Fill(iCutFlow); // cutflow: primary vertex |zvtx| < 80 cm
        HistPrimaryVertexAbsPosZPostSelection->Fill(abs(upcEvt->getVertex(0)->getPosZ()));      

		vector <StUPCTrack*> tracksWithTofHit;
        for (Int_t i = 0; i<upcEvt->getNumberOfTracks(); i++)
		{
			if(upcEvt->getTrack(i)->getFlag(StUPCTrack::kTof))
			{
                tracksWithTofHit.push_back(upcEvt->getTrack(i));
			}
		}       		

		// two TOF-matched primary tracks
        HistNumTofMatchedTracks->Fill(tracksWithTofHit.size());
		bool isValidNumberOfTofMatchedTracks = 0;
		if (tracksWithTofHit.size() == 4)
		{
			isValidNumberOfTofMatchedTracks = 1;
		}
		if (isValidNumberOfTofMatchedTracks == 0) 
        {
            tracksWithTofHit.clear();
            continue;
        } 
        iCutFlow+=1;
		HistCutFlow->Fill(iCutFlow); // cutflow: two TOF matched tracks        
        HistNumTofMatchedTracksPostSelection->Fill(tracksWithTofHit.size());

		// TOF-matched tracks are of opposied signs		
        int sumCharge = 0;
        for (unsigned long int i = 0; i < tracksWithTofHit.size(); i++)
        {
            HistTofMatchedTracksCharge->Fill(tracksWithTofHit[i]->getCharge());
            sumCharge+=tracksWithTofHit[i]->getCharge();
        
        }
        if (sumCharge != 0)  
        {
            tracksWithTofHit.clear();
            continue;
        }
        iCutFlow+=1;
		HistCutFlow->Fill(iCutFlow);  // cutflow: sum of the charge is zero
    
        sumCharge = 0;
        for (unsigned long int i = 0; i < tracksWithTofHit.size(); i++)
        {
            HistTofMatchedTracksChargePostSelection->Fill(tracksWithTofHit[i]->getCharge());        
        }

        bool isEtaValid = 1;
        for (unsigned long int i = 0; i < tracksWithTofHit.size(); i++)
        {	
            HistTofMatchedTracksAbsEta->Fill(abs(tracksWithTofHit[i]->getEta()));
            if (abs(tracksWithTofHit[i]->getEta()) >= 0.7)
            {
                isEtaValid = 0;
            }
        }
        if (isEtaValid == 0) {continue;}
        iCutFlow+=1;
		HistCutFlow->Fill(iCutFlow);  // cutflow: TOF-matched tracks within appropraite pseudorapidity range
      
        for (unsigned long int i = 0; i < tracksWithTofHit.size(); i++)
        {	
            HistTofMatchedTracksAbsEtaPostSelection->Fill(abs(tracksWithTofHit[i]->getEta()));
        }
        
		//  associated global tracks match well to the primary vertex |DCA(z)| < 2.0 // ROZSZERZENIE ZAKRESU Z 1.0		
        bool isDcazSmall = 1;
        for (unsigned long int i = 0; i < tracksWithTofHit.size(); i++)
        {	
            HistTofMatchedTracksAbsDcaZ->Fill(abs(tracksWithTofHit[i]->getDcaZ()));	
            if (abs(tracksWithTofHit[i]->getDcaZ())>=2.0)
            {
                isDcazSmall = 0;
            }
        }
        if (isDcazSmall == 0) {continue;}
        iCutFlow+=1;
		HistCutFlow->Fill(iCutFlow);  // cutflow: TOF-matched tracks match primary vertex

        for (unsigned long int i = 0; i < tracksWithTofHit.size(); i++)
        {	
            HistTofMatchedTracksAbsDcaZPostSelection->Fill(abs(tracksWithTofHit[i]->getDcaZ()));	
        }


		//  associated global tracks match well to the primary vertex DCA(R) < 3.0 //ROZSZERZENIE ZAKRESU Z 1.5
        bool isDcaxySmall = 1;
        for (unsigned long int i = 0; i < tracksWithTofHit.size(); i++)
        {	
            HistTofMatchedTracksDcaXY->Fill(tracksWithTofHit[i]->getDcaXY());
            if (tracksWithTofHit[i]->getDcaXY() >= 3.0)
            {
                isDcaxySmall = 0;
            }
        }
        if (isDcaxySmall == 0) {continue;}
        iCutFlow+=1;
		HistCutFlow->Fill(iCutFlow);   //  cutflow: TOF-matched tracks match primary vertex

        for (unsigned long int i = 0; i < tracksWithTofHit.size(); i++)
        {	
            HistTofMatchedTracksDcaXYPostSelection->Fill(tracksWithTofHit[i]->getDcaXY());
        }

		// associated global tracks satisfy quality criteria N_{fit} >= 25
        bool isNfitAtLeast25 = 1;
        for (unsigned long int i = 0; i < tracksWithTofHit.size(); i++)
        {		
            HistTofMatchedTracksNfit->Fill(tracksWithTofHit[i]->getNhitsFit());
            if (tracksWithTofHit[i]->getNhitsFit() < 25)
            {
                isNfitAtLeast25 = 0;
            }
        }
        if (isNfitAtLeast25 == 0) {continue;}
        iCutFlow+=1;
		HistCutFlow->Fill(iCutFlow);    // cutflow: two tracks quality N_{fit} >= 25

        for (unsigned long int i = 0; i < tracksWithTofHit.size(); i++)
        {		
            HistTofMatchedTracksNfitPostSelection->Fill(tracksWithTofHit[i]->getNhitsFit());
        }


		// associated global tracks satisfy quality criteria N_{dE/dx} >= 15 
        bool isNdEdxAtLeast15 = 1;
        for (unsigned long int i = 0; i < tracksWithTofHit.size(); i++)
        {		
            HistTofMatchedTracksNdEdx->Fill(tracksWithTofHit[i]->getNhitsDEdx());
            if (tracksWithTofHit[i]->getNhitsDEdx() < 15)
            {
                isNdEdxAtLeast15 = 0;
            }
        }
        if (isNdEdxAtLeast15 == 0) {continue;}
        iCutFlow+=1;
		HistCutFlow->Fill(iCutFlow);  // cutflow: two tracks quality N_{dE/dx} >= 15 

        for (unsigned long int i = 0; i < tracksWithTofHit.size(); i++)
        {		
            HistTofMatchedTracksNdEdxPostSelection->Fill(tracksWithTofHit[i]->getNhitsDEdx());
        }
	    
        for (unsigned long int i = 0; i < tracksWithTofHit.size(); i++)
        {
            tracksWithTofHit[i]->getLorentzVector(trackVector, massPion);
            if( tracksWithTofHit[i]->getCharge()==1)
            {
                posPion.push_back(trackVector);
            }
            else if (tracksWithTofHit[i]->getCharge()==-1)
            {
                negPion.push_back(trackVector);
            }
        }
          
        if(posPion.size() == 0 or negPion.size() == 0)
        {
                posPion.clear();
                negPion.clear();
                continue;
        }

        TLorentzVector lVecKaonComb1a = {0, 0, 0, 0};
        TLorentzVector lVecKaonComb1b = {0, 0, 0, 0};

        TLorentzVector lVecKaonComb2a = {0, 0, 0, 0};
        TLorentzVector lVecKaonComb2b = {0, 0, 0, 0};

        // two combinations
        lVecKaonComb1a = posPion[0] + negPion[0];
        lVecKaonComb1b = posPion[1] + negPion[1];

        lVecKaonComb2a = posPion[0] + negPion[1];
        lVecKaonComb2b = posPion[1] + negPion[0];

        double kaonMassWindowLow = 0.46;
        double kaonMassWindowHigh = 0.52;

        bool arePionPairsInMassWindowFirstCombination = 0;
        bool isPionPairInMassWindowFirstCombinationA = 0;
        bool isPionPairInMassWindowFirstCombinationB = 0;

        bool arePionPairsInMassWindowSecondCombination = 0;
        bool isPionPairInMassWindowSecondCombinationA = 0;
        bool isPionPairInMassWindowSecondCombinationB = 0;

        //check which paris are within mass window for the first combination
        if ((lVecKaonComb1a.M() > kaonMassWindowLow and lVecKaonComb1a.M() < kaonMassWindowHigh) and  (lVecKaonComb1b.M() > kaonMassWindowLow and lVecKaonComb1b.M() < kaonMassWindowHigh))
        {
            arePionPairsInMassWindowFirstCombination = 1;
        }

        else if (lVecKaonComb1a.M() > kaonMassWindowLow and lVecKaonComb1a.M() < kaonMassWindowHigh)
        {
            isPionPairInMassWindowFirstCombinationA = 1;
        }

        else if (lVecKaonComb1b.M() > kaonMassWindowLow and lVecKaonComb1b.M() < kaonMassWindowHigh)
        {
            isPionPairInMassWindowFirstCombinationB = 1;
        }

        //check which paris are within mass window for the second combination  
        if ((lVecKaonComb2a.M() > kaonMassWindowLow and lVecKaonComb2a.M() < kaonMassWindowHigh) and  (lVecKaonComb2b.M() > kaonMassWindowLow and lVecKaonComb2b.M() < kaonMassWindowHigh))
        {
            arePionPairsInMassWindowSecondCombination = 1;
        }

        else if (lVecKaonComb2a.M() > kaonMassWindowLow and lVecKaonComb2a.M() < kaonMassWindowHigh)
        {
            isPionPairInMassWindowSecondCombinationA = 1;
        }

        else if (lVecKaonComb2b.M() > kaonMassWindowLow and lVecKaonComb2b.M() < kaonMassWindowHigh)
        {
            isPionPairInMassWindowSecondCombinationB = 1;
        }

        vector <TLorentzVector> lVecCombinations;
        vector <double> pTCombinations;
      
        int maxPtIndex;
        TLorentzVector conditioningKaon;
        TLorentzVector complementaryKaon;

        if (arePionPairsInMassWindowFirstCombination and arePionPairsInMassWindowSecondCombination)
        {
  
            lVecCombinations =  {lVecKaonComb1a, lVecKaonComb1b, lVecKaonComb2a, lVecKaonComb2b}; 
            pTCombinations = {lVecKaonComb1a.Pt(), lVecKaonComb1b.Pt(), lVecKaonComb2a.Pt(), lVecKaonComb2b.Pt()}; 
            maxPtIndex = distance(pTCombinations.begin(), max_element(pTCombinations.begin(), pTCombinations.end()));

            if (maxPtIndex <= 1)
            {
                conditioningKaon = lVecCombinations[maxPtIndex]; 
                complementaryKaon = lVecCombinations[abs(maxPtIndex-1)];
            }

            else
            {  
                conditioningKaon = lVecCombinations[maxPtIndex]; 
                complementaryKaon = lVecCombinations[2 + abs(3-maxPtIndex)]; 
            }
            HistInvMassPiPi->Fill(conditioningKaon.M());
            HistInvMassPiPiPeak->Fill(conditioningKaon.M());
            HistInvMassPiPi2->Fill(complementaryKaon.M());
            HistInvMassPiPiPeak2->Fill(complementaryKaon.M());

        }

        else if (arePionPairsInMassWindowFirstCombination)
        {

            if (lVecKaonComb1a.Pt() > lVecKaonComb1b.Pt())
            {
                conditioningKaon = lVecKaonComb1a;
                complementaryKaon = lVecKaonComb1b;
            }

            else
            {
                conditioningKaon = lVecKaonComb1b;
                complementaryKaon = lVecKaonComb1a;
            }
            HistInvMassPiPi->Fill(conditioningKaon.M());
            HistInvMassPiPiPeak->Fill(conditioningKaon.M());
            HistInvMassPiPi2->Fill(complementaryKaon.M());
            HistInvMassPiPiPeak2->Fill(complementaryKaon.M());
        }


        else if (arePionPairsInMassWindowSecondCombination)
        {

            if (lVecKaonComb2a.Pt() > lVecKaonComb2b.Pt())
            {
                conditioningKaon = lVecKaonComb2a;
                complementaryKaon = lVecKaonComb2b;
            }

            else
            {
                conditioningKaon = lVecKaonComb2b;
                complementaryKaon = lVecKaonComb2a;
            }
            HistInvMassPiPi->Fill(conditioningKaon.M());
            HistInvMassPiPiPeak->Fill(conditioningKaon.M());
            HistInvMassPiPi2->Fill(complementaryKaon.M());
            HistInvMassPiPiPeak2->Fill(complementaryKaon.M());

        }       

        else 
        {

            lVecCombinations =  {lVecKaonComb1a, lVecKaonComb1b, lVecKaonComb2a, lVecKaonComb2b}; 
            pTCombinations = {lVecKaonComb1a.Pt(), lVecKaonComb1b.Pt(), lVecKaonComb2a.Pt(), lVecKaonComb2b.Pt()}; 
            maxPtIndex = distance(pTCombinations.begin(), max_element(pTCombinations.begin(), pTCombinations.end()));

            if (maxPtIndex <= 1)
            {
                conditioningKaon = lVecCombinations[maxPtIndex]; 
                complementaryKaon = lVecCombinations[abs(maxPtIndex-1)];
            }

            else
            {  
                conditioningKaon = lVecCombinations[maxPtIndex]; 
                complementaryKaon = lVecCombinations[2 + abs(3-maxPtIndex)]; 
            }

            HistInvMassPiPi->Fill(conditioningKaon.M());
            HistInvMassPiPiPeak->Fill(conditioningKaon.M());
            HistInvMassPiPi2->Fill(complementaryKaon.M());
            HistInvMassPiPiPeak2->Fill(complementaryKaon.M());

            posPion.clear();
            negPion.clear();  
            pTCombinations.clear();
            lVecCombinations.clear();
            continue; 
        }
    
        iCutFlow+=1;
		HistCutFlow->Fill(iCutFlow);  

        HistInvMassPiPiPostSelection->Fill(conditioningKaon.M());
        HistInvMassPiPiPeakPostSelection->Fill(conditioningKaon.M());
        HistInvMassPiPi2PostSelection->Fill(complementaryKaon.M());
        HistInvMassPiPiPeak2PostSelection->Fill(complementaryKaon.M());

        HistPtProtonPtKaonsX->Fill(conditioningKaon.Px()+complementaryKaon.Px(), proton1.X() +  proton2.X());
        HistPtProtonPtKaonsY->Fill(conditioningKaon.Py()+complementaryKaon.Py(), proton1.Y() +  proton2.Y());
        HistPtProtonPtKaonsX1d->Fill(conditioningKaon.Px()+complementaryKaon.Px() + proton1.X() +  proton2.X());
        HistPtProtonPtKaonsY1d->Fill(conditioningKaon.Py()+complementaryKaon.Py() +proton1.Y() +  proton2.Y());

        if (protonSignsY < 0.0)
        {
            HistPtProtonPtKaonsXOppY->Fill(conditioningKaon.Px()+complementaryKaon.Px(), proton1.X() +  proton2.X());
            HistPtProtonPtKaonsYOppY->Fill(conditioningKaon.Py()+complementaryKaon.Py(), proton1.Y() +  proton2.Y());
            HistPtProtonPtKaonsX1dOppY->Fill(conditioningKaon.Px()+complementaryKaon.Px() + proton1.X() +  proton2.X());
            HistPtProtonPtKaonsY1dOppY->Fill(conditioningKaon.Py()+complementaryKaon.Py() +proton1.Y() +  proton2.Y());
        }

        else
        {
            HistPtProtonPtKaonsXSameY->Fill(conditioningKaon.Px()+complementaryKaon.Px(), proton1.X() +  proton2.X());
            HistPtProtonPtKaonsYSameY->Fill(conditioningKaon.Py()+complementaryKaon.Py(), proton1.Y() +  proton2.Y());
            HistPtProtonPtKaonsX1dSameY->Fill(conditioningKaon.Px()+complementaryKaon.Px() + proton1.X() +  proton2.X());
            HistPtProtonPtKaonsY1dSameY->Fill(conditioningKaon.Py()+complementaryKaon.Py() +proton1.Y() +  proton2.Y());
        }


        posPion.clear();
        negPion.clear();  
        pTCombinations.clear();
        lVecCombinations.clear();



        if (isMC == 1) // begin Monte Carlo efficiency
        {
            TLorentzVector lorentzVectorPosPion, lorentzVectorNegPion;
            TLorentzVector lorentzVectorNeutKaon;
            TLorentzVector trackVector;
            double truthVertexR, truthVertexZ, truthEta, truthPt;
            double detVertexR, detVertexZ, detEta, detPt;

            
            // extract all Pi+, Pi- and diffractive protons
            for (int i = 0; i < upcEvt->getNumberOfMCParticles(); i++)
            {
                particle = upcEvt->getMCParticle(i);
            
                if (particle->GetPDG()->PdgCode() == 211)
                {
                    vPosPionsMC.push_back(particle);
                }

                else if (particle->GetPDG()->PdgCode() == -211)
                {
                    vNegPionsMC.push_back(particle);
                }
            }

            //only paired pions - based on the GetFirstMother() method
            for (int i = 0; i < vPosPionsMC.size(); i++)
            {
                for (int j = 0; j < vNegPionsMC.size(); j++)
                {
                    if (vPosPionsMC[i]->GetFirstMother() == vNegPionsMC[j]->GetFirstMother())
                    {
                        vPosPionsMCPaired.push_back(vPosPionsMC[i]);
                        vNegPionsMCPaired.push_back(vNegPionsMC[j]);
                    }
                }
            }

            if (vPosPionsMCPaired.size() == 0 or vNegPionsMCPaired.size() == 0 ) 
            {
                vPosPionsMC.clear();
                vNegPionsMC.clear();
                vPosPionsMCPaired.clear();
                vNegPionsMCPaired.clear();
                vProtons.clear();
                continue;
            }
            
            for (int i = 0; i < vPosPionsMCPaired.size(); i++)
            {   
                vector <StUPCTrack*> vTpcTracks;
                int num = 0;

                lorentzVectorPosPion.SetPxPyPzE(vPosPionsMCPaired[i]->Px(), vPosPionsMCPaired[i]->Py(), vPosPionsMCPaired[i]->Pz(), vPosPionsMCPaired[i]->Energy());
                lorentzVectorNegPion.SetPxPyPzE(vNegPionsMCPaired[i]->Px(), vNegPionsMCPaired[i]->Py(), vNegPionsMCPaired[i]->Pz(), vNegPionsMCPaired[i]->Energy());

                // if at least one pion is outside of the acceptance then proceed to the next pair...
                if ( lorentzVectorPosPion.Pt() <= 0.1 or abs(lorentzVectorPosPion.Eta()) >= 1 )
                {
                    vPosPionsMC.clear();
                    vNegPionsMC.clear();
                    vPosPionsMCPaired.clear();
                    vNegPionsMCPaired.clear();
                    vProtons.clear();
                    continue;
                }
                
                if ( lorentzVectorNegPion.Pt() <= 0.1 or abs(lorentzVectorNegPion.Eta()) >= 1 )
                {
                    vPosPionsMC.clear();
                    vNegPionsMC.clear();
                    vPosPionsMCPaired.clear();
                    vNegPionsMCPaired.clear();
                    vProtons.clear();
                    continue;
                }

                lorentzVectorNeutKaon = lorentzVectorPosPion + lorentzVectorNegPion;

                TLorentzVector productionVertex;
                vPosPionsMCPaired[i]->ProductionVertex(productionVertex);

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
                            vTpcTracks.push_back(upcEvt->getTrack(j));
                            num+=1;
                        }

                        if ( abs(lorentzVectorNegPion.Eta() - upcEvt->getTrack(j)->getEta()) < 0.1  and  abs(lorentzVectorNegPion.Phi() - upcEvt->getTrack(j)->getPhi()) < 0.1  and upcEvt->getTrack(i)->getFlag(StUPCTrack::kTof)) 
                        {
                            vTpcTracks.push_back(upcEvt->getTrack(j));
                            num+=1;
                        }
                    }         

                    // phi - eta consistency with truth level for exacly two tracks                               
                    if (num != 2) 
                    {       
                        vPosPionsMC.clear();
                        vNegPionsMC.clear();
                        vPosPionsMCPaired.clear();
                        vNegPionsMCPaired.clear();
                        vProtons.clear();
                        vTpcTracks.clear();
                        continue;
                    }
                }

                else 
                {       
                    vPosPionsMC.clear();
                    vNegPionsMC.clear();
                    vPosPionsMCPaired.clear();
                    vNegPionsMCPaired.clear();
                    vProtons.clear();
                    vTpcTracks.clear();
                    continue;
                }
                
                TLorentzVector tpcKaon = {0,0,0,0};

                int c = 0;
                double vertexTpcR = 0;
                double vertexTpcZ = 0;

                for (int i = 0; i < vTpcTracks.size(); i++)
                {
                    vTpcTracks[i]->getLorentzVector(trackVector, massPion);
                    if ( vTpcTracks[i]->getCharge() < 0.0)
                    {
                        tpcKaon+=trackVector;
                        c+=1;
                    }

                    else
                    {
                        tpcKaon+=trackVector;
                        c+=1;
                    }
                    vertexTpcR+=sqrt( pow(vTpcTracks[i]->getVertex()->getPosX() ,2)+ pow(vTpcTracks[i]->getVertex()->getPosY() ,2) )/2.0;
                    vertexTpcZ+=vTpcTracks[i]->getVertex()->getPosZ()/2.0;
                }    
    
                if (c!=2) 
                {
                    vPosPionsMC.clear();
                    vNegPionsMC.clear();
                    vPosPionsMCPaired.clear();
                    vNegPionsMCPaired.clear();
                    vTpcTracks.clear();
                    vProtons.clear();
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

            vPosPionsMC.clear();
            vNegPionsMC.clear();
            vPosPionsMCPaired.clear();
            vNegPionsMCPaired.clear();
            vProtons.clear();   

        }   // End Monte Carlo Efficiency
    }

    TFile *outfile = TFile::Open(argv[2], "recreate"); 
 
    HistCutFlow->Write(); 

    ProfileDistVertexBeamX->Write();
    ProfileDistVertexBeamY->Write(); 
    HistR->Write();
   	HistSig->Write();    

    HistVtxX->Write();
    HistVtxY->Write();
    HistVtxZ->Write();
       
    HistVtxErrX->Write();
    HistVtxErrY->Write();
    HistVtxErrZ->Write();

    HistNumPrimaryVertices->Write(); 
    HistPrimaryVertexAbsPosZ->Write();
    HistNumTofMatchedTracks->Write();
    HistTofMatchedTracksCharge->Write();
    HistTofMatchedTracksAbsEta->Write();
    HistTofMatchedTracksAbsDcaZ->Write();
    HistTofMatchedTracksDcaXY->Write();
    HistTofMatchedTracksNfit->Write();
    HistTofMatchedTracksNdEdx->Write();

    HistNumPrimaryVerticesPostSelection->Write();
    HistPrimaryVertexAbsPosZPostSelection->Write();
    HistNumTofMatchedTracksPostSelection->Write();
    HistTofMatchedTracksChargePostSelection->Write();
    HistTofMatchedTracksAbsEtaPostSelection->Write();
    HistTofMatchedTracksAbsDcaZPostSelection->Write();
    HistTofMatchedTracksDcaXYPostSelection->Write();
    HistTofMatchedTracksNfitPostSelection->Write();
    HistTofMatchedTracksNdEdxPostSelection->Write();

    HistPtProtonPtKaonsX->Write();
    HistPtProtonPtKaonsY->Write();
    HistPtProtonPtKaonsX1d->Write();
    HistPtProtonPtKaonsY1d->Write();
    HistPtProtonPtKaonsXOppY->Write();
    HistPtProtonPtKaonsYOppY->Write();
    HistPtProtonPtKaonsX1dOppY->Write();
    HistPtProtonPtKaonsY1dOppY->Write();
    HistPtProtonPtKaonsXSameY->Write();
    HistPtProtonPtKaonsYSameY->Write();
    HistPtProtonPtKaonsX1dSameY->Write();
    HistPtProtonPtKaonsY1dSameY->Write();

    HistInvMassPiPi->Write();
    HistInvMassPiPiPeak->Write();
    HistInvMassPiPi2->Write();
    HistInvMassPiPiPeak2->Write();

    HistInvMassPiPiPostSelection->Write();
    HistInvMassPiPiPeakPostSelection->Write();
    HistInvMassPiPi2PostSelection->Write();
    HistInvMassPiPiPeak2PostSelection->Write();   

    if (isMC == 1)
    {
        HistKaonPtTruth->Write();
        HistKaonEtaTruth->Write();
        HistKaonVtxRTruth->Write();
        HistKaonVtxZTruth->Write();
        
        HistKaonPtDet->Write();
        HistKaonEtaDet->Write();
        HistKaonVtxRDet->Write();
        HistKaonVtxZDet->Write();
    }

    outfile->Close();

    return 0;
}
