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
    while (std::getline(inputFilePathList, inputFileName))
    {
        chain->Add(inputFileName.c_str());
    }

    inputFilePathList.close();

    // to be changed... there is no warning in case of the invalid input file
	vector <vector<double>> fillNumberWithPosition = ReadFillPositionData("InputData/Run7PolarizationWithPosition.csv");

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

    //cut histograms 
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

    TH2D* HistPtProtonPtKaonsX = new TH2D("HistPtProtonPtKaonsX", ";[p_{K^{0}p_{K^{0}}]_{x}; [p_{p'}^{W}p_{p'}^{E}]_{x} ",  100, -2.0, 2.0, 100, -2.0, 2.0);
    TH2D* HistPtProtonPtKaonsY = new TH2D("HistPtProtonPtKaonsY", ";[p_{K^{0}p_{K^{0}}]_{y}; [p_{p'}^{W}p_{p'}^{E}]_{y} ",  100, -2.0, 2.0, 100, -2.0, 2.0);
    TH1D* HistPtProtonPtKaonsX1d = new TH1D("HistPtProtonPtKaonsX1d", ";[p_{K^{0}p_{K^{0}} + [p_{p'}^{W}p_{p'}^{E}]_{x}; events ",  200, -5.0, 5.0);
    TH1D* HistPtProtonPtKaonsY1d = new TH1D("HistPtProtonPtKaonsY1d", ";[p_{K^{0}p_{K^{0}} +  [p_{p'}^{W}p_{p'}^{E}]_{y}; events ",  200, -5.0, 5.0);
    TH2D* HistPtProtonPtKaonsXOppY = new TH2D("HistPtProtonPtKaonsXOppY", ";[p_{K^{0}p_{K^{0}}]_{x}; [p_{p'}^{W}p_{p'}^{E}]_{x} ",  100, -2.0, 2.0, 100, -2.0, 2.0);
    TH2D* HistPtProtonPtKaonsYOppY = new TH2D("HistPtProtonPtKaonsYOppY", ";[p_{K^{0}p_{K^{0}}]_{y}; [p_{p'}^{W}p_{p'}^{E}]_{y} ",  100, -2.0, 2.0, 100, -2.0, 2.0);
    TH1D* HistPtProtonPtKaonsX1dOppY = new TH1D("HistPtProtonPtKaonsX1dOppY", ";[p_{K^{0}p_{K^{0}} + [p_{p'}^{W}p_{p'}^{E}]_{x}; events ",  200, -5.0, 5.0);
    TH1D* HistPtProtonPtKaonsY1dOppY = new TH1D("HistPtProtonPtKaonsY1dOppY", ";[p_{K^{0}p_{K^{0}} +  [p_{p'}^{W}p_{p'}^{E}]_{y}; events ",  200, -5.0, 5.0);
    TH2D* HistPtProtonPtKaonsXSameY = new TH2D("HistPtProtonPtKaonsXSameY", ";[p_{K^{0}p_{K^{0}}]_{x}; [p_{p'}^{W}p_{p'}^{E}]_{x} ",  100, -2.0, 2.0, 100, -2.0, 2.0);
    TH2D* HistPtProtonPtKaonsYSameY = new TH2D("HistPtProtonPtKaonsYSameY", ";[p_{K^{0}p_{K^{0}}]_{y}; [p_{p'}^{W}p_{p'}^{E}]_{y} ",  100, -2.0, 2.0, 100, -2.0, 2.0);
    TH1D* HistPtProtonPtKaonsX1dSameY = new TH1D("HistPtProtonPtKaonsX1dSameY", ";[p_{K^{0}p_{K^{0}} + [p_{p'}^{W}p_{p'}^{E}]_{x}; events ",  200, -5.0, 5.0);
    TH1D* HistPtProtonPtKaonsY1dSameY = new TH1D("HistPtProtonPtKaonsY1dSameY", ";[p_{K^{0}p_{K^{0}} +  [p_{p'}^{W}p_{p'}^{E}]_{y}; events ",  200, -5.0, 5.0);

    TH1D* HistInvMassPiPi = new TH1D("HistInvMassPiPi", "; m_{#pi0#pi0} [GeV]; # events", 250 ,0, 3);
    TH1D* HistInvMassPiPiPeak = new TH1D("HistInvMassPiPiPeak", "; m_{#pi0#pi0} [GeV]; # events", 100 ,0.4, 0.6);
    TH1D* HistInvMassPiPi2 = new TH1D("HistInvMassPiPi2", "; m_{#pi0#pi0} [GeV]; # events", 250 ,0, 3);
    TH1D* HistInvMassPiPiPeak2 = new TH1D("HistInvMassPiPiPeak2", "; m_{#pi0#pi0} [GeV]; # events", 100 ,0.4, 0.6);

    TH1D* HistInvMassPiPiPostSelection = new TH1D("HistInvMassPiPiPs", "; m_{#pi0#pi0} [GeV]; # events", 250 ,0, 3);
    TH1D* HistInvMassPiPiPeakPostSelection = new TH1D("HistInvMassPiPiPeakPs", "; m_{#pi0#pi0} [GeV]; # events", 100 ,0.4, 0.6);
    TH1D* HistInvMassPiPi2PostSelection = new TH1D("HistInvMassPiPi2Ps", "; m_{#pi0#pi0} [GeV]; # events", 250 ,0, 3);
    TH1D* HistInvMassPiPiPeak2PostSelection = new TH1D("HistInvMassPiPiPeak2Ps", "; m_{#pi0#pi0} [GeV]; # events", 100 ,0.4, 0.6);



    static StUPCEvent *upcEvt = 0x0;
    static StRPEvent *rpEvt = 0x0;
    chain->SetBranchAddress("mUPCEvent", &upcEvt);
   // chain->SetBranchAddress("mRPEvent", &rpEvt);

    vector <TLorentzVector> posPion;
    vector <TLorentzVector> negPion;
    TLorentzVector trackVector;
	TVector3 pVector;
	TVector3 proton1, proton2;

    double beamPositionX, beamPositionY;
    double primVertexPosX, primVertexPosY, primVertexPosZ, primVertexPosErrX, primVertexPosErrY, primVertexPosErrZ;
    double distVertexBeamX, distVertexBeamY, distR, significance;
    int iCutFlow;

    double protonSignsY;
    TParticle* particle;
    vector <TParticle*> Protons;
    const vector <int> triggerID = {570209, 570219, 570229, 570701, 570702, 570703, 570704, 570705, 570709, 570711, 570712, 570719, 590701, 590703, 590705, 590708, 590709};
    const vector <int> triggerCEP = { 570701, 570705, 570711, 590701, 590705, 590708};

    for (Long64_t i = 0; i < chain->GetEntries(); ++i) 
    {

        iCutFlow = 0;

        chain->GetEntry(i);
        HistCutFlow->Fill(iCutFlow);  // cutflow: all data

        if (i%100000 == 0)  { cout << i << "/" <<  chain->GetEntries()  << endl;}

        for (int i = 0; i < upcEvt->getNumberOfMCParticles(); i++)
        {
            particle = upcEvt->getMCParticle(i);
            if (particle->GetPDG()->PdgCode() == 2212 and particle->GetFirstMother() == 1)
            {
                Protons.push_back(particle);
            }
        }

        proton1.SetX(Protons[0]->Px());
        proton1.SetY(Protons[0]->Py());
        proton1.SetZ(Protons[0]->Pz());


        proton2.SetX(Protons[1]->Px());
        proton2.SetY(Protons[1]->Py());
        proton2.SetZ(Protons[1]->Pz());
        cout << Protons[0]->Py() << ", " << Protons[1]->Py() << endl;
        protonSignsY = 1.0;
        protonSignsY = 1.0*Protons[0]->Py()*Protons[1]->Py();

        Protons.clear();





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

		//primary vertex is placed within |zvtx| < 80 cm
        HistPrimaryVertexAbsPosZ->Fill(abs(upcEvt->getVertex(0)->getPosZ()));
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
		HistCutFlow->Fill(iCutFlow);  // cutflow: two tracks quality N_{dE/dx} >= 15 

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


        // fill number and beam position
        int nFillNumber = upcEvt->getFillNumber();
		beamPositionX = FindPosition(nFillNumber, fillNumberWithPosition[0], fillNumberWithPosition[1], fillNumberWithPosition[2], fillNumberWithPosition[3], fillNumberWithPosition[4])[0];
		beamPositionY = FindPosition(nFillNumber, fillNumberWithPosition[0], fillNumberWithPosition[1], fillNumberWithPosition[2], fillNumberWithPosition[3], fillNumberWithPosition[4])[1];			  

        // profile - distance between vertex and beam in x-y plane, parameter R and significance
        primVertexPosX = upcEvt->getVertex(0)->getPosX();
        primVertexPosY = upcEvt->getVertex(0)->getPosY();
        primVertexPosZ = upcEvt->getVertex(0)->getPosZ();

        primVertexPosErrX = upcEvt->getVertex(0)->getErrX();
        primVertexPosErrY = upcEvt->getVertex(0)->getErrY();
        primVertexPosErrZ = upcEvt->getVertex(0)->getErrZ();

        distVertexBeamX = primVertexPosX - 0.0;
        distVertexBeamY = primVertexPosY - 0.0;
        
        distR = sqrt(pow(distVertexBeamX,2) + pow(distVertexBeamY,2));
   		significance = sqrt(pow((distVertexBeamX/primVertexPosErrX),2) + pow((distVertexBeamY/primVertexPosErrY),2));       
                     		
        ProfileDistVertexBeamX->Fill(nFillNumber, distVertexBeamX); 
        ProfileDistVertexBeamY->Fill(nFillNumber, distVertexBeamY); 
        HistR->Fill(distR);
   		HistSig->Fill(significance);    

        HistVtxX->Fill(primVertexPosX);
        HistVtxY->Fill(primVertexPosY);
        HistVtxZ->Fill(primVertexPosZ);
        
        HistVtxErrX->Fill(primVertexPosErrX);
        HistVtxErrY->Fill(primVertexPosErrY);
        HistVtxErrZ->Fill(primVertexPosErrZ);


        posPion.clear();
        negPion.clear();  
        pTCombinations.clear();
        lVecCombinations.clear();
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


    outfile->Close();

    return 0;
}
