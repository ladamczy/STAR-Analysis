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
    TH1D* HistInvMassPiPi = new TH1D("HistInvMassPiPi", "; m_{#pi0#pi0} [GeV]; # events", 250 ,0, 3);
    TH1D* HistInvMassPiPiPeak = new TH1D("HistInvMassPiPiPeak", "; m_{#pi0#pi0} [GeV]; # events", 100 ,0.4, 0.6);
    TH1D* HistInvMassPiPi2 = new TH1D("HistInvMassPiPi2", "; m_{#pi0#pi0} [GeV]; # events", 250 ,0, 3);
    TH1D* HistInvMassPiPiPeak2 = new TH1D("HistInvMassPiPiPeak2", "; m_{#pi0#pi0} [GeV]; # events", 100 ,0.4, 0.6);
    TH1D* HistInvMassPiPiMassWindow = new TH1D("HistInvMassPiPiMassWindow", "; m_{#pi0#pi0} [GeV]; # events", 100, 0.4, 0.6);
    TProfile* ProfileDistVertexBeamX = new TProfile("ProfileDistVertexBeamX","",639, 20511.5, 21150.5);
    TProfile* ProfileDistVertexBeamY = new TProfile("ProfileDistVertexBeamY","",639, 20511.5, 21150.5);
    TH1D* HistR = new TH1D("HistR", ";R [cm]; # events", 150, 0, 0.3);
    TH1D* HistSig = new TH1D("HistSig"," ; significance; # events", 250, 0.0, 10);

    //cut histograms 
    TH1D* HistTriggers = new TH1D("1_HistTriggers", ";trigger number; count", 17, -0.5, 16.5);
    TH1D* HistRpTrackBranch = new TH1D("2_HistRpTrackBranch", ";RP track branch; count", 4, -0.5, 3.5);
    TH1D* HistRpTrackPointPlanesUsed = new TH1D("3_HistRpTrackPointPlanesUsed", ";Number Of Planes per TrackPoint; count", 5, -0.5, 4.5);
    TH2D* HistRpTrackPxPy= new TH2D("4_HistRpTrackPxPy", ";p_{x} [GeV]; p_{y} [GeV]", 100, -1.5, 1.5, 100, -1.5, 1.5);	
    TH1D* HistNumPrimaryVertices = new TH1D("5_HistNumPrimaryVertices", " ;Number of primary vertices; count", 11, -0.5, 10.5); 
    TH1D* HistPrimaryVertexAbsPosZ = new TH1D("6_HistPrimaryVertexAbsPosZ", ";primary vertex z [cm]; count", 100, 0, 200);
    TH1D* HistNumTofMatchedTracks = new TH1D("7_HistNumTofMatchedTracks", ";number of TOF matched tracks; count", 16, -0.5, 15.5);    
    TH1D* HistTofMatchedTracksCharge = new TH1D("8_HistTofMatchedTracksCharge", ";charge [e]; count", 3, -1.5, 1.5);
    TH1D* HistTofMatchedTracksAbsEta = new TH1D("9_HistTofMatchedTracksAbsEta", ";#eta; count", 50, 0, 2);
    TH1D* HistTofMatchedTracksAbsDcaZ = new TH1D("10_HistTofMatchedTracksAbsDcaZ", ";DCA_{Z} [cm]; count ", 60, 0, 4);
    TH1D* HistTofMatchedTracksDcaXY = new TH1D("11_HistTofMatchedTracksDcaXY", ";DCA_{XY} [cm]; count", 60, 0, 4);
    TH1D* HistTofMatchedTracksNfit = new TH1D("12_HistTofMatchedTracksNfit", ";N_{fit}; count", 61, -0.5, 60.5);
    TH1D* HistTofMatchedTracksNdEdx = new TH1D("13_HistTofMatchedTracksNdEdx", ";N_{dE/dx}; count", 61, -0.5, 60.5);
    
    TH1D* HistTriggersPostSelection = new TH1D("1_HistTriggersPs", ";trigger number; count", 17, -0.5, 16.5);
    TH1D* HistRpTrackBranchPostSelection = new TH1D("2_HistRpTrackBranchPs", ";RP track branch; count", 4, -0.5, 3.5);
    TH1D* HistRpTrackPointPlanesUsedPostSelection = new TH1D("3_HistRpTrackPointPlanesUsedPs", ";Number Of Planes per TrackPoint; count", 5, -0.5, 4.5);
    TH2D* HistRpTrackPxPyPostSelection = new TH2D("4_HistRpTrackPxPyPs", ";p_{x} [GeV]; p_{y} [GeV]", 100, -1.5, 1.5, 100, -1.5, 1.5);	
    TH1D* HistNumPrimaryVerticesPostSelection = new TH1D("5_HistNumPrimaryVerticesPs", " ;Number of primary vertices; count", 11, -0.5, 10.5); 
    TH1D* HistPrimaryVertexAbsPosZPostSelection = new TH1D("6_HistPrimaryVertexAbsPosZPs", ";primary vertex z [cm]; count", 100, 0, 200);
    TH1D* HistNumTofMatchedTracksPostSelection = new TH1D("7_HistNumTofMatchedTracksPs", ";number of TOF matched tracks; count", 16, -0.5, 15.5);    
    TH1D* HistTofMatchedTracksChargePostSelection = new TH1D("8_HistTofMatchedTracksChargePs", ";charge [e]; count", 3, -1.5, 1.5);
    TH1D* HistTofMatchedTracksAbsEtaPostSelection = new TH1D("9_HistTofMatchedTracksAbsEtaPs", ";#eta; count", 50, 0, 2);
    TH1D* HistTofMatchedTracksAbsDcaZPostSelection = new TH1D("10_HistTofMatchedTracksAbsDcaZPs", ";DCA_{Z} [cm]; count ", 60, 0, 4);
    TH1D* HistTofMatchedTracksDcaXYPostSelection = new TH1D("11_HistTofMatchedTracksDcaXYPs", ";DCA_{XY} [cm]; count", 60, 0, 4);
    TH1D* HistTofMatchedTracksNfitPostSelection = new TH1D("12_HistTofMatchedTracksNfitPs", ";N_{fit}; count", 61, -0.5, 60.5);
    TH1D* HistTofMatchedTracksNdEdxPostSelection = new TH1D("13_HistTofMatchedTracksNdEdxPs", ";N_{dE/dx}; count", 61, -0.5, 60.5);

    // kaon histograms
    TH1D* HistKaonE = new TH1D("HistKaonE", "; E_{K0} [GeV]; # events", 100 ,0, 3);
    TH1D* HistKaonPt = new TH1D("HistKaonPt", "; pT_{K0} [GeV]; # events", 100 ,0, 3);
    TH1D* HistKaonEta = new TH1D("HistKaonEta", "; #eta_{K0}; # events", 101, -2.5, 2.5);
    TH1D* HistKaonTheta = new TH1D("HistKaonTheta", "; #theta_{K0} [rad]; # events", 101, 0, 2*M_PI);
    TH1D* HistKaonPhi = new TH1D("HistKaonPhi", "; #phi; # events;", 101, -M_PI, M_PI);
    TH2D* HistKaonThetaEta = new TH2D("HistKaonThetaEta", "; #theta_{K0} [rad]; #eta_{K0}",101, 0, 2*M_PI, 101, -2.5, 2.5);

    TH1D* HistKaonE2 = new TH1D("HistKaonE2", "; E_{K0} [GeV]; # events", 100 ,0, 3);
    TH1D* HistKaonPt2 = new TH1D("HistKaonPt2", "; pT_{K0} [GeV]; # events", 100 ,0, 3);
    TH1D* HistKaonEta2 = new TH1D("HistKaonEta2", "; #eta_{K0}; # events", 101, -2.5, 2.5);
    TH1D* HistKaonTheta2 = new TH1D("HistKaonTheta2", "; #theta_{K0} [rad]; # events", 101, 0, 2*M_PI);
    TH1D* HistKaonPhi2 = new TH1D("HistKaonPhi2", "; #phi; # events;", 101, -M_PI, M_PI);
    TH2D* HistKaonThetaEta2 = new TH2D("HistKaonThetaEta2", "; #theta_{K0} [rad]; #eta_{K0}",101, 0, 2*M_PI, 101, -2.5, 2.5);

    // additional
    TH1D* HistNumTracksIf3TofTracks = new TH1D("HistNumTracksIf3TofTracks", ";number tracks; count", 16, -0.5, 15.5);
    TH1D* HistPionSelection = new TH1D("HistPionSelection", ";;count", 6, -0.5, 5.5);
    TH2D* HistPtProtonPtKaonsX = new TH2D("HistPtProtonPtKaonsX", ";[p_{K^{0}p_{K^{0}}]_{x}; [p_{p'}^{W}p_{p'}^{E}]_{x} ",  100, -2.0, 2.0, 100, -2.0, 2.0);
    TH2D* HistPtProtonPtKaonsY = new TH2D("HistPtProtonPtKaonsY", ";[p_{K^{0}p_{K^{0}}]_{y}; [p_{p'}^{W}p_{p'}^{E}]_{y} ",  100, -2.0, 2.0, 100, -2.0, 2.0);

    static StUPCEvent *upcEvt = 0x0;
    static StRPEvent *rpEvt = 0x0;
    chain->SetBranchAddress("mUPCEvent", &upcEvt);
    chain->SetBranchAddress("mRPEvent", &rpEvt);

    vector<TLorentzVector> posPion;
    vector<TLorentzVector> negPion;
    TLorentzVector trackVector;
	TVector3 pVector;
	TVector3 proton1, proton2;

    double beamPositionX, beamPositionY;
    double primVertexPosX, primVertexPosY, primVertexPosErrX, primVertexPosErrY;
    double distVertexBeamX, distVertexBeamY, distR, significance;
    int iCutFlow;

    const vector <int> triggerID = {570209, 570219, 570229, 570701, 570702, 570703, 570704, 570705, 570709, 570711, 570712, 570719, 590701, 590703, 590705, 590708, 590709};
    const vector <int> triggerCEP = { 570701, 570705, 570711, 590701, 590705, 590708};

    for (Long64_t i = 0; i < chain->GetEntries(); ++i) 
    {

        iCutFlow = 0;

        chain->GetEntry(i);
        HistCutFlow->Fill(iCutFlow);  // cutflow: all data

        if (i%100000 == 0)  { cout << i << "/" <<  chain->GetEntries()  << endl;}

        // trigger histogram - event can have several triggers
        for (unsigned long int i = 0; i < triggerID.size(); i++)
        {
            if (upcEvt->isTrigger(triggerID[i]))
            {
                HistTriggers->Fill(i);
            }
        }
        
        // check CEP triggers
        bool isCepTrigger = 0;
        for (unsigned long int i = 0; i < triggerCEP.size(); i++)
        {
            if (upcEvt->isTrigger(triggerCEP[i]))
            {
                isCepTrigger = 1;
            }
        }
        if (isCepTrigger == 0) {continue;} 
        iCutFlow+=1;
        HistCutFlow->Fill(iCutFlow); // CEP trigger

        for (unsigned long int i = 0; i < triggerID.size(); i++)
        {
            if (upcEvt->isTrigger(triggerID[i]))
            {
                HistTriggersPostSelection->Fill(i);
            }
        }        

        // two RP tracks on the opposite sides of the detector
        bool hasTwoRpTracksOppositeSide = 0;
		int numberOfTracksEast = 0;
		int numberOfTracksWest = 0;
		int indexWest = 0;
		int indexEast = 0;
		for(UInt_t i = 0; i < rpEvt->getNumberOfTracks(); i++)
		{
			int iBranch = rpEvt->getTrack(i)->branch();
            HistRpTrackBranch->Fill(iBranch);
			if (iBranch < 2) {numberOfTracksEast+=1; indexEast = i;}
		    else {numberOfTracksWest+=1; indexWest  = i;}
	    }
		if (numberOfTracksWest == 1 and numberOfTracksEast == 1)
		{
			hasTwoRpTracksOppositeSide = 1;
		}
		if (hasTwoRpTracksOppositeSide == 0) {continue;}	
        iCutFlow+=1;
  		HistCutFlow->Fill(iCutFlow); // cutflow: two RP tracks on the opposite sides of the detector

		for(UInt_t i = 0; i < rpEvt->getNumberOfTracks(); i++)
        {
			int iBranch = rpEvt->getTrack(i)->branch();
            HistRpTrackBranchPostSelection->Fill(iBranch);
        }

        // RP trackpoints planes used
    	bool useAtLeastThreePlanes = 1;
		for(UInt_t i = 0; i < rpEvt->getNumberOfTrackPoints(); i++)
		{
            HistRpTrackPointPlanesUsed->Fill(rpEvt->getTrackPoint(i)->planesUsed());
			if(rpEvt->getTrackPoint(i)->planesUsed() < 3)
			{
				useAtLeastThreePlanes = 0;		
				break;
			}
		}		
		if (useAtLeastThreePlanes == 0) {continue;}	
        iCutFlow+=1;
  		HistCutFlow->Fill(iCutFlow); //  cutflow 3 out of 4 silicon planes were used for the reconstruction

		for(UInt_t i = 0; i < rpEvt->getNumberOfTrackPoints(); i++)
		{
            HistRpTrackPointPlanesUsedPostSelection->Fill(rpEvt->getTrackPoint(i)->planesUsed());
        }


        // fiducial region
    	bool withinFiducialRegion = 1;
	    for(unsigned int k = 0; k < rpEvt->getNumberOfTracks(); ++k)
	    {		
			StUPCRpsTrack *trk = rpEvt->getTrack(k);
		    trk->setEvent(rpEvt);
		    pVector = trk->pVec();
            if (k == 0)
            {
		        proton1 = trk->pVec();
            }
            else if (k == 1)
            {
		        proton2 = trk->pVec();
            }

			HistRpTrackPxPy->Fill(pVector.X(), pVector.Y());

			if ((abs(pVector.Y()) >= 0.8 and abs(pVector.Y()) <= 0.4) or (pVector.X() <= -0.27) or (pow(pVector.X()+0.6,2) + pow(pVector.Y(), 2)  >=  pow(1.25,2)))
			{
				withinFiducialRegion = 0;
			}  	
		}
		if (withinFiducialRegion == 0) {continue;}    
        iCutFlow+=1;
		HistCutFlow->Fill(iCutFlow); // cutflow: fiducial region

	    for(unsigned int k = 0; k < rpEvt->getNumberOfTracks(); ++k)
	    {		
			StUPCRpsTrack *trk = rpEvt->getTrack(k);
		    trk->setEvent(rpEvt);
		    pVector = trk->pVec();
            if (k == 0)
            {
		        proton1 = trk->pVec();
            }
            else if (k == 1)
            {
		        proton2 = trk->pVec();
            }

			HistRpTrackPxPyPostSelection->Fill(pVector.X(), pVector.Y());
        }



        // number of primary vertices - fos single Kaon we demand exacly one primary vertex
        Int_t numberOfPrimaryVertices = upcEvt->getNumberOfVertices();       	
        HistNumPrimaryVertices->Fill(numberOfPrimaryVertices);		
        if(numberOfPrimaryVertices!=1) {continue;}		
        iCutFlow+=1;
		HistCutFlow->Fill(iCutFlow); // cutflow: one vertex
        HistNumPrimaryVerticesPostSelection->Fill(numberOfPrimaryVertices);		
        // fill number and beam position
        int nFillNumber = upcEvt->getFillNumber();
		beamPositionX = FindPosition(nFillNumber, fillNumberWithPosition[0], fillNumberWithPosition[1], fillNumberWithPosition[2], fillNumberWithPosition[3], fillNumberWithPosition[4])[0];
		beamPositionY = FindPosition(nFillNumber, fillNumberWithPosition[0], fillNumberWithPosition[1], fillNumberWithPosition[2], fillNumberWithPosition[3], fillNumberWithPosition[4])[1];			  

        // profile - distance between vertex and beam in x-y plane, parameter R and significance
        primVertexPosX = upcEvt->getVertex(0)->getPosX();
        primVertexPosY = upcEvt->getVertex(0)->getPosY();
        primVertexPosErrX = upcEvt->getVertex(0)->getErrX();
        primVertexPosErrY = upcEvt->getVertex(0)->getErrY();
        distVertexBeamX = primVertexPosX - beamPositionX;
        distVertexBeamY = primVertexPosY - beamPositionY;
        distR = sqrt(pow(distVertexBeamX,2) + pow(distVertexBeamY,2));
   		significance = sqrt(pow((distVertexBeamX/primVertexPosErrX),2) + pow((distVertexBeamY/primVertexPosErrY),2));            		
        ProfileDistVertexBeamX->Fill(nFillNumber, distVertexBeamX); 
        ProfileDistVertexBeamY->Fill(nFillNumber, distVertexBeamY); 
        HistR->Fill(distR);
   		HistSig->Fill(significance);    

		//primary vertex is placed within |zvtx| < 80 cm
        HistPrimaryVertexAbsPosZ->Fill(abs(upcEvt->getVertex(0)->getPosZ()));
		if(abs(upcEvt->getVertex(0)->getPosZ())>=80)
		{
			continue;
        }   
        iCutFlow+=1;
		HistCutFlow->Fill(iCutFlow); // cutflow: primary vertex |zvtx| < 80 cm
        HistPrimaryVertexAbsPosZPostSelection->Fill(abs(upcEvt->getVertex(0)->getPosZ()));      

		// vector for TOF matched tracks
		vector <StUPCTrack*> tracksWithTofHit;
        for (Int_t i = 0; i<upcEvt->getNumberOfTracks(); i++)
		{
			if(upcEvt->getTrack(i)->getFlag(StUPCTrack::kTof))
			{
                tracksWithTofHit.push_back(upcEvt->getTrack(i));
			}
		}       		

        //additional
		if (tracksWithTofHit.size() == 3)
		{
			HistNumTracksIf3TofTracks->Fill(upcEvt->getNumberOfTracks()-3);
		}

		// two TOF-matched primary tracks
        HistNumTofMatchedTracks->Fill(tracksWithTofHit.size());
		bool isValidNumberOfTofMatchedTracks = 0;
		if (tracksWithTofHit.size() == 4)
		{
			isValidNumberOfTofMatchedTracks = 1;
		}
		if (isValidNumberOfTofMatchedTracks == 0) {continue;} 
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
        if (sumCharge != 0)    {continue;}
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

        double kaonMassWindowLow = 0.485;
        double kaonMassWindowHigh = 0.509;

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
            HistPionSelection->Fill(0); // 4
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
        }

        else if (arePionPairsInMassWindowFirstCombination and isPionPairInMassWindowSecondCombinationA)
        {
            HistPionSelection->Fill(1); // 2+1
            lVecCombinations = {lVecKaonComb1a, lVecKaonComb1b, lVecKaonComb2a}; 
            pTCombinations = {lVecKaonComb1a.Pt(), lVecKaonComb1b.Pt(), lVecKaonComb2a.Pt()}; 
            maxPtIndex = distance(pTCombinations.begin(),  max_element(pTCombinations.begin(), pTCombinations.end()));

            if (maxPtIndex <= 1)
            {
                conditioningKaon = lVecCombinations[maxPtIndex]; 
                complementaryKaon = lVecCombinations[abs(maxPtIndex-1)];
            }

            else
            {
                conditioningKaon = lVecKaonComb2a;
                complementaryKaon = lVecKaonComb2b;
            }
        }

        else if (arePionPairsInMassWindowFirstCombination and isPionPairInMassWindowSecondCombinationB)
        {
            HistPionSelection->Fill(1); // 2+1
            lVecCombinations = {lVecKaonComb1a, lVecKaonComb1b, lVecKaonComb2b}; 
            pTCombinations = {lVecKaonComb1a.Pt(), lVecKaonComb1b.Pt(), lVecKaonComb2b.Pt()}; 
            maxPtIndex = distance(pTCombinations.begin(), max_element(pTCombinations.begin(), pTCombinations.end()));

            if (maxPtIndex <= 1)
            {
                conditioningKaon = lVecCombinations[maxPtIndex]; 
                complementaryKaon = lVecCombinations[abs(maxPtIndex-1)];
            }

            else
            {
                conditioningKaon = lVecKaonComb2b;
                complementaryKaon = lVecKaonComb2a;
            }
        }

        else if (arePionPairsInMassWindowSecondCombination and isPionPairInMassWindowFirstCombinationA)
        {
            HistPionSelection->Fill(1); // 2+1
            lVecCombinations = {lVecKaonComb2a, lVecKaonComb2b, lVecKaonComb1a}; 
            pTCombinations = {lVecKaonComb2a.Pt(), lVecKaonComb2b.Pt(), lVecKaonComb1a.Pt()}; 
            maxPtIndex = distance(pTCombinations.begin(), max_element(pTCombinations.begin(), pTCombinations.end()));

            if (maxPtIndex <= 1)
            {
                conditioningKaon = lVecCombinations[maxPtIndex]; 
                complementaryKaon = lVecCombinations[abs(maxPtIndex-1)];
            }

            else
            {
                conditioningKaon = lVecKaonComb1a;
                complementaryKaon = lVecKaonComb1b;
            }
        }

        else if (arePionPairsInMassWindowSecondCombination and isPionPairInMassWindowFirstCombinationB)
        {
            HistPionSelection->Fill(1); // 2+1
            lVecCombinations = {lVecKaonComb2a, lVecKaonComb2b, lVecKaonComb1b}; 
            pTCombinations = {lVecKaonComb2a.Pt(), lVecKaonComb2b.Pt(), lVecKaonComb1b.Pt()}; 
            maxPtIndex = distance(pTCombinations.begin(), max_element(pTCombinations.begin(), pTCombinations.end()));

            if (maxPtIndex <= 1)
            {
                conditioningKaon = lVecCombinations[maxPtIndex]; 
                complementaryKaon = lVecCombinations[abs(maxPtIndex-1)];
            }

            else
            {
                conditioningKaon = lVecKaonComb1b;
                complementaryKaon = lVecKaonComb1a;
            }
        }


        else if (arePionPairsInMassWindowFirstCombination)
        {
            HistPionSelection->Fill(2); // 2
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
        }


        else if (arePionPairsInMassWindowSecondCombination)
        {
            HistPionSelection->Fill(2); // 2
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
        }       

        else if (isPionPairInMassWindowFirstCombinationA and isPionPairInMassWindowSecondCombinationA)
        {
            HistPionSelection->Fill(3); // 1+1
            if (lVecKaonComb1a.Pt() > lVecKaonComb2a.Pt())
            {
                conditioningKaon = lVecKaonComb1a;
                complementaryKaon = lVecKaonComb1b;
            }

            else
            {
                conditioningKaon = lVecKaonComb2a;
                complementaryKaon = lVecKaonComb2b;              
            }
        }

        else if (isPionPairInMassWindowFirstCombinationA and isPionPairInMassWindowSecondCombinationB)
        {
            HistPionSelection->Fill(3); // 1+1
            if (lVecKaonComb1a.Pt() > lVecKaonComb2b.Pt())
            {
                conditioningKaon = lVecKaonComb1a;
                complementaryKaon = lVecKaonComb1b;
            }

            else
            {
                conditioningKaon = lVecKaonComb2b;
                complementaryKaon = lVecKaonComb2a;              
            }  
        }

        else if (isPionPairInMassWindowFirstCombinationB and isPionPairInMassWindowSecondCombinationA)
        {
            HistPionSelection->Fill(3); // 1+1
            if (lVecKaonComb1b.Pt() > lVecKaonComb2a.Pt())
            {
                conditioningKaon = lVecKaonComb1b;
                complementaryKaon = lVecKaonComb1a;
            }

            else
            {
                conditioningKaon = lVecKaonComb2a;
                complementaryKaon = lVecKaonComb2b;              
            }
        }

        else if (isPionPairInMassWindowFirstCombinationB and isPionPairInMassWindowSecondCombinationB)
        {
            HistPionSelection->Fill(3); // 1+1
            if (lVecKaonComb1b.Pt() > lVecKaonComb2b.Pt())
            {
                conditioningKaon = lVecKaonComb1b;
                complementaryKaon = lVecKaonComb1a;
            }

            else
            {
                conditioningKaon = lVecKaonComb2b;
                complementaryKaon = lVecKaonComb2a;              
            }   
        }

        else if (isPionPairInMassWindowFirstCombinationA)
        {
            HistPionSelection->Fill(4); // 1
            conditioningKaon = lVecKaonComb1a;
            complementaryKaon = lVecKaonComb1b;
        }

        else if (isPionPairInMassWindowFirstCombinationB)
        {
            HistPionSelection->Fill(4); // 1
            conditioningKaon = lVecKaonComb1b;
            complementaryKaon = lVecKaonComb1a;
        }

        else if (isPionPairInMassWindowSecondCombinationA)
        {
            HistPionSelection->Fill(4); // 1
            conditioningKaon = lVecKaonComb2a;
            complementaryKaon = lVecKaonComb2b;
        }

        else if (isPionPairInMassWindowSecondCombinationB)
        {
            HistPionSelection->Fill(4); // 1
            conditioningKaon = lVecKaonComb2b;
            complementaryKaon = lVecKaonComb2a;
        }

        else 
        {
            HistPionSelection->Fill(5); //0
            
            posPion.clear();
            negPion.clear();  
            pTCombinations.clear();
            continue; 
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
        }

        HistInvMassPiPi->Fill(conditioningKaon.M());
        HistInvMassPiPiPeak->Fill(conditioningKaon.M());
        HistKaonE->Fill(conditioningKaon.E());
        HistKaonPt->Fill(conditioningKaon.Pt());
        HistKaonEta->Fill(conditioningKaon.Eta());
        HistKaonTheta->Fill(conditioningKaon.Theta());
        HistKaonPhi->Fill(conditioningKaon.Phi());
        HistKaonThetaEta->Fill(conditioningKaon.Theta(), conditioningKaon.Eta());

        HistInvMassPiPi2->Fill(complementaryKaon.M());
        HistInvMassPiPiPeak2->Fill(complementaryKaon.M());
        HistKaonE2->Fill(complementaryKaon.E());
        HistKaonPt2->Fill(complementaryKaon.Pt());
        HistKaonEta2->Fill(complementaryKaon.Eta());
        HistKaonTheta2->Fill(complementaryKaon.Theta());
        HistKaonPhi2->Fill(complementaryKaon.Phi());
        HistKaonThetaEta2->Fill(complementaryKaon.Theta(), complementaryKaon.Eta());

        HistPtProtonPtKaonsX->Fill(conditioningKaon.Px()+complementaryKaon.Px(), proton1.X() +  proton2.X());
        HistPtProtonPtKaonsY->Fill(conditioningKaon.Py()+complementaryKaon.Py(), proton1.Y() +  proton2.Y());

        posPion.clear();
        negPion.clear();   
        pTCombinations.clear();              
    }

    TFile *outfile = TFile::Open(argv[2], "recreate"); 
 
    HistCutFlow->Write(); 

    ProfileDistVertexBeamX->Write();
    ProfileDistVertexBeamY->Write();
    HistR->Write();
    HistSig->Write();
    HistTriggers->Write();
    HistRpTrackBranch->Write();
    HistRpTrackPointPlanesUsed->Write();
    HistRpTrackPxPy->Write();
    HistNumPrimaryVertices->Write();
    HistPrimaryVertexAbsPosZ->Write();
    HistNumTofMatchedTracks->Write();
    HistTofMatchedTracksCharge->Write();
    HistTofMatchedTracksAbsEta->Write();
    HistTofMatchedTracksAbsDcaZ->Write();
    HistTofMatchedTracksDcaXY->Write();
    HistTofMatchedTracksNfit->Write();
    HistTofMatchedTracksNdEdx->Write(); 

    HistTriggersPostSelection->Write();
    HistRpTrackBranchPostSelection->Write();
    HistRpTrackPointPlanesUsedPostSelection->Write();
    HistRpTrackPxPyPostSelection->Write();
    HistNumPrimaryVerticesPostSelection->Write();
    HistPrimaryVertexAbsPosZPostSelection->Write();
    HistNumTofMatchedTracksPostSelection->Write();
    HistTofMatchedTracksChargePostSelection->Write();
    HistTofMatchedTracksAbsEtaPostSelection->Write();
    HistTofMatchedTracksAbsDcaZPostSelection->Write();
    HistTofMatchedTracksDcaXYPostSelection->Write();
    HistTofMatchedTracksNfitPostSelection->Write();
    HistTofMatchedTracksNdEdxPostSelection->Write(); 

    HistInvMassPiPi->Write();
    HistInvMassPiPiPeak->Write();
    HistInvMassPiPi2->Write();
    HistInvMassPiPiPeak2->Write();

    HistKaonE->Write();
    HistKaonPt->Write();
    HistKaonEta->Write();
    HistKaonTheta->Write();
    HistKaonPhi->Write();
    HistKaonThetaEta->Write();

    HistKaonE2->Write();
    HistKaonPt2->Write();
    HistKaonEta2->Write();
    HistKaonTheta2->Write();
    HistKaonPhi2->Write();
    HistKaonThetaEta2->Write();

    HistNumTracksIf3TofTracks->Write();
    HistPionSelection->Write();
    HistPtProtonPtKaonsX->Write();
    HistPtProtonPtKaonsY->Write();
    outfile->Close();

    return 0;
}
