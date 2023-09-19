#include "Includes.h"

using namespace std;

int main(int argc, char** argv)  
{
    double massPion = 0.13957; 
    const vector <int> triggerID = {570209, 570219, 570229, 570701, 570702, 570703, 570704, 570705, 570709, 570711, 570712, 570719, 590701, 590703, 590705, 590708, 590709};
    const vector <int> triggerCEP = { 570701, 570705, 570711, 590701, 590705, 590708};

    ifstream file(argv[1]);
    if (!file) 
    {
        cerr << "Failed to open lis.list file." << std::endl;
        return 1;
    }

    TChain chain("mUPCTree"); 

    string line;
    vector <string> rootFiles;
    while (std::getline(file, line))
    {
        chain.Add(line.c_str());
    }
    file.close();

	vector <vector<double>> fillNumberWithPosition = ReadFillPositionData("Run7PolarizationWithPosition.csv");

    TH1D* HistInvMassPiPi = new TH1D("HistInvMassPiPi", "", 100, 0.42, 0.56);
    TProfile* ProfileDistVertexBeamX = new TProfile("ProfileDistVertexBeamX","",639, 20511.5, 21150.5);
    TProfile* ProfileDistVertexBeamY = new TProfile("ProfileDistVertexBeamY","",639, 20511.5, 21150.5);
    TH1D* HistR = new TH1D("HistR", "", 150, 0, 0.3);
    TH1D* HistSig = new TH1D("HistSig","", 100, 0.0, 2);
    TH1D* HistCutFlow = new TH1D("HistCutFlow", "", 16, -0.5, 15.5);

    //cut histograms 
    TH1D* HistTriggers = new TH1D("1_HistTriggers", ";trigger number; count", 17, -0.5, 16.5);
    TH1D* HistNumPrimaryVertices = new TH1D("2_HistNumPrimaryVertices", " ;Number of primary vertices; count", 11, -0.5, 10.5); 
    TH1D* HistPrimaryVertexAbsPosZ = new TH1D("3_HistPrimaryVertexAbsPosZ", ";primary vertex z [cm];count", 100, 0, 200);
    TH1D* HistNumTofMatchedTracks = new TH1D("4_HistNumTofMatchedTracks", ";number of TOF matched tracks;count", 16, -0.5, 15.5);
    TH1D* HistTwoTofMatchedTracksCharge = new TH1D("5_HistTwoTofMatchedTracksCharge", ";charge [e];count", 3, -1.5, 1.5);
    TH1D* HistTwoTofMatchedTracksAbsEta = new TH1D("6_HistTwoTofMatchedTracksAbsEta", ";#eta;count", 50, 0, 2);
    TH1D* HistTwoTofMatchedTracksPt = new TH1D("6_HistTwoTofMatchedTracksPt", ";p_{T} [GeV];count", 50, 0, 3.5);
    TH1D* HistTwoTofMatchedTracksNfit = new TH1D("7_HistTwoTofMatchedTracksNfit", ";N_{fit};count", 61, -0.5, 60.5);
    TH1D* HistTwoTofMatchedTracksNdEdx = new TH1D("7_HistTwoTofMatchedTracksNdEdx", ";N_{dE/dx}; count", 61, -0.5, 60.5);
    TH1D* HistTwoTofMatchedTracksDcaXY = new TH1D("8_HistTwoTofMatchedTracksDcaXY", ";DCA_{XY} [cm]; count", 60, 0, 4);
    TH1D* HistTwoTofMatchedTracksAbsDcaZ = new TH1D("8_HistTwoTofMatchedTracksAbsDcaZ", ";DCA_{Z} [cm]; count ", 60, 0, 4);
    TH1D* HistAbsDeltaDcaZTwoTofMatched  = new TH1D("9_HistAbsDeltaDcaZTwoTofMatched", ";#Delta DCA_{Z} [cm]; count", 60, 0, 3);
    TH1D* HistRpTrackPointPlanesUsed = new TH1D("10_HistRpTrackPointPlanesUsed", ";Number Of Planes per TrackPoint;count", 5, -0.5, 4.5);
    TH1D* HistRpTrackBranch = new TH1D("11_HistRpTrackBranch", ";RP track branch;count", 4, -0.5, 3.5);
    TH1D* HistRpTrackPointStation = new TH1D("12_HistRpTrackPointStation", ";RP track station; count", 8, -0.5, 7.5);
    TH1D* HistThetaRpX =  new TH1D("13_HistThetaRpX", ";#theta_{x}^{RP}; count ", 200, -0.05, 0.05);
    TH1D* HistThetaRpY =  new TH1D("13_HistThetaRpY", ";#theta_{y}^{RP}; count", 200, -0.05, 0.05);
    TH1D* HistThetaX =  new TH1D("13_HistThetaX", ";#theta_x;count", 200, -0.05, 0.05);
    TH1D* HistThetaY =  new TH1D("13_HistThetaY", ";#theta_y; count", 200, -0.05, 0.05);
    TH1D* HistDiffThetaX =  new TH1D("13_HistDiffThetaX", ";#theta_{x}^{RP} - #theta_{x};count", 200,-0.01, 0.01);
    TH1D* HistDiffThetaY =  new TH1D("13_HistDiffThetaY", ";#theta_{y}^{RP} - #theta_{y};count", 200,-0.0002, 0.0002);
    TH2D* HistRpTrackPxPy= new TH2D("14_HistRpTrackPxPy", ";p_{x} [GeV]; p_{y} [GeV]", 100, -1.5, 1.5, 100, -1.5, 1.5);	


	ROOT::TThreadedObject<TH2D> DoublePxPyEastBeforeCuts;
	new(&DoublePxPyEastBeforeCuts)ROOT::TThreadedObject<TH2D>("DoublePxPyEastBeforeCuts","DoublePxPyEastBeforeCuts", 100, -1.5, 1.5, 100, -1.5, 1.5);	

    static StUPCEvent *upcEvt = 0x0;
    static StRPEvent *rpEvt = 0x0;
    chain.SetBranchAddress("mUPCEvent", &upcEvt);
    chain.SetBranchAddress("mRPEvent", &rpEvt);

    vector<TLorentzVector> posPion;
    vector<TLorentzVector> negPion;
    TLorentzVector trackVector;
	TVector3 pVector;
	
    double beamPositionX, beamPositionY;
    double primVertexPosX, primVertexPosY, primVertexPosErrX, primVertexPosErrY;
    double distVertexBeamX, distVertexBeamY, distR, significance;

    for (Long64_t i = 0; i < chain.GetEntries(); ++i) 
    {
        // cout << i << "/" <<  chain.GetEntries()  << endl;
        chain.GetEntry(i);
        HistCutFlow->Fill(0); //all data
        if (i%100000 == 0)
        {
            cout << i << "/" <<  chain.GetEntries()  << endl;
        }

        // check CEP trigger
        for (unsigned long int i = 0; i < triggerID.size(); i++)
        {
            if (upcEvt->isTrigger(triggerID[i]))
            {
                HistTriggers->Fill(i);
            }
        }
        
        bool cepTrigger = 0;
        for (unsigned long int i = 0; i < triggerCEP.size(); i++)
        {
            if (upcEvt->isTrigger(triggerCEP[i]))
            {
                cepTrigger = 1;
            }
        }
        if (cepTrigger == 0) {continue;} 
        HistCutFlow->Fill(1); //check CEP triggers

        int nFillNumber = upcEvt->getFillNumber();
		beamPositionX = FindPosition(nFillNumber, fillNumberWithPosition[0], fillNumberWithPosition[1], fillNumberWithPosition[2], fillNumberWithPosition[3], fillNumberWithPosition[4])[0];
		beamPositionY = FindPosition(nFillNumber, fillNumberWithPosition[0], fillNumberWithPosition[1], fillNumberWithPosition[2], fillNumberWithPosition[3], fillNumberWithPosition[4])[1];			  
    
    	//SC1: Exactly one primary vertex with TPC track(s) matched with hits in TOF is found
		vector<Int_t> primaryVertices;
		Int_t VertexId;
        for (Int_t i = 0; i < upcEvt->getNumberOfTracks(); i++)
        {
            VertexId = upcEvt->getTrack(i)->getVertexId();
            if(upcEvt->getTrack(i)->getFlag(StUPCTrack::kTof) && upcEvt->getTrack(i)->getFlag(StUPCTrack::kPrimary) && !(find(primaryVertices.begin(), primaryVertices.end(), VertexId)!=primaryVertices.end()))
            {
                primaryVertices.push_back(VertexId);
            }
        }
        HistNumPrimaryVertices->Fill(primaryVertices.size());
            				
        if(primaryVertices.size()!=1) //if the number of primary vertices greater than one - continue
        {
            continue;
        }		
        HistCutFlow->Fill(2); //one vertex
           		
 		VertexId = primaryVertices[0];
        primVertexPosX = upcEvt->getVertexId(VertexId)->getPosX();
        primVertexPosY = upcEvt->getVertexId(VertexId)->getPosY();
        primVertexPosErrX = upcEvt->getVertexId(VertexId)->getErrX();
        primVertexPosErrY = upcEvt->getVertexId(VertexId)->getErrY();
        distVertexBeamX = primVertexPosX - beamPositionX;
        distVertexBeamY = primVertexPosY - beamPositionY;
        distR = sqrt(pow(distVertexBeamX,2) + pow(distVertexBeamY,2));
   		significance = sqrt(pow(distVertexBeamX,2)/primVertexPosErrX + pow(distVertexBeamY,2)/primVertexPosErrY);
	            		
        ProfileDistVertexBeamX->Fill(nFillNumber, distVertexBeamX); 
        ProfileDistVertexBeamY->Fill(nFillNumber, distVertexBeamY); 
        HistR->Fill(distR);
   		HistSig->Fill(significance);     

		//SC2: TPC vertex from SC1 is placed within |zvtx| < 80 cm
        HistPrimaryVertexAbsPosZ->Fill(abs(upcEvt->getVertexId(VertexId)->getPosZ()));
		if(abs(upcEvt->getVertexId(VertexId)->getPosZ())>=80)
		{
			continue;
        }   
        HistCutFlow->Fill(3); //80 cm
            			
		vector <StUPCTrack*> tracksWithTofHit;
        for (Int_t i = 0; i<upcEvt->getNumberOfTracks(); i++)
		{
			if(upcEvt->getTrack(i)->getFlag(StUPCTrack::kTof))
			{
                tracksWithTofHit.push_back(upcEvt->getTrack(i));
			}
		}       		
          	
		//SC3.1:  Exactly two TOF-matched primary tracks
        HistNumTofMatchedTracks->Fill(tracksWithTofHit.size());
		bool SC3_1 = 0;
		if (tracksWithTofHit.size() == 2)
		{
			SC3_1 = 1;
		}
		if (SC3_1 == 0) {continue;} 
		HistCutFlow->Fill(4); //two tracks

		//SC3.2:  Tracks are of opposied signs		
        int sumCharge = 0;
        for (unsigned long int i = 0; i < tracksWithTofHit.size(); i++)
        {
            HistTwoTofMatchedTracksCharge->Fill(tracksWithTofHit[i]->getCharge());
            sumCharge+=tracksWithTofHit[i]->getCharge();
        }
        if (sumCharge != 0)    {continue;}
		HistCutFlow->Fill(5); 
			
		//SC3.3: Both tracks are contained within the kinematic range	
        bool SC3_3 = 1;
        for (unsigned long int i = 0; i < tracksWithTofHit.size(); i++)
        {	
            HistTwoTofMatchedTracksAbsEta->Fill(abs(tracksWithTofHit[i]->getEta()));
            HistTwoTofMatchedTracksPt->Fill(tracksWithTofHit[i]->getPt()); 

            if (abs(tracksWithTofHit[i]->getEta()) >= 0.7 or tracksWithTofHit[i]->getPt() <= 0.2)
            {
                SC3_3 = 0;
            }
        }
        if (SC3_3 == 0) {continue;}
		HistCutFlow->Fill(6); //two tracks kinematic range

		//SC3.4: Associated global tracks satisfy quality criteria: N_{fit} >= 25 and N_{dE/dx} >= 15 
        bool SC3_4 = 1;
        for (unsigned long int i = 0; i < tracksWithTofHit.size(); i++)
        {		
            HistTwoTofMatchedTracksNfit->Fill(tracksWithTofHit[i]->getNhitsFit());
            HistTwoTofMatchedTracksNdEdx->Fill(tracksWithTofHit[i]->getNhitsDEdx());
            if (tracksWithTofHit[i]->getNhitsFit() < 25 or tracksWithTofHit[i]->getNhitsDEdx() < 15)
            {
                SC3_4 = 0;
            }
        }
        if (SC3_4 == 0) {continue;}
		HistCutFlow->Fill(7); //two tracks quality

		//SC3.5 Associated global tracks match well to the primary vertex DCA(R) < 1.5 |DCA(z)| < 1.0			
        bool SC3_5 = 1;
        for (unsigned long int i = 0; i < tracksWithTofHit.size(); i++)
        {	
            HistTwoTofMatchedTracksDcaXY->Fill(tracksWithTofHit[i]->getDcaXY());
            HistTwoTofMatchedTracksAbsDcaZ->Fill(abs(tracksWithTofHit[i]->getDcaZ()));	
            if (tracksWithTofHit[i]->getDcaXY() >= 1.5 or abs(tracksWithTofHit[i]->getDcaZ())>=1.0)
            {
                SC3_5 = 0;
            }
        }
        if (SC3_5 == 0) {continue;}
		HistCutFlow->Fill(8); //two tracks match primary vertex
			
		//SC3.6 Associated global tracks are close at the beamline: |dz0| < 2 I usec dist in DCA
        bool SC3_6 = 1;
        for (unsigned long int i = 0; i < tracksWithTofHit.size(); i++)
        {	
            for (unsigned long int j = i+1; j < tracksWithTofHit.size(); j++)
            {	
                if (j!=i)
                {
                    HistAbsDeltaDcaZTwoTofMatched->Fill(abs(tracksWithTofHit[i]->getDcaZ() - tracksWithTofHit[j]->getDcaZ()));
                    if (abs(tracksWithTofHit[i]->getDcaZ() - tracksWithTofHit[j]->getDcaZ()) >=2.0 )
                    {
                        SC3_6 = 0;
                    }
                }
            }
        }
        if (SC3_6 == 0) {continue;}
		HistCutFlow->Fill(9); 


		bool SC4_1 = 1;
		for(unsigned int i = 0; i < rpEvt->getNumberOfTrackPoints(); i++)
		{
            HistRpTrackPointPlanesUsed->Fill(rpEvt->getTrackPoint(i)->planesUsed());
			if(rpEvt->getTrackPoint(i)->planesUsed() <3)
			{
				SC4_1 = 0;		
				break;
			}
		}		
		if ( SC4_1 == 0) {continue;}	
		HistCutFlow->Fill(10); //at least three out of the four plane

        bool SC4_3 = 0;
		int numberOfTracksEast = 0;
		int numberOfTracksWest = 0;
		int indexWest = 0;
		int indexEast = 0;
		for(unsigned int i = 0; i < rpEvt->getNumberOfTracks(); i++)
		{
			int iBranch = rpEvt->getTrack(i)->branch();
            HistRpTrackBranch->Fill(iBranch);
			if (iBranch < 2) {numberOfTracksEast+=1; indexEast = i;}
		    else {numberOfTracksWest+=1; indexWest  = i;}
	    }
		if (numberOfTracksWest==1 and numberOfTracksEast == 1)
		{
			SC4_3 = 1;
		}
		if (SC4_3 == 0) {continue;}	
		HistCutFlow->Fill(11); 

        bool SC4_3_1 = 0;
		int numberOfTrackPointsEast = 0;
		int numberOfTrackPointsWest = 0;
		int indexStationWest = 0;
		int indexStationEast = 0;
        
		for(unsigned int i = 0; i < rpEvt->getNumberOfTrackPoints(); i++)
		{
            int iStation = static_cast<int>(rpEvt->getTrackPoint(i)->rpId());
            HistRpTrackPointStation->Fill(iStation);
            if (iStation < 4) {numberOfTrackPointsEast+=1; indexStationEast = i;}
		    else {numberOfTrackPointsWest+=1; indexStationWest  = i;}
        } 
		if (numberOfTrackPointsWest == 2 and numberOfTrackPointsEast == 2)
		{
			SC4_3_1 = 1;
		}
		if (SC4_3_1 == 0) {continue;}	
		HistCutFlow->Fill(12); 

        bool SC4_2 = 1;
	    for(unsigned int k = 0; k < rpEvt->getNumberOfTracks(); ++k)
		{
		    StUPCRpsTrack *trk = rpEvt->getTrack(k);
			trk->setEvent(rpEvt);
		    double thetaRPx = rpEvt->getTrack(k)->thetaRp(0);
		    double thetaRPy = rpEvt->getTrack(k)->thetaRp(1);  
		    double thetax = rpEvt->getTrack(k)->theta(0);
		    double thetay = rpEvt->getTrack(k)->theta(1);	
            HistThetaRpX->Fill(thetaRPx);
            HistThetaRpY->Fill(thetaRPy);
            HistThetaX->Fill(thetax);
            HistThetaY->Fill(thetay);
            HistDiffThetaX->Fill(thetaRPx-thetax);
            HistDiffThetaY->Fill(thetaRPy-thetay);

			if ((thetaRPx-thetax) <= -2/10e3 or (thetaRPx-thetax) >= 4/10e3 or (thetaRPy-thetay) <= -2/10e3 or (thetaRPy-thetay) >= 2/10e3)
			{
				SC4_2 = 0;	 				
	      	}
        }		
		if (SC4_2 == 0) {continue;}
		HistCutFlow->Fill(13); // difference between measured and recondstructed momenta
				
    	bool SC4_4 = 1;
	    for(unsigned int k = 0; k < rpEvt->getNumberOfTracks(); ++k)
	    {			
			StUPCRpsTrack *trk = rpEvt->getTrack(k);
		    trk->setEvent(rpEvt);
		    pVector = trk->pVec();
			HistRpTrackPxPy->Fill(pVector.X(), pVector.Y());

			if (abs(pVector.Y()) >= 0.8 or abs(pVector.Y()) <= 0.4 or pVector.X() <= -0.27 or (pow(pVector.X()+0.6,2) + pow(pVector.Y(), 2)  >=  pow(1.1,2)))
			{
				SC4_4 = 0;
			}  	
		}
		if (SC4_4 == 0) {continue;}    
		HistCutFlow->Fill(14); //momenta
        

        for (unsigned long int i = 0; i < tracksWithTofHit.size(); i++)
        {
            if (abs(tracksWithTofHit[i]->getNSigmasTPCPion())<3)
            {
                tracksWithTofHit[i]->getLorentzVector(trackVector, massPion);
            }
            else
            {
                break;
            }

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

	    HistCutFlow->Fill(15); //two tracks are good candidate for pion

        TLorentzVector tempSum = {0, 0, 0, 0};
        tempSum = posPion[0] + negPion[0];
        HistInvMassPiPi->Fill(tempSum.M());

        posPion.clear();
        negPion.clear();           
    }

    for (int i = 0; i < 17; i++)
    {
        char constNumberStr[20];
        sprintf(constNumberStr, "%d", triggerID[i]);
        HistTriggers->GetXaxis()->SetBinLabel(i+1,constNumberStr);
    }

    TFile *outfile = TFile::Open(argv[2], "recreate"); 
    HistInvMassPiPi->Write();
    ProfileDistVertexBeamX->Write();
    HistR->Write();
    HistSig->Write();
    HistCutFlow->Write();
    ProfileDistVertexBeamY->Write();
    HistTriggers->Write();
    HistNumPrimaryVertices->Write();
    HistPrimaryVertexAbsPosZ->Write();
    HistNumTofMatchedTracks->Write();
    HistTwoTofMatchedTracksCharge->Write();
    HistTwoTofMatchedTracksAbsEta->Write();
    HistTwoTofMatchedTracksPt->Write();
    HistTwoTofMatchedTracksNfit->Write();
    HistTwoTofMatchedTracksNdEdx->Write();
    HistTwoTofMatchedTracksDcaXY->Write();
    HistTwoTofMatchedTracksAbsDcaZ->Write();
    HistAbsDeltaDcaZTwoTofMatched->Write();
    HistRpTrackPointPlanesUsed->Write();
    HistRpTrackBranch->Write();
    HistRpTrackPointStation->Write();
    HistThetaRpX->Write();
    HistThetaRpY->Write();
    HistThetaX->Write();
    HistThetaY->Write();
    HistDiffThetaX->Write();
    HistDiffThetaY->Write();
    HistRpTrackPxPy->Write();
    outfile->Close();

    return 0;
}
