#include "Includes.h"

using namespace std;

int main(int argc, char** argv)  
{
    // histograms 
    TH1D* HistCutFlow = new TH1D("HistCutFlow", ";;count", 16, -0.5, 15.5);

    TH1D* HistTriggers = new TH1D("1_HistTriggers", ";trigger number; count", 17, -0.5, 16.5);
    TH1D* HistRpTrackBranch = new TH1D("2_HistRpTrackBranch", ";RP track branch; count", 4, -0.5, 3.5);
    TH1D* HistRpTrackPointPlanesUsed = new TH1D("3_HistRpTrackPointPlanesUsed", ";Number Of Planes per TrackPoint; count", 5, -0.5, 4.5);
    TH2D* HistRpTrackPxPy= new TH2D("4_HistRpTrackPxPy", ";p_{x} [GeV]; p_{y} [GeV]", 100, -1.5, 1.5, 100, -1.5, 1.5);	
    TH1D* HistNumTofMatchedTracks = new TH1D("5_HistNumTofMatchedTracks", ";number of TOF matched tracks; count", 16, -0.5, 15.5);    

    TH1D* HistTriggersPostSelection = new TH1D("1_HistTriggersPs", ";trigger number; count", 17, -0.5, 16.5);
    TH1D* HistRpTrackBranchPostSelection = new TH1D("2_HistRpTrackBranchPs", ";RP track branch; count", 4, -0.5, 3.5);
    TH1D* HistRpTrackPointPlanesUsedPostSelection = new TH1D("3_HistRpTrackPointPlanesUsedPs", ";Number Of Planes per TrackPoint; count", 5, -0.5, 4.5);
    TH2D* HistRpTrackPxPyPostSelection = new TH2D("4_HistRpTrackPxPyPs", ";p_{x} [GeV]; p_{y} [GeV]", 100, -1.5, 1.5, 100, -1.5, 1.5);	
    TH1D* HistNumTofMatchedTracksPostSelection = new TH1D("5_HistNumTofMatchedTracksPs", ";number of TOF matched tracks; count", 16, -0.5, 15.5);    

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

    static StUPCEvent *upcEvt = 0x0;
    static StRPEvent *rpEvt = 0x0;
    chain->SetBranchAddress("mUPCEvent", &upcEvt);
    chain->SetBranchAddress("mRPEvent", &rpEvt);

    TTree* mUPCTree = new TTree("mUPCTree", "mUPCTree");
    mUPCTree->Branch("mUPCEvent", &upcEvt);
    mUPCTree->Branch("mRPEvent", &rpEvt);

	TVector3 pVector;

    int iCutFlow;

    const vector <int> triggerID = {570209, 570219, 570229, 570701, 570702, 570703, 570704, 570705, 570709, 570711, 570712, 570719, 590701, 590703, 590705, 590708, 590709};
    const vector <int> triggerCEP = {570701, 570705, 570711, 590701, 590705, 590708};

    for (Long64_t i = 0; i < chain->GetEntries(); ++i) 
    {

        iCutFlow = 0;

        chain->GetEntry(i);
        HistCutFlow->Fill(iCutFlow);  // cutflow: all data

        if (i%1000000 == 0)  { cout << i << "/" <<  chain->GetEntries()  << endl;}

        // triggers - one event can have several triggers
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
        HistCutFlow->Fill(iCutFlow); // cutflow: CEP trigger

        // post-selection triggers
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
  		HistCutFlow->Fill(iCutFlow); // cutflow: two RP tracks on opposite sides of the detector

        // post-selection RP branches
		for(UInt_t i = 0; i < rpEvt->getNumberOfTracks(); i++)
        {
			int iBranch = rpEvt->getTrack(i)->branch();
            HistRpTrackBranchPostSelection->Fill(iBranch);
        }

        // RP trackpoints: # of used planes - at least 3 out of 4
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
  		HistCutFlow->Fill(iCutFlow); // cutflow: 3 out of 4 silicon planes were used for the reconstruction

        // post-selection: # of planes per tracpoint
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
			HistRpTrackPxPy->Fill(pVector.X(), pVector.Y());

			if ((abs(pVector.Y()) >= 0.8 or abs(pVector.Y()) <= 0.4) or (pVector.X() <= -0.27) or (pow(pVector.X()+0.6,2) + pow(pVector.Y(), 2)  >=  pow(1.25,2)))
			{
				withinFiducialRegion = 0;
			}  	
		}
		if (withinFiducialRegion == 0) {continue;}    
        iCutFlow+=1;
		HistCutFlow->Fill(iCutFlow); // cutflow: fiducial region

        // post-seletcion: fiducial region
	    for(unsigned int k = 0; k < rpEvt->getNumberOfTracks(); ++k)
	    {		
			StUPCRpsTrack *trk = rpEvt->getTrack(k);
		    trk->setEvent(rpEvt);
		    pVector = trk->pVec();
			HistRpTrackPxPyPostSelection->Fill(pVector.X(), pVector.Y());
        }
 
		// vector for TOF matched tracks
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
		if (tracksWithTofHit.size() == 4 or tracksWithTofHit.size() == 3)
		{
			isValidNumberOfTofMatchedTracks = 1;
		}
		if (isValidNumberOfTofMatchedTracks == 0) {continue;} 
        iCutFlow+=1;
		HistCutFlow->Fill(iCutFlow); // cutflow: 3 or 4 tpc tracks matched w/ TOF       
        HistNumTofMatchedTracksPostSelection->Fill(tracksWithTofHit.size()); // post-selection: # events with 3 or 4 TOF-matched tracks
    
        // set branch adress
        mUPCTree->SetBranchAddress("mUPCEvent", &upcEvt);
        mUPCTree->SetBranchAddress("mRPEvent", &rpEvt);
        
        // fill tree
        mUPCTree->Fill();
    }

    // save root files that passed selection
    TFile outputFile(argv[3], "RECREATE");
    mUPCTree->Write();
    outputFile.Close();

    // save histograms
    TFile *outfile = TFile::Open(argv[2], "recreate"); 
    HistCutFlow->Write(); 
    HistTriggers->Write();
    HistRpTrackBranch->Write();
    HistRpTrackPointPlanesUsed->Write();
    HistRpTrackPxPy->Write();
    HistNumTofMatchedTracks->Write();
    HistTriggersPostSelection->Write();
    HistRpTrackBranchPostSelection->Write();
    HistRpTrackPointPlanesUsedPostSelection->Write();
    HistRpTrackPxPyPostSelection->Write();
    HistNumTofMatchedTracksPostSelection->Write();
    outfile->Close();

    return 0;
}
