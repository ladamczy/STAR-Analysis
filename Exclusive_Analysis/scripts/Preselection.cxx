#include <iostream>
#include <string>    
#include <utility>
#include <sstream> 
#include <algorithm> 
#include <stdio.h> 
#include <stdlib.h> 
#include <vector> 
#include <fstream> 
#include <cmath> 
#include <cstdlib>
#include <sys/stat.h>
#include <iterator>
#include <ostream>
#include <iomanip>
#include <stdexcept>
#include <limits>
#include "TROOT.h"
#include "TSystem.h"
#include "TThread.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include "TProfile.h"
#include <TH2.h> 
#include <TF1.h> 
#include <TF2.h> 
#include <THStack.h> 
#include <TStyle.h> 
#include <TGraph.h> 
#include <TGraph2D.h> 
#include <TGraphErrors.h> 
#include <TCanvas.h> 
#include <TLegend.h> 
#include <TGaxis.h> 
#include <TString.h> 
#include <TColor.h> 
#include <TLine.h> 
#include <TExec.h> 
#include <TFitResultPtr.h> 
#include <TFitResult.h> 
#include <TLatex.h> 
#include <TMath.h>
#include <TLorentzVector.h>
#include <ROOT/TThreadedObject.hxx>
#include <TTreeReader.h>
#include <ROOT/TTreeProcessorMT.hxx>
#include "StRPEvent.h"
#include "StUPCRpsTrack.h"
#include "StUPCRpsTrackPoint.h"
#include "StUPCEvent.h"
#include "StUPCTrack.h"
#include "StUPCBemcCluster.h"
#include "StUPCVertex.h"
#include "StUPCTofHit.h"
#include <Afterburner.h>
#include "StUPCV0.h"
#include <StRPEvent.h>
#include <map>

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

    TH1D* HistTriggersPostSelection = new TH1D("1_HistTriggersPostSelection", ";trigger number; count", 17, -0.5, 16.5);
    TH1D* HistRpTrackBranchPostSelection = new TH1D("2_HistRpTrackBranchPostSelection", ";RP track branch; count", 4, -0.5, 3.5);
    TH1D* HistRpTrackPointPlanesUsedPostSelection = new TH1D("3_HistRpTrackPointPlanesUsedPostSelection", ";Number Of Planes per TrackPoint; count", 5, -0.5, 4.5);
    TH2D* HistRpTrackPxPyPostSelection = new TH2D("4_HistRpTrackPxPyPostSelection", ";p_{x} [GeV]; p_{y} [GeV]", 100, -1.5, 1.5, 100, -1.5, 1.5);	
    TH1D* HistNumTofMatchedTracksPostSelection = new TH1D("5_HistNumTofMatchedTracksPostSelection", ";number of TOF matched tracks; count", 16, -0.5, 15.5);    

    TH1D* HistTriggersN1 = new TH1D("1_HistTriggersN1", ";trigger number; count", 17, -0.5, 16.5);
    TH1D* HistRpTrackBranchN1 = new TH1D("2_HistRpTrackBranchN1", ";RP track branch; count", 4, -0.5, 3.5);
    TH1D* HistRpTrackPointPlanesUsedN1 = new TH1D("3_HistRpTrackPointPlanesUsedN1", ";Number Of Planes per TrackPoint; count", 5, -0.5, 4.5);
    TH2D* HistRpTrackPxPyN1 = new TH2D("4_HistRpTrackPxPyN1", ";p_{x} [GeV]; p_{y} [GeV]", 100, -1.5, 1.5, 100, -1.5, 1.5);	
   

    ifstream inputFilePathList(argv[1]);
    if (!inputFilePathList) 
    {
        cerr << "Failed to open input file." << std::endl;
        return 1;
    }

    TChain *chain = new TChain("mUPCTree"); 

    string inputFileName;

    std::string baseDir = "/run/media/pysia/One_Touch/";
    std::vector<std::string> subDirs = {"star_data_1", "star_data_2", "star_data_3", "star_data_4", "star_data_5/star_data_6"};  // Zakładając, że znamy nazwy podkatalogów
std::string filePath;
    while (std::getline(inputFilePathList, inputFileName))
    {    
        
        std::string targetFile = inputFileName;

        for (const auto& subDir : subDirs) {
         filePath = baseDir + subDir + "/" + targetFile;
        TFile* file = TFile::Open(filePath.c_str());
        if (file && file->IsOpen()) {
            std::cout << "Pomyślnie otworzono plik: " << filePath << std::endl;
     
            file->Close();
            delete file;
            break;  // Zakładamy, że plik jest tylko w jednym podkatalogu
        } else {
            std::cerr << "Nie można otworzyć pliku: " << filePath << std::endl;
        }
    }
        chain->Add(filePath.c_str());
    }

    inputFilePathList.close();

    static StUPCEvent *upcEvt = 0x0;
    static StRPEvent *rpEvt = 0x0;
    static StRPEvent  *correctedRpEvent = 0x0;

    chain->SetBranchAddress("mUPCEvent", &upcEvt);
    chain->SetBranchAddress("mRPEvent", &rpEvt);

    LoadOffsetFile("../share/OffSetsCorrectionsRun17.list", mCorrection);

    TTree* mUPCTree = new TTree("mUPCTree", "mUPCTree");
    mUPCTree->Branch("mUPCEvent", &upcEvt);
    mUPCTree->Branch("correctedRpEvent", &correctedRpEvent);

	TVector3 pVector;

    int iCutFlow;

    const vector <int> triggerID = {570209, 570219, 570229, 570701, 570702, 570703, 570704, 570705, 570709, 570711, 570712, 570719, 590701, 590703, 590705, 590708, 590709};
    const vector <int> triggerCEP = {570701, 570705, 570711};

    for (Long64_t i = 0; i < chain->GetEntries(); ++i) 
    {
        if (i%1000000 == 0)  { cout << i << "/" <<  chain->GetEntries()  << endl;}
        chain->GetEntry(i);

     
        correctedRpEvent = new StRPEvent(*rpEvt);
        correctedRpEvent->clearEvent();
        runAfterburner(rpEvt, correctedRpEvent, upcEvt->getRunNumber());

        // selection
        // CEP Triggers
        // two RP tracks
        // 3 out of 4 planes
        // fiducial region

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

        // two RP tracks on the opposite sides of the detector
        bool hasTwoRpTracksOppositeSide = 0;
		int numberOfTracksEast = 0;
		int numberOfTracksWest = 0;
		int indexWest = 0;
		int indexEast = 0;

		for(UInt_t i = 0; i < correctedRpEvent->getNumberOfTracks(); i++)
		{
			int iBranch = correctedRpEvent->getTrack(i)->branch();
            HistRpTrackBranch->Fill(iBranch);
			if (iBranch < 2) {numberOfTracksEast+=1; indexEast = i;}
		    else {numberOfTracksWest+=1; indexWest  = i;}
	    }
		if (numberOfTracksWest == 1 and numberOfTracksEast == 1)
		{
			hasTwoRpTracksOppositeSide = 1;

		}

        // RP trackpoints: # of used planes - at least 3 out of 4
    	bool useAtLeastThreePlanes = 1;
		for(UInt_t i = 0; i < correctedRpEvent->getNumberOfTrackPoints(); i++)
		{
            HistRpTrackPointPlanesUsed->Fill(correctedRpEvent->getTrackPoint(i)->planesUsed());
			if(correctedRpEvent->getTrackPoint(i)->planesUsed() < 3)
			{
				useAtLeastThreePlanes = 0;		
				break;
			}
		}		

        // fiducial region
    	bool withinFiducialRegion = 1;
	    for(unsigned int k = 0; k < correctedRpEvent->getNumberOfTracks(); ++k)
	    {		
			StUPCRpsTrack *trk = correctedRpEvent->getTrack(k);
		    trk->setEvent(correctedRpEvent);
            double thetaX = trk->thetaRp(0);
            double thetaY = trk->thetaRp(1);
            double thetaZ = trk->thetaRp(2);
            double px = 254.867*thetaX;
            double py = 254.867*thetaY;
            double pz = 254.867*thetaZ;     

			HistRpTrackPxPy->Fill(px, py);

			if ((abs(py) >= 0.8 or abs(py) <= 0.4) or (px <= -0.27) or (pow(px+0.6,2) + pow(py, 2)  >=  (1.25)))
			{
				withinFiducialRegion = 0;
			}  	
		}


		// vector for TOF matched tracks
		vector <StUPCTrack*> tracksWithTofHit;
        vector <StUPCTrack*> tracksWithTofHistClassA;
        vector <StUPCTrack*> tracksWithTofHistClassB;
              
        for (Int_t i = 0; i<upcEvt->getNumberOfTracks(); i++)
		{
			if(upcEvt->getTrack(i)->getFlag(StUPCTrack::kTof))
			{
                tracksWithTofHit.push_back(upcEvt->getTrack(i));
			}

			if(upcEvt->getTrack(i)->getFlag(StUPCTrack::kTof) and !(upcEvt->getTrack(i)->getFlag(StUPCTrack::kV0)) and !(upcEvt->getTrack(i)->getFlag(StUPCTrack::kCEP))  )
			{
                tracksWithTofHistClassA.push_back(upcEvt->getTrack(i));
			}
            
			if(upcEvt->getTrack(i)->getFlag(StUPCTrack::kTof) and upcEvt->getTrack(i)->getFlag(StUPCTrack::kV0))
			{
                tracksWithTofHistClassB.push_back(upcEvt->getTrack(i));
			}
		}       		

		// two TOF-matched primary tracks
        HistNumTofMatchedTracks->Fill(tracksWithTofHit.size());
		bool isValidNumberOfTofMatchedTracks = 0;
		bool isValidNumberOfTofMatchedTracksClassA = 0;
        bool isValidNumberOfTofMatchedTracksClassB = 0;
        
              
        
		if (tracksWithTofHit.size() >=2 )
		{
			isValidNumberOfTofMatchedTracks = 1;
		}

        
		if (tracksWithTofHistClassA.size() >=2 )
		{
			isValidNumberOfTofMatchedTracksClassA = 1;
		}

        
		if (tracksWithTofHistClassB.size() >=2 )
		{
			isValidNumberOfTofMatchedTracksClassB = 1;
		}




        if (hasTwoRpTracksOppositeSide and useAtLeastThreePlanes and withinFiducialRegion)
        {
            for (unsigned long int i = 0; i < triggerID.size(); i++)
            {
                if (upcEvt->isTrigger(triggerID[i]))
                {
                    HistTriggersN1->Fill(i);
                }
            }
        }

        if (isCepTrigger and useAtLeastThreePlanes and withinFiducialRegion)
        {
            for(UInt_t i = 0; i < correctedRpEvent->getNumberOfTracks(); i++)
            {
                int iBranch = correctedRpEvent->getTrack(i)->branch();
                HistRpTrackBranchN1->Fill(iBranch);
            }
        }

        if (isCepTrigger and hasTwoRpTracksOppositeSide and withinFiducialRegion)
        {
            for(UInt_t i = 0; i < correctedRpEvent->getNumberOfTrackPoints(); i++)
            {
                HistRpTrackPointPlanesUsedN1->Fill(correctedRpEvent->getTrackPoint(i)->planesUsed());
            }	
        }

        if (isCepTrigger and hasTwoRpTracksOppositeSide and useAtLeastThreePlanes)
        {

            for(unsigned int k = 0; k < correctedRpEvent->getNumberOfTracks(); ++k)
            {		
                StUPCRpsTrack *trk = correctedRpEvent->getTrack(k);
                trk->setEvent(correctedRpEvent);
                double thetaX = trk->thetaRp(0);
                double thetaY = trk->thetaRp(1);
                double thetaZ = trk->thetaRp(2);
                double px = 254.867*thetaX;
                double py = 254.867*thetaY;
                double pz = 254.867*thetaZ;     
                HistRpTrackPxPyN1->Fill(px, py);
            }

        }


        if (isCepTrigger and hasTwoRpTracksOppositeSide and useAtLeastThreePlanes and withinFiducialRegion and (isValidNumberOfTofMatchedTracksClassA or isValidNumberOfTofMatchedTracksClassB) )
        {
            // set branch adress
            mUPCTree->SetBranchAddress("mUPCEvent", &upcEvt);
            mUPCTree->SetBranchAddress("correctedRpEvent", &correctedRpEvent);
            
            // fill tree
            mUPCTree->Fill();

            for (unsigned long int i = 0; i < triggerID.size(); i++)
            {
                if (upcEvt->isTrigger(triggerID[i]))
                {
                    HistTriggersPostSelection->Fill(i);
                }
            }
            
            for(UInt_t i = 0; i < correctedRpEvent->getNumberOfTracks(); i++)
            {
                int iBranch = correctedRpEvent->getTrack(i)->branch();
                HistRpTrackBranchPostSelection->Fill(iBranch);
            }

            for(UInt_t i = 0; i < correctedRpEvent->getNumberOfTrackPoints(); i++)
            {
                HistRpTrackPointPlanesUsedPostSelection->Fill(correctedRpEvent->getTrackPoint(i)->planesUsed());
            }	

            for(unsigned int k = 0; k < correctedRpEvent->getNumberOfTracks(); ++k)
            {		
                StUPCRpsTrack *trk = correctedRpEvent->getTrack(k);
                trk->setEvent(correctedRpEvent);
                double thetaX = trk->thetaRp(0);
                double thetaY = trk->thetaRp(1);
                double thetaZ = trk->thetaRp(2);
                double px = 254.867*thetaX;
                double py = 254.867*thetaY;
                double pz = 254.867*thetaZ;     
                HistRpTrackPxPyPostSelection->Fill(px, py);
            }
            HistNumTofMatchedTracksPostSelection->Fill(tracksWithTofHit.size());
        }

        delete correctedRpEvent;
    
    }
    // save root files that passed selection
    TFile outputFile(argv[3], "RECREATE");
    mUPCTree->Write();
    outputFile.Close();

    // save histograms
    TFile *outfile = TFile::Open(argv[2], "recreate"); 

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

    HistTriggersN1->Write();
    HistRpTrackBranchN1->Write();
    HistRpTrackPointPlanesUsedN1->Write();
    HistRpTrackPxPyN1->Write();

    outfile->Close();

    return 0;
}
