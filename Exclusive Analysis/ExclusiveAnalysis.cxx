// c++ headers
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

// ROOT headers
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

// picoDst headers
#include "StRPEvent.h"
#include "StUPCRpsTrack.h"
#include "StUPCRpsTrackPoint.h"
#include "StUPCEvent.h"
#include "StUPCTrack.h"
#include "StUPCBemcCluster.h"
#include "StUPCVertex.h"
#include "StUPCTofHit.h"

#include "MatchFillPosition.h"
#include "ReadFillPositionFile.h"

using namespace std;

enum {kAll = 1, kCPT,  kRP, kOneVertex, kTPCTOF, kTotQ, kMax};
enum SIDE {E = 0, East = 0, W = 1, West = 1, nSides};
enum PARTICLES {Pion = 0, Kaon = 1, Proton = 2, nParticles};
enum BRANCH_ID {EU, ED, WU, WD, nBranches };
enum RP_ID {E1U, E1D, E2U, E2D, W1U, W1D, W2U, W2D, nRomanPots};

const double particleMass[nParticles] = {0.13957, 0.497611, 0.93827}; // pion, kaon, proton in GeV /c^2 
const int nTriggers = 17;
const int triggerID[] = {570209, 570219, 570229, 570701, 570702, 570703, 570704, 570705, 570709, 570711, 570712, 570719, 590701, 590703, 590705, 590708, 590709};
const int CEPtriggers[] = {570701, 570705, 570711, 590701, 590705, 590708};

bool ConnectInput(int argc, char** argv, TChain* fileChain);
bool CheckTriggers(StUPCEvent* localupcEvt);
long long GetFileSize(string filename);

int main(int argc, char** argv) 
{
	int nthreads = 2;
	if(argc == 4)
	{
		nthreads = atoi(argv[3]);
	}
	cout << "Program is running on " << nthreads << " threads" << endl;
	ROOT::EnableThreadSafety();
	ROOT::EnableImplicitMT(nthreads); //turn on multicore processing
	ROOT::EnableThreadSafety();
	   
	TChain* upcChain = new TChain("mUPCTree");    //chain with files to iterate through
	if(!ConnectInput(argc, argv, upcChain))
	{
	        cout << "Wrong input parameters..." << endl; 
	        return 1;
	}
	const string& outputFolder = argv[2];
	
	//Data reading
	string file = "Run7PolarizationWithPosition.csv";
	vector <vector<double>> Data = ReadFillPositionData(file);
	
	ROOT::TThreadedObject<TH1D> InvMassHist;
        new(&InvMassHist)ROOT::TThreadedObject<TH1D>("InvMassHist","InvMassHist", 100, 0.42, 0.56);
        
	ROOT::TThreadedObject<TH1D> InvMassHistR0_01;
        new(&InvMassHistR0_01)ROOT::TThreadedObject<TH1D>("InvMassHistR0_01","InvMassHistR0_01", 100, 0.42, 0.56);

	ROOT::TThreadedObject<TH1D> InvMassHistR01_025;
        new(&InvMassHistR01_025)ROOT::TThreadedObject<TH1D>("InvMassHistR01_025","InvMassHistR01_025", 100, 0.42, 0.56);
        
	ROOT::TThreadedObject<TH1D> InvMassHistR025;
        new(&InvMassHistR025)ROOT::TThreadedObject<TH1D>("InvMassHistR025","InvMassHistR025", 100, 0.42, 0.56);        
        
        
        

	ROOT::TThreadedObject<TProfile> DifferenceVertexPositionBeamLineX;
	new(&DifferenceVertexPositionBeamLineX)ROOT::TThreadedObject<TProfile>("DifferenceVertexPositionBeamLineX","DifferenceVertexPositionBeamLineX",639, 20511.5, 21150.5);

	ROOT::TThreadedObject<TProfile> DifferenceVertexPositionBeamLineY;
	new(&DifferenceVertexPositionBeamLineY)ROOT::TThreadedObject<TProfile>("DifferenceVertexPositionBeamLineY","DifferenceVertexPositionBeamLineY",639, 20511.5, 21150.5);

	ROOT::TThreadedObject<TH1D> PerpendicularDistanceVertexPositionBeamLine; new(&PerpendicularDistanceVertexPositionBeamLine)ROOT::TThreadedObject<TH1D>("PerpendicularDistanceVertexPositionBeamLine","PerpendicularDistanceVertexPositionBeamLine", 150, 0.0, 0.3);
	
	ROOT::TThreadedObject<TH1D> Significance;
	new(&Significance)ROOT::TThreadedObject<TH1D>("Significance","Significance", 30000, 0.0, 30);

	ROOT::TThreadedObject<TH1D> histL;
	new(&histL)ROOT::TThreadedObject<TH1D>("histL","histL", 30000, -40.0, 40.0);

	ROOT::TThreadedObject<TH1D> histL0_01;
	new(&histL0_01)ROOT::TThreadedObject<TH1D>("histL0_01","histL0_01", 30000, -40.0, 40.0);	

	ROOT::TThreadedObject<TH1D> histL01_025;
	new(&histL01_025)ROOT::TThreadedObject<TH1D>("histL01_025","histL01_025", 30000, -40.0, 40.0);	

	ROOT::TThreadedObject<TH1D> histL025;
	new(&histL025)ROOT::TThreadedObject<TH1D>("histL025","histL025", 30000, -40.0, 40.0);		
	

	auto myFunction = [&](TFile* myFile) 
	{

		TFile* tempFile = new TFile(myFile->GetTitle());
		TTree* tempTree = (TTree*)tempFile->Get("mUPCTree");
		
       	 	if(tempTree->GetEntries()==0)
       	 	{
			delete tempFile;
			return 0;
        	}	

		TTreeReader myReader(tempTree);
		TTreeReaderValue<StUPCEvent> StUPCEventInstance(myReader, "mUPCEvent");
		TTreeReaderValue<StRPEvent> StRPEventInstance(myReader, "mRPEvent");
		
	        myReader.Next();
	        
		std::shared_ptr<TH1D> InvMassHistLocal;
		InvMassHistLocal = InvMassHist.Get();

		std::shared_ptr<TH1D> InvMassHistR0_01Local;
		InvMassHistR0_01Local = InvMassHistR0_01.Get();

		std::shared_ptr<TH1D> InvMassHistR01_025Local;
		InvMassHistR01_025Local = InvMassHistR01_025.Get();

		std::shared_ptr<TH1D> InvMassHistR025Local;
		InvMassHistR025Local = InvMassHistR025.Get();
		
		std::shared_ptr<TProfile> DifferenceVertexPositionBeamLineXLocal;
		DifferenceVertexPositionBeamLineXLocal = DifferenceVertexPositionBeamLineX.Get();
		
		std::shared_ptr<TProfile> DifferenceVertexPositionBeamLineYLocal;
		DifferenceVertexPositionBeamLineYLocal = DifferenceVertexPositionBeamLineY.Get();
		
		std::shared_ptr<TH1D> PerpendicularDistanceVertexPositionBeamLineLocal;
		PerpendicularDistanceVertexPositionBeamLineLocal = PerpendicularDistanceVertexPositionBeamLine.Get();	
		
		std::shared_ptr<TH1D> SignificanceLocal;
		SignificanceLocal = Significance.Get();		
		
		std::shared_ptr<TH1D> histLLocal;		
		histLLocal = histL.Get();

		std::shared_ptr<TH1D> histLLocal0_01;		
		histLLocal0_01 = histL0_01.Get();		

		std::shared_ptr<TH1D> histLLocal01_025;		
		histLLocal01_025 = histL01_025.Get();	

		std::shared_ptr<TH1D> histLLocal025;		
		histLLocal025 = histL025.Get();	
								
	
		int filtered_entries = 0;
	
		StUPCEvent* tempUPCevent = StUPCEventInstance.Get();
		StRPEvent* tempRPevent = StRPEventInstance.Get();

		vector<TLorentzVector> posPion;
		vector<TLorentzVector> negPion;
		vector<TLorentzVector> posPionSuspicious;
		vector<TLorentzVector> negPionSuspicious;
		StUPCTrack* trackInstanceLocalUPC;
		TLorentzVector trackVector;
		TVector3 pVector;
	
		do{
	
			tempUPCevent = StUPCEventInstance.Get();
			tempRPevent = StRPEventInstance.Get();
		   	
		   	if(!CheckTriggers(StUPCEventInstance.Get())) //check triggers
		    	{
		    	    continue;
		    	}
   
			//Match Fill Number with Beam Position (to be converted into a function)
			int nFillNumber = tempUPCevent->getFillNumber();
	
			double beamPositionX = FindPosition(nFillNumber, Data[0], Data[1], Data[2], Data[3], Data[4])[0];
			double beamPositionY = FindPosition(nFillNumber, Data[0], Data[1], Data[2], Data[3], Data[4])[1];			 

			if (isnan(beamPositionX)) {continue;}
			
			
			//SC1: Exactly one primary vertex with TPC track(s) matched with hits in TOF is found
			vector<Int_t> primaryVertices;
		    	Int_t VertexId;
           		for (Int_t i = 0; i < tempUPCevent->getNumberOfTracks(); i++)
           		{
                		VertexId = tempUPCevent->getTrack(i)->getVertexId();
                		if(tempUPCevent->getTrack(i)->getFlag(StUPCTrack::kTof) && tempUPCevent->getTrack(i)->getFlag(StUPCTrack::kPrimary) && !(find(primaryVertices.begin(), primaryVertices.end(), VertexId)!=primaryVertices.end()))
                		{
                    			primaryVertices.push_back(VertexId);
                		}
            		}
            		
            		if(primaryVertices.size()!=1) //if the number of primary vertices greater than one - continue
            		{
                		continue;
            		}
     
  		
		   	VertexId = primaryVertices[0];
            		
            		Double_t  primaryVertexPositionX = tempUPCevent->getVertexId(VertexId)->getPosX();
            		Double_t  primaryVertexPositionY = tempUPCevent->getVertexId(VertexId)->getPosY();
  
              		Double_t  primaryVertexPositionXError = tempUPCevent->getVertexId(VertexId)->getErrX();
            		Double_t  primaryVertexPositionYError = tempUPCevent->getVertexId(VertexId)->getErrY();
            		        		
            		Double_t differenceVertexBeamX = primaryVertexPositionX - beamPositionX;
            		Double_t differenceVertexBeamY = primaryVertexPositionY - beamPositionY;
            		
            		Double_t R = sqrt(pow(differenceVertexBeamX,2) + pow(differenceVertexBeamY,2));
            		
            		DifferenceVertexPositionBeamLineXLocal->Fill(nFillNumber, differenceVertexBeamX); //Fill TProfile - distance: vertex - beam in xplane
            		DifferenceVertexPositionBeamLineYLocal->Fill(nFillNumber, differenceVertexBeamY); // Fill TPorgile - distnace: vertex - beam in yplane
            		PerpendicularDistanceVertexPositionBeamLineLocal->Fill(R); //Fill histogram R = sqrt(dx^2 + dy^2)
   
   
   			Double_t significance = sqrt(pow(differenceVertexBeamX,2)/primaryVertexPositionXError + pow(differenceVertexBeamY,2)/primaryVertexPositionYError);
   			SignificanceLocal->Fill(significance);   
   
           		//SC2: TPC vertex from SC1 is placed within |zvtx| < 80 cm
           		if(abs(tempUPCevent->getVertexId(VertexId)->getPosZ())>=80)
           		{
                		continue;
           		}

	

			int flagTwoTofMatchedPrimaryTracks = 0;
			for (Int_t i = 0; i<tempUPCevent->getNumberOfTracks(); i++)
			{
				if(tempUPCevent->getTrack(i)->getFlag(StUPCTrack::kTof))
				{
					flagTwoTofMatchedPrimaryTracks+=1;
				}
			}       		
          

			StUPCTrack *firstTofMatchedTrack = tempUPCevent->getTrack(0);
			StUPCTrack *secondTofMatchedTrack = tempUPCevent->getTrack(1);

		
			//SC3.1:  Exactly two TOF-matched primary tracks
			bool SC3_1 = 0;
			if (flagTwoTofMatchedPrimaryTracks == 2)
			{
				SC3_1 = 1;
			}
			if (SC3_1 == 0) {continue;} 
			
			//SC3.2:  Tracks are of opposied signs				
			bool SC3_2 = 0;
			if (firstTofMatchedTrack->getEta()*secondTofMatchedTrack->getEta() < 0.0) 
			{
				SC3_2 = 1;
			}
			if (SC3_2 == 0) {continue;}
			
			
			//SC3.3: Both tracks are contained within the kinematic range			
			bool SC3_3 = 0;	
			if (abs(firstTofMatchedTrack->getEta()) < 0.7 and abs(secondTofMatchedTrack->getEta()) < 0.7 and  firstTofMatchedTrack->getPt() > 0.2 and  secondTofMatchedTrack->getPt() > 0.2)   
			
			if (firstTofMatchedTrack->getEta()*secondTofMatchedTrack->getEta() < 0.0) 
			{
				SC3_3 = 1;
			}
			if (SC3_3 == 0) {continue;}
			
			//SC3.4: Associated global tracks satisfy quality criteria: N_{fit} >= 25 and N_{dE/dx} >= 15 and d0 > 1.5
			bool SC3_4 = 0;
			if (firstTofMatchedTrack->getNhitsFit() >= 25 and secondTofMatchedTrack->getNhitsFit() >= 25 and firstTofMatchedTrack->getNhitsDEdx() >= 15 
			and secondTofMatchedTrack->getNhitsDEdx() >= 15 and firstTofMatchedTrack->getDcaXY() < 1.5 and  secondTofMatchedTrack->getDcaXY() < 1.5)
			{
				SC3_4 = 1;
			}
			if (SC3_4 == 0) {continue;}		
			
			//SC3.5 Associated global tracks match well to the primary vertex DCA(R) < 1.5 |DCA(z)| < 1.0			
			bool SC3_5 = 0;
			 if (sqrt(pow(firstTofMatchedTrack->getDcaXY(),2) + pow(firstTofMatchedTrack->getDcaZ(),2)) < 1.5 and  sqrt(pow(secondTofMatchedTrack->getDcaXY(),2) + pow(secondTofMatchedTrack->getDcaZ(),2)) < 1.5  and abs(firstTofMatchedTrack->getDcaZ()) < 1.0 and abs(secondTofMatchedTrack->getDcaZ()) < 1.0)
			{
				SC3_5 = 1;
			}
			if (SC3_5 == 0) {continue;}	
		
			//SC3.6 Associated global tracks are close at the beamline: |dz0| < 2 I usec difference in DCA
			bool SC3_6 = 0;
			if (abs(firstTofMatchedTrack->getDcaZ() - secondTofMatchedTrack->getDcaZ()) < 2.0 )
			{
				SC3_6 = 1;
			}
			
			if (SC3_6 == 0) {continue;}
	
			bool SC4_1 = 1;
			for(unsigned int k = 0; k < tempRPevent->getNumberOfTrackPoints(); ++k)
			{
				if(tempRPevent->getTrackPoint(k)->planesUsed() <3 )
				{
					SC4_1 = 0;		
					break;
				}
				
			}	
			if ( SC4_1 == 0) {continue;}	
			
			bool SC4_2 = 1;
            		for(unsigned int k = 0; k < tempRPevent->getNumberOfTracks(); ++k)
            		{
            			StUPCRpsTrack *trk = tempRPevent->getTrack(k);
				trk->setEvent(tempRPevent);
               			
                		double thetaRPx = tempRPevent->getTrack(k)->thetaRp(0);
                 		double thetaRPy = tempRPevent->getTrack(k)->thetaRp(1);  
                		double thetax = tempRPevent->getTrack(k)->theta(0);
                 		double thetay = tempRPevent->getTrack(k)->theta(1);
        			
        			if ((thetaRPx-thetax) <= -2/10e3 or (thetaRPx-thetax) >= 4/10e3 or (thetaRPy-thetay) <= -2/10e3 or (thetaRPy-thetay) >= 2/10e3)
        			{
        				SC4_2 = 0;	 				
      				}
            		}		
			if (SC4_2 == 0) {continue;}
			
				
					
			bool SC4_3 = 0;

			int numberOfTracksEast = 0;
			int numberOfTracksWest = 0;
			for(unsigned int k = 0; k < tempRPevent->getNumberOfTracks(); ++k)
            		{
				int iBranch = tempRPevent->getTrack(k)->branch();
				if (iBranch < 2) {numberOfTracksEast+=1;}
				else numberOfTracksWest+=1;
            		}

			if (numberOfTracksWest==1 and numberOfTracksEast == 1)
			{
				SC4_3 = 1;
			}
			if (SC4_3 == 0) {continue;}
		
			
				
			bool SC4_4 = 0;
            		for(unsigned int k = 0; k < tempRPevent->getNumberOfTracks(); ++k)
            		{		
			
               			StUPCRpsTrack *trk = tempRPevent->getTrack(k);
                		trk->setEvent(tempRPevent);
                		pVector = trk->pVec();
			
				if ( abs(pVector.Y()) >= 0.4 or abs(pVector.Y()) <= 0.2 or pVector.X() <= -0.2 or (  pow(pVector.X()+0.3,2) + pow(pVector.Y(), 2)  >=  pow(0.5,2)))
				{
					SC4_4 = 1;
				}  	
			}
			if (SC4_4 == 0) {continue;} 
			
			for (Int_t i = 0; i < tempUPCevent->getNumberOfTracks(); i++) 			
 			{
		    		if (abs(tempUPCevent->getTrack(i)->getNSigmasTPCPion())<3)
		    		{
		    		
		    			trackInstanceLocalUPC = tempUPCevent->getTrack(i);
					trackInstanceLocalUPC->getLorentzVector(trackVector, particleMass[Pion]);
                   			if(tempUPCevent->getTrack(i)->getCharge()==1)
                   			{
                        			posPion.push_back(trackVector);
                    			}
           				else if (tempUPCevent->getTrack(i)->getCharge()==-1)
           				{
           					negPion.push_back(trackVector);
                   			} 	 
				}							
			}	
			TLorentzVector tempSum = {0, 0, 0, 0};
			for (long unsigned int posi = 0; posi < posPion.size(); posi++)
			{
				for (long unsigned int negi = 0; negi < negPion.size(); negi++)
				{
					tempSum = posPion[posi] + negPion[negi];
					InvMassHistLocal->Fill(tempSum.M());
					double px = tempSum.Px();
					double py = tempSum.Py();
					double pz = tempSum.Pz();
					double p = sqrt(pow(px,2) + pow(py,2) + pow(pz,2));
					double sintheta = py/p;
					double l = R/sintheta;
					
					histL->Fill(l);
					
					if (R < 0.1)
					{
						InvMassHistR0_01Local->Fill(tempSum.M());
						histLLocal0_01->Fill(l);
					}
					
					else if (R >= 0.1 and R < 0.25)
					{
						InvMassHistR01_025Local->Fill(tempSum.M());
						histLLocal01_025->Fill(l);
					}
					
					else if (R >= 0.25)			
					{
						InvMassHistR025Local->Fill(tempSum.M());
						histLLocal025->Fill(l);
					}	
				}
			}
   
   
   
   
		filtered_entries++;
		

		}while (myReader.Next());
            posPion.clear();
            negPion.clear();
            posPionSuspicious.clear();
            negPionSuspicious.clear();
        //waiting for file opening to check if there were any filtered entries
        if(filtered_entries==0)
        {
            cout<<"Finished operation on file "<<myFile->GetTitle()<<endl;
            cout<<"Analyzed "<<tempTree->GetEntries()<<" entries"<<endl;
            cout<<"There were 0 filtered entries"<<endl;

            delete tempFile;
            return 0;
        }

        cout<<"Finished operation on file "<<myFile->GetTitle()<<endl;
        cout<<"Analyzed "<<tempTree->GetEntries()<<" entries"<<endl;
        cout<<"Filtered "<<filtered_entries<<" entries"<<endl;

        delete tempFile;

        return filtered_entries;
	};

	auto redFunction = [](const std::vector<int> &mapV)
	{
		return std::accumulate(mapV.begin(), mapV.end(), 0);
	};

	int filtered_entries = 0;
	vector<TFile*> listOfFiles;

	TObjArray* tempList = upcChain->GetListOfFiles();
	for (int i = 0; i < tempList->GetEntries(); i++)
	{
		listOfFiles.push_back((TFile*)(tempList->At(i)));
	}

	ROOT::TThreadExecutor TreeProcessor(nthreads);
	filtered_entries = TreeProcessor.MapReduce(myFunction, listOfFiles, redFunction);

	std::shared_ptr<TH1D> InvMassHistFinal;
	InvMassHistFinal = InvMassHist.Merge();
	
	std::shared_ptr<TH1D> InvMassHistR0_01Final;	
	InvMassHistR0_01Final = InvMassHistR0_01.Merge();

	std::shared_ptr<TH1D> InvMassHistR01_025Final;	
	InvMassHistR01_025Final = InvMassHistR01_025.Merge();
		
	std::shared_ptr<TH1D> InvMassHistR025Final;	
	InvMassHistR025Final = InvMassHistR025.Merge();

	std::shared_ptr<TProfile> DifferenceVertexPositionBeamLineXFinal;
	DifferenceVertexPositionBeamLineXFinal = DifferenceVertexPositionBeamLineX.Merge();
	std::shared_ptr<TProfile> DifferenceVertexPositionBeamLineYFinal;
	DifferenceVertexPositionBeamLineYFinal = DifferenceVertexPositionBeamLineY.Merge();

	std::shared_ptr<TH1D> PerpendicularDistanceVertexPositionBeamLineFinal;
	PerpendicularDistanceVertexPositionBeamLineFinal = PerpendicularDistanceVertexPositionBeamLine.Merge();

	std::shared_ptr<TH1D> SignificanceFinal;
	SignificanceFinal = Significance.Merge();
	
	std::shared_ptr<TH1D> histLFinal;
	histLFinal = histL.Merge();

	std::shared_ptr<TH1D> histLFinal0_01;
	histLFinal0_01 = histL0_01.Merge();
	
	std::shared_ptr<TH1D> histLFinal01_025;
	histLFinal01_025 = histL01_025.Merge();

	std::shared_ptr<TH1D> histLFinal025;
	histLFinal025 = histL025.Merge();
	string outfileName = outputFolder + "new.root";
	cout<<"Created output file "<<outfileName<<endl;
	
	TFile* outputFileHist = TFile::Open(outfileName.c_str(), "recreate");

	outputFileHist->cd();
	   
	InvMassHistFinal->Write();
	DifferenceVertexPositionBeamLineXFinal->Write();
	DifferenceVertexPositionBeamLineYFinal->Write();   
	PerpendicularDistanceVertexPositionBeamLineFinal->Write(); 
	InvMassHistR0_01Final->Write();
	InvMassHistR01_025Final->Write();
	InvMassHistR025Final->Write();
	SignificanceFinal->Write();
	histLFinal->Write();
	histLFinal0_01->Write();	
	histLFinal01_025->Write();
	histLFinal025->Write();
	outputFileHist->Close();

	return 0;
}





bool ConnectInput(int argc, char** argv, TChain* fileChain) 
{
	int fileId = -1;
	string line;
	int lineId=0;

	const string& input = argv[1];
	cout<<"Using list "<<input<<endl;
	if(input.find(".list") != string::npos )
	{
		cout << "Input from chain" << endl;
		ifstream instr(input.c_str());
        	if (!instr.is_open())
        	{
        		cout<< "Couldn't open: "<<input.c_str()<<endl;
            		return false;
        	}	
        	TFile* infile;
        	while(getline(instr, line)) 
        	{
            		if(fileId==lineId || fileId== -1)
            		{
		        	fileChain->AddFile(line.c_str());
		        	infile = TFile::Open(line.c_str(), "read");
		        	if(!infile)
		        	{
		            		cout<< "Couldn't open: "<<line.c_str()<<endl;
		            		return false;
		        	}
				infile->Close();
        		}
            		lineId++;
        	}	
		instr.close();
	}
	return true;
}

bool CheckTriggers(StUPCEvent* localupcEvt)
{
	bool CPTtrigger = false;
	for(int var = 0; var < nTriggers; ++var)
	{
		if(localupcEvt->isTrigger(triggerID[var]))
        	{
            	for (int i = 0; i < *(&CEPtriggers + 1) - CEPtriggers; ++i)
                	if(triggerID[var] == CEPtriggers[i])
                CPTtrigger=true;
        	}
	}
	return CPTtrigger;
}

long long GetFileSize(string filename)
{
	struct stat64 stat_buf;
	int rc = stat64(filename.c_str(), &stat_buf);
	return rc==0 ? stat_buf.st_size : -1;
}
