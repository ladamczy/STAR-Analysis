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
#include <TParticle.h>
#include <TParticlePDG.h>
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
#include "StUPCV0.h"


using namespace std;

int main(int argc, char** argv)  
{
//    ifstream inputFilePathList(argv[1]);
//    if (!inputFilePathList) 
//
//    {
//        cerr << "Failed to open input file." << std::endl;
//        return 1;
//    }

    TChain *chain = new TChain("mUPCTree"); 

    string inputFileName;
    string outputFileName;

    std::string run = argv[1];
    inputFileName="/data2/sd_star_2017/presel_"+run+".root";      
    outputFileName="Run/"+run+".root";       

    chain->AddFile(inputFileName.c_str());

    static StUPCEvent * upcEvt = 0x0;
    static StRPEvent  * rpEvt = 0x0;

    chain->SetBranchAddress("mUPCEvent", &upcEvt);
    chain->SetBranchAddress("correctedRpEvent",  &rpEvt);

    TH1D* HistNumOfPrimaryTracksToF = new TH1D("HistNumOfPrimaryTracksToF", "; Num of tracks; # events", 50 ,0, 50);
    TH1D* HistXiProtonWest = new TH1D("HistXiProtonWest", "; xi of the proton on west; # events", 110 ,-0.1,1);
    TH1D* HistXiProtonEast = new TH1D("HistXiProtonEast", "; xi of the proton on east; # events", 110 ,-0.1,1);

    for (Long64_t i = 0; i < chain->GetEntries(); ++i) 
    {
       chain->GetEntry(i);
       if( upcEvt->getNumberOfVertices() != 1) continue;
       if (rpEvt->getNumberOfTracks() != 1 ) continue;

//       cout << upcEvt->getEventNumber() << " " 
//            << upcEvt->getNPrimVertices() <<  endl;

       int NumOfPrimaryTracksToF = 0;

       for (int j = 0; j < upcEvt->getNumberOfTracks(); j++) {
              if(   upcEvt->getTrack(j)->getNhits()>15  
                 && upcEvt->getTrack(j)->getFlag(StUPCTrack::kTof) 
                 && upcEvt->getTrack(j)->getFlag(StUPCTrack::kPrimary)
                 && upcEvt->getTrack(j)->getPt()>0.15 
                 && abs(upcEvt->getTrack(j)->getEta())<1.0 ) {

              NumOfPrimaryTracksToF++;
              }
       }
       HistNumOfPrimaryTracksToF->Fill(NumOfPrimaryTracksToF);

       StUPCRpsTrack *proton = rpEvt->getTrack(0);
       if (proton->branch() < 2)  HistXiProtonEast->Fill(proton->xi(254.867));
       if (proton->branch() > 1)  HistXiProtonWest->Fill(proton->xi(254.867));
     }

    
    TFile *outfile = TFile::Open(outputFileName.c_str(), "recreate"); 

    HistNumOfPrimaryTracksToF->Write();
    HistXiProtonEast->Write();
    HistXiProtonWest->Write();

    outfile->Close();
    return 0;
}
 
