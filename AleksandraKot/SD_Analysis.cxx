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
    inputFileName="/data1/sd_star_2017/presel_"+run+".root";      
    outputFileName="Run/"+run+".root";       

    chain->AddFile(inputFileName.c_str());

    static StUPCEvent * upcEvt = 0x0;
    static StRPEvent  * rpEvt = 0x0;

    chain->SetBranchAddress("mUPCEvent", &upcEvt);
    chain->SetBranchAddress("correctedRpEvent",  &rpEvt);


    TFile *outfile = TFile::Open(outputFileName.c_str(), "recreate"); 
    TDirectory* dirNSigma = outfile->mkdir("dEdx_Nsigma");

    outfile->cd();
    TH1D* HistNumOfPrimaryTracksToF = new TH1D("HistNumOfPrimaryTracksToF", "; Num of tracks; # events", 50 ,0, 50);
    TH1D* HistNumOfPToFWest = new TH1D("HistNumOfPToFWest", "; Num of tracks when proton on west; # events", 50 ,0, 50);
    TH1D* HistNumOfPToFEast = new TH1D("HistNumOfPToFEast", "; Num of tracks when proton on east; # events", 50 ,0, 50);
    TH1D* HistXiProtonWest = new TH1D("HistXiProtonWest", "xi of the proton on west;xi; # events", 100, -0.1, 1);
    TH1D* HistXiProtonEast = new TH1D("HistXiProtonEast", "xi of the proton on east;xi; # events", 100, -0.1, 1);
    TH1D* HistLogXiProtonWest = new TH1D("HistLogXiProtonWest", "logxi of the proton on west; logxi; # events", 100, -6, 0);
    TH1D* HistLogXiProtonEast = new TH1D("HistLogXiProtonEast", "logxi of the proton on east; logxi; # events", 100, -6, 0);
    TH1D* HistPtProtonWest = new TH1D("HistPtProtonWest", "pT of the proton on west; pT [GeV]; # events", 100, 0, 4);
    TH1D* HistPtProtonEast = new TH1D("HistPtProtonEast", "pT of the proton on west; pT [GeV]; pT [GeV]; # events", 100, 0, 4);
    TH1D* HistEtaProtonWest = new TH1D("HistEtaProtonWest", "Eta of the proton on west; #eta; # events", 100, -10, 10);
    TH1D* HistEtaProtonEast = new TH1D("HistEtaProtonEast", "Eta of the proton on east; #eta; # events", 100, -10, 10);
    TH1D* HistPtTracksWest = new TH1D("HistPtTracksWest", "pT of the reconstructed tracks when proton on west; pT [GeV]; # events", 100, 0, 4);
    TH1D* HistPtTracksEast = new TH1D("HistPtTracksEast", "pT of the reconstructed tracks when proton on east; pT [GeV]; # events", 100, 0, 4);
    TH1D* HistEtaTracksWest = new TH1D("HistEtaTracksWest", "Eta of the reconstructed tracks when proton on west; #eta; # events", 100, -3, 3);
    TH1D* HistEtaTracksEast = new TH1D("HistEtaTracksEast", "Eta of the reconstructed tracks when proton on east; #eta; # events", 100, -3, 3);
    TH1D* HistPtTracksWestCut = new TH1D("HistPtTracksWestCut", "pT of the reconstructed tracks when proton on west; pT [GeV]; # events", 100, 0, 4);
    TH1D* HistPtTracksEastCut = new TH1D("HistPtTracksEastCut", "pT of the reconstructed tracks when proton on east; pT [GeV]; # events", 100, 0, 4);
    TH1D* HistEtaTracksWestCut = new TH1D("HistEtaTracksWestCut", "Eta of the reconstructed tracks when proton on west; #eta; # events", 100, -3, 3);
    TH1D* HistEtaTracksEastCut = new TH1D("HistEtaTracksEastCut", "Eta of the reconstructed tracks when proton on east; #eta; # events", 100, -3, 3);
    
    TH2F* hNSigmaPiPlus  = new TH2F("hNSigmaPiPlus",  ";p_{T} [GeV/c];n#sigma_{#pi}^{TPC} (+)", 300, 0.0, 3.0, 200, -50.0, 50.0); //300 bins 0.0-3.0 for pT and 200 bins -50.0-50.0 for n_sigma
    TH2F* hNSigmaPiMinus = new TH2F("hNSigmaPiMinus", ";p_{T} [GeV/c];n#sigma_{#pi}^{TPC} (-)", 300, 0.0, 3.0, 200, -50.0, 50.0);
    TH2F* hNSigmaKPlus   = new TH2F("hNSigmaKPlus",   ";p_{T} [GeV/c];n#sigma_{K}^{TPC} (+)", 300, 0.0, 3.0, 200, -50.0, 50.0);
    TH2F* hNSigmaKMinus  = new TH2F("hNSigmaKMinus",  ";p_{T} [GeV/c];n#sigma_{K}^{TPC} (-)", 300, 0.0, 3.0, 200, -50.0, 50.0);
    TH2F* hNSigmaPPlus   = new TH2F("hNSigmaPPlus",   ";p_{T} [GeV/c];n#sigma_{p}^{TPC} (+)", 300, 0.0, 3.0, 200, -50.0, 50.0);
    TH2F* hNSigmaPMinus  = new TH2F("hNSigmaPMinus",  ";p_{T} [GeV/c];n#sigma_{p}^{TPC} (-)", 300, 0.0, 3.0, 200, -50.0, 50.0);

    TH2F* hdEdx = new TH2F("hdEdx",  ";q #times p [GeV/c]; dE/dx [keV/cm]", 600, -5.0, 5.0, 600, 0, 100);


       for (Long64_t i = 0; i < chain->GetEntries(); ++i) {
       chain->GetEntry(i);
       if( upcEvt->getNumberOfVertices() != 1) continue;
       if (rpEvt->getNumberOfTracks() != 1 ) continue;

//       cout << upcEvt->getEventNumber() << " " 
//            << upcEvt->getNPrimVertices() <<  endl;

       StUPCRpsTrack *proton = rpEvt->getTrack(0);

       bool isEast = (proton->branch() < 2);
       bool isWest = (proton->branch() > 1);

       int NumOfPrimaryTracksToF = 0;

       for (int j = 0; j < upcEvt->getNumberOfTracks(); j++) {
              if(   upcEvt->getTrack(j)->getNhits()>15  
                 && upcEvt->getTrack(j)->getFlag(StUPCTrack::kTof) 
                 && upcEvt->getTrack(j)->getFlag(StUPCTrack::kPrimary)) {
                     if (upcEvt->getTrack(j)->getPt()>0.2  && abs(upcEvt->getTrack(j)->getEta())<0.9){
                            NumOfPrimaryTracksToF++;
                           if (isEast) {
                            HistPtTracksEastCut->Fill(upcEvt->getTrack(j)->getPt());
                            HistEtaTracksEastCut->Fill(upcEvt->getTrack(j)->getEta());
                           }
                           if (isWest) {
                            HistPtTracksWestCut->Fill(upcEvt->getTrack(j)->getPt());
                            HistEtaTracksWestCut->Fill(upcEvt->getTrack(j)->getEta());
                           }
                     }
                     if (abs(upcEvt->getTrack(j)->getEta())<0.9) {
                           if (isEast) HistPtTracksEast->Fill(upcEvt->getTrack(j)->getPt());
                           if (isWest) HistPtTracksWest->Fill(upcEvt->getTrack(j)->getPt()); 
                           if(upcEvt->getTrack(j)->getCharge()>0) {
                                   hNSigmaPiPlus->Fill(upcEvt->getTrack(j)->getPt(), upcEvt->getTrack(j)->getNSigmasTPCPion());
                                   hNSigmaKPlus->Fill(upcEvt->getTrack(j)->getPt(), upcEvt->getTrack(j)->getNSigmasTPCKaon());
                                   hNSigmaPPlus->Fill(upcEvt->getTrack(j)->getPt(), upcEvt->getTrack(j)->getNSigmasTPCProton());
                           } else {
                                   hNSigmaPiMinus->Fill(upcEvt->getTrack(j)->getPt(), upcEvt->getTrack(j)->getNSigmasTPCPion());
                                   hNSigmaKMinus->Fill(upcEvt->getTrack(j)->getPt(), upcEvt->getTrack(j)->getNSigmasTPCKaon());
                                   hNSigmaPMinus->Fill(upcEvt->getTrack(j)->getPt(), upcEvt->getTrack(j)->getNSigmasTPCProton());
                           }
                           double pq = (upcEvt->getTrack(j)->getPt())*(upcEvt->getTrack(j)->getCharge());
                           hdEdx->Fill(pq, upcEvt->getTrack(j)->getDEdxSignal()*1e6); //GeV originally, keV now
                     } 
                     if (upcEvt->getTrack(j)->getPt()>0.2) {
                           if (isEast) HistEtaTracksEast->Fill(upcEvt->getTrack(j)->getEta());
                           if (isWest) HistEtaTracksWest->Fill(upcEvt->getTrack(j)->getEta()); 
                     }
              }
       }
       HistNumOfPrimaryTracksToF->Fill(NumOfPrimaryTracksToF);
       if (isEast) HistNumOfPToFEast->Fill(NumOfPrimaryTracksToF);
       if (isWest) HistNumOfPToFWest->Fill(NumOfPrimaryTracksToF);

       if (proton->branch() < 2)  {
       HistXiProtonEast->Fill(proton->xi(254.867));
       HistLogXiProtonEast->Fill(log10(proton->xi(254.867)));
       HistPtProtonEast->Fill(proton->pt());
       HistEtaProtonEast->Fill(proton->eta());
       }
       if (proton->branch() > 1)  {
       HistXiProtonWest->Fill(proton->xi(254.867));
       HistLogXiProtonWest->Fill(log10(proton->xi(254.867)));
       HistPtProtonWest->Fill(proton->pt());
       HistEtaProtonWest->Fill(proton->eta());
       }

     }
    
    HistNumOfPrimaryTracksToF->Write();
    HistNumOfPToFWest->Write();
    HistNumOfPToFEast->Write();

    HistXiProtonWest->Write();
    HistXiProtonEast->Write();

    HistLogXiProtonWest->Write();
    HistLogXiProtonEast->Write();

    HistPtProtonWest->Write();
    HistPtProtonEast->Write();

    HistEtaProtonWest->Write();
    HistEtaProtonEast->Write();

    HistPtTracksWest->Write();
    HistPtTracksEast->Write();

    HistEtaTracksWest->Write();
    HistEtaTracksEast->Write();

    HistPtTracksWestCut->Write();
    HistPtTracksEastCut->Write();

    HistEtaTracksWestCut->Write();
    HistEtaTracksEastCut->Write();

    dirNSigma->cd();

    hNSigmaPiPlus->Write();
    hNSigmaPiMinus->Write();
    hNSigmaKPlus->Write();
    hNSigmaKMinus->Write();
    hNSigmaPPlus->Write();
    hNSigmaPMinus->Write();

    hdEdx->Write();

    outfile->Close();
    return 0;
}
 
