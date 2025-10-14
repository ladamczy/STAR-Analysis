// Run by: ./Ana file.list output.root
// e.g. ./Ana /star/u/truhlar/star-upcDst/build/run17.list ./output.root


// Table of RP indecies and names
// RP_ID   0,    1,    2,   3,   4,   5,   6, 7
// RP_name E1U, E1D, E2U, E2D, W1U, W1D, W2U, W2D

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
#include <filesystem>

// ROOT headers
#include "TROOT.h"
#include "TSystem.h"
#include "TThread.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
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

//repeating things
#include "UsefulThings.h"

using namespace std;

// enums are very usefull 
enum { kAll = 1, kCPT,  kRP, kOneVertex, kTPCTOF, 
    kTotQ, kMax};
enum SIDE {E = 0, East = 0, W = 1, West = 1, nSides};
enum PARTICLES {Pion = 0, Kaon = 1, Proton = 2, nParticles};
enum BRANCH_ID { EU, ED, WU, WD, nBranches };
enum RP_ID {E1U, E1D, E2U, E2D, W1U, W1D, W2U, W2D, nRomanPots};

const double particleMass[nParticles] = { 0.13957, 0.493677, 0.93827 }; // pion, kaon, proton in GeV /c^2 

//_____________________________________________________________________________
int main(int argc, char** argv) 
{
    int nthreads = 2;
    if(argc==4){
        nthreads = atoi(argv[3]);
    }
    cout<<"Program is running on "<<nthreads<<" threads"<<endl;
    ROOT::EnableThreadSafety();
    ROOT::EnableImplicitMT(nthreads); //turn on multicore processing
    ROOT::EnableThreadSafety();
    
    TChain* upcChain = new TChain("mUPCTree");    //chain with files to iterate through

    if(!ConnectInput(argc, argv, upcChain)){
        cout << "Wrong input parameters..." << endl; 
        return 1;
    }

    const string& outputFolder = argv[2];

    //HISTOGRAMS
    ROOT::TThreadedObject<TH1D> invmasshistbefore = ROOT::TThreadedObject<TH1D>("invmasshistbefore", "Inv. mass of pions before cut;m_{inv} [GeV];events", 100, 0.42, 0.56);
    ROOT::TThreadedObject<TH1D> invmasshistafter = ROOT::TThreadedObject<TH1D>("invmasshistafter", "Inv. mass of pions after cut;m_{inv} [GeV];events", 100, 0.42, 0.56);

    ROOT::TThreadedObject<TH1D> phianglehistbefore = ROOT::TThreadedObject<TH1D>("phianglehistbefore", "Azimuthal angle of K_{0}^{S} before cut;m_{inv} [GeV];events", 100, -3.14159, 3.14159);
    ROOT::TThreadedObject<TH1D> phianglehistafter = ROOT::TThreadedObject<TH1D>("phianglehistafter", "Azimuthal angle of K_{0}^{S} after cut;m_{inv} [GeV];events", 100, -3.14159, 3.14159);

    ROOT::TThreadedObject<TH1D> pipairanglehist = ROOT::TThreadedObject<TH1D>("pipairanglehist", "Angle between pions in detector FoR;#phi angle [rad];events", 30, 0, 3.1416);

    int nEvents = 7;
    ROOT::TThreadedObject<TH1D> cutflow("cutflow", "Cutflow;;Tracks", nEvents, 0, nEvents);
    cutflow->GetXaxis()->SetBinLabel(1, "Before tests");
    cutflow->GetXaxis()->SetBinLabel(2, "Good PV position");
    cutflow->GetXaxis()->SetBinLabel(3, "N hits");
    cutflow->GetXaxis()->SetBinLabel(4, "p_{T} & #eta");
    cutflow->GetXaxis()->SetBinLabel(5, "Not p & not K");
    cutflow->GetXaxis()->SetBinLabel(6, "Probably #pi");
    cutflow->GetXaxis()->SetBinLabel(7, "Check DCA");
    
    // Define the function that will process a subrange of the tree.
    // The function must receive only one parameter, a TTreeReader,
    // and it must be thread safe. To enforce the latter requirement,
    // TThreadedObject histograms will be used.
    //but maybe later
    auto myFunction = [&](TFile* myFile) {
        //test if tree is not empty
        TFile* tempFile = new TFile(myFile->GetTitle());
        TTree* tempTree = (TTree*)tempFile->Get("mUPCTree");
        if(tempTree->GetEntries()==0){
            delete tempFile;
            return 0;
        }
        //creating a reader and all stuff
        TTreeReader myReader(tempTree);
        TTreeReaderValue<StUPCEvent> StUPCEventInstance(myReader, "mUPCEvent");
        //unavailable
        // TTreeReaderValue<StRPEvent> StRPEventInstance(myReader, "mRPEvent");
        //variable initialization
        myReader.Next();
        //histograms
        std::shared_ptr<TH1D> invmasshistbeforeLocal;
        std::shared_ptr<TH1D> invmasshistafterLocal;
        std::shared_ptr<TH1D> phianglehistbeforeLocal;
        std::shared_ptr<TH1D> phianglehistafterLocal;
        std::shared_ptr<TH1D> pipairanglehistLocal;
        std::shared_ptr<TH1D> cutflowLocal;

        invmasshistbeforeLocal = invmasshistbefore.Get();
        invmasshistafterLocal = invmasshistafter.Get();
        phianglehistbeforeLocal = phianglehistbefore.Get();
        phianglehistafterLocal = phianglehistafter.Get();
        pipairanglehistLocal = pipairanglehist.Get();
        cutflowLocal = cutflow.Get();

        int filtered_entries = 0;
        //for changing branch address
        StUPCEvent* tempUPCpointer = StUPCEventInstance.Get();
        //unavailable
        // StRPEvent* tempRPpointer = StRPEventInstance.Get();

        TLorentzVector trackVector;
        //unavailable
        // Double_t verDeltaZ;
        TVector3 pVector;
        vector<TLorentzVector> posPion;
        vector<TLorentzVector> negPion;
        do{
            tempUPCpointer = StUPCEventInstance.Get();
            //unavailable
            // tempRPpointer = StRPEventInstance.Get();

            //TESTS & HISTOGRAMS
            //triggers
            //unavailable
            // if(!CheckTriggers(StUPCEventInstance.Get())){
            //     continue;
            // }
            
            //one proton each side
            //unavailable
            // int numberOfTracksPerSide[nSides] = {0, 0};
            // for(unsigned int k = 0; k < tempRPpointer->getNumberOfTracks(); ++k)
            // {
            //     // Get pointer to k-th track in Roman Pot data collection
            //     StUPCRpsTrack *trk = tempRPpointer->getTrack(k);
            //     trk->setEvent(tempRPpointer);
            //     // Get ID of a branch in which this k-th track was reconstructed
            //     int j = trk->branch();
            //     int side = j<2 ? E : W;
            //     numberOfTracksPerSide[side]++;
            // }
            // if(numberOfTracksPerSide[0]!=1 || numberOfTracksPerSide[1]!=1){
            //     continue;
            // }
            
            //four planes on each TrackPoint
            //unavailable
            // bool areAllPlanesPresent = true;
            // for(unsigned int k = 0; k < tempRPpointer->getNumberOfTracks(); ++k)
            // {
            //     // Get pointer to k-th track in Roman Pot data collection
            //     StUPCRpsTrack *trk = tempRPpointer->getTrack(k);
            //     trk->setEvent(tempRPpointer);
            //     // test if all 8 planes are present
            //     if(trk->planesUsed()!=8){
            //         areAllPlanesPresent = false;
            //         break;
            //     }
            // }
            // if(!areAllPlanesPresent){
            //     continue;
            // }

            //exactly one primary vertex to which TOF tracks match to
            //unavailable
            //replaced below
            // vector<Int_t> primaryVertices;
            // Int_t VertexId;
            // for (Int_t i = 0; i < tempUPCpointer->getNumberOfTracks(); i++){
            //     VertexId = tempUPCpointer->getTrack(i)->getVertexId();
            //     if(tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kTof) && tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kPrimary) && !(find(primaryVertices.begin(), primaryVertices.end(), VertexId)!=primaryVertices.end())){
            //         primaryVertices.push_back(VertexId);
            //     }
            // }
            // if(primaryVertices.size()!=1){
            //     continue;
            // }
            // VertexId = primaryVertices[0];

            //exactly one primary vertex
            //replacement
            vector<Int_t> primaryVertices;
            Int_t VertexId;
            //cutflow
            cutflowLocal->Fill(0.0, tempUPCpointer->getNumberOfTracks());
            for (Int_t i = 0; i < tempUPCpointer->getNumberOfTracks(); i++){
                VertexId = tempUPCpointer->getTrack(i)->getVertexId();
                if(tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kPrimary) && !(find(primaryVertices.begin(), primaryVertices.end(), VertexId)!=primaryVertices.end())){
                    primaryVertices.push_back(VertexId);
                }
            }
            if(primaryVertices.size()!=1){
                continue;
            }
            VertexId = primaryVertices[0];

            //primary vertex is placed within 80cm of the centre
            if(abs(tempUPCpointer->getVertexId(VertexId)->getPosZ())>=80){
                continue;
            }

            //vertexes from TPC and RP are not too far away (36cm)
            //unavailable
            // bool areBothRPTracksWithTime = true;
            // Double_t time = 0;
            // for(unsigned int k = 0; k < tempRPpointer->getNumberOfTracks(); ++k)
            // {
            //     // Get pointer to k-th track in Roman Pot data collection
            //     StUPCRpsTrack *trk = tempRPpointer->getTrack(k);
            //     trk->setEvent(tempRPpointer);
            //     // Get ID of a branch in which this k-th track was reconstructed
            //     int j = trk->branch();
            //     //deletes those tracks which do not have confirmed time
            //     if(trk->time()<0){
            //         areBothRPTracksWithTime = false;
            //         break;
            //     }
            //     int side = j<2 ? E : W;
            //     //to not make the loop go twice, also west side is the "+" one
            //     if(side==E){
            //         time -= trk->time();
            //     }
            //     if(side==W){
            //         time += trk->time();
            //     }
            // }
            // if(!areBothRPTracksWithTime){
            //     continue;
            // }
            // //*100 at the end it to convert from m to cm
            // verDeltaZ = tempUPCpointer->getVertexId(VertexId)->getPosZ()-299792458.0*time/2*100;
            // if(abs(verDeltaZ)>36){
            //     continue;
            // }

            //fiducial cuts
            //unavailable
            // bool passedFiducialCuts = true;
            // //px barrier
            // Double_t pxbarrier = -0.27;
            // //py low nand high barrier
            // Double_t pylowbarrier = 0.4;
            // Double_t pyhighbarrier = 0.8;
            // //p circle parameters
            // Double_t pxcenter = -0.6;
            // Double_t pycenter = 0;
            // Double_t pradius = 1.1;
            // for(unsigned int k = 0; k < tempRPpointer->getNumberOfTracks(); ++k){
            //     // Get pointer to k-th track in Roman Pot data collection
            //     StUPCRpsTrack *trk = tempRPpointer->getTrack(k);
            //     trk->setEvent(tempRPpointer);
            //     pVector = trk->pVec();
            //     bool f1 = pow(pVector.X()-pxcenter, 2) + pow(pVector.Y()-pycenter, 2) < pradius*pradius;
            //     bool f2 = pylowbarrier < abs(pVector.Y()) && abs(pVector.Y()) < pyhighbarrier;
            //     bool f3 = pxbarrier < pVector.X();
            //     if(!(f1 && f2 && f3)){
            //         passedFiducialCuts = false;
            //         break;
            //     }
            // }            
            // if (!passedFiducialCuts){
            //     continue;
            // }

            //elastic cuts
            //unavailable
            // TLorentzVector RPsum = {0,0,0,0};
            // //checking total momentum
            // for(unsigned int k = 0; k < tempRPpointer->getNumberOfTracks(); ++k){
            //     // Get pointer to k-th track in Roman Pot data collection
            //     StUPCRpsTrack *trk = tempRPpointer->getTrack(k);
            //     trk->setEvent(tempRPpointer);
            //     trackVector.SetVectM(trk->pVec(), particleMass[Proton]);
            //     RPsum += trackVector;
            // }
            // //xposition of the gauss
            // Double_t xGauss = -0.035;
            // //yposition of the gauss
            // Double_t yGauss = 0.01;
            // //radius of the gauss
            // //was 0.07, 0.1 FOR TEST
            // Double_t rGauss = 0.07;
            // if (pow(RPsum.X()-xGauss, 2) + pow(RPsum.Y()-yGauss, 2) < rGauss*rGauss){
            //     continue;
            // }

            //quality test of central tracks
            //1. Number of hits
            //2. DCA in R and z (the second one also doubles as spread limiter)
            //3. kinematic range for optimal measurement
            for (Int_t i = 0; i < tempUPCpointer->getNumberOfTracks(); i++){
                cutflowLocal->Fill(1);
                //quality
                if(tempUPCpointer->getTrack(i)->getNhitsDEdx()<15 || tempUPCpointer->getTrack(i)->getNhitsFit()<25){
                    continue;
                }
                cutflowLocal->Fill(2);
                if(tempUPCpointer->getTrack(i)->getPt()<=0.2 || abs(tempUPCpointer->getTrack(i)->getEta())>=0.7){
                    continue;
                }
                cutflowLocal->Fill(3);
                if(abs(tempUPCpointer->getTrack(i)->getNSigmasTPCProton())<3 || abs(tempUPCpointer->getTrack(i)->getNSigmasTPCKaon())<3){
                    continue;
                }
                cutflowLocal->Fill(4);
                if (abs(tempUPCpointer->getTrack(i)->getNSigmasTPCPion())<3){
                    tempUPCpointer->getTrack(i)->getLorentzVector(trackVector, particleMass[Pion]);
                    cutflowLocal->Fill(5);
                }else{
                    continue;
                }
                //DCA only active for pion candidates, this way it allows more not-pions and thus changes the amounts of additional particles
                if((pow(tempUPCpointer->getTrack(i)->getDcaXY(), 2) + pow(tempUPCpointer->getTrack(i)->getDcaZ(), 2))<1 && tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kPrimary)){
                    continue;
                }
                cutflowLocal->Fill(6);
                //basic test of ToF
                //unavailable
                // if(!tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kTof) || tempUPCpointer->getTrack(i)->getTofPathLength()<=0 || tempUPCpointer->getTrack(i)->getTofTime()<=0){
                //     continue;
                // }
                //test for charge
                if(tempUPCpointer->getTrack(i)->getCharge()==1){
                    posPion.push_back(trackVector);
                }else if (tempUPCpointer->getTrack(i)->getCharge()==-1){
                    negPion.push_back(trackVector);
                }
            }

            if(posPion.size()==0 || negPion.size()==0){
                continue;
            }

            //for mass
            TLorentzVector tempSum = {0, 0, 0, 0};
            //invariant mass, pi+ pi-, legit ones
            for (long unsigned int posi = 0; posi < posPion.size(); posi++){
                for (long unsigned int negi = 0; negi < negPion.size(); negi++)
                {
                    tempSum = posPion[posi] + negPion[negi];
                    pipairanglehistLocal->Fill(posPion[posi].Angle(negPion[negi].Vect()));
                    invmasshistbeforeLocal->Fill(tempSum.M());
                    phianglehistbeforeLocal->Fill(tempSum.Phi());
                    if(posPion[posi].Angle(negPion[negi].Vect())>1.8){
                        continue;
                    }
                    invmasshistafterLocal->Fill(tempSum.M());
                    phianglehistafterLocal->Fill(tempSum.Phi());
                }
            }

            filtered_entries++;

            posPion.clear();
            negPion.clear();
        }while (myReader.Next());

        //waiting for file opening to check if there were any filtered entries
        if(filtered_entries==0){
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
    //creating the list of TFile*
    TObjArray* tempList = upcChain->GetListOfFiles();
    for (int i = 0; i < tempList->GetEntries(); i++)
    {
        listOfFiles.push_back((TFile*)(tempList->At(i)));
    }

    // Create a TreeProcessor: specify the file and the tree in it
    ROOT::TThreadExecutor TreeProcessor(nthreads);
    // Launch the parallel processing of the tree
    filtered_entries = TreeProcessor.MapReduce(myFunction, listOfFiles, redFunction);
    // Use the TThreadedObject::Merge method to merge the thread private tree
    // into the final result
    std::shared_ptr<TH1D> invmasshistbeforeFinal;
    std::shared_ptr<TH1D> invmasshistafterFinal;
    std::shared_ptr<TH1D> phianglehistbeforeFinal;
    std::shared_ptr<TH1D> phianglehistafterFinal;
    std::shared_ptr<TH1D> pipairanglehistFinal;
    std::shared_ptr<TH1D> cutflowFinal;

    invmasshistbeforeFinal = invmasshistbefore.Merge();
    invmasshistafterFinal = invmasshistafter.Merge();
    phianglehistbeforeFinal = phianglehistbefore.Merge();
    phianglehistafterFinal = phianglehistafter.Merge();
    pipairanglehistFinal = pipairanglehist.Merge();
    cutflowFinal = cutflow.Merge();

    invmasshistbeforeFinal->SetMinimum(0);
    invmasshistafterFinal->SetMinimum(0);
    phianglehistbeforeFinal->SetMinimum(0);
    phianglehistafterFinal->SetMinimum(0);
    pipairanglehistFinal->SetMinimum(0);

    //setting up a tree & output file
    string path = string(argv[0]);
    string outfileName = outputFolder + "AnaOutput_" + path.substr(path.find_last_of("/\\") + 1) + ".root";
    cout<<"Created output file "<<outfileName<<endl;
    TFile* outputFileHist = TFile::Open(outfileName.c_str(), "recreate");

    outputFileHist->cd();
    invmasshistbeforeFinal->Write();
    invmasshistafterFinal->Write();
    phianglehistbeforeFinal->Write();
    phianglehistafterFinal->Write();
    pipairanglehistFinal->Write();
    cutflowFinal->Write();
    outputFileHist->Close();

    cout<<"Finished processing "<<endl;
    cout<<"Analyzed total "<<upcChain->GetEntries()<<" entries"<<endl;
    cout<<"Filtered total "<<filtered_entries<<" entries"<<endl;
    cout<<"Ending Analysis... GOOD BYE!"<<endl;
    return 0;
}//main