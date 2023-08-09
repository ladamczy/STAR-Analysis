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

using namespace std;

// enums are very usefull 
enum { kAll = 1, kCPT,  kRP, kOneVertex, kTPCTOF, 
    kTotQ, kMax};
enum SIDE {E = 0, East = 0, W = 1, West = 1, nSides};
enum PARTICLES {Pion = 0, Kaon = 1, Proton = 2, nParticles};
enum BRANCH_ID { EU, ED, WU, WD, nBranches };
enum RP_ID {E1U, E1D, E2U, E2D, W1U, W1D, W2U, W2D, nRomanPots};

const double particleMass[nParticles] = { 0.13957, 0.497611, 0.93827}; // pion, kaon, proton in GeV /c^2 
const int nTriggers = 17;
const int triggerID[] = { 570209, 570219, 570229, 570701, 570702, 570703, 570704, 570705, 
                  570709, 570711, 570712, 570719, 590701, 590703, 590705, 590708, 590709};
// 570702 RP_UPC // 570712 RP_UPC // 570703 RP_SDT // 570709 RP_ET // 570719 RP_ET // 570701 RP_CPT2 // 570711 RP_CPT2 // 570705 RP_CPT2noBBCL // 570704 RP_Zerobias // 590703 RP_SDT // 590709 RP_ET // 590701 RP_CPT2 // 590705 RP_CPT2noBBCL // 590708 RP_CPTnoBBCL // 570209 JPsi*HTTP // 570219 JPsi*HTTP // 570229 JPsi*HTTP
const int CEPtriggers[] = { 570701, 570705, 570711, 590701, 590705, 590708};

bool ConnectInput(int argc, char** argv, TChain* fileChain);
bool CheckTriggers(StUPCEvent* localupcEvt);
long long GetFileSize(string filename);

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

    if(!ConnectInput(argc, argv, upcChain))
    {
        cout << "Wrong input parameters..." << endl; 
        return 1;
    }

    const string& outputFolder = argv[2];

    //HISTOGRAMS
    //afterthesis analysis
    string ToF_no_pion[2] = {"Exactly 1 other track with ToF", "More than one track with ToF"};
    string ToF_pion[3] = {"Pion pairs without ToF", "Pion pairs with 1 ToF confirmed", "Pion pairs with 2 ToF confirmed"};
    string Hist_type[3] = {"Invariant mass", "Vertices vs invariant mass", "Azimuthal angle"};
    string fixedDAC = "DAC only affects pion candidates";

    ROOT::TThreadedObject<TH1D> invmasshisttab[6];
    for (int i = 0; i < 6; i++){
        string name = "invmasshist" + to_string(i);
        string title = Hist_type[0] + ", " + ToF_no_pion[i/3] + ", " + ToF_pion[i%3] + ", " + fixedDAC + ";inv. mass [GeV];events";
        new(&invmasshisttab[i])ROOT::TThreadedObject<TH1D>(name.c_str(), title.c_str(), 100, 0.42, 0.56);
    }

    ROOT::TThreadedObject<TH2D> invmassvsverticeshisttab[6];
    for (int i = 0; i < 6; i++){
        string name = "invmassvsverticeshist" + to_string(i);
        string title = Hist_type[1] + ", " + ToF_no_pion[i/3] + ", " + ToF_pion[i%3] + ", " + fixedDAC + ";inv. mass [GeV];vertices";
        new(&invmassvsverticeshisttab[i])ROOT::TThreadedObject<TH2D>(name.c_str(), title.c_str(), 100, 0.42, 0.56, 10, 0, 10);
    }

    ROOT::TThreadedObject<TH1D> phianglehisttab[6];
    for (int i = 0; i < 6; i++){
        string name = "phianglehist" + to_string(i);
        string title = Hist_type[2] + ", " + ToF_no_pion[i/3] + ", " + ToF_pion[i%3] + ", " + fixedDAC + ";#phi angle [rad];events";
        new(&phianglehisttab[i])ROOT::TThreadedObject<TH1D>(name.c_str(), title.c_str(), 100, -3.14159, 3.14159);
    }

    ROOT::TThreadedObject<TH1D> pipairanglehist = ROOT::TThreadedObject<TH1D>("pipairanglehist", "Angle between pions in detector FoR;#phi angle [rad];events", 30, 0, 3.1416);

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
        TTreeReaderValue<StRPEvent> StRPEventInstance(myReader, "mRPEvent");
        //variable initialization
        myReader.Next();
        //histograms
        std::shared_ptr<TH1D> invmasshisttabLocal[6];
        std::shared_ptr<TH2D> invmassvsverticeshisttabLocal[6];
        std::shared_ptr<TH1D> phianglehisttabLocal[6];
        std::shared_ptr<TH1D> pipairanglehistLocal;
        for (int i = 0; i < 6; i++){
            invmasshisttabLocal[i] = invmasshisttab[i].Get();
            invmassvsverticeshisttabLocal[i] = invmassvsverticeshisttab[i].Get();
            phianglehisttabLocal[i] = phianglehisttab[i].Get();
        }
        pipairanglehistLocal = pipairanglehist.Get();

        int filtered_entries = 0;
        //for changing branch address
        StUPCEvent* tempUPCpointer = StUPCEventInstance.Get();
        StRPEvent* tempRPpointer = StRPEventInstance.Get();

        TLorentzVector trackVector;
        Double_t verDeltaZ;
        TVector3 pVector;
        vector<TLorentzVector> posPion;
        vector<TLorentzVector> negPion;
        vector<TLorentzVector> posPionSuspicious;
        vector<TLorentzVector> negPionSuspicious;
        do{
            tempUPCpointer = StUPCEventInstance.Get();
            tempRPpointer = StRPEventInstance.Get();

            //TESTS & HISTOGRAMS
            //triggers
            if(!CheckTriggers(StUPCEventInstance.Get())){
                continue;
            }
            
            //one proton each side
            int numberOfTracksPerSide[nSides] = {0, 0};
            for(unsigned int k = 0; k < tempRPpointer->getNumberOfTracks(); ++k)
            {
                // Get pointer to k-th track in Roman Pot data collection
                StUPCRpsTrack *trk = tempRPpointer->getTrack(k);
                trk->setEvent(tempRPpointer);
                // Get ID of a branch in which this k-th track was reconstructed
                int j = trk->branch();
                int side = j<2 ? E : W;
                numberOfTracksPerSide[side]++;
            }
            if(numberOfTracksPerSide[0]!=1 || numberOfTracksPerSide[1]!=1){
                continue;
            }
            
            //four planes on each TrackPoint
            bool areAllPlanesPresent = true;
            for(unsigned int k = 0; k < tempRPpointer->getNumberOfTracks(); ++k)
            {
                // Get pointer to k-th track in Roman Pot data collection
                StUPCRpsTrack *trk = tempRPpointer->getTrack(k);
                trk->setEvent(tempRPpointer);
                // test if all 8 planes are present
                if(trk->planesUsed()!=8){
                    areAllPlanesPresent = false;
                    break;
                }
            }
            if(!areAllPlanesPresent){
                continue;
            }

            //exactly one primary vertex to which TOF tracks match to
            vector<Int_t> primaryVertices;
            Int_t VertexId;
            for (Int_t i = 0; i < tempUPCpointer->getNumberOfTracks(); i++){
                VertexId = tempUPCpointer->getTrack(i)->getVertexId();
                if(tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kTof) && tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kPrimary) && !(find(primaryVertices.begin(), primaryVertices.end(), VertexId)!=primaryVertices.end())){
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
            bool areBothRPTracksWithTime = true;
            Double_t time = 0;
            for(unsigned int k = 0; k < tempRPpointer->getNumberOfTracks(); ++k)
            {
                // Get pointer to k-th track in Roman Pot data collection
                StUPCRpsTrack *trk = tempRPpointer->getTrack(k);
                trk->setEvent(tempRPpointer);
                // Get ID of a branch in which this k-th track was reconstructed
                int j = trk->branch();
                //deletes those tracks which do not have confirmed time
                if(trk->time()<0){
                    areBothRPTracksWithTime = false;
                    break;
                }
                int side = j<2 ? E : W;
                //to not make the loop go twice, also west side is the "+" one
                if(side==E){
                    time -= trk->time();
                }
                if(side==W){
                    time += trk->time();
                }
            }
            if(!areBothRPTracksWithTime){
                continue;
            }
            //*100 at the end it to convert from m to cm
            verDeltaZ = tempUPCpointer->getVertexId(VertexId)->getPosZ()-299792458.0*time/2*100;
            if(abs(verDeltaZ)>36){
                continue;
            }

            //fiducial cuts
            bool passedFiducialCuts = true;
            //px barrier
            Double_t pxbarrier = -0.27;
            //py low nand high barrier
            Double_t pylowbarrier = 0.4;
            Double_t pyhighbarrier = 0.8;
            //p circle parameters
            Double_t pxcenter = -0.6;
            Double_t pycenter = 0;
            Double_t pradius = 1.1;
            for(unsigned int k = 0; k < tempRPpointer->getNumberOfTracks(); ++k){
                // Get pointer to k-th track in Roman Pot data collection
                StUPCRpsTrack *trk = tempRPpointer->getTrack(k);
                trk->setEvent(tempRPpointer);
                pVector = trk->pVec();
                bool f1 = pow(pVector.X()-pxcenter, 2) + pow(pVector.Y()-pycenter, 2) < pradius*pradius;
                bool f2 = pylowbarrier < abs(pVector.Y()) && abs(pVector.Y()) < pyhighbarrier;
                bool f3 = pxbarrier < pVector.X();
                if(!(f1 && f2 && f3)){
                    passedFiducialCuts = false;
                    break;
                }
            }            
            if (!passedFiducialCuts){
                continue;
            }

            //elastic cuts
            TLorentzVector RPsum = {0,0,0,0};
            //checking total momentum
            for(unsigned int k = 0; k < tempRPpointer->getNumberOfTracks(); ++k){
                // Get pointer to k-th track in Roman Pot data collection
                StUPCRpsTrack *trk = tempRPpointer->getTrack(k);
                trk->setEvent(tempRPpointer);
                trackVector.SetVectM(trk->pVec(), particleMass[Proton]);
                RPsum += trackVector;
            }
            //xposition of the gauss
            Double_t xGauss = -0.035;
            //yposition of the gauss
            Double_t yGauss = 0.01;
            //radius of the gauss
            //was 0.07, 0.1 FOR TEST
            Double_t rGauss = 0.07;
            if (pow(RPsum.X()-xGauss, 2) + pow(RPsum.Y()-yGauss, 2) < rGauss*rGauss){
                continue;
            }

            //quality test of central tracks
            //1. Number of hits
            //2. DCA in R and z (the second one also doubles as spread limiter)
            //3. kinematic range for optimal measurement
            int numberOfTOFnotPions = 0;
            for (Int_t i = 0; i < tempUPCpointer->getNumberOfTracks(); i++){
                //quality
                if(tempUPCpointer->getTrack(i)->getNhitsDEdx()<15 || tempUPCpointer->getTrack(i)->getNhitsFit()<25){
                    continue;
                }
                if(tempUPCpointer->getTrack(i)->getPt()<=0.2 || abs(tempUPCpointer->getTrack(i)->getEta())>=0.7){
                    continue;
                }
                if(abs(tempUPCpointer->getTrack(i)->getNSigmasTPCProton())<3 || abs(tempUPCpointer->getTrack(i)->getNSigmasTPCKaon())<3){
                    if(tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kTof) && tempUPCpointer->getTrack(i)->getTofPathLength()>0 && tempUPCpointer->getTrack(i)->getTofTime()>0){
                        numberOfTOFnotPions++;
                    }
                    continue;
                }else if (abs(tempUPCpointer->getTrack(i)->getNSigmasTPCPion())<3){
                    tempUPCpointer->getTrack(i)->getLorentzVector(trackVector, particleMass[Pion]);
                }else{
                    continue;
                }
                //DCA only active for pion candidates, this way it allows more not-pions and thus changes the amounts of additional particles
                if((pow(tempUPCpointer->getTrack(i)->getDcaXY(), 2) + pow(tempUPCpointer->getTrack(i)->getDcaZ(), 2))<1 && tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kPrimary)){
                    continue;
                }
                //final checks; for charge and ToF
                if(tempUPCpointer->getTrack(i)->getCharge()==1){
                    //basic test of ToF
                    if(tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kTof) && tempUPCpointer->getTrack(i)->getTofPathLength()>0 && tempUPCpointer->getTrack(i)->getTofTime()>0){
                        posPion.push_back(trackVector);
                    }else{
                        posPionSuspicious.push_back(trackVector);
                    }
                }else if (tempUPCpointer->getTrack(i)->getCharge()==-1){
                    //basic test of ToF
                    if(tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kTof) && tempUPCpointer->getTrack(i)->getTofPathLength()>0 && tempUPCpointer->getTrack(i)->getTofTime()>0){
                        negPion.push_back(trackVector);
                    }else{
                        negPionSuspicious.push_back(trackVector);
                    }
                }
            }

            if(posPion.size()==0 && negPion.size()==0 && posPionSuspicious.size()==0 && negPionSuspicious.size()==0){
                continue;
            }
            // GUIDE
            // string ToF_no_pion[2] = {"Exactly 1 other track with ToF", "More than one track with ToF"};
            // string ToF_pion[3] = {"Pion pairs without ToF", "Pion pairs with 1 ToF confirmed", "Pion pairs with 2 ToF confirmed"};
            // string Hist_type[2] = {"Invariant mass", "Vertices vs invariant mass"};
            // std::shared_ptr<TH1D> invmasshisttabLocal[6];
            // std::shared_ptr<TH2D> invmassvsverticeshisttabLocal[6];
            //for mass
            TLorentzVector tempSum = {0, 0, 0, 0};
            //invariant mass, pi+ pi-, legit ones
            for (long unsigned int posi = 0; posi < posPion.size(); posi++){
                for (long unsigned int negi = 0; negi < negPion.size(); negi++)
                {
                    tempSum = posPion[posi] + negPion[negi];
                    if(numberOfTOFnotPions==1){
                        invmasshisttabLocal[2]->Fill(tempSum.M());
                        invmassvsverticeshisttabLocal[2]->Fill(tempSum.M(), tempUPCpointer->getNumberOfVertices());
                        phianglehisttabLocal[2]->Fill(tempSum.Phi());
                    }else if (numberOfTOFnotPions>1){
                        invmasshisttabLocal[5]->Fill(tempSum.M());
                        invmassvsverticeshisttabLocal[5]->Fill(tempSum.M(), tempUPCpointer->getNumberOfVertices());
                        phianglehisttabLocal[5]->Fill(tempSum.Phi());
                    }                   
                    pipairanglehistLocal->Fill(posPion[posi].Angle(negPion[negi].Vect()));
                }
            }
            //invariant mass, pi+ pi-, mixed ones
            for (long unsigned int posi = 0; posi < posPionSuspicious.size(); posi++){
                for (long unsigned int negi = 0; negi < negPion.size(); negi++)
                {
                    tempSum = posPionSuspicious[posi] + negPion[negi];
                    if(numberOfTOFnotPions==1){
                        invmasshisttabLocal[1]->Fill(tempSum.M());
                        invmassvsverticeshisttabLocal[1]->Fill(tempSum.M(), tempUPCpointer->getNumberOfVertices());
                        phianglehisttabLocal[1]->Fill(tempSum.Phi());
                    }else if (numberOfTOFnotPions>1){
                        invmasshisttabLocal[4]->Fill(tempSum.M());
                        invmassvsverticeshisttabLocal[4]->Fill(tempSum.M(), tempUPCpointer->getNumberOfVertices());
                        phianglehisttabLocal[4]->Fill(tempSum.Phi());
                    }   
                }
            }
            for (long unsigned int posi = 0; posi < posPion.size(); posi++){
                for (long unsigned int negi = 0; negi < negPionSuspicious.size(); negi++)
                {
                    tempSum = posPion[posi] + negPionSuspicious[negi];
                    if(numberOfTOFnotPions==1){
                        invmasshisttabLocal[1]->Fill(tempSum.M());
                        invmassvsverticeshisttabLocal[1]->Fill(tempSum.M(), tempUPCpointer->getNumberOfVertices());
                        phianglehisttabLocal[1]->Fill(tempSum.Phi());
                    }else if (numberOfTOFnotPions>1){
                        invmasshisttabLocal[4]->Fill(tempSum.M());
                        invmassvsverticeshisttabLocal[4]->Fill(tempSum.M(), tempUPCpointer->getNumberOfVertices());
                        phianglehisttabLocal[4]->Fill(tempSum.Phi());
                    }   
                }
            }
            //invariant mass, pi+ pi-, suspicious ones
            for (long unsigned int posi = 0; posi < posPionSuspicious.size(); posi++){
                for (long unsigned int negi = 0; negi < negPionSuspicious.size(); negi++)
                {
                    tempSum = posPionSuspicious[posi] + negPionSuspicious[negi];
                    if(numberOfTOFnotPions==1){
                        invmasshisttabLocal[0]->Fill(tempSum.M());
                        invmassvsverticeshisttabLocal[0]->Fill(tempSum.M(), tempUPCpointer->getNumberOfVertices());
                        phianglehisttabLocal[0]->Fill(tempSum.Phi());
                    }else if (numberOfTOFnotPions>1){
                        invmasshisttabLocal[3]->Fill(tempSum.M());
                        invmassvsverticeshisttabLocal[3]->Fill(tempSum.M(), tempUPCpointer->getNumberOfVertices());
                        phianglehisttabLocal[3]->Fill(tempSum.Phi());
                    }   
                }
            }

            filtered_entries++;

            posPion.clear();
            negPion.clear();
            posPionSuspicious.clear();
            negPionSuspicious.clear();
        }while (myReader.Next());

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
    std::shared_ptr<TH1D> invmasshisttabFinal[6];
    std::shared_ptr<TH2D> invmassvsverticeshisttabFinal[6];
    std::shared_ptr<TH1D> phianglehisttabFinal[6];
    std::shared_ptr<TH1D> pipairanglehistFinal;
    for (int i = 0; i < 6; i++){
        invmasshisttabFinal[i] = invmasshisttab[i].Merge();
        invmassvsverticeshisttabFinal[i] = invmassvsverticeshisttab[i].Merge();
        phianglehisttabFinal[i] = phianglehisttab[i].Merge();
    }
    pipairanglehistFinal = pipairanglehist.Merge();

    //setting up a tree & output file
    string outfileName = outputFolder + "AnaOutput_NewAnalysis_FixedDCA.root";
    cout<<"Created output file "<<outfileName<<endl;
    TFile* outputFileHist = TFile::Open(outfileName.c_str(), "recreate");

    outputFileHist->cd();
    for (int i = 0; i < 6; i++){
        invmasshisttabFinal[i]->Write();
        invmassvsverticeshisttabFinal[i]->Write();
        phianglehisttabFinal[i]->Write();
    }
    pipairanglehistFinal->Write();
    outputFileHist->Close();

    cout<<"Finished processing "<<endl;
    cout<<"Analyzed total "<<upcChain->GetEntries()<<" entries"<<endl;
    cout<<"Filtered total "<<filtered_entries<<" entries"<<endl;
    cout<<"Ending Analysis... GOOD BYE!"<<endl;
    return 0;
}//main

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
        //for testing if file exists
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
}//ConnectInput

bool CheckTriggers(StUPCEvent* localupcEvt)
{

    bool CPTtrigger = false;
    for(int var = 0; var < nTriggers; ++var)
    {
        if(localupcEvt->isTrigger(triggerID[var]))
        {
            //Checked if it is CPT trigger
            for (int i = 0; i < *(&CEPtriggers + 1) - CEPtriggers; ++i)
                if(triggerID[var] == CEPtriggers[i])
                    CPTtrigger=true;
        }
    }

    return CPTtrigger;
}

long long GetFileSize(string filename){
    struct stat64 stat_buf;
    int rc = stat64(filename.c_str(), &stat_buf);
    return rc==0 ? stat_buf.st_size : -1;
}