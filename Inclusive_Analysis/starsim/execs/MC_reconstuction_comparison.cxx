#include <iostream>

//ROOT headers
#include <ROOT/TThreadedObject.hxx>
#include <TTreeReader.h>
#include <ROOT/TTreeProcessorMT.hxx>
#include <TGraph.h>
#include <TCanvas.h>
#include <TParticlePDG.h>

//STAR headers
#include <StarGenEvent.h>
#include <StarGenParticle.h>

// picoDst headers
#include "StRPEvent.h"
#include "StUPCRpsTrack.h"
#include "StUPCRpsTrackPoint.h"
#include "StUPCEvent.h"
#include "StUPCTrack.h"
#include "StUPCBemcCluster.h"
#include "StUPCVertex.h"
#include "StUPCTofHit.h"

//my headers
#include <ProcessingInsideLoop.h>
#include <ProcessingOutsideLoop.h>
#include <UsefulThings.h>

int main(int argc, char* argv[]){

    //argv[1] - input file
    //argv[2] - output folder
    //argv[3] - #of cores

    int nthreads = 1;
    if(argc==4){
        nthreads = atoi(argv[3]);
    }

    cout<<"Program is running on "<<nthreads<<" threads"<<endl;
    ROOT::EnableThreadSafety();
    //actually i'm not sure if it's needed here
    // ROOT::EnableImplicitMT(nthreads); //turn on multicore processing

    //preparing input & output
    TChain* eventFiles = new TChain("mUPCTree");
    if(ConnectInput(argc, argv, eventFiles)){
        cout<<"All files connected"<<endl;
    }
    const string& outputFolder = argv[2];

    //Useful IDs
    const int K0sPDGid = 310;
    const int K0sbarPDGid = -310;
    const int LambdaPDGid = 3122;
    const int LambdabarPDGid = -3122;
    const int phiPDGid = 333;
    const int phibarPDGid = -333;
    const int KstarPDGid = 313;
    const int KstarbarPDGid = -313;
    const int piplusPDGid = 211;
    const int piminusPDGid = -211;
    const int KplusPDGid = 321;
    const int KminusPDGid = -321;
    const int pplusPDGid = 2212;
    const int pminusPDGid = -2212;

    //histograms
    // ProcessingOutsideLoop outsideprocessing;
    //mixing TOF between events proved to be a failure
    //mass histograms with TOF first and mixing event pairs after
    // outsideprocessing.AddHistogram(TH1D("Name", "Name of simulated particles;id;Number of particles", 1, 0, 1));
    // outsideprocessing.AddHistogram(TH2D("etapTK0S", "K^{0}_{S} number;eta;p_{T}", 10, -1, 1, 10, 0, 2.5));
    // outsideprocessing.AddHistogram(TH2D("etapTLambda", "#Lambda^{0} number;eta;p_{T}", 10, -1, 1, 10, 0, 2.5));
    // outsideprocessing.AddHistogram(TH2D("etapTLambdabar", "#bar{#Lambda}^{0} number;eta;p_{T}", 10, -1, 1, 10, 0, 2.5));
    // outsideprocessing.AddHistogram(TH2D("etapTKstar", "K^{*}(892) number;eta;p_{T}", 10, -1, 1, 10, 0, 2.5));
    // outsideprocessing.AddHistogram(TH2D("etapTKstarbar", "#bar{K}^{*}(892) number;eta;p_{T}", 10, -1, 1, 10, 0, 2.5));
    // outsideprocessing.AddHistogram(TH2D("etapTphi", "#varphi(1020) number;eta;p_{T}", 10, -1, 1, 10, 0, 2.5));

    TH1D MCParticles("MCParticles", "MC number of particles;particles;events", 25, 0, 25);
    TH1D TPCParticles("TPCParticles", "TPC number of particles;particles;events", 25, 0, 25);
    TH1D DifferenceParticles("DifferenceParticles", "(MC - TPC) difference in number of particles;difference;events", 20, -10, 10);
    TH1D Distance("Distance", "Distance MC-TPC in #eta-#phi space;distance;track pairs", 100, 0, 5);

    //processing
    //defining TreeProcessor
    // ROOT::TTreeProcessorMT TreeProc(*eventFiles, nthreads);

    int tempCounter = -1;

    //defining processing function
    // auto myFunction = [&](TTreeReader& myReader){

    //setting up TTreeReader without multithread processing

    TTreeReader myReader(eventFiles);

    //getting values from TChain, in-loop histogram initialization
    TTreeReaderValue<StUPCEvent> StUPCEventInstance(myReader, "mUPCEvent");
    TTreeReaderValue<StRPEvent> StRPEventInstance(myReader, "mRPEvent");
    // ProcessingInsideLoop insideprocessing;
    StUPCEvent* tempUPCpointer;
    StRPEvent* tempRPpointer;
    // insideprocessing.GetLocalHistograms(&outsideprocessing);

    std::vector<TParticle*> positiveMC;
    std::vector<TParticle*> negativeMC;
    std::vector<StUPCTrack*> positiveTrack;
    std::vector<StUPCTrack*> negativeTrack;

    while(myReader.Next()){
        tempUPCpointer = StUPCEventInstance.Get();
        tempRPpointer = StRPEventInstance.Get();
        tempCounter++;
        //cleaning
        positiveMC.clear();
        negativeMC.clear();
        positiveTrack.clear();
        negativeTrack.clear();
        //below is the loop

        printf("Event %d:\n", tempCounter);
        //filtering MC particles
        printf("Number of MC particles (excluding protons) before filter: %d\n", tempUPCpointer->getNumberOfMCParticles()-2);
        for(size_t i = 0; i<tempUPCpointer->getNumberOfMCParticles(); i++){
            TParticle* temp = tempUPCpointer->getMCParticle(i);
            //filtering
            if(fabs(temp->Eta())>0.9||temp->Pt()<0.2){
                continue;
            }
            if(temp->GetPDG()->Charge()==0||abs(temp->GetPdgCode())<10){
                //if not charged or gluon/quark
                continue;
            }
            //filling
            if(temp->GetPDG()->Charge()>0){
                positiveMC.push_back(temp);
            } else{
                negativeMC.push_back(temp);
            }
        }
        printf("Number of MC particles (excluding protons) after filter: %d\n", positiveMC.size()+negativeMC.size());

        //filtering tracks
        printf("Number of tracks (excluding protons) before filter: %d\n", tempUPCpointer->getNumberOfTracks());
        for(int i = 0; i<tempUPCpointer->getNumberOfTracks(); i++){
            StUPCTrack* tempTrack = tempUPCpointer->getTrack(i);
            if(!tempTrack->getFlag(StUPCTrack::kTof)){
                continue;
            }
            // if(!tempTrack->getFlag(StUPCTrack::kPrimary)){
            //     continue;
            // }
            //filtering
            if(fabs(tempTrack->getEta())>0.9||tempTrack->getPt()<0.2){
                continue;
            }
            //filling
            if(tempTrack->getCharge()>0){
                positiveTrack.push_back(tempTrack);
            } else{
                negativeTrack.push_back(tempTrack);
            }
        }
        printf("Number of tracks (excluding protons) after filter: %d\n", positiveTrack.size()+negativeTrack.size());

        //histograms
        MCParticles.Fill(positiveMC.size()+negativeMC.size());
        TPCParticles.Fill(positiveTrack.size()+negativeTrack.size());
        DifferenceParticles.Fill(int(positiveMC.size()+negativeMC.size())-int(positiveTrack.size()+negativeTrack.size()));
        if(positiveMC.size()>=positiveTrack.size()&&negativeMC.size()>=negativeTrack.size()){
            for(size_t i = 0; i<positiveMC.size(); i++){
                for(size_t j = 0; j<positiveTrack.size(); j++){
                    //getting distance in eta-phi space
                    TVector3 vecMC, vecTPC;
                    vecMC.SetXYZ(positiveMC[i]->Px(), positiveMC[i]->Py(), positiveMC[i]->Pz());
                    positiveTrack[j]->getMomentum(vecTPC);
                    Distance.Fill(vecMC.DrEtaPhi(vecTPC));
                }
            }
        }





        //special part where one event is drawn
        if(tempCounter==nthreads){
            TCanvas c1("c1", "c1", 1200, 800);

            //MC particles graph
            TGraph ParticlesMC(0);
            ParticlesMC.SetNameTitle("MC", "MC;#eta;#phi");
            ParticlesMC.SetMarkerStyle(20);
            ParticlesMC.SetMarkerSize(2);
            ParticlesMC.SetMarkerColor(4);
            printf("MC positive:\n");
            for(size_t particle_index = 0; particle_index<positiveMC.size(); particle_index++){
                TParticle* temp = positiveMC[particle_index];
                ParticlesMC.AddPoint(temp->Eta(), temp->Phi());
                printf("Eta:\t%f,\tPhi:\t%f\n", temp->Eta(), temp->Phi());
            }
            printf("MC negative:\n");
            for(size_t particle_index = 0; particle_index<negativeMC.size(); particle_index++){
                TParticle* temp = negativeMC[particle_index];
                ParticlesMC.AddPoint(temp->Eta(), temp->Phi());
                printf("Eta:\t%f,\tPhi:\t%f\n", temp->Eta(), temp->Phi());
            }

            //TPC particles graph
            TGraph ParticlesTPC(0);
            ParticlesTPC.SetNameTitle("TPC", "TPC;#eta;#phi");
            ParticlesTPC.SetMarkerStyle(21);
            ParticlesTPC.SetMarkerSize(1.5);
            ParticlesTPC.SetMarkerColor(2);
            printf("Track positive:\n");
            for(int particle_index = 0; particle_index<positiveTrack.size(); particle_index++){
                StUPCTrack* tempTrack = positiveTrack[particle_index];
                double corrected_phi = tempTrack->getPhi()<0 ? tempTrack->getPhi()+2*TMath::Pi() : tempTrack->getPhi();
                ParticlesTPC.AddPoint(tempTrack->getEta(), corrected_phi);
                printf("Eta:\t%f,\tPhi:\t%f\n", tempTrack->getEta(), corrected_phi);
            }
            printf("Track negative:\n");
            for(int particle_index = 0; particle_index<negativeTrack.size(); particle_index++){
                StUPCTrack* tempTrack = negativeTrack[particle_index];
                double corrected_phi = tempTrack->getPhi()<0 ? tempTrack->getPhi()+2*TMath::Pi() : tempTrack->getPhi();
                ParticlesTPC.AddPoint(tempTrack->getEta(), corrected_phi);
                printf("Eta:\t%f,\tPhi:\t%f\n", tempTrack->getEta(), corrected_phi);
            }

            //drawing and saving canvas (with protection against empty graphs)
            bool drawnMC = false;
            if(ParticlesMC.GetN()){
                ParticlesMC.Draw("ap");
                ParticlesMC.GetXaxis()->SetLimits(-1.0, 1.0);
                ParticlesMC.GetHistogram()->SetMinimum(0.0);
                ParticlesMC.GetHistogram()->SetMaximum(2*TMath::Pi());
                drawnMC = true;
            }
            if(ParticlesTPC.GetN()){
                if(drawnMC){
                    ParticlesTPC.Draw("same p");
                } else{
                    ParticlesTPC.Draw("ap");
                    ParticlesTPC.GetXaxis()->SetLimits(-1.0, 1.0);
                    ParticlesTPC.GetHistogram()->SetMinimum(0.0);
                    ParticlesTPC.GetHistogram()->SetMaximum(2*TMath::Pi());
                }
            }
            c1.BuildLegend();
            c1.SetTitle("Matching test");
            c1.SaveAs("Matching.png");
        }

    }
    //event loop finish

//     //lambda finish
//     return 0;
// };

// TreeProc.Process(myFunction);

// outsideprocessing.Merge();
// outsideprocessing.GetPointerAfterMerge1D("Name")->LabelsDeflate();
// outsideprocessing.GetPointerAfterMerge1D("Name")->SetMinimum(0);
// outsideprocessing.GetPointerAfterMerge1D("Name")->LabelsOption("a", "X");

    //setting up a tree & output file
    string path = string(argv[0]);
    string outfileName;
    if(outputFolder.find(".root")!=std::string::npos){
        outfileName = outputFolder;
    } else{
        outfileName = outputFolder+"SimOutput_"+path.substr(path.find_last_of("/\\")+1)+".root";
    }
    cout<<"Created output file "<<outfileName<<endl;

    TFile* outputFileHist = TFile::Open(outfileName.c_str(), "recreate");
    // outsideprocessing.SaveToFile(outputFileHist);
    MCParticles.Write();
    TPCParticles.Write();
    DifferenceParticles.Write();
    Distance.Write();
    outputFileHist->Close();

    return 0;
}
