#include <iostream>

//ROOT headers
#include <ROOT/TThreadedObject.hxx>
#include <TTreeReader.h>
#include <ROOT/TTreeProcessorMT.hxx>
#include <TGraph.h>
#include <TCanvas.h>
#include <TParticlePDG.h>
#include <TApplication.h>

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
#include <UsefulThings.h>

int main(int argc, char* argv[]){

    TApplication a("a", 0, 0);
    //because it was spamming the stdout with warnings
    gErrorIgnoreLevel = kError;

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

    //setting up input trees
    TFile pythia_input("/home/adam/STAR-Analysis/Inclusive_Analysis/starsim/execs/centralDiffractive0pythia.root");
    TFile MCafterGeant_input("/home/adam/STAR-Analysis/Inclusive_Analysis/starsim/execs/centralDiffractive0upcDst.root");
    // TTree* pythia_tree = pythia_input.Get<TTree>("genevents");
    // TTree* MCafterGeant_tree = MCafterGeant_input.Get<TTree>("mUPCTree");
    TTree* pythia_tree = (TTree*)pythia_input.Get("genevents");
    TTree* MCafterGeant_tree = (TTree*)MCafterGeant_input.Get("mUPCTree");
    StarGenEvent* pythia_event = new StarGenEvent();
    StUPCEvent* MCafterGeant_event = new StUPCEvent();
    pythia_tree->SetBranchAddress("primaryEvent", &pythia_event);
    MCafterGeant_tree->SetBranchAddress("mUPCEvent", &MCafterGeant_event);

    //setting up canvas
    TCanvas c1("c1", "c1", 1200, 800);


    //loop to check different events
    while(true){
        int pythia_event_number, MCafterGeant_event_number;
        printf("Write event numbers; PYTHIA, then Geant\n");
        scanf("%d %d", &pythia_event_number, &MCafterGeant_event_number);
        pythia_tree->GetEntry(pythia_event_number);
        MCafterGeant_tree->GetEntry(MCafterGeant_event_number);
        printf("PYTHIA particles:\t%d\nGeant particles:\t%d\n\n\n\n\n", pythia_event->GetNumberOfParticles(), MCafterGeant_event->getNumberOfMCParticles());
        c1.Clear();

        //Pythia particles graph
        TGraph ParticlesPYTHIA(0);
        ParticlesPYTHIA.SetNameTitle("PYTHIA", "PYTHIA;#eta;#phi");
        ParticlesPYTHIA.SetMarkerStyle(20);
        ParticlesPYTHIA.SetMarkerSize(2);
        ParticlesPYTHIA.SetMarkerColor(4);
        for(size_t particle_index = 0; particle_index<pythia_event->GetNumberOfParticles(); particle_index++){
            StarGenParticle* temp = (*pythia_event)[particle_index];
            TVector3 tempVec(temp->GetPx(), temp->GetPy(), temp->GetPz());
            ParticlesPYTHIA.AddPoint(tempVec.Eta(), tempVec.Phi());
        }
        //Geant particles graph
        TGraph ParticlesGeant(0);
        ParticlesGeant.SetNameTitle("Geant", "Geant;#eta;#phi");
        ParticlesGeant.SetMarkerStyle(21);
        ParticlesGeant.SetMarkerSize(1.5);
        ParticlesGeant.SetMarkerColor(2);
        for(int particle_index = 0; particle_index<MCafterGeant_event->getNumberOfMCParticles(); particle_index++){
            TParticle* tempPart = MCafterGeant_event->getMCParticle(particle_index);
            tempPart->Print();
            TVector3 tempVec(tempPart->Px(), tempPart->Py(), tempPart->Pz());
            ParticlesGeant.AddPoint(tempVec.Eta(), tempVec.Phi());
        }

        //drawing and saving canvas
        ParticlesPYTHIA.RemovePoint(0);
        bool drawnMC = false;
        if(ParticlesPYTHIA.GetN()){
            ParticlesPYTHIA.Draw("ap");
            ParticlesPYTHIA.GetXaxis()->SetLimits(-1.0, 1.0);
            ParticlesPYTHIA.GetHistogram()->SetMinimum(-TMath::Pi());
            ParticlesPYTHIA.GetHistogram()->SetMaximum(TMath::Pi());
            drawnMC = true;
        }
        if(ParticlesGeant.GetN()){
            if(drawnMC){
                ParticlesGeant.Draw("same p");
            } else{
                ParticlesGeant.Draw("ap");
                ParticlesGeant.GetXaxis()->SetLimits(-1.0, 1.0);
                ParticlesGeant.GetHistogram()->SetMinimum(-TMath::Pi());
                ParticlesGeant.GetHistogram()->SetMaximum(TMath::Pi());
            }
        }
        c1.BuildLegend();
        gPad->Modified();
        gPad->Update();

        printf("Full event:\n");
        pythia_event->Print();
        printf("\nPart passed to Geant:\n");
        pythia_event->Print("simu");
        printf("\n\n\n\nEvent number:\t%d\n", MCafterGeant_event->getEventNumber());

    }

    a.Run();
    return 0;
}
