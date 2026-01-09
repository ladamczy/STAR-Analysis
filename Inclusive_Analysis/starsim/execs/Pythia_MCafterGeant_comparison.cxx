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
    TFile pythia_input("/home/adam/STAR-Analysis/Inclusive_Analysis/starsim/execs/centralDiffractive.root");
    TFile MCafterGeant_input("/home/adam/STAR-Analysis/Inclusive_Analysis/starsim/execs/centralDiffractive.upcDst.root");
    // TTree* pythia_tree = pythia_input.Get<TTree>("genevents");
    // TTree* MCafterGeant_tree = MCafterGeant_input.Get<TTree>("mUPCTree");
    TTree* pythia_tree = (TTree*)pythia_input.Get("genevents");
    TTree* MCafterGeant_tree = (TTree*)MCafterGeant_input.Get("mUPCTree");
    StarGenEvent* pythia_event = new StarGenEvent();
    StUPCEvent* MCafterGeant_event = new StUPCEvent();
    pythia_tree->SetBranchAddress("primaryEvent", &pythia_event);
    MCafterGeant_tree->SetBranchAddress("mUPCEvent", &MCafterGeant_event);

    //setting up graph histograms
    TH1D pythia_given_particles("pythia_given_particles", "Number of particles given to GEANT for further processing", 50, 0, 50);
    TH1D GEANT_taken_particles("GEANT_taken_particles", "Number of particles first in GEANT event tree", 50, 0, 50);

    //setting up canvas
    TCanvas c1("c1", "c1", 1200, 800);

    //loop to check all the events
    int max_events = max(pythia_tree->GetEntries(), MCafterGeant_tree->GetEntries());
    int offset = 0;
    int making_sure_all_events_were_checked = 0;
    std::vector<int> pythia_pgd_ids, geant_pdg_ids;
    // for(size_t event_number = 0; event_number<99; event_number++){
    for(size_t event_number = 0; event_number<max_events; event_number++){
        pythia_tree->GetEntry(event_number);
        MCafterGeant_tree->GetEntry(event_number-offset);
        pythia_pgd_ids.clear();
        geant_pdg_ids.clear();
        making_sure_all_events_were_checked++;

        //part that synchronises events
        //we loop after pythia events, as there will be more of them surely
        //if pythia event is smaller, that means there is 
        while(pythia_event->GetEventNumber()<MCafterGeant_event->getEventNumber()){
            printf("Offset increase at position %d\n", event_number);
            offset++;
            event_number++;
            pythia_tree->GetEntry(event_number);
            MCafterGeant_tree->GetEntry(event_number-offset);
        }

        //pythia
        int part_pythia = 0;
        for(size_t part_number = 0; part_number<pythia_event->GetNumberOfParticles(); part_number++){
            StarGenParticle* part = (*pythia_event)[part_number];
            if(part->Simulate()&&part->GetStack()!=-1){
                part_pythia++;
                pythia_pgd_ids.push_back(part->GetId());
            }
        }
        pythia_given_particles.Fill(part_pythia);
        //Geant
        int part_geant = 0;
        for(size_t part_number = 0; part_number<MCafterGeant_event->getNumberOfMCParticles(); part_number++){
            TParticle* part = MCafterGeant_event->getMCParticle(part_number);
            if(part->GetFirstMother()==1){
                part_geant++;
                geant_pdg_ids.push_back(part->GetPdgCode());
            }
        }
        GEANT_taken_particles.Fill(part_geant);

        //other
        if(part_geant>part_pythia){
            printf("More particles given to Geant than implied it should be on event %d\n", event_number);
        } else if(part_geant<part_pythia){
            printf("Less particles given to Geant than implied it should be on event %d\n", event_number);
        } else{
            for(size_t i = 0; i<pythia_pgd_ids.size(); i++){
                if(pythia_pgd_ids[i]!=geant_pdg_ids[i]){
                    printf("MISMATCHED IDS at event %d\n", event_number);
                }
            }
        }
    }

    //summary
    printf("Events checked: %d\n", making_sure_all_events_were_checked);

    //saving histograms
    TFile output("~/output.root", "RECREATE");
    output.cd();
    pythia_given_particles.Write();
    GEANT_taken_particles.Write();
    output.Close();

    //loop to check different events
    while(true){
        int pythia_event_number, MCafterGeant_event_number;
        printf("Write event numbers; PYTHIA, then Geant\n\n");
        scanf("%d %d", &pythia_event_number, &MCafterGeant_event_number);
        pythia_tree->GetEntry(pythia_event_number);
        MCafterGeant_tree->GetEntry(MCafterGeant_event_number);
        printf("\nPYTHIA particles:\t%d\nGeant particles:\t%d\n\n", pythia_event->GetNumberOfParticles(), MCafterGeant_event->getNumberOfMCParticles());
        c1.Clear();

        //Pythia particles graph
        TGraph ParticlesPYTHIA(0);
        ParticlesPYTHIA.SetNameTitle("PYTHIA", "PYTHIA;#eta;#phi");
        ParticlesPYTHIA.SetMarkerStyle(20);
        ParticlesPYTHIA.SetMarkerSize(2);
        ParticlesPYTHIA.SetMarkerColor(4);
        for(size_t particle_index = 0; particle_index<pythia_event->GetNumberOfParticles(); particle_index++){
            StarGenParticle* temp = (*pythia_event)[particle_index];
            if(!temp->Simulate()){
                continue;
            }
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
            TVector3 tempVec(tempPart->Px(), tempPart->Py(), tempPart->Pz());
            ParticlesGeant.AddPoint(tempVec.Eta(), tempVec.Phi());
        }

        //drawing and saving canvas
        ParticlesPYTHIA.RemovePoint(0);
        bool drawnMC = false;
        if(ParticlesPYTHIA.GetN()){
            ParticlesPYTHIA.Draw("ap");
            ParticlesPYTHIA.GetXaxis()->SetLimits(-6.0, 6.0);
            ParticlesPYTHIA.GetHistogram()->SetMinimum(-TMath::Pi());
            ParticlesPYTHIA.GetHistogram()->SetMaximum(TMath::Pi());
            drawnMC = true;
        }
        if(ParticlesGeant.GetN()){
            if(drawnMC){
                ParticlesGeant.Draw("same p");
            } else{
                ParticlesGeant.Draw("ap");
                ParticlesGeant.GetXaxis()->SetLimits(-6.0, 6.0);
                ParticlesGeant.GetHistogram()->SetMinimum(-TMath::Pi());
                ParticlesGeant.GetHistogram()->SetMaximum(TMath::Pi());
            }
        }
        c1.BuildLegend();
        gPad->Modified();
        gPad->Update();

        //writing event
        //Pythia part
        printf("Full Pythia event:\n");
        pythia_event->Print();
        printf("\nPart passed from Pythia to Geant:\n");
        pythia_event->Print("simu");
        //Geant part
        printf("\nFull Geant event:\n");
        printf("Geant event number:\t%d\n", MCafterGeant_event->getEventNumber());
        for(int particle_index = 0; particle_index<MCafterGeant_event->getNumberOfMCParticles(); particle_index++){
            TParticle* tempPart = MCafterGeant_event->getMCParticle(particle_index);
            tempPart->Print();
        }
    }

    a.Run();
    return 0;
}
