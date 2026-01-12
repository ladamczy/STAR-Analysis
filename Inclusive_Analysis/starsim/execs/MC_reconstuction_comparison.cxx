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

    TH1D FlowOfEvents("FlowOfEvents", "Independent event checks (of MC particles);;events", 1, 0, 1);
    TH1D K0SIndex("K0SIndex", "Index of K^{0}_{S};index;particles", 30, 0, 30);
    TH3D NumberOfPions("NumberOfPions", "Number of pions in event;positive;negative;neutral", 5, 0, 5, 5, 0, 5, 5, 0, 5);
    TH1D MCParticles("MCParticles", "MC number of particles;particles;events", 25, 0, 25);
    TH1D TPCParticles("TPCParticles", "TPC number of particles;particles;events", 25, 0, 25);
    TH1D DifferenceParticles("DifferenceParticles", "(MC - TPC) difference in number of particles;difference;events", 20, -10, 10);
    TH1D Distance("Distance", "Distance MC-TPC in #eta-#phi space;distance;track pairs", 120, 0, 6);
    TH1D DistanceCloser("DistanceCloser", "Distance MC-TPC in #eta-#phi space;distance;track pairs", 100, 0, 0.5);
    TH1D MpipiTPC("MpipiTPC", "#pi^{+}#pi^{-} pair mass;m_{#pi^{+}#pi^{-}} [GeV];pairs", 40, 0.4, 0.6);
    TH1D MpipiTPCExtremelyWide("MpipiTPCExtremelyWide", "#pi^{+}#pi^{-} pair mass;m_{#pi^{+}#pi^{-}} [GeV];pairs", 300, 0.0, 3.0);
    TH1D MpipiMC("MpipiMC", "#pi^{+}#pi^{-} pair mass (only K^{0}_{S} decay products);m_{#pi^{+}#pi^{-}} [GeV];pairs", 40, 0.4, 0.6);
    TH1D MpipiMCExtremelyWide("MpipiMCExtremelyWide", "#pi^{+}#pi^{-} pair mass (only K^{0}_{S} decay products);m_{#pi^{+}#pi^{-}} [GeV];pairs", 300, 0.0, 3.0);
    TH1D MpipiFlow("MpipiFlow", "#pi^{+}#pi^{-} pairs after MC cuts;;pairs", 1, 0, 1);
    TH2D MpipiPairs("MpipiPairs", "particle pairs;positive;negative", 1, 0, 1, 1, 0, 1);
    TH2D MpipiMothers("MpipiMothers", "particle number;positive;negative", 20, 0, 20, 20, 0, 20);
    TH1D MpipiMotherName("MpipiMotherName", "#pi^{+}#pi^{-} pair mothers;;mothers", 1, 0, 1);
    TH1D MpipiDeltaT("MpipiDeltaT", "Time between pion decay;#Delta t_{0} [ns];events", 100, -5, 5);
    TH1D MpipiDeltaTOnlyTwoTracksDetected("MpipiDeltaTOnlyTwoTracksDetected", "Time between pion decay;#Delta t_{0} [ns];events", 100, -5, 5);
    TH1D MpipiDeltaTMoreThanTwoTracksDetected("MpipiDeltaTMoreThanTwoTracksDetected", "Time between pion decay;#Delta t_{0} [ns];events", 100, -5, 5);
    TH1D MpipiAfterDeltaT("MpipiAfterDeltaT", "#pi^{+}#pi^{-} pair mass after #Delta t_{0} cut;m_{#pi^{+}#pi^{-}} [GeV];pairs", 40, 0.4, 0.6);
    TH1D MpipiAfterDeltaTExtremelyWide("MpipiAfterDeltaTExtremelyWide", "#pi^{+}#pi^{-} pair mass after #Delta t_{0} cut;m_{#pi^{+}#pi^{-}} [GeV];pairs", 300, 0.0, 3.0);
    TH2D MpipiAfterDeltaTMass("MpipiAfterDeltaTMass", "#pi^{+}#pi^{-} pair mass after #Delta t_{0} cut;m_{#pi^{+}#pi^{-}} [GeV];pairs", 40, 0.4, 0.6, 100, -5, 5);
    TH1D MpipiAfterDeltaTNotPassed("MpipiAfterDeltaTNotPassed", "#pi^{+}#pi^{-} pair mass after failed #Delta t_{0} cut;m_{#pi^{+}#pi^{-}} [GeV];pairs", 40, 0.4, 0.6);
    TH1D MpipiAfterDeltaTNotPassedExtremelyWide("MpipiAfterDeltaTNotPassedExtremelyWide", "#pi^{+}#pi^{-} pair mass after failed #Delta t_{0} cut;m_{#pi^{+}#pi^{-}} [GeV];pairs", 300, 0.0, 3.0);
    TH2D MpipiAfterDeltaTNotPassedMass("MpipiAfterDeltaTNotPassedMass", "#pi^{+}#pi^{-} pair mass after failed #Delta t_{0} cut;m_{#pi^{+}#pi^{-}} [GeV];pairs", 40, 0.4, 0.6, 50, -5, 5);
    TH1D MpipiMCnotK0SMother("MpipiMCnotK0SMother", "#pi^{+}#pi^{-} pair mass (everything except K^{0}_{S} mother verification);m_{#pi^{+}#pi^{-}} [GeV];pairs", 40, 0.4, 0.6);
    TH1D MpipiMCnotK0SMotherExtremelyWide("MpipiMCnotK0SMotherExtremelyWide", "#pi^{+}#pi^{-} pair mass (everything except K^{0}_{S} mother verification);m_{#pi^{+}#pi^{-}} [GeV];pairs", 300, 0.0, 3.0);

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
    // TTreeReaderValue<StRPEvent> StRPEventInstance(myReader, "mRPEvent");
    // ProcessingInsideLoop insideprocessing;
    StUPCEvent* tempUPCpointer;
    // StRPEvent* tempRPpointer;
    // insideprocessing.GetLocalHistograms(&outsideprocessing);

    std::vector<TParticle*> positiveMC;
    std::vector<TParticle*> negativeMC;
    std::vector<StUPCTrack*> positiveTrack;
    std::vector<StUPCTrack*> negativeTrack;
    std::vector<int> chosen_MC_list_positive;
    std::vector<int> chosen_MC_list_negative;
    std::vector<int> chosen_MC_list_copy;

    while(myReader.Next()){
        tempUPCpointer = StUPCEventInstance.Get();
        // tempRPpointer = StRPEventInstance.Get();
        tempCounter++;
        //cleaning
        positiveMC.clear();
        negativeMC.clear();
        positiveTrack.clear();
        negativeTrack.clear();
        chosen_MC_list_positive.clear();
        chosen_MC_list_negative.clear();
        //below is the loop

        //filtering MC particles
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

        //filtering tracks
        for(int i = 0; i<tempUPCpointer->getNumberOfTracks(); i++){
            StUPCTrack* tempTrack = tempUPCpointer->getTrack(i);
            //TODO check with & without TOF
            if(!tempTrack->getFlag(StUPCTrack::kTof)){
                continue;
            }
            // ALL TRACKS ARE PRIMARY
            // if(tempTrack->getFlag(StUPCTrack::kPrimary)){
            //     continue;
            // }
            //filtering
            if(fabs(tempTrack->getEta())>0.9||tempTrack->getPt()<0.2){
                continue;
            }
            if(tempTrack->getNhitsFit()<=20||tempTrack->getNhitsDEdx()<15){
                continue;
            }
            //filling
            if(tempTrack->getCharge()>0){
                positiveTrack.push_back(tempTrack);
            } else{
                negativeTrack.push_back(tempTrack);
            }
        }

        bool duplicated_matches_exist = false;

        //histograms
        //flow
        FlowOfEvents.Fill("All events", 1.0);
        //for proper order
        FlowOfEvents.Fill("2x K^{0}_{S}", 0.0);
        FlowOfEvents.Fill("4x#pi^{#pm/0}", 0.0);
        FlowOfEvents.Fill("2x#pi^{+}, 2x#pi^{-}", 0.0);
        int n_of_K0S = 0;
        int n_of_piplus = 0;
        int n_of_piminus = 0;
        int n_of_pizero = 0;
        for(size_t i = 0; i<tempUPCpointer->getNumberOfMCParticles(); i++){
            if(abs(tempUPCpointer->getMCParticle(i)->GetPdgCode())==310){
                n_of_K0S++;
                K0SIndex.Fill(i);
            }
            if(tempUPCpointer->getMCParticle(i)->GetPdgCode()==piplusPDGid)
                n_of_piplus++;
            if(tempUPCpointer->getMCParticle(i)->GetPdgCode()==piminusPDGid)
                n_of_piminus++;
            if(abs(tempUPCpointer->getMCParticle(i)->GetPdgCode())==111)
                n_of_pizero++;
        }
        NumberOfPions.Fill(n_of_piplus, n_of_piminus, n_of_pizero);
        if(n_of_K0S==2)
            FlowOfEvents.Fill("2x K^{0}_{S}", 1.0);
        if(n_of_piplus==2&&n_of_piminus==2)
            FlowOfEvents.Fill("2x#pi^{+}, 2x#pi^{-}", 1.0);
        if(n_of_piplus+n_of_piminus+n_of_pizero==4)
            FlowOfEvents.Fill("4x#pi^{#pm/0}", 1.0);
        //the rest
        MCParticles.Fill(positiveMC.size()+negativeMC.size());
        TPCParticles.Fill(positiveTrack.size()+negativeTrack.size());
        DifferenceParticles.Fill(int(positiveMC.size()+negativeMC.size())-int(positiveTrack.size()+negativeTrack.size()));
        //loops for fitting the tracks
        //positive
        for(size_t i = 0; i<positiveTrack.size(); i++){
            double min_distance = 100;
            int min_distance_id = -1;
            for(size_t j = 0; j<positiveMC.size(); j++){
                //getting distance in eta-phi space
                TVector3 vecMC, vecTPC;
                vecMC.SetXYZ(positiveMC[j]->Px(), positiveMC[j]->Py(), positiveMC[j]->Pz());
                positiveTrack[i]->getMomentum(vecTPC);
                double etaphi_distance = vecMC.DrEtaPhi(vecTPC);
                //distance cutoff
                if(etaphi_distance>0.1){
                    continue;
                }
                Distance.Fill(etaphi_distance);
                DistanceCloser.Fill(etaphi_distance);
                //checking if the distance to j-th MC track is the smallest
                min_distance = min(min_distance, etaphi_distance);
                if(min_distance==etaphi_distance){
                    min_distance_id = j;
                }
            }
            //checking for which MC track we have the nearest match
            chosen_MC_list_positive.push_back(min_distance_id);
        }
        chosen_MC_list_copy = chosen_MC_list_positive;
        //removing all the -1
        chosen_MC_list_copy.erase(std::remove(chosen_MC_list_copy.begin(), chosen_MC_list_copy.end(), -1), chosen_MC_list_copy.end());
        std::sort(chosen_MC_list_copy.begin(), chosen_MC_list_copy.end());
        const auto duplicate_positive = std::adjacent_find(chosen_MC_list_copy.begin(), chosen_MC_list_copy.end());
        if(duplicate_positive!=chosen_MC_list_copy.end()){
            printf("Event with conflict (positive) = %d\n", tempCounter);
            for(size_t index = 0; index<chosen_MC_list_positive.size(); index++){
                printf("%d, ", chosen_MC_list_positive[index]);
            }
            printf("\n");
            duplicated_matches_exist = true;
        }
        //negative
        for(size_t i = 0; i<negativeTrack.size(); i++){
            double min_distance = 100;
            int min_distance_id = -1;
            for(size_t j = 0; j<negativeMC.size(); j++){
                //getting distance in eta-phi space
                TVector3 vecMC, vecTPC;
                vecMC.SetXYZ(negativeMC[j]->Px(), negativeMC[j]->Py(), negativeMC[j]->Pz());
                negativeTrack[i]->getMomentum(vecTPC);
                double etaphi_distance = vecMC.DrEtaPhi(vecTPC);
                //distance cutoff
                if(etaphi_distance>0.1){
                    continue;
                }
                Distance.Fill(etaphi_distance);
                DistanceCloser.Fill(etaphi_distance);
                //checking if the distance to j-th MC track is the smallest
                min_distance = min(min_distance, etaphi_distance);
                if(min_distance==etaphi_distance){
                    min_distance_id = j;
                }
            }
            //checking for which MC track we have the nearest match
            chosen_MC_list_negative.push_back(min_distance_id);
        }
        chosen_MC_list_copy = chosen_MC_list_negative;
        //removing all the -1
        chosen_MC_list_copy.erase(std::remove(chosen_MC_list_copy.begin(), chosen_MC_list_copy.end(), -1), chosen_MC_list_copy.end());
        std::sort(chosen_MC_list_copy.begin(), chosen_MC_list_copy.end());
        const auto duplicate_negative = std::adjacent_find(chosen_MC_list_copy.begin(), chosen_MC_list_copy.end());
        if(duplicate_negative!=chosen_MC_list_copy.end()){
            printf("Event with conflict (negative) = %d\n", tempCounter);
            for(size_t index = 0; index<chosen_MC_list_negative.size(); index++){
                printf("%d, ", chosen_MC_list_negative[index]);
            }
            printf("\n");
            duplicated_matches_exist = true;
        }

        //checking the mass-fit of pairs when no duplicates
        if(!duplicated_matches_exist){
            for(size_t i = 0; i<positiveTrack.size(); i++){
                for(size_t j = 0; j<negativeTrack.size(); j++){
                    //check if the track is associated with anything
                    if(chosen_MC_list_positive[i]<0||chosen_MC_list_negative[j]<0){
                        continue;
                    }
                    //tracks
                    TLorentzVector posTrack, negTrack;
                    positiveTrack[i]->getLorentzVector(posTrack, 0.13957);
                    negativeTrack[j]->getLorentzVector(negTrack, 0.13957);
                    MpipiTPC.Fill((posTrack+negTrack).M());
                    MpipiTPCExtremelyWide.Fill((posTrack+negTrack).M());
                    //MC particles
                    //choosing MC particle associated with these particular tracks
                    int posPDG = positiveMC[chosen_MC_list_positive[i]]->GetPdgCode();
                    int posProductionVertex = positiveMC[chosen_MC_list_positive[i]]->GetFirstMother();
                    int negPDG = negativeMC[chosen_MC_list_negative[j]]->GetPdgCode();
                    int negProductionVertex = negativeMC[chosen_MC_list_negative[j]]->GetFirstMother();
                    //loop for finding mother particle
                    int posMotherPDG;
                    for(size_t MCindex = 0; MCindex<tempUPCpointer->getNumberOfMCParticles(); MCindex++){
                        if(tempUPCpointer->getMCParticle(MCindex)->GetFirstDaughter()==posProductionVertex){
                            posMotherPDG = tempUPCpointer->getMCParticle(MCindex)->GetPdgCode();
                            break;
                        }
                    }
                    //both pions, same mother, mother is K0S
                    MpipiFlow.Fill("TPC", 1.0);
                    MpipiPairs.Fill(positiveMC[chosen_MC_list_positive[i]]->GetPDG()->GetName(), negativeMC[chosen_MC_list_negative[j]]->GetPDG()->GetName(), 1.0);
                    if(posPDG==211&&negPDG==-211){
                        MpipiFlow.Fill("Pion pair", 1.0);
                        MpipiMothers.Fill(posProductionVertex, negProductionVertex);
                        if(posProductionVertex==negProductionVertex){
                            MpipiFlow.Fill("Same mother", 1.0);
                            for(size_t MCindex = 0; MCindex<tempUPCpointer->getNumberOfMCParticles(); MCindex++){
                                if(tempUPCpointer->getMCParticle(MCindex)->GetFirstDaughter()==posProductionVertex){
                                    MpipiMotherName.Fill(tempUPCpointer->getMCParticle(MCindex)->GetPDG()->GetName(), 1.0);
                                    break;
                                }
                            }
                            if(abs(posMotherPDG)==310){
                                MpipiFlow.Fill("K^{0}_{S} mother", 1.0);
                                MpipiMC.Fill((posTrack+negTrack).M());
                                MpipiMCExtremelyWide.Fill((posTrack+negTrack).M());
                                double deltaT = DeltaT0(positiveTrack[i], negativeTrack[j], 0.13957, 0.13957);
                                MpipiDeltaT.Fill(deltaT);
                                if(positiveTrack.size()==1&&negativeTrack.size()==1){
                                    MpipiDeltaTOnlyTwoTracksDetected.Fill(deltaT);
                                } else{
                                    //if no tracks either positive or negative
                                    //then no loops, so the only alternative to 1&1
                                    //is more tracks, so no checking necessary
                                    MpipiDeltaTMoreThanTwoTracksDetected.Fill(deltaT);
                                }
                                //value nicked from .txt file with actual data
                                if(fabs(deltaT)<3*0.13124272289253383){
                                    MpipiAfterDeltaT.Fill((posTrack+negTrack).M());
                                    MpipiAfterDeltaTExtremelyWide.Fill((posTrack+negTrack).M());
                                    MpipiAfterDeltaTMass.Fill((posTrack+negTrack).M(), deltaT);
                                } else{
                                    MpipiAfterDeltaTNotPassed.Fill((posTrack+negTrack).M());
                                    MpipiAfterDeltaTNotPassedExtremelyWide.Fill((posTrack+negTrack).M());
                                    MpipiAfterDeltaTNotPassedMass.Fill((posTrack+negTrack).M(), deltaT);
                                }
                            } else{
                                MpipiMCnotK0SMother.Fill((posTrack+negTrack).M());
                                MpipiMCnotK0SMotherExtremelyWide.Fill((posTrack+negTrack).M());
                            }
                        }
                    }
                }
            }
        }


        //special part where one event is drawn
        if(tempCounter==nthreads){
            //statistics
            printf("Event %d:\n", tempCounter);
            printf("Number of MC particles (excluding protons) before filter: %d\n", tempUPCpointer->getNumberOfMCParticles()-2);
            printf("Number of MC particles (excluding protons) after filter: %d\n", positiveMC.size()+negativeMC.size());
            printf("Number of tracks (excluding protons) before filter: %d\n", tempUPCpointer->getNumberOfTracks());
            printf("Number of tracks (excluding protons) after filter: %d\n", positiveTrack.size()+negativeTrack.size());

            //drawing
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
                TVector3 tempVec(temp->Px(), temp->Py(), temp->Pz());
                ParticlesMC.AddPoint(tempVec.Eta(), tempVec.Phi());
                printf("Eta:\t%f,\tPhi:\t%f\n", tempVec.Eta(), tempVec.Phi());
            }
            printf("MC negative:\n");
            for(size_t particle_index = 0; particle_index<negativeMC.size(); particle_index++){
                TParticle* temp = negativeMC[particle_index];
                TVector3 tempVec(temp->Px(), temp->Py(), temp->Pz());
                ParticlesMC.AddPoint(tempVec.Eta(), tempVec.Phi());
                printf("Eta:\t%f,\tPhi:\t%f\n", tempVec.Eta(), tempVec.Phi());
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
                TVector3 tempVec;
                tempTrack->getMomentum(tempVec);
                ParticlesTPC.AddPoint(tempVec.Eta(), tempVec.Phi());
                printf("Eta:\t%f,\tPhi:\t%f\n", tempVec.Eta(), tempVec.Phi());
            }
            printf("Track negative:\n");
            for(int particle_index = 0; particle_index<negativeTrack.size(); particle_index++){
                StUPCTrack* tempTrack = negativeTrack[particle_index];
                TVector3 tempVec;
                tempTrack->getMomentum(tempVec);
                ParticlesTPC.AddPoint(tempVec.Eta(), tempVec.Phi());
                printf("Eta:\t%f,\tPhi:\t%f\n", tempVec.Eta(), tempVec.Phi());
            }

            //drawing and saving canvas (with protection against empty graphs)
            bool drawnMC = false;
            if(ParticlesMC.GetN()){
                ParticlesMC.Draw("ap");
                ParticlesMC.GetXaxis()->SetLimits(-1.0, 1.0);
                ParticlesMC.GetHistogram()->SetMinimum(-TMath::Pi());
                ParticlesMC.GetHistogram()->SetMaximum(TMath::Pi());
                drawnMC = true;
            }
            if(ParticlesTPC.GetN()){
                if(drawnMC){
                    ParticlesTPC.Draw("same p");
                } else{
                    ParticlesTPC.Draw("ap");
                    ParticlesTPC.GetXaxis()->SetLimits(-1.0, 1.0);
                    ParticlesTPC.GetHistogram()->SetMinimum(-TMath::Pi());
                    ParticlesTPC.GetHistogram()->SetMaximum(TMath::Pi());
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
    FlowOfEvents.LabelsDeflate();
    FlowOfEvents.Write();
    K0SIndex.Write();
    NumberOfPions.Write();
    MCParticles.Write();
    TPCParticles.Write();
    DifferenceParticles.Write();
    Distance.Write();
    DistanceCloser.Write();
    MpipiTPC.Write();
    MpipiTPCExtremelyWide.Write();
    MpipiMC.Write();
    MpipiMCExtremelyWide.Write();
    MpipiFlow.LabelsDeflate();
    MpipiFlow.Write();
    MpipiPairs.LabelsDeflate("X");
    MpipiPairs.LabelsDeflate("Y");
    MpipiPairs.Write();
    MpipiMothers.Write();
    MpipiMotherName.LabelsDeflate();
    MpipiMotherName.Write();
    MpipiDeltaT.Write();
    MpipiDeltaTOnlyTwoTracksDetected.Write();
    MpipiDeltaTMoreThanTwoTracksDetected.Write();
    MpipiAfterDeltaT.Write();
    MpipiAfterDeltaTExtremelyWide.Write();
    MpipiAfterDeltaTMass.Write();
    MpipiAfterDeltaTNotPassed.Write();
    MpipiAfterDeltaTNotPassedExtremelyWide.Write();
    MpipiAfterDeltaTNotPassedMass.Write();
    MpipiMCnotK0SMother.Write();
    MpipiMCnotK0SMotherExtremelyWide.Write();
    outputFileHist->Close();

    return 0;
}
