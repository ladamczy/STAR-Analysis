#include <iostream>

//ROOT headers
#include <ROOT/TThreadedObject.hxx>
#include <TTreeReader.h>
#include <ROOT/TTreeProcessorMT.hxx>
#include <TGraph.h>
#include <TCanvas.h>
#include <TParticlePDG.h>
#include <TLine.h>
#include <TEfficiency.h>

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
#include "StPicoPhysicalHelix.h"

//my headers
#include <ProcessingInsideLoop.h>
#include <ProcessingOutsideLoop.h>
#include <UsefulThings.h>

void PrintBigger(TParticle* input, std::string additional_stuff = "");

int main(int argc, char* argv[]){

    //argv[1] - input file
    //argv[2] - output folder/file
    //argv[3] - PDG id/codename of particle tested
    //argv[4] - event to print
    if(argc>5||argc<4){
        printf("Invalid number of arguments. Proper argument usage:\n");
        printf("argv[1] - input file\n");
        printf("argv[2] - output folder/file\n");
        printf("argv[3] - PDG id/codename of particle tested\n");
        printf("argv[4] - event to print (optional)\n");
        printf("\n");
        printf("Particle\tPDG code\tcodename\n");
        printf("K0S\t\t310\t\tK0S\n");
        printf("Lambda0\t\t3122\t\tLambda0\n");
        printf("Lambda0bar\t-3122\t\tLambda0bar\n");
        printf("K*(892)\t\t313\t\tKstar\n");
        printf("K*(892)bar\t-313\t\tKstarbar\n");
        printf("phi(1020)\t333\t\tphi\n");
        return 0;
    }

    //Useful IDs
    const int K0SPDGid = 310;
    const int LambdaPDGid = 3122;
    const int LambdabarPDGid = -3122;
    const int KstarPDGid = 313;
    const int KstarbarPDGid = -313;
    const int phiPDGid = 333;
    const int piplusPDGid = 211;
    const int piminusPDGid = -211;
    const int KplusPDGid = 321;
    const int KminusPDGid = -321;
    const int pplusPDGid = 2212;
    const int pminusPDGid = -2212;
    //masses
    std::map<int, double> massmap = { {K0SPDGid, 0.497611},
                                    {LambdaPDGid, 1.115683},
                                    {LambdabarPDGid, 1.115683},
                                    {KstarPDGid, 0.89167},
                                    {KstarbarPDGid, 0.89167},
                                    {phiPDGid, 1.019461},
                                    {piplusPDGid, 0.139570},
                                    {piminusPDGid, 0.139570},
                                    {KplusPDGid, 0.493677},
                                    {KminusPDGid, 0.493677},
                                    {pplusPDGid, 0.938272},
                                    {pminusPDGid, 0.938272} };

    int PDGmain, eventToPrint = -1;
    int PDGpositive, PDGnegative;
    if(argc==5){
        eventToPrint = atoi(argv[4]);
    }
    //atoi(not number) returns 0
    //so there is a nice check if argv[3] is a PDG id or a codename
    PDGmain = atoi(argv[3]);
    //consult the map for the key of particle codename
    //and fill PDGmain with proper id
    if(PDGmain==0){
        //why is there no switch for strings, I will never fathom
        //instead i have std::map with custom comparison function for c strings
        auto PDGmap = std::map<const char*, int, std::function<bool(const char*, const char*)>>{
            [](const char* a, const char* b){
                return strcmp(a,b)<0;
            }
        };
        PDGmap = { {"K0S", K0SPDGid},
            {"Lambda0", LambdaPDGid},
            {"Lambda0bar", LambdabarPDGid},
            {"Kstar", KstarPDGid},
            {"Kstatbar", KstarbarPDGid},
            {"phi", phiPDGid} };
        if(PDGmap.count(argv[3])==0){
            printf("This codename (%s) is not implemented/invalid\n", argv[3]);
            return 1;
        } else{
            PDGmain = PDGmap[argv[3]];
        }
    }

    printf("Chosen main particle PDG id: %d\n", PDGmain);

    //setting particle PDG number and decay products
    switch(PDGmain){
    case K0SPDGid:
        PDGpositive = piplusPDGid;
        PDGnegative = piminusPDGid;
        break;
    case LambdaPDGid:
        PDGpositive = pplusPDGid;
        PDGnegative = piminusPDGid;
        break;
    case LambdabarPDGid:
        PDGpositive = piplusPDGid;
        PDGnegative = pminusPDGid;
        break;
    case KstarPDGid:
        PDGpositive = KplusPDGid;
        PDGnegative = piminusPDGid;
        break;
    case KstarbarPDGid:
        PDGpositive = piplusPDGid;
        PDGnegative = KminusPDGid;
        break;
    case phiPDGid:
        PDGpositive = KplusPDGid;
        PDGnegative = KminusPDGid;
        break;
    default:
        printf("This PDG number is not implemented\n");
        return 1;
    }

    //preparing input & output
    TChain* eventFiles = new TChain("mUPCTree");
    if(ConnectInput(argc, argv, eventFiles)){
        cout<<"All files connected"<<endl;
    }
    const string& outputFolder = argv[2];

    //histograms
    TH1D FlowChart("FlowChart", "Number of #phi(1020) passing criteria", 1, 0, 1);
    TH1D Chi2Signal("Chi2Signal", "#chi^{2} characteristic of the signal", 100, 0, 100);
    TH1D Chi2All("Chi2All", "#chi^{2} characteristic of the signal+background", 100, 0, 100);


    //setting up TTreeReader without multithread processing
    TTreeReader myReader(eventFiles);
    int tempCounter = -1;

    //getting values from TChain, in-loop histogram initialization
    TTreeReaderValue<StUPCEvent> StUPCEventInstance(myReader, "mUPCEvent");
    StUPCEvent* tempUPCpointer;

    //TODO change for proper stuff
    int searchedPDGid = 333;

    while(myReader.Next()){
        tempUPCpointer = StUPCEventInstance.Get();
        tempCounter++;
        if(tempCounter%10000==0){
            printf("Event %d analysed\n", tempCounter);
        }
        //cleaning

        //loop for tracks with confirmed proper parent
        for(size_t MC_main_particle_index = 0; MC_main_particle_index<tempUPCpointer->getNumberOfMCParticles(); MC_main_particle_index++){
            // looking for a decayed particle
            if(tempUPCpointer->getMCParticle(MC_main_particle_index)->GetPdgCode()!=searchedPDGid){
                continue;
            }
            FlowChart.Fill("All #phi(1020)", 1.0);


            //looking for decay products correlated to TPC tracks
            std::vector<int> MC_decay;
            MC_decay.clear();
            for(size_t MC_product_particle_index = 0; MC_product_particle_index<tempUPCpointer->getNumberOfMCParticles(); MC_product_particle_index++){
                if(tempUPCpointer->getMCParticle(MC_main_particle_index)->GetFirstDaughter()==tempUPCpointer->getMCParticle(MC_product_particle_index)->GetFirstMother()){
                    MC_decay.push_back(MC_product_particle_index);
                }
            }
            std::vector<StUPCTrack*> TPC_tracks;
            TPC_tracks.clear();
            for(size_t TPC_track_index = 0; TPC_track_index<tempUPCpointer->getNumberOfTracks(); TPC_track_index++){
                if(std::find(MC_decay.begin(), MC_decay.end(), tempUPCpointer->getTrack(TPC_track_index)->getIdTruth()-1)!=MC_decay.end()){
                    TPC_tracks.push_back(tempUPCpointer->getTrack(TPC_track_index));
                }
            }
            if(TPC_tracks.size()!=2){
                continue;
            }
            FlowChart.Fill("Decay products detected in TPC", 1.0);


            //quality of tracks
            bool allTracksFine = true;
            for(size_t TPC_track_index = 0; TPC_track_index<TPC_tracks.size(); TPC_track_index++){
                if(TPC_tracks[TPC_track_index]->getNhitsFit()<=20||TPC_tracks[TPC_track_index]->getNhitsDEdx()<15){
                    allTracksFine = false;
                    break;
                }
            }
            if(!allTracksFine){
                continue;
            }
            FlowChart.Fill("Decay products properly reconstructed", 1.0);


            //fiducial cuts
            bool allTracksWithinFiducial = true;
            for(size_t TPC_track_index = 0; TPC_track_index<TPC_tracks.size(); TPC_track_index++){
                if(fabs(TPC_tracks[TPC_track_index]->getEta())>0.9||TPC_tracks[TPC_track_index]->getPt()<0.2){
                    allTracksWithinFiducial = false;
                    break;
                }
            }
            if(!allTracksWithinFiducial){
                continue;
            }
            FlowChart.Fill("Decay products inside TPC fiducial region", 1.0);


            //proper TOF flags
            bool allTracksWithProperTOFQualities = true;
            for(size_t TPC_track_index = 0; TPC_track_index<TPC_tracks.size(); TPC_track_index++){
                if(!TPC_tracks[TPC_track_index]->getFlag(StUPCTrack::kTof)){
                    allTracksWithProperTOFQualities = false;
                    break;
                }
                if(TPC_tracks[TPC_track_index]->getTofPathLength()<=0){
                    allTracksWithProperTOFQualities = false;
                    break;
                }
                if(TPC_tracks[TPC_track_index]->getTofTime()<=0){
                    allTracksWithProperTOFQualities = false;
                    break;
                }
            }
            if(!allTracksWithProperTOFQualities){
                continue;
            }
            FlowChart.Fill("Decay products properly detected in TOF", 1.0);


            //deltaT section
            //TODO fix that
            double deltaT = DeltaT0(TPC_tracks[0], TPC_tracks[1], massmap[KplusPDGid], massmap[KplusPDGid]);
            //fixing the +-1ns peaks
            if(fabs(deltaT-1)<deltaT){
                deltaT -= 1;
            } else if(fabs(deltaT+1)<deltaT){
                deltaT += 1;
            }
            //TODO change for better
            double sigmaT = 0.13555533188873461;
            double sigmaParticle1 = TPC_tracks[0]->getNSigmasTPCKaon();
            double sigmaParticle2 = TPC_tracks[1]->getNSigmasTPCKaon();
            double Chi2 = pow(deltaT/sigmaT, 2)+sigmaParticle1*sigmaParticle1+sigmaParticle2*sigmaParticle2;
            Chi2Signal.Fill(Chi2);
        }


        //loop for all track pairs
        for(int i = 0; i<tempUPCpointer->getNumberOfTracks()-1; i++){
            //check if track okay
            StUPCTrack* tempTrack1 = tempUPCpointer->getTrack(i);
            if(!tempTrack1->getFlag(StUPCTrack::kTof)){
                continue;
            }
            if(tempTrack1->getTofPathLength()<=0){
                continue;
            }
            if(tempTrack1->getTofTime()<=0){
                continue;
            }
            if(fabs(tempTrack1->getEta())>0.9||tempTrack1->getPt()<0.2){
                continue;
            }
            if(tempTrack1->getNhitsFit()<=20||tempTrack1->getNhitsDEdx()<15){
                continue;
            }
            for(int j = i+1; j<tempUPCpointer->getNumberOfTracks(); j++){
                //check if track okay
                StUPCTrack* tempTrack2 = tempUPCpointer->getTrack(j);
                if(!tempTrack2->getFlag(StUPCTrack::kTof)){
                    continue;
                }
                if(tempTrack2->getTofPathLength()<=0){
                    continue;
                }
                if(tempTrack2->getTofTime()<=0){
                    continue;
                }
                if(fabs(tempTrack2->getEta())>0.9||tempTrack2->getPt()<0.2){
                    continue;
                }
                if(tempTrack2->getNhitsFit()<=20||tempTrack2->getNhitsDEdx()<15){
                    continue;
                }

                //actual filling
                //deltaT section
                //TODO fix that
                double deltaT = DeltaT0(tempTrack1, tempTrack2, massmap[KplusPDGid], massmap[KplusPDGid]);
                //fixing the +-1ns peaks
                if(fabs(deltaT-1)<deltaT){
                    deltaT -= 1;
                } else if(fabs(deltaT+1)<deltaT){
                    deltaT += 1;
                }
                //TODO change for better
                double sigmaT = 0.13555533188873461;
                double sigmaParticle1 = tempTrack1->getNSigmasTPCKaon();
                double sigmaParticle2 = tempTrack2->getNSigmasTPCKaon();
                double Chi2 = pow(deltaT/sigmaT, 2)+sigmaParticle1*sigmaParticle1+sigmaParticle2*sigmaParticle2;
                Chi2All.Fill(Chi2);
            }
        }








        //special part where one event is drawn
        // if(tempCounter==eventToPrint){
        //     //statistics
        //     printf("Number of MC particles (excluding protons) before filter: %d\n", tempUPCpointer->getNumberOfMCParticles()-2);
        //     printf("Number of MC particles (excluding protons) after filter: %d\n", positiveMC.size()+negativeMC.size());
        //     printf("Number of tracks (excluding protons) before filter: %d\n", tempUPCpointer->getNumberOfTracks());
        //     printf("Number of tracks (excluding protons) after filter: %d\n", positiveTrack.size()+negativeTrack.size());
        //     for(size_t MCindex = 0; MCindex<tempUPCpointer->getNumberOfMCParticles(); MCindex++){
        //         //we only know IdTruth because there is great care taken not to change the order of particles
        //         //it is NOT written into upcDst file!!!
        //         PrintBigger(tempUPCpointer->getMCParticle(MCindex), "\tIdTruth:\t"+to_string(MCindex+1));
        //     }


        //     //drawing
        //     TCanvas c1("c1", "c1", 1200, 800);

        //     //MC particles graph
        //     TGraph ParticlesMC(0);
        //     ParticlesMC.SetNameTitle("MC", "MC;#eta;#phi");
        //     ParticlesMC.SetMarkerStyle(20);
        //     ParticlesMC.SetMarkerSize(2);
        //     ParticlesMC.SetMarkerColor(4);
        //     printf("MC positive:\n");
        //     for(size_t particle_index = 0; particle_index<positiveMC.size(); particle_index++){
        //         TParticle* temp = positiveMC[particle_index];
        //         TVector3 tempVec(temp->Px(), temp->Py(), temp->Pz());
        //         ParticlesMC.AddPoint(tempVec.Eta(), tempVec.Phi());
        //         printf("Eta:\t%f,\tPhi:\t%f\n", tempVec.Eta(), tempVec.Phi());
        //     }
        //     printf("MC negative:\n");
        //     for(size_t particle_index = 0; particle_index<negativeMC.size(); particle_index++){
        //         TParticle* temp = negativeMC[particle_index];
        //         TVector3 tempVec(temp->Px(), temp->Py(), temp->Pz());
        //         ParticlesMC.AddPoint(tempVec.Eta(), tempVec.Phi());
        //         printf("Eta:\t%f,\tPhi:\t%f\n", tempVec.Eta(), tempVec.Phi());
        //     }

        //     //TPC particles graph
        //     TGraph ParticlesTPC(0);
        //     ParticlesTPC.SetNameTitle("TPC", "TPC;#eta;#phi");
        //     ParticlesTPC.SetMarkerStyle(21);
        //     ParticlesTPC.SetMarkerSize(1.5);
        //     ParticlesTPC.SetMarkerColor(2);
        //     printf("Track positive:\n");
        //     for(int particle_index = 0; particle_index<positiveTrack.size(); particle_index++){
        //         StUPCTrack* tempTrack = positiveTrack[particle_index];
        //         TVector3 tempVec;
        //         tempTrack->getMomentum(tempVec);
        //         ParticlesTPC.AddPoint(tempVec.Eta(), tempVec.Phi());
        //         printf("Eta:\t%f,\tPhi:\t%f\n", tempVec.Eta(), tempVec.Phi());
        //     }
        //     printf("Track negative:\n");
        //     for(int particle_index = 0; particle_index<negativeTrack.size(); particle_index++){
        //         StUPCTrack* tempTrack = negativeTrack[particle_index];
        //         TVector3 tempVec;
        //         tempTrack->getMomentum(tempVec);
        //         ParticlesTPC.AddPoint(tempVec.Eta(), tempVec.Phi());
        //         printf("Eta:\t%f,\tPhi:\t%f\n", tempVec.Eta(), tempVec.Phi());
        //     }

        //     //drawing and saving canvas (with protection against empty graphs)
        //     bool drawnMC = false;
        //     if(ParticlesMC.GetN()){
        //         ParticlesMC.Draw("ap");
        //         ParticlesMC.GetXaxis()->SetLimits(-1.0, 1.0);
        //         ParticlesMC.GetHistogram()->SetMinimum(-TMath::Pi());
        //         ParticlesMC.GetHistogram()->SetMaximum(TMath::Pi());
        //         drawnMC = true;
        //     }
        //     if(ParticlesTPC.GetN()){
        //         if(drawnMC){
        //             ParticlesTPC.Draw("same p");
        //         } else{
        //             ParticlesTPC.Draw("ap");
        //             ParticlesTPC.GetXaxis()->SetLimits(-1.0, 1.0);
        //             ParticlesTPC.GetHistogram()->SetMinimum(-TMath::Pi());
        //             ParticlesTPC.GetHistogram()->SetMaximum(TMath::Pi());
        //         }
        //     }
        //     c1.BuildLegend();
        //     c1.SetTitle("Matching test");
        //     c1.SaveAs("Matching.png");
        // }

    }
    //event loop finish

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

    //saving histograms
    FlowChart.LabelsDeflate();
    FlowChart.Write();
    Chi2Signal.Write();
    Chi2All.Write();

    outputFileHist->Close();

    return 0;
}

void PrintBigger(TParticle* input, std::string additional_stuff){
    Printf("TParticle: %-13s  p: %8f %8f %8f \tVertex: %8e %8e %8e \tProd. Vertex: %5d %5d \tDecay Vertex: %5d \tTOF tray:%5d \tTOF module:%5d%s",
        input->GetName(), input->Px(), input->Py(), input->Pz(), input->Vx(), input->Vy(), input->Vz(),
        input->GetFirstMother(), input->GetSecondMother(), input->GetFirstDaughter(),
        input->GetLastDaughter()/100, input->GetLastDaughter()%100, additional_stuff.c_str());
}