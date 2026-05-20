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
    //histogram names
    std::map<int, std::string> namemap = { {K0SPDGid, "K^{0}_{S}"},
                                    {LambdaPDGid, "#Lambda^0"},
                                    {LambdabarPDGid, "#bar{#Lambda}^{0}"},
                                    {KstarPDGid, "K^{*}(892)"},
                                    {KstarbarPDGid, "#bar{K}^{*}(892)"},
                                    {phiPDGid, "#phi(1020)"},
                                    {piplusPDGid, "#pi^{+}"},
                                    {piminusPDGid, "#pi^{-}"},
                                    {KplusPDGid, "K^{+}"},
                                    {KminusPDGid, "K^{-}"},
                                    {pplusPDGid, "p^{+}"},
                                    {pminusPDGid, "p^{-}"} };
    //deltaT width times
    std::map<int, double> sigmaTmap = { {K0SPDGid, 0.13124272289253383},
                                {LambdaPDGid, 0.16761990972301405},
                                {LambdabarPDGid, 0.15911655471031255},
                                {KstarPDGid, 0.1251159662169776},
                                {KstarbarPDGid, 0.12956977399721734},
                                {phiPDGid, 0.13707309430776241} };

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
    TH1D FlowChart("FlowChart", ("Number of "+namemap[PDGmain]+" passing criteria;;particles").c_str(), 1, 0, 1);
    TH1D Chi2Signal("Chi2Signal", "#chi^{2} characteristic of the signal;#chi^{2};pairs of tracks", 100, 0, 100);
    TH1D Chi2All("Chi2All", "#chi^{2} characteristic of the signal+background;#chi^{2};pairs of tracks", 100, 0, 100);
    std::vector<TH1D*> histograms;
    for(auto&& sign:{ 1, -1 }){
        for(auto&& MCparticle:{ 0,1,2,3,4 }){
            for(auto&& particle:{ 0,1,2,3 }){
                std::string mainTitle = "NSigmaTestMC";
                std::string axisTitle = "n#sigma_{";
                switch(MCparticle){
                case 0:
                    mainTitle += "Electron";
                    break;
                case 1:
                    mainTitle += "Pion";
                    break;
                case 2:
                    mainTitle += "Kaon";
                    break;
                case 3:
                    mainTitle += "Proton";
                    break;
                default:
                    mainTitle += "Unidentified";
                    break;
                }
                mainTitle += "TPCPID";
                switch(particle){
                case 0:
                    axisTitle += "e";
                    mainTitle += "Electron";
                    break;
                case 1:
                    axisTitle += "#pi";
                    mainTitle += "Pion";
                    break;
                case 2:
                    axisTitle += "K";
                    mainTitle += "Kaon";
                    break;
                case 3:
                    axisTitle += "p";
                    mainTitle += "Proton";
                    break;
                default:
                    axisTitle += "X";
                    mainTitle += "Mysterious";
                    break;
                }
                axisTitle += "^{";
                switch(sign){
                case 1:
                    axisTitle += "+";
                    mainTitle += "Positive";
                    break;
                case -1:
                    axisTitle += "-";
                    mainTitle += "Negative";
                    break;
                default:
                    axisTitle += "#pm";
                    mainTitle += "Unknown";
                    break;
                }
                axisTitle += "}}";
                histograms.push_back(new TH1D(mainTitle.c_str(), (axisTitle+" histogram for tests;"+axisTitle+";tracks").c_str(), 320, -40, 40));
            }
        }
    }
    TH1D UnidentifiedPositive("UnidentifiedPositive", "Names of unidentified particles;;counts", 1, 0, 1);
    TH1D UnidentifiedNegative("UnidentifiedNegative", "Names of unidentified particles;;counts", 1, 0, 1);

    //setting up TTreeReader without multithread processing
    TTreeReader myReader(eventFiles);
    int tempCounter = -1;

    //getting values from TChain, in-loop histogram initialization
    TTreeReaderValue<StUPCEvent> StUPCEventInstance(myReader, "mUPCEvent");
    StUPCEvent* tempUPCpointer;

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
            if(tempUPCpointer->getMCParticle(MC_main_particle_index)->GetPdgCode()!=PDGmain){
                continue;
            }
            FlowChart.Fill(("All "+namemap[PDGmain]).c_str(), 1.0);


            //looking for decay products correlated to TPC tracks
            std::vector<int> MC_decay;
            MC_decay.clear();
            //loop through all MC particles in search of main particle daughters
            for(size_t MC_product_particle_index = 0; MC_product_particle_index<tempUPCpointer->getNumberOfMCParticles(); MC_product_particle_index++){
                //mother and daughter here refer to vertex numbers!!!
                int motherDecayVertex = tempUPCpointer->getMCParticle(MC_main_particle_index)->GetFirstDaughter();
                int daugterProductionVertex = tempUPCpointer->getMCParticle(MC_product_particle_index)->GetFirstMother();
                //if we found a daughter (a particle whose production vertex is its mothers' decay vertex) we save its MC container number
                if(motherDecayVertex==daugterProductionVertex){
                    MC_decay.push_back(MC_product_particle_index);
                }
            }
            std::vector<StUPCTrack*> TPC_tracks;
            TPC_tracks.clear();
            TPC_tracks.assign(MC_decay.size(), nullptr);
            for(size_t TPC_track_index = 0; TPC_track_index<tempUPCpointer->getNumberOfTracks(); TPC_track_index++){
                //check which index MC track correlated with this tpc track has in MC_decay table
                //if it has one, we write StUPCTrack pointer at the same location 
                int MC_decay_index = std::find(MC_decay.begin(), MC_decay.end(), tempUPCpointer->getTrack(TPC_track_index)->getIdTruth()-1)-MC_decay.begin();
                if(MC_decay_index<MC_decay.size()){
                    TPC_tracks[MC_decay_index] = (tempUPCpointer->getTrack(TPC_track_index));
                }
            }
            //if we find that there exist MC track that has no TPC track assigned (nullptr) we skip this one
            if(std::find(TPC_tracks.begin(), TPC_tracks.end(), nullptr)!=TPC_tracks.end()){
                continue;
            }
            FlowChart.Fill("#splitline{Decay products detected in TPC}{(MC track has associated StUPCTrack)}", 1.0);


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
            FlowChart.Fill("#splitline{Decay products properly reconstructed}{(NhitsFit>20 & NhitsDEdx#geq15)}", 1.0);


            //fiducial cuts
            bool allTracksWithinFiducial = true;
            for(size_t TPC_track_index = 0; TPC_track_index<TPC_tracks.size(); TPC_track_index++){
                if(fabs(TPC_tracks[TPC_track_index]->getEta())>0.9||TPC_tracks[TPC_track_index]->getPt()<0.2||TPC_tracks[TPC_track_index]->getCharge()==0){
                    allTracksWithinFiducial = false;
                    break;
                }
            }
            if(!allTracksWithinFiducial){
                continue;
            }
            FlowChart.Fill("#splitline{Decay products inside TPC fiducial region}{(|#eta|<0.9 & p_{T}>0.2 & q #neq 0)}", 1.0);


            //TOF flag
            bool allTracksWithTOFFlag = true;
            for(size_t TPC_track_index = 0; TPC_track_index<TPC_tracks.size(); TPC_track_index++){
                if(!TPC_tracks[TPC_track_index]->getFlag(StUPCTrack::kTof)){
                    allTracksWithTOFFlag = false;
                    break;
                }
            }
            if(!allTracksWithTOFFlag){
                continue;
            }
            FlowChart.Fill("#splitline{Decay products detected in TOF}{(kToF flag present)}", 1.0);


            //proper TOF measurements
            bool allTracksWithProperTOFMeasurements = true;
            for(size_t TPC_track_index = 0; TPC_track_index<TPC_tracks.size(); TPC_track_index++){
                if(TPC_tracks[TPC_track_index]->getTofPathLength()<=0){
                    allTracksWithProperTOFMeasurements = false;
                    break;
                }
                if(TPC_tracks[TPC_track_index]->getTofTime()<=0){
                    allTracksWithProperTOFMeasurements = false;
                    break;
                }
            }
            if(!allTracksWithProperTOFMeasurements){
                continue;
            }
            FlowChart.Fill("#splitline{Decay products properly detected in TOF}{(kToF flag present, time & distance >0)}", 1.0);


            //deltaT section
            //only for decay into 2 tracks
            if(TPC_tracks.size()!=2){
                continue;
            }
            StUPCTrack* positiveTrack;
            StUPCTrack* negativeTrack;
            TParticle* positiveMCParticle;
            TParticle* negativeMCParticle;
            if(TPC_tracks[0]->getCharge()>0){
                positiveTrack = TPC_tracks[0];
                positiveMCParticle = tempUPCpointer->getMCParticle(MC_decay[0]);
                negativeTrack = TPC_tracks[1];
                negativeMCParticle = tempUPCpointer->getMCParticle(MC_decay[1]);
            } else{
                positiveTrack = TPC_tracks[1];
                positiveMCParticle = tempUPCpointer->getMCParticle(MC_decay[1]);
                negativeTrack = TPC_tracks[0];
                negativeMCParticle = tempUPCpointer->getMCParticle(MC_decay[0]);
            }
            //checking if particles are correct
            if(positiveMCParticle->GetPdgCode()!=PDGpositive||negativeMCParticle->GetPdgCode()!=PDGnegative){
                continue;
            }
            //the rest of deltaT0 calculations
            double deltaT = DeltaT0(positiveTrack, negativeTrack, massmap[PDGpositive], massmap[PDGnegative]);
            //fixing the +-1ns peaks
            if(fabs(deltaT-1)<fabs(deltaT)){
                deltaT -= 1;
            } else if(fabs(deltaT+1)<fabs(deltaT)){
                deltaT += 1;
            }
            double sigmaT = sigmaTmap[PDGmain];
            double sigmaParticle1, sigmaParticle2;
            switch(PDGpositive){
            case piplusPDGid:
                sigmaParticle1 = positiveTrack->getNSigmasTPCPion();
                break;
            case KplusPDGid:
                sigmaParticle1 = positiveTrack->getNSigmasTPCKaon();
                break;
            case pplusPDGid:
                sigmaParticle1 = positiveTrack->getNSigmasTPCProton();
                break;
            default:
                sigmaParticle1 = 0;
                break;
            }
            switch(PDGnegative){
            case piminusPDGid:
                sigmaParticle2 = negativeTrack->getNSigmasTPCPion();
                break;
            case KminusPDGid:
                sigmaParticle2 = negativeTrack->getNSigmasTPCKaon();
                break;
            case pminusPDGid:
                sigmaParticle2 = negativeTrack->getNSigmasTPCProton();
                break;
            default:
                sigmaParticle2 = 0;
                break;
            }

            //chi2 histogram
            double Chi2 = pow(deltaT/sigmaT, 2)+pow(sigmaParticle1, 2)+pow(sigmaParticle2, 2);
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
                //check if the pair is of opposite charges
                if(tempTrack1->getCharge()*tempTrack2->getCharge()>=0){
                    continue;
                }

                //deltaT section
                StUPCTrack* positiveTrack;
                StUPCTrack* negativeTrack;
                if(tempTrack1->getCharge()>0){
                    positiveTrack = tempTrack1;
                    negativeTrack = tempTrack2;
                } else{
                    positiveTrack = tempTrack2;
                    negativeTrack = tempTrack1;
                }
                double deltaT = DeltaT0(positiveTrack, negativeTrack, massmap[PDGpositive], massmap[PDGnegative]);
                //fixing the +-1ns peaks
                if(fabs(deltaT-1)<fabs(deltaT)){
                    deltaT -= 1;
                } else if(fabs(deltaT+1)<fabs(deltaT)){
                    deltaT += 1;
                }
                double sigmaT = sigmaTmap[PDGmain];
                double sigmaParticle1, sigmaParticle2;
                switch(PDGpositive){
                case piplusPDGid:
                    sigmaParticle1 = positiveTrack->getNSigmasTPCPion();
                    break;
                case KplusPDGid:
                    sigmaParticle1 = positiveTrack->getNSigmasTPCKaon();
                    break;
                case pplusPDGid:
                    sigmaParticle1 = positiveTrack->getNSigmasTPCProton();
                    break;
                default:
                    sigmaParticle1 = 0;
                    break;
                }
                switch(PDGnegative){
                case piminusPDGid:
                    sigmaParticle2 = negativeTrack->getNSigmasTPCPion();
                    break;
                case KminusPDGid:
                    sigmaParticle2 = negativeTrack->getNSigmasTPCKaon();
                    break;
                case pminusPDGid:
                    sigmaParticle2 = negativeTrack->getNSigmasTPCProton();
                    break;
                default:
                    sigmaParticle2 = 0;
                    break;
                }
                //chi2 histogram
                double Chi2 = pow(deltaT/sigmaT, 2)+pow(sigmaParticle1, 2)+pow(sigmaParticle2, 2);
                Chi2All.Fill(Chi2);
            }
        }

        //loop for all tracks
        for(int i = 0; i<tempUPCpointer->getNumberOfTracks(); i++){
            //check if track okay
            StUPCTrack* tempTrack = tempUPCpointer->getTrack(i);
            if(!tempTrack->getFlag(StUPCTrack::kTof)){
                continue;
            }
            if(tempTrack->getTofPathLength()<=0){
                continue;
            }
            if(tempTrack->getTofTime()<=0){
                continue;
            }
            if(fabs(tempTrack->getEta())>0.9||tempTrack->getPt()<0.2){
                continue;
            }
            if(tempTrack->getNhitsFit()<=20||tempTrack->getNhitsDEdx()<15){
                continue;
            }
            //fill all the nsigma histograms
            int addDependingOnSign = tempTrack->getCharge()>0 ? 0 : 1;
            int MCparticleId = 4;
            //fill particle id in a way that will make histograms fillable
            for(size_t MC_particle_id = 0; MC_particle_id<tempUPCpointer->getNumberOfMCParticles(); MC_particle_id++){
                if((tempTrack->getIdTruth()-1)==MC_particle_id){
                    switch(abs(tempUPCpointer->getMCParticle(MC_particle_id)->GetPdgCode())){
                    case 11://electron
                        MCparticleId = 0;
                        break;
                    case piplusPDGid:
                        MCparticleId = 1;
                        break;
                    case KplusPDGid:
                        MCparticleId = 2;
                        break;
                    case pplusPDGid:
                        MCparticleId = 3;
                        break;
                    default:
                        MCparticleId = 4;
                        if(tempTrack->getCharge()>0){
                            UnidentifiedPositive.Fill(tempUPCpointer->getMCParticle(MC_particle_id)->GetName(), 1.0);
                        } else{
                            UnidentifiedNegative.Fill(tempUPCpointer->getMCParticle(MC_particle_id)->GetName(), 1.0);
                        }
                        break;
                    }
                    //we identified the particle, we no longer need the loop
                    break;
                }
            }

            for(auto&& particle:{ 0,1,2,3 }){
                histograms[particle+4*MCparticleId+20*addDependingOnSign]->Fill(tempTrack->getNSigmasTPC(static_cast<StUPCTrack::Part>(particle)));
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
    FlowChart.GetXaxis()->SetLabelSize(0.08);
    FlowChart.Write();
    Chi2Signal.Write();
    Chi2All.Write();
    for(size_t i = 0; i<histograms.size(); i++){
        histograms[i]->Write();
    }
    UnidentifiedPositive.LabelsDeflate();
    UnidentifiedPositive.LabelsOption("a");
    UnidentifiedPositive.Write();
    UnidentifiedNegative.LabelsDeflate();
    UnidentifiedNegative.LabelsOption("a");
    UnidentifiedNegative.Write();

    outputFileHist->Close();

    return 0;
}

void PrintBigger(TParticle* input, std::string additional_stuff){
    Printf("TParticle: %-13s  p: %8f %8f %8f \tVertex: %8e %8e %8e \tProd. Vertex: %5d %5d \tDecay Vertex: %5d \tTOF tray:%5d \tTOF module:%5d%s",
        input->GetName(), input->Px(), input->Py(), input->Pz(), input->Vx(), input->Vy(), input->Vz(),
        input->GetFirstMother(), input->GetSecondMother(), input->GetFirstDaughter(),
        input->GetLastDaughter()/100, input->GetLastDaughter()%100, additional_stuff.c_str());
}