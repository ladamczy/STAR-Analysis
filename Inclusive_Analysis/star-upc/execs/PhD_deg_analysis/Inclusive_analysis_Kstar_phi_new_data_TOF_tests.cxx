//cpp headers

//ROOT headers
#include <ROOT/TThreadedObject.hxx>
#include <TTreeReader.h>
#include <ROOT/TTreeProcessorMT.hxx>
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "TStyle.h"
#include "THStack.h"

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
#include "UsefulThings.h"
#include "ProcessingInsideLoop.h"
#include "ProcessingOutsideLoop.h"
#include "MyStyles.h"

enum{
    kAll = 1, kCPT, kRP, kOneVertex, kTPCTOF,
    kTotQ, kMax
};
enum SIDE{ E = 0, East = 0, W = 1, West = 1, nSides };
enum PARTICLES{ Pion = 0, Kaon = 1, Proton = 2, nParticles };
const double particleMass[nParticles] = { 0.13957, 0.493677, 0.93827 }; // pion, kaon, proton in GeV /c^2 
enum BRANCH_ID{ EU, ED, WU, WD, nBranches };
enum RP_ID{ E1U, E1D, E2U, E2D, W1U, W1D, W2U, W2D, nRomanPots };
enum SUSPECTED_PARTICLES{ K0S, Lambda, Kstar, Phi };

double drawFit(TH1D* hist, string outfileFolder, string outfileName, double mint0 = -3., double maxt0 = 3., double* params = nullptr, string histName = "");

int main(int argc, char** argv){

    int nthreads = 1;
    std::string drawing_folder = "";
    switch(argc){
    case 5:
        nthreads = atoi(argv[3]);
        drawing_folder = argv[4];
        break;
    case 4:
        //if it is the path, then we have path to efficiency corrections
        //if it is the number, we have number of cores
        if(atoi(argv[3])==0){
            drawing_folder = argv[3];
        } else{
            nthreads = atoi(argv[3]);
        }
        break;
    default:
        printf("Invalid number of arguments. Proper argument usage:\n");
        printf("argv[1] - .root file or .list list\n");
        printf("argv[2] - output folder or file\n");
        printf("argv[3] - (optional) number of cores (default 1)\n");
        printf("argv[4] - (optional) folder to save plots in\n");
        return 1;
        break;
    }
    //summary
    printf("Output folder: %s\n", argv[2]);
    printf("Program is running on %d threads\n", nthreads);
    printf("Histogram plot folder: %s\n", drawing_folder.length()!=0 ? drawing_folder.c_str() : "Not used");

    ROOT::EnableThreadSafety();
    //actually i'm not sure if it's needed here
    ROOT::EnableImplicitMT(nthreads); //turn on multicore processing

    //preparing input & output
    TChain* upcChain = new TChain("mUPCTree");
    if(ConnectInput(argc, argv, upcChain)){
        std::cout<<"All files connected"<<std::endl;
    }
    const string& outputFolder = argv[2];

    //histograms
    ProcessingOutsideLoop outsideprocessing;
    string particleNicks[EXTENDED_PARTICLES::nParticlesExtended] = { "e", "pi", "K", "p" };
    //deltaT0
    //old data, new data, and total
    for(std::string&& histtype:{ "old", "new", "total" }){
        for(size_t i = 0; i<nParticlesExtended; i++){
            for(size_t j = 0; j<nParticlesExtended; j++){
                outsideprocessing.AddHistogram(TH1D(("deltaT0_"+histtype+"_"+particleNicks[i]+"_"+particleNicks[j]+"_Wide").c_str(), ";t_{0}^{+}-t_{0}^{-} [ns];pair count", 4000, -2000, 2000));
                outsideprocessing.AddHistogram(TH1D(("deltaT0_"+histtype+"_"+particleNicks[i]+"_"+particleNicks[j]).c_str(), ";t_{0}^{+}-t_{0}^{-} [ns];pair count", 4000, -20, 20));
                outsideprocessing.AddHistogram(TH1D(("deltaT0_"+histtype+"_"+particleNicks[i]+"_"+particleNicks[j]+"_Narrow").c_str(), ";t_{0}^{+}-t_{0}^{-} [ns];pair count", 200, -1, 1));
            }
        }
    }
    //a special one to prove that electrons are difficult to fit
    outsideprocessing.AddHistogram(TH1D("deltaT0p0Narrow", ";t_{0}^{+}-t_{0}^{-} [ns];pair count", 100, -5, 5));
    //choice
    TH1D temp("choice", ";choice", 16, 0, 16);
    // int pairchoice = ppiPair*8+pipiPair*4+KpiPair*2+KKPair;
    temp.GetXaxis()->SetBinLabel(1, "Nothing");
    temp.GetXaxis()->SetBinLabel(2, "KK");
    temp.GetXaxis()->SetBinLabel(3, "K#pi");
    temp.GetXaxis()->SetBinLabel(4, "K#pi+KK");
    temp.GetXaxis()->SetBinLabel(5, "#pi#pi");
    temp.GetXaxis()->SetBinLabel(6, "#pi#pi+KK");
    temp.GetXaxis()->SetBinLabel(7, "#pi#pi+Kpi");
    temp.GetXaxis()->SetBinLabel(8, "#pi#pi+Kpi+KK");
    temp.GetXaxis()->SetBinLabel(9, "p#pi");
    temp.GetXaxis()->SetBinLabel(11, "p#pi+K#pi");
    temp.GetXaxis()->SetBinLabel(15, "p#pi+KK+K#pi");
    temp.GetXaxis()->SetBinLabel(16, "All");
    outsideprocessing.AddHistogram(temp);
    //the rest
    outsideprocessing.AddHistogram(TH1D("MKKNarrowChoice", ";m_{K^{#pm}K^{#pm}} [GeV/c^{2}];Number of pairs", 100, 0.9, 1.7));
    outsideprocessing.AddHistogram(TH1D("MKpiNarrowChoice", ";m_{K^{#pm}#pi^{#mp}} [GeV/c^{2}];Number of pairs", 25, 0.8, 1.0));
    outsideprocessing.AddHistogram(TH1D("MppiNarrowChoice", ";m_{K^{#pm}#pi^{#mp}} [GeV/c^{2}];Number of pairs", 75, 0.9, 1.5));
    outsideprocessing.AddHistogram(TH1D("MpipiChoice", ";m_{#pi^{#pm}#pi^{#mp}} [GeV/c^{2}];Number of pairs", 100, 0., 2.));
    outsideprocessing.AddHistogram(TH1D("MpipiNoChoice", ";m_{#pi^{#pm}#pi^{#mp}} [GeV/c^{2}];Number of pairs", 100, 0., 2.));
    outsideprocessing.AddHistogram(TH1D("M2TOFpipiKK", ";m_{TOF}^{2} [GeV^{2}/c^{4}];Number of pairs", 100, -0.25, 0.75));
    outsideprocessing.AddHistogram(TH2D("chi2pipiKK", ";n#sigma_{#pi};n#sigma_{K}", 50, 0, 10, 50, 0, 10));

    int events_processed = 0;

    //processing
    //defining TreeProcessor
    ROOT::TTreeProcessorMT TreeProc(*upcChain, nthreads);

    //defining processing function
    auto myFunction = [&](TTreeReader& myReader){
        //getting values from TChain, in-loop histogram initialization
        TTreeReaderValue<StUPCEvent> StUPCEventInstance(myReader, "mUPCEvent");
        TTreeReaderValue<StRPEvent> StRPEventInstance(myReader, "mRPEvent");
        ProcessingInsideLoop insideprocessing;
        StUPCEvent* tempUPCpointer;
        StRPEvent* tempRPpointer;
        insideprocessing.GetLocalHistograms(&outsideprocessing);

        //helpful variables
        std::vector<StUPCTrack*> vector_Track_positive_old;
        std::vector<StUPCTrack*> vector_Track_negative_old;
        std::vector<StUPCTrack*> vector_Track_positive_new;
        std::vector<StUPCTrack*> vector_Track_negative_new;
        StUPCTrack* tempTrack;
        TLorentzVector positive_track;
        TLorentzVector negative_track;
        double mass;

        //actual loop
        while(myReader.Next()){
            //in a TTree, it *would* be constant, in TChain however not necessarily
            tempUPCpointer = StUPCEventInstance.Get();
            tempRPpointer = StRPEventInstance.Get();

            //cleaning the loop
            vector_Track_positive_old.clear();
            vector_Track_negative_old.clear();
            vector_Track_positive_new.clear();
            vector_Track_negative_new.clear();

            //cause I want to see what's going on
            if(events_processed%10000==0){
                std::cout<<"Processed "<<events_processed<<" events"<<endl;
            }
            events_processed++;

            //additional cuts that normally are used in only-RP-cuts examples
            //cause the data i used isn't properly filtered
            //at least 2 good tracks, either from new or old data
            int nOfGoodTracksOld = 0;
            int nOfGoodTracksNew = 0;
            for(int i = 0; i<tempUPCpointer->getNumberOfTracks(); i++){
                StUPCTrack* tmptrk = tempUPCpointer->getTrack(i);
                if(!tmptrk->getFlag(StUPCTrack::kTof)||std::abs(tmptrk->getEta())>0.9||tmptrk->getPt()<0.2){
                    continue;
                }
                if(tmptrk->getFlag(StUPCTrack::kPrimary)){
                    nOfGoodTracksOld++;
                } else if(tmptrk->getFlag(StUPCTrack::kV0)){
                    nOfGoodTracksNew++;
                }
            }
            if(nOfGoodTracksOld<2&&nOfGoodTracksNew<2){
                continue;
            }
            //exactly one vertex
            if(tempUPCpointer->getNumberOfVertices()!=1){
                continue;
            }

            //cuts & histogram filling
            //selecting tracks matching criteria:
            //TOF
            //pt & eta 
            //Nhits
            //kPrimary or kV0 flag
            for(int i = 0; i<tempUPCpointer->getNumberOfTracks(); i++){
                tempTrack = tempUPCpointer->getTrack(i);
                if(!tempTrack->getFlag(StUPCTrack::kTof)){
                    continue;
                }
                if(tempTrack->getPt()<=0.2 or abs(tempTrack->getEta())>=0.9){
                    continue;
                }
                if(tempTrack->getNhits()<=20){
                    continue;
                }
                //old vs new separation
                if(tempTrack->getFlag(StUPCTrack::kPrimary)){
                    if(tempTrack->getCharge()>0){
                        vector_Track_positive_old.push_back(tempTrack);
                    } else{
                        vector_Track_negative_old.push_back(tempTrack);
                    }
                } else if(tempTrack->getFlag(StUPCTrack::kV0)){
                    if(tempTrack->getCharge()>0){
                        vector_Track_positive_new.push_back(tempTrack);
                    } else{
                        vector_Track_negative_new.push_back(tempTrack);
                    }
                }
            }

            //loop through particles, old (primary)
            for(long unsigned int i = 0; i<vector_Track_positive_old.size(); i++){
                for(long unsigned int j = 0; j<vector_Track_negative_old.size(); j++){
                    //test if the pair contains good particles (in TOF sense)
                    if(vector_Track_positive_old[i]->getTofPathLength()<=0 or
                        vector_Track_positive_old[i]->getTofTime()<=0 or
                        vector_Track_negative_old[j]->getTofPathLength()<=0 or
                        vector_Track_negative_old[j]->getTofTime()<=0){
                        continue;
                    }

                    //filling the deltaT0 histograms
                    for(size_t pos = 0; pos<nParticlesExtended; pos++){
                        for(size_t neg = 0; neg<nParticlesExtended; neg++){
                            insideprocessing.Fill(("deltaT0_old_"+particleNicks[pos]+"_"+particleNicks[neg]+"_Wide").c_str(), DeltaT0(vector_Track_positive_old[i], vector_Track_negative_old[j], particleMassExtended[pos], particleMassExtended[neg]));
                            insideprocessing.Fill(("deltaT0_old_"+particleNicks[pos]+"_"+particleNicks[neg]).c_str(), DeltaT0(vector_Track_positive_old[i], vector_Track_negative_old[j], particleMassExtended[pos], particleMassExtended[neg]));
                            insideprocessing.Fill(("deltaT0_old_"+particleNicks[pos]+"_"+particleNicks[neg]+"_Narrow").c_str(), DeltaT0(vector_Track_positive_old[i], vector_Track_negative_old[j], particleMassExtended[pos], particleMassExtended[neg]));
                            insideprocessing.Fill(("deltaT0_total_"+particleNicks[pos]+"_"+particleNicks[neg]+"_Wide").c_str(), DeltaT0(vector_Track_positive_old[i], vector_Track_negative_old[j], particleMassExtended[pos], particleMassExtended[neg]));
                            insideprocessing.Fill(("deltaT0_total_"+particleNicks[pos]+"_"+particleNicks[neg]).c_str(), DeltaT0(vector_Track_positive_old[i], vector_Track_negative_old[j], particleMassExtended[pos], particleMassExtended[neg]));
                            insideprocessing.Fill(("deltaT0_total_"+particleNicks[pos]+"_"+particleNicks[neg]+"_Narrow").c_str(), DeltaT0(vector_Track_positive_old[i], vector_Track_negative_old[j], particleMassExtended[pos], particleMassExtended[neg]));
                        }
                    }

                    //test histogram for peaks close to 0
                    insideprocessing.Fill("deltaT0p0Narrow", DeltaT0(vector_Track_positive_old[i], vector_Track_negative_old[j], particleMassExtended[ExtProton], 0.));
                    double t0cutoff = 0.6;
                    bool pipiPair = abs(DeltaT0(vector_Track_positive_old[i], vector_Track_negative_old[j], particleMass[Pion], particleMass[Pion]))<t0cutoff;
                    bool KpiPair = abs(DeltaT0(vector_Track_positive_old[i], vector_Track_negative_old[j], particleMass[Pion], particleMass[Kaon]))<t0cutoff or abs(DeltaT0(vector_Track_positive_old[i], vector_Track_negative_old[j], particleMass[Kaon], particleMass[Pion]))<t0cutoff;
                    bool KKPair = abs(DeltaT0(vector_Track_positive_old[i], vector_Track_negative_old[j], particleMass[Kaon], particleMass[Kaon]))<t0cutoff;
                    bool ppiPair = abs(DeltaT0(vector_Track_positive_old[i], vector_Track_negative_old[j], particleMass[Pion], particleMass[Proton]))<t0cutoff or abs(DeltaT0(vector_Track_positive_old[i], vector_Track_negative_old[j], particleMass[Proton], particleMass[Pion]))<t0cutoff;
                    int pairchoice = ppiPair*8+pipiPair*4+KpiPair*2+KKPair;
                    insideprocessing.Fill("choice", pairchoice);

                    switch(pairchoice){
                        //check for every pair
                    case 5:
                    {
                        insideprocessing.Fill("M2TOFpipiKK", M2TOF(vector_Track_positive_old[i], vector_Track_negative_old[j]));
                        double chi2pion = pow(vector_Track_positive_old[i]->getNSigmasTPCPion(), 2)+pow(vector_Track_negative_old[j]->getNSigmasTPCPion(), 2);
                        double chi2kaon = pow(vector_Track_positive_old[i]->getNSigmasTPCKaon(), 2)+pow(vector_Track_negative_old[j]->getNSigmasTPCKaon(), 2);
                        insideprocessing.Fill("chi2pipiKK", sqrt(chi2pion), sqrt(chi2kaon));
                        //mutually exclusive
                    }
                    case 8:
                        if(abs(DeltaT0(vector_Track_positive_old[i], vector_Track_negative_old[j], particleMass[Pion], particleMass[Proton]))<abs(DeltaT0(vector_Track_positive_old[i], vector_Track_negative_old[j], particleMass[Proton], particleMass[Pion]))){
                            vector_Track_positive_old[i]->getLorentzVector(positive_track, particleMass[Pion]);
                            vector_Track_negative_old[j]->getLorentzVector(negative_track, particleMass[Proton]);
                        } else{
                            vector_Track_positive_old[i]->getLorentzVector(positive_track, particleMass[Proton]);
                            vector_Track_negative_old[j]->getLorentzVector(negative_track, particleMass[Pion]);
                        }
                        mass = (positive_track+negative_track).M();
                        insideprocessing.Fill("MppiNarrowChoice", mass);
                        break;
                    case 4:
                        vector_Track_positive_old[i]->getLorentzVector(positive_track, particleMass[Pion]);
                        vector_Track_negative_old[j]->getLorentzVector(negative_track, particleMass[Pion]);
                        mass = (positive_track+negative_track).M();
                        insideprocessing.Fill("MpipiChoice", mass);
                        break;
                    case 2:
                        if(abs(DeltaT0(vector_Track_positive_old[i], vector_Track_negative_old[j], particleMass[Pion], particleMass[Kaon]))<abs(DeltaT0(vector_Track_positive_old[i], vector_Track_negative_old[j], particleMass[Kaon], particleMass[Pion]))){
                            vector_Track_positive_old[i]->getLorentzVector(positive_track, particleMass[Pion]);
                            vector_Track_negative_old[j]->getLorentzVector(negative_track, particleMass[Kaon]);
                        } else{
                            vector_Track_positive_old[i]->getLorentzVector(positive_track, particleMass[Kaon]);
                            vector_Track_negative_old[j]->getLorentzVector(negative_track, particleMass[Pion]);
                        }
                        mass = (positive_track+negative_track).M();
                        insideprocessing.Fill("MKpiNarrowChoice", mass);
                        break;
                    case 1:
                        vector_Track_positive_old[i]->getLorentzVector(positive_track, particleMass[Kaon]);
                        vector_Track_negative_old[j]->getLorentzVector(negative_track, particleMass[Kaon]);
                        mass = (positive_track+negative_track).M();
                        insideprocessing.Fill("MKKNarrowChoice", mass);
                        break;
                    case 0:
                        vector_Track_positive_old[i]->getLorentzVector(positive_track, particleMass[Pion]);
                        vector_Track_negative_old[j]->getLorentzVector(negative_track, particleMass[Pion]);
                        mass = (positive_track+negative_track).M();
                        insideprocessing.Fill("MpipiNoChoice", mass);
                        break;
                    default:
                        break;
                    }
                }
            }
            //loop through new particles (new, kV0)
            for(long unsigned int i = 0; i<vector_Track_positive_new.size(); i++){
                for(long unsigned int j = 0; j<vector_Track_negative_new.size(); j++){
                    //test if the pair contains good particles (in TOF sense)
                    if(vector_Track_positive_new[i]->getTofPathLength()<=0 or
                        vector_Track_positive_new[i]->getTofTime()<=0 or
                        vector_Track_negative_new[j]->getTofPathLength()<=0 or
                        vector_Track_negative_new[j]->getTofTime()<=0){
                        continue;
                    }

                    //filling the deltaT0 histograms
                    for(size_t pos = 0; pos<nParticlesExtended; pos++){
                        for(size_t neg = 0; neg<nParticlesExtended; neg++){
                            insideprocessing.Fill(("deltaT0_new_"+particleNicks[pos]+"_"+particleNicks[neg]+"_Wide").c_str(), DeltaT0(vector_Track_positive_new[i], vector_Track_negative_new[j], particleMassExtended[pos], particleMassExtended[neg]));
                            insideprocessing.Fill(("deltaT0_new_"+particleNicks[pos]+"_"+particleNicks[neg]).c_str(), DeltaT0(vector_Track_positive_new[i], vector_Track_negative_new[j], particleMassExtended[pos], particleMassExtended[neg]));
                            insideprocessing.Fill(("deltaT0_new_"+particleNicks[pos]+"_"+particleNicks[neg]+"_Narrow").c_str(), DeltaT0(vector_Track_positive_new[i], vector_Track_negative_new[j], particleMassExtended[pos], particleMassExtended[neg]));
                            insideprocessing.Fill(("deltaT0_total_"+particleNicks[pos]+"_"+particleNicks[neg]+"_Wide").c_str(), DeltaT0(vector_Track_positive_new[i], vector_Track_negative_new[j], particleMassExtended[pos], particleMassExtended[neg]));
                            insideprocessing.Fill(("deltaT0_total_"+particleNicks[pos]+"_"+particleNicks[neg]).c_str(), DeltaT0(vector_Track_positive_new[i], vector_Track_negative_new[j], particleMassExtended[pos], particleMassExtended[neg]));
                            insideprocessing.Fill(("deltaT0_total_"+particleNicks[pos]+"_"+particleNicks[neg]+"_Narrow").c_str(), DeltaT0(vector_Track_positive_new[i], vector_Track_negative_new[j], particleMassExtended[pos], particleMassExtended[neg]));
                        }
                    }
                }
            }



            //test of normalising to earliest occurence
            //didn't work
            //test of checking mass as calculated by TOF
            //kinda worked for KK pair
            //test of mass 


            //lambda finish
        }
        return 0;
    };

    TreeProc.Process(myFunction);

    outsideprocessing.Merge();

    //setting up output file name and folder and stuff
    string path = string(argv[0]);
    string outfileFullPath, outfileFolder, outfileName;
    if(outputFolder.find(".root")!=std::string::npos){
        outfileFullPath = outputFolder;
    } else{
        outfileFullPath = outputFolder+"AnaOutput_"+path.substr(path.find_last_of("/\\")+1)+".root";
    }
    outfileFolder = outfileFullPath.substr(0, outfileFullPath.find_last_of("/\\")+1);
    outfileName = outfileFullPath.substr(outfileFullPath.find_last_of("/\\")+1);

    //fitting the gauss+pol2 and drawing histograms
    std::vector<string> sigmaNames;
    std::vector<double> sigmaValues;
    string temptitle, tempname;
    double tempsigma;
    double lower_fitting_edge = -0.6;
    double upper_fitting_edge = 0.6;
    for(size_t i = 0; i<nParticlesExtended; i++){
        for(size_t j = 0; j<nParticlesExtended; j++){
            temptitle = "#Delta t_{0}: "+((i==1) ? "#pi" : particleNicks[i])+"^{+}"+((j==1) ? "#pi" : particleNicks[j])+"^{-}";
            tempname = "deltaT0_total_"+particleNicks[i]+"_"+particleNicks[j]+"_Narrow";
            //getting the sigma
            tempsigma = drawFit(outsideprocessing.GetPointerAfterMerge1D(tempname.c_str()).get(), outfileFolder, outfileName, lower_fitting_edge, upper_fitting_edge, nullptr, temptitle);
            sigmaNames.push_back(particleNicks[i]+"_"+particleNicks[j]);
            sigmaValues.push_back(tempsigma);
        }
    }
    drawFit(outsideprocessing.GetPointerAfterMerge1D("deltaT0p0Narrow").get(), outfileFolder, outfileName, lower_fitting_edge, upper_fitting_edge, nullptr, "#Delta t_{0}: p\"#gamma\"");
    //fitting the same, but not drawing - only taking data
    std::vector<double> sigmaValuesOld, sigmaValuesNew;
    for(size_t i = 0; i<nParticlesExtended; i++){
        for(size_t j = 0; j<nParticlesExtended; j++){
            //old
            tempname = "deltaT0_old_"+particleNicks[i]+"_"+particleNicks[j]+"_Narrow";
            //getting the sigma
            tempsigma = drawFit(outsideprocessing.GetPointerAfterMerge1D(tempname.c_str()).get(), "", "", lower_fitting_edge, upper_fitting_edge, nullptr, "");
            sigmaValuesOld.push_back(tempsigma);

            //new
            tempname = "deltaT0_new_"+particleNicks[i]+"_"+particleNicks[j]+"_Narrow";
            //getting the sigma
            tempsigma = drawFit(outsideprocessing.GetPointerAfterMerge1D(tempname.c_str()).get(), "", "", lower_fitting_edge, upper_fitting_edge, nullptr, "");
            sigmaValuesNew.push_back(tempsigma);
        }
    }
    //making bar histogram
    THStack histogramSigmaAll("histogramSigmaAll", "");
    TH1D histogramSigmaOld("histogramSigmaOld", "", 16, 0, 16);
    TH1D histogramSigmaNew("histogramSigmaNew", "", 16, 0, 16);
    TH1D histogramSigmaTotal("histogramSigmaTotal", "", 16, 0, 16);
    histogramSigmaOld.SetFillColor(kRed);
    histogramSigmaNew.SetFillColor(kGreen);
    histogramSigmaTotal.SetFillColor(kBlue);
    for(size_t i = 0; i<nParticlesExtended*nParticlesExtended; i++){
        //labels
        histogramSigmaOld.GetXaxis()->SetBinLabel(i+1, sigmaNames[i].c_str());
        histogramSigmaNew.GetXaxis()->SetBinLabel(i+1, sigmaNames[i].c_str());
        histogramSigmaTotal.GetXaxis()->SetBinLabel(i+1, sigmaNames[i].c_str());
        //values
        histogramSigmaOld.SetBinContent(i+1, sigmaValuesOld[i]);
        histogramSigmaNew.SetBinContent(i+1, sigmaValuesNew[i]);
        histogramSigmaTotal.SetBinContent(i+1, sigmaValuesOld[i]);
    }
    histogramSigmaAll.Add(&histogramSigmaOld);
    histogramSigmaAll.Add(&histogramSigmaNew);
    histogramSigmaAll.Add(&histogramSigmaTotal);
    histogramSigmaAll.SetDrawOption("nostackb");

    //actully making the output file
    std::cout<<"Created output file "<<outfileFullPath<<endl;
    TFile* outputFileHist = TFile::Open(outfileFullPath.c_str(), "recreate");
    outsideprocessing.SaveToFile(outputFileHist);
    histogramSigmaAll.Write();
    outputFileHist->Close();

    //writing calculated sigma values to a txt file
    string sigmaOutputFullPath = outfileFullPath.substr(0, outfileFullPath.find_last_of("."))+"_sigmaValues.txt";
    string tempString;
    TFile* sigmaOutputFile = TFile::Open((sigmaOutputFullPath+"?filetype=raw").c_str(), "recreate");
    for(size_t i = 0; i<sigmaNames.size(); i++){
        printf("%s:\t%lf\n", sigmaNames[i].c_str(), sigmaValues[i]);
        tempString = sigmaNames[i]+"\t\t"+sigmaValues[i]+"\n";
        sigmaOutputFile->WriteBuffer(tempString.c_str(), tempString.length());
        sigmaOutputFile->Flush();
    }
    sigmaOutputFile->Close();

    return 0;
}

double drawFit(TH1D* hist, string outfileFolder, string outfileName, double mint0, double maxt0, double* params, string histName){
    TCanvas* result = new TCanvas("result", "result", 1800, 1600);
    TF1* GfitK = new TF1("GfitK", "gausn(0) + pol2(3)");
    GfitK->SetRange(mint0, maxt0);
    GfitK->SetParNames("Constant", "Mean", "Sigma", "c", "b", "a");
    if(params!=nullptr){
        GfitK->SetParameters(params);
    } else{
        GfitK->SetParameters(1000, 0.0, 0.2, 0, 0, 0);
    }
    GfitK->SetParLimits(2, 0., 0.5);
    hist->SetMinimum(0);
    hist->SetMarkerStyle(kFullCircle);
    hist->Fit(GfitK, "R0");
    if(histName.length()!=0){
        hist->SetTitle(histName.c_str());
    }
    hist->Draw("E");
    Double_t paramsK[6];
    GfitK->GetParameters(paramsK);
    GfitK->SetNpx(1000);
    GfitK->Draw("CSAME");
    result->UseCurrentStyle();

    if(outfileFolder.size()!=0&&outfileName.size()!=0){
        TStyle mystyle = MyStyles::Hist2DQuarterSize(true);
        mystyle.cd();
        gROOT->ForceStyle();
        result->UseCurrentStyle();
        string outfileFullPath = outfileFolder+outfileName;
        string output = outfileFullPath.insert(outfileFullPath.find_last_of("."), "_"+string(hist->GetName())).substr(0, outfileFullPath.find_last_of("."))+".pdf";
        result->SaveAs(output.c_str());
    }
    gStyle->SetOptStat(1);
    return paramsK[2];
}