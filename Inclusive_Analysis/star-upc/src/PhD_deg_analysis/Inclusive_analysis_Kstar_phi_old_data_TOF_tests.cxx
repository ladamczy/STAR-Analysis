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
enum EXTENDED_PARTICLES{ ExtElectron = 0, ExtPion = 1, ExtKaon = 2, ExtProton = 3, nParticlesExtended };
const double particleMass[nParticles] = { 0.13957, 0.493677, 0.93827 }; // pion, kaon, proton in GeV /c^2 
const double particleMassExtended[nParticlesExtended] = { 0.000510999, 0.13957, 0.493677, 0.93827 }; // electron, pion, kaon, proton in GeV /c^2 
enum BRANCH_ID{ EU, ED, WU, WD, nBranches };
enum RP_ID{ E1U, E1D, E2U, E2D, W1U, W1D, W2U, W2D, nRomanPots };
enum SUSPECTED_PARTICLES{ K0S, Lambda, Kstar, Phi };

double drawFit(TH1D* hist, string outfileName, double mint0 = -3., double maxt0 = 3., double* params = nullptr, string histName = "");
int main(int argc, char** argv){

    int nthreads = 1;
    if(argc==4){
        nthreads = atoi(argv[3]);
    }

    cout<<"Program is running on "<<nthreads<<" threads"<<endl;
    ROOT::EnableThreadSafety();
    //actually i'm not sure if it's needed here
    // ROOT::EnableImplicitMT(nthreads); //turn on multicore processing

    //preparing input & output
    TChain* upcChain = new TChain("mUPCTree");
    if(ConnectInput(argc, argv, upcChain)){
        cout<<"All files connected"<<endl;
    }
    const string& outputFolder = argv[2];

    //histograms
    ProcessingOutsideLoop outsideprocessing;
    //deltaT0 - all the combinations (+narrow versions)
    string particleNicks[nParticlesExtended] = { "e", "pi", "K", "p" };
    for(size_t i = 0; i<nParticlesExtended; i++){
        for(size_t j = 0; j<nParticlesExtended; j++){
            outsideprocessing.AddHistogram(TH1D(("deltaT0"+particleNicks[i]+particleNicks[j]).c_str(), ";t_{0}^{+}-t_{0}^{-} [ns];pair count", 200, -20, 20));
        }
    }
    for(size_t i = 0; i<nParticlesExtended; i++){
        for(size_t j = 0; j<nParticlesExtended; j++){
            outsideprocessing.AddHistogram(TH1D(("deltaT0"+particleNicks[i]+particleNicks[j]+"Narrow").c_str(), ";t_{0}^{+}-t_{0}^{-} [ns];pair count", 100, -5, 5));
        }
    }
    outsideprocessing.AddHistogram(TH1D("deltaT0p0Narrow", ";t_{0}^{+}-t_{0}^{-} [ns];pair count", 100, -5, 5));
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
    outsideprocessing.AddHistogram(TH1D("MKKNarrowChoice", ";m_{K^{#pm}K^{#pm}} [GeV];Number of pairs", 100, 0.9, 1.7));
    outsideprocessing.AddHistogram(TH1D("MKpiNarrowChoice", ";m_{K^{#pm}#pi^{#mp}} [GeV];Number of pairs", 25, 0.8, 1.0));
    outsideprocessing.AddHistogram(TH1D("MppiNarrowChoice", ";m_{K^{#pm}#pi^{#mp}} [GeV];Number of pairs", 75, 0.9, 1.5));
    outsideprocessing.AddHistogram(TH1D("MpipiChoice", ";m_{#pi^{#pm}#pi^{#mp}} [GeV];Number of pairs", 100, 0., 2.));
    outsideprocessing.AddHistogram(TH1D("MpipiNoChoice", ";m_{#pi^{#pm}#pi^{#mp}} [GeV];Number of pairs", 100, 0., 2.));
    outsideprocessing.AddHistogram(TH1D("M2TOFpipiKK", ";m_{TOF}^{2} [GeV^{2}];Number of pairs", 100, -0.25, 0.75));
    outsideprocessing.AddHistogram(TH2D("chi2pipiKK", ";n#sigma_{#pi};n#sigma_{K}", 50, 0, 10, 50, 0, 10));

    //processing
    //defining TreeProcessor
    ROOT::TTreeProcessorMT TreeProc(*upcChain, nthreads);

    //other things
    int eventsProcessed = 0;
    double K0WindowLow = 0.48;
    double K0WindowHigh = 0.51;
    double Lambda0WindowLow = 1.05;
    double Lambda0WindowHigh = 1.15;

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
        std::vector<StUPCTrack*> vector_Track_positive;
        std::vector<StUPCTrack*> vector_Track_negative;
        StUPCTrack* tempTrack;
        TLorentzVector positive_track;
        TLorentzVector negative_track;
        double mass;
        std::vector<double> vector_t0_all_pions;
        double min_t0;

        //actual loop
        while(myReader.Next()){
            //in a TTree, it *would* be constant, in TChain however not necessarily
            tempUPCpointer = StUPCEventInstance.Get();
            tempRPpointer = StRPEventInstance.Get();

            //cleaning the loop
            vector_Track_positive.clear();
            vector_Track_negative.clear();

            //cause I want to see what's going on
            if(eventsProcessed%10000==0){
                cout<<"Processed "<<eventsProcessed<<" events"<<endl;
            }
            eventsProcessed++;

            //additional cuts that normally are used in only-RP-cuts examples
            //cause the data i used isn't properly filtered
            //at least 2 good tracks
            int nOfGoodTracks = 0;
            for(int i = 0; i<tempUPCpointer->getNumberOfTracks(); i++){
                StUPCTrack* tmptrk = tempUPCpointer->getTrack(i);
                if(tmptrk->getFlag(StUPCTrack::kTof)&&abs(tmptrk->getEta())<0.9&&tmptrk->getPt()>0.2){
                    nOfGoodTracks++;
                }
            }
            if(nOfGoodTracks<2){
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
            for(int i = 0; i<tempUPCpointer->getNumberOfTracks(); i++){
                tempTrack = tempUPCpointer->getTrack(i);
                if(!tempTrack->getFlag(StUPCTrack::kTof)){
                    continue;
                }
                if(!tempTrack->getFlag(StUPCTrack::kPrimary)){
                    continue;
                }
                if(tempTrack->getPt()<=0.2 or abs(tempTrack->getEta())>=0.9){
                    continue;
                }
                if(tempTrack->getNhits()<=20){
                    continue;
                }
                if(tempTrack->getCharge()>0){
                    vector_Track_positive.push_back(tempTrack);
                } else{
                    vector_Track_negative.push_back(tempTrack);
                }
            }

            //loop through particles
            for(long unsigned int i = 0; i<vector_Track_positive.size(); i++){
                for(long unsigned int j = 0; j<vector_Track_negative.size(); j++){
                    //test if the pair contains good particles (in TOF sense)
                    if(vector_Track_positive[i]->getTofPathLength()<=0 or
                        vector_Track_positive[i]->getTofTime()<=0 or
                        vector_Track_negative[j]->getTofPathLength()<=0 or
                        vector_Track_negative[j]->getTofTime()<=0){
                        continue;
                    }
                    //filling the histograms
                    for(size_t pos = 0; pos<nParticlesExtended; pos++){
                        for(size_t neg = 0; neg<nParticlesExtended; neg++){
                            insideprocessing.Fill(("deltaT0"+particleNicks[pos]+particleNicks[neg]).c_str(), DeltaT0(vector_Track_positive[i], vector_Track_negative[j], particleMassExtended[pos], particleMassExtended[neg]));
                            insideprocessing.Fill(("deltaT0"+particleNicks[pos]+particleNicks[neg]+"Narrow").c_str(), DeltaT0(vector_Track_positive[i], vector_Track_negative[j], particleMassExtended[pos], particleMassExtended[neg]));
                        }
                    }
                    //test histogram for peaks close to 0
                    insideprocessing.Fill("deltaT0p0Narrow", DeltaT0(vector_Track_positive[i], vector_Track_negative[j], particleMassExtended[ExtProton], 0.));
                    double t0cutoff = 0.6;
                    bool pipiPair = abs(DeltaT0(vector_Track_positive[i], vector_Track_negative[j], particleMass[Pion], particleMass[Pion]))<t0cutoff;
                    bool KpiPair = abs(DeltaT0(vector_Track_positive[i], vector_Track_negative[j], particleMass[Pion], particleMass[Kaon]))<t0cutoff or abs(DeltaT0(vector_Track_positive[i], vector_Track_negative[j], particleMass[Kaon], particleMass[Pion]))<t0cutoff;
                    bool KKPair = abs(DeltaT0(vector_Track_positive[i], vector_Track_negative[j], particleMass[Kaon], particleMass[Kaon]))<t0cutoff;
                    bool ppiPair = abs(DeltaT0(vector_Track_positive[i], vector_Track_negative[j], particleMass[Pion], particleMass[Proton]))<t0cutoff or abs(DeltaT0(vector_Track_positive[i], vector_Track_negative[j], particleMass[Proton], particleMass[Pion]))<t0cutoff;
                    int pairchoice = ppiPair*8+pipiPair*4+KpiPair*2+KKPair;
                    if(pairchoice==5){
                        insideprocessing.Fill("M2TOFpipiKK", M2TOF(vector_Track_positive[i], vector_Track_negative[j]));
                        double chi2pion = pow(vector_Track_positive[i]->getNSigmasTPCPion(), 2)+pow(vector_Track_negative[j]->getNSigmasTPCPion(), 2);
                        double chi2kaon = pow(vector_Track_positive[i]->getNSigmasTPCKaon(), 2)+pow(vector_Track_negative[j]->getNSigmasTPCKaon(), 2);
                        insideprocessing.Fill("chi2pipiKK", sqrt(chi2pion), sqrt(chi2kaon));
                    }
                    insideprocessing.Fill("choice", pairchoice);
                    if(pairchoice==8){
                        if(abs(DeltaT0(vector_Track_positive[i], vector_Track_negative[j], particleMass[Pion], particleMass[Proton]))<abs(DeltaT0(vector_Track_positive[i], vector_Track_negative[j], particleMass[Proton], particleMass[Pion]))){
                            vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                            vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Proton]);
                        } else{
                            vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Proton]);
                            vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Pion]);
                        }
                        mass = (positive_track+negative_track).M();
                        insideprocessing.Fill("MppiNarrowChoice", mass);
                    } else if(pairchoice==4){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                        vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Pion]);
                        mass = (positive_track+negative_track).M();
                        insideprocessing.Fill("MpipiChoice", mass);
                    } else if(pairchoice==2){
                        if(abs(DeltaT0(vector_Track_positive[i], vector_Track_negative[j], particleMass[Pion], particleMass[Kaon]))<abs(DeltaT0(vector_Track_positive[i], vector_Track_negative[j], particleMass[Kaon], particleMass[Pion]))){
                            vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                            vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Kaon]);
                        } else{
                            vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Kaon]);
                            vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Pion]);
                        }
                        mass = (positive_track+negative_track).M();
                        insideprocessing.Fill("MKpiNarrowChoice", mass);
                    } else if(pairchoice==1){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Kaon]);
                        vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Kaon]);
                        mass = (positive_track+negative_track).M();
                        insideprocessing.Fill("MKKNarrowChoice", mass);
                    } else if(pairchoice==0){
                        vector_Track_positive[i]->getLorentzVector(positive_track, particleMass[Pion]);
                        vector_Track_negative[j]->getLorentzVector(negative_track, particleMass[Pion]);
                        mass = (positive_track+negative_track).M();
                        insideprocessing.Fill("MpipiNoChoice", mass);
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

    //setting up output file name
    string path = string(argv[0]);
    string outfileName;
    if(outputFolder.find(".root")!=std::string::npos){
        outfileName = outputFolder;
    } else{
        outfileName = outputFolder+"AnaOutput_"+path.substr(path.find_last_of("/\\")+1)+".root";
    }
    //fitting the gauss+pol2
    std::vector<string> sigmaNames;
    std::vector<double> sigmaValues;
    string temptitle;
    double tempsigma;
    for(size_t i = 0; i<nParticlesExtended; i++){
        for(size_t j = 0; j<nParticlesExtended; j++){
            temptitle = "#Delta t_{0}: "+((i==1) ? "#pi" : particleNicks[i])+"^{+}"+((j==1) ? "#pi" : particleNicks[j])+"^{-}";
            tempsigma = drawFit(outsideprocessing.GetPointerAfterMerge1D(("deltaT0"+particleNicks[i]+particleNicks[j]+"Narrow").c_str()).get(), outfileName, -1.0, 1.0, nullptr, temptitle);
            sigmaNames.push_back(particleNicks[i]+"_"+particleNicks[j]);
            sigmaValues.push_back(tempsigma);
        }
    }
    drawFit(outsideprocessing.GetPointerAfterMerge1D("deltaT0p0Narrow").get(), outfileName, -1.0, 1.0, nullptr, "#Delta t_{0}: p\"#gamma\"");
    //actully making the output file
    cout<<"Created output file "<<outfileName<<endl;
    TFile* outputFileHist = TFile::Open(outfileName.c_str(), "recreate");
    outsideprocessing.SaveToFile(outputFileHist);
    outputFileHist->Close();

    //writing calculated sigma values to a txt file
    string sigmaOutput = outfileName.substr(0, outfileName.find_last_of("."))+"_sigmaValues.txt";
    string tempString;
    TFile* sigmaOutputFile = TFile::Open((sigmaOutput+"?filetype=raw").c_str(), "recreate");
    for(size_t i = 0; i<sigmaNames.size(); i++){
        printf("%s:\t%lf\n", sigmaNames[i].c_str(), sigmaValues[i]);
        tempString = sigmaNames[i]+"\t\t"+sigmaValues[i]+"\n";
        sigmaOutputFile->WriteBuffer(tempString.c_str(), tempString.length());
        sigmaOutputFile->Flush();
    }
    sigmaOutputFile->Close();

    return 0;
}

double drawFit(TH1D* hist, string outfileName, double mint0, double maxt0, double* params, string histName){
    MyStyles styleLibrary;
    TStyle mystyle = styleLibrary.Hist2DQuarterSize(true);
    mystyle.cd();
    gROOT->ForceStyle();
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

    string output = outfileName.insert(outfileName.find_last_of("."), "_"+string(hist->GetName())).substr(0, outfileName.find_last_of("."))+".pdf";
    result->SaveAs(output.c_str());
    gStyle->SetOptStat(1);
    return paramsK[2];
}