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
#include "TH3D.h"
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
    ROOT::TThreadedObject<TH1D> metricbetween = ROOT::TThreadedObject<TH1D>("metricbetween", "Distance between simulation and event in #phi-#eta space;;events", 18, 0, TMath::Pi());
    ROOT::TThreadedObject<TH1D> metricbetweenprecise = ROOT::TThreadedObject<TH1D>("metricbetweenprecise", "Distance between simulation and event in #phi-#eta space;;events", 20, 0, 0.1);
    ROOT::TThreadedObject<TH1D> nongeneratorparticles = ROOT::TThreadedObject<TH1D>("nongeneratorparticles", "Particles not created by generator;;particles", 3, 0, 3);
    ROOT::TThreadedObject<TH1D> xvertexdifference = ROOT::TThreadedObject<TH1D>("xvertexdifference", "Difference in X position between vertex and its reconstruction;#Delta x [cm];events", 40, -20, 20);
    ROOT::TThreadedObject<TH1D> yvertexdifference = ROOT::TThreadedObject<TH1D>("yvertexdifference", "Difference in Y position between vertex and its reconstruction;#Delta y [cm];events", 40, -20, 20);
    ROOT::TThreadedObject<TH1D> zvertexdifference = ROOT::TThreadedObject<TH1D>("zvertexdifference", "Difference in Z position between vertex and its reconstruction;#Delta z [cm];events", 40, -20, 20);
    ROOT::TThreadedObject<TH1D> vertexdifference = ROOT::TThreadedObject<TH1D>("vertexdifference", "Difference in position between vertex and its reconstruction;#Delta d [cm];events", 20, 0, 20);
    ROOT::TThreadedObject<TH2D> massvsxvertexdifference = ROOT::TThreadedObject<TH2D>("massvsxvertexdifference", "Reconstructed kaon mass vs difference in X position between vertex and its reconstruction;K^{0}_{S} mass [GeV];#Delta x [cm]", 100, 0.42, 0.56, 40, -20, 20);
    ROOT::TThreadedObject<TH2D> massvsyvertexdifference = ROOT::TThreadedObject<TH2D>("massvsyvertexdifference", "Reconstructed kaon mass vs difference in Y position between vertex and its reconstruction;K^{0}_{S} mass [GeV];#Delta y [cm]", 100, 0.42, 0.56, 40, -20, 20);
    ROOT::TThreadedObject<TH2D> massvszvertexdifference = ROOT::TThreadedObject<TH2D>("massvszvertexdifference", "Reconstructed kaon mass vs difference in Z position between vertex and its reconstruction;K^{0}_{S} mass [GeV];#Delta z [cm]", 100, 0.42, 0.56, 40, -20, 20);
    ROOT::TThreadedObject<TH2D> massvsvertexdifference = ROOT::TThreadedObject<TH2D>("massvsvertexdifference", "Reconstructed kaon mass vs difference in position between vertex and its reconstruction;K^{0}_{S} mass [GeV];#Delta d [cm]", 100, 0.42, 0.56, 20, 0, 20);
    ROOT::TThreadedObject<TH1D> invmass = ROOT::TThreadedObject<TH1D>("invmass", "K^{0}_{S} invariant mass;m [GeV];entries", 100, 0.42, 0.56);
    //dE/dx type, particle  name, value
    ROOT::TThreadedObject<TH3D> dEdxsigma = ROOT::TThreadedObject<TH3D>("dEdxsigma", "dE/dx for different particles", 4, 0, 4, 4, 0, 4, 40, -5, 5);
    ROOT::TThreadedObject<TH2D> massvsadditionalparticles = ROOT::TThreadedObject<TH2D>("massvsadditionalparticles", "Reconstructed kaon mass vs number of additional particles beside 4 pions;K^{0}_{S} mass [GeV];tracks", 100, 0.42, 0.56, 30, 0, 30);    

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
        std::shared_ptr<TH1D> metricbetweenLocal;
        std::shared_ptr<TH1D> metricbetweenpreciseLocal;
        std::shared_ptr<TH1D> nongeneratorparticlesLocal;
        std::shared_ptr<TH1D> xvertexdifferenceLocal;
        std::shared_ptr<TH1D> yvertexdifferenceLocal;
        std::shared_ptr<TH1D> zvertexdifferenceLocal;
        std::shared_ptr<TH1D> vertexdifferenceLocal;
        std::shared_ptr<TH2D> massvsxvertexdifferenceLocal;
        std::shared_ptr<TH2D> massvsyvertexdifferenceLocal;
        std::shared_ptr<TH2D> massvszvertexdifferenceLocal;
        std::shared_ptr<TH2D> massvsvertexdifferenceLocal;
        std::shared_ptr<TH1D> invmassLocal;
        std::shared_ptr<TH3D> dEdxsigmaLocal;
        std::shared_ptr<TH2D> massvsadditionalparticlesLocal;

        metricbetweenLocal = metricbetween.Get();
        metricbetweenpreciseLocal = metricbetweenprecise.Get();
        nongeneratorparticlesLocal = nongeneratorparticles.Get();
        xvertexdifferenceLocal = xvertexdifference.Get();
        yvertexdifferenceLocal = yvertexdifference.Get();
        zvertexdifferenceLocal = zvertexdifference.Get();
        vertexdifferenceLocal = vertexdifference.Get();
        massvsxvertexdifferenceLocal = massvsxvertexdifference.Get();
        massvsyvertexdifferenceLocal = massvsyvertexdifference.Get();
        massvszvertexdifferenceLocal = massvszvertexdifference.Get();
        massvsvertexdifferenceLocal = massvsvertexdifference.Get();
        invmassLocal = invmass.Get();
        dEdxsigmaLocal = dEdxsigma.Get();
        massvsadditionalparticlesLocal = massvsadditionalparticles.Get();

        int filtered_entries = 0;
        //for changing branch address
        StUPCEvent* tempUPCpointer = StUPCEventInstance.Get();
        //unavailable
        // StRPEvent* tempRPpointer = StRPEventInstance.Get();

        int total_charged_MC_particles = 0;
        do{
            tempUPCpointer = StUPCEventInstance.Get();
            vector<TParticle*> true_particles;
            vector<StUPCTrack*> detected_particles;
            vector<double> metric_difference;
            vector<bool> was_used(tempUPCpointer->getNumberOfMCParticles(), false);

            //assigning and identifying
            for (int i = 0; i < tempUPCpointer->getNumberOfTracks(); i++){
                //processing
                int best_choice = -1;
                double best_metric = 100.0; //absurdly high
                for (int j = 0; j < tempUPCpointer->getNumberOfMCParticles(); j++){
                    //counting particles present and eliminating neutral ones
                    if (int(round(tempUPCpointer->getMCParticle(j)->GetPDG()->Charge()))!=0){
                        total_charged_MC_particles++;
                    }else{
                        was_used[j] = true;
                        continue;
                    }
                    //eliminating protons diffractively scattered
                    if(tempUPCpointer->getMCParticle(j)->GetPdgCode()==2212 && abs(tempUPCpointer->getMCParticle(j)->Eta())>5){
                        was_used[j] = true;
                        continue;
                    }
                    //check metric and match
                    double temp_metric = MetricCheck(tempUPCpointer->getMCParticle(j), tempUPCpointer->getTrack(i));
                    if(was_used[j] || temp_metric<0){
                        continue;
                    }else if(temp_metric<best_metric){
                        best_choice = j;
                        best_metric = temp_metric;
                        was_used[j] = true;
                    }
                }
                if(best_metric<100.0 && best_choice>=0){
                    true_particles.push_back(tempUPCpointer->getMCParticle(best_choice));
                    detected_particles.push_back(tempUPCpointer->getTrack(i));
                    metric_difference.push_back(best_metric);
                }
            }

            //things with MC particles
            //test of non-generator particles
            int piplusomitted = 0;
            int piminusomitted = 0;
            TLorentzVector vertexposition;
            bool hasvertexposition = false;
            for (int j = 0; j < tempUPCpointer->getNumberOfMCParticles(); j++){
                //omitting generator particles

                //eliminating protons diffractively scattered
                if(tempUPCpointer->getMCParticle(j)->GetPdgCode()==2212 && abs(tempUPCpointer->getMCParticle(j)->Eta())>5){
                    if(!hasvertexposition){
                        tempUPCpointer->getMCParticle(j)->ProductionVertex(vertexposition);
                    }
                    continue;
                }
                //K0S
                if(tempUPCpointer->getMCParticle(j)->GetPdgCode()==310){
                    continue;
                }
                //pi+- from K0S
                //because i cant identify mothers, there will be workaround
                //cause in 5th there is a particle created between pions and kaons
                if(tempUPCpointer->getMCParticle(j)->GetPdgCode()==211 && piplusomitted<2){
                    piplusomitted++;
                    continue;
                }
                if(tempUPCpointer->getMCParticle(j)->GetPdgCode()==-211 && piminusomitted<2){
                    piminusomitted++;
                    continue;
                }
                nongeneratorparticlesLocal->Fill(tempUPCpointer->getMCParticle(j)->GetName(), 1.);
            }

            //things to do with identified particles
            bool hasdetectedvertexposition = false;
            TVector3 detectedvertexposition;
            //test of pairing
            for(size_t i = 0; i<metric_difference.size(); i++){
                if(!hasdetectedvertexposition && detected_particles[i]->getFlag(StUPCTrack::kPrimary)){
                    hasdetectedvertexposition = true;
                    detectedvertexposition = {detected_particles[i]->getVertex()->getPosX(), detected_particles[i]->getVertex()->getPosY(), detected_particles[i]->getVertex()->getPosZ()};
                }
                metricbetweenLocal->Fill(metric_difference[i]);
                metricbetweenpreciseLocal->Fill(metric_difference[i]);
                for (int j = 0; j < 4; j++){
                    //enum Part{kElectron=0, kPion, kKaon, kProton};
                    //dE/dx type, particle  name, value
                    dEdxsigmaLocal->Fill(j, true_particles[i]->GetName(), detected_particles[i]->getNSigmasTPC(static_cast<StUPCTrack::Part>(j)), 1);
                }

                filtered_entries++;
            }

            double deltaX = vertexposition.X() - detectedvertexposition.X();
            double deltaY = vertexposition.Y() - detectedvertexposition.Y();
            double deltaZ = vertexposition.Z() - detectedvertexposition.Z();
            double delta = sqrt(deltaX*deltaX+deltaY*deltaY+deltaZ*deltaZ);
            xvertexdifferenceLocal->Fill(deltaX);
            yvertexdifferenceLocal->Fill(deltaY);
            zvertexdifferenceLocal->Fill(deltaZ);
            vertexdifferenceLocal->Fill(delta);

            //loop to get pions from K0Ss
            TLorentzVector pion1;
            TLorentzVector pion2;
            for(size_t i = 0; i<metric_difference.size(); i++){
                if(true_particles[i]->GetPdgCode()==211){
                    for(size_t j = 0; j<metric_difference.size(); j++){
                        if(true_particles[j]->GetPdgCode()==-211 && true_particles[i]->GetMother(0) == true_particles[j]->GetMother(0)){
                            detected_particles[i]->getLorentzVector(pion1, particleMass[Pion]);
                            detected_particles[j]->getLorentzVector(pion2, particleMass[Pion]);
                            double pionmass = (pion1+pion2).M();
                            massvsxvertexdifferenceLocal->Fill(pionmass, deltaX);
                            massvsyvertexdifferenceLocal->Fill(pionmass, deltaY);
                            massvszvertexdifferenceLocal->Fill(pionmass, deltaZ);
                            massvsvertexdifferenceLocal->Fill(pionmass, delta);
                            invmassLocal->Fill(pionmass);
                            // 2p, 2K0, 4pi means 8 there will always be >=8 particles
                            massvsadditionalparticlesLocal->Fill(pionmass, tempUPCpointer->getNumberOfMCParticles()-8);
                            continue;
                        }
                    }
                }
            }
            
            
            
            true_particles.clear();
            detected_particles.clear();
            metric_difference.clear();
        }while (myReader.Next());

        //deleting some no-name entries
        nongeneratorparticlesLocal->LabelsDeflate();
        dEdxsigmaLocal->LabelsDeflate("Y");

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
        cout<<"Matched "<<filtered_entries<<" out of "<<total_charged_MC_particles<<endl;

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
    std::shared_ptr<TH1D> metricbetweenFinal;
    std::shared_ptr<TH1D> metricbetweenpreciseFinal;
    std::shared_ptr<TH1D> nongeneratorparticlesFinal;
    std::shared_ptr<TH1D> xvertexdifferenceFinal;
    std::shared_ptr<TH1D> yvertexdifferenceFinal;
    std::shared_ptr<TH1D> zvertexdifferenceFinal;
    std::shared_ptr<TH1D> vertexdifferenceFinal;
    std::shared_ptr<TH2D> massvsxvertexdifferenceFinal;
    std::shared_ptr<TH2D> massvsyvertexdifferenceFinal;
    std::shared_ptr<TH2D> massvszvertexdifferenceFinal;
    std::shared_ptr<TH2D> massvsvertexdifferenceFinal;
    std::shared_ptr<TH1D> invmassFinal;
    std::shared_ptr<TH3D> dEdxsigmaFinal;
    std::shared_ptr<TH2D> massvsadditionalparticlesFinal;

    metricbetweenFinal = metricbetween.Merge();
    metricbetweenpreciseFinal = metricbetweenprecise.Merge();
    nongeneratorparticlesFinal = nongeneratorparticles.Merge();
    xvertexdifferenceFinal = xvertexdifference.Merge();
    yvertexdifferenceFinal = yvertexdifference.Merge();
    zvertexdifferenceFinal = zvertexdifference.Merge();
    vertexdifferenceFinal = vertexdifference.Merge();
    massvsxvertexdifferenceFinal = massvsxvertexdifference.Merge();
    massvsyvertexdifferenceFinal = massvsyvertexdifference.Merge();
    massvszvertexdifferenceFinal = massvszvertexdifference.Merge();
    massvsvertexdifferenceFinal = massvsvertexdifference.Merge();
    invmassFinal = invmass.Merge();
    dEdxsigmaFinal = dEdxsigma.Merge();
    massvsadditionalparticlesFinal = massvsadditionalparticles.Merge();

    metricbetweenFinal->SetMinimum(0);
    metricbetweenpreciseFinal->SetMinimum(0);
    nongeneratorparticlesFinal->SetMinimum(0);
    xvertexdifferenceFinal->SetMinimum(0);
    yvertexdifferenceFinal->SetMinimum(0);
    zvertexdifferenceFinal->SetMinimum(0);
    vertexdifferenceFinal->SetMinimum(0);
    massvsxvertexdifferenceFinal->SetMinimum(0);
    massvsyvertexdifferenceFinal->SetMinimum(0);
    massvszvertexdifferenceFinal->SetMinimum(0);
    massvsvertexdifferenceFinal->SetMinimum(0);
    invmassFinal->SetMinimum(0);
    // dEdxsigmaFinal->SetMinimum(0);
    massvsadditionalparticlesFinal->SetMinimum(0);

    //setting up a tree & output file
    string path = string(argv[0]);
    string outfileName = outputFolder + "AnaOutput_" + path.substr(path.find_last_of("/\\") + 1) + ".root";
    cout<<"Created output file "<<outfileName<<endl;
    TFile* outputFileHist = TFile::Open(outfileName.c_str(), "recreate");

    outputFileHist->cd();
    metricbetweenFinal->Write();
    metricbetweenpreciseFinal->Write();
    nongeneratorparticlesFinal->Write();
    xvertexdifferenceFinal->Write();
    yvertexdifferenceFinal->Write();
    zvertexdifferenceFinal->Write();
    vertexdifferenceFinal->Write();
    massvsxvertexdifferenceFinal->Write();
    massvsyvertexdifferenceFinal->Write();
    massvszvertexdifferenceFinal->Write();
    massvsvertexdifferenceFinal->Write();
    invmassFinal->Write();
    // dEdxsigmaFinal->Write();
    //dE/dx type, particle  name, value
    TH3D* tempcopy;
    TH2D* tempproj;
    string dEdxtypes[] = {"Electron", "Pion", "Kaon", "Proton"};
    for (int i = 1; i <= 4; i++){
        tempcopy = (TH3D*)dEdxsigmaFinal->Clone();
        tempcopy->GetXaxis()->SetRange(i, i);
        tempproj = (TH2D*)tempcopy->Project3D("zy");
        tempproj->SetTitle(string("dE/dx for different particles, type " + dEdxtypes[i-1]).c_str());
        tempproj->GetXaxis()->SetTitle("Particle type");
        tempproj->GetYaxis()->SetTitle("#frac{dE}{dx} #sigma value");
        tempproj->Write();
    }
    massvsadditionalparticlesFinal->Write();
    outputFileHist->Close();

    cout<<"Finished processing "<<endl;
    cout<<"Analyzed total "<<upcChain->GetEntries()<<" entries"<<endl;
    cout<<"Filtered total "<<filtered_entries<<" entries"<<endl;
    cout<<"Ending Analysis... GOOD BYE!"<<endl;
    return 0;
}//main