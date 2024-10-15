//cpp headers

//ROOT headers
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
#include "StUPCV0.h"
#include "BeamPosition.h"
#include "StPicoPhysicalHelix.h"

//my headers
#include "UsefulThings.h"
#include "ProcessingInsideLoop.h"
#include "ProcessingOutsideLoop.h"

enum{
    kAll = 1, kCPT, kRP, kOneVertex, kTPCTOF,
    kTotQ, kMax
};
enum SIDE{ E = 0, East = 0, W = 1, West = 1, nSides };
enum PARTICLES{ Pion = 0, Kaon = 1, Proton = 2, nParticles };
const double particleMass[nParticles] = { 0.13957, 0.493677, 0.93827 }; // pion, kaon, proton in GeV /c^2 
enum BRANCH_ID{ EU, ED, WU, WD, nBranches };
enum RP_ID{ E1U, E1D, E2U, E2D, W1U, W1D, W2U, W2D, nRomanPots };

int main(int argc, char **argv){

    int nthreads = 1;
    if(argc==4){
        nthreads = atoi(argv[3]);
    }

    cout<<"Program is running on "<<nthreads<<" threads"<<endl;
    ROOT::EnableThreadSafety();
    //actually i'm not sure if it's needed here
    // ROOT::EnableImplicitMT(nthreads); //turn on multicore processing

    //preparing input & output
    TChain *upcChain = new TChain("mUPCTree");
    if(ConnectInput(argc, argv, upcChain)){
        cout<<"All files connected"<<endl;
    }
    const string &outputFolder = argv[2];

    //useful constants
    vector<vector<double>> beamData = ReadFillPositionData("STAR-Analysis/share/Run7PolarizationWithPosition.csv");

    //histograms
    ProcessingOutsideLoop outsideprocessing;
    outsideprocessing.AddHistogram(TH1D("DataFigure3_4a", ";Number of vertices;Number of events", 10, 0., 10));
    outsideprocessing.AddHistogram(TH1D("DataFigure3_4b", ";Vertex z coordinate [cm];Number of events", 80, -200, 200));
    outsideprocessing.AddHistogram(TH1D("DataFigure3_5a", ";N^{fit}_{hits};number of tracks", 50, 0, 50));
    outsideprocessing.AddHistogram(TH1D("DataFigure3_5b", ";N^{fit}_{dE/dx};number of tracks", 50, 0, 50));
    outsideprocessing.AddHistogram(TH1D("DataFigure3_5c", ";DCA_{xy} [cm];number of tracks", 50, 0, 5));
    outsideprocessing.AddHistogram(TH1D("DataFigure3_5d", ";DCA_{z} [cm];number of tracks", 100, -5, 5));
    // outsideprocessing.AddHistogram(TH1D("DataFigure3_5e", ";d_{0} [cm];number of tracks", 50, 0, 5));
    outsideprocessing.AddHistogram(TH1D("DataFigure3_6ab", ";#eta;number of tracks", 100, -2, 2));
    outsideprocessing.AddHistogram(TH1D("DataFigure3_6c", ";p_{T} [GeV/c];number of tracks", 50, 0, 5));
    outsideprocessing.AddHistogram(TH1D("DataFigure3_6d", ";#Phi [rad];number of tracks", 100, -TMath::Pi(), TMath::Pi()));
    outsideprocessing.AddHistogram(TH1D("DataFigure3_18g", ";#xi_{1}#xi_{2};#frac{1}{N}#frac{dN}{d(#xi_{1}#xi_{2})}", 1000, 0, 0.04));
    outsideprocessing.AddHistogram(TH1D("DataFigure3_19a", ";n_{sel};#frac{1}{N}#frac{dN}{dn_{sel}}", 10, 2, 12));
    outsideprocessing.AddHistogram(TH1D("DataFigure3_20", ";p_{T} [GeV/c];#frac{1}{N}#frac{dN}{dp_{T}}", 50, 0, 5));
    outsideprocessing.AddHistogram(TH1D("DataFigure3_21", ";#eta;#frac{1}{N}#frac{dN}{d#eta}", 100, -2, 2));
    outsideprocessing.AddHistogram(TH2D("DataFigure3_60", ";q #times p [GeV/c];dE/dx [GeV/cm]", 500, -5, 5, 500, 0, 1e-4));
    outsideprocessing.AddHistogram(TH1D("DataFigure3_62a", ";n#sigma^{#pi};#frac{1}{N}#frac{1}{p_{T}}#frac{dN}{d(n#sigma^{#pi})}", 100, -50, 50));
    outsideprocessing.AddHistogram(TH1D("DataFigure3_62b", ";n#sigma^{K};#frac{1}{N}#frac{1}{p_{T}}#frac{dN}{d(n#sigma^{K})}", 100, -50, 50));
    outsideprocessing.AddHistogram(TH1D("DataFigure3_62c", ";n#sigma^{p};#frac{1}{N}#frac{1}{p_{T}}#frac{dN}{d(n#sigma^{p})}", 100, -50, 50));
    outsideprocessing.AddHistogram(TH1D("DataFigure3_62d", ";n#sigma^{e};#frac{1}{N}#frac{1}{p_{T}}#frac{dN}{d(n#sigma^{e})}", 100, -50, 50));

    //processing
    //defining TreeProcessor
    ROOT::TTreeProcessorMT TreeProc(*upcChain, nthreads);

    int eventsProcessed = 0;

    //defining processing function
    auto myFunction = [&](TTreeReader &myReader){
        //getting values from TChain, in-loop histogram initialization
        TTreeReaderValue<StUPCEvent> StUPCEventInstance(myReader, "mUPCEvent");
        TTreeReaderValue<StRPEvent> StRPEventInstance(myReader, "mRPEvent");
        ProcessingInsideLoop insideprocessing;
        StUPCEvent *tempUPCpointer;
        StRPEvent *tempRPpointer;
        insideprocessing.GetLocalHistograms(&outsideprocessing);

        //helpful variables
        // TVector3 vertexPrimary;
        // std::vector<double> tempBeamVector;
        // double beamValues[4];
        StUPCTrack *tempTrack;
        // StPicoPhysicalHelix *beamLine;
        // StPicoPhysicalHelix *trackLine;
        int numberOfSelectedTracks = 0;
        TVector3 tempVector;

        //actual loop
        while(myReader.Next()){
            //in a TTree, it *would* be constant, in TChain however not necessarily
            tempUPCpointer = StUPCEventInstance.Get();
            tempRPpointer = StRPEventInstance.Get();

            //cause I want to see what's going on
            if(eventsProcessed%10000==0){
                cout<<"Processed "<<eventsProcessed<<" events"<<endl;
            }
            eventsProcessed++;

            // vertexPrimary = { tempUPCpointer->getVertex(0)->getPosX(), tempUPCpointer->getVertex(0)->getPosY(), tempUPCpointer->getVertex(0)->getPosZ() };
            // tempBeamVector = FindPosition(tempUPCpointer->getFillNumber(), vertexPrimary.Z(), beamData[0], beamData[1], beamData[2], beamData[3], beamData[4], beamData[5], beamData[6], beamData[7], beamData[8]);
            // beamValues[0] = tempBeamVector[0];
            // beamValues[1] = tempBeamVector[1];
            // beamValues[2] = tempBeamVector[2];
            // beamValues[3] = tempBeamVector[3];

            //histogram filling
            //possibly make it better with StPicoPhysicslHelix? TODO
            insideprocessing.Fill("DataFigure3_4a", tempUPCpointer->getNumberOfVertices());
            if(tempUPCpointer->getNumberOfVertices()==1){
                insideprocessing.Fill("DataFigure3_4b", tempUPCpointer->getVertex(0)->getPosZ());
            } else{
                continue;
            }
            numberOfSelectedTracks = 0;
            for(int i = 0; i<tempUPCpointer->getNumberOfTracks(); i++){
                //used ToF, should I? TODO
                tempTrack = tempUPCpointer->getTrack(i);
                if(!tempTrack->getFlag(StUPCTrack::kTof)){
                    continue;
                }
                //3.5
                insideprocessing.Fill("DataFigure3_5a", tempTrack->getNhitsFit());
                insideprocessing.Fill("DataFigure3_5b", tempTrack->getNhitsDEdx());
                //dca to the vertex
                insideprocessing.Fill("DataFigure3_5c", tempTrack->getDcaXY());
                insideprocessing.Fill("DataFigure3_5d", tempTrack->getDcaZ());
                //dca to the beamline - TODO
                // insideprocessing.Fill("DataFigure3_5e", 0);
                //3.6
                insideprocessing.Fill("DataFigure3_6ab", tempTrack->getEta());
                insideprocessing.Fill("DataFigure3_6c", tempTrack->getPt());
                insideprocessing.Fill("DataFigure3_6d", tempTrack->getPhi());
                if(tempTrack->getPt()>0.2&&abs(tempTrack->getEta())<0.7&&tempTrack->getNhitsFit()>20){
                    numberOfSelectedTracks++;
                }
            }
            //histograms of normalised events
            if(numberOfSelectedTracks>=2){
                //3.18g
                insideprocessing.Fill("DataFigure3_18g", tempRPpointer->getTrack(0)->xi(beamMomentum)*tempRPpointer->getTrack(1)->xi(beamMomentum));
                //3.19a
                insideprocessing.Fill("DataFigure3_19a", numberOfSelectedTracks);
                //tracks
                for(int i = 0; i<tempUPCpointer->getNumberOfTracks(); i++){
                    //used ToF, should I? TODO
                    tempTrack = tempUPCpointer->getTrack(i);
                    if(!tempTrack->getFlag(StUPCTrack::kTof)){
                        continue;
                    }
                    //3.20
                    insideprocessing.Fill("DataFigure3_20", tempTrack->getPt());
                    //3.21
                    insideprocessing.Fill("DataFigure3_21", tempTrack->getEta());
                    //3.60
                    tempTrack->getMomentum(tempVector);
                    insideprocessing.Fill("DataFigure3_60", tempVector.Mag()*tempTrack->getCharge(), tempTrack->getDEdxSignal());
                    //3.62
                    insideprocessing.Fill("DataFigure3_62a", tempTrack->getNSigmasTPCPion(), 1/tempTrack->getPt());
                    insideprocessing.Fill("DataFigure3_62b", tempTrack->getNSigmasTPCKaon(), 1/tempTrack->getPt());
                    insideprocessing.Fill("DataFigure3_62c", tempTrack->getNSigmasTPCProton(), 1/tempTrack->getPt());
                    insideprocessing.Fill("DataFigure3_62d", tempTrack->getNSigmasTPCElectron(), 1/tempTrack->getPt());
                }
            }
        }
        return 0;
        };

    TreeProc.Process(myFunction);

    outsideprocessing.Merge();

    //setting up a tree & output file
    string path = string(argv[0]);
    string outfileName;
    if(outputFolder.find(".root")!=std::string::npos){
        outfileName = outputFolder;
    } else{
        outfileName = outputFolder+"AnaOutput_"+path.substr(path.find_last_of("/\\")+1)+".root";
    }
    cout<<"Created output file "<<outfileName<<endl;
    TFile *outputFileHist = TFile::Open(outfileName.c_str(), "recreate");
    outsideprocessing.SaveToFile(outputFileHist);
    outputFileHist->Close();

    return 0;
}