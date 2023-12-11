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
const double particleMass[nParticles] = { 0.13957, 0.497611, 0.93827 }; // pion, kaon, proton in GeV /c^2 
enum BRANCH_ID{ EU, ED, WU, WD, nBranches };
enum RP_ID{ E1U, E1D, E2U, E2D, W1U, W1D, W2U, W2D, nRomanPots };

int main(int argc, char **argv){

    int nthreads = 2;
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
    double kaonMassWindowNarrowLow = 0.48;
    double kaonMassWindowNarrowHigh = 0.51;
    double kaonMassWindowWideLow = 0.44;
    double kaonMassWindowWideHigh = 0.54;
    double kaonMassWindowPresentationLow = 0.46;
    double kaonMassWindowPresentationHigh = 0.53;
    vector<vector<double>> beamData = ReadFillPositionData("../../../../share/Run7PolarizationWithPosition.csv");

    //histograms
    ProcessingOutsideLoop outsideprocessing;
    outsideprocessing.AddHistogram(TH1D("Mpipibefore", "K^{0}_{S} mass;m_{#pi^{+}#pi^{-}} [GeV];Number of pairs", 100, kaonMassWindowWideLow, kaonMassWindowWideHigh));
    outsideprocessing.AddHistogram(TH1D("DCApipiK0", "DCA between #pi^{#pm} from K0 K^{0}_{S};DCA_{#pi^{+}#pi^{-}-K^{0}_{S}};Number of pairs", 50, 0, 5));
    outsideprocessing.AddHistogram(TH1D("DCApipiPV", "DCA between #pi^{#pm} from vertex when K^{0}_{S} in narrow mass window;DCA_{#pi^{+}#pi^{-}-PV};Number of pairs", 50, 0, 5));
    outsideprocessing.AddHistogram(TH1D("DCAK0PV", "DCA between K^{0}_{S} and vertex;DCA_{#pi^{+}#pi^{-}-K^{0}_{S}};Number of pairs", 50, 0, 5));
    outsideprocessing.AddHistogram(TH1D("LogProtons", "log(#xi_{E}*#xi_{W});log(#xi_{E}*#xi_{W});events", 100, -10, 0));
    outsideprocessing.AddHistogram(TH1D("DivProtons", "ln(#xi_{E}/#xi_{W});ln(#xi_{E}/#xi_{W});events", 100, -10, 10));
    outsideprocessing.AddHistogram(TH2D("Log2DProtons", "log#xi_{W} vs log#xi_{E};log#xi_{E};log#xi_{W}", 60, -5, 1, 60, -5, 1));
    outsideprocessing.AddHistogram(TH1D("LogEproton", "log#xi_{E};log#xi_{E};events", 60, -5, 1));
    outsideprocessing.AddHistogram(TH1D("LogWproton", "log#xi_{W};log#xi_{W};events", 60, -5, 1));

    // int triggers[] = { 570701, 570705, 570711, 590701, 590705, 590708 };
    // outsideprocessing.AddHistogram(TH1D("triggerHist", "Data triggers;Trigger ID;Number of events", 6, 0, 6));
    // for(int i = 0;i<6;i++){
    //     outsideprocessing.GetPointer1D(1)->GetXaxis()->SetBinLabel(i+1, to_string(triggers[i]).c_str());
    // }

    //processing
    //defining TreeProcessor
    ROOT::TTreeProcessorMT TreeProc(*upcChain, nthreads);

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
        std::vector<StUPCTrack *> vector_Track;
        int K0_pair_indices[] = { -1,-1 };
        int vertex_pair_indices[] = { -1,-1 };
        StUPCV0 *K0_pair = nullptr;
        StUPCV0 *vertex_pair = nullptr;
        StUPCV0 *tempParticle = nullptr;
        vector<double> tempBeamVector;
        double beamValues[4];
        TVector3 vertexPrimary;
        StUPCRpsTrack *eastTrack;
        StUPCRpsTrack *westTrack;

        //actual loop
        while(myReader.Next()){
            //in a TTree, it *would* be constant, in TChain however not necessarily
            tempUPCpointer = StUPCEventInstance.Get();
            tempRPpointer = StRPEventInstance.Get();

            //cleaning the loop
            //these deletes were unnecessary
            // delete K0_pair;
            // delete vertex_pair;
            K0_pair_indices[0] = -1;
            K0_pair_indices[1] = -1;
            vertex_pair_indices[0] = -1;
            vertex_pair_indices[1] = -1;
            vector_Track.clear();

            //cuts & histogram filling
            //selecting TOF tracks
            for(int i = 0; i<tempUPCpointer->getNumberOfTracks(); i++){
                if(tempUPCpointer->getTrack(i)->getFlag(StUPCTrack::kTof)){
                    vector_Track.push_back(tempUPCpointer->getTrack(i));
                }
            }
            //selecting two highest kaons
            vertexPrimary = { tempUPCpointer->getVertex(0)->getPosX(), tempUPCpointer->getVertex(0)->getPosY(), tempUPCpointer->getVertex(0)->getPosZ() };
            tempBeamVector = FindPosition(tempUPCpointer->getFillNumber(), vertexPrimary.Z(), beamData[0], beamData[1], beamData[2], beamData[3], beamData[4], beamData[5], beamData[6], beamData[7], beamData[8]);
            beamValues[0] = tempBeamVector[0];
            beamValues[1] = tempBeamVector[1];
            beamValues[2] = tempBeamVector[2];
            beamValues[3] = tempBeamVector[3];
            //actual loop
            for(long unsigned int i = 0; i<vector_Track.size()-1; i++){
                for(long unsigned int j = i+1; j<vector_Track.size(); j++){
                    tempParticle = new StUPCV0(vector_Track[i], vector_Track[j], particleMass[0], particleMass[0], 1, 1, vertexPrimary, beamValues, tempUPCpointer->getMagneticField(), true);
                    //tests if accept the particle
                    bool K0test1 = vector_Track[i]->getCharge()*vector_Track[j]->getCharge()<0;
                    bool K0test2 = tempParticle->dcaDaughters()<1.5;
                    bool K0test3 = tempParticle->pointingAngleHypo()>0.925;
                    bool K0test4 = tempParticle->DCABeamLine()<1.5;
                    bool K0test5 = tempParticle->m()>kaonMassWindowPresentationLow||tempParticle->m()<kaonMassWindowPresentationHigh;
                    if(!(K0test1&&K0test2&&K0test3&&K0test4&&K0test5)){
                        continue;
                    }
                    //filling
                    if(K0_pair_indices[0]<0){
                        K0_pair_indices[0] = i;
                        K0_pair_indices[1] = j;
                        K0_pair = new StUPCV0(vector_Track[i], vector_Track[j], particleMass[0], particleMass[0], 1, 1, vertexPrimary, beamValues, tempUPCpointer->getMagneticField(), true);
                    } else if(abs(tempParticle->m()-particleMass[0])<abs(K0_pair->m()-particleMass[0])){
                        K0_pair_indices[0] = i;
                        K0_pair_indices[1] = j;
                        delete K0_pair;
                        K0_pair = new StUPCV0(vector_Track[i], vector_Track[j], particleMass[0], particleMass[0], 1, 1, vertexPrimary, beamValues, tempUPCpointer->getMagneticField(), true);
                    }
                    //finishing
                    delete tempParticle;
                }
            }
            //check if any K0 particle was found
            if(K0_pair_indices[0]<0){
                continue;
            }
            //setting second vertex
            //under assumption there are 4 particles
            bool isvertexPresent = false;
            vector<int> temp = { 0, 1, 2, 3 };
            temp.erase(std::find(temp.begin(), temp.end(), K0_pair_indices[0]));
            temp.erase(std::find(temp.begin(), temp.end(), K0_pair_indices[1]));
            vertex_pair_indices[0] = temp[0];
            vertex_pair_indices[1] = temp[1];
            vertex_pair = new StUPCV0(vector_Track[vertex_pair_indices[0]], vector_Track[vertex_pair_indices[1]], particleMass[0], particleMass[0], 1, 1, vertexPrimary, beamValues, tempUPCpointer->getMagneticField(), true);
            //tests to see if vertex suggestion is okay
            TVector3 correctedVertex;
            bool PVtest1 = vertex_pair->dcaDaughters()<1.5;
            bool PVtest2 = vertex_pair->DCABeamLine()<1.5;
            if(PVtest1&&PVtest2){
                correctedVertex = vertex_pair->decayVertex();
                isvertexPresent = true;
            } else{
                correctedVertex = { beamValues[0], beamValues[1], tempUPCpointer->getVertex(0)->getPosZ() };
            }

            // TH1D("Mpipibefore", "K^{0}_{S} mass;m_{#pi^{+}#pi^{-}} [GeV];Number of pairs", 100, kaonMassWindowWideLow, kaonMassWindowWideHigh));
            // TH1D("DCApipiK0", "DCA between #pi^{#pm} from K0 K^{0}_{S};DCA_{#pi^{+}#pi^{-}-K^{0}_{S}};Number of pairs", 50, 0, 5));
            // TH1D("DCApipiPV", "DCA between #pi^{#pm} from vertex K^{0}_{S} (PV) when in narrow mass window;DCA_{#pi^{+}#pi^{-}-PV};Number of pairs", 50, 0, 5));
            // TH1D("DCAK0PV", "DCA between K0 K^{0}_{S} and vertex K^{0}_{S};DCA_{#pi^{+}#pi^{-}-K^{0}_{S}};Number of pairs", 50, 0, 5));
            // TH1D("LogProtons", "log(#xi_{E}*#xi_{W});log(#xi_{E}*#xi_{W});events", 100, -10, 0));
            // TH1D("DivProtons", "ln(#xi_{E}/#xi_{W});ln(#xi_{E}/#xi_{W});events", 100, -10, 10));
            // TH2D("Log2DProtons", "log#xi_{W} vs log#xi_{E};log#xi_{E};log#xi_{W}", 60, -5, 1, 60, -5, 1));
            // TH1D("LogEproton", "log#xi_{E};log#xi_{E};events", 60, -5, 1));
            // TH1D("LogWproton", "log#xi_{W};log#xi_{W};events", 60, -5, 1));

            insideprocessing.Fill(0, K0_pair->m());
            insideprocessing.Fill(1, K0_pair->dcaDaughters());
            if(isvertexPresent&&K0_pair->m()>kaonMassWindowNarrowLow&&K0_pair->m()<kaonMassWindowNarrowHigh){
                insideprocessing.Fill(2, vertex_pair->dcaDaughters());
                insideprocessing.Fill(3, (K0_pair->decayVertex()-vertex_pair->decayVertex()).Mag());
            }

            //filter to filter out badly reconstructed protons
            //as in with xi>1, cause those with xi<0 will get log(xi)=NaN and get registered as overflow
            if(tempRPpointer->getTrack(0)->xi(255.0)>1||tempRPpointer->getTrack(1)->xi(255.0)>1){
                cout<<log10(tempRPpointer->getTrack(0)->xi(255.0))<<" "<<log10(tempRPpointer->getTrack(1)->xi(255.0))<<endl;
                continue;
            }

            //histograms
            insideprocessing.Fill(4, log10(tempRPpointer->getTrack(0)->xi(255.0)*tempRPpointer->getTrack(1)->xi(255.0)));
            //0th track is east if branch <2
            if(tempRPpointer->getTrack(0)->branch()<2){
                eastTrack = tempRPpointer->getTrack(0);
                westTrack = tempRPpointer->getTrack(1);
            } else{
                eastTrack = tempRPpointer->getTrack(1);
                westTrack = tempRPpointer->getTrack(0);
            }
            insideprocessing.Fill(5, log(eastTrack->xi(255.0)/westTrack->xi(255.0)));
            insideprocessing.Fill(6, log10(eastTrack->xi(255.0)), log10(westTrack->xi(255.0)));
            insideprocessing.Fill(7, log10(eastTrack->xi(255.0)));
            insideprocessing.Fill(8, log10(westTrack->xi(255.0)));
        }
        return 0;
        };

    TreeProc.Process(myFunction);

    outsideprocessing.Merge();

    //setting up a tree & output file
    string path = string(argv[0]);
    string outfileName = outputFolder+"AnaOutput_"+path.substr(path.find_last_of("/\\")+1)+".root";
    cout<<"Created output file "<<outfileName<<endl;
    TFile *outputFileHist = TFile::Open(outfileName.c_str(), "recreate");
    outsideprocessing.SaveToFile(outputFileHist);
    outputFileHist->Close();

    return 0;
}
