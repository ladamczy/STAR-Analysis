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

const double particleMass[3] = { 0.13957, 0.497611, 0.93827 }; // pion, kaon, proton in GeV /c^2 

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
    vector<vector<double>> beamData = ReadFillPositionData("../../../../share/Run7PolarizationWithPosition.csv");

    //histograms
    ProcessingOutsideLoop outsideprocessing;
    outsideprocessing.AddHistogram(TH1D("MpipiforMatchHist", "K^{0}_{S} mass;m_{#pi^{+}#pi^{-}} [GeV];Number of pairs", 100, kaonMassWindowWideLow, kaonMassWindowWideHigh));
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
        int leading_pair_indices[] = { -1,-1 };
        int subleading_pair_indices[] = { -1,-1 };
        StUPCV0 *leading_particle = nullptr;
        StUPCV0 *subleading_particle = nullptr;
        StUPCV0 *tempParticle = nullptr;
        vector<double> tempBeamVector;
        double beamValues[4];
        TVector3 vertexPrimary;

        //actual loop
        while(myReader.Next()){
            //in a TTree, it *would* be constant, in TChain however not necessarily
            tempUPCpointer = StUPCEventInstance.Get();
            tempRPpointer = StRPEventInstance.Get();

            //cleaning the loop
            //these deletes were unnecessary
            // delete leading_particle;
            // delete subleading_particle;
            leading_pair_indices[0] = -1;
            leading_pair_indices[1] = -1;
            subleading_pair_indices[0] = -1;
            subleading_pair_indices[1] = -1;
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
                    tempParticle = new StUPCV0(vector_Track[i], vector_Track[j], particleMass[0], particleMass[0], 1, 1, vertexPrimary, beamValues, tempUPCpointer->getMagneticField(), false, true);
                    if(tempParticle->m()<kaonMassWindowWideLow||tempParticle->m()>kaonMassWindowWideHigh||vector_Track[i]->getCharge()*vector_Track[j]->getCharge()>0){
                        continue;
                    }
                    //testing the place it should belong
                    if(leading_pair_indices[0]<0){
                        leading_pair_indices[0] = i;
                        leading_pair_indices[1] = j;
                        leading_particle = new StUPCV0(vector_Track[i], vector_Track[j], particleMass[0], particleMass[0], 1, 1, vertexPrimary, beamValues, tempUPCpointer->getMagneticField(), false, true);
                    } else if(abs(tempParticle->m()-particleMass[0])<abs(leading_particle->m()-particleMass[0])){
                        subleading_pair_indices[0] = leading_pair_indices[0];
                        subleading_pair_indices[1] = leading_pair_indices[1];
                        leading_pair_indices[0] = i;
                        leading_pair_indices[1] = j;
                        delete subleading_particle;
                        subleading_particle = new StUPCV0(vector_Track[subleading_pair_indices[0]], vector_Track[subleading_pair_indices[1]], particleMass[0], particleMass[0], 1, 1, vertexPrimary, beamValues, tempUPCpointer->getMagneticField(), false, true);
                        delete leading_particle;
                        leading_particle = new StUPCV0(vector_Track[i], vector_Track[j], particleMass[0], particleMass[0], 1, 1, vertexPrimary, beamValues, tempUPCpointer->getMagneticField(), false, true);
                    } else if(abs(tempParticle->m()-particleMass[0])>abs(leading_particle->m()-particleMass[0])&&subleading_pair_indices[0]<0){
                        subleading_pair_indices[0] = i;
                        subleading_pair_indices[1] = j;
                        delete subleading_particle;
                        subleading_particle = new StUPCV0(vector_Track[i], vector_Track[j], particleMass[0], particleMass[0], 1, 1, vertexPrimary, beamValues, tempUPCpointer->getMagneticField(), false, true);
                    } else if((abs(tempParticle->m()-particleMass[0])>abs(leading_particle->m()-particleMass[0]))&&(abs(tempParticle->m()-particleMass[0])<abs(subleading_particle->m()-particleMass[0]))){
                        subleading_pair_indices[0] = i;
                        subleading_pair_indices[1] = j;
                        delete subleading_particle;
                        subleading_particle = new StUPCV0(vector_Track[i], vector_Track[j], particleMass[0], particleMass[0], 1, 1, vertexPrimary, beamValues, tempUPCpointer->getMagneticField(), false, true);
                    }
                    delete tempParticle;
                }
            }
            //check if any leading particle was found
            if(leading_pair_indices[0]<0){
                continue;
            }
            //setting second vertex
            TVector3 correctedVertex;
            if(subleading_pair_indices[0]<0){
                correctedVertex = { beamValues[0], beamValues[1], tempUPCpointer->getVertex(0)->getPosZ() };
            } else{
                correctedVertex = subleading_particle->decayVertex();
            }

            insideprocessing.Fill(0, leading_particle->m());

            // //used triggers histograms
            // for(int i = 0;i<6;i++){
            //     if(tempUPCpointer->isTrigger(triggers[i])){
            //         insideprocessing.Fill(1, i+1);
            //     }
            // }
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
