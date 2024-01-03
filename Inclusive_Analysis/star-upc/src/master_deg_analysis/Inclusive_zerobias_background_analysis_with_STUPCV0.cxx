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

bool IsInXiElasticSpot(StUPCRpsTrack *, StUPCRpsTrack *);
bool IsInMomElasticSpot(StUPCRpsTrack *, StUPCRpsTrack *);

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
    vector<vector<double>> beamData = ReadFillPositionData("STAR-Analysis/share/Run7PolarizationWithPosition.csv");

    //histograms
    ProcessingOutsideLoop outsideprocessing;
    outsideprocessing.AddHistogram(TH1D("LogProtons", "log(#xi_{E}*#xi_{W});log(#xi_{E}*#xi_{W});events", 100, -10, 0));
    outsideprocessing.AddHistogram(TH1D("DivProtons", "ln(#xi_{E}/#xi_{W});ln(#xi_{E}/#xi_{W});events", 100, -10, 10));
    outsideprocessing.AddHistogram(TH2D("Log2DProtons", "log#xi_{W} vs log#xi_{E};log#xi_{E};log#xi_{W}", 60, -5, 1, 60, -5, 1));
    outsideprocessing.AddHistogram(TH1D("LogEproton", "log#xi_{E};log#xi_{E};events", 60, -5, 1));
    outsideprocessing.AddHistogram(TH1D("LogWproton", "log#xi_{W};log#xi_{W};events", 60, -5, 1));
    outsideprocessing.AddHistogram(TH1D("XiEproton", "#xi_{E};#xi_{E};events", 84, -0.05, 1));
    outsideprocessing.AddHistogram(TH1D("XiWproton", "#xi_{W};#xi_{W};events", 84, -0.05, 1));
    outsideprocessing.AddHistogram(TH1D("XiEprotoncloser", "#xi_{E};#xi_{E};events", 400, -0.05, 0.15));
    outsideprocessing.AddHistogram(TH1D("XiWprotoncloser", "#xi_{W};#xi_{W};events", 400, -0.05, 0.15));
    outsideprocessing.AddHistogram(TH2D("Xi2DProtons", "#xi_{W} vs #xi_{E};#xi_{E};#xi_{W}", 60, -0.01, 0.05, 60, -0.01, 0.05));
    outsideprocessing.AddHistogram(TH2D("sumTheta2DProtons", "Sum of #theta_{x} and #theta_{y} of protons;#Delta#theta_{x};#Delta#theta_{y}", 100, -3e-3, 3e-3, 100, -2e-3, 2e-3));
    outsideprocessing.AddHistogram(TH2D("sump2DProtons", "Sum of p_{x} and p_{y} of protons;#Sigmap_{x};#Sigmap_{y}", 100, -2, 2, 100, -2, 2));
    outsideprocessing.AddHistogram(TH2D("sump2DProtonsExact", "Sum of p_{x} and p_{y} of protons;#Sigmap_{x};#Sigmap_{y}", 80, -0.6, 1., 50, -0.5, 0.5));

    outsideprocessing.AddHistogram(TH1D("XiEprotoncloserAfterElasticCut", "#xi_{E};#xi_{E};events", 400, -0.05, 0.15));
    outsideprocessing.AddHistogram(TH1D("XiWprotoncloserAfterElasticCut", "#xi_{W};#xi_{W};events", 400, -0.05, 0.15));

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
        // StUPCEvent *tempUPCpointer;
        StRPEvent *tempRPpointer;
        insideprocessing.GetLocalHistograms(&outsideprocessing);

        //helpful variables
        vector<double> tempBeamVector;
        StUPCRpsTrack *eastTrack;
        StUPCRpsTrack *westTrack;

        //actual loop
        while(myReader.Next()){
            //in a TTree, it *would* be constant, in TChain however not necessarily
            // tempUPCpointer = StUPCEventInstance.Get();
            tempRPpointer = StRPEventInstance.Get();

            //0
            // TH1D("Mpipibefore", "K^{0}_{S} mass;m_{#pi^{+}#pi^{-}} [GeV];Number of pairs", 100, kaonMassWindowWideLow, kaonMassWindowWideHigh));
            // TH1D("DCApipiK0", "DCA between #pi^{#pm} from K0 K^{0}_{S};DCA_{#pi^{+}#pi^{-}-K^{0}_{S}};Number of pairs", 50, 0, 5));
            // TH1D("DCApipiPV", "DCA between #pi^{#pm} from vertex K^{0}_{S} (PV) when in narrow mass window;DCA_{#pi^{+}#pi^{-}-PV};Number of pairs", 50, 0, 5));
            // TH1D("DCAK0PV", "DCA between K0 K^{0}_{S} and vertex K^{0}_{S};DCA_{#pi^{+}#pi^{-}-K^{0}_{S}};Number of pairs", 50, 0, 5));
            // TH1D("LogProtons", "log(#xi_{E}*#xi_{W});log(#xi_{E}*#xi_{W});events", 100, -10, 0));
            //5
            // TH1D("DivProtons", "ln(#xi_{E}/#xi_{W});ln(#xi_{E}/#xi_{W});events", 100, -10, 10));
            // TH2D("Log2DProtons", "log#xi_{W} vs log#xi_{E};log#xi_{E};log#xi_{W}", 60, -5, 1, 60, -5, 1));
            // TH1D("LogEproton", "log#xi_{E};log#xi_{E};events", 60, -5, 1));
            // TH1D("LogWproton", "log#xi_{W};log#xi_{W};events", 60, -5, 1));
            // TH2D("dcaDaughtersvsMass", ";m_{#pi^{+}#pi^{-}} [GeV];dcaDaughters", 100, kaonMassWindowWideLow, kaonMassWindowWideHigh, 20, 0, 10));
            //10
            // TH2D("pointingAngleHypovsMass", ";m_{#pi^{+}#pi^{-}} [GeV];pointingAngleHypo", 100, kaonMassWindowWideLow, kaonMassWindowWideHigh, 20, -1, 1));
            // TH2D("DCABeamLinevsMass", ";m_{#pi^{+}#pi^{-}} [GeV];DCABeamLine", 100, kaonMassWindowWideLow, kaonMassWindowWideHigh, 20, 0, 10));
            // TH1D("XiEproton", "#xi_{E};#xi_{E};events", 84, -0.05, 1));
            // TH1D("XiWproton", "#xi_{W};#xi_{W};events", 84, -0.05, 1));
            // TH1D("XiEprotoncloser", "#xi_{E};#xi_{E};events", 40, -0.05, 0.15));
            //15
            // TH1D("XiWprotoncloser", "#xi_{W};#xi_{W};events", 40, -0.05, 0.15));
            // TH1D("vertex_pair_dcaDaughters", "Vertex pair dcaDaughters();dcaDaughters;events", 60, 0, 3));
            // TH1D("vertex_pair_DCABeamLine", "Vertex pair DCABeamLine();DCABeamLine;events", 60, 0, 3));
            // TH1D("K0decayLengthHypo", "decayLengthHypo() of K^{0}_{S} pair;length [cm];events", 60, 0, 3));
            // TH1D("vertexdecayLengthHypo", "decayLengthHypo() of vertex pair;length [cm];events", 60, 0, 3));
            //20
            // TH2D("decayLengthHypovsMass", "decayLengthHypo ();m_{#pi^{+}#pi^{-}} [GeV];decayLengthHypo() [cm]", 60, 0, 3, 70, kaonMassWindowPresentationLow, kaonMassWindowPresentationHigh));
            // TH1D("decayvertexZdifference", "#Delta Z of K^{0}_{S} and PV decay vertices;#Delta Z [cm];events", 200, -100, 100));
            // TH2D("Xi2DProtons", "#xi_{W} vs #xi_{E};#xi_{E};#xi_{W}", 400, -0.05, 0.15, 400, -0.05, 0.15));
            // TH2D("deltaTheta2DProtons", "Difference in #theta_{x} and #theta_{y} of protons;#Delta#theta_{x};#Delta#theta_{y}", 200, -0.2, 0.2, 200, -0.2, 0.2));
            // TH2D("sump2DProtons", "Sum of p_{x} and p_{y} of protons;#Sigmap_{x};#Sigmap_{y}", 100, -1, 1, 100, -1, 1));
            //25
            // TH2D("sump2DProtonsExact", "Sum of p_{x} and p_{y} of protons;#Sigmap_{x};#Sigmap_{y}", 80, -0.6, 1., 50, -0.5, 0.5));
            // TH2D("etavsK0Mass", "#eta in function of K^{0}_{S} candidate mass;m_{#pi^{+}#pi^{-}} [GeV];#eta", 70, kaonMassWindowPresentationLow, kaonMassWindowPresentationHigh, 60, -3, 3));
            // TH1D("XiEprotoncloserAfterElasticCut", "#xi_{E};#xi_{E};events", 400, -0.05, 0.15));
            // TH1D("XiWprotoncloserAfterElasticCut", "#xi_{W};#xi_{W};events", 400, -0.05, 0.15));
            // TH1D("MpipiAfterElasticCut", "K^{0}_{S} mass;m_{#pi^{+}#pi^{-}} [GeV];Number of pairs", 70, kaonMassWindowPresentationLow, kaonMassWindowPresentationHigh));

            //filter to filter out badly reconstructed protons
            //as in with xi>1, cause those with xi<0 will get log(xi)=NaN and get registered as overflow
            if(tempRPpointer->getTrack(0)->xi(255.0)>1||tempRPpointer->getTrack(1)->xi(255.0)>1){
                cout<<log10(tempRPpointer->getTrack(0)->xi(255.0))<<" "<<log10(tempRPpointer->getTrack(1)->xi(255.0))<<endl;
                continue;
            }

            //0th track is east if branch <2
            if(tempRPpointer->getTrack(0)->branch()<2){
                eastTrack = tempRPpointer->getTrack(0);
                westTrack = tempRPpointer->getTrack(1);
            } else{
                eastTrack = tempRPpointer->getTrack(1);
                westTrack = tempRPpointer->getTrack(0);
            }
            //histograms
            if(eastTrack->xi(255.0)<0.005||westTrack->xi(255.0)<0.005){
                insideprocessing.Fill("LogProtons", 2*log10(0.005));
            } else{
                insideprocessing.Fill("LogProtons", log10(tempRPpointer->getTrack(0)->xi(255.0)*tempRPpointer->getTrack(1)->xi(255.0)));
            }
            if(eastTrack->xi(255.0)<0.005||westTrack->xi(255.0)<0.005){
                insideprocessing.Fill("LogProtons", 0);
            } else{
                insideprocessing.Fill("DivProtons", log(eastTrack->xi(255.0)/westTrack->xi(255.0)));
            }
            insideprocessing.Fill("Log2DProtons", log10(eastTrack->xi(255.0)), log10(westTrack->xi(255.0)));
            insideprocessing.Fill("LogEproton", log10(eastTrack->xi(255.0)));
            insideprocessing.Fill("LogWproton", log10(westTrack->xi(255.0)));
            insideprocessing.Fill("XiEproton", eastTrack->xi(255.0));
            insideprocessing.Fill("XiWproton", westTrack->xi(255.0));
            insideprocessing.Fill("XiEprotoncloser", eastTrack->xi(255.0));
            insideprocessing.Fill("XiWprotoncloser", westTrack->xi(255.0));
            insideprocessing.Fill("Xi2DProtons", eastTrack->xi(255.0), westTrack->xi(255.0));
            insideprocessing.Fill("sumTheta2DProtons", eastTrack->theta(StUPCRpsTrack::rpsAngleThetaX)+westTrack->theta(StUPCRpsTrack::rpsAngleThetaX), eastTrack->theta(StUPCRpsTrack::rpsAngleThetaY)+westTrack->theta(StUPCRpsTrack::rpsAngleThetaY));
            insideprocessing.Fill("sump2DProtons", eastTrack->pVec().X()+westTrack->pVec().X(), eastTrack->pVec().Y()+westTrack->pVec().Y());
            insideprocessing.Fill("sump2DProtonsExact", eastTrack->pVec().X()+westTrack->pVec().X(), eastTrack->pVec().Y()+westTrack->pVec().Y());

            //part dedicated to check the elastic cut
            if(!IsInXiElasticSpot(eastTrack, westTrack)&&!IsInMomElasticSpot(eastTrack, westTrack)){
                insideprocessing.Fill("XiEprotoncloserAfterElasticCut", eastTrack->xi(255.0));
                insideprocessing.Fill("XiWprotoncloserAfterElasticCut", westTrack->xi(255.0));
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

bool IsInXiElasticSpot(StUPCRpsTrack *east, StUPCRpsTrack *west){
    double x_0 = 4.30588e-03;
    double sigma_x = 2.02340e-03;
    double y_0 = 1.72097e-03;
    double sigma_y = 2.26638e-03;
    return pow((east->xi(255.0)-x_0)/sigma_x, 2)+pow((west->xi(255.0)-y_0)/sigma_y, 2)<3*3;
}

bool IsInMomElasticSpot(StUPCRpsTrack *east, StUPCRpsTrack *west){
    double x_0 = -3.82151e-02;
    double sigma_x = 3.67545e-02;
    double y_0 = 1.98348e-03;
    double sigma_y = 3.40440e-02;
    double x = east->pVec().X()+west->pVec().X();
    double y = east->pVec().Y()+west->pVec().Y();
    return pow((x-x_0)/sigma_x, 2)+pow((y-y_0)/sigma_y, 2)<3*3;
}