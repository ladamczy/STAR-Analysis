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
enum SUSPECTED_PARTICLES{ K0S, Lambda, Kstar, Phi };

//cuts in here (except trigger ones) taken from paper (both Drupal link and PDF link):
//https://drupal.star.bnl.gov/STAR/pwg/lfs-upc/CEP-510-GeV/PWG-Paper-Proposal
//https://drupal.star.bnl.gov/STAR/system/files/CEP510_proposal_0.pdf
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
    //other
    outsideprocessing.AddHistogram(TH1D("MissingpT", ";p_{T} [GeV];events", 100, 0., 1.));
    outsideprocessing.AddHistogram(TH1D("MissingpTWide", "MissingpT", 250, 0., 5.));
    outsideprocessing.AddHistogram(TH1D("MKKExtraNarrowNoVeto", ";m_{K^{+}K^{-}} [GeV];Number of pairs", 140, 0.98, 1.05));
    outsideprocessing.AddHistogram(TH1D("MKKWideNoVeto", ";m_{K^{+}K^{-}} [GeV];Number of pairs", 500, 0., 5.));
    outsideprocessing.AddHistogram(TH1D("MpipiExtraNarrowNoVeto", ";m_{#pi^{+}#pi^{-}} [GeV];Number of pairs", 100, 0.3, 0.7));
    outsideprocessing.AddHistogram(TH1D("MKKLikeInPaper", ";m_{K^{+}K^{-}} [GeV];Number of pairs", 30, 0.9, 2.4));
    outsideprocessing.AddHistogram(TH1D("M2TOFLikeInPaper", ";m^{2}_{TOF} [GeV^{2}];Number of pairs", 200, -0.5, 1.5));
    outsideprocessing.AddHistogram(TH1D("M2TOFpipiLikeInPaper", ";m^{2}_{TOF} [GeV^{2}];Number of pairs", 200, -0.5, 1.5));
    outsideprocessing.AddHistogram(TH1D("M2TOFKKLikeInPaper", ";m^{2}_{TOF} [GeV^{2}];Number of pairs", 200, -0.5, 1.5));
    outsideprocessing.AddHistogram(TH1D("M2TOFppLikeInPaper", ";m^{2}_{TOF} [GeV^{2}];Number of pairs", 200, -0.5, 1.5));
    outsideprocessing.AddHistogram(TH1D("pT", ";p_{T} [GeV];tracks", 100, 0., 1.));
    outsideprocessing.AddHistogram(TH1D("eta", ";#eta;tracks", 220, -1.1, 1.1));
    outsideprocessing.AddHistogram(TH1D("vertices", ";vertices;events", 10, 0, 10));
    outsideprocessing.AddHistogram(TH2D("xi", ";#xi_{E};#xi_{W}", 150, -0.05, 0.25, 150, -0.05, 0.25));

    //processing
    //defining TreeProcessor
    ROOT::TTreeProcessorMT TreeProc(*upcChain, nthreads);

    //other things
    int eventsProcessed = 0;

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
        TVector3 tempMomentum;
        TVector3 totalMomentum;
        TLorentzVector temp4Vector1, temp4Vector2;
        bool goodQuality;
        double firstBranch, secondBranch;
        bool f1, f2, f3;
        double px, py;
        double m2TOF;
        double chi2proton, chi2kaon, chi2pion;

        //actual loop
        while(myReader.Next()){
            //in a TTree, it *would* be constant, in TChain however not necessarily
            tempUPCpointer = StUPCEventInstance.Get();
            tempRPpointer = StRPEventInstance.Get();

            //cleaning the loop
            vector_Track_positive.clear();
            vector_Track_negative.clear();
            totalMomentum = { 0, 0, 0 };
            goodQuality = true;

            //cause I want to see what's going on
            if(eventsProcessed%10000==0){
                cout<<"Processed "<<eventsProcessed<<" events"<<endl;
            }
            eventsProcessed++;

            //tests
            //trigger RP_CPT2noBBCL
            if(!tempUPCpointer->isTrigger(570705)){
                continue;
            }
            //2 tracks
            if(tempRPpointer->getNumberOfTracks()!=2){
                continue;
            }
            //1 track east, 1 track west (neat trick - assigning negative to east by
            //substracting 1.5, and if after multiplying they are <0, they are from opposite sides
            firstBranch = tempRPpointer->getTrack(0)->branch();
            secondBranch = tempRPpointer->getTrack(1)->branch();
            if((firstBranch-1.5)*(secondBranch-1.5)>0){
                continue;
            }
            //at least 3 out of 4 planes on both and both should have both RPs hit
            //also fiducial
            for(unsigned int k = 0; k<tempRPpointer->getNumberOfTracks(); ++k){
                // Get pointer to k-th track in Roman Pot data collection
                StUPCRpsTrack* trk = tempRPpointer->getTrack(k);
                trk->setEvent(tempRPpointer);
                //there were problems with apparently not having track point like, entirely???
                //so the first is check point if they do have them
                //and then if points are of good quality
                if(trk->getTrackPoint(0)==nullptr||trk->getTrackPoint(1)==nullptr){
                    goodQuality = false;
                    break;
                }
                //check if track has at least 3 of 4 RP planes used
                if(trk->getTrackPoint(0)->planesUsed()<3||trk->getTrackPoint(1)->planesUsed()<3){
                    goodQuality = false;
                    break;
                }

                //fiducial
                px = trk->pVec().X();
                py = trk->pVec().Y();
                f1 = (0.4<abs(py)&&abs(py)<0.8);
                f2 = (-0.27<px);
                f3 = (pow(px+0.6, 2)+pow(py, 2)<1.25);
                if(!(f1&&f2&&f3)){
                    goodQuality = false;
                    break;
                }
                totalMomentum += trk->pVec();
            }
            if(!goodQuality){ continue; }
            //tracks checks
            for(int i = 0; i<tempUPCpointer->getNumberOfTracks(); i++){
                tempTrack = tempUPCpointer->getTrack(i);
                if(!tempTrack->getFlag(StUPCTrack::kTof) or !tempTrack->getFlag(StUPCTrack::kPrimary)){
                    continue;
                }
                if(tempTrack->getNhitsFit()<20 or tempTrack->getNhitsDEdx()<15){
                    continue;
                }
                if(tempTrack->getDcaXY()>=1.5 or abs(tempTrack->getDcaZ())>=1.0){
                    continue;
                }
                if(abs(tempTrack->getVertex()->getPosZ())>=100 or abs(tempTrack->getEta())>(1.0-1./250*abs(tempTrack->getVertex()->getPosZ()))){
                    continue;
                }
                if(tempTrack->getCharge()>0){
                    vector_Track_positive.push_back(tempTrack);
                } else{
                    vector_Track_negative.push_back(tempTrack);
                }
            }
            //looking for a pair with the same vertex
            if(vector_Track_positive.size()!=1 or vector_Track_negative.size()!=1){
                continue;
            }
            if(vector_Track_positive[0]->getVertexId()!=vector_Track_negative[0]->getVertexId()){
                continue;
            }
            vector_Track_positive[0]->getMomentum(tempMomentum);
            totalMomentum += tempMomentum;
            vector_Track_negative[0]->getMomentum(tempMomentum);
            totalMomentum += tempMomentum;
            insideprocessing.Fill("MissingpT", totalMomentum.Pt());
            insideprocessing.Fill("MissingpTWide", totalMomentum.Pt());
            // //pT cut at 120 MeV
            // if(totalMomentum.Pt()>=0.12){
            //     continue;
            // }
            //test pT cut at 60 MeV
            if(totalMomentum.Pt()>=0.06){
                continue;
            }
            vector_Track_positive[0]->getLorentzVector(temp4Vector1, particleMass[Kaon]);
            vector_Track_negative[0]->getLorentzVector(temp4Vector2, particleMass[Kaon]);
            insideprocessing.Fill("MKKExtraNarrowNoVeto", (temp4Vector1+temp4Vector2).M());
            insideprocessing.Fill("MKKWideNoVeto", (temp4Vector1+temp4Vector2).M());
            vector_Track_positive[0]->getLorentzVector(temp4Vector1, particleMass[Pion]);
            vector_Track_negative[0]->getLorentzVector(temp4Vector2, particleMass[Pion]);
            insideprocessing.Fill("MpipiExtraNarrowNoVeto", (temp4Vector1+temp4Vector2).M());
            insideprocessing.Fill("pT", vector_Track_positive[0]->getPt());
            insideprocessing.Fill("pT", vector_Track_negative[0]->getPt());
            insideprocessing.Fill("eta", vector_Track_positive[0]->getEta());
            insideprocessing.Fill("eta", vector_Track_negative[0]->getEta());
            insideprocessing.Fill("vertices", tempUPCpointer->getNumberOfVertices());
            if(tempRPpointer->getTrack(0)->pVec().Z()>0){
                insideprocessing.Fill("xi", tempRPpointer->getTrack(1)->xi(254.867), tempRPpointer->getTrack(0)->xi(254.867));
            } else{
                insideprocessing.Fill("xi", tempRPpointer->getTrack(0)->xi(254.867), tempRPpointer->getTrack(1)->xi(254.867));
            }
            //making cuts like in the paper
            m2TOF = M2TOF(vector_Track_positive[0], vector_Track_negative[0]);
            chi2pion = pow(vector_Track_positive[0]->getNSigmasTPCPion(), 2)+pow(vector_Track_negative[0]->getNSigmasTPCPion(), 2);
            chi2kaon = pow(vector_Track_positive[0]->getNSigmasTPCKaon(), 2)+pow(vector_Track_negative[0]->getNSigmasTPCKaon(), 2);
            chi2proton = pow(vector_Track_positive[0]->getNSigmasTPCProton(), 2)+pow(vector_Track_negative[0]->getNSigmasTPCProton(), 2);
            insideprocessing.Fill("M2TOFLikeInPaper", m2TOF);
            if(chi2pion>9&&chi2kaon>9&&chi2proton<9){
                insideprocessing.Fill("M2TOFppLikeInPaper", m2TOF);
            } else if(chi2pion>9&&chi2kaon<9&&chi2proton>9){
                insideprocessing.Fill("M2TOFKKLikeInPaper", m2TOF);
                if(m2TOF<0.6&&m2TOF>0.15){
                    vector_Track_positive[0]->getLorentzVector(temp4Vector1, particleMass[Kaon]);
                    vector_Track_negative[0]->getLorentzVector(temp4Vector2, particleMass[Kaon]);
                    insideprocessing.Fill("MKKLikeInPaper", (temp4Vector1+temp4Vector2).M());
                }
            } else if(chi2pion<12){
                insideprocessing.Fill("M2TOFpipiLikeInPaper", m2TOF);
            }

            //lambda finish
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
    TFile* outputFileHist = TFile::Open(outfileName.c_str(), "recreate");
    outsideprocessing.SaveToFile(outputFileHist);
    outputFileHist->Close();

    return 0;
}