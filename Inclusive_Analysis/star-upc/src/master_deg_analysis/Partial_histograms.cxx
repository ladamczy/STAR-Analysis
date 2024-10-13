//cpp headers
#include <algorithm>

//ROOT headers
#include <ROOT/TThreadedObject.hxx>
#include <TTreeReader.h>
#include <ROOT/TTreeProcessorMT.hxx>
#include <TH2D.h>

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
#include "Afterburner.h"

//my headers
#include "UsefulThings.h"
#include "PartialHistograms.h"

enum{
    kAll = 1, kCPT, kRP, kOneVertex, kTPCTOF,
    kTotQ, kMax
};
// enum SIDE{ E = 0, East = 0, W = 1, West = 1, nSides };
// enum PARTICLES{ Pion = 0, Kaon = 1, Proton = 2, nParticles };
const double particleMass[nParticles] = { 0.13957, 0.493677, 0.93827 }; // pion, kaon, proton in GeV /c^2 
const double K0mass = 0.497611;
// enum BRANCH_ID{ EU, ED, WU, WD, nBranches };
// enum RP_ID{ E1U, E1D, E2U, E2D, W1U, W1D, W2U, W2D, nRomanPots };

//necesssary cuts
bool eventAndRunCut(StRPEvent *, StUPCEvent *, std::vector<int> *, std::vector<int> *);         //ok
bool protonCuts(StRPEvent *, StUPCEvent *, std::vector<int> *, std::vector<int> *);             //ok
//histogram cuts
bool enoughClusters(StRPEvent *, StUPCEvent *, TH1 *, std::vector<int> *, std::vector<int> *);
bool oneVertex(StRPEvent *, StUPCEvent *, TH1 *, std::vector<int> *, std::vector<int> *);       //ok
bool NTOFTracks(StRPEvent *, StUPCEvent *, TH1 *, std::vector<int> *, std::vector<int> *);   //ok
bool OppositeCharges(StRPEvent *, StUPCEvent *, TH1 *, std::vector<int> *, std::vector<int> *);        //ok
bool KinematicCut(StRPEvent *, StUPCEvent *, TH1 *, std::vector<int> *, std::vector<int> *);    //consider it ok
bool QualityCut(StRPEvent *, StUPCEvent *, TH1 *, std::vector<int> *, std::vector<int> *);      //ok
bool K0massTest(StRPEvent *, StUPCEvent *, TH1 *, std::vector<int> *, std::vector<int> *);      //ok
bool LambdamassTest(StRPEvent *, StUPCEvent *, TH1 *, std::vector<int> *, std::vector<int> *);      //ok
//a special one
bool CounterCut(StRPEvent *, StUPCEvent *, TH1 *, std::vector<int> *, std::vector<int> *);

//global variables
double kaonMassWindowPresentationLow = 0.46;
double kaonMassWindowPresentationHigh = 0.53;
double lambdaMassWindowPresentationLow = 1.080;
double lambdaMassWindowPresentationHigh = 1.150;
double Lambdamass = 1.115;
//Afterburner, yet again
vector<vector<double>> beamData = ReadFillPositionData("STAR-Analysis/share/Run7PolarizationWithPosition.csv");

int main(int argc, char **argv){
    //preparing input & output
    TChain *upcChain = new TChain("mUPCTree");
    if(ConnectInput(argc, argv, upcChain)){
        cout<<"All files connected"<<endl;
    }
    const string &outputFolder = argv[2];

    TTreeReader myReader(upcChain);
    TTreeReaderValue<StUPCEvent> StUPCEventInstance(myReader, "mUPCEvent");
    TTreeReaderValue<StRPEvent> StRPEventInstance(myReader, "mRPEvent");
    StUPCEvent *tempUPCpointer;
    StRPEvent *tempRPpointer;
    PartialHistograms histObject;

    //Afterburner things
    LoadOffsetFile("STAR-Analysis/share/OffSetsCorrectionsRun17.list", mCorrection);

    //histograms
    TH1D VertexNumber("VertexNumber", "Number of primary vertices", 10, 0, 10);
    TH1D TOFhitsNumber("TOFhitsNumber", "Number of hits in TOF", 100, 0, 100);
    TH2D ptEtaHist("ptEtaHist", "p_{T} vs #eta", 40, -2, 2, 40, 0, 2);
    TH2D QualityCheck("QualityCheck", "NdE/dx vs Nhits", 50, 0, 50, 50, 0, 50);
    TH1D ClusterNumber("ClusterNumber", "Cluster hit number", 50, 0, 50);
    TH1D K0massCheck("K0massCheck", "K0 mass check", 70, kaonMassWindowPresentationLow, kaonMassWindowPresentationHigh);
    TH1D LambdamassCheck("LambdamassCheck", "Lambda mass check", 70, lambdaMassWindowPresentationLow, lambdaMassWindowPresentationHigh);

    //just in case
    TH1D dummyHist("dummyHist", "dummy histogram", 1, 0, 1);
    TH1D counterHist("counterHist", "counter", 1, 0, 1);

    //necessary cuts
    // histObject.AddNecessaryCut(eventAndRunCut, true);
    histObject.AddNecessaryCut(protonCuts, true);

    //histogram cuts
    histObject.AddHistogramCut(oneVertex, &VertexNumber, true);
    histObject.AddHistogramCut(NTOFTracks, &TOFhitsNumber, true);
    histObject.AddHistogramCut(OppositeCharges, &dummyHist, true);
    histObject.AddHistogramCut(KinematicCut, &ptEtaHist, true);
    histObject.AddHistogramCut(QualityCut, &QualityCheck, true);
    histObject.AddHistogramCut(enoughClusters, &ClusterNumber, true);
    histObject.AddHistogramCut(K0massTest, &K0massCheck, true);
    histObject.AddHistogramCut(LambdamassTest, &LambdamassCheck, true);
    // histObject.AddHistogramCut(CounterCut, &counterHist, true);

    int counter = 0;
    while(myReader.Next()){
        counter++;
        if(counter%100000==0){
            cout<<"Analysed entry no "<<counter<<endl;
        }
        //in a TTree, it *would* be constant, in TChain however not necessarily
        tempUPCpointer = StUPCEventInstance.Get();
        tempRPpointer = StRPEventInstance.Get();
        //AFTERBURNER ALREADY APPLIED
        //TURN BACK BEFORE REAL STUFF
        //modified for afterburner
        // tempRPpointer = new StRPEvent(*StRPEventInstance.Get());
        // tempRPpointer->clearEvent();
        // runAfterburner(StRPEventInstance.Get(), tempRPpointer, tempUPCpointer->getRunNumber());

        //actual  processing
        histObject.SetEventPointers(tempRPpointer, tempUPCpointer);
        histObject.ProcessEvent();

        //AFTERBURNER ALREADY APPLIED
        //TURN BACK BEFORE REAL STUFF
        //final in-loop cleaning after Afterburner
        // delete tempRPpointer;
    }

    //setting up a tree & output file
    string path = string(argv[0]);
    string outfileName = outputFolder+"AnaOutput_"+path.substr(path.find_last_of("/\\")+1)+".root";
    cout<<"Created output file "<<outfileName<<endl;
    TFile *outputFileHist = TFile::Open(outfileName.c_str(), "recreate");

    //histograms
    outputFileHist->cd();
    ClusterNumber.Write();
    VertexNumber.Write();
    TOFhitsNumber.Write();
    ptEtaHist.Write();
    QualityCheck.Write();
    K0massCheck.Write();
    LambdamassCheck.Write();
    // counterHist.Write();

    outputFileHist->Close();
    return 0;
}






///////////////////
//cuts
///////////////////





//necessary cuts
bool eventAndRunCut(StRPEvent *RPEvent, StUPCEvent *UPCEvent, std::vector<int> *RPTrackIDs, std::vector<int> *UPCTrackIDs){
    // if(UPCEvent->getRunNumber()<=18083025&&UPCEvent->isTrigger(570704)){
    //     return true;
    // }
    // return false;
    // 570701, 570705, 570711
    if(UPCEvent->isTrigger(570701)||UPCEvent->isTrigger(570705)||UPCEvent->isTrigger(570711)){
        return true;
    }
    return false;
}

bool protonCuts(StRPEvent *RPEvent, StUPCEvent *UPCEvent, std::vector<int> *RPTrackIDs, std::vector<int> *UPCTrackIDs){
    //2 tracks
    if(RPTrackIDs->size()!=2){
        return false;
    }
    //1 track east, 1 track west (neat trick - assigning negative to east by
    //substracting 1.5, and if after multiplying they are <0, they are from opposite sides
    double firstBranch = RPEvent->getTrack((*RPTrackIDs)[0])->branch();
    double secondBranch = RPEvent->getTrack((*RPTrackIDs)[1])->branch();
    if((firstBranch-1.5)*(secondBranch-1.5)>0){
        return false;
    }
    //at least 3 out of 4 planes on both and both should have both RPs hit
    //also fiducial
    for(unsigned int k = 0; k<RPTrackIDs->size(); ++k){
        // Get pointer to k-th track in Roman Pot data collection
        StUPCRpsTrack *trk = RPEvent->getTrack((*RPTrackIDs)[k]);
        trk->setEvent(RPEvent);
        //there were problems with apparently not having track point like, entirely???
        //so the first is check point if they do have them
        //and then if points are of good quality
        if(trk->getTrackPoint(0)==nullptr||trk->getTrackPoint(1)==nullptr){
            return false;
        }
        //check if track has at least 3 of 4 RP planes used
        if(trk->getTrackPoint(0)->planesUsed()<3||trk->getTrackPoint(1)->planesUsed()<3){
            return false;
        }

        //fiducial
        double px = trk->pVec().X();
        double py = trk->pVec().Y();
        bool f1 = (0.4<abs(py)&&abs(py)<0.8);
        bool f2 = (-0.27<px);
        bool f3 = (pow(px+0.6, 2)+pow(py, 2)<1.25);
        if(!(f1&&f2&&f3)){
            return false;
        }
    }
    return true;
}



////////////////
//histogram cuts
////////////////



bool oneVertex(StRPEvent *RPEvent, StUPCEvent *UPCEvent, TH1 *hist, std::vector<int> *RPTrackIDs, std::vector<int> *UPCTrackIDs){
    //histogram
    if(hist!=nullptr){
        hist->Fill(UPCEvent->getNPrimVertices());
    }
    //cut
    if(UPCEvent->getNPrimVertices()==0){
        return false;
    }
    return true;
}

bool NTOFTracks(StRPEvent *RPEvent, StUPCEvent *UPCEvent, TH1 *hist, std::vector<int> *RPTrackIDs, std::vector<int> *UPCTrackIDs){
    //histogram
    if(hist!=nullptr){
        int nOfTOFTracks = 0;
        for(long unsigned int i = 0; i<UPCTrackIDs->size(); i++){
            int trackID = (*UPCTrackIDs)[i];
            if(UPCEvent->getTrack(trackID)->getFlag(StUPCTrack::kTof)){
                nOfTOFTracks++;
            }
        }
        return false;
    }
    //cut
    int nOfTOFTracks = 0;
    for(long unsigned int i = 0; i<UPCTrackIDs->size(); i++){
        int trackID = (*UPCTrackIDs)[i];
        if(UPCEvent->getTrack(trackID)->getFlag(StUPCTrack::kTof)){
            nOfTOFTracks++;
        } else{
            UPCTrackIDs->erase(std::find(UPCTrackIDs->begin(), UPCTrackIDs->end(), trackID));
            //because if it erases ith element, it's place gets (i+1)th, and then after i++
            //you get i+2nd, so  the i+1st is not checked
            i--;
        }
    }
    if(nOfTOFTracks<2){
        return false;
    }
    return true;
}

bool OppositeCharges(StRPEvent *RPEvent, StUPCEvent *UPCEvent, TH1 *hist, std::vector<int> *RPTrackIDs, std::vector<int> *UPCTrackIDs){
    int charge = 0;
    for(size_t i = 0; i<UPCTrackIDs->size(); i++){
        charge += UPCEvent->getTrack((*UPCTrackIDs)[i])->getCharge();
    }
    return abs(charge)!=UPCTrackIDs->size();
}

bool enoughClusters(StRPEvent *RPEvent, StUPCEvent *UPCEvent, TH1 *hist, std::vector<int> *RPTrackIDs, std::vector<int> *UPCTrackIDs){
    //TODO selekcja Patrycji
    // SELECTION: number of cluster 
    int isNumberOfTofClusterSmall = 0;
    Int_t nTofHits = UPCEvent->getNumberOfHits();
    vector <Int_t> vTray;
    vector <Int_t> vTrayUniqueVector;
    vector <Int_t> vMmodule;

    //filling with trays and modules (a pair for every hit)
    for(long unsigned int i = 0; i<nTofHits; i++){
        vTray.push_back(Int_t(UPCEvent->getHit(i)->getTray()));
        vTrayUniqueVector.push_back(Int_t(UPCEvent->getHit(i)->getTray()));
        vMmodule.push_back(Int_t(UPCEvent->getHit(i)->getModule()));
    }

    //making vTrayUniqueVector unique 
    sort(vTrayUniqueVector.begin(), vTrayUniqueVector.end());
    auto last = unique(vTrayUniqueVector.begin(), vTrayUniqueVector.end());
    vTrayUniqueVector.erase(last, vTrayUniqueVector.end());


    vector <vector <Int_t>> vModuleUniqueTray;
    //for every unique tray we add modules from hits with this unique tray
    for(long unsigned int i = 0; i<vTrayUniqueVector.size(); i++){
        vector <Int_t> vModuleUnique;
        for(long unsigned int j = 0; j<nTofHits; j++){
            if(vTrayUniqueVector[i]==vTray[j]){
                vModuleUnique.push_back(vMmodule[j]);
            }
        }
        vModuleUniqueTray.push_back(vModuleUnique);
        vModuleUnique.clear();
    }

    int totalCluster = 0;
    for(long unsigned int i = 0; i<vModuleUniqueTray.size(); i++){
        //from modules of unique trays we make unique modules of unique trays
        vector <Int_t> vec = vModuleUniqueTray[i];
        sort(vec.begin(), vec.end());
        auto last = unique(vec.begin(), vec.end());
        vec.erase(last, vec.end());

        //if there is one unique module for that unique tray for that hit we increase counter by 1
        if(vec.size()==1){
            totalCluster += 1;
        }

        //if there is more than one unique module per hit per tray we
        //
        for(long unsigned int j = 0; j<vec.size()-1; j++){
            Int_t modNum = vec[j];
            int diff = 1;
            int num = 0;
            for(long unsigned int z = j+1; z<vec.size(); z++){
                //if the module is the next module by diff
                //you increase the difference and check the next pair
                if(modNum+diff==vec[z]){
                    num += 0;
                    diff += 1;

                    if(z==j+1){
                        num += 1;
                    }
                } else if(j==(vec.size()-2)&&vec[vec.size()-2]+1!=vec[vec.size()-1]){
                    num += 2;
                    continue;
                } else{
                    num += 1;
                }
                j += (diff);
            }
            totalCluster += num;
        }
    }

    if(hist!=nullptr){
        hist->Fill(totalCluster);
        return false;
    } else if(totalCluster<=9){
        return true;
    }
    return false;
}



bool KinematicCut(StRPEvent *RPEvent, StUPCEvent *UPCEvent, TH1 *hist, std::vector<int> *RPTrackIDs, std::vector<int> *UPCTrackIDs){
    //histogram
    if(hist!=nullptr){
        for(size_t i = 0; i<UPCTrackIDs->size(); i++){
            hist->Fill(UPCEvent->getTrack((*UPCTrackIDs)[i])->getEta(), UPCEvent->getTrack((*UPCTrackIDs)[i])->getPt());
        }
    }

    //good version
    // for(long unsigned int i = 0; i<UPCTrackIDs->size(); i++){
    //     if(!(abs(UPCEvent->getTrack((*UPCTrackIDs)[i])->getEta())<0.9&&UPCEvent->getTrack((*UPCTrackIDs)[i])->getPt()>0.2)){
    //         UPCTrackIDs->erase(std::find(UPCTrackIDs->begin(), UPCTrackIDs->end(), (*UPCTrackIDs)[i]));
    //         //because if it erases ith element, it's place gets (i+1)th, and then after i++
    //         //you get i+2nd, so  the i+1st is not checked
    //         i--;
    //     }
    // }
    // if(UPCTrackIDs->size()==0){
    //     return false;
    // }
    // return true;

    //bad version
    for(long unsigned int i = 0; i<UPCTrackIDs->size(); i++){
        if(!((abs(UPCEvent->getTrack((*UPCTrackIDs)[i])->getEta())<0.9)&&(UPCEvent->getTrack((*UPCTrackIDs)[i])->getPt()>0.2))){
            return false;
        }
    }
    return true;
}

bool QualityCut(StRPEvent *RPEvent, StUPCEvent *UPCEvent, TH1 *hist, std::vector<int> *RPTrackIDs, std::vector<int> *UPCTrackIDs){
    //histogram
    if(hist!=nullptr){
        for(long unsigned int i = 0; i<UPCTrackIDs->size(); i++){
            hist->Fill(UPCEvent->getTrack((*UPCTrackIDs)[i])->getNhitsFit(), UPCEvent->getTrack((*UPCTrackIDs)[i])->getNhitsDEdx());
        }
        return false;
    }
    //bad version
    for(long unsigned int i = 0; i<UPCTrackIDs->size(); i++){
        if(!(UPCEvent->getTrack((*UPCTrackIDs)[i])->getNhitsFit()>25&&UPCEvent->getTrack((*UPCTrackIDs)[i])->getNhitsDEdx()>15)){
            return false;
        }
    }
    return true;
}

bool K0massTest(StRPEvent *RPEvent, StUPCEvent *UPCEvent, TH1 *hist, std::vector<int> *RPTrackIDs, std::vector<int> *UPCTrackIDs){
    //histogram
    if(hist!=nullptr){
        //helpful variables
        double beamValues[4];
        TVector3 vertexPrimary;
        vector<double> tempBeamVector;

        //processing
        vertexPrimary = { UPCEvent->getVertex(0)->getPosX(), UPCEvent->getVertex(0)->getPosY(), UPCEvent->getVertex(0)->getPosZ() };
        tempBeamVector = FindPosition(UPCEvent->getFillNumber(), vertexPrimary.Z(), beamData[0], beamData[1], beamData[2], beamData[3], beamData[4], beamData[5], beamData[6], beamData[7], beamData[8]);
        beamValues[0] = tempBeamVector[0];
        beamValues[1] = tempBeamVector[1];
        beamValues[2] = tempBeamVector[2];
        beamValues[3] = tempBeamVector[3];
        int first_K0_pion = -1;
        int second_K0_pion = -1;
        int first_vertex_pion = -1;
        int second_vertex_pion = -1;
        StUPCTrack *first_track;
        StUPCTrack *second_track;
        StUPCV0 *tempParticle;
        StUPCV0 *K0_pair;
        StUPCV0 *vertex_pair;
        //actual loop
        for(long unsigned int i = 0; i<UPCTrackIDs->size()-1; i++){
            for(long unsigned int j = i+1; j<UPCTrackIDs->size(); j++){
                first_track = UPCEvent->getTrack((*UPCTrackIDs)[i]);
                second_track = UPCEvent->getTrack((*UPCTrackIDs)[j]);
                tempParticle = new StUPCV0(first_track, second_track, particleMass[0], particleMass[0], 1, 1, vertexPrimary, beamValues, UPCEvent->getMagneticField(), true);
                //tests if accept the particle
                bool K0test1 = first_track->getCharge()*second_track->getCharge()<0;
                bool K0test2 = tempParticle->dcaDaughters()<1.5;
                bool K0test3 = tempParticle->pointingAngleHypo()>0.925;
                bool K0test4 = tempParticle->DCABeamLine()<1.5;
                bool K0test5 = tempParticle->m()>kaonMassWindowPresentationLow&&tempParticle->m()<kaonMassWindowPresentationHigh;
                if(!(K0test1&&K0test2&&K0test3&&K0test4&&K0test5)){
                    continue;
                }
                //filling
                if(first_K0_pion<0){
                    first_K0_pion = i;
                    second_K0_pion = j;
                    K0_pair = new StUPCV0(first_track, second_track, particleMass[0], particleMass[0], 1, 1, vertexPrimary, beamValues, UPCEvent->getMagneticField(), true);
                } else if(abs(tempParticle->m()-K0mass)<abs(K0_pair->m()-K0mass)){
                    first_K0_pion = i;
                    second_K0_pion = j;
                    delete K0_pair;
                    K0_pair = new StUPCV0(first_track, second_track, particleMass[0], particleMass[0], 1, 1, vertexPrimary, beamValues, UPCEvent->getMagneticField(), true);
                }
                //finishing
                delete tempParticle;
            }
        }
        for(long unsigned int i = 0; i<UPCTrackIDs->size()-1; i++){
            if(i!=first_K0_pion&&i!=second_K0_pion&&first_vertex_pion>0){
                second_vertex_pion = i;
            } else if(i!=first_K0_pion&&i!=second_K0_pion){
                first_vertex_pion = i;
            }
        }
        vertex_pair = new StUPCV0(UPCEvent->getTrack((*UPCTrackIDs)[first_vertex_pion]), UPCEvent->getTrack((*UPCTrackIDs)[second_vertex_pion]), particleMass[0], particleMass[0], 1, 1, vertexPrimary, beamValues, UPCEvent->getMagneticField(), true);
        //tests to see if vertex suggestion is okay
        bool PVtest1 = vertex_pair->dcaDaughters()<1.5;
        bool PVtest2 = vertex_pair->DCABeamLine()<1.5;
        if(PVtest1&&PVtest2&&first_K0_pion>0){
            hist->Fill(K0_pair->m());
        }
    }
    //cut
    return true;
}

bool LambdamassTest(StRPEvent *RPEvent, StUPCEvent *UPCEvent, TH1 *hist, std::vector<int> *RPTrackIDs, std::vector<int> *UPCTrackIDs){
    //histogram
    if(hist!=nullptr){
        //helpful variables
        double beamValues[4];
        TVector3 vertexPrimary;
        vector<double> tempBeamVector;

        //processing
        vertexPrimary = { UPCEvent->getVertex(0)->getPosX(), UPCEvent->getVertex(0)->getPosY(), UPCEvent->getVertex(0)->getPosZ() };
        tempBeamVector = FindPosition(UPCEvent->getFillNumber(), vertexPrimary.Z(), beamData[0], beamData[1], beamData[2], beamData[3], beamData[4], beamData[5], beamData[6], beamData[7], beamData[8]);
        beamValues[0] = tempBeamVector[0];
        beamValues[1] = tempBeamVector[1];
        beamValues[2] = tempBeamVector[2];
        beamValues[3] = tempBeamVector[3];
        int first_K0_pion = -1;
        int second_K0_pion = -1;
        int first_vertex_pion = -1;
        int second_vertex_pion = -1;
        StUPCTrack *first_track;
        StUPCTrack *second_track;
        StUPCV0 *tempParticle;
        StUPCV0 *K0_pair;
        StUPCV0 *vertex_pair;
        //actual loop
        for(long unsigned int i = 0; i<UPCTrackIDs->size()-1; i++){
            for(long unsigned int j = i+1; j<UPCTrackIDs->size(); j++){
                first_track = UPCEvent->getTrack((*UPCTrackIDs)[i]);
                second_track = UPCEvent->getTrack((*UPCTrackIDs)[j]);
                tempParticle = new StUPCV0(first_track, second_track, particleMass[0], particleMass[0], 1, 1, vertexPrimary, beamValues, UPCEvent->getMagneticField(), true);
                //tests if accept the particle
                bool K0test1 = first_track->getCharge()*second_track->getCharge()<0;
                bool K0test2 = tempParticle->dcaDaughters()<1.5;
                bool K0test3 = tempParticle->pointingAngleHypo()>0.925;
                bool K0test4 = tempParticle->DCABeamLine()<1.5;
                bool K0test5 = tempParticle->m()>kaonMassWindowPresentationLow&&tempParticle->m()<kaonMassWindowPresentationHigh;
                if(!(K0test1&&K0test2&&K0test3&&K0test4&&K0test5)){
                    continue;
                }
                //filling
                if(first_K0_pion<0){
                    first_K0_pion = i;
                    second_K0_pion = j;
                    K0_pair = new StUPCV0(first_track, second_track, particleMass[2], particleMass[0], 1, 1, vertexPrimary, beamValues, UPCEvent->getMagneticField(), true);
                } else if(abs(tempParticle->m()-Lambdamass)<abs(K0_pair->m()-Lambdamass)){
                    first_K0_pion = i;
                    second_K0_pion = j;
                    delete K0_pair;
                    K0_pair = new StUPCV0(first_track, second_track, particleMass[2], particleMass[0], 1, 1, vertexPrimary, beamValues, UPCEvent->getMagneticField(), true);
                }
                //finishing
                delete tempParticle;
            }
        }
        for(long unsigned int i = 0; i<UPCTrackIDs->size()-1; i++){
            if(i!=first_K0_pion&&i!=second_K0_pion&&first_vertex_pion>0){
                second_vertex_pion = i;
            } else if(i!=first_K0_pion&&i!=second_K0_pion){
                first_vertex_pion = i;
            }
        }
        vertex_pair = new StUPCV0(UPCEvent->getTrack((*UPCTrackIDs)[first_vertex_pion]), UPCEvent->getTrack((*UPCTrackIDs)[second_vertex_pion]), particleMass[0], particleMass[0], 1, 1, vertexPrimary, beamValues, UPCEvent->getMagneticField(), true);
        //tests to see if vertex suggestion is okay
        bool PVtest1 = vertex_pair->dcaDaughters()<1.5;
        bool PVtest2 = vertex_pair->DCABeamLine()<1.5;
        if(PVtest1&&PVtest2&&first_K0_pion>0){
            hist->Fill(K0_pair->m());
        }
    }
    //cut
    return true;
}





//just in case

bool CounterCut(StRPEvent *RPEvent, StUPCEvent *UPCEvent, TH1 *hist, std::vector<int> *RPTrackIDs, std::vector<int> *UPCTrackIDs){
    //histogram
    if(hist!=nullptr){
        hist->Fill(0);
    }
    //cut
    return true;
}

// bool a(StRPEvent *RPEvent, StUPCEvent *UPCEvent, std::vector<int> *RPTrackIDs, std::vector<int> *UPCTrackIDs){

// }