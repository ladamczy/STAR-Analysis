#include <iostream>

//ROOT headers
#include <ROOT/TThreadedObject.hxx>
#include <TTreeReader.h>
#include <ROOT/TTreeProcessorMT.hxx>

//STAR headers
#include <StarGenEvent.h>
#include <StarGenParticle.h>

//my headers
#include <ProcessingInsideLoop.h>
#include <ProcessingOutsideLoop.h>
#include <UsefulThings.h>

int main(int argc, char* argv[]){

    //argv[1] - input file
    //argv[2] - output folder
    //argv[3] - #of cores

    int nthreads = 1;
    if(argc==4){
        nthreads = atoi(argv[3]);
    }

    cout<<"Program is running on "<<nthreads<<" threads"<<endl;
    ROOT::EnableThreadSafety();
    //actually i'm not sure if it's needed here
    // ROOT::EnableImplicitMT(nthreads); //turn on multicore processing

    //preparing input & output
    TChain* eventFiles = new TChain("genevents");
    if(ConnectInput(argc, argv, eventFiles)){
        cout<<"All files connected"<<endl;
    }
    const string& outputFolder = argv[2];

    //Useful IDs
    const int K0sPDGid = 310;
    const int K0sbarPDGid = -310;
    const int LambdaPDGid = 3122;
    const int LambdabarPDGid = -3122;
    const int phiPDGid = 333;
    const int phibarPDGid = -333;
    const int KstarPDGid = 313;
    const int KstarbarPDGid = -313;
    const int piplusPDGid = 211;
    const int piminusPDGid = -211;
    const int KplusPDGid = 321;
    const int KminusPDGid = -321;
    const int pplusPDGid = 2212;
    const int pminusPDGid = -2212;

    //histograms
    ProcessingOutsideLoop outsideprocessing;
    //mixing TOF between events proved to be a failure
    //mass histograms with TOF first and mixing event pairs after
    outsideprocessing.AddHistogram(TH1D("Name", "Name of simulated particles;id;Number of particles", 1, 0, 1));
    outsideprocessing.AddHistogram(TH2D("etapTK0S", "K^{0}_{S} number;eta;p_{T}", 10, -1, 1, 10, 0, 2.5));
    outsideprocessing.AddHistogram(TH2D("etapTLambda", "#Lambda^{0} number;eta;p_{T}", 10, -1, 1, 10, 0, 2.5));
    outsideprocessing.AddHistogram(TH2D("etapTLambdabar", "#bar{#Lambda}^{0} number;eta;p_{T}", 10, -1, 1, 10, 0, 2.5));
    outsideprocessing.AddHistogram(TH2D("etapTKstar", "K^{*}(892) number;eta;p_{T}", 10, -1, 1, 10, 0, 2.5));
    outsideprocessing.AddHistogram(TH2D("etapTKstarbar", "#bar{K}^{*}(892) number;eta;p_{T}", 10, -1, 1, 10, 0, 2.5));
    outsideprocessing.AddHistogram(TH2D("etapTphi", "#varphi(1020) number;eta;p_{T}", 10, -1, 1, 10, 0, 2.5));

    //processing
    //defining TreeProcessor
    ROOT::TTreeProcessorMT TreeProc(*eventFiles, nthreads);

    //defining processing function
    auto myFunction = [&](TTreeReader& myReader){
        //getting values from TChain, in-loop histogram initialization
        TTreeReaderValue<StarGenEvent> tempStarGenEventInstance(myReader, "primaryEvent");
        ProcessingInsideLoop insideprocessing;
        StarGenEvent* tempStarGenEventpointer;
        insideprocessing.GetLocalHistograms(&outsideprocessing);

        while(myReader.Next()){
            tempStarGenEventpointer = tempStarGenEventInstance.Get();
            //below is the loop

            //checking created particles - how many of them existed
            for(size_t particle_index = 0; particle_index<tempStarGenEventpointer->GetNumberOfParticles(); particle_index++){
                //putting result (detection or not) into histogram
                TLorentzVector particle;
                particle.SetXYZM((*tempStarGenEventpointer)[particle_index]->GetPx(), (*tempStarGenEventpointer)[particle_index]->GetPy(), (*tempStarGenEventpointer)[particle_index]->GetPz(), (*tempStarGenEventpointer)[particle_index]->GetMass());
                switch((*tempStarGenEventpointer)[particle_index]->GetId()){
                case K0sPDGid:
                case K0sbarPDGid:
                    insideprocessing.Fill("Name", "K^{0}_{S}", 1.);
                    insideprocessing.Fill("etapTK0S", particle.Eta(), particle.Pt());
                    break;
                case LambdaPDGid:
                    insideprocessing.Fill("Name", "#Lambda^{0}", 1.);
                    insideprocessing.Fill("etapTLambda", particle.Eta(), particle.Pt());
                    break;
                case LambdabarPDGid:
                    insideprocessing.Fill("Name", "#bar{#Lambda}^{0}", 1.);
                    insideprocessing.Fill("etapTLambdabar", particle.Eta(), particle.Pt());
                    break;
                case KstarPDGid:
                    insideprocessing.Fill("Name", "K^{*}", 1.);
                    insideprocessing.Fill("etapTKstar", particle.Eta(), particle.Pt());
                    break;
                case KstarbarPDGid:
                    insideprocessing.Fill("Name", "#bar{K}^{*}", 1.);
                    insideprocessing.Fill("etapTKstarbar", particle.Eta(), particle.Pt());
                    break;
                case phiPDGid:
                case phibarPDGid:
                    insideprocessing.Fill("Name", "#phi", 1.);
                    insideprocessing.Fill("etapTphi", particle.Eta(), particle.Pt());
                    break;
                default:
                    break;
                }
            }

            //event loop finish
        }

        //lambda finish
        return 0;
    };

    TreeProc.Process(myFunction);

    outsideprocessing.Merge();
    outsideprocessing.GetPointerAfterMerge1D("Name")->LabelsDeflate();
    outsideprocessing.GetPointerAfterMerge1D("Name")->SetMinimum(0);
    outsideprocessing.GetPointerAfterMerge1D("Name")->LabelsOption("a", "X");

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
    outsideprocessing.SaveToFile(outputFileHist);
    outputFileHist->Close();

    return 0;
}
