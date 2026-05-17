class StChain;
class StMuDstMaker;
class StMuDst;
class StMuEvent;
class StMuBTofPidTraits;


// for test root4star -b -q -l 'preselectionTest.C(10)'. //for 10 events to test

void preselectionMuDst(int nEvents = -1) {  // -1 : run for all events

    gROOT->Macro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
    StMuDebug::setLevel(0);  // too many STAR verbose so removing it

    // input and out
    const char *inFile  = "./centralDiffractive.MuDst.root";
    const char *outFile = "./centralDiffractivePreselected.MuDst.root";

    
    TFile *ftmp = TFile::Open(inFile);
    TTree *ttmp = (TTree*)ftmp->Get("MuDst");
    Long64_t totalEvents = ttmp->GetEntries();
    ftmp->Close(); delete ftmp;

    if (nEvents < 0 || nEvents > totalEvents) nEvents = totalEvents;
    printf("[INFO] Input:  %s\n",   inFile);
    printf("[INFO] Output: %s\n",   outFile);
    printf("[INFO] Processing %d / %lld events\n\n", nEvents, totalEvents);

    
    // PASS 1

    StChain      *chain   = new StChain("preselChain");
    StMuDstMaker *inMaker = new StMuDstMaker(0, 0, "", inFile, "MuDst", 1000, "inMaker");
    inMaker->SetStatus("*", 1);
    chain->Init();

    TEntryList *elist = new TEntryList("passing", "Events passing preselection");

    int processed = 0, selected = 0;
    for (int i = 0; i < nEvents; i++) {
        if (chain->Make() != kStOK) break;
        processed++;

        StMuDst   *mMuDst = inMaker->muDst();    if (!mMuDst) { chain->Clear(); continue; }
        StMuEvent *ev     = mMuDst->event();     if (!ev)     { chain->Clear(); continue; }

       
        // YOUR PRESELECTION CUT HERE
        
        bool pass = false;
        bool oneTrackAreadyFound = false;
        //vertex test
        if(mMuDst->numberOfPrimaryVertices()==1){
            //static call to set current primary vertex
            StMuDst::setVertexIndex(0);

            //get array of primary tracks
            TObjArray* trkArray = mMuDst->primaryTracks();
            if(!trkArray) continue;

            //tracks loop
            for(Int_t itrk = 0; itrk<trkArray->GetEntriesFast(); itrk++){
                StMuTrack* track = dynamic_cast<StMuTrack*>(trkArray->At(itrk));
                if(!track) continue;

                //TOF matching
                const StMuBTofPidTraits& tofPid = track->btofPidTraits();
                if(tofPid.matchFlag()!=0&&track->pt()>0.2&&fabs(track->eta())<0.9){
                    if(oneTrackAreadyFound){
                        pass=true;
                        break;
                    }else{
                        oneTrackAreadyFound=true;
                    }
                }
            }//tracks loop
        }
        

        if (pass) {
            elist->Enter(i);
            selected++;
        }

        if (processed % 1000 == 0)
            printf("  Processed %d / %d  (selected so far: %d)\n",
                   processed, nEvents, selected);

        chain->Clear();
    }

    printf("\nPass 1 done: processed=%d  selected=%d\n\n", processed, selected);
    chain->Finish();
    delete chain;

    if (selected == 0) {
        printf("No events passed the cut. No output written.\n");
        gSystem->Exit(0);
    }

    // PASS 2

    TFile *fin  = TFile::Open(inFile,  "READ");
    TFile *fout = TFile::Open(outFile, "RECREATE");

    TTree *tin = (TTree*)fin->Get("MuDst");
    tin->SetEntryList(elist);

    fout->cd();
    gErrorIgnoreLevel = kWarning;  //silencing TTree mode verbose since no data has been lost both files are created withing same schema
    TTree *tout = tin->CopyTree("");
    printf("Wrote %lld events to %s  (%.1f KB)\n",
           tout->GetEntries(), outFile, fout->GetSize()/1024.0);

    fout->Write();
    fout->Close();
    fin->Close();
    delete fout;
    delete fin;

    //return the number of entries copied so that the loop will know if it should continue
    gSystem->Exit(tout->GetEntries());
}