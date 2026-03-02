//dummy class declarations so that I can actually load them later
class StMuMcVertex;
class StMuMcTrack;
class StarGenEvent;
class StarGenParticle;

bool stringEndsWith(const char* str, const char* suffix){
    if(!str||!suffix)
        return 0;
    size_t lenstr = strlen(str);
    size_t lensuffix = strlen(suffix);
    if(lensuffix>lenstr)
        return 0;
    return strncmp(str+lenstr-lensuffix, suffix, lensuffix)==0;
}

void copyStMuMcVertex(struct g2t_vertex_st* copy, StMuMcVertex* original){
    copy->id = original->Id();
    copy->n_daughter = original->NoDaughters();
    copy->parent_p = original->IdParTrk();
    copy->is_itrmd = original->IsIntermedate();
    copy->ge_tof = original->Time();
    copy->ge_x[0] = original->XyzV().x();
    copy->ge_x[1] = original->XyzV().y();
    copy->ge_x[2] = original->XyzV().z();
}

void copyStMuMcTrack(struct g2t_track_st* copy, StMuMcTrack* original){
    copy->ge_pid = original->GePid();
    copy->id = original->Id();
    copy->is_shower = original->IsShower();
    copy->itrmd_vertex_p = original->ItrmdVertex();
    copy->start_vertex_p = original->IdVx();
    copy->stop_vertex_p = original->IdVxEnd();
    copy->charge = original->Charge();
    copy->e = original->E();
    copy->eta = original->Eta();
    copy->p[0] = original->Pxyz().x();
    copy->p[1] = original->Pxyz().y();
    copy->p[2] = original->Pxyz().z();
    copy->pt = original->pT();
    copy->ptot = original->Ptot();
    copy->rapidity = original->Rapidity();
    copy->n_ctb_hit = original->No_ctb_hit();  /* Nhits in ctb */
    copy->n_eem_hit = original->No_eem_hit();  /* Nhits in eem (endcap em cal) */
    copy->n_emc_hit = original->No_emc_hit();  /* Nhits in emc */
    copy->n_esm_hit = original->No_esm_hit();  /* Nhits in esm (endcap shower max) */
    copy->n_ftp_hit = original->No_ftp_hit();  /* Nhits in forward tpc */
    copy->n_gem_hit = original->No_gem_hit();  /* Nhits in gem barrel */
    copy->n_hpd_hit = original->No_hpd_hit();  /* Nhits in hpd */
    copy->n_ist_hit = original->No_ist_hit();  /* Nhits in ist */
    copy->n_igt_hit = original->No_igt_hit();  /* Nhits in igt */
    copy->n_fst_hit = original->No_fst_hit();  /* Nhits in fst */
    copy->n_fgt_hit = original->No_fgt_hit();  /* Nhits in fgt */
    copy->n_fpd_hit = original->No_fpd_hit();  /* Nhits in fpd */
    copy->n_mwc_hit = original->No_mwc_hit();  /* Nhits in mwc */
    copy->n_pgc_hit = original->No_pgc_hit();  /* Nhits in pgc  ???  */
    copy->n_pmd_hit = original->No_pmd_hit();  /* Nhits in pmd (PMD) */
    copy->n_smd_hit = original->No_smd_hit();  /* number of hits in shower max */
    copy->n_ssd_hit = original->No_ssd_hit();  /* Nhits in ssd */
    copy->n_svt_hit = original->No_svt_hit();  /* Nhits in svt */
    copy->n_pix_hit = original->No_pix_hit();  /* Nhits in pix */
    copy->n_tof_hit = original->No_tof_hit();  /* Nhits in tof */
    copy->n_tpc_hit = original->No_tpc_hit();  /* Nhits in tpc */
    copy->n_vpd_hit = original->No_vpd_hit();  /* Nhits in vpd */
    copy->n_etr_hit = original->No_etr_hit();  /* Nhits in etr */
    copy->n_hca_hit = original->No_hca_hit();  /* Nhits in hca */
    copy->n_fts_hit = original->No_fts_hit();  /* Nhits in fts */
    copy->n_eto_hit = original->No_eto_hit();  /* Nhits in eto */
    copy->n_stg_hit = original->No_stg_hit();  /* Nhits in stg */
    copy->n_wca_hit = original->No_wca_hit();  /* Nhits in wca */
    copy->n_pre_hit = original->No_pre_hit();  /* Nhits in pre */
    copy->n_epd_hit = original->No_epd_hit();  /* Nhits in epd */
}

void createTrackFromPythia(struct g2t_track_st* copy, StarGenParticle* original){
    //should be set outside because it depends on particle position in the event array
    copy->id = 0;
    copy->stop_vertex_p = 0;

    //copyable/can be set now

    copy->ge_pid = original->GetId();

    TDatabasePDG* pdgdat = TDatabasePDG::Instance();
    //charge in units of |e|/3, needs conversion
    copy->charge = TMath::Nint(pdgdat->GetParticle(original->GetId())->Charge()/3.);
    delete pdgdat;

    TLorentzVector fourmom = original->momentum();
    copy->e = fourmom.E();
    copy->p[0] = fourmom.X();
    copy->p[1] = fourmom.Y();
    copy->p[2] = fourmom.Z();
    copy->ptot = fourmom.P();

    //set to predetermined values, the same as other particles
    copy->start_vertex_p = 1;
    copy->is_shower = 0;
    copy->itrmd_vertex_p = 0;
    copy->eta = -999.;
    copy->pt = -999.;
    copy->rapidity = -999.;
    copy->n_ctb_hit = 0;  /* Nhits in ctb */
    copy->n_eem_hit = 0;  /* Nhits in eem (endcap em cal) */
    copy->n_emc_hit = 0;  /* Nhits in emc */
    copy->n_esm_hit = 0;  /* Nhits in esm (endcap shower max) */
    copy->n_ftp_hit = 0;  /* Nhits in forward tpc */
    copy->n_gem_hit = 0;  /* Nhits in gem barrel */
    copy->n_hpd_hit = 0;  /* Nhits in hpd */
    copy->n_ist_hit = 0;  /* Nhits in ist */
    copy->n_igt_hit = 0;  /* Nhits in igt */
    copy->n_fst_hit = 0;  /* Nhits in fst */
    copy->n_fgt_hit = 0;  /* Nhits in fgt */
    copy->n_fpd_hit = 0;  /* Nhits in fpd */
    copy->n_mwc_hit = 0;  /* Nhits in mwc */
    copy->n_pgc_hit = 0;  /* Nhits in pgc  ???  */
    copy->n_pmd_hit = 0;  /* Nhits in pmd (PMD) */
    copy->n_smd_hit = 0;  /* number of hits in shower max */
    copy->n_ssd_hit = 0;  /* Nhits in ssd */
    copy->n_svt_hit = 0;  /* Nhits in svt */
    copy->n_pix_hit = 0;  /* Nhits in pix */
    copy->n_tof_hit = 0;  /* Nhits in tof */
    copy->n_tpc_hit = 0;  /* Nhits in tpc */
    copy->n_vpd_hit = 0;  /* Nhits in vpd */
    copy->n_etr_hit = 0;  /* Nhits in etr */
    copy->n_hca_hit = 0;  /* Nhits in hca */
    copy->n_fts_hit = 0;  /* Nhits in fts */
    copy->n_eto_hit = 0;  /* Nhits in eto */
    copy->n_stg_hit = 0;  /* Nhits in stg */
    copy->n_wca_hit = 0;  /* Nhits in wca */
    copy->n_pre_hit = 0;  /* Nhits in pre */
    copy->n_epd_hit = 0;  /* Nhits in epd */
}

bool isParticleAddable(StarGenEvent* pythia_event, int i){
    //333 = phi(1020)
    //313 = K*(892)
    switch(((*pythia_event)[i])->GetId()){
    case 333:
    case 313:
    case -313:
        return true;
        break;

    default:
        return false;
        break;
    }
}

int addPythiaParticles(char const* firstFile, char const* secondFile, int run_number = 0){
    //load all the libraries needed
    gROOT->Macro("loadMuDst.C");
    gSystem->Load("libStarGeneratorUtil.so");
    gSystem->Load("libStarGeneratorEvent.so");


    //opening MuDst file and Pythia files (with correct extensions)
    TFile* MuFile, * PythiaFile;
    if(stringEndsWith(firstFile, ".MuDst.root")&&stringEndsWith(secondFile, ".root")){
        MuFile = new TFile(firstFile, "update");
        PythiaFile = new TFile(secondFile, "read");
    } else if(stringEndsWith(firstFile, ".root")&&stringEndsWith(secondFile, ".MuDst.root")){
        MuFile = new TFile(secondFile, "update");
        PythiaFile = new TFile(firstFile, "read");
    } else{
        printf("There is no files with suitable extensions\n");
        return -1;
    }
    //checking if they contain correct trees
    TTree* MuTree, * PythiaTree;
    if(MuFile->Get("MuDst")==0){
        printf("MuDst TTree not found\n");
        return -1;
    } else{
        MuTree = (TTree*)MuFile->Get("MuDst");
    }
    if(PythiaFile->Get("genevents")==0){
        printf("genevents TTree not found\n");
        return -1;
    } else{
        PythiaTree = (TTree*)PythiaFile->Get("genevents");
    }
    //making them friends to iterate over event number as an index
    //normally not needed, but might be when some filters between Pythia and Geant cut number of events
    //event number copied to upcDst is in MuEvent.mEventInfo.mId branch
    //event taken from Pythia is in mEventNumber, with mFilterResult=1 when passed filter and mFilterResult=2 when not
    PythiaTree->BuildIndex("mEventNumber", "mFilterResult");


    //TClonesArrays and other things needed for MuDst file manipulation
    TClonesArray* MuEvent;
    TClonesArray* MuMcVtx;
    TClonesArray* MuMcTracks;
    StarGenEvent* PythiaEvent;
    //branch connecting
    MuTree->SetBranchAddress("MuEvent", &MuEvent);
    MuTree->SetBranchAddress("StMuMcVertex", &MuMcVtx);
    MuTree->SetBranchAddress("StMuMcTrack", &MuMcTracks);
    PythiaTree->SetBranchAddress("primaryEvent", &PythiaEvent);
    //creating a new tree for copying the events (it has already all the branches connected and all)
    //(0 is there because we want no entries copied)
    gROOT->cd(0); // tells ROOT the next TTree is "memory-resident" (whatever that means), otherwise we get problems with writing
    TTree* NewMuTree = MuTree->CloneTree(0);
    NewMuTree->SetName("MuDstTemp");


    //loop for copying and inserting the phi and K* particles
    for(size_t entry = 0; entry<MuTree->GetEntries(); entry++){
        //setting trees to the same event
        MuTree->GetEntry(entry);
        StMuEvent* MuEventInstance = (StMuEvent*)MuEvent->UncheckedAt(0);
        PythiaTree->GetEntryWithIndex(MuEventInstance->eventNumber(), 1);

        //setting proper run number (default is 18091010)
        if(run_number==0){
            run_number = 18091010;
        }
        StEventInfo& eventInfo = MuEventInstance->eventInfo();
        eventInfo.setRunId(run_number);

        //test printout
        if(entry==0){
            MuMcVtx->Print();
            MuMcTracks->Print();
            PythiaEvent->Print();
        }


        //loop for catching add-able particles
        for(Int_t particle = 0; particle<PythiaEvent->GetNumberOfParticles(); particle++){
            if(isParticleAddable(PythiaEvent, particle)){
                //REMINDER: 
                //PARTICLES AND VERTICES AFTER GEANT ARE INDEXED FROM 1 TO n
                //SO THERE IS A NEED TO +-1 IN CERTAIN PLACES
                //BUT THOSE IN PYTHIA ARE NOT


                //adding a vertex (in a weird way, since "AddLast" and similar methods are forbidden in TClonesArray)
                //filling only members of the temporary struct used later in upcDst conversion, and NoTk
                struct g2t_vertex_st dummy_vertex;
                dummy_vertex.id = MuMcVtx->GetEntries()+1;
                dummy_vertex.n_daughter = ((*PythiaEvent)[particle])->GetLastDaughter()-((*PythiaEvent)[particle])->GetFirstDaughter()+1;
                //we're getting decay vertex of the particle by checking the production vertex of its daughters
                //and the time of the vertex by t=s/v
                //vertex will probably be in almost the same position with mean life of 10^-20 s
                //but no one will tell me i'm cutting corners
                int daughterIndex = ((*PythiaEvent)[particle])->GetFirstDaughter();
                StarGenParticle* daughter = (*PythiaEvent)[daughterIndex];
                //vertex position
                dummy_vertex.ge_x[0] = daughter->GetVx();
                dummy_vertex.ge_x[1] = daughter->GetVy();
                dummy_vertex.ge_x[2] = daughter->GetVz();
                //vertex time
                TVector3 productionVertex(((*PythiaEvent)[particle])->GetVx(), ((*PythiaEvent)[particle])->GetVy(), ((*PythiaEvent)[particle])->GetVz());
                TVector3 decayVertex(daughter->GetVx(), daughter->GetVy(), daughter->GetVz());
                dummy_vertex.ge_tof = (decayVertex-productionVertex).Mag()/((*PythiaEvent)[particle])->momentum().BoostVector().Mag()/299.792458; //mm/c to ns conversion
                //creating a new vertex
                MuMcVtx->ExpandCreate(MuMcVtx->GetEntries()+1);
                new(MuMcVtx->At(MuMcVtx->GetEntries()-1)) StMuMcVertex(dummy_vertex);


                //changing production vertices of tracks from 1 to the freshly added one
                for(Int_t decayPythiaIndex = ((*PythiaEvent)[particle])->GetFirstDaughter(); decayPythiaIndex<=((*PythiaEvent)[particle])->GetLastDaughter(); decayPythiaIndex++){
                    //getting Geant array number = mStack-1
                    //NOT Stack number!!!
                    int decayGeantIndex = ((*PythiaEvent)[decayPythiaIndex])->GetStack()-1;
                    //switch made specifically to counter K0/K0_bar (not passed to geant) turning into K0S/K0L (which is passed into Geant)
                    //we just take the particle that K0/K0_bar "decays" into
                    if(fabs(((*PythiaEvent)[decayPythiaIndex])->GetId())==311){
                        int K0DaughterIndex = ((*PythiaEvent)[decayPythiaIndex])->GetFirstDaughter();
                        decayGeantIndex = ((*PythiaEvent)[K0DaughterIndex])->GetStack()-1;
                    }
                    struct g2t_track_st dummy_track;
                    //modifying positive one (had to cast from TObject*, because ROOT is stupid like that)
                    copyStMuMcTrack(&dummy_track, static_cast<StMuMcTrack*>(MuMcTracks->At(decayGeantIndex)));
                    dummy_track.start_vertex_p = MuMcVtx->GetEntries();//vertices are numbered from 1 to n, unlike arrays, and we take the last one
                    new(MuMcTracks->At(decayGeantIndex)) StMuMcTrack(dummy_track);
                }


                //adding the decayed track from pythia in the last place                
                struct g2t_track_st dummy_decaying_track;
                createTrackFromPythia(&dummy_decaying_track, (*PythiaEvent)[particle]);
                dummy_decaying_track.stop_vertex_p = MuMcVtx->GetEntries();//vertices are numbered from 1 to n, unlike arrays, and we take the last one, freshly created
                dummy_decaying_track.id = MuMcTracks->GetEntries()+1;
                MuMcTracks->ExpandCreate(MuMcTracks->GetEntries()+1);
                new(MuMcTracks->At(MuMcTracks->GetEntries()-1)) StMuMcTrack(dummy_decaying_track);


                //removing correct number of tracks from the first vertex (gains 1, loses however many particle decayed into)
                //so the difference is equal to GetLastDaughter()-GetFirstDaughter()
                struct g2t_vertex_st dummy_first_vertex;
                copyStMuMcVertex(&dummy_first_vertex, static_cast<StMuMcVertex*>(MuMcVtx->At(0)));
                dummy_first_vertex.n_daughter -= ((*PythiaEvent)[particle])->GetLastDaughter()-((*PythiaEvent)[particle])->GetFirstDaughter();
                new(MuMcVtx->At(0)) StMuMcVertex(dummy_first_vertex);
            }
        }

        if(entry==0){
            MuMcVtx->Print();
            MuMcTracks->Print();
        }

        //writing to a tree part
        NewMuTree->Fill();
    }

    //saving the new tree and deleting the old one
    //file grows in side because it needs to hold the second TTree, and there is no shrinking of thr .root files
    //I can mitigate that by making a brand new file I guess
    MuFile->cd();
    NewMuTree->AutoSave();
    gDirectory->Delete("MuDst;*");
    NewMuTree->SetName("MuDst");
    NewMuTree->Write();
    gDirectory->Delete("MuDstTemp;*");
    MuFile->Write();

    return 0;
}