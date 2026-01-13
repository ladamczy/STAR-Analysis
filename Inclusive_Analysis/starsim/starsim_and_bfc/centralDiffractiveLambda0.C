// macro to instantiate the Geant3 from within
// STAR  C++  framework and get the starsim prompt
// To use it do
//  root4star starsim.C

class St_geant_Maker;
St_geant_Maker* geant_maker = 0;

class StarGenEvent;
StarGenEvent* event = 0;

class StarPrimaryMaker;
StarPrimaryMaker* _primary = 0;

class StarFilterMaker;
StarFilterMaker* filter = 0;

std::map<int, int> date;
std::map<int, int> time;

// ----------------------------------------------------------------------------
bool loadTimestampFile(){
    ifstream fp;
    fp.open("./StRoot/run17_date_time.txt");
    if(!fp.is_open()){
        printf("Run timestamp file could not be loaded\n");
        return false;
    }
    int runNumb;
    int DateSt;
    int timeSt;
    int rowNumber = 0;
    while(fp.good()){
        std::string line;
        fp>>DateSt>>timeSt>>runNumb;
        date[runNumb] = DateSt;
        time[runNumb] = timeSt;
        rowNumber++;
    }
    printf("Imported %d rows\n", rowNumber);

    return true;
}
// ----------------------------------------------------------------------------
bool getBeamlineParameters(int runNumber, double* beamline){
    loadTimestampFile();
    //checking if run actually exists on a list
    if(date.find(runNumber)==date.end()||time.find(runNumber)==time.end()){
        printf("Run does not exist in the list, possibly it isn't from Run17?\n");
        return false;
    }
    //database maker setup
    printf("Database setup start; run %d\n", runNumber);
    St_db_Maker* dbMk = new St_db_Maker("dbMaker", "MySQL:StarDb", "$STAR/StarDb");
    dbMk->SetDebug();
    int dt = date.find(runNumber)->second;
    int tt = time.find(runNumber)->second;
    dbMk->SetDateTime(dt, tt);
    dbMk->SetFlavor("ofl");
    dbMk->Init();
    dbMk->Make();
    printf("Database setup finished\n");

    TDataSet* DB = 0;
    DB = dbMk->GetDataBase("Calibrations/rhic/vertexSeed");
    if(!DB){
        printf("ERROR: no table found in db, or malformed local db config\n");
        return false;
    }

    St_vertexSeed* dataset = 0;
    dataset = (St_vertexSeed*)DB->Find("vertexSeed");
    Int_t rows = dataset->GetNRows();
    if(rows>1)
        printf("INFO: found INDEXED table with "<<rows<<" rows\n");

    if(dataset){
        vertexSeed_st* table = dataset->GetTable();
        beamline[0] = table->x0;
        beamline[1] = table->y0;
        beamline[2] = table->dxdz;
        beamline[3] = table->dydz;
    } else{
        printf("ERROR: dataset does not contain requested table\n");
        return false;
    }

    printf("Database set up successfully\n");
    return true;
}
// ----------------------------------------------------------------------------
void geometry(TString tag, Bool_t agml = true){
    TString cmd = "DETP GEOM "; cmd += tag;
    if(!geant_maker) geant_maker = (St_geant_Maker*)chain->GetMaker("geant");
    geant_maker->LoadGeometry(cmd);
    //  if ( agml ) command("gexec $STAR_LIB/libxgeometry.so");
}
// ----------------------------------------------------------------------------
void command(TString cmd){
    if(!geant_maker) geant_maker = (St_geant_Maker*)chain->GetMaker("geant");
    geant_maker->Do(cmd);
}
// ----------------------------------------------------------------------------
// trig()  -- generates one event
// trig(n) -- generates n+1 events.
//
// NOTE:  last event generated will be corrupt in the FZD file
//
void trig(Int_t n = 1){
    chain->EventLoop(n);
    _primary->event()->Print();
}
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
void Pythia8(){
    // Create the pythia 8 event generator and add it to 
    // the primary generator
    StarPythia8* pythia8 = new StarPythia8("pythia8");

    pythia8->SetFrame("CMS", 510.0);
    pythia8->SetBlue("proton");
    pythia8->SetYell("proton");

    pythia8->Set("SoftQCD:centralDiffractive = on");
    pythia8->Set("SigmaTotal:zeroAXB = off");

    //setting only charged decay products
    // //K0S
    // pythia8->Set("310:onMode=0");
    // pythia8->Set("310:OnIfMatch=211 -211");
    // pythia8->Set("-310:onMode=0");
    // pythia8->Set("-310:OnIfMatch=-211 211");
    //Lambda0
    pythia8->Set("3122:onMode=0");
    pythia8->Set("3122:OnIfMatch=2212 -211");
    pythia8->Set("-3122:onMode=0");
    pythia8->Set("-3122:OnIfMatch=-2212 211");
    // //K*(892)
    // pythia8->Set("313:onMode=0");
    // pythia8->Set("313:OnIfMatch=321 -211");
    // pythia8->Set("-313:onMode=0");
    // pythia8->Set("-313:OnIfMatch=-321 211");
    // //phi
    // pythia8->Set("333:onMode=0");
    // pythia8->Set("333:OnIfMatch=321 -321");
    // pythia8->Set("-333:onMode=0");
    // pythia8->Set("-333:OnIfMatch=-321 321");

    _primary->AddGenerator(pythia8);
}
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
void centralDiffractiveLambda0(Int_t nevents = 10, Int_t runNumber = 18091010, Int_t rngSeed = 1234){
    gROOT->ProcessLine(".L bfc.C");
    {
        TString simple = "y2017a geant gstar usexgeom agml";
        bfc(0, simple);
    }

    gSystem->Load("libVMC.so");
    gSystem->Load("StarGeneratorUtil.so");
    gSystem->Load("StarGeneratorEvent.so");
    gSystem->Load("StarGeneratorBase.so");
    gSystem->Load("./StRoot/StarGeneratorFilt.so");
    gSystem->Load("Pythia8_3_03.so");
    gSystem->Load("libMathMore.so");

    // Force loading of xgeometry
    gSystem->Load("xgeometry.so");

    // Setup RNG seed and map all ROOT TRandom here
    StarRandom::seed(rngSeed);
    StarRandom::capture();

    // Create the primary event generator and insert it
    // before the geant maker
    _primary = new StarPrimaryMaker();
    _primary->SetFileName("results/centralDiffractive.root");
    chain->AddBefore("geant", _primary);

    //addig filter to primary event maker
    filter = new ParticleLambda0Filter();
    filter->SetAttr(".Privilege", 1);
    _primary->AddFilter(filter);
    //not needed, set to 1 if you want to keep all despite filter working
    _primary->SetAttr("FilterKeepAll", int(0));
    //not needed, set to 1 if you want 1000 tries instead of 1000 events generated
    _primary->SetAttr("FilterSkipRejects", int(0));
    //ACTUALLY NEEDED, set to 1 if you want to keep rejected events
    _primary->SetAttr("FilterKeepHeader", int(0));
    //needed for ability to skip events
    _primary->SetAttr(".Privilege", 1);

    Pythia8();

    // Setup cuts on which particles get passed to geant for
    //   simulation.  (To run generator in standalone mode,
    //   set ptmin=1.0E9.)
    _primary->SetPtRange(0.0, -1.0);         // GeV
    _primary->SetEtaRange(-3.0, +3.0);
    _primary->SetPhiRange(0., TMath::TwoPi());

    // Setup a realistic vertex distribution:
    //z parameter fit from data
    //TODO: actually change it to a reasonable value
    _primary->SetSigma(0.015, 0.015, 50.);

    //settings taken from the database for the run given
    //default run 18091010
    //RTS Start Time 2017-04-01 08:36:11 GMT
    double beamline[4] = { 0,0,0,0 };
    if(runNumber==0)
        runNumber = 18091010;
    getBeamlineParameters(runNumber, beamline);
    printf("Beamline parameters:\nx0:\t%lf\ny0:\t%lf\ndxdz:\t%lf\ndydz:\t%lf\n", beamline[0], beamline[1], beamline[2], beamline[3]);
    _primary->SetVertex(beamline[0], beamline[1], 0.);
    _primary->SetSlope(beamline[2], beamline[3]);

    _primary->Init();

    // Setup geometry
    geometry("y2017a field=-5.0");
    command("gkine -4 0");
    command("gfile o results/centralDiffractive.fzd");

    trig(nevents);

    chain->Finish();

    command("call agexit");  // Make sure that STARSIM exits properly

}
