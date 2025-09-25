// macro to instantiate the Geant3 from within
// STAR  C++  framework and get the starsim prompt
// To use it do
//  root4star starsim.C

class St_geant_Maker;
St_geant_Maker *geant_maker = 0;

class StarGenEvent;
StarGenEvent   *event       = 0;

class StarPrimaryMaker;
StarPrimaryMaker *_primary = 0;

// ----------------------------------------------------------------------------
void geometry( TString tag, Bool_t agml=true ){
    TString cmd = "DETP GEOM "; cmd += tag;
    if ( !geant_maker ) geant_maker = (St_geant_Maker *)chain->GetMaker("geant");
    geant_maker -> LoadGeometry(cmd);
    //  if ( agml ) command("gexec $STAR_LIB/libxgeometry.so");
}
// ----------------------------------------------------------------------------
void command( TString cmd ){
    if ( !geant_maker ) geant_maker = (St_geant_Maker *)chain->GetMaker("geant");
    geant_maker -> Do( cmd );
}
// ----------------------------------------------------------------------------
// trig()  -- generates one event
// trig(n) -- generates n+1 events.
//
// NOTE:  last event generated will be corrupt in the FZD file
//
void trig( Int_t n=1 ){
    chain->EventLoop(n);
    _primary->event()->Print();
}
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
void Pythia8(){
    // Create the pythia 8 event generator and add it to 
    // the primary generator
    StarPythia8 *pythia8 = new StarPythia8();  

    pythia8->SetFrame("CMS", 510.0);
    pythia8->SetBlue("proton");
    pythia8->SetYell("proton");
    
    pythia8->Set("SoftQCD:elastic = on");
    pythia8->Set("SigmaTotal:zeroAXB = off");
    //K0S exclusive
    // pythia8->Set("310:onMode=0");
    // pythia8->Set("310:OnIfMatch=211 -211");
    // pythia8->Set("-310:onMode=0");
    // pythia8->Set("-310:OnIfMatch=-211 211");
        
    _primary -> AddGenerator( pythia8 );
}
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
void elastic( Int_t nevents=10, Int_t rngSeed=1234 ){ 
    gROOT->ProcessLine(".L bfc.C");
    {
        TString simple = "y2017a geant gstar usexgeom agml";
    bfc(0, simple );
    }

    gSystem->Load("libVMC.so");
    gSystem->Load("StarGeneratorUtil.so");
    gSystem->Load("StarGeneratorEvent.so");
    gSystem->Load("StarGeneratorBase.so");
    gSystem->Load("Pythia8_3_03.so");
    gSystem->Load("libMathMore.so");  

    // Force loading of xgeometry
    gSystem->Load("xgeometry.so");

    // Setup RNG seed and map all ROOT TRandom here
    StarRandom::seed( rngSeed );
    StarRandom::capture();

    // Create the primary event generator and insert it
    // before the geant maker
    _primary = new StarPrimaryMaker();
    _primary -> SetFileName( "elastic.root");
    _primary -> SetVertex( 0.1, -0.1, 0.0 );
    _primary->SetSigma(0.015, 0.015, 50.);
    chain -> AddBefore( "geant", _primary );

    Pythia8();

    // Setup cuts on which particles get passed to geant for
    //   simulation.  (To run generator in standalone mode,
    //   set ptmin=1.0E9.)
    _primary->SetPtRange  (0.0,  -1.0);         // GeV
    _primary->SetEtaRange ( -3.0, +3.0 );
    _primary->SetPhiRange ( 0., TMath::TwoPi() );

    // Setup a realistic z-vertex distribution:
    //   x = 0 gauss width = 1mm
    //   y = 0 gauss width = 1mm
    //   z = 0 gauss width = 30cm
    _primary->SetVertex( 0., 0., 0. );
    _primary->SetSigma(0.015, 0.015, 50.);

    _primary -> Init();

    // Setup geometry
    geometry("y2017a field=-5.0");
    command("gkine -4 0");
    command("gfile o elastic.fzd");

    trig( nevents );

    chain->Finish();

    command("call agexit");  // Make sure that STARSIM exits properly

}