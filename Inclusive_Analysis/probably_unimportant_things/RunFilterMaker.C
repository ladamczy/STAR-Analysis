
#include "TArgumentParser.h"

class StUPCFilterMaker;
StUPCFilterMaker *anaMaker;
void AddTrigger(UInt_t id, Int_t rmin, Int_t rmax);

//_____________________________________________________________________________
void RunFilterMaker(string filelist, Int_t nFiles, string outfile, string config) {

  //print the full list of files including their size
  PrintFilelist(filelist);

  //load libraries to work with muDst
  gROOT->Macro("loadMuDst.C");

  // Load St_db_Maker and co
  gSystem->Load("StDbLib.so");
  gSystem->Load("StDbBroker.so");
  gSystem->Load("St_db_Maker");

  // Load Emc libraries
  gSystem->Load("StDaqLib");
  gSystem->Load("StEmcRawMaker");
  gSystem->Load("StEmcADCtoEMaker");
  gSystem->Load("StPreEclMaker");
  gSystem->Load("StEpcMaker");

  //load the analysis maker compiled before with cons
  gSystem->Load("StUPCFilterMaker.so");
  //TOF calib maker
  gSystem->Load("StBTofCalibMaker.so");

  //create chain directory-like structure for maker
  //top level
  StChain *chain = new StChain;
  //maker to access muDST data
  StMuDstMaker *maker = new StMuDstMaker(0, 0, "", filelist.c_str(), "", nFiles);

  //St_db_Maker for Emc calibration
  St_db_Maker *db1 = new St_db_Maker("db","$HOME/StarDb","MySQL:StarDb","$STAR/StarDb");
  // Maker to apply calibration
  StEmcADCtoEMaker *adc_to_e = new StEmcADCtoEMaker();
  adc_to_e->setPrint(kFALSE);
  // Makers for clusterfinding
  StPreEclMaker *pre_ecl = new StPreEclMaker();
  pre_ecl->setPrint(kFALSE);
  StEpcMaker *epc = new StEpcMaker();
  epc->setPrint(kFALSE);
  //TOF calib
  StBTofCalibMaker *tofCalib = new StBTofCalibMaker();
  tofCalib->setMuDstIn(kTRUE);
  //analysis maker
  anaMaker = new StUPCFilterMaker(maker, outfile); //maker for muDst passed to the constructor

  //----------------------------------------------------------------------
  // analysis maker configuration, do not change default values here,
  // they are supposed to be loaded from configuration file

  //data/MC
  Int_t isMC = 0; // 0 - data,  1 - starsim MC,  2 - embedding MC

  // use BEMC cluster conditions below if set to true
  Bool_t useClusterParam = kFALSE;
  Int_t sizeMax = 4;
  Float_t energySeed = 0.4;
  Float_t energyAdd = 0.001;
  Float_t energyThresholdAll = 0.4;

  //use TOF start time override to zero
  Bool_t tof_start_zero = kFALSE;

  //write Roman Pot data
  Bool_t make_RP_event = kFALSE;

  //load values from config file
  TArgumentParser parser;
  parser.AddInt("is_mc", &isMC);
  parser.AddBool("bemc_cluster_param", &useClusterParam);
  parser.AddInt("bemc_size_max", &sizeMax);
  parser.AddFloat("bemc_energy_seed", &energySeed);
  parser.AddFloat("bemc_energy_add", &energyAdd);
  parser.AddFloat("bemc_energy_threshold_all", &energyThresholdAll);
  parser.AddBool("tof_start_zero", &tof_start_zero);
  parser.AddFuncInt3("add_trigger", AddTrigger);
  parser.AddBool("make_RP_event", &make_RP_event);

  cout << "RunFilterMaker, using config from: " << config << endl;
  parser.Parse(config);

  //----------------------------------------------------------------------

  //no debug printouts
  StMuDebug::setLevel(0);

  //show input and output for the maker
  cout << "RunFilterMaker, filelist: " << filelist << endl;
  cout << "RunFilterMaker, nFiles:   " << nFiles << endl;
  cout << "RunFilterMaker, outfile:  " << outfile << endl;

  //apply data/mc selection
  cout << "RunFilterMaker, isMC: " << isMC << endl;
  anaMaker->setIsMC(isMC);

  //apply TOF start time override
  cout << "RunFilterMaker, tof_start_zero: " << tof_start_zero << endl;
  if( tof_start_zero == kTRUE ) {
    tofCalib->forceTStartZero();
  }

  //Roman Pot data
  cout << "RunFilterMaker, make_RP_event: " << make_RP_event << endl;
  anaMaker->setMakeRPEvent(make_RP_event);

  Int_t nevt = maker->chain()->GetEntries();
  cout << "Number of events: " << nevt << endl;

  //initialize the makers
  chain->Init();

  //apply BEMC clustering parameters if requested
  if( useClusterParam ) {
    //call need to happen after the StPreEclMaker has been initialized
    pre_ecl->SetClusterConditions("bemc", sizeMax, energySeed, energyAdd, energyThresholdAll, kFALSE);
  }
  cout << "-------------------------------------------" << endl;
  cout << "StEmcOldFinder cluster parameters for BEMC:" << endl;
  StEmcOldFinder *finder = dynamic_cast<StEmcOldFinder*>(pre_ecl->finder());
  //bemc has index 1, according to StRoot/StEmcUtil/others/emcDetectorName.h and StPreEclMaker.cxx
  cout << "  sizeMax: " << finder->sizeMax(1) << endl;
  cout << "  energySeed: " << finder->energySeed(1) << endl;
  cout << "  energyAdd: " << finder->energyAdd(1) << endl;
  cout << "  energyThresholdAll: " << finder->energyThresholdAll(1) << endl;
  cout << "-------------------------------------------" << endl;

  //loop over events
  chain->EventLoop(nevt);
  chain->Finish();

  //release allocated memory
  delete chain;







}//RunFilterMaker

//_____________________________________________________________________________
void AddTrigger(UInt_t id, Int_t rmin, Int_t rmax) {

  //wrapper function to call maker::addTriggerId

  anaMaker->addTriggerId(id, rmin, rmax);

}//AddTrigger

//_____________________________________________________________________________
void PrintFilelist(const string& filelist) {

  //print the full list of files including their size

  cout << "------- In RunFilterMaker.C, the list of files is: -------" << endl;

  ifstream in(filelist.c_str());

  string line;
  while(in.good()) {
    getline(in, line);
    if(line.empty()) continue;

    string cmd = "ls " + line + " -alh";

    gSystem->Exec(cmd.c_str());
  }

  cout << "------------------ End of list of files ------------------" << endl;

  in.close();

}//PrintFilelist






