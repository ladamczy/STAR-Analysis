// #ifdef __CLING__
#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class StarGenerator+;
#pragma link C++ class StarGenPPEvent+;
#pragma link C++ class StarGenEvent+;
// #pragma link C++ class StarGenEventReader+;
#pragma link C++ class StarGenParticle+;
#pragma link C++ class StarGenStats+;
#pragma link C++ class StarParticleData+;
#pragma link C++ class StTestMaker+;
// #pragma link C++ class TColumnView+;
#pragma link C++ class TDataSet+;
// #pragma link C++ class StCloseFileOnTerminate+;
// #pragma link C++ class StNotifyOnTerminate+;
#pragma link C++ class TDataSetIter+;
// #pragma link C++ class TChair+;
#pragma link C++ class StEvtHddr+;
// #pragma link C++ class TAttr+;
#pragma link C++ class TObjectSet+;
// #pragma link C++ class TUnixTime+;
//those classes have a custom streamer
//so rootcint shouldn't generate a custom one
#pragma link C++ class StMaker-;
#pragma link C++ class StMemStat-;
// #pragma link C++ class TTableDescriptor-;
#pragma link C++ class TTableMap-;
// #pragma link C++ class StChain-;
#pragma link C++ class TTable-;

#pragma link C++ class StUPCEvent+;
#pragma link C++ class StUPCTrack+;
#pragma link C++ class StUPCBemcCluster+;
#pragma link C++ class StUPCVertex+;
#pragma link C++ class StRPEvent+;
#pragma link C++ class StUPCRpsTrack+;
#pragma link C++ class StUPCRpsTrackPoint+;
#pragma link C++ class StUPCRpsCluster+;
#pragma link C++ class StUPCTofHit+;
#pragma link C++ class StUPCV0+;
#pragma link C++ class StPicoPhysicalHelix+;
#pragma link C++ class StPicoHelix+;
#endif

