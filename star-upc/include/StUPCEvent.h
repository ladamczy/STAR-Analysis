#ifndef StUPCEvent_h
#define StUPCEvent_h

//_____________________________________________________________________________
//    Class for UPC data
//    Author: Jaroslav Adam
//_____________________________________________________________________________

class StUPCTrack;
class TClonesArray;
class StUPCBemcCluster;
class StUPCTofHit;
class TIterator;
class StUPCVertex;
class TParticle;

#include "TArrayI.h"

class StUPCEvent
{

public:

  StUPCEvent();
  virtual ~StUPCEvent();

  void clearEvent();

  //setters
  void addTriggerId(Int_t id);
  void setRunNumber(Int_t run) { mRunNum = run; }
  void setEventNumber(Int_t num) { mEvtNum = num; }
  void setFillNumber(Int_t num) { mFillNum = num; }
  void setBunchCrossId(UInt_t id) { mbCrossId = id; }
  void setBunchCrossId7bit(UInt_t id) { mbCrossId7bit = id; }
  void setMagneticField(Double_t mag) { mMagField = mag; }

  void setZDCEastRate(Double_t rate) { mZdcEastRate = rate; }
  void setZDCWestRate(Double_t rate) { mZdcWestRate = rate; }
  void setZDCCoincRate(Double_t rate) { mZdcCoincRate = rate; }

  void setLastDSM0(UShort_t dsm) { mLastDSM0 = dsm; }
  void setLastDSM1(UShort_t dsm) { mLastDSM1 = dsm; }
  void setLastDSM3(UShort_t dsm) { mLastDSM3 = dsm; }

  void setZDCUnAttEast(UShort_t signal) { mZdcEastUA = signal; }
  void setZDCUnAttWest(UShort_t signal) { mZdcWestUA = signal; }
  void setZDCEastADC(UShort_t signal, Int_t pmt);
  void setZDCWestADC(UShort_t signal, Int_t pmt);
  void setZDCEastTDC(UShort_t time) { mZdcEastTDC = time; }
  void setZDCWestTDC(UShort_t time) { mZdcWestTDC = time; }
  void setZDCTimeDiff(UShort_t deltaT) { mZdcTimeDiff = deltaT; }
  void setZdcVertexZ(Float_t vtx) { mZdcVertexZ = vtx; }

  void setZdcSMDEast(Int_t verthori, Int_t strip, UShort_t smd);
  void setZdcSMDWest(Int_t verthori, Int_t strip, UShort_t smd);

  void setBBCSmallEast(UInt_t sum) { mBBCSmallEast = sum; }
  void setBBCSmallWest(UInt_t sum) { mBBCSmallWest = sum; }
  void setBBCLargeEast(UInt_t sum) { mBBCLargeEast = sum; }
  void setBBCLargeWest(UInt_t sum) { mBBCLargeWest = sum; }

  void setVPDSumEast(UShort_t sum) { mVPDSumEast = sum; }
  void setVPDSumWest(UShort_t sum) { mVPDSumWest = sum; }
  void setVPDTimeDiff(UShort_t deltaT) { mVPDTimeDiff = deltaT; }

  void setTOFMultiplicity(UShort_t tof) { mTofMult = tof; }
  void setBEMCMultiplicity(UInt_t be) { mBemcMult = be; }

  void setNGlobTracks(Int_t ntrk) { mNGlobTracks = ntrk; }
  void setNPrimTracks(Int_t ntrk) { mNPrimTracks = ntrk; }
  void setNPrimVertices(UInt_t nvtx) { mNPrimVertices = nvtx; }

  StUPCTrack *addTrack();

  StUPCBemcCluster *addCluster();

  StUPCTofHit *addHit();

  StUPCVertex *addVertex();

  void setIsMC(Bool_t mc=kTRUE);
  TParticle *addMCParticle();

  //getters
  Bool_t isTrigger(Int_t id) const;
  Int_t getRunNumber() const { return mRunNum; }
  Int_t getEventNumber() const { return mEvtNum; }
  Int_t getFillNumber() const { return mFillNum; }
  UInt_t getBunchCrossId() const { return mbCrossId; }
  UInt_t getBunchCrossId7bit() const { return mbCrossId7bit; }
  Double_t getMagneticField() const { return mMagField; }

  Double_t getZDCEastRate() const { return mZdcEastRate; }
  Double_t getZDCWestRate() const { return mZdcWestRate; }
  Double_t getZDCCoincRate() const { return mZdcCoincRate; }

  UShort_t getLastDSM0() const { return mLastDSM0; }
  UShort_t getLastDSM1() const { return mLastDSM1; }
  UShort_t getLastDSM3() const { return mLastDSM3; }

  UShort_t getZDCUnAttEast() const { return mZdcEastUA; }
  UShort_t getZDCUnAttWest() const { return mZdcWestUA; }
  UShort_t getZDCEastADC(Int_t pmt) const { return mZdcEastADC[pmt-1]; }
  UShort_t getZDCWestADC(Int_t pmt) const { return mZdcWestADC[pmt-1]; }
  UShort_t getZDCEastTDC() const { return mZdcEastTDC; }
  UShort_t getZDCWestTDC() const { return mZdcWestTDC; }
  UShort_t getZDCTimeDiff() const { return mZdcTimeDiff; }
  Float_t getZdcVertexZ() const { return mZdcVertexZ/10.; } // convert from mm to cm

  UShort_t getZdcSMDEast(Int_t verthori, Int_t strip) const { return mZdcSmdEast[strip-1 + verthori*8]; }
  UShort_t getZdcSMDWest(Int_t verthori, Int_t strip) const { return mZdcSmdWest[strip-1 + verthori*8]; }

  UInt_t getBBCSmallEast() const { return mBBCSmallEast; }
  UInt_t getBBCSmallWest() const { return mBBCSmallWest; }
  UInt_t getBBCLargeEast() const { return mBBCLargeEast; }
  UInt_t getBBCLargeWest() const { return mBBCLargeWest; }

  UShort_t getVPDSumEast() const { return mVPDSumEast; }
  UShort_t getVPDSumWest() const { return mVPDSumWest; }
  UShort_t getVPDTimeDiff() const { return mVPDTimeDiff; }

  UShort_t getTOFMultiplicity() const { return mTofMult; }
  UInt_t getBEMCMultiplicity() const { return mBemcMult; }

  Int_t getNGlobTracks() const { return mNGlobTracks; }
  Int_t getNPrimTracks() const { return mNPrimTracks; }
  UInt_t getNPrimVertices() const { return mNPrimVertices; }

  Int_t getNumberOfTracks() const;
  StUPCTrack *getTrack(Int_t iTrack) const;

  Int_t getNumberOfClusters() const;
  StUPCBemcCluster *getCluster(Int_t iCls) const;
  StUPCBemcCluster *getClusterId(UInt_t clsId) const;
  TIterator *makeClustersIter() const;

  Int_t getNumberOfHits() const;
  StUPCTofHit *getHit(Int_t iHit) const;

  Int_t getNumberOfVertices() const;
  StUPCVertex *getVertex(Int_t iVtx) const;
  StUPCVertex *getVertexId(UInt_t vtxId) const;
  TIterator *makeVerticesIter() const;

  Bool_t getIsMC() const { return mMCParticles != NULL ? kTRUE : kFALSE; }
  Int_t getNumberOfMCParticles() const;
  TParticle *getMCParticle(Int_t iMC) const;

private:

  StUPCEvent(const StUPCEvent &o); //not implemented
  StUPCEvent &operator=(const StUPCEvent &o); //not implemented

  TArrayI mTrgIDs; // fired trigger IDs

  Int_t mRunNum; // number of current run
  Int_t mEvtNum; // event number
  Int_t mFillNum; // beam fill number
  UInt_t mbCrossId; // bunch crossing ID
  UInt_t mbCrossId7bit; // bunch crossing ID 7bit
  Double32_t mMagField; // magnetic field

  Double32_t mZdcEastRate; // ZDC rate, east
  Double32_t mZdcWestRate; // ZDC rate, west
  Double32_t mZdcCoincRate; // ZDC rate, east&&west coincidence

  UShort_t mLastDSM0; // TCU bits 0-15 from StTriggerData::lastDSM(0)
  UShort_t mLastDSM1; // TCU bits 16-31
  UShort_t mLastDSM3; // TCU bits 58-62

  UShort_t mZdcEastUA; // ZDC unattenuated signal, east
  UShort_t mZdcWestUA; // ZDC unattenuated signal, west
  static const Int_t mNZdcPmt = 3;
  UShort_t mZdcEastADC[mNZdcPmt]; // ZDC 3 PMT ADCs, east
  UShort_t mZdcWestADC[mNZdcPmt]; // ZDC 3 PMT ADCs, west
  UShort_t mZdcEastTDC; // ZDC TDC, east
  UShort_t mZdcWestTDC; // ZDC TDC, west
  UShort_t mZdcTimeDiff; // ZDC time difference, east-west
  Float_t mZdcVertexZ; // ZDC vertex z position, mm
  static const Int_t mNZdcSmd = 16; // number of ZDC SMD channels, east and west
  UShort_t mZdcSmdEast[mNZdcSmd]; // ZDC SMD data east
  UShort_t mZdcSmdWest[mNZdcSmd]; // ZDC SMD data west

  UInt_t mBBCSmallEast; // BBC truncated sum, small tiles, east
  UInt_t mBBCSmallWest; // BBC truncated sum, small tiles, west
  UInt_t mBBCLargeEast; // BBC truncated sum, large tiles, east
  UInt_t mBBCLargeWest; // BBC truncated sum, large tiles, west

  UShort_t mVPDSumEast; // VPD ADC sum east
  UShort_t mVPDSumWest; // VPD ADC sum west
  UShort_t mVPDTimeDiff; // VPD time difference

  UShort_t mTofMult; // TOF multiplicity from StTriggerData
  UInt_t mBemcMult; // BEMC multiplicity, number of all BEMC clusters in event

  Int_t mNGlobTracks; // number of global tracks
  Int_t mNPrimTracks; // number of primary tracks
  UInt_t mNPrimVertices; // number of primary vertices

  static TClonesArray *mgUPCTracks; // array of upc tracks
  TClonesArray *mUPCTracks; //-> array of upc tracks
  Int_t mNtracks; //! number of upc tracks in event, local use when filling

  static TClonesArray *mgUPCBemcClusters; // array of BEMC clusters
  TClonesArray *mUPCBemcClusters; //-> array of BEMC clusters
  Int_t mNclusters; //! number of BEMC clusters written in event, local use when filling

  static TClonesArray *mgUPCTOFHits; // array of upc TOF hits
  TClonesArray *mUPCTOFHits; //-> array of upc TOF hits
  Int_t mNhits; //! number of TOF hits written in event, local use when filling

  static TClonesArray *mgUPCVertices; // array of UPC vertices
  TClonesArray *mUPCVertices; //-> array of UPC vertices
  Int_t mNvertices; //! number of vertices written in event, local use when filling

  static TClonesArray *mgMCParticles; // array of MC particles
  TClonesArray *mMCParticles; // array of MC particles
  Int_t mNmc; //! number of mc particles in event, local use when filling

  ClassDef(StUPCEvent, 3);
};

#endif













