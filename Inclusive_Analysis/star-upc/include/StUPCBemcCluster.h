#ifndef StUPCBemcCluster_h
#define StUPCBemcCluster_h

//_____________________________________________________________________________
//    Class for UPC data
//    Author: Jaroslav Adam
//_____________________________________________________________________________

#include "TObject.h"

class StUPCBemcCluster: public TObject
{

public:

  StUPCBemcCluster();
  ~StUPCBemcCluster() {}

  //setters
  void setEta(Float_t eta) { mEta = eta; }
  void setPhi(Float_t phi) { mPhi = phi; }
  void setSigmaEta(Float_t se) { mSigmaEta = se; }
  void setSigmaPhi(Float_t sp) { mSigmaPhi = sp; }
  void setEnergy(Float_t en) { mEnergy = en; }
  void setId(UInt_t id) { mId = id; }

  void setHTEnergy(Float_t en) { mHTEnergy = en; }
  void setHTsoftID(Int_t id) { mHTsoftID = id; }

  //getters
  Float_t getEta() const { return mEta; }
  Float_t getPhi() const { return mPhi; }
  Float_t getSigmaEta() const { return mSigmaEta; }
  Float_t getSigmaPhi() const { return mSigmaPhi; }
  Float_t getEnergy() const { return mEnergy; }
  UInt_t getId() const { return mId; }

  Float_t getHTEnergy() const { return mHTEnergy; }
  Int_t getHTsoftID() const { return mHTsoftID; }

private:

  StUPCBemcCluster(const StUPCBemcCluster &o); //not implemented
  StUPCBemcCluster &operator=(const StUPCBemcCluster &o); //not implemented

  Float_t mEta; // eta of cluster
  Float_t mPhi; // phi of cluster
  Float_t mSigmaEta; // sigma of eta
  Float_t mSigmaPhi; // sigma of phi
  Float_t mEnergy; // energy of cluster;
  UInt_t mId; // original ID

  Float_t mHTEnergy; // energy of highest energy cell in cluster (high tower)
  Int_t mHTsoftID; // softID of highest energy cell






  ClassDef(StUPCBemcCluster, 2)

};

#endif






































