
//_____________________________________________________________________________
//    Class for UPC data
//    Author: Jaroslav Adam
//
//    contains parameters of BEMC cluster relevant for UPC analysis
//_____________________________________________________________________________


//local headers
#include "StUPCBemcCluster.h"

ClassImp(StUPCBemcCluster);

//_____________________________________________________________________________
StUPCBemcCluster::StUPCBemcCluster(): TObject(),
  mEta(0), mPhi(0), mSigmaEta(0), mSigmaPhi(0), mEnergy(0), mId(0),
  mHTEnergy(0), mHTsoftID(0)
{
  //default constructor

}//StUPCBemcCluster
