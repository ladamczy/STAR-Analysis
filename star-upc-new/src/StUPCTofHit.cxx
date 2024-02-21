
//_____________________________________________________________________________
//    Class for UPC data
//    Author: Tomas Truhlar
//
//    contains parameters of TOF hit relevant for UPC+RP analysis
//_____________________________________________________________________________

//local headers
#include "StUPCTofHit.h"

ClassImp(StUPCTofHit);

//_____________________________________________________________________________
StUPCTofHit::StUPCTofHit(): TObject(),
  mTray(0), mModule(0), mCell(0), mLeadingEdgeTime(0), mTrailingEdgeTime(0)
{
  //default constructor

}//StUPCTofHit
