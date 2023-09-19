#ifndef StUPCTofHit_h
#define StUPCTofHit_h

//_____________________________________________________________________________
//    Class for UPC data
//    Author: Tomas Truhlar
//
//    contains parameters of TOF hit relevant for UPC+RP analysis
//_____________________________________________________________________________
#include "TObject.h"

class StUPCTofHit: public TObject
{

public:

  StUPCTofHit();
  ~StUPCTofHit() {}

  //setters
  void setTray(UChar_t tray)  { mTray = tray; }
  void setModule(UChar_t module)  { mModule = module; }
  void setCell(UChar_t cell)  { mCell = cell; }
  void setLeadingEdgeTime(Double_t leTime)  { mLeadingEdgeTime = leTime; }
  void setTrailingEdgeTime(Double_t teTime)  { mTrailingEdgeTime = teTime; }

  //getters
  UChar_t getTray() const { return mTray; }
  UChar_t getModule() const { return mModule; }
  UChar_t getCell() const { return mCell; }
  Double_t getLeadingEdgeTime() const { return mLeadingEdgeTime; }
  Double_t getTrailingEdgeTime() const { return mTrailingEdgeTime; }

private:

  StUPCTofHit(const StUPCTofHit &o); //not implemented
  StUPCTofHit &operator=(const StUPCTofHit &o); //not implemented

  UChar_t mTray; // Tray of TOF hit
  UChar_t mModule; // Module of TOF hit
  UChar_t mCell; // Cell of TOF hit
  Double_t mLeadingEdgeTime; // Leading Edge Time of TOF hit
  Double_t mTrailingEdgeTime; // Trailing Edge Time of TOF hit

  ClassDef(StUPCTofHit, 2)

};

#endif






































