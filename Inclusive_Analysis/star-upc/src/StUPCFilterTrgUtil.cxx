
//_____________________________________________________________________________
//
//    Utility class for trigger, BBC and ZDC in UPC filter maker
//    Author: Jaroslav Adam
//
//_____________________________________________________________________________

//StRoot headers
#include "StEvent/StTriggerData.h"

//local headers
#include "StUPCEvent.h"

#include "StUPCFilterTrgUtil.h"

//_____________________________________________________________________________
StUPCFilterTrgUtil::StUPCFilterTrgUtil() {
  //constructor

}//StUPCFilterTrgUtil

//_____________________________________________________________________________
void StUPCFilterTrgUtil::processEvent(const StTriggerData *trgdat, StUPCEvent *upcEvt) {

  //TCU bits
  upcEvt->setLastDSM0( trgdat->lastDSM(0) );
  upcEvt->setLastDSM1( trgdat->lastDSM(1) );
  upcEvt->setLastDSM3( trgdat->lastDSM(3) );

  //ZDC
  runZDC(trgdat, upcEvt);
  //VPD
  runVPD(trgdat, upcEvt);
  //BBC
  runBBC(trgdat, upcEvt);
  //TOF multiplicity
  upcEvt->setTOFMultiplicity( trgdat->tofMultiplicity() );

}//processEvent

//_____________________________________________________________________________
void StUPCFilterTrgUtil::runZDC(const StTriggerData *trgdat, StUPCEvent *upcEvt) {

  //function to write ZDC data to output UPC event

  upcEvt->setZDCUnAttEast( trgdat->zdcUnAttenuated(east) );
  upcEvt->setZDCUnAttWest( trgdat->zdcUnAttenuated(west) );
  for (Int_t ipmt=1; ipmt<=3; ipmt++) {
    upcEvt->setZDCEastADC( trgdat->zdcADC(east, ipmt), ipmt );
    upcEvt->setZDCWestADC( trgdat->zdcADC(west, ipmt), ipmt );
  }
  upcEvt->setZDCEastTDC( trgdat->zdcTDC(east) );
  upcEvt->setZDCWestTDC( trgdat->zdcTDC(west) );
  upcEvt->setZDCTimeDiff( trgdat->zdcTimeDifference() );
  upcEvt->setZdcVertexZ( trgdat->zdcVertexZ() );

  //ZDC SMD
  //vertical - horizontal loop
  for(Int_t verthori=0; verthori<2; verthori++) {
    //strip loop
    for(Int_t strip=1; strip<9; strip++) {
      upcEvt->setZdcSMDEast(verthori, strip, trgdat->zdcSMD(east, verthori, strip));
      upcEvt->setZdcSMDWest(verthori, strip, trgdat->zdcSMD(west, verthori, strip));
    }//strip loop
  }//vertical - horizontal loop

}//runZDC

//_____________________________________________________________________________
void StUPCFilterTrgUtil::runVPD(const StTriggerData *trgdat, StUPCEvent *upcEvt) {

  //function to write VPD data to output UPC event

  UShort_t sum_east=0, sum_west=0;
  //pmt loop
  for(Int_t ipmt=1; ipmt<=16; ipmt++) {
    sum_east += trgdat->vpdADC(east, ipmt);
    sum_west += trgdat->vpdADC(west, ipmt);
  }//pmt loop
  upcEvt->setVPDSumEast( sum_east );
  upcEvt->setVPDSumWest( sum_west );

  upcEvt->setVPDTimeDiff( trgdat->vpdTimeDifference() );

}//runVPD

//_____________________________________________________________________________
void StUPCFilterTrgUtil::runBBC(const StTriggerData *trgdat, StUPCEvent *upcEvt) {

  //funcion to write BBC data to output UPC event

  //BBC small and large tiles truncated sum
  UInt_t bbcSmallE = 0, bbcSmallW = 0, bbcLargeE = 0, bbcLargeW = 0;

  //BBC tile loop
  for(Int_t ipmt=1; ipmt<=24; ipmt++) {

    //east adc and tdc
    UShort_t adcE = trgdat->bbcADC(east, ipmt);
    UShort_t tdcE = trgdat->bbcTDC(east, ipmt);

    //west adc and tdc
    UShort_t adcW = trgdat->bbcADC(west, ipmt);
    UShort_t tdcW = trgdat->bbcTDC(west, ipmt);

    //small/large sum after ADC threshold and TDC window
    if( ipmt <= mIPmtLastSmall ) {
      if( bbcLimits(adcE, tdcE) ) bbcSmallE += adcE;
      if( bbcLimits(adcW, tdcW) ) bbcSmallW += adcW;
    }
    if( ipmt >  mIPmtLastSmall ) {
      if( bbcLimits(adcE, tdcE) ) bbcLargeE += adcE;
      if( bbcLimits(adcW, tdcW) ) bbcLargeW += adcW;
    }

  }//BBC tile loop

  //set BBC truncated sum
  upcEvt->setBBCSmallEast(bbcSmallE);
  upcEvt->setBBCSmallWest(bbcSmallW);
  upcEvt->setBBCLargeEast(bbcLargeE);
  upcEvt->setBBCLargeWest(bbcLargeW);

}//runBBC

//_____________________________________________________________________________
Bool_t StUPCFilterTrgUtil::bbcLimits(UShort_t adc, UShort_t tdc) const
{
  //function to compare BBC ADC and TDC limits

  //ADC threshold
  if( adc < mBBCADCmin ) return kFALSE;

  //select TDC window
  if( tdc < mBBCTDCmin || tdc > mBBCTDCmax ) return kFALSE;

  return kTRUE;

}//bbcLimits























