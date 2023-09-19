#ifndef StUPCFilterTrgUtil_h
#define StUPCFilterTrgUtil_h

class StUPCFilterTrgUtil {

public:

  StUPCFilterTrgUtil();
  virtual ~StUPCFilterTrgUtil() {}

  void processEvent(const StTriggerData *trgdat, StUPCEvent *upcEvt);

private:

  StUPCFilterTrgUtil(const StUPCFilterTrgUtil &o); // not implemented
  StUPCFilterTrgUtil &operator=(const StUPCFilterTrgUtil &o); // not implemented

  void runZDC(const StTriggerData *trgdat, StUPCEvent *upcEvt);
  void runVPD(const StTriggerData *trgdat, StUPCEvent *upcEvt);
  void runBBC(const StTriggerData *trgdat, StUPCEvent *upcEvt);

  Bool_t bbcLimits(UShort_t adc, UShort_t tdc) const;

  //BBC tile TDC window and ADC threshold
  static const UShort_t mBBCTDCmin = 100;
  static const UShort_t mBBCTDCmax = 2400;
  static const UShort_t mBBCADCmin = 20;
  static const Int_t mIPmtLastSmall = 16; // ipmt for last small BBC tile



};

#endif

