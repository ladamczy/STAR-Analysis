#ifndef StUPCFilterRPUtil_h
#define StUPCFilterRPUtil_h

//_____________________________________________________________________________
//    Class for making picoDst RP data 2019
//    Author: Truhlar Tomas
//_____________________________________________________________________________

class StRPEvent;
class StMuDst;
class StMuRpsTrackPoint;

class StUPCFilterRPUtil {

public:

  StUPCFilterRPUtil();
  virtual ~StUPCFilterRPUtil() {}

  void clear();

  void processEvent(StRPEvent *rpEvt, StMuDst *mMuDstw);

private:

	vector<const StMuRpsTrackPoint*>  mMuTrackPoints;

   StUPCFilterRPUtil(const StUPCFilterRPUtil &o); // not implemented
   StUPCFilterRPUtil &operator=(const StUPCFilterRPUtil &o); // not implemented

};

#endif

