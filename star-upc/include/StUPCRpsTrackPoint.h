#ifndef StUPCRpsTrackPoint_hh
#define StUPCRpsTrackPoint_hh

//_____________________________________________________________________________
//    Class for making picoDst RP data 2019
//    Author: Truhlar Tomas
//_____________________________________________________________________________


#include "TObject.h" 
#include "TVector3.h" 


class StUPCRpsTrackPoint : public TObject {
public:
	StUPCRpsTrackPoint();
	virtual ~StUPCRpsTrackPoint();

	enum StUPCRpsTrackPointQuality {rpsNormal, rpsGolden, rpsNotSet}; 

	void Clear(Option_t * /*option*/ ="");
	Int_t planesUsed() const;


	static const Int_t mNumberOfPmtsInRp = 2;
	static const Int_t mNumberOfPlanesInRp = 4;

	TVector3 positionVec() const { return mPosition; }
	UChar_t rpId() const { return mRpId; }
	Int_t clusterId(Int_t planeId ) const
	{
	    return planeId<mNumberOfPlanesInRp ? mClusterId[planeId] : -1;
	}
	Float_t time(Int_t pmtId) const
	{
	    return pmtId<mNumberOfPmtsInRp ? mTime[pmtId] : -1;
	}
	StUPCRpsTrackPointQuality quality() const { return mQuality; }
	Float_t x() const { return mPosition.x(); }
	Float_t y() const { return mPosition.y(); }
	Float_t z() const { return mPosition.z(); }

	void setPosition(const TVector3& position)
	{
	    mPosition = position;
	}
	void setRpId(UChar_t rpId) { mRpId = rpId; }
	void setClusterId(Int_t clusterId, Int_t planeId)
	{
	    if( planeId<mNumberOfPlanesInRp )
	        mClusterId[planeId] = clusterId;
	}
	void setTime(Float_t timeVal, UInt_t pmtId)
	{
	    if( pmtId<mNumberOfPmtsInRp ) mTime[pmtId] = timeVal;
	}
	void setQuality(StUPCRpsTrackPointQuality quality )
	{
	    mQuality = quality;
	}
    
private:
	StUPCRpsTrackPoint(const StUPCRpsTrackPoint &o); //not implemented
	StUPCRpsTrackPoint& operator=(const StUPCRpsTrackPoint&); 
	
	Float_t mTime[mNumberOfPmtsInRp];
	UChar_t mRpId; // Id of Roman Pot, E1U = 0, E1D = 1, E2U = 2, E2D, W1U, W1D, W2U, W2D = 7
	Int_t mClusterId[mNumberOfPlanesInRp];

	TVector3         mPosition;
	StUPCRpsTrackPointQuality mQuality;

	ClassDef( StUPCRpsTrackPoint, 1 )
};

#endif
