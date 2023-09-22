#ifndef StUPCRpsTrack_hh
#define StUPCRpsTrack_hh

//_____________________________________________________________________________
//    Class for making picoDst RP data 2019
//    Author: Truhlar Tomas
//_____________________________________________________________________________

#include "TObject.h"
#include "TRef.h"
#include "TClonesArray.h"
#include "TVector3.h"

using namespace std;

class StUPCRpsTrackPoint;
class StRPEvent;

class StUPCRpsTrack : public TObject {
public:
	StUPCRpsTrack();
	virtual ~StUPCRpsTrack();
	void Clear(Option_t * /*option*/ ="");

	enum StRpsTrackType { rpsLocal, rpsGlobal, rpsUndefined };
	enum StRpsAngles { rpsAngleThetaX, rpsAngleThetaY, rpsAngleTheta, mNumberOfAngleTypes };

	//getters
	unsigned int planesUsed() const;

	double theta(unsigned int = rpsAngleTheta) const;
	double thetaRp(unsigned int = rpsAngleTheta) const;
	double phiRp() const;
	double time() const;


	static const Int_t mNumberOfStationsInBranch = 2;

	Int_t getFirstTrackPointId() { return mFirstTrackPointId; }
	Int_t getSecondTrackPointId() { return mSecondTrackPointId; }
	TVector3 pVec() const { return mP; }
	Int_t branch() const { return mBranch; }
	StRpsTrackType type() const { return mType; }
	double phi() const { return mP.Phi(); } 
	double t(Float_t beamMomentum) const
	{
	  return -2*beamMomentum*beamMomentum*(1-xi(beamMomentum))*(1-cos(theta(rpsAngleTheta)));
	}
	double xi(Float_t beamMomentum) const
	{
	    return (beamMomentum - mP.Mag())/beamMomentum; 
	}
	double p() const { return mP.Mag(); }
	double pt() const { return mP.Perp(); } 
	double eta() const { return mP.PseudoRapidity(); } 
	StRPEvent *getEvent() const { return mEvt; }
	StUPCRpsTrackPoint* getTrackPoint(UInt_t  station) const;

	// setters
	void setFirstTrackPointId(Int_t Id) { mFirstTrackPointId = Id; }
	void setSecondTrackPointId(Int_t Id) { mSecondTrackPointId = Id; }
	void setP(const TVector3& P) { mP = P; }
	void setBranch(Int_t branch) { mBranch = branch; }
	void setType(StRpsTrackType type) { mType = type; }
	void setEvent(StRPEvent *evt) { mEvt = evt; }

private:	
	StUPCRpsTrack(const StUPCRpsTrack &o); //not implemented
	StUPCRpsTrack& operator=(const StUPCRpsTrack&); 

	// pointers to trackPoints (local tracks)
	Int_t mFirstTrackPointId;
	Int_t mSecondTrackPointId;
	
	TVector3 mP;				// three-vector with reconstructed track momentum
	StRpsTrackType mType;			// type of the track
	Int_t          mBranch;			// detectors branch, EU=0, ED=1, WU=2, WD=3 

	StRPEvent *mEvt; //! pointer to current event, local use only

	ClassDef(StUPCRpsTrack, 1)
};

#endif
