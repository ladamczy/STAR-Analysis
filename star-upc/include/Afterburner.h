#ifndef Afterburner_h
#define Afterburner_h

#include <StUPCRpsCluster.h>
#include <Util.h>

//added includes because of compilation issues
#include <map>
#include "StRPEvent.h"
#include "StUPCRpsTrackPoint.h"
#include "StUPCRpsTrack.h"
#include <fstream>
#include <sstream>

const int nPlanesUsed = 3;


map<unsigned int, TVector3> mCorrection[nRomanPots];

bool LoadOffsetFile(TString fileName, map<unsigned int, TVector3> (&offsets)[nRomanPots]);
void runAfterburner(StRPEvent *event, StRPEvent *newRpEvent, int RunNumber);
TVector3 CalculateMomentumVector(double TPx, double TPy, double TPz, StUPCRpsTrack *track);

bool LoadOffsetFile(TString fileName, map<unsigned int, TVector3> (&offsets)[nRomanPots])
{
   string line;
   ifstream file( fileName );
   if (!file.is_open() )
   {
      cout << "\n ERROR in LoadOffsets(): Problems with opening a file: "<< fileName<< endl;
      return false;
   }   
   int RUNNUMBER;
   double var[2*nRomanPots]; // placeholder for offsets or offsets correction in x and y
   while ( getline(file, line) )
   {
      stringstream ss(line);
      ss >> RUNNUMBER;
      for (int i = 0; i < 2*nRomanPots; ++i)
         ss >> var[i]; 

      for (int iRP = 0; iRP < nRomanPots; ++iRP)
         offsets[iRP][RUNNUMBER].SetXYZ(var[2*iRP], var[2*iRP + 1], 0.0); 
   }

   file.close();
   return true;
}

void runAfterburner(StRPEvent *event, StRPEvent *newRpEvent, int RunNumber)
{
   for(unsigned int iCl = 0; iCl < event->getNumberOfClusters(); ++iCl){
      StUPCRpsCluster *oldCl = event->getCluster(iCl);
      float position = oldCl->position();
      int RP = oldCl->romanPotId();
      int plane = oldCl->planeId();
      // Set new position
      position -= mCorrection[RP][RunNumber][plane%2 ? 0 : 1];
      // Create new TrackPoint with new position
      StUPCRpsCluster *rpCluster = newRpEvent->addCluster();
      rpCluster->setPosition(position);
      rpCluster->setPositionRMS(oldCl->positionRMS()); 
      rpCluster->setLength(oldCl->length()); 
      rpCluster->setEnergy(oldCl->energy()); 
      rpCluster->setXY(oldCl->xy()); 
      rpCluster->setQuality(oldCl->quality()); 
      rpCluster->setPlaneId(plane);
      rpCluster->setRomanPotId(RP);
   }

   for(unsigned int iTP = 0; iTP < event->getNumberOfTrackPoints(); ++iTP){
      StUPCRpsTrackPoint *oldTP = event->getTrackPoint(iTP);
      TVector3 position = oldTP->positionVec();
      // Set new position
      int RP = oldTP->rpId();
      position -= mCorrection[RP][RunNumber];
      // Create new TrackPoint with new position
      StUPCRpsTrackPoint *rpTrackPoint = newRpEvent->addTrackPoint();
      rpTrackPoint->setPosition(position);
      rpTrackPoint->setRpId(RP);
      rpTrackPoint->setTime(oldTP->time(0), 0); // pmtId = 0
      rpTrackPoint->setTime(oldTP->time(1), 1); // pmtId = 1
      rpTrackPoint->setQuality(oldTP->quality());    
      for(UInt_t iPlane=0; iPlane < nPlanes; ++iPlane) // 4 planes in a RP
         rpTrackPoint->setClusterId(oldTP->clusterId(iPlane), iPlane);
   }

   for(unsigned int k = 0; k < event->getNumberOfTracks(); ++k)
   {// Get pointer to k-th track in Roman Pot data collection
      StUPCRpsTrack *track = event->getTrack(k);
      track->setEvent(event);

      if( !(track->getTrackPoint(RP1) ? track->getTrackPoint(RP1)->planesUsed()>=nPlanesUsed : false) ||
      !(track->getTrackPoint(RP2) ? track->getTrackPoint(RP2)->planesUsed()>=nPlanesUsed : false))
         continue;

      StUPCRpsTrack *rpTrack = newRpEvent->addTrack();
      rpTrack->setEvent(newRpEvent);
      rpTrack->setFirstTrackPointId(track->getFirstTrackPointId());
      rpTrack->setSecondTrackPointId(track->getSecondTrackPointId());

      StUPCRpsTrackPoint *newTP = track->getTrackPoint(RP1);
      TVector3 position = newTP->positionVec();

      rpTrack->setType(track->type());
      rpTrack->setBranch(track->branch());
      rpTrack->setP( CalculateMomentumVector(position.X(), position.Y(), position.Z(), rpTrack)); // setting the momentum vector
   }
}

TVector3 CalculateMomentumVector(double TPx, double TPy, double TPz, StUPCRpsTrack *track)
{
   // default values for pp run17 
   double mXYZ_IP[] = { -0.0011, 0, 0}; /* collision coordinates at the IP; 0=X, 1=Y, 2=Z */
   double mThetaXY_tilt[] = { -0.000070, 0}; /* tilt angles of the beam at collision; 0=X, 1=Y */
   double mDistanceFromIPtoDX[] = { 9.8, 9.8}; /* distance from the IP to the DX magnet in the East and West; 0=E, 1=W */
   double mLDX[] = { 3.7, 3.7};     /* length of DX in the East and West; 0=E, 1=W */
   double mBendingAngle[] = { 0.018832292, 0.018826657};     /* DX bending angles in the East and West; 0=E, 1=W */
   double beamMomenta[2] = { 254.867  , 254.867 };

   int sign = (TPz < 0 ? -1 : 1 );
   int iSide = (TPz < 0 ? E : W );
   // below calculating momentum vector

   if( track->type() == StUPCRpsTrack::rpsLocal)
   {
      double x_BCS = TPx - mXYZ_IP[X] - sin(mThetaXY_tilt[X])*( TPz - mXYZ_IP[2] ); // x_RP in beam coordinate system
      double y_BCS = TPy - mXYZ_IP[Y] - sin(mThetaXY_tilt[Y])*( TPz - mXYZ_IP[2] ); // y_RP in beam coordinate system
      double localThetaX = x_BCS / abs( TPz);
      double localThetaY = y_BCS / abs( TPz);
      double momentumValue = beamMomenta[iSide];

      TVector3 momentumVector( 0, 0, sign*momentumValue );
      //cout<<"Angle: "<<localThetaY<<endl;
      momentumVector.RotateX( -sign*localThetaY );
      momentumVector.RotateY( sign*localThetaX );
      return momentumVector;
   }

   double localThetaX = track->thetaRp( 0 ) - sign*mThetaXY_tilt[0]; // adding tilt of the beam (set to 0)
   double localThetaY = track->thetaRp( 1 ) - sign*mThetaXY_tilt[1]; 
   double x_BCS =  TPx - mXYZ_IP[0] - sin(mThetaXY_tilt[0])*( TPz - mXYZ_IP[2] ); // x_RP1 in beam coordinate system
   double d2 = abs( TPz ) - mLDX[iSide] - mDistanceFromIPtoDX[iSide]; // distance from DX magnet exit to first RP station
   double thetaX_IP = ( x_BCS - (d2 + 0.5*mLDX[iSide])*localThetaX ) / ( mDistanceFromIPtoDX[iSide]  + 0.5*mLDX[iSide] );
   double xi = 1. / ( 1 + (mBendingAngle[iSide]*(mDistanceFromIPtoDX[iSide] + 0.5*mLDX[iSide])) / ( localThetaX*abs( TPz ) - x_BCS ) );
   double momentumValue = beamMomenta[iSide] * (1.-xi);

   TVector3 momentumVector( 0, 0, sign*momentumValue );
   //cout<<"Angle: "<<localThetaY<<endl;
   momentumVector.RotateX( -sign*localThetaY );
   momentumVector.RotateY( sign*thetaX_IP );
   return momentumVector;
   
}

#endif
