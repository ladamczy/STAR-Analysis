#include "MatchFillPosition.h"
#include <iostream>
#include <string>    
#include <sstream> 
#include <vector> 
#include <cmath>
#include <fstream> 
#include <ostream>

using namespace std;
vector <double> FindPosition(int nFillNumber, double zPos, vector <double> &vFillNumber,  vector <double> &vXPosition, vector <double> &vYPosition,   vector <double> &vX2Position, vector <double> &vY2Position, vector <double> &vXSlope, vector <double> &vYSlope,   vector <double> &vX2Slope, vector <double> &vY2Slope)
{

	bool found = false;
	int indexFill;

	vector <double> vFinalPositions;
	
	for (unsigned int i = 0; i < vFillNumber.size(); i++)
	{
		if (nFillNumber == vFillNumber[i] )
		{
			indexFill = i;	
			found = true;
		}
	}
			
	if (found == false)
	{
		vFinalPositions.push_back(nan(""));
		vFinalPositions.push_back(nan(""));
		return vFinalPositions;
		
	}
						
	double posX, posY, posX2, posY2, beamPositionX, beamPositionY;
	double sX, sY, sX2, sY2;
	posX = vXPosition[indexFill];
	posY = vYPosition[indexFill];
	posX2 = vX2Position[indexFill];
	posY2 = vY2Position[indexFill];		 

	sX = vXSlope[indexFill];
	sY = vYSlope[indexFill];
	sX2 = vX2Slope[indexFill];
	sY2 = vY2Slope[indexFill];		 

	if (isnan(posX) or isnan(posY))
	{
		vFinalPositions.push_back(nan(""));
		vFinalPositions.push_back(nan(""));
		return vFinalPositions;
	}
		
	if (isnan(posX2) or isnan(posY2))
	{
		beamPositionX = posX + zPos*sX;
		beamPositionY = posY + zPos*sY;	
	}
			
	else
	{
		beamPositionX = ( ( posX + zPos*sX )+ ( posX2 + zPos*sX2) )/2.0;
		beamPositionY = ( ( posY + zPos*sY )+ ( posY2 + zPos*sY2) )/2.0;			
	}
			

	vFinalPositions.push_back(beamPositionX);
	vFinalPositions.push_back(beamPositionY);
	
	return vFinalPositions;

}
