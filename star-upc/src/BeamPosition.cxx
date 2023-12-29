#include "BeamPosition.h"
#include <iostream>
#include <string>    
#include <sstream> 
#include <vector> 
#include <cmath>
#include <fstream> 
#include <ostream>

using namespace std;

vector <vector<double>> ReadFillPositionData(string file)
{
	ifstream FileName(file);
	string line;
	vector <vector<string>> Data;
	while (getline(FileName, line))
	{
		vector <string> data;
		string str;
		stringstream lines(line);
		while(getline(lines,str,','))
		{		
			data.push_back(str);			
		}	
		Data.push_back(data);
	}
	
	vector <double> vFillNumber;
	vector <double> vXPosition;
	vector <double> vYPosition;
	vector <double> vX2Position;
	vector <double> vY2Position; 

	vector <double> vXSlope;
	vector <double> vYSlope;
	vector <double> vX2Slope;
	vector <double> vY2Slope; 

	int nFill;
	double dXPosition, dYPosition, dX2Position, dY2Position;	
	double dXSlope, dYSlope, dX2Slope, dY2Slope;	

	for (unsigned long int i = 0; i < Data.size(); i ++)
	{
		stringstream stringSteramFill, stringSteramX, stringSteramY, stringSteramX2, stringSteramY2,  stringSteramXslope, stringSteramYslope, stringSteramX2slope, stringSteramY2slope;
		if (Data[i][1] != "")
		{
			stringSteramFill << Data[i][1]; 
			stringSteramFill >> nFill;
		}
	
		if (Data[i][2] != "")
		{
			stringSteramX << Data[i][2];
			stringSteramX >> dXPosition;
		}
		else dXPosition = nan("");

		if (Data[i][3] != "")
		{		
			stringSteramY << Data[i][3];
			stringSteramY >> dYPosition;	
		}
		else dYPosition = nan("");
		
		if (Data[i][4] != "")
		{
			stringSteramX2 << Data[i][4];
			stringSteramX2 >> dX2Position;
		}
		else dX2Position = nan("");

		if (Data[i][5] != "")
		{		
			stringSteramY2 << Data[i][5];
			stringSteramY2 >> dY2Position;	
		}
		else dY2Position = nan("");
	
		if (Data[i][6] != "")
		{		
			stringSteramXslope << Data[i][6];
			stringSteramXslope >> dXSlope;	
		}
		else dXSlope = nan("");
	
		if (Data[i][7] != "")
		{		
			stringSteramYslope << Data[i][7];
			stringSteramYslope >> dYSlope;	
		}
		else dYSlope = nan("");
	
	
		if (Data[i][8] != "")
		{		
			stringSteramX2slope << Data[i][8];
			stringSteramX2slope >> dX2Slope;	
		}
		else dX2Slope = nan("");
	
		if (Data[i][9] != "")
		{		
			stringSteramY2slope << Data[i][9];
			stringSteramY2slope >> dY2Slope;	
		}
		else dY2Slope = nan("");
		
		vFillNumber.push_back(nFill);
		vXPosition.push_back(dXPosition);
		vYPosition.push_back(dYPosition);
		vX2Position.push_back(dX2Position);
		vY2Position.push_back(dY2Position);

		vXSlope.push_back(dXSlope);
		vYSlope.push_back(dYSlope);
		vX2Slope.push_back(dX2Slope);
		vY2Slope.push_back(dY2Slope);
	}			

	vector <vector<double>> AllData;
	AllData.push_back(vFillNumber);
 	AllData.push_back(vXPosition);
 	AllData.push_back(vYPosition);
	AllData.push_back(vX2Position);
 	AllData.push_back(vY2Position);

 	AllData.push_back(vXSlope);
 	AllData.push_back(vYSlope);
	AllData.push_back(vX2Slope);
 	AllData.push_back(vY2Slope);

 	return 	AllData;	
}


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
		vFinalPositions.push_back(nan(""));
		vFinalPositions.push_back(nan(""));
		return vFinalPositions;
		
	}
						
	double posX, posY, posX2, posY2, beamPositionX, beamPositionY, slopeX, slopeY;
	double sX, sY, sX2, sY2;
	posX = vXPosition[indexFill];
	posY = vYPosition[indexFill];
	posX2 = vX2Position[indexFill];
	posY2 = vY2Position[indexFill];		 

	sX = vXSlope[indexFill];
	sY = vYSlope[indexFill];
	sX2 = vX2Slope[indexFill];
	sY2 = vY2Slope[indexFill];		 

	if (std::isnan(posX) or std::isnan(posY))
	{
		vFinalPositions.push_back(nan(""));
		vFinalPositions.push_back(nan(""));
		vFinalPositions.push_back(nan(""));
		vFinalPositions.push_back(nan(""));
		return vFinalPositions;
	}
		
	if (std::isnan(posX2) or std::isnan(posY2))
	{
		beamPositionX = posX + zPos*sX;
		beamPositionY = posY+zPos*sY;
		slopeX = sX;
		slopeY = sY;
	}
			
	else
	{
		beamPositionX = ( ( posX + zPos*sX )+ ( posX2 + zPos*sX2) )/2.0;
		beamPositionY = ((posY+zPos*sY)+(posY2+zPos*sY2))/2.0;
		slopeX = (sX+sX2)/2.0;
		slopeY = (sY+sY2)/2.0;
	}
			

	vFinalPositions.push_back(beamPositionX);
	vFinalPositions.push_back(beamPositionY);
	vFinalPositions.push_back(slopeX);
	vFinalPositions.push_back(slopeY);
	
	return vFinalPositions;

}



vector <double> GetBeamParameters(int nFillNumber, vector <double> &vFillNumber,  vector <double> &vXPosition, vector <double> &vYPosition,   vector <double> &vX2Position, vector <double> &vY2Position, vector <double> &vXSlope, vector <double> &vYSlope,   vector <double> &vX2Slope, vector <double> &vY2Slope)
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
						
	double posX, posY, posX2, posY2;
	double sX, sY, sX2, sY2;
	posX = vXPosition[indexFill];
	posY = vYPosition[indexFill];
	posX2 = vX2Position[indexFill];
	posY2 = vY2Position[indexFill];		 

	sX = vXSlope[indexFill];
	sY = vYSlope[indexFill];
	sX2 = vX2Slope[indexFill];
	sY2 = vY2Slope[indexFill];		 

	if (std::isnan(posX) or std::isnan(posY))
	{
		vFinalPositions.push_back(nan(""));
		vFinalPositions.push_back(nan(""));
		return vFinalPositions;
	}
		
	if (std::isnan(posX2) or std::isnan(posY2))
	{
		vFinalPositions.push_back(posX);
		vFinalPositions.push_back(posY);
		vFinalPositions.push_back(sX);
		vFinalPositions.push_back(sY);
	}
			
	else
	{
		vFinalPositions.push_back((posX+posX2)/2.0);
		vFinalPositions.push_back((posY+posY2)/2.0);
		vFinalPositions.push_back((sX+sX2)/2.0);
		vFinalPositions.push_back((sY+sY2)/2.0);
	}
	
	

	
	return vFinalPositions;

}