#include "MatchFillPosition.h"
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
