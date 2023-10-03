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
	
	int nFill;
	double dXPosition, dYPosition, dX2Position, dY2Position;	
	
	for (unsigned long int i = 0; i < Data.size(); i ++)
	{
		stringstream stringSteramFill, stringSteramX, stringSteramY, stringSteramX2, stringSteramY2;
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
	
		vFillNumber.push_back(nFill);
		vXPosition.push_back(dXPosition);
		vYPosition.push_back(dYPosition);
		vX2Position.push_back(dX2Position);
		vY2Position.push_back(dY2Position);
	}			

	vector <vector<double>> AllData;
	AllData.push_back(vFillNumber);
 	AllData.push_back(vXPosition);
 	AllData.push_back(vYPosition);
	AllData.push_back(vX2Position);
 	AllData.push_back(vY2Position);
 	
 	return 	AllData;	
}
