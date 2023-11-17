#ifndef BeamPosition_h
#define BeamPosition_h

#include <string>
#include <vector> 
using namespace std;

vector <vector<double>> ReadFillPositionData(string file);
vector <double> FindPosition(int nFillNumber, double zPos, vector <double> &vFillNumber,  vector <double> &vXPosition, vector <double> &vYPosition,  vector <double> &vX2Position, vector <double> &vY2Position, vector <double> &vXSlope, vector <double> &vYSlope,   vector <double> &vX2Slope, vector <double> &vY2Slope);

#endif

