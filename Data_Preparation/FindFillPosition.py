import numpy as np 
import pandas as pd
import datetime
import time
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

pd.set_option("display.max_rows", None, "display.max_columns", None)
pd.set_option('display.precision', 10)

BeamLineData = pd.read_csv("beamLine.txt", header = None, delim_whitespace=True)
Run7PolarizationData = pd.read_csv("Run7Polarization.txt", delim_whitespace=True)

BeamLineData["EpochDate"] = BeamLineData[0].astype(str) + " " +  BeamLineData[1]
BeamLineData["UnixTime"] = BeamLineData["EpochDate"].apply(lambda x: pd.Timestamp(x).timestamp()) 

BeamLineData["UnixTime"] = pd.to_numeric(BeamLineData["UnixTime"])
Run7PolarizationData["Start"] = pd.to_numeric(Run7PolarizationData["Start"])
Run7PolarizationData["Stop"] = pd.to_numeric(Run7PolarizationData["Stop"])

Run7PolarizationData["DateTimeStart"] = Run7PolarizationData["Start"].apply(lambda x: datetime.datetime.utcfromtimestamp(x))
Run7PolarizationData["DateTimeStop"] = Run7PolarizationData["Stop"].apply(lambda x: datetime.datetime.utcfromtimestamp(x))

beamLineUnixDate = BeamLineData["UnixTime"]
beamLineDate = BeamLineData["EpochDate"]

x = []
y = []
x2 = []
y2 = []

for ind in Run7PolarizationData.index:
				
	start = Run7PolarizationData["Start"][ind]
	stop = Run7PolarizationData["Stop"][ind]
				
	diffStart = (beamLineUnixDate - start)
	diffStop  = (beamLineUnixDate - stop)
				
	diffStart = diffStart.index[diffStart <= 0].tolist()
	diffStop = diffStop.index[diffStop >= 0].tolist()

	if len(diffStart)== 0 or len(diffStart)==0:
		x.append(np.nan) 
		y.append(np.nan)
		x2.append(np.nan)
		y2.append(np.nan)	
	
	else:
		indexStart = diffStart[-1]
		indexStop = diffStop[0]

		x.append(BeamLineData[2][indexStart]) 
		y.append(BeamLineData[4][indexStart])
		
	
		if (indexStop-indexStart) == 2:
			x2.append(BeamLineData[2][indexStart+1])
			y2.append(BeamLineData[4][indexStart+1])
		else:
			x2.append(np.nan)
			y2.append(np.nan)			


Run7PolarizationData["x"]  = x
Run7PolarizationData["y"]  = y
Run7PolarizationData["x2"]  = x2
Run7PolarizationData["y2"]  = y2

Run7PolarizationData2 = Run7PolarizationData[[ "Fill","x", "y", "x2", "y2", "Start", "Stop", "DateTimeStart",  "DateTimeStop"]].copy()
Run7PolarizationData2.to_csv("Run7PolarizationWithPosition.csv")	


