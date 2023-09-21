import pandas as pd
import numpy as np
import os

currentDirectory = os.getcwd()
parentDirectory = os.path.dirname(currentDirectory)

lumFileName = "lum_perrun_RP_CPT2.txt"
inputFilesList = "starlist.list"

lumFilePath = parentDirectory + "/share/" + lumFileName
inputFilesListPath = parentDirectory + "/share/" + inputFilesList

columnNames = ["run", "start time", "end time", "fill", "lum", "prescale", "trigger lifetime", \
               "base lum trigger", "base lum trigger lifetime", "L times pblue2pluspyellow2/2", "L times pblue2pellow2", \
               "number of events", "number of events times p4", "number of events times pblue2pluspyellow2/2", "L times pblue2"]

data = pd.read_csv(lumFilePath, delimiter = " ", header = None)
data.dropna(axis = 1, how = "all", inplace = True)
data.columns = columnNames
runsFromFile = data["run"]

inputFiles =  pd.read_csv(inputFilesListPath, delimiter = " ", header = None)
inputFiles[0] = inputFiles[0].str.replace(r'\.root$', '', regex=True)

inputFiles[0] = inputFiles[0].astype(int)

dataWithUsedRuns =  data[data["run"].isin(inputFiles[0])]
integratedLuminosity = dataWithUsedRuns["lum"].sum()
print("intergrated luminosity: ", integratedLuminosity)


