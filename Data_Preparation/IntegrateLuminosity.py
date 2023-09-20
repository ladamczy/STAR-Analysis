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

integratedLuminosity = data["lum"].sum()
runsFromFile = data["run"]

print("intergrated luminosity: ", integratedLuminosity)


inputFiles =  pd.read_csv(inputFilesListPath, delimiter = " ", header = None)
inputFiles[0] = inputFiles[0].str.replace(r'\.root$', '', regex=True)

inputFiles[0] = inputFiles[0].astype(int)
elementsNotInInputFiles =  data["run"][~data["run"].isin( inputFiles[0])]
elementsNotInLumiFile =  inputFiles[0][~inputFiles[0].isin( data["run"])]

print(elementsNotInInputFiles)
print(elementsNotInLumiFile)
