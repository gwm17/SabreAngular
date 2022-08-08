#!/usr/bin/env python

'''
GenerateSimFile

Use this script to generate a text file of simulated angular distribution data of the proper format for performing
angular distribution fits to extract branching ratios. Currently, just enter your simulated counts into binnedCounts,
the cosThetaCM coordinates of the bins into binCoords, and the total simulated samples into totalSamples. Modify filename as needed. 

Intended for use with SPS-SABRE data at FSU, but should be broadly applicable.

GWM - Aug 2022
'''

import numpy as np

totalSamples = 1.0e6
binnedCounts = np.array([41227.0, 32970.0, 27351.0, 21461.0])
binCoords = np.array([0.95, 0.85, 0.75, 0.65])
countsErr = np.zeros(len(binnedCounts))

efficiency = binnedCounts/totalSamples

for i in range(len(binnedCounts)):
    countsErr[i] = efficiency[i]*np.sqrt(1.0/binnedCounts[i] + 1.0/totalSamples)


filename = "../data/simulation/simulationData_7Begs_08082022.txt"
with open(filename, "w") as outputfile:
    line = "Cos(ThetaCM)\tEfficiency\tErr\n"
    outputfile.write(line)
    for i in range(len(binnedCounts)):
        line = str(binCoords[i]) + "\t" + str(efficiency[i]) + "\t" + str(countsErr[i]) + "\n"
        outputfile.write(line)
    outputfile.close()