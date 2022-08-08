#!/usr/bin/env python

'''
GenerateExpFile

Use this script to generate a text file of experimental angular distribution data of the proper format for performing
angular distribution fits to extract branching ratios. Currently, just enter your experimental counts into binnedCounts,
the cosThetaCM coordinates of the bins into binCoords, and the total counts into totalCounts. Modify filename as needed. 

Intended for use with SPS-SABRE data at FSU, but should be broadly applicable.

GWM - Aug 2022
'''


import numpy as np

totalCounts = 129985.7
binnedCounts = np.array([168.859, 158.625, 220.097, 122.958])
binCoords = np.array([-0.95, -0.85, -0.75, -0.65])
countsErr = np.zeros(len(binnedCounts))
normedCounts = binnedCounts * 1.0/totalCounts

for i in range(len(binnedCounts)):
    countsErr[i] = normedCounts[i]*np.sqrt(1.0/binnedCounts[i] + 1.0/totalCounts)

filename = "../data/experiment/experimentalData_7Begs_08082022.txt"
with open(filename, "w") as outputfile:
    line = "Cos(ThetaCM)\tNormedCounts\tErr\n"
    outputfile.write(line)

    for i in range(len(normedCounts)):
        line = str(binCoords[i]) + "\t" + str(normedCounts[i]) + "\t" + str(countsErr[i]) + "\n"
        outputfile.write(line)
    outputfile.close()