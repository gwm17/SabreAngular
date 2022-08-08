#!/usr/bin/env python

from ast import Param
import numpy as np
from lmfit import Parameters, minimize, fit_report
import scipy.special as funcs
from matplotlib import pyplot as plt

def AngularDistribution(params : Parameters, x : np.ndarray, L : int=0, data : np.ndarray=None, err : np.ndarray=None) -> np.ndarray:
    vals = params.valuesdict()
    result = np.zeros(len(x))
    paramStr = ""
    for i in range(L+1):
        paramStr = "a" + str(i)
        result += vals[paramStr] * funcs.eval_legendre(i*2, x)

    if data is None:
        return result
    elif err is None:
        return result - data
    else:
        return (result-data)/err

def ReadDatafile(filename : str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    xPoints = np.empty(0)
    yPoints = np.empty(0)
    yErr = np.empty(0)
    with open(filename, "r") as inputfile:
        inputfile.readline()
        for line in inputfile.readlines():
            entries = line.split("\t")
            xPoints = np.append(xPoints, float(entries[0]))
            yPoints = np.append(yPoints, float(entries[1]))
            yErr = np.append(yErr, float(entries[2]))
    return xPoints, yPoints, yErr

def GenerateDistribution(experimentCounts : np.ndarray, simEfficiency : np.ndarray , expErr : np.ndarray,
                         simErr : np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    diffCrossSection = np.zeros(len(experimentCounts))
    diffCrossSectionErr = np.zeros(len(experimentCounts))
    for i in range(len(experimentCounts)):
        diffCrossSection[i] = experimentCounts[i] * 1.0/(simEfficiency[i]*4.0*np.pi)
        diffCrossSectionErr[i] = diffCrossSection[i] * np.sqrt((expErr[i]/experimentCounts[i])**2.0 + (simErr[i]/simEfficiency[i])**2.0 + 0.075**2.0)
    return diffCrossSection, diffCrossSectionErr

def GenerateParameters(L : int) -> Parameters:
    params = Parameters()
    paramStr = ""
    for i in range(L+1):
        paramStr = "a" + str(i)
        params.add(paramStr, value=1.0, min=0.0, max=10.0)
    return params

def RunFit(expfile : str, simfile : str, L : int):
    x, normedCounts, countsErr = ReadDatafile(expfile)
    junk, efficiency, effErr = ReadDatafile(simfile)

    diffCS, diffCSErr = GenerateDistribution(normedCounts, efficiency, countsErr, effErr)

    params = GenerateParameters(L)

    out = minimize(AngularDistribution, params, args=(x,), kws={"L": L, "data": diffCS, "err": diffCSErr})
    print(fit_report(out))
    branchingRatio = out.params["a0"] * 4.0 * np.pi
    print("Branching Ratio: ", branchingRatio)

    plt.errorbar(x, diffCS, yerr = diffCSErr, marker="o",label="Data")
    plt.plot(x, AngularDistribution(out.params, x=x), linestyle="-", label="best fit")
    plt.xlabel(r"cos($\theta_{CM}$)", size=12)
    plt.ylabel(r"$\frac{d\sigma}{d\Omega_{CM}}$", size=14)
    plt.legend()
    plt.tight_layout()
    plt.show()

datafile = "../data/experiment/experimentalData_7Begs_08082022.txt"
simfile = "../data/simulation/simulationData_7Begs_08082022.txt"
RunFit(datafile, simfile, L=0)
