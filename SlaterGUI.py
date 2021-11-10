"""A simple GUI for displaying Slater electron densities.
   By Neal Coleman
"""
from easygui import *
import sys
import Slater
import numpy
import inputFunctions
import matplotlib.pyplot as plt


labelList = ["cumulative density", "1s subshell", "2s&p subshell", "3s&p subshell", "3d subshell", "4s&p subshell", "4d subshell", "4f subshell", "5s&p subshell", "5d subshell"]

while 1:
    # msgbox("A simple GUI for displaying Slater electron densities\n" + "By Neal Coleman", "Introduction")

    atomicNumber = inputFunctions.getAtomicNumber()
    plotType = inputFunctions.getPlotType()
    derivativeNumber = inputFunctions.chooseDerivativeOptions()
    scaleType = inputFunctions.getScaleType()
    arrayX = inputFunctions.getArrayX(scaleType)
    orbitalConfigList = inputFunctions.getElectronConfigInput()

    yListMaster = []

    dty, components = Slater.density(arrayX, atomicNumber, orbitalConfigList)
    dty = 4 * numpy.pi * arrayX ** 2 * dty
    yList = []

    if plotType == "cumulative" or plotType == "both":  # Plots the density of the entire system as one plot.
        for i in range(len(arrayX)):
            yList.append(dty[derivativeNumber][i])
        yListMaster.append(yList)

    if plotType == "components" or plotType == "both":  # Plots the density of each individual shell without considering screening from other shells.
        for index in range(len(components)):
            components[index][derivativeNumber] = 4 * numpy.pi * arrayX ** 2 * components[index][derivativeNumber]
            yListMaster.append(components[index][derivativeNumber])

    for i in range(len(yListMaster)):
        if plotType == "cumulative" or plotType == "both":
            plt.plot(arrayX, yListMaster[i], label=labelList[i])
        else:
            plt.plot(arrayX, yListMaster[i], label=labelList[i + 1])  # This skips the first label in the list of label handles so the legend correctly labels the plots when the cumulative plot is not present.

    if derivativeNumber != 0:
        plt.title("Plot of charge density (" + inputFunctions.derivativeOptions[derivativeNumber] + ") vs. radius for atomic number " + str(atomicNumber) + "\nScale type: " + scaleType)
    else:
        plt.title("Plot of charge density vs. radius for atomic number " + str(atomicNumber) + "\nScale type: " + scaleType)
    plt.xlabel("radius")
    plt.ylabel("charge density 4pi r ^2")
    plt.legend()
    plt.grid()
    plt.show()

    if ccbox("There was your density.  Shall we do it again?", "Finale"):  # show a Continue/Cancel dialog
        pass  # user chose Continue
    else:
        sys.exit(0)  # user chose Cancel
