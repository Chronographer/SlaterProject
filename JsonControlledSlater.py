"""A means of accessing the functionality of SlaterGUI.py without a GUI.
   By Daniel Isenberg
"""
import Slater
import inputFunctions
import numpy
import matplotlib.pyplot as plt
import FileIO


json = FileIO.openJsonFile()
serializedJson = FileIO.serializeJson(json)


labelList = ["cumulative density", "1s subshell", "2s&p subshell", "3s&p subshell", "3d subshell", "4s&p subshell", "4d subshell", "4f subshell", "5s&p subshell", "5d subshell"]
run = True

while run:
    plotType = serializedJson["plotType"]
    derivativeNumber = serializedJson["derivativeNumber"]
    scaleType = serializedJson["scaleType"]
    arrayX = inputFunctions.getArrayXFromJSON(scaleType, serializedJson["plotRadius"])
    atomicNumber = serializedJson["atoms"][0]["atomicNumber"]
    orbitalConfigList = inputFunctions.getElectronConfigFromJson(serializedJson["atoms"][0]["shellOccupation"])

    dty, components = Slater.density(arrayX, atomicNumber, orbitalConfigList)
    # dty = Slater.grlaglll(arrayX, atomicNumber, orbitalConfigList)
    dty = 4 * numpy.pi * arrayX ** 2 * dty
    yList = []
    yListMaster = []

    if plotType == "cumulative" or plotType == "both":  # Plots the density of the entire system as a single data set.
        for i in range(len(arrayX)):
            yList.append(dty[derivativeNumber][i])
        yListMaster.append(yList)

    if plotType == "components" or plotType == "both":  # Plots the density of each individual shell without considering screening from other shells, each as its own data set.
        for index in range(len(components)):
            components[index][derivativeNumber] = 4 * numpy.pi * arrayX ** 2 * components[index][derivativeNumber]
            yListMaster.append(components[index][derivativeNumber])

    for i in range(len(yListMaster)):  # This ensures that labels in the legend are applied correctly by skipping the first label in the list of label handles if the cumulative plot is not present.
        if plotType == "cumulative" or plotType == "both":
            plt.plot(arrayX, yListMaster[i], label=labelList[i])
        else:
            plt.plot(arrayX, yListMaster[i], label=labelList[i + 1])

    if derivativeNumber != 0:  # makes the title reflect whether you are plotting just the density or one of it's derivatives.
        plt.title("Plot of charge density (" + inputFunctions.derivativeOptions[derivativeNumber] + ") vs. radius for atomic number " + str(atomicNumber) + "\nScale type: " + scaleType)
    else:
        plt.title("Plot of charge density vs. radius for atomic number " + str(atomicNumber) + "\nScale type: " + scaleType)
    plt.xlabel("Distance from atomic center (Bohr radii)")
    plt.ylabel("Radial electron density")
    plt.legend()
    plt.grid()
    plt.show()
    encodedNumpyData = FileIO.saveJsonTest(dty)
    FileIO.saveToJson(encodedNumpyData)

    run = inputFunctions.askToRepeat()
