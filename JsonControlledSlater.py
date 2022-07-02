"""A means of accessing the functionality of SlaterGUI.py without a GUI.
   By Daniel Isenberg
   Stuff to do: Consider refactoring the input orbital configuration for atoms so that they are integer lists instead of strings to save an extra step in converting it. (It only took string input originally because of the manual input system I started with)
                Consider refactoring Slater.density() so that it only returns the component contributions, which can then be added together later so it only has to return one thing instead of a tuple.
                Refactor pretty much everything here so it works with all the changes made to regular Slater, ie usage of atom object, etc.
"""
import Slater
import inputFunctions
import numpy
import matplotlib.pyplot as plt
import FileIO


json = FileIO.openJsonFile("C:/Users/Daniel/Desktop/inputJson.json")
serializedJson = FileIO.serializeJson(json)


labelList = ["cumulative density", "1s subshell", "2s&p subshell", "3s&p subshell", "3d subshell", "4s&p subshell", "4d subshell", "4f subshell", "5s&p subshell", "5d subshell"]
run = True

atomList = []
atomicNumberList = []
shellOccupationList = []
dtyList = []

for atom in range(len(serializedJson["atoms"])):
    plotType = serializedJson["plotType"]
    derivativeNumber = serializedJson["derivativeNumber"]
    scaleType = serializedJson["scaleType"]
    arrayX = inputFunctions.getArrayXFromJSON(scaleType, serializedJson["plotRadius"])
    atomicNumber = serializedJson["atoms"][atom]["atomicNumber"]
    print("Atom index is " + str(atom))
    print("Atomic Number for current atom is " + str(atomicNumber))
    orbitalConfigList = inputFunctions.getElectronConfigFromJson(serializedJson["atoms"][atom]["shellOccupation"])

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

    atomicNumberList.append(atomicNumber)
    shellOccupationList.append(orbitalConfigList)
    dtyList.append(dty)
for i in range(len(atomicNumberList)):
    atomList.append(
        {
            "atomicNumber": atomicNumberList[i],
            "shellOccupation": shellOccupationList[i],
            "dtyRaw": dtyList[i].tolist()
        }
    )


encodedNumpyData = FileIO.saveToJson(atomList)

outputListTest = []
outTest = FileIO.openJsonFile("C:/Users/Daniel/Desktop/outputJson.json")
serializedOutTest = FileIO.serializeJson(outTest)
for i in range(len(serializedOutTest["atoms"])):
    testListY = serializedOutTest["atoms"][i]["dtyRaw"][0]
    testListX = inputFunctions.getArrayXFromJSON(scaleType, serializedJson["plotRadius"])
    plt.plot(testListX, testListY, label="cumulative density")
    plt.title("Cumulative density for atomic number " + str(serializedOutTest["atoms"][i]["atomicNumber"]) + " retrieved from saved json file")
    plt.xlabel("Distance from atomic center (Bohr radii)")
    plt.ylabel("Radial electron density")
    plt.legend()
    plt.grid()
    plt.show()
