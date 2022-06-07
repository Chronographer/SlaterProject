"""A simple GUI for displaying Slater electron densities.
   By Neal Coleman
   Modified and expanded by Daniel Isenberg
"""
import Slater
import inputFunctions
import numpy
import matplotlib.pyplot as plt
import Atoms

# need to do: add documentation about how to use stuff and what it means, look into pypi stuff
labelList = ["cumulative density", "1s subshell", "2s&p subshell", "3s&p subshell", "3d subshell", "4s&p subshell", "4d subshell", "4f subshell", "5s&p subshell", "5d subshell"]  # this list of strings is used in the legend of matplotlib plots.
run = True


while run:
    target = inputFunctions.getElementNameInput()

    if target == "manual":
        atomicNumber = inputFunctions.getAtomicNumber()
    else:
        atomicNumber = Atoms.periodictable[target].Z

        #  make an actual atom object here
        elementName = target
        atomicNumber = Atoms.periodictable[target].Z
        electronOccupancy = Atoms.periodictable[target].occupancy
        atom = Atoms.NewAtom(atomicNumber, elementName, electronOccupancy)
        lists = Slater.newShieldingConstantComputer(atom)
        print("output of new computeShieldingConstants is: " + str(lists))

    plotType = inputFunctions.getPlotType()
    derivativeNumber = inputFunctions.chooseDerivativeOptions()
    scaleType = inputFunctions.getScaleType()
    arrayX = inputFunctions.getArrayX(scaleType)

    if target == "manual":
        orbitalConfigList = inputFunctions.getElectronConfigInput()
        dty, components = Slater.density(arrayX, atomicNumber, orbitalConfigList)
    else:
        orbitalConfigList = Atoms.periodictable[target].occupancy
        dty, components = Slater.newDensity(arrayX, atom)


    # dty, components = Slater.newDensity(atom)
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

    run = inputFunctions.askToRepeat()
