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
run = True


while run:
    target = inputFunctions.getElementNameInput()

    if target == "manual":
        atomicNumber = inputFunctions.getAtomicNumber()
        electronOccupancy = inputFunctions.getElectronConfigInput()
        elementName = "placeholder element name"
        atom = Atoms.NewAtom(atomicNumber, elementName, electronOccupancy)
    else:
        elementName = target
        atomicNumber = Atoms.atomData[target]['atomicNumber']
        electronOccupancy = Atoms.atomData[target]['occupancy']
        atom = Atoms.NewAtom(atomicNumber, elementName, electronOccupancy)

    plotType = inputFunctions.getPlotType()
    derivativeNumber = inputFunctions.chooseDerivativeOptions()
    scaleType = inputFunctions.getScaleType()
    arrayX = inputFunctions.getArrayX(scaleType)
    inputFunctions.showEnergy(atom.totalEnergy)

    if target == "manual":
        dty, components = Slater.newDensity(arrayX, atom)
    else:
        dty, components = Slater.newDensity(arrayX, atom)

    dty = 4 * numpy.pi * arrayX ** 2 * dty
    yList = []
    yListMaster = []

    if plotType == "cumulative" or plotType == "both":  # Plots the combined density of every shell in the system as a single data set, without showing the contributions of each shell individual .
        for i in range(len(arrayX)):
            yList.append(dty[derivativeNumber][i])
        yListMaster.append(yList)

    if plotType == "components" or plotType == "both":  # Plots the contribution to the density from each individual shell as its own data set.
        for index in range(len(components)):
            components[index][derivativeNumber] = 4 * numpy.pi * arrayX ** 2 * components[index][derivativeNumber]
            yListMaster.append(components[index][derivativeNumber])

    for i in range(len(yListMaster)):  # This loop adds the desired plots to the output graph, and ensures they have the appropriate label in the legend.
        if plotType == "cumulative" or plotType == "both":
            if i == 0:
                plt.plot(arrayX, yListMaster[i], label="cumulative density")  # Since the cumulative density is the first thing to be plotted, we want that to be the first label. We insert it manually as a special case because that string can't be constructed from stuff already in the atom object.
            else:
                legendLabel = str(atom.principalQuantumNumberLabelList[i-1]) + atom.azimuthalQuantumNumberLabelList[i-1] + " subshell"  # After the first label, we can construct the "<PrincipalQuantumNumber><AzimuthalQuantumNumber> subshell" string/label from data stored in the atom itself. (ie, "2sp subshell", "4d subshell", etc) but we have to use [i-1] because the first label is taken by the cumulative plot.
                orbitalIsEmptyFlag = True
                for index in range(len(yListMaster[i])):
                    if not yListMaster[i][index] == 0:
                        orbitalIsEmptyFlag = False
                        break
                if orbitalIsEmptyFlag:
                    legendLabel = legendLabel + " (empty)"
                plt.plot(arrayX, yListMaster[i], label=legendLabel)
        else:
            legendLabel = str(atom.principalQuantumNumberLabelList[i]) + atom.azimuthalQuantumNumberLabelList[i] + " subshell"  # This does not need the [i-1] part, because the cumulative density is not being plotted, and so is not being labeled. Thus we can generate ALL labels from the atom object itself.
            orbitalIsEmptyFlag = True
            for index in range(len(yListMaster[i])):
                if not yListMaster[i][index] == 0:
                    orbitalIsEmptyFlag = False
                    break
            if orbitalIsEmptyFlag:
                legendLabel = legendLabel + " (empty)"
            plt.plot(arrayX, yListMaster[i], label=legendLabel)

    if derivativeNumber != 0:  # This makes the title reflect whether you are plotting just the density or one of it's derivatives.
        plt.title("Plot of charge density (" + inputFunctions.derivativeOptions[derivativeNumber] + ") vs. radius for atomic number " + str(atomicNumber) + "\nScale type: " + scaleType)
    else:
        plt.title("Plot of charge density vs. radius for atomic number " + str(atomicNumber) + "\nScale type: " + scaleType)
    plt.xlabel("Distance from atomic center (Bohr radii)")
    plt.ylabel("Radial electron density")
    plt.legend()
    plt.grid()
    plt.show()

    run = inputFunctions.askToRepeat()
