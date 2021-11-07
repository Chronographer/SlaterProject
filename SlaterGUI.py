"""A simple GUI for displaying Slater electron densities.
   By Neal Coleman
"""
from easygui import *
import sys
import Slater
import numpy
import Routines
import matplotlib.pyplot as plt

labelList = ["cumulative density", "1s subshell", "2s&p subshell", "3s&p subshell", "3d subshell", "4s&p subshell", "4d subshell", "4f subshell", "5s&p subshell", "5d subshell"]

while 1:
    # msgbox("A simple GUI for displaying Slater electron densities\n" + "By Neal Coleman", "Introduction")

    Z = integerbox("Enter the atomic number of the element.", "Element Selection.")

    listLength = integerbox("Enter (integer) length of grid.", "Grid Length")

    exponentialOrUniform = ["Exponential", "Uniform"]

    scaleType = buttonbox("Would you like an exponential or uniform grid? (They produce identical plots for now)", "Grid Type", exponentialOrUniform)

    plotType = enterbox("What would you like to plot? (type 'cumulative' for cumulative only, 'components' for subshell components only, or 'both' for both.)")
    while not(plotType == "components" or plotType == "cumulative" or plotType == "both"):
        msgbox("'" + plotType + "' is not a valid plot type. Please enter only one of the specified options.")
        plotType = enterbox("What would you like to plot? (type 'cumulative' for cumulative only, 'components' for subshell components only, or 'both' for both.)")

    if scaleType == "Exponential":
        arrayx = Routines.ExpGridStretch2(numpy.arange(0.01, 1.0 * listLength, 0.01))
    else:
        arrayx = numpy.arange(0.01, 1.0 * listLength, 0.01)

    subshell = ["1s", "2s&p", "3s&p", "3d", "4s&p", "4d", "4f", "5s&p", "5d"]

    yListMaster = []
    orbitalConfigList = []
    while len(orbitalConfigList) != 9:
        orbitalConfigString = enterbox("Enter the orbital configuration of each shell as a period separated list.\n\nThere should be no more than 9 elements.\n\nTrailing zero's can be omitted and will be automatically added to the end of the list as needed.")
        if not orbitalConfigString[len(orbitalConfigString)-1].isdigit():
            print("\nWARNING: the last character of your input string was '" + orbitalConfigString[len(orbitalConfigString)-1] + "', which is not a digit. You might have accidentally left out a number or entered it incorrectly. The input you provided was: '" + orbitalConfigString + "'\nThe anomalous character has been removed in an attempt to prevent a runtime error.\nBe cautious; even if a runtime error is not thrown, your results might not reflect the system you intended to model!\n")
            orbitalConfigString = orbitalConfigString[:-1]  # this removes the last character from the input string; the idea is that if you accidentally end with a '.', this will remove it for you. Will produce incorrect results if the user left out a non-zero number after the last character.
        orbitalConfigList = orbitalConfigString.split(".")
        if len(orbitalConfigList) > 9:
            msgbox("There are " + str(len(orbitalConfigList)) + " elements in orbitalConfigList, there should be no more than 9.\nRe-enter the orbital configuration of each shell as a period separated list.")
        elif len(orbitalConfigList) < 9:
            while len(orbitalConfigList) < 9:  # Allows user to omit all trailing 0's by adding as many as are needed to generate a list of correct length, thus saving time.
                orbitalConfigList.append(0)
    for element in range(0, len(orbitalConfigList)):
        orbitalConfigList[element] = int(orbitalConfigList[element])

    dty, components = Slater.density(arrayx, Z, orbitalConfigList)
    dty = 4 * numpy.pi * arrayx**2 * dty
    yList = []

    if plotType == "cumulative" or plotType == "both":
        for i in range(len(arrayx)):
            yList.append(dty[i])
        yListMaster.append(yList)

    if plotType == "components" or plotType == "both":
        for index in range(len(components)):
            components[index] = 4 * numpy.pi * arrayx**2 * components[index]
            yListMaster.append(components[index])

    for i in range(len(yListMaster)):
        if plotType == "combined" or plotType == "both":
            plt.plot(arrayx, yListMaster[i], label=labelList[i])
        else:
            plt.plot(arrayx, yListMaster[i], label=labelList[i+1])

    plt.title("Plot of <something> vs. radius for atomic number " + str(Z))
    plt.xlabel("radius")
    plt.ylabel("4pi r ^2")
    plt.legend()
    plt.show()

    if ccbox("There was your density.  Shall we do it again?", "Finale"):  # show a Continue/Cancel dialog
        pass  # user chose Continue
    else:
        sys.exit(0)  # user chose Cancel
