"""A simple GUI for displaying Slater electron densities.
   By Neal Coleman
"""
import easygui
from easygui import *
import sys
import Slater
import numpy
import Routines
import matplotlib.pyplot as plt

listx = [0, 1, 2, 3, 4, 5, 6]
listy = [1, 5, 6, 2, 3, 8, 7]
list1x = [0, 1, 2, 3, 4, 5, 6]
list1y = [4, 5, 5, 2, 1, 7, 8]
list2x = [0, 1, 2, 3, 4, 5, 6]
list2y = [6, 4, 1, 1, 5, 9, 9]

mlistx = []
mlisty = []

mlistx.append(listx)
mlisty.append(listy)
mlistx.append(list1x)
mlisty.append(list1y)
mlistx.append(list2x)
mlisty.append(list2y)

while 1:
    # msgbox("A simple GUI for displaying Slater electron densities\n" + "By Neal Coleman", "Introduction")

    Z = integerbox("Enter the atomic number of the element.", "Element Selection.")

    listLength = integerbox("Enter (integer) length of grid.", "Grid Length")

    exponentialOrUniform = ["Exponential", "Uniform"]

    scaleType = buttonbox("Would you like an exponential or uniform grid? (They produce identical plots for now)", "Grid Type", exponentialOrUniform)

    plotType = enterbox("What would you like to plot? (type 'combined' for combined only, 'components' for subshell components only, or 'both' for both.)")
    while not(plotType == "components" or plotType == "combined" or plotType == "both"):
        msgbox("'" + plotType + "' is not a valid plot type. Please enter only one of the specified options.")
        plotType = enterbox("What would you like to plot? (type 'combined' for combined only, 'components' for subshell components only, or 'both' for both.)")

    """legendLabelList = [] # for eventual automatic creation of legend labels.
    for legendIndex in range(independentOrbitalNumber):
        legendBox = easygui.textbox("Enter the legend label for the plot number" + str(legendIndex))
        legendLabelList.append(legendBox)"""

    if scaleType == "Exponential":
        arrayx = Routines.ExpGridStretch2(numpy.arange(0.01, 1.0 * listLength, 0.01))
    else:
        arrayx = numpy.arange(0.01, 1.0 * listLength, 0.01)

    subshell = ["1s", "2s&p", "3s&p", "3d", "4s&p", "4d", "4f", "5s&p", "5d"]

    xListMaster = []
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
    xList = []
    yList = []

    if plotType == "combined" or plotType == "both":
        for i in range(len(arrayx)):
            xList.append(arrayx[i])
            yList.append(dty[i])
        xListMaster.append(xList)
        yListMaster.append(yList)

    if plotType == "components" or plotType == "both":
        for index in range(len(components)):
            components[index] = 4 * numpy.pi * arrayx**2 * components[index]
            yListMaster.append(components[index])

    for i in range(len(yListMaster)):
        plt.plot(arrayx, yListMaster[i], label="plot number "+str(i))

    plt.title("Plot of <something> vs. radius for atomic number " + str(Z))
    plt.xlabel("radius")
    plt.ylabel("4pi r ^2")
    plt.legend()
    plt.show()

    if ccbox("There was your density.  Shall we do it again?", "Finale"):  # show a Continue/Cancel dialog
        pass  # user chose Continue
    else:
        sys.exit(0)  # user chose Cancel
