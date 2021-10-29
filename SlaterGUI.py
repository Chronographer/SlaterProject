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


while 1:
    # msgbox("A simple GUI for displaying Slater electron densities\n" + "By Neal Coleman", "Introduction")

    Z = integerbox("Enter the atomic number of the element.", "Element Selection.")

    listLength = integerbox("Enter (integer) length of grid.", "Grid Length")

    exponentialOrUniform = ["Exponential", "Uniform"]

    scaleType = buttonbox("Would you like an exponential or uniform grid? (They produce identical plots for now)", "Grid Type", exponentialOrUniform)

    independentOrbitalNumber = integerbox("How many orbitals would you like to plot independently?")

    """legendLabelList = [] # for eventual automatic creation of legend labels.
    for legendIndex in range(independentOrbitalNumber):
        legendBox = easygui.textbox("Enter the legend label for the plot number" + str(legendIndex))
        legendLabelList.append(legendBox)"""


    if scaleType == "Exponential":
        arrayx = Routines.ExpGridStretch2(numpy.arange(0.01, 1.0 * listLength, 0.01))
    else:
        arrayx = numpy.arange(0.01, 1.0 * listLength, 0.01)

    subshell = ["1s", "2s&p", "3s&p", "3d", "4s&p", "4d", "4f", "5s&p", "5d"]

    occupancy = []
    xListMaster = []
    yListMaster = []
    orbitalConfig = ""
    orbitalConfig = enterbox("Enter the orbital configuration of each shell as a comma separated list.")
    print(orbitalConfig)
    orbitalConfigList = orbitalConfig.split(",")
    print(orbitalConfigList)
    for element in range(0, len(orbitalConfigList)):
        orbitalConfigList[element] = int(orbitalConfigList[element])
    print(orbitalConfigList)
    for orbital in range(independentOrbitalNumber):
        msgbox("Enter the Orbital Configuration for plot number " + str(orbital) + ".")
        for i in range(9):
            occupancy.append(integerbox("Enter the occupancy of the " + subshell[i] + " subshell.", "Subshell Selection"))
        dty = 4 * numpy.pi * arrayx**2 * Slater.density(arrayx, Z, occupancy)[0]
        xList = []
        yList = []
        for i in range(len(arrayx)):
            xList.append(arrayx[i])
            yList.append(dty[i])
        xListMaster.append(xList)
        yListMaster.append(yList)
        occupancy.clear()
    for plotNumber in range(len(xListMaster)):
        plt.plot(xListMaster[plotNumber], yListMaster[plotNumber])
    plt.title("Plot of <something> vs. radius for atomic number " + str(Z))
    plt.xlabel("radius")
    plt.ylabel("4pi r ^2")
    plt.plot()
    plt.legend()
    plt.show()

    if ccbox("There was your density.  Shall we do it again?", "Finale"):  # show a Continue/Cancel dialog
        pass  # user chose Continue
    else:
        sys.exit(0)  # user chose Cancel
