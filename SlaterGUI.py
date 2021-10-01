"""A simple GUI for displaying Slater electron densities.
   By Neal Coleman
"""

from easygui import *
import sys
import Slater
import numpy
import Routines
import matplotlib.pyplot as plt


while 1:
    msgbox("A simple GUI for displaying Slater electron densities\n" + "By Neal Coleman", "Introduction")

    Z = integerbox("Enter the atomic number of the element.", "Element Selection.")

    listLength = integerbox("Enter (integer) length of grid.", "Grid Length")

    choices1 = ["Exponential", "Uniform"]

    scaleType = buttonbox("Would you like an exponential or uniform grid?", "Grid Type", choices1)

    if scaleType == "Exponential":
        arrayx = Routines.ExpGridStretch2(numpy.arange(0.01, 1.0 * listLength, 0.01))
    else:
        arrayx = numpy.arange(0.01, 1.0 * listLength, 0.01)

    subshell = ["1s", "2s&p", "3s&p", "3d", "4s&p", "4d", "4f", "5s&p", "5d"]

    occupancy = []

    for i in range(9):
        occupancy.append(integerbox("Enter the occupancy of the " + subshell[i] + " subshell.", "Subshell Selection"))
    dty = 4 * numpy.pi * arrayx**2 * Slater.density(arrayx, Z, occupancy)[0]  # what is 'dty'? Density?
    xList = []
    yList = []
    for i in range(len(arrayx)):
        xList.append(arrayx[i])
        yList.append(dty[i])
    plt.plot(xList, yList)
    plt.title("Plot of <something> vs. <something> for atomic number " + str(Z))
    plt.xlabel("xLable")
    plt.ylabel("ylable")
    plt.show()

    if ccbox("There was your density.  Shall we do it again?", "Finale"):  # show a Continue/Cancel dialog
        pass  # user chose Continue
    else:
        sys.exit(0)  # user chose Cancel
