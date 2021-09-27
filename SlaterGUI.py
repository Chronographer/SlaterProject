"""A simple GUI for displaying Slater electron densities.
   By Neal Coleman
"""

from easygui import *
import sys
import Slater
import numpy
import Routines
import Gnuplot
import Gnuplot.funcutils

while 1:
    msgbox("A simple GUI for displaying Slater electron densities\n" + "By Neal Coleman", "Introduction")

    Z = integerbox("Enter the atomic number of the element.", "Element Selection.")

    listlength = integerbox("Enter (integer) length of grid.", "Grid Length")

    choices1 = ["Exponential", "Uniform"]

    a = buttonbox("Would you like an exponential or uniform grid?", "Grid Type", choices1)

    if a == "Exponential":
        arrayx = Routines.ExpGridStretch2(numpy.arange(0.01, 1.0 * listlength, 0.01))
    else:
        arrayx = numpy.arange(0.01, 1.0 * listlength, 0.01)

    subshell = ["1s", "2s&p", "3s&p", "3d", "4s&p", "4d", "4f", "5s&p", "5d"]

    occ = []

    for i in range(9):
        occ.append(integerbox("Enter the occupancy of the "+subshell[i]+" subshell.", "Subshell Selection"))

    g = Gnuplot.Gnuplot()
    g.title("tmp")
    g('set data style linespoints')
    dty = 4 * numpy.pi * arrayx**2 * Slater.density(arrayx, Z, occ)[0]
    printout = []
    for i in range(len(arrayx)):
        printout.append([arrayx[i],dty[i]])
    # print printout
    g.plot(printout)
    # g.hardcopy('tmp.gif',enhanced=1,color=1)
    g.reset()

    if ccbox("There was your density.  Shall we do it again?", "Finale"):     # show a Continue/Cancel dialog
        pass  # user chose Continue
    else:
        sys.exit(0)           # user chose Cancel
