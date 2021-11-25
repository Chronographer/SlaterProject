from easygui import *
import numpy
import Routines

derivativeOptions = ["None", "First derivative", "Second derivative", "Third Derivative", "Fourth derivative"]

#  msgbox("A simple GUI for displaying Slater electron densities\n" + "By Neal Coleman", "Introduction")


def getElectronConfigInput():
    """Handles user input of electron shell population.
    \nTakes input as a string of numbers separated by periods, then parses and converts it to a list of integers.\n
    The returned list will always have 9 elements, and trailing zeros omitted by the user during input will
    be automatically added by this function.\n
    It also checks for and attempts to repair an input not ending in
    a number, warning the user of possible unintended consequences in the process."""
    orbitalConfigList = []
    while len(orbitalConfigList) != 9:
        orbitalConfigString = enterbox("Enter the occupation of each Slater shell group as a period separated list.\n(The Slater shells combine certain orbitals as follows: 1s, 2sp, 3sp, 3d, 4sp, 4d, 4f, 5sp, 5d, 6s)\n\nThere should be no more than 9 elements.\n\nTrailing zero's can be omitted and will be automatically added to the end of the list as needed.")
        if not orbitalConfigString[len(orbitalConfigString) - 1].isdigit():
            print("\nWARNING: the last character of your input string was '" + orbitalConfigString[len(orbitalConfigString) - 1] + "', which is not a digit. You might have accidentally left out a number or entered it incorrectly. The input you provided was: '" + orbitalConfigString + "'\nThe anomalous character has been removed in an attempt to prevent a runtime error.\nBe cautious; even if a runtime error is not thrown, your results might not reflect the system you intended to model!\n")
            orbitalConfigString = orbitalConfigString[:-1]  # this removes the last character from the input string; the idea is that if you accidentally end with a '.', this will remove it for you. Will produce incorrect results if the user left out a non-zero number after the last character.
        orbitalConfigList = orbitalConfigString.split(".")
        if len(orbitalConfigList) > 9:
            msgbox("There are " + str(len(orbitalConfigList)) + " elements in orbitalConfigList, there should be no more than 9.\nRe-enter the orbital configuration of each shell as a period separated list.")
        elif len(orbitalConfigList) < 9:
            while len(orbitalConfigList) < 9:  # Allows user to omit all trailing 0's by adding as many as are needed to generate a list of correct length, thus saving time.
                orbitalConfigList.append(0)
    for element in range(0, len(orbitalConfigList)):
        orbitalConfigList[element] = int(orbitalConfigList[element])
    return orbitalConfigList


def chooseDerivativeOptions():
    """Asks the user to choose between plotting the density or any of its first 4 derivatives."""
    derivativeType = buttonbox("Which derivative (if any) do you want to plot?", "Derivative Type", derivativeOptions)
    derivativeNumber = -1  # By initializing this to an unacceptable/nonsense value, if the code following this were to somehow break in a spectacular manner and not assign a value to derivativeNumber, it would prevent it from being uninitialized and make it easier to trace the problem back to here. I don't think its possible for this to occur, however. Think of this as looking both ways before crossing a one-way street.
    for i in range(len(derivativeOptions)):  # This lets you choose what derivative (if any) you want to plot and makes the plot title update automatically to reflect this.
        if derivativeType == derivativeOptions[i]:
            derivativeNumber = i
            break
    return derivativeNumber


def getArrayX(scaleType):
    """Takes a string value for scaleType and user input for the length of the x axis.\n
    Returns arrayX, appropriately scaled to either an exponential or uniform scale factor."""
    listLength = integerbox("Enter maximum radius for the figure (Bohr units).", "Grid Length")
    if scaleType == "Exponential":
        arrayX = Routines.ExpGridStretch2(numpy.arange(0.01, 1.0 * listLength, 0.01))
    else:
        arrayX = numpy.arange(0.01, 1.0 * listLength, 0.01)
    return arrayX


def getScaleType():
    """Gives the user the choice between scaling the x-axis of the plot in either an exponential or uniform manner.\n
    The returned string 'scaleType' is used in getArrayX()."""
    exponentialOrUniform = ["Exponential", "Uniform"]
    scaleType = buttonbox("Would you like an exponential or uniform grid?", "Grid Type", exponentialOrUniform)
    return scaleType


def getPlotType():
    """Asks the user if they want to plot the cumulative density of the
    system, (factoring in all shells and their screening), just the component
    density for each shell (not considering screening from other shells) or both."""
    plotTypeOptions = ["cumulative", "components", "both"]
    plotType = buttonbox("What would you like to plot?", "Plot Type", plotTypeOptions)
    return plotType


def getAtomicNumber():
    """Asks the user to enter the atomic number of the element they wish to model."""
    atomicNumber = integerbox("Enter the atomic number of the element.", "Element Selection.")
    return atomicNumber


def askToRepeat():
    """Asks the user if they would like to run the program again after the plot window has been closed."""
    reRunProgram = boolbox("There was your density. Would you like to run the program again?", "End of program")
    return reRunProgram
