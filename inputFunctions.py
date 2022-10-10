from easygui import *
import numpy
import GridStretch
import Atoms

derivativeOptions = ["None", "First derivative", "Second derivative", "Third Derivative", "Fourth derivative"]

#  msgbox("A simple GUI for displaying Slater electron densities\n" + "By Neal Coleman", "Introduction")


def getElementNameInput():
    """Handles user input of desired element to model.\n
    Takes a string as input which should either be 'manual' or a case sensitive element abbreviation.\n
    Returns a string containing either a valid element abbreviation or the string 'manual'.\n
    It also handles invalid selections by clearing the user input and trying again if the input cannot be matched to an element from the 'periodictable{}' dictionary in Atoms.py"""
    elementInput = enterbox("Enter the abbreviation for the element you wish to model, or type 'manual' to manually specify the atomic number and Slater electron configuration.", "Element Selection (automatic)")
    elementCodeIsValid = False
    while not elementCodeIsValid:
        if elementInput == "manual":
            break  # if you want manual control, we do not need to waste time checking every key in the dictionary.
        for key in Atoms.atomData:  # We want to make sure that there is a dictionary key matching the input before we pass it back to SlaterGUI.py
            if key == elementInput:
                elementCodeIsValid = True
                break  # once we find a key matching the user input, we do not need to check the rest.
        if not elementCodeIsValid:
            msgbox("'" + elementInput + "' is not a recognized element!\n\nElement abbreviations are case sensitive. Check to make sure you spelled it correctly.\n\nNot all elements are in the database. For a list of all recognized elements, see dictionary 'periodictable{}' in script 'Atoms.py'", "Input Error: Unrecognized Element")
            elementInput = enterbox("Enter the abbreviation for the element you wish to model, or type 'manual' to manually specify the atomic number and Slater electron configuration.", "Element Selection (automatic)")

    return elementInput


def getElectronConfigInput():
    """Handles user input of electron shell population.
    \nTakes input as a string of numbers separated by periods, then parses and converts it to a list of integers.\n
    It also checks for and attempts to automatically repair an input not ending in
    a number, warning the user of possible unintended consequences in the process.\n
    Additionally, it handles empty string inputs and inputs which contain characters which are neither numbers or
    periods by clearing the user input and asking them to try again."""
    orbitalConfigList = []
    while len(orbitalConfigList) == 0:
        orbitalConfigString = enterbox("Enter the occupation of each Slater orbital as a period separated list. (Ie, combine the s and p orbitals)\n\nThe general input format should thus be as follows:\n1s, 2sp, 3sp, 3d, 4sp, 4d, 4f, 5sp, 5d, 5f, 5g, 6sp, 6d, 6f, 6g, 6h, 7sp, etc.\n\nYour input MUST include zeros for empty orbitals which are inside of populated orbitals. Trailing zeros may be omitted.\n\nFor example, Potassium (19) populates the 4sp orbital before 3d, and thus it should be formatted as '2.8.8.0.1'.", "Electron Configuration Entry")
        if orbitalConfigString == "":
            msgbox("There are no characters in 'orbitalConfigString'. (You provided an empty string)\n\nRe-enter the orbital configuration of each shell as a period separated list.", "Input Error: Empty string")
        elif not orbitalConfigString[len(orbitalConfigString) - 1].isdigit():
            print("\nWARNING: The last character of your input string was '" + orbitalConfigString[len(orbitalConfigString) - 1] + "', which is not a digit. You might have accidentally left out a number or entered it incorrectly. The input you provided was: '" + orbitalConfigString + "'\nThe anomalous character has been removed in an attempt to prevent a runtime error.\nBe cautious; even if a runtime error is not thrown, your results might not reflect the system you intended to model!\n")
            orbitalConfigString = orbitalConfigString[:-1]  # this removes the last character from the input string; the idea is that if you accidentally end with a '.', this will remove it for you. Will produce incorrect results if the user left out a non-zero number after the last character.
        orbitalConfigList = orbitalConfigString.split(".")
        for i in range(0, len(orbitalConfigList)):
            if not orbitalConfigList[i].isdigit():
                if not orbitalConfigList[i] == "":  # we already gave the user a message about providing an empty string above, so we don't want to give them a different message here for the same problem. We DO however still want to clear orbitalConfigList so it won't throw an error further on.
                    msgbox("There are non-number elements in 'orbitalConfigList' after parsing to remove the periods. It should (at this point) contain only numbers (stored as strings). You most likely hit an errant key by accident.\n\nCurrently, 'orbitalConfigList' is " + str(orbitalConfigList) + ". The non-number element is '" + str(orbitalConfigList[i]) + "'.\n\nRe-enter the orbital configuration of each shell as a period separated list.", "Parsing Error: Non-Number Element")
                elif len(orbitalConfigList) > 1 and orbitalConfigList[i] == "":
                    msgbox("There is at least one empty string in 'orbitalConfigList' after parsing to remove the periods. It should (at this point) contain only numbers (stored as strings). You mostly likely either hit the period key twice or omitted a number in between periods.\n\nCurrently, 'orbitalConfigList' is " + str(orbitalConfigList) + ".\n\nRe-enter the orbital configuration of each shell as a period separated list.", "Parsing Error: Empty String")
                orbitalConfigList.clear()
                break  # clear the list so we can start over, then break out of this loop so we don't try to check the next element in it, which no longer exists.

    for element in range(0, len(orbitalConfigList)):  # This converts the list of numbers stored as strings into a list of integers.
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
        arrayX = GridStretch.ExpGridStretch2(numpy.arange(0.01, 1.0 * listLength, 0.01))
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
    density for each shell (just showing the contributions from each individual shell) or both."""
    plotTypeOptions = ["cumulative", "components", "both"]
    plotType = buttonbox("What would you like to plot?", "Plot Type", plotTypeOptions)
    return plotType


def getAtomicNumber():
    """Asks the user to enter the atomic number of the element they wish to model."""
    atomicNumber = integerbox("Enter the atomic number of the element.", "Element Selection (manual)")
    return atomicNumber


def showEnergy(energy):
    msgbox("The total energy is " + str(energy) + " Hartrees.", "Total Energy")


def askToRepeat():
    """Asks the user if they would like to run the program again after the plot window has been closed."""
    reRunProgram = boolbox("There was your density. Would you like to run the program again?", "End of program")
    return reRunProgram


def getElectronConfigFromJson(orbitalConfigString):
    """Converts a string of shell occupancies from a JSON file into something usable by Slater.py. This is ONLY used by the JSON controlled Slater program.
        \nTakes input as a string of numbers separated by periods, then parses and converts it to a list of integers.\n
        It also checks for and attempts to repair an input not ending in
        a number, warning the user of possible unintended consequences in the process."""

    if not orbitalConfigString[len(orbitalConfigString) - 1].isdigit():
        print("\nWARNING: the last character of your input string was '" + orbitalConfigString[len(orbitalConfigString) - 1] + "', which is not a digit. You might have accidentally left out a number or entered it incorrectly. The input you provided was: '" + orbitalConfigString + "'\nThe anomalous character has been removed in an attempt to prevent a runtime error.\nBe cautious; even if a runtime error is not thrown, your results might not reflect the system you intended to model!\n")
        orbitalConfigString = orbitalConfigString[:-1]  # this removes the last character from the input string; the idea is that if you accidentally end with a '.', this will remove it for you. Will produce incorrect results if the user left out a non-zero number after the last character.
    orbitalConfigList = orbitalConfigString.split(".")
    if len(orbitalConfigList) > 10:
        msgbox("There are " + str(len(orbitalConfigList)) + " elements in orbitalConfigList, there should be no more than 10.\nRe-enter the orbital configuration of each shell as a period separated list.\n\nThis problem almost certainly started in your input .json file.")
    for element in range(0, len(orbitalConfigList)):
        orbitalConfigList[element] = int(orbitalConfigList[element])
    return orbitalConfigList


def getArrayXFromJSON(scaleType, listLength):
    """Takes a string value for scaleType and an integer extracted from a JSON file for the length of the x axis.\n
        Returns arrayX, appropriately scaled to either an exponential or uniform scale factor."""
    if scaleType == "Exponential":
        arrayX = GridStretch.ExpGridStretch2(numpy.arange(0.01, 1.0 * listLength, 0.01))
    else:
        arrayX = numpy.arange(0.01, 1.0 * listLength, 0.01)
    return arrayX
