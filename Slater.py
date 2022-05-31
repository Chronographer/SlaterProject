"""The following are the Slater atomic density approximation routines.

Implements charge densities and derivatives using the Slater model and
Slater-type orbitals, from

J. C. Slater, "Atomic Shielding Constants", Phys. Rev. 36, 57-64, (1930)

The working end of this is the charmingly named function "grlaglll" but there
are plenty more goodies to look at.

__author__ == "Neal Coleman"
__version__ == "Stable"
__documenter__ == "Antonio Cancio"
"""

import numpy
import Gamma

# List to extract from occupancy list the energy
energy = [1, 2, 3, 3, 4, 4, 4, 5, 5, 6]

# The effective energy level dict
nStar = {1: 1.0, 2: 2.0, 3: 3.0, 4: 3.7, 5: 4.0, 6: 4.2}

# The shielding level dict
nShield = {"1s_valence": 0.3, "Valence": 0.35, "sp_semicore": 0.85, "Core": 1.0}


# No shielding, 1s2 toy!
# nShield = {"1s_valence":0.0, "Valence":1.0, "sp_semicore":1.0, "Core":1.0}


def newShieldingConstantComputer(atom):
    shellLength = len(atom.occupancy)
    shellOccupancy = atom.occupancy
    principalQuantumNumber = atom.principalQuantumNumberLabelList
    angularMomentumLabel = atom.magneticQuantumNumberLabelList  # l[i] is the magnetic quantum number of an electron, ie 'sp' 'd' 'f' etc
    shielding = 0.0
    lists = []
    print("empty lists: " + str(lists))
    for i in range(0, shellLength):
        shielding = 0.0
        if i == 0:  # this is the 1s shell.
            shielding = shielding + 0.3 * (shellOccupancy[i]-1)
        else:
            for j in range(0, i+1):
                if j == i:
                    shielding = shielding + 0.35 * (shellOccupancy[j]-1)
                elif angularMomentumLabel[i] == 'sp' and principalQuantumNumber[j] == principalQuantumNumber[i]-1:
                    shielding = shielding + 0.85 * (shellOccupancy[j])
                else:
                    shielding = shielding + (shellOccupancy[j])
        lists.append(shielding)
    return lists  # this should be (Z effective) 


def ComputeShieldingConstants(electronConfigList):
    """Extract shielding constant from occupancy list - through 5d
    These are Ansatz functions, so there's no general formula."""
    shells = len(electronConfigList)
    lists = []
    for i in range(0, shells):  # this is really 2 nested loops
        if i == 0:  # 1s
            lists.append((electronConfigList[0] - 1) * nShield["1s_valence"])
        elif i == 1:  # 2sp
            lists.append((electronConfigList[1] - 1) * nShield["Valence"] + electronConfigList[0] * nShield["sp_semicore"])
        elif i == 2:  # 3sp
            lists.append((electronConfigList[2] - 1) * nShield["Valence"] + electronConfigList[1] * nShield["sp_semicore"] + electronConfigList[0] * nShield["Core"])
        elif i == 3:  # 3d
            lists.append((electronConfigList[3] - 1) * nShield["Valence"] + (electronConfigList[2] + electronConfigList[1] + electronConfigList[0]) * nShield["Core"])
        elif i == 4:  # 4sp
            lists.append((electronConfigList[4] - 1) * nShield["Valence"] + (electronConfigList[3] + electronConfigList[2]) * nShield["sp_semicore"] + (electronConfigList[1] + electronConfigList[0]) * nShield["Core"])
        elif i == 5:  # 4d
            lists.append((electronConfigList[5] - 1) * nShield["Valence"] + (electronConfigList[4] + electronConfigList[3] + electronConfigList[2] + electronConfigList[1] + electronConfigList[0]) * nShield["Core"])
        elif i == 6:  # 4f
            lists.append((electronConfigList[6] - 1) * nShield["Valence"] + (electronConfigList[5] + electronConfigList[4] + electronConfigList[3] + electronConfigList[2] + electronConfigList[1] + electronConfigList[0]) * nShield["Core"])
        elif i == 7:  # 5sp
            lists.append((electronConfigList[7] - 1) * nShield["Valence"] + (electronConfigList[6] + electronConfigList[5] + electronConfigList[4]) * nShield["sp_semicore"] + (electronConfigList[3] + electronConfigList[2] + electronConfigList[1] + electronConfigList[0]) * nShield["Core"])
        elif i == 8:  # 5d
            lists.append((electronConfigList[8] - 1) * nShield["Valence"] + (electronConfigList[7] + electronConfigList[6] + electronConfigList[5] + electronConfigList[4] + electronConfigList[3] + electronConfigList[2] + electronConfigList[1] + electronConfigList[0]) * nShield["Core"])
        elif i == 9:  # 6s
            lists.append((electronConfigList[9] - 1) * nShield["Valence"] + (electronConfigList[8] + electronConfigList[7] + electronConfigList[6] + electronConfigList[5] + electronConfigList[4] + electronConfigList[3] + electronConfigList[2] + electronConfigList[1] + electronConfigList[0]) * nShield["Core"])
        #elif i == 10:  # *** WARNING: Work in progress. THIS IS NOT CORRECT ***
         #   lists.append((electronConfigList[10] - 1) * nShield["Valence"] + (electronConfigList[9] + electronConfigList[9]))  # *** THIS IS NOT CORRECT ***
    return lists


""" e -> energyQuantumNumber  # This was put here to remind myself (Daniel) what the original variable names were, in case I missed something somewhere so I can be consistent with how I rename them.
    s -> shielding
    Z -> netCharge """


def A(shielding, energyQuantumNumber, netCharge):
    """Normalization constant for ith shell"""
    nx = nStar[energyQuantumNumber]
    # print("Energy quantum number, shielding, netCharge:  ", energyQuantumNumber, shielding, netCharge)
    if (netCharge - shielding) < 0:
        print("Fatal Error: UNBOUND ATOM")
        print("Shielding charge {0} exceeds nuclear charge {1} ".format(shielding, netCharge))
        print("EXITING SOON")
    normalizationConstant = numpy.sqrt((2.0 * (netCharge - shielding)) ** (2.0 * nx + 1) / (4.0 * numpy.pi * nx ** (2.0 * nx + 1) * Gamma.gamma(2 * nx + 1.0)))
    return normalizationConstant


def phi(shielding, energyQuantumNumber, netCharge, arrayX):
    """wavefunction for ith shell"""
    nx = nStar[energyQuantumNumber]
    arrayY = A(shielding, energyQuantumNumber, netCharge) * arrayX ** (nx - 1) * numpy.exp(-(netCharge - shielding) * arrayX / nx)
    return arrayY


def n(shielding, energyQuantumNumber, netCharge, arrayX, N):
    """density for ith shell; returns density and four derivatives, calculated recursively
    shielding -- shielding charge s
    energyQuantumNumber -- principal quantum number n
    netCharge -- nuclear charge Z
    arrayX -- grid of radial positions
    N -- occupancy number of shell
    Returns:  shell wavefunction squared * N, first, second, third,  fourth derivatives of nphi"""
    if N == 0:
        array0 = 0.0 * arrayX
        return array0, array0, array0, array0, array0
    else:
        nx = nStar[energyQuantumNumber]
        nphi = numpy.absolute(phi(shielding, energyQuantumNumber, netCharge, arrayX)) ** 2   # probability  density of orbital phi.
        arrayY = N * nphi
        ayp = 2 * ((nx - 1) / arrayX - (netCharge - shielding) / nx) * arrayY
        ayp2 = 2 * (nx - 1) * (-1 / arrayX ** 2) * arrayY + 2 * ((nx - 1) / arrayX - (netCharge - shielding) / nx) * ayp
        ayp3 = 2 * (nx - 1) * (2 / arrayX ** 3) * arrayY + 2 * 2 * (nx - 1) * (-1 / arrayX ** 2) * ayp + 2 * ((nx - 1) / arrayX - (netCharge - shielding) / nx) * ayp2
        ayp4 = 2 * (nx - 1) * (-6 / arrayX ** 4) * arrayY + 3 * 2 * (nx - 1) * (2 / arrayX ** 3) * ayp + 3 * 2 * (nx - 1) * (-1 / arrayX ** 2) * ayp2 + 2 * ((nx - 1) / arrayX - (netCharge - shielding) / nx) * ayp3
        return arrayY, ayp, ayp2, ayp3, ayp4


def shellDensities(arrayX, Z, listN):  # check to see if i made this
    """ gives the densities of each shell """
    shieldingValues = ComputeShieldingConstants(listN)
    returnList = []
    for j in range(len(listN)):
        shieldingConstant = shieldingValues[j]
        e = energy[j]
        N = listN[j]
        # only want the density of each shell, and want list over shells
        returnList.append(n(shieldingConstant, e, Z, arrayX, N)[0])
    return returnList


'''def newDensity(atom):
    """gets total density and its derivatives, summing over shells \n
    arrayX -- array of radial positions \n
    listN  -- list of shell occupancy numbers \n
    netCharge -- net charge"""
    shieldingValues = ComputeShieldingConstants(atom)
    length = len(arrayX)
    final = numpy.zeros(length)
    finalP1 = numpy.zeros(length)
    finalP2 = numpy.zeros(length)
    finalP3 = numpy.zeros(length)
    finalP4 = numpy.zeros(length)
    componentList = []
    for j in range(len(listN)):
        shieldingConstant = shieldingValues[j]
        e = energy[j]
        N = listN[j]
        if N > 0:
            dens, densP1, densP2, densP3, densP4 = n(shieldingConstant, e, netCharge, arrayX, N)
            final = final + dens
            finalP1 = finalP1 + densP1
            finalP2 = finalP2 + densP2
            finalP3 = finalP3 + densP3
            finalP4 = finalP4 + densP4
            densDerivativeList = [dens, densP1, densP2, densP3, densP4]
            componentList.append(densDerivativeList)
    finalDerivativeList = [final, finalP1, finalP2, finalP3, finalP4]
    return finalDerivativeList, componentList'''


def density(arrayX, netCharge, listN):
    """gets total density and its derivatives, summing over shells \n
    arrayX -- array of radial positions \n
    listN  -- list of shell occupancy numbers \n
    netCharge -- net charge"""
    shieldingValues = ComputeShieldingConstants(listN)
    print("output of old computeShieldingConstants is: " + str(shieldingValues))
    length = len(arrayX)
    final = numpy.zeros(length)
    finalP1 = numpy.zeros(length)
    finalP2 = numpy.zeros(length)
    finalP3 = numpy.zeros(length)
    finalP4 = numpy.zeros(length)
    componentList = []
    for j in range(len(listN)):
        shieldingConstant = shieldingValues[j]
        e = energy[j]
        N = listN[j]
        if N > 0:
            dens, densP1, densP2, densP3, densP4 = n(shieldingConstant, e, netCharge, arrayX, N)
            final = final + dens
            finalP1 = finalP1 + densP1
            finalP2 = finalP2 + densP2
            finalP3 = finalP3 + densP3
            finalP4 = finalP4 + densP4
            densDerivativeList = [dens, densP1, densP2, densP3, densP4]
            componentList.append(densDerivativeList)
    finalDerivativeList = [final, finalP1, finalP2, finalP3, finalP4]
    return finalDerivativeList, componentList


def grlaglll(arrayX, netCharge, listN):
    """Returns the RADIAL density and its gradient, laplacian, grad(lapl), lapl(lapl)\n
    arrayX -- array of radial positions \n
    listN  -- list of shell occupancy numbers \n
    netCharge      -- net charge"""
    N = len(arrayX)

    # print("grlaglll: listN", listN)
    d0, d1, d2, d3, d4 = density(arrayX, netCharge, listN)[0]

    grad = d1
    lapl = d2 + 2 * d1 / arrayX
    glap = d3 + 2 * d2 / arrayX - 2 * d1 / arrayX ** 2
    llap = d4 + 4 * d3 / arrayX

    return d0, grad, lapl, glap, llap


# The following are some miscellaneous analytic density routines.


def H(listX):
    """Most likely for debugging purposes"""
    y = numpy.exp(-2 * listX) / numpy.pi
    yPrime = -2 * numpy.exp(-2 * listX) / numpy.pi
    return y, yPrime


def Ne(listX):
    """Most likely for debugging purposes"""
    # N.B.: Not a real wavefunction!
    norm1s = numpy.sqrt(27.)
    y = numpy.exp(-listX) + norm1s * numpy.exp(-3 * listX)
    yPrime = -numpy.exp(-listX) - 3.0 * norm1s * numpy.exp(-3 * listX)
    y2Prime = numpy.exp(-listX) + 3.0 ** 2 * norm1s * numpy.exp(-3 * listX)
    y3Prime = -numpy.exp(-listX) - 3.0 ** 3 * norm1s * numpy.exp(-3 * listX)
    y4Prime = numpy.exp(-listX) + 3.0 ** 4 * norm1s * numpy.exp(-3 * listX)
    return y, yPrime, y2Prime, y3Prime, y4Prime
