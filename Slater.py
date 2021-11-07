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
import Routines

# List to extract from occupancy list the energy
energy = [1, 2, 3, 3, 4, 4, 4, 5, 5, 6]

# The effective energy level dict
nstar = {1: 1.0, 2: 2.0, 3: 3.0, 4: 3.7, 5: 4.0, 6: 4.2}

# The shielding level dict
nshield = {"1s_valence": 0.3, "Valence": 0.35, "sp_semicore": 0.85, "Core": 1.0}


# No shielding, 1s2 toy!
# nshield = {"1s_valence":0.0, "Valence":1.0, "sp_semicore":1.0, "Core":1.0}

def s(listN):
    """Extract shielding constant from occupancy list - through 5d
    These are ansatz functions, so there's no general formula."""
    lists = [(listN[0] - 1) * nshield["1s_valence"],
             (listN[1] - 1) * nshield["Valence"] + listN[0] * nshield["sp_semicore"],
             (listN[2] - 1) * nshield["Valence"] + listN[1] * nshield["sp_semicore"] + listN[0] * nshield["Core"],
             (listN[3] - 1) * nshield["Valence"] + (listN[2] + listN[1] + listN[0]) * nshield["Core"],
             (listN[4] - 1) * nshield["Valence"] + (listN[3] + listN[2]) * nshield["sp_semicore"] + (listN[1] + listN[0]) * nshield["Core"],
             (listN[5] - 1) * nshield["Valence"] + (listN[4] + listN[3] + listN[2] + listN[1] + listN[0]) * nshield["Core"],
             (listN[6] - 1) * nshield["Valence"] + (listN[5] + listN[4] + listN[3] + listN[2] + listN[1] + listN[0]) * nshield["Core"],
             (listN[7] - 1) * nshield["Valence"] + (listN[6] + listN[5] + listN[4]) * nshield["sp_semicore"] + (listN[3] + listN[2] + listN[1] + listN[0]) * nshield["Core"],
             (listN[8] - 1) * nshield["Valence"] + (listN[7] + listN[6] + listN[5] + listN[4] + listN[3] + listN[2] + listN[1] + listN[0]) * nshield["Core"]]#,
             #(listN[9] - 1) * nshield["Valence"] + (listN[8] + listN[7] + listN[6] + listN[5] + listN[4] + listN[3] + listN[2] + listN[1] + listN[0]) * nshield["Core"]]
    # 1s
    # 2sp
    # 3sp
    # 3d
    # 4sp
    # 4d
    # 4f
    # 5sp
    # 5d
    # 6s
    return lists


""" e -> energyQuantumNumber
    s -> shielding
    Z -> netCharge """


def A(shielding, energyQuantumNumber, netCharge):
    # Normalization constant for ith shell
    nx = nstar[energyQuantumNumber]
    print("Energy quantum number, shielding, netCharge:  ", energyQuantumNumber, shielding, netCharge)
    if (netCharge - shielding) < 0:
        print("Fatal Error: UNBOUND ATOM")
        print("Shielding charge {0} exceeds nuclear charge {1} ".format(shielding, netCharge))
        print("EXITING SOON")
    normalizationConstant = numpy.sqrt((2.0 * (netCharge - shielding)) ** (2.0 * nx + 1) / (4.0 * numpy.pi * nx ** (2.0 * nx + 1) * Gamma.gamma(2 * nx + 1.0)))
    return normalizationConstant


def phi(shielding, energyQuantumNumber, netCharge, arrayX):
    # wavefunction for ith shell
    nx = nstar[energyQuantumNumber]
    arrayY = A(shielding, energyQuantumNumber, netCharge) * arrayX ** (nx - 1) * numpy.exp(-(netCharge - shielding) * arrayX / nx)
    return arrayY


def n(shielding, energyQuantumNumber, netCharge, arrayX, N):
    # density for ith shell; returns density and four derivatives, calculated recursively
    if N == 0:
        array0 = 0.0 * arrayX
        return array0, array0, array0, array0, array0
    else:
        nx = nstar[energyQuantumNumber]
        mphi = numpy.absolute(phi(shielding, energyQuantumNumber, netCharge, arrayX)) ** 2
        arrayY = N * mphi
        ayp = 2 * ((nx - 1) / arrayX - (netCharge - shielding) / nx) * arrayY
        aypp = 2 * (nx - 1) * (-1 / arrayX ** 2) * arrayY + 2 * ((nx - 1) / arrayX - (netCharge - shielding) / nx) * ayp
        ayp3 = 2 * (nx - 1) * (2 / arrayX ** 3) * arrayY + 2 * 2 * (nx - 1) * (-1 / arrayX ** 2) * ayp + 2 * ((nx - 1) / arrayX - (netCharge - shielding) / nx) * aypp
        ayp4 = 2 * (nx - 1) * (-6 / arrayX ** 4) * arrayY + 3 * 2 * (nx - 1) * (2 / arrayX ** 3) * ayp + 3 * 2 * (nx - 1) * (-1 / arrayX ** 2) * aypp + 2 * ((nx - 1) / arrayX - (netCharge - shielding) / nx) * ayp3
        return arrayY, ayp, aypp, ayp3, ayp4


def shellDensities(arrayX, Z, listN):
    # gives the densities of each shell
    sValues = s(listN)
    returnList = []
    for j in range(len(listN)):
        sConstant = sValues[j]
        e = energy[j]
        N = listN[j]
        # only want the density of each shell, and want list over shells
        returnList.append(n(sConstant, e, Z, arrayX, N)[0])
    return returnList


def density(arrayX, netCharge, listN):
    """gets total density and its derivatives, summing over shells

    arrayX -- array of radial positions
    listN  -- list of shell occupancy numbers
    netCharge      -- net charge """
    sValues = s(listN)
    # print("sValues", sValues)
    length = len(arrayX)
    # print("length", length)
    final = numpy.zeros(length)
    finalp = numpy.zeros(length)
    finalpp = numpy.zeros(length)
    finalp3 = numpy.zeros(length)
    finalp4 = numpy.zeros(length)
    componentList = []
    returnList = []

    for j in range(len(listN)):
        sConstant = sValues[j]
        e = energy[j]
        N = listN[j]
        if N > 0:
            dens, densp, denspp, densp3, densp4 = n(sConstant, e, netCharge, arrayX, N)  # here
            final = final + dens
            finalp = finalp + densp
            finalpp = finalpp + denspp
            finalp3 = finalp3 + densp3
            finalp4 = finalp4 + densp4
            componentList.append(dens)
    returnList.append(final)
    return final, componentList


def grlaglll(arrayX, netCharge, listN):
    """Returns the RADIAL density and its gradient, laplacian, grad(lapl), lapl(lapl)
    arrayX -- array of radial positions
    listN  -- list of shell occupancy numbers
    netCharge      -- net charge """
    N = len(arrayX)

    # print("grlaglll: listN", listN)
    d0, d1, d2, d3, d4 = density(arrayX, netCharge, listN)

    grad = d1
    lapl = d2 + 2 * d1 / arrayX
    glap = d3 + 2 * d2 / arrayX - 2 * d1 / arrayX ** 2
    llap = d4 + 4 * d3 / arrayX

    return d0, grad, lapl, glap, llap


# The following are some miscellaneous analytic density routines.


def H(listX):
    y = numpy.exp(-2 * listX) / numpy.pi
    yprime = -2 * numpy.exp(-2 * listX) / numpy.pi
    return y, yprime


def Ne(listX):
    # N.B.: Not a real wavefunction!
    norm1s = numpy.sqrt(27.)
    y = numpy.exp(-listX) + norm1s * numpy.exp(-3 * listX)
    yprime = -numpy.exp(-listX) - 3.0 * norm1s * numpy.exp(-3 * listX)
    y2prime = numpy.exp(-listX) + 3.0 ** 2 * norm1s * numpy.exp(-3 * listX)
    y3prime = -numpy.exp(-listX) - 3.0 ** 3 * norm1s * numpy.exp(-3 * listX)
    y4prime = numpy.exp(-listX) + 3.0 ** 4 * norm1s * numpy.exp(-3 * listX)
    return y, yprime, y2prime, y3prime, y4prime
