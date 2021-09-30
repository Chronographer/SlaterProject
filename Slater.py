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

#List to extract from occupancy list the energy
energy = [1,2,3,3,4,4,4,5,5,6]

#The effective energy level dict
nstar = {1:1.0, 2:2.0, 3:3.0, 4:3.7, 5:4.0, 6:4.2}

#The shielding level dict
nshield = {"1s_valence":0.3, "Valence":0.35, "sp_semicore":0.85, "Core":1.0}
#No shielding, 1s2 toy!
#nshield = {"1s_valence":0.0, "Valence":1.0, "sp_semicore":1.0, "Core":1.0}

def s(listN):
    """Extract shielding constant from occupancy list - through 5d
    These are ansatz functions, so there's no general formula."""
    lists = []
    #1s
    lists.append((listN[0]-1)*nshield["1s_valence"])
    #2sp
    lists.append((listN[1]-1)*nshield["Valence"]+listN[0]*nshield["sp_semicore"])
    #3sp
    lists.append((listN[2]-1)*nshield["Valence"]+listN[1]*nshield["sp_semicore"]+listN[0]*nshield["Core"])
    #3d
    lists.append((listN[3]-1)*nshield["Valence"]+(listN[2]+listN[1]+listN[0])*nshield["Core"])
    #4sp
    lists.append((listN[4]-1)*nshield["Valence"]+(listN[3]+listN[2])*nshield["sp_semicore"]+(listN[1]+listN[0])*nshield["Core"])
    #4d
    lists.append((listN[5]-1)*nshield["Valence"]+(listN[4]+listN[3]+listN[2]+listN[1]+listN[0])*nshield["Core"])
    #4f
    lists.append((listN[6]-1)*nshield["Valence"]+(listN[5]+listN[4]+listN[3]+listN[2]+listN[1]+listN[0])*nshield["Core"])
    #5sp
    lists.append((listN[7]-1)*nshield["Valence"]+(listN[6]+listN[5]+listN[4])*nshield["sp_semicore"]+(listN[3]+listN[2]+listN[1]+listN[0])*nshield["Core"])
    #5d
    lists.append((listN[8]-1)*nshield["Valence"]+(listN[7]+listN[6]+listN[5]+listN[4]+listN[3]+listN[2]+listN[1]+listN[0])*nshield["Core"])
    #6s
    lists.append((listN[9]-1)*nshield["Valence"]+(listN[8]+listN[7]+listN[6]+listN[5]+listN[4]+listN[3]+listN[2]+listN[1]+listN[0])*nshield["Core"])
    return lists

def A(s,e,Z):
    """Normalization constant for ith shell"""
    nx=nstar[e]
    print "Energy quantum number, shielding, Z:  ", e, s, Z
    if (Z - s) < 0:
        print "Fatal Error: UNBOUND ATOM"
        print "Shielding charge {0} exceeds nuclear charge {1}".format(s,Z)
        print "EXITING SOON"
    a = numpy.sqrt( (2.0*(Z-s))**(2.0*nx+1) / \
        (4.0*numpy.pi*nx**(2.0*nx+1) * Gamma.gamma(2*nx+1.0)) )
    return a

def phi(s,e,Z,arrayx):
    """wavefunction for ith shell"""
    nx = nstar[e]
    arrayy = A(s,e,Z)*arrayx**(nx-1) * numpy.exp( -(Z-s)*arrayx/nx )
    return arrayy

def n(s,e,Z,arrayx,N):
    """density for ith shell; returns density and four derivatives,
    calculated recursively"""
    if ( N == 0 ):
        array0 = 0.0*arrayx
        return array0, array0, array0, array0, array0
    else:
        nx = nstar[e]
        mphi = numpy.absolute(phi(s,e,Z,arrayx))**2
        arrayy = N*mphi
        ayp = 2*( (nx - 1)/arrayx - (Z-s)/nx )*arrayy
        aypp = 2*(nx-1)*(-1/arrayx**2)*arrayy + 2*( (nx - 1)/arrayx - \
            (Z-s)/nx )*ayp
        ayp3 = 2*(nx-1)*(2/arrayx**3)*arrayy + 2*2*(nx-1)*(-1/arrayx**2)\
            *ayp + 2*( (nx - 1)/arrayx - (Z-s)/nx )*aypp
        ayp4 = 2*(nx-1)*(-6/arrayx**4)*arrayy + 3*2*(nx-1)*(2/arrayx**3)\
            *ayp + 3*2*(nx-1)*(-1/arrayx**2)*aypp + 2*( (nx - 1)/arrayx -\
            (Z-s)/nx )*ayp3
        return arrayy, ayp, aypp, ayp3, ayp4

def shelldensities(arrayx,Z,listN):
    """gives the densities of each shell"""
    svalues = s(listN)
    returnlist = []
    for j in range(len(listN)):
        sconst = svalues[j]
        e = energy[j]
        N = listN[j]
        #only want the density of each shell, and want list over shells
        returnlist.append(n(sconst,e,Z,arrayx,N)[0])
    return returnlist

def density(arrayx,Z,listN):
    """gets total density and its derivatives, summing over shells

    arrayx -- array of radial positions
    listN  -- list of shell occupancy numbers
    Z      -- net charge
    """
    svalues = s(listN)
    #print "svalues", svalues
    length = len(arrayx)
    #print "length", length
    final = numpy.zeros(length)
    finalp = numpy.zeros(length)
    finalpp = numpy.zeros(length)
    finalp3 = numpy.zeros(length)
    finalp4 = numpy.zeros(length)

    #print "density: listN", listN
    for j in range(len(listN)):
        #print "j, svalues[j]", j, svalues[j]
        sconst = svalues[j]
        e = energy[j]
        N = listN[j]
        if (N > 0):
            dens, densp, denspp, densp3, densp4 = n(sconst,e,Z,arrayx,N)
            final = final + dens
            finalp = finalp + densp
            finalpp = finalpp + denspp
            finalp3 = finalp3 + densp3
            finalp4 = finalp4 + densp4

    return final, finalp, finalpp, finalp3, finalp4

def grlaglll(arrayx, Z, listN):
    """Returns the RADIAL density and its gradient, laplacian, grad(lapl), lapl(lapl)
    arrayx -- array of radial positions
    listN  -- list of shell occupancy numbers
    Z      -- net charge
    """
    N = len(arrayx)

    #print "grlaglll: listN", listN
    d0, d1, d2, d3, d4 = density(arrayx, Z, listN)

    grad = d1
    lapl = d2 + 2 * d1 / arrayx
    glap = d3 + 2 * d2 / arrayx - 2 * d1 / arrayx**2
    llap = d4 + 4 * d3 / arrayx

    return d0, grad, lapl, glap, llap

"""The following are some miscellaneous analytic density routines."""

def H(listx):
    y = numpy.exp(-2*listx)/numpy.pi
    yprime = -2*numpy.exp(-2*listx)/numpy.pi
    return y, yprime

def Ne(listx):
    """N.B.: Not a real wavefunction!"""
    norm1s = numpy.sqrt(27.)
    y = numpy.exp(-listx)+norm1s*numpy.exp(-3*listx)
    yprime = -numpy.exp(-listx)-3.0*norm1s*numpy.exp(-3*listx)
    y2prime = numpy.exp(-listx)+3.0**2*norm1s*numpy.exp(-3*listx)
    y3prime = -numpy.exp(-listx)-3.0**3*norm1s*numpy.exp(-3*listx)
    y4prime = numpy.exp(-listx)+3.0**4*norm1s*numpy.exp(-3*listx)
    return y,yprime,y2prime,y3prime,y4prime
