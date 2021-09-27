"""This is a set of routines to approximate first and second derivatives,
integration, and error-checking routines."""

import numpy


def ArrayDotProduct(vector1array, vector2array):
    """Takes two arrays of 3-vectors, calculates the 3-vector dot product"""
    x1, y1, z1 = vector1array[:, 0], vector1array[:, 1], vector1array[:, 2]
    x2, y2, z2 = vector2array[:, 0], vector2array[:, 1], vector2array[:, 2]
    return x1 * x2 + y1 * y2 + z1 * z2


def ArrayApproxDerivative(arrayx, arrayf):
    """Takes a numpy arrray of x values and an array of f values and finds an
    approximate centered first derivative."""
    arrayfprime = 0.0 * arrayf
    # first point indexed by 0
    arrayfprime[0] = (arrayf[1] - arrayf[0]) / (arrayx[1] - arrayx[0])
    # middle points done with slice
    arrayfprime[1:-1] = (arrayf[2:] - arrayf[0:-2]) / (arrayx[2:] - arrayx[0:-2])
    # last point indexed by -1
    arrayfprime[-1] = (arrayf[-1] - arrayf[-2]) / (arrayx[-1] - arrayx[-2])
    return arrayfprime


def ApproximateDerivative(listx, listf):
    """Takes a list of x values and a list of f values and finds an
    approximate centered first derivative."""
    N = len(listx)
    listfprime = []
    listfprime.append((listf[1] - listf[0]) / (listx[1] - listx[0]))
    for i in range(1, N - 1, 1):
        fprime = (listf[i + 1] - listf[i - 1]) / (listx[i + 1] - listx[i - 1])
        listfprime.append(fprime)
    listfprime.append((listf[-1] - listf[-2]) / (listx[-1] - listx[-2]))
    return listfprime


def CrudelyApproximateSecondDerivative(listx, listf):
    """Takes a list of x values and a list of f values and finds an
    approximate centered second derivative."""
    N = len(listx)
    listfdoubleprime = []
    fprimeprime = (listf[2] - listf[0]) / ((listx[2] - listx[0]) * \
                                           (listx[1] - listx[0])) - (listf[1] - listf[0]) / (listx[1] - listx[0]) ** 2
    listfdoubleprime.append(fprimeprime)
    for i in range(1, N - 1, 1):
        fprimeprime = ((listf[i + 1] - listf[i]) / (listx[i + 1] - listx[i]) \
                       - (listf[i] - listf[i - 1]) / (listx[i] - listx[i - 1])) / ((listx[i + 1] - listx[i - 1]) / 2)
        listfdoubleprime.append(fprimeprime)
    fprimeprime = (listf[-1] - listf[-2]) / (listx[-1] - listx[-2]) ** 2 \
                  - (listf[-1] - listf[-3]) / ((listx[-1] - listx[-3]) * (listx[-1] - listx[-2]))
    listfdoubleprime.append(fprimeprime)
    return listfdoubleprime


def ArrayIntegrate(arrayx, arrayf):
    """Takes a array of f values and a grid of x values and integrates
    f with respect to f."""
    area = (arrayx[1:] - arrayx[0:-1]) * (arrayf[1:] + arrayf[0:-1]) * 0.5
    # add all values of area and reduce to scalar quantity
    return numpy.add.reduce(area)


def Integrate(listx, listf):
    """Takes a list of f values and a grid of x values and integrates
    f with respect to f."""
    rsum = 0
    for i in range(len(listx) - 1):
        rsum = rsum + (listx[i + 1] - listx[i]) * (listf[i + 1] + listf[i]) / 2
    return rsum


def RMSError(listapprox, listexact):
    """Takes a list of approximate function values and a list of exact
    function values and returns a root-mean-square error."""
    ngrid = len(listapprox)
    if len(listexact) == ngrid:
        sum = numpy.add.reduce((listapprox - listexact) ** 2)
        chi = numpy.sqrt(sum / ngrid)
        return chi
    else:
        return "Make sure the lists have the same length."


def MARE(weight, listapprox, listexact):
    """Returns a mean absolute relative error for two arrays."""
    numerator = numpy.add.reduce(weight * numpy.abs(listapprox - listexact))
    denominator = numpy.add.reduce(weight)
    mare = numerator / denominator
    return mare


def radialMARE(grid, dens, listapprox, listexact):
    """Returns a mean absolute relative error for two lists."""
    weight = grid ** 2 * dens
    return MARE(weight, listapprox, listexact)


def PctgError(exact, approx):
    """Returns percentage error at each point."""
    final = numpy.zeros(len(exact))
    final = numpy.abs(exact - approx) / exact
    return final


def WgtRMSError(weight, listapprox, listexact):
    """Returns a weighted RMS error"""
    ngrid = len(listapprox)
    if len(weight) == ngrid & len(listexact) == ngrid:
        sum = numpy.add.reduce(weight * (listapprox - listexact) ** 2)
        chi = numpy.sqrt(sum / numpy.add.reduce(weight))
        return chi
    else:
        return "Make sure all lists have same length."


def FurthestIndex(grid, array, a):
    """Finds the index and concurrent grid value beyond which all points have magnitude less than a."""
    N = len(grid)
    i = N - 1
    while (array[i] > a) and i > -1:
        i = i - 1
    return i - 1, grid[i - 1]


def ChainDerivPQ(fp, fq, fpp, fpq, fqq, fqqq, gf, gff, gfff):
    """Find derivatives of G(f(p,q)) in terms of those of f(p,q) and G(f)"""

    gp = gf * fp
    gq = gf * fq
    gpp = gff * fp * fp + gf * fpp
    gpq = gff * fp * fq + gf * fpq
    gqq = gff * fq * fq + gf * fqq
    gqqq = gfff * fq * fq * fq + 3.0 * gff * fqq * fq + gf * fqqq
    # unused in Laplacian DFT
    # gppp = gfff*fp*fp*fp + 3.0*gff*fpp*fp + gf*fppp
    # gppq = gfff*fp*fp*fq + gff*(2.0*fpp*fq + fpq*fp) + gf*fppq
    # gpqq = gfff*fp*fq*fq + gff*(2.0*fqq*fp + fpq*fq) + gf*fpqq

    return gp, gq, gpp, gpq, gqq, gqqq


def ChainDerivQ(fq, fqq, fqqq, gf, gff, gfff):
    """Find derivatives of G(f(q)) in terms of those of f(q) and G(f)"""

    gq = gf * fq
    gqq = gff * fq * fq + gf * fqq
    gqqq = gfff * fq * fq * fq + 3.0 * gff * fqq * fq + gf * fqqq

    return gq, gqq, gqqq


# def ChainDerivP(fp,fpp,gf,gff):
#    """Find derivatives of G(f(q)) in terms of those of f(q) and G(f)"""
#
#    gp = gf*fp
#    gpp = gff*fp*fp + gf*fpp
#
#    return gp, gpp, gpp

"""The following are the two exponential grid routines."""


def ExpGridStretch1(N, x0, a):
    """An exponential grid-creating routine.  Takes length in Bohr radii (N),
    minimum length (x0), and factor (a).  Returns array of grid points."""
    stretchedlistx = [x0]
    i = x0
    while i < N:
        i = i * a
        stretchedlistx.append(i)
    return numpy.array(stretchedlistx)


def ExpGridStretch2(listx):
    """Takes a (uniform) grid and returns a new grid, stretched exponentially,
    that has the same length (N), starting point (listx[0]), and ending point
    (listx[N-1])."""
    N = len(listx) - 1
    x0 = listx[0]
    xN = listx[N]
    stretchedlistx = x0 * numpy.exp(numpy.arange(N + 1) * numpy.log(xN / x0) / N)
    return stretchedlistx
