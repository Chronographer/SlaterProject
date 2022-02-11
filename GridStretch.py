# This is a set of routines to approximate first and second derivatives, integration, and error-checking routines.

import numpy


def ArrayDotProduct(vector1array, vector2array):
    """Takes two arrays of 3-vectors, calculates the 3-vector dot product\n
     # numpy.dot(vector1array, vector2array)  # This *should* do the same thing (?)"""
    x1, y1, z1 = vector1array[:, 0], vector1array[:, 1], vector1array[:, 2]
    x2, y2, z2 = vector2array[:, 0], vector2array[:, 1], vector2array[:, 2]
    return x1 * x2 + y1 * y2 + z1 * z2


def ArrayApproxDerivative(arrayX, arrayF):
    """Takes a numpy array of x values and an array of f values and finds an approximate centered first derivative."""
    arrayfprime = 0.0 * arrayF
    # first point indexed by 0
    arrayfprime[0] = (arrayF[1] - arrayF[0]) / (arrayX[1] - arrayX[0])
    # middle points done with slice
    arrayfprime[1:-1] = (arrayF[2:] - arrayF[0:-2]) / (arrayX[2:] - arrayX[0:-2])
    # last point indexed by -1
    arrayfprime[-1] = (arrayF[-1] - arrayF[-2]) / (arrayX[-1] - arrayX[-2])
    return arrayfprime


def ApproximateDerivative(listX, listF):
    """Takes a list of x values and a list of f values and finds an approximate centered first derivative."""
    N = len(listX)
    listfprime = [(listF[1] - listF[0]) / (listX[1] - listX[0])]
    for i in range(1, N - 1, 1):
        fprime = (listF[i + 1] - listF[i - 1]) / (listX[i + 1] - listX[i - 1])
        listfprime.append(fprime)
    listfprime.append((listF[-1] - listF[-2]) / (listX[-1] - listX[-2]))
    return listfprime


def CrudelyApproximateSecondDerivative(listX, listF):
    """Takes a list of x values and a list of f values and finds an approximate centered second derivative."""
    N = len(listX)
    listfdoubleprime = []
    fprimeprime = (listF[2] - listF[0]) / ((listX[2] - listX[0]) * (listX[1] - listX[0])) - (listF[1] - listF[0]) / (listX[1] - listX[0]) ** 2
    listfdoubleprime.append(fprimeprime)
    for i in range(1, N - 1, 1):
        fprimeprime = ((listF[i + 1] - listF[i]) / (listX[i + 1] - listX[i]) - (listF[i] - listF[i - 1]) / (listX[i] - listX[i - 1])) / ((listX[i + 1] - listX[i - 1]) / 2)
        listfdoubleprime.append(fprimeprime)
    fprimeprime = (listF[-1] - listF[-2]) / (listX[-1] - listX[-2]) ** 2 - (listF[-1] - listF[-3]) / ((listX[-1] - listX[-3]) * (listX[-1] - listX[-2]))
    listfdoubleprime.append(fprimeprime)
    return listfdoubleprime


def ArrayIntegrate(arrayX, arrayF):
    """Takes a array of f values and a grid of x values and integrates f with respect to f."""
    area = (arrayX[1:] - arrayX[0:-1]) * (arrayF[1:] + arrayF[0:-1]) * 0.5
    # add all values of area and reduce to scalar quantity
    return numpy.add.reduce(area)


def Integrate(listX, listF):
    """Takes a list of f values and a grid of x values and integrates f with respect to f."""
    rSum = 0
    for i in range(len(listX) - 1):
        rSum = rSum + (listX[i + 1] - listX[i]) * (listF[i + 1] + listF[i]) / 2
    return rSum


def RMSError(listApproximate, listExact):
    """Takes a list of approximate function values and a list of exact function values and returns a root-mean-square error."""
    nGrid = len(listApproximate)
    if len(listExact) == nGrid:
        sum = numpy.add.reduce((listApproximate - listExact) ** 2)
        chi = numpy.sqrt(sum / nGrid)
        return chi
    else:
        return "Make sure the lists have the same length."


def MARE(weight, listApproximate, listExact):
    """Returns a Mean Absolute Relative Error for two arrays."""
    numerator = numpy.add.reduce(weight * numpy.abs(listApproximate - listExact))
    denominator = numpy.add.reduce(weight)
    mare = numerator / denominator
    return mare


def radialMARE(grid, dens, listApproximate, listExact):
    """Returns a Mean Absolute Relative Error for two lists."""
    weight = grid ** 2 * dens
    return MARE(weight, listApproximate, listExact)


def PctgError(exact, approximate):
    """Returns percentage error at each point."""
    final = numpy.zeros(len(exact))
    final = numpy.abs(exact - approximate) / exact
    return final


def WeightedRMSError(weight, listApproximate, listExact):
    """Returns a weighted RMS error"""
    nGrid = len(listApproximate)
    if len(weight) == nGrid & len(listExact) == nGrid:
        sum = numpy.add.reduce(weight * (listApproximate - listExact) ** 2)
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


def ChainDerivativePQ(fp, fq, fpp, fpq, fqq, fqqq, gf, gff, gfff):
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


def ChainDerivativeQ(fq, fqq, fqqq, gf, gff, gfff):
    """Find derivatives of G(f(q)) in terms of those of f(q) and G(f)"""

    gq = gf * fq
    gqq = gff * fq * fq + gf * fqq
    gqqq = gfff * fq * fq * fq + 3.0 * gff * fqq * fq + gf * fqqq

    return gq, gqq, gqqq


""" def ChainDerivP(fp,fpp,gf,gff):
    # Find derivatives of G(f(q)) in terms of those of f(q) and G(f)

    gp = gf*fp
    gpp = gff*fp*fp + gf*fpp

    return gp, gpp, gpp """


# The following are the two exponential grid routines.


def ExpGridStretch1(N, minimumLength, a):
    """An exponential grid-creating routine.  Takes length in Bohr radii (N),
    minimum length (minimumLength), and factor (a).  Returns array of grid points."""
    stretchedListX = [minimumLength]
    i = minimumLength
    while i < N:
        i = i * a
        stretchedListX.append(i)
    return numpy.array(stretchedListX)


def ExpGridStretch2(listX):
    """Takes a (uniform) grid and returns a new grid, stretched exponentially,
    that has the same length (N), starting point (listX[0]), and ending point
    (listX[N-1])."""
    N = len(listX) - 1
    x0 = listX[0]
    xN = listX[N]
    stretchedListX = x0 * numpy.exp(numpy.arange(N + 1) * numpy.log(xN / x0) / N)
    return stretchedListX
