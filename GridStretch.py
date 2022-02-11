import numpy


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
